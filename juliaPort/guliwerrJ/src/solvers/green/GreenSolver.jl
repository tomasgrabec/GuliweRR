# =========================================================================================
# GreenSolver.jl -- boundary-condition assembly, the velocity loop, the +-k assembly and the
# optional full output of `guliwerr/solvers/green/@GreenSolver/calc.m`.
#
# As in MATLAB's calc.m:
#   * ONE solve per point: only k1 = +1 is solved; G_IJ(-k) = G_JI(+k)
#     by elastodynamic reciprocity, so GFull_minus is GFull transposed in (I, J) and
#     `computeMinus` only decides whether the minus/eff fields are filled (it costs nothing);
#   * G_eff is the parity-signed mean (G_+ + (-1)^n G_-)/2, n = number of (I, J) equal to 1
#     -- e.g. G13_eff = (G13_+ - G31_+)/2, G33_eff = G33_+;
#   * `eta` keyword (calc.m 'Eta'), default LIMITING_ABSORPTION_ETA = 3e-4;
#   * a non-finite B counts as singular (rcond NaN -> 0).
#
# PERFORMANCE DESIGN (Julia-side; results stay one-to-one with MATLAB):
#   * a per-point kernel (`_green_point!`) on preallocated chunk buffers, no allocation in
#     the velocity loop: Aberth roots warm-started from the previous point, null vectors by
#     cross products, only the substrate waves that are actually used get a polarisation,
#     B assembled in place, one partial-pivoting LU per point giving the inverse and an
#     exact 1-norm rcond;
#   * velocity chunks of `GREEN_CHUNK` points, each cold-started (companion-matrix roots at
#     its first point), so the result is bitwise independent of the thread count;
#   * THREADING POLICY: `calc` spreads the chunks over Julia's task scheduler
#     (`threaded = true`, default). Julia's scheduler is composable, so a caller that
#     already threads over directions (SpotGreen, sweeps) may keep the default -- nested
#     tasks are queued, never OS-oversubscribed; pass `threaded = false` there to save the
#     (small) task overhead on short axes.
# =========================================================================================

using LinearAlgebra

using ...CoreTypes: Structure, nlay
import ..Ritz: calc          # ADD a method to Ritz's generic `calc`; never a new binding

# -----------------------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------------------

"""Depth-axis resolution of the half-space `fullOutput` displacement profile (`calc.m`)."""
const DISP_FIELD_NDEPTH = 6003

"""Depth-axis extent of the half-space `fullOutput` profile, normalised units (30 wavelengths)."""
const DISP_FIELD_MAX_DEPTH = 2 * pi * 30

"""
Velocity points per kernel chunk. Each chunk cold-starts its root iteration, so results do
not depend on how chunks are distributed over threads; 64 keeps the cold-start share
(one companion-matrix eigensolve per layer per chunk) small.
"""
const GREEN_CHUNK = 64

# -----------------------------------------------------------------------------------------
# Small dense linear algebra, allocation-free
# -----------------------------------------------------------------------------------------

@inline _cabs1(z::Complex) = abs(real(z)) + abs(imag(z))

"""
    _lu_inv!(X, A, ipiv, n) -> Bool

Overwrite `A[1:n, 1:n]` by its LU factors (partial pivoting on `|re| + |im|`, LAPACK's
`izamax` rule) and write `inv(A)` into `X[1:n, 1:n]`. A zero pivot (exactly singular)
gives an all-`Inf` `X`, MATLAB's `inv` behaviour, and returns `false`.
"""
function _lu_inv!(X::AbstractMatrix{ComplexF64}, A::AbstractMatrix{ComplexF64},
                  ipiv::AbstractVector{Int}, n::Int)
    @inbounds begin
        for k in 1:n
            p = k
            amax = _cabs1(A[k, k])
            for i in (k + 1):n
                a = _cabs1(A[i, k])
                if a > amax
                    amax = a
                    p = i
                end
            end
            ipiv[k] = p
            if iszero(A[p, k])
                for j in 1:n, i in 1:n
                    X[i, j] = ComplexF64(Inf)
                end
                return false
            end
            if p != k
                for j in 1:n
                    A[k, j], A[p, j] = A[p, j], A[k, j]
                end
            end
            piv = inv(A[k, k])
            for i in (k + 1):n
                A[i, k] *= piv
            end
            for j in (k + 1):n
                akj = A[k, j]
                iszero(akj) && continue
                for i in (k + 1):n
                    A[i, j] -= A[i, k] * akj
                end
            end
        end
        # X = P * I, then forward (unit L) and back (U) substitution, column by column
        for j in 1:n, i in 1:n
            X[i, j] = i == j ? one(ComplexF64) : zero(ComplexF64)
        end
        for k in 1:n
            p = ipiv[k]
            if p != k
                for j in 1:n
                    X[k, j], X[p, j] = X[p, j], X[k, j]
                end
            end
        end
        for j in 1:n
            for k in 1:n
                xk = X[k, j]
                iszero(xk) && continue
                for i in (k + 1):n
                    X[i, j] -= A[i, k] * xk
                end
            end
            for k in n:-1:1
                X[k, j] /= A[k, k]
                xk = X[k, j]
                for i in 1:(k - 1)
                    X[i, j] -= A[i, k] * xk
                end
            end
        end
    end
    return true
end

"1-norm (max absolute column sum) of `A[1:n, 1:n]`."
function _norm1(A::AbstractMatrix, n::Int)
    m = 0.0
    @inbounds for j in 1:n
        s = 0.0
        for i in 1:n
            s += abs(A[i, j])
        end
        m = max(m, s)
    end
    return m
end

"""
    _rcond_from_inv(B, X, n) -> Float64

`1 / (norm(B, 1) * norm(inv(B), 1))`, the exact reciprocal 1-norm condition number.
MATLAB's `rcond` is LAPACK's `zgecon` *estimate* of the same quantity (equal on small
well-behaved matrices; the estimate can only be larger). A non-finite B or inverse gives
`0` -- `calc.m`'s `if isnan(rc), rc = 0`.
"""
function _rcond_from_inv(B::AbstractMatrix, X::AbstractMatrix, n::Int)
    a = _norm1(B, n)
    x = _norm1(X, n)
    (isfinite(a) && isfinite(x)) || return 0.0
    (a == 0 || x == 0) && return 0.0
    return inv(a * x)
end

"""
    _pow2_reciprocal(x) -> Float64

The power of two nearest below `1/x` (`1.0` for a non-finite or non-positive `x`), clamped
against over/underflow. Multiplying by a power of two is exact in IEEE arithmetic.
"""
function _pow2_reciprocal(x::Real)
    (isfinite(x) && x > 0) || return 1.0
    _, e = frexp(Float64(x))
    return ldexp(1.0, -clamp(e, -1000, 1000))
end

"""
    _equilibrate_pow2!(r, c, B, n)

Row then column power-of-two scale vectors making `Diagonal(r) * B * Diagonal(c)` have
every row/column infinity-norm within a factor of two of 1.

**Deliberate divergence from MATLAB:** `calc.m` uses
`equilibrate(B)` (MC64 permutation + scaling), which Julia does not have. Both give the
identical true inverse (`invB = C*inv(R*P*B*C)*R*P` for any invertible `P, R, C`) and the
branch only runs where `rcond(B) < nWaves*eps`, i.e. on points `calc` already reports as
near-singular; measured on a relDepth-160 case, G agrees to 1.5e-16 relative. The
only observable difference is whether such a point is counted in the warning.
"""
function _equilibrate_pow2!(r::AbstractVector{Float64}, c::AbstractVector{Float64},
                            B::AbstractMatrix, n::Int)
    @inbounds for i in 1:n
        mx = 0.0
        for j in 1:n
            mx = max(mx, abs(B[i, j]))
        end
        r[i] = _pow2_reciprocal(mx)
    end
    @inbounds for j in 1:n
        mx = 0.0
        for i in 1:n
            mx = max(mx, r[i] * abs(B[i, j]))
        end
        c[j] = _pow2_reciprocal(mx)
    end
    return r, c
end

"Allocating form: `(r, c)`."
function _equilibrate_pow2(B::AbstractMatrix)
    n = size(B, 1)
    r = ones(n); c = ones(n)
    _equilibrate_pow2!(r, c, B, n)
    return r, c
end

# -----------------------------------------------------------------------------------------
# Reference (allocating) forms of calc.m's local functions -- kept for diagnostics and MATLAB comparison
# -----------------------------------------------------------------------------------------

"""
    calc_traction_cols(C3, U, k1, k2, k3, omega) -> Matrix{ComplexF64}      # 3 x nW

`T[l, n] = sum_{p,q} C_{3,l,p,q} U[p, n] k_q[n] / omega` (`calc.m` `calcTractionCols`),
`C3 = reshape(elTen[3, :, :, :], 3, 9)`.
"""
function calc_traction_cols(C3::AbstractMatrix{<:Real}, U::AbstractMatrix{<:Complex},
                             k1::Real, k2::Real, k3::AbstractVector{<:Complex},
                             omega::Complex)
    S = kvec_cols(k1, k2, k3) ./ omega
    return C3 * outer_cols(U, S)
end

"""
    calc_gfull_from_partial_waves(invB, U, phase, omega) -> Matrix{ComplexF64}   # 3 x 3

`G[I, J] = (i/omega) * sum_n invB[n, I] * U[J, n] * phase[n]` (`calc.m`
`calcGFullFromPartialWaves`): force index `I` reads invB's COLUMN, displacement `J` U's ROW.
"""
function calc_gfull_from_partial_waves(invB::AbstractMatrix{<:Complex},
                                        U::AbstractMatrix{<:Complex},
                                        phase, omega::Complex)
    nW = size(U, 2)
    scalarPhase = phase isa Number
    scalarPhase || length(phase) == nW || throw(DimensionMismatch(
        "calc_gfull_from_partial_waves: phase must be scalar or have $nW entries"))
    G = Matrix{ComplexF64}(undef, 3, 3)
    fac = im / omega
    @inbounds for J in 1:3, I in 1:3
        s = zero(ComplexF64)
        for n in 1:nW
            ph = scalarPhase ? phase : phase[n]
            s += (invB[n, I] * U[J, n]) * ph
        end
        G[I, J] = fac * s
    end
    return G
end

"""
    build_bc_matrix(C3lay, Ulay, k3lay, k1, k2, omega, nFilm, cumDepthFromBelow, hScale)
        -> (B, surfPhase)

`(6*nFilm + 3)^2` boundary-condition matrix (Every 2002; `calc.m` `buildBCMatrix`), internal
(surface-first) layer order; `surfPhase` = film 1's free-surface phase factors.
"""
function build_bc_matrix(C3lay::AbstractVector{<:AbstractMatrix},
                          Ulay::AbstractVector{<:AbstractMatrix},
                          k3lay::AbstractVector{<:AbstractVector},
                          k1::Real, k2::Real, omega::Complex, nFilm::Integer,
                          cumDepthFromBelow::AbstractVector{<:Real}, hScale::Real)
    nWaves = 6 * nFilm + 3
    B = zeros(ComplexF64, nWaves, nWaves)
    Tlay = [calc_traction_cols(C3lay[i], Ulay[i], k1, k2, k3lay[i], omega) for i in 1:(nFilm + 1)]
    hEff1 = cumDepthFromBelow[1] * hScale
    surfPhase = exp.((-im * hEff1) .* k3lay[1])
    @inbounds for c in 1:6, r in 1:3
        B[r, c] = Tlay[1][r, c] * surfPhase[c]
    end
    @inbounds for iLay in 1:nFilm
        rD = 3 + 6 * (iLay - 1) + 1
        rT = rD + 3
        dIface = iLay < nFilm ? cumDepthFromBelow[iLay + 1] * hScale : 0.0
        cS = 6 * (iLay - 1) + 1
        ph = exp.((-im * dIface) .* k3lay[iLay])
        for c in 1:6, r in 1:3
            B[rD + r - 1, cS + c - 1] = Ulay[iLay][r, c] * ph[c]
            B[rT + r - 1, cS + c - 1] = Tlay[iLay][r, c] * ph[c]
        end
        if iLay < nFilm
            cN = 6 * iLay + 1
            phN = exp.((-im * dIface) .* k3lay[iLay + 1])
            for c in 1:6, r in 1:3
                B[rD + r - 1, cN + c - 1] = -Ulay[iLay + 1][r, c] * phN[c]
                B[rT + r - 1, cN + c - 1] = -Tlay[iLay + 1][r, c] * phN[c]
            end
        else
            cSub = 6 * nFilm + 1
            for c in 1:3, r in 1:3
                B[rD + r - 1, cSub + c - 1] = -Ulay[nFilm + 1][r, c]
                B[rT + r - 1, cSub + c - 1] = -Tlay[nFilm + 1][r, c]
            end
        end
    end
    return B, surfPhase
end

# -----------------------------------------------------------------------------------------
# The per-point kernel
# -----------------------------------------------------------------------------------------

"Per-call constant context of the kernel."
struct _GreenCtx
    layers::Vector{GreenLayer}          # internal order: 1 = surface ... nLay = substrate
    vel::Vector{Float64}
    wAllD::Vector{ComplexF64}
    velRangeD::Vector{ComplexF64}
    isoFreq::Bool
    omega0::Float64
    cumDepthFromBelow::Vector{Float64}  # nFilm entries (empty for a half-space)
    rcondWarnThresh::Float64
end

"Per-task scratch of the kernel."
struct _GreenBuf
    lb::Vector{LayerPointBuffers}
    B::Matrix{ComplexF64}
    A::Matrix{ComplexF64}
    X::Matrix{ComplexF64}
    ipiv::Vector{Int}
    r::Vector{Float64}
    c::Vector{Float64}
    phase::Vector{ComplexF64}
    G::Matrix{ComplexF64}
end
function _GreenBuf(nLay::Int)
    nW = 6 * (nLay - 1) + 3
    return _GreenBuf([LayerPointBuffers() for _ in 1:nLay],
                     zeros(ComplexF64, nW, nW), zeros(ComplexF64, nW, nW),
                     zeros(ComplexF64, nW, nW), zeros(Int, nW), ones(nW), ones(nW),
                     zeros(ComplexF64, 6), zeros(ComplexF64, 3, 3))
end

"Optional full-output storage (`fullOutput = true`); internal layer order for the layered group."
struct _GreenFull
    dispField::Array{ComplexF64,4}
    xx3::Vector{Float64}
    k3Lay::Vector{Matrix{ComplexF64}}
    ULay::Vector{Array{ComplexF64,3}}
    ALay::Vector{Array{ComplexF64,3}}
    hScale::Vector{Float64}             # per velocity (the object keeps the last)
end

"Traction of one partial wave (unit k1, k2 = 0): `T[l] = sum_p (Ct1[l,p] + k3 Ct3[l,p]) U[p] / omega`."
@inline function _traction(L::GreenLayer, U::AbstractMatrix{ComplexF64}, n::Int,
                           k3::ComplexF64, invOmega::ComplexF64)
    @inbounds begin
        u1 = U[1, n]; u2 = U[2, n]; u3 = U[3, n]
        s1 = invOmega; s3 = k3 * invOmega
        t(l) = (L.Ct1[l, 1] * u1 + L.Ct1[l, 2] * u2 + L.Ct1[l, 3] * u3) * s1 +
               (L.Ct3[l, 1] * u1 + L.Ct3[l, 2] * u2 + L.Ct3[l, 3] * u3) * s3
        return (t(1), t(2), t(3))
    end
end

"`need` predicate for the substrate: only the depth-decaying waves get a polarisation."
_need_adm(k3, n) = imag(k3[n]) > 0

"""
    _green_point!(buf, ctx, iv, warm, full) -> (rc, usedFallback)

G_IJ(+k) (all nine pairs) at velocity `iv` into `buf.G`; returns the point's rcond and
whether the substrate admissibility needed the largest-Im(k3) fallback.
"""
function _green_point!(buf::_GreenBuf, ctx::_GreenCtx, iv::Int, warm::Bool, full)
    layers = ctx.layers
    nLay = length(layers)
    nFilm = nLay - 1
    w = ctx.wAllD[iv]
    omegaD = ctx.velRangeD[iv]
    v = ctx.vel[iv]
    invOmega = inv(omegaD)

    # --- partial waves, every layer (substrate: only the admissible ones need U)
    for iLay in 1:nLay
        need = (iLay == nLay && full === nothing) ? _need_adm : _need_all
        _layer_point!(buf.lb[iLay], layers[iLay], w, omegaD, warm, need; v = v, iLay = iLay)
    end
    bs = buf.lb[nLay]
    if count(z -> imag(z) > 0, bs.k3) != 3 && full === nothing
        # fallback picks may include a non-admissible wave: make sure all have a U
        for n in 1:6
            u = _null_vector(layers[nLay], bs.k3[n], layers[nLay].rho * omegaD^2)
            bs.U[1, n] = u[1]; bs.U[2, n] = u[2]; bs.U[3, n] = u[3]
        end
        bs.ident[1] != 0 && _degenerate_polarisations!(bs.U, layers[nLay], bs.k3, bs.ident, omegaD)
    end
    inds, usedFallback = substrate_wave_indices(bs.k3)

    B, A, X = buf.B, buf.A, buf.X
    G = buf.G
    fac = im * invOmega
    if nFilm == 0
        # --- half-space: 3x3 B (Every 2015 Eqs. 5-6)
        L = layers[1]
        @inbounds for m in 1:3
            n = inds[m]
            t1, t2, t3 = _traction(L, bs.U, n, bs.k3[n], invOmega)
            B[1, m] = t1; B[2, m] = t2; B[3, m] = t3
        end
        copyto!(A, B)
        _lu_inv!(X, A, buf.ipiv, 3)
        rc = _rcond_from_inv(B, X, 3)
        @inbounds for J in 1:3, I in 1:3
            s = zero(ComplexF64)
            for m in 1:3
                s += X[m, I] * bs.U[J, inds[m]]
            end
            G[I, J] = fac * s
        end
        if full !== nothing
            _halfspace_dispfield!(full, iv, X, bs, inds, fac)
        end
        return rc, usedFallback
    end

    # --- layered: (6 nFilm + 3)^2 B (Every 2002), calc.m buildBCMatrix in place
    nW = 6 * nFilm + 3
    hScale = ctx.isoFreq ? ctx.omega0 / v : 1.0
    cum = ctx.cumDepthFromBelow
    fill!(B, zero(ComplexF64))
    b1 = buf.lb[1]
    hEff1 = cum[1] * hScale
    @inbounds for c in 1:6
        buf.phase[c] = exp(-im * b1.k3[c] * hEff1)
        t1, t2, t3 = _traction(layers[1], b1.U, c, b1.k3[c], invOmega)
        ph = buf.phase[c]
        B[1, c] = t1 * ph; B[2, c] = t2 * ph; B[3, c] = t3 * ph
    end
    @inbounds for iLay in 1:nFilm
        rD = 3 + 6 * (iLay - 1) + 1
        rT = rD + 3
        dIface = iLay < nFilm ? cum[iLay + 1] * hScale : 0.0
        cS = 6 * (iLay - 1) + 1
        bl = buf.lb[iLay]
        for c in 1:6
            ph = exp(-im * bl.k3[c] * dIface)
            t1, t2, t3 = _traction(layers[iLay], bl.U, c, bl.k3[c], invOmega)
            col = cS + c - 1
            B[rD, col] = bl.U[1, c] * ph; B[rD + 1, col] = bl.U[2, c] * ph
            B[rD + 2, col] = bl.U[3, c] * ph
            B[rT, col] = t1 * ph; B[rT + 1, col] = t2 * ph; B[rT + 2, col] = t3 * ph
        end
        if iLay < nFilm
            cN = 6 * iLay + 1
            bn = buf.lb[iLay + 1]
            for c in 1:6
                ph = exp(-im * bn.k3[c] * dIface)
                t1, t2, t3 = _traction(layers[iLay + 1], bn.U, c, bn.k3[c], invOmega)
                col = cN + c - 1
                B[rD, col] = -bn.U[1, c] * ph; B[rD + 1, col] = -bn.U[2, c] * ph
                B[rD + 2, col] = -bn.U[3, c] * ph
                B[rT, col] = -t1 * ph; B[rT + 1, col] = -t2 * ph; B[rT + 2, col] = -t3 * ph
            end
        else
            cSub = 6 * nFilm + 1
            for m in 1:3
                n = inds[m]
                t1, t2, t3 = _traction(layers[nLay], bs.U, n, bs.k3[n], invOmega)
                col = cSub + m - 1
                B[rD, col] = -bs.U[1, n]; B[rD + 1, col] = -bs.U[2, n]
                B[rD + 2, col] = -bs.U[3, n]
                B[rT, col] = -t1; B[rT + 1, col] = -t2; B[rT + 2, col] = -t3
            end
        end
    end

    copyto!(A, B)
    _lu_inv!(X, A, buf.ipiv, nW)
    rc = _rcond_from_inv(B, X, nW)
    # EQUILIBRATION (calc.m): a thick film spreads B's column norms over many
    # orders of magnitude; scaling (exact algebra) fixes the reported conditioning. Gated
    # on the rcond check, skipped for a non-finite B (data lost at source).
    if rc < ctx.rcondWarnThresh && all(isfinite, B)
        r, c = buf.r, buf.c
        _equilibrate_pow2!(r, c, B, nW)
        @inbounds for j in 1:nW, i in 1:nW
            A[i, j] = r[i] * B[i, j] * c[j]
        end
        Beq = buf.B              # B no longer needed: reuse it for the scaled matrix
        copyto!(Beq, A)
        Xeq = similar(X)         # rare branch: allocation is fine here
        _lu_inv!(Xeq, A, buf.ipiv, nW)
        rcEq = _rcond_from_inv(Beq, Xeq, nW)
        if rcEq > rc
            rc = rcEq
            @inbounds for j in 1:nW, i in 1:nW
                X[i, j] = c[i] * Xeq[i, j] * r[j]
            end
        end
    end

    # G from film 1's six waves (rows 1:6 of invB) with the free-surface phase
    @inbounds for J in 1:3, I in 1:3
        s = zero(ComplexF64)
        for n in 1:6
            s += (X[n, I] * b1.U[J, n]) * buf.phase[n]
        end
        G[I, J] = fac * s
    end
    if full !== nothing
        _layered_store!(full, iv, X, buf, inds, fac, nFilm, hScale)
    end
    return rc, usedFallback
end

"Half-space depth profile `u(x3, ci, cj) = sum_n A[n, cj] U[ci, n] exp(i k3[n] x3)`, `A = (i/omega) invB`."
function _halfspace_dispfield!(full::_GreenFull, iv::Int, X, bs::LayerPointBuffers,
                               inds, fac::ComplexF64)
    xx3 = full.xx3
    D = full.dispField
    @inbounds for cj in 1:3, ci in 1:3
        for d in eachindex(xx3)
            D[d, ci, cj, iv] = 0
        end
        for m in 1:3
            n = inds[m]
            coef = (fac * X[m, cj]) * bs.U[ci, n]
            k3n = bs.k3[n]
            for d in eachindex(xx3)
                D[d, ci, cj, iv] += coef * exp(im * k3n * xx3[d])
            end
        end
    end
    return nothing
end

"Layered full output: every layer's k3 / U / amplitudes `(i/omega) invB[rows, 1:3]`."
function _layered_store!(full::_GreenFull, iv::Int, X, buf::_GreenBuf, inds, fac::ComplexF64,
                         nFilm::Int, hScale::Float64)
    nLay = nFilm + 1
    full.hScale[iv] = hScale
    @inbounds for iLay in 1:nLay
        bl = buf.lb[iLay]
        if iLay <= nFilm
            r0 = 6 * (iLay - 1)
            for n in 1:6
                full.k3Lay[iLay][n, iv] = bl.k3[n]
                for p in 1:3
                    full.ULay[iLay][p, n, iv] = bl.U[p, n]
                    full.ALay[iLay][n, p, iv] = fac * X[r0 + n, p]
                end
            end
        else
            r0 = 6 * nFilm
            for m in 1:3
                n = inds[m]
                full.k3Lay[iLay][m, iv] = bl.k3[n]
                for p in 1:3
                    full.ULay[iLay][p, m, iv] = bl.U[p, n]
                    full.ALay[iLay][m, p, iv] = fac * X[r0 + m, p]
                end
            end
        end
    end
    return nothing
end

"Conditioning / fallback tallies of one chunk."
struct _ChunkStats
    nearSing::Int
    vMin::Float64
    vMax::Float64
    worst::Float64
    laFallback::Int
end

"""
    _green_chunk!(store, ctx, ivs, full) -> _ChunkStats

Run the kernel over the velocity indices `ivs` (cold start at the first, warm after),
handing each point's 3x3 `G` to `store(iv, G)`.
"""
function _green_chunk!(store::S, ctx::_GreenCtx, ivs::UnitRange{Int}, full) where {S}
    buf = _GreenBuf(length(ctx.layers))
    nearSing = 0; vMin = Inf; vMax = -Inf; worst = Inf; laFallback = 0
    warm = false
    for iv in ivs
        rc, fb = _green_point!(buf, ctx, iv, warm, full)
        warm = true
        fb && (laFallback += 1)
        store(iv, buf.G)
        if rc < ctx.rcondWarnThresh
            v = ctx.vel[iv]
            nearSing += 1
            vMin = min(vMin, v); vMax = max(vMax, v); worst = min(worst, rc)
        end
    end
    return _ChunkStats(nearSing, vMin, vMax, worst, laFallback)
end

"""
    _green_run!(store, structure, angle, vel; eta, scanMode, omega0, threaded, full)
        -> nothing

Validate, set up the layers and the damped axis, run all chunks (threaded or not), and
issue `calc.m`'s two consolidated warnings.
"""
function _green_run!(store::S, structure::Structure, angle::Real, vel::Vector{Float64};
                     eta::Real, scanMode::Symbol, omega0, threaded::Bool,
                     full = nothing) where {S}
    isfinite(angle) || throw(ArgumentError("GreenSolver.calc: angle must be finite"))
    scanMode in (:isoWavevec, :isoFreq) || throw(ArgumentError(
        "GreenSolver:calc:badScanMode -- scanMode must be :isoWavevec or :isoFreq; got :$scanMode"))
    (isfinite(eta) && eta > 0) || throw(ArgumentError(
        "GreenSolver.calc: eta must be finite and positive; got $eta"))
    isempty(vel) && throw(ArgumentError("GreenSolver.calc: velRange must be non-empty"))
    (all(isfinite, vel) && all(>(0), vel)) || throw(ArgumentError(
        "GreenSolver.calc: velRange must be finite and positive"))
    om0 = omega0 === nothing ? maximum(vel) : Float64(omega0)
    (isfinite(om0) && om0 > 0) || throw(ArgumentError("GreenSolver.calc: omega0 must be finite and positive"))
    nLay = nlay(structure)
    nLay >= 1 || throw(ArgumentError(
        "GreenSolver:calc:emptyStructure -- Structure must contain at least one layer (substrate)."))

    layers = green_layers(structure, angle)
    nFilm = nLay - 1
    cum = Float64[]
    if nFilm > 0
        h = [2 * pi * structure.calcParams[nLay - iLay + 1].relDepth for iLay in 1:nFilm]
        cum = sum(h) .- vcat(0.0, cumsum(h[1:(end - 1)]))
    end
    wAllD, velRangeD = _damped_axis(vel, eta)
    ctx = _GreenCtx(layers, vel, wAllD, velRangeD, scanMode === :isoFreq, om0, cum,
                    (3 + 6 * nFilm) * eps(Float64))

    nVel = length(vel)
    chunks = [i:min(i + GREEN_CHUNK - 1, nVel) for i in 1:GREEN_CHUNK:nVel]
    stats = Vector{_ChunkStats}(undef, length(chunks))
    if threaded && length(chunks) > 1 && Threads.nthreads() > 1
        tasks = [Threads.@spawn _green_chunk!(store, ctx, ch, full) for ch in chunks]
        for (k, t) in enumerate(tasks)
            stats[k] = fetch(t)
        end
    else
        for (k, ch) in enumerate(chunks)
            stats[k] = _green_chunk!(store, ctx, ch, full)
        end
    end

    laFallback = sum(s.laFallback for s in stats)
    nearSing = sum(s.nearSing for s in stats)
    if laFallback > 0
        @warn string("GreenSolver:calc:absorptionFallback -- $laFallback/$nVel ",
                     "substrate selections did not yield exactly 3 admissible partial waves ",
                     "under the limiting-absorption test and fell back to a largest-Im(k3) ",
                     "tiebreak.")
    end
    if nearSing > 0
        vMin = minimum(s.vMin for s in stats)
        vMax = maximum(s.vMax for s in stats)
        worst = minimum(s.worst for s in stats)
        @warn string("GreenSolver:calc:nearSingular -- $nearSing/$nVel velocity points ",
                     "(v = $(round(vMin; sigdigits = 4))-$(round(vMax; sigdigits = 4)) km/s) ",
                     "had a near-singular boundary-condition matrix (worst RCOND = ",
                     "$(round(worst; sigdigits = 4)) < threshold ",
                     "$(round(ctx.rcondWarnThresh; sigdigits = 4))) -- G values there may be ",
                     "unreliable. See docs/reference/architecture/green.md, \"Conditioning\".")
    end
    return ctx
end

_vel_vector(velRange::AbstractVector{<:Real}) = Float64.(collect(velRange))
_vel_vector(velRange::Real) = [Float64(velRange)]

"+1 for directions 2, 3; -1 for 1 (the propagation index): (-1)^n of the +-k parity."
@inline _refl_sign(k::Integer) = k == 1 ? -1.0 : 1.0

"Store sink: the full 3x3 tensor per point."
struct _StoreFull
    GFull::Array{ComplexF64,3}
end
@inline function (s::_StoreFull)(iv::Int, G::Matrix{ComplexF64})
    @inbounds for J in 1:3, I in 1:3
        s.GFull[I, J, iv] = G[I, J]
    end
end

"Store sink: one scalar per point -- G_IJ(+k), G_IJ(-k) = G_JI(+k), or the parity-signed mean."
struct _StoreComponent
    out::Vector{ComplexF64}
    I::Int
    J::Int
    side::Symbol
end
@inline function (s::_StoreComponent)(iv::Int, G::Matrix{ComplexF64})
    @inbounds begin
        gp = G[s.I, s.J]
        gm = G[s.J, s.I]
        s.out[iv] = s.side === :plus ? gp : s.side === :minus ? gm :
                    (gp + _refl_sign(s.I) * _refl_sign(s.J) * gm) / 2
    end
end

# -----------------------------------------------------------------------------------------
# Public entry points
# -----------------------------------------------------------------------------------------

"""
    calc(GreenSolver, structure, angle, velRange = range(0.1, 10, 5000);
         forceDir = 1, displDir = 3, scanMode = :isoWavevec, omega0 = nothing,
         computeMinus = true, eta = LIMITING_ABSORPTION_ETA, fullOutput = false,
         threaded = true) -> GreenSolver

Surface dynamic Green's function `G_IJ(k_par, omega)` of a [`Structure`](@ref), half-space
or layered (Every et al. 1997, Every 2002, Maznev et al. 1999). Port of
`guliwerr/solvers/green/@GreenSolver/calc.m` (MATLAB's name-value options `'ScanMode'`,
`'Omega0'`, `'ComputeMinus'`, `'Eta'` are keywords here, as are `forceDir`/`displDir`).

- `angle` -- in-plane propagation angle [deg]. Callers holding an already-rotated
  structure pass `0` (the double-rotation rule).
- `velRange` -- phase-velocity axis [km/s] (a scalar is accepted).
- `scanMode` -- `:isoWavevec` (|k_par| = 1, omega = v) or `:isoFreq` (omega = `omega0`,
  |k_par| = omega0/v; `omega0` defaults to `maximum(velRange)`).
- `eta` -- limiting-absorption damping, `omega^2 -> omega^2 (1 + i eta)` everywhere;
  larger values emulate a finite measurement linewidth (the scores use `SawSpectrum.ETA`).
- `computeMinus` -- fill `G_minus`, `G_eff`, `GFull_minus`, `GFull_eff` (default `true`).
  Free either way: only k1 = +1 is solved and `G_IJ(-k) = G_JI(+k)` (reciprocity, exact
  for any tensor, cut, stack, damped or not).
- `fullOutput` -- MATLAB's `nargout > 1`: half-space `dispField`/`depthAxis`; layered
  `k3OutLayered`/`polarizationULayered`/`partialAmplitudeLayered`/`hScaleOut` (Structure
  order). Intended for a scalar velocity.
- `threaded` -- spread the velocity chunks over Julia's threads (results are bitwise
  independent of it).

`G_eff = (G_+ + (-1)^n G_-)/2`, `n` = how many of `(forceDir, displDir)` equal 1: e.g.
`G13_eff = (G13_+ - G31_+)/2`; for a diagonal pair `G_eff = G_+` exactly.

For the scores' hot loops, [`green_component`](@ref) returns one complex component per
velocity without building the object or the 3x3xnVel tensors.
"""
function calc(::Type{GreenSolver}, structure::Structure, angle::Real,
               velRange::Union{Real,AbstractVector{<:Real}} = range(0.1, 10.0, length = 5000);
               forceDir::Integer = 1, displDir::Integer = 3,
               scanMode::Symbol = :isoWavevec, omega0::Union{Real,Nothing} = nothing,
               computeMinus::Bool = true, eta::Real = LIMITING_ABSORPTION_ETA,
               fullOutput::Bool = false, threaded::Bool = true)
    1 <= forceDir <= 3 || throw(ArgumentError(
        "GreenSolver.calc: forceDir must be 1, 2 or 3; got $forceDir"))
    1 <= displDir <= 3 || throw(ArgumentError(
        "GreenSolver.calc: displDir must be 1, 2 or 3; got $displDir"))
    vel = _vel_vector(velRange)
    nVel = length(vel)
    nLay = nlay(structure)
    nFilm = nLay - 1

    GFull = fill(ComplexF64(NaN), 3, 3, max(nVel, 0))
    full = nothing
    if fullOutput && nVel > 0
        wantHS = nFilm == 0
        nW(iLay) = iLay <= nFilm ? 6 : 3
        full = _GreenFull(
            wantHS ? zeros(ComplexF64, DISP_FIELD_NDEPTH, 3, 3, nVel) :
                     Array{ComplexF64,4}(undef, 0, 0, 0, 0),
            wantHS ? collect(range(0.0, DISP_FIELD_MAX_DEPTH, length = DISP_FIELD_NDEPTH)) :
                     Float64[],
            wantHS ? Matrix{ComplexF64}[] : [fill(ComplexF64(NaN), nW(i), nVel) for i in 1:nLay],
            wantHS ? Array{ComplexF64,3}[] : [fill(ComplexF64(NaN), 3, nW(i), nVel) for i in 1:nLay],
            wantHS ? Array{ComplexF64,3}[] : [fill(ComplexF64(NaN), nW(i), 3, nVel) for i in 1:nLay],
            fill(NaN, nVel))
    end
    _green_run!(_StoreFull(GFull), structure, angle, vel; eta = eta, scanMode = scanMode,
                omega0 = omega0, threaded = threaded, full = full)

    G = GFull[forceDir, displDir, :]
    if computeMinus
        GFullMinus = permutedims(GFull, (2, 1, 3))
        Gminus = GFull[displDir, forceDir, :]
        sgn = _refl_sign(forceDir) * _refl_sign(displDir)
        Geff = (G .+ sgn .* Gminus) ./ 2
        signMat = [_refl_sign(I) * _refl_sign(J) for I in 1:3, J in 1:3]
        GFullEff = (GFull .+ signMat .* GFullMinus) ./ 2
    else
        GFullMinus = Array{ComplexF64,3}(undef, 0, 0, 0)
        GFullEff = Array{ComplexF64,3}(undef, 0, 0, 0)
        Gminus = ComplexF64[]
        Geff = ComplexF64[]
    end

    wantHS = full !== nothing && nFilm == 0
    wantLay = full !== nothing && nFilm > 0
    return GreenSolver(;
        forceDir = Int(forceDir), displDir = Int(displDir),
        velRange = vel, computeMinus = computeMinus,
        G = G, GFull = GFull,
        G_minus = Gminus, G_eff = Geff, GFull_minus = GFullMinus, GFull_eff = GFullEff,
        dispField = wantHS ? full.dispField : Array{ComplexF64,4}(undef, 0, 0, 0, 0),
        depthAxis = wantHS ? full.xx3 : Float64[],
        # internal (surface-first) -> Structure order (1 = substrate ... nLay = top film)
        k3OutLayered = wantLay ? reverse(full.k3Lay) : Matrix{ComplexF64}[],
        polarizationULayered = wantLay ? reverse(full.ULay) : Array{ComplexF64,3}[],
        partialAmplitudeLayered = wantLay ? reverse(full.ALay) : Array{ComplexF64,3}[],
        hScaleOut = wantLay ? full.hScale[end] : NaN)
end

"""
    green_component(structure, angle, velRange, I = 1, J = 3; side = :eff,
                    eta = LIMITING_ABSORPTION_ETA, scanMode = :isoWavevec,
                    omega0 = nothing, threaded = true) -> Vector{ComplexF64}

One complex Green component per velocity, straight from the kernel (no `GreenSolver`
object, no 3x3xnVel tensors): `side = :plus` -> `G_IJ(+k)` (MATLAB `obj.G`), `:minus` ->
`G_IJ(-k) = G_JI(+k)` (`obj.G_minus`), `:eff` -> the parity-signed mean (`obj.G_eff`;
`(1, 3)` gives G13_eff = (G13_+ - G31_+)/2). Bitwise identical to the corresponding
field of `calc(GreenSolver, ...)`. This is the call for the classification scores
(`CalcCommon.greenScoreSpectrum`, `GreenScore` probes, `SpotGreen.solveOne`,
`GreenSpotScore.solveOne`).
"""
function green_component(structure::Structure, angle::Real,
                         velRange::Union{Real,AbstractVector{<:Real}},
                         I::Integer = 1, J::Integer = 3; side::Symbol = :eff,
                         eta::Real = LIMITING_ABSORPTION_ETA, scanMode::Symbol = :isoWavevec,
                         omega0::Union{Real,Nothing} = nothing, threaded::Bool = true)
    (1 <= I <= 3 && 1 <= J <= 3) || throw(ArgumentError(
        "green_component: I and J must be in 1:3; got ($I, $J)"))
    side in (:plus, :minus, :eff) || throw(ArgumentError(
        "green_component: side must be :plus, :minus or :eff; got :$side"))
    vel = _vel_vector(velRange)
    out = fill(ComplexF64(NaN), length(vel))
    _green_run!(_StoreComponent(out, Int(I), Int(J), side), structure, angle, vel;
                eta = eta, scanMode = scanMode, omega0 = omega0, threaded = threaded)
    return out
end

"""
    green_gfull(structure, angle, velRange; eta, scanMode, omega0, threaded)
        -> Array{ComplexF64,3}                                     # 3 x 3 x nVel

`G_IJ(+k)` for all nine pairs (MATLAB `obj.GFull`), without the object; `G_IJ(-k)` is its
transpose in `(I, J)`.
"""
function green_gfull(structure::Structure, angle::Real,
                     velRange::Union{Real,AbstractVector{<:Real}};
                     eta::Real = LIMITING_ABSORPTION_ETA, scanMode::Symbol = :isoWavevec,
                     omega0::Union{Real,Nothing} = nothing, threaded::Bool = true)
    vel = _vel_vector(velRange)
    GFull = fill(ComplexF64(NaN), 3, 3, length(vel))
    _green_run!(_StoreFull(GFull), structure, angle, vel; eta = eta, scanMode = scanMode,
                omega0 = omega0, threaded = threaded)
    return GFull
end

export calc_traction_cols, calc_gfull_from_partial_waves, build_bc_matrix,
       green_component, green_gfull, DISP_FIELD_NDEPTH, DISP_FIELD_MAX_DEPTH, GREEN_CHUNK
