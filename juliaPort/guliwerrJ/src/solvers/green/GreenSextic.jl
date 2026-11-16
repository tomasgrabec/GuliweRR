# =========================================================================================
# GreenSextic.jl -- the boundary-condition-independent half of
# `guliwerr/solvers/green/@GreenSolver/calc.m`: layer setup, the sextic k3 polynomial, its
# roots, the degeneracy scan and the partial-wave polarisations. The sextic is that of the
# partial-wave method (Farnell 1970; Every et al. 1997, Eq. 10), solved as a polynomial in
# k3 -- not as the 6x6 eigenvalue problem of the Stroh (1958) / Barnett-Lothe (1973)
# formalism, which holds the same roots.
#
# As in MATLAB's calc.m:
#   * only the k1 = +1 branch is ever solved (G(-k) = G(+k)^T by reciprocity), so there is
#     no k-sign dimension anywhere;
#   * the damping eta is an argument (calc.m 'Eta'), LIMITING_ABSORPTION_ETA its default;
#   * the degenerate-eigenvalue match compares against real(omega), and the "eigenvalue not
#     found" fallback stays orthogonal to an assigned partner;
#   * the near-isotropy screen reads ElastConsts' stored shearAnisotropyRatio when present,
#     scanning only when it is missing or NaN (calc.m).
#
# STRUCTURE (Julia-side, not a MATLAB mirror). MATLAB batches the whole velocity axis per
# layer (`roots` loop, `pagesvd`) because MATLAB's per-call overhead dominates. Julia has no
# such overhead, so the hot path here is a PER-POINT kernel, `_layer_point!`, working on a
# handful of preallocated 6-/3x6-element buffers with no allocation:
#   1. sextic coefficients at this w (7 complex numbers, from the 7x4 DetW);
#   2. roots by Aberth-Ehrlich iteration warm-started from the previous velocity's roots
#      (cold start and any non-convergence: `matlab_roots`, i.e. the companion-matrix
#      eigenvalues MATLAB's `roots` uses);
#   3. canonical (real, imag) order -- the port's long-standing deliberate divergence
#      (`canonicalize_k3`), which also makes the root-finder's own order irrelevant;
#   4. the degeneracy scan of `CheckIdenticalRecords` for one column;
#   5. polarisations: the null vector of the 3x3 Christoffel matrix as the largest of the
#      three row cross products (exact for a rank-2 matrix; MATLAB's SVD gives the same ray),
#      the closed-form SH/SV basis / `eigen` branch only at a degenerate pair.
# The whole-axis functions (`solve_k3_allvel`, `check_identical_records_allvel`,
# `calc_u_allvel`, `green_sextic_prepass`) remain as thin wrappers -- they are what is
# compared against MATLAB's own pre-pass, and a diagnostic API.
# =========================================================================================

using LinearAlgebra

using ...CoreTypes: Structure, Material, nlay, SetDirMat, GetRotElastTensor,
                    GetStandardDirMat, matlab_roots
using ..Christoffel: calc_anisotropy_ratio

# -----------------------------------------------------------------------------------------
# Numerical constants (calc.m)
# -----------------------------------------------------------------------------------------

"""Frequency-match window for degenerate-k3 detection. `calc.m` `OMEGA_TOLERANCE`."""
const OMEGA_TOLERANCE = 5e-3

"""k3 degeneracy gap threshold of the degeneracy scan. `calc.m` `ZERO_TOLERANCE`."""
const ZERO_TOLERANCE = 1e-4

"""
Default limiting-absorption damping: the whole solve runs at complex
`w = omega^2 * (1 + i*eta)` (retarded Green's function, Every et al. 1997, p. 1347). `calc.m`
`LIMITING_ABSORPTION_ETA` = 3e-4, the default of `calc`'s `eta` keyword. Larger values
(e.g. the classification scores' `SawSpectrum.ETA`) emulate a finite measurement linewidth.
"""
const LIMITING_ABSORPTION_ETA = 3e-4

"""Coarse angular scan step [deg] for the near-isotropy screen. `calc.m`."""
const NEAR_ISO_ANGLE_STEP = 15.0

"""`(fastest shear / slowest shear)^2` threshold below which a layer is treated as isotropic. `calc.m`."""
const NEAR_ISO_WARN_THRESH = 1.01

# -----------------------------------------------------------------------------------------
# Small shared helpers -- calc.m's kVecCols / outerCols
# -----------------------------------------------------------------------------------------

"""
    kvec_cols(k1, k2, k3) -> Matrix{ComplexF64}

`3 x nW` wave-vector columns `[k1; k2; k3[n]]`. Port of `calc.m` `kVecCols`.
"""
function kvec_cols(k1::Number, k2::Number, k3::AbstractArray{<:Number})
    n = length(k3)
    K = Matrix{ComplexF64}(undef, 3, n)
    @inbounds for j in 1:n
        K[1, j] = k1
        K[2, j] = k2
        K[3, j] = k3[j]
    end
    return K
end

"""
    outer_cols(A, B) -> Matrix

Column-wise 3x3 outer products flattened column-major: `P[(q-1)*3 + p, n] == A[p, n] *
B[q, n]`. No adjoint anywhere (complex inputs). Port of `calc.m` `outerCols`.
"""
function outer_cols(A::AbstractMatrix, B::AbstractMatrix)
    n = size(A, 2)
    size(B, 2) == n ||
        throw(DimensionMismatch("outer_cols: A and B must have the same number of columns"))
    return reshape(reshape(A, 3, 1, n) .* reshape(B, 1, 3, n), 9, n)
end

"Full 2-D convolution (MATLAB `conv2(A, B)`, shape 'full') for the tiny DetW products."
function _conv2_full(A::AbstractMatrix{TA}, B::AbstractMatrix{TB}) where {TA<:Number,TB<:Number}
    ma, na = size(A)
    mb, nb = size(B)
    C = zeros(promote_type(TA, TB), ma + mb - 1, na + nb - 1)
    @inbounds for jb in 1:nb, ib in 1:mb
        b = B[ib, jb]
        iszero(b) && continue
        for ja in 1:na, ia in 1:ma
            C[ia + ib - 1, ja + jb - 1] += A[ia, ja] * b
        end
    end
    return C
end

# -----------------------------------------------------------------------------------------
# Layer setup + near-isotropy screen -- calc.m
# -----------------------------------------------------------------------------------------

"Stored shear anisotropy ratio of a layer's constants if the port has one, else NaN."
_stored_aniso_ratio(ec) = hasproperty(ec, :shearAnisotropyRatio) ?
    Float64(getproperty(ec, :shearAnisotropyRatio)) : NaN

"""
    green_layer_setup(structure, angle) -> (elTen, density, isIsotropic)

Rotated elastic tensors, densities and per-layer near-isotropy flags in GreenSolver's
**internal** layer order (`iLay = 1` = surface, `iLay = nLay` = substrate; the reverse of
`Structure`'s). Each tensor is in the standard frame (`x1` = propagation, `x3` = surface
normal) rotated in-plane to `angle`.

`isIsotropic[iLay]`: the layer's shear anisotropy ratio is below `NEAR_ISO_WARN_THRESH`.
The ratio is the stored `ElastConsts.shearAnisotropyRatio` when finite (exact for isotropic
and cubic constants), otherwise a `NEAR_ISO_ANGLE_STEP` scan of the rotated tensor --
`calc.m`. The notice is logged once per session (`maxlog = 1`), the analogue of
MATLAB's warning-state dedup.
"""
function green_layer_setup(structure::Structure, angle::Real)
    nLay = nlay(structure)
    nLay >= 1 || throw(ArgumentError(
        "GreenSolver:calc:emptyStructure -- Structure must contain at least one layer (substrate)."))

    elTen = Vector{Array{Float64,4}}(undef, nLay)
    density = Vector{Float64}(undef, nLay)
    isIsotropic = Vector{Bool}(undef, nLay)
    for iLay in 1:nLay
        vIdx = nLay - iLay + 1                 # index into structure.mater (reversed)
        mat0 = structure.mater[vIdx]
        matR = SetDirMat(mat0, GetStandardDirMat(mat0.dirMat, angle))
        elTen[iLay] = GetRotElastTensor(matR, 0, 0)
        density[iLay] = mat0.rho

        anisoRatio = _stored_aniso_ratio(mat0.elConsts)
        isfinite(anisoRatio) || (anisoRatio = calc_anisotropy_ratio(elTen[iLay], NEAR_ISO_ANGLE_STEP))
        isIsotropic[iLay] = anisoRatio < NEAR_ISO_WARN_THRESH
        if isIsotropic[iLay]
            @warn string("GreenSolver:calc:nearIsotropic -- Layer $iLay/$nLay (internal order, ",
                         "surface-first) has shear anisotropy ratio ",
                         round(anisoRatio; digits = 6), " (< $NEAR_ISO_WARN_THRESH threshold) ",
                         "-- GreenSolver uses a closed-form SH/SV basis for the degenerate shear ",
                         "branch instead of eigen(). Shown once per session.") maxlog = 1
        end
    end
    return elTen, density, isIsotropic
end

# -----------------------------------------------------------------------------------------
# Calc_DetW_byChristoffel -- calc.m
# -----------------------------------------------------------------------------------------

"""
    calc_detw_by_christoffel(elTen, rho, k1, k2) -> Matrix{ComplexF64}   # 7 x 4

Velocity-independent part of the sextic k3-finder: `DetW[i, j]` is the coefficient of
`P^(7-i) * w^(j-1)` (`P = -i*k3`, `w = omega^2`), so the polynomial at one velocity is
`DetW * [1; w; w^2; w^3]`. Exact, via bivariate `conv2` products of the Christoffel
quadratics. Port of `calc.m` `Calc_DetW_byChristoffel`; `calc` only builds `k1 = +1`.
"""
function calc_detw_by_christoffel(elTen::AbstractArray{<:Real,4}, rho::Real,
                                   k1::Real, k2::Real)
    d2 = Matrix{Matrix{ComplexF64}}(undef, 3, 3)
    @inbounds for ll in 1:3, n in 1:3
        a_ln = -elTen[ll, 1, n, 1] * k1 * k1 -
               (elTen[ll, 1, n, 2] + elTen[ll, 2, n, 1]) * k1 * k2 -
               elTen[ll, 2, n, 2] * k2 * k2
        b_ln = -im * (elTen[ll, 1, n, 3] + elTen[ll, 3, n, 1]) * k1 -
               im * (elTen[ll, 2, n, 3] + elTen[ll, 3, n, 2]) * k2
        c_ln = complex(elTen[ll, 3, n, 3])
        M = zeros(ComplexF64, 3, 2)
        M[1, 1] = c_ln
        M[2, 1] = b_ln
        M[3, 1] = a_ln
        M[3, 2] = rho * (ll == n)
        d2[ll, n] = M
    end
    tp(x, y, z) = _conv2_full(x, _conv2_full(y, z))
    return tp(d2[3, 3], d2[1, 1], d2[2, 2]) .+
           tp(d2[3, 1], d2[1, 2], d2[2, 3]) .+
           tp(d2[2, 1], d2[1, 3], d2[3, 2]) .-
           tp(d2[3, 1], d2[1, 3], d2[2, 2]) .-
           tp(d2[3, 3], d2[1, 2], d2[2, 1]) .-
           tp(d2[1, 1], d2[2, 3], d2[3, 2])
end

# -----------------------------------------------------------------------------------------
# Per-layer constant data for the per-point kernel
# -----------------------------------------------------------------------------------------

"""
    GreenLayer

Everything the per-point kernel needs about one layer, built once per `calc` call
(`green_layers`). With `k1 = 1`, `k2 = 0`:

* Christoffel matrix `Gamma(k3)[l, j] = A0[l, j] + k3*A1[l, j] + k3^2*A2[l, j]`
  (`A0 = C[l,1,j,1]`, `A1 = C[l,1,j,3] + C[l,3,j,1]`, `A2 = C[l,3,j,3]`);
* traction of a partial wave `T[l] = sum_p (Ct1[l, p] + k3*Ct3[l, p]) * U[p] / omega`
  (`Ct1 = C[3,l,p,1]`, `Ct3 = C[3,l,p,3]`) -- `calcTractionCols` with the slowness
  `[1, 0, k3]/omega` written out;
* `DetW`, the 7x4 sextic coefficient table.
"""
struct GreenLayer
    DetW::Matrix{ComplexF64}
    A0::Matrix{Float64}
    A1::Matrix{Float64}
    A2::Matrix{Float64}
    Ct1::Matrix{Float64}
    Ct3::Matrix{Float64}
    elTen::Array{Float64,4}
    rho::Float64
    isIsotropic::Bool
end

function GreenLayer(elTen::Array{Float64,4}, rho::Real, isIsotropic::Bool)
    A0 = [elTen[l, 1, j, 1] for l in 1:3, j in 1:3]
    A1 = [elTen[l, 1, j, 3] + elTen[l, 3, j, 1] for l in 1:3, j in 1:3]
    A2 = [elTen[l, 3, j, 3] for l in 1:3, j in 1:3]
    Ct1 = [elTen[3, l, p, 1] for l in 1:3, p in 1:3]
    Ct3 = [elTen[3, l, p, 3] for l in 1:3, p in 1:3]
    DetW = calc_detw_by_christoffel(elTen, rho, 1.0, 0.0)
    return GreenLayer(DetW, A0, A1, A2, Ct1, Ct3, elTen, Float64(rho), isIsotropic)
end

"""
    green_layers(structure, angle) -> Vector{GreenLayer}

Per-layer kernel data in internal (surface-first) order; see [`green_layer_setup`](@ref).
"""
function green_layers(structure::Structure, angle::Real)
    elTen, density, isIso = green_layer_setup(structure, angle)
    return [GreenLayer(elTen[i], density[i], isIso[i]) for i in eachindex(density)]
end

"""
    _damped_axis(vel, eta) -> (wAllD, velRangeD)

`wAllD = v.^2 * (1 + i*eta)`, `velRangeD = sqrt.(wAllD)` (`calc.m`). Every quantity
of the solve uses these, never the real `v`.
"""
function _damped_axis(vel::AbstractVector{<:Real}, eta::Real)
    wAllD = [v^2 * (1 + im * eta) for v in vel]
    return wAllD, sqrt.(wAllD)
end

# -----------------------------------------------------------------------------------------
# Sextic roots: Aberth-Ehrlich with warm start, matlab_roots as cold start / fallback
# -----------------------------------------------------------------------------------------

"Sextic coefficients (descending in P) at damped `w`: `DetW * [1; w; w^2; w^3]`."
@inline function _sextic_coeffs!(c::AbstractVector{ComplexF64}, DetW::Matrix{ComplexF64},
                                 w::ComplexF64)
    w2 = w * w
    w3 = w2 * w
    @inbounds for k in 1:7
        c[k] = DetW[k, 1] + DetW[k, 2] * w + DetW[k, 3] * w2 + DetW[k, 4] * w3
    end
    return c
end

"""
    _aberth!(z, c; maxit = 50) -> Bool

Refine the six roots `z` of the sextic `c` (descending coefficients) in place by
Aberth-Ehrlich iteration (Gauss-Seidel form). `z` must hold six distinct starting values --
here the previous velocity's roots, which are within one axis step of the new ones, so 2-4
iterations reach full precision. Stops once every relative correction is below `1e-12` and
then applies one more polishing sweep (the cubic convergence takes that to rounding level).
Returns `false` -- leaving `z` unusable -- if it does not converge, or meets a coincident
pair or a non-finite correction; the caller then falls back to `matlab_roots`.
"""
function _aberth!(z::AbstractVector{ComplexF64}, c::AbstractVector{ComplexF64};
                  maxit::Int = 50)
    converged = false
    @inbounds for it in 1:maxit
        maxrel = 0.0
        for i in 1:6
            zi = z[i]
            p = c[1]
            dp = zero(ComplexF64)
            for k in 2:7
                dp = muladd(dp, zi, p)
                p = muladd(p, zi, c[k])
            end
            iszero(p) && continue                 # exact root
            ratio = p / dp
            s = zero(ComplexF64)
            for j in 1:6
                j == i && continue
                d = zi - z[j]
                iszero(d) && return false
                s += inv(d)
            end
            wi = ratio / (1 - ratio * s)
            isfinite(wi) || return false
            z[i] = zi - wi
            maxrel = max(maxrel, abs(wi) / abs(zi))
        end
        converged && return true                  # this was the polishing sweep
        maxrel < 1e-12 && (converged = true)
    end
    return false
end

"Stable insertion sort of `k3` (and the companion permutation `perm`) by `(real, imag)`."
@inline function _canonical_sort!(k3::AbstractVector{ComplexF64})
    @inbounds for i in 2:6
        x = k3[i]
        j = i - 1
        while j >= 1 && (real(k3[j]) > real(x) ||
                         (real(k3[j]) == real(x) && imag(k3[j]) > imag(x)))
            k3[j + 1] = k3[j]
            j -= 1
        end
        k3[j + 1] = x
    end
    return k3
end

"""
    _ident_col!(ident, pos, key, k3, zeroTolerance) -> ident

`CheckIdenticalRecords_allVel` for one velocity column (`calc.m`): sort the six
roots by real part (by imaginary part when `mean(|real|) <= tol`), stably; a consecutive
gap `|k3[pos[ii+1]] - k3[pos[ii]]| < tol` records the pair; at most two pairs, in gap
order. `ident` (6 Ints) gets the indices, zeros elsewhere; `pos`/`key` are scratch.
"""
function _ident_col!(ident::AbstractVector{Int}, pos::AbstractVector{Int},
                     key::AbstractVector{Float64}, k3::AbstractVector{ComplexF64},
                     zeroTolerance::Real)
    @inbounds begin
        s = 0.0
        for n in 1:6
            s += abs(real(k3[n]))
        end
        useReal = s / 6 > zeroTolerance
        for n in 1:6
            key[n] = useReal ? real(k3[n]) : imag(k3[n])
            pos[n] = n
            ident[n] = 0
        end
        for i in 2:6                               # stable insertion sort of pos by key
            p = pos[i]
            kp = key[p]
            j = i - 1
            while j >= 1 && key[pos[j]] > kp
                pos[j + 1] = pos[j]
                j -= 1
            end
            pos[j + 1] = p
        end
        cnt = 1
        for ii in 1:5
            if abs(k3[pos[ii + 1]] - k3[pos[ii]]) < zeroTolerance && cnt <= 4
                ident[cnt] = pos[ii]
                ident[cnt + 1] = pos[ii + 1]
                cnt += 2
            end
        end
    end
    return ident
end

# -----------------------------------------------------------------------------------------
# Polarisations
# -----------------------------------------------------------------------------------------

"""
    _null_vector(L, k3, wr) -> NTuple{3,ComplexF64}

Unit null vector of `Gamma(k3) - wr*I` (`wr = rho*omegaD^2`): the largest of the three
row cross products. For a rank-2 matrix every row pair's cross product (bilinear, no
conjugation) is orthogonal to all rows and hence spans the null space; taking the largest
avoids a nearly parallel pair, as MATLAB calc.m (Calc_U_allVel). Its phase is arbitrary in both languages and cancels in every G (each U
column scales B's column and invB's row oppositely).
"""
@inline function _null_vector(L::GreenLayer, k3::ComplexF64, wr::ComplexF64)
    A0, A1, A2 = L.A0, L.A1, L.A2
    k32 = k3 * k3
    @inbounds begin
        g(l, j) = A0[l, j] + k3 * A1[l, j] + k32 * A2[l, j] - (l == j ? wr : zero(wr))
        a11 = g(1, 1); a12 = g(1, 2); a13 = g(1, 3)
        a21 = g(2, 1); a22 = g(2, 2); a23 = g(2, 3)
        a31 = g(3, 1); a32 = g(3, 2); a33 = g(3, 3)
    end
    # r1 x r2, r1 x r3, r2 x r3
    x1 = a12 * a23 - a13 * a22; y1 = a13 * a21 - a11 * a23; z1 = a11 * a22 - a12 * a21
    x2 = a12 * a33 - a13 * a32; y2 = a13 * a31 - a11 * a33; z2 = a11 * a32 - a12 * a31
    x3 = a22 * a33 - a23 * a32; y3 = a23 * a31 - a21 * a33; z3 = a21 * a32 - a22 * a31
    n1 = abs2(x1) + abs2(y1) + abs2(z1)
    n2 = abs2(x2) + abs2(y2) + abs2(z2)
    n3 = abs2(x3) + abs2(y3) + abs2(z3)
    # Numerically rank-1 matrix (largest cross product < 1e-6 of the matrix's own scale): a
    # double eigenvalue not paired as degenerate, e.g. both shear waves at k3 ~ 0 on a
    # symmetry axis with eta ~ 1e-12. The cross product is then noise and the pair would get
    # the same vector (singular B); use the SVD's null vector, as MATLAB's calc.m does.
    # Never taken at the default eta.
    nb = max(n1, n2, n3)
    frob2 = abs2(a11) + abs2(a12) + abs2(a13) + abs2(a21) + abs2(a22) + abs2(a23) +
            abs2(a31) + abs2(a32) + abs2(a33)
    if !(nb > 1e-12 * frob2^2) || !isfinite(nb)
        return _svd_null_vector(a11, a12, a13, a21, a22, a23, a31, a32, a33)
    end
    if n1 >= n2 && n1 >= n3
        s = inv(sqrt(n1)); return (x1 * s, y1 * s, z1 * s)
    elseif n2 >= n3
        s = inv(sqrt(n2)); return (x2 * s, y2 * s, z2 * s)
    else
        s = inv(sqrt(n3)); return (x3 * s, y3 * s, z3 * s)
    end
end

"Last right singular vector of the 3x3 matrix (MATLAB `pagesvd` `V(:,3)`): the rank-1 fallback of `_null_vector`."
@noinline function _svd_null_vector(a11, a12, a13, a21, a22, a23, a31, a32, a33)
    v = svd(ComplexF64[a11 a12 a13; a21 a22 a23; a31 a32 a33]).V
    return (v[1, 3], v[2, 3], v[3, 3])
end

"Christoffel matrix `Gamma(k3)` (k1 = 1, k2 = 0) as a dense 3x3 -- the degenerate branch only."
_gamma_matrix(L::GreenLayer, k3::ComplexF64) =
    ComplexF64[L.A0[l, j] + k3 * L.A1[l, j] + k3 * k3 * L.A2[l, j] for l in 1:3, j in 1:3]

"""
    _degenerate_polarisations!(U, L, k3, ident, omega; omegaTolerance, zeroTolerance)

Fill `U[:, ii]` for every degenerate index `ii` of this column (ascending order -- the
generic branch reads its partner's already-assigned column). Port of the degenerate part
of `calc.m` `Calc_U_allVel`:

* isotropic / near-isotropic layer: closed-form SH (`[0,1,0]`, first of the pair) / SV
  (`k x eSH`, normalised) basis;
* generic anisotropic: `eigen` of `Gamma` (unsorted, as MATLAB's `eig`), eigenvalues
  matched against `real(omega)` (`omega` is the damped frequency, and its
  imaginary part must not make a genuine match fail), orthogonal to an assigned partner;
  if fewer than two eigenvalues match, warn and still pick the eigenvector most orthogonal
  to an assigned partner.

The Christoffel matrix is complex SYMMETRIC (not Hermitian) for complex k3 -- it is
deliberately not symmetrised.
"""
function _degenerate_polarisations!(U::AbstractMatrix{ComplexF64}, L::GreenLayer,
                                    k3::AbstractVector{ComplexF64},
                                    ident::AbstractVector{Int}, omega::ComplexF64;
                                    omegaTolerance::Real = OMEGA_TOLERANCE,
                                    zeroTolerance::Real = ZERO_TOLERANCE)
    for ii in 1:6
        any(==(ii), ident) || continue
        pairPos = findfirst(==(ii), ident)
        if L.isIsotropic
            eSH = ComplexF64[0, 1, 0]
            eSV = cross(ComplexF64[1, 0, k3[ii]], eSH)
            if pairPos !== nothing && isodd(pairPos)
                U[:, ii] .= eSH
            else
                U[:, ii] .= eSV ./ norm(eSV)
            end
            continue
        end
        F = eigen(_gamma_matrix(L, k3[ii]), sortby = nothing)
        VV = F.vectors
        freq = sqrt.(abs.(F.values) ./ L.rho)
        degenLocation = zeros(Int, 3)
        count = 1
        for jj in 1:3
            if abs(freq[jj] - real(omega)) < omegaTolerance && count <= 2
                degenLocation[count] = jj
                count += 1
            end
        end
        pair2 = pairPos === nothing ? 0 :
                (isodd(pairPos) ? ident[pairPos + 1] : ident[pairPos - 1])
        partnerAssigned = pair2 > 0 && any(!iszero, @view U[:, pair2])
        if degenLocation[1] == 0 || degenLocation[2] == 0
            @warn "GreenSolver:calc:degenEig -- Cannot find correct degenerate eigenvalue for k3($ii); result may be inaccurate."
            if partnerAssigned
                best = argmin([abs(dot(VV[:, m], U[:, pair2])) for m in 1:3])
                U[:, ii] .= VV[:, best]
            else
                U[:, ii] .= VV[:, 1]
            end
            continue
        end
        if pairPos === nothing
            U[:, ii] .= VV[:, degenLocation[1]]
            continue
        end
        if partnerAssigned
            if abs(dot(U[:, pair2], VV[:, degenLocation[1]])) < zeroTolerance
                U[:, ii] .= VV[:, degenLocation[1]]
            elseif abs(dot(U[:, pair2], VV[:, degenLocation[2]])) < zeroTolerance
                U[:, ii] .= VV[:, degenLocation[2]]
            elseif U[:, pair2] != VV[:, degenLocation[1]]
                U[:, ii] .= VV[:, degenLocation[1]]
            else
                U[:, ii] .= VV[:, degenLocation[2]]
            end
        else
            U[:, ii] .= VV[:, degenLocation[1]]
        end
    end
    return U
end

# -----------------------------------------------------------------------------------------
# The per-point kernel
# -----------------------------------------------------------------------------------------

"""
    LayerPointBuffers

Scratch for [`_layer_point!`](@ref), one per layer per task: the sextic coefficients `c`,
the running root estimates `z` (`P = -i k3`, carried to the next velocity as the warm
start), the sorted `k3`, the degeneracy record `ident`, sort scratch, and `U` (3 x 6).
"""
struct LayerPointBuffers
    c::Vector{ComplexF64}
    z::Vector{ComplexF64}
    k3::Vector{ComplexF64}
    ident::Vector{Int}
    pos::Vector{Int}
    key::Vector{Float64}
    U::Matrix{ComplexF64}
end
LayerPointBuffers() = LayerPointBuffers(zeros(ComplexF64, 7), zeros(ComplexF64, 6),
                                        zeros(ComplexF64, 6), zeros(Int, 6), zeros(Int, 6),
                                        zeros(6), zeros(ComplexF64, 3, 6))

"""
    _layer_point!(b, L, w, omegaD, warm, need) -> Bool

One layer at one velocity: roots, canonical order, degeneracy record and the
polarisations of the waves selected by `need(k3, n)` (all six where the column has a
degenerate pair). `warm`: `b.z` holds the previous velocity's roots. Returns whether the
column is degenerate. Errors `GreenSolver:calc:rootCount` like `calc.m` if the sextic
degenerates (vanishing leading/trailing coefficient).
"""
function _layer_point!(b::LayerPointBuffers, L::GreenLayer, w::ComplexF64,
                       omegaD::ComplexF64, warm::Bool, need::F;
                       v::Real = NaN, iLay::Int = 1) where {F}
    c, z, k3 = b.c, b.z, b.k3
    _sextic_coeffs!(c, L.DetW, w)
    ok = warm && !L.isIsotropic && !iszero(c[1]) && !iszero(c[7]) && _aberth!(z, c)
    if !ok
        r = matlab_roots(c)
        length(r) == 6 || error(
            "GreenSolver:calc:rootCount -- Christoffel polynomial gave $(length(r)) roots ",
            "(expected 6) at v = $v km/s, layer $iLay -- leading/trailing coefficient vanished.")
        copyto!(z, r)
    end
    @inbounds for n in 1:6
        k3[n] = im * z[n]                       # P = -i k3
    end
    _canonical_sort!(k3)
    _ident_col!(b.ident, b.pos, b.key, k3, ZERO_TOLERANCE)
    degenerate = b.ident[1] != 0

    wr = L.rho * omegaD^2                       # velRangeD^2, exactly as calc.m
    U = b.U
    @inbounds for n in 1:6
        if degenerate && any(==(n), b.ident)
            U[1, n] = 0; U[2, n] = 0; U[3, n] = 0
        elseif degenerate || need(k3, n)
            u1, u2, u3 = _null_vector(L, k3[n], wr)
            U[1, n] = u1; U[2, n] = u2; U[3, n] = u3
        end
    end
    degenerate && _degenerate_polarisations!(U, L, k3, b.ident, omegaD)
    return degenerate
end

_need_all(k3, n) = true

# -----------------------------------------------------------------------------------------
# Whole-axis wrappers (diagnostics and MATLAB comparison of the pre-pass)
# -----------------------------------------------------------------------------------------

"""
    solve_k3_allvel(D, velRange, iLay = 1) -> Matrix{ComplexF64}        # 6 x nVel

Reference root-find, the literal port of `calc.m`'s `roots` loop: `k3 = i * matlab_roots`
of each coefficient column, raw (unsorted) order. Not used by `calc` (which uses the
warm-started Aberth iteration of `_layer_point!`); kept as the MATLAB-equivalent reference
that iteration is compared against.
"""
function solve_k3_allvel(D::AbstractMatrix{<:Complex}, velRange::AbstractVector{<:Real},
                          iLay::Integer = 1)
    nV = size(D, 2)
    k3M = Matrix{ComplexF64}(undef, 6, nV)
    @inbounds for iv in 1:nV
        r = matlab_roots(@view D[:, iv])
        length(r) == 6 || error(
            "GreenSolver:calc:rootCount -- Christoffel polynomial gave $(length(r)) roots ",
            "(expected 6) at v = $(velRange[iv]) km/s, layer $iLay -- leading/trailing ",
            "coefficient vanished.")
        for n in 1:6
            k3M[n, iv] = im * r[n]
        end
    end
    return k3M
end

"""
    canonicalize_k3(k3M) -> (k3sorted, perms)

Sort each column's six k3 lexicographically by `(real, imag)` (stable). **Deliberate
divergence from MATLAB**, which keeps the vendor-LAPACK-dependent `roots()` order: the
port's partial-wave labelling is canonical, every summed quantity (G, dispField) is
permutation-invariant. Smallest gap between canonically ordered roots on the cases
checked against MATLAB: 0.031, far above the 1e-14 cross-language value noise.
"""
function canonicalize_k3(k3M::AbstractMatrix{<:Complex})
    nV = size(k3M, 2)
    out = similar(k3M, ComplexF64)
    perms = Matrix{Int}(undef, size(k3M, 1), nV)
    @inbounds for iv in 1:nV
        p = sortperm(@view(k3M[:, iv]); by = z -> (real(z), imag(z)), alg = MergeSort)
        perms[:, iv] = p
        for n in eachindex(p)
            out[n, iv] = k3M[p[n], iv]
        end
    end
    return out, perms
end

"""
    check_identical_records_allvel(k3M, zeroTolerance = ZERO_TOLERANCE) -> Matrix{Int}

`calc.m` `CheckIdenticalRecords_allVel`: per column, the indices of up to two degenerate
k3 pairs (zeros elsewhere). Column-wise [`_ident_col!`](@ref).
"""
function check_identical_records_allvel(k3M::AbstractMatrix{<:Complex},
                                         zeroTolerance::Real = ZERO_TOLERANCE)
    size(k3M, 1) == 6 || throw(DimensionMismatch(
        "check_identical_records_allvel: k3M must have 6 rows, got $(size(k3M, 1))"))
    nV = size(k3M, 2)
    identM = zeros(Int, 6, nV)
    pos = zeros(Int, 6); key = zeros(6); col = zeros(ComplexF64, 6)
    for iv in 1:nV
        col .= @view k3M[:, iv]
        _ident_col!(@view(identM[:, iv]), pos, key, col, zeroTolerance)
    end
    return identM
end

"""
    calc_u_allvel(elTenPerm, rho, k1, k2, k3M, velRangeD, identM;
                  omegaTolerance, zeroTolerance, isIsotropic) -> Array{ComplexF64,3}

Polarisations of all six waves at every velocity (`calc.m` `Calc_U_allVel`), for a given
root matrix. Non-degenerate: null vector of `Gamma - rho*omegaD^2*I` (the ray MATLAB's
SVD gives); degenerate: [`_degenerate_polarisations!`](@ref). `elTenPerm` is the
`reshape(permutedims(elTen, (1,3,2,4)), 9, 9)` of `calc.m`; only `k1 = 1, k2 = 0` is
supported (the only branch `calc` solves).
"""
function calc_u_allvel(elTenPerm::AbstractMatrix{<:Real}, rho::Real, k1::Real, k2::Real,
                        k3M::AbstractMatrix{<:Complex}, velRangeD::AbstractVector{<:Complex},
                        identM::AbstractMatrix{<:Integer};
                        omegaTolerance::Real = OMEGA_TOLERANCE,
                        zeroTolerance::Real = ZERO_TOLERANCE,
                        isIsotropic::Bool = false)
    (k1 == 1 && k2 == 0) || throw(ArgumentError("calc_u_allvel: only k1 = 1, k2 = 0"))
    elTen = reshape(permutedims(reshape(elTenPerm, 3, 3, 3, 3), (1, 3, 2, 4)), 3, 3, 3, 3)
    L = GreenLayer(Array{Float64,4}(elTen), rho, isIsotropic)
    nV = size(k3M, 2)
    UAll = zeros(ComplexF64, 3, 6, nV)
    k3 = zeros(ComplexF64, 6)
    ident = zeros(Int, 6)
    for iv in 1:nV
        k3 .= @view k3M[:, iv]
        ident .= @view identM[:, iv]
        wr = L.rho * velRangeD[iv]^2
        U = @view UAll[:, :, iv]
        for n in 1:6
            any(==(n), ident) && continue
            u = _null_vector(L, k3[n], ComplexF64(wr))
            U[1, n] = u[1]; U[2, n] = u[2]; U[3, n] = u[3]
        end
        ident[1] != 0 && _degenerate_polarisations!(U, L, k3, ident, ComplexF64(velRangeD[iv]);
                                                    omegaTolerance, zeroTolerance)
    end
    return UAll
end

"""
    GreenSexticPrepass

The boundary-condition-independent per-point data of one `calc` call, over the whole axis,
in internal (surface-first) layer order: damped axis, per-layer canonical `k3` (`6 x nVel`),
polarisations `U` (`3 x 6 x nVel`, all six waves), degeneracy records `ident`, and the
substrate admissibility mask `adm = imag(k3) .> 0`. Built by [`green_sextic_prepass`](@ref)
with the same kernel `calc` uses; a diagnostic API, not on `calc`'s hot path.
"""
struct GreenSexticPrepass
    velRange::Vector{Float64}
    velRangeD::Vector{ComplexF64}
    eta::Float64
    elTen::Vector{Array{Float64,4}}
    density::Vector{Float64}
    isIsotropic::Vector{Bool}
    k3::Vector{Matrix{ComplexF64}}
    U::Vector{Array{ComplexF64,3}}
    ident::Vector{Matrix{Int}}
    adm::BitMatrix
end

nvel(p::GreenSexticPrepass) = length(p.velRange)
nlayers(p::GreenSexticPrepass) = length(p.density)

"""
    green_sextic_prepass(structure, angle, velRange; eta = LIMITING_ABSORPTION_ETA)
        -> GreenSexticPrepass

Run the per-point kernel of `calc` over the whole axis for every layer and keep its
k3 / U / degeneracy / admissibility output (`k1 = +1` only, as `calc.m`).
"""
function green_sextic_prepass(structure::Structure, angle::Real,
                               velRange::AbstractVector{<:Real};
                               eta::Real = LIMITING_ABSORPTION_ETA)
    (all(isfinite, velRange) && all(>(0), velRange)) || throw(ArgumentError(
        "GreenSolver.calc: velRange must be finite and positive"))
    elTen, density, isIso = green_layer_setup(structure, angle)
    layers = [GreenLayer(elTen[i], density[i], isIso[i]) for i in eachindex(density)]
    vel = Float64.(collect(velRange))
    nV = length(vel)
    nLay = length(layers)
    wAllD, velRangeD = _damped_axis(vel, eta)
    k3 = [zeros(ComplexF64, 6, nV) for _ in 1:nLay]
    U = [zeros(ComplexF64, 3, 6, nV) for _ in 1:nLay]
    ident = [zeros(Int, 6, nV) for _ in 1:nLay]
    for iLay in 1:nLay
        b = LayerPointBuffers()
        for iv in 1:nV
            _layer_point!(b, layers[iLay], wAllD[iv], velRangeD[iv], iv > 1, _need_all;
                          v = vel[iv], iLay = iLay)
            k3[iLay][:, iv] .= b.k3
            U[iLay][:, :, iv] .= b.U
            ident[iLay][:, iv] .= b.ident
        end
    end
    adm = imag.(k3[nLay]) .> 0
    return GreenSexticPrepass(vel, velRangeD, Float64(eta), elTen, density, isIso,
                              k3, U, ident, BitMatrix(adm))
end

"""
    substrate_wave_indices(k3col) -> (inds::NTuple{3,Int}, usedFallback)

The three admissible substrate waves (`imag(k3) > 0`, exactly three by the 3/3
split of the damped sextic, Every et al. 1997, pp. 1346-1347), as ascending indices;
otherwise the three with the largest `imag(k3)` (stable descending sort, as `calc.m`)
and `usedFallback = true`.
"""
function substrate_wave_indices(k3col::AbstractVector{ComplexF64})
    nadm = 0
    i1 = i2 = i3 = 0
    @inbounds for n in 1:6
        if imag(k3col[n]) > 0
            nadm += 1
            nadm == 1 ? (i1 = n) : nadm == 2 ? (i2 = n) : nadm == 3 ? (i3 = n) : nothing
        end
    end
    nadm == 3 && return (i1, i2, i3), false
    ord = sortperm(imag.(k3col); rev = true, alg = MergeSort)
    s = sort(ord[1:3])
    return (s[1], s[2], s[3]), true
end

"Prepass form: column `iv` of the substrate's k3."
substrate_wave_indices(p::GreenSexticPrepass, iv::Integer) =
    substrate_wave_indices(p.k3[nlayers(p)][:, iv])

export GreenSexticPrepass, green_sextic_prepass, green_layer_setup, green_layers,
       GreenLayer, calc_detw_by_christoffel, solve_k3_allvel, canonicalize_k3,
       check_identical_records_allvel,
       calc_u_allvel, substrate_wave_indices, kvec_cols, outer_cols,
       OMEGA_TOLERANCE, ZERO_TOLERANCE, LIMITING_ABSORPTION_ETA,
       NEAR_ISO_ANGLE_STEP, NEAR_ISO_WARN_THRESH
