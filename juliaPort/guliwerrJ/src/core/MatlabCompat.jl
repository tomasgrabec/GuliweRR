# src/core/MatlabCompat.jl
#
# Six MATLAB-builtin compatibility shims used by the solvers: `matlab_roots`,
# `movmedian`, `matchpairs0`, `rcond_matlab`, `fminbnd_brent`, `findpeaks_matlab`.
#
# Every recipe here was empirically verified against MATLAB R2025b — the approach, the exact
# tolerances, and the documented gotchas (companion-matrix element type, `strict = false`, the
# `'1'` `gecon!` norm type, …) are recorded in the comments below. Do not re-derive these
# from scratch; if a comparison against MATLAB disagrees with the tolerance measured here,
# re-check the fixture before loosening the tolerance.
#
# Plain top-level definitions, not `module MatlabCompat` — `CoreTypes.jl` `include`s this
# file directly, one module per top-level folder rather than one per file (a
# submodule named after its own main export collides with that export's name once
# `using`d from the same scope).

using LinearAlgebra
using Statistics: median
using Peaks: argmaxima, peakproms, peakwidths
using Optim: Optim, optimize, Brent
using Hungarian: hungarian

export matlab_roots, movmedian, matchpairs0, rcond_matlab, fminbnd_brent, findpeaks_matlab

# ---------------------------------------------------------------------------------------
# 1. matlab_roots
# ---------------------------------------------------------------------------------------

"""
    matlab_roots(c::AbstractVector{<:Number}) -> Vector{<:Complex}

Port of MATLAB's `roots(c)`: eigenvalues of the coefficient vector `c`'s companion matrix,
in LAPACK's raw (unsorted) order, with one zero root prepended per stripped trailing
coefficient.

**Type preservation is the whole finding**: the companion
matrix is built with `eltype(cc)`, NOT forced to `ComplexF64`. Forcing complex on real
coefficients routes LAPACK through `zgeev` instead of `dgeev` and reorders conjugate root
pairs (measured: same-order error goes from 5.6e-16 to 1.30 on a real fixture). `eigvals`
is called with `sortby = nothing` because Julia's `eigvals` sorts by default whereas
MATLAB's `eig`/`roots` do not, and `GreenSolver.calc`'s degeneracy/partial-wave
logic depends on the exact k3 ordering this produces.

Residuals of 1e-5 to 1e-8 on exactly- or near-degenerate roots are root *conditioning*
(`eps^(1/m)` for an m-fold root) — both MATLAB and Julia are equally imprecise there, and
it is safely below `GreenSolver.calc`'s own `ZERO_TOLERANCE = 1e-4` degeneracy gap
(and the limiting-absorption damping lifts exact degeneracies anyway), so a degeneracy
verdict cannot flip.
"""
function matlab_roots(c::AbstractVector{T}) where {T<:Number}
    nz = findall(!iszero, c)
    isempty(nz) && return Complex{real(float(T))}[]
    ntrail = length(c) - nz[end]
    cc = float.(c[nz[1]:nz[end]])
    n = length(cc)
    RT = Complex{real(eltype(cc))}
    n == 1 && return zeros(RT, ntrail)
    A = zeros(eltype(cc), n - 1, n - 1)   # NOTE: eltype(cc), NOT ComplexF64
    for i in 2:(n-1)
        A[i, i-1] = one(eltype(cc))
    end
    A[1, :] .= -cc[2:n] ./ cc[1]
    return vcat(zeros(RT, ntrail), RT.(eigvals(A, sortby = nothing)))
end

# ---------------------------------------------------------------------------------------
# 2. movmedian
# ---------------------------------------------------------------------------------------

"""
    movmedian(v::AbstractVector, w::Integer) -> Vector

Port of MATLAB's `movmedian(v, w)`: a centred moving median that **shrinks** at the ends
(it never pads or wraps). For even `w` the window spans `[i - w/2, i + w/2 - 1]`, matching
MATLAB's convention exactly (verified `max|Δ| = 0.0` for `w = 5` and `w = 11` against a
50-point signal with a step change). `RollingFunctions.jl` is
not used — its default endpoint handling differs, and this is six lines.
"""
function movmedian(v::AbstractVector, w::Integer)
    n = length(v)
    out = similar(float.(v))
    back = w ÷ 2
    fwd = iseven(w) ? w ÷ 2 - 1 : w ÷ 2
    for i in 1:n
        out[i] = median(@view v[max(1, i - back):min(n, i + fwd)])
    end
    return out
end

# ---------------------------------------------------------------------------------------
# 3. matchpairs0
# ---------------------------------------------------------------------------------------

"""
    matchpairs0(C::AbstractMatrix{<:Real}) -> Vector{Tuple{Int,Int}}

Port of MATLAB's `matchpairs(-subScore, 0)` as used by `ModeAssociation.m`:
rectangular linear assignment with cost of non-assignment `0`. With `costUnmatched = 0`
the objective is `Σ_matched C(i,j)`, so any pair with `C ≥ 0` never lowers the total and is
never taken — the reduction below (discard non-negative entries, run Hungarian on what is
left) is only valid for exactly this `costUnmatched = 0` case; a future caller passing a
non-zero value needs the general `(m+n)×(n+m)` padding reduction instead.

Verified exactly (including the degenerate all-zero case, where MATLAB returns zero pairs)
against 4 MATLAB fixtures: a 3×3 with one weak row, a 3×2 rectangular case, a near-tie, and
an all-zero matrix.
"""
function matchpairs0(C::AbstractMatrix{<:Real})
    W = Matrix{Union{Missing,Float64}}(missing, size(C))
    for i in eachindex(C)
        C[i] < 0 && (W[i] = Float64(C[i]))
    end
    all(ismissing, W) && return Tuple{Int,Int}[]
    asg, _ = hungarian(W)
    return [(i, asg[i]) for i in eachindex(asg) if asg[i] != 0]
end

# ---------------------------------------------------------------------------------------
# 4. rcond_matlab
# ---------------------------------------------------------------------------------------

"""
    rcond_matlab(A::AbstractMatrix) -> Float64

Port of MATLAB's `rcond(A)`: the reciprocal 1-norm condition number estimate, via LAPACK's
`gecon`. The `'1'` norm-type spelling is required — Julia's `gecon!` rejects LAPACK's `'O'`
spelling that the raw Fortran call accepts.

Measured against 4 MATLAB fixtures: exact agreement on `randn(5)`, `eye(4)*3`, and an
*exactly* singular 2×2 (`[1 2; 2 4+eps]` — `4+eps` rounds straight back to `4.0` in IEEE
double precision, since `eps(1.0) ≈ 2.22e-16` is smaller than `eps(4.0) ≈ 8.88e-16`; MATLAB
and Julia agree that this matrix is exactly rank-1); 1.64e-8 relative difference on
`hilb(8)`, which only appears on an input already far past any useful conditioning
threshold and does not matter for a warning gate (the 5 call sites in `GreenSolver.calc`
use `rcond` only to decide whether to warn).

**Deviation from the naive `lu(A)` + `gecon!` recipe: MATLAB's `rcond` returns `0.0` gracefully on an
exactly singular matrix, but Julia's `lu(A)` — checked, the default — *throws*
`SingularException` on the same input.** `lu(A; check = false)` plus an `issuccess` guard
reproduces MATLAB's graceful-zero behaviour; every GULIWERR call site only ever compares
`rcond` against a threshold to decide whether to warn, so returning exactly `0.0` here is
the right answer, not a special case to special-case away.

**Second deviation**: on a matrix with a non-finite entry
(reachable in practice — `GreenSolver.calc`'s thick-film/near-singular path can produce one
after equilibration), `opnorm(A, 1)` is itself non-finite, and Julia's `LAPACK.gecon!`
*throws* `ArgumentError("invalid argument #5 to LAPACK call")` on a non-finite `ANORM`.
MATLAB's `rcond` instead returns `NaN` gracefully. Guarded below the same way as the
exact-singularity case, for the same reason: every call site only ever compares against a
threshold to decide whether to warn, and MATLAB's own `min` (unlike Julia's) already ignores
`NaN`, so callers written against MATLAB's semantics expect a `NaN` here to fall out of a
`min`/threshold comparison silently, not to propagate or throw.
"""
function rcond_matlab(A::AbstractMatrix)
    isfinite(opnorm(A, 1)) || return NaN
    F = lu(A; check = false)
    issuccess(F) || return 0.0
    return LAPACK.gecon!('1', F.factors, opnorm(A, 1))
end

# ---------------------------------------------------------------------------------------
# 5. fminbnd_brent
# ---------------------------------------------------------------------------------------

"""
    fminbnd_brent(f, lo::Real, hi::Real; abs_tol::Real = 1e-12) -> (x_min, f_min)

Port of MATLAB's `fminbnd(f, lo, hi, optimset('TolX', tol))`, as used by
`LimitingVels.refineRoot`. Both `fminbnd` and `Optim.Brent()` implement Brent's method
(golden-section + parabolic interpolation) and converge to the **same minimum**, but their
internal bookkeeping differs slightly, so they do **not** take the same iterate path
(measured: `Δx` from 0 to 1.8e-9, `Δf` from 0 to 1.6e-10 across 3 test objectives).

**Consequence for every caller**: refined limiting velocities carry `rtol = 1e-8`, not the
tighter tolerance the Ritz eigenproblem gets. Anything that *thresholds* on a limiting
velocity should assert the resulting mask/label, which is stable, rather than the boundary
value itself.
"""
function fminbnd_brent(f, lo::Real, hi::Real; abs_tol::Real = 1e-12)
    res = optimize(f, lo, hi, Brent(); abs_tol = abs_tol)
    return Optim.minimizer(res), Optim.minimum(res)
end

# ---------------------------------------------------------------------------------------
# 6. findpeaks_matlab
# ---------------------------------------------------------------------------------------

"""
    findpeaks_matlab(y, x; minheight = -Inf, relheight = 0.5) -> NamedTuple

Port of MATLAB's
`findpeaks(y, x, 'MinPeakHeight', minheight, 'WidthReference', 'halfprom')`, as used by
`ModeClassification.m` (mode classification reads directly off peak
location/height/width/prominence, so an approximate match would be a physics change, not
merely a numerical one).

Returns `(loc, height, width, prominence)`, each a vector over the surviving peaks.

Three gotchas are handled by this exact recipe:
- `strict = false` is required on **all three** Peaks.jl calls. MATLAB's `findpeaks`
  treats the array ends as valid minima when computing prominence/width; Peaks.jl's
  default `strict = true` instead returns `prom = 0` / `width = NaN` for a peak whose
  reference interval is truncated by an edge. The three `strict` flags must all match, or
  `peakwidths` throws.
- MATLAB's `'WidthReference','halfprom'` is Peaks.jl's `relheight = 0.5` (the default
  here, kept as a keyword since some MATLAB callers could in principle pass another
  reference level).
- MATLAB's `findpeaks` never reports the first or last sample of the array as a peak — a
  peak needs a *strictly lower* neighbour on both sides, and a boundary sample (or a
  plateau whose flat run touches either array edge) cannot demonstrate a neighbour past the
  edge. `strict = false` on `argmaxima` fixes the *prominence/width reference-window*
  truncation described above, but it does **not** itself exclude a maximum whose plateau
  touches the array boundary — confirmed with a fixture whose last two samples are a tied
  maximum plateau: `argmaxima(y, 1; strict = false)` flags it as a 6th peak, MATLAB's
  `findpeaks` on the identical input does not. The explicit boundary walk below reproduces
  MATLAB's exclusion; skipping it would let a spurious trailing-edge "peak" through into
  `ModeClassification`'s mode identification.

Verified against a 17-point synthetic signal with a flat plateau, a value tie between two
separate peaks, and a plateau at each array edge (4 real peaks) and
against the real Cu(001) 22° `|G_eff|` curve (800 points).
"""
function findpeaks_matlab(y, x; minheight = -Inf, relheight = 0.5)
    n = length(y)
    idx0 = argmaxima(y, 1; strict = false)   # plateaus allowed

    # Boundary exclusion (see docstring):
    # MATLAB never treats the first/last sample as a peak. `strict = false` above lets a
    # plateau be reported by its first index even when the plateau's *far* end touches an
    # array edge, so walk each candidate's flat run explicitly and require a strictly
    # lower sample immediately outside it on BOTH sides.
    valid = falses(length(idx0))
    for (k, i) in pairs(idx0)
        v = y[i]
        okleft = i > 1 && y[i-1] < v
        j = i
        while j < n && y[j+1] == v
            j += 1
        end
        okright = j < n && y[j+1] < v
        valid[k] = okleft && okright
    end
    idx0 = idx0[valid]

    idx, proms = peakproms(idx0, y; strict = false)               # ends count as minima
    idx, w, lo, hi = peakwidths(idx, y, proms; relheight, strict = false)
    dx = x[2] - x[1]
    keep = findall(k -> y[idx[k]] >= minheight, eachindex(idx))
    return (loc = x[idx[keep]], height = y[idx[keep]],
            width = (hi[keep] .- lo[keep]) .* dx, prominence = proms[keep])
end
