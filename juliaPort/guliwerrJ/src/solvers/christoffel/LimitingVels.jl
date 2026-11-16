# =========================================================================================
# LimitingVels.jl -- `include`d directly into `module Christoffel` (see Christoffel.jl's
# header for the "one module per folder" rationale).
#
# Port of guliwerr/solvers/christoffel/LimitingVels.m: limiting bulk-wave velocities -- of
# the bulk waves with wavevector in the sagittal plane whose group velocity there has no
# surface-normal component (vg2 = 0): the energy flux runs parallel to the surface, and the
# limiting velocity is its along-surface component vg1, equal to the phase velocity along the
# surface, v/cos(theta*). The lowest one bounds the existence domain of SAW. `findRoots` (limiting velocities for one pre-rotated tensor;
# the rotation to the standard frame is the caller's, e.g. CalcCommon.limitingAndBulk) and
# `refineRoot` (re-find ONE already-located root in a small bracket).
#
# Behaviour worth knowing (as in MATLAB):
#   - findRoots' fminbnd runs at TolX 1e-12 (MATLAB's default is 1e-4): at a steep
#     slow-shear crossing the default leaves a residual above the accept test and rejects a
#     real root (synthetic cubic of Zener ratio 20, (112) cut, 118 deg). Both methods
#     share one refine-and-accept helper (`refineAt`).
#   - input validation at both entries is raised as ArgumentErrors with MATLAB ids.
#   - performance: the 9x9 permuted tensor is hoisted once per call (`_kronperm`), and the
#     coarse scan's eigenvectors are not recomputed.
#
# LimitingRoot STRUCT: `findRoots` in MATLAB returns three PARALLEL cell arrays --
# `vLim_sA{iwt}`/`eigenVecs_sA{iwt}`/`thetaRoots_sA{iwt}` -- that must be kept in step by
# convention alone (nothing stops a caller from indexing `vLim_sA{2}(3)` against
# `eigenVecs_sA{2}(5,:)` by a typo). Bundling `(iwt, vel, theta, pol)` into one struct per
# accepted root removes that whole hazard structurally: a `Vector{LimitingRoot}` can never
# have its three "columns" drift out of sync, because there is only one column. Storing
# `iwt` ON the struct (rather than leaving it implicit in "which of the 3 output vectors
# this came from") gives the "branch label" MATLAB has to reconstruct by hand downstream
# (`CalcCommon.limitingAndBulk`'s `shearLimIwt`, `ModeAssociation.limitingRootScore`) --
# with it on the struct, flattening the fast+slow shear branches together
# (`vcat(roots[2], roots[3])`) carries its own label for free.
#
# =========================================================================================

"""
    LimitingRoot

One accepted sagittal-plane limiting-velocity root. See this file's header for why
`findRoots` bundles what MATLAB keeps as three parallel cell arrays into one struct per
root.

Fields:
- `iwt::Int`             -- wave-type branch, 1 = longitudinal, 2 = fast shear,
                            3 = slow shear (Christoffel's descending-eigenvalue convention)
- `vel::Float64`         -- limiting velocity `vg1` at the accepted root
- `theta::Float64`       -- sagittal-plane angle theta* [rad] of the root
- `pol::Vector{Float64}` -- polarisation eigenvector at theta*
"""
struct LimitingRoot
    iwt::Int
    vel::Float64
    theta::Float64
    pol::Vector{Float64}
end

# -----------------------------------------------------------------------------------------
# islocalmin -- MATLAB builtin, restricted to what findRoots needs (real vector, default
# options). NOT one of MatlabCompat.jl's general-purpose shims -- used only here, so it
# stays local.
# -----------------------------------------------------------------------------------------

"""
    islocalmin(v) -> BitVector

Minimal port of MATLAB's `islocalmin(v)` for a real vector with default options (no
`MinProminence`, no `SamplePoints`).

Verified against MATLAB R2025b: MATLAB's default `islocalmin` never flags an endpoint
(`v[1]`/`v[end]` can never be `true`), and for a maximal flat run `v[i0] == v[i0+1] == ... == v[i1]`
sandwiched between two STRICTLY larger neighbours, it flags exactly the run's own midpoint
(`(i0+i1) ÷ 2`, rounding DOWN for an even-length run -- e.g. a run at positions 6:7 is
flagged at 6, not 7) rather than every point in the run.

On data with no exact ties -- the generic case for a smooth Christoffel `|vg2|` scan -- this
reduces to the textbook "strictly lower than both neighbours" test (confirmed identical
index sets on a real Cu(001) case). The flat-run handling is included anyway: near-ties are
a real possibility for a near-isotropic material, and getting this wrong would silently
change the root COUNT `findRoots` returns.
"""
function islocalmin(v::AbstractVector{<:Real})
    n = length(v)
    tf = falses(n)
    n < 3 && return tf
    i = 1
    while i <= n
        j = i
        while j < n && v[j + 1] == v[i]
            j += 1
        end
        if i > 1 && j < n && v[i - 1] > v[i] && v[j + 1] > v[j]
            tf[(i + j) ÷ 2] = true
        end
        i = j + 1
    end
    return tf
end


# -----------------------------------------------------------------------------------------
# Private helpers -- MATLAB `methods (Static, Access = private)`
# -----------------------------------------------------------------------------------------

"""fminbnd tolerance of every root refinement (MATLAB `LimitingVels.FMIN_OPTS`, TolX 1e-12)."""
const FMIN_TOLX = 1e-12

"""Domain of the sagittal-plane scan, also the clamp of `refineRoot`'s bracket."""
const THETA_LO = -pi / 2 + 1e-4
const THETA_HI =  pi / 2 - 1e-4

"""
    DEGENERATE_REL_TOL, DEGENERATE_THETA_TOL

Two roots of DIFFERENT wave types are one degenerate root (an acoustic axis in the sagittal
plane, e.g. cubic (001) at 0 deg, or any direction of an isotropic medium) when their
velocities agree to `DEGENERATE_REL_TOL` (relative) and their theta* to
`DEGENERATE_THETA_TOL` [rad]. Set from the root finder's own precision: Brent/fminbnd stop
at `sqrt(eps)|theta|` ~ 1.5e-8 rad, which moves the velocity by up to ~3e-8 relative at the
steepest crossings seen; 1e-6 is 30x above that and far below any split the eigensolver
resolves (polarisations defined to ~eps/gap ~ 1e-10 at a 1e-6 gap). MATLAB
`LimitingVels.DEGENERATE_REL_TOL` / `DEGENERATE_THETA_TOL`.
"""
const DEGENERATE_REL_TOL = 1e-6
const DEGENERATE_THETA_TOL = 1e-6

"""
    canonicalBasis(E) -> Matrix{Float64}

Basis-independent orthonormal basis of the span of the orthonormal columns `E` (MATLAB
`LimitingVels.canonicalBasis`): from the projector `P = E E'` (the same matrix for ANY
orthonormal basis of the subspace, to round-off) take `P`'s longest column (ties within
1e-10 relative -> the lowest index), normalise it, deflate `P` by it, repeat. Its own
entry `P[c, c] > 0` fixes the sign. Cu(001) at 0 deg (eigenspace y-z) gives (0,1,0), then
(0,0,1).
"""
function canonicalBasis(E::AbstractMatrix{<:Real})
    P = E * E'
    P = (P + P') / 2
    g = size(E, 2)
    B = zeros(3, g)
    for k in 1:g
        nrm = [sqrt(sum(abs2, view(P, :, j))) for j in 1:3]
        c = findfirst(>=(maximum(nrm) * (1 - 1e-10)), nrm)
        b = P[:, c] ./ nrm[c]
        B[:, k] = b
        P = P - b * b'
    end
    return B
end

"""
    canonicaliseDegenerate!(M9, rootsByType) -> rootsByType

Give the roots of a degenerate group a canonical basis of their eigenspace, in place (MATLAB
`LimitingVels.canonicaliseDegenerate`). A group: roots of different wave types at one
velocity and one theta* (`DEGENERATE_REL_TOL` / `DEGENERATE_THETA_TOL`) whose phase
velocities coincide there too. Their eigenvalue is multiple, so each root's own eigenvector
is an arbitrary vector of the eigenspace, decided by round-off (MATLAB returned (0,1,0) for
BOTH shear roots of Cu(001) at 0 deg; MATLAB and Julia disagree); the group's polarisations
are replaced by [`canonicalBasis`](@ref) of the eigenspace at the first member's theta*, in
wave-type order. Non-degenerate roots are untouched.
"""
function canonicaliseDegenerate!(M9::AbstractMatrix{Float64}, rootsByType)
    idx = [(iwt, q) for iwt in 1:3 for q in eachindex(rootsByType[iwt])]
    root(k) = rootsByType[idx[k][1]][idx[k][2]]
    done = falses(length(idx))
    for a in eachindex(idx)
        done[a] && continue
        ra = root(a)
        grp = [a]
        for b in a+1:length(idx)
            rb = root(b)
            if !done[b] && !any(root(m).iwt == rb.iwt for m in grp) &&
               abs(rb.vel - ra.vel) <= DEGENERATE_REL_TOL * abs(ra.vel) &&
               abs(rb.theta - ra.theta) <= DEGENERATE_THETA_TOL
                push!(grp, b)
            end
        end
        length(grp) < 2 && continue
        _, vPh, _, E = vg2All(M9, ra.theta)
        iw = [root(m).iwt for m in grp]
        maximum(vPh[iw]) - minimum(vPh[iw]) > DEGENERATE_REL_TOL * maximum(vPh[iw]) && continue
        done[grp] .= true
        B = canonicalBasis(E[:, iw])
        for (k, m) in enumerate(grp)
            root(m).pol .= view(B, :, k)
        end
    end
    return rootsByType
end

"""
    checkTensorRho(rotElTensor, rho, methodName)

A real, finite 3x3x3x3 tensor and a positive finite scalar density (MATLAB
`LimitingVels.checkTensorRho`; `validateattributes` errors become `ArgumentError`s carrying
the MATLAB method name).
"""
function checkTensorRho(rotElTensor::AbstractArray{<:Real,4}, rho::Real, methodName::String)
    (size(rotElTensor) == (3, 3, 3, 3) && all(isfinite, rotElTensor)) ||
        throw(ArgumentError("LimitingVels.$methodName: rotElTensor must be a finite 3x3x3x3 tensor."))
    (isfinite(rho) && rho > 0) ||
        throw(ArgumentError("LimitingVels.$methodName: rho must be a positive finite scalar."))
    return nothing
end

"""
    vg2All(M9, thetaRad) -> (vg2, vels, dVels, eigVecs)

Out-of-plane (x3) group-velocity component for all three wave types at one sagittal-plane
angle, `vg2(theta) = v sin(theta) + dv/dtheta cos(theta)`, on a hoisted density-normalised
[`_kronperm`](@ref) matrix `M9`. Also returns what it was built from, so a caller that needs
`vg1` or the polarisation at the same angle does not re-solve the eigenproblem.
"""
function vg2All(M9::AbstractMatrix{Float64}, thetaRad::Real)
    v, dv, eigVecs = _pvad_perm(M9, thetaRad, 1.0)
    s, c = sin(thetaRad), cos(thetaRad)
    return v .* s .+ dv .* c, v, dv, eigVecs
end

"""fminbnd objective of `refineAt`: `|vg2|` for one wave type."""
vg2Abs(M9::AbstractMatrix{Float64}, thetaRad::Real, iwt::Integer) =
    abs(vg2All(M9, thetaRad)[1][iwt])

"""
    refineAt(M9, lo, hi, iwt, relTol) -> (vLim, thetaRoot, residual, ok, pol)

Refine one root of `|vg2|` in `[lo, hi]` (Brent, `abs_tol = FMIN_TOLX`) and apply the accept
test, shared by [`findRoots`](@ref) and [`refineRoot`](@ref) so the two cannot drift apart
(MATLAB `LimitingVels.refineAt`). `vLim = vg1 = v cos(theta) - dv/dtheta sin(theta)` at the
root; accepted when `|vg2| < relTol |vg1|` and `vg1 > 0`. `pol` is the polarisation at the
root (MATLAB's `polAt`, here taken from the same eigendecomposition as `vg1`).
"""
function refineAt(M9::AbstractMatrix{Float64}, lo::Real, hi::Real, iwt::Integer, relTol::Real)
    thetaRoot, residual = fminbnd_brent(th -> vg2Abs(M9, th, iwt), lo, hi; abs_tol = FMIN_TOLX)
    _, vRoot, dvRoot, eigVecs = vg2All(M9, thetaRoot)
    vLim = vRoot[iwt] * cos(thetaRoot) - dvRoot[iwt] * sin(thetaRoot)
    ok = residual < relTol * abs(vLim) && vLim > 0
    return vLim, thetaRoot, residual, ok, eigVecs[:, iwt]
end

# -----------------------------------------------------------------------------------------
# findRoots
# -----------------------------------------------------------------------------------------

"""
    findRoots(rotElTensor, rho=1.0, alphaPrecision=1e-2, relTol=1e-3) -> NTuple{3,Vector{LimitingRoot}}

Limiting velocities for a single pre-rotated elastic tensor (rotated to the standard frame
by the caller, e.g. `CalcCommon.limitingAndBulk`).

INPUT: `rotElTensor` -- 3x3x3x3 rotated elastic tensor; `rho` -- density (default `1`);
`alphaPrecision` -- sagittal-plane scan step [rad] (default `1e-2`, positive, `< 1`);
`relTol` -- accept threshold, relative to the in-plane group-velocity magnitude at the
candidate root (default `1e-3`).
OUTPUT: a 3-tuple, indexed by wave type (`1` = longitudinal, `2` = fast shear, `3` = slow
shear per Christoffel's descending-eigenvalue convention), each a `Vector{LimitingRoot}` of
the accepted roots for that branch, in scan order. Each root carries its `theta` (MATLAB's
third output `thetaRoots_sA`), so a caller can anchor it and re-find it with
[`refineRoot`](@ref).

Coarse `islocalmin` scan over `|vg2|`, then a bounded Brent refinement of every candidate
at TolX 1e-12 before the accept/reject test. `islocalmin` (rather than a sign change) also
catches tangential touches -- bulk waves grazing near-zero out-of-plane group velocity
without crossing sign, a real surface-skimming feature. Algorithm:
`docs/reference/architecture/christoffel.md`.

Why TolX 1e-12, not fminbnd's default 1e-4 rad: at a steep sign crossing the residual
`|vg2|` left by a theta error is slope x error. A strongly anisotropic cubic can cross at
`|dvg2/dtheta|` of tens of km/s/rad, where the default leaves a residual above the accept
tolerance `relTol vg1` and a genuine root is rejected (synthetic cubic of Zener ratio 20,
(112) cut, 118 deg). The theta error would also move the accepted velocity, as `vg1` is not stationary
at the root.

**Downstream, `min(shearLimVels)` taken across ALL entries is the correct Rayleigh existence
boundary, not a branch-mixing bug.** This function itself makes no such reduction; the note
is for callers.

Port of `LimitingVels.findRoots`.
"""
function findRoots(rotElTensor::AbstractArray{<:Real,4}, rho::Real = 1.0,
                    alphaPrecision::Real = 1e-2, relTol::Real = 1e-3)
    checkTensorRho(rotElTensor, rho, "findRoots")
    (isfinite(alphaPrecision) && 0 < alphaPrecision < 1) ||
        throw(ArgumentError("LimitingVels.findRoots: alphaPrecision must be a positive scalar < 1."))
    (isfinite(relTol) && relTol > 0) ||
        throw(ArgumentError("LimitingVels.findRoots: relTol must be a positive scalar."))

    M9 = _kronperm(rotElTensor ./ rho)

    # coarse scan, radians throughout (phaseVelAngleDeriv is radian-native)
    theta = collect(THETA_LO:alphaPrecision:THETA_HI)
    nTheta = length(theta)
    absVg2 = Matrix{Float64}(undef, nTheta, 3)       # column per wave type: contiguous rows
    for it in 1:nTheta
        vg2 = vg2All(M9, theta[it])[1]
        @inbounds for iwt in 1:3
            absVg2[it, iwt] = abs(vg2[iwt])
        end
    end

    rootsByType = ntuple(_ -> LimitingRoot[], 3)
    for iwt in 1:3
        for ic in findall(islocalmin(@view absVg2[:, iwt]))
            lo = theta[max(ic - 1, 1)]
            hi = theta[min(ic + 1, nTheta)]
            vg1Root, thetaRoot, _, ok, pol = refineAt(M9, lo, hi, iwt, relTol)
            ok && push!(rootsByType[iwt], LimitingRoot(iwt, vg1Root, thetaRoot, pol))
        end
    end
    return canonicaliseDegenerate!(M9, rootsByType)
end

# -----------------------------------------------------------------------------------------
# refineRoot
# -----------------------------------------------------------------------------------------

"""
    refineRoot(rotElTensor, rho, thetaGuess, iwt, halfWidth=0.05, relTol=1e-3)
        -> (vLim, thetaRoot, pol, residual, ok)

Re-find ONE already-located sagittal-plane root in a small bracket.

INPUT: `rotElTensor` -- 3x3x3x3 rotated elastic tensor; `rho` -- density; `thetaGuess` --
sagittal-plane angle [rad] of the root being tracked; `iwt` -- wave-type branch index
(1 = longitudinal, 2 = fast shear, 3 = slow shear); `halfWidth` -- bracket half-width [rad]
(default `0.05`, five steps of `findRoots`' coarse grid); `relTol` -- accept threshold
(default `1e-3`, the same test `findRoots` applies).

OUTPUT: `vLim` -- limiting velocity at the refined root; `thetaRoot` -- refined theta*
[rad]; `pol` -- `[3]` polarisation at the refined root; `residual` -- `|vg2|` at the refined
root; `ok` -- accept test passed. `false` is the root-DISAPPEARANCE detector: a limiting
root can vanish outright under a small change in the elastic constants, and then the
minimum of `|vg2|` inside the bracket stops being near zero. **A caller MUST check it and
fall back to a full `findRoots` rescan; `vLim` is meaningless when `false`.**

Use this, not `findRoots`, inside a finite-difference derivative: `findRoots` brackets each
root from whichever coarse-grid cell `islocalmin` flagged, which makes the returned
velocity a non-smooth function of the elastic constants at the 1e-4 GPa scale a central
difference perturbs them by (a tiny change can move the flagged cell by one step and jump
the bracket). Here the bracket is fixed by the CALLER (identical for the `+h` and `-h`
evaluations), which makes the value smooth and the difference valid.

Freezing theta* instead and differentiating only the explicit `c`-dependence is NOT a valid
shortcut: `v_lim = v_g1(theta*(c), c)` and `v_g1` is not stationary in theta at the root
(only `v_g2` is), so the `dtheta*/dc` term is not small.

`pol` is always returned (MATLAB gates it behind `nargout > 2`; here it comes from the
eigendecomposition `vLim` already needs, so it costs nothing).

Port of `LimitingVels.refineRoot`.
"""
function refineRoot(rotElTensor::AbstractArray{<:Real,4}, rho::Real, thetaGuess::Real,
                     iwt::Integer, halfWidth::Real = 0.05, relTol::Real = 1e-3)
    checkTensorRho(rotElTensor, rho, "refineRoot")
    isfinite(thetaGuess) ||
        throw(ArgumentError("LimitingVels.refineRoot: thetaGuess must be a finite scalar."))
    1 <= iwt <= 3 ||
        throw(ArgumentError("LimitingVels.refineRoot: iwt must be 1, 2 or 3."))
    (isfinite(halfWidth) && halfWidth > 0) ||
        throw(ArgumentError("LimitingVels.refineRoot: halfWidth must be a positive scalar."))
    (isfinite(relTol) && relTol > 0) ||
        throw(ArgumentError("LimitingVels.refineRoot: relTol must be a positive scalar."))

    M9 = _kronperm(rotElTensor ./ rho)
    lo = max(thetaGuess - halfWidth, THETA_LO)
    hi = min(thetaGuess + halfWidth, THETA_HI)
    vLim, thetaRoot, residual, ok, pol = refineAt(M9, lo, hi, iwt, relTol)
    return vLim, thetaRoot, pol, residual, ok
end

export LimitingRoot, findRoots, refineRoot
