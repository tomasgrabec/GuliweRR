# =========================================================================================
# CalcWavevecAngSweep2Plate.jl -- port of guliwerr/forwardCalcs/CalcWavevecAngSweep2Plate.m:
# 2D sweep (propagation angle AND real wavevector) for the free-plate (Lamb-wave) forward calculation.
# Built on CalcAngSweepPlate: constructs one CalcAngSweepPlate per wavevector value (at
# one FIXED realThickness) and drives it through its own `calc`. Mirrors
# CalcWavevecAngSweep2Layered.jl's per-wavelength composition, simplified in two ways that
# both follow from simpler physics (see CalcSinglePlate.jl):
#   1. No `Structure.ApplyWavelength` step -- `RitzSolverLamb`/`CalcAngSweepPlate` accept
#      `(realThickness, realWavVec)` directly and rescale `relDepth` internally per point.
#   2. No Green/Christoffel phase barrier -- every (angle, wavevector) point is an
#      independent Ritz solve, so there is no `flattenSweep`/`calcFlat` counterpart to the
#      Layered class.
#
# NO per-file module (see Calc.jl's header for why).
#
# MUTABLE: a `mutable struct`, like CalcSingleHS/Layered and their sweeps
# (MATLAB handle semantics): every step mutates the object and returns it, so both
# `calc(x)` and `x = calc(x)` work; `copy(x)` is a shallow copy.
#
# PARALLELISM: this family has no per-wavevector barrier, so `calc` builds every
# (wavevector, angle) point up front and runs ONE flat `forEachElement` over the whole grid,
# then re-assembles one `CalcAngSweepPlate` per wavevector (`sweeps`, unchanged shape).
# MATLAB keeps the outer loop serial (no nested parfor) and parallelises the angles only; the
# per-point results are the same. One `nullConstrMat` serves every point (it depends on
# `calcParams` alone).
#
# RAM: `retentionLevel` DEFAULTS TO `:none`, same reasoning as the Layered sibling
# (eigenvector RAM multiplies across the added wavevector axis) -- but note the ALLOWED SET
# is `(:full, :none)` here, NOT `(:full, :saw_only, :none)`: there is no SAW/limiting-mode
# concept for a free plate (`CalcAngSweepPlate`/`CalcSinglePlate` never had a
# `:saw_only` tier either -- see those files' own headers for why).
#
# BRANCH TRACKING: `DispersionBranches.track` (mass-weighted MAC, one-to-one
# assignment) -- it reads `sweep.angles` and
# `sweep.sweeps[iw].cs[angleIdx].ritzSolver`, which this struct provides.
#
# `confidenceThresh = 0.3` is unvalidated, as in MATLAB: see `classifyBranches`.
# =========================================================================================

using ..CoreTypes: isMirrorSymmetric, isSagittalMirrorSymmetric

"""
    CalcWavevecAngSweep2Plate

2D sweep (propagation angle AND real wavevector) for the free-plate (Lamb-wave) forward
calculation, at one FIXED `realThickness`.

Port of `guliwerr/forwardCalcs/CalcWavevecAngSweep2Plate.m`. **Mutable struct** (MATLAB handle semantics). Built on
[`CalcAngSweepPlate`](@ref): `calc` constructs one
`CalcAngSweepPlate` per wavevector and drives it through its own `calc`. No cross-angle
barrier (unlike [`CalcWavevecAngSweep2Layered`](@ref)) -- every (angle, wavevector) point is
an independent Ritz solve.

# Fields
- `structure::Structure` -- the UNROTATED base structure; `nlay >= 1`.
- `angles::Vector{Float64}` -- propagation angles [deg]; default `0:2:180`.
- `wavVecs::Vector{Float64}` -- positive real wavevectors, strictly increasing -- REQUIRED,
  no default.
- `realThickness::Vector{Float64}` -- physical plate thickness, FIXED across the whole
  sweep -- REQUIRED (without it `wavVecs` has no effect on `relDepth`, making a wavevector
  sweep meaningless).
- `retentionLevel::Symbol` -- `:full` | `:none` (**default `:none`**) -- no `:saw_only` tier.
- `useParallel::Symbol` -- `:auto` (default) | `:always` | `:never` -- governs BOTH the outer
  wavevector loop and each wavevector's own inner angle loop.
- `sweeps::Vector{CalcAngSweepPlate}` -- one per wavevector; filled by [`calc`](@ref) only.
- `velsMap3D`/`velsGroupMap3D`/`velsGroupTransMap3D`/`modeTypeConfidenceMap3D`/
  `freqsPhysicalMap3D`::`Array{Float64,3}`
  -- `[nModes x nAngles x nWav]`, dense (angle-independent mode count, as in
  `CalcAngSweepPlate`, carried one dimension further).
- `modesTypeMap3D::Array{Int,3}` -- `[nModes x nAngles x nWav]`, 1=Asymm/2=Symm/3=TH-Asymm/
  4=TH-Symm.
"""
mutable struct CalcWavevecAngSweep2Plate
    structure::Structure
    angles::Vector{Float64}
    wavVecs::Vector{Float64}
    realThickness::Vector{Float64}

    retentionLevel::Symbol
    useParallel::Symbol

    sweeps::Vector{CalcAngSweepPlate}

    velsMap3D::Array{Float64,3}
    velsGroupMap3D::Array{Float64,3}
    velsGroupTransMap3D::Array{Float64,3}
    modesTypeMap3D::Array{Int,3}
    modeTypeConfidenceMap3D::Array{Float64,3}
    freqsPhysicalMap3D::Array{Float64,3}
end

_calc2plate_isempty(x) = x === nothing || (!(x isa Real) && isempty(x))
_calc2plate_tovec(x) = x isa Real ? [Float64(x)] : Float64.(vec(collect(x)))

"""
    CalcWavevecAngSweep2Plate(; structure = nothing, angles = nothing, wavVecs = nothing,
                                 realThickness = nothing, retentionLevel = :none,
                                 useParallel = :auto) -> CalcWavevecAngSweep2Plate

Construct `CalcWavevecAngSweep2Plate` for a 2D (angle x wavevector) sweep.

Port of `CalcWavevecAngSweep2Plate.m`. `structure === nothing` (MATLAB's
`nargin == 0`) returns a well-formed empty object (`angles` still defaults to `0:2:180`,
`retentionLevel` still defaults to `:none`).

Otherwise: `wavVecs` required (no default), real/finite/positive/strictly increasing, else `ArgumentError`
(`"CalcWavevecAngSweep2Plate:constructor:missingWavVecs"`); `realThickness` required, else
`ArgumentError` (`"CalcWavevecAngSweep2Plate:constructor:missingRealThickness"`) -- without
it `wavVecs` has no effect on `relDepth` (`RitzSolverLamb`'s own constructor contract), making
a wavevector sweep meaningless.

`useParallel` is validated eagerly (mirrors MATLAB's `mustBeMember` property validator).

Keyword-only, not positional -- see `CalcSingleHS`'s constructor for why.
"""
function CalcWavevecAngSweep2Plate(; structure::Union{Structure,Nothing} = nothing,
                                      angles = nothing,
                                      wavVecs = nothing,
                                      realThickness = nothing,
                                      retentionLevel::Symbol = :none,
                                      useParallel::Symbol = :auto)
    _validate_use_parallel(useParallel, "CalcWavevecAngSweep2Plate:constructor")
    validateRetentionLevel(retentionLevel, "CalcWavevecAngSweep2Plate:constructor", (:full, :none))

    angs = (angles === nothing || isempty(angles)) ? _ANGSWEEP_DEFAULT_ANGLES :
           Float64.(vec(collect(angles)))
    all(isfinite, angs) ||
        throw(ArgumentError("CalcWavevecAngSweep2Plate: angles must be a real, finite vector"))

    if structure === nothing
        return CalcWavevecAngSweep2Plate(Structure(), angs, Float64[], Float64[],
            retentionLevel, useParallel, CalcAngSweepPlate[],
            zeros(Float64, 0, 0, 0), zeros(Float64, 0, 0, 0), zeros(Float64, 0, 0, 0),
            zeros(Int, 0, 0, 0),
            zeros(Float64, 0, 0, 0), zeros(Float64, 0, 0, 0))
    end

    _calc2plate_isempty(wavVecs) &&
        throw(ArgumentError("CalcWavevecAngSweep2Plate:constructor:missingWavVecs -- " *
            "wavVecs is required (no default) -- pass a vector of positive real wavevectors."))
    wv = _calc2plate_tovec(wavVecs)
    (all(isfinite, wv) && all(>(0), wv)) ||
        throw(ArgumentError("CalcWavevecAngSweep2Plate: wavVecs must be a real, finite, " *
            "positive vector"))
    # locateZGVs brackets each ZGV between neighbouring grid points
    all(>(0), diff(wv)) ||
        throw(ArgumentError("CalcWavevecAngSweep2Plate:constructor:wavVecsNotIncreasing -- " *
            "wavVecs must be strictly increasing (sort them and drop duplicates)."))

    _calc2plate_isempty(realThickness) &&
        throw(ArgumentError("CalcWavevecAngSweep2Plate:constructor:missingRealThickness -- " *
            "realThickness is required (no default) -- without it, wavVecs has no effect on " *
            "relDepth (see RitzSolverLamb's own constructor), making a wavevector sweep " *
            "meaningless."))
    rt = _calc2plate_tovec(realThickness)

    return CalcWavevecAngSweep2Plate(structure, angs, wv, rt, retentionLevel, useParallel,
        CalcAngSweepPlate[], zeros(Float64, 0, 0, 0), zeros(Float64, 0, 0, 0),
        zeros(Float64, 0, 0, 0), zeros(Int, 0, 0, 0), zeros(Float64, 0, 0, 0),
        zeros(Float64, 0, 0, 0))
end

"""
    _rebuild(sweep::CalcWavevecAngSweep2Plate; <field> = sweep.<field>, ...) ->
        CalcWavevecAngSweep2Plate

In-place field update (the Calc types are mutable, MATLAB handle semantics). Internal; not exported.
"""
function _rebuild(x::CalcWavevecAngSweep2Plate; kw...)
    for (k, v) in kw
        setproperty!(x, k, v)
    end
    return x
end

"""
    copy(x::CalcWavevecAngSweep2Plate)

Shallow copy (MATLAB `copy` of a handle object, one level).
"""
Base.copy(x::CalcWavevecAngSweep2Plate) = CalcWavevecAngSweep2Plate((getfield(x, k) for k in fieldnames(CalcWavevecAngSweep2Plate))...)

# -----------------------------------------------------------------------------------------
# calc -- CalcWavevecAngSweep2Plate.m
# -----------------------------------------------------------------------------------------

"""
    calc(sweep::CalcWavevecAngSweep2Plate, retentionLevel = nothing) -> CalcWavevecAngSweep2Plate

Run the 2D sweep: one [`CalcAngSweepPlate`](@ref) per wavevector, BOTH axes parallelised
directly (see this file's header), then assemble the 3D convenience maps.

Port of `CalcWavevecAngSweep2Plate.calc`.
`retentionLevel` defaults to `sweep.retentionLevel`; validated against `(:full, :none)`.
"""
function calc(sweep::CalcWavevecAngSweep2Plate, retentionLevel::Union{Symbol,Nothing} = nothing)
    rl = retentionLevel === nothing ? sweep.retentionLevel : retentionLevel
    validateRetentionLevel(rl, "CalcWavevecAngSweep2Plate:calc", (:full, :none))
    sweep = _rebuild(sweep; retentionLevel = rl)
    nWav = length(sweep.wavVecs)
    nAng = length(sweep.angles)

    # nullConstrMat and rawByDepthCell depend on calcParams' legOrder alone (not on the
    # wavevector, which only rescales relDepth inside RitzSolverLamb, nor on the angle): built
    # ONCE for the whole grid. (MATLAB builds rawByDepthCell once per wavevector.)
    precomp = (; nullConstrMat = buildConstraintMatrix(sweep.structure.calcParams, LambBoundary()),
                 rawByDepthCell = RawCellFor(sweep.structure))
    pending = Vector{CalcAngSweepPlate}(undef, nWav)
    for iWav in 1:nWav
        pending[iWav] = CalcAngSweepPlate(; structure = sweep.structure, angles = sweep.angles,
                                             realThickness = sweep.realThickness,
                                             realWavVec = sweep.wavVecs[iWav],
                                             retentionLevel = rl, useParallel = sweep.useParallel,
                                             precompStruct = precomp)
    end

    # ONE flat parallel loop over the whole (wavevector, angle) grid -- every point is an
    # independent Lamb solve, so flattening balances the load across threads better than
    # MATLAB's serial outer loop (no nested parfor there) or nested per-wavevector loops.
    # Results are the same per point (each solve is deterministic and BLAS runs single-
    # threaded inside the parallel loop either way).
    flat = Vector{CalcSinglePlate}(undef, nWav * nAng)
    for iWav in 1:nWav, ia in 1:nAng
        flat[(iWav - 1) * nAng + ia] = pending[iWav].cs[ia]
    end
    useP = resolveParallel(sweep.useParallel, length(flat))
    flat = forEachElement(flat, useP, (c, _) -> calcAll(c, rl, precomp))

    sweeps = Vector{CalcAngSweepPlate}(undef, nWav)
    for iWav in 1:nWav
        csW = flat[(iWav - 1) * nAng + 1:iWav * nAng]
        sweeps[iWav] = assembleMaps(_rebuild(pending[iWav]; cs = csW, retentionLevel = rl))
    end
    sweep = _rebuild(sweep; sweeps = sweeps)
    return assembleMaps(sweep)
end

"""
    assembleMaps(sweep::CalcWavevecAngSweep2Plate) -> CalcWavevecAngSweep2Plate

3D re-indexing of each wavevector-slice's already-computed 2D maps -- no new classification
logic.

Port of `CalcWavevecAngSweep2Plate.assembleMaps`.
`freqsPhysicalMap3D` is always fully populated here (unlike `CalcAngSweepPlate`'s own
possibly-empty field) because `wavVecs` is validated positive at construction, so every
`sweeps[iWav]`'s own `realWavVec != 0`.
"""
function assembleMaps(sweep::CalcWavevecAngSweep2Plate)
    nWav = length(sweep.sweeps)
    nAngles = length(sweep.angles)
    nModes = size(sweep.sweeps[1].velsMap, 1)

    velsMap3D = fill(NaN, nModes, nAngles, nWav)
    velsGroupMap3D = fill(NaN, nModes, nAngles, nWav)
    velsGroupTransMap3D = fill(NaN, nModes, nAngles, nWav)
    modesTypeMap3D = zeros(Int, nModes, nAngles, nWav)
    modeTypeConfidenceMap3D = fill(NaN, nModes, nAngles, nWav)
    freqsPhysicalMap3D = fill(NaN, nModes, nAngles, nWav)

    for iWav in 1:nWav
        velsMap3D[:, :, iWav] = sweep.sweeps[iWav].velsMap
        velsGroupMap3D[:, :, iWav] = sweep.sweeps[iWav].velsGroupMap
        velsGroupTransMap3D[:, :, iWav] = sweep.sweeps[iWav].velsGroupTransMap
        modesTypeMap3D[:, :, iWav] = sweep.sweeps[iWav].modesTypeMap
        modeTypeConfidenceMap3D[:, :, iWav] = sweep.sweeps[iWav].modeTypeConfidenceMap
        freqsPhysicalMap3D[:, :, iWav] = sweep.sweeps[iWav].freqsPhysicalMap
    end

    return _rebuild(sweep; velsMap3D, velsGroupMap3D, velsGroupTransMap3D, modesTypeMap3D,
                            modeTypeConfidenceMap3D, freqsPhysicalMap3D)
end

# -----------------------------------------------------------------------------------------
# reindexByBranch -- CalcWavevecAngSweep2Plate.m (MATLAB `Static, Access = private`;
# not exported)
# -----------------------------------------------------------------------------------------

"""
    reindexByBranch(map3D::AbstractArray{T,3}, angleIdx::Integer,
                     branchModeNoMap::AbstractMatrix{<:Integer}) where {T<:Real} -> Matrix{T}

Re-index one `[nModes x nAngles x nWav]` map from RAW Ritz mode order into TRACKED BRANCH
order, at one fixed angle.

Port of `CalcWavevecAngSweep2Plate.reindexByBranch`.
`byBranch[iw, ib]` is that branch's value at wavevector step `iw`. The indirection is not
optional: the RAW Ritz mode index at a given wavevector step is a numerical ordering by
eigenvalue, not a physical identity, and it shifts between steps as branches approach and
cross -- `branchModeNoMap` (from `DispersionBranches.track`) is the MAC-tracked answer
to "which raw index continues this branch at this step".

Generic over `T` (unlike MATLAB's blanket `double`): called with `Int`-eltype
`modesTypeMap3D` and `Float64`-eltype `velsGroupMap3D`/`modeTypeConfidenceMap3D` alike, and
every output element is unconditionally overwritten inside the loop below (MATLAB's own
`NaN(nWav, nBranch)` preallocation is likewise never actually read as NaN in the returned
result -- both languages fully overwrite every entry), so `Matrix{T}(undef, ...)` needs no
NaN sentinel, avoiding the "NaN literal in an `Int` array" type mismatch a direct transliteration
would hit.
"""
function reindexByBranch(map3D::AbstractArray{T,3}, angleIdx::Integer,
                          branchModeNoMap::AbstractMatrix{<:Integer}) where {T<:Real}
    nWav, nBranch = size(branchModeNoMap)
    byBranch = Matrix{T}(undef, nWav, nBranch)
    for iw in 1:nWav
        rawSlice = @view map3D[:, angleIdx, iw]
        byBranch[iw, :] = rawSlice[branchModeNoMap[iw, :]]
    end
    return byBranch
end

# -----------------------------------------------------------------------------------------
# locateZGVs -- CalcWavevecAngSweep2Plate.locateZGVs / refineZGV / nearestBranch
# -----------------------------------------------------------------------------------------

"""
    ZGVPoint

One entry of [`locateZGVs`](@ref)'s output: a point of zero group velocity along k.

# Fields
- `branchNo::Int` -- tracked-branch index (column of `DispersionBranches.track`'s output,
  NOT a raw Ritz mode number).
- `wavVecZGV::Float64` -- refined wavevector.
- `freqZGV::Float64` -- refined physical frequency.
- `angle::Float64` -- the angle [deg] this point was searched at.
- `velGroupTrans::Float64` -- transverse group velocity `c_gtheta` there (perpendicular to k,
  positive towards increasing angle) [same unit as the phase velocity].
- `isTrueZGV::Bool` -- `|velGroupTrans| <= ZGV_VG_REL_TOL * vp`: the whole group-velocity
  vector vanishes (a true ZGV point); `false` = a TGV wave.
"""
struct ZGVPoint
    branchNo::Int
    wavVecZGV::Float64
    freqZGV::Float64
    angle::Float64
    velGroupTrans::Float64
    isTrueZGV::Bool
end

"""
Largest `|vg|/vp` accepted as a ZGV, checked on an independent Ritz solve at the Newton point
(a genuine root sits at `|vg|/vp` ~1e-12 or below); the same bound on `|c_gtheta|/vp` decides
`isTrueZGV` (a TGV wave sits at ~1e-2, Si(001) at 18.4 deg). Near a symmetry direction
`|c_gtheta|/vp` grows linearly with the angle offset, so a point a little off it still counts as
a true ZGV: on Si(001) the slope is 1e-3 per degree off [100] and 7e-4 off [110], i.e. within
0.10 deg / 0.14 deg. MATLAB `CalcWavevecAngSweep2Plate.ZGV_VG_REL_TOL`.
"""
const ZGV_VG_REL_TOL = 1e-4
"Two Newton results closer than this in k and in f (relative) are one point."
const ZGV_DUP_REL_TOL = 1e-8
"Largest relative frequency mismatch between the Newton point and its matched Ritz mode."
const ZGV_FREQ_MATCH_REL_TOL = 1e-6
"""
A step where `|vg|` is below half its value at both neighbours and below this fraction of
`vp`, with no sign change next to it, also seeds Newton (two ZGV points between two grid steps;
[`_zgvCandidates`](@ref)). Chosen from one case (the dip there sits
at `|vg|/vp` = 1.1e-3); revisit when another case misses a pair of points. MATLAB
`CalcWavevecAngSweep2Plate.ZGV_DIP_REL_TOL`.
"""
const ZGV_DIP_REL_TOL = 0.01

"""
    locateZGVs(sweep::CalcWavevecAngSweep2Plate, angleIdx::Integer = 1,
               maxCount::Real = Inf) -> Vector{ZGVPoint}

Locate zero-group-velocity (ZGV) points along the wavevector axis at one fixed angle:
candidate scan, Newton refinement, ZGV/TGV flag. Port of `CalcWavevecAngSweep2Plate.locateZGVs`.

A located point has zero group velocity ALONG k (`velsGroup`). It is a true ZGV point (a
minimum or saddle of the dispersion surface `omega(kx, ky)`) when the component perpendicular
to k, `velsGroupTrans`, vanishes too: on a direction with a sagittal mirror (e.g. <100>, <110>
of a cubic (001) plate), and generally not elsewhere, where the point is a TGV wave that
carries energy sideways (Kiefer, Mezil & Prada 2023, Sci. Adv. 9, eadk6846). Both kinds are
returned; `isTrueZGV` tells them apart. Requires `retentionLevel = :full`
(`DispersionBranches.track`'s own requirement).

STAGE 1 (candidate scan): a sign change of a tracked branch's `velsGroupMap3D` between two
wavevector steps -- no extra solves.

STAGE 1 also seeds from a dip of `|vg|` without a sign change (two ZGV points between two grid
steps; [`_zgvCandidates`](@ref)).

STAGE 2 (Newton refinement, `_refineZGV`): the first seed is the linear interpolation of the
sign change (k0 and f0 along the tracked branch), then the bracket ends (with the sweep's stored
eigenvectors) when Newton does not converge or lands outside the widened bracket;
[`ZgvNewton`](@ref)`.solve` then solves the plate
eigenproblem written as a polynomial in k ([`kPolynomial`](@ref), built once per call) with the
zero-group-velocity condition, for k, omega and the mode shape at once (adapted from GEWtool's
ZGVNewtonBeta, Kiefer et al. 2023, JASA 153, 1386). No raw mode number is held fixed, so a
branch crossing inside the bracket does not matter. The candidates are refined independently
(threads, per `useParallel`).

ACCEPTANCE (as MATLAB): Newton converged, inside the swept k range, and an independent
`RitzSolverLamb` solve there has a mode at that frequency (`ZGV_FREQ_MATCH_REL_TOL`) with
`|vg| <= ZGV_VG_REL_TOL * vp`. Duplicates (two seeds onto one point, `ZGV_DUP_REL_TOL`) are
reported once. `branchNo`: the candidate's branch when the point lies inside its bracket, else
the tracked branch nearest in frequency in the step containing it. `maxCount` counts accepted
points, in candidate order (branch, then wavevector step).

Off-axis stationary points of `omega(kx, ky)` (a true ZGV off the symmetry directions of a
low-symmetry plate) need a 2-D search in (k, angle), not done here.
"""
function locateZGVs(sweep::CalcWavevecAngSweep2Plate, angleIdx::Integer = 1,
                     maxCount::Real = Inf)
    nAngles = length(sweep.angles)
    (1 <= angleIdx <= nAngles) ||
        throw(ArgumentError("CalcWavevecAngSweep2Plate.locateZGVs: angleIdx must be an " *
            "integer in 1:$nAngles, got $angleIdx"))
    branchModeNoMap = DispersionBranches.track(sweep, angleIdx)
    nCompModes = size(branchModeNoMap, 2)
    velsGroupByBranch = reindexByBranch(sweep.velsGroupMap3D, angleIdx, branchModeNoMap)
    freqByBranch = reindexByBranch(sweep.freqsPhysicalMap3D, angleIdx, branchModeNoMap)
    ang = sweep.angles[angleIdx]
    maxCount <= 0 && return ZGVPoint[]

    pc0 = sweep.sweeps[1].precompStructForRitz
    precomp = (; nullConstrMat = pc0.nullConstrMat,
                 rawByDepthCell = precomp_field(pc0, :rawByDepthCell,
                                                () -> RawCellFor(sweep.structure)))
    poly = kPolynomial(sweep.structure, ang, sweep.realThickness, precomp)
    k = sweep.wavVecs

    # STAGE 1: the candidates, in MATLAB's order (branch, then wavevector step), with seeds
    velsByBranch = reindexByBranch(sweep.velsMap3D, angleIdx, branchModeNoMap)
    cands = Tuple{Int,NamedTuple}[]
    for ib in 1:nCompModes
        for c in _zgvCandidates(k, velsGroupByBranch[:, ib], freqByBranch[:, ib],
                                velsByBranch[:, ib], branchModeNoMap[:, ib])
            push!(cands, (ib, c))
        end
    end

    # STAGE 2: refine every candidate (independent -> threads), then dedupe and cap in
    # candidate order -- the same list as MATLAB's sequential loop with its early break
    useP = resolveParallel(sweep.useParallel, length(cands))
    refined = forEachElement(Any[c for c in cands], useP,
        (c, _) -> _refineZGV(sweep, angleIdx, poly, c[2].seeds, c[2].kWin, precomp))
    out = ZGVPoint[]
    for (c, r) in zip(cands, refined)
        length(out) >= maxCount && break
        ib, cand = c
        wavVecZGV, freqZGV, isRoot, velGroupTrans, velPhase = r
        isRoot || continue
        any(p -> abs(p.wavVecZGV - wavVecZGV) <= ZGV_DUP_REL_TOL * wavVecZGV &&
                 abs(p.freqZGV - freqZGV) <= ZGV_DUP_REL_TOL * freqZGV, out) && continue
        branchNo = (cand.kBracket[1] <= wavVecZGV <= cand.kBracket[2]) ? ib :
                   _nearestBranch(k, freqByBranch, wavVecZGV, freqZGV)
        push!(out, ZGVPoint(branchNo, wavVecZGV, freqZGV, ang, velGroupTrans,
                            abs(velGroupTrans) <= ZGV_VG_REL_TOL * abs(velPhase)))
    end
    return out
end

"""
    _zgvCandidates(k, vg, f, vp, rawMode) -> Vector{NamedTuple}

Newton seeds for the ZGV points of one tracked branch, in wavevector order; each candidate is
`(; iw, seeds, kWin, kBracket)` with `seeds` an `nSeed x 4` matrix of rows `(k0, f0, iw,
rawMode)` (`iw = 0`: interpolated, no stored eigenvector). Port of MATLAB's private
`CalcWavevecAngSweep2Plate.zgvCandidates`:
(1) a sign change of `vg` between steps `iw0` and `iw0+1`: the linear interpolation of the zero,
then the two ends; window = the bracket widened by one step;
(2) a dip: an interior step where `|vg|` is below half its value at both neighbours and below
`ZGV_DIP_REL_TOL * vp`, with no sign change next to it (two ZGV points of one branch between two
grid steps leave `vg` with one sign at both ends): three single-seed candidates, from the dip
step and its two neighbours.
"""
function _zgvCandidates(k::AbstractVector, vg::AbstractVector, f::AbstractVector,
                        vp::AbstractVector, rawMode::AbstractVector)
    nW = length(k)
    cands = NamedTuple[]
    sc = [sign(vg[i]) * sign(vg[i + 1]) < 0 for i in 1:nW-1]
    for iw0 in findall(sc)
        t  = vg[iw0] / (vg[iw0] - vg[iw0 + 1])
        k0 = k[iw0] + t * (k[iw0 + 1] - k[iw0])
        f0 = f[iw0] + t * (f[iw0 + 1] - f[iw0])
        push!(cands, (; iw = Float64(iw0),
            seeds = [k0 f0 0 0; k[iw0] f[iw0] iw0 rawMode[iw0]; k[iw0 + 1] f[iw0 + 1] iw0 + 1 rawMode[iw0 + 1]],
            kWin = (k[max(iw0 - 1, 1)], k[min(iw0 + 2, nW)]), kBracket = (k[iw0], k[iw0 + 1])))
    end
    a = abs.(vg)
    for iw in 2:nW-1
        isDip = a[iw] < 0.5 * min(a[iw - 1], a[iw + 1]) && !sc[iw - 1] && !sc[iw] &&
                a[iw] <= ZGV_DIP_REL_TOL * abs(vp[iw])
        isDip || continue
        for j in (iw - 1, iw, iw + 1)
            push!(cands, (; iw = iw + 0.1 * (j - iw), seeds = [k[j] f[j] j rawMode[j]],
                kWin = (k[iw - 1], k[iw + 1]), kBracket = (k[iw - 1], k[iw + 1])))
        end
    end
    return cands[sortperm([c.iw for c in cands])]
end

"""
    _refineZGV(sweep, angleIdx, poly, seeds, kWin, precomp)
        -> (wavVecZGV, freqZGV, isRoot, velGroupTrans, velPhase)

One Newton ZGV refinement from up to three seeds, then its check (see [`locateZGVs`](@ref)).
Port of MATLAB's private `CalcWavevecAngSweep2Plate.refineZGV`: the seeds are tried in order; a
converged point inside `kWin` is taken at once, else the first converged point inside the swept
range. A seed row with `iw > 0` starts from the sweep's stored eigenvector of its raw mode.
"""
function _refineZGV(sweep::CalcWavevecAngSweep2Plate, angleIdx::Integer, poly::NamedTuple,
                    seeds::AbstractMatrix, kWin, precomp::NamedTuple)
    kTol = 1e-12 * maximum(sweep.wavVecs)
    kZ, fZ = NaN, NaN
    for is in 1:size(seeds, 1)
        iwS = Int(seeds[is, 3])
        u0 = iwS > 0 ? sweep.sweeps[iwS].cs[angleIdx].ritzSolver.eigVecs_con[:, Int(seeds[is, 4])] : nothing
        kN, omegaN, _, ok, _ = ZgvNewton.solve(poly.L2, poly.L1, poly.L0, poly.M, seeds[is, 1],
                                               2π * seeds[is, 2], u0)
        (!ok || kN < minimum(sweep.wavVecs) - kTol || kN > maximum(sweep.wavVecs) + kTol) && continue
        isnan(kZ) && ((kZ, fZ) = (kN, omegaN / (2π)))
        if kWin[1] <= kN <= kWin[2]
            kZ, fZ = kN, omegaN / (2π)
            break
        end
    end
    isnan(kZ) && return kZ, fZ, false, NaN, NaN
    rs = calc(RitzSolverLamb(; structure = sweep.structure, angleDeg = sweep.angles[angleIdx],
                               realThickness = sweep.realThickness, realWavVec = kZ),
              false, precomp)
    dev, m = findmin(abs.(rs.freqsPhysical .- fZ))
    dev > ZGV_FREQ_MATCH_REL_TOL * fZ && return kZ, fZ, false, NaN, NaN
    isRoot = abs(rs.velsGroup[m]) <= ZGV_VG_REL_TOL * abs(rs.vels[m])
    return isRoot ? (kZ, fZ, true, rs.velsGroupTrans[m], rs.vels[m]) : (kZ, fZ, false, NaN, NaN)
end

"""
    _nearestBranch(wavVecs, freqByBranch, wavVecZGV, freqZGV) -> Int

Tracked branch nearest a point in frequency, in the wavevector step containing it (linear
interpolation; the end steps on the range boundary). Port of MATLAB's private
`CalcWavevecAngSweep2Plate.nearestBranch`.
"""
function _nearestBranch(wavVecs::AbstractVector, freqByBranch::AbstractMatrix,
                        wavVecZGV::Real, freqZGV::Real)
    iwF = findlast(<=(wavVecZGV), wavVecs)
    iw  = clamp(iwF === nothing ? 1 : iwF, 1, length(wavVecs) - 1)
    t   = (wavVecZGV - wavVecs[iw]) / (wavVecs[iw + 1] - wavVecs[iw])
    fI  = freqByBranch[iw, :] .+ t .* (freqByBranch[iw + 1, :] .- freqByBranch[iw, :])
    return argmin(abs.(fI .- freqZGV))
end

# -----------------------------------------------------------------------------------------
# classifyBranches -- CalcWavevecAngSweep2Plate.m
# -----------------------------------------------------------------------------------------

const _PLATE_MODE_TYPE_NAMES = ("Asymm", "Symm", "TH-Asymm", "TH-Symm")
const _PLATE_PARITY_NAMES = ("Asymm-class", "Symm-class", "Asymm-class", "Symm-class")

"""
    ClassifyBranchesResult

Output of [`classifyBranches`](@ref).

# Fields
- `isExactMirrorSymmetric::Bool` -- `sweep.structure`'s WHOLE-STRUCTURE mid-plane mirror flag
  ([`isMirrorSymmetric`](@ref)): every mode's parity is exact.
- `isExactSagittalSymmetric::Bool` -- [`isSagittalMirrorSymmetric`](@ref) at this angle: the
  Lamb/TH split is exact. Both `true`: every label is rigorous regardless of confidence.
- `confidenceThresh::Float64` -- echoes the input.
- `branchModeNoMap::Matrix{Int}` -- `[nWav x nCompModes]`, from
  `DispersionBranches.track`.
- `modesTypeByBranch::Matrix{Int}` -- `[nWav x nCompModes]`, raw `modesType`
  (1=Asymm, 2=Symm, 3=TH-Asymm, 4=TH-Symm) re-indexed into branch order.
- `confidenceByBranch::Matrix{Float64}` -- `[nWav x nCompModes]`, likewise for
  `modeTypeConfidence`.
- `trustLabelByBranch::Matrix{Bool}` -- `[nWav x nCompModes]`: (both exact flags)
  `|| confidence >= confidenceThresh`, evaluated PER (step, branch).
- `displayLabel::Matrix{String}` -- the type name where `trustLabelByBranch` is `true`; else
  `"Asymm-class"`/`"Symm-class"` where only the parity is exact; else `"branch N"`.
"""
struct ClassifyBranchesResult
    isExactMirrorSymmetric::Bool
    isExactSagittalSymmetric::Bool
    confidenceThresh::Float64
    branchModeNoMap::Matrix{Int}
    modesTypeByBranch::Matrix{Int}
    confidenceByBranch::Matrix{Float64}
    trustLabelByBranch::Matrix{Bool}
    displayLabel::Matrix{String}
end

"""
    classifyBranches(sweep::CalcWavevecAngSweep2Plate, angleIdx::Integer = 1,
                      confidenceThresh::Real = 0.3) -> ClassifyBranchesResult

Multi-layer-aware mode-type classification with a confidence-driven fallback.

Port of `CalcWavevecAngSweep2Plate.classifyBranches`.
A per-point confidence score rather than a static geometric threshold, because how much a
multi-layer stack's asymmetry matters depends on `kd` (wavevector x thickness), not on
geometry alone, and a fixed cutoff cannot track that. Requires `retentionLevel = :full`,
inherited from `DispersionBranches.track` (throws `ArgumentError` otherwise, as in MATLAB).

# `confidenceThresh` default (0.3): unvalidated

`confidence` = top parity fraction minus the runner-up (`modeTypeConfidence`), so `0` is
exactly tied and `1` is a pure type. The default is NOT tuned against measured multi-layer data
or the literature. Treat it as a starting point, not a derived constant.
"""
function classifyBranches(sweep::CalcWavevecAngSweep2Plate, angleIdx::Integer = 1,
                           confidenceThresh::Real = 0.3)
    branchModeNoMap = DispersionBranches.track(sweep, angleIdx)
    nWav, nCompModes = size(branchModeNoMap)

    isExactSymmetric = isMirrorSymmetric(sweep.structure)
    isExactSagittal = isSagittalMirrorSymmetric(sweep.structure, sweep.angles[angleIdx])

    modesTypeByBranch = reindexByBranch(sweep.modesTypeMap3D, angleIdx, branchModeNoMap)
    confidenceByBranch = reindexByBranch(sweep.modeTypeConfidenceMap3D, angleIdx, branchModeNoMap)

    # Per (step, branch), not once globally -- confidence varies along a branch. Mid-plane
    # AND sagittal mirrors make every label rigorous; the mid-plane mirror alone, the parity.
    trustLabelByBranch = (isExactSymmetric && isExactSagittal) .|
        (confidenceByBranch .>= confidenceThresh)

    displayLabel = Matrix{String}(undef, nWav, nCompModes)
    for iw in 1:nWav, ib in 1:nCompModes
        displayLabel[iw, ib] =
            trustLabelByBranch[iw, ib] ? _PLATE_MODE_TYPE_NAMES[modesTypeByBranch[iw, ib]] :
            isExactSymmetric ? _PLATE_PARITY_NAMES[modesTypeByBranch[iw, ib]] : "branch $ib"
    end

    return ClassifyBranchesResult(isExactSymmetric, isExactSagittal, Float64(confidenceThresh),
        branchModeNoMap, modesTypeByBranch, confidenceByBranch, trustLabelByBranch, displayLabel)
end
