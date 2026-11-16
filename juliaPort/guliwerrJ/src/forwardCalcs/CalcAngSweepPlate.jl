# =========================================================================================
# CalcAngSweepPlate.jl -- port of guliwerr/forwardCalcs/CalcAngSweepPlate.m: angular sweep of
# CalcSinglePlate objects for the free-plate (Lamb-wave) forward calculation, at one FIXED
# thickness/wavevector. Mirrors CalcAngSweepHS/CalcAngSweepLayered's per-angle composition,
# simplified the same way CalcSinglePlate simplifies CalcSingleHS/CalcSingleLayered: no
# Christoffel/Green, and no "which mode is physically real" question.
#
# NO per-file module (see Calc.jl's header for why).
#
# MUTABLE: a `mutable struct`, like CalcSingleHS/Layered and their sweeps
# (MATLAB handle semantics): every step mutates the object and returns it, so both
# `calc(x)` and `x = calc(x)` work; `copy(x)` is a shallow copy.
#
# DENSE, NOT RAGGED (as `CalcAngSweepPlate.m` explains): every angle shares the SAME number of computed
# modes -- `nCompModes` depends only on `calcParams` (legOrder/layer count), never on angle or
# material rotation -- unlike `CalcAngSweepHS`/`CalcAngSweepLayered`'s ragged per-angle SAW-mode count (different angles
# CAN classify different numbers of significant modes there, since that count comes from
# ModeClassification's Green-peak matching, not from a fixed Ritz basis size). The convenience
# maps below are therefore DENSE `[nCompModes x nAngles]` matrices (`Matrix{Float64}`/
# `Matrix{Int}`), never `Vector{Vector{...}}` -- do not force the ragged HS/layered shape onto
# this class.
#
# SINGLE FLAT LOOP, NO PHASE SPLIT: every angle's Lamb solve
# is fully independent (no Green/Christoffel shared-axis barrier exists here, unlike the HS/
# layered sweeps) -- `calc()` is ONE `forEachElement` call over `cs[i].calcAll(retentionLevel)`,
# not a `calcFullPhases`-style multi-step orchestration. `forEachElement` already pins
# BLAS to 1 thread for its own parallel branch -- no extra wrapping needed here.
# =========================================================================================

"""
    CalcAngSweepPlate

Angular sweep of [`CalcSinglePlate`](@ref) objects for the free-plate (Lamb-wave) forward
calculation, at one fixed thickness/wavevector.

Port of `guliwerr/forwardCalcs/CalcAngSweepPlate.m`. **Mutable struct** (MATLAB handle semantics; see this file's
header). Single flat loop (no phase split -- every angle's Lamb solve is independent), and
DENSE (not ragged) convenience maps, since `nCompModes` is angle-independent.

# Fields
- `structure::Structure` -- the UNROTATED base structure; `nlay >= 1`.
- `angles::Vector{Float64}` -- propagation angles [deg]; default `0:2:180`.
- `realThickness::Vector{Float64}` -- physical plate thickness, FIXED across the whole sweep
  (same meaning as `CalcSinglePlate`'s own field).
- `realWavVec::Float64` -- physical wavevector `k`; `0.0` means "no physical-frequency
  conversion".
- `retentionLevel::Symbol` -- `:full` (default) | `:none` -- no `:saw_only` tier (see
  `CalcSinglePlate.jl`'s header for why).
- `useParallel::Symbol` -- `:auto` (default) | `:always` | `:never`.
- `cs::Vector{CalcSinglePlate}` -- one per angle.
- `velsMap`/`velsGroupMap`/`modeTypeConfidenceMap`::`Matrix{Float64}` -- `[nCompModes x
  nAngles]` phase/group velocities and mode-type-recognition confidence.
- `modesTypeMap::Matrix{Int}` -- `[nCompModes x nAngles]`, 1=Asymm/2=Symm/3=TH-Asymm/
  4=TH-Symm (matches `RitzSolverLamb.modesType`'s own numeric encoding).
- `freqsPhysicalMap::Matrix{Float64}` -- `[nCompModes x nAngles]`, only when `realWavVec !=
  0`; `zeros(Float64, 0, 0)` (MATLAB's `[]`) otherwise.
- `precompStructForRitz::NamedTuple` -- `nullConstrMat` and `rawByDepthCell` (both depend on
  `calcParams`' legOrder alone; `gammaByDepthCell` bakes in a SPECIFIC rotated angle and is
  therefore never shareable across angles -- each `cs[i]` builds its own).
"""
mutable struct CalcAngSweepPlate
    structure::Structure
    angles::Vector{Float64}

    realThickness::Vector{Float64}
    realWavVec::Float64

    retentionLevel::Symbol
    useParallel::Symbol

    cs::Vector{CalcSinglePlate}

    velsMap::Matrix{Float64}
    velsGroupMap::Matrix{Float64}
    velsGroupTransMap::Matrix{Float64}
    modesTypeMap::Matrix{Int}
    modeTypeConfidenceMap::Matrix{Float64}
    freqsPhysicalMap::Matrix{Float64}

    precompStructForRitz::NamedTuple
end

"""
    CalcAngSweepPlate(; structure = nothing, angles = nothing, realThickness = nothing,
                        realWavVec = nothing, retentionLevel = :full,
                        useParallel = :auto, precompStruct = NamedTuple()) -> CalcAngSweepPlate

Construct `CalcAngSweepPlate` for an angular sweep at fixed thickness/wavevector: one
`CalcSinglePlate` per angle, all sharing the SAME `nullConstrMat`.

Port of `CalcAngSweepPlate.m`. `structure === nothing` (MATLAB's `nargin == 0`) returns
a well-formed empty object (`angles` still defaults to `0:2:180`, `realWavVec` stays `0.0`).

Otherwise: `nullConstrMat` is built once via [`buildConstraintMatrix`](@ref)`(structure.
calcParams, LambBoundary())` and shared by every `cs[i] = CalcSinglePlate(; structure, angle =
angles[i], realThickness, realWavVec, precompStruct = precompStructForRitz, retentionLevel)`.

`useParallel` is validated eagerly (mirrors MATLAB's `mustBeMember` property validator), and
so is `retentionLevel` (MATLAB validates a given one in the constructor).

`precompStruct` (Julia-only keyword): a caller that already holds this `calcParams`'
`nullConstrMat` (the 2D [`CalcWavevecAngSweep2Plate`](@ref), one matrix for every
wavevector) passes it in and it is not rebuilt. It depends on `calcParams` alone (MATLAB's
`locateZGVs` relies on the same fact).

Keyword-only, not positional -- see `CalcSingleHS`'s constructor for why.
"""
function CalcAngSweepPlate(; structure::Union{Structure,Nothing} = nothing,
                              angles = nothing,
                              realThickness = nothing,
                              realWavVec::Union{Real,Nothing} = nothing,
                              retentionLevel::Symbol = :full,
                              useParallel::Symbol = :auto,
                              precompStruct::NamedTuple = NamedTuple())
    _validate_use_parallel(useParallel, "CalcAngSweepPlate:constructor")
    validateRetentionLevel(retentionLevel, "CalcAngSweepPlate:constructor", (:full, :none))

    angs = (angles === nothing || isempty(angles)) ? _ANGSWEEP_DEFAULT_ANGLES :
           Float64.(vec(collect(angles)))
    rt = realThickness === nothing ? Float64[] :
         (realThickness isa Real ? [Float64(realThickness)] : Float64.(collect(realThickness)))
    rwv = realWavVec === nothing ? 0.0 : Float64(realWavVec)

    if structure === nothing
        return CalcAngSweepPlate(Structure(), angs, rt, rwv, retentionLevel, useParallel,
                                  CalcSinglePlate[],
                                  zeros(Float64, 0, 0), zeros(Float64, 0, 0), zeros(Float64, 0, 0),
                                  zeros(Int, 0, 0), zeros(Float64, 0, 0),
                                  zeros(Float64, 0, 0), NamedTuple())
    end

    nullConstrMat = precomp_field(precompStruct, :nullConstrMat,
        () -> buildConstraintMatrix(structure.calcParams, LambBoundary()))
    rawByDepthCell = precomp_field(precompStruct, :rawByDepthCell, () -> RawCellFor(structure))
    precomp = (; nullConstrMat = nullConstrMat, rawByDepthCell = rawByDepthCell)

    nAngles = length(angs)
    cs = Vector{CalcSinglePlate}(undef, nAngles)
    for iAngle in 1:nAngles
        cs[iAngle] = CalcSinglePlate(; structure = structure, angle = angs[iAngle],
                                        realThickness = rt, realWavVec = rwv,
                                        precompStruct = precomp, retentionLevel = retentionLevel)
    end

    return CalcAngSweepPlate(structure, angs, rt, rwv, retentionLevel, useParallel, cs,
                              zeros(Float64, 0, 0), zeros(Float64, 0, 0), zeros(Float64, 0, 0),
                              zeros(Int, 0, 0), zeros(Float64, 0, 0),
                              zeros(Float64, 0, 0), precomp)
end

"""
    _rebuild(sweep::CalcAngSweepPlate; <field> = sweep.<field>, ...) -> CalcAngSweepPlate

In-place field update (the Calc types are mutable, MATLAB handle semantics). Internal; not exported.
"""
function _rebuild(x::CalcAngSweepPlate; kw...)
    for (k, v) in kw
        setproperty!(x, k, v)
    end
    return x
end

"""
    copy(x::CalcAngSweepPlate)

Shallow copy (MATLAB `copy` of a handle object, one level).
"""
Base.copy(x::CalcAngSweepPlate) = CalcAngSweepPlate((getfield(x, k) for k in fieldnames(CalcAngSweepPlate))...)

# -----------------------------------------------------------------------------------------
# calc -- CalcAngSweepPlate.m
# -----------------------------------------------------------------------------------------

"""
    calc(sweep::CalcAngSweepPlate, retentionLevel = nothing) -> CalcAngSweepPlate

Run the sweep (Ritz only, every angle -- single flat loop, no phase split), then assemble the
DENSE convenience maps.

Port of `CalcAngSweepPlate.calc`. `retentionLevel` defaults to
`sweep.retentionLevel`; validated against `(:full, :none)` (no `:saw_only` tier).
"""
function calc(sweep::CalcAngSweepPlate, retentionLevel::Union{Symbol,Nothing} = nothing)
    rl = retentionLevel === nothing ? sweep.retentionLevel : retentionLevel
    validateRetentionLevel(rl, "CalcAngSweepPlate:calc", (:full, :none))

    sweep = _rebuild(sweep; retentionLevel = rl)

    # Single flat loop, no phase split: every angle's Lamb solve is fully independent (no
    # Green/Christoffel shared-axis barrier exists here), sharing only the precompute, which
    # goes in explicitly (as in MATLAB).
    precomp = sweep.precompStructForRitz
    useP = resolveParallel(sweep.useParallel, length(sweep.cs))
    cs = forEachElement(sweep.cs, useP, (cs0, _) -> calcAll(cs0, rl, precomp))
    sweep = _rebuild(sweep; cs = cs)

    return assembleMaps(sweep)
end

# -----------------------------------------------------------------------------------------
# assembleMaps -- CalcAngSweepPlate.m (MATLAB `Access = private`; not exported)
# -----------------------------------------------------------------------------------------

"""
    assembleMaps(sweep::CalcAngSweepPlate) -> CalcAngSweepPlate

Extract the sweep-level convenience maps from `sweep.cs` -- DENSE `[nCompModes x nAngles]`
matrices, not ragged (see this file's header).

Port of `CalcAngSweepPlate.assembleMaps`. `nCompModes` is read
from `cs[1]` and assumed identical for every angle (the whole point of this class being safe
to store densely -- see this file's header).
"""
function assembleMaps(sweep::CalcAngSweepPlate)
    nAngles = length(sweep.cs)
    nModes = sweep.cs[1].ritzSolver.nCompModes

    velsMap = zeros(Float64, nModes, nAngles)
    velsGroupMap = zeros(Float64, nModes, nAngles)
    velsGroupTransMap = zeros(Float64, nModes, nAngles)
    modesTypeMap = zeros(Int, nModes, nAngles)
    modeTypeConfidenceMap = zeros(Float64, nModes, nAngles)
    hasPhysicalFreqs = sweep.realWavVec != 0
    freqsPhysicalMap = hasPhysicalFreqs ? zeros(Float64, nModes, nAngles) : zeros(Float64, 0, 0)

    for iAngle in 1:nAngles
        rs = sweep.cs[iAngle].ritzSolver
        velsMap[:, iAngle] = rs.vels[1:nModes]
        velsGroupMap[:, iAngle] = rs.velsGroup[1:nModes]
        velsGroupTransMap[:, iAngle] = rs.velsGroupTrans[1:nModes]
        modesTypeMap[:, iAngle] = rs.modesType[1:nModes]
        modeTypeConfidenceMap[:, iAngle] = rs.modeTypeConfidence[1:nModes]
        if hasPhysicalFreqs
            freqsPhysicalMap[:, iAngle] = rs.freqsPhysical[1:nModes]
        end
    end

    return _rebuild(sweep; velsMap, velsGroupMap, velsGroupTransMap, modesTypeMap,
                            modeTypeConfidenceMap, freqsPhysicalMap)
end
