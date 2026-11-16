# =========================================================================================
# CalcSinglePlate.jl -- port of guliwerr/forwardCalcs/CalcSinglePlate.m:
# single-angle container for the free-plate (Lamb-wave) forward calculation (RitzSolverLamb
# only).
#
# NO per-file module (see Calc.jl's header for why).
#
# MUTABLE: a `mutable struct`, like CalcSingleHS/Layered and their sweeps
# (MATLAB handle semantics): every step mutates the object and returns it, so both
# `calc(x)` and `x = calc(x)` work; `copy(x)` is a shallow copy.
#
# LEANER than CalcSingleHS/CalcSingleLayered -- no `calcLim`, no `calcGreen`, no
# `identifyModes`/`classification` field. Both reasons are physical, per MATLAB's own
# docblock: (1) there is no Green's-function analogue for a free plate in this codebase
# (GreenSolver models a half-space/layered-substrate surface response to a point force); (2)
# a free plate's boundary conditions are the TRUE physical ones, not an artificial truncation
# of an infinite domain the way the half-space's rigid-bottom BC is, so every converged Ritz
# mode is a genuine Lamb branch -- there is no "which numerical mode is the real one"
# question for ModeClassification to answer. Mode TYPE (Asymm/Symm/TH-Asymm/TH-Symm) is a
# separate question, answered inline by `RitzSolverLamb.modeRecognition`
# (`ritzSolver.modesType`/`.modeTypeConfidence`).
#
# `applyRetention`/`saw_only` REJECTION: this file
# does its OWN two-way (`:full`/`:none`) retention switch and deliberately does NOT route
# through `CalcCommon.pruneRitzTo` -- that helper ACCEPTS `:saw_only` (the SAW family's own
# three-tier contract), and a plate has no "significant subset" concept to prune to (every
# mode is physically real, per the docblock above) -- routing a plate through `pruneRitzTo`
# would make it silently ACCEPT `:saw_only`.
# `applyRetention`'s `else` branch below calls [`validateRetentionLevel`](@ref) with the
# allowed set `(:full, :none)`, which THROWS for `:saw_only` (or anything else).
#
# DOUBLE-ROTATION INVARIANT -- as for `CalcSingleHS.jl`/`CalcSingleLayered.jl`: `structure` is
# rotated ONCE at construction and the solver is then called at angle 0. Checked against
# MATLAB with a rotated plate at angle = 25.
# =========================================================================================

"""
    CalcSinglePlate

Single-angle container for the free-plate (Lamb-wave) forward calculation.

Port of `guliwerr/forwardCalcs/CalcSinglePlate.m`. **Mutable struct** (MATLAB handle semantics; see this file's
header). Given one `Structure` (free plate, `nlay >= 1`) and one propagation angle, runs
`RitzSolverLamb` for that direction. Build one with the keyword constructor, then thread it
through [`calcRitz`](@ref)/[`applyRetention`](@ref) individually, or run the whole pipeline
at once with [`calcAll`](@ref).

# Fields
- `structure::Structure` -- the ROTATED structure for this angle.
- `angle::Float64` -- propagation angle [deg], reference only.
- `realThickness::Vector{Float64}` -- physical plate thickness, one value per layer or a
  single value broadcast to all layers (empty = normalised `relDepth` units only, matching
  `RitzSolverLamb`'s own default).
- `realWavVec::Float64` -- physical wavevector `k`; `0.0` means "no physical-frequency
  conversion" (matches `RitzSolverLamb`'s own default).
- `ritzSolver::RitzSolverLamb` -- [`calcRitz`](@ref)'s result; mode-type classification runs
  inline (`withModeTypes = true`), so `modesType`/`modeTypeConfidence`/`freqsPhysical` are
  populated as soon as `calcRitz` returns.
- `retentionLevel::Symbol` -- `:full` (default) | `:none` -- **no `:saw_only` tier** (see this
  file's header).
- `precompStructForRitz::NamedTuple` -- `nullConstrMat` only (unconditionally buildable from
  `calcParams` alone); `gammaByDepthCell` is deliberately NEVER auto-built here (it bakes in
  the ROTATED elastic tensor at a SPECIFIC angle, reusable only across repeated wavevector
  evaluations at one fixed angle -- a caller in that position may still pass one in, but this
  constructor cannot invent one, having no angle to invent it for).
"""
mutable struct CalcSinglePlate
    structure::Structure
    angle::Float64

    realThickness::Vector{Float64}
    realWavVec::Float64

    ritzSolver::RitzSolverLamb

    retentionLevel::Symbol

    precompStructForRitz::NamedTuple
end

"""
    CalcSinglePlate(; structure = nothing, angle = nothing, realThickness = nothing,
                      realWavVec = nothing, precompStruct = NamedTuple(),
                      retentionLevel = :full) -> CalcSinglePlate

Construct `CalcSinglePlate` for one propagation angle.

Port of `CalcSinglePlate.m`. `structure === nothing` (MATLAB's `nargin == 0`) returns
a well-formed empty object with no rotation/precompute work attempted.

Otherwise: `structure` is rotated ONCE (see `CalcSingleHS.jl`'s header for the invariant this
maintains), and `precompStruct.nullConstrMat` is filled in via
[`buildConstraintMatrix`](@ref)`(structRot.calcParams, LambBoundary())` if missing -- any
OTHER field the caller supplied (e.g. `gammaByDepthCell`) survives untouched.

Keyword-only, not positional -- see `CalcSingleHS`'s constructor for why.
"""
function CalcSinglePlate(; structure::Union{Structure,Nothing} = nothing,
                            angle::Union{Real,Nothing} = nothing,
                            realThickness = nothing,
                            realWavVec::Union{Real,Nothing} = nothing,
                            precompStruct::NamedTuple = NamedTuple(),
                            retentionLevel::Symbol = :full)
    if structure === nothing
        return CalcSinglePlate(Structure(), 0.0, Float64[], 0.0, RitzSolverLamb(),
                                :full, NamedTuple())
    end

    # MATLAB validates a given retentionLevel already in the constructor
    validateRetentionLevel(retentionLevel, "CalcSinglePlate:constructor", (:full, :none))

    ang = angle === nothing ? 0.0 : Float64(angle)
    structRot = rotateStructure(structure, ang)

    rt = realThickness === nothing ? Float64[] :
         (realThickness isa Real ? [Float64(realThickness)] : Float64.(collect(realThickness)))
    rwv = realWavVec === nothing ? 0.0 : Float64(realWavVec)

    nullConstrMat = precomp_field(precompStruct, :nullConstrMat,
        () -> buildConstraintMatrix(structRot.calcParams, LambBoundary()))
    precomp = merge(precompStruct, (; nullConstrMat = nullConstrMat))

    return CalcSinglePlate(structRot, ang, rt, rwv, RitzSolverLamb(), retentionLevel, precomp)
end

"""
    _rebuild(cs::CalcSinglePlate; <field> = cs.<field>, ...) -> CalcSinglePlate

In-place field update (the Calc types are mutable, MATLAB handle semantics). Internal; not exported.
"""
function _rebuild(x::CalcSinglePlate; kw...)
    for (k, v) in kw
        setproperty!(x, k, v)
    end
    return x
end

"""
    copy(x::CalcSinglePlate)

Shallow copy (MATLAB `copy` of a handle object, one level).
"""
Base.copy(x::CalcSinglePlate) = CalcSinglePlate((getfield(x, k) for k in fieldnames(CalcSinglePlate))...)

# -----------------------------------------------------------------------------------------
# calcRitz -- CalcSinglePlate.m
# -----------------------------------------------------------------------------------------

"""
    calcRitz(cs::CalcSinglePlate, precomp = nothing) -> CalcSinglePlate

Run `RitzSolverLamb` for this angle, with mode-type classification inline
(`withModeTypes = true`).

Port of `CalcSinglePlate.calcRitz`. `precomp` is the sweep's
shared Ritz precompute (`nullConstrMat`); `nothing` = this object's own
`precompStructForRitz`. A precompute without `nullConstrMat` gets one built (MATLAB does the
same, with a worker warning that has no Julia analogue: Julia threads share memory, so a
sweep's precompute always arrives). `cs.structure` is already rotated, so `angleDeg = 0.0`
is passed here, never `cs.angle`.
"""
function calcRitz(cs::CalcSinglePlate, precomp::Union{NamedTuple,Nothing} = nothing)
    pc = precomp === nothing ? cs.precompStructForRitz : precomp
    if !haskey(pc, :nullConstrMat) || isempty(pc.nullConstrMat)
        pc = merge(pc, (; nullConstrMat = buildConstraintMatrix(cs.structure.calcParams,
                                                                LambBoundary())))
    end
    rs = RitzSolverLamb(; structure = cs.structure, angleDeg = 0.0,
                           realThickness = cs.realThickness, realWavVec = cs.realWavVec)
    rs = calc(rs, true, pc)
    return _rebuild(cs; ritzSolver = rs)
end

# -----------------------------------------------------------------------------------------
# applyRetention -- CalcSinglePlate.m
# -----------------------------------------------------------------------------------------

"""
    applyRetention(cs::CalcSinglePlate, retentionLevel) -> CalcSinglePlate

Apply the eigenvector retention policy to `cs.ritzSolver`.

Port of `CalcSinglePlate.applyRetention`. `:full` keeps
`eigVecs_con` untouched; `:none` clears it entirely (`vels`/`modesType`/
`modeTypeConfidence`/`freqsPhysical`/`velsGroup` all stay -- cheap scalars, not the large
`[strBasisVecL x nModes]` matrix). Anything else (including `:saw_only`) is rejected via
[`validateRetentionLevel`](@ref), which always throws for a value outside `(:full, :none)` --
this file's own two-way switch, deliberately NOT routed through [`pruneRitzTo`](@ref) (see
this file's header).
"""
function applyRetention(cs::CalcSinglePlate, retentionLevel)
    if retentionLevel === :full
        return cs
    elseif retentionLevel === :none
        return _rebuild(cs; ritzSolver = pruneEigVecs(cs.ritzSolver, nothing))
    else
        validateRetentionLevel(retentionLevel, "CalcSinglePlate:applyRetention", (:full, :none))
        return cs   # unreachable: validateRetentionLevel always throws in this branch
    end
end

# -----------------------------------------------------------------------------------------
# calcAll -- CalcSinglePlate.m
# -----------------------------------------------------------------------------------------

"""
    calcAll(cs::CalcSinglePlate, retentionLevel = nothing, precomp = nothing) -> CalcSinglePlate

Run the (only) pipeline for this angle, then prune.

Port of `CalcSinglePlate.calcAll`. `precomp` (the sweep's shared
Ritz precompute) goes to [`calcRitz`](@ref). `calcRitz` -> `applyRetention`
-- no `mode` argument, unlike `CalcSingleHS`/`CalcSingleLayered`'s `calcAll(mode,
retentionLevel)`: there is only one pipeline here.
"""
function calcAll(cs::CalcSinglePlate, retentionLevel::Union{Symbol,Nothing} = nothing,
                 precomp::Union{NamedTuple,Nothing} = nothing)
    rl = retentionLevel === nothing ? cs.retentionLevel : retentionLevel
    validateRetentionLevel(rl, "CalcSinglePlate:calcAll", (:full, :none))

    cs = calcRitz(cs, precomp)
    cs = applyRetention(cs, rl)
    cs = _rebuild(cs; retentionLevel = rl)
    return cs
end
