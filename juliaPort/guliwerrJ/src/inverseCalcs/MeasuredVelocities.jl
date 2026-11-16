# src/inverseCalcs/MeasuredVelocities.jl
#
# EXPERIMENTAL (as in MATLAB): validated on synthetic data only; the interface may change.
#
# Port of `guliwerr/inverseCalcs/MeasuredVelocities.m` -- the measured dataset from one TGS
# angular sweep, and nothing else. Immutable, calculation-free: what was measured, with what
# uncertainty, in what units. No mode number, no eigenvector, no calculated velocity -- that
# binding lives in `ModeAssociation`.
#
# MATLAB is a *value* class (`addPhaseVelocity`/`setBulkModulus` return the modified object) --
# so it becomes an immutable `struct` plus setter *functions* returning a new value (the
# convention used throughout `core/`, e.g. `CalcParams.jl`). `MeasuredVelocities.m`'s own
# name matches this file's main exported type, so (as everywhere in this port) the file
# carries NO `module` wrapper of its own and is `include`d as plain top-level definitions
# into `InverseCalcs` (see `InverseCalcs.jl`'s header).
#
# ERROR-ID CONVENTION (established by `Structure.jl`'s `ApplyWavelength`, e.g.
# `"Structure:ApplyWavelength:expectedPositive"`): MATLAB's hand-written `error(id, ...)`
# calls (constructor `sizeMismatch`/`badWaveType`, `addPhaseVelocity:badWaveType`,
# `fromTable:missingColumns`) keep their id VERBATIM. MATLAB's `validateattributes`-GENERATED
# ids come in two kinds:
#   - pure TYPE checks (e.g. `{'double'}` class failing on a char/string input) are NOT
#     replicated with a bespoke id -- Julia's own type system already rejects these via
#     multiple dispatch (`MethodError`), which is the Julia-idiomatic equivalent (a MATLAB
#     `invalidType` case becomes a `MethodError`).
#   - genuine ATTRIBUTE checks Julia's type system cannot express structurally (`'positive'`,
#     `'finite'`, `numel == 3`) ARE ported, as `ArgumentError`s whose message is prefixed
#     `"ClassName:methodName:attributeCheck"` (dropping MATLAB's own `"MATLAB:"` prefix and
#     dot separator, keeping the semantic suffix -- e.g. `expectedPositive`), matching
#     `Structure:ApplyWavelength:expectedPositive`'s own translation of a `'positive'` check.

using LinearAlgebra: norm

export WAVE_TYPES, MeasuredVelocities, PhaseVelocityEntry, addPhaseVelocity, setBulkModulus,
       fromTable

# -----------------------------------------------------------------------------------------
# WAVE_TYPES -- MeasuredVelocities.m (`properties (Constant)`)
# -----------------------------------------------------------------------------------------

"""
    WAVE_TYPES

Fixed wave-type vocabulary, matching `ModeClassification`'s `type` names (plus
`"sezawa"`/`"love"` for layered systems and `"unknown"` as an escape hatch).

Port of `MeasuredVelocities.WAVE_TYPES`. A `Tuple` of `String`,
not a `Vector`, since it is a fixed compile-time constant, never mutated (MATLAB's
`properties (Constant)` cell array).
"""
const WAVE_TYPES = ("longitudinal", "fastShear", "slowShear", "rayleigh",
                     "pseudoRayleigh", "sezawa", "love", "unknown")

# -----------------------------------------------------------------------------------------
# PhaseVelocityEntry -- one entry of MeasuredVelocities.m's `phaseVelocities` struct array
# -----------------------------------------------------------------------------------------

"""
    PhaseVelocityEntry

One independently measured pulse-echo bulk phase velocity (`MeasuredVelocities.m`'s
`phaseVelocities` struct-array entry, added by [`addPhaseVelocity`](@ref)).

# Fields
- `propagationDir::Vector{Float64}` -- unit propagation direction, `[3]`.
- `waveType::String` -- member of [`WAVE_TYPES`](@ref).
- `velocity::Float64` -- measured bulk phase velocity [km/s].
"""
struct PhaseVelocityEntry
    propagationDir::Vector{Float64}
    waveType::String
    velocity::Float64
end

# -----------------------------------------------------------------------------------------
# The type
# -----------------------------------------------------------------------------------------

"""
    MeasuredVelocities

The measured dataset from one TGS angular sweep. Port of `MeasuredVelocities.m`.

# Fields
- `angleDeg::Vector{Float64}` -- `[nPoints]` propagation angle [deg].
- `velocity::Vector{Float64}` -- `[nPoints]` measured phase velocity [km/s].
- `waveType::Vector{String}` -- `[nPoints]`, each a member of [`WAVE_TYPES`](@ref).
- `velocityStd::Vector{Float64}` -- `[nPoints]` measurement std [km/s]; `NaN` = unknown.
- `phaseVelocities::Vector{PhaseVelocityEntry}` -- optional pulse-echo bulk channel.
- `bulkModulus::Float64` -- soft-fixed bulk modulus [GPa]; `NaN` = not set.
- `nPoints::Int`.

See also the class docblock's note on JOINT MULTI-DATASET INVERSION: for several
crystallographic cuts of the same crystal, an `AbstractVector{MeasuredVelocities}` (one per
dataset) is passed directly to `SensitivityAnalysis.calc`/`ElasticityInversion` -- unchanged
from this type's own per-dataset logic.
"""
struct MeasuredVelocities
    angleDeg::Vector{Float64}
    velocity::Vector{Float64}
    waveType::Vector{String}
    velocityStd::Vector{Float64}
    phaseVelocities::Vector{PhaseVelocityEntry}
    bulkModulus::Float64
    nPoints::Int
end

"""
    MeasuredVelocities() -> MeasuredVelocities

Empty dataset (MATLAB's `nargin == 0` early return, `MeasuredVelocities.m`).
"""
MeasuredVelocities() = MeasuredVelocities(Float64[], Float64[], String[], Float64[],
                                           PhaseVelocityEntry[], NaN, 0)

"""
    MeasuredVelocities(angleDeg, velocity, waveType=nothing, velocityStd=nothing)

Construct the measured-velocity dataset. Port of `MeasuredVelocities.m`.

INPUT: `angleDeg`/`velocity` -- real, finite vectors of the same length `n` (`velocity` must
additionally be positive); `waveType` (optional) -- a single string applied to every point,
or a vector of `n` strings, each a member of [`WAVE_TYPES`](@ref); default `"unknown"` for
every point. `velocityStd` (optional) -- a scalar applied to every point, or a vector of `n`;
default `NaN` (unknown) for every point.

`angleDeg`/`velocity` are typed `AbstractVector{<:Real}` so that a non-numeric-vector
argument (MATLAB's `validateattributes(..., {'double'}, ...)` class check, generating
`MATLAB:MeasuredVelocities:invalidType`) is rejected by Julia's own dispatch (`MethodError`)
instead -- see this file's header note on the type-vs-attribute convention. MATLAB's
row/column distinction has no Julia analogue: a
`Vector` is already one-dimensional.
"""
function MeasuredVelocities(angleDeg::AbstractVector{<:Real}, velocity::AbstractVector{<:Real},
                             waveTypeIn::Union{Nothing,AbstractString,AbstractVector} = nothing,
                             velocityStdIn::Union{Nothing,Real,AbstractVector} = nothing)
    all(isfinite, angleDeg) || throw(ArgumentError(
        "MeasuredVelocities:constructor:expectedFinite -- angleDeg must be finite."))
    all(isfinite, velocity) || throw(ArgumentError(
        "MeasuredVelocities:constructor:expectedFinite -- velocity must be finite."))
    all(>(0), velocity) || throw(ArgumentError(
        "MeasuredVelocities:constructor:expectedPositive -- velocity must be positive."))

    angleDegV = Float64.(collect(angleDeg))
    velocityV = Float64.(collect(velocity))
    n = length(angleDegV)
    length(velocityV) == n || throw(ArgumentError(
        "MeasuredVelocities:constructor:sizeMismatch -- angleDeg ($n) and velocity " *
        "($(length(velocityV))) must have the same length."))

    waveType = _normalize_wave_type(waveTypeIn, n, "MeasuredVelocities:constructor")
    velocityStd = _normalize_velocity_std(velocityStdIn, n, "MeasuredVelocities:constructor")

    return MeasuredVelocities(angleDegV, velocityV, waveType, velocityStd,
                               PhaseVelocityEntry[], NaN, n)
end

# -----------------------------------------------------------------------------------------
# Constructor helpers -- MeasuredVelocities.m
# -----------------------------------------------------------------------------------------

function _normalize_wave_type(waveTypeIn, n::Int, errPrefix::AbstractString)
    if waveTypeIn === nothing || (waveTypeIn isa AbstractVector && isempty(waveTypeIn))
        wt = fill("unknown", n)
    elseif waveTypeIn isa AbstractString
        wt = fill(String(waveTypeIn), n)
    else
        wt = String.(collect(waveTypeIn))
        length(wt) == n || throw(ArgumentError(
            "$errPrefix:sizeMismatch -- waveType ($(length(wt))) must have the same length " *
            "as angleDeg ($n)."))
    end
    badIdx = findfirst(w -> !(w in WAVE_TYPES), wt)
    if badIdx !== nothing
        throw(ArgumentError("$errPrefix:badWaveType -- waveType entry $badIdx " *
            "(\"$(wt[badIdx])\") is not one of: $(join(WAVE_TYPES, ", "))."))
    end
    return wt
end

function _normalize_velocity_std(velocityStdIn, n::Int, errPrefix::AbstractString)
    if velocityStdIn === nothing || (velocityStdIn isa AbstractVector && isempty(velocityStdIn))
        return fill(NaN, n)
    elseif velocityStdIn isa Real
        return fill(Float64(velocityStdIn), n)
    else
        vs = Float64.(collect(velocityStdIn))
        length(vs) == n || throw(ArgumentError(
            "$errPrefix:sizeMismatch -- velocityStd ($(length(vs))) must have length 1 " *
            "or n ($n)."))
        return vs
    end
end

# -----------------------------------------------------------------------------------------
# addPhaseVelocity -- MeasuredVelocities.m
# -----------------------------------------------------------------------------------------

"""
    addPhaseVelocity(obj::MeasuredVelocities, propagationDir, waveType, velocity)
        -> MeasuredVelocities

Append one independently measured pulse-echo bulk phase velocity. Port of
`MeasuredVelocities.addPhaseVelocity`. Value semantics: `obj`
is untouched, a new object is returned. `propagationDir` need not be unit length -- normalised
on use. `propagationDir` is in the inverted layer's CRYSTAL frame, not the dirMat/lab frame;
`waveType` is `"longitudinal"`, `"fastShear"` or `"slowShear"`.
"""
function addPhaseVelocity(obj::MeasuredVelocities, propagationDir::AbstractVector{<:Real},
                           waveType::AbstractString, velocity::Real)
    length(propagationDir) == 3 || throw(ArgumentError(
        "MeasuredVelocities:addPhaseVelocity:incorrectNumel -- propagationDir must have " *
        "exactly 3 elements."))
    all(isfinite, propagationDir) || throw(ArgumentError(
        "MeasuredVelocities:addPhaseVelocity:expectedFinite -- propagationDir must be finite."))
    wt = String(waveType)
    # a bulk wave: the Christoffel branch the channel evaluates (a surface-wave label is
    # rejected here, not later inside SensitivityAnalysis)
    bulkTypes = ("longitudinal", "fastShear", "slowShear")
    wt in bulkTypes || throw(ArgumentError(
        "MeasuredVelocities:addPhaseVelocity:badWaveType -- \"$wt\" is not a bulk wave " *
        "type (one of: $(join(bulkTypes, ", ")))."))
    isfinite(velocity) && velocity > 0 || throw(ArgumentError(
        "MeasuredVelocities:addPhaseVelocity:expectedPositive -- velocity must be positive."))

    dir = Float64.(collect(propagationDir))
    newEntry = PhaseVelocityEntry(dir ./ norm(dir), wt, Float64(velocity))
    return MeasuredVelocities(obj.angleDeg, obj.velocity, obj.waveType, obj.velocityStd,
                               vcat(obj.phaseVelocities, [newEntry]), obj.bulkModulus, obj.nPoints)
end

# -----------------------------------------------------------------------------------------
# setBulkModulus -- MeasuredVelocities.m
# -----------------------------------------------------------------------------------------

"""
    setBulkModulus(obj::MeasuredVelocities, bulkModulusGPa::Real) -> MeasuredVelocities

Set the soft-fixed bulk modulus datum. Port of `MeasuredVelocities.setBulkModulus`
(`MeasuredVelocities.m`). Value semantics: `obj` is untouched.
"""
function setBulkModulus(obj::MeasuredVelocities, bulkModulusGPa::Real)
    isfinite(bulkModulusGPa) && bulkModulusGPa > 0 || throw(ArgumentError(
        "MeasuredVelocities:setBulkModulus:expectedPositive -- bulkModulusGPa must be positive."))
    return MeasuredVelocities(obj.angleDeg, obj.velocity, obj.waveType, obj.velocityStd,
                               obj.phaseVelocities, Float64(bulkModulusGPa), obj.nPoints)
end

# -----------------------------------------------------------------------------------------
# fromTable -- MeasuredVelocities.m
# -----------------------------------------------------------------------------------------

"""
    fromTable(cols::NamedTuple) -> MeasuredVelocities

Build a `MeasuredVelocities` from a plain columns-container, the port of
`MeasuredVelocities.fromTable`. MATLAB's `table` is the ONLY
use of that type anywhere in the ported codebase -- rather than add a `DataFrames.jl`
dependency for this single call site, `cols` is any plain Julia `NamedTuple` exposing the same columns a table would (`cols.angleDeg`,
`cols.velocity`, and optionally `cols.waveType`/`cols.velocityStd`), e.g.
`fromTable((angleDeg = [0;30], velocity = [3.1;3.0]))`. Typing the argument `NamedTuple`
means a non-table-like input (MATLAB's `validateattributes(tbl, {'table'}, ...)` class check,
`MATLAB:MeasuredVelocities:fromTable:invalidType`) is rejected by Julia's own dispatch
(`MethodError`) instead, per this file's header note.
"""
function fromTable(cols::NamedTuple)
    (haskey(cols, :angleDeg) && haskey(cols, :velocity)) || throw(ArgumentError(
        "MeasuredVelocities:fromTable:missingColumns -- cols must have fields \"angleDeg\" " *
        "and \"velocity\" at minimum."))
    waveType = get(cols, :waveType, nothing)
    velocityStd = get(cols, :velocityStd, nothing)
    return MeasuredVelocities(cols.angleDeg, cols.velocity, waveType, velocityStd)
end
