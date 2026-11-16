# src/inverseCalcs/ModeAssociation.jl -- the binding table
#
# EXPERIMENTAL (as in MATLAB): validated on synthetic data only; the interface may change.
#
# Port of `guliwerr/inverseCalcs/ModeAssociation.m` -- the binding table: which measured
# velocity point corresponds to which calculated object (a Ritz eigenvector for
# "ritzMode"/"pseudoSAW" points, a Christoffel limiting root for "limiting" points). See
# `ModeAssociation.m`'s own class docblock for the full BINDING/ANCHOR/RE-IDENTIFICATION
# vocabulary -- reused verbatim in this file's comments.
#
# `include`d as plain top-level definitions into `module InverseCalcs`, exactly like
# `MeasuredVelocities.jl`/`SensitivityAnalysis.jl` -- no `module ModeAssociation` wrapper
# (`ModeAssociation.m`'s own name matches this file's main exported type, which a wrapping
# module of the same name would collide with).
#
# VALUE TYPE. `ModeAssociation.m`'s own docblock states it: "Value class:
# rebind/rebindLimiting/addPoint/setPinned return the modified object. That is what lets
# ElasticityInversion re-bind into a SCRATCH copy every objective evaluation and promote it
# only when the cost improved, so a rejected trial step can never move an anchor." The
# fields are only ever mutated inside `ModeAssociation.m`'s own methods, never from outside,
# the same discipline as `MeasuredVelocities`/`SensitivityAnalysis`. `ModeAssociation`
# therefore ports to an IMMUTABLE `struct` wrapping `points::Vector{AssocPoint}` (itself an
# immutable `struct`), and every mutating-looking method below (`rebind`, `setPinned`,
# `addPoint`, ...) follows the copy-with-replacements idiom: COPY `obj.points` into a fresh
# `Vector` before writing any element, then wrap the fresh vector in a NEW `ModeAssociation`.
# This is not cosmetic -- `Vector` is a mutable container even when held by an immutable
# struct field, so skipping the copy (e.g. `push!(obj.points, x)`) would silently mutate every
# other reference to the SAME points vector (the exact "current best" corruption the MATLAB
# value-class discipline exists to prevent). `AssocPoint` itself is immutable for the same
# reason one level down: `obj.points(k).pinned = true`-style single-field writes become
# `_updatePoint(p; pinned = true)`, which reconstructs a whole new `AssocPoint` rather than
# mutating a field in place.
#
# THE MATLAB null()/nullspace() BASIS TRAP. The null-space basis N of the constraint matrix is
# not reproducible across implementations: an eigenvector expressed in MATLAB's null-space
# basis cannot be combined with a mass matrix `M_con` expressed in Julia's, because two
# different bases of the same subspace multiplied together give wrong numbers, not merely
# noisy ones. `ModeAssociation.rebind`'s MAC formula has exactly that shape
# (`eigVecs_con(:,m)' * M_con * anchorEigVec_con`), but the trap does not arise here, for two
# independent reasons:
#   1. Production use is self-consistent. `fromSweep` extracts `anchorEigVec_con` straight
#      from THE SAME sweep's `eigVecs_con` (`sweep.cs[iAngle].ritzSolver.
#      eigVecs_con[:, best.modeNo]`), and every later `rebind` call is fed a freshly computed
#      `eigVecs_con`/`M_con` from the same RitzSolver at the SAME angle -- and the null-space
#      basis N is a property of the (elasticity-independent) constraint matrix, so the SAME N
#      is reused across every iterate within one run. Anchor and candidate are therefore
#      always expressed in the same reduced basis, by construction.
#   2. MAC is basis- AND phase-invariant even if it weren't. `alpha_m^H M_con alpha_anchor`
#      with `M_con = N^T M N` equals the FULL-SPACE inner product `(N alpha_m)^H M (N
#      alpha_anchor)` for ANY orthonormal N reconstructing the SAME full-space eigenvectors,
#      and the `abs(...)` around it (`ModeAssociation.m`) additionally kills any residual
#      global phase. Only the raw eigenvector is basis- and phase-dependent; quantities such
#      as MAC and `confidence()` / `limitingRootScore()` are physical and agree between MATLAB
#      and Julia whichever basis each side computed them in.
#
# DUCK-TYPED `pointsAtAngle`. `SensitivityAnalysis.jl` takes `association` as a
# "ModeAssociation-shaped value" that it deliberately does NOT type-check. `pointsAtAngle`
# below is therefore left UNTYPED on its first argument (it reads `assoc.points` generically,
# exactly like `ModeAssociation.pointsAtAngle`'s MATLAB body does on `obj.points`), so that
# `SensitivityAnalysis.jl` can call this same function for any value exposing `.points`, each
# with `.bindingClass`/`.angleDeg`. A real `ModeAssociation` naturally satisfies the same duck
# type: one function, two shapes of caller, not two competing implementations.
#
# FIELD COUNT. `AssocPoint` has 22 fields: the number of name/value pairs in
# `ModeAssociation.m`'s own `struct(...)` literal, all reproduced in `emptyPoint`.
#
# MATLAB-COMMA-LIST-EXPANSION SITES (13 in this file) all become comprehensions below, e.g.
# `{obj.points.bindingClass}` -> `[p.bindingClass for p in obj.points]`,
# `[obj.points.angleDeg]` -> `Float64[p.angleDeg for p in obj.points]`.
#
# `unique` WITH THREE OUTPUTS (`ModeAssociation.m`) -- `removeDuplicatePoints` needs
# MATLAB's `[uniqueKeys, ~, groupId] = unique(keyOf)` group-index output, which Julia's
# `unique` does not provide. Ported as manual `Dict{K,Vector{Int}}` grouping (`get!` per key),
# with a PLAIN TUPLE key (`(angle, ritzModeNo)` or `(angle, waveType)`) rather than a
# `sprintf`-formatted string -- Julia's `Dict` needs no string-hashing workaround the way
# MATLAB's `unique` on a cellstr does, and this also avoids adding a `Printf` dependency for
# a key that is never displayed. Grouping by an UNORDERED `Dict` (instead of MATLAB's SORTED
# `uniqueKeys` loop order) gives an IDENTICAL `keep` mask regardless of iteration order, since
# which duplicate wins is decided by minimum error within each group, not by group visitation
# order -- and the final result preserves the ORIGINAL point order via boolean masking either
# way (`ModeAssociation.m`'s own `obj.points = obj.points(keep)`).
#
# `matchpairs(-subScore, 0)` (`ModeAssociation.m`) -> `matchpairs0(-subScore)` from
# `MatlabCompat` (Hungarian assignment restricted to the negative-cost entries, as
# `matchpairs(C, 0)` does).
#
# ACCESS-MODIFIER MAPPING -- MATLAB's `methods (Static)` (public) vs `methods (Static, Access
# = private)` is preserved at this file's OWN export boundary (Julia has no enforced privacy,
# but the package's public surface should still reflect it): `fromSweep`, `defaultTrackOpts`,
# `limitingRootScore` are exported; `emptyPoint`/`mergeTrackOpts`/`assignWithPinnedPriority`/
# `classifyWaveType` (MATLAB `Access = private`) are defined but NOT exported -- still directly
# reachable by fully-qualified name (Julia modules have no real privacy).
#
# WARNING TRANSLATION -- the 4 suppressible-identifier warnings this file owns
# (`removeDuplicatePoints:dropped`, `fromSweep:noLimitingRoot`, `fromSweep:pointsDropped`,
# `fromSweep:noGreenPeak`) become plain `@warn "ModeAssociation:method:id --
# message"`, WITHOUT `maxlog` -- unlike `GreenSextic.jl`'s `nearIsotropic`/`degenEig` warnings
# (which fire inside a tight per-(layer,velocity) loop and need `maxlog = 1` to avoid spamming
# thousands of times per solve), MATLAB's own `ModeAssociation.m` source never self-silences
# these 4 (no internal `warning('off', ...)` call at any of the 4 sites), so each one
# legitimately fires once per `fromSweep`/`removeDuplicatePoints` call, matching
# `Material.jl`'s simpler (no-`maxlog`) `@warn` convention for the same reason.

# -----------------------------------------------------------------------------------------
# Constants -- ModeAssociation.m (`properties (Constant)`)
# -----------------------------------------------------------------------------------------

const LIM_VELOCITY_WINDOW = 0.15
const LIM_THETA_WINDOW    = 0.15
const LIM_BRANCH_PENALTY  = 0.5

# The binding classes tracked through the Ritz eigenproblem (rebind), as opposed to
# "limiting" (rebindLimiting). Every per-class selection in ModeAssociation /
# SensitivityAnalysis / ElasticityInversion reads this set (or activeRitzClasses): the
# sensitivity rows and the residual rows must select the same points. (MATLAB RITZ_CLASSES)
const RITZ_CLASSES = ("ritzMode", "pseudoSAW")

# Default re-identification thresholds (defaultTrackOpts) and rebind window, the one
# definition ElasticityInversion's option defaults also read. (MATLAB DEFAULT_TRACK_OPTS /
# DEFAULT_REBIND_WINDOW)
const DEFAULT_TRACK_OPTS = (macWarnThreshold = 0.7, limMatchScoreThresh = 0.5,
                            useGreenSteering = false, greenWindowRel = 0.02,
                            greenFilterRel = 0.04)
const DEFAULT_REBIND_WINDOW = 1.0   # [km/s]

"""
    activeRitzClasses(includePseudoSAW::Bool) -> Tuple

The Ritz binding classes that contribute data rows: `("ritzMode",)`, plus `"pseudoSAW"` when
`includePseudoSAW` (the option of `SensitivityAnalysis` `calc` and `ElasticityInversion`).
Port of `ModeAssociation.activeRitzClasses`.
"""
activeRitzClasses(includePseudoSAW::Bool) = includePseudoSAW ? RITZ_CLASSES : (RITZ_CLASSES[1],)

export LIM_VELOCITY_WINDOW, LIM_THETA_WINDOW, LIM_BRANCH_PENALTY, RITZ_CLASSES,
       DEFAULT_TRACK_OPTS, DEFAULT_REBIND_WINDOW, activeRitzClasses

# -----------------------------------------------------------------------------------------
# AssocPoint -- one `points` entry. ModeAssociation.m (properties block)
# (emptyPoint, the neutral-default constructor every construction path must route through).
# -----------------------------------------------------------------------------------------

"""
    AssocPoint

One `(measured point -> calculated counterpart) `binding, exactly one MATLAB `points` struct
array entry (`ModeAssociation.m`). See that file's own class docblock for the
BINDING/ANCHOR/RE-IDENTIFICATION vocabulary these fields implement.

Immutable (see this file's header) -- built via the keyword constructor (`Base.@kwdef`
defaults reproduce `ModeAssociation.emptyPoint`'s neutral values exactly), and "mutated" via
[`_updatePoint`](@ref), which returns a new `AssocPoint` with the named fields replaced.

# Fields
- `angleDeg::Float64` -- propagation angle [deg].
- `measuredVelocity::Float64` -- [km/s].
- `waveType::String` -- `MeasuredVelocities.WAVE_TYPES` member.
- `ritzModeNo::Float64` -- column index into `eigVecs_con` (stored as `Float64`, matching
  MATLAB's `double`, so the "unset" state on a `"limiting"` point can be `NaN`).
- `calculatedVelocity::Float64` -- [km/s], current best match.
- `anchorEigVec_con::Vector{ComplexF64}` -- mode-identity eigenvector; empty for `"limiting"`.
- `bindingClass::String` -- `"ritzMode"` | `"limiting"` | `"pseudoSAW"`.
- `macValue::Float64` -- achieved MAC at last rebind; `NaN` until rebound.
- `modeCollision::Bool` -- true if the last rebind could not uniquely match.
- `limBranchIwt::Float64` -- 1 = longitudinal, 2 = fast shear, 3 = slow shear (`NaN` off the
  limiting class); PART OF THE ANCHOR. `Float64` for the same NaN-sentinel reason as
  `ritzModeNo` -- convert with `Int(...)` at any call site requiring a true `Integer`
  (e.g. `LimitingVels.refineRoot`'s `iwt` argument).
- `limAnchorPol::Vector{Float64}` -- `[3]` polarisation at the tracked root; PART OF THE
  ANCHOR; empty off the limiting class.
- `limRootTheta::Float64` -- sagittal root angle theta* [rad]; NOT part of the anchor -- live
  state, refreshed on every promoted evaluation.
- `limMatchScore::Float64` -- confidence at the last re-identification, `[0,1]`.
- `limRootLost::Bool` -- no candidate root cleared the threshold.
- `anchorGreenPeakVel::Float64` -- Green peak velocity this point was bound against at
  CONSTRUCTION [km/s]; `NaN` if unavailable.
- `greenPeakVel::Float64` -- LIVE Green peak velocity at the last rescan [km/s]; `NaN`
  disables Green steering for this point.
- `greenAgreement::Float64` -- `|v_matched - greenPeakVel| / greenPeakVel` at the last
  re-identification; `NaN` if no live peak.
- `pinned::Bool` -- human-confirmed ANCHOR.
- `lowConfidenceCount::Int` -- how many PROMOTED evaluations this point spent below its
  class's confidence threshold.
- `activeThisEvaluation::Bool` -- `false` => this evaluation's residual AND Jacobian row are
  both zeroed.
"""
Base.@kwdef struct AssocPoint
    angleDeg::Float64 = NaN
    measuredVelocity::Float64 = NaN
    waveType::String = "unknown"
    ritzModeNo::Float64 = NaN
    calculatedVelocity::Float64 = NaN
    anchorEigVec_con::Vector{ComplexF64} = ComplexF64[]
    bindingClass::String = "ritzMode"
    macValue::Float64 = NaN
    modeCollision::Bool = false
    limBranchIwt::Float64 = NaN
    limAnchorPol::Vector{Float64} = Float64[]
    limRootTheta::Float64 = NaN
    limMatchScore::Float64 = NaN
    limRootLost::Bool = false
    anchorGreenPeakVel::Float64 = NaN
    greenPeakVel::Float64 = NaN
    greenAgreement::Float64 = NaN
    pinned::Bool = false
    lowConfidenceCount::Int = 0
    activeThisEvaluation::Bool = true
end

"""
    emptyPoint() -> AssocPoint

One fully-populated `AssocPoint` with every field at its neutral default. Port of
`ModeAssociation.emptyPoint` (`ModeAssociation.m`, `Access = private`) -- kept as a
thin alias (rather than inlining `AssocPoint()` everywhere) purely so this file reads
line-for-line against the MATLAB source. NOT exported (MATLAB `Access = private`).
"""
emptyPoint() = AssocPoint()

"""
    _updatePoint(p::AssocPoint; kwargs...) -> AssocPoint

Copy-with-replacements for one point (the `AssocPoint` analogue of the `_rebuild` pattern
used for the other immutable structs). Internal; not exported. `merge` first, then a
SINGLE keyword splat -- splatting two separately-keyed collections with overlapping keys
directly (`f(; a..., b...)`) is a Julia keyword-call error, not a silent override.
"""
function _updatePoint(p::AssocPoint; kwargs...)
    base = NamedTuple{fieldnames(AssocPoint)}(getfield.(Ref(p), fieldnames(AssocPoint)))
    return AssocPoint(; merge(base, NamedTuple(kwargs))...)
end

# -----------------------------------------------------------------------------------------
# ModeAssociation -- ModeAssociation.m (classdef header, properties, empty constructor)
# -----------------------------------------------------------------------------------------

"""
    ModeAssociation

The binding table: which measured velocity point corresponds to which calculated object. See
`ModeAssociation.m`'s own class docblock for the full account (reproduced in this file's
header). Immutable -- see this file's header for why, and for the confirmed value-class
grep evidence.

# Fields
- `points::Vector{AssocPoint}`
"""
struct ModeAssociation
    points::Vector{AssocPoint}
end

"""
    ModeAssociation() -> ModeAssociation

Empty `ModeAssociation`. Port of `ModeAssociation.m`.
"""
ModeAssociation() = ModeAssociation(AssocPoint[])

export ModeAssociation, AssocPoint

# -----------------------------------------------------------------------------------------
# nPoints / uniqueRitzAngles / uniqueLimitingAngles -- ModeAssociation.m
# -----------------------------------------------------------------------------------------

"""
    nPoints(obj::ModeAssociation) -> Int

Number of bound points. Port of `ModeAssociation.nPoints` (`:137-141`).
"""
nPoints(obj::ModeAssociation) = length(obj.points)

"""
    uniqueRitzAngles(obj::ModeAssociation) -> Vector{Float64}

Unique propagation angles among Ritz-mode-bound points (`"ritzMode"`/`"pseudoSAW"`), i.e. the
angles at which the inversion loop must solve an eigenproblem. Port of
`ModeAssociation.uniqueRitzAngles` (`:143-149`). MATLAB's `unique` SORTS its output ascending
by default -- `sort ∘ unique` reproduces that (Julia's `unique` preserves first-occurrence
order instead).
"""
function uniqueRitzAngles(obj::ModeAssociation)
    idx = findall(p -> p.bindingClass in RITZ_CLASSES, obj.points)
    return sort(unique(Float64[obj.points[i].angleDeg for i in idx]))
end

"""
    uniqueLimitingAngles(obj::ModeAssociation) -> Vector{Float64}

Unique propagation angles among `"limiting"`-bound points. Port of
`ModeAssociation.uniqueLimitingAngles` (`:151-161`).
"""
function uniqueLimitingAngles(obj::ModeAssociation)
    isempty(obj.points) && return Float64[]
    idx = findall(p -> p.bindingClass == "limiting", obj.points)
    return sort(unique(Float64[obj.points[i].angleDeg for i in idx]))
end

# -----------------------------------------------------------------------------------------
# pointsAtAngle
# -----------------------------------------------------------------------------------------

"""
    pointsAtAngle(assoc, angleDeg::Real, bindingClasses = nothing) -> Vector{Int}

Indices of points at a given angle, optionally filtered by `bindingClass`. Port of
`ModeAssociation.pointsAtAngle` (`:163-178`).

`assoc` is deliberately UNTYPED (duck-typed on `.points`, each a value exposing `.bindingClass`
/`.angleDeg`) -- see this file's header on why: `SensitivityAnalysis.jl` reuses this exact
function for any `ModeAssociation`-shaped value.
"""
function pointsAtAngle(assoc, angleDeg::Real,
                        bindingClasses::Union{Nothing,AbstractVector,Tuple} = nothing)
    pts = assoc.points
    isempty(pts) && return Int[]
    mask = bindingClasses === nothing ? trues(length(pts)) :
           [p.bindingClass in bindingClasses for p in pts]
    return findall(i -> abs(pts[i].angleDeg - angleDeg) < 1e-9 && mask[i], eachindex(pts))
end

# -----------------------------------------------------------------------------------------
# confidence / pinnedPoints / setPinned -- ModeAssociation.m
# -----------------------------------------------------------------------------------------

"""
    confidence(obj::ModeAssociation) -> Vector{Float64}

The class-appropriate confidence signal for every point, as one vector: `macValue` for
`"ritzMode"`/`"pseudoSAW"`, `limMatchScore` for `"limiting"`, `NaN` where never re-identified.
Port of `ModeAssociation.confidence` (`:205-224`).
"""
function confidence(obj::ModeAssociation)
    n = length(obj.points)
    c = fill(NaN, n)
    for k in 1:n
        p = obj.points[k]
        c[k] = p.bindingClass == "limiting" ? p.limMatchScore : p.macValue
    end
    return c
end

"""
    pinnedPoints(obj::ModeAssociation) -> Vector{Int}

Indices of human-confirmed (pinned) points, for reporting. Port of
`ModeAssociation.pinnedPoints` (`:226-234`).
"""
function pinnedPoints(obj::ModeAssociation)
    isempty(obj.points) && return Int[]
    return findall(p -> p.pinned, obj.points)
end

"""
    setPinned(obj::ModeAssociation, idx, tf::Bool = true) -> ModeAssociation

Mark one or more points as human-confirmed (or un-confirm them). Port of
`ModeAssociation.setPinned` (`:236-259`). `idx` -- an `Integer` or a vector thereof, 1-based
(MATLAB and Julia agree on 1-based indexing here, no shift needed).
"""
function setPinned(obj::ModeAssociation, idx, tf::Bool = true)
    idxV = idx isa AbstractVector ? idx : [idx]
    n = length(obj.points)
    all(i -> i isa Integer && 1 <= i <= n, idxV) || throw(ArgumentError(
        "ModeAssociation:setPinned:badIdx -- idx must be a vector of positive integers <= $n."))
    points = copy(obj.points)
    for k in idxV
        points[k] = _updatePoint(points[k]; pinned = tf)
    end
    return ModeAssociation(points)
end

# -----------------------------------------------------------------------------------------
# setLimitingAnchor / setLimitingRefined / setLimitingLost / setGreenPeak
# -----------------------------------------------------------------------------------------

"""
    setLimitingAnchor(obj::ModeAssociation, idx::Integer, iwt::Integer,
                       pol::AbstractVector{<:Real}, thetaRoot::Real) -> ModeAssociation

Set one limiting point's ROOT ANCHOR by hand. Port of `ModeAssociation.setLimitingAnchor`
(`:270-301`).
"""
function setLimitingAnchor(obj::ModeAssociation, idx::Integer, iwt::Integer,
                            pol::AbstractVector{<:Real}, thetaRoot::Real)
    (1 <= iwt <= 3) || throw(ArgumentError(
        "ModeAssociation:setLimitingAnchor:badIwt -- iwt must be an integer in 1:3, got $iwt."))
    (length(pol) == 3 && all(isfinite, pol)) || throw(ArgumentError(
        "ModeAssociation:setLimitingAnchor:badPol -- pol must be a finite 3-vector."))
    isfinite(thetaRoot) || throw(ArgumentError(
        "ModeAssociation:setLimitingAnchor:badTheta -- thetaRoot must be finite."))
    p = obj.points[idx]
    p.bindingClass == "limiting" || throw(ArgumentError(
        "ModeAssociation:setLimitingAnchor:wrongClass -- Point $idx has bindingClass " *
        "\"$(p.bindingClass)\"; a root anchor only means anything on \"limiting\"."))
    points = copy(obj.points)
    points[idx] = _updatePoint(p; limBranchIwt = Float64(iwt), limAnchorPol = Float64.(vec(pol)),
                                   limRootTheta = Float64(thetaRoot), limRootLost = false)
    return ModeAssociation(points)
end

"""
    setLimitingRefined(obj::ModeAssociation, idx::Integer, velocity::Real, thetaRoot::Real)
        -> ModeAssociation

Advance ONE limiting point from a successful anchored local refinement, without a full
re-identification. Port of `ModeAssociation.setLimitingRefined` (`:303-320`).
"""
function setLimitingRefined(obj::ModeAssociation, idx::Integer, velocity::Real, thetaRoot::Real)
    points = copy(obj.points)
    points[idx] = _updatePoint(points[idx]; calculatedVelocity = Float64(velocity),
                                             limRootTheta = Float64(thetaRoot),
                                             limRootLost = false, activeThisEvaluation = true)
    return ModeAssociation(points)
end

"""
    setLimitingLost(obj::ModeAssociation, idx::Integer) -> ModeAssociation

Flag ONE limiting point's tracked root as lost at this evaluation. Port of
`ModeAssociation.setLimitingLost`: written by ElasticityInversion when a perturbed evaluation
of the point's finite-difference Jacobian row lost the root (`channelDerivative`'s
`ok = false`); the point is excluded from this evaluation exactly as `rebindLimiting`
excludes a vanished root.
"""
function setLimitingLost(obj::ModeAssociation, idx::Integer)
    points = copy(obj.points)
    p = points[idx]
    points[idx] = _updatePoint(p; limRootLost = true, activeThisEvaluation = false,
                                  lowConfidenceCount = p.lowConfidenceCount + 1)
    return ModeAssociation(points)
end

"""
    setGreenPeak(obj::ModeAssociation, idx::Integer, peakVel::Real) -> ModeAssociation

Store a freshly measured live Green peak velocity on one point. Port of
`ModeAssociation.setGreenPeak` (`:322-337`).
"""
function setGreenPeak(obj::ModeAssociation, idx::Integer, peakVel::Real)
    points = copy(obj.points)
    points[idx] = _updatePoint(points[idx]; greenPeakVel = Float64(peakVel))
    return ModeAssociation(points)
end

# -----------------------------------------------------------------------------------------
# addPoint / removeDuplicatePoints
# -----------------------------------------------------------------------------------------

"""
    addPoint(obj::ModeAssociation, angleDeg, ritzModeNo, measuredVelocity, calculatedVelocity,
              anchorEigVec_con, bindingClass = "ritzMode", waveType = "unknown",
              pinned = false) -> ModeAssociation

Add one measured-point <-> calculated-object binding. Port of `ModeAssociation.addPoint`
(`:339-365`).
"""
function addPoint(obj::ModeAssociation, angleDeg::Real, ritzModeNo::Real, measuredVelocity::Real,
                   calculatedVelocity::Real, anchorEigVec_con::AbstractVector,
                   bindingClass::AbstractString = "ritzMode",
                   waveType::AbstractString = "unknown", pinned::Bool = false)
    bindingClass in (RITZ_CLASSES..., "limiting") || throw(ArgumentError(
        "ModeAssociation:addPoint:badBindingClass -- bindingClass must be 'ritzMode', " *
        "'pseudoSAW' or 'limiting'; got '$bindingClass'."))
    waveType in WAVE_TYPES || throw(ArgumentError("ModeAssociation:addPoint:badWaveType -- " *
        "\"$waveType\" is not one of: $(join(WAVE_TYPES, ", "))."))
    newPt = AssocPoint(; angleDeg = Float64(angleDeg), measuredVelocity = Float64(measuredVelocity),
                          waveType = String(waveType), ritzModeNo = Float64(ritzModeNo),
                          calculatedVelocity = Float64(calculatedVelocity),
                          anchorEigVec_con = ComplexF64.(anchorEigVec_con),
                          bindingClass = String(bindingClass), pinned = pinned)
    return ModeAssociation(vcat(obj.points, [newPt]))
end

"""
    removeDuplicatePoints(obj::ModeAssociation) -> ModeAssociation

Keep at most one binding per `(angle, ritzModeNo)` for `"ritzMode"`/`"pseudoSAW"` points, and
at most one per `(angle, waveType)` for `"limiting"` points -- the closest-matching duplicate
wins, the rest are dropped with a warning. Port of `ModeAssociation.removeDuplicatePoints`
(`:373-412`). `unique` with 3 outputs (MATLAB `:391`) -> manual `Dict`-based grouping by a
plain tuple key -- see this file's header.
"""
function removeDuplicatePoints(obj::ModeAssociation)
    n = length(obj.points)
    n < 2 && return obj

    groupOf = Dict{Tuple{Float64,String},Vector{Int}}()
    for k in 1:n
        p = obj.points[k]
        key = p.bindingClass == "limiting" ? (round(p.angleDeg; digits = 6), p.waveType) :
                                              (round(p.angleDeg; digits = 6), string(Int(p.ritzModeNo)))
        push!(get!(groupOf, key, Int[]), k)
    end

    keep = trues(n)
    nDropped = 0
    for members in values(groupOf)
        length(members) < 2 && continue
        errs = [abs(obj.points[m].calculatedVelocity - obj.points[m].measuredVelocity) for m in members]
        iBest = argmin(errs)
        for (j, m) in enumerate(members)
            if j != iBest
                keep[m] = false
                nDropped += 1
            end
        end
    end

    if nDropped > 0
        @warn "ModeAssociation:removeDuplicatePoints:dropped -- $nDropped duplicate point(s) " *
              "removed (same angle+mode bound more than once)."
    end
    return ModeAssociation(obj.points[keep])
end

# -----------------------------------------------------------------------------------------
# rebind -- ModeAssociation.m
# -----------------------------------------------------------------------------------------

"""
    rebind(obj::ModeAssociation, angleDeg::Real, eigVecs_con::AbstractMatrix{<:Number},
           M_con::AbstractMatrix{<:Number}, vels::AbstractVector{<:Real},
           rebindVelocityWindow::Real = 1.0, updateAnchor::Bool = false,
           trackOpts::Union{Nothing,NamedTuple} = nothing) -> ModeAssociation

Re-establish the RITZ-class mode binding at one angle from a freshly solved eigenproblem.
Port of `ModeAssociation.rebind` (`:428-607`) -- see that method's own extensive MATLAB
comments (M-weighted MAC, Green steering, pinned-priority) for the physics/design reasoning,
reproduced faithfully below without re-deriving it here.

`M_con` may be real or complex (`AbstractMatrix{<:Number}`) -- MATLAB's own `'` on
`eigVecs_con` is a genuine conjugate transpose (the eigenvectors are complex); Julia's
`adjoint` (`'`) does the same and additionally reduces to a plain transpose when the operand
is real, matching MATLAB exactly either way.
"""
function rebind(obj::ModeAssociation, angleDeg::Real, eigVecs_con::AbstractMatrix{<:Number},
                 M_con::AbstractMatrix{<:Number}, vels::AbstractVector{<:Real},
                 rebindVelocityWindow::Real = DEFAULT_REBIND_WINDOW, updateAnchor::Bool = false,
                 trackOpts::Union{Nothing,NamedTuple} = nothing)
    opts = mergeTrackOpts(trackOpts)

    idxPts = pointsAtAngle(obj, angleDeg, RITZ_CLASSES)
    isempty(idxPts) && return obj

    nModes = length(vels)
    nPts = length(idxPts)
    MAC = zeros(nModes, nPts)
    score = zeros(nModes, nPts)
    steered = falses(nPts)

    for k in 1:nPts
        p = obj.points[idxPts[k]]
        useGreen = opts.useGreenSteering && p.bindingClass == "pseudoSAW" &&
                   isfinite(p.greenPeakVel) && p.greenPeakVel > 0
        steered[k] = useGreen
        centre = useGreen ? p.greenPeakVel : p.calculatedVelocity
        inWindow = abs.(vels .- centre) .< rebindVelocityWindow
        if useGreen
            # HARD candidate filter at the GROSS-ERROR scale, MAC free to choose inside it --
            # see ModeAssociation.m's own comment for the measured reasoning.
            inWindow = inWindow .& (abs.(vels .- p.greenPeakVel) .<= opts.greenFilterRel * p.greenPeakVel)
        end
        !any(inWindow) && continue

        # M * anchor first: O(n^2) per point instead of forming V' * M (O(n^2 k))
        macCol = abs.(view(eigVecs_con, :, inWindow)' * (M_con * p.anchorEigVec_con))
        MAC[inWindow, k] .= macCol
        if useGreen
            # inside the filter, the Gaussian only breaks ties
            w = opts.greenWindowRel * p.greenPeakVel
            score[inWindow, k] .= macCol .* exp.(-((vels[inWindow] .- p.greenPeakVel) ./ w) .^ 2)
        else
            score[inWindow, k] .= macCol
        end
    end

    pinnedMask = [obj.points[i].pinned for i in idxPts]
    assignedMode = assignWithPinnedPriority(score, pinnedMask)

    points = copy(obj.points)
    for k in 1:nPts
        gi = idxPts[k]
        p = points[gi]
        m = assignedMode[k]

        matched = !isnan(m) && score[Int(m), k] > 0

        if !matched
            # ritzModeNo/calculatedVelocity deliberately LEFT at their previous values (no
            # better estimate); the row is zeroed via activeThisEvaluation so a stale velocity
            # can never pair with a current-eigenvector Jacobian row.
            p = _updatePoint(p; macValue = 0.0, modeCollision = true,
                                 lowConfidenceCount = p.lowConfidenceCount + 1,
                                 activeThisEvaluation = p.pinned)
        else
            mi = Int(m)
            greenAgreement = p.greenAgreement
            if isfinite(p.greenPeakVel) && p.greenPeakVel > 0
                greenAgreement = abs(vels[mi] - p.greenPeakVel) / p.greenPeakVel
            end
            # A STEERED point's MAC is a diagnostic, not the gate (ModeAssociation.m):
            # greenAgreement is already bounded by the hard filter above.
            lowConf = steered[k] ? false : (MAC[mi, k] < opts.macWarnThreshold)
            newAnchor = (updateAnchor && !p.pinned && !lowConf) ? eigVecs_con[:, mi] : p.anchorEigVec_con
            p = _updatePoint(p; ritzModeNo = Float64(mi), calculatedVelocity = vels[mi],
                                 macValue = MAC[mi, k], modeCollision = false,
                                 greenAgreement = greenAgreement,
                                 lowConfidenceCount = lowConf ? p.lowConfidenceCount + 1 :
                                                                 p.lowConfidenceCount,
                                 activeThisEvaluation = lowConf ? p.pinned : true,
                                 anchorEigVec_con = newAnchor)
        end
        points[gi] = p
    end
    return ModeAssociation(points)
end

# -----------------------------------------------------------------------------------------
# rebindLimiting -- ModeAssociation.m
# -----------------------------------------------------------------------------------------

"""
    rebindLimiting(obj::ModeAssociation, angleDeg::Real, roots::AbstractVector{LimitingRoot},
                    updateAnchor::Bool = false,
                    trackOpts::Union{Nothing,NamedTuple} = nothing) -> ModeAssociation

Re-establish the LIMITING-class root binding at one angle from a freshly computed Christoffel
root set. Port of `ModeAssociation.rebindLimiting` (`:609-692`) -- a deliberate sibling of
[`rebind`](@ref) (same shape, same `matchpairs`/pinned-priority pass), scoring candidates by
[`limitingRootScore`](@ref) instead of MAC.
"""
function rebindLimiting(obj::ModeAssociation, angleDeg::Real, roots::AbstractVector{LimitingRoot},
                         updateAnchor::Bool = false,
                         trackOpts::Union{Nothing,NamedTuple} = nothing)
    opts = mergeTrackOpts(trackOpts)

    idxPts = pointsAtAngle(obj, angleDeg, ["limiting"])
    isempty(idxPts) && return obj

    nRoots = length(roots)
    nPts = length(idxPts)
    score = zeros(nRoots, nPts)
    for k in 1:nPts
        if nRoots > 0
            score[:, k] = limitingRootScore(obj.points[idxPts[k]], roots)
        end
    end

    pinnedMask = [obj.points[i].pinned for i in idxPts]
    assignedRoot = assignWithPinnedPriority(score, pinnedMask)

    points = copy(obj.points)
    for k in 1:nPts
        gi = idxPts[k]
        p = points[gi]
        q = assignedRoot[k]

        if isnan(q) || score[Int(q), k] <= opts.limMatchScoreThresh
            newMatchScore = isnan(q) ? 0.0 : score[Int(q), k]
            # calculatedVelocity/limRootTheta/limBranchIwt deliberately left at their previous
            # values -- a vanished root must NOT silently take the best remaining one.
            p = _updatePoint(p; limMatchScore = newMatchScore, limRootLost = true,
                                 lowConfidenceCount = p.lowConfidenceCount + 1,
                                 activeThisEvaluation = p.pinned)
        else
            qi = Int(q)
            newAnchorPol = (updateAnchor && !p.pinned) ? Float64.(vec(roots[qi].pol)) : p.limAnchorPol
            p = _updatePoint(p; calculatedVelocity = roots[qi].vel, limMatchScore = score[qi, k],
                                 limRootTheta = roots[qi].theta, limBranchIwt = Float64(roots[qi].iwt),
                                 limRootLost = false, activeThisEvaluation = true,
                                 limAnchorPol = newAnchorPol)
        end
        points[gi] = p
    end
    return ModeAssociation(points)
end

# -----------------------------------------------------------------------------------------
# fromSweep -- ModeAssociation.m
# -----------------------------------------------------------------------------------------

"""
    fromSweep(sweep, measured::MeasuredVelocities, velocityTolerance::Real = 0.3,
              angleTolerance::Real = 1e-6, invertedLayerNo::Integer = 1) -> ModeAssociation

Build the initial mode association from a forward sweep at the initial elasticity guess. Port
of `ModeAssociation.fromSweep` (`:696-841`).

`sweep` is deliberately UNTYPED (duck-typed: `CalcAngSweepHS`/`CalcAngSweepLayered`-shaped,
matching MATLAB's own docstring) -- both share `.structure`, `.angles`, `.significantModesMap`,
`.cs[i].ritzSolver.eigVecs_con`.

Wave-type routing: see [`classifyWaveType`](@ref). The two classes bind DIFFERENTLY -- Ritz
classes filter `sweep.significantModesMap` by the accepted `ModeClassification` types and take
the entry minimising `|vel - measured|`; the limiting class computes the Christoffel root set
at that angle directly, filters by branch, and takes the nearest by velocity (never routed
through `significantModesMap`, which would let every limiting branch compete against every
'limiting' entry and apply `velocityTolerance` to a nearby Ritz mode instead of the fitted
root).
"""
function fromSweep(sweep, measured::MeasuredVelocities, velocityTolerance::Real = 0.3,
                    angleTolerance::Real = 1e-6, invertedLayerNo::Integer = 1)
    nDropped = 0
    nNoGreen = 0   # pseudoSAW points bound without a Green peak
    pts = AssocPoint[]

    matLay = sweep.structure.mater[invertedLayerNo]
    symmNo = matLay.elConsts.symmClassNo
    cVec0 = Float64.(vec(matLay.elConsts.vector))
    limRootsCache = Dict{Int,Vector{LimitingRoot}}()

    for i in 1:measured.nPoints
        angleDeg = measured.angleDeg[i]
        dAngle, iAngle = findmin(abs.(sweep.angles .- angleDeg))
        if dAngle > angleTolerance
            throw(ArgumentError("ModeAssociation:fromSweep:angleNotFound -- Measured angle " *
                "$(angleDeg) deg has no sweep angle within tolerance $(angleTolerance) deg " *
                "(closest sweep angle: $(sweep.angles[iAngle]) deg)."))
        end

        wtMeas = measured.waveType[i]
        acceptedTypes, bindingClass = classifyWaveType(wtMeas)

        if bindingClass == "limiting"
            roots = get!(limRootsCache, iAngle) do
                limitingRootsAt(matLay, symmNo, cVec0, sweep.angles[iAngle])
            end
            iwtWanted = waveTypeIndex(wtMeas)
            cand = filter(r -> r.iwt == iwtWanted, roots)
            if isempty(cand)
                # The branch can have no root at this angle at all; never fall back to a root
                # of the OTHER branch.
                @warn "ModeAssociation:fromSweep:noLimitingRoot -- No $wtMeas limiting root " *
                      "exists at $angleDeg deg; the point was dropped rather than bound to " *
                      "another branch."
                nDropped += 1
                continue
            end
            errs = [abs(c.vel - measured.velocity[i]) for c in cand]
            dv, iBest = findmin(errs)
            if dv > velocityTolerance
                nDropped += 1
                continue
            end
            best = cand[iBest]
            newPt = AssocPoint(; angleDeg = sweep.angles[iAngle], measuredVelocity = measured.velocity[i],
                                  waveType = wtMeas, bindingClass = bindingClass,
                                  calculatedVelocity = best.vel, ritzModeNo = NaN,
                                  anchorEigVec_con = ComplexF64[],
                                  limBranchIwt = Float64(best.iwt), limAnchorPol = Float64.(vec(best.pol)),
                                  limRootTheta = best.theta, limMatchScore = NaN, limRootLost = false)
        else
            candidates = sweep.significantModesMap[iAngle]
            if !isempty(acceptedTypes)
                candidates = filter(c -> c.type in acceptedTypes, candidates)
            end
            if isempty(candidates)
                nDropped += 1
                continue
            end
            # modeVel: the bound Ritz mode's eigenvelocity, what the inversion fits (an
            # entry's vel is the Green peak)
            errs = [abs(c.modeVel - measured.velocity[i]) for c in candidates]
            dv, iBest = findmin(errs)
            if dv > velocityTolerance
                nDropped += 1
                continue
            end
            best = candidates[iBest]

            greenAgreement = NaN
            if isfinite(best.greenPeakVel) && best.greenPeakVel > 0
                greenAgreement = abs(best.modeVel - best.greenPeakVel) / best.greenPeakVel
            elseif bindingClass == "pseudoSAW"
                nNoGreen += 1
            end

            # An INDEPENDENT, Green-derived estimate of where the physical resonance sits,
            # not subject to the near-degenerate-cluster problem the MAC has.
            newPt = AssocPoint(; angleDeg = sweep.angles[iAngle], measuredVelocity = measured.velocity[i],
                                  waveType = wtMeas, bindingClass = bindingClass,
                                  calculatedVelocity = best.modeVel, ritzModeNo = Float64(best.modeNo),
                                  anchorEigVec_con = sweep.cs[iAngle].ritzSolver.eigVecs_con[:, best.modeNo],
                                  anchorGreenPeakVel = best.greenPeakVel,
                                  greenPeakVel = best.greenPeakVel,
                                  greenAgreement = greenAgreement)
        end

        push!(pts, newPt)
    end

    if nDropped > 0
        @warn "ModeAssociation:fromSweep:pointsDropped -- $nDropped of $(measured.nPoints) " *
              "measured point(s) had no matching mode within tolerance and were dropped."
    end
    if nNoGreen > 0
        # A sweep classified without Green (the 'ritz' method) has no peak to seed the
        # steering with: the first rescan would then start from the point's own Ritz
        # velocity, the mirror of the Ritz tracking that steering exists to avoid.
        @warn "ModeAssociation:fromSweep:noGreenPeak -- $nNoGreen pseudoSAW point(s) were " *
              "bound to entries without a Green peak (a sweep classified without Green, e.g. " *
              "the 'ritz' method): Green steering has no seed for them. Classify the sweep " *
              "with 'green' for a pseudo-SAW inversion."
    end

    return removeDuplicatePoints(ModeAssociation(pts))
end

# -----------------------------------------------------------------------------------------
# defaultTrackOpts / mergeTrackOpts -- ModeAssociation.m
# -----------------------------------------------------------------------------------------

"""
    defaultTrackOpts() -> NamedTuple

Default re-identification thresholds shared by [`rebind`](@ref)/[`rebindLimiting`](@ref).
Port of `ModeAssociation.defaultTrackOpts` (`:843-880`) -- a `NamedTuple` rather than a
mutable MATLAB struct; override one field with `merge(defaultTrackOpts(), (useGreenSteering =
true,))` rather than MATLAB's `opts.useGreenSteering = true` dot-assignment.

# Fields
- `macWarnThreshold::Float64 = 0.7`
- `limMatchScoreThresh::Float64 = 0.5`
- `useGreenSteering::Bool = false`
- `greenWindowRel::Float64 = 0.02`
- `greenFilterRel::Float64 = 0.04`
"""
defaultTrackOpts() = DEFAULT_TRACK_OPTS

"""
    mergeTrackOpts(userOpts) -> NamedTuple

Fill any field the caller omitted from [`defaultTrackOpts`](@ref); errors on an unrecognised
field name. Port of `ModeAssociation.mergeTrackOpts` (`:966-981`, `Access = private`; NOT
exported). `userOpts` -- `nothing` or a `NamedTuple` (partial or full).
"""
function mergeTrackOpts(userOpts)
    opts = defaultTrackOpts()
    userOpts === nothing && return opts
    for (k, v) in pairs(userOpts)
        haskey(opts, k) || throw(ArgumentError("ModeAssociation:mergeTrackOpts:unknownField -- " *
            "\"$k\" is not a recognised tracking option."))
        opts = merge(opts, NamedTuple{(k,)}((v,)))
    end
    return opts
end

# -----------------------------------------------------------------------------------------
# limitingRootScore -- ModeAssociation.m
# -----------------------------------------------------------------------------------------

"""
    limitingRootScore(pt, roots::AbstractVector{LimitingRoot}) -> Vector{Float64}

Re-identification score of every candidate root against ONE point's anchor, `[0,1]`. Port of
`ModeAssociation.limitingRootScore` (`:882-928`):

    score = |pol_a' pol_k| / (|pol_a| |pol_k|)          polarisation overlap
          * exp(-((v_k - v_a) / (LIM_VELOCITY_WINDOW*v_a))^2)
          * exp(-((theta_k - theta_a) / LIM_THETA_WINDOW)^2)
          * (iwt_k == iwt_a ? 1 : LIM_BRANCH_PENALTY)

`pt` is duck-typed on `.limAnchorPol`/`.calculatedVelocity`/`.limRootTheta`/`.limBranchIwt` --
both a real `AssocPoint` and `SensitivityAnalysis.jl`'s `limitingVelocityTracked` fallback
pass a value exposing exactly these fields.
"""
function limitingRootScore(pt, roots::AbstractVector{LimitingRoot})
    n = length(roots)
    score = zeros(n)
    polA = pt.limAnchorPol
    nA = norm(polA)
    (nA == 0 || isempty(polA)) && return score
    vA = pt.calculatedVelocity
    thA = pt.limRootTheta
    iwtA = pt.limBranchIwt
    wV = LIM_VELOCITY_WINDOW * abs(vA)
    for q in 1:n
        polK = roots[q].pol
        nK = norm(polK)
        nK == 0 && continue
        overlap = abs(polA' * polK) / (nA * nK)
        fVel = exp(-((roots[q].vel - vA) / wV)^2)
        fTheta = exp(-((roots[q].theta - thA) / LIM_THETA_WINDOW)^2)
        fBranch = roots[q].iwt == iwtA ? 1.0 : LIM_BRANCH_PENALTY
        score[q] = overlap * fVel * fTheta * fBranch
    end
    return score
end

# -----------------------------------------------------------------------------------------
# assignWithPinnedPriority -- ModeAssociation.m (Access = private; NOT exported)
# -----------------------------------------------------------------------------------------

"""
    assignWithPinnedPriority(score::AbstractMatrix{<:Real}, pinnedMask::AbstractVector{Bool})
        -> Vector{Float64}

One-to-one candidate assignment, pinned points claiming first. Port of
`ModeAssociation.assignWithPinnedPriority` (`:983-1025`, `Access = private`; NOT exported --
see this file's header). `score` -- `[nCandidates x nPoints]`, `>= 0`, `0` = ineligible.
`assigned` -- `[nPoints]`, candidate row index (as `Float64`, MATLAB's `NaN`-sentinel
convention) or `NaN`.

Pinned points are removed from the assignment problem first, by their own argmax over the
still-free candidates -- without this, an unpinned point with a marginally higher score could
displace a human's confirmed binding onto its second choice. The remainder is solved by
[`matchpairs0`](@ref) (`MatlabCompat`) on the NEGATED score
(`matchpairs0` minimises; maximising total score <=> minimising `-score`) -- not
reimplemented here.
"""
function assignWithPinnedPriority(score::AbstractMatrix{<:Real}, pinnedMask::AbstractVector{Bool})
    nCand, nPts = size(score)
    assigned = fill(NaN, nPts)
    (nCand == 0 || nPts == 0) && return assigned

    freeRows = trues(nCand)
    for k in findall(pinnedMask)
        col = copy(score[:, k])
        col[.!freeRows] .= 0
        m = argmax(col)
        if col[m] > 0
            assigned[k] = m
            freeRows[m] = false
        end
    end

    remaining = findall(.!pinnedMask)
    (isempty(remaining) || !any(freeRows)) && return assigned

    rowMap = findall(freeRows)
    subScore = score[freeRows, remaining]
    assignment = matchpairs0(-subScore)   # maximise total score <=> minimise -score
    for (rowIdxInSub, colIdxInSub) in assignment
        assigned[remaining[colIdxInSub]] = rowMap[rowIdxInSub]
    end
    return assigned
end

# -----------------------------------------------------------------------------------------
# classifyWaveType -- ModeAssociation.m (Access = private; NOT exported)
# -----------------------------------------------------------------------------------------

"""
    classifyWaveType(waveType::AbstractString) -> (Vector{Symbol}, String)

Map a `MeasuredVelocities` wave-type label onto `ModeClassification` `type` names and a
`ModeAssociation` `bindingClass`. Port of `ModeAssociation.classifyWaveType` (`:1027-1056`,
`Access = private`; NOT exported).

`acceptedTypes` is `Vector{Symbol}` (not `String`), matching the port's own established
`SignificantMode.type::Symbol` convention (`:Rayleigh`/`:SAW`/`:pseudoRayleigh`/... --
`ModeClassification.jl`), meaningful for the Ritz classes only: the limiting class does not
route through `significantModesMap` at all (see [`fromSweep`](@ref)), so its `[:limiting]`
entry is kept for documentation value rather than being used as a filter.
"""
function classifyWaveType(waveType::AbstractString)
    if waveType == "rayleigh"
        return Symbol[:Rayleigh, :SAW], "ritzMode"
    elseif waveType in ("sezawa", "love")
        return Symbol[:SAW], "ritzMode"
    elseif waveType == "pseudoRayleigh"
        return Symbol[:pseudoRayleigh], "pseudoSAW"
    elseif waveType in ("longitudinal", "fastShear", "slowShear")
        return Symbol[:limiting], "limiting"
    elseif waveType == "unknown"
        return Symbol[], "ritzMode"
    else
        throw(ArgumentError("ModeAssociation:classifyWaveType:badWaveType -- \"$waveType\" " *
            "is not a recognised wave type."))
    end
end

export nPoints, uniqueRitzAngles, uniqueLimitingAngles, pointsAtAngle, confidence,
       pinnedPoints, setPinned, setLimitingAnchor, setLimitingRefined, setLimitingLost,
       setGreenPeak, addPoint, removeDuplicatePoints, rebind, rebindLimiting,
       fromSweep, defaultTrackOpts, limitingRootScore
