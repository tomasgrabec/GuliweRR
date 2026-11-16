# =========================================================================================
# CalcIsofreqAngSweepLayered.jl -- port of guliwerr/forwardCalcs/CalcIsofreqAngSweepLayered.m:
# angular sweep of CalcIsofreqSingleLayered objects at FIXED reference frequency omega0.
# Rotates elasticity per angle and runs one CalcIsofreqSingleLayered (fixed omega0, internal
# sweep over k_par) at each, delivering raw per-angle isoFreq spectra -- NOT peak-tracked into connected contours.
#
# NO per-file module (see Calc.jl's header for why).
#
# MUTABLE: a `mutable struct`, like CalcSingleHS/Layered and their sweeps
# (MATLAB handle semantics): every step mutates the object and returns it, so both
# `calc(x)` and `x = calc(x)` work; `copy(x)` is a shallow copy.
#
# NO PHASE SPLIT (as in CalcAngSweepPlate): no
# Ritz and no GammaMatrix here at all, so there is no analogue of a shared cross-angle
# reduction -- every angle's `CalcIsofreqSingleLayered.calc()` is fully independent (fixed
# `omega0`/`kParRange`, shared across every angle by construction, not derived from any
# phase's results). `calc()` is ONE flat `forEachElement` call.
# =========================================================================================

"""
    CalcIsofreqAngSweepLayered

Angular sweep of [`CalcIsofreqSingleLayered`](@ref) objects at a fixed reference frequency
`omega0`.

Port of `guliwerr/forwardCalcs/CalcIsofreqAngSweepLayered.m`. **Mutable struct** (MATLAB handle semantics). Single flat loop, no phase split, no
retention policy -- every angle's calculation is fully independent.

# Fields
- `structure::Structure` -- the UNROTATED base structure; any `nlay >= 1` (`GreenSolver` is
  layered-generic).
- `angles::Vector{Float64}` -- propagation angles [deg]; default `0:2:180`.
- `omega0::Float64` -- fixed reference frequency; `NaN` (unset) only on the empty/default
  object -- REQUIRED at construction otherwise.
- `kParRange::Vector{Float64}` -- wavevector-magnitude axis to sweep, shared across every
  angle; REQUIRED at construction otherwise.
- `useParallel::Symbol` -- `:auto` (default) | `:always` | `:never`.
- `cs::Vector{CalcIsofreqSingleLayered}` -- one per angle.
- `absGMap::Matrix{Float64}` -- `[nK x nAngles]`, filled by [`calc`](@ref).
"""
mutable struct CalcIsofreqAngSweepLayered
    structure::Structure
    angles::Vector{Float64}
    omega0::Float64
    kParRange::Vector{Float64}
    useParallel::Symbol

    cs::Vector{CalcIsofreqSingleLayered}

    absGMap::Matrix{Float64}
end

"""
    CalcIsofreqAngSweepLayered(; structure = nothing, angles = nothing, omega0 = nothing,
                                  kParRange = nothing, useParallel = :auto) ->
        CalcIsofreqAngSweepLayered

Construct `CalcIsofreqAngSweepLayered` for an angular sweep at fixed `omega0`: one
`CalcIsofreqSingleLayered` per angle, all sharing the SAME `omega0`/`kParRange`.

Port of `CalcIsofreqAngSweepLayered.m`. `structure === nothing` (MATLAB's `nargin == 0`)
returns a well-formed empty object (`angles` still defaults to `0:2:180`, `omega0` stays
`NaN`, `kParRange`/`cs`/`absGMap` stay empty).

Otherwise `omega0` and `kParRange` are BOTH required, with no default (mirrors
`CalcIsofreqSingleLayered`'s own constructor contract exactly):
- `omega0 === nothing` throws `ArgumentError`
  (`"CalcIsofreqAngSweepLayered:constructor:missingOmega0"`); otherwise must be real, finite,
  positive.
- `kParRange === nothing` or empty throws `ArgumentError`
  (`"CalcIsofreqAngSweepLayered:constructor:missingKParRange"`); otherwise must be a real,
  finite, positive vector.

Each `cs[i]` is built by `CalcIsofreqSingleLayered`'s own constructor, which does the
per-angle rotation itself -- this class never rotates anything itself, and has no Ritz/
GammaMatrix precompute to share (unlike the HS, layered and plate angle sweeps).

`useParallel` is validated eagerly (mirrors MATLAB's `mustBeMember` property validator).

Keyword-only, not positional -- see `CalcSingleHS`'s constructor for why.
"""
function CalcIsofreqAngSweepLayered(; structure::Union{Structure,Nothing} = nothing,
                                       angles = nothing,
                                       omega0::Union{Real,Nothing} = nothing,
                                       kParRange = nothing,
                                       useParallel::Symbol = :auto)
    _validate_use_parallel(useParallel, "CalcIsofreqAngSweepLayered:constructor")

    angs = (angles === nothing || isempty(angles)) ? _ANGSWEEP_DEFAULT_ANGLES :
           Float64.(vec(collect(angles)))

    if structure === nothing
        return CalcIsofreqAngSweepLayered(Structure(), angs, NaN, Float64[], useParallel,
                                           CalcIsofreqSingleLayered[], zeros(Float64, 0, 0))
    end

    omega0 === nothing &&
        throw(ArgumentError("CalcIsofreqAngSweepLayered:constructor:missingOmega0 -- omega0 " *
            "is required (no default) -- a CalcIsofreqAngSweepLayered with no frequency " *
            "reference is meaningless."))
    (isfinite(omega0) && omega0 > 0) ||
        throw(ArgumentError("CalcIsofreqAngSweepLayered: omega0 must be a real, finite, " *
            "positive scalar; got $omega0"))

    (kParRange === nothing || isempty(kParRange)) &&
        throw(ArgumentError("CalcIsofreqAngSweepLayered:constructor:missingKParRange -- " *
            "kParRange is required (no default) -- the wavevector-magnitude axis to sweep, " *
            "shared across every angle."))
    kpr = Float64.(vec(collect(kParRange)))
    (all(isfinite, kpr) && all(>(0), kpr)) ||
        throw(ArgumentError("CalcIsofreqAngSweepLayered: kParRange must be a real, finite, " *
            "positive vector"))

    nAngles = length(angs)
    cs = Vector{CalcIsofreqSingleLayered}(undef, nAngles)
    for iAngle in 1:nAngles
        cs[iAngle] = CalcIsofreqSingleLayered(; structure = structure, angle = angs[iAngle],
                                                 omega0 = Float64(omega0), kParRange = kpr)
    end

    return CalcIsofreqAngSweepLayered(structure, angs, Float64(omega0), kpr, useParallel,
                                       cs, zeros(Float64, 0, 0))
end

"""
    _rebuild(sweep::CalcIsofreqAngSweepLayered; <field> = sweep.<field>, ...) ->
        CalcIsofreqAngSweepLayered

In-place field update (the Calc types are mutable, MATLAB handle semantics). Internal; not exported.
"""
function _rebuild(x::CalcIsofreqAngSweepLayered; kw...)
    for (k, v) in kw
        setproperty!(x, k, v)
    end
    return x
end

"""
    copy(x::CalcIsofreqAngSweepLayered)

Shallow copy (MATLAB `copy` of a handle object, one level).
"""
Base.copy(x::CalcIsofreqAngSweepLayered) = CalcIsofreqAngSweepLayered((getfield(x, k) for k in fieldnames(CalcIsofreqAngSweepLayered))...)

# -----------------------------------------------------------------------------------------
# calc -- CalcIsofreqAngSweepLayered.m
# -----------------------------------------------------------------------------------------

"""
    calc(sweep::CalcIsofreqAngSweepLayered) -> CalcIsofreqAngSweepLayered

Run the isoFreq sweep: one `CalcIsofreqSingleLayered.calc()` per angle, then assemble the
convenience map.

Port of `CalcIsofreqAngSweepLayered.calc`. No phase
split, no shared-axis preamble, no retention policy -- every angle's calculation is fully
independent, so this is a single flat loop.
"""
function calc(sweep::CalcIsofreqAngSweepLayered)
    useP = resolveParallel(sweep.useParallel, length(sweep.cs))
    cs = forEachElement(sweep.cs, useP, (cs0, _) -> calc(cs0))
    sweep = _rebuild(sweep; cs = cs)
    return assembleMaps(sweep)
end

# -----------------------------------------------------------------------------------------
# assembleMaps -- CalcIsofreqAngSweepLayered.m (MATLAB `Access = private`; not
# exported)
# -----------------------------------------------------------------------------------------

"""
    assembleMaps(sweep::CalcIsofreqAngSweepLayered) -> CalcIsofreqAngSweepLayered

Extract the sweep-level convenience map from `sweep.cs`.

Port of `CalcIsofreqAngSweepLayered.assembleMaps`.
"""
function assembleMaps(sweep::CalcIsofreqAngSweepLayered)
    # collect(...) first, not a bare generator: reduce(hcat, <length-1 generator>) returns the
    # lone element UNCHANGED (a Vector, not a 1-column Matrix), which would break a
    # single-angle sweep (as in a single-angle CalcFreqAngSweep2Layered).
    absGMap = reduce(hcat, [c.absG for c in sweep.cs])
    return _rebuild(sweep; absGMap = absGMap)
end
