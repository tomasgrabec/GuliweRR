# =========================================================================================
# CalcFreqAngSweep2Layered.jl -- port of guliwerr/forwardCalcs/CalcFreqAngSweep2Layered.m:
# 2D sweep (propagation angle AND real reference frequency omega0) for the layered forward
# calculation, isoFreq side. The counterpart to CalcWavevecAngSweep2Layered.jl: constructs
# one CalcIsofreqAngSweepLayered per omega0 value and drives it through its own `calc`.
#
# NO per-file module (see Calc.jl's header for why).
#
# MUTABLE: a `mutable struct`, like CalcSingleHS/Layered and their sweeps
# (MATLAB handle semantics): every step mutates the object and returns it, so both
# `calc(x)` and `x = calc(x)` work; `copy(x)` is a shallow copy.
#
# DELETED, NOT PORTED: `flattenSweep` / `calcFlat`. Same reasoning as
# `CalcWavevecAngSweep2Layered.jl`'s header, and if anything an EASIER case: this family has
# NO cross-angle dependency at all (no phase split, no shared axis -- every
# `CalcIsofreqSingleLayered.calc()` at a given (angle, omega0) is fully independent, per
# `CalcIsofreqAngSweepLayered.jl`'s own header), so the MATLAB source's own comment
# calls its `flattenSweep=true` path "mechanical" (`CalcFreqAngSweep2Layered.m`: "This
# family has NO cross-angle dependency ..., so flattening is exact here") -- there is not
# even a barrier to preserve across the deletion, unlike the Layered/wavevector sibling.
# Genuine two-axis parallelism (outer omega0 loop nested with each `CalcIsofreqAngSweepLayered`'s
# own inner angle loop, both via `forEachElement`/`Threads.@threads`) replaces BOTH MATLAB
# code paths (default-nested-serial-outer AND flattenSweep=true) with the one Julia can do
# safely and directly -- see `CalcWavevecAngSweep2Layered.jl`'s header for the full
# nested-`Threads.@threads`/BLAS-pinning safety argument, unchanged here.
#
# NO `calcMode`/`retentionLevel`/`nVelGreen`/`vMax` field at all (unlike the wavevector
# sibling) -- matches `CalcIsofreqAngSweepLayered` itself having none: no Ritz/
# GammaMatrix phase split, no retention policy, only a Green isoFreq spectrum scan.
#
# `sweeps` is populated ONLY by `calc`, NOT the constructor -- same convention as the
# wavevector sibling (see that file's header).
#
# No `assembleMaps` step at all -- matches the MATLAB source exactly (`calc` just builds
# `obj.sweeps`, no top-level aggregation; `docs/reference/architecture/forwardCalcs.md`'s
# "Precompute sharing differs between the two" and "`retentionLevel='none'` default, and no
# top-level aggregation" sections, referenced in the MATLAB docblock, apply here unchanged). Callers read `sweep.sweeps[iOmega].absGMap` directly, exactly as
# `CalcWavevecAngSweep2Layered`'s own per-wavelength Green maps are read via
# `sweep.sweeps[iWav].greenAbsGMap` rather than a top-level aggregate.
# =========================================================================================

"""
    CalcFreqAngSweep2Layered

2D sweep (propagation angle AND real reference frequency `omega0`) for the layered forward
calculation, isoFreq side.

Port of `guliwerr/forwardCalcs/CalcFreqAngSweep2Layered.m`. **Mutable struct** (MATLAB handle semantics). Built on [`CalcIsofreqAngSweepLayered`](@ref): `calc` constructs one per
`omega0` value and drives it through its own `calc`. No cross-angle dependency at all (no
phase split, no shared axis) -- every point is fully independent.

**No `flattenSweep` field** (see this file's header) -- both sweep axes parallelise directly
and simultaneously whenever `useParallel` resolves to `true` on each.

# Fields
- `structure::Structure` -- the UNROTATED base structure; any `nlay >= 1`.
- `angles::Vector{Float64}` -- propagation angles [deg]; default `0:2:180`.
- `omega0List::Vector{Float64}` -- fixed reference frequencies -- REQUIRED, no default.
- `kParRange::Vector{Float64}` -- wavevector-magnitude axis, shared across every (angle,
  omega0) point -- REQUIRED, no default.
- `useParallel::Symbol` -- `:auto` (default) | `:always` | `:never` -- governs BOTH the outer
  omega0 loop and each omega0's own inner angle loop.
- `sweeps::Vector{CalcIsofreqAngSweepLayered}` -- one per omega0; filled by [`calc`](@ref)
  only (empty before the first `calc` call). No top-level aggregation -- read
  `sweeps[iOmega].absGMap` directly.
"""
mutable struct CalcFreqAngSweep2Layered
    structure::Structure
    angles::Vector{Float64}
    omega0List::Vector{Float64}
    kParRange::Vector{Float64}
    useParallel::Symbol

    sweeps::Vector{CalcIsofreqAngSweepLayered}
end

_calc2freq_isempty(x) = x === nothing || (!(x isa Real) && isempty(x))
_calc2freq_tovec(x) = x isa Real ? [Float64(x)] : Float64.(vec(collect(x)))

"""
    CalcFreqAngSweep2Layered(; structure = nothing, angles = nothing, omega0List = nothing,
                                kParRange = nothing, useParallel = :auto) ->
        CalcFreqAngSweep2Layered

Construct `CalcFreqAngSweep2Layered` for a 2D (angle x omega0) sweep.

Port of `CalcFreqAngSweep2Layered.m`. `structure === nothing` (MATLAB's `nargin == 0`)
returns a well-formed empty object (`angles` still defaults to `0:2:180`).

Otherwise: `omega0List` required (no default), real/finite/positive, else `ArgumentError`
(`"CalcFreqAngSweep2Layered:constructor:missingOmega0List"`); `kParRange` required, real/
finite/positive, else `ArgumentError`
(`"CalcFreqAngSweep2Layered:constructor:missingKParRange"`).

`useParallel` is validated eagerly (mirrors MATLAB's `mustBeMember` property validator).

Keyword-only, not positional -- see `CalcSingleHS`'s constructor for why.
"""
function CalcFreqAngSweep2Layered(; structure::Union{Structure,Nothing} = nothing,
                                     angles = nothing,
                                     omega0List = nothing,
                                     kParRange = nothing,
                                     useParallel::Symbol = :auto)
    _validate_use_parallel(useParallel, "CalcFreqAngSweep2Layered:constructor")

    angs = (angles === nothing || isempty(angles)) ? _ANGSWEEP_DEFAULT_ANGLES :
           Float64.(vec(collect(angles)))

    if structure === nothing
        return CalcFreqAngSweep2Layered(Structure(), angs, Float64[], Float64[], useParallel,
                                         CalcIsofreqAngSweepLayered[])
    end

    _calc2freq_isempty(omega0List) &&
        throw(ArgumentError("CalcFreqAngSweep2Layered:constructor:missingOmega0List -- " *
            "omega0List is required (no default) -- a CalcFreqAngSweep2Layered with no " *
            "frequency reference is meaningless."))
    omL = _calc2freq_tovec(omega0List)
    (all(isfinite, omL) && all(>(0), omL)) ||
        throw(ArgumentError("CalcFreqAngSweep2Layered: omega0List must be a real, finite, " *
            "positive vector"))

    _calc2freq_isempty(kParRange) &&
        throw(ArgumentError("CalcFreqAngSweep2Layered:constructor:missingKParRange -- " *
            "kParRange is required (no default) -- the wavevector-magnitude axis to sweep, " *
            "shared across every (angle, omega0) point."))
    kpr = _calc2freq_tovec(kParRange)
    (all(isfinite, kpr) && all(>(0), kpr)) ||
        throw(ArgumentError("CalcFreqAngSweep2Layered: kParRange must be a real, finite, " *
            "positive vector"))

    return CalcFreqAngSweep2Layered(structure, angs, omL, kpr, useParallel,
                                     CalcIsofreqAngSweepLayered[])
end

"""
    _rebuild(sweep::CalcFreqAngSweep2Layered; <field> = sweep.<field>, ...) ->
        CalcFreqAngSweep2Layered

In-place field update (the Calc types are mutable, MATLAB handle semantics). Internal; not exported.
"""
function _rebuild(x::CalcFreqAngSweep2Layered; kw...)
    for (k, v) in kw
        setproperty!(x, k, v)
    end
    return x
end

"""
    copy(x::CalcFreqAngSweep2Layered)

Shallow copy (MATLAB `copy` of a handle object, one level).
"""
Base.copy(x::CalcFreqAngSweep2Layered) = CalcFreqAngSweep2Layered((getfield(x, k) for k in fieldnames(CalcFreqAngSweep2Layered))...)

# -----------------------------------------------------------------------------------------
# calc -- CalcFreqAngSweep2Layered.m (the DEFAULT, `flattenSweep = false` path -- the
# ONLY path ported, see this file's header)
# -----------------------------------------------------------------------------------------

"""
    calc(sweep::CalcFreqAngSweep2Layered) -> CalcFreqAngSweep2Layered

Run the 2D sweep: one [`CalcIsofreqAngSweepLayered`](@ref) per omega0, BOTH axes
parallelised directly (see this file's header), leaving the per-omega0 spectra distributed
with no top-level aggregation.

Port of `CalcFreqAngSweep2Layered.calc`. Plain nested
loop, no shared precompute -- each `sweeps[iOmega]` is constructed and driven through its
own `calc`.
"""
function calc(sweep::CalcFreqAngSweep2Layered)
    nOmega = length(sweep.omega0List)
    useP, innerP = nestedParallel(sweep.useParallel, nOmega, length(sweep.angles))
    pending = Vector{CalcIsofreqAngSweepLayered}(undef, nOmega)
    for iOmega in 1:nOmega
        pending[iOmega] = CalcIsofreqAngSweepLayered(; structure = sweep.structure,
                                                         angles = sweep.angles,
                                                         omega0 = sweep.omega0List[iOmega],
                                                         kParRange = sweep.kParRange,
                                                         useParallel = innerP)
    end

    # OUTER axis (omega0), genuinely thread-parallel -- see this file's header (and
    # CalcWavevecAngSweep2Layered.jl's, for the full nested-Threads.@threads safety
    # argument, unchanged here). No barrier at all in this family, so nesting is even more
    # straightforward than the wavevector sibling's.
    sweeps = forEachElement(pending, useP, (sw0, _) -> calc(sw0))
    for sw in sweeps                       # the stored setting is the user's, not the run's
        sw.useParallel = sweep.useParallel
    end

    return _rebuild(sweep; sweeps = sweeps)
end
