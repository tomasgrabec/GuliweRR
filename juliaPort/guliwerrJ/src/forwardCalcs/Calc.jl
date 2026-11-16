# =========================================================================================
# Calc.jl -- the `Calc.*` launchers
#
# Port of guliwerr/forwardCalcs/@Calc/ (Calc.m itself + the one-line launchers).
# `classdef Calc` is a MATLAB STATIC NAMESPACE (no properties, `methods (Static)` only, never
# instantiated), so it has no Julia type: the launchers are plain top-level functions in
# `module ForwardCalcs`, exactly like `CalcCommon.jl`'s own functions. (A per-file module named
# after the file's main type would collide with that type name, so the port uses none.)
#
# THINNEST POSSIBLE LAUNCHER, not a positional MATLAB-order re-implementation: every MATLAB
# `@Calc/*.m` file is a 6-8 line `obj = CalcXxx(varargin{:})` passthrough (see
# `guliwerr/forwardCalcs/@Calc/angSweepHS.m` etc.) -- literally a renamed call, no logic of its
# own, no argument-order translation to preserve. Every calculator struct in this port is
# KEYWORD-ONLY (see `CalcSingleHS.jl`'s constructor docstring for why: MATLAB's own
# positional-with-`nargin`-guards style breaks `Pkg.precompile()` for a positional outer
# constructor). Reconstructing MATLAB's positional call signature here
# (`Calc.angSweepHS(structure, angles, calcMode, retentionLevel, nVelGreen, useParallel)`)
# would mean re-deriving and hand-maintaining each target struct's own argument order a SECOND
# time in a place with no compiler check that the two stay in sync -- exactly the kind of
# drift-prone duplication the "thin launcher" MATLAB design was avoiding in the first place.
# Each launcher below is therefore a bare alias, forwarding ALL keyword arguments verbatim to
# the struct's own constructor -- e.g. `angSweepHS(; kwargs...) = CalcAngSweepHS(; kwargs...)`
# -- which also means `angSweepHS()` (zero args) works exactly like `Calc.angSweepHS()` does in
# MATLAB, with no special-casing needed here.
#
#   launchers:  singleHS, angSweepHS, singleLayered, angSweepLayered, isofreqSingleLayered,
#               isofreqAngSweepLayered, singlePlate, angSweepPlate, modeClassification,
#               freqToOmega0, wavevecAngSweep2Layered, wavevecAngSweep2Plate,
#               freqAngSweep2Layered
# The three 2D-sweep launchers (the last three structs) are the entry points MATLAB's
# `StructureControlGUI.onRunClicked` calls for its "2-axis sweep" calcType; see
# `StructureControlJob.jl`.
#
# `modeClassification(; kw...)` -> `ModeClassification(; kw...)`: the settings struct;
# `calc(mc, ritz, ...)` is MATLAB's `.calc(...)` (multiple dispatch on the first argument).
#
# Not ported: MATLAB's `Calc.saveResult` / `Calc.loadResult` (MATLAB .mat files with
# provenance metadata) -- reading and writing MATLAB-saved results is not available in the
# Julia port.
#
# `freqToOmega0` is NOT a struct-constructor forward (there is no `CalcFreqToOmega0` type) --
# it is ported directly, keeping MATLAB's two required POSITIONAL arguments (`freqHz`,
# `lambdaRefNm`) and its `validateattributes` input checks (the numeric formula and the two
# rejection paths).
# =========================================================================================

singleHS(; kwargs...) = CalcSingleHS(; kwargs...)
angSweepHS(; kwargs...) = CalcAngSweepHS(; kwargs...)

singleLayered(; kwargs...) = CalcSingleLayered(; kwargs...)
angSweepLayered(; kwargs...) = CalcAngSweepLayered(; kwargs...)

isofreqSingleLayered(; kwargs...) = CalcIsofreqSingleLayered(; kwargs...)
isofreqAngSweepLayered(; kwargs...) = CalcIsofreqAngSweepLayered(; kwargs...)

singlePlate(; kwargs...) = CalcSinglePlate(; kwargs...)
angSweepPlate(; kwargs...) = CalcAngSweepPlate(; kwargs...)

wavevecAngSweep2Layered(; kwargs...) = CalcWavevecAngSweep2Layered(; kwargs...)

wavevecAngSweep2Plate(; kwargs...) = CalcWavevecAngSweep2Plate(; kwargs...)

freqAngSweep2Layered(; kwargs...) = CalcFreqAngSweep2Layered(; kwargs...)

# Help text: `?angSweepHS` shows the launcher line followed by the type's own docstring, so
# the launchers (the names users call) carry the full documentation of what they build.
for (f, T) in ((:singleHS, :CalcSingleHS), (:angSweepHS, :CalcAngSweepHS),
               (:singleLayered, :CalcSingleLayered), (:angSweepLayered, :CalcAngSweepLayered),
               (:isofreqSingleLayered, :CalcIsofreqSingleLayered),
               (:isofreqAngSweepLayered, :CalcIsofreqAngSweepLayered),
               (:singlePlate, :CalcSinglePlate), (:angSweepPlate, :CalcAngSweepPlate),
               (:wavevecAngSweep2Layered, :CalcWavevecAngSweep2Layered),
               (:wavevecAngSweep2Plate, :CalcWavevecAngSweep2Plate),
               (:freqAngSweep2Layered, :CalcFreqAngSweep2Layered))
    head = "    $f(; kwargs...) -> $T\n\nLauncher: builds a `$T` with the same keywords; run it with `calc`.\n\n"
    md = Base.Docs.meta(@__MODULE__)[Base.Docs.Binding(@__MODULE__, T)]
    body = join((join(string.(d.text)) for d in values(md.docs)), "\n\n")
    @eval @doc $(head * body) $f
end

"""
    modeClassification(; kwargs...) -> ModeClassification

Thin `Calc` namespace launcher for [`ModeClassification`](@ref) (settings as keywords);
run it with [`calc`](@ref) (`calc(Calc.modeClassification(), ritz, structure, ...)`).

Port of `guliwerr/forwardCalcs/@Calc/modeClassification.m`.
"""
modeClassification(; kwargs...) = ModeClassification(; kwargs...)

"""
    freqToOmega0(freqHz, lambdaRefNm::Real) -> omega0

Convert a real frequency to `GreenSolver`'s normalised isoFreq `omega0`.

Port of `guliwerr/forwardCalcs/@Calc/freqToOmega0.m`. `freqHz` -- real frequency [Hz], scalar
or array (broadcasts elementwise, matching MATLAB's `.*`); `lambdaRefNm` -- reference
wavelength [nm], a SCALAR -- MUST be the same wavelength passed to `ApplyWavelength` when
building `relDepth` for the isoFreq structure (that call is what fixes what `|k_par| = 1`
means in real units; a different reference here and `omega0`/`relDepth` disagree about the
length scale).

`omega0` is the phase velocity that would give `|k_par| = 1` at the chosen `lambdaRefNm` --
the standard `v = f*lambda` relation at that one reference wavelength, converted from
`[Hz]*[nm]` to `[km/s]` by `1e-12` (`freqToOmega0.m`, reproduced here in the SAME
operation order for bit-exact agreement).

`lambdaRefNm` non-scalar is rejected by Julia's own dispatch (`MethodError`; no
`validateattributes`-equivalent is needed for a case the type system already prevents
structurally) rather than a hand-written runtime check.

Throws `ArgumentError` if `freqHz` is not real/finite/positive (elementwise), or if
`lambdaRefNm` is not real/finite/positive.
"""
function freqToOmega0(freqHz, lambdaRefNm::Real)
    freqHzVec = freqHz isa Real ? (freqHz,) : freqHz
    all(x -> isfinite(x) && x > 0, freqHzVec) ||
        throw(ArgumentError("freqToOmega0: freqHz must be real, finite, and positive"))
    (isfinite(lambdaRefNm) && lambdaRefNm > 0) ||
        throw(ArgumentError("freqToOmega0: lambdaRefNm must be a real, finite, positive scalar"))
    return freqHz .* lambdaRefNm .* 1e-12
end
