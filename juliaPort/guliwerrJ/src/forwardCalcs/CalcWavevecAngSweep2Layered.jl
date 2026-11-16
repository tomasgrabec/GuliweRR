# =========================================================================================
# CalcWavevecAngSweep2Layered.jl -- 2D sweep (propagation angle AND real wavelength) for the
# layered forward calculation. Port of guliwerr/forwardCalcs/CalcWavevecAngSweep2Layered.m.
# Built on CalcAngSweepLayered: one per wavelength (via ApplyWavelength), each driven through
# its own `calc`.
#
# NO per-file module (see Calc.jl's header for why).
#
# MUTABLE: a `mutable struct`, like CalcSingleHS/Layered and their sweeps
# (MATLAB handle semantics): every step mutates the object and returns it, so both
# `calc(x)` and `x = calc(x)` work; `copy(x)` is a shallow copy.
#
# Settings (as in MATLAB):
# - `vMax` is optional (NaN = each angle's band end, 1.05 x the fastest longitudinal limit);
#   there is no `missingVMax` error. MATLAB's 7th constructor argument is unused (it must be
#   empty) and has no Julia keyword: the Green axis is each angle's own.
# - `sawClassMethod` / `sawThreshold` are forwarded to every wavelength's sweep, which
#   resolves the method per run (unset: 'green' on :full, 'ritz' on :ritzOnly);
#   `classify(sweep, method, threshold)` / `reclassify(sweep, threshold)` switch the active
#   classification of every wavelength without recomputing.
# - No reference is computed per wavelength: every grid point is scored against its own band
#   / floor.
#
# DELIBERATE DIVERGENCE: no `flattenSweep` / `calcFlat`.
# MATLAB's flat path exists only because MATLAB forbids a parfor inside a parfor; its default
# path keeps the outer wavelength loop serial. Julia's scheduler is composable: the OUTER
# wavelength loop here and the INNER angle loop of every CalcAngSweepLayered.calc are both
# `forEachElement` (Threads.@threads :dynamic) and nest without restriction. BLAS is pinned
# to one thread by the outermost `with_blas_threads(1)`, and every nested call reads and
# restores that same value. Each wavelength keeps its own per-wavelength Green-axis barrier
# (CalcAngSweepLayered.commonGreenAxis over THAT wavelength's angles), exactly MATLAB's.
#
# `rawByDepthCell` (GammaMatrix.RawCellFor) depends on legOrder only, so it is built once in
# the constructor and shared by every wavelength and every angle.
#
# WARNINGS (carried from the MATLAB source): legOrder must stay constant across the
# wavelength sweep (rawByDepthCell is shared); the substrate (layer 1) is not swept -- choose
# its relDepth deep relative to the smallest wavelength.
#
# retentionLevel DEFAULTS TO :none (eigenvector RAM multiplies across the wavelength axis);
# :full is needed by DispersionBranches.track / trackSeededBranch.
# =========================================================================================

using ..CoreTypes: ApplyWavelength
using ..Solvers: RawCellFor

"""
    CalcWavevecAngSweep2Layered

2D sweep (propagation angle AND real wavelength), layered. Fields (MATLAB names):
`structure` (UNROTATED reference; every film layer's `realThickness` set), `angles`,
`wavelengths`, `vMax` (NaN = the band end), `calcMode` (`:full` | `:ritzOnly`),
`retentionLevel` (`:full` | `:saw_only` | `:none`, default `:none`), `useParallel`,
`sawClassMethod`, `sawThreshold` (forwarded to every wavelength), `rawByDepthCell`,
`sweeps::Vector{CalcAngSweepLayered}` (filled by [`calc`](@ref)), and the `[nWav x nAngles]`
maps `sawVelMap2D`, `sawModeNoMap2D`, `sawTypeMap2D`, `significantModesMap2D`.
"""
mutable struct CalcWavevecAngSweep2Layered
    structure::Structure
    angles::Vector{Float64}
    wavelengths::Vector{Float64}
    vMax::Float64
    calcMode::Symbol
    retentionLevel::Symbol
    useParallel::Symbol
    sawClassMethod::Union{Nothing,Symbol}
    sawThreshold::Union{Nothing,Symbol,Float64}
    rawByDepthCell::Vector{Array{ComplexF64,4}}
    sweeps::Vector{CalcAngSweepLayered}
    sawVelMap2D::Matrix{Vector{Float64}}
    sawModeNoMap2D::Matrix{Vector{Int}}
    sawTypeMap2D::Matrix{Vector{Symbol}}
    significantModesMap2D::Matrix{Vector{SignificantMode}}
end

"""
    CalcWavevecAngSweep2Layered(; structure = nothing, angles = 0:2:180, wavelengths = nothing,
                                  vMax = nothing, calcMode = nothing, retentionLevel = nothing,
                                  useParallel = :auto, sawClassMethod = nothing,
                                  sawThreshold = nothing)

Construct the 2D sweep (MATLAB constructor, `CalcWavevecAngSweep2Layered.m`). No
structure: an empty object (`angles` 0:2:180, `retentionLevel` `:none`). Otherwise, in
MATLAB's order: nLay >= 2 (`:needsFilmAndSubstrate`), calcMode / retentionLevel validated,
`wavelengths` required and positive (`:missingWavelengths`), every film layer's
`realThickness` set (`:missingRealThickness`). `rawByDepthCell` is built here, once.
`sawClassMethod` / `sawThreshold` are MATLAB's public properties (settable on the object too).
"""
function CalcWavevecAngSweep2Layered(; structure::Union{Structure,Nothing} = nothing,
                                     angles = nothing, wavelengths = nothing, vMax = nothing,
                                     calcMode = nothing, retentionLevel = nothing,
                                     useParallel = :auto, sawClassMethod = nothing,
                                     sawThreshold = nothing)
    up = _validate_use_parallel(useParallel, "CalcWavevecAngSweep2Layered:constructor")
    angs = (angles === nothing || isempty(angles)) ? copy(_ANGSWEEP_DEFAULT_ANGLES) :
           Float64.(vec(collect(angles)))
    sw = CalcWavevecAngSweep2Layered(Structure(), angs, Float64[], NaN, :full, :none, up,
        sawClassMethod === nothing ? nothing : Symbol(sawClassMethod),
        normaliseThreshold(sawThreshold),
        Array{ComplexF64,4}[], CalcAngSweepLayered[],
        Matrix{Vector{Float64}}(undef, 0, 0), Matrix{Vector{Int}}(undef, 0, 0),
        Matrix{Vector{Symbol}}(undef, 0, 0), Matrix{Vector{SignificantMode}}(undef, 0, 0))
    structure === nothing && return sw

    nlay(structure) >= 2 ||
        throw(ArgumentError("CalcWavevecAngSweep2Layered:constructor:needsFilmAndSubstrate -- " *
            "CalcWavevecAngSweep2Layered needs at least one film plus a substrate (nLay >= 2); " *
            "use CalcAngSweepHS for a homogeneous half-space."))
    all(isfinite, angs) || throw(ArgumentError("CalcWavevecAngSweep2Layered: angles must be finite"))
    validateConstructorFlags(calcMode, retentionLevel, "CalcWavevecAngSweep2Layered:constructor",
                             (:full, :saw_only, :none))
    (wavelengths === nothing || isempty(wavelengths)) &&
        throw(ArgumentError("CalcWavevecAngSweep2Layered:constructor:missingWavelengths -- " *
            "wavelengths is required (no default) -- pass a vector of positive real wavelengths."))
    wls = Float64.(vec(collect(wavelengths)))
    (all(isfinite, wls) && all(>(0), wls)) ||
        throw(ArgumentError("CalcWavevecAngSweep2Layered: wavelengths must be a real, finite, " *
            "positive vector"))
    for ilay in 2:nlay(structure)
        isnan(structure.calcParams[ilay].realThickness) &&
            throw(ArgumentError("CalcWavevecAngSweep2Layered:constructor:missingRealThickness -- " *
                "Layer $ilay has no realThickness set -- cannot sweep wavelength."))
    end
    sw.structure = structure
    sw.wavelengths = wls
    vMax === nothing || (sw.vMax = Float64(vMax))
    calcMode === nothing || (sw.calcMode = calcMode)
    retentionLevel === nothing || (sw.retentionLevel = retentionLevel)
    sw.rawByDepthCell = RawCellFor(structure)
    return sw
end

function _rebuild(x::CalcWavevecAngSweep2Layered; kw...)
    for (k, v) in kw
        setproperty!(x, k, v)
    end
    return x
end

"""
    copy(x::CalcWavevecAngSweep2Layered)

Shallow copy (MATLAB `copy` of a handle object, one level).
"""
Base.copy(x::CalcWavevecAngSweep2Layered) =
    CalcWavevecAngSweep2Layered((getfield(x, k) for k in fieldnames(CalcWavevecAngSweep2Layered))...)

"""
    calc(sweep::CalcWavevecAngSweep2Layered, mode = nothing, retentionLevel = nothing) -> sweep

Run the 2D sweep (`CalcWavevecAngSweep2Layered.calc`): validate, resolve the run's
classification method (errors as MATLAB: numeric threshold without a method, a non-'ritz'
method on :ritzOnly), then one [`CalcAngSweepLayered`](@ref) per wavelength (sharing
`rawByDepthCell`, carrying `vMax`, `sawClassMethod`, `sawThreshold`), each run through its
own `calc`, wavelengths in parallel (see this file's header). Then the 2D maps.
"""
function calc(sweep::CalcWavevecAngSweep2Layered, mode = nothing, retentionLevel = nothing)
    m = mode === nothing ? sweep.calcMode : mode
    rl = retentionLevel === nothing ? sweep.retentionLevel : retentionLevel
    validateCalcMode(m, "CalcWavevecAngSweep2Layered:calc")
    validateRetentionLevel(rl, "CalcWavevecAngSweep2Layered:calc", (:full, :saw_only, :none))
    resolveSawClassMethod(sweep.sawClassMethod, sweep.sawThreshold, m, true,
                          "CalcWavevecAngSweep2Layered:calc")   # fail early, before any solve
    sweep.calcMode = m
    sweep.retentionLevel = rl
    isempty(sweep.rawByDepthCell) && (sweep.rawByDepthCell = RawCellFor(sweep.structure))

    nWav = length(sweep.wavelengths)
    useP, innerP = nestedParallel(sweep.useParallel, nWav, length(sweep.angles))
    pending = Vector{CalcAngSweepLayered}(undef, nWav)
    for iWav in 1:nWav
        structAtWav = ApplyWavelength(sweep.structure, sweep.wavelengths[iWav])
        sw = CalcAngSweepLayered(; structure = structAtWav, angles = sweep.angles,
                                 vMax = isnan(sweep.vMax) ? nothing : sweep.vMax,
                                 calcMode = m, retentionLevel = rl,
                                 rawByDepthCell = sweep.rawByDepthCell,
                                 useParallel = innerP)
        sw.sawClassMethod = sweep.sawClassMethod
        sw.sawThreshold = sweep.sawThreshold
        pending[iWav] = sw
    end
    sweep.sweeps = forEachElement(pending, useP, (sw0, _) -> calc(sw0, m, rl))
    for sw in sweep.sweeps                 # the stored setting is the user's, not the run's
        sw.useParallel = sweep.useParallel
    end
    return assembleMaps(sweep)
end

"""
    classify(sweep::CalcWavevecAngSweep2Layered, method = nothing, threshold = nothing) -> sweep

Set the active classification of every wavelength's sweep (`CalcAngSweepLayered` `classify`:
'green' | 'ritz', threshold 'low' (default) | 'high' | numeric), then reassemble the 2D maps.
Throws `"CalcWavevecAngSweep2Layered:classify:notCalculated"` before `calc`.
"""
function classify(sweep::CalcWavevecAngSweep2Layered, method = nothing, threshold = nothing)
    isempty(sweep.sweeps) &&
        throw(ArgumentError("CalcWavevecAngSweep2Layered:classify:notCalculated -- calc() must run first."))
    for sw in sweep.sweeps
        classify(sw, method, threshold)
    end
    return assembleMaps(sweep)
end

"""
    reclassify(sweep::CalcWavevecAngSweep2Layered, threshold) -> sweep

Re-threshold the active classification of every wavelength (no recomputation).
"""
reclassify(sweep::CalcWavevecAngSweep2Layered, threshold) = classify(sweep, nothing, threshold)

isStale(sweep::CalcWavevecAngSweep2Layered) = any(isStale, sweep.sweeps)

"""
    assembleMaps(sweep::CalcWavevecAngSweep2Layered) -> sweep

`[nWav x nAngles]` re-indexing of every wavelength's 1D maps -- no new classification logic.
"""
function assembleMaps(sweep::CalcWavevecAngSweep2Layered)
    nWav = length(sweep.sweeps)
    nAngles = length(sweep.angles)
    sweep.sawVelMap2D = Matrix{Vector{Float64}}(undef, nWav, nAngles)
    sweep.sawModeNoMap2D = Matrix{Vector{Int}}(undef, nWav, nAngles)
    sweep.sawTypeMap2D = Matrix{Vector{Symbol}}(undef, nWav, nAngles)
    sweep.significantModesMap2D = Matrix{Vector{SignificantMode}}(undef, nWav, nAngles)
    for iWav in 1:nWav
        sweep.sawVelMap2D[iWav, :] = sweep.sweeps[iWav].sawVelMap
        sweep.sawModeNoMap2D[iWav, :] = sweep.sweeps[iWav].sawModeNoMap
        sweep.sawTypeMap2D[iWav, :] = sweep.sweeps[iWav].sawTypeMap
        sweep.significantModesMap2D[iWav, :] = sweep.sweeps[iWav].significantModesMap
    end
    return sweep
end
