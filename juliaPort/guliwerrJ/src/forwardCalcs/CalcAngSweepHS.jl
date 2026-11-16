# =========================================================================================
# CalcAngSweepHS.jl -- port of guliwerr/forwardCalcs/CalcAngSweepHS.m.
#
# Angular sweep of the homogeneous half-space forward calculation: one CalcSingleHS per
# angle, one shared Ritz precompute, one common Green axis, one GreenScore reference.
#
# Julia execution of 'full' (same results as MATLAB's phase order, less memory and fewer
# barriers): calcLim of every angle first (milliseconds; the common axis needs the largest
# longLimVel), the reference, then ONE threaded loop over the angles that runs, per angle,
# calcRitz -> greenStep (Green solve, scores, classify, prune). MATLAB runs every angle's
# Ritz solve before any Green solve and so holds all angles' eigenvectors at once; here at
# most one per running task. Each task writes only its own angle object (disjoint slots,
# so no locking is needed); the shared Ritz precompute and reference are read-only.
#
# calcSpotGreenMap (the experimental SpotGreen picture) is
# ported: `spotGreenMap`, `spotGreenVelAxis`, `spotGreenInfo`, kept by calc()/classify() (they
# need only the structure and the angles). Not ported: saveResult/loadResult.
# =========================================================================================

"""
    CalcAngSweepHS

Half-space angular sweep (MATLAB `CalcAngSweepHS`): `structure` (unrotated), `angles`,
`calcMode`, `retentionLevel`, `nVelGreen` (NaN: use `greenVelStep`), `useParallel`
(`:auto` | `:always` | `:never`), `cs::Vector{CalcSingleHS}`, `greenVelAxis` (the common
axis), `longVelMap`, `bulkVelMap` [3 x nAngles], `greenAbsG13EffMap` [nVel x nAngles]
Float32, `sawVelMap`, `sawModeNoMap`, `sawTypeMap`, `significantModesMap` (per angle),
`greenScoreRef`, `refinedSawRef`, the spot-averaged Green map of [`calcSpotGreenMap`](@ref)
(EXPERIMENTAL) `spotGreenMap` [nVel x nAngles] Float32, `spotGreenVelAxis`, `spotGreenInfo`
(`nothing` until then), `activeMethod`, `activeThreshold`, `sawClassMethod`,
`sawThreshold`, `greenVelStep`, `precompStructForRitz`.
"""
mutable struct CalcAngSweepHS
    structure::Structure
    angles::Vector{Float64}
    calcMode::Symbol
    retentionLevel::Symbol
    nVelGreen::Float64
    useParallel::Symbol
    cs::Vector{CalcSingleHS}
    greenVelAxis::Vector{Float64}
    longVelMap::Vector{Float64}
    bulkVelMap::Matrix{Float64}
    greenAbsG13EffMap::Matrix{Float32}
    sawVelMap::Vector{Vector{Float64}}
    sawModeNoMap::Vector{Vector{Int}}
    sawTypeMap::Vector{Vector{Symbol}}
    significantModesMap::Vector{Vector{SignificantMode}}
    greenScoreRef::Union{Nothing,GreenScore.Reference}
    refinedSawRef::Union{Nothing,GreenSpotScore.Reference}
    spotGreenMap::Matrix{Float32}
    spotGreenVelAxis::Vector{Float64}
    spotGreenInfo::Union{Nothing,NamedTuple}
    activeMethod::Union{Nothing,Symbol}
    activeThreshold::Union{Nothing,Symbol,Float64}
    sawClassMethod::Union{Nothing,Symbol}
    sawThreshold::Union{Nothing,Symbol,Float64}
    greenVelStep::Float64
    precompStructForRitz::NamedTuple
end

const _ANGSWEEP_DEFAULT_ANGLES = collect(0.0:2.0:180.0)

function _validate_use_parallel(useParallel, errPrefix)
    Symbol(useParallel) in (:auto, :always, :never) || throw(ArgumentError(
        "$errPrefix:badUseParallel -- useParallel must be 'auto', 'always' or 'never'; got '$useParallel'."))
    return Symbol(useParallel)
end

"""
    CalcAngSweepHS(; structure = nothing, angles = 0:2:180, calcMode = nothing,
                   retentionLevel = nothing, nVelGreen = nothing, useParallel = :auto)

Build the sweep (one CalcSingleHS per angle, sharing one Ritz precompute). Run it with
[`calc`](@ref).
"""
function CalcAngSweepHS(; structure::Union{Structure,Nothing} = nothing,
                        angles = _ANGSWEEP_DEFAULT_ANGLES, calcMode = nothing,
                        retentionLevel = nothing, nVelGreen = nothing, useParallel = :auto)
    up = _validate_use_parallel(useParallel, "CalcAngSweepHS:constructor")
    sw = CalcAngSweepHS(Structure(), copy(_ANGSWEEP_DEFAULT_ANGLES), :full, :full, NaN, up,
                        CalcSingleHS[], Float64[], Float64[], zeros(3, 0), zeros(Float32, 0, 0),
                        Vector{Float64}[], Vector{Int}[], Vector{Symbol}[], Vector{SignificantMode}[],
                        nothing, nothing, zeros(Float32, 0, 0), Float64[], nothing,
                        nothing, nothing, nothing, nothing,
                        GreenScore.DEFAULT_STEP, NamedTuple())
    structure === nothing && return sw
    all(isfinite, angles) || throw(ArgumentError("CalcAngSweepHS: angles must be finite"))
    validateConstructorFlags(calcMode, retentionLevel, "CalcAngSweepHS:constructor", (:full, :saw_only, :none))
    sw.structure = structure
    sw.angles = collect(Float64, vec(angles))
    calcMode === nothing || (sw.calcMode = calcMode)
    retentionLevel === nothing || (sw.retentionLevel = retentionLevel)
    if nVelGreen !== nothing
        (nVelGreen isa Real && nVelGreen >= 2) ||
            throw(ArgumentError("CalcAngSweepHS: nVelGreen must be a real scalar >= 2; got $nVelGreen"))
        sw.nVelGreen = Float64(nVelGreen)
    end
    sw.precompStructForRitz = ritzSAWPrecompute(structure)
    sw.cs = [CalcSingleHS(; structure = structure, angle = a, precompStruct = sw.precompStructForRitz,
                          calcMode = sw.calcMode, retentionLevel = sw.retentionLevel)
             for a in sw.angles]
    return sw
end

"""
    ritzPrecompute(sweep) -> NamedTuple
"""
function ritzPrecompute(sw::CalcAngSweepHS)
    sw.precompStructForRitz = ritzSAWPrecompute(sw.structure, sw.precompStructForRitz)
    return sw.precompStructForRitz
end

"""
    calc(sweep::CalcAngSweepHS, mode = nothing, retentionLevel = nothing) -> sweep

Run every angle. `:ritzOnly`: per angle `calcAll(:ritzOnly)`. `:full`: per angle Ritz +
Green + scores, classified 'green' (then the chosen method, if another). Every angle carries
the run's method and standard.
"""
function calc(sw::CalcAngSweepHS, mode = nothing, retentionLevel = nothing)
    mode === nothing && (mode = sw.calcMode)
    retentionLevel === nothing && (retentionLevel = sw.retentionLevel)
    validateCalcMode(mode, "CalcAngSweepHS:calc")
    validateRetentionLevel(retentionLevel, "CalcAngSweepHS:calc", (:full, :saw_only, :none))
    method = resolveSawClassMethod(sw.sawClassMethod, sw.sawThreshold, mode, false, "CalcAngSweepHS:calc")
    sw.calcMode = mode
    sw.retentionLevel = retentionLevel
    sw.refinedSawRef = nothing
    for c in sw.cs
        c.sawClassMethod = method
        c.sawThreshold = sw.sawThreshold
    end
    precomp = ritzPrecompute(sw)
    if mode === :ritzOnly
        useP = resolveParallel(sw.useParallel, length(sw.cs))
        thr = sw.sawThreshold
        sw.cs = forEachElement(sw.cs, useP, (c, _) -> ritzOnlyStep(c, retentionLevel, thr, precomp))
    else
        calcFullPhases(sw, precomp, method)
        if method !== :green
            sw.activeMethod = :green
            return classify(sw, method, sw.sawThreshold)
        end
    end
    return assembleMaps(sw)
end

"""
    calcFullPhases(sweep, precomp, method) -> sweep

The 'full' phases (see the file header for the Julia execution order).
"""
function calcFullPhases(sw::CalcAngSweepHS, precomp, method)
    useP = resolveParallel(sw.useParallel, length(sw.cs))
    # Christoffel limiting velocities: independent per angle, so on the threads too (a
    # serial loop here was ~12 % of a serial Fe3Al(316) 'full' sweep -- an Amdahl cap on
    # the parallel one). Needed before the loop below for the common Green axis.
    sw.cs = forEachElement(sw.cs, useP, (c, _) -> calcLim(c))
    vBound = SawSpectrum.AXIS_HI_X_VL * maximum(c.longLimVel for c in sw.cs)
    step = sw.greenVelStep
    isfinite(sw.nVelGreen) &&
        (step = (vBound + SawSpectrum.EDGE_TAPER - SawSpectrum.VEL_MIN) / (sw.nVelGreen - 1))
    step = GreenScore.resolveStep(step)
    sw.greenVelAxis = GreenScore.gridAxis(vBound, step)
    for c in sw.cs
        c.greenVelStep = step
        c.nVelGreen = NaN
    end
    sw.greenScoreRef = GreenScore.reference(sw.structure; useParallel = sw.useParallel)
    ref = sw.greenScoreRef; axis = sw.greenVelAxis; rl = sw.retentionLevel
    thr = method === :green ? sw.sawThreshold : nothing
    sw.cs = forEachElement(sw.cs, useP, (c, _) -> begin
        calcRitz(c, precomp)
        greenStep(c, axis, rl, ref, thr, precomp)
    end)
    return sw
end

"""
    classify(sweep::CalcAngSweepHS, method = nothing, threshold = nothing) -> sweep

Switch every angle's classification (`:green` | `:refined` | `:ritz`; 'refined' solves the
spot directions of all angles, plus the reference's on a structure's first use, in one
parallel batch).
"""
function classify(sw::CalcAngSweepHS, method = nothing, threshold = nothing)
    threshold = normaliseThreshold(threshold)
    method === nothing && (method = sw.activeMethod)
    method === nothing && (method = sw.sawClassMethod)
    method === nothing && (method = resolveSawClassMethod(nothing, threshold, sw.calcMode, false, "CalcAngSweepHS:classify"))
    method = validateSawClassMethod(method, false, "CalcAngSweepHS:classify")
    threshold === nothing && (threshold = defaultThreshold(method, false))
    validateSawThreshold(threshold, "CalcAngSweepHS:classify")
    refuseStale(sw, "CalcAngSweepHS:classify")
    (isempty(sw.cs) || sw.cs[1].classification === nothing) &&
        throw(ArgumentError("CalcAngSweepHS:classify:notCalculated -- calc() must run first."))
    method !== :ritz && sw.calcMode !== :full && throw(ArgumentError("CalcAngSweepHS:classify:needsFull -- " *
        "Method '$method' needs calcMode=:full; this sweep is '$(sw.calcMode)'."))
    if method === :refined && any(c -> c.refinedSaw === nothing, sw.cs)
        D = GreenSpotScore.directionSet(sw.angles)
        if memoGet(GreenSpotScore._REF_MEMO, structureKey(sw.structure)) === nothing
            D = vcat(D, GreenSpotScore.directionSet(GreenSpotScore.REF_ANGLES))
        end
        dc = GreenSpotScore.solveDirections(sw.structure, D; useParallel = sw.useParallel)
        ref = GreenSpotScore.reference(sw.structure; dc = dc, useParallel = sw.useParallel)
        sw.refinedSawRef = ref
        for c in sw.cs
            calcRefined(c, dc, ref)
        end
    end
    useP = method !== :refined && resolveParallel(sw.useParallel, length(sw.cs))
    sw.cs = forEachElement(sw.cs, useP, (c, _) -> classify(c, method, threshold))
    return assembleMaps(sw)
end

reclassify(sw::CalcAngSweepHS, threshold) = classify(sw, sw.activeMethod, threshold)

"""
    refinedSawAgreement(sweep) -> Vector{NamedTuple}

Every angle's `refinedSawAgreement` rows, concatenated.
"""
function refinedSawAgreement(sw::CalcAngSweepHS)
    (isempty(sw.cs) || any(c -> c.refinedSaw === nothing, sw.cs)) && throw(ArgumentError(
        "CalcAngSweepHS:refinedSawAgreement:notRefined -- classify(:refined) must run first."))
    return reduce(vcat, [refinedSawAgreement(c) for c in sw.cs])
end

"""
    greenComponentMap(sweep, component) -> Matrix{Float32}

Any Green component over the sweep on greenVelAxis (`"G13eff"` stored, others solved now).
"""
greenComponentMap(sw::CalcAngSweepHS, component) = greenComponentMapOf(sw, component)

"""
    assembleMaps(sweep) -> sweep

The per-angle SAW maps, the active method, and (full mode) longVelMap, bulkVelMap and the
|G13_eff| map.
"""
function assembleMaps(sw::CalcAngSweepHS)
    sw.sawVelMap, sw.sawModeNoMap, sw.sawTypeMap, sw.significantModesMap =
        sawMapsFrom(sw.cs, (:Rayleigh, :pseudoRayleigh))
    if !isempty(sw.cs) && sw.cs[1].activeMethod !== nothing
        sw.activeMethod = sw.cs[1].activeMethod
        sw.activeThreshold = sw.cs[1].activeThreshold
    end
    if sw.calcMode === :full
        sw.longVelMap = [c.longLimVel for c in sw.cs]
        sw.bulkVelMap = reduce(hcat, [c.bulkVels for c in sw.cs])
        sw.greenAbsG13EffMap = greenG13EffMapFrom(sw.cs)
    end
    return sw
end

"""
    calcSpotGreenMap(sweep, kernel::ProbeKernel; DirStep, RelKStep, Eta, VelStep,
                     useParallel = sweep.useParallel) -> sweep

EXPERIMENTAL (as in MATLAB; tried on Fe3Al(316) only, not validated against measurements).
The |G| map as a finite-pattern measurement sees it: every displayed angle is the coherent
sum of the complex Green's function over the directions and |k| of a [`ProbeKernel`](@ref)
(`SpotGreen.compute`): a steep branch broadens and drops, a flat one keeps its height. Stores `spotGreenMap` (`[nVel x nAngles]`
`Float32`, |X|, same size as `greenAbsG13EffMap` on the 4 m/s axis), `spotGreenVelAxis` [km/s]
and `spotGreenInfo` (kernel, steps, eta, `vEnd`, `nDirections`, `nNodes`, per-angle
`dirConvergence` / `relKConvergence`). Options as [`calcSpotGreen`](@ref): `DirStep` [deg],
`RelKStep`, `Eta`, `VelStep` [km/s]; `useParallel` (`:auto` | `:always` | `:never`, default
the sweep's) is the thread pool over the direction solves (the result does not depend on it).
One Green solve per distinct direction, shared between angles (with integer-degree angles and
the default direction steps, a 0:1:180 map needs 180 / `DirStep` solves). Any calcMode;
half-space only (a layered structure throws an `ArgumentError`
`CalcAngSweepHS:calcSpotGreenMap:layered`). Plot with `plotGreenMap(sweep; component = "spot")`.
A warning (`SpotGreen:notConverged`) reports a quadrature step too coarse for the branch
slopes (the "comb").
"""
function calcSpotGreenMap(sw::CalcAngSweepHS, kernel::ProbeKernel; useParallel = nothing, kw...)
    nlay(sw.structure) == 1 || throw(ArgumentError("CalcAngSweepHS:calcSpotGreenMap:layered -- " *
        "calcSpotGreenMap is for a half-space (nLay = 1) only."))
    opts = SpotGreen.options(kernel; kw...)
    up = useParallel === nothing ? sw.useParallel :
         _validate_use_parallel(useParallel, "CalcAngSweepHS:calcSpotGreenMap")
    sw.spotGreenMap, sw.spotGreenVelAxis, sw.spotGreenInfo =
        SpotGreen.compute(sw.structure, sw.angles, kernel, opts; useParallel = up)
    return sw
end
