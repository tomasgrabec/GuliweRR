# =========================================================================================
# CalcAngSweepLayered.jl -- port of guliwerr/forwardCalcs/CalcAngSweepLayered.m.
#
# Angular sweep of a layered stack: one CalcSingleLayered per angle, one shared Ritz
# precompute, one common Green axis (commonGreenAxis: the 4 m/s lattice spanning every
# angle's own axis). 'full' runs calcLim of all angles first, then one threaded loop of
# calcRitz -> greenStep per angle (see CalcAngSweepHS.jl's header).
# =========================================================================================

"""
    CalcAngSweepLayered

Layered angular sweep (MATLAB `CalcAngSweepLayered`): `structure` (unrotated, nLay >= 2),
`angles`, `calcMode`, `retentionLevel`, `useParallel`, `vMax`, `cs::Vector{CalcSingleLayered}`,
`subLongVelMap`, `subBulkVelMap`, `greenVelAxis`, `greenAbsG13EffMap`, `sawVelMap`,
`sawModeNoMap`, `sawTypeMap`, `significantModesMap`, `activeMethod`, `activeThreshold`,
`sawClassMethod`, `sawThreshold`, `precompStructForRitz`.
"""
mutable struct CalcAngSweepLayered
    structure::Structure
    angles::Vector{Float64}
    calcMode::Symbol
    retentionLevel::Symbol
    useParallel::Symbol
    vMax::Float64
    cs::Vector{CalcSingleLayered}
    subLongVelMap::Vector{Float64}
    subBulkVelMap::Matrix{Float64}
    greenVelAxis::Vector{Float64}
    greenAbsG13EffMap::Matrix{Float32}
    sawVelMap::Vector{Vector{Float64}}
    sawModeNoMap::Vector{Vector{Int}}
    sawTypeMap::Vector{Vector{Symbol}}
    significantModesMap::Vector{Vector{SignificantMode}}
    activeMethod::Union{Nothing,Symbol}
    activeThreshold::Union{Nothing,Symbol,Float64}
    sawClassMethod::Union{Nothing,Symbol}
    sawThreshold::Union{Nothing,Symbol,Float64}
    precompStructForRitz::NamedTuple
end

"""
    CalcAngSweepLayered(; structure = nothing, angles = 0:2:180, vMax = nothing,
                        calcMode = nothing, retentionLevel = nothing, rawByDepthCell = nothing,
                        useParallel = :auto)

Build the sweep. `rawByDepthCell`: a pre-built GammaMatrix.rawByDepth cell (a 2D sweep
shares it across wavelengths). MATLAB's unused 6th argument (formerly nVelGreen) has no
Julia counterpart.
"""
function CalcAngSweepLayered(; structure::Union{Structure,Nothing} = nothing,
                             angles = _ANGSWEEP_DEFAULT_ANGLES, vMax = nothing, calcMode = nothing,
                             retentionLevel = nothing, rawByDepthCell = nothing, useParallel = :auto)
    up = _validate_use_parallel(useParallel, "CalcAngSweepLayered:constructor")
    sw = CalcAngSweepLayered(Structure(), copy(_ANGSWEEP_DEFAULT_ANGLES), :full, :full, up, NaN,
                             CalcSingleLayered[], Float64[], zeros(3, 0), Float64[],
                             zeros(Float32, 0, 0), Vector{Float64}[], Vector{Int}[],
                             Vector{Symbol}[], Vector{SignificantMode}[], nothing, nothing,
                             nothing, nothing, NamedTuple())
    structure === nothing && return sw
    nlay(structure) < 2 && throw(ArgumentError("CalcAngSweepLayered:constructor:needsFilmAndSubstrate -- " *
        "CalcAngSweepLayered needs at least one film plus a substrate (nLay >= 2); use " *
        "CalcAngSweepHS for a homogeneous half-space."))
    all(isfinite, angles) || throw(ArgumentError("CalcAngSweepLayered: angles must be finite"))
    validateConstructorFlags(calcMode, retentionLevel, "CalcAngSweepLayered:constructor", (:full, :saw_only, :none))
    sw.structure = structure
    sw.angles = collect(Float64, vec(angles))
    vMax === nothing || (sw.vMax = Float64(vMax))
    calcMode === nothing || (sw.calcMode = calcMode)
    retentionLevel === nothing || (sw.retentionLevel = retentionLevel)
    sw.precompStructForRitz = ritzSAWPrecompute(structure, NamedTuple(), rawByDepthCell)
    sw.cs = [CalcSingleLayered(; structure = structure, angle = a, vMax = sw.vMax,
                               precompStruct = sw.precompStructForRitz, calcMode = sw.calcMode,
                               retentionLevel = sw.retentionLevel) for a in sw.angles]
    return sw
end

function ritzPrecompute(sw::CalcAngSweepLayered)
    sw.precompStructForRitz = ritzSAWPrecompute(sw.structure, sw.precompStructForRitz)
    return sw.precompStructForRitz
end

"""
    setElement(sweep, iAngle, cs) -> sweep

Replace one angle's object (the 2D sweeps' flat dispatch).
"""
function setElement(sw::CalcAngSweepLayered, iAngle::Integer, csElem::CalcSingleLayered)
    sw.cs[iAngle] = csElem
    return sw
end

function setGreenVelAxis(sw::CalcAngSweepLayered, axis)
    sw.greenVelAxis = collect(Float64, axis)
    return sw
end

"""
    commonGreenAxis(sweep) -> Vector{Float64}

The 4 m/s lattice from min over angles of 0.3 bandVSmin to max over angles of
max(vMax, 1.05 bandVLmax) (+ 0.1 km/s). Needs every angle's calcLim.
"""
function commonGreenAxis(sw::CalcAngSweepLayered)
    lo = minimum(GreenScoreLayered.BAND_LO_X_VSMIN * c.bandVSmin for c in sw.cs)
    hi = maximum(scoreAxisBound(c) for c in sw.cs)
    return GreenScoreLayered.gridAxis(lo, hi)
end

finalizeMaps(sw::CalcAngSweepLayered) = assembleMaps(sw)

"""
    calc(sweep::CalcAngSweepLayered, mode = nothing, retentionLevel = nothing) -> sweep
"""
function calc(sw::CalcAngSweepLayered, mode = nothing, retentionLevel = nothing)
    mode === nothing && (mode = sw.calcMode)
    retentionLevel === nothing && (retentionLevel = sw.retentionLevel)
    validateCalcMode(mode, "CalcAngSweepLayered:calc")
    validateRetentionLevel(retentionLevel, "CalcAngSweepLayered:calc", (:full, :saw_only, :none))
    method = resolveSawClassMethod(sw.sawClassMethod, sw.sawThreshold, mode, true, "CalcAngSweepLayered:calc")
    sw.calcMode = mode
    sw.retentionLevel = retentionLevel
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

function calcFullPhases(sw::CalcAngSweepLayered, precomp, method)
    useP = resolveParallel(sw.useParallel, length(sw.cs))
    sw.cs = forEachElement(sw.cs, useP, (c, _) -> begin    # independent per angle
        c.vMax = sw.vMax
        calcLim(c)
    end)
    sw.greenVelAxis = commonGreenAxis(sw)
    axis = sw.greenVelAxis; rl = sw.retentionLevel
    thr = method === :green ? sw.sawThreshold : nothing
    sw.cs = forEachElement(sw.cs, useP, (c, _) -> begin
        calcRitz(c, precomp)
        greenStep(c, axis, rl, nothing, thr, precomp)
    end)
    return sw
end

"""
    classify(sweep::CalcAngSweepLayered, method = nothing, threshold = nothing) -> sweep
"""
function classify(sw::CalcAngSweepLayered, method = nothing, threshold = nothing)
    threshold = normaliseThreshold(threshold)
    method === nothing && (method = sw.activeMethod)
    method === nothing && (method = sw.sawClassMethod)
    method === nothing && (method = resolveSawClassMethod(nothing, threshold, sw.calcMode, true, "CalcAngSweepLayered:classify"))
    method = validateSawClassMethod(method, true, "CalcAngSweepLayered:classify")
    threshold === nothing && (threshold = defaultThreshold(method, true))
    validateSawThreshold(threshold, "CalcAngSweepLayered:classify")
    refuseStale(sw, "CalcAngSweepLayered:classify")
    (isempty(sw.cs) || sw.cs[1].classification === nothing) &&
        throw(ArgumentError("CalcAngSweepLayered:classify:notCalculated -- calc() must run first."))
    method === :green && sw.calcMode !== :full && throw(ArgumentError("CalcAngSweepLayered:classify:needsFull -- " *
        "Method 'green' needs calcMode=:full; this sweep is '$(sw.calcMode)'."))
    useP = resolveParallel(sw.useParallel, length(sw.cs))
    sw.cs = forEachElement(sw.cs, useP, (c, _) -> classify(c, method, threshold))
    return assembleMaps(sw)
end

reclassify(sw::CalcAngSweepLayered, threshold) = classify(sw, sw.activeMethod, threshold)

greenComponentMap(sw::CalcAngSweepLayered, component) = greenComponentMapOf(sw, component)

function assembleMaps(sw::CalcAngSweepLayered)
    sw.sawVelMap, sw.sawModeNoMap, sw.sawTypeMap, sw.significantModesMap = sawMapsFrom(sw.cs, (:SAW,))
    if !isempty(sw.cs) && sw.cs[1].activeMethod !== nothing
        sw.activeMethod = sw.cs[1].activeMethod
        sw.activeThreshold = sw.cs[1].activeThreshold
    end
    if sw.calcMode === :full
        sw.subLongVelMap = [c.subLongLimVel for c in sw.cs]
        sw.subBulkVelMap = reduce(hcat, [c.subBulkVels for c in sw.cs])
        sw.greenAbsG13EffMap = greenG13EffMapFrom(sw.cs)
    end
    return sw
end
