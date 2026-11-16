# =========================================================================================
# CalcSingleLayered.jl -- port of guliwerr/forwardCalcs/CalcSingleLayered.m.
#
# Single-angle container for a layered film-on-substrate structure (nLay >= 2): RitzSolverSAW,
# the substrate's Christoffel reference, the whole-stack band, the one Green solve and the
# existence scores (GreenScoreLayered on 'full', RitzScore layered on 'ritzOnly'). Mutable,
# steps mutate and return `cs` (see CalcSingleHS.jl's header).
# =========================================================================================

"""
    CalcSingleLayered

Layered single-angle calculation (MATLAB `CalcSingleLayered`; fields as MATLAB):
`structure` (rotated), `angle`, `structureBase`, `calcMode`, `ritzSolver`; substrate-only
Christoffel reference `subLongLimVel`, `subLongLimPol`, `subShearLimVels`,
`subShearLimPols`, `subBulkVels`, `subBulkPols`; whole-stack band `bandVSmin`, `bandVLmax`;
Green `greenVelAxis`, `greenAbsG13Eff`, `greenScoreVelAxis`, `greenScoreG13Eff`,
`greenScorePeaks` (GreenScoreLayered); Ritz score `ritzScoreNum`, `ritzScoreVelAxis`,
`ritzScoreRef`, `ritzScoreCands`; classification `classification`, `classificationVersion`,
`activeMethod`, `activeThreshold`, `limitingEntries`, `limitingDone`; settings
`retentionLevel`, `sawClassMethod`, `sawThreshold`, `vMax` (NaN: the band's own end);
`precompStructForRitz`.
"""
mutable struct CalcSingleLayered
    structure::Structure
    angle::Float64
    structureBase::Structure
    calcMode::Symbol
    ritzSolver::RitzSolverSAW
    subLongLimVel::Float64
    subLongLimPol::Vector{Float64}
    subShearLimVels::Vector{Float64}
    subShearLimPols::Matrix{Float64}
    subBulkVels::Vector{Float64}
    subBulkPols::Matrix{Float64}
    bandVSmin::Float64
    bandVLmax::Float64
    greenVelAxis::Vector{Float64}
    greenAbsG13Eff::Vector{Float32}
    greenScoreVelAxis::Vector{Float64}
    greenScoreG13Eff::Vector{ComplexF64}
    greenScorePeaks::Vector{GreenScoreLayered.Peak}
    greenScoreDone::Bool
    ritzScoreNum::Vector{ComplexF64}
    ritzScoreVelAxis::Vector{Float64}
    ritzScoreRef::Union{Nothing,RitzScore.OwnReference}
    ritzScoreCands::Vector{RitzScore.Peak}
    classification::Union{Nothing,ModeClassification}
    classificationVersion::Int
    activeMethod::Union{Nothing,Symbol}
    activeThreshold::Union{Nothing,Symbol,Float64}
    limitingEntries::Vector{SignificantMode}
    limitingDone::Bool
    retentionLevel::Symbol
    sawClassMethod::Union{Nothing,Symbol}
    sawThreshold::Union{Nothing,Symbol,Float64}
    vMax::Float64
    precompStructForRitz::NamedTuple
end

Base.copy(cs::CalcSingleLayered) =
    CalcSingleLayered((getfield(cs, f) for f in fieldnames(CalcSingleLayered))...)

function _emptyCalcSingleLayered(structure, angle, structureBase, precomp)
    return CalcSingleLayered(structure, angle, structureBase, :full, RitzSolverSAW(),
                             NaN, Float64[], Float64[], zeros(3, 0), Float64[], zeros(3, 0),
                             NaN, NaN, Float64[], Float32[], Float64[], ComplexF64[],
                             GreenScoreLayered.Peak[], false, ComplexF64[], Float64[], nothing,
                             RitzScore.Peak[], nothing, 0, nothing, nothing, SignificantMode[],
                             false, :full, nothing, nothing, NaN, precomp)
end

"""
    CalcSingleLayered(; structure = nothing, angle = 0, vMax = nothing,
                      precompStruct = NamedTuple(), calcMode = nothing, retentionLevel = nothing)

One propagation angle of an UNROTATED layered structure (`nLay >= 2`; rotated here, once).
`vMax` [km/s]: extends the own Green axis beyond 1.05 x the fastest longitudinal limit of
any layer (NaN / nothing: no extension). No structure: an empty object.
"""
function CalcSingleLayered(; structure::Union{Structure,Nothing} = nothing, angle::Real = 0.0,
                           vMax = nothing, precompStruct::NamedTuple = NamedTuple(),
                           calcMode = nothing, retentionLevel = nothing)
    structure === nothing &&
        return _emptyCalcSingleLayered(Structure(), 0.0, Structure(), NamedTuple())
    nlay(structure) < 2 && throw(ArgumentError("CalcSingleLayered:constructor:needsFilmAndSubstrate -- " *
        "CalcSingleLayered needs at least one film plus a substrate (nLay >= 2); use CalcSingleHS " *
        "for a homogeneous half-space."))
    isfinite(angle) || throw(ArgumentError("CalcSingleLayered: angle must be finite; got $angle"))
    validateConstructorFlags(calcMode, retentionLevel, "CalcSingleLayered:constructor", (:full, :saw_only, :none))
    sr = rotateStructure(structure, angle)
    cs = _emptyCalcSingleLayered(sr, Float64(angle), structure, ritzSAWPrecompute(sr, precompStruct))
    vMax === nothing || (cs.vMax = Float64(vMax))
    calcMode === nothing || (cs.calcMode = calcMode)
    retentionLevel === nothing || (cs.retentionLevel = retentionLevel)
    return cs
end

"""
    calcRitz(cs, precomp = nothing) -> cs

RitzSolverSAW of the stack at this angle (angle 0 on the rotated structure), diagnostics on.
"""
function calcRitz(cs::CalcSingleLayered, precomp = nothing)
    p = precomp === nothing ? ritzPrecompute(cs) : completeRitzPrecompute(cs.structure, precomp)
    cs.ritzSolver = calc(RitzSolverSAW(; structure = cs.structure, angleDeg = 0.0), true, p)
    cs.limitingEntries = SignificantMode[]; cs.limitingDone = false
    return cs
end

"""
    calcLim(cs) -> cs

The SUBSTRATE's (`mater[1]`) Christoffel limits (a soft reference: the leaky-mode label and
the limiting entries, not a system-wide bound) and the whole-stack band
(`GreenScoreLayered.bandLimits`: slowest shear and fastest longitudinal limit of any layer).
"""
function calcLim(cs::CalcSingleLayered)
    # every layer's limits once: the substrate's (layer 1) are also one of the band's
    Ls  = [limitingAndBulk(m) for m in cs.structure.mater]
    lim = Ls[1]
    cs.subLongLimVel = lim.longLimVel; cs.subLongLimPol = lim.longLimPol
    cs.subShearLimVels = lim.shearLimVels; cs.subShearLimPols = lim.shearLimPols
    cs.subBulkVels = lim.bulkVels; cs.subBulkPols = lim.bulkPols
    cs.bandVSmin, cs.bandVLmax = GreenScoreLayered.bandLimits(Ls)
    return cs
end

scoreAxisBound(cs::CalcSingleLayered) =
    isnan(cs.vMax) ? GreenScoreLayered.BAND_HI_X_VLMAX * cs.bandVLmax :
                     max(cs.vMax, GreenScoreLayered.BAND_HI_X_VLMAX * cs.bandVLmax)

"""
    greenAxis(cs) -> Vector{Float64}

This angle's own Green axis: the 4 m/s lattice from the band's lower end
(0.3 bandVSmin) to max(vMax, 1.05 bandVLmax) + 0.1 km/s. Needs calcLim.
"""
function greenAxis(cs::CalcSingleLayered)
    isnan(cs.bandVSmin) && throw(ArgumentError("CalcSingleLayered:greenAxis:noLim -- calcLim() " *
        "must run first: the axis runs to 1.05 x the fastest longitudinal limit of any layer."))
    vLo, _ = GreenScoreLayered.band(cs.bandVSmin, cs.bandVLmax)
    return GreenScoreLayered.gridAxis(vLo, scoreAxisBound(cs))
end

"""
    calcGreen(cs, velAxis = nothing; threaded = true) -> cs

The one Green solve (G13_eff at SawSpectrum.ETA), keeping the band
`[0.3 bandVSmin, 1.05 bandVLmax] (+- 1.5 steps)` for the score when the axis covers it.
"""
function calcGreen(cs::CalcSingleLayered, velAxis = nothing; threaded::Bool = true)
    v0 = velAxis === nothing ? greenAxis(cs) : velAxis
    g = greenScoreSpectrum(cs.structure, v0, "CalcSingleLayered.calcGreen"; threaded = threaded)
    cs.greenVelAxis = g.velAxis
    cs.greenAbsG13Eff = g.absG13Eff
    cs.greenScoreVelAxis = Float64[]; cs.greenScoreG13Eff = ComplexF64[]
    cs.greenScorePeaks = GreenScoreLayered.Peak[]; cs.greenScoreDone = false
    isnan(cs.bandVSmin) && return cs
    vLo, vHi = GreenScoreLayered.band(cs.bandVSmin, cs.bandVLmax)
    v = g.velAxis; dv = v[2] - v[1]
    b = findall(x -> x >= vLo - 1.5 * dv && x <= vHi + 1.5 * dv, v)
    covers = v[1] <= vLo && v[end] >= vHi
    if covers && length(b) >= 5 && maximum(abs.(diff(v[b]) .- dv)) < 1e-9
        cs.greenScoreVelAxis = v[b]
        cs.greenScoreG13Eff = g.G13Eff[b]
    end
    return cs
end

greenComponent(cs::CalcSingleLayered, component) = greenComponentOf(cs, component)

"""
    calcGreenScore(cs) -> cs

GreenScoreLayered candidates (probe slopes, polarisation, scores; no reference). Solves the
band itself if calcGreen left none.
"""
function calcGreenScore(cs::CalcSingleLayered)
    isnan(cs.bandVSmin) && throw(ArgumentError("CalcSingleLayered:calcGreenScore:noLim -- calcLim() must run first."))
    if isempty(cs.greenScoreG13Eff)
        vLo, vHi = GreenScoreLayered.band(cs.bandVSmin, cs.bandVLmax)
        v = GreenScoreLayered.gridAxis(vLo, vHi); dv = v[2] - v[1]
        v = v[(v .>= vLo - 1.5 * dv) .& (v .<= vHi + 1.5 * dv)]
        g = greenScoreSpectrum(cs.structure, v, "CalcSingleLayered.calcGreenScore")
        cs.greenScoreVelAxis = g.velAxis; cs.greenScoreG13Eff = g.G13Eff
    end
    cs.greenScorePeaks = GreenScoreLayered.direction(cs.structureBase, cs.angle, cs.greenScoreG13Eff,
                                                     cs.greenScoreVelAxis, cs.bandVSmin, cs.bandVLmax)
    cs.greenScoreDone = true
    return cs
end

"""
    calcRitzScore(cs, precomp = nothing) -> cs

Layered Ritz-only existence input: modal-sum weights of the modes up to 1.2 bandVLmax, the
band axis, the own floor and the scored candidates (RitzScore.directionLayered).
"""
function calcRitzScore(cs::CalcSingleLayered, precomp = nothing)
    isnan(cs.bandVLmax) && throw(ArgumentError("CalcSingleLayered:calcRitzScore:noLim -- calcLim() must run first."))
    precomp === nothing && (precomp = ritzPrecompute(cs))
    dKdC = precomp_field(precomp, :dKdC_cell)
    vLo, vHi = GreenScoreLayered.band(cs.bandVSmin, cs.bandVLmax)
    prof = RitzScore.modeProfiles(cs.ritzSolver, cs.bandVLmax)
    cs.ritzScoreNum = prof.num
    cs.ritzScoreVelAxis = RitzScore.layeredAxis(vLo, vHi)
    cs.ritzScoreRef = RitzScore.ownReference(cs.structure, cs.ritzSolver.vels[1:length(prof.num)],
                                             prof.num, cs.ritzScoreVelAxis, vLo, vHi)
    cs.ritzScoreCands = RitzScore.directionLayered(cs.ritzSolver, prof, cs.ritzScoreVelAxis,
                                                   cs.ritzScoreRef, vLo, vHi,
                                                   minimum(cs.subShearLimVels), dKdC)
    return cs
end

"""
    classify(cs, method = nothing, threshold = nothing) -> cs

`:green` | `:ritz` (no 'refined' for layered stacks) at `:high` | `:low` | a number (default
'low'); re-thresholds the stored candidates.
"""
function classify(cs::CalcSingleLayered, method = nothing, threshold = nothing)
    threshold = normaliseThreshold(threshold)
    method === nothing && (method = cs.activeMethod)
    method === nothing && (method = cs.sawClassMethod)
    method === nothing && (method = resolveSawClassMethod(nothing, threshold, cs.calcMode, true, "CalcSingleLayered:classify"))
    method = validateSawClassMethod(method, true, "CalcSingleLayered:classify")
    threshold === nothing && (threshold = defaultThreshold(method, true))
    validateSawThreshold(threshold, "CalcSingleLayered:classify")
    refuseStale(cs, "CalcSingleLayered:classify")
    method === :green && cs.calcMode !== :full && throw(ArgumentError("CalcSingleLayered:classify:needsFull -- " *
        "Method 'green' needs calcMode=:full; this result is '$(cs.calcMode)'."))
    return identifyModes(cs, method, threshold)
end

reclassify(cs::CalcSingleLayered, threshold) = classify(cs, cs.activeMethod, threshold)

"""
    identifyModes(cs, method = nothing, thr = nothing) -> cs

Classify `:green` (GreenScoreLayered candidates, substrate limiting entries) or `:ritz`.
"""
function identifyModes(cs::CalcSingleLayered, method = nothing, thr = nothing)
    thr = thr === nothing ? cs.sawThreshold : normaliseThreshold(thr)
    method === nothing && (method = cs.sawClassMethod)
    method === nothing && (method = resolveSawClassMethod(nothing, thr, cs.calcMode, true, "CalcSingleLayered:identifyModes"))
    method = validateSawClassMethod(method, true, "CalcSingleLayered:identifyModes")
    thr === nothing && (thr = defaultThreshold(method, true))
    if method === :ritz
        (isempty(cs.ritzScoreVelAxis) || isnan(cs.subLongLimVel)) && throw(ArgumentError(
            "CalcSingleLayered:identifyModes:noRitzScore -- Ritz classification needs calcLim() " *
            "and calcRitzScore() first (calcAll runs both)."))
        cs.classification = calc(ModeClassification(threshold = thr), cs.ritzSolver, cs.structure,
                                 nothing, nothing, (cands = cs.ritzScoreCands,))
    else
        cs.greenScoreDone || throw(ArgumentError("CalcSingleLayered:identifyModes:noScore -- " *
            "calcGreenScore() must run first (calcAll(:full) runs it)."))
        cr = (subLongLimVel = cs.subLongLimVel, subLongLimPol = cs.subLongLimPol,
              subShearLimVels = cs.subShearLimVels, subShearLimPols = cs.subShearLimPols)
        gr = (greenScoreLayered = (peaks = cs.greenScorePeaks,),)
        cs.classification, cs.limitingEntries = classifyGreenPath(cs, thr, cr, gr)
        cs.limitingDone = true
    end
    cs.activeMethod = method
    cs.activeThreshold = thr
    cs.classificationVersion = CLASSIFICATION_VERSION
    return cs
end

"""
    calcAll(cs, mode = nothing, retentionLevel = nothing, precomp = nothing) -> cs

As CalcSingleHS.calcAll; `:full` runs calcGreenScore without a reference (GreenScoreLayered
needs none).
"""
function calcAll(cs::CalcSingleLayered, mode = nothing, retentionLevel = nothing, precomp = nothing)
    mode === nothing && (mode = cs.calcMode)
    retentionLevel === nothing && (retentionLevel = cs.retentionLevel)
    validateCalcMode(mode, "CalcSingleLayered:calcAll")
    validateRetentionLevel(retentionLevel, "CalcSingleLayered:calcAll", (:full, :saw_only, :none))
    target = resolveSawClassMethod(cs.sawClassMethod, cs.sawThreshold, mode, true, "CalcSingleLayered:calcAll")
    targetThr = cs.sawThreshold
    cs.calcMode = mode
    calcRitz(cs, precomp)
    cs.classificationVersion = 0
    calcLim(cs)
    if mode === :full
        calcGreen(cs)
        calcGreenScore(cs)
        calcRitzScore(cs, precomp)
        identifyModes(cs, :green, target === :green ? targetThr : nothing)
        target === :green || classify(cs, target, targetThr)
    else
        calcRitzScore(cs, precomp)
        identifyModes(cs, target, targetThr)
    end
    applyRetention(cs, retentionLevel)
    cs.retentionLevel = retentionLevel
    return cs
end

function applyRetention(cs::CalcSingleLayered, retentionLevel)
    cs.ritzSolver = pruneRitzTo(cs.ritzSolver, cs.classification, retentionLevel, "CalcSingleLayered:applyRetention")
    return cs
end

function setRunMode(cs::CalcSingleLayered, mode)
    validateCalcMode(mode, "CalcSingleLayered:setRunMode")
    cs.calcMode = mode
    return cs
end
