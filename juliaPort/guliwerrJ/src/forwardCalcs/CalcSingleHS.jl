# =========================================================================================
# CalcSingleHS.jl -- port of guliwerr/forwardCalcs/CalcSingleHS.m.
#
# Single-angle container for the homogeneous half-space forward calculation (RitzSolverSAW +
# Christoffel + the one Green solve + the existence scores). The per-angle building block of
# CalcAngSweepHS.
#
# MUTABLE struct: MATLAB's class is a handle whose steps fill fields
# one by one (calcRitz -> calcLim -> calcGreen -> calcGreenScore -> calcRitzScore ->
# identifyModes / classify -> applyRetention), and the sweeps set per-angle properties in
# place. Every step mutates `cs` and returns it, so both `calcAll(cs)` and
# `cs = calcAll(cs)` work, as in MATLAB. A sweep's angles are distinct objects, so a
# threaded loop over them writes disjoint memory. `copy(cs)` is a shallow copy (MATLAB
# copyElement; every field is a value or immutable).
#
# Representation: calcMode / retentionLevel / methods are Symbols; MATLAB's '' / [] is
# `nothing` (sawClassMethod, sawThreshold, activeMethod, activeThreshold, the references).
#
# calcSpotGreen (the experimental SpotGreen picture) is ported:
# `spotGreen`, `spotGreenVelAxis`, `spotGreenInfo` are filled on request, and kept by
# calcAll/classify (they need only the structure). Not ported: the .mat I/O.
# =========================================================================================

"""
    RefinedSaw

`calcRefined` result (MATLAB `refinedSaw` struct): `threshold` (the standard the peaks were
last accepted at), `ref` (GreenSpotScore.Reference), `peaks` (every candidate after
`GreenSpotScore.identify`, `accepted` filled), `info` (`nVel`, `velEnd`, `minShearSpot`).
"""
struct RefinedSaw
    threshold::Union{Symbol,Float64}
    ref::GreenSpotScore.Reference
    peaks::Vector{GreenSpotScore.Peak}
    info::NamedTuple{(:nVel, :velEnd, :minShearSpot),Tuple{Int,Float64,Float64}}
end

"""
    CalcSingleHS

Half-space single-angle calculation (MATLAB `CalcSingleHS`; fields as MATLAB):
`structure` (rotated), `angle`, `structureBase` (unrotated), `calcMode`, `ritzSolver`,
Christoffel `longLimVel`, `longLimPol`, `shearLimVels`, `shearLimPols`, `shearLimThetas`,
`bulkVels`, `bulkPols`; Green `greenVelAxis`, `greenAbsG13Eff` (Float32); GreenScore input and
output `greenScoreVelAxis`, `greenScoreG13Eff`, `greenScorePeaks`, `greenScoreRef`; Ritz-only
score `ritzScoreNum`, `ritzScoreVelAxis`, `ritzScoreRef`, `ritzScoreCands`; classification
`classification`, `classificationVersion`, `activeMethod`, `activeThreshold`,
`limitingEntries`, `limitingDone`, `refinedSaw`; the spot-averaged Green curve of
[`calcSpotGreen`](@ref) (EXPERIMENTAL) `spotGreen` (Float32), `spotGreenVelAxis`,
`spotGreenInfo` (`nothing` until then); settings `retentionLevel`, `sawClassMethod`,
`sawThreshold`, `greenVelStep`, `nVelGreen`; `precompStructForRitz`.
"""
mutable struct CalcSingleHS
    structure::Structure
    angle::Float64
    structureBase::Structure
    calcMode::Symbol
    ritzSolver::RitzSolverSAW
    longLimVel::Float64
    longLimPol::Vector{Float64}
    shearLimVels::Vector{Float64}
    shearLimPols::Matrix{Float64}
    shearLimThetas::Vector{Float64}
    bulkVels::Vector{Float64}
    bulkPols::Matrix{Float64}
    greenVelAxis::Vector{Float64}
    greenAbsG13Eff::Vector{Float32}
    greenScoreVelAxis::Vector{Float64}
    greenScoreG13Eff::Vector{ComplexF64}
    greenScorePeaks::Vector{GreenScore.Peak}
    greenScoreRef::Union{Nothing,GreenScore.Reference}
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
    refinedSaw::Union{Nothing,RefinedSaw}
    spotGreen::Vector{Float32}
    spotGreenVelAxis::Vector{Float64}
    spotGreenInfo::Union{Nothing,NamedTuple}
    retentionLevel::Symbol
    sawClassMethod::Union{Nothing,Symbol}
    sawThreshold::Union{Nothing,Symbol,Float64}
    greenVelStep::Float64
    nVelGreen::Float64
    precompStructForRitz::NamedTuple
end

Base.copy(cs::CalcSingleHS) = CalcSingleHS((getfield(cs, f) for f in fieldnames(CalcSingleHS))...)

function _emptyCalcSingleHS(structure, angle, structureBase, precomp)
    return CalcSingleHS(structure, angle, structureBase, :full, RitzSolverSAW(),
                        NaN, Float64[], Float64[], zeros(3, 0), Float64[], Float64[], zeros(3, 0),
                        Float64[], Float32[], Float64[], ComplexF64[], GreenScore.Peak[], nothing,
                        ComplexF64[], Float64[], nothing, RitzScore.Peak[],
                        nothing, 0, nothing, nothing, SignificantMode[], false, nothing,
                        Float32[], Float64[], nothing,
                        :full, nothing, nothing, GreenScore.DEFAULT_STEP, NaN, precomp)
end

"""
    CalcSingleHS(; structure = nothing, angle = 0, precompStruct = NamedTuple(),
                 calcMode = nothing, retentionLevel = nothing, nVelGreen = nothing)

One propagation angle of an UNROTATED half-space structure (rotated here, once). The Ritz
precompute is completed from `precompStruct` (a sweep passes its shared one). `nVelGreen`:
number of points of this angle's Green axis (sets the step, capped at GreenScore.MAX_STEP);
default: the `greenVelStep` field (GreenScore.DEFAULT_STEP). No structure: an empty object.
"""
function CalcSingleHS(; structure::Union{Structure,Nothing} = nothing, angle::Real = 0.0,
                      precompStruct::NamedTuple = NamedTuple(), calcMode = nothing,
                      retentionLevel = nothing, nVelGreen = nothing)
    structure === nothing && return _emptyCalcSingleHS(Structure(), 0.0, Structure(), NamedTuple())
    isfinite(angle) || throw(ArgumentError("CalcSingleHS: angle must be finite; got $angle"))
    validateConstructorFlags(calcMode, retentionLevel, "CalcSingleHS:constructor", (:full, :saw_only, :none))
    sr = rotateStructure(structure, angle)
    cs = _emptyCalcSingleHS(sr, Float64(angle), structure, ritzSAWPrecompute(sr, precompStruct))
    calcMode === nothing || (cs.calcMode = calcMode)
    retentionLevel === nothing || (cs.retentionLevel = retentionLevel)
    if nVelGreen !== nothing
        (nVelGreen isa Real && nVelGreen >= 2) ||
            throw(ArgumentError("CalcSingleHS: nVelGreen must be a real scalar >= 2; got $nVelGreen"))
        cs.nVelGreen = Float64(nVelGreen)
    end
    return cs
end

"""
    calcRitz(cs, precomp = nothing) -> cs

RitzSolverSAW for this angle (angle 0 on the rotated structure), with the per-mode
diagnostics. `precomp`: a sweep's shared precompute (completed if partial); default the
object's own. Clears the refined candidates and the limiting entries (new eigenvectors).
"""
function calcRitz(cs::CalcSingleHS, precomp = nothing)
    p = precomp === nothing ? ritzPrecompute(cs) : completeRitzPrecompute(cs.structure, precomp)
    cs.ritzSolver = calc(RitzSolverSAW(; structure = cs.structure, angleDeg = 0.0), true, p)
    cs.refinedSaw = nothing
    cs.limitingEntries = SignificantMode[]; cs.limitingDone = false
    return cs
end

"""
    ritzPrecompute(cs) -> NamedTuple

The object's Ritz precompute, completed (and stored back) if partial.
"""
function ritzPrecompute(cs)
    cs.precompStructForRitz = completeRitzPrecompute(cs.structure, cs.precompStructForRitz)
    return cs.precompStructForRitz
end

"""
    calcLim(cs) -> cs

Christoffel limiting velocities of the medium (the surface layer, `mater[nLay]`: for a
half-space the whole medium) and the bulk velocities.
"""
function calcLim(cs::CalcSingleHS)
    lim = limitingAndBulk(cs.structure.mater[nlay(cs.structure)])
    cs.longLimVel = lim.longLimVel; cs.longLimPol = lim.longLimPol
    cs.shearLimVels = lim.shearLimVels; cs.shearLimPols = lim.shearLimPols
    cs.shearLimThetas = lim.shearLimTheta
    cs.bulkVels = lim.bulkVels; cs.bulkPols = lim.bulkPols
    return cs
end

shearRange(cs::CalcSingleHS) = (minimum(cs.shearLimVels), maximum(cs.shearLimVels))
scoreAxisBound(cs::CalcSingleHS) = SawSpectrum.AXIS_HI_X_VL * cs.longLimVel

"""
    greenAxis(cs) -> Vector{Float64}

This angle's own Green axis: 0.3 km/s to 1.05 longLimVel + 0.1 km/s in `greenVelStep` steps
(or the step giving `nVelGreen` points), capped at GreenScore.MAX_STEP. Needs calcLim.
"""
function greenAxis(cs::CalcSingleHS)
    isnan(cs.longLimVel) && throw(ArgumentError("CalcSingleHS:greenAxis:noLim -- calcLim() " *
        "must run first: the axis runs to 1.05 x longLimVel."))
    vBound = scoreAxisBound(cs)
    step = cs.greenVelStep
    isfinite(cs.nVelGreen) &&
        (step = (vBound + SawSpectrum.EDGE_TAPER - SawSpectrum.VEL_MIN) / (cs.nVelGreen - 1))
    return GreenScore.gridAxis(vBound, GreenScore.resolveStep(step))
end

"""
    calcGreen(cs, velAxis = nothing; threaded = true) -> cs

The one Green solve of this angle (G13_eff at SawSpectrum.ETA): stores `greenVelAxis`,
`greenAbsG13Eff` and, on a uniform axis covering it, the band
`[BAND_LO vS, BAND_SOLVE_HI vTf] (+- 1.5 steps)` the existence score reads.
"""
function calcGreen(cs::CalcSingleHS, velAxis = nothing; threaded::Bool = true)
    v0 = velAxis === nothing ? greenAxis(cs) : velAxis
    g = greenScoreSpectrum(cs.structure, v0, "CalcSingleHS.calcGreen"; threaded = threaded)
    cs.refinedSaw = nothing
    cs.greenVelAxis = g.velAxis
    cs.greenAbsG13Eff = g.absG13Eff
    cs.greenScoreVelAxis = Float64[]; cs.greenScoreG13Eff = ComplexF64[]
    isnan(cs.longLimVel) && return cs
    v = g.velAxis; dv = v[2] - v[1]
    vS, vTf = shearRange(cs)
    b = findall(x -> x >= GreenScore.BAND_LO * vS - 1.5 * dv &&
                     x <= GreenScore.BAND_SOLVE_HI * vTf + 1.5 * dv, v)
    if length(b) >= 5 && maximum(abs.(diff(v[b]) .- dv)) < 1e-9
        cs.greenScoreVelAxis = v[b]
        cs.greenScoreG13Eff = g.G13Eff[b]
    end
    return cs
end

"""
    calcGreenScore(cs, ref = nothing) -> cs

GreenScore candidates of this angle (probe slopes, scores, the 'low' reference). `ref`: a
sweep's shared GreenScore.Reference; default that of `structureBase` (memoised).
"""
function calcGreenScore(cs::CalcSingleHS, ref = nothing)
    isnan(cs.longLimVel) && throw(ArgumentError("CalcSingleHS:calcGreenScore:noLim -- calcLim() must run first."))
    ref === nothing && (ref = GreenScore.reference(cs.structureBase))
    cs.greenScoreRef = ref
    vS, vTf = shearRange(cs)
    if isempty(cs.greenScoreG13Eff)
        v = GreenScore.bandAxis(vS, vTf, GreenScore.resolveStep(cs.greenVelStep))
        g = greenScoreSpectrum(cs.structure, v, "CalcSingleHS.calcGreenScore")
        cs.greenScoreVelAxis = g.velAxis; cs.greenScoreG13Eff = g.G13Eff
    end
    cs.greenScorePeaks = GreenScore.direction(cs.structureBase, cs.angle, cs.greenScoreG13Eff,
                                              cs.greenScoreVelAxis, vS, vTf, ref.vRef, ref)
    return cs
end

"""
    greenComponent(cs, component) -> Vector{Float32}

|G| of any component (`"G13eff"`, `"G33"`, `"G13"`, ...) on greenVelAxis; only G13eff is
stored, any other is one solve now.
"""
greenComponent(cs::CalcSingleHS, component) = greenComponentOf(cs, component)

"""
    calcRitzScore(cs, precomp = nothing) -> cs

Ritz-only existence input: modal-sum weights, the uniform score axis, this direction's own
floor (+ slab resolution), and the scored candidates (RitzScore.direction). Needs calcRitz
(eigenvectors present) and calcLim.
"""
function calcRitzScore(cs::CalcSingleHS, precomp = nothing)
    isnan(cs.longLimVel) && throw(ArgumentError("CalcSingleHS:calcRitzScore:noLim -- calcLim() must run first."))
    precomp === nothing && (precomp = ritzPrecompute(cs))
    dKdC = precomp_field(precomp, :dKdC_cell)
    vS, vTf = shearRange(cs); vL = cs.longLimVel
    RitzScore.checkBasisHS(cs.structure.calcParams[1])
    prof = RitzScore.modeProfiles(cs.ritzSolver, vL)
    cs.ritzScoreNum = prof.num
    cs.ritzScoreVelAxis = SawSpectrum.velocityAxis(SawSpectrum.AXIS_HI_X_VL * vL, RitzScore.VEL_STEP)
    ref = RitzScore.ownReference(cs.structure, cs.ritzSolver.vels[1:length(prof.num)], prof.num,
                                 cs.ritzScoreVelAxis, RitzScore.HS_FLOOR_LO * vS, vTf)
    iS = argmin(cs.shearLimVels)
    ref = RitzScore.withSlabRes(ref, RitzScore.slabResolution(cs.structure.mater[1], vS,
                                cs.shearLimThetas[iS], cs.structure.calcParams[1].relDepth))
    cs.ritzScoreRef = ref
    cs.ritzScoreCands = RitzScore.direction(cs.ritzSolver, prof, cs.ritzScoreVelAxis, ref, vS, vL, dKdC)
    return cs
end

"""
    classify(cs, method = nothing, threshold = nothing) -> cs

Set the active classification (re-thresholds the stored candidates, no recomputation except
'refined' on first use). `method`: `:green` | `:refined` | `:ritz` (default: the active one);
`threshold`: `:high` | `:low` | a number (default: the method's default standard). Stale
results are refused.
"""
function classify(cs::CalcSingleHS, method = nothing, threshold = nothing)
    threshold = normaliseThreshold(threshold)
    method === nothing && (method = cs.activeMethod)
    method === nothing && (method = cs.sawClassMethod)
    method === nothing && (method = resolveSawClassMethod(nothing, threshold, cs.calcMode, false, "CalcSingleHS:classify"))
    method = validateSawClassMethod(method, false, "CalcSingleHS:classify")
    threshold === nothing && (threshold = defaultThreshold(method, false))
    validateSawThreshold(threshold, "CalcSingleHS:classify")
    refuseStale(cs, "CalcSingleHS:classify")
    method !== :ritz && cs.calcMode !== :full && throw(ArgumentError("CalcSingleHS:classify:needsFull -- " *
        "Method '$method' needs calcMode=:full; this result is '$(cs.calcMode)'."))
    method === :refined && cs.refinedSaw === nothing && calcRefined(cs)
    return identifyModes(cs, method, threshold)
end

"""
    reclassify(cs, threshold) -> cs

`classify(cs, cs.activeMethod, threshold)`.
"""
reclassify(cs::CalcSingleHS, threshold) = classify(cs, cs.activeMethod, threshold)

"""
    calcRefined(cs, dirCache = nothing, ref = nothing) -> cs

Refined (GreenSpotScore) candidates in `refinedSaw`: the 11 kernel directions are solved
(or read from a sweep's `dirCache`), averaged, scored, identified with the Ritz matcher and
accepted at 'high'. Identity: each spot peak via its nearest single-direction pole
(`greenScorePeaks`), so `calcAll(:full)` must have run. 'full' half-space only.
"""
function calcRefined(cs::CalcSingleHS, dirCache = nothing, ref = nothing)
    cs.calcMode === :full || throw(ArgumentError("CalcSingleHS:calcRefined:notFull -- The refined " *
        "classification needs calcMode :full; this result is '$(cs.calcMode)'."))
    nlay(cs.structure) == 1 || throw(ArgumentError("CalcSingleHS:calcRefined:layered -- The " *
        "refined classification is for a half-space (nLay = 1) only."))
    (isnan(cs.longLimVel) || isempty(cs.ritzSolver.vels)) && throw(ArgumentError(
        "CalcSingleHS:calcRefined:notCalculated -- calcAll(:full) must run first."))
    (isempty(cs.greenScoreVelAxis) || cs.greenScoreRef === nothing) && throw(ArgumentError(
        "CalcSingleHS:calcRefined:noGreenScore -- calcGreen() and calcGreenScore() must run first " *
        "(calcAll(:full) runs both): the spot peaks are identified via the single-direction poles."))
    ref === nothing && (ref = GreenSpotScore.reference(cs.structureBase))
    if dirCache === nothing
        o, _ = GreenSpotScore.kernel()
        dirCache = GreenSpotScore.solveDirections(cs.structureBase, cs.angle .+ o; useParallel = :never)
    end
    X, v, minShearSpot = GreenSpotScore.average(dirCache, cs.angle)
    vS, vTf = shearRange(cs)
    pk = GreenSpotScore.direction(X, v, vS, vTf, cs.longLimVel, ref.vRef, ref)
    mc = ModeClassification()
    poleVels = [q.vel for q in cs.greenScorePeaks]   # the single-direction poles of this angle
    pk = GreenSpotScore.identify(pk, vel -> matchRitzToGreenPeak(cs.ritzSolver, vel, mc, false),
                                 vel -> matchSpotPeak(cs.ritzSolver, vel, mc), minShearSpot, poleVels)
    acc = GreenSpotScore.accepted(pk, :high)
    pk = GreenSpotScore.Peak[_with(pk[j]; accepted = acc[j]) for j in eachindex(pk)]
    cs.refinedSaw = RefinedSaw(:high, ref, pk, (nVel = length(v), velEnd = Float64(v[end]),
                                               minShearSpot = Float64(minShearSpot)))
    return cs
end

"""
    refinedSawAgreement(cs) -> Vector{NamedTuple}

Default ('green') vs refined surface-wave entries of this angle, paired one to one (nearest
first, within max(1.2 % v, 30 m/s)): rows `(angle, defaultVel, refinedVel, dv, defaultType,
refinedType, status)`, status `:both` | `:defaultOnly` | `:refinedOnly`.
"""
function refinedSawAgreement(cs::CalcSingleHS)
    cs.refinedSaw === nothing && throw(ArgumentError("CalcSingleHS:refinedSawAgreement:notRefined -- " *
        "classify(:refined) must run first."))
    thr = cs.refinedSaw.threshold
    thr isa Symbol || (thr = nothing)
    cp = classify(copy(cs), :green, thr)
    sm = [m for m in cp.classification.significantModes if m.type in (:Rayleigh, :pseudoRayleigh)]
    pk = [p for p in cs.refinedSaw.peaks if p.accepted]
    return pairEntries(cs.angle, [m.greenPeakVel for m in sm], [m.type for m in sm],
                       [p.vel for p in pk], [Symbol(p.type) for p in pk])
end

function pairEntries(angle, vA, tA, vB, tB)
    D = [abs(vB[b] - vA[a]) for a in eachindex(vA), b in eachindex(vB)]
    for a in eachindex(vA), b in eachindex(vB)
        D[a, b] > max(GreenSpotScore.POLE_WIN_REL * vA[a], GreenSpotScore.POLE_WIN_ABS) && (D[a, b] = Inf)
    end
    iA = Int[]; iB = Int[]
    while any(isfinite, D)
        k = argmin(D); a, b = Tuple(k)            # column-major first minimum, as MATLAB min(D(:))
        push!(iA, a); push!(iB, b)
        D[a, :] .= Inf; D[:, b] .= Inf
    end
    rows = NamedTuple{(:angle, :defaultVel, :refinedVel, :dv, :defaultType, :refinedType, :status),
                      Tuple{Float64,Float64,Float64,Float64,Symbol,Symbol,Symbol}}[]
    for (a, b) in zip(iA, iB)
        push!(rows, (angle = angle, defaultVel = vA[a], refinedVel = vB[b], dv = vB[b] - vA[a],
                     defaultType = tA[a], refinedType = tB[b], status = :both))
    end
    for a in setdiff(eachindex(vA), iA)
        push!(rows, (angle = angle, defaultVel = vA[a], refinedVel = NaN, dv = NaN,
                     defaultType = tA[a], refinedType = :none, status = :defaultOnly))
    end
    for b in setdiff(eachindex(vB), iB)
        push!(rows, (angle = angle, defaultVel = NaN, refinedVel = vB[b], dv = NaN,
                     defaultType = :none, refinedType = tB[b], status = :refinedOnly))
    end
    return rows
end

"""
    identifyModes(cs, method = nothing, thr = nothing) -> cs

Classify with a method at a standard (`ModeClassification`): `:green` (needs
calcGreenScore), `:refined` (needs calcRefined), `:ritz` (needs calcLim + calcRitzScore).
The limiting entries are associated on the first Green-path classification after a
calculation and carried by every later one. Stamps `classificationVersion`.
"""
function identifyModes(cs::CalcSingleHS, method = nothing, thr = nothing)
    thr = thr === nothing ? cs.sawThreshold : normaliseThreshold(thr)
    method === nothing && (method = cs.sawClassMethod)
    method === nothing && (method = resolveSawClassMethod(nothing, thr, cs.calcMode, false, "CalcSingleHS:identifyModes"))
    method = validateSawClassMethod(method, false, "CalcSingleHS:identifyModes")
    thr === nothing && (thr = defaultThreshold(method, false))
    if method === :ritz
        (isempty(cs.ritzScoreVelAxis) || isnan(cs.longLimVel)) && throw(ArgumentError(
            "CalcSingleHS:identifyModes:noRitzScore -- Ritz classification needs calcLim() and " *
            "calcRitzScore() first (calcAll runs both)."))
        cs.classification = calc(ModeClassification(threshold = thr), cs.ritzSolver, cs.structure,
                                 nothing, nothing, (cands = cs.ritzScoreCands,))
    else
        cr = (longLimVel = cs.longLimVel, longLimPol = cs.longLimPol,
              shearLimVels = cs.shearLimVels, shearLimPols = cs.shearLimPols)
        if method === :green
            cs.greenScoreRef === nothing && throw(ArgumentError("CalcSingleHS:identifyModes:noScore -- " *
                "calcGreenScore() must run first (calcAll(:full) runs it)."))
            gr = (greenScore = (peaks = cs.greenScorePeaks,),)
        else
            cs.refinedSaw === nothing && throw(ArgumentError("CalcSingleHS:identifyModes:noRefined -- " *
                "calcRefined() must run first (classify(:refined) runs it)."))
            gr = (refined = (peaks = cs.refinedSaw.peaks,),)
        end
        cs.classification, cs.limitingEntries = classifyGreenPath(cs, thr, cr, gr)
        cs.limitingDone = true
        if method === :refined
            acc = GreenSpotScore.accepted(cs.refinedSaw.peaks, thr)
            r = cs.refinedSaw
            cs.refinedSaw = RefinedSaw(thr, r.ref, GreenSpotScore.Peak[_with(r.peaks[j]; accepted = acc[j])
                                                                      for j in eachindex(r.peaks)], r.info)
        end
    end
    cs.activeMethod = method
    cs.activeThreshold = thr
    cs.classificationVersion = CLASSIFICATION_VERSION
    return cs
end

"""
    calcAll(cs, mode = nothing, retentionLevel = nothing, precomp = nothing) -> cs

`:full`: calcRitz -> calcLim -> calcGreen -> calcGreenScore -> calcRitzScore ->
identifyModes(:green) (then the chosen method, if another) -> prune. `:ritzOnly`: calcRitz ->
calcLim -> calcRitzScore -> identifyModes(:ritz) -> prune. An unset `sawClassMethod` runs
the mode's own method. Pruning always comes after the classification.
"""
function calcAll(cs::CalcSingleHS, mode = nothing, retentionLevel = nothing, precomp = nothing)
    mode === nothing && (mode = cs.calcMode)
    retentionLevel === nothing && (retentionLevel = cs.retentionLevel)
    validateCalcMode(mode, "CalcSingleHS:calcAll")
    validateRetentionLevel(retentionLevel, "CalcSingleHS:calcAll", (:full, :saw_only, :none))
    target = resolveSawClassMethod(cs.sawClassMethod, cs.sawThreshold, mode, false, "CalcSingleHS:calcAll")
    targetThr = cs.sawThreshold
    cs.calcMode = mode
    calcRitz(cs, precomp)
    cs.classificationVersion = 0
    calcLim(cs)
    if mode === :full
        calcGreen(cs)
        calcGreenScore(cs, cs.greenScoreRef)
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

"""
    applyRetention(cs, retentionLevel) -> cs

Eigenvector retention (`pruneRitzTo`): `:full` | `:saw_only` | `:none`.
"""
function applyRetention(cs::CalcSingleHS, retentionLevel)
    cs.ritzSolver = pruneRitzTo(cs.ritzSolver, cs.classification, retentionLevel, "CalcSingleHS:applyRetention")
    return cs
end

"""
    setRunMode(cs, mode) -> cs

Stamp the calcMode of a sweep's run (its 'full' phases call the steps, not calcAll).
"""
function setRunMode(cs::CalcSingleHS, mode)
    validateCalcMode(mode, "CalcSingleHS:setRunMode")
    cs.calcMode = mode
    return cs
end

"""
    calcSpotGreen(cs, kernel::ProbeKernel; DirStep, RelKStep, Eta, VelStep, useParallel = :auto)
        -> cs

EXPERIMENTAL (as in MATLAB; tried on Fe3Al(316) only, not validated against measurements).
|G| as a finite-pattern measurement sees it: the coherent sum of the complex Green's function
over the wavevectors a [`ProbeKernel`](@ref) couples to (`SpotGreen`), for this angle. Stores
`spotGreen` (`[nVel]` `Float32`, |X|), `spotGreenVelAxis` [km/s] and `spotGreenInfo` (kernel, steps, eta,
`vEnd`, `nDirections`, `nNodes`, per-angle quadrature checks `dirConvergence` /
`relKConvergence`). A picture only; no classification reads it. Options: `DirStep` [deg]
(default `kernel.dirStep`), `RelKStep` (default `min(kernel.relKStep, Eta / 4)`), `Eta`
(Green damping, default `SawSpectrum.ETA`), `VelStep` [km/s] (solve/output axis step, default
0.004). One Green solve per kernel direction; the |k| spread costs none. Equal, bit for bit,
to this angle's column of a sweep's [`calcSpotGreenMap`](@ref) with the same kernel and
options. Any calcMode; half-space only (a layered structure throws an `ArgumentError`
`CalcSingleHS:calcSpotGreen:layered`). Plot with `plotGreenCurve(cs; component = "spot")`.
"""
function calcSpotGreen(cs::CalcSingleHS, kernel::ProbeKernel; useParallel = :auto, kw...)
    nlay(cs.structure) == 1 || throw(ArgumentError("CalcSingleHS:calcSpotGreen:layered -- " *
        "calcSpotGreen is for a half-space (nLay = 1) only."))
    opts = SpotGreen.options(kernel; kw...)
    absX, vOut, info = SpotGreen.compute(cs.structureBase, [cs.angle], kernel, opts; useParallel)
    cs.spotGreen = vec(absX)
    cs.spotGreenVelAxis = vOut
    cs.spotGreenInfo = info
    return cs
end
