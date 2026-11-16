## Tutorial 09 (code) -- inversion: recover Cu constants from synthetic SAW/pSAW/limiting velocities
# Narrated version: juliaPort/tutorialsJ/tutorial_09_inversion_synthetic.jl.
# Run: at the prompt of `julia --project=. -t auto` started in your own environment with GuliwerrJ and GLMakie
#   (juliaPort/docsJ/gettingStartedJ.md, section 2): include("/path/to/GuliweRR/juliaPort/tutorialsCodeJ/tutorial_09_inversion_synthetic.jl")

using GuliwerrJ
using GLMakie
using Printf
using Random

## Ground truth, deliberately wrong initial guess, and a Cu {001} plane
dirMat        = [1 0 0; 0 1 0]     # normal [100], reference [010] (cubic: same as the default)
matTrue, _    = MaterialFromLibrary("Cu"; dirMat = dirMat)   # "true" Cu from the materials library
cuConstsTrue  = matTrue.elConsts.vector    # [C11, C12, C44] GPa
density       = matTrue.rho                # g/cm^3
cuConstsGuess = cuConstsTrue .* [1.10, 0.90, 1.08]   # +10%, -10%, +8%
symmNo        = 2                  # cubic

## Discretisation -- "easy" tier (Zener ratio ~3.2, principal (001) plane)
relDepth   = 20.0
legOrder   = 50
calcParams = CalcParams(; relDepth = relDepth, legOrder = legOrder)
anglesDeg  = collect(0.0:2.0:180.0)
nMaxModes  = 500

## Step 1 -- forward sweep of the true material
calcMode       = :full             # Green needed for pseudoRayleigh points
retentionLevel = :full             # association() needs every eigenvector
strTrue   = Structure(; materIn = matTrue, calcParamsIn = calcParams, name = "cuTrue",
                        nMaxModes = nMaxModes)
sweepTrue = calc(angSweepHS(; structure = strTrue, angles = anglesDeg, calcMode = calcMode,
                              retentionLevel = retentionLevel))

## Step 2a -- harvest Rayleigh and the most prominent pseudoRayleigh per angle
# modeVel: the bound Ritz mode's eigenvelocity, what the inversion fits (vel = the Green peak)
rayAngle = Float64[];    rayVel = Float64[]
pseudoAngle = Float64[]; pseudoVel = Float64[]
for iAngle in eachindex(anglesDeg)
    sigModes = sweepTrue.significantModesMap[iAngle]
    idxRay = findfirst(m -> m.type == :Rayleigh, sigModes)
    if idxRay !== nothing
        push!(rayAngle, sweepTrue.angles[iAngle]); push!(rayVel, sigModes[idxRay].modeVel)
    end
    idxPseudo = findall(m -> m.type == :pseudoRayleigh, sigModes)
    if !isempty(idxPseudo)
        iBest = idxPseudo[argmax([sigModes[i].relPeakHeight for i in idxPseudo])]   # tallest Green peak
        push!(pseudoAngle, sweepTrue.angles[iAngle]); push!(pseudoVel, sigModes[iBest].modeVel)
    end
end

## Step 2b -- Cu(001) Rayleigh-observability gap, longitudinal points, dropout, noise
gapAngleLow  = 36.0                # deg -- no Rayleigh velocity kept here (pseudoRayleigh range)
gapAngleHigh = 54.0                # deg
keepRay  = .!((rayAngle .>= gapAngleLow) .& (rayAngle .<= gapAngleHigh))
rayAngle = rayAngle[keepRay]
rayVel   = rayVel[keepRay]

longAngleStep = 10.0               # deg
longAngle     = anglesDeg[mod.(anglesDeg, longAngleStep) .== 0]
longIdx       = Int.(indexin(longAngle, sweepTrue.angles))
longVel       = sweepTrue.longVelMap[longIdx]

allAngle = vcat(rayAngle, pseudoAngle, longAngle)
allVel   = vcat(rayVel, pseudoVel, longVel)
allType  = vcat(fill("rayleigh", length(rayAngle)), fill("pseudoRayleigh", length(pseudoAngle)),
                fill("longitudinal", length(longAngle)))

randomSeed      = 20260911         # Julia's RNG differs from MATLAB's: same seed, different draws
dropoutFraction = 0.15
noiseSigma      = 0.005            # km/s, ~ real TGS RMSE
Random.seed!(randomSeed)
nAll     = length(allAngle)
keepMask = trues(nAll)
keepMask[randperm(nAll)[1:round(Int, dropoutFraction * nAll)]] .= false
allAngle = allAngle[keepMask]
allVel   = allVel[keepMask] .+ noiseSigma .* randn(count(keepMask))
allType  = allType[keepMask]

measuredData = measured(allAngle, allVel, allType, noiseSigma)   # waveType list: WAVE_TYPES
@printf("%d measured points\n", measuredData.nPoints)

## Step 3 -- forward sweep at the guess, association built against IT
ecGuess    = ElastConsts(cuConstsGuess, symmNo)
matGuess   = Material(; elConsts = ecGuess, rho = density, dirMat = dirMat)
strGuess   = Structure(; materIn = matGuess, calcParamsIn = calcParams, name = "cuGuess",
                         nMaxModes = nMaxModes)
sweepGuess = calc(angSweepHS(; structure = strGuess, angles = anglesDeg, calcMode = calcMode,
                               retentionLevel = retentionLevel))

velocityTolerance = 0.5            # km/s, widened from the 0.3 default for a far guess
angleTolerance    = 1e-6           # deg
assocGuess = association(sweepGuess, measuredData, velocityTolerance, angleTolerance)

## Step 4 -- sensitivity: which constant combinations the data determine
invertedLayerNo = 1
sens = sensitivity(strGuess, assocGuess, invertedLayerNo, measuredData)
@printf("%d of %d combinations well-determined\n", sens.nOptimisedConstants, sens.nIndepConsts)
axsSens = plotSensitivitySpectrum(sens)   # two Axis in one figure
display(axsSens[1].parent)                 # plot functions only build figures; display shows them

## Step 5 -- inversion (keywords: stepwise, anchorPolicy, includePseudoSAW, pseudoSAWRefinement)
stepwise  = true                   # true: grow the optimised set one at a time | false: all at once
inversion = elasticity(strGuess, measuredData, assocGuess, sens, invertedLayerNo; stepwise = stepwise)
inversion = calc!(inversion)

constNames = ("c11", "c12", "c44")
for iConst in eachindex(constNames)
    @printf("%s: true %.2f, guess %.2f, recovered %.2f +/- %.2f GPa\n", constNames[iConst],
            cuConstsTrue[iConst], cuConstsGuess[iConst],
            inversion.result.voigtElasticVec[iConst], inversion.result.voigtElasticVecErrors[iConst])
end
rmse = sqrt(inversion.result.costFunctionValue / inversion.result.nDataPoints)
@printf("RMSE %.4f km/s (noise %.4f km/s)\n", rmse, noiseSigma)

## Step 5b -- second pass from the recovered constants (fixes points matched to the wrong branch from the far guess)
constsPass1 = copy(inversion.result.voigtElasticVec)
matPass2    = Material(; elConsts = ElastConsts(constsPass1, symmNo), rho = density, dirMat = dirMat)
strPass2    = Structure(; materIn = matPass2, calcParamsIn = calcParams, name = "cuPass2",
                          nMaxModes = nMaxModes)
sweepPass2  = calc(angSweepHS(; structure = strPass2, angles = anglesDeg, calcMode = calcMode,
                                retentionLevel = retentionLevel))
assocPass2 = association(sweepPass2, measuredData, velocityTolerance, angleTolerance)
sensPass2  = sensitivity(strPass2, assocPass2, invertedLayerNo, measuredData)
inversion  = elasticity(strPass2, measuredData, assocPass2, sensPass2, invertedLayerNo; stepwise = stepwise)
inversion  = calc!(inversion)
for iConst in eachindex(constNames)
    @printf("%s: true %.2f, pass 1 %.2f, pass 2 %.2f +/- %.2f GPa\n", constNames[iConst],
            cuConstsTrue[iConst], constsPass1[iConst],
            inversion.result.voigtElasticVec[iConst], inversion.result.voigtElasticVecErrors[iConst])
end
rmse = sqrt(inversion.result.costFunctionValue / inversion.result.nDataPoints)
@printf("RMSE %.4f km/s (noise %.4f km/s)\n", rmse, noiseSigma)

## Step 6 -- fit quality and mode-tracking diagnostic
assocFinal = inversion.association[1]   # one element per crystal cut; a single cut here
axFit = plotInversionFit(measuredData, assocFinal, assocGuess)   # assocGuess = initial-guess snapshot
display(axFit.parent)
axAssoc = plotAssociationDiagnostic(assocFinal)
display(axAssoc.parent)

## Step 7 -- Green maps: true vs. recovered, with the association overlaid
sweepRecovered = recalcForwardAtResult(inversion; angles = anglesDeg, calcMode = calcMode,
                                       retentionLevel = retentionLevel)
overlayMode  = true                # markers only, keep the Green colormap
figCompare   = Figure(; size = (1100, 850))
sweepsToShow = (sweepTrue, sweepRecovered)
panelTitles  = ("True material", "Recovered material")
for iPanel in 1:2
    axTop, _ = plotGreenMap(sweepsToShow[iPanel]; ax = Axis(figCompare[1, iPanel]))
    axTop.title = panelTitles[iPanel]
    axBottom, _ = plotGreenMap(sweepsToShow[iPanel]; ax = Axis(figCompare[2, iPanel]))
    plotAssociationDiagnostic(assocFinal, overlayMode; ax = axBottom)
    axBottom.title = panelTitles[iPanel] * " + association"
end
display(figCompare)
