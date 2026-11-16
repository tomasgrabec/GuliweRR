## Tutorial: elastic-constant inversion from synthetic SAW/limiting velocities (Julia port)
#
# Same physics as examples/tutorials/tutorial_09_inversion_synthetic.m: the INVERSE problem --
# given measured surface-wave velocities (here, synthesised from a known material so the answer
# is known in advance), recover the elastic constants that produced them. Copper, cubic, on a
# {001} plane (surface normal [100], propagation reference [010]; for a cubic crystal equivalent
# to Material's default dirMat). The "measurement" is Cu's own Rayleigh, pseudo-Rayleigh and
# longitudinal limiting velocities with noise added; with real TGS (transient grating
# spectroscopy) data, built with `measured` or `fromTable` (MATLAB: MeasuredVelocities.fromTable),
# everything downstream is the same.
#
# See tutorial_01_isotropic_halfspace.jl's header for how to run a script in this series.
#
# Julia concepts introduced in this tutorial:
#   - `push!`, Julia's replacement for MATLAB's `array(end+1) = x` growing idiom -- the FIRST
#     genuinely mutating function this series has you call directly (as opposed to the
#     package's own internal use of the pattern); the trailing `!` is a NAMING CONVENTION
#     (not special syntax) meaning "this function modifies its argument in place"
#   - `ElasticityInversion`, a mutable struct like all `Calc*` results, and its `calc!`
#     (the `!` marks that it mutates its argument in place and also returns it)
#   - a real naming collision to avoid: several GuliwerrJ functions share their MATLAB
#     `Inv.xxx` name exactly (`measured`, `association`, `sensitivity`) with no namespace
#     prefix to disambiguate -- naming the RESULT of such a call the same as the function
#     itself would silently shadow that function for the rest of the session
#   - `Random.seed!`/`randperm`/`randn` -- Julia's RNG is a DIFFERENT algorithm from MATLAB's,
#     so "the same seed" does not reproduce the same numbers across languages (only within one)

using GuliwerrJ
using GLMakie   # opens an interactive window; `using CairoMakie` instead renders the same
                 # Figure/Axis objects headlessly (e.g. for save("fig.png", fig)) -- see
                 # tutorial_02_angular_map.jl for the full explanation
using Printf
using Random

## Ground truth and initial (deliberately wrong) guess
dirMat = [1 0 0;    # surface normal [100]
          0 1 0]    # propagation reference [010]
# A matrix literal with MORE THAN ONE ROW is written with SPACES between elements and `;`
# between rows -- `[1,3,7; -3,1,0]` (commas instead of spaces) is a Julia SYNTAX ERROR, not
# merely unconventional: commas already mean "build a plain Vector" on their own, and mixing
# that reading with a row-separating `;` in the same brackets does not parse. Space-separated
# elements is the one spelling that means "matrix row" in Julia, exactly as in MATLAB -- so the
# safest habit coming from MATLAB is to keep using spaces (never commas) inside a multi-row
# literal, which also happens to be identical to what MATLAB itself would accept.

# The "true" Cu comes from the materials library; the guess is the truth scaled by
# +10%, -10%, +8% -- a genuinely wrong starting point (starting from the truth would
# make this exercise vacuous).
matTrue, _ = MaterialFromLibrary("Cu"; dirMat = dirMat)
cuVecTrue  = matTrue.elConsts.vector              # [c11 c12 c44] GPa
cuRho      = matTrue.rho                          # g/cm^3
cuVecGuess = cuVecTrue .* [1.10, 0.90, 1.08]      # +10%, -10%, +8%
symmNo     = 2                                    # cubic

relDepth = 20   # Legendre-domain truncation depth [wavelengths]
legOrder = 50   # number of Legendre basis functions spanning the depth -- must stay
                 # >= 2x relDepth (docs/pitfalls.md MTH1, "easy" tier: Cu's Zener ratio
                 # 2*c44/(c11-c12) ~ 3.2 at this principal plane)
calcParams = CalcParams(; relDepth = relDepth, legOrder = legOrder)
angleGrid  = collect(0.0:2.0:180.0)

## Step 1 -- forward-calculate the "true" material, to synthesise a measurement from
strTrue = Structure(; materIn = matTrue,
                       calcParamsIn = calcParams, name = "cuTrue", nMaxModes = 500)
calcMode       = :full   # Green is needed: pseudoRayleigh points come from Green peaks
retentionLevel = :full   # association() needs every Ritz eigenvector
sweepTrue = angSweepHS(; structure = strTrue, angles = angleGrid, calcMode = calcMode,
                          retentionLevel = retentionLevel)
sweepTrue = calc(sweepTrue)

## Step 2 -- synthesise the "measured" dataset: Rayleigh AND pseudo-Rayleigh
# velocities wherever each is classified, longitudinal limiting velocities
# every 10 deg, a Cu(001) observability gap (no Rayleigh velocity kept at 36-54 deg, an illustrative range
# where the pseudo-Rayleigh wave takes over in the calculation), a further random 15% dropout, and
# 0.005 km/s Gaussian noise (an assumed value for this example).
#
# The pseudo-Rayleigh points matter: over the 36-54 deg gap (and its mirror, 126-144 deg) they
# are the only surface-wave data, which is why includePseudoSAW defaults to true in
# sensitivity/elasticity. Where one branch ends and the other starts is a threshold on a
# continuum (docs/pitfalls.md PHY1).
# modeVel: the Ritz eigenvelocity of each entry's mode, the quantity the inversion fits (an
# entry's vel is its Green peak velocity since CLASSIFICATION_VERSION 3).
rayAngle = Float64[]; rayVel = Float64[]
pseudoAngle = Float64[]; pseudoVel = Float64[]
for ia in eachindex(angleGrid)
    sm = sweepTrue.significantModesMap[ia]

    idxR = findfirst(m -> m.type == :Rayleigh, sm)
    if idxR !== nothing
        push!(rayAngle, sweepTrue.angles[ia])
        push!(rayVel, sm[idxR].modeVel)
    end
    # `push!(vec, x)` appends `x` to the end of `vec`, MUTATING `vec` in place -- the Julia
    # replacement for MATLAB's `rayAngle(end+1,1) = ...` growing-by-assignment idiom. The
    # trailing `!` is how Julia spells "mutates its argument" by CONVENTION (not a special
    # operator): any function ending in `!` should be assumed to modify at least one of its
    # arguments, exactly like `push!` here.

    # Several pseudoRayleigh entries can coexist at one angle (a genuine forward-model feature,
    # not noise) -- take the most prominent one, the same choice a real experimenter would make
    # when reading an angular map by eye.
    idxP = findall(m -> m.type == :pseudoRayleigh, sm)
    if !isempty(idxP)
        iBest = idxP[argmax([sm[i].relPeakHeight for i in idxP])]
        push!(pseudoAngle, sweepTrue.angles[ia])
        push!(pseudoVel, sm[iBest].modeVel)
    end
end

gapAngleLow  = 36.0   # Cu(001) Rayleigh-observability gap [deg]: no Rayleigh
gapAngleHigh = 54.0   # velocity is kept here
keepMask = .!((rayAngle .>= gapAngleLow) .& (rayAngle .<= gapAngleHigh))
rayAngle = rayAngle[keepMask]
rayVel   = rayVel[keepMask]

longitudinalAngleStep = 10.0   # [deg] spacing of the synthetic longitudinal-limiting points
longAngle = angleGrid[mod.(angleGrid, longitudinalAngleStep) .== 0]
longIdx = Int.(indexin(longAngle, sweepTrue.angles))
longVel = sweepTrue.longVelMap[longIdx]

allAngle = vcat(rayAngle, pseudoAngle, longAngle)
allVel   = vcat(rayVel, pseudoVel, longVel)
allType  = vcat(fill("rayleigh", length(rayAngle)),
                 fill("pseudoRayleigh", length(pseudoAngle)),
                 fill("longitudinal", length(longAngle)))
nCandidateBeforeDropout = length(allAngle)   # after the Cu(001) gap, before the random dropout

randomSeed            = 20260911   # fixed for a reproducible tutorial run -- but see this
                                     # file's header note on Julia vs. MATLAB random streams
randomDropoutFraction = 0.15        # fraction of the (already gap-reduced) points dropped
                                      # at random, mimicking further gaps in the data
noiseSigma = 0.005   # [km/s] additive Gaussian noise (assumed value) -- reused below for measuredData's velocityStd and
                      # for judging the fit's residual against the noise floor

Random.seed!(randomSeed)
# Fixing the seed makes THIS script reproducible from run to run, exactly as MATLAB's own
# `rng(randomSeed, 'twister')` does for the MATLAB tutorial -- but Julia's default random-number
# generator is a DIFFERENT algorithm from MATLAB's Mersenne Twister, so the same seed value does
# NOT give the same sequence of "random" numbers across the two languages. The dropped points
# and added noise will differ between this script and tutorial_09_inversion_synthetic.m even
# though both are individually reproducible -- expected, and irrelevant to what the tutorial
# demonstrates (recovering the true elastic constants within the noise budget).
nAll = length(allAngle)
nDrop = round(Int, randomDropoutFraction * nAll)
dropIdx = randperm(nAll)[1:nDrop]   # a uniformly random nDrop-subset of 1:nAll, Julia's
                                      # analogue of MATLAB's `randperm(nAll, nDrop)`
keepMask2 = trues(nAll)
keepMask2[dropIdx] .= false
allAngle = allAngle[keepMask2]
allVel   = allVel[keepMask2] .+ noiseSigma .* randn(count(keepMask2))
allType  = allType[keepMask2]

measuredData = measured(allAngle, allVel, allType, noiseSigma)
# Named `measuredData`, NOT `measured` (unlike the MATLAB tutorial's own
# `measured = Inv.measured(...)`): MATLAB's `Inv.` namespace prefix lets the function and the
# resulting variable share the exact same name with no ambiguity. GuliwerrJ's `measured` is a
# plain top-level function with no namespace to qualify it (see this file's header) -- reusing
# that name for a local variable would work for this ONE assignment (Julia evaluates the
# right-hand side, `measured(...)`, before binding the left-hand name), but would then leave
# `measured` referring to the DATA, not the function, for the rest of the session; the same
# renaming is applied to `association`/`sensitivity` below for the identical reason.
@printf("Synthesised %d measured points (%d dropped at random from the %d that survived the Cu(001) gap):\n",
        measuredData.nPoints, nCandidateBeforeDropout - measuredData.nPoints, nCandidateBeforeDropout)
@printf("  rayleigh %d, pseudoRayleigh %d, longitudinal %d\n",
        count(==("rayleigh"), measuredData.waveType), count(==("pseudoRayleigh"), measuredData.waveType),
        count(==("longitudinal"), measuredData.waveType))

## Step 3 -- forward-calculate at the (wrong) initial guess, and build the
# initial mode association against IT, not against the truth -- associating
# against the truth would smuggle the answer into the input.
strGuess = Structure(; materIn = Material(; elConsts = ElastConsts(cuVecGuess, symmNo),
                                             rho = cuRho, dirMat = dirMat),
                        calcParamsIn = calcParams, name = "cuGuess", nMaxModes = 500)
sweepGuess = angSweepHS(; structure = strGuess, angles = angleGrid, calcMode = calcMode,
                           retentionLevel = retentionLevel)
sweepGuess = calc(sweepGuess)

velocityTolerance = 0.5    # [km/s] max |v_calc - v_meas| to accept a candidate mode match.
                             # Widened from the association default (0.3 km/s): the guess is
                             # far enough from the truth that 0.3 would drop most Rayleigh
                             # points
angleTolerance = 1e-6      # [deg] max mismatch to treat a measured angle as "the same" sweep
                             # angle -- tight, since both are drawn from the same angleGrid here
assocGuess = association(sweepGuess, measuredData, velocityTolerance, angleTolerance)

## Step 4 -- sensitivity analysis: which combinations of elastic constants
# can this dataset actually determine?
invertedLayerNo = 1   # which layer's constants to recover (1 = the half-space itself)
sens = sensitivity(strGuess, assocGuess, invertedLayerNo, measuredData)
@printf("Sensitivity: %d of %d transformed combinations judged well-determined ",
        sens.nOptimisedConstants, sens.nIndepConsts)
@printf("(cutoff = 1/%.0f error amplification).\n", 1 / sens.sensitivityCutoff)

axsSens = plotSensitivitySpectrum(sens)   # a tuple of two Axis (spectrum, transform matrix)
display(axsSens[1].parent)

## Step 5 -- run the inversion
# Keyword options (all optional, lower-case first letter in this port): stepwise (true: grow
# the optimised set one combination at a time | false: all at once), anchorPolicy (:chained |
# :fixed), includePseudoSAW (true), pseudoSAWRefinement (:green | :none).
inversion = elasticity(strGuess, measuredData, assocGuess, sens, invertedLayerNo)
inversion = calc!(inversion)
# `calc!(inversion)`, not `inversion.calc()` (MATLAB) and not `calc(inversion)` (every other
# `calc` so far in this series). The `!` marks a function that changes its argument: `calc!`
# edits `inversion`'s own fields as the optimiser iterates, and also returns it (so `inversion = calc!(inversion)` reads the same
# as MATLAB's `inversion = inversion.calc()`, even though the reassignment is not strictly
# needed here -- `inversion` was already changed in place).

@printf("\n%8s %12s %12s %12s\n", "", "true", "guess", "recovered +/- error")
labels = ("c11", "c12", "c44")
for i in 1:3
    @printf("%8s %12.3f %12.3f %8.3f +/- %.3f GPa\n", labels[i], cuVecTrue[i], cuVecGuess[i],
            inversion.result.voigtElasticVec[i], inversion.result.voigtElasticVecErrors[i])
end
@printf("\nFinal cost function value: %.4g (RMSE %.4f km/s, noise floor %.4f km/s)\n",
        inversion.result.costFunctionValue,
        sqrt(inversion.result.costFunctionValue / inversion.result.nDataPoints), noiseSigma)

## Step 5b -- second pass, restarted from the recovered constants
# From a far guess, Step 3 matches some measured points to the wrong branch, and the fit then
# stops short of the true minimum (RMSE well above the noise). Re-running the forward sweep and
# the association at the first result, then inverting again, corrects those matches. Do the
# same with real data: repeat until the constants stop changing and the RMSE reaches the noise
# level.
vecPass1 = copy(inversion.result.voigtElasticVec)
strPass2 = Structure(; materIn = Material(; elConsts = ElastConsts(vecPass1, symmNo),
                                             rho = cuRho, dirMat = dirMat),
                        calcParamsIn = calcParams, name = "cuPass2", nMaxModes = 500)
sweepPass2 = calc(angSweepHS(; structure = strPass2, angles = angleGrid, calcMode = calcMode,
                               retentionLevel = retentionLevel))
assocPass2 = association(sweepPass2, measuredData, velocityTolerance, angleTolerance)
sensPass2  = sensitivity(strPass2, assocPass2, invertedLayerNo, measuredData)
inversion  = elasticity(strPass2, measuredData, assocPass2, sensPass2, invertedLayerNo)
inversion  = calc!(inversion)

@printf("\nSecond pass:\n%8s %12s %12s %12s\n", "", "true", "pass 1", "recovered +/- error")
for i in 1:3
    @printf("%8s %12.3f %12.3f %8.3f +/- %.3f GPa\n", labels[i], cuVecTrue[i], vecPass1[i],
            inversion.result.voigtElasticVec[i], inversion.result.voigtElasticVecErrors[i])
end
@printf("\nFinal cost function value: %.4g (RMSE %.4f km/s, noise floor %.4f km/s)\n",
        inversion.result.costFunctionValue,
        sqrt(inversion.result.costFunctionValue / inversion.result.nDataPoints), noiseSigma)

# Mode-tracking health, per binding class -- worth checking explicitly whenever pseudoRayleigh
# data is included, since a mis-tracked pseudo-SAW point can read as perfectly confident on its
# Ritz-side MAC alone even when locked onto the wrong mode. ElasticityInversion's Green's-
# function-based cross-check is what greenAgreement* below reports on.
trk = inversion.result.tracking
@printf("\nMode-tracking health (%s pseudo-SAW refinement):\n", String(trk.pseudoSAWRefinement))
@printf("  %-12s %8s %14s %10s\n", "class", "nPoints", "minConfidence", "nBelowThr")
for cls in (:ritzMode, :pseudoSAW, :limiting)
    trk.nPointsPerClass[cls] == 0 && continue
    @printf("  %-12s %8d %14.4f %10d\n", String(cls), trk.nPointsPerClass[cls],
            trk.minConfidence[cls], trk.nBelowThresholdPerClass[cls])
end
# `NamedTuple` fields index by Symbol with square brackets (`trk.nPointsPerClass[cls]`) exactly
# as MATLAB's dynamic field access `trk.nPointsPerClass.(c)` does with a string -- one of the
# few places a MATLAB dynamic-field idiom translates almost verbatim.
@printf("  collisions=%d, limitingRootsLost=%d, excludedAtFinalEval=%d\n",
        trk.nModeCollisions, trk.nLimitingRootsLost, trk.nExcludedAtFinalEval)
if trk.pseudoSAWRefinement === :green
    @printf("  Green rescans=%d, pseudo-SAW/Green agreement: max %.4f, median %.4f\n",
            trk.nGreenRescans, trk.greenAgreementMax, trk.greenAgreementMedian)
end

## Step 6 -- results analysis: fit quality and mode-tracking diagnostic
# One plot, every stage: the full measured dataset (dropped points marked with a grey x --
# association is a value type, so the ORIGINAL `assocGuess` variable from Step 3 is untouched by
# calc!(inversion) and still holds the initial-guess fit; only inversion.association is the
# final, converged one).
axFit = plotInversionFit(measuredData, inversion.association[1], assocGuess)
axFit.title = "Measured, initial guess, and recovered fit -- all in one"
display(axFit.parent)   # a plot function only builds the figure; `display` shows it (tutorial 2)

axAssoc = plotAssociationDiagnostic(inversion.association[1])
display(axAssoc.parent)
# `inversion.association` is a `Vector{ModeAssociation}` (one entry per jointly-inverted
# dataset, see ElasticityInversion.jl -- this tutorial uses exactly one dataset, so `[1]` is the
# only entry that exists), unlike MATLAB's own `inversion.association`, which for a
# single-dataset inversion IS the one `ModeAssociation` object directly.

## Step 7 -- Green's-function maps: true vs. recovered, with tracking overlaid
# Top row: the |G13_eff| map of the true and of the recovered material -- nearly identical if
# the inversion converged. Bottom row: the same maps with every associated point overlaid
# (overlay mode: markers only, the map keeps its colormap), coloured by tracking confidence;
# grey diamonds are limiting (Christoffel) points. Comparing the two bottom panels shows
# whether every point, pseudo-SAW included, sits on the same ridge in both.
#
# recalcForwardAtResult runs a fresh forward sweep at the inverted elasticity -- deliberately
# NOT done inside the optimisation loop itself, so it costs one extra sweep here. (vMax, which
# a layered stack may pass to extend the Green axis, is optional and not needed here.)
sweepRecovered = recalcForwardAtResult(inversion; angles = angleGrid, calcMode = calcMode,
                                                    retentionLevel = retentionLevel)

figCompare = Figure(; size = (1100, 850))

axTrueTop, _ = plotGreenMap(sweepTrue; ax = Axis(figCompare[1, 1]))
axTrueTop.title = "True material"

axRecTop, _ = plotGreenMap(sweepRecovered; ax = Axis(figCompare[1, 2]))
axRecTop.title = "Recovered material"

axTrueBottom, _ = plotGreenMap(sweepTrue; ax = Axis(figCompare[2, 1]))
plotAssociationDiagnostic(inversion.association[1], true; ax = axTrueBottom)
axTrueBottom.title = "True material + association overlay"

axRecBottom, _ = plotGreenMap(sweepRecovered; ax = Axis(figCompare[2, 2]))
plotAssociationDiagnostic(inversion.association[1], true; ax = axRecBottom)
axRecBottom.title = "Recovered material + association overlay"

Label(figCompare[0, :], "Green's-function maps: true vs. recovered, with mode tracking overlaid";
      fontsize = 16)
display(figCompare)
# `plotGreenMap(sweepTrue; ax = Axis(figCompare[1, 1]))` -- passing an explicit `ax` keyword
# into a caller-owned Figure cell replaces MATLAB's `subplot(2,2,i)` + passing that axes handle
# as plotGreenMap's own 5th positional argument; the panel LAYOUT itself is expressed the same
# way tutorial 3 used for its own multi-panel Label (`figCompare[row, col]`).

## What to try next
# - Loosen/tighten cuVecGuess to see how far a starting guess can be from the truth before mode
#   association or convergence starts failing.
# - Try a general (non-principal) crystallographic plane: build dirMat from an arbitrary surface
#   normal and CalcParams(; relDepth = 30, legOrder = 70) instead. The recovered crystal-frame
#   constants should agree with this script's, independent of which plane the "measurement"
#   was cut on.
# - Add addPhaseVelocity(...)/setBulkModulus(...) to measuredData for the two auxiliary data
#   channels, and re-run from Step 3.
