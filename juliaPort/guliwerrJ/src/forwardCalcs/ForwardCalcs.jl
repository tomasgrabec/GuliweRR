# =========================================================================================
# ForwardCalcs.jl -- the forward-calculation module (`GuliwerrJ.ForwardCalcs`), one module for
# the whole `forwardCalcs/` folder; the files below are `include`d as plain top-level
# definitions (the score namespaces are submodules).
#
# Relative imports: `ForwardCalcs` sits one level below `GuliwerrJ`, so its siblings are
# reached with two dots (`using ..CoreTypes`, `using ..Solvers`).
#
# Shared `calc` generic: the forward types add methods to the ONE `calc` generic defined in
# `Solvers.Ritz`, which needs `import ..Solvers: calc` (a plain `using` would not allow the
# extension, and a fresh `calc` here would collide with Solvers' in GuliwerrJ).
#
# INCLUDE ORDER: CalcCommon.jl first (the score modules import its helpers),
# then the score submodules, then ModeClassification.jl (its struct types name the score
# modules' Peak types), then the Calc types.
# =========================================================================================

module ForwardCalcs

using LinearAlgebra: norm, dot, qr
using Statistics: median
using Interpolations: interpolate, Gridded, Constant

using ..CoreTypes: Structure, Material, nlay, ModifyLayer, SetDirMat, GetRotDirMat,
                    GetStandardDirMat, GetRotElastTensor,
                    movmedian, findpeaks_matlab,
                    should_use_parallel, threads_available, with_blas_threads, with_task_blas_single,
                    ElastConsts, CalcParams, AddLayer, IndepVoigtPositions, ApplyWavelength,
                    HillAverages, ProbeKernel

import ..Solvers: calc
using ..Solvers: RitzSolverSAW, RitzSolverLamb, calcDisplacement,
                  precomp_field,
                  GreenSolver, matrix, findRoots, green_component, green_gfull,
                  buildConstraintMatrix, SAWBoundary, LambBoundary, derivsByDepth,
                  pruneEigVecs, ritzStructurePrecompute,
                  ritzGramBlocks, gaugeBackMap, RawCellFor, kPolynomial, ZgvNewton
import ..Solvers: RitzWorkspace          # + a constructor from a RitzWavelengthSweep (CalcCommon.jl)
using ..Solvers.Ritz: _contract_depth_rho!, _structure_from_gram!, _check_workspace
using ..Solvers.Christoffel: DEGENERATE_REL_TOL

include("CalcCommon.jl")
# --- scores: the SAW-classification building blocks, one submodule each -------------
include("CalcCommonScores.jl")
include("SawSpectrum.jl")          # module SawSpectrum + module GreenPeaks
include("GreenScore.jl")
include("GreenScoreLayered.jl")
include("RitzScore.jl")
include("GreenSpotScore.jl")
include("SpotGreen.jl")
include("DispersionBranches.jl")
include("ModeClassification.jl")
include("CalcSingleHS.jl")
include("CalcSingleLayered.jl")
include("CalcSinglePlate.jl")
include("CalcIsofreqSingleLayered.jl")
include("CalcAngSweepHS.jl")
include("CalcAngSweepLayered.jl")
include("CalcAngSweepPlate.jl")
include("CalcIsofreqAngSweepLayered.jl")
include("CalcWavevecAngSweep2Layered.jl")
include("CalcWavevecAngSweep2Plate.jl")
include("CalcFreqAngSweep2Layered.jl")
include("Calc.jl")
include("StructureControlJob.jl")

export ModeClassification, SignificantMode, RAYLEIGH_LABEL_TOL, CLASSIFICATION_VERSION,
       matchRitzToGreenPeak, matchSpotPeak, matchScoreOf, GREEN_TIE_REL, associateLimitingMode,
       associateLimitingGroup, degenerateLimitGroups, LIMIT_ORTH_TOL, identityWeight
export LimitingAndBulk, rotateStructure, resolveParallel, forEachElement, forEachElementPooled,
       validateCalcMode, validateRetentionLevel, validateConstructorFlags, validateSawThreshold,
       resolveSawClassMethod, modeSawClassMethod, classificationMethods, validateSawClassMethod,
       defaultThreshold, normaliseThreshold, orList,
       ritzSAWPrecompute, ritzSAWPrecompute!, completeRitzPrecompute, RitzWavelengthSweep, ritzWavelengthSweep, limitingAndBulk, greenScoreSpectrum,
       greenComponentSpectrum, parseGreenComponent, greenSpectrum, buildAdaptiveVelAxis,
       pruneRitzTo, isStale, refuseStale, greenComponentOf, greenComponentMapOf,
       classifyGreenPath, greenStep, ritzOnlyStep, sawMapsFrom, greenG13EffMapFrom
export CalcSingleHS, CalcSingleLayered, CalcSinglePlate, CalcIsofreqSingleLayered, RefinedSaw,
       calcRitz, calcLim, calcAll, calcGreen, calcGreenScore, calcRitzScore, calcRefined,
       calcSpotGreen, identifyModes, applyRetention, classify, reclassify, greenAxis, greenComponent,
       greenComponentMap, ritzPrecompute, refinedSawAgreement, setRunMode, commonGreenAxis,
       setElement, setGreenVelAxis, finalizeMaps
export CalcAngSweepHS, CalcAngSweepLayered, CalcAngSweepPlate, CalcIsofreqAngSweepLayered,
       calcFullPhases, assembleMaps, calcSpotGreenMap
export singleHS, angSweepHS, singleLayered, angSweepLayered,
       isofreqSingleLayered, isofreqAngSweepLayered, singlePlate, angSweepPlate,
       modeClassification, freqToOmega0,
       wavevecAngSweep2Layered, wavevecAngSweep2Plate, freqAngSweep2Layered
export CalcWavevecAngSweep2Layered, CalcWavevecAngSweep2Plate, CalcFreqAngSweep2Layered,
       ZGVPoint, ClassifyBranchesResult, locateZGVs, classifyBranches, reindexByBranch
export SawSpectrum, GreenPeaks, GreenScore, GreenScoreLayered, RitzScore, GreenSpotScore,
       SpotGreen, DispersionBranches, thresholdValue, surfaceVRef, shearRange, scoreGreenSolve
export StructureControlSpec, LayerSpec, build_job, build_structure, build_hs_structure

end # module ForwardCalcs
