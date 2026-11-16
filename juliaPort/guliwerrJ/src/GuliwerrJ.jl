"""
    GuliwerrJ

Julia port of GULIWERR (guided-wave forward/inverse elasticity solver, originally MATLAB).
See `juliaPort/docsJ/quickStartJ.md` for how to use it.

This top-level module does two things only: `include` every source file in dependency
order, and re-export the public names each submodule defines. No numerical code lives
here directly.

Submodules: `CoreTypes` (Legendre basis, `CalcParams`, `ElastConsts`, `Material`,
`Structure`, `ParallelPolicy`), `Solvers` (Christoffel bulk-wave and limiting velocities,
Ritz, Green), `ForwardCalcs`, `InverseCalcs` and `Visual`.
"""
module GuliwerrJ

include("core/CoreTypes.jl")
using .CoreTypes
include("solvers/Solvers.jl")
using .Solvers
include("forwardCalcs/ForwardCalcs.jl")
using .ForwardCalcs
include("inverseCalcs/InverseCalcs.jl")
using .InverseCalcs
include("visual/Visual.jl")
using .Visual

# Re-export the public surface of every submodule.
export matlab_roots, movmedian, matchpairs0, rcond_matlab, fminbnd_brent, findpeaks_matlab
export legendre_norm, basis_decomposition_2d, legendre_matrix_diag_nondiag
export CalcParams, set_rel_depth, set_real_thickness, set_leg_order
export ElastConsts, Voigt2Vector, Vector2Voigt, Any2Tensor, Tensor2Voigt, Tensor2Vector, ZenerRatio
export SymmClassLibrary, FindSymmClass_input, IndepVoigtPositions, GetIndOf2DReduct_15of21,
       CubicConversion, CubicDeconversion, HexagonalConversion, HexagonalDeconversion,
       TetragonalConversion, TetragonalDeconversion
export RotateTensor, ElastRotReductMat, RightHandFrame, Tensor2Youngs, Voigt2Bulk,
       HillAverages, VoigtAverages, ReussAverages
export Material, SetElastConsts, SetDirMat, SetDensity, SetName,
       GetRotMat, GetRotElastTensor, GetRotElastVoigt, GetRotDirMat, GetStandardElastTensor
export Euler2DCM, DCM2Euler, Euler2DirMat, DirMat2Euler, Miller2DirMat, DirMat2Miller,
       CompleteDirMat, GetStandardDirMat, RotDirMatFromAngle, CalcRotMat
export MaterialFromLibrary, ListMaterialLibrary
export Structure, nlay, str_basis_vec_l, AddLayer, RemoveLayer, ModifyLayer,
       ApplyWavelength, isMirrorSymmetric, isSagittalMirrorSymmetric
export AUTO_THRESHOLD, should_use_parallel, threads_available, with_blas_threads
export blasBackend, with_task_blas_single, BLAS_ENV_VAR
export guliwerrVersion
export ProbeKernel, sample, probe_point, probe_gaussian, probe_tgs, probe_sras, probe_from_real_space, probe_validate_component, probe_default_dir_step, probe_default_relk_step
export matrix, GammaFromKron, phaseVels, phaseVelScan
export phaseVelAngleDeriv, groupVelVectorsForDir
export calc_anisotropy_ratio
export LimitingRoot, findRoots, refineRoot
export GammaMatrix, SAWBoundary, LambBoundary, calc,
       rawByDepth, derivsByDepth, elastByDepth, dummyForLayer,
       RawCellFor, LayerRotMat, RedElastCellFor, LayerNormFac, ContractDepth
export RitzSolver, solveConstrainedEig, realGauge, gaugeVector, isGaugeReal, phaseFix,
       REAL_GAUGE_TOL, precomp_field,
       buildConstraintMatrix, solveEigenvalueProblem,
       calcGroupVelocity, calcDisplacement, angleSlope, pruneEigVecs
export RitzSolverSAW, RitzDiagnostics, modeDiagnostics, calcWithGamma, ritzSensitivities,
       ritzSensitivities!, RitzWorkspace
export RitzSolverLamb, modeRecognition, parityFractions, fastCalcGroupVel, locateZGV,
       calcTransGroupVelocity, kPolynomial, ZgvNewton
export GreenSolver, nvel, has_minus, has_half_space_full_output, has_layered_full_output
export GreenSexticPrepass, green_sextic_prepass, green_layer_setup, green_layers, GreenLayer,
       calc_detw_by_christoffel, solve_k3_allvel, canonicalize_k3,
       check_identical_records_allvel,
       calc_u_allvel, substrate_wave_indices, kvec_cols, outer_cols,
       OMEGA_TOLERANCE, ZERO_TOLERANCE, LIMITING_ABSORPTION_ETA,
       NEAR_ISO_ANGLE_STEP, NEAR_ISO_WARN_THRESH
export calc_traction_cols, calc_gfull_from_partial_waves, build_bc_matrix,
       green_component, green_gfull, DISP_FIELD_NDEPTH, DISP_FIELD_MAX_DEPTH, GREEN_CHUNK
export SawSpectrum, GreenPeaks, GreenScore, GreenScoreLayered, RitzScore, GreenSpotScore,
       SpotGreen, DispersionBranches, thresholdValue, surfaceVRef, shearRange, scoreGreenSolve
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
export StructureControlSpec, LayerSpec, build_job, build_structure, build_hs_structure
export WAVE_TYPES, MeasuredVelocities, PhaseVelocityEntry, addPhaseVelocity, setBulkModulus,
       fromTable
export ModeAssociation, AssocPoint, LIM_VELOCITY_WINDOW, LIM_THETA_WINDOW, LIM_BRANCH_PENALTY,
       RITZ_CLASSES, DEFAULT_TRACK_OPTS, DEFAULT_REBIND_WINDOW, activeRitzClasses,
       nPoints, uniqueRitzAngles, uniqueLimitingAngles, pointsAtAngle, confidence,
       pinnedPoints, setPinned, setLimitingAnchor, setLimitingRefined, setLimitingLost,
       setGreenPeak, addPoint, removeDuplicatePoints, rebind, rebindLimiting,
       fromSweep, defaultTrackOpts, limitingRootScore
export SensitivityAnalysis, DEFAULT_SENSITIVITY_CUTOFF, waveTypeIndex, phaseVelocityAt,
       bulkVelocityAt, limitingStep, phaseBulkStep, limitingRootsAt, limitingVelocityNear,
       limitingVelocityTracked, rotatedTensorFor, channelDerivative, channelValue,
       ritzVelocityDerivative, checkDatasets
export ElasticityInversion, calc!, recalcForwardAtResult, pseudoSAWGreenPeaks
export measured, association, sensitivity, elasticity
export ColMapGrayOrange, ColormapRWB, ColormapRYG, plateModeTypeStyle
export visDispField2D, visDecayCutoffIdx, visDepthMagnitude, visSphereDirGrid,
       visFlattenSignificantModes, visClassificationColorVals, visClassificationSizes,
       visIsLayered, visHasRitzSpectrum, visRitzSpectrumData, visRitzCriterion,
       visRitzCriterionNames, visValidateGreenComponent, visPolarisationColours, visIsRgbDirs,
       visPlateRitzModeData, visPropOrEmpty, visValidateChoice, visDispersionCurveTypes,
       visClassificationTypesAllowed, visStandardFrameTensor,
       visPlateLimVels, visPlateLongLimVel
export calcDispField2D, calcGreenDispField2D, calcGreenDispField2DLayered, extractPlanarCut,
       scaleGreenData, calcSubstrateCutoffDepth, calcModalShapeRowHeights, pcolorCenteredGrid,
       greenComponentData, ritzModeScores, classificationScoreScale, calcPolarisationRgb
export velSurfaceData, groupVelSurfaceData, youngsSurfaceData, velCurveData, youngsCurveData,
       ritzSpectrumMapData, ritzSpectrumCurveData, dispersionCurvesData,
       wavevecGreenBackgroundData
export plotVelSurface, plotGroupVelSurface, plotYoungsSurface, plotVelCurve, plotYoungsCurve,
       plotGreenCurve, plotRitzSpectrumCurve, plotRitzModeCurve, plotRitzModeCurvePlate,
       plotGreenMap, plotDispersionMap, plotPolarHeatmap, plotAngularMap, plotRitzModeMap,
       plotRitzModeMapPlate, plotRitzSpectrumMap, plotPolarisationKey, plotClassificationCurve,
       plotClassificationDiagnostic, plotDispersionCurves, plotDispersionCurvesPlate,
       plotModalShape, plotGreenModalShape, plotGreenModalShapeLayered,
       plotAssociationDiagnostic, plotInversionFit, plotSensitivitySpectrum
export structureControl, singleAngleViewer, angSweepViewer, modalShapeViewer, combinedViewer,
       doubleSweepViewer, dispersionViewer, setAxisRange!

end # module GuliwerrJ
