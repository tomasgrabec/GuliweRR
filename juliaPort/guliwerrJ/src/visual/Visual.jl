# =========================================================================================
# Visual.jl -- plotting-data preparation and the plotting entry points (no Makie here)
#
# `Visual` is the module behind MATLAB's `@Vis` static namespace: the numeric helpers that
# prepare arrays for plots (`VisData.jl`, `VisResults.jl`, `PlotData.jl`), the colormaps
# (`Colormaps.jl`), and the plotting/GUI function stubs that the Makie extensions fill in
# (`PlotStubs.jl`, `GuiStubs.jl`). Nothing in this module draws; it only returns finished
# arrays, so it loads without Makie. `Vis` itself names no Julia type -- its methods are
# this module's top-level functions.
#
# Dependencies: `CoreTypes`, `Solvers` and `ForwardCalcs` (sweep-result and classification
# types). `InverseCalcs` is deliberately not imported: no file here references one of its
# types or functions, and importing it unused would only be dead weight.
#
# `Visual` is included directly into `GuliwerrJ`, one lexical level below it, so sibling
# modules are reached with two dots (`using ..CoreTypes`).
#
# `VisResults.jl` holds the helpers that read the Calc results of the classification system
# (greenComponentData, the Ritz spectrum, ritzModeScores, classificationScoreScale, the
# Ritz-mode criteria, classification marker sizes) and the polarisation-RGB view
# (calcPolarisationRgb).
#
# Include order: `VisData.jl` before `PlotData.jl`, so a reader meets each helper before its
# first caller. Julia would resolve the cross-file calls at call time in any order.
# `Colormaps.jl` depends on neither.
# =========================================================================================

module Visual

using LinearAlgebra: norm, dot, cross
import Contour
using Interpolations: interpolate, extrapolate, Gridded, Linear

using ..CoreTypes: Structure, Material, nlay, SetDirMat, GetStandardDirMat, GetRotElastTensor,
                    GetStandardElastTensor, GetRotDirMat, GetRotMat, Tensor2Youngs

using ..Solvers: RitzSolver, RitzSolverSAW, RitzSolverLamb, calcDisplacement,
                  GreenSolver, findRoots, phaseVels, phaseVelScan, groupVelVectorsForDir,
                  calc

using ..ForwardCalcs: orList, CalcSingleHS, CalcSingleLayered, CalcAngSweepHS,
                       CalcAngSweepLayered, CalcWavevecAngSweep2Layered, SignificantMode,
                       DispersionBranches, RitzScore, GreenScore, GreenScoreLayered,
                       GreenSpotScore, SawSpectrum, classificationMethods, parseGreenComponent,
                       greenComponent, greenComponentMap

include("Colormaps.jl")
include("VisData.jl")
include("VisResults.jl")
include("PlotData.jl")
include("PlotStubs.jl")
include("GuiStubs.jl")

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

end # module Visual
