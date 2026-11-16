# =========================================================================================
# PlotStubs.jl -- drawing-function stubs for the Makie plotting extension
#
# Empty function stubs for every `guliwerr/visual/@Vis/plot*.m` drawing function (26 total).
# Each stub's real implementation lives in `ext/GuliwerrJMakieExt.jl` (a package extension,
# per Project.toml's `[extensions]` table, loaded automatically whenever the end user does
# `using GuliwerrJ, GLMakie` or `using GuliwerrJ, WGLMakie` -- the extension triggers on
# `Makie` being loaded, not a specific backend, so either works with zero code change).
#
# This is the same "stub now, extend later" pattern `CalcCommon.jl` uses for
# `calcGreen`/`identifyModes`/`applyRetention`: the extension's plotting sub-files add
# methods to these generic functions, and nothing in `Visual.jl`/`GuliwerrJ.jl` changes when
# plot methods are added.
#
# Calling one of these WITHOUT Makie loaded throws Julia's ordinary `MethodError: no method
# matching plotGreenMap(...)` -- there is no dispatch-table entry until the extension adds
# one. This is the expected failure mode for a package extension with no fallback bodies;
# no custom error message is needed.
#
# Every stub takes the ALREADY-COMPUTED data (a `PlotData.jl`/`VisData.jl` output, or a
# solver/result type directly) -- these functions only draw, they do no physics computation
# of their own; the compute half lives in `PlotData.jl`/`VisData.jl`.
# =========================================================================================

"""
    plotVelSurface(mater, surfaceType, nTheta, nPhi; fig=nothing) -> Figure

3D phase-velocity surface over all propagation directions. Port of `Vis.plotVelSurface`.
Implemented in the Makie extension; see `PlotData.velSurfaceData` for the numeric half.
"""
function plotVelSurface end

"""
    plotGroupVelSurface(mater, nTheta, nPhi; fig=nothing) -> Figure

3D group-velocity surface. Port of `Vis.plotGroupVelSurface`.
Implemented in the Makie extension; see `PlotData.groupVelSurfaceData`.
"""
function plotGroupVelSurface end

"""
    plotYoungsSurface(mater, nTheta, nPhi; fig=nothing) -> Figure

3D Young's-modulus surface. Port of `Vis.plotYoungsSurface`.
Implemented in the Makie extension; see `PlotData.youngsSurfaceData`.
"""
function plotYoungsSurface end

"""
    plotVelCurve(mater, angleDeg; ax=nothing) -> Axis

Phase-velocity vs. angle curve (bulk modes) in one plane. Port of `Vis.plotVelCurve`.
Implemented in the Makie extension; see `PlotData.velCurveData`.
"""
function plotVelCurve end

"""
    plotYoungsCurve(mater, angleDeg; ax=nothing) -> Axis

Young's-modulus vs. angle curve. Port of `Vis.plotYoungsCurve`.
Implemented in the Makie extension; see `PlotData.youngsCurveData`.
"""
function plotYoungsCurve end

"""
    plotGreenCurve(singleResult, scaleMode; ax=nothing) -> Axis

`|G_33|` vs. velocity for one angle. Port of `Vis.plotGreenCurve`.
Implemented in the Makie extension.
"""
function plotGreenCurve end

"""
    plotRitzSpectrumCurve(singleResult, scaleMode = :log; ax = nothing, showCands = false)
        -> (Axis, spec)

The Ritz spectrum `|S(v)|` of one angle, rebuilt from the stored modal numerators (no
solve); optional overlay of the Ritz candidates. Port of `Vis.plotRitzSpectrumCurve`.
Implemented in the Makie extension; see
`PlotData.ritzSpectrumCurveData`.
"""
function plotRitzSpectrumCurve end

"""
    plotRitzModeCurve(singleResult, criterion, velCutoffFactor; ax=nothing) -> Axis

Ritz per-mode criterion (`:ritzScore`/`:freqDisp`/`:depthInsensitivity`) vs. velocity for one angle,
half-space/layered case. Port of `Vis.plotRitzModeCurve`. Implemented in the Makie extension.
"""
function plotRitzModeCurve end

"""
    plotRitzModeCurvePlate(singleResult, criterion, velCutoffFactor; ax=nothing) -> Axis

Plate-geometry counterpart of [`plotRitzModeCurve`](@ref). Port of
`Vis.plotRitzModeCurvePlate`. Implemented in the Makie extension.
"""
function plotRitzModeCurvePlate end

"""
    plotGreenMap(sweepResult, scaleMode, climRange, plotType; ax=nothing) -> (Axis, Colorbar)

Green's-function angular map -- `|G_+|` (or the requested `plotType` component) heatmap,
angle x velocity, via [`pcolorCenteredGrid`](@ref)'s edge convention. Port of
`Vis.plotGreenMap`. Implemented in the Makie extension.
"""
function plotGreenMap end

"""
    plotDispersionMap(sweepResult, omegaIdx, scaleMode, climRange, plotType; ax=nothing) -> Axis

Dispersion map (2D wavevec/angle or freq/angle sweep) at one outer-loop index. Port of
`Vis.plotDispersionMap`. Implemented in the Makie extension; see
`PlotData.wavevecGreenBackgroundData`.
"""
function plotDispersionMap end

"""
    plotPolarHeatmap(angles, radius, C; ax=nothing) -> Axis

Generic polar heatmap helper shared by several angular-map plots. Port of
`Vis.plotPolarHeatmap`. Implemented in the Makie extension.
"""
function plotPolarHeatmap end

"""
    plotAngularMap(sweepResult, plotType; ax=nothing) -> Axis

SAW-velocity (or other scalar) polar map over a full angular sweep. Port of
`Vis.plotAngularMap`. Implemented in the Makie extension.
"""
function plotAngularMap end

"""
    plotRitzModeMap(sweepResult, criterion, velCutoffFactor; ax=nothing) -> Axis

[`plotRitzModeCurve`](@ref)'s map (all angles) counterpart. Port of `Vis.plotRitzModeMap`.
Implemented in the Makie extension.
"""
function plotRitzModeMap end

"""
    plotRitzModeMapPlate(sweepResult, criterion, velCutoffFactor; ax=nothing) -> Axis

Plate-geometry counterpart of [`plotRitzModeMap`](@ref). Port of `Vis.plotRitzModeMapPlate`.
Implemented in the Makie extension; consumes [`visPlateLimVels`](@ref) internally, per
`Vis.plotRitzModeMapPlate`'s own MATLAB source.
"""
function plotRitzModeMapPlate end

"""
    plotRitzSpectrumMap(sweepResult, scaleMode = :log, climRange = (-4, -1); ax = nothing,
                        showCands = false) -> (Axis, Colorbar)

Angle x velocity map of the Ritz spectrum `|S(v)|` (stored numerators, no solve; the default
background of a `:ritzOnly` result). Port of `Vis.plotRitzSpectrumMap`.
Implemented in the Makie extension; see
`PlotData.ritzSpectrumMapData`.
"""
function plotRitzSpectrumMap end

"""
    plotPolarisationKey(; ax = nothing) -> Axis

Colour key of the polarisation-RGB view (`dirs = :rgb`): the positive octant of the unit
sphere of directions coloured as [`calcPolarisationRgb`](@ref) colours them. Port of
`Vis.plotPolarisationKey`. Implemented in the Makie extension.
"""
function plotPolarisationKey end

"""
    plotClassificationCurve(singleResult, colorBy, showGreenBackground, refLevel; ax=nothing) -> Axis

Classified mode points overlaid on a Green's-function curve, for one angle. Port of
`Vis.plotClassificationCurve`. Implemented in the Makie extension.
"""
function plotClassificationCurve end

"""
    plotClassificationDiagnostic(sweepResult, colorBy, showGreenBackground, refLevel; ax=nothing) -> Axis

[`plotClassificationCurve`](@ref)'s map (all angles) counterpart, classified points over the
Green's-function background. Port of `Vis.plotClassificationDiagnostic`. Implemented in the
Makie extension.
"""
function plotClassificationDiagnostic end

"""
    plotDispersionCurves(sweepResult, curveType, angleIdx, showRitzBranches, showGreenBackground, scaleMode, climRange; ax=nothing) -> Axis

Tracked dispersion branches (from `trackDispersionBranches`) over a Green's-function
background. Port of `Vis.plotDispersionCurves`. Implemented in the Makie extension.
"""
function plotDispersionCurves end

"""
    plotDispersionCurvesPlate(sweepResult, curveType, angleIdx, colorMode, confidenceThresh, showChristoffelBounds; ax=nothing) -> Axis

Plate-geometry counterpart of [`plotDispersionCurves`](@ref), with Christoffel bulk-velocity
bound overlays. Port of `Vis.plotDispersionCurvesPlate`. Implemented in the Makie extension.
`colorMode = :type` colours a branch by its mode type where [`classifyBranches`](@ref) trusts the
label, by its parity class ([`plateModeTypeStyle`](@ref)'s `parityColors`, dashed) where only
the parity is exact, and by branch identity otherwise; `:branch` always by branch identity.
"""
function plotDispersionCurvesPlate end

"""
    plotModalShape(ritzSolver, modeNo, substrateCutoff, cutoffParam, dirs, meshSpacing; fig=nothing) -> Vector{Axis}

Signed 2D (depth x lateral) displacement-field plot for one Ritz mode, one or more
directions. Port of `Vis.plotModalShape`. Implemented in the Makie extension; see
`calcDispField2D` for the numeric field.
"""
function plotModalShape end

"""
    plotGreenModalShape(greenObj, forceDir, cutoffParam, dirs, meshSpacing; fig=nothing) -> Vector{Axis}

Signed 2D displacement field from a half-space Green's-function response. Port of
`Vis.plotGreenModalShape`. Implemented in the Makie extension; see `calcGreenDispField2D`.
"""
function plotGreenModalShape end

"""
    plotGreenModalShapeLayered(greenObj, structure, forceDir, substrateCutoff, cutoffParam, dirs, meshSpacing; fig=nothing) -> Vector{Axis}

Layered-structure counterpart of [`plotGreenModalShape`](@ref). Port of
`Vis.plotGreenModalShapeLayered`. Implemented in the Makie extension; see
`calcGreenDispField2DLayered`.
"""
function plotGreenModalShapeLayered end

"""
    plotAssociationDiagnostic(association, overlayMode; ax=nothing) -> Axis

`ModeAssociation` points (measured vs. calculated velocity, binding class) diagnostic plot.
Port of `Vis.plotAssociationDiagnostic`. Implemented in the Makie extension.
"""
function plotAssociationDiagnostic end

"""
    plotInversionFit(measured, associationFinal, associationInitial; ax=nothing) -> Axis

Measured-vs-fitted velocity comparison, before/after inversion. Port of
`Vis.plotInversionFit`. Implemented in the Makie extension.
"""
function plotInversionFit end

"""
    plotSensitivitySpectrum(sensitivityAnalysis; fig=nothing) -> (Axis, Axis)

Sensitivity eigenvalue spectrum (bar chart) plus transform-matrix heatmap. Port of
`Vis.plotSensitivitySpectrum`. Implemented in the Makie extension.
"""
function plotSensitivitySpectrum end

export plotVelSurface, plotGroupVelSurface, plotYoungsSurface, plotVelCurve, plotYoungsCurve,
       plotGreenCurve, plotRitzSpectrumCurve, plotRitzModeCurve, plotRitzModeCurvePlate,
       plotGreenMap, plotDispersionMap, plotPolarHeatmap, plotAngularMap, plotRitzModeMap,
       plotRitzModeMapPlate, plotRitzSpectrumMap, plotPolarisationKey, plotClassificationCurve,
       plotClassificationDiagnostic, plotDispersionCurves, plotDispersionCurvesPlate,
       plotModalShape, plotGreenModalShape, plotGreenModalShapeLayered,
       plotAssociationDiagnostic, plotInversionFit, plotSensitivitySpectrum
