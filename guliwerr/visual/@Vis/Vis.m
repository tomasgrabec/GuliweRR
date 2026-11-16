classdef Vis
%%% Namespace for GuliweRR's plotting and visualisation functions.
%%% All methods are static -- do not instantiate this class.
%%%
%%% Usage:  Vis.plotModalShape(ritzSolver, modeNo)
%%%         Vis.plotAngularMap(sweepResult)

    methods (Static)
        %%% Angular / velocity-surface plots
        [hSaw, hAx]        = plotAngularMap(sweepResult, plotType, hAx)
        [hOut, hAx]        = plotVelSurface(mater, surfaceType, nTheta, nPhi, hFig)
        [hOut, hAx]        = plotGroupVelSurface(mater, nTheta, nPhi, hFig)
        [hOut, hAx]        = plotYoungsSurface(mater, nTheta, nPhi, hFig)
        [hOut, hAx]        = plotVelCurve(mater, angleDeg, hAx)
        [hOut, hAx]        = plotYoungsCurve(mater, angleDeg, hAx)
        [hOut, hAx]        = plotRitzModeMap(sweepResult, criterion, velCutoffFactor, hAx)
        [hOut, hAx]        = plotRitzModeMapPlate(sweepResult, criterion, velCutoffFactor, hAx)
        [hPoints, hAx, hBg, hLim] = plotClassificationDiagnostic(sweepResult, colorBy, showGreenBackground, refLevel, hAx, component)

        %%% Ritz spectrum |S(v)| (RitzScore.spectrum of the stored modal numerators),
        %%% log-normalised like plotGreenMap/plotGreenCurve; no solve
        [hOut, hAx, hCb]   = plotRitzSpectrumMap(sweepResult, scaleMode, climRange, plotType, hAx, showCands)

        %%% Single-angle curve plots (1D analogues of the maps above, for
        %%% CalcSingleHS/CalcSingleLayered results -- no angle axis to sweep)
        [hOut, hAx]        = plotRitzModeCurve(singleResult, criterion, velCutoffFactor, hAx)
        [hOut, hAx]        = plotRitzModeCurvePlate(singleResult, criterion, velCutoffFactor, hAx)
        [hOut, hAx]        = plotGreenCurve(singleResult, scaleMode, hAx, component)
        [hPoints, hAx, hBg, hLim] = plotClassificationCurve(singleResult, colorBy, showGreenBackground, refLevel, hAx, component)
        [hOut, hAx, spec]  = plotRitzSpectrumCurve(singleResult, scaleMode, hAx, showCands)

        %%% Green's-function / dispersion maps
        [hOut, hAx, hCb]   = plotGreenMap(sweepResult, scaleMode, climRange, plotType, hAx, component)
        [hOut, hAx]        = plotDispersionMap(sweepResult, omegaIdx, scaleMode, climRange, plotType, hAx)
        [hLines, hAx, hBg] = plotDispersionCurves(sweepResult, curveType, angleIdx, showRitzBranches, showGreenBackground, scaleMode, climRange, hAx)
        [hLines, hAx, hChristLines] = plotDispersionCurvesPlate(sweepResult, curveType, angleIdx, colorMode, confidenceThresh, showChristoffelBounds, hAx)
        [hOut, hAx]        = plotPolarHeatmap(angles, radius, C, hAx)

        %%% Grid helpers
        [xEdges, yEdges, CPadded] = pcolorCenteredGrid(x, y, C)

        %%% Modal-shape (displacement field) plots
        [hOut, hAxCell]    = plotModalShape(ritzSolver, modeNo, substrateCutoff, cutoffParam, dirs, meshSpacing, hFig, hAxCellIn)
        [hOut, hAxCell]    = plotGreenModalShape(greenObj, forceDir, cutoffParam, dirs, meshSpacing, hFig, hAxCellIn)
        [hOut, hAxCell]    = plotGreenModalShapeLayered(greenObj, structure, forceDir, substrateCutoff, cutoffParam, dirs, meshSpacing, hFig, hAxCellIn)
        [hPatch, hAx]      = plotPolarisationKey(hAx, nSub)

        %%% Field-extraction / geometry helpers
        dispField          = calcDispField2D(ritzSolver, modeNo, meshSpacing)
        dispField          = calcGreenDispField2D(greenObj, forceDir, meshSpacing)
        dispField          = calcGreenDispField2DLayered(greenObj, structure, forceDir, meshSpacing)
        [rgbCell, ampCell] = calcPolarisationRgb(dispField)
        cutoffRelDepth      = calcSubstrateCutoffDepth(substrateCutoff, cutoffParam, relDepths, substrateDispMag)
        [uCoord, wCoord, ptsXYZ] = extractPlanarCut(X, Y, Z, theta, phi, planeNormal, inPlaneRef)
        frac               = calcModalShapeRowHeights(relDepths, minFrac)

        %%% Display conventions shared with callers outside visual/ -- a GUI
        %%% overlaying its own annotation on one of these plots must use the same
        %%% scaling and vocabulary. Conventions used only INSIDE visual/ live in
        %%% @Vis/private/ and are not exposed here.
        [plotData, dataLabel]   = scaleGreenData(absG, scaleMode, componentLabel)
        [absG, componentLabel, velAxis] = greenComponentData(result, component)
        [typeNames, typeColors, parityNames, parityColors] = plateModeTypeStyle()
        s                       = classificationScoreScale(method, isLayered)
        scores                  = ritzModeScores(singleResult)

        %%% Colormaps
        cMap = ColMapGrayOrange()
        cMap = ColormapRWB()
        cMap = ColormapRYG()

        %%% Inverse-characterisation diagnostics (guliwerr/inverseCalcs/)
        [hBar, hAx, hHeat]  = plotSensitivitySpectrum(sensitivityAnalysis, hFig)
        [hOut, hAx]         = plotInversionFit(measured, associationFinal, associationInitial, hAx)
        [hOut, hAx]         = plotAssociationDiagnostic(association, hAx, overlayMode)
    end
end
