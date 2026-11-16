classdef Gui
%%% Namespace for GuliweRR's GUI viewer windows.
%%% Static methods are thin launchers over the underlying ...GUI handle
%%% classes in guliwerr/gui/ -- the same @ClassName/ pattern as
%%% Christoffel/GammaMatrix/GreenSolver/Calc/Vis.
%%% Do not instantiate this class.
%%%
%%% Usage:  win = Gui.angSweepViewer(sweepResult)
%%%         win = Gui.viewSaved('results/<file>.mat')   % result saved by Calc.saveResult

    methods (Static)
        obj = angSweepViewer(sweepResult, parentContainer, classifyTarget)
        obj = angSweepViewerPlate(sweepResult, parentContainer)
        obj = combinedViewer(sweepResult, windowMode)
        obj = combinedViewerPlate(sweepResult, windowMode)
        obj = dispersionViewer(sweepResult, windowMode, angleIdx, parentContainer)
        obj = dispersionViewerPlate(sweepResult, windowMode, angleIdx, parentContainer)
        obj = doubleSweepViewer(sweepResult)
        obj = doubleSweepViewerPlate(sweepResult)
        obj = freqAngSweep2LayeredViewer(sweepResult)
        obj = isofreqAngSweepViewer(sweepResult, parentContainer)
        obj = isofreqCombinedViewer(sweepResult, windowMode)
        obj = isofreqModalShapeViewer(sweepResult, windowMode, controlsParent, shapeParent)
        obj = isofreqSingleLayeredViewer(singleResult, parentContainer)
        obj = modalShapeViewer(sweepResult, windowMode, controlsParent, shapeParent)
        obj = modalShapeViewerPlate(sweepResult, windowMode, controlsParent, shapeParent)
        obj = singleAngleViewer(singleResult, parentContainer)
        obj = singleAngleViewerPlate(singleResult, parentContainer)
        obj = structureControl(initialStructure)

        %%% Any Calc* result -> its default viewer; viewSaved loads one from file first
        win = viewResult(calcResult, windowMode, showNow)
        [win, calcResult, meta] = viewSaved(filePath, windowMode, showNow)
    end
end
