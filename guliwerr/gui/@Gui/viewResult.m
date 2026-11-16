function win = viewResult(calcResult, windowMode, showNow)
%%% Open the default viewer for any Calc* result, picked from its class
%%% >> win = Gui.viewResult(calcResult, windowMode, showNow)
%%%
%%% INPUT:  1. calcResult .. any CalcCommon.resultClasses() object, already
%%%                          calc()'d / calcAll()'d
%%%  (opt)  2. windowMode .. passed on to viewers that take one ('combined' |
%%%                          'split' | 'fullSplit', as the chosen viewer allows);
%%%                          default = [] (the viewer's own default)
%%%  (opt)  3. showNow    .. true (default) to make the window(s) visible; false
%%%                          builds them invisible (headless/testing use)
%%% OUTPUT: win .. the viewer object
%%%
%%% Result class -> viewer:
%%%   CalcSingleHS/AngSweepHS/SingleLayered/AngSweepLayered -> Gui.combinedViewer
%%%   CalcSinglePlate/AngSweepPlate                         -> Gui.combinedViewerPlate
%%%   CalcIsofreqSingleLayered/IsofreqAngSweepLayered       -> Gui.isofreqCombinedViewer
%%%   CalcFreqAngSweep2Layered                              -> Gui.freqAngSweep2LayeredViewer
%%%   CalcWavevecAngSweep2Layered/Plate, by shape:
%%%     >1 angle AND >1 wavelength/wavevector -> Gui.doubleSweepViewer(Plate)
%%%     1 angle, >1 wavelength/wavevector     -> Gui.dispersionViewer(Plate)
%%%     1 wavelength/wavevector               -> Gui.combinedViewer(Plate) on sweeps{1}
%%%
%%% The table lives here rather than in GuiCommon: choosing what a window shows
%%% is a viewer decision, which GuiCommon deliberately never makes.

    if nargin < 2; windowMode = []; end
    if nargin < 3 || isempty(showNow); showNow = true; end
    CalcCommon.validateResultObject(calcResult, 'Gui:viewResult');

    switch class(calcResult)
        case {'CalcSingleHS', 'CalcAngSweepHS', 'CalcSingleLayered', 'CalcAngSweepLayered'}
            win = Gui.combinedViewer(calcResult, windowMode);
        case {'CalcSinglePlate', 'CalcAngSweepPlate'}
            win = Gui.combinedViewerPlate(calcResult, windowMode);
        case {'CalcIsofreqSingleLayered', 'CalcIsofreqAngSweepLayered'}
            win = Gui.isofreqCombinedViewer(calcResult, windowMode);
        case 'CalcFreqAngSweep2Layered'
            requireNoWindowMode(windowMode, 'FreqAngSweep2LayeredViewerGUI');
            win = Gui.freqAngSweep2LayeredViewer(calcResult);
        case 'CalcWavevecAngSweep2Layered'
            win = viewTwoAxisSweep(calcResult, windowMode, numel(calcResult.wavelengths), false);
        case 'CalcWavevecAngSweep2Plate'
            win = viewTwoAxisSweep(calcResult, windowMode, numel(calcResult.wavVecs), true);
        otherwise
            %%% unreachable while this switch and CalcCommon.resultClasses() agree
            error('Gui:viewResult:noViewer', 'No default viewer for class %s.', class(calcResult));
    end

    if showNow
        win.show();
    end
end

function win = viewTwoAxisSweep(calcResult, windowMode, nOuter, isPlate)
%%% Dispatch a wavelength/wavevector x angle sweep by its actual shape
    if isempty(calcResult.sweeps)
        error('Gui:viewResult:notCalculated', ...
            'calcResult has not been calc()''d yet -- nothing to view.');
    end
    nAngles = numel(calcResult.angles);
    if nOuter == 1
        %%% one wavelength: really an angular sweep
        if isPlate
            win = Gui.combinedViewerPlate(calcResult.sweeps{1}, windowMode);
        else
            win = Gui.combinedViewer(calcResult.sweeps{1}, windowMode);
        end
    elseif nAngles == 1
        if isPlate
            win = Gui.dispersionViewerPlate(calcResult, windowMode);
        else
            win = Gui.dispersionViewer(calcResult, windowMode);
        end
    else
        if isPlate
            requireNoWindowMode(windowMode, 'DoubleSweepViewerGUIPlate');
            win = Gui.doubleSweepViewerPlate(calcResult);
        else
            requireNoWindowMode(windowMode, 'DoubleSweepViewerGUI');
            win = Gui.doubleSweepViewer(calcResult);
        end
    end
end

function requireNoWindowMode(windowMode, viewerName)
%%% The chosen viewer has no windowMode argument -- refuse one rather than
%%% silently ignore it
    if ~isempty(windowMode)
        error('Gui:viewResult:windowModeUnsupported', ...
            '%s (the viewer for this result) takes no windowMode; got ''%s''.', ...
            viewerName, char(windowMode));
    end
end
