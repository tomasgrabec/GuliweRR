classdef FreqAngSweep2LayeredViewerGUI < handle
%%% FreqAngSweep2LayeredViewerGUI -- interactive viewer for a
%%% CalcFreqAngSweep2Layered result (isoFreq 2D sweep: angle x real
%%% reference frequency omega0, each omega0 its own angle x kParRange Green
%%% |G_33| heatmap slice -- see CalcFreqAngSweep2Layered's own docblock).
%%% An omega0 stepper rebuilds IsofreqModalShapeViewerGUI +
%%% IsofreqAngSweepViewerGUI against a different sweepResult.sweeps{iOmega}
%%% element on every step.
%%%
%%% Does not embed IsofreqCombinedViewerGUI, for the same reason
%%% DoubleSweepViewerGUI does not embed CombinedViewerGUI: the composite
%%% class has no parentContainer support, only its two children do. This
%%% class instead directly builds+wires IsofreqModalShapeViewerGUI/
%%% IsofreqAngSweepViewerGUI into its own layout on every omega0 step, a
%%% direct copy of DoubleSweepViewerGUI.buildAngularContent/onWavChanged
%%% substituting sweeps{currentOmegaNo} for sweeps{currentWavNo}.
%%%
%%% UNLIKE DoubleSweepViewerGUI, there is no second mode to switch to (no
%%% dispersion-curve view exists for isoFreq -- no Ritz, no branch
%%% tracking): every omega0 step rebuilds the same one view, always
%%% resetting to the map. There is no "drilled-in" state worth carrying
%%% over either, since the underlying CalcIsofreqSingleLayered object at a
%%% given angle is a genuinely different object at a different omega0.
%%% Single window only.
%%%
%%% See docs/reference/architecture/gui_viewers.md for the full isoFreq viewer family.
%%%
%%% GuliweRR v6.0

    properties
        fig
        statusLine
    end

    properties (SetAccess = private)
        sweepResult    % CalcFreqAngSweep2Layered, already calc()'d
        currentOmegaNo % index into sweepResult.omega0List
    end

    properties (Access = {?FreqAngSweep2LayeredViewerGUI, ?test_FreqAngSweep2LayeredViewerGUI})
        contentParent
        omegaDropdown           % kept in sync by onOmegaChanged, so << / >> update it too
        modalShapeViewer = []   % IsofreqModalShapeViewerGUI, current omega0 slice
        angSweepViewer   = []   % IsofreqAngSweepViewerGUI, current omega0 slice
        savedAxisRange   = []   % the map's range-field values, carried across an
                                 % omega0 step -- see onOmegaChanged
        %%% Friend access to the test class -- matches
        %%% DoubleSweepViewerGUI.angularModalShapeViewer's own precedent, so
        %%% tests can confirm WHICH sweeps{i} element a step actually rebuilt
        %%% against without those needing to be part of the public surface.
    end

    methods
        function obj = FreqAngSweep2LayeredViewerGUI(sweepResult)
        %%% Construct the isoFreq frequency-sweep viewer
        %%% >> obj = FreqAngSweep2LayeredViewerGUI(sweepResult)
        %%%
        %%% INPUT: 1. sweepResult .. CalcFreqAngSweep2Layered, already calc()'d
            validateattributes(sweepResult, {'CalcFreqAngSweep2Layered'}, ...
                {'scalar'}, 'FreqAngSweep2LayeredViewerGUI', 'sweepResult');
            if isempty(sweepResult.sweeps)
                error('FreqAngSweep2LayeredViewerGUI:constructor:notCalculated', ...
                    'sweepResult has not been calc()''d yet -- nothing to view.');
            end

            obj.sweepResult    = sweepResult;
            obj.currentOmegaNo = 1;

            obj.buildFigure();
            obj.buildContent();
        end

        function show(obj)
        %%% Make the window visible (constructed invisible so headless/test use is safe)
        %%% >> show(obj)
            obj.fig.Visible = 'on';
        end

        function setAxisRange(obj, axisName, lims)
        %%% Set one axis range of the bottom map/curve viewer from code
        %%% >> setAxisRange(obj, axisName, lims)    e.g. setAxisRange(obj, 'c', [-3 0])
        %%%
        %%% Forwarded to angSweepViewer.setAxisRange -- see GuiCommon.setAxisRange
            obj.angSweepViewer.setAxisRange(axisName, lims);
        end
    end

    %%% ===================================================================
    %%% Figure construction
    %%% ===================================================================
    methods (Access = private)
        function buildFigure(obj)
            obj.fig = uifigure('Name', 'GuliweRR -- isoFreq Frequency Sweep Viewer', ...
                'Position', [80 60 1200 900], 'Visible', 'off');
            mainLayout = uigridlayout(obj.fig, [3 1]);
            mainLayout.RowHeight = {36, '1x', 26};
            mainLayout.Padding = [4 4 4 4];

            omegaStrs = arrayfun(@(w) sprintf('%.4g', w), obj.sweepResult.omega0List, 'UniformOutput', false);
            obj.omegaDropdown = GuiCommon.dropdownStepperRow(mainLayout, 1, 'omega0:', omegaStrs, ...
                {80, 160, 36, 36}, @(src,~) obj.onOmegaChanged(src.Value), ...
                @(~,~) obj.onPrevOmega(), @(~,~) obj.onNextOmega());
            obj.omegaDropdown.Value = obj.currentOmegaNo;

            %%% stable content slot -- buildContent() swaps what's INSIDE it
            %%% (never the slot itself), same pattern
            %%% DoubleSweepViewerGUI/AngSweepViewerGUI use
            obj.contentParent = uigridlayout(mainLayout, [1 1]);
            obj.contentParent.Layout.Row = 2;
            obj.contentParent.Padding = [0 0 0 0];

            obj.statusLine = GuiCommon.statusLabel(mainLayout, 3);
        end

        function buildContent(obj)
        %%% (Re)build IsofreqModalShapeViewerGUI + IsofreqAngSweepViewerGUI on
        %%% sweepResult.sweeps{currentOmegaNo}, wired exactly as
        %%% IsofreqCombinedViewerGUI wires them -- direct copy of
        %%% DoubleSweepViewerGUI.buildAngularContent, see class docblock.
        %%% >> buildContent(obj)
            delete(obj.contentParent.Children);

            %%% same 34/66 top/bottom split + thin divider
            %%% IsofreqCombinedViewerGUI's own 'combined' windowMode uses
            [shapeParent, controlsParent, bottomParent] = ...
                GuiCommon.dividedTopBottom(obj.contentParent, {'34x', 3, '66x'}, [0 0 0 0]);

            sliceResult = obj.sweepResult.sweeps{obj.currentOmegaNo};

            obj.modalShapeViewer = IsofreqModalShapeViewerGUI(sliceResult, [], controlsParent, shapeParent);
            obj.angSweepViewer   = IsofreqAngSweepViewerGUI(sliceResult, bottomParent);
            GuiCommon.restoreAxisRange(obj.angSweepViewer, obj.savedAxisRange);

            obj.angSweepViewer.OnPointSelected = ...
                @(angle, kParVal) obj.modalShapeViewer.goToPoint(angle, kParVal);
            obj.modalShapeViewer.OnSelectionChanged = ...
                @(angle, kParVal) obj.angSweepViewer.showCrosshairAt(angle, kParVal);

            GuiCommon.setStatusOk(obj.statusLine, ...
                sprintf('omega0 = %.4g', obj.sweepResult.omega0List(obj.currentOmegaNo)));
        end
    end

    %%% ===================================================================
    %%% Public entry point + callbacks -- public (not Access=private),
    %%% matching DoubleSweepViewerGUI's own convention so tests can call
    %%% these directly.
    %%% ===================================================================
    methods
        function onOmegaChanged(obj, newOmegaNo)
        %%% omega0 dropdown callback -- rebuilds against a different
        %%% sweeps{iOmega} element (see buildContent).
        %%% >> onOmegaChanged(obj, newOmegaNo)
            %%% keep the map's current ranges across the slice rebuild, typed or
            %%% set by setAxisRange (the fields are dead while drilled in -- the
            %%% last saved range is then kept instead); same as
            %%% DoubleSweepViewerGUI.onWavChanged
            if ~obj.angSweepViewer.isDrilledIn
                obj.savedAxisRange = GuiCommon.captureAxisRange(obj.angSweepViewer);
            end
            obj.currentOmegaNo = newOmegaNo;
            obj.omegaDropdown.Value = newOmegaNo;
            obj.buildContent();
        end

        function onPrevOmega(obj)
            obj.onOmegaChanged(max(obj.currentOmegaNo - 1, 1));
        end

        function onNextOmega(obj)
            obj.onOmegaChanged(min(obj.currentOmegaNo + 1, numel(obj.sweepResult.omega0List)));
        end
    end
end
