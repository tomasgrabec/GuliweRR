classdef SingleAngleViewerGUIPlate < handle
%%% SingleAngleViewerGUIPlate -- single-angle "curve" viewer window for a
%%% CalcSinglePlate result. 1D analogue of AngSweepViewerGUIPlate, the Plate
%%% counterpart of SingleAngleViewerGUI.
%%%
%%% A DELIBERATELY SEPARATE class from SingleAngleViewerGUI, not a
%%% class-dispatch branch inside it (same "separate parallel sibling"
%%% convention as ModalShapeViewerGUIPlate -- see its own docblock).
%%% Differences from SingleAngleViewerGUI, all following from
%%% CalcSinglePlate.m's own docblock:
%%%   - No graphTypeDropdown at all -- ABSENT, not disabled
%%%     (DispersionViewerGUIPlate's own precedent). There is exactly one
%%%     graph type here (the Ritz-mode curve); no Green's-function curve
%%%     (no analogue exists), no classification curve (no
%%%     significant-vs-insignificant distinction to classify INTO -- every
%%%     Ritz mode already IS a genuine physical branch).
%%%   - criterionDropdown offers 'Mode type' | 'Confidence'
%%%     (Vis.plotRitzModeCurvePlate's own 'modeType'/'confidence'), replacing
%%%     'Ritz score'/'freqDisp'/'depthInsensitivity' (all SAW-specific, no Plate
%%%     analogue).
%%%   - No 'Show types:' per-type checkbox filter. plotRitzModeCurvePlate
%%%     returns ONE scatter handle (colour-coded by type via a colour
%%%     matrix), not one handle per type the way plotClassificationCurve
%%%     does -- there is nothing to individually show/hide without splitting
%%%     that function into 4 series.
%%%   - onCurveClicked's info panel reports mode type + confidence
%%%     (RitzSolverLamb.modesType/modeTypeConfidence) instead of a
%%%     classification.significantModes lookup.
%%%
%%% Implements the same small duck-typed interface AngSweepViewerGUIPlate
%%% does (.fig, .show(), .OnModeSelected property, .showCrosshairAt(angle,
%%% vel)) so CombinedViewerGUIPlate/ModalShapeViewerGUIPlate's wiring works
%%% unmodified against either class.
%%%
%%% GuliweRR v6.0

    properties
        fig
        curveAxes
        criterionLabel
        criterionDropdown
        xMinField, xMaxField, yMinField, yMaxField
        infoPanel
        statusLine

        OnModeSelected = []   % function handle @(angle, modeNo, velSnapped, velRaw), for a linked modal-shape viewer
    end

    properties (SetAccess = private)
        singleResult   % CalcSinglePlate, already calc()'d
    end

    properties (Access = private)
        hCrosshair = gobjects(0)
        hCurvePlot = gobjects(0)
    end

    methods
        function obj = SingleAngleViewerGUIPlate(singleResult, parentContainer)
        %%% Construct the Plate single-angle curve viewer window
        %%% >> obj = SingleAngleViewerGUIPlate(singleResult, parentContainer)
        %%%
        %%% INPUT: 1. singleResult    .. CalcSinglePlate, already calc()'d
        %%%  (opt) 2. parentContainer .. a uigridlayout region to build into
        %%%                              instead of a fresh own uifigure --
        %%%                              used by CombinedViewerGUIPlate
            validateattributes(singleResult, {'CalcSinglePlate'}, ...
                {'scalar'}, 'SingleAngleViewerGUIPlate', 'singleResult');
            if isempty(singleResult.ritzSolver) || singleResult.ritzSolver.nCompModes == 0
                error('SingleAngleViewerGUIPlate:constructor:notCalculated', ...
                    'singleResult has not been calc()''d yet -- nothing to view.');
            end
            obj.singleResult = singleResult;
            if nargin < 2; parentContainer = []; end
            obj.buildFigure(parentContainer);
            obj.onGraphTypeChanged();
        end

        function show(obj)
        %%% Make the window visible (constructed invisible so headless/test use is safe)
        %%% >> show(obj)
            obj.fig.Visible = 'on';
        end
    end

    %%% ===================================================================
    %%% Figure construction
    %%% ===================================================================
    methods (Access = private)
        function buildFigure(obj, parentContainer)
            [obj.fig, mainLayout] = GuiCommon.hostLayout(parentContainer, ...
                'GuliweRR -- Single-Angle Viewer (Plate)', [100 100 1000 700], [3 1]);
            mainLayout.RowHeight = {40, '1x', 26};

            %%% top bar -- criterion selector (no graph-type dropdown, see class docblock)
            topBar = uigridlayout(mainLayout, [1 2]);
            topBar.Layout.Row = 1;
            topBar.ColumnWidth = {80, 160};
            topBar.Padding = [0 0 0 0];

            obj.criterionLabel = uilabel(topBar, 'Text', 'Criterion:');
            obj.criterionDropdown = uidropdown(topBar, ...
                'Items', {'Mode type', 'Confidence'}, 'ItemsData', {'modeType', 'confidence'}, ...
                'ValueChangedFcn', @(~,~) obj.onGraphTypeChanged());

            %%% curve + right-hand side panel (axis-range controls, info panel)
            [obj.curveAxes, sideGl, content] = GuiCommon.plotWithSidePanel(mainLayout, ...
                {18, 26, 18, 26, '1x'}, @(src, event) obj.onCurveClicked(src, event));
            content.Layout.Row = 2;

            ritzVels = obj.singleResult.ritzSolver.vels(1:obj.singleResult.ritzSolver.nCompModes);
            onRangeEdit = @(~,~) obj.applyAxisRange();
            [obj.xMinField, obj.xMaxField] = GuiCommon.rangeFieldPair(sideGl, ...
                'Velocity [km/s] min/max', num2str(min(ritzVels)), num2str(max(ritzVels)), onRangeEdit);
            [obj.yMinField, obj.yMaxField] = GuiCommon.rangeFieldPair(sideGl, ...
                'Value [Y] min/max', '', '', onRangeEdit);

            obj.infoPanel = uilabel(sideGl, 'Text', 'Click the curve to select a mode.', ...
                'VerticalAlignment', 'top');
            obj.infoPanel.Layout.Row = 5;

            obj.statusLine = GuiCommon.statusLabel(mainLayout, 3);
        end
    end

    %%% ===================================================================
    %%% Callbacks -- public (not Access=private), matching
    %%% SingleAngleViewerGUI's own convention.
    %%% ===================================================================
    methods
        function onGraphTypeChanged(obj)
        %%% Redraw the curve for whichever criterion is now selected
        %%% >> onGraphTypeChanged(obj)
            cla(obj.curveAxes);
            obj.hCrosshair = gobjects(0);

            try
                [hOut, ~] = Vis.plotRitzModeCurvePlate(obj.singleResult, obj.criterionDropdown.Value, [], obj.curveAxes);
                obj.hCurvePlot = hOut;
                set(obj.hCurvePlot, 'ButtonDownFcn', @(src, event) obj.onCurveClicked(src, event));
                GuiCommon.setStatusOk(obj.statusLine, '');
            catch err
                obj.hCurvePlot = gobjects(0);
                GuiCommon.setStatusError(obj.statusLine, err.message);
            end

            obj.applyAxisRange();
        end

        function applyAxisRange(obj)
        %%% Apply the velocity/value range fields to the curve axes
        %%% >> applyAxisRange(obj)
            GuiCommon.applyFieldLimit(obj.curveAxes, 'x', obj.xMinField, obj.xMaxField);
            GuiCommon.applyFieldLimit(obj.curveAxes, 'y', obj.yMinField, obj.yMaxField);
        end

        function setAxisRange(obj, axisName, lims)
        %%% Set one axis range from code, as if typed into its min/max fields
        %%% >> setAxisRange(obj, axisName, lims)    e.g. setAxisRange(obj, 'y', [NaN 7])
        %%%
        %%% axisName 'x' | 'y' (no colour axis here); lims [lo hi], a NaN
        %%% entry keeps that limit, [] keeps both -- see GuiCommon.setAxisRange
            GuiCommon.setAxisRange(obj, obj.curveAxes, axisName, lims);
        end

        function onCurveClicked(obj, ~, event)
        %%% Snap a curve click to the nearest computed Ritz mode by velocity.
        %%% Fires OnModeSelected.
        %%% >> onCurveClicked(obj, src, event)
            pt = GuiCommon.clickPoint(obj.curveAxes, event);
            velPos = pt(1);

            ritz  = obj.singleResult.ritzSolver;
            nM    = ritz.nCompModes;
            vels  = ritz.vels(1:nM);
            [~, modeNo] = min(abs(vels - velPos));
            velSel   = vels(modeNo);
            angleVal = obj.singleResult.angle;

            obj.infoPanel.Text = { ...
                sprintf('angle = %.1f deg', angleVal), ...
                sprintf('mode # = %d', modeNo), ...
                sprintf('velocity = %.4f km/s', velSel), ...
                sprintf('type = %s', GuiCommon.plateModeTypeText(ritz, modeNo))};

            obj.showCrosshairAt(angleVal, velSel);

            if ~isempty(obj.OnModeSelected)
                obj.OnModeSelected(angleVal, modeNo, velSel, velPos);
            end
        end

        function showCrosshairAt(obj, ~, velVal)
        %%% Move (or create) the crosshair marker at velVal, WITHOUT
        %%% triggering any click/mode-lookup logic. First argument (angle)
        %%% accepted-but-unused, matching AngSweepViewerGUIPlate's own
        %%% showCrosshairAt(angle, vel) signature.
        %%% >> showCrosshairAt(obj, angleVal, velVal)
            yPos = obj.valueAtVelocity(velVal);
            obj.hCrosshair = GuiCommon.markCrosshair(obj.curveAxes, obj.hCrosshair, velVal, yPos);
        end
    end

    %%% ===================================================================
    %%% Run-time helpers
    %%% ===================================================================
    methods (Access = private)
        function yPos = valueAtVelocity(obj, velVal)
        %%% Y-value the currently-displayed curve would show at velVal --
        %%% nearest mode's own criterion value (the curve is a discrete Ritz
        %%% spectrum, not a continuous function, so "nearest" is the only
        %%% sensible interpolation, same as SingleAngleViewerGUI's own
        %%% Ritz-mode-curve case).
        %%% >> yPos = valueAtVelocity(obj, velVal)
            ritz = obj.singleResult.ritzSolver;
            nM   = ritz.nCompModes;
            vels = ritz.vels(1:nM);
            [~, idx] = min(abs(vels - velVal));
            if strcmp(obj.criterionDropdown.Value, 'modeType')
                yPos = ritz.modesType(idx);
            else
                yPos = ritz.modeTypeConfidence(idx);
            end
        end
    end
end
