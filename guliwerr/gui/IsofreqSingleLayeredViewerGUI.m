classdef IsofreqSingleLayeredViewerGUI < handle
%%% IsofreqSingleLayeredViewerGUI -- single-angle, single-omega0 isoFreq
%%% "curve" viewer window, the isoFreq counterpart of SingleAngleViewerGUI
%%% for CalcIsofreqSingleLayered results.
%%%
%%% Unlike SingleAngleViewerGUI (three graph types: Ritz-mode/Green's-
%%% function/classification curve), this family has exactly ONE possible
%%% graph -- the Green |G_33| field vs. kParRange -- since CalcIsofreqSingleLayered
%%% never runs a Ritz solve at all (fixed omega0, k-swept Green's-function
%%% scan only; see docs/reference/architecture/forwardCalcs.md, "isoWavevec vs.
%%% isoFreq"). So there is no graphTypeDropdown/criterionDropdown at all --
%%% absent, not disabled, matching the Plate viewer family's own precedent.
%%% Rendering reuses Vis.plotDispersionMap directly (it already dispatches
%%% on CalcIsofreqSingleLayered for its 1D line case), not a new dedicated
%%% plotting function.
%%%
%%% Implements a renamed duck-typed interface: OnPointSelected/
%%% showCrosshairAt(angle, kParVal), NOT SingleAngleViewerGUI's OnModeSelected/
%%% showCrosshairAt(angle, modeNo, velSnapped, velRaw) -- kParRange is a
%%% continuous, physically-meaningful axis here (CalcIsofreqSingleLayered's
%%% own docblock), not a discrete Ritz mode index, so the callback signature
%%% says so honestly instead of reusing a shape built for a different
%%% quantity. angle is accepted only for interface symmetry with
%%% IsofreqAngSweepViewerGUI (always obj.singleResult.angle here).
%%%
%%% See docs/reference/architecture/gui_viewers.md for the full isoFreq viewer family.
%%%
%%% GuliweRR v6.0

    properties
        fig
        curveAxes
        scaleModeDropdown
        xMinField, xMaxField, yMinField, yMaxField
        infoPanel
        statusLine

        OnPointSelected = []   % function handle @(angle, kParVal), set by IsofreqModalShapeViewerGUI
    end

    properties (SetAccess = private)
        singleResult   % CalcIsofreqSingleLayered, already calc()'d
    end

    properties (Access = private)
        hCrosshair = gobjects(0)
        hCurvePlot = gobjects(0)   % handle the active plot returned
    end

    methods
        function obj = IsofreqSingleLayeredViewerGUI(singleResult, parentContainer)
        %%% Construct the isoFreq single-angle curve viewer window
        %%% >> obj = IsofreqSingleLayeredViewerGUI(singleResult, parentContainer)
        %%%
        %%% INPUT: 1. singleResult    .. CalcIsofreqSingleLayered, already
        %%%                              calc()'d (obj.lastResult off
        %%%                              StructureControlGUI's Run button, typically)
        %%%  (opt) 2. parentContainer .. a uigridlayout region (or any valid
        %%%                              uigridlayout parent) to build into instead
        %%%                              of a fresh own uifigure -- used by
        %%%                              IsofreqCombinedViewerGUI to embed this
        %%%                              window's content into a shared figure.
        %%%                              Omitted (default): standalone, owns its
        %%%                              own uifigure.
            validateattributes(singleResult, {'CalcIsofreqSingleLayered'}, ...
                {'scalar'}, 'IsofreqSingleLayeredViewerGUI', 'singleResult');
            if isempty(singleResult.absG)
                error('IsofreqSingleLayeredViewerGUI:constructor:notCalculated', ...
                    'singleResult has not been calc()''d yet -- nothing to view.');
            end
            obj.singleResult = singleResult;
            if nargin < 2; parentContainer = []; end
            obj.buildFigure(parentContainer);
            obj.onScaleModeChanged();
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
                'GuliweRR -- isoFreq Single-Angle Viewer', [100 100 1000 700], [3 1]);
            mainLayout.RowHeight = {40, '1x', 26};

            %%% top bar -- scale-mode selector only (no graph-type/criterion --
            %%% see class docblock)
            topBar = uigridlayout(mainLayout, [1 2]);
            topBar.Layout.Row = 1;
            topBar.ColumnWidth = {80, 100};
            topBar.Padding = [0 0 0 0];

            uilabel(topBar, 'Text', 'Scale:');
            obj.scaleModeDropdown = uidropdown(topBar, 'Items', {'log', 'linear'}, ...
                'ValueChangedFcn', @(~,~) obj.onScaleModeChanged());

            %%% curve + right-hand side panel (axis-range controls, info panel)
            [obj.curveAxes, sideGl, content] = GuiCommon.plotWithSidePanel(mainLayout, ...
                {18, 26, 18, 26, '1x'}, @(src, event) obj.onCurveClicked(src, event));
            content.Layout.Row = 2;

            onRangeEdit = @(~,~) obj.applyAxisRange();
            [obj.xMinField, obj.xMaxField] = GuiCommon.rangeFieldPair(sideGl, ...
                'kParRange min/max', num2str(min(obj.singleResult.kParRange)), ...
                num2str(max(obj.singleResult.kParRange)), onRangeEdit);

            %%% Value (y) range -- defaults to Vis.plotDispersionMap's own
            %%% 1D-line default ([-5,0], top at the normalised peak) since 'log' is the initial scale
            %%% mode; onScaleModeChanged re-defaults these on every scale switch,
            %%% same "re-default per active mode" mechanic AngSweepViewerGUI's
            %%% colour-range fields use for its own three graph types.
            [obj.yMinField, obj.yMaxField] = GuiCommon.rangeFieldPair(sideGl, ...
                'Value [Y] min/max', '-5', '0', onRangeEdit);

            obj.infoPanel = uilabel(sideGl, 'Text', 'Click the curve to select a point.', ...
                'VerticalAlignment', 'top');
            obj.infoPanel.Layout.Row = 5;

            obj.statusLine = GuiCommon.statusLabel(mainLayout, 3);
        end
    end

    %%% ===================================================================
    %%% Callbacks -- public (not Access=private): same convention every
    %%% other viewer in this project uses, so tests can call these directly.
    %%% ===================================================================
    methods
        function onScaleModeChanged(obj)
        %%% Redraw the curve for the currently selected scale mode
        %%% >> onScaleModeChanged(obj)
            cla(obj.curveAxes);
            obj.hCrosshair = gobjects(0);

            isLog = strcmp(obj.scaleModeDropdown.Value, 'log');
            if isLog
                obj.yMinField.Value = '-5';
                obj.yMaxField.Value = '0';
            else
                obj.yMinField.Value = '';
                obj.yMaxField.Value = '';
            end

            try
                [hOut, ~] = Vis.plotDispersionMap(obj.singleResult, [], obj.scaleModeDropdown.Value, ...
                    [], 'cartesian', obj.curveAxes);
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
        %%% Apply the kParRange/value range fields to the curve axes
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
        %%% Snap a curve click to the nearest computed kParRange sample --
        %%% AngSweepViewerGUI's onMapClicked with the nearest-angle step
        %%% dropped (only one angle exists here) and "nearest mode" replaced
        %%% by "nearest kParRange sample" (continuous axis, no discrete mode
        %%% list). Fires OnPointSelected for IsofreqModalShapeViewerGUI.
        %%% >> onCurveClicked(obj, src, event)
            pt = GuiCommon.clickPoint(obj.curveAxes, event);
            kPos = pt(1);

            kRange = obj.singleResult.kParRange;
            [~, idx] = min(abs(kRange - kPos));
            kParVal  = kRange(idx);
            angleVal = obj.singleResult.angle;
            velVal   = obj.singleResult.omega0 / kParVal;

            obj.infoPanel.Text = { ...
                sprintf('angle = %.1f deg', angleVal), ...
                sprintf('kParRange = %.4g', kParVal), ...
                sprintf('velocity = %.4f km/s', velVal), ...
                sprintf('|G_{33}| = %.4g', obj.singleResult.absG(idx))};

            obj.showCrosshairAt(angleVal, kParVal);

            if ~isempty(obj.OnPointSelected)
                obj.OnPointSelected(angleVal, kParVal);
            end
        end

        function showCrosshairAt(obj, ~, kParVal)
        %%% Move (or create) the crosshair marker at kParVal on the curve,
        %%% WITHOUT triggering any click/point-lookup logic -- the external
        %%% entry point a linked IsofreqModalShapeViewerGUI (or any other
        %%% caller) uses to mark its own current selection. Public for
        %%% exactly that reason. First argument (angle) accepted-but-unused --
        %%% kept only so this method's call signature matches
        %%% IsofreqAngSweepViewerGUI's own showCrosshairAt(angle, kParVal),
        %%% letting IsofreqCombinedViewerGUI wire either class identically;
        %%% there is only ever one valid angle here (obj.singleResult.angle).
        %%% >> showCrosshairAt(obj, angleVal, kParVal)
            yPos = obj.valueAtKPar(kParVal);
            obj.hCrosshair = GuiCommon.markCrosshair(obj.curveAxes, obj.hCrosshair, kParVal, yPos);
        end
    end

    %%% ===================================================================
    %%% Run-time helpers
    %%% ===================================================================
    methods (Access = private)
        function yPos = valueAtKPar(obj, kParVal)
        %%% Y-value the currently-displayed curve shows at kParVal -- linear
        %%% interpolation onto the same (log-normalised or raw) |G_33| curve
        %%% Vis.plotDispersionMap itself draws.
        %%%
        %%% Vis.scaleGreenData is the ONE definition of that log/linear
        %%% scaling -- shared with the plot function itself, so the crosshair
        %%% cannot drift off the curve it is meant to sit on.
        %%% >> yPos = valueAtKPar(obj, kParVal)
            yData = Vis.scaleGreenData(obj.singleResult.absG, obj.scaleModeDropdown.Value, 'G_{33}');
            yPos = interp1(obj.singleResult.kParRange, yData, kParVal, 'linear', 'extrap');
        end
    end
end
