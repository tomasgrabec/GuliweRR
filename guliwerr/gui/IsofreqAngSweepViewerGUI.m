classdef IsofreqAngSweepViewerGUI < handle
%%% IsofreqAngSweepViewerGUI -- isoFreq angular-sweep viewer window, the
%%% isoFreq counterpart of AngSweepViewerGUI for CalcIsofreqAngSweepLayered
%%% results (angle x kParRange Green |G_33| heatmap at one fixed omega0).
%%%
%%% Unlike AngSweepViewerGUI (three graph types: Ritz-mode/Green's-function
%%% map/classification diagnostic), this family has exactly ONE possible
%%% graph -- the Green |G_33| field -- since CalcIsofreqAngSweepLayered never
%%% runs a Ritz solve at all (see docs/reference/architecture/forwardCalcs.md,
%%% "isoWavevec vs. isoFreq"). So there is no graphTypeDropdown/
%%% criterionDropdown at all -- absent, not disabled, matching the Plate
%%% viewer family's own precedent. Rendering reuses Vis.plotDispersionMap
%%% directly (it already dispatches on CalcIsofreqAngSweepLayered for its 2D
%%% heatmap case, plotType cartesian/polar included), not a new dedicated
%%% plotting function.
%%%
%%% In 'polar' view the map is drawn in Cartesian X/Y (pol2cart of angle and
%%% kParRange): a click is converted back with cart2pol, the crosshair is placed
%%% with pol2cart, and the angle/kParRange range fields do not apply (they are
%%% disabled; the colour range still applies).
%%%
%%% Click-to-select fires OnPointSelected(angle, kParVal) directly from the
%%% clicked point -- simpler than AngSweepViewerGUI's own nearest-angle-
%%% then-nearest-mode search, since there is no discrete mode list to snap
%%% to here (kParRange is a continuous, physically-meaningful axis --
%%% CalcIsofreqSingleLayered's own docblock). Retains the "View single
%%% angle" drill-down to IsofreqSingleLayeredViewerGUI, same
%%% isDrilledIn/returnToSweepView mechanism AngSweepViewerGUI itself uses.
%%%
%%% See docs/reference/architecture/gui_viewers.md for the full isoFreq viewer family.
%%%
%%% GuliweRR v6.0

    properties
        fig
        mapAxes
        scaleModeDropdown
        plotTypeDropdown
        xMinField, xMaxField, yMinField, yMaxField
        cMinField, cMaxField
        infoPanel
        statusLine
        viewSingleAngleBtn

        OnPointSelected = []   % function handle @(angle, kParVal), set by IsofreqModalShapeViewerGUI
    end

    properties (SetAccess = private)
        sweepResult   % CalcIsofreqAngSweepLayered, already calc()'d
        isDrilledIn = false       % true while showing an embedded
                                   % IsofreqSingleLayeredViewerGUI instead of the map --
                                   % same readable-externally reasoning as
                                   % AngSweepViewerGUI.isDrilledIn
    end

    properties (Access = private)
        hCrosshair = gobjects(0)
        hMapPlot   = gobjects(0)   % handle the active plot returned
        currentAngNo = 1          % last-clicked (or default first) angle index --
                                   % drill-down target for viewSingleAngleBtn
        contentParent             % stable uigridlayout slot the map/drill-down
                                   % content is (re)built into, as in AngSweepViewerGUI
        singleAngleViewer = []    % embedded IsofreqSingleLayeredViewerGUI while drilled in, else []
        lastCrosshairAngle = []   % (angleVal, kParVal) last passed to showCrosshairAt,
        lastCrosshairKPar  = []   % kept across drill-down rebuilds, as in AngSweepViewerGUI
        savedAxisRange = []       % struct(xMin,xMax,yMin,yMax,cMin,cMax) field VALUES
                                   % captured right before drilling in
    end

    methods
        function obj = IsofreqAngSweepViewerGUI(sweepResult, parentContainer)
        %%% Construct the isoFreq angular-sweep viewer window
        %%% >> obj = IsofreqAngSweepViewerGUI(sweepResult, parentContainer)
        %%%
        %%% INPUT: 1. sweepResult     .. CalcIsofreqAngSweepLayered, already
        %%%                              calc()'d (obj.lastResult off
        %%%                              StructureControlGUI's Run button, typically)
        %%%  (opt) 2. parentContainer .. a uigridlayout region (or any valid
        %%%                              uigridlayout parent) to build into instead
        %%%                              of a fresh own uifigure -- used by
        %%%                              IsofreqCombinedViewerGUI/
        %%%                              FreqAngSweep2LayeredViewerGUI to embed
        %%%                              this window's content into a shared
        %%%                              figure. Omitted (default): standalone,
        %%%                              owns its own uifigure.
            validateattributes(sweepResult, {'CalcIsofreqAngSweepLayered'}, ...
                {'scalar'}, 'IsofreqAngSweepViewerGUI', 'sweepResult');
            if isempty(sweepResult.absGMap)
                error('IsofreqAngSweepViewerGUI:constructor:notCalculated', ...
                    'sweepResult has not been calc()''d yet -- nothing to view.');
            end
            obj.sweepResult = sweepResult;
            if nargin < 2; parentContainer = []; end
            obj.buildFigure(parentContainer);
            obj.onGraphChanged();
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
                'GuliweRR -- isoFreq Angular Sweep Viewer', [100 100 1000 700], [3 1]);
            mainLayout.RowHeight = {40, '1x', 26};

            %%% top bar -- scale-mode/plot-type selectors + single-angle
            %%% drill-down toggle. Stays put across drill-down/back (unlike
            %%% obj.contentParent's children).
            topBar = uigridlayout(mainLayout, [1 5]);
            topBar.Layout.Row = 1;
            topBar.ColumnWidth = {60, 90, 80, 90, 170};
            topBar.Padding = [0 0 0 0];

            uilabel(topBar, 'Text', 'Scale:');
            obj.scaleModeDropdown = uidropdown(topBar, 'Items', {'log', 'linear'}, ...
                'ValueChangedFcn', @(~,~) obj.onGraphChanged());
            uilabel(topBar, 'Text', 'Plot type:');
            obj.plotTypeDropdown = uidropdown(topBar, 'Items', {'cartesian', 'polar'}, ...
                'ValueChangedFcn', @(~,~) obj.onGraphChanged());

            obj.viewSingleAngleBtn = uibutton(topBar, 'Text', 'View single angle ->', ...
                'ButtonPushedFcn', @(~,~) obj.onToggleSingleAngleView());

            %%% stable content slot -- buildMapContent()/drillDownToAngle() swap
            %%% what's INSIDE it (never the slot itself), same pattern
            %%% AngSweepViewerGUI's own contentParent uses
            obj.contentParent = uigridlayout(mainLayout, [1 1]);
            obj.contentParent.Layout.Row = 2;
            obj.contentParent.Padding = [0 0 0 0];

            obj.buildMapContent();

            obj.statusLine = GuiCommon.statusLabel(mainLayout, 3);
        end

        function buildMapContent(obj)
        %%% (Re)build the normal map + right-hand side panel (axis-range
        %%% controls, info panel) inside obj.contentParent -- called once at
        %%% construction and again by returnToSweepView() after a drill-down/
        %%% back cycle has destroyed the previous content.
        %%% >> buildMapContent(obj)
            [obj.mapAxes, sideGl] = GuiCommon.plotWithSidePanel(obj.contentParent, ...
                {18, 26, 18, 26, 18, 26, '1x'}, ...
                @(src, event) obj.onMapClicked(src, event));
            obj.hCrosshair = gobjects(0);   % previous handle, if any, died with the old axes

            onRangeEdit = @(~,~) obj.applyAxisRange();

            [obj.xMinField, obj.xMaxField] = GuiCommon.rangeFieldPair(sideGl, ...
                'Angle [deg] min/max', num2str(min(obj.sweepResult.angles)), ...
                num2str(max(obj.sweepResult.angles)), onRangeEdit);

            [obj.yMinField, obj.yMaxField] = GuiCommon.rangeFieldPair(sideGl, ...
                'kParRange min/max', num2str(min(obj.sweepResult.kParRange)), ...
                num2str(max(obj.sweepResult.kParRange)), onRangeEdit);

            %%% Colour min/max -- fields start empty; populated from the
            %%% first draw's own CLim by onGraphChanged, same "re-default
            %%% from the freshly-drawn CLim" mechanic AngSweepViewerGUI's own
            %%% colour-range control uses.
            [obj.cMinField, obj.cMaxField] = GuiCommon.rangeFieldPair(sideGl, ...
                'Colour min/max', '', '', onRangeEdit);

            obj.infoPanel = uilabel(sideGl, 'Text', 'Click the map to select a point.', ...
                'VerticalAlignment', 'top');
            obj.infoPanel.Layout.Row = 7;
        end

        function tf = isPolar(obj)
            tf = strcmp(obj.plotTypeDropdown.Value, 'polar');
        end

        function [x, y] = mapXY(obj, angleVal, kParVal)
        %%% Map-axes coordinates of (angle [deg], kPar): pol2cart in polar view
            if obj.isPolar()
                [x, y] = pol2cart(deg2rad(angleVal), kParVal);
            else
                x = angleVal; y = kParVal;
            end
        end
    end

    %%% ===================================================================
    %%% Callbacks -- public (not Access=private): same convention every
    %%% other viewer in this project uses, so tests can call these directly.
    %%% ===================================================================
    methods
        function onGraphChanged(obj)
        %%% Redraw the map for the currently selected scale/plot-type
        %%% >> onGraphChanged(obj)
            cla(obj.mapAxes);
            %%% cla() clears CHILDREN only, not axis-limit modes -- reset
            %%% explicitly so a later draw doesn't inherit a stale CLim from
            %%% a previous log-mode draw (the gotcha AngSweepViewerGUI.
            %%% onGraphTypeChanged documents, recurring here across
            %%% scale-mode switches).
            obj.mapAxes.XLimMode = 'auto';
            obj.mapAxes.YLimMode = 'auto';
            obj.mapAxes.CLimMode = 'auto';
            %%% nor what plotPolarHeatmap's 'axis equal'/'axis off' set
            axis(obj.mapAxes, 'normal');
            obj.mapAxes.Visible = 'on';
            obj.hCrosshair = gobjects(0);

            try
                [hOut, ~] = Vis.plotDispersionMap(obj.sweepResult, [], obj.scaleModeDropdown.Value, ...
                    [], obj.plotTypeDropdown.Value, obj.mapAxes);
                obj.hMapPlot = hOut;
                set(obj.hMapPlot, 'ButtonDownFcn', @(src, event) obj.onMapClicked(src, event));
                obj.cMinField.Value = num2str(obj.mapAxes.CLim(1));
                obj.cMaxField.Value = num2str(obj.mapAxes.CLim(2));
                GuiCommon.setStatusOk(obj.statusLine, '');
            catch err
                obj.hMapPlot = gobjects(0);
                GuiCommon.setStatusError(obj.statusLine, err.message);
            end

            set([obj.xMinField, obj.xMaxField, obj.yMinField, obj.yMaxField], ...
                'Enable', matlab.lang.OnOffSwitchState(~obj.isPolar()));
            obj.applyAxisRange();
        end

        function applyAxisRange(obj)
        %%% Apply the angle/kParRange/colour range fields to the map axes
        %%% >> applyAxisRange(obj)
        %%% In polar view only the colour range: the axes hold pol2cart X/Y, not
        %%% (angle, kParRange).
            if ~obj.isPolar()
                GuiCommon.applyFieldLimit(obj.mapAxes, 'x', obj.xMinField, obj.xMaxField);
                GuiCommon.applyFieldLimit(obj.mapAxes, 'y', obj.yMinField, obj.yMaxField);
            end
            GuiCommon.applyFieldLimit(obj.mapAxes, 'c', obj.cMinField, obj.cMaxField);
        end

        function setAxisRange(obj, axisName, lims)
        %%% Set one axis range from code, as if typed into its min/max fields
        %%% >> setAxisRange(obj, axisName, lims)    e.g. setAxisRange(obj, 'y', [NaN 7])
        %%%
        %%% axisName 'x' | 'y' | 'c'; lims [lo hi], a NaN
        %%% entry keeps that limit, [] keeps both -- see GuiCommon.setAxisRange
            GuiCommon.setAxisRange(obj, obj.mapAxes, axisName, lims);
        end

        function onMapClicked(obj, ~, event)
        %%% Snap a map click to the nearest computed (angle, kParRange) grid
        %%% point -- AngSweepViewerGUI's onMapClicked with "nearest mode
        %%% velocity" replaced by "nearest kParRange sample" (continuous
        %%% axis, no discrete mode list to snap to). Fires OnPointSelected
        %%% for IsofreqModalShapeViewerGUI to attach to.
        %%% >> onMapClicked(obj, src, event)
            pt = GuiCommon.clickPoint(obj.mapAxes, event);
            if obj.isPolar()
                [th, kPos] = cart2pol(pt(1), pt(2));
                %%% nearest angle on the circle: cart2pol returns (-180, 180]
                [~, angNo] = min(abs(mod(obj.sweepResult.angles - rad2deg(th) + 180, 360) - 180));
            else
                kPos = pt(2);
                [~, angNo] = min(abs(obj.sweepResult.angles - pt(1)));
            end
            obj.currentAngNo = angNo;   % drillDownToAngle's target if the user hits viewSingleAngleBtn next
            [~, kIdx] = min(abs(obj.sweepResult.kParRange - kPos));
            kParVal  = obj.sweepResult.kParRange(kIdx);
            angleVal = obj.sweepResult.angles(angNo);
            absGVal  = obj.sweepResult.absGMap(kIdx, angNo);

            obj.infoPanel.Text = { ...
                sprintf('angle = %.1f deg', angleVal), ...
                sprintf('kParRange = %.4g', kParVal), ...
                sprintf('|G_{33}| = %.4g', absGVal)};

            obj.showCrosshairAt(angleVal, kParVal);

            if ~isempty(obj.OnPointSelected)
                obj.OnPointSelected(angleVal, kParVal);
            end
        end

        function showCrosshairAt(obj, angleVal, kParVal)
        %%% Move (or create) the crosshair marker at (angleVal, kParVal)
        %%% WITHOUT triggering any click/point-lookup logic -- the external
        %%% entry point a linked IsofreqModalShapeViewerGUI (or any other
        %%% caller) uses to mark its own current selection on this map.
        %%% Public for exactly that reason.
        %%%
        %%% While drilled into a single angle (isDrilledIn), forwards to the
        %%% nested IsofreqSingleLayeredViewerGUI instead of drawing on the map
        %%% -- same re-drill-if-mismatched-angle behaviour
        %%% AngSweepViewerGUI.showCrosshairAt uses.
        %%% >> showCrosshairAt(obj, angleVal, kParVal)
            obj.lastCrosshairAngle = angleVal;
            obj.lastCrosshairKPar  = kParVal;
            if obj.isDrilledIn
                [~, angNo] = min(abs(obj.sweepResult.angles - angleVal));
                if angNo ~= obj.currentAngNo
                    obj.drillDownToAngle(angNo);
                end
                obj.singleAngleViewer.showCrosshairAt(angleVal, kParVal);
                return
            end
            [x, y] = obj.mapXY(angleVal, kParVal);
            obj.hCrosshair = GuiCommon.markCrosshair(obj.mapAxes, obj.hCrosshair, x, y);
        end

        function onToggleSingleAngleView(obj)
        %%% Toggle between the normal 2D map and a drilled-down
        %%% IsofreqSingleLayeredViewerGUI at obj.currentAngNo (the
        %%% last-clicked angle, or the first angle if nothing has been
        %%% clicked yet).
        %%% >> onToggleSingleAngleView(obj)
            if obj.isDrilledIn
                obj.returnToSweepView();
            else
                obj.savedAxisRange = GuiCommon.captureAxisRange(obj);
                obj.drillDownToAngle(obj.currentAngNo);
            end
        end

        function drillDownToAngle(obj, angNo)
        %%% Replace the map content with an embedded
        %%% IsofreqSingleLayeredViewerGUI for sweepResult.cs{angNo} --
        %%% already calc()'d as part of the sweep itself, so no
        %%% recomputation is needed.
        %%% >> drillDownToAngle(obj, angNo)
            delete(obj.contentParent.Children);
            obj.hCrosshair = gobjects(0);

            singleResult = obj.sweepResult.cs{angNo};
            obj.singleAngleViewer = IsofreqSingleLayeredViewerGUI(singleResult, obj.contentParent);
            obj.singleAngleViewer.OnPointSelected = ...
                @(angle, kParVal) obj.dispatchPointSelected(angle, kParVal);
            obj.singleAngleViewer.show();

            %%% restore the crosshair if the last one shown belongs to THIS
            %%% angle -- see AngSweepViewerGUI.drillDownToAngle's own comment
            if ~isempty(obj.lastCrosshairAngle)
                [~, crosshairAngNo] = min(abs(obj.sweepResult.angles - obj.lastCrosshairAngle));
                if crosshairAngNo == angNo
                    obj.singleAngleViewer.showCrosshairAt(obj.lastCrosshairAngle, obj.lastCrosshairKPar);
                end
            end

            obj.scaleModeDropdown.Enable = 'off';
            obj.plotTypeDropdown.Enable  = 'off';
            obj.viewSingleAngleBtn.Text = '<< Back to angular sweep';
            obj.currentAngNo = angNo;
            obj.isDrilledIn = true;
        end

        function returnToSweepView(obj)
        %%% Tear down the embedded IsofreqSingleLayeredViewerGUI and rebuild
        %%% the normal map content in its place -- restoring the axis-range
        %%% fields onToggleSingleAngleView saved right before drilling in.
        %%% >> returnToSweepView(obj)
            obj.singleAngleViewer = [];
            delete(obj.contentParent.Children);
            obj.buildMapContent();
            obj.onGraphChanged();
            GuiCommon.restoreAxisRange(obj, obj.savedAxisRange);

            %%% restore the crosshair -- see AngSweepViewerGUI.returnToSweepView's
            %%% own comment
            if ~isempty(obj.lastCrosshairAngle)
                [x, y] = obj.mapXY(obj.lastCrosshairAngle, obj.lastCrosshairKPar);
                obj.hCrosshair = GuiCommon.markCrosshair(obj.mapAxes, obj.hCrosshair, x, y);
            end

            obj.scaleModeDropdown.Enable = 'on';
            obj.plotTypeDropdown.Enable  = 'on';
            obj.viewSingleAngleBtn.Text = 'View single angle ->';
            obj.isDrilledIn = false;
        end

        function dispatchPointSelected(obj, angle, kParVal)
        %%% Forward a click inside the drilled-down
        %%% IsofreqSingleLayeredViewerGUI out through THIS object's own
        %%% OnPointSelected -- same guard pattern onMapClicked itself uses.
        %%% >> dispatchPointSelected(obj, angle, kParVal)
            if ~isempty(obj.OnPointSelected)
                obj.OnPointSelected(angle, kParVal);
            end
        end
    end
end
