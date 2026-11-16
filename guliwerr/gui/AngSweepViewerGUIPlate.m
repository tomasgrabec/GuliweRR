classdef AngSweepViewerGUIPlate < handle
%%% AngSweepViewerGUIPlate -- angular-sweep viewer window for a
%%% CalcAngSweepPlate result. Plate counterpart of AngSweepViewerGUI.
%%%
%%% A SEPARATE class from AngSweepViewerGUI, not a class-dispatch branch inside
%%% it — the same parallel-sibling convention every Plate GUI class uses.
%%% Differences from AngSweepViewerGUI, all following from the plate physics:
%%%   - No graphTypeDropdown at all — ABSENT, not disabled. There is exactly one
%%%     map here, the Ritz-mode map: no Green's-function map (no analogue
%%%     exists) and no classification diagnostic (no significant-vs-
%%%     insignificant distinction to classify INTO).
%%%   - criterionDropdown offers 'Mode type' | 'Confidence'
%%%     (Vis.plotRitzModeMapPlate's 'modeType'/'confidence'), in place of
%%%     'Ritz score'/'freqDisp'/'depthInsensitivity' (all SAW-specific).
%%%   - Colour-range fields (cMinField/cMaxField) are meaningful only for
%%%     criterion='confidence', a genuine continuous scale whose CLim is
%%%     re-defaulted from the just-drawn axes. They are disabled, not hidden,
%%%     for criterion='modeType', whose scatter points carry explicit RGB
%%%     colour — a 4-way categorical legend, not a CLim-mapped scale.
%%%   - No 'Show types:' per-type checkbox filter: plotRitzModeMapPlate returns
%%%     ONE scatter handle for 'modeType', colour-coded via a colour matrix
%%%     rather than one handle per type, so there is nothing to show/hide
%%%     individually.
%%%   - onMapClicked's info panel reports mode type and confidence directly off
%%%     RitzSolverLamb.modesType/modeTypeConfidence, not a significantModesMap
%%%     lookup.
%%%
%%% Single-angle drill-down (viewSingleAngleBtn) works exactly as
%%% AngSweepViewerGUI's does: sweepResult.cs{angNo} is already a fully calc()'d
%%% CalcSinglePlate, so drilling into an embedded SingleAngleViewerGUIPlate
%%% needs no recomputation.
%%%
%%% GuliweRR v6.0

    properties
        fig
        mapAxes
        criterionLabel
        criterionDropdown
        xMinField, xMaxField, yMinField, yMaxField
        cMinField, cMaxField
        infoPanel
        statusLine
        viewSingleAngleBtn

        OnModeSelected = []   % function handle @(angle, modeNo, velSnapped, velRaw), for a linked modal-shape viewer
    end

    properties (SetAccess = private)
        sweepResult   % CalcAngSweepPlate, already calc()'d
        isDrilledIn = false
    end

    properties (Access = private)
        hCrosshair = gobjects(0)
        hMapPlot   = gobjects(0)
        currentAngNo = 1
        contentParent
        singleAngleViewer = []
        lastCrosshairAngle = []   % (angleVal, velVal) last passed to showCrosshairAt,
        lastCrosshairVel   = []   % kept across drill-down rebuilds, as in AngSweepViewerGUI
        savedAxisRange = []
    end

    methods
        function obj = AngSweepViewerGUIPlate(sweepResult, parentContainer)
        %%% Construct the Plate angular-sweep viewer window
        %%% >> obj = AngSweepViewerGUIPlate(sweepResult, parentContainer)
        %%%
        %%% INPUT: 1. sweepResult     .. CalcAngSweepPlate, already calc()'d
        %%%  (opt) 2. parentContainer .. a uigridlayout region to build into
        %%%                              instead of a fresh own uifigure --
        %%%                              used by CombinedViewerGUIPlate
            validateattributes(sweepResult, {'CalcAngSweepPlate'}, ...
                {'scalar'}, 'AngSweepViewerGUIPlate', 'sweepResult');
            %%% .cs is populated by the CONSTRUCTOR itself (uncalculated
            %%% CalcSinglePlate objects, one per angle) -- not a valid "has
            %%% calc() run" signal, same pitfall AngSweepViewerGUI's
            %%% constructor flags for significantModesMap. velsMap stays []
            %%% until assembleMaps() runs at the end of calc().
            if isempty(sweepResult.velsMap)
                error('AngSweepViewerGUIPlate:constructor:notCalculated', ...
                    'sweepResult has not been calc()''d yet -- nothing to view.');
            end
            obj.sweepResult = sweepResult;
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
                'GuliweRR -- Angular Sweep Viewer (Plate)', [100 100 1000 700], [3 1]);
            mainLayout.RowHeight = {40, '1x', 26};

            %%% top bar -- criterion selector + single-angle drill-down toggle
            %%% (no graph-type dropdown, see class docblock)
            topBar = uigridlayout(mainLayout, [1 3]);
            topBar.Layout.Row = 1;
            topBar.ColumnWidth = {80, 160, 170};
            topBar.Padding = [0 0 0 0];

            obj.criterionLabel = uilabel(topBar, 'Text', 'Criterion:');
            obj.criterionDropdown = uidropdown(topBar, ...
                'Items', {'Mode type', 'Confidence'}, 'ItemsData', {'modeType', 'confidence'}, ...
                'ValueChangedFcn', @(~,~) obj.onGraphTypeChanged());

            obj.viewSingleAngleBtn = uibutton(topBar, 'Text', 'View single angle ->', ...
                'ButtonPushedFcn', @(~,~) obj.onToggleSingleAngleView());

            obj.contentParent = uigridlayout(mainLayout, [1 1]);
            obj.contentParent.Layout.Row = 2;
            obj.contentParent.Padding = [0 0 0 0];

            obj.buildMapContent();

            obj.statusLine = GuiCommon.statusLabel(mainLayout, 3);
        end

        function buildMapContent(obj)
        %%% (Re)build the normal map + right-hand side panel inside
        %%% obj.contentParent -- see AngSweepViewerGUI.buildMapContent's own
        %%% docblock for why this is a separate rebuildable method.
        %%% >> buildMapContent(obj)
            [obj.mapAxes, sideGl] = GuiCommon.plotWithSidePanel(obj.contentParent, ...
                {18, 26, 18, 26, 18, 26, '1x'}, ...
                @(src, event) obj.onMapClicked(src, event));
            obj.hCrosshair = gobjects(0);

            onRangeEdit = @(~,~) obj.applyAxisRange();

            [obj.xMinField, obj.xMaxField] = GuiCommon.rangeFieldPair(sideGl, ...
                'Angle [deg] min/max', num2str(min(obj.sweepResult.angles)), ...
                num2str(max(obj.sweepResult.angles)), onRangeEdit);

            %%% velocity range left blank here and filled from the drawn map
            %%% by onGraphTypeChanged -- no Green velocity axis to default to
            %%% (see class docblock), and the whole Ritz spectrum reaches
            %%% hundreds of km/s, far above the map's own velocity cutoff
            [obj.yMinField, obj.yMaxField] = GuiCommon.rangeFieldPair(sideGl, ...
                'Velocity [km/s] min/max', '', '', onRangeEdit);

            %%% colour min/max -- only meaningful for criterion='confidence',
            %%% see class docblock; re-defaulted from the just-drawn CLim by
            %%% onGraphTypeChanged, same mechanic AngSweepViewerGUI uses
            [obj.cMinField, obj.cMaxField] = GuiCommon.rangeFieldPair(sideGl, ...
                'Colour min/max', '', '', onRangeEdit);

            obj.infoPanel = uilabel(sideGl, 'Text', 'Click the map to select a mode.', ...
                'VerticalAlignment', 'top');
            obj.infoPanel.Layout.Row = 7;
        end
    end

    %%% ===================================================================
    %%% Callbacks -- public (not Access=private), matching
    %%% AngSweepViewerGUI's own convention.
    %%% ===================================================================
    methods
        function onGraphTypeChanged(obj)
        %%% Redraw the map for whichever criterion is now selected
        %%% >> onGraphTypeChanged(obj)
            cla(obj.mapAxes);
            obj.mapAxes.Color = [1 1 1];
            obj.mapAxes.CLimMode = 'auto';
            obj.hCrosshair = gobjects(0);

            isConfidence = strcmp(obj.criterionDropdown.Value, 'confidence');
            try
                [hOut, ~] = Vis.plotRitzModeMapPlate(obj.sweepResult, obj.criterionDropdown.Value, [], obj.mapAxes);
                obj.hMapPlot = hOut;
                set(obj.hMapPlot, 'ButtonDownFcn', @(src, event) obj.onMapClicked(src, event));
                if isempty(obj.yMinField.Value)   % first draw: the map's own range
                    obj.yMinField.Value = num2str(obj.mapAxes.YLim(1));
                    obj.yMaxField.Value = num2str(obj.mapAxes.YLim(2));
                end
                if isConfidence
                    obj.cMinField.Value = num2str(obj.mapAxes.CLim(1));
                    obj.cMaxField.Value = num2str(obj.mapAxes.CLim(2));
                end
                obj.cMinField.Enable = isConfidence;
                obj.cMaxField.Enable = isConfidence;
                GuiCommon.setStatusOk(obj.statusLine, '');
            catch err
                obj.hMapPlot = gobjects(0);
                GuiCommon.setStatusError(obj.statusLine, err.message);
            end

            obj.applyAxisRange();
        end

        function applyAxisRange(obj)
        %%% Apply the angle/velocity/colour range fields to the map axes --
        %%% the colour pair only while it is enabled at all (criterion=
        %%% 'confidence'; see class docblock), unlike AngSweepViewerGUI where
        %%% every graph type has a genuine colour scale to control.
        %%% >> applyAxisRange(obj)
            GuiCommon.applyFieldLimit(obj.mapAxes, 'x', obj.xMinField, obj.xMaxField);
            GuiCommon.applyFieldLimit(obj.mapAxes, 'y', obj.yMinField, obj.yMaxField);
            if strcmp(obj.cMinField.Enable, 'on')
                GuiCommon.applyFieldLimit(obj.mapAxes, 'c', obj.cMinField, obj.cMaxField);
            end
        end

        function setAxisRange(obj, axisName, lims)
        %%% Set one axis range from code, as if typed into its min/max fields
        %%% >> setAxisRange(obj, axisName, lims)    e.g. setAxisRange(obj, 'y', [NaN 7])
        %%%
        %%% axisName 'x' | 'y' | 'c' (only while the colour pair is enabled); lims [lo hi], a NaN
        %%% entry keeps that limit, [] keeps both -- see GuiCommon.setAxisRange
            GuiCommon.setAxisRange(obj, obj.mapAxes, axisName, lims);
        end

        function onMapClicked(obj, ~, event)
        %%% Snap a map click to the nearest computed Ritz mode (nearest angle,
        %%% then nearest mode velocity at that angle). Fires OnModeSelected.
        %%% >> onMapClicked(obj, src, event)
            pt = GuiCommon.clickPoint(obj.mapAxes, event);
            angPos = pt(1); velPos = pt(2);

            [~, angNo] = min(abs(obj.sweepResult.angles - angPos));
            obj.currentAngNo = angNo;
            ritz  = obj.sweepResult.cs{angNo}.ritzSolver;
            nM    = ritz.nCompModes;
            vels  = ritz.vels(1:nM);
            [~, modeNo] = min(abs(vels - velPos));
            velSel = vels(modeNo);
            angleVal = obj.sweepResult.angles(angNo);

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

        function showCrosshairAt(obj, angleVal, velVal)
        %%% Move (or create) the crosshair marker at (angleVal, velVal)
        %%% WITHOUT triggering any click/mode-lookup logic -- same
        %%% drilled-in re-drill behaviour as AngSweepViewerGUI's own
        %%% showCrosshairAt (see that method's own docblock).
        %%% >> showCrosshairAt(obj, angleVal, velVal)
            obj.lastCrosshairAngle = angleVal;
            obj.lastCrosshairVel   = velVal;
            if obj.isDrilledIn
                [~, angNo] = min(abs(obj.sweepResult.angles - angleVal));
                if angNo ~= obj.currentAngNo
                    obj.drillDownToAngle(angNo);
                end
                obj.singleAngleViewer.showCrosshairAt(angleVal, velVal);
                return
            end
            obj.hCrosshair = GuiCommon.markCrosshair(obj.mapAxes, obj.hCrosshair, angleVal, velVal);
        end

        function onToggleSingleAngleView(obj)
        %%% Toggle between the normal 2D map and a drilled-down
        %%% SingleAngleViewerGUIPlate at obj.currentAngNo.
        %%% >> onToggleSingleAngleView(obj)
            if obj.isDrilledIn
                obj.returnToSweepView();
            else
                obj.savedAxisRange = GuiCommon.captureAxisRange(obj);
                obj.drillDownToAngle(obj.currentAngNo);
            end
        end

        function drillDownToAngle(obj, angNo)
        %%% Replace the map content with an embedded SingleAngleViewerGUIPlate
        %%% for sweepResult.cs{angNo} -- already calc()'d as part of the
        %%% sweep itself, so no recomputation is needed.
        %%% >> drillDownToAngle(obj, angNo)
            delete(obj.contentParent.Children);
            obj.hCrosshair = gobjects(0);

            singleResult = obj.sweepResult.cs{angNo};
            obj.singleAngleViewer = SingleAngleViewerGUIPlate(singleResult, obj.contentParent);
            obj.singleAngleViewer.OnModeSelected = ...
                @(angle, modeNo, velSnap, velRaw) obj.dispatchModeSelected(angle, modeNo, velSnap, velRaw);
            obj.singleAngleViewer.show();

            %%% restore the crosshair if the last one shown belongs to THIS
            %%% angle -- see AngSweepViewerGUI.drillDownToAngle's own comment
            if ~isempty(obj.lastCrosshairAngle)
                [~, crosshairAngNo] = min(abs(obj.sweepResult.angles - obj.lastCrosshairAngle));
                if crosshairAngNo == angNo
                    obj.singleAngleViewer.showCrosshairAt(obj.lastCrosshairAngle, obj.lastCrosshairVel);
                end
            end

            obj.criterionDropdown.Enable = 'off';
            obj.viewSingleAngleBtn.Text = '<< Back to angular sweep';
            obj.currentAngNo = angNo;
            obj.isDrilledIn = true;
        end

        function returnToSweepView(obj)
        %%% Tear down the embedded SingleAngleViewerGUIPlate and rebuild the
        %%% normal map content in its place, restoring the axis-range fields
        %%% saved right before drilling in.
        %%% >> returnToSweepView(obj)
            obj.singleAngleViewer = [];
            delete(obj.contentParent.Children);
            obj.buildMapContent();
            obj.onGraphTypeChanged();
            GuiCommon.restoreAxisRange(obj, obj.savedAxisRange);

            %%% restore the crosshair -- see AngSweepViewerGUI.returnToSweepView's
            %%% own comment
            if ~isempty(obj.lastCrosshairAngle)
                obj.hCrosshair = GuiCommon.markCrosshair(obj.mapAxes, obj.hCrosshair, ...
                    obj.lastCrosshairAngle, obj.lastCrosshairVel);
            end

            obj.criterionDropdown.Enable = 'on';
            obj.viewSingleAngleBtn.Text = 'View single angle ->';
            obj.isDrilledIn = false;
        end

        function dispatchModeSelected(obj, angle, modeNo, velSnap, velRaw)
        %%% Forward a click inside the drilled-down SingleAngleViewerGUIPlate
        %%% out through THIS object's own OnModeSelected.
        %%% >> dispatchModeSelected(obj, angle, modeNo, velSnap, velRaw)
            if ~isempty(obj.OnModeSelected)
                obj.OnModeSelected(angle, modeNo, velSnap, velRaw);
            end
        end
    end
end
