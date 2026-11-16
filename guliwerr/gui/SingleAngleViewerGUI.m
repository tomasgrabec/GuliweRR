classdef SingleAngleViewerGUI < handle
%%% SingleAngleViewerGUI -- single-angle "curve" viewer window, the 1D analogue
%%% of AngSweepViewerGUI for CalcSingleHS/CalcSingleLayered results, which have
%%% no angle axis to sweep and so no 2D map to draw.
%%%
%%% One uifigure window: a curve (classification curve / Green's-function curve /
%%% Ritz-mode curve / Ritz spectrum curve), a velocity/value axis-range
%%% control, and click-to-select-a-mode interactivity — structurally the role
%%% AngSweepViewerGUI plays for an angular sweep, driven by
%%% Vis.plotClassificationCurve/plotRitzModeCurve/plotGreenCurve/
%%% plotRitzSpectrumCurve instead of their *Map/*Diagnostic counterparts.
%%%
%%% CLASSIFICATION BAR (second top-bar row, a ClassificationBar): method, standard
%%% and Reclassify for this result -- see ClassificationBar. After a Reclassify the
%%% curve is redrawn (refreshClassification) and OnClassificationChanged fires, so a
%%% host (CombinedViewerGUI) can refresh a linked modal-shape viewer. Hidden when
%%% AngSweepViewerGUI embeds this class as its drill-down: the sweep's own bar then
%%% classifies every angle at once, and reclassifying one angle alone would leave the
%%% sweep's maps out of step.
%%%
%%% HOST CONTROLS (4th constructor argument): in that drill-down the sweep viewer's
%%% own Graph type / Criterion dropdowns drive this curve (setView), so this class
%%% hides its own two -- one set of controls on screen, kept in step across the
%%% drill-down and back. The dropdowns still exist (setView sets them), only hidden.
%%%
%%% EXPORT FIGURE (button beside 'Show types:'): GuiCommon.exportAxes copies the
%%% curve as shown into an ordinary figure window, to save in any format.
%%%
%%% The Ritz-mode curve's criterion dropdown: 'Ritz score' (each mode's score of
%%% the Ritz candidate it supports, Vis.ritzModeScores; modes supporting none drawn
%%% faint), 'freqDisp' (the half-space veto quantity, with the veto level
%%% RitzScore.MIN_FREQDISP drawn; a click below it reads "vetoed") and
%%% 'depthInsensitivity'. On a layered result freqDisp is greyed in the list and the
%%% status line notes it is not used on layered stacks.
%%%
%%% The Ritz spectrum curve (Vis.plotRitzSpectrumCurve: |S(v)| from the stored modal
%%% numerators) is the Ritz-side response next to the Green's-function curve, and the
%%% classification curve's background on a 'ritzOnly' result, which has no Green data.
%%%
%%% CLASSIFICATION OVERLAY and SELECTION READOUT: as in AngSweepViewerGUI (see
%%% its docblock). On the Ritz-mode, Green's-function and Ritz spectrum curves the
%%% overlaid classified points are placed at the DISPLAYED curve's own value at
%%% each mode's velocity (valueAtVelocity -- the crosshair's own lookup, which
%%% reads the plotted line itself), not at plotClassificationCurve's projection, so
%%% they sit on whichever curve is shown.
%%%
%%% The Green's-function curve (and the classification curve's background) is
%%% |G13_eff|, the +-k symmetrised TGS response -- the one Green component a
%%% result stores and the one the existence score reads. Other
%%% components are not offered here: Vis.plotGreenCurve(singleResult, [], [],
%%% 'G33') etc. recomputes any of them on request, as in AngSweepViewerGUI.
%%% A 'ritzOnly' result has no Green data; the Green curve then reports the
%%% plotGreenCurve error on the status line, and the classification curve is drawn
%%% on the Ritz spectrum instead.
%%%
%%% Implements the SAME duck-typed interface AngSweepViewerGUI does (.fig,
%%% .show(), .OnModeSelected, .OnClassificationChanged, .showCrosshairAt(angle,
%%% vel)), so
%%% CombinedViewerGUI/ModalShapeViewerGUI's wiring works against either class
%%% without knowing which it has. Angle is always the fixed
%%% obj.singleResult.angle here, accepted as an argument for interface symmetry.
%%%
%%% CROSSHAIR PERSISTENCE: the top-half mode's crosshair (last position passed
%%% to showCrosshairAt, kept in lastCrosshairVel) is redrawn after every
%%% graph-type/criterion change rather than lost to onGraphTypeChanged's cla()
%%% -- the mechanic AngSweepViewerGUI uses across its drill-down toggle, applied
%%% here to the curve redraw. modeCrosshairCheckbox (on by default) lets the
%%% user switch it off, e.g. when the marker sits awkwardly on a busy curve --
%%% see drawOrHideCrosshair.
%%%
%%% Only TWO axis-range fields (velocity, value), not three: none of this
%%% class's curve types use colour as a data dimension by default. The Ritz-mode,
%%% Green's-function and Ritz spectrum curves are plain 2D plots, and the
%%% classification curve's default colorBy='type' uses fixed per-type swatches
%%% rather than a continuous scale, so a colour-range control would control
%%% nothing.
%%%
%%% GuliweRR v6.0

    properties
        fig
        curveAxes
        graphTypeLabel
        graphTypeDropdown
        criterionLabel
        criterionDropdown
        modeCrosshairCheckbox
        xMinField, xMaxField, yMinField, yMaxField
        statusLine
        classBar = []         % ClassificationBar, [] when hidden (drill-down of AngSweepViewerGUI)
        exportBtn             % 'Export figure' -- see class docblock

        OnModeSelected = []   % function handle @(angle, modeNo, velSnapped, velRaw), for a linked modal-shape viewer
        OnClassificationChanged = []   % function handle @(), fired after this viewer's own
                                       % Reclassify, for a linked modal-shape viewer
        showSelectionInStatus   % true: a click's selectionInfo is also written to the status
                                % line. Defaults to standalone (no parentContainer);
                                % AngSweepViewerGUI's drill-down passes its own value on
    end

    properties (SetAccess = private)
        singleResult   % CalcSingleHS or CalcSingleLayered, already calcAll()'d
        hostControls = false   % true: a host drives graph type / criterion (setView)
                               % and this viewer's own two dropdowns are hidden
        selectionInfo = {}   % {angle, mode #, velocity, type, ...} lines of the last click, {}
                             % before any; lines 5 onwards (classification readout,
                             % criterion value) are also shown when embedded
    end

    properties (Access = private)
        hCrosshair = gobjects(0)
        hCurvePlot = gobjects(0)   % handle(s) the active plot function returned
        hCurveBg   = gobjects(0)   % the classification curve's background line (Green,
                                    % or the Ritz spectrum on 'ritzOnly'); gobjects(0)
                                    % for the other graph types
        typeCheckPanel             % 'Show types:' panel, built once, contents rebuilt
                                    % per draw -- same mechanic as AngSweepViewerGUI's own
        typeCheckboxes = gobjects(0)
        diagnosticTypeVisible = struct()   % per-type checkbox states, field = type name --
        overlayTypeVisible    = struct()   % see AngSweepViewerGUI's CLASSIFICATION OVERLAY note
        lastCrosshairVel = []      % velVal last passed to showCrosshairAt, kept across
                                    % onGraphTypeChanged's cla() so drawOrHideCrosshair can
                                    % redraw it on the fresh axes -- [] before any selection
    end

    methods
        function obj = SingleAngleViewerGUI(singleResult, parentContainer, showClassBar, hostControls)
        %%% Construct the single-angle curve viewer window
        %%% >> obj = SingleAngleViewerGUI(singleResult, parentContainer, showClassBar, hostControls)
        %%%
        %%% INPUT: 1. singleResult    .. CalcSingleHS or CalcSingleLayered,
        %%%                              already calcAll()'d (obj.lastResult off
        %%%                              StructureControlGUI's Run button, typically)
        %%%  (opt) 2. parentContainer .. a uigridlayout region (or any valid
        %%%                              uigridlayout parent) to build into instead
        %%%                              of a fresh own uifigure -- used by
        %%%                              CombinedViewerGUI to embed this window's
        %%%                              content into a shared figure. Omitted
        %%%                              (default): standalone, owns its own uifigure.
        %%%  (opt) 3. showClassBar    .. true (default): build the classification bar;
        %%%                              false: AngSweepViewerGUI's drill-down (see
        %%%                              class docblock)
        %%%  (opt) 4. hostControls    .. false (default); true: the host drives the
        %%%                              graph type and criterion through setView and
        %%%                              this viewer hides its own dropdowns (see
        %%%                              class docblock)
            validateattributes(singleResult, {'CalcSingleHS', 'CalcSingleLayered'}, ...
                {'scalar'}, 'SingleAngleViewerGUI', 'singleResult');
            %%% classification stays [] (CalcSingleHS/Layered's own default) until
            %%% identifyModes() runs, called at the end of calcAll() -- same
            %%% "has this actually run" signal AngSweepViewerGUI uses via
            %%% significantModesMap, one level down (per-object instead of
            %%% per-sweep).
            if isempty(singleResult.classification)
                error('SingleAngleViewerGUI:constructor:notCalculated', ...
                    'singleResult has not been calcAll()''d yet -- nothing to view.');
            end
            obj.singleResult = singleResult;
            if nargin < 2; parentContainer = []; end
            if nargin < 3 || isempty(showClassBar); showClassBar = true; end
            if nargin < 4 || isempty(hostControls); hostControls = false; end
            obj.hostControls = logical(hostControls);
            obj.showSelectionInStatus = isempty(parentContainer);
            obj.buildFigure(parentContainer, showClassBar);
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
        function buildFigure(obj, parentContainer, showClassBar)
            [obj.fig, mainLayout] = GuiCommon.hostLayout(parentContainer, ...
                'GuliweRR -- Single-Angle Viewer', [100 100 1000 730], [4 1]);
            if showClassBar
                mainLayout.RowHeight = {40, 30, '1x', 26};
            else
                mainLayout.RowHeight = {40, 0, '1x', 26};
            end
            if obj.hostControls
                %%% only the crosshair toggle is left in the top bar -- see HOST CONTROLS
                mainLayout.RowHeight{1} = 24;
            end

            %%% top bar -- graph-type + criterion selectors, crosshair toggle
            %%% columns: graph-type label/dropdown | criterion label/dropdown |
            %%% spacer | crosshair toggle, kept right-aligned -- the same column
            %%% layout as AngSweepViewerGUI's top bar
            topBar = uigridlayout(mainLayout, [1 6]);
            topBar.Layout.Row = 1;
            topBar.ColumnWidth = {80, 190, 80, 160, '1x', 170};
            topBar.Padding = [0 0 0 0];

            obj.graphTypeLabel = uilabel(topBar, 'Text', 'Graph type:');
            obj.graphTypeLabel.Layout.Column = 1;
            obj.graphTypeDropdown = uidropdown(topBar, ...
                'Items', {'Classification curve', 'Green''s-function curve', 'Ritz-mode curve', ...
                          'Ritz spectrum curve'}, ...
                'ValueChangedFcn', @(~,~) obj.onGraphTypeChanged());
            obj.graphTypeDropdown.Layout.Column = 2;
            if obj.hostControls   % see class docblock's HOST CONTROLS note
                obj.graphTypeLabel.Visible = 'off';
                obj.graphTypeDropdown.Visible = 'off';
            end

            obj.criterionLabel = uilabel(topBar, 'Text', 'Criterion:', 'Visible', 'off');
            obj.criterionLabel.Layout.Column = 3;
            obj.criterionDropdown = GuiCommon.ritzCriterionDropdown(topBar, ...
                isa(obj.singleResult, 'CalcSingleLayered'), @(~,~) obj.onGraphTypeChanged());
            obj.criterionDropdown.Layout.Column = 4;

            %%% on by default, like DispersionViewerGUI's "Show ..." checkboxes --
            %%% see class docblock's CROSSHAIR PERSISTENCE note
            obj.modeCrosshairCheckbox = uicheckbox(topBar, 'Text', 'Mode crosshair', ...
                'Value', true, 'ValueChangedFcn', @(~,~) obj.drawOrHideCrosshair());
            obj.modeCrosshairCheckbox.Layout.Column = 6;

            %%% classification bar -- see class docblock
            if showClassBar
                obj.classBar = ClassificationBar(mainLayout, obj.singleResult, 2);
                obj.classBar.setReclassifyHandler(@() obj.onReclassified());
            end

            %%% curve + right-hand side panel (axis-range controls, type checkboxes)
            [obj.curveAxes, sideGl, content] = GuiCommon.plotWithSidePanel(mainLayout, ...
                {18, 26, 18, 26, 10, 24, '1x'}, ...
                @(src, event) obj.onCurveClicked(src, event));
            content.Layout.Row = 3;

            %%% default to the Green velocity axis' own range where available
            %%% (calcMode='full'), but starting at 0.9 x the slowest Ritz mode: the
            %%% axis itself starts at 0.3 km/s (SawSpectrum.VEL_MIN), far below any
            %%% surface wave. Falls back to the Ritz spectrum's own vel range under
            %%% 'ritzOnly' (no greenVelAxis at all then).
            if ~isempty(obj.singleResult.greenVelAxis)
                xMinDefault = num2str(max(min(obj.singleResult.greenVelAxis), ...
                    0.9 * min(obj.singleResult.ritzSolver.vels)));
                xMaxDefault = num2str(max(obj.singleResult.greenVelAxis));
            else
                ritzVels = obj.singleResult.ritzSolver.vels(1:obj.singleResult.ritzSolver.nCompModes);
                xMinDefault = num2str(min(ritzVels));
                xMaxDefault = num2str(max(ritzVels));
            end
            onRangeEdit = @(~,~) obj.applyAxisRange();
            [obj.xMinField, obj.xMaxField] = GuiCommon.rangeFieldPair(sideGl, ...
                'Velocity [km/s] min/max', xMinDefault, xMaxDefault, onRangeEdit);

            %%% Value (y) range -- left blank/auto by default, unlike the velocity
            %%% range above: the y-axis meaning changes per graph type (criterion
            %%% value / |G| / log-normalised |G|), so no single sensible
            %%% default range covers all three the way the velocity axis does.
            [obj.yMinField, obj.yMaxField] = GuiCommon.rangeFieldPair(sideGl, ...
                'Value [Y] min/max', '', '', onRangeEdit);

            %%% 'Show types:' header with the Export figure button, then the
            %%% checkbox panel -- built here empty, contents filled in by
            %%% GuiCommon.rebuildTypeCheckboxes on every onGraphTypeChanged draw
            %%% (see class docblock's CLASSIFICATION OVERLAY note). Takes all the
            %%% remaining side-panel height.
            [~, obj.exportBtn] = GuiCommon.typesHeaderRow(sideGl, 6, @() obj.onExportFigure());
            obj.typeCheckPanel = uipanel(sideGl, 'BorderType', 'none');
            obj.typeCheckPanel.Layout.Row = 7;

            obj.statusLine = GuiCommon.statusLabel(mainLayout, 4);
        end
    end

    %%% ===================================================================
    %%% Callbacks -- public (not Access=private): same convention
    %%% AngSweepViewerGUI/StructureControlGUI use, so tests can call these
    %%% directly without going through actual UI events.
    %%% ===================================================================
    methods
        function onGraphTypeChanged(obj)
        %%% Redraw the curve for whichever graph type is now selected
        %%% >> onGraphTypeChanged(obj)
            isRitz = strcmp(obj.graphTypeDropdown.Value, 'Ritz-mode curve');
            obj.criterionLabel.Visible = isRitz && ~obj.hostControls;
            obj.criterionDropdown.Visible = isRitz && ~obj.hostControls;

            GuiCommon.clearAxes(obj.curveAxes);
            obj.hCrosshair = gobjects(0);
            obj.hCurveBg = gobjects(0);

            try
                hBg = gobjects(0);   % background object(s) drawn UNDER hOut, if any -- see below
                switch obj.graphTypeDropdown.Value
                    case 'Ritz-mode curve'
                        [hOut, ~] = Vis.plotRitzModeCurve(obj.singleResult, obj.criterionDropdown.Value, [], obj.curveAxes);
                    case 'Green''s-function curve'
                        %%% 'ritzOnly' (no Green data): plotGreenCurve errors into
                        %%% the status line
                        [hOut, ~] = Vis.plotGreenCurve(obj.singleResult, [], obj.curveAxes);
                    case 'Ritz spectrum curve'
                        [hOut, ~] = Vis.plotRitzSpectrumCurve(obj.singleResult, [], obj.curveAxes, true);
                    otherwise   % 'Classification curve'
                        %%% hBg is the background curve (Green, or the Ritz spectrum on
                        %%% 'ritzOnly'), BEHIND hOut's scatter -- it gets the same
                        %%% ButtonDownFcn below so a click landing on the curve itself
                        %%% (not a marker) still registers, as on AngSweepViewerGUI's
                        %%% map (where the background is a full-axes pcolor)
                        [hOut, ~, hBg] = Vis.plotClassificationCurve(obj.singleResult, [], [], [], obj.curveAxes);
                        if ~isempty(hBg); obj.hCurveBg = hBg; end
                end
                obj.hCurvePlot = hOut;
                set([obj.hCurvePlot(:); hBg(:)], 'ButtonDownFcn', @(src, event) obj.onCurveClicked(src, event));
                %%% type checkboxes -- see class docblock's CLASSIFICATION OVERLAY note
                obj.typeCheckboxes = GuiCommon.classificationTypeCheckboxes(obj.typeCheckPanel, ...
                    strcmp(obj.graphTypeDropdown.Value, 'Classification curve'), hOut, ...
                    @() obj.drawClassificationOverlay(), obj.diagnosticTypeVisible, ...
                    obj.overlayTypeVisible, @(stateProp, typeName, isOn) obj.rememberTypeVisible(stateProp, typeName, isOn));
                GuiCommon.setStatusOk(obj.statusLine, obj.idleStatusText());
            catch err
                obj.hCurvePlot = gobjects(0);
                obj.hCurveBg = gobjects(0);
                obj.typeCheckboxes = GuiCommon.rebuildTypeCheckboxes(obj.typeCheckPanel, gobjects(0));
                GuiCommon.setStatusError(obj.statusLine, err.message);
            end

            obj.applyAxisRange();
            %%% redraw the persisted crosshair on the fresh axes -- last, since
            %%% valueAtVelocity reads the lines drawn above (see class docblock's
            %%% CROSSHAIR PERSISTENCE note)
            obj.drawOrHideCrosshair();
        end

        function setView(obj, graphType, criterion)
        %%% Show graphType (one of graphTypeDropdown.Items) with criterion (the
        %%% Ritz-mode curve's; [] keeps the current one) and redraw -- the entry
        %%% point of a host that drives this viewer (see HOST CONTROLS)
        %%% >> setView(obj, graphType, criterion)
            obj.graphTypeDropdown.Value = graphType;
            if nargin > 2 && ~isempty(criterion)
                obj.criterionDropdown.Value = criterion;
            end
            obj.onGraphTypeChanged();
        end

        function fig = onExportFigure(obj)
        %%% Export figure button: the curve as shown, in an ordinary figure window
        %%% >> fig = onExportFigure(obj)
            fig = GuiCommon.exportAxes(obj.curveAxes);
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
        %%% Snap a curve click to the nearest computed Ritz mode by velocity --
        %%% AngSweepViewerGUI's onMapClicked with the nearest-angle step dropped
        %%% (only one angle exists here). Fires OnModeSelected.
        %%% >> onCurveClicked(obj, src, event)
            pt = GuiCommon.clickPoint(obj.curveAxes, event);
            velPos = pt(1);

            ritz  = obj.singleResult.ritzSolver;
            nM    = ritz.nCompModes;
            vels  = ritz.vels(1:nM);
            [~, modeNo] = min(abs(vels - velPos));
            velSel   = vels(modeNo);
            angleVal = obj.singleResult.angle;

            sm = obj.singleResult.classification.significantModes;
            typeStr = GuiCommon.significantModeType(sm, modeNo);

            obj.selectionInfo = [GuiCommon.selectionLines(angleVal, modeNo, velSel, typeStr), ...
                GuiCommon.selectionExtras(obj.singleResult, sm, modeNo, obj.activeCriterion())];
            GuiCommon.setStatusOk(obj.statusLine, obj.idleStatusText());

            obj.showCrosshairAt(angleVal, velSel);

            if ~isempty(obj.OnModeSelected)
                obj.OnModeSelected(angleVal, modeNo, velSel, velPos);
            end
        end

        function showCrosshairAt(obj, ~, velVal)
        %%% Move (or create) the crosshair marker at velVal on whichever curve is
        %%% currently displayed, WITHOUT triggering any click/mode-lookup logic --
        %%% the entry point a linked modal-shape viewer (or any other caller)
        %%% uses to mark its own selection, and public for that reason. The angle
        %%% argument is unused, kept so the signature matches
        %%% AngSweepViewerGUI.showCrosshairAt(angle, vel). velVal is stored in
        %%% lastCrosshairVel so later redraws restore it (see class docblock).
        %%% >> showCrosshairAt(obj, angleVal, velVal)
            obj.lastCrosshairVel = velVal;
            obj.drawOrHideCrosshair();
        end

        function rememberTypeVisible(obj, stateProp, typeName, isOn)
        %%% A type checkbox's onToggle -- store its new state so the next draw
        %%% of the same context (classification curve or overlay) restores it
        %%% >> rememberTypeVisible(obj, stateProp, typeName, isOn)
            obj.(stateProp).(typeName) = isOn;
        end

        function refreshClassification(obj)
        %%% Redraw after the result's classification changed (a Reclassify here or in
        %%% a host): the curve and its overlay, the bar, and the selection readout
        %%% >> refreshClassification(obj)
            if ~isempty(obj.classBar)
                obj.classBar.refresh();
            end
            obj.redrawClassified();
        end

        function onReclassified(obj)
        %%% The classification bar's reclassify handler: redraw, then notify a host
        %%% >> onReclassified(obj)
            obj.redrawClassified();
            if ~isempty(obj.OnClassificationChanged)
                obj.OnClassificationChanged();
            end
        end

        function drawOrHideCrosshair(obj)
        %%% Draw the crosshair at lastCrosshairVel on the CURRENTLY displayed
        %%% curve, or remove it -- whichever modeCrosshairCheckbox and
        %%% lastCrosshairVel currently call for. The single place both
        %%% showCrosshairAt (new selection) and onGraphTypeChanged (redraw
        %%% after a graph-type/criterion switch) go through, so the two stay
        %%% in sync without duplicating the on/off logic.
        %%% >> drawOrHideCrosshair(obj)
            %%% isempty(hCurvePlot) also covers a failed draw (onGraphTypeChanged's
            %%% catch sets it to gobjects(0)) -- nothing to place the marker AGAINST
            %%% then, e.g. a Green curve request under calcMode='ritzOnly'
            if isempty(obj.lastCrosshairVel) || ~obj.modeCrosshairCheckbox.Value || isempty(obj.hCurvePlot)
                if ~isempty(obj.hCrosshair) && isgraphics(obj.hCrosshair)
                    delete(obj.hCrosshair);
                end
                obj.hCrosshair = gobjects(0);
                return
            end
            yPos = obj.valueAtVelocity(obj.lastCrosshairVel);
            obj.hCrosshair = GuiCommon.markCrosshair(obj.curveAxes, obj.hCrosshair, obj.lastCrosshairVel, yPos);
        end
    end

    %%% ===================================================================
    %%% Run-time helpers
    %%% ===================================================================
    methods (Access = private)
        function redrawClassified(obj)
        %%% Redraw the curve and its overlay with the selection readout cleared --
        %%% the shared part of refreshClassification and onReclassified
        %%% >> redrawClassified(obj)
            obj.selectionInfo = {};
            obj.onGraphTypeChanged();
        end

        function hOverlay = drawClassificationOverlay(obj)
        %%% Draw the classified points on top of the Ritz-mode/Green's-function
        %%% curve already on curveAxes -- plotClassificationCurve itself, without
        %%% its Green background, restyled as an overlay, with each point's y
        %%% moved onto the DISPLAYED curve (see class docblock)
        %%% >> hOverlay = drawClassificationOverlay(obj)
            savedLabels = GuiCommon.captureAxesLabels(obj.curveAxes);
            [hOverlay, ~, ~, hLim] = Vis.plotClassificationCurve( ...
                obj.singleResult, 'type', false, [], obj.curveAxes);
            for k = 1:numel(hOverlay)
                hOverlay(k).YData = arrayfun(@(v) obj.valueAtVelocity(v), hOverlay(k).XData);
            end
            GuiCommon.styleClassificationOverlay(obj.curveAxes, hOverlay, hLim, savedLabels);
        end

        function msg = idleStatusText(obj)
        %%% Status-line text after a successful draw or click (GuiCommon.statusText)
        %%% >> msg = idleStatusText(obj)
            msg = GuiCommon.statusText(obj.selectionInfo, obj.showSelectionInStatus, ...
                obj.activeCriterion(), isa(obj.singleResult, 'CalcSingleLayered'));
        end

        function crit = activeCriterion(obj)
        %%% The Ritz-mode criterion while the Ritz-mode curve is shown, else ''
        %%% >> crit = activeCriterion(obj)
            crit = '';
            if strcmp(obj.graphTypeDropdown.Value, 'Ritz-mode curve')
                crit = obj.criterionDropdown.Value;
            end
        end

        function yPos = valueAtVelocity(obj, velVal)
        %%% Y-value the currently-displayed curve shows at velVal, read off the
        %%% plotted line itself, so the marker cannot drift off the curve it is
        %%% meant to sit on --
        %%% Ritz-mode curve: the nearest PLOTTED mode's own value (the curve is a
        %%% discrete Ritz spectrum, not a continuous function, so "nearest" is the
        %%% only sensible interpolation; read off every plotted series, so the
        %%% 'ritzScore' faint baseline is matched too); Ritz spectrum curve: the
        %%% nearest sample; Green's-function curve: linear interpolation along the
        %%% line; classification curve: linear interpolation along its background
        %%% line (Green, or the Ritz spectrum on a 'ritzOnly' result). A
        %%% classification curve without a background has nothing to place the
        %%% marker on: SingleAngleViewerGUI:valueAtVelocity:noCurve.
        %%% >> yPos = valueAtVelocity(obj, velVal)
            switch obj.graphTypeDropdown.Value
                case {'Ritz-mode curve', 'Ritz spectrum curve'}
                    [~, idx, ys] = nearestSample(obj.hCurvePlot, velVal);
                    yPos = ys(idx);
                case 'Green''s-function curve'
                    yPos = interp1(obj.hCurvePlot.XData, obj.hCurvePlot.YData, velVal, 'linear', 'extrap');
                otherwise   % 'Classification curve'
                    if isempty(obj.hCurveBg) || ~isgraphics(obj.hCurveBg)
                        error('SingleAngleViewerGUI:valueAtVelocity:noCurve', ...
                            'The classification curve has no background curve to place a marker on.');
                    end
                    yPos = interp1(obj.hCurveBg.XData, obj.hCurveBg.YData, velVal, 'linear', 'extrap');
            end
        end
    end
end

function [xs, idx, ys] = nearestSample(hLines, velVal)
%%% All plotted samples of hLines (one or more line handles) and the index of the
%%% one nearest velVal in x
    xs = []; ys = [];
    for h = hLines(:)'
        xs = [xs; h.XData(:)]; ys = [ys; h.YData(:)]; %#ok<AGROW> two handles at most
    end
    [~, idx] = min(abs(xs - velVal));
end
