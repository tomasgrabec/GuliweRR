classdef AngSweepViewerGUI < handle
%%% AngSweepViewerGUI -- angular-sweep viewer window
%%%
%%% One uifigure window: a map (classification diagnostic / Green's-function map /
%%% Ritz-mode map / Ritz spectrum map), angle/velocity axis-range controls,
%%% and click-to-select-a-mode interactivity. Built on CalcAngSweepHS/
%%% CalcAngSweepLayered results and the Vis plotting functions
%%% (plotClassificationDiagnostic/plotRitzModeMap/plotGreenMap/plotRitzSpectrumMap),
%%% reused directly rather than re-implemented.
%%%
%%% CLASSIFICATION BAR (second top-bar row, a ClassificationBar): method (Green |
%%% Refined | Ritz; layered Green | Ritz), standard and Reclassify for every angle of
%%% the sweep -- see ClassificationBar. After a Reclassify the bar calls
%%% onReclassified, which redraws the map (or the drilled-in single-angle curve)
%%% with the axis ranges kept (redrawClassified) and then fires
%%% OnClassificationChanged, so a host (CombinedViewerGUI, DoubleSweepViewerGUI)
%%% refreshes its linked ModalShapeViewerGUI. refreshClassification is the entry
%%% point for a host that reclassified the sweep elsewhere: it refreshes the bar
%%% and redraws the same way. The nested single-angle viewer has no bar of its
%%% own (see SingleAngleViewerGUI).
%%%
%%% The Ritz-mode map's criterion dropdown: 'Ritz score' (each mode coloured by
%%% the score of the Ritz candidate it supports, Vis.ritzModeScores; modes
%%% supporting none drawn faint), 'freqDisp' (the half-space veto quantity; the
%%% colour scale marks the veto level RitzScore.MIN_FREQDISP and a click below it
%%% reads "vetoed") and 'depthInsensitivity'. On a layered result freqDisp is
%%% greyed in the list and the status line notes it is not used on layered stacks.
%%%
%%% The Ritz spectrum map (Vis.plotRitzSpectrumMap: |S(v)| from the stored modal
%%% numerators, no solve) is the Ritz-side response next to the Green's-function map,
%%% and the classification diagnostic's background on a 'ritzOnly' result.
%%%
%%% CLASSIFICATION OVERLAY: the 'Show types:' checkboxes apply to every graph
%%% type. On the classification diagnostic they filter its own per-type series
%%% (all on by default); on the Ritz-mode and Green's-function maps they toggle
%%% the same classified points drawn on top as hollow rings
%%% (GuiCommon.styleClassificationOverlay), off by default so those maps open
%%% unchanged. Each of the two contexts remembers its own per-type state across
%%% graph-type switches (diagnosticTypeVisible/overlayTypeVisible).
%%%
%%% SELECTION READOUT: a map click's angle/mode #/velocity/type is kept in
%%% selectionInfo, followed by the classification readout of a classified mode
%%% (method, score in dB, grade, supporting mode) and the shown criterion's value.
%%% Standalone, all of it is shown on the status line; embedded (CombinedViewerGUI,
%%% DoubleSweepViewerGUI) only the classification and criterion lines, since the
%%% modal-shape viewer above already shows the selection (showSelectionInStatus).
%%%
%%% OnModeSelected below is the hand-off point ModalShapeViewerGUI attaches to;
%%% see docs/reference/architecture/gui_viewers.md for the full viewer family.
%%%
%%% SCOPE: a single angular sweep only. The double-sweep dispersion classes
%%% have their own viewers (DispersionViewerGUI, IsofreqAngSweepViewerGUI,
%%% FreqAngSweep2LayeredViewerGUI).
%%%
%%% COLOUR SCALE: the graph types colour their data on fundamentally different
%%% scales (Green / Ritz spectrum map: log10-normalised |G| / |S|; Ritz-mode map: a
%%% criterion value; diagnostic: a categorical type or a score). One shared
%%% cMinField/cMaxField pair stays safe across them by RE-DEFAULTING on every
%%% graph-type change — onGraphTypeChanged reads the freshly-drawn axes' own CLim
%%% right after drawing, before applyAxisRange re-applies it — rather than by
%%% picking one fixed range that would suit only one of them.
%%%
%%% SINGLE-ANGLE DRILL-DOWN (viewSingleAngleBtn): toggles the map content
%%% between the 2D angle x velocity view and an embedded SingleAngleViewerGUI at
%%% obj.currentAngNo, the last-clicked angle or the first if none was clicked.
%%% sweepResult.cs{angNo} is already a fully calcAll()'d CalcSingleHS/
%%% CalcSingleLayered, so drilling down needs no recomputation. The map's
%%% graphTypeDropdown/criterionDropdown stay the one set of controls while drilled
%%% in: the nested viewer is built with hostControls=true (its own two dropdowns
%%% hidden) and each graph type drives the curve of the same position in the
%%% lists (Classification diagnostic -> curve, Green's-function map -> curve,
%%% ...; the two Items lists are kept in the same order for this), through
%%% SingleAngleViewerGUI.setView, so map and curve stay in step both ways; and
%%% the nested viewer's OnModeSelected is forwarded through THIS object's
%%% (dispatchModeSelected), so a caller that wired this object once keeps
%%% working whether the user is looking at the map or at one angle. The
%%% crosshair (lastCrosshairAngle/lastCrosshairVel) is re-drawn on the fresh
%%% axes in both directions of the toggle rather than lost with the old ones.
%%%
%%% EXPORT FIGURE (button beside 'Show types:'): GuiCommon.exportAxes copies the
%%% map as shown into an ordinary figure window, to save in any format. While
%%% drilled in, the nested viewer's own button exports the curve.
%%%
%%% GuliweRR v6.0

    properties
        fig
        mapAxes
        graphTypeDropdown
        criterionLabel
        criterionDropdown
        xMinField, xMaxField, yMinField, yMaxField
        cMinField, cMaxField
        statusLine
        viewSingleAngleBtn
        classBar              % ClassificationBar -- see class docblock
        exportBtn             % 'Export figure' of the map -- see class docblock

        OnModeSelected = []   % function handle @(angle, modeNo, velSnapped, velRaw), for a linked modal-shape viewer
        OnClassificationChanged = []   % function handle @(), fired after a Reclassify
                                       % (refreshClassification chain), for a linked
                                       % modal-shape viewer
        showSelectionInStatus   % true: a map click's selectionInfo is also written to the
                                % status line. Defaults to standalone (no parentContainer)
                                % -- see class docblock's SELECTION READOUT note
    end

    properties (SetAccess = private)
        sweepResult   % CalcAngSweepHS or CalcAngSweepLayered, already calc()'d
        selectionInfo = {}   % {angle, mode #, velocity, type, ...} lines of the last map
                             % click, {} before any -- see SELECTION READOUT note
        isDrilledIn = false       % true while showing an embedded SingleAngleViewerGUI
                                   % instead of the map. Readable externally so a
                                   % composite host (DoubleSweepViewerGUI) can tell
                                   % whether xMinField/yMinField/cMinField are still
                                   % valid graphics handles before touching them: they
                                   % are destroyed by
                                   % delete(obj.contentParent.Children) while drilled
                                   % in, and nothing resets the property references,
                                   % so reading .xMinField.Value then throws "Invalid
                                   % or deleted object"
    end

    properties (Access = private)
        hCrosshair = gobjects(0)
        hMapPlot   = gobjects(0)   % handle(s) the active plot function returned
        currentAngNo = 1          % last-clicked (or default first) angle index --
                                   % drill-down target for viewSingleAngleBtn
        contentParent             % stable uigridlayout slot the map/drill-down
                                   % content is (re)built into -- itself never
                                   % destroyed, see buildMapContent/drillDownToAngle
        lastCrosshairAngle = []   % (angleVal, velVal) last passed to showCrosshairAt,
        lastCrosshairVel   = []   % kept so the marker survives drill-down/return rebuilds
        savedAxisRange = []       % struct(xMin,xMax,yMin,yMax,cMin,cMax) field VALUES
                                   % captured right before drilling in, restored by
                                   % returnToSweepView instead of re-defaulting --
                                   % see onToggleSingleAngleView/returnToSweepView
        typeCheckPanel            % 'Show types:' panel, built once, contents rebuilt
                                   % per draw -- see GuiCommon.rebuildTypeCheckboxes
        typeCheckboxes = gobjects(0)
        diagnosticTypeVisible = struct()   % per-type checkbox states, field = type name --
        overlayTypeVisible    = struct()   % see class docblock's CLASSIFICATION OVERLAY note
    end

    properties (Access = {?AngSweepViewerGUI, ?test_AngSweepViewerGUI})
        singleAngleViewer = []    % embedded SingleAngleViewerGUI while drilled in, else [].
                                   % Friend access to the test class (ModalShapeViewerGUI's
                                   % fullRitzCache precedent), so the refresh chain into the
                                   % drill-down is testable directly
    end

    methods
        function obj = AngSweepViewerGUI(sweepResult, parentContainer, classifyTarget)
        %%% Construct the angular-sweep viewer window
        %%% >> obj = AngSweepViewerGUI(sweepResult, parentContainer, classifyTarget)
        %%%
        %%% INPUT: 1. sweepResult     .. CalcAngSweepHS or CalcAngSweepLayered,
        %%%                              already calc()'d (obj.lastResult off
        %%%                              StructureControlGUI's Run button, typically)
        %%%  (opt) 2. parentContainer .. a uigridlayout region to build into instead
        %%%                              of a fresh own uifigure — used by
        %%%                              CombinedViewerGUI to embed this window's
        %%%                              content in a shared figure. Omitted: this
        %%%                              window owns its own uifigure
        %%%  (opt) 3. classifyTarget  .. the object the classification bar's Reclassify
        %%%                              classifies (default: sweepResult) --
        %%%                              DoubleSweepViewerGUI passes the whole double
        %%%                              sweep, of which sweepResult is one wavelength
        %%%                              (see ClassificationBar)
            validateattributes(sweepResult, {'CalcAngSweepHS', 'CalcAngSweepLayered'}, ...
                {'scalar'}, 'AngSweepViewerGUI', 'sweepResult');
            %%% .cs is populated by the CONSTRUCTOR itself (uncalculated CalcSingleHS/
            %%% Layered objects, one per angle) -- not a valid "has calc() run" signal.
            %%% significantModesMap stays {} until assembleMaps() runs at the end of
            %%% calc().
            if isempty(sweepResult.significantModesMap)
                error('AngSweepViewerGUI:constructor:notCalculated', ...
                    'sweepResult has not been calc()''d yet -- nothing to view.');
            end
            obj.sweepResult = sweepResult;
            if nargin < 2; parentContainer = []; end
            if nargin < 3; classifyTarget = []; end
            obj.showSelectionInStatus = isempty(parentContainer);
            obj.buildFigure(parentContainer, classifyTarget);
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
        function buildFigure(obj, parentContainer, classifyTarget)
            [obj.fig, mainLayout] = GuiCommon.hostLayout(parentContainer, ...
                'GuliweRR -- Angular Sweep Viewer', [100 100 1000 730], [4 1]);
            mainLayout.RowHeight = {40, 30, '1x', 26};

            %%% top bar -- graph-type + criterion selectors + single-angle
            %%% drill-down toggle. Stays put across drill-down/back (unlike
            %%% obj.contentParent's children) -- viewSingleAngleBtn is always
            %%% reachable regardless of which content is currently showing.
            %%% columns: graph-type label/dropdown | criterion label/dropdown |
            %%% spacer | drill-down toggle, kept right-aligned
            topBar = uigridlayout(mainLayout, [1 6]);
            topBar.Layout.Row = 1;
            topBar.ColumnWidth = {80, 190, 80, 160, '1x', 170};
            topBar.Padding = [0 0 0 0];

            lbl = uilabel(topBar, 'Text', 'Graph type:');
            lbl.Layout.Column = 1;
            obj.graphTypeDropdown = uidropdown(topBar, ...
                'Items', {'Classification diagnostic', 'Green''s-function map', 'Ritz-mode map', ...
                          'Ritz spectrum map'}, ...
                'ValueChangedFcn', @(~,~) obj.onGraphTypeChanged());
            obj.graphTypeDropdown.Layout.Column = 2;

            obj.criterionLabel = uilabel(topBar, 'Text', 'Criterion:', 'Visible', 'off');
            obj.criterionLabel.Layout.Column = 3;
            obj.criterionDropdown = GuiCommon.ritzCriterionDropdown(topBar, ...
                isa(obj.sweepResult, 'CalcAngSweepLayered'), @(~,~) obj.onGraphTypeChanged());
            obj.criterionDropdown.Layout.Column = 4;

            obj.viewSingleAngleBtn = uibutton(topBar, 'Text', 'View single angle ->', ...
                'ButtonPushedFcn', @(~,~) obj.onToggleSingleAngleView());
            obj.viewSingleAngleBtn.Layout.Column = 6;

            %%% classification bar -- see class docblock
            obj.classBar = ClassificationBar(mainLayout, obj.sweepResult, 2, classifyTarget);
            obj.classBar.setReclassifyHandler(@() obj.onReclassified());

            %%% stable content slot -- buildMapContent()/drillDownToAngle() swap
            %%% what's INSIDE it (never the slot itself), so this uigridlayout
            %%% is created exactly once regardless of how many times the user
            %%% toggles the drill-down button
            obj.contentParent = uigridlayout(mainLayout, [1 1]);
            obj.contentParent.Layout.Row = 3;
            obj.contentParent.Padding = [0 0 0 0];

            obj.buildMapContent();

            obj.statusLine = GuiCommon.statusLabel(mainLayout, 4);
        end

        function buildMapContent(obj)
        %%% (Re)build the normal map + right-hand side panel (axis-range
        %%% controls, type checkboxes) inside obj.contentParent -- called once at
        %%% construction (from buildFigure) and again by returnToSweepView()
        %%% after a drill-down/back cycle has destroyed the previous content.
        %%% >> buildMapContent(obj)
            [obj.mapAxes, sideGl] = GuiCommon.plotWithSidePanel(obj.contentParent, ...
                {18, 26, 18, 26, 18, 26, 10, 24, '1x'}, ...
                @(src, event) obj.onMapClicked(src, event));
            obj.hCrosshair = gobjects(0);   % previous handle, if any, died with the old axes

            onRangeEdit = @(~,~) obj.applyAxisRange();

            [obj.xMinField, obj.xMaxField] = GuiCommon.rangeFieldPair(sideGl, ...
                'Angle [deg] min/max', num2str(min(obj.sweepResult.angles)), ...
                num2str(max(obj.sweepResult.angles)), onRangeEdit);

            %%% default to the Green velocity axis' own range where available
            %%% (calcMode='full'), up to its end (1.05 x the largest longitudinal
            %%% limiting velocity + 0.1 km/s for a half-space) but starting at 0.9 x
            %%% the slowest Ritz mode of any angle: the axis itself starts at 0.3 km/s
            %%% (SawSpectrum.VEL_MIN), far below any surface wave. Without it
            %%% ('ritzOnly'), the Ritz modes' own range over all angles, as
            %%% SingleAngleViewerGUI uses for one angle.
            if ~isempty(obj.sweepResult.greenVelAxis)
                ritzVelMin  = min(cellfun(@(c) min(c.ritzSolver.vels), obj.sweepResult.cs));
                yMinDefault = num2str(max(min(obj.sweepResult.greenVelAxis), 0.9 * ritzVelMin));
                yMaxDefault = num2str(max(obj.sweepResult.greenVelAxis));
            else
                compVels = @(c) c.ritzSolver.vels(1:c.ritzSolver.nCompModes);
                yMinDefault = num2str(min(cellfun(@(c) min(compVels(c)), obj.sweepResult.cs)));
                yMaxDefault = num2str(max(cellfun(@(c) max(compVels(c)), obj.sweepResult.cs)));
            end
            [obj.yMinField, obj.yMaxField] = GuiCommon.rangeFieldPair(sideGl, ...
                'Velocity [km/s] min/max', yMinDefault, yMaxDefault, onRangeEdit);

            %%% Colour min/max -- start empty; onGraphTypeChanged (called right
            %%% after this method returns) fills them from the first draw's own
            %%% CLim, the same re-default-per-type mechanic as every later
            %%% graph-type switch (see class docblock).
            [obj.cMinField, obj.cMaxField] = GuiCommon.rangeFieldPair(sideGl, ...
                'Colour min/max', '', '', onRangeEdit);

            %%% 'Show types:' header with the Export figure button, then the
            %%% checkbox panel -- built here empty, contents filled in by
            %%% GuiCommon.rebuildTypeCheckboxes on every onGraphTypeChanged draw
            %%% (see class docblock's CLASSIFICATION OVERLAY note). Takes all the
            %%% remaining side-panel height.
            [~, obj.exportBtn] = GuiCommon.typesHeaderRow(sideGl, 8, @() obj.onExportFigure());
            obj.typeCheckPanel = uipanel(sideGl, 'BorderType', 'none');
            obj.typeCheckPanel.Layout.Row = 9;
        end
    end

    %%% ===================================================================
    %%% Callbacks -- public (not Access=private), StructureControlGUI's
    %%% convention, so tests can call them directly without going through
    %%% actual UI events.
    %%% ===================================================================
    methods
        function onGraphTypeChanged(obj)
        %%% Redraw the map for whichever graph type is now selected -- or, while
        %%% drilled in, the single-angle curve of the same graph type (see class
        %%% docblock)
        %%% >> onGraphTypeChanged(obj)
            isRitz = strcmp(obj.graphTypeDropdown.Value, 'Ritz-mode map');
            obj.criterionLabel.Visible = isRitz;
            obj.criterionDropdown.Visible = isRitz;
            if obj.isDrilledIn
                obj.syncDrilledInView();
                return
            end

            GuiCommon.clearAxes(obj.mapAxes);
            obj.mapAxes.Color = [1 1 1];   % undo plotRitzModeMap's black background, if it was active
            %%% cla() clears CHILDREN only, not axis-limit modes -- CLimMode
            %%% stays 'manual' from a previous plotGreenMap draw (the only
            %%% one of the three that calls clim() itself) otherwise, so a
            %%% later switch to Ritz-mode map/classification diagnostic
            %%% (neither sets CLim itself) would silently inherit a stale
            %%% range instead of auto-scaling to its own new data, breaking
            %%% the colour-range fields' re-default per type.
            obj.mapAxes.CLimMode = 'auto';
            obj.hCrosshair = gobjects(0);

            try
                hBg = gobjects(0);   % background object(s) drawn UNDER hOut, if any -- see below
                switch obj.graphTypeDropdown.Value
                    case 'Ritz-mode map'
                        [hOut, ~] = Vis.plotRitzModeMap(obj.sweepResult, obj.criterionDropdown.Value, [], obj.mapAxes);
                    case 'Green''s-function map'
                        %%% 'ritzOnly' (no Green data): plotGreenMap errors into
                        %%% the status line
                        [hOut, ~] = Vis.plotGreenMap(obj.sweepResult, [], [], [], obj.mapAxes);
                    case 'Ritz spectrum map'
                        [hOut, ~] = Vis.plotRitzSpectrumMap(obj.sweepResult, [], [], [], obj.mapAxes);
                    otherwise   % 'Classification diagnostic'
                        %%% hBg is the pcolor background (Green, or the Ritz spectrum
                        %%% on 'ritzOnly'), BEHIND hOut's scatter and covering the
                        %%% WHOLE axes. It gets the same ButtonDownFcn below, or this
                        %%% opaque surface would swallow every click that misses a
                        %%% scatter marker (only significant points would then reach
                        %%% onMapClicked).
                        [hOut, ~, hBg] = Vis.plotClassificationDiagnostic(obj.sweepResult, [], [], [], obj.mapAxes);
                end
                obj.hMapPlot = hOut;
                set([obj.hMapPlot(:); hBg(:)], 'ButtonDownFcn', @(src, event) obj.onMapClicked(src, event));
                %%% re-default the colour-range fields to THIS graph type's
                %%% own just-drawn CLim -- see class docblock: this is what
                %%% makes one shared C min/max control safe across three
                %%% different-scale graph types. Read before any overlay is
                %%% drawn (its RGB-coloured rings leave CLim alone anyway).
                obj.cMinField.Value = num2str(obj.mapAxes.CLim(1));
                obj.cMaxField.Value = num2str(obj.mapAxes.CLim(2));
                %%% type checkboxes -- see class docblock's CLASSIFICATION OVERLAY note
                obj.typeCheckboxes = GuiCommon.classificationTypeCheckboxes(obj.typeCheckPanel, ...
                    strcmp(obj.graphTypeDropdown.Value, 'Classification diagnostic'), hOut, ...
                    @() obj.drawClassificationOverlay(), obj.diagnosticTypeVisible, ...
                    obj.overlayTypeVisible, @(stateProp, typeName, isOn) obj.rememberTypeVisible(stateProp, typeName, isOn));
                GuiCommon.setStatusOk(obj.statusLine, obj.idleStatusText());
            catch err
                obj.hMapPlot = gobjects(0);
                obj.typeCheckboxes = GuiCommon.rebuildTypeCheckboxes(obj.typeCheckPanel, gobjects(0));
                GuiCommon.setStatusError(obj.statusLine, err.message);
            end

            obj.applyAxisRange();
        end

        function fig = onExportFigure(obj)
        %%% Export figure button: the map as shown, in an ordinary figure window
        %%% >> fig = onExportFigure(obj)
            fig = GuiCommon.exportAxes(obj.mapAxes);
        end

        function applyAxisRange(obj)
        %%% Apply the angle/velocity/colour range fields to the map axes
        %%% >> applyAxisRange(obj)
            GuiCommon.applyFieldLimit(obj.mapAxes, 'x', obj.xMinField, obj.xMaxField);
            GuiCommon.applyFieldLimit(obj.mapAxes, 'y', obj.yMinField, obj.yMaxField);
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
        %%% Snap a map click to the nearest computed Ritz mode -- nearest angle,
        %%% then nearest mode velocity at that angle -- and fire OnModeSelected.
        %%% >> onMapClicked(obj, src, event)
            pt = GuiCommon.clickPoint(obj.mapAxes, event);
            angPos = pt(1); velPos = pt(2);

            [~, angNo] = min(abs(obj.sweepResult.angles - angPos));
            obj.currentAngNo = angNo;   % drillDownToAngle's target if the user hits viewSingleAngleBtn next
            ritz  = obj.sweepResult.cs{angNo}.ritzSolver;
            nM    = ritz.nCompModes;
            vels  = ritz.vels(1:nM);
            [~, modeNo] = min(abs(vels - velPos));
            velSel = vels(modeNo);
            angleVal = obj.sweepResult.angles(angNo);

            sm = obj.sweepResult.significantModesMap{angNo};
            typeStr = GuiCommon.significantModeType(sm, modeNo);

            obj.selectionInfo = [GuiCommon.selectionLines(angleVal, modeNo, velSel, typeStr), ...
                GuiCommon.selectionExtras(obj.sweepResult.cs{angNo}, sm, modeNo, obj.activeCriterion())];
            GuiCommon.setStatusOk(obj.statusLine, obj.idleStatusText());

            obj.showCrosshairAt(angleVal, velSel);

            if ~isempty(obj.OnModeSelected)
                %%% velPos, the RAW clicked velocity before snapping, is passed
                %%% alongside the snapped (modeNo, velSel) pair. The listening
                %%% modal-shape viewer picks according to its own source: modeNo for
                %%% a Ritz modal shape, velPos directly for a Green's response,
                %%% which has no discrete mode to snap to. This class stays unaware
                %%% of that state.
                obj.OnModeSelected(angleVal, modeNo, velSel, velPos);
            end
        end

        function showCrosshairAt(obj, angleVal, velVal)
        %%% Move (or create) the crosshair marker at (angleVal, velVal) WITHOUT
        %%% triggering any click/mode-lookup logic -- the entry point a linked
        %%% modal-shape viewer uses to mark its own selection on this map, and
        %%% public for that reason.
        %%%
        %%% While drilled into a single angle (isDrilledIn) this forwards to the
        %%% nested SingleAngleViewerGUI instead of drawing on the map, which does
        %%% not exist in the content slot then. If angleVal is a DIFFERENT angle
        %%% from the one drilled into -- a linked viewer's angle stepper having
        %%% moved on -- it re-drills first, rather than leaving the single-angle
        %%% curve showing a stale angle with no indication. Nearest-angle match,
        %%% as onMapClicked uses, not an exact match against sweepResult.angles.
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
        %%% SingleAngleViewerGUI at obj.currentAngNo (the last-clicked angle,
        %%% or the first angle if nothing has been clicked yet).
        %%% >> onToggleSingleAngleView(obj)
            if obj.isDrilledIn
                obj.returnToSweepView();
            else
                %%% capture the CURRENT axis-range field values before
                %%% drillDownToAngle destroys them -- returnToSweepView restores
                %%% them rather than re-defaulting, so a custom colour/angle/
                %%% velocity range survives the round trip
                obj.savedAxisRange = GuiCommon.captureAxisRange(obj);
                obj.drillDownToAngle(obj.currentAngNo);
            end
        end

        function drillDownToAngle(obj, angNo)
        %%% Replace the map content with an embedded SingleAngleViewerGUI for
        %%% sweepResult.cs{angNo} -- already calcAll()'d as part of the sweep
        %%% itself (CalcAngSweepHS/Layered's own per-angle building block), so
        %%% no recomputation is needed. The nested viewer's OnModeSelected is
        %%% forwarded through dispatchModeSelected so CombinedViewerGUI's
        %%% existing wiring (set on THIS object, not the nested one) keeps
        %%% driving the modal-shape viewer above while drilled in.
        %%% >> drillDownToAngle(obj, angNo)
            delete(obj.contentParent.Children);
            obj.hCrosshair = gobjects(0);

            singleResult = obj.sweepResult.cs{angNo};
            %%% no bar of its own: the sweep's bar classifies every angle; no graph-type
            %%% or criterion dropdown of its own: this viewer's drive it
            obj.singleAngleViewer = SingleAngleViewerGUI(singleResult, obj.contentParent, false, true);
            %%% embedded, so off by the nested viewer's own default -- but a
            %%% STANDALONE map window has no shape viewer to show the selection
            obj.singleAngleViewer.showSelectionInStatus = obj.showSelectionInStatus;
            obj.singleAngleViewer.OnModeSelected = ...
                @(angle, modeNo, velSnap, velRaw) obj.dispatchModeSelected(angle, modeNo, velSnap, velRaw);
            obj.singleAngleViewer.show();
            obj.currentAngNo = angNo;
            obj.isDrilledIn = true;
            obj.syncDrilledInView();

            %%% restore the crosshair if the last one shown belongs to THIS
            %%% angle (see class docblock)
            if ~isempty(obj.lastCrosshairAngle)
                [~, crosshairAngNo] = min(abs(obj.sweepResult.angles - obj.lastCrosshairAngle));
                if crosshairAngNo == angNo
                    obj.singleAngleViewer.showCrosshairAt(obj.lastCrosshairAngle, obj.lastCrosshairVel);
                end
            end

            obj.viewSingleAngleBtn.Text = '<< Back to angular sweep';
        end

        function returnToSweepView(obj)
        %%% Tear down the embedded SingleAngleViewerGUI and rebuild the normal
        %%% map content in its place, restoring the axis-range fields
        %%% onToggleSingleAngleView saved right before drilling in rather than
        %%% keeping the rebuild's fresh defaults.
        %%% >> returnToSweepView(obj)
            obj.singleAngleViewer = [];
            obj.isDrilledIn = false;   % before the redraw, which forwards while drilled in
            delete(obj.contentParent.Children);
            obj.buildMapContent();
            obj.onGraphTypeChanged();
            GuiCommon.restoreAxisRange(obj, obj.savedAxisRange);

            %%% re-draw the crosshair: buildMapContent/onGraphTypeChanged both
            %%% reset hCrosshair to gobjects(0) on the fresh axes
            if ~isempty(obj.lastCrosshairAngle)
                obj.hCrosshair = GuiCommon.markCrosshair(obj.mapAxes, obj.hCrosshair, ...
                    obj.lastCrosshairAngle, obj.lastCrosshairVel);
            end

            obj.viewSingleAngleBtn.Text = 'View single angle ->';
        end

        function rememberTypeVisible(obj, stateProp, typeName, isOn)
        %%% A type checkbox's onToggle -- store its new state so the next draw
        %%% of the same context (diagnostic or overlay) restores it
        %%% >> rememberTypeVisible(obj, stateProp, typeName, isOn)
            obj.(stateProp).(typeName) = isOn;
        end

        function refreshClassification(obj)
        %%% Redraw after the sweep's classification changed (a Reclassify here or in
        %%% a host): the bar, then the map with its axis ranges and crosshair kept,
        %%% or the drilled-in single-angle curve
        %%% >> refreshClassification(obj)
            obj.classBar.refresh();
            obj.redrawClassified();
        end

        function onReclassified(obj)
        %%% The classification bar's reclassify handler: redraw, then notify a host
        %%% (CombinedViewerGUI: the linked modal-shape viewer)
        %%% >> onReclassified(obj)
            obj.redrawClassified();
            if ~isempty(obj.OnClassificationChanged)
                obj.OnClassificationChanged();
            end
        end

        function dispatchModeSelected(obj, angle, modeNo, velSnap, velRaw)
        %%% Forward a click inside the drilled-down SingleAngleViewerGUI out
        %%% through THIS object's own OnModeSelected -- same guard pattern
        %%% onMapClicked itself already uses.
        %%% >> dispatchModeSelected(obj, angle, modeNo, velSnap, velRaw)
            if ~isempty(obj.OnModeSelected)
                obj.OnModeSelected(angle, modeNo, velSnap, velRaw);
            end
        end
    end

    %%% ===================================================================
    %%% Run-time helpers
    %%% ===================================================================
    methods (Access = private)
        function syncDrilledInView(obj)
        %%% Show on the drilled-in curve the graph type at the same list position
        %%% as the map's, with the map's criterion (see class docblock)
        %%% >> syncDrilledInView(obj)
            idx = find(strcmp(obj.graphTypeDropdown.Items, obj.graphTypeDropdown.Value), 1);
            curveTypes = obj.singleAngleViewer.graphTypeDropdown.Items;
            obj.singleAngleViewer.setView(curveTypes{idx}, obj.criterionDropdown.Value);
        end

        function hOverlay = drawClassificationOverlay(obj)
        %%% Draw the classified points on top of the Ritz-mode/Green's-function
        %%% map already on mapAxes -- plotClassificationDiagnostic itself, without
        %%% its Green background, restyled as an overlay (see class docblock)
        %%% >> hOverlay = drawClassificationOverlay(obj)
            savedLabels = GuiCommon.captureAxesLabels(obj.mapAxes);
            [hOverlay, ~, ~, hLim] = Vis.plotClassificationDiagnostic( ...
                obj.sweepResult, 'type', false, [], obj.mapAxes);
            GuiCommon.styleClassificationOverlay(obj.mapAxes, hOverlay, hLim, savedLabels);
        end

        function redrawClassified(obj)
        %%% Redraw the classified content: the drilled-in curve, or the map with the
        %%% axis-range fields and the crosshair kept (onGraphTypeChanged alone would
        %%% re-default the colour range and drop the crosshair)
        %%% >> redrawClassified(obj)
            obj.selectionInfo = {};
            if obj.isDrilledIn
                obj.singleAngleViewer.refreshClassification();
                return
            end
            saved = GuiCommon.captureAxisRange(obj);
            obj.onGraphTypeChanged();
            GuiCommon.restoreAxisRange(obj, saved);
            if ~isempty(obj.lastCrosshairAngle)
                obj.hCrosshair = GuiCommon.markCrosshair(obj.mapAxes, obj.hCrosshair, ...
                    obj.lastCrosshairAngle, obj.lastCrosshairVel);
            end
        end

        function msg = idleStatusText(obj)
        %%% Status-line text after a successful draw or click (GuiCommon.statusText)
        %%% >> msg = idleStatusText(obj)
            msg = GuiCommon.statusText(obj.selectionInfo, obj.showSelectionInStatus, ...
                obj.activeCriterion(), isa(obj.sweepResult, 'CalcAngSweepLayered'));
        end

        function crit = activeCriterion(obj)
        %%% The Ritz-mode criterion while the Ritz-mode map is shown, else ''
        %%% >> crit = activeCriterion(obj)
            crit = '';
            if strcmp(obj.graphTypeDropdown.Value, 'Ritz-mode map')
                crit = obj.criterionDropdown.Value;
            end
        end
    end
end
