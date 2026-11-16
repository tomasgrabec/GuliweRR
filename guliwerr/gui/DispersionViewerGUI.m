classdef DispersionViewerGUI < handle
%%% DispersionViewerGUI -- fixed-angle wavelength/frequency dispersion
%%% viewer. See docs/reference/architecture/gui_viewers.md.
%%%
%%% Modal shapes plus controls on one side, a dispersion map on the other. ONE
%%% class rather than a wired-together pair, because switching wavelength step
%%% here switches WHICH CalcAngSweepLayered/ritzSolver the modal shape comes
%%% from (sweepResult.sweeps{iw}), not just which mode within a fixed one.
%%% ModalShapeViewerGUI's angle-oriented API does not fit that, so this class
%%% renders modal shapes directly via Vis.plotModalShape, driven by its own
%%% (wavelength step, branch) state.
%%%
%%% windowMode: 'combined' (default, ONE window, self-building the same
%%% shapes-left/controls-right-over-map layout CombinedViewerGUI uses) |
%%% 'split' (TWO windows, a modal-shape figure and a map figure).
%%%
%%% sweepResult MAY carry more than one angle: angleIdx (3rd constructor
%%% argument, default 1) picks which one THIS view is fixed at. It adds no
%%% angle stepper — this is always a SINGLE fixed-angle view — it only lets a
%%% genuine multi-angle CalcWavevecAngSweep2Layered be sliced at a chosen
%%% angle. DoubleSweepViewerGUI is what offers angle switching, by re-embedding
%%% a fresh DispersionViewerGUI at a new angleIdx each time.
%%%
%%% EMBEDDING: parentContainer (4th constructor argument) hosts this class's
%%% content inside a caller-owned uigridlayout region instead of a fresh
%%% uifigure — the mechanism DoubleSweepViewerGUI uses. One region is enough,
%%% unlike ModalShapeViewerGUI's two-argument controlsParent/shapeParent,
%%% because this class is already a self-contained composite and builds its own
%%% internal shapes/controls/map split inside whatever it is given. windowMode
%%% must be 'combined' when parentContainer is given: 'split' needs two
%%% independent windows, which is meaningless once this class owns none.
%%%
%%% CLASSIFICATION BAR (above the dispersion map, a ClassificationBar): method
%%% (Green | Ritz) and standard for the WHOLE wavelength sweep -- the bar reads its
%%% state off sweeps{1} and calls sweepResult.classify, which reclassifies every
%%% wavelength and reassembles the maps. refreshClassification then redraws the
%%% branches (Vis.plotDispersionCurves highlights the classified ones, read live
%%% from significantModesMap) and the info panel's type line. The tracked branches
%%% themselves do not depend on the classification.
%%%
%%% No Green's-response toggle (CalcWavevecAngSweep2Layered is layered-only,
%%% with no half-space equivalent) and no SAW/limiting override panel. Branch
%%% "prev/next" stepping needs no "follow physical wave" logic the way
%%% ModalShapeViewerGUI's angle-stepping does: continuity across wavelength
%%% steps is already baked into branchModeNoMap
%%% (DispersionBranches.track), so stepping a branch is just
%%% incrementing a column index.
%%%
%%% GuliweRR v6.0

    properties
        fig       % hosts the modal-shape half (shapes uiaxes grid + controls)
        mapFig    % hosts the dispersion-curve map + its controls; == fig in
                  % 'combined' windowMode, a separate uifigure in 'split'

        wavDropdown, prevWavBtn, nextWavBtn
        branchField, prevBranchBtn, nextBranchBtn
        velLabel
        displayDropdown     % 'components' | 'rgb' -- GuiCommon.shapeDisplayDropdown
        shapeStatusLine

        mapAxes
        curveTypeDropdown
        showRitzBranchesCheckbox
        showGreenBgCheckbox
        xRangeLabel, yRangeLabel
        xMinField, xMaxField, yMinField, yMaxField
        infoPanel
        mapStatusLine
        classBar            % ClassificationBar -- see class docblock
    end

    properties (SetAccess = private)
        sweepResult       % CalcWavevecAngSweep2Layered, already calc()'d
        windowMode        % 'combined' (default) | 'split' -- see class docblock
        angleIdx          % index into sweepResult.angles this view is fixed at (default 1) -- see class docblock
        branchModeNoMap   % [nWav x nCompModes], DispersionBranches.track's own output
        currentWavNo      % index into sweepResult.wavelengths
        currentBranchIdx  % column index into branchModeNoMap (a tracked physical branch)
    end

    properties (Access = private)
        shapeAxesCell = {}   % {3 x nLay} pre-built uiaxes, same pattern as ModalShapeViewerGUI's own
        rgbAxesCell   = {}   % {1 x nLay} uiaxes of the polarisation-RGB display
        keyAx         = []   % its colour-key uiaxes
        hLinePlots = gobjects(0)   % current map's line handles, one per branch
        hCrosshair = gobjects(0)
        hasGreenData = false   % sweeps{1}.greenVelAxis non-empty, i.e. sweepResult was calc()'d with
                                 % calcMode='full' -- computed once at construction, gates
                                 % showGreenBgCheckbox.Enable alongside curveType (see onCurveTypeChanged)
        drawnCurveType = ''    % curveType of the last redrawMap -- range fields reset only when it changes
    end

    properties (Constant, Access = private)
        VPHASE_CURVE_TYPES = {'wavelength-vphase', 'k-vphase'}   % the only two curveTypes
                                                                     % Vis.plotDispersionCurves'
                                                                     % Green background supports
                                                                     % (velocity y-axis) -- mirrored
                                                                     % here so the checkbox can be
                                                                     % greyed out instead of ever
                                                                     % hitting that function's error
    end

    methods
        function obj = DispersionViewerGUI(sweepResult, windowMode, angleIdx, parentContainer)
        %%% Construct the dispersion viewer
        %%% >> obj = DispersionViewerGUI(sweepResult, windowMode, angleIdx, parentContainer)
        %%%
        %%% INPUT: 1. sweepResult .. CalcWavevecAngSweep2Layered, already
        %%%                          calc()'d, with retentionLevel='full'
        %%%                          (needed for
        %%%                          DispersionBranches.track
        %%%                          -- its own missingEigVecs error surfaces
        %%%                          directly if this isn't satisfied)
        %%%  (opt) 2. windowMode  .. 'combined' (default) | 'split' -- see
        %%%                          class docblock. Must be 'combined' (or
        %%%                          omitted) when parentContainer is given.
        %%%  (opt) 3. angleIdx    .. which angle (default 1) -- see class
        %%%                          docblock's "sweepResult MAY carry more
        %%%                          than one angle" paragraph
        %%%  (opt) 4. parentContainer .. a uigridlayout region to build into
        %%%                          instead of a fresh own uifigure -- see
        %%%                          class docblock's "Embedding" paragraph
            validateattributes(sweepResult, {'CalcWavevecAngSweep2Layered'}, ...
                {'scalar'}, 'DispersionViewerGUI', 'sweepResult');
            if nargin < 2 || isempty(windowMode); windowMode = 'combined'; end
            GuiCommon.validateWindowMode('DispersionViewerGUI', windowMode, {'combined', 'split'});
            %%% angleIdx bounds check -- doesn't need calc() to have run
            %%% (sweepResult.angles is set at construction), so an
            %%% out-of-range angleIdx gets its own specific error even
            %%% before it's calc()'d, same "cheap checks first" ordering the
            %%% old single-angle check used
            if nargin < 3 || isempty(angleIdx); angleIdx = 1; end
            validateattributes(angleIdx, {'double'}, ...
                {'scalar', 'integer', 'positive', '<=', numel(sweepResult.angles)}, ...
                'DispersionViewerGUI', 'angleIdx');
            if nargin < 4; parentContainer = []; end
            GuiCommon.requireStandaloneForSplit('DispersionViewerGUI', windowMode, ...
                parentContainer, 'parentContainer');
            if isempty(sweepResult.sweeps)
                error('DispersionViewerGUI:constructor:notCalculated', ...
                    'sweepResult has not been calc()''d yet -- nothing to view.');
            end

            obj.sweepResult     = sweepResult;
            obj.windowMode      = windowMode;
            obj.angleIdx        = angleIdx;
            obj.branchModeNoMap = DispersionBranches.track(sweepResult, angleIdx);
            obj.hasGreenData    = ~isempty(sweepResult.sweeps{1}.greenVelAxis);

            obj.buildFigures(parentContainer);
            obj.goTo(1, 1);
        end

        function show(obj)
        %%% Make every window visible (constructed invisible so headless/
        %%% test use is safe) -- a no-op on the mapFig side beyond the first
        %%% call in 'combined' windowMode, since .fig and .mapFig are then
        %%% the same figure.
        %%% >> show(obj)
            obj.fig.Visible = 'on';
            obj.mapFig.Visible = 'on';
        end
    end

    %%% ===================================================================
    %%% Figure construction
    %%% ===================================================================
    methods (Access = private)
        function buildFigures(obj, parentContainer)
            if isempty(parentContainer) && strcmp(obj.windowMode, 'split')
                %%% 'split' -- two independent windows, so the map cannot
                %%% share a 2-row layout with the shapes
                obj.fig = uifigure('Name', 'GuliweRR -- Dispersion Viewer: Modal Shapes', ...
                    'Position', [80 460 1400 440], 'Visible', 'off');
                [shapeParent, controlsParent] = GuiCommon.shapeControlsSplit(obj.fig, [4 4 4 4]);

                obj.mapFig = uifigure('Name', 'GuliweRR -- Dispersion Viewer: Map', ...
                    'Position', [80 60 1000 350], 'Visible', 'off');
                mapParent = uigridlayout(obj.mapFig, [1 1]);
                mapParent.Padding = [4 4 4 4];
            else
                %%% 'combined' standalone AND embedded (windowMode is
                %%% guaranteed 'combined' when parentContainer is given -- the
                %%% constructor rejects 'split' with it) build the SAME
                %%% shapes+controls-over-map 2-row layout; hostLayout is the
                %%% only difference between owning a uifigure and living in
                %%% the caller's own region (see class docblock, "Embedding").
                [obj.fig, mainLayout] = GuiCommon.hostLayout(parentContainer, ...
                    'GuliweRR -- Dispersion Viewer', [80 60 1200 900], [2 1]);
                obj.mapFig = obj.fig;
                mainLayout.RowHeight = {'44x', '56x'};
                if ~isempty(parentContainer)
                    mainLayout.Padding = [0 0 0 0];   % the caller already padded its own region
                end

                [shapeParent, controlsParent, topLayout] = ...
                    GuiCommon.shapeControlsSplit(mainLayout, [0 0 0 0]);
                topLayout.Layout.Row = 1;

                mapParent = uigridlayout(mainLayout, [1 1]);
                mapParent.Layout.Row = 2; mapParent.Padding = [0 0 0 0];
            end

            obj.buildShapeAxesGrid(shapeParent);
            obj.buildShapeControls(controlsParent);
            obj.buildMapPanel(mapParent);
        end

        function buildShapeAxesGrid(obj, shapeParent)
        %%% Pre-built {3 x nLay} uiaxes grid -- geometry owned by
        %%% GuiCommon.shapeAxesGrid, shared with the ModalShape viewer family
        %%% (see that helper's docblock for why the axes are manually
        %%% positioned inside a plain uipanel rather than a uigridlayout).
            relDepths = [obj.sweepResult.structure.calcParams.relDepth];
            [obj.shapeAxesCell, obj.rgbAxesCell, obj.keyAx] = GuiCommon.shapeAxesGrid(shapeParent, relDepths);
        end

        function buildShapeControls(obj, controlsParent)
            mainLayout = uigridlayout(controlsParent, [5 1]);
            mainLayout.RowHeight = {36, 36, 36, '1x', 26};
            mainLayout.Padding = [4 4 4 4];

            wavStrs = arrayfun(@(w) sprintf('%.4g', w), obj.sweepResult.wavelengths, 'UniformOutput', false);
            [obj.wavDropdown, obj.prevWavBtn, obj.nextWavBtn] = ...
                GuiCommon.dropdownStepperRow(mainLayout, 1, 'Wavelength:', wavStrs, {70, 120, 36, 36}, ...
                    @(~,~) obj.goTo(obj.wavDropdown.Value, obj.currentBranchIdx), ...
                    @(~,~) obj.onPrevWav(), @(~,~) obj.onNextWav());

            [obj.branchField, obj.prevBranchBtn, obj.nextBranchBtn, obj.velLabel] = ...
                GuiCommon.indexStepperRow(mainLayout, 2, 'Branch:', {60, 60, 36, 36, '1x'}, ...
                    @(~,~) obj.goTo(obj.currentWavNo, obj.branchField.Value), ...
                    @(~,~) obj.onPrevBranch(), @(~,~) obj.onNextBranch());

            %%% row 3 ('1x') is the shape display itself in the caller's own
            %%% layout -- nothing to build here, matches ModalShapeViewerGUI's
            %%% own row-allocation style (an empty/flexible spacer row)
            %%% row 3 -- display: signed components / polarisation RGB
            displayRow = uigridlayout(mainLayout, [1 2]);
            displayRow.Layout.Row = 3;
            displayRow.ColumnWidth = {60, 200};
            displayRow.Padding = [0 0 0 0];
            uilabel(displayRow, 'Text', 'Display:');
            obj.displayDropdown = GuiCommon.shapeDisplayDropdown(displayRow, @(~,~) obj.onDisplayChanged());

            spacer = uipanel(mainLayout, 'BorderType', 'none');
            spacer.Layout.Row = 4;

            obj.shapeStatusLine = GuiCommon.statusLabel(mainLayout, 5);
        end

        function buildMapPanel(obj, mapParent)
            mainLayout = uigridlayout(mapParent, [4 1]);
            mainLayout.RowHeight = {36, 30, '1x', 26};
            mainLayout.Padding = [4 4 4 4];

            topBar = uigridlayout(mainLayout, [1 5]);
            topBar.Layout.Row = 1;
            topBar.ColumnWidth = {90, 220, 170, 190, 140};
            topBar.Padding = [0 0 0 0];
            uilabel(topBar, 'Text', 'Curve type:');
            obj.curveTypeDropdown = uidropdown(topBar, ...
                'Items', {'wavelength-freq', 'wavelength-vphase', 'k-freq', 'k-vphase', 'freq-vphase', 'freq-vgroup'}, ...
                'ValueChangedFcn', @(~,~) obj.onCurveTypeChanged());

            %%% Show Ritz branches -- unlike the Green checkbox below, available for
            %%% every curveType and regardless of calcMode (the branch lines exist
            %%% either way), so it needs no Enable-gating.
            obj.showRitzBranchesCheckbox = uicheckbox(topBar, 'Text', 'Show Ritz branches', ...
                'Value', true, 'ValueChangedFcn', @(~,~) obj.onCurveTypeChanged());

            %%% Show Green background -- meaningful only for the 2 v_phase curveTypes
            %%% and only if sweepResult was calc()'d with calcMode='full'
            %%% (obj.hasGreenData). Greyed out otherwise rather than left enabled to
            %%% hit Vis.plotDispersionCurves' hard error; Enable is (re)set by
            %%% onCurveTypeChanged alongside the redraw.
            obj.showGreenBgCheckbox = uicheckbox(topBar, 'Text', 'Show Green background', ...
                'Value', false, ...
                'Enable', obj.hasGreenData && ismember(obj.curveTypeDropdown.Value, obj.VPHASE_CURVE_TYPES), ...
                'ValueChangedFcn', @(~,~) obj.onCurveTypeChanged());

            %%% Read-only: angleIdx is fixed for this whole view, not a stepper.
            %%% Informative only once sweepResult carries more than one angle, but
            %%% shown unconditionally.
            uilabel(topBar, 'Text', sprintf('Angle: %.4g deg', obj.sweepResult.angles(obj.angleIdx)));

            %%% Axis-range side panel -- the shared plot-plus-side-panel body. Unlike
            %%% AngSweepViewerGUI, where angle/velocity never change meaning, this
            %%% class's X/Y meaning changes with curveType, so both the range-field
            %%% DEFAULTS and their LABELS are refreshed on every redraw, not set here.
            %%% classification bar of the whole wavelength sweep -- see class docblock
            obj.classBar = ClassificationBar(mainLayout, obj.sweepResult.sweeps{1}, 2, obj.sweepResult);
            obj.classBar.setReclassifyHandler(@() obj.refreshClassification());

            [obj.mapAxes, sideGl, content] = GuiCommon.plotWithSidePanel(mainLayout, ...
                {18, 26, 18, 26, 10, '1x'}, @(src, event) obj.onMapClicked(src, event));
            content.Layout.Row = 3;

            onRangeEdit = @(~,~) obj.applyAxisRange();
            [obj.xMinField, obj.xMaxField, obj.xRangeLabel] = ...
                GuiCommon.rangeFieldPair(sideGl, 'X min/max', '', '', onRangeEdit);
            [obj.yMinField, obj.yMaxField, obj.yRangeLabel] = ...
                GuiCommon.rangeFieldPair(sideGl, 'Y min/max', '', '', onRangeEdit);

            obj.infoPanel = uilabel(sideGl, 'Text', 'Click the map to select a branch.', ...
                'VerticalAlignment', 'top');
            obj.infoPanel.Layout.Row = 6;

            obj.mapStatusLine = GuiCommon.statusLabel(mainLayout, 4);

            obj.redrawMap();
        end
    end

    %%% ===================================================================
    %%% Public entry point + callbacks -- public (not Access=private),
    %%% matching ModalShapeViewerGUI/AngSweepViewerGUI's own convention so
    %%% tests can call these directly.
    %%% ===================================================================
    methods
        function goTo(obj, wavNo, branchIdx)
        %%% Jump directly to a (wavelength step, branch) selection
        %%% >> goTo(obj, wavNo, branchIdx)
            nWav    = size(obj.branchModeNoMap, 1);
            nBranch = size(obj.branchModeNoMap, 2);
            obj.currentWavNo     = min(max(round(wavNo), 1), nWav);
            obj.currentBranchIdx = min(max(round(branchIdx), 1), nBranch);
            obj.syncControlsAndRedraw();
        end

        function onPrevWav(obj)
            obj.goTo(obj.currentWavNo - 1, obj.currentBranchIdx);
        end

        function onNextWav(obj)
            obj.goTo(obj.currentWavNo + 1, obj.currentBranchIdx);
        end

        function onPrevBranch(obj)
            obj.goTo(obj.currentWavNo, obj.currentBranchIdx - 1);
        end

        function onNextBranch(obj)
            obj.goTo(obj.currentWavNo, obj.currentBranchIdx + 1);
        end

        function onDisplayChanged(obj)
        %%% Switch between the signed component panels and the polarisation-RGB
        %%% panel -- as ModalShapeViewerGUI.onDisplayChanged
        %%% >> onDisplayChanged(obj)
            GuiCommon.setShapeDisplayMode(obj.shapeAxesCell, obj.rgbAxesCell, obj.keyAx, ...
                strcmp(obj.displayDropdown.Value, 'rgb'));
            obj.redrawShape();
        end

        function onCurveTypeChanged(obj)
        %%% Also fired by showGreenBgCheckbox's own ValueChangedFcn (not just
        %%% curveTypeDropdown's) -- both need the exact same redraw +
        %%% Enable-gating response, so one shared handler covers both widgets
        %%% rather than duplicating the logic.
            obj.showGreenBgCheckbox.Enable = obj.hasGreenData && ...
                ismember(obj.curveTypeDropdown.Value, obj.VPHASE_CURVE_TYPES);
            obj.redrawMap();
            obj.syncCrosshair();
        end

        function applyAxisRange(obj)
        %%% Apply the X/Y range fields to the map axes -- same guarded
        %%% parse-and-apply pattern as AngSweepViewerGUI.applyAxisRange.
        %%% Public (not Access=private): same reasoning as every other
        %%% callback in this block, so tests can call it directly.
        %%% >> applyAxisRange(obj)
            GuiCommon.applyFieldLimit(obj.mapAxes, 'x', obj.xMinField, obj.xMaxField);
            GuiCommon.applyFieldLimit(obj.mapAxes, 'y', obj.yMinField, obj.yMaxField);
        end

        function setAxisRange(obj, axisName, lims)
        %%% Set one axis range from code, as if typed into its min/max fields
        %%% >> setAxisRange(obj, axisName, lims)    e.g. setAxisRange(obj, 'y', [NaN 7])
        %%%
        %%% axisName 'x' | 'y' (no colour axis here); lims [lo hi], a NaN
        %%% entry keeps that limit, [] keeps both -- see GuiCommon.setAxisRange
            GuiCommon.setAxisRange(obj, obj.mapAxes, axisName, lims);
        end

        function refreshClassification(obj)
        %%% Redraw after the sweep's classification changed: the bar, the branch
        %%% highlighting (with the axis ranges kept) and the selection's type line
        %%% >> refreshClassification(obj)
            obj.classBar.refresh();
            saved = GuiCommon.captureAxisRange(obj);
            obj.redrawMap();
            GuiCommon.restoreAxisRange(obj, saved);
            obj.syncControlsAndRedraw();
        end

        function onMapClicked(obj, ~, event)
        %%% Snap a map click to the nearest plotted (wavelength step,
        %%% branch) point -- nearest across every branch's line, in
        %%% axis-range-normalised distance (so wildly different x/y scales,
        %%% e.g. wavelength in microns vs. velocity in km/s, don't silently
        %%% bias the match toward whichever axis happens to span more raw
        %%% units) -- GuiCommon.nearestLinePoint owns that search, including
        %%% the "skip hidden branches" rule showRitzBranchesCheckbox needs.
        %%% Same isstruct/isprop IntersectionPoint handling as every other
        %%% viewer, via GuiCommon.clickPoint (tests pass a plain struct
        %%% instead of a real HitEvent object).
        %%% >> onMapClicked(obj, src, event)
            pt = GuiCommon.clickPoint(obj.mapAxes, event);
            [bestWav, bestBranch] = GuiCommon.nearestLinePoint(obj.mapAxes, obj.hLinePlots, pt, ...
                obj.currentWavNo, obj.currentBranchIdx);
            obj.goTo(bestWav, bestBranch);
        end
    end

    %%% ===================================================================
    %%% Run-time helpers
    %%% ===================================================================
    methods (Access = private)
        function syncControlsAndRedraw(obj)
        %%% Shared tail every state-changing method (goTo) funnels through --
        %%% sync every widget, redraw the shape, refresh the info panel and
        %%% crosshair.
        %%% >> syncControlsAndRedraw(obj)
            obj.wavDropdown.Value  = obj.currentWavNo;
            obj.branchField.Value  = obj.currentBranchIdx;

            obj.redrawShape();
            obj.syncCrosshair();
        end

        function redrawShape(obj)
        %%% >> redrawShape(obj)
            try
                ritz   = obj.sweepResult.sweeps{obj.currentWavNo}.cs{obj.angleIdx}.ritzSolver;
                modeNo = obj.branchModeNoMap(obj.currentWavNo, obj.currentBranchIdx);
                dirs = GuiCommon.shapeDirs(obj.displayDropdown);
                if ischar(dirs); axCell = obj.rgbAxesCell; else; axCell = obj.shapeAxesCell; end
                Vis.plotModalShape(ritz, modeNo, 'boosted', [], dirs, [], [], axCell);

                wav  = obj.sweepResult.wavelengths(obj.currentWavNo);
                vel  = ritz.vels(modeNo);
                obj.velLabel.Text = sprintf('vel: %.4f km/s, freq: %.4g', vel, vel/wav);

                %%% classification type of the raw mode index at this
                %%% wavelength step, if it is a significant mode there at all
                typeStr = GuiCommon.significantModeType( ...
                    obj.sweepResult.sweeps{obj.currentWavNo}.significantModesMap{obj.angleIdx}, modeNo);
                obj.infoPanel.Text = { ...
                    sprintf('wavelength = %.4g', wav), ...
                    sprintf('branch # = %d (mode %d)', obj.currentBranchIdx, modeNo), ...
                    sprintf('velocity = %.4f km/s', vel), ...
                    sprintf('type = %s', typeStr)};

                GuiCommon.setStatusOk(obj.shapeStatusLine, '');
            catch err
                GuiCommon.setStatusError(obj.shapeStatusLine, err.message);
            end
        end

        function redrawMap(obj)
        %%% Re-plot the dispersion-curve map for the current curveType
        %%% >> redrawMap(obj)
            cla(obj.mapAxes);
            %%% cla() clears CHILDREN only, not axis-limit modes -- a
            %%% previous call's own applyAxisRange() (below) sets
            %%% XLimMode/YLimMode to 'manual', so without resetting them
            %%% here, reading xlim()/ylim() back after THIS redraw would
            %%% silently return the PREVIOUS curveType's stale range instead
            %%% of the new data's actual extent (AngSweepViewerGUI.onGraphTypeChanged
            %%% resets CLimMode for the same reason).
            obj.mapAxes.XLimMode = 'auto';
            obj.mapAxes.YLimMode = 'auto';
            try
                curveType = obj.curveTypeDropdown.Value;
                showBranches = obj.showRitzBranchesCheckbox.Value;
                %%% Recomputed here from scratch (not read off showGreenBgCheckbox.Enable)
                %%% so this call is correct regardless of whether Enable has been synced
                %%% for the CURRENT curveType yet -- e.g. buildMapPanel's own first
                %%% redrawMap() call happens before onCurveTypeChanged has ever run.
                showBg = obj.hasGreenData && obj.showGreenBgCheckbox.Value && ...
                    ismember(curveType, obj.VPHASE_CURVE_TYPES);
                hLines = Vis.plotDispersionCurves(obj.sweepResult, curveType, obj.angleIdx, showBranches, showBg, [], [], obj.mapAxes);
                for im = 1:numel(hLines)
                    hLines(im).ButtonDownFcn = @(src, event) obj.onMapClicked(src, event);
                end
                obj.hLinePlots = hLines;

                %%% X/Y meaning changes with curveType (unlike
                %%% AngSweepViewerGUI's fixed angle/velocity axes) -- refresh
                %%% both the range-field LABELS and their DEFAULT values from
                %%% this fresh draw's own axis labels/limits. Only when the
                %%% curveType changed: a checkbox toggle (Ritz branches, Green
                %%% background) keeps the same axes, so the user's range stays.
                obj.xRangeLabel.Text = obj.mapAxes.XLabel.String;
                obj.yRangeLabel.Text = obj.mapAxes.YLabel.String;
                if ~strcmp(curveType, obj.drawnCurveType)
                    xl = xlim(obj.mapAxes); yl = ylim(obj.mapAxes);
                    obj.xMinField.Value = num2str(xl(1)); obj.xMaxField.Value = num2str(xl(2));
                    obj.yMinField.Value = num2str(yl(1)); obj.yMaxField.Value = num2str(yl(2));
                    obj.drawnCurveType = curveType;
                end

                GuiCommon.setStatusOk(obj.mapStatusLine, '');
            catch err
                obj.hLinePlots = gobjects(0);
                GuiCommon.setStatusError(obj.mapStatusLine, err.message);
            end
            obj.applyAxisRange();
        end

        function syncCrosshair(obj)
        %%% Mark the current (wavelength, branch) selection on the map --
        %%% read directly off the already-plotted line data (already in
        %%% whatever x/y the active curveType projects to), rather than
        %%% recomputing frequency/velocity a second time.
        %%% >> syncCrosshair(obj)
            im = obj.currentBranchIdx;
            if im < 1 || im > numel(obj.hLinePlots) || ~isgraphics(obj.hLinePlots(im))
                return
            end
            xVal = obj.hLinePlots(im).XData(obj.currentWavNo);
            yVal = obj.hLinePlots(im).YData(obj.currentWavNo);
            obj.hCrosshair = GuiCommon.markCrosshair(obj.mapAxes, obj.hCrosshair, xVal, yVal);
        end
    end
end
