classdef DispersionViewerGUIPlate < handle
%%% DispersionViewerGUIPlate -- fixed-angle wavevector/frequency dispersion
%%% viewer for a free-plate (Lamb-wave) 2D sweep. See
%%% docs/reference/architecture/gui_viewers.md and docs/reference/architecture/plate.md.
%%%
%%% A SEPARATE class from DispersionViewerGUI (the CalcWavevecAngSweep2Layered
%%% analogue), not a class(sweepResult) dispatch branch inside it — the same
%%% parallel-sibling split Plate and Layered/HS take everywhere their data
%%% models genuinely differ. Differences, all following from the plate physics
%%% (no Green's function, every mode real):
%%%   - No Green's-response toggle at all: absent, not gated, since no analogue
%%%     exists for a free plate.
%%%   - No significantModesMap-based branch highlighting (Plate has no
%%%     "significant subset" concept) — uses CalcWavevecAngSweep2Plate's
%%%     classifyBranches instead, via Vis.plotDispersionCurvesPlate.
%%%   - Adds a Christoffel-bounds checkbox and a colour-mode selector
%%%     ('type'/'branch'), both capabilities Vis.plotDispersionCurvesPlate has
%%%     and Vis.plotDispersionCurves does not.
%%%   - The wavevector stepper reads sweepResult.wavVecs directly, not a real
%%%     wavelength.
%%%   - The modal-shape frequency label reads ritz.freqsPhysical directly, not a
%%%     re-derived vel/wavelength ratio.
%%%
%%% Everything else — windowMode ('combined'/'split'), embedding via
%%% parentContainer, the {3 x nLay} shape-axes grid, axis-range fields,
%%% crosshair/click-to-select, branch stepping via branchModeNoMap — follows
%%% DispersionViewerGUI's mechanisms, none of which are Layered-specific.
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
        colorModeDropdown
        showChristoffelBoundsCheckbox
        xRangeLabel, yRangeLabel
        xMinField, xMaxField, yMinField, yMaxField
        infoPanel
        mapStatusLine
    end

    properties (SetAccess = private)
        sweepResult       % CalcWavevecAngSweep2Plate, already calc()'d
        windowMode        % 'combined' (default) | 'split' -- see class docblock
        angleIdx          % index into sweepResult.angles this view is fixed at (default 1)
        branchModeNoMap   % [nWav x nCompModes], DispersionBranches.track's own output
        currentWavNo      % index into sweepResult.wavVecs
        currentBranchIdx  % column index into branchModeNoMap (a tracked physical branch)
    end

    properties (Access = private)
        shapeAxesCell = {}   % {3 x nLay} pre-built uiaxes, same pattern as DispersionViewerGUI's own
        rgbAxesCell   = {}   % {1 x nLay} uiaxes of the polarisation-RGB display
        keyAx         = []   % its colour-key uiaxes
        hLinePlots = gobjects(0)   % current map's line handles, one per branch
        drawnCurveType = ''    % curveType of the last redrawMap -- range fields reset only when it changes
        hCrosshair = gobjects(0)
    end

    methods
        function obj = DispersionViewerGUIPlate(sweepResult, windowMode, angleIdx, parentContainer)
        %%% Construct the Plate dispersion viewer
        %%% >> obj = DispersionViewerGUIPlate(sweepResult, windowMode, angleIdx, parentContainer)
        %%%
        %%% INPUT: 1. sweepResult .. CalcWavevecAngSweep2Plate, already
        %%%                          calc()'d, with retentionLevel='full'
        %%%                          (needed for
        %%%                          DispersionBranches.track
        %%%                          / classifyBranches -- their own
        %%%                          missingEigVecs error surfaces directly if
        %%%                          this isn't satisfied)
        %%%  (opt) 2. windowMode  .. 'combined' (default) | 'split'
        %%%  (opt) 3. angleIdx    .. which angle (default 1)
        %%%  (opt) 4. parentContainer .. a uigridlayout region to build into
        %%%                          instead of a fresh own uifigure
            validateattributes(sweepResult, {'CalcWavevecAngSweep2Plate'}, ...
                {'scalar'}, 'DispersionViewerGUIPlate', 'sweepResult');
            if nargin < 2 || isempty(windowMode); windowMode = 'combined'; end
            GuiCommon.validateWindowMode('DispersionViewerGUIPlate', windowMode, {'combined', 'split'});
            if nargin < 3 || isempty(angleIdx); angleIdx = 1; end
            validateattributes(angleIdx, {'double'}, ...
                {'scalar', 'integer', 'positive', '<=', numel(sweepResult.angles)}, ...
                'DispersionViewerGUIPlate', 'angleIdx');
            if nargin < 4; parentContainer = []; end
            GuiCommon.requireStandaloneForSplit('DispersionViewerGUIPlate', windowMode, ...
                parentContainer, 'parentContainer');
            if isempty(sweepResult.sweeps)
                error('DispersionViewerGUIPlate:constructor:notCalculated', ...
                    'sweepResult has not been calc()''d yet -- nothing to view.');
            end

            obj.sweepResult     = sweepResult;
            obj.windowMode      = windowMode;
            obj.angleIdx        = angleIdx;
            obj.branchModeNoMap = DispersionBranches.track(sweepResult, angleIdx);

            obj.buildFigures(parentContainer);
            obj.goTo(1, 1);
        end

        function show(obj)
        %%% Make every window visible (constructed invisible so headless/
        %%% test use is safe).
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
                obj.fig = uifigure('Name', 'GuliweRR -- Plate Dispersion Viewer: Modal Shapes', ...
                    'Position', [80 460 1400 440], 'Visible', 'off');
                [shapeParent, controlsParent] = GuiCommon.shapeControlsSplit(obj.fig, [4 4 4 4]);

                obj.mapFig = uifigure('Name', 'GuliweRR -- Plate Dispersion Viewer: Map', ...
                    'Position', [80 60 1000 350], 'Visible', 'off');
                mapParent = uigridlayout(obj.mapFig, [1 1]);
                mapParent.Padding = [4 4 4 4];
            else
                %%% 'combined' standalone AND embedded build the SAME
                %%% shapes+controls-over-map 2-row layout -- see
                %%% DispersionViewerGUI.buildFigures' own comment.
                [obj.fig, mainLayout] = GuiCommon.hostLayout(parentContainer, ...
                    'GuliweRR -- Plate Dispersion Viewer', [80 60 1200 900], [2 1]);
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
        %%% GuiCommon.shapeAxesGrid (see that helper's docblock for why a
        %%% plain uipanel, not a uigridlayout).
            relDepths = [obj.sweepResult.structure.calcParams.relDepth];
            [obj.shapeAxesCell, obj.rgbAxesCell, obj.keyAx] = GuiCommon.shapeAxesGrid(shapeParent, relDepths);
        end

        function buildShapeControls(obj, controlsParent)
            mainLayout = uigridlayout(controlsParent, [5 1]);
            mainLayout.RowHeight = {36, 36, 36, '1x', 26};
            mainLayout.Padding = [4 4 4 4];

            wavStrs = arrayfun(@(w) sprintf('%.4g', w), obj.sweepResult.wavVecs, 'UniformOutput', false);
            [obj.wavDropdown, obj.prevWavBtn, obj.nextWavBtn] = ...
                GuiCommon.dropdownStepperRow(mainLayout, 1, 'Wavevector:', wavStrs, {80, 120, 36, 36}, ...
                    @(~,~) obj.goTo(obj.wavDropdown.Value, obj.currentBranchIdx), ...
                    @(~,~) obj.onPrevWav(), @(~,~) obj.onNextWav());

            [obj.branchField, obj.prevBranchBtn, obj.nextBranchBtn, obj.velLabel] = ...
                GuiCommon.indexStepperRow(mainLayout, 2, 'Branch:', {60, 60, 36, 36, '1x'}, ...
                    @(~,~) obj.goTo(obj.currentWavNo, obj.branchField.Value), ...
                    @(~,~) obj.onPrevBranch(), @(~,~) obj.onNextBranch());

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
            mainLayout = uigridlayout(mapParent, [3 1]);
            mainLayout.RowHeight = {36, '1x', 26};
            mainLayout.Padding = [4 4 4 4];

            topBar = uigridlayout(mainLayout, [1 7]);
            topBar.Layout.Row = 1;
            topBar.ColumnWidth = {90, 200, 90, 130, 180, 20, 140};
            topBar.Padding = [0 0 0 0];
            uilabel(topBar, 'Text', 'Curve type:');
            obj.curveTypeDropdown = uidropdown(topBar, ...
                'Items', {'wavelength-freq', 'wavelength-vphase', 'k-freq', 'k-vphase', 'freq-vphase', 'freq-vgroup'}, ...
                'ValueChangedFcn', @(~,~) obj.onCurveTypeChanged());

            uilabel(topBar, 'Text', 'Colour by:');
            obj.colorModeDropdown = uidropdown(topBar, 'Items', {'type', 'branch'}, ...
                'ValueChangedFcn', @(~,~) obj.onCurveTypeChanged());

            %%% Christoffel bounds (docs/reference/architecture/plate.md: "Christoffel
            %%% limiting velocities are still relevant" despite no Green
            %%% data) -- unconditionally available (no Enable-gating the way
            %%% the Layered viewer's Green checkbox needs, since Christoffel
            %%% data costs nothing extra to compute and is not curveType- or
            %%% calcMode-dependent the way Green data is); Vis.plotDispersionCurvesPlate
            %%% itself silently skips drawing on a non-velocity-y-axis
            %%% curveType rather than erroring, so no gating is needed here either.
            obj.showChristoffelBoundsCheckbox = uicheckbox(topBar, 'Text', 'Show Christoffel bounds', ...
                'Value', true, 'ValueChangedFcn', @(~,~) obj.onCurveTypeChanged());

            uilabel(topBar, 'Text', sprintf('Angle: %.4g deg', obj.sweepResult.angles(obj.angleIdx)));

            [obj.mapAxes, sideGl, content] = GuiCommon.plotWithSidePanel(mainLayout, ...
                {18, 26, 18, 26, 10, '1x'}, @(src, event) obj.onMapClicked(src, event));
            content.Layout.Row = 2;

            onRangeEdit = @(~,~) obj.applyAxisRange();
            [obj.xMinField, obj.xMaxField, obj.xRangeLabel] = ...
                GuiCommon.rangeFieldPair(sideGl, 'X min/max', '', '', onRangeEdit);
            [obj.yMinField, obj.yMaxField, obj.yRangeLabel] = ...
                GuiCommon.rangeFieldPair(sideGl, 'Y min/max', '', '', onRangeEdit);

            obj.infoPanel = uilabel(sideGl, 'Text', 'Click the map to select a branch.', ...
                'VerticalAlignment', 'top');
            obj.infoPanel.Layout.Row = 6;

            obj.mapStatusLine = GuiCommon.statusLabel(mainLayout, 3);

            obj.redrawMap();
        end
    end

    %%% ===================================================================
    %%% Public entry point + callbacks
    %%% ===================================================================
    methods
        function goTo(obj, wavNo, branchIdx)
        %%% Jump directly to a (wavevector step, branch) selection
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
        %%% Also fired by colorModeDropdown/showChristoffelBoundsCheckbox's
        %%% own ValueChangedFcn -- all three need the same redraw response.
            obj.redrawMap();
            obj.syncCrosshair();
        end

        function applyAxisRange(obj)
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

        function onMapClicked(obj, ~, event)
        %%% Snap a map click to the nearest plotted (wavevector step,
        %%% branch) point -- identical mechanism to DispersionViewerGUI's own,
        %%% both through GuiCommon.clickPoint/nearestLinePoint.
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

                wavVec = obj.sweepResult.wavVecs(obj.currentWavNo);
                vel    = ritz.vels(modeNo);
                freq   = ritz.freqsPhysical(modeNo);   % RitzSolverLamb's own convention -- NOT
                                                         % re-derived here, see class docblock
                obj.velLabel.Text = sprintf('vel: %.4f km/s, freq: %.4g', vel, freq);

                %%% names from Vis.plateModeTypeStyle -- the one source for the
                %%% plate mode-type vocabulary, shared with the plot this label
                %%% sits next to (was a second hand-maintained copy here)
                typeNames = Vis.plateModeTypeStyle();
                typeStr = sprintf('%s (confidence=%.2f)', ...
                    typeNames{ritz.modesType(modeNo)}, ritz.modeTypeConfidence(modeNo));
                obj.infoPanel.Text = { ...
                    sprintf('wavevector = %.4g', wavVec), ...
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
            obj.mapAxes.XLimMode = 'auto';
            obj.mapAxes.YLimMode = 'auto';
            try
                curveType   = obj.curveTypeDropdown.Value;
                colorMode   = obj.colorModeDropdown.Value;
                showChrist  = obj.showChristoffelBoundsCheckbox.Value;
                hLines = Vis.plotDispersionCurvesPlate(obj.sweepResult, curveType, obj.angleIdx, ...
                    colorMode, [], showChrist, obj.mapAxes);
                for im = 1:numel(hLines)
                    hLines(im).ButtonDownFcn = @(src, event) obj.onMapClicked(src, event);
                end
                obj.hLinePlots = hLines;

                obj.xRangeLabel.Text = obj.mapAxes.XLabel.String;
                obj.yRangeLabel.Text = obj.mapAxes.YLabel.String;
                %%% Reset the range fields only on a curveType change; the colour
                %%% mode and Christoffel-bounds toggles keep the user's range.
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
