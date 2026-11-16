classdef IsofreqModalShapeViewerGUI < handle
%%% IsofreqModalShapeViewerGUI -- on-demand Green's-response depth-resolved
%%% displacement field viewer for the isoFreq Layered family, the isoFreq
%%% counterpart of ModalShapeViewerGUI's Green's-response source.
%%%
%%% CalcIsofreqSingleLayered/CalcIsofreqAngSweepLayered never run a Ritz
%%% solve (see docs/reference/architecture/forwardCalcs.md, "isoWavevec vs. isoFreq"),
%%% so unlike ModalShapeViewerGUI there is no Ritz modal shape to show and no
%%% sourceDropdown -- Green's response is the ONLY source, so the dropdown is
%%% dropped entirely (absent, not disabled), matching the Plate viewer
%%% family's own precedent. Everything Ritz-specific in ModalShapeViewerGUI
%%% (limModePanel, session-local SAW/limiting overrides, "follow the physical
%%% wave" angle stepping, retention-aware recompute-on-demand) has no
%%% isoFreq analogue and is dropped along with it.
%%%
%%% The underlying mechanism is reused directly: GreenSolver.calc with
%%% nargout>1 (full output), then Vis.calcGreenDispField2D/
%%% Vis.plotGreenModalShape -- exactly what ModalShapeViewerGUI's own
%%% Green's-response mode does, with scanMode='isoFreq', the right omega0,
%%% and kParRange converted to velocity via v = omega0/kParRange (the same
%%% conversion CalcIsofreqSingleLayered.calc() performs). A layered structure
%%% (nLay>=2, the usual isoFreq case) goes through
%%% Vis.plotGreenModalShapeLayered instead, one row per layer, with the same
%%% 'boosted' row heights the Ritz viewers use (see buildShapeAxesGrid). The
%%% isoFreq depth scaling needs nothing extra here: GreenSolver.calc stores
%%% hScaleOut = omega0/v = |k_parallel|, so every depth is already in units of
%%% the selected wave's own wavelength (docs/reference/architecture/green.md).
%%% Any solver error is caught onto the status line, the same try/catch
%%% pattern ModalShapeViewerGUI.redrawShape uses.
%%%
%%% "kParRange" is a continuous, physically-meaningful axis (not a discrete
%%% Ritz mode index), so this class's state is (currentAngNo, currentKParVal,
%%% forceDir) rather than ModalShapeViewerGUI's (currentAngNo, currentModeNo/
%%% currentVel) dual-meaning pair -- there is nothing to switch between.
%%% Angle stepping keeps the SAME kParRange across the step, the direct
%%% analogue of ModalShapeViewerGUI's "same velocity, different angle"
%%% behaviour in its own Green's-response mode.
%%%
%%% ACCEPTS EITHER CalcIsofreqAngSweepLayered OR CalcIsofreqSingleLayered
%%% directly, same "adapt the single-point case into a degenerate 1-angle
%%% struct" pattern ModalShapeViewerGUI uses for CalcSingleHS/
%%% CalcSingleLayered. IMPORTANT DIFFERENCE: CalcIsofreqSingleLayered.structure
%%% is ALREADY ROTATED to its own angle (constructor: CalcCommon.rotateStructure(structure,
%%% angle)), whereas CalcIsofreqAngSweepLayered.structure is the UNROTATED base
%%% (each element rotates its own internal copy). GreenSolver.calc rotates its
%%% `structure` argument by whatever `angle` it is given, UNCONDITIONALLY --
%%% this is the rotation rule documented in docs/reference/architecture/
%%% forwardCalcs.md, "Rotation convention: pass 0 to solvers after
%%% construction", and the reason
%%% CalcIsofreqSingleLayered.calc() itself passes angle=0, not obj.angle, to
%%% GreenSolver.calc. So redrawShape below passes angle=0 when built from a
%%% single-point adapter (already rotated) and the real angleVal when built
%%% from a genuine sweep (unrotated), tracked via the private
%%% structureAlreadyRotated flag set at construction.
%%%
%%% See docs/reference/architecture/gui_viewers.md for the full isoFreq viewer family.
%%%
%%% GuliweRR v6.0

    properties
        fig
        shapeFig

        angleDropdown, prevAngleBtn, nextAngleBtn
        kParField, prevKParBtn, nextKParBtn, velLabel
        forceDirLabel, forceDirDropdown
        displayDropdown     % 'components' | 'rgb' -- GuiCommon.shapeDisplayDropdown
        statusLine

        OnSelectionChanged = []   % function handle @(angle, kParVal), set by the map/curve viewer
    end

    properties (SetAccess = private)
        sweepResult       % struct(angles, omega0, structure, kParRange) -- either a
                           % genuine CalcIsofreqAngSweepLayered or a degenerate 1-angle
                           % adapter built from a CalcIsofreqSingleLayered (see class docblock)
        currentAngNo
        currentKParVal
    end

    properties (Access = {?IsofreqModalShapeViewerGUI, ?test_IsofreqModalShapeViewerGUI})
        embedded  = false     % true when rendering via the uiaxes grid (real embedding
                               % OR standalone 'combined'); false only for standalone 'split'
        shapeAxesCell = {}    % {3 x nLay} pre-built uiaxes, whenever obj.embedded is true
                               % (one row per layer, layer 1 = substrate at the bottom)
        rgbAxesCell   = {}    % {1 x nLay} uiaxes of the polarisation-RGB display (embedded only)
        keyAx         = []    % its colour-key uiaxes
        structureAlreadyRotated = false   % see class docblock's "IMPORTANT DIFFERENCE" note
        %%% Friend access to the test class -- matches ModalShapeViewerGUI.
        %%% fullRitzCache's own precedent -- so a test can confirm the
        %%% double-rotation-avoidance logic directly (structureAlreadyRotated)
        %%% and read the rendered field back (shapeAxesCell) without either
        %%% needing to be part of this class' own public surface.
    end

    methods
        function obj = IsofreqModalShapeViewerGUI(sweepResult, windowMode, controlsParent, shapeParent)
        %%% Construct the isoFreq modal-shape/Green's-response viewer window
        %%% >> obj = IsofreqModalShapeViewerGUI(sweepResult, windowMode, controlsParent, shapeParent)
        %%%
        %%% INPUT: 1. sweepResult     .. CalcIsofreqAngSweepLayered, already
        %%%                              calc()'d -- OR CalcIsofreqSingleLayered,
        %%%                              already calc()'d (packaged internally as
        %%%                              a degenerate 1-angle "sweep", see class
        %%%                              docblock)
        %%%  (opt) 2. windowMode      .. 'combined' (default, one window) or
        %%%                              'split' (two windows), same meaning as
        %%%                              ModalShapeViewerGUI's own. Ignored when embedding.
        %%%  (opt) 3. controlsParent  .. a uigridlayout region to build the
        %%%                              control panel into instead of a fresh
        %%%                              own uifigure
        %%%  (opt) 4. shapeParent     .. a uigridlayout region to build the
        %%%                              shape uiaxes grid into. Required together
        %%%                              with controlsParent (both or neither).
            validateattributes(sweepResult, ...
                {'CalcIsofreqAngSweepLayered', 'CalcIsofreqSingleLayered'}, ...
                {'scalar'}, 'IsofreqModalShapeViewerGUI', 'sweepResult');
            if isa(sweepResult, 'CalcIsofreqSingleLayered')
                if isempty(sweepResult.absG)
                    error('IsofreqModalShapeViewerGUI:constructor:notCalculated', ...
                        'sweepResult has not been calc()''d yet -- nothing to view.');
                end
                singleResult = sweepResult;
                sweepResult = struct( ...
                    'angles', singleResult.angle, ...
                    'omega0', singleResult.omega0, ...
                    'structure', singleResult.structure, ...
                    'kParRange', singleResult.kParRange);
                obj.structureAlreadyRotated = true;
            else
                if isempty(sweepResult.absGMap)
                    error('IsofreqModalShapeViewerGUI:constructor:notCalculated', ...
                        'sweepResult has not been calc()''d yet -- nothing to view.');
                end
                obj.structureAlreadyRotated = false;
            end
            if nargin < 2 || isempty(windowMode); windowMode = 'combined'; end
            if nargin < 3; controlsParent = []; end
            if nargin < 4; shapeParent = []; end
            GuiCommon.validateWindowMode('IsofreqModalShapeViewerGUI', windowMode, {'combined', 'split'});
            GuiCommon.validateEmbedParents('IsofreqModalShapeViewerGUI', windowMode, controlsParent, shapeParent);
            obj.sweepResult = sweepResult;
            obj.embedded    = ~isempty(controlsParent);

            if isempty(controlsParent) && strcmp(windowMode, 'combined')
                ownFig = uifigure('Name', 'GuliweRR -- isoFreq Modal Shape Viewer', ...
                    'Position', [80 60 1200 360], 'Visible', 'off');
                [shapeParent, controlsParent] = GuiCommon.shapeControlsSplit(ownFig, [4 4 4 4]);
                obj.embedded = true;
            end

            if obj.embedded
                obj.buildShapeAxesGrid(shapeParent);
                obj.shapeFig = ancestor(shapeParent, 'figure');
            else
                obj.shapeFig = figure('Visible', 'off', ...
                    'Name', 'GuliweRR -- isoFreq Green''s Response');
            end
            obj.buildFigure(controlsParent);
            obj.goToPoint(sweepResult.angles(1), sweepResult.kParRange(1));
        end

        function show(obj)
        %%% Make both windows visible (constructed invisible so headless/test
        %%% use is safe) -- a no-op on the shapeFig side beyond the first call
        %%% when embedded, since .fig and .shapeFig are then the same figure.
        %%% >> show(obj)
            obj.fig.Visible = 'on';
            obj.shapeFig.Visible = 'on';
        end
    end

    %%% ===================================================================
    %%% Figure construction
    %%% ===================================================================
    methods (Access = private)
        function buildShapeAxesGrid(obj, shapeParent)
        %%% Pre-build the {3 x nLay} uiaxes grid Vis.plotGreenModalShapeLayered
        %%% (or, for a half-space, Vis.plotGreenModalShape) renders into via
        %%% hAxCellIn, through the same GuiCommon.shapeAxesGrid every other
        %%% shape viewer uses. Its row heights are the 'boosted' weights, so
        %%% redrawShape must render with substrateCutoff='boosted' as well. The
        %%% structure never changes for the life of the viewer, so the grid is
        %%% built once.
        %%% >> buildShapeAxesGrid(obj, shapeParent)
            relDepths = [obj.sweepResult.structure.calcParams.relDepth];
            [obj.shapeAxesCell, obj.rgbAxesCell, obj.keyAx] = GuiCommon.shapeAxesGrid(shapeParent, relDepths);
        end

        function buildFigure(obj, controlsParent)
            [obj.fig, mainLayout] = GuiCommon.hostLayout(controlsParent, ...
                'GuliweRR -- isoFreq Modal Shape Viewer', [100 100 620 300], [5 1]);
            mainLayout.RowHeight = {36, 36, 36, 36, 26};

            %%% row 1 -- angle controls
            angleStrs = arrayfun(@(a) sprintf('%.1f', a), obj.sweepResult.angles, 'UniformOutput', false);
            [obj.angleDropdown, obj.prevAngleBtn, obj.nextAngleBtn] = ...
                GuiCommon.dropdownStepperRow(mainLayout, 1, 'Angle:', angleStrs, {50, 140, 36, 36}, ...
                    @(~,~) obj.stepToAngle(obj.angleDropdown.Value), ...
                    @(~,~) obj.onPrevAngle(), @(~,~) obj.onNextAngle());

            %%% row 2 -- kParRange field (the sole continuous state -- no
            %%% Ritz/velocity dual meaning to switch between, unlike
            %%% ModalShapeViewerGUI's re-purposed modeField). << / >> step to
            %%% the neighbouring kParRange sample (stepKPar); the trailing label
            %%% shows the matching phase velocity v = omega0/kParRange.
            [obj.kParField, obj.prevKParBtn, obj.nextKParBtn, obj.velLabel] = ...
                GuiCommon.indexStepperRow(mainLayout, 2, 'kParRange:', {70, 70, 36, 36, '1x'}, ...
                    @(~,~) obj.onKParFieldChanged(), ...
                    @(~,~) obj.stepKPar(-1), @(~,~) obj.stepKPar(+1));
            obj.kParField.Limits = [0 Inf];
            obj.kParField.RoundFractionalValues = 'off';
            obj.kParField.Value = obj.sweepResult.kParRange(1);
            obj.velLabel.FontWeight = 'bold';

            %%% row 3 -- force direction (always visible -- Green's response
            %%% is the only source, see class docblock)
            forceRow = uigridlayout(mainLayout, [1 2]);
            forceRow.Layout.Row = 3;
            forceRow.ColumnWidth = {70, 160};
            forceRow.Padding = [0 0 0 0];
            obj.forceDirLabel = uilabel(forceRow, 'Text', 'Force dir:');
            obj.forceDirDropdown = uidropdown(forceRow, ...
                'Items', {'x (propag. dir.)', 'y (in-plane)', 'z (out-of-plane)'}, ...
                'ItemsData', 1:3, 'Value', 1, ...
                'ValueChangedFcn', @(~,~) obj.onForceDirChanged());

            %%% row 4 -- display: signed components / polarisation RGB
            displayRow = uigridlayout(mainLayout, [1 2]);
            displayRow.Layout.Row = 4;
            displayRow.ColumnWidth = {70, 200};
            displayRow.Padding = [0 0 0 0];
            uilabel(displayRow, 'Text', 'Display:');
            obj.displayDropdown = GuiCommon.shapeDisplayDropdown(displayRow, @(~,~) obj.onDisplayChanged());

            obj.statusLine = GuiCommon.statusLabel(mainLayout, 5);
        end
    end

    %%% ===================================================================
    %%% Public entry point + callbacks -- public (not Access=private):
    %%% StructureControlGUI's own convention, so tests can call these
    %%% directly, and so IsofreqAngSweepViewerGUI/IsofreqSingleLayeredViewerGUI's
    %%% OnPointSelected can be wired straight to goToPoint.
    %%% ===================================================================
    methods
        function goToPoint(obj, angleVal, kParVal)
        %%% Jump directly to a specific (angle, kParRange) point
        %%% >> goToPoint(obj, angleVal, kParVal)
            [~, angNo] = min(abs(obj.sweepResult.angles - angleVal));
            obj.currentAngNo   = angNo;
            obj.currentKParVal = kParVal;
            obj.syncControlsAndRedraw();
        end

        function stepToAngle(obj, newAngNo)
        %%% Keeps the SAME kParRange across the angle step -- see class docblock
        %%% >> stepToAngle(obj, newAngNo)
            obj.goToPoint(obj.sweepResult.angles(newAngNo), obj.currentKParVal);
        end

        function onPrevAngle(obj)
            obj.stepToAngle(max(obj.currentAngNo - 1, 1));
        end

        function onNextAngle(obj)
            obj.stepToAngle(min(obj.currentAngNo + 1, numel(obj.sweepResult.angles)));
        end

        function stepKPar(obj, direction)
        %%% Step to the neighbouring kParRange sample (direction -1 / +1),
        %%% clamped at the ends -- the isoFreq counterpart of
        %%% ModalShapeViewerGUI's fixed VEL_STEP in Green's-response mode, on the
        %%% map's own k grid. A value typed between two samples steps to the
        %%% nearer side's next sample.
        %%% >> stepKPar(obj, direction)
            grid = sort(obj.sweepResult.kParRange(:));
            tol  = 1e-9 * max(abs(grid));
            if direction > 0
                k = find(grid > obj.currentKParVal + tol, 1, 'first');
                if isempty(k), k = numel(grid); end
            else
                k = find(grid < obj.currentKParVal - tol, 1, 'last');
                if isempty(k), k = 1; end
            end
            obj.goToPoint(obj.sweepResult.angles(obj.currentAngNo), grid(k));
        end

        function onKParFieldChanged(obj)
            obj.goToPoint(obj.sweepResult.angles(obj.currentAngNo), obj.kParField.Value);
        end

        function onForceDirChanged(obj)
            obj.redrawShape();
        end

        function onDisplayChanged(obj)
        %%% Switch between the signed component panels and the polarisation-RGB
        %%% panel -- as ModalShapeViewerGUI.onDisplayChanged
        %%% >> onDisplayChanged(obj)
            if obj.embedded
                GuiCommon.setShapeDisplayMode(obj.shapeAxesCell, obj.rgbAxesCell, obj.keyAx, ...
                    strcmp(obj.displayDropdown.Value, 'rgb'));
            end
            obj.redrawShape();
        end
    end

    %%% ===================================================================
    %%% Run-time helpers
    %%% ===================================================================
    methods (Access = private)
        function syncControlsAndRedraw(obj)
        %%% Shared tail every state-changing method (goToPoint) funnels
        %%% through: sync every widget to (currentAngNo, currentKParVal),
        %%% redraw the shape, notify OnSelectionChanged.
        %%% >> syncControlsAndRedraw(obj)
            obj.angleDropdown.Value = obj.currentAngNo;
            obj.kParField.Value = obj.currentKParVal;
            obj.velLabel.Text = sprintf('%.3f km/s', obj.sweepResult.omega0 / obj.currentKParVal);

            obj.redrawShape();

            if ~isempty(obj.OnSelectionChanged)
                obj.OnSelectionChanged(obj.sweepResult.angles(obj.currentAngNo), obj.currentKParVal);
            end
        end

        function redrawShape(obj)
        %%% Re-render the shape display for the current angle/kParRange/
        %%% force-dir selection -- into shapeAxesCell (obj.embedded) or the
        %%% classic shapeFig (windowMode='split'). See class docblock for
        %%% the angleForGreen/structureAlreadyRotated reasoning.
        %%% >> redrawShape(obj)
            try
                angleVal = obj.sweepResult.angles(obj.currentAngNo);
                forceDir = obj.forceDirDropdown.Value;
                velForGreen = obj.sweepResult.omega0 / obj.currentKParVal;
                if obj.structureAlreadyRotated
                    angleForGreen = 0;
                else
                    angleForGreen = angleVal;
                end
                [greenObj, ~] = GreenSolver.calc(obj.sweepResult.structure, angleForGreen, ...
                    velForGreen, forceDir, 3, 'ScanMode', 'isoFreq', 'Omega0', obj.sweepResult.omega0, ...
                    'ComputeMinus', false);
                dirs = GuiCommon.shapeDirs(obj.displayDropdown);
                if ischar(dirs); axCell = obj.rgbAxesCell; else; axCell = obj.shapeAxesCell; end
                if ~obj.embedded; axCell = []; end
                str = obj.sweepResult.structure;
                if str.nLay > 1
                    Vis.plotGreenModalShapeLayered(greenObj, str, forceDir, 'boosted', [], dirs, [], ...
                        obj.shapeFig, axCell);
                elseif obj.embedded
                    Vis.plotGreenModalShape(greenObj, forceDir, [], dirs, [], [], axCell);
                else
                    Vis.plotGreenModalShape(greenObj, forceDir, [], dirs, [], obj.shapeFig);
                end
                obj.shapeFig.Visible = obj.fig.Visible;   % follow the control window's own visibility
                GuiCommon.setStatusOk(obj.statusLine, '');
            catch err
                GuiCommon.setStatusError(obj.statusLine, err.message);
            end
        end
    end
end
