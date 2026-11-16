classdef ModalShapeViewerGUIPlate < handle
%%% ModalShapeViewerGUIPlate -- modal-shape viewer window for the Plate
%%% (Lamb-wave) family. Ritz-mode-shape counterpart of ModalShapeViewerGUI.
%%%
%%% A DELIBERATELY SEPARATE class from ModalShapeViewerGUI, not a
%%% class-dispatch branch inside it -- same "separate parallel sibling"
%%% convention as every other Plate GUI class (see DispersionViewerGUIPlate's
%%% own docblock). Differences from ModalShapeViewerGUI, all following
%%% directly from CalcSinglePlate.m's own docblock (no Green analogue exists
%%% for a free plate; every Ritz mode is a genuine physical branch, not a
%%% "significant subset" of one):
%%%   - No source toggle, no force-direction control, no continuous-velocity
%%%     Green's-response mode -- ABSENT, not disabled (DispersionViewerGUIPlate's
%%%     own precedent). Every control here drives the one and only source,
%%%     the Ritz modal shape.
%%%   - No significantModesMap/"SAW vs limiting" concept, so no session-local
%%%     "mark as shown" override panel with per-entry Go-To buttons -- Plate
%%%     has nothing analogous to classify a mode INTO (every mode already IS
%%%     a genuine physical branch). Replaced by a single status label showing
%%%     the CURRENTLY DISPLAYED mode's type/confidence
%%%     (RitzSolverLamb.modesType/modeTypeConfidence), the same text
%%%     DispersionViewerGUIPlate.redrawShape already renders.
%%%   - No followPhysicalWave/classifyCurrentMode angle-stepping logic.
%%%     ModalShapeViewerGUI's version specifically follows a SAW/limiting
%%%     TYPE through significantModesMap across an angle step -- Plate has no
%%%     such classification to follow. Angle stepping here keeps the SAME RAW
%%%     MODE INDEX instead (clamped to the new angle's nCompModes) -- a plain
%%%     default rather than an unrequested mode-type-tracking heuristic; if a
%%%     mode's physical identity needs following across angle (mode ordering
%%%     by velocity can shuffle, same as SAW), the user re-selects by mode
%%%     type/confidence directly.
%%%
%%% Retention-aware recompute-on-demand (mirrors ModalShapeViewerGUI's own
%%% resolveDisplayRitz exactly, RitzSolverLamb in place of RitzSolverSAW):
%%% CalcSinglePlate's retentionLevel is 'full'/'none' only (no partial-prune
%%% 'saw_only' tier to worry about column-index mismatches for) -- 'none'
%%% clears eigVecs_con entirely via pruneEigVecs([]), so this class
%%% recomputes and caches per angle whenever retentionLevel~='full', exactly
%%% as ModalShapeViewerGUI already does.
%%%
%%% GuliweRR v6.0

    properties
        fig
        shapeFig

        angleDropdown, prevAngleBtn, nextAngleBtn
        modeField, prevModeBtn, nextModeBtn, velLabel
        modeTypeLabel
        displayDropdown     % 'components' | 'rgb' -- GuiCommon.shapeDisplayDropdown
        statusLine

        OnSelectionChanged = []   % function handle @(angle, vel), fired on every state change
    end

    properties (SetAccess = private)
        sweepResult      % struct(angles, cs, structure) -- .cs{iAngle} are CalcSinglePlate
                          % objects, built directly (CalcAngSweepPlate already exposes exactly
                          % this shape) or adapted from a bare CalcSinglePlate (see constructor)
        windowMode        % 'combined' (default, one window) | 'split' (two windows)
        currentAngNo      % index into sweepResult.angles
        currentModeNo     % active Ritz mode index
        currentVel        % ritz.vels(currentModeNo)
    end

    properties (Access = private)
        embedded  = false
        shapeAxesCell = {}
        rgbAxesCell   = {}   % {1 x nLay} uiaxes of the polarisation-RGB display (embedded only)
        keyAx         = []   % its colour-key uiaxes
    end

    properties (Access = {?ModalShapeViewerGUIPlate, ?test_ModalShapeViewerGUIPlate})
        fullRitzCache = {}   % {1 x nAngles}, lazily filled -- see resolveDisplayRitz
    end

    methods
        function obj = ModalShapeViewerGUIPlate(sweepResult, windowMode, controlsParent, shapeParent)
        %%% Construct the Plate modal-shape viewer window
        %%% >> obj = ModalShapeViewerGUIPlate(sweepResult, windowMode, controlsParent, shapeParent)
        %%%
        %%% INPUT: 1. sweepResult     .. CalcAngSweepPlate, already calc()'d --
        %%%                              OR CalcSinglePlate, already calc()'d
        %%%                              (packaged internally as a degenerate
        %%%                              1-angle "sweep", see below)
        %%%  (opt) 2. windowMode      .. 'combined' (default, one window) or
        %%%                              'split' (two windows)
        %%%  (opt) 3. controlsParent  .. a uigridlayout region to build the
        %%%                              control panel into instead of a fresh
        %%%                              own uifigure
        %%%  (opt) 4. shapeParent     .. a uigridlayout region to build the
        %%%                              modal-shape uiaxes grid into. Required
        %%%                              together with controlsParent.
            validateattributes(sweepResult, {'CalcAngSweepPlate', 'CalcSinglePlate'}, ...
                {'scalar'}, 'ModalShapeViewerGUIPlate', 'sweepResult');
            if isa(sweepResult, 'CalcSinglePlate')
                if isempty(sweepResult.ritzSolver) || sweepResult.ritzSolver.nCompModes == 0
                    error('ModalShapeViewerGUIPlate:constructor:notCalculated', ...
                        'sweepResult has not been calc()''d yet -- nothing to view.');
                end
                sweepResult = struct( ...
                    'angles', sweepResult.angle, ...
                    'cs', {{sweepResult}}, ...
                    'structure', sweepResult.structure);
            elseif isempty(sweepResult.velsMap)
                %%% .cs is populated by the CONSTRUCTOR itself (uncalculated
                %%% CalcSinglePlate objects) -- not a valid "has calc() run"
                %%% signal; velsMap stays [] until assembleMaps() runs at the
                %%% end of calc() (same pitfall AngSweepViewerGUIPlate's own
                %%% constructor avoids the same way).
                error('ModalShapeViewerGUIPlate:constructor:notCalculated', ...
                    'sweepResult has not been calc()''d yet -- nothing to view.');
            end
            if nargin < 2 || isempty(windowMode); windowMode = 'combined'; end
            if nargin < 3; controlsParent = []; end
            if nargin < 4; shapeParent = []; end
            GuiCommon.validateWindowMode('ModalShapeViewerGUIPlate', windowMode, {'combined', 'split'});
            GuiCommon.validateEmbedParents('ModalShapeViewerGUIPlate', windowMode, controlsParent, shapeParent);
            obj.sweepResult   = sweepResult;
            obj.windowMode    = windowMode;
            obj.embedded      = ~isempty(controlsParent);
            obj.fullRitzCache = cell(1, numel(sweepResult.angles));

            if isempty(controlsParent) && strcmp(windowMode, 'combined')
                ownFig = uifigure('Name', 'GuliweRR -- Modal Shape Viewer (Plate)', ...
                    'Position', [80 60 1400 640], 'Visible', 'off');
                [shapeParent, controlsParent] = GuiCommon.shapeControlsSplit(ownFig, [4 4 4 4]);
                obj.embedded = true;
            end

            if obj.embedded
                obj.buildShapeAxesGrid(shapeParent);
                obj.shapeFig = ancestor(shapeParent, 'figure');
            else
                obj.shapeFig = figure('Visible', 'off', ...
                    'Name', 'GuliweRR -- Modal Shape Viewer (Plate)');
            end
            obj.buildFigure(controlsParent);
            obj.goToMode(sweepResult.angles(1), 1);
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
        %%% Pre-build the {3 x nLay} uiaxes grid Vis.plotModalShape renders
        %%% into via hAxCellIn -- geometry owned by GuiCommon.shapeAxesGrid,
        %%% see its own docblock for the manually-positioned-uiaxes rationale.
        %%% >> buildShapeAxesGrid(obj, shapeParent)
            relDepths = [obj.sweepResult.structure.calcParams.relDepth];
            [obj.shapeAxesCell, obj.rgbAxesCell, obj.keyAx] = GuiCommon.shapeAxesGrid(shapeParent, relDepths);
        end

        function buildFigure(obj, controlsParent)
            [obj.fig, mainLayout] = GuiCommon.hostLayout(controlsParent, ...
                'GuliweRR -- Modal Shape Viewer (Plate)', [100 100 560 360], [5 1]);
            mainLayout.RowHeight = {36, 36, 36, 36, 26};

            %%% row 1 -- angle controls
            angleStrs = arrayfun(@(a) sprintf('%.1f', a), obj.sweepResult.angles, 'UniformOutput', false);
            [obj.angleDropdown, obj.prevAngleBtn, obj.nextAngleBtn] = ...
                GuiCommon.dropdownStepperRow(mainLayout, 1, 'Angle:', angleStrs, {50, 140, 36, 36}, ...
                    @(~,~) obj.stepToAngle(obj.angleDropdown.Value), ...
                    @(~,~) obj.onPrevAngle(), @(~,~) obj.onNextAngle());

            %%% row 2 -- mode controls
            [obj.modeField, obj.prevModeBtn, obj.nextModeBtn, obj.velLabel] = ...
                GuiCommon.indexStepperRow(mainLayout, 2, 'Mode:', {50, 70, 36, 36, '1x'}, ...
                    @(~,~) obj.onModeFieldChanged(), ...
                    @(~,~) obj.onPrevMode(), @(~,~) obj.onNextMode());
            obj.velLabel.FontWeight = 'bold';   % the mode's velocity, as ModalShapeViewerGUI

            %%% row 3 -- current mode's type/confidence (replaces
            %%% ModalShapeViewerGUI's SAW/limiting list -- see class docblock)
            obj.modeTypeLabel = uilabel(mainLayout, 'Text', '');
            obj.modeTypeLabel.Layout.Row = 3;

            %%% row 4 -- display: signed components / polarisation RGB
            displayRow = uigridlayout(mainLayout, [1 2]);
            displayRow.Layout.Row = 4;
            displayRow.ColumnWidth = {60, 200};
            displayRow.Padding = [0 0 0 0];
            uilabel(displayRow, 'Text', 'Display:');
            obj.displayDropdown = GuiCommon.shapeDisplayDropdown(displayRow, @(~,~) obj.onDisplayChanged());

            obj.statusLine = GuiCommon.statusLabel(mainLayout, 5);
        end
    end

    %%% ===================================================================
    %%% Public entry point + callbacks
    %%% ===================================================================
    methods
        function goToMode(obj, angleVal, modeNo)
        %%% Jump directly to a specific (angle, Ritz mode).
        %%% >> goToMode(obj, angleVal, modeNo)
            [~, angNo] = min(abs(obj.sweepResult.angles - angleVal));
            nM = obj.sweepResult.cs{angNo}.ritzSolver.nCompModes;
            modeNo = min(max(round(modeNo), 1), nM);

            obj.currentAngNo  = angNo;
            obj.currentModeNo = modeNo;
            ritz = obj.sweepResult.cs{angNo}.ritzSolver;
            obj.currentVel = ritz.vels(modeNo);

            obj.syncControlsAndRedraw();
        end

        function goToSelection(obj, angleVal, modeNo, ~)
        %%% External entry point for AngSweepViewerGUIPlate.OnModeSelected --
        %%% same duck-typed 4-argument signature CombinedViewerGUI's wiring
        %%% uses ((angle, modeNo, velSnap, velRaw)), 4th argument unused since
        %%% Plate has no continuous-velocity source to fall back to.
        %%% >> goToSelection(obj, angleVal, modeNo, velRaw)
            obj.goToMode(angleVal, modeNo);
        end

        function stepToAngle(obj, newAngNo)
        %%% Keeps the same raw mode index across the angle step -- see class
        %%% docblock for why no type-following heuristic is applied here.
        %%% >> stepToAngle(obj, newAngNo)
            obj.goToMode(obj.sweepResult.angles(newAngNo), obj.currentModeNo);
        end

        function onPrevAngle(obj)
            obj.stepToAngle(max(obj.currentAngNo - 1, 1));
        end

        function onNextAngle(obj)
            obj.stepToAngle(min(obj.currentAngNo + 1, numel(obj.sweepResult.angles)));
        end

        function onModeFieldChanged(obj)
            obj.goToMode(obj.sweepResult.angles(obj.currentAngNo), obj.modeField.Value);
        end

        function onPrevMode(obj)
            obj.goToMode(obj.sweepResult.angles(obj.currentAngNo), obj.currentModeNo - 1);
        end

        function onNextMode(obj)
            obj.goToMode(obj.sweepResult.angles(obj.currentAngNo), obj.currentModeNo + 1);
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
        %%% >> syncControlsAndRedraw(obj)
            obj.angleDropdown.Value = obj.currentAngNo;
            obj.modeField.Value = obj.currentModeNo;
            ritz = obj.sweepResult.cs{obj.currentAngNo}.ritzSolver;
            obj.velLabel.Text = sprintf('vel: %.2f km/s', ritz.vels(obj.currentModeNo));
            obj.modeTypeLabel.Text = sprintf('Type: %s', ...
                GuiCommon.plateModeTypeText(ritz, obj.currentModeNo));

            obj.redrawShape();

            if ~isempty(obj.OnSelectionChanged)
                obj.OnSelectionChanged(obj.sweepResult.angles(obj.currentAngNo), obj.currentVel);
            end
        end

        function redrawShape(obj)
        %%% >> redrawShape(obj)
            angNo = obj.currentAngNo;
            dirs  = GuiCommon.shapeDirs(obj.displayDropdown);
            if ischar(dirs); axCell = obj.rgbAxesCell; else; axCell = obj.shapeAxesCell; end
            try
                ritz = obj.resolveDisplayRitz(angNo);
                if obj.embedded
                    Vis.plotModalShape(ritz, obj.currentModeNo, 'boosted', [], dirs, [], [], axCell);
                else
                    Vis.plotModalShape(ritz, obj.currentModeNo, 'boosted', [], dirs, [], obj.shapeFig);
                end
                obj.shapeFig.Visible = obj.fig.Visible;
                GuiCommon.setStatusOk(obj.statusLine, '');
            catch err
                GuiCommon.setStatusError(obj.statusLine, err.message);
            end
        end

        function ritz = resolveDisplayRitz(obj, angNo)
        %%% Return a RitzSolverLamb with a FULL eigVecs_con for angNo,
        %%% recomputing on demand when retentionLevel~='full' pruned it away
        %%% -- direct mirror of ModalShapeViewerGUI.resolveDisplayRitz,
        %%% RitzSolverLamb in place of RitzSolverSAW (see class docblock).
        %%% >> ritz = resolveDisplayRitz(obj, angNo)
            cs = obj.sweepResult.cs{angNo};
            if strcmp(cs.retentionLevel, 'full')
                ritz = cs.ritzSolver;
                return
            end
            if ~isempty(obj.fullRitzCache{angNo})
                ritz = obj.fullRitzCache{angNo};
                return
            end
            %%% angle 0 -- cs.structure is already rotated for this angle,
            %%% passing obj.angle here would double-rotate it (same
            %%% reasoning as ModalShapeViewerGUI's own resolveDisplayRitz)
            rs   = RitzSolverLamb(cs.structure, 0, cs.realThickness, cs.realWavVec);
            ritz = rs.calc(true, cs.ritzPrecompute());   % rebuilt if empty, e.g. after load()
            obj.fullRitzCache{angNo} = ritz;
        end
    end
end
