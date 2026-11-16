classdef ModalShapeViewerGUI < handle
%%% ModalShapeViewerGUI -- modal-shape / Green's-response viewer window
%%%
%%% windowMode, the 2nd constructor argument, controls how many physical windows
%%% this class opens:
%%%
%%%   ModalShapeViewerGUI(sweepResult)                 -- ONE window (default,
%%%       'combined'): self-builds a single uifigure with the SAME {'3x','1x'}
%%%       shapes-left/controls-right layout CombinedViewerGUI uses for this
%%%       class's top half, then renders into a pre-built {3 x nLay} uiaxes grid
%%%       (buildShapeAxesGrid) via plotModalShape/plotGreenModalShape's
%%%       hAxCellIn parameter.
%%%   ModalShapeViewerGUI(sweepResult, 'split')        -- TWO windows: a uifigure
%%%       control panel (this class's own .fig — angle/mode stepping, source
%%%       toggle, limiting-wave/SAW panel) plus a separate CLASSIC figure
%%%       (.shapeFig) that plotModalShape/plotGreenModalShape render into with
%%%       their classic axes() grid. Classic axes cannot be parented inside a
%%%       uifigure at all, which is why this mode needs its own figure.
%%%   ModalShapeViewerGUI(sweepResult, [], controlsParent, shapeParent)
%%%       -- EMBEDDED: controlsParent/shapeParent are two uigridlayout regions
%%%       of a caller-owned SHARED uifigure, and windowMode is ignored —
%%%       embedding is the same uiaxes-grid mechanism 'combined' uses, hosted by
%%%       the caller. Passing 'split' together with controlsParent is an error
%%%       (ModalShapeViewerGUI:constructor:splitRequiresStandalone): splitting only makes
%%%       sense when this class owns its own windows.
%%%
%%% obj.embedded (internal) is true whenever rendering goes through the
%%% uiaxes-grid path — both real embedding and standalone 'combined' — and false
%%% only for 'split' mode's classic-figure path.
%%%
%%% See docs/reference/architecture/gui_viewers.md for the full design rationale.
%%% Standalone constructible (like AngSweepViewerGUI): takes a sweepResult directly, not a parent
%%% controller object. Link two standalone windows manually via:
%%%   win1 = ModalShapeViewerGUI(sweep); win2 = AngSweepViewerGUI(sweep);
%%%   win2.OnModeSelected    = @(angle, modeNo, ~, rawVel) win1.goToSelection(angle, modeNo, rawVel);
%%%   win1.OnSelectionChanged = @(angle, vel) win2.showCrosshairAt(angle, vel);
%%% CombinedViewerGUI wires exactly this automatically.
%%%
%%% ANGLE STEPPING with source=Ritz modal shape follows the SAME PHYSICAL WAVE
%%% (SAW/limiting-type identity), not the same raw mode index: Ritz mode
%%% ordering by velocity can shuffle with angle, so "stay on this index" would
%%% jump to the wrong branch. With source=Green's response it instead keeps the
%%% same VELOCITY — there is no discrete mode or type to follow for an arbitrary
%%% continuous velocity, so "same velocity, different angle" is the analogue.
%%%
%%% CONTINUOUS-VELOCITY GREEN'S MODE: a map click on AngSweepViewerGUI with
%%% Green's response active has no discrete Ritz mode to snap to, so
%%% goToVelocity/currentVel decouple the displayed velocity from any Ritz mode
%%% index. The "Mode:" field is re-purposed as a direct velocity entry in this
%%% state — a relabelled meaning, not a new widget. currentModeNo is still kept
%%% pointed at the nearest Ritz mode so a switch back has a sensible starting
%%% point, but never drives the Green's-response calculation.
%%%
%%% RETENTION-AWARE RECOMPUTE: under retentionLevel~='full',
%%% cs{angNo}.ritzSolver.eigVecs_con may have been pruned to the significant
%%% modes or cleared entirely. resolveDisplayRitz re-solves the missing
%%% eigenvectors on demand, cached per angle, rather than freezing on stale data
%%% or requiring retentionLevel='full' to browse modal shapes at all.
%%%
%%% "Mark current mode as SAW/limiting" is a manual QC override of the automatic
%%% classification, SESSION-LOCAL ONLY: stored on this object's own overrides
%%% list, never written back into the result's classification, which is
%%% SetAccess=protected. See docs/reference/architecture/gui_viewers.md,
%%% "Session-local classification overrides".
%%%
%%% LIVE CLASSIFICATION: the classified modes are read through the result handles
%%% on every use (significantModesAt: cs{angNo}.classification.significantModes),
%%% never copied at construction, so a Reclassify elsewhere (the classification bar
%%% of the linked sweep / single-angle viewer) is seen here.
%%% refreshClassification rebuilds the SAW / limiting panel for it; hosts
%%% (CombinedViewerGUI, DoubleSweepViewerGUI) call it from the sweep viewer's
%%% OnClassificationChanged.
%%%
%%% GuliweRR v6.0

    properties
        fig
        shapeFig

        angleDropdown, prevAngleBtn, nextAngleBtn
        modeField, prevModeBtn, nextModeBtn, velLabel, modeLabel
        sourceDropdown, forceDirLabel, forceDirDropdown
        displayDropdown     % 'components' | 'rgb' -- GuiCommon.shapeDisplayDropdown
        limModePanel        % fixed outer panel/title -- 'Mark as...' buttons pinned at its
                              % bottom (built once), a scrollable sub-panel above them
        limModeScrollPanel   % the scrollable sub-panel refreshLimModePanel rebuilds per angle
        statusLine

        OnSelectionChanged = []   % function handle @(angle, vel), fired on every state change
    end

    properties (SetAccess = private)
        sweepResult      % CalcAngSweepHS or CalcAngSweepLayered, already calc()'d -- OR a
                          % plain struct(angles, cs, structure) adapter built by the
                          % constructor from a CalcSingleHS/CalcSingleLayered input (see
                          % constructor docblock); cs{1} is that result's own handle
        windowMode       % 'combined' (default, one window) | 'split' (two windows) --
                          % see class docblock. Always 'combined' when embedded.
        currentAngNo     % index into sweepResult.angles
        currentModeNo    % nearest/active Ritz mode index -- always kept in sync,
                          % even in Green's-response mode (see class docblock)
        currentVel       % active velocity -- ritz.vels(currentModeNo) in Ritz mode,
                          % an arbitrary continuous value in Green's-response mode
    end

    properties (Access = private)
        overrides = struct('angNo', {}, 'category', {}, 'modeNo', {})   % session-local "Mark as..." overrides
        embedded  = false    % true when rendering via the uiaxes grid (real embedding
                              % OR standalone windowMode='combined'); false only for
                              % standalone windowMode='split' (classic shapeFig)
        shapeAxesCell = {}   % {3 x nLay} pre-built uiaxes, whenever obj.embedded is true
        rgbAxesCell   = {}   % {1 x nLay} uiaxes of the polarisation-RGB display (embedded only)
        controlsLayout       % the controls' 5-row layout; row 4 (force direction) is
                              % shown only for Green's response (onSourceChanged)
        keyAx         = []   % its colour-key uiaxes
        VEL_STEP = 0.05      % km/s per << >> step in Green's-response mode
    end

    properties (Access = {?ModalShapeViewerGUI, ?test_ModalShapeViewerGUI})
        fullRitzCache = {}   % {1 x nAngles}, lazily filled -- see resolveDisplayRitz.
                              % Friend access to the test class (matches
                              % ModeClassification.m's own precedent) so the
                              % per-angle caching behaviour is directly
                              % testable without a flaky timing check.
    end

    properties (Constant, Access = private)
        SAW_TYPES          = {'Rayleigh', 'pseudoRayleigh', 'SAW'}   % 'SAW' is the layered
                                                                       % path's single type
                                                                       % (recognizeLayered)
    end

    methods
        function obj = ModalShapeViewerGUI(sweepResult, windowMode, controlsParent, shapeParent)
        %%% Construct the modal-shape/Green's-response viewer window
        %%% >> obj = ModalShapeViewerGUI(sweepResult, windowMode, controlsParent, shapeParent)
        %%%
        %%% INPUT: 1. sweepResult     .. CalcAngSweepHS or CalcAngSweepLayered,
        %%%                              already calc()'d -- OR CalcSingleHS/
        %%%                              CalcSingleLayered, already calcAll()'d
        %%%                              (packaged internally as a degenerate
        %%%                              1-angle "sweep", see below)
        %%%  (opt) 2. windowMode      .. 'combined' (default, one window) or
        %%%                              'split' (two windows) -- see class
        %%%                              docblock. Only meaningful for standalone
        %%%                              construction (no controlsParent/shapeParent);
        %%%                              ignored (must not be 'split') when embedding.
        %%%  (opt) 3. controlsParent  .. a uigridlayout region to build the
        %%%                              control panel into instead of a fresh
        %%%                              own uifigure
        %%%  (opt) 4. shapeParent     .. a uigridlayout region to build the
        %%%                              modal-shape uiaxes grid into. Required
        %%%                              together with controlsParent (both or
        %%%                              neither) -- used by CombinedViewerGUI.
        %%%                              Omitted (default): standalone, governed
        %%%                              by windowMode instead.
            validateattributes(sweepResult, ...
                {'CalcAngSweepHS', 'CalcAngSweepLayered', 'CalcSingleHS', 'CalcSingleLayered'}, ...
                {'scalar'}, 'ModalShapeViewerGUI', 'sweepResult');
            if isa(sweepResult, 'CalcSingleHS') || isa(sweepResult, 'CalcSingleLayered')
                %%% Single-angle result: package it as a degenerate 1-angle "sweep" --
                %%% CalcAngSweepHS/Layered's own cs{iAngle} IS a CalcSingleHS/
                %%% CalcSingleLayered instance, so a single result already has the
                %%% shape of one sweep element. This class reads only three fields off
                %%% sweepResult (.angles/.cs/.structure), all satisfied by a plain
                %%% struct; the classification is read through cs{1}, the result's
                %%% own handle (see class docblock, LIVE CLASSIFICATION). Angle
                %%% prev/next degrade to no-ops since onPrevAngle/onNextAngle clamp
                %%% against numel(angles)==1.
                if isempty(sweepResult.classification)
                    error('ModalShapeViewerGUI:constructor:notCalculated', ...
                        'sweepResult has not been calcAll()''d yet -- nothing to view.');
                end
                singleResult = sweepResult;
                sweepResult = struct( ...
                    'angles', singleResult.angle, ...
                    'cs', {{singleResult}}, ...
                    'structure', singleResult.structure);
            elseif isempty(sweepResult.significantModesMap)
                error('ModalShapeViewerGUI:constructor:notCalculated', ...
                    'sweepResult has not been calc()''d yet -- nothing to view.');
            end
            if nargin < 2 || isempty(windowMode); windowMode = 'combined'; end
            if nargin < 3; controlsParent = []; end
            if nargin < 4; shapeParent = []; end
            GuiCommon.validateWindowMode('ModalShapeViewerGUI', windowMode, {'combined', 'split'});
            GuiCommon.validateEmbedParents('ModalShapeViewerGUI', windowMode, controlsParent, shapeParent);
            obj.sweepResult   = sweepResult;
            obj.windowMode    = windowMode;
            obj.embedded      = ~isempty(controlsParent);
            obj.fullRitzCache = cell(1, numel(sweepResult.angles));

            if isempty(controlsParent) && strcmp(windowMode, 'combined')
                %%% standalone 'combined': self-build the SAME top-half layout
                %%% CombinedViewerGUI itself builds for this class (shapes left/
                %%% controls right, {'3x','1x'} column split) inside one owned
                %%% uifigure -- one window instead of two, via the identical
                %%% uiaxes-grid rendering path real embedding already uses.
                %%% Wider than CombinedViewerGUI's own split (1400 vs. ~1100 total)
                %%% since this window doesn't share vertical room with a bottom
                %%% angular-sweep panel, and the control rows below (e.g. the
                %%% ~440px-wide source-toggle row) need a reasonably wide 1x column.
                ownFig = uifigure('Name', 'GuliweRR -- Modal Shape Viewer', ...
                    'Position', [80 60 1400 640], 'Visible', 'off');
                [shapeParent, controlsParent] = GuiCommon.shapeControlsSplit(ownFig, [4 4 4 4]);
                obj.embedded = true;
            end

            if obj.embedded
                obj.buildShapeAxesGrid(shapeParent);
                obj.shapeFig = ancestor(shapeParent, 'figure');
            else
                obj.shapeFig = figure('Visible', 'off', ...
                    'Name', 'GuliweRR -- Modal Shape / Green''s Response');
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
        %%% Pre-build the {3 x nLay} uiaxes grid plotModalShape/
        %%% plotGreenModalShape render into via hAxCellIn, whenever obj.embedded
        %%% is true (real embedding, or standalone windowMode='combined').
        %%% GuiCommon.shapeAxesGrid owns the geometry (row weights, margins,
        %%% the manually-positioned-uiaxes-in-a-uipanel mechanism) -- see its
        %%% own docblock for the full reasoning behind every constant there.
        %%%
        %%% Green's response fills the same grid: Vis.plotGreenModalShapeLayered
        %%% for a layered structure (rendered 'boosted', like the Ritz shape, so
        %%% the drawn rows match these weights), Vis.plotGreenModalShape into its
        %%% single row for a half-space.
        %%% >> buildShapeAxesGrid(obj, shapeParent)
            relDepths = [obj.sweepResult.structure.calcParams.relDepth];
            [obj.shapeAxesCell, obj.rgbAxesCell, obj.keyAx] = GuiCommon.shapeAxesGrid(shapeParent, relDepths);
        end

        function buildFigure(obj, controlsParent)
            [obj.fig, mainLayout] = GuiCommon.hostLayout(controlsParent, ...
                'GuliweRR -- Modal Shape Viewer', [100 100 620 520], [5 1]);
            %%% Angle/Mode/Source rows get just what a label+field+button row
            %%% needs, leaving the rest to row 5's SAW/limiting-mode panel. Row 4
            %%% (force direction) has height 0 until the source is Green's response
            %%% (onSourceChanged). No trailing status row: statusLine lives inside
            %%% that panel (above the mark buttons), since the shape-plot side has no
            %%% matching row and a dedicated one here would leave the panel's bottom
            %%% edge above the plot's.
            mainLayout.RowHeight = {30, 30, 30, 0, '1x'};
            obj.controlsLayout = mainLayout;

            %%% row 1 -- angle controls
            angleStrs = arrayfun(@(a) sprintf('%.1f', a), obj.sweepResult.angles, 'UniformOutput', false);
            [obj.angleDropdown, obj.prevAngleBtn, obj.nextAngleBtn] = ...
                GuiCommon.dropdownStepperRow(mainLayout, 1, 'Angle:', angleStrs, {50, 140, 36, 36}, ...
                    @(~,~) obj.stepToAngle(obj.angleDropdown.Value), ...
                    @(~,~) obj.onPrevAngle(), @(~,~) obj.onNextAngle());

            %%% row 2 -- mode controls (re-purposed as a velocity field in
            %%% Green's-response mode -- see onSourceChanged, which relaxes the
            %%% shared row's Limits/RoundFractionalValues at run time)
            [obj.modeField, obj.prevModeBtn, obj.nextModeBtn, obj.velLabel, obj.modeLabel] = ...
                GuiCommon.indexStepperRow(mainLayout, 2, 'Mode:', {50, 70, 36, 36, '1x'}, ...
                    @(~,~) obj.onModeFieldChanged(), ...
                    @(~,~) obj.onPrevMode(), @(~,~) obj.onNextMode());
            %%% the selected mode's velocity, short and bold so it reads at a glance
            obj.velLabel.FontWeight = 'bold';

            %%% row 3 -- source toggle (Ritz modal shape / Green's response) and,
            %%% beside it without a label of its own, the display (signed components /
            %%% polarisation RGB): one row, to leave the height to the mode list
            sourceRow = uigridlayout(mainLayout, [1 3]);
            sourceRow.Layout.Row = 3;
            sourceRow.ColumnWidth = {50, '1x', '1x'};
            sourceRow.Padding = [0 0 0 0];
            uilabel(sourceRow, 'Text', 'Source:');
            obj.sourceDropdown = uidropdown(sourceRow, ...
                'Items', {'Ritz modal shape', 'Green''s response'}, ...
                'ValueChangedFcn', @(~,~) obj.onSourceChanged());
            obj.displayDropdown = GuiCommon.shapeDisplayDropdown(sourceRow, @(~,~) obj.onDisplayChanged());

            %%% row 4 -- force direction, Green's response only (height 0 otherwise)
            forceRow = uigridlayout(mainLayout, [1 2]);
            forceRow.Layout.Row = 4;
            forceRow.ColumnWidth = {70, 140};
            forceRow.Padding = [0 0 0 0];
            obj.forceDirLabel = uilabel(forceRow, 'Text', 'Force dir:', 'Visible', 'off');
            obj.forceDirDropdown = uidropdown(forceRow, ...
                'Items', {'x (propag. dir.)', 'y (in-plane)', 'z (out-of-plane)'}, ...
                'ItemsData', 1:3, 'Value', 1, 'Visible', 'off', ...
                'ValueChangedFcn', @(~,~) obj.onForceDirChanged());

            %%% row 5 -- SAW / limiting-mode panel. The "Mark as..." buttons are
            %%% pinned in a fixed row at the bottom, built once here -- the
            %%% variable-length mode list (rebuilt per angle by
            %%% refreshLimModePanel) lives in a NESTED scrollable sub-panel above
            %%% them, so scrolling the list can never push the buttons out of view.
            %%% The list's row count is data-driven (however many SAW/limiting
            %%% modes this angle has) and can always exceed the space available,
            %%% hence a scrollbar rather than trying to guarantee a fit.
            obj.limModePanel = uipanel(mainLayout, 'Title', 'SAW / limiting modes');
            obj.limModePanel.Layout.Row = 5;
            limModeGl = uigridlayout(obj.limModePanel, [3 1]);
            %%% row 1: scrollable list; row 2: the status line (see RowHeight comment
            %%% above); row 3, the panel's bottom line: the inline "Mark current mode
            %%% as: [SAW] [lim]" row
            limModeGl.RowHeight = {'1x', 16, 30};
            limModeGl.Padding = [0 0 0 0];

            %%% AutoResizeChildren='off' is required alongside Scrollable='on' here --
            %%% a uipanel's default AutoResizeChildren='on' forcibly rescales its
            %%% child grid to the panel's own box on every layout pass, which
            %%% silently defeats Scrollable's "let content exceed the box, add a
            %%% scrollbar" behaviour (uigridlayout has no AutoResizeChildren
            %%% property, so this doesn't apply to StructureControlGUI's Scrollable
            %%% grid). The panel's Scrollable alone is not enough: refreshLimModePanel
            %%% also sets Scrollable on the grid it rebuilds inside, so the grid
            %%% itself opts out of being force-fit too.
            obj.limModeScrollPanel = uipanel(limModeGl, 'BorderType', 'none', ...
                'Scrollable', 'on', 'AutoResizeChildren', 'off');
            obj.limModeScrollPanel.Layout.Row = 1;

            %%% one inline row -- label + two compact buttons, the label right-aligned
            %%% against the buttons -- to leave vertical space for the list above
            markGl = uigridlayout(limModeGl, [1 3]);
            markGl.Layout.Row = 3;
            markGl.ColumnWidth = {'1x', 50, 46};
            markGl.Padding = [0 0 0 0];
            markLbl = uilabel(markGl, 'Text', 'Mark current mode as:', 'HorizontalAlignment', 'right');
            markLbl.Layout.Column = 1;
            markSaw = uibutton(markGl, 'Text', 'SAW', 'Tooltip', 'Mark current mode as SAW', ...
                'ButtonPushedFcn', @(~,~) obj.onMarkAsShown('SAW'));
            markSaw.Layout.Column = 2;
            markLim = uibutton(markGl, 'Text', 'lim', 'Tooltip', 'Mark current mode as limiting', ...
                'ButtonPushedFcn', @(~,~) obj.onMarkAsShown('limiting'));
            markLim.Layout.Column = 3;

            obj.statusLine = GuiCommon.statusLabel(limModeGl, 2);
        end
    end

    %%% ===================================================================
    %%% Public entry point + callbacks -- public (not Access=private):
    %%% StructureControlGUI's own convention, so tests can call these
    %%% directly, and so AngSweepViewerGUI.OnModeSelected can be wired
    %%% straight to goToSelection.
    %%% ===================================================================
    methods
        function goToMode(obj, angleVal, modeNo)
        %%% Jump directly to a specific (angle, Ritz mode) -- switches source
        %%% to 'Ritz modal shape' if it wasn't already.
        %%% >> goToMode(obj, angleVal, modeNo)
            [~, angNo] = min(abs(obj.sweepResult.angles - angleVal));
            nM = obj.sweepResult.cs{angNo}.ritzSolver.nCompModes;
            modeNo = min(max(round(modeNo), 1), nM);

            obj.currentAngNo  = angNo;
            obj.currentModeNo = modeNo;
            ritz = obj.sweepResult.cs{angNo}.ritzSolver;
            obj.currentVel = ritz.vels(modeNo);   % kept in sync so a later source switch is continuous

            obj.sourceDropdown.Value = 'Ritz modal shape';
            obj.syncControlsAndRedraw();
        end

        function goToVelocity(obj, angleVal, velVal)
        %%% Jump to a specific (angle, velocity) in Green's-response mode --
        %%% the continuous-velocity counterpart to goToMode's discrete Ritz
        %%% mode index. Switches source to 'Green's response' if it wasn't
        %%% already; there is no "arbitrary velocity" concept for a discrete
        %%% Ritz mode.
        %%% >> goToVelocity(obj, angleVal, velVal)
            [~, angNo] = min(abs(obj.sweepResult.angles - angleVal));
            obj.currentAngNo = angNo;
            obj.currentVel   = velVal;

            %%% keep currentModeNo pointed at the nearest Ritz mode purely so
            %%% a later switch BACK to Ritz modal shape starts somewhere sensible
            ritz = obj.sweepResult.cs{angNo}.ritzSolver;
            nM   = ritz.nCompModes;
            [~, obj.currentModeNo] = min(abs(ritz.vels(1:nM) - velVal));

            obj.sourceDropdown.Value = 'Green''s response';
            obj.syncControlsAndRedraw();
        end

        function goToSelection(obj, angleVal, modeNo, rawVel)
        %%% External entry point for AngSweepViewerGUI.OnModeSelected -- picks
        %%% goToMode or goToVelocity based on THIS window's own current
        %%% source, since only this window knows which one currently makes
        %%% sense (AngSweepViewerGUI deliberately stays unaware of it).
        %%% >> goToSelection(obj, angleVal, modeNo, rawVel)
            if strcmp(obj.sourceDropdown.Value, 'Ritz modal shape')
                obj.goToMode(angleVal, modeNo);
            else
                obj.goToVelocity(angleVal, rawVel);
            end
        end

        function stepToAngle(obj, newAngNo)
        %%% >> stepToAngle(obj, newAngNo)
            if strcmp(obj.sourceDropdown.Value, 'Ritz modal shape')
                newModeNo = obj.followPhysicalWave(obj.currentAngNo, obj.currentModeNo, newAngNo);
                obj.goToMode(obj.sweepResult.angles(newAngNo), newModeNo);
            else
                %%% no discrete mode/type to "follow" for an arbitrary
                %%% velocity -- same velocity, different angle instead
                obj.goToVelocity(obj.sweepResult.angles(newAngNo), obj.currentVel);
            end
        end

        function onPrevAngle(obj)
            obj.stepToAngle(max(obj.currentAngNo - 1, 1));
        end

        function onNextAngle(obj)
            obj.stepToAngle(min(obj.currentAngNo + 1, numel(obj.sweepResult.angles)));
        end

        function onModeFieldChanged(obj)
        %%% modeField's MEANING depends on the active source -- see
        %%% onSourceChanged/syncControlsAndRedraw.
        %%% >> onModeFieldChanged(obj)
            if strcmp(obj.sourceDropdown.Value, 'Ritz modal shape')
                obj.goToMode(obj.sweepResult.angles(obj.currentAngNo), obj.modeField.Value);
            else
                obj.goToVelocity(obj.sweepResult.angles(obj.currentAngNo), obj.modeField.Value);
            end
        end

        function onPrevMode(obj)
            if strcmp(obj.sourceDropdown.Value, 'Ritz modal shape')
                obj.goToMode(obj.sweepResult.angles(obj.currentAngNo), obj.currentModeNo - 1);
            else
                obj.goToVelocity(obj.sweepResult.angles(obj.currentAngNo), obj.currentVel - obj.VEL_STEP);
            end
        end

        function onNextMode(obj)
            if strcmp(obj.sourceDropdown.Value, 'Ritz modal shape')
                obj.goToMode(obj.sweepResult.angles(obj.currentAngNo), obj.currentModeNo + 1);
            else
                obj.goToVelocity(obj.sweepResult.angles(obj.currentAngNo), obj.currentVel + obj.VEL_STEP);
            end
        end

        function onSourceChanged(obj)
        %%% Source dropdown callback; the source-dependent widgets are set up
        %%% in syncControlsAndRedraw (applySourceWidgets)
        %%% >> onSourceChanged(obj)
            obj.syncControlsAndRedraw();
        end

        function onForceDirChanged(obj)
            obj.redrawShape();
        end

        function onDisplayChanged(obj)
        %%% Switch between the signed component panels and the polarisation-RGB
        %%% panel ('split' mode redraws its classic figure from scratch anyway)
        %%% >> onDisplayChanged(obj)
            if obj.embedded
                GuiCommon.setShapeDisplayMode(obj.shapeAxesCell, obj.rgbAxesCell, obj.keyAx, ...
                    strcmp(obj.displayDropdown.Value, 'rgb'));
            end
            obj.redrawShape();
        end

        function onMarkAsShown(obj, category)
        %%% Register the CURRENTLY DISPLAYED mode as the shown SAW/limiting
        %%% mode at the current angle (session-local, see class docblock).
        %%% >> onMarkAsShown(obj, category)
            keep = true(1, numel(obj.overrides));
            for k = 1:numel(obj.overrides)
                if obj.overrides(k).angNo == obj.currentAngNo && strcmp(obj.overrides(k).category, category)
                    keep(k) = false;
                end
            end
            obj.overrides = obj.overrides(keep);
            obj.overrides(end+1) = struct('angNo', obj.currentAngNo, ...
                'category', category, 'modeNo', obj.currentModeNo);

            GuiCommon.setStatusOk(obj.statusLine, sprintf(...
                'Marked mode %d as the shown %s at angle %.1f deg (this session only).', ...
                obj.currentModeNo, category, obj.sweepResult.angles(obj.currentAngNo)));
            obj.refreshLimModePanel();
        end

        function refreshClassification(obj)
        %%% Rebuild what depends on the classification after it changed (a
        %%% Reclassify in a linked viewer): the SAW / limiting panel of the current
        %%% angle. The shown mode and shape stay as they are.
        %%% >> refreshClassification(obj)
            obj.refreshLimModePanel();
        end

        function category = classifyCurrentMode(obj, angNo, modeNo)
        %%% 'SAW' | 'limiting' | 'none' -- session-local overrides take
        %%% priority over the automatic classification. Public (not
        %%% Access=private): this is the core "follow the same physical
        %%% wave" logic, precise enough to warrant its own direct tests
        %%% rather than only verifying it through emergent side effects.
        %%% >> category = classifyCurrentMode(obj, angNo, modeNo)
            for k = 1:numel(obj.overrides)
                if obj.overrides(k).angNo == angNo && obj.overrides(k).modeNo == modeNo
                    category = obj.overrides(k).category;
                    return
                end
            end
            sm  = obj.significantModesAt(angNo);
            idx = find([sm.modeNo] == modeNo, 1);
            if isempty(idx)
                category = 'none';
            elseif ismember(sm(idx).type, obj.SAW_TYPES)
                category = 'SAW';
            elseif strcmp(sm(idx).type, 'limiting')
                category = 'limiting';
            else
                category = 'none';
            end
        end

        function newModeNo = followPhysicalWave(obj, oldAngNo, oldModeNo, newAngNo)
        %%% Try to keep showing the SAME PHYSICAL WAVE across an angle step --
        %%% see the class docblock. Public for the same reason as
        %%% classifyCurrentMode above.
        %%% >> newModeNo = followPhysicalWave(obj, oldAngNo, oldModeNo, newAngNo)
            category  = obj.classifyCurrentMode(oldAngNo, oldModeNo);
            newModeNo = oldModeNo;   % default: keep the same raw index

            smNew = obj.significantModesAt(newAngNo);
            switch category
                case 'SAW'
                    override = obj.findOverride(newAngNo, 'SAW');
                    isSAW = ismember({smNew.type}, obj.SAW_TYPES);
                    if ~isempty(override)
                        newModeNo = override;
                    elseif any(isSAW)
                        entries = smNew(isSAW);
                        newModeNo = entries(1).modeNo;
                    end
                case 'limiting'
                    override = obj.findOverride(newAngNo, 'limiting');
                    isLim = strcmp({smNew.type}, 'limiting');
                    if ~isempty(override)
                        newModeNo = override;
                    elseif any(isLim)
                        oldVel  = obj.sweepResult.cs{oldAngNo}.ritzSolver.vels(oldModeNo);
                        entries = smNew(isLim);
                        [~, closest] = min(abs([entries.vel] - oldVel));   % a limiting entry's vel is its mode's eigenvelocity
                        newModeNo = entries(closest).modeNo;
                    end
            end

            nM = obj.sweepResult.cs{newAngNo}.ritzSolver.nCompModes;
            newModeNo = min(max(newModeNo, 1), nM);
        end
    end

    %%% ===================================================================
    %%% Run-time helpers
    %%% ===================================================================
    methods (Access = private)
        function applySourceWidgets(obj, isGreen)
        %%% Set the widgets whose meaning depends on the source: the force-
        %%% direction row and the mode field (integer Ritz mode vs. continuous
        %%% velocity). Called from syncControlsAndRedraw, NOT only from the
        %%% dropdown callback: goToMode/goToVelocity switch the source from
        %%% code, and a programmatic .Value change fires no ValueChangedFcn --
        %%% the mode field then kept rounding a velocity to an integer.
        %%% >> applySourceWidgets(obj, isGreen)
            obj.forceDirLabel.Visible    = isGreen;
            obj.forceDirDropdown.Visible = isGreen;
            obj.controlsLayout.RowHeight{4} = 30 * isGreen;   % force-direction row
            if isGreen
                obj.modeField.Limits = [0 Inf];
                obj.modeField.RoundFractionalValues = 'off';
                obj.modeLabel.Text = 'Velocity:';
            else
                obj.modeField.Limits = [1 Inf];
                obj.modeField.RoundFractionalValues = 'on';
                obj.modeLabel.Text = 'Mode:';
            end
        end

        function syncControlsAndRedraw(obj)
        %%% Shared tail every state-changing method (goToMode/goToVelocity)
        %%% funnels through: sync every widget to (currentAngNo, currentModeNo/
        %%% currentVel, source), redraw the shape, rebuild the limiting-mode
        %%% panel, and notify OnSelectionChanged.
        %%% >> syncControlsAndRedraw(obj)
            obj.angleDropdown.Value = obj.currentAngNo;
            isGreen = strcmp(obj.sourceDropdown.Value, 'Green''s response');
            obj.applySourceWidgets(isGreen);
            if isGreen
                obj.modeField.Value = obj.currentVel;
                obj.velLabel.Text   = '';   % redundant with the field itself now
            else
                obj.modeField.Value = obj.currentModeNo;
                ritz = obj.sweepResult.cs{obj.currentAngNo}.ritzSolver;
                obj.velLabel.Text = sprintf('vel: %.2f km/s', ritz.vels(obj.currentModeNo));
            end

            obj.redrawShape();
            obj.refreshLimModePanel();

            if ~isempty(obj.OnSelectionChanged)
                obj.OnSelectionChanged(obj.sweepResult.angles(obj.currentAngNo), obj.currentVel);
            end
        end

        function redrawShape(obj)
        %%% Re-render the shape display for the current angle/velocity/
        %%% source/forceDir selection -- into shapeAxesCell (obj.embedded) or
        %%% the classic shapeFig (windowMode='split').
        %%% >> redrawShape(obj)
            angNo = obj.currentAngNo;
            dirs  = GuiCommon.shapeDirs(obj.displayDropdown);
            if ischar(dirs)
                axCell = obj.rgbAxesCell;     % {1 x nLay}
            else
                axCell = obj.shapeAxesCell;   % {3 x nLay}
            end
            try
                if strcmp(obj.sourceDropdown.Value, 'Ritz modal shape')
                    ritz = obj.resolveDisplayRitz(angNo);
                    %%% 'boosted' unconditionally -- see buildShapeAxesGrid; safe
                    %%% for nLay==1 too, where calcModalShapeRowHeights returns 1
                    if obj.embedded
                        Vis.plotModalShape(ritz, obj.currentModeNo, 'boosted', [], dirs, [], [], axCell);
                    else
                        Vis.plotModalShape(ritz, obj.currentModeNo, 'boosted', [], dirs, [], obj.shapeFig);
                    end
                else
                    forceDir = obj.forceDirDropdown.Value;
                    cs       = obj.sweepResult.cs{angNo};
                    %%% cs.structure with angle 0, NOT sweepResult.structure with
                    %%% angleVal. cs.structure is always already rotated for this
                    %%% angle -- the same object resolveDisplayRitz uses above --
                    %%% whereas sweepResult.structure is unrotated only for a real
                    %%% multi-angle sweep. For a single-angle CalcSingleHS/
                    %%% CalcSingleLayered input, packaged as a degenerate 1-angle
                    %%% sweep by the constructor, sweepResult.structure IS
                    %%% cs.structure, so passing angleVal would double-rotate it.
                    [greenObj, ~] = GreenSolver.calc(cs.structure, 0, obj.currentVel);
                    if cs.structure.nLay > 1
                        %%% one row per layer, 'boosted' as the grid (see buildShapeAxesGrid)
                        if ~obj.embedded; axCell = []; end
                        Vis.plotGreenModalShapeLayered(greenObj, cs.structure, forceDir, 'boosted', [], ...
                            dirs, [], obj.shapeFig, axCell);
                    elseif obj.embedded
                        Vis.plotGreenModalShape(greenObj, forceDir, [], dirs, [], [], axCell(:, 1));
                    else
                        Vis.plotGreenModalShape(greenObj, forceDir, [], dirs, [], obj.shapeFig);
                    end
                end
                obj.shapeFig.Visible = obj.fig.Visible;   % follow the control window's own visibility
                GuiCommon.setStatusOk(obj.statusLine, '');
            catch err
                GuiCommon.setStatusError(obj.statusLine, err.message);
            end
        end

        function ritz = resolveDisplayRitz(obj, angNo)
        %%% Return a RitzSolver with a FULL eigVecs_con for angNo, recomputing on
        %%% demand when cs{angNo}.retentionLevel pruned it away. plotModalShape
        %%% deliberately errors rather than guessing; recomputing here, one layer
        %%% up, keeps that contract while making the GUI transparent to
        %%% retentionLevel.
        %%%
        %%% ALWAYS recompute under a non-'full' retentionLevel, never merely when
        %%% eigVecs_con looks "too small" for modeNo: applyRetention('saw_only')
        %%% calls pruneEigVecs(keepCols), which COMPACTS eigVecs_con down to
        %%% numel(keepCols) columns rather than zeroing the dropped ones in place,
        %%% so column j of the pruned matrix is original mode keepCols(j), not
        %%% mode j. A naive "nColsAvail >= modeNo" check would silently plot the
        %%% WRONG physical mode whenever modeNo <= numel(keepCols).
        %%%
        %%% Cached per angle, session-local: the recompute is a real Ritz solve,
        %%% so repeated mode-stepping at one angle must not re-solve every redraw.
        %%% NOT written back into obj.sweepResult, which would defeat
        %%% retentionLevel's whole purpose.
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
            %%% same call sequence CalcSingleHS/Layered.calcRitz() itself uses
            %%% (angle 0 -- cs.structure is already rotated for this angle,
            %%% passing obj.angle here would double-rotate it)
            rs   = RitzSolverSAW(cs.structure, 0);
            ritz = rs.calc(true, cs.ritzPrecompute());   % rebuilt if empty, e.g. after load()
            obj.fullRitzCache{angNo} = ritz;
        end

        function sm = significantModesAt(obj, angNo)
        %%% The classified modes of one angle, read through the result handle (see
        %%% class docblock, LIVE CLASSIFICATION)
        %%% >> sm = significantModesAt(obj, angNo)
            sm = obj.sweepResult.cs{angNo}.classification.significantModes;
        end

        function modeNo = findOverride(obj, angNo, category)
        %%% >> modeNo = findOverride(obj, angNo, category)  -- [] if none
            modeNo = [];
            for k = 1:numel(obj.overrides)
                if obj.overrides(k).angNo == angNo && strcmp(obj.overrides(k).category, category)
                    modeNo = obj.overrides(k).modeNo;
                    return
                end
            end
        end

        function refreshLimModePanel(obj)
        %%% Rebuild the SAW/limiting-mode list for the current angle --
        %%% variable row count, so the whole gridlayout is reconstructed
        %%% rather than resized (as AngSweepViewerGUI does for variable-shape
        %%% content). Only limModeScrollPanel is rebuilt; the "Mark as..."
        %%% buttons are built once in buildFigure, outside the scrollable region.
        %%% >> refreshLimModePanel(obj)
            delete(obj.limModeScrollPanel.Children);

            sm = obj.significantModesAt(obj.currentAngNo);
            isSAW = ismember({sm.type}, obj.SAW_TYPES);
            isLim = strcmp({sm.type}, 'limiting');
            sawEntries = sm(isSAW);
            limEntries = sm(isLim);

            nRows = 2 + numel(sawEntries) + numel(limEntries);
            gl = uigridlayout(obj.limModeScrollPanel, [nRows, 2]);
            gl.ColumnWidth = {'1x', 70};
            %%% 'fit' rows, so a long mode-type/velocity label grows its own row
            %%% instead of overflowing (as in StructureControlGUI's calc-fields box)
            gl.RowHeight = repmat({'fit'}, 1, nRows);
            %%% Scrollable on the grid as well as the panel -- the panel's alone
            %%% leaves this grid force-fit and clipped (see buildFigure)
            gl.Scrollable = 'on';

            r = 1;
            hdr1 = uilabel(gl, 'Text', 'SAW-like modes:', 'FontWeight', 'bold');
            hdr1.Layout.Row = r; hdr1.Layout.Column = [1 2];
            r = r + 1;
            for k = 1:numel(sawEntries)
                e = sawEntries(k);
                txt = uilabel(gl, 'Text', sprintf('%s: mode %d, %.4f km/s, %s %.1f dB (%s)', ...
                    e.type, e.modeNo, e.vel, e.method, e.scoreDb, e.grade));
                txt.Layout.Row = r; txt.Layout.Column = 1;
                btn = uibutton(gl, 'Text', 'Go To', ...
                    'ButtonPushedFcn', @(~,~) obj.goToMode(obj.sweepResult.angles(obj.currentAngNo), e.modeNo));
                btn.Layout.Row = r; btn.Layout.Column = 2;
                r = r + 1;
            end

            hdr2 = uilabel(gl, 'Text', 'Limiting-velocity modes:', 'FontWeight', 'bold');
            hdr2.Layout.Row = r; hdr2.Layout.Column = [1 2];
            r = r + 1;
            for k = 1:numel(limEntries)
                e = limEntries(k);
                txt = uilabel(gl, 'Text', sprintf('mode %d, %.4f km/s', e.modeNo, e.vel));
                txt.Layout.Row = r; txt.Layout.Column = 1;
                btn = uibutton(gl, 'Text', 'Go To', ...
                    'ButtonPushedFcn', @(~,~) obj.goToMode(obj.sweepResult.angles(obj.currentAngNo), e.modeNo));
                btn.Layout.Row = r; btn.Layout.Column = 2;
                r = r + 1;
            end
        end
    end
end
