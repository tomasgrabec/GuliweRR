classdef DoubleSweepViewerGUI < handle
%%% DoubleSweepViewerGUI -- interactive viewer for a genuine DOUBLE sweep
%%% (CalcWavevecAngSweep2Layered with BOTH numel(angles)>1 AND
%%% numel(wavelengths)>1), switching between an angular-map view (at a
%%% chosen fixed wavelength) and a dispersion-curve view (at a chosen fixed
%%% angle), remembering axis-range settings on each side of the switch.
%%%
%%% A separate class rather than a CombinedViewerGUI mode. CombinedViewerGUI
%%% already has a "map <-> single-angle curve" swap, but that lives entirely
%%% inside AngSweepViewerGUI and works bottom-only, because both sides pair
%%% with the SAME top (ModalShapeViewerGUI, angle-oriented). Angular-map and
%%% dispersion-curve mode do not: DispersionViewerGUI exists as its own class
%%% precisely because ModalShapeViewerGUI's angle-oriented API does not fit
%%% wavelength-stepping, so it renders modal shapes itself. This class
%%% therefore swaps a whole TOP+BOTTOM PAIR at once
%%% (ModalShapeViewerGUI+AngSweepViewerGUI <-> DispersionViewerGUI's
%%% self-contained shape+curve content), by the same mechanics
%%% AngSweepViewerGUI's own swap uses — stable content slot, external toggle
%%% button, saved-and-restored axis ranges — one level up.
%%%
%%% ANGULAR MODE: ModalShapeViewerGUI + AngSweepViewerGUI, wired as
%%% CombinedViewerGUI wires them, pointed at sweepResult.sweeps{currentWavNo} —
%%% already a complete, independently-calc()'d CalcAngSweepLayered, so
%%% switching wavelength needs no recomputation, just rebuilding both viewers
%%% against a different element of sweepResult.sweeps.
%%%
%%% DISPERSION MODE: DispersionViewerGUI embedded via its parentContainer
%%% support, pointed at sweepResult with angleIdx=currentAngNo — again no
%%% recomputation, just re-embedding at a different angle slice.
%%%
%%% RECLASSIFY: both modes carry a classification bar (AngSweepViewerGUI's,
%%% DispersionViewerGUI's), and both classify the WHOLE double sweep
%%% (CalcWavevecAngSweep2Layered.classify: every wavelength, 2D maps reassembled),
%%% never one wavelength's slice alone: the angular sweep viewer is built with
%%% classifyTarget = sweepResult. After a Reclassify in angular mode it redraws and
%%% the modal-shape viewer refreshes (the CombinedViewerGUI chain); the dispersion
%%% viewer redraws its own branches. A mode switch rebuilds the other side from the
%%% live result, so it shows the new classification too.
%%%
%%% currentWavNo/currentAngNo are tracked independently and persist across mode
%%% switches: entering angular mode uses whatever wavelength was last active
%%% (from a previous visit or the dispersion view's wavelength stepper), and
%%% entering dispersion mode uses whatever angle was last active (a previous
%%% visit, a map click, or ModalShapeViewerGUI's angle-stepping).
%%%
%%% RAM WARNING: dispersion mode needs sweepResult calc()'d with
%%% retentionLevel='full' (DispersionBranches.track's own requirement). This is
%%% NOT validated at construction, since angular mode alone does not need it
%%% and requiring it unconditionally would block viewing a large double sweep
%%% computed cheaply for its angular side. Switching to dispersion mode without
%%% it surfaces DispersionBranches.track's missingEigVecs error instead. Full
%%% retention across BOTH swept axes of a real-sized double sweep is a genuinely
%%% large RAM cost — see docs/pitfalls.md NUM4.
%%%
%%% Single window only, unlike CombinedViewerGUI's combined/split/fullSplit
%%% toggle.
%%%
%%% GuliweRR v6.0

    properties
        fig
        switchModeBtn
        statusLine
    end

    properties (SetAccess = private)
        sweepResult    % CalcWavevecAngSweep2Layered, calc()'d, numel(angles)>1 AND numel(wavelengths)>1
        mode           % 'angular' | 'dispersion'
        currentWavNo   % index into sweepResult.wavelengths -- used by angular mode
        currentAngNo   % index into sweepResult.angles -- used by dispersion mode
        angularSweepViewer = []   % AngSweepViewerGUI, while mode=='angular' -- readable like
                                  % CombinedViewerGUI.angSweepViewer, for scripted set-up
        dispersionViewer   = []   % DispersionViewerGUI, while mode=='dispersion' (same)
    end

    properties (Access = {?DoubleSweepViewerGUI, ?test_DoubleSweepViewerGUI})
        contentParent
        angularModalShapeViewer = []   % ModalShapeViewerGUI, while mode=='angular'
        savedAxisRangeAngular    = []
        savedAxisRangeDispersion = []
        %%% Friend access to the test class -- matches ModalShapeViewerGUI.
        %%% fullRitzCache's own precedent -- so tests can confirm WHICH
        %%% wavelength/angle slice a mode actually rebuilt against, and
        %%% inspect saved axis-range state directly, without those needing
        %%% to be part of this class' own public surface.
    end

    methods
        function obj = DoubleSweepViewerGUI(sweepResult)
        %%% Construct the double-sweep viewer
        %%% >> obj = DoubleSweepViewerGUI(sweepResult)
        %%%
        %%% INPUT: 1. sweepResult .. CalcWavevecAngSweep2Layered, already
        %%%                          calc()'d, with numel(angles)>1 AND
        %%%                          numel(wavelengths)>1 (a genuine double
        %%%                          sweep -- use CombinedViewerGUI directly
        %%%                          for an angle-only sweep, or
        %%%                          DispersionViewerGUI directly for a
        %%%                          wavelength-only one)
            validateattributes(sweepResult, {'CalcWavevecAngSweep2Layered'}, ...
                {'scalar'}, 'DoubleSweepViewerGUI', 'sweepResult');
            if numel(sweepResult.angles) < 2 || numel(sweepResult.wavelengths) < 2
                error('DoubleSweepViewerGUI:constructor:notADoubleSweep', ...
                    ['DoubleSweepViewerGUI needs a genuine double sweep (numel(angles)>1 AND ' ...
                     'numel(wavelengths)>1); got %d angle(s), %d wavelength(s) -- use ' ...
                     'CombinedViewerGUI directly for an angle-only sweep, or DispersionViewerGUI ' ...
                     'directly for a wavelength-only one.'], ...
                    numel(sweepResult.angles), numel(sweepResult.wavelengths));
            end
            if isempty(sweepResult.sweeps)
                error('DoubleSweepViewerGUI:constructor:notCalculated', ...
                    'sweepResult has not been calc()''d yet -- nothing to view.');
            end

            obj.sweepResult  = sweepResult;
            obj.currentWavNo = 1;
            obj.currentAngNo = 1;

            obj.buildFigure();
            obj.buildAngularContent();
        end

        function show(obj)
        %%% Make the window visible, maximised (constructed invisible so headless/test
        %%% use is safe)
        %%% >> show(obj)
            GuiCommon.showMaximised(obj.fig);
        end

        function setAxisRange(obj, axisName, lims)
        %%% Set one axis range of the CURRENT mode's view from code
        %%% >> setAxisRange(obj, axisName, lims)    e.g. setAxisRange(obj, 'y', [NaN 7])
        %%%
        %%% Forwarded to the angular map viewer or, in dispersion mode, to the
        %%% dispersion viewer (no colour axis there: 'c' errors) -- see
        %%% GuiCommon.setAxisRange. The other mode keeps its own saved range.
            if strcmp(obj.mode, 'angular')
                obj.angularSweepViewer.setAxisRange(axisName, lims);
            else
                obj.dispersionViewer.setAxisRange(axisName, lims);
            end
        end
    end

    %%% ===================================================================
    %%% Figure construction
    %%% ===================================================================
    methods (Access = private)
        function buildFigure(obj)
            obj.fig = uifigure('Name', 'GuliweRR -- Double Sweep Viewer', ...
                'Position', [80 60 1200 900], 'Visible', 'off');
            mainLayout = uigridlayout(obj.fig, [3 1]);
            mainLayout.RowHeight = {36, '1x', 26};
            mainLayout.Padding = [4 4 4 4];

            topBar = uigridlayout(mainLayout, [1 2]);
            topBar.Layout.Row = 1;
            topBar.ColumnWidth = {260, '1x'};
            topBar.Padding = [0 0 0 0];
            obj.switchModeBtn = uibutton(topBar, 'Text', 'Switch to dispersion view ->', ...
                'ButtonPushedFcn', @(~,~) obj.onSwitchMode());

            %%% stable content slot -- buildAngularContent()/buildDispersionContent()
            %%% swap what's INSIDE it (never the slot itself), same pattern
            %%% AngSweepViewerGUI's own contentParent uses for its drill-down
            obj.contentParent = uigridlayout(mainLayout, [1 1]);
            obj.contentParent.Layout.Row = 2;
            obj.contentParent.Padding = [0 0 0 0];

            obj.statusLine = GuiCommon.statusLabel(mainLayout, 3);
        end

        function buildAngularContent(obj)
        %%% (Re)build angular-map mode: ModalShapeViewerGUI + AngSweepViewerGUI
        %%% on sweepResult.sweeps{currentWavNo}, wired exactly as
        %%% CombinedViewerGUI wires them -- see class docblock.
        %%% >> buildAngularContent(obj)
            delete(obj.contentParent.Children);

            innerLayout = uigridlayout(obj.contentParent, [2 1]);
            innerLayout.RowHeight = {36, '1x'};
            innerLayout.Padding = [0 0 0 0];

            wavRow = uigridlayout(innerLayout, [1 2]);
            wavRow.Layout.Row = 1;
            wavRow.ColumnWidth = {90, 160};
            wavRow.Padding = [0 0 0 0];
            uilabel(wavRow, 'Text', 'Wavelength:');
            wavStrs = arrayfun(@(w) sprintf('%.4g', w), obj.sweepResult.wavelengths, 'UniformOutput', false);
            uidropdown(wavRow, 'Items', wavStrs, 'ItemsData', 1:numel(obj.sweepResult.wavelengths), ...
                'Value', obj.currentWavNo, 'ValueChangedFcn', @(src,~) obj.onWavChanged(src.Value));

            %%% same 44/56 top/bottom split + thin divider CombinedViewerGUI's
            %%% own 'combined' windowMode uses
            [shapeParent, controlsParent, bottomParent, splitLayout] = ...
                GuiCommon.dividedTopBottom(innerLayout, {'44x', 3, '56x'}, [0 0 0 0]);
            splitLayout.Layout.Row = 2;   % below the wavelength stepper row

            sliceResult = obj.sweepResult.sweeps{obj.currentWavNo};

            obj.angularModalShapeViewer = ModalShapeViewerGUI(sliceResult, [], controlsParent, shapeParent);
            %%% the bar classifies the whole double sweep, not this one slice
            %%% (see class docblock)
            obj.angularSweepViewer      = AngSweepViewerGUI(sliceResult, bottomParent, obj.sweepResult);

            %%% the same bidirectional wiring CombinedViewerGUI performs, plus
            %%% tracking the selected angle as obj.currentAngNo
            obj.angularSweepViewer.OnModeSelected = ...
                @(angle, modeNo, velSnap, velRaw) obj.onAngularModeSelected(angle, modeNo, velRaw);
            obj.angularModalShapeViewer.OnSelectionChanged = ...
                @(angle, vel) obj.onAngularSelectionChanged(angle, vel);
            %%% after a Reclassify: the shape viewer's refresh (see class docblock)
            obj.angularSweepViewer.OnClassificationChanged = ...
                @() obj.angularModalShapeViewer.refreshClassification();

            GuiCommon.restoreAxisRange(obj.angularSweepViewer, obj.savedAxisRangeAngular);

            obj.switchModeBtn.Text = 'Switch to dispersion view ->';
            obj.mode = 'angular';
            GuiCommon.setStatusOk(obj.statusLine, sprintf( ...
                'Angular map at wavelength = %.4g', obj.sweepResult.wavelengths(obj.currentWavNo)));
        end

        function buildDispersionContent(obj)
        %%% (Re)build dispersion-curve mode: DispersionViewerGUI embedded via its
        %%% parentContainer support, angleIdx=currentAngNo.
        %%% >> buildDispersionContent(obj)
            delete(obj.contentParent.Children);

            obj.dispersionViewer = DispersionViewerGUI(obj.sweepResult, 'combined', obj.currentAngNo, obj.contentParent);
            %%% DispersionViewerGUI's constructor always goTo(1,1), so override
            %%% the wavelength back to whatever was last active — otherwise a
            %%% mode switch silently resets it. Branch index has no cross-mode
            %%% equivalent to restore, so it stays at 1.
            obj.dispersionViewer.goTo(obj.currentWavNo, 1);
            GuiCommon.restoreAxisRange(obj.dispersionViewer, obj.savedAxisRangeDispersion);

            obj.switchModeBtn.Text = '<- Switch to angular view';
            obj.mode = 'dispersion';
            GuiCommon.setStatusOk(obj.statusLine, sprintf( ...
                'Dispersion curve at angle = %.4g deg', obj.sweepResult.angles(obj.currentAngNo)));
        end
    end

    %%% ===================================================================
    %%% Public entry point + callbacks -- public (not Access=private),
    %%% matching AngSweepViewerGUI/DispersionViewerGUI's own convention so
    %%% tests can call these directly.
    %%% ===================================================================
    methods
        function onSwitchMode(obj)
        %%% Switch between angular-map and dispersion-curve mode, carrying over
        %%% the wavelength/angle active in whichever mode is being left, and
        %%% saving/restoring each mode's own axis-range fields.
        %%% >> onSwitchMode(obj)
            if strcmp(obj.mode, 'angular')
                if obj.angularSweepViewer.isDrilledIn
                    %%% xMinField/yMinField/cMinField are DESTROYED while
                    %%% AngSweepViewerGUI is drilled into its own single-angle
                    %%% view (delete(contentParent.Children) in drillDownToAngle,
                    %%% with nothing resetting the property references), so
                    %%% reading .Value on them below would throw "Invalid or
                    %%% deleted object". Return to the normal sweep view first,
                    %%% as clicking "<< Back to angular sweep" would.
                    %%% obj.currentAngNo is THIS class's own and is updated on
                    %%% every map click and shape-stepper move regardless of
                    %%% drill state, so no angle is lost by undrilling; only the
                    %%% drilled-in curve's own axis range goes, and neither of
                    %%% its axes match the 2D map's angle/velocity meaning anyway.
                    obj.angularSweepViewer.returnToSweepView();
                end
                obj.savedAxisRangeAngular = GuiCommon.captureAxisRange(obj.angularSweepViewer);
                obj.buildDispersionContent();
            else
                %%% no cMin/cMax here -- DispersionViewerGUI has no colour axis
                %%% to control, which captureAxisRange detects for itself
                obj.savedAxisRangeDispersion = GuiCommon.captureAxisRange(obj.dispersionViewer);
                obj.currentWavNo = obj.dispersionViewer.currentWavNo;
                obj.buildAngularContent();
            end
        end

        function onWavChanged(obj, newWavNo)
        %%% Angular-mode wavelength dropdown callback -- rebuilds against a
        %%% different sweeps{iWav} element (see buildAngularContent).
        %%% >> onWavChanged(obj, newWavNo)
            %%% keep the map's current ranges across the slice rebuild, typed or
            %%% set by setAxisRange (the fields are dead while drilled in -- the
            %%% last saved range is then kept instead)
            if ~obj.angularSweepViewer.isDrilledIn
                obj.savedAxisRangeAngular = GuiCommon.captureAxisRange(obj.angularSweepViewer);
            end
            obj.currentWavNo = newWavNo;
            obj.buildAngularContent();
        end

        function onAngularModeSelected(obj, angle, modeNo, velRaw)
        %%% AngSweepViewerGUI.OnModeSelected handler -- forwards to the modal-shape
        %%% viewer AND tracks the clicked angle as obj.currentAngNo, so a later
        %%% switch to dispersion mode lands on the selected angle.
        %%% >> onAngularModeSelected(obj, angle, modeNo, velRaw)
            [~, obj.currentAngNo] = min(abs(obj.sweepResult.angles - angle));
            obj.angularModalShapeViewer.goToSelection(angle, modeNo, velRaw);
        end

        function onAngularSelectionChanged(obj, angle, vel)
        %%% ModalShapeViewerGUI.OnSelectionChanged handler -- forwards to the
        %%% map's crosshair (CombinedViewerGUI's own wiring) AND tracks the
        %%% angle as obj.currentAngNo, same reasoning as
        %%% onAngularModeSelected (angle stepping via the shape viewer's own
        %%% prev/next buttons counts as "selecting" an angle too).
        %%% >> onAngularSelectionChanged(obj, angle, vel)
            [~, obj.currentAngNo] = min(abs(obj.sweepResult.angles - angle));
            obj.angularSweepViewer.showCrosshairAt(angle, vel);
        end
    end
end
