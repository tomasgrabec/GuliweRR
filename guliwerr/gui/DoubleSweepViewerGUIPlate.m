classdef DoubleSweepViewerGUIPlate < handle
%%% DoubleSweepViewerGUIPlate -- interactive viewer for a genuine DOUBLE
%%% sweep (CalcWavevecAngSweep2Plate with BOTH numel(angles)>1 AND
%%% numel(wavVecs)>1), switching between an angular-map view (at a chosen
%%% fixed wavevector) and a dispersion-curve view (at a chosen fixed angle).
%%% Plate counterpart of DoubleSweepViewerGUI -- see that class for the design
%%% rationale (the angular/dispersion mode split, axis-range save/restore).
%%% This follows it mechanically, against the Plate sibling classes and Plate's
%%% `wavVecs` sweep axis in place of Layered's `wavelengths`.
%%%
%%% A SEPARATE class from DoubleSweepViewerGUI, not a class-dispatch branch
%%% inside it — the same parallel-sibling convention every Plate GUI class uses.
%%%
%%% Angular mode: ModalShapeViewerGUIPlate + AngSweepViewerGUIPlate, wired as
%%% CombinedViewerGUIPlate wires them, pointed at
%%% sweepResult.sweeps{currentWavNo} — already a complete, independently
%%% calc()'d CalcAngSweepPlate.
%%%
%%% Dispersion mode: DispersionViewerGUIPlate embedded via its own
%%% parentContainer support, pointed at sweepResult directly with
%%% angleIdx=currentAngNo.
%%%
%%% RAM warning: dispersion mode needs sweepResult calc()'d with
%%% retentionLevel='full' (DispersionBranches.track's requirement, identical for
%%% Plate). NOT validated at construction, for the same reason
%%% DoubleSweepViewerGUI gives.
%%%
%%% Single window only, as DoubleSweepViewerGUI is.
%%%
%%% GuliweRR v6.0

    properties
        fig
        switchModeBtn
        statusLine
    end

    properties (SetAccess = private)
        sweepResult    % CalcWavevecAngSweep2Plate, calc()'d, numel(angles)>1 AND numel(wavVecs)>1
        mode           % 'angular' | 'dispersion'
        currentWavNo   % index into sweepResult.wavVecs -- used by angular mode
        currentAngNo   % index into sweepResult.angles -- used by dispersion mode
    end

    properties (Access = {?DoubleSweepViewerGUIPlate, ?test_DoubleSweepViewerGUIPlate})
        contentParent
        angularModalShapeViewer = []   % ModalShapeViewerGUIPlate, while mode=='angular'
        angularSweepViewer      = []   % AngSweepViewerGUIPlate, while mode=='angular'
        dispersionViewer        = []   % DispersionViewerGUIPlate, while mode=='dispersion'
        savedAxisRangeAngular    = []
        savedAxisRangeDispersion = []
    end

    methods
        function obj = DoubleSweepViewerGUIPlate(sweepResult)
        %%% Construct the Plate double-sweep viewer
        %%% >> obj = DoubleSweepViewerGUIPlate(sweepResult)
        %%%
        %%% INPUT: 1. sweepResult .. CalcWavevecAngSweep2Plate, already
        %%%                          calc()'d, with numel(angles)>1 AND
        %%%                          numel(wavVecs)>1 (a genuine double
        %%%                          sweep -- use CombinedViewerGUIPlate
        %%%                          directly for an angle-only sweep, or
        %%%                          DispersionViewerGUIPlate directly for a
        %%%                          wavevector-only one)
            validateattributes(sweepResult, {'CalcWavevecAngSweep2Plate'}, ...
                {'scalar'}, 'DoubleSweepViewerGUIPlate', 'sweepResult');
            if numel(sweepResult.angles) < 2 || numel(sweepResult.wavVecs) < 2
                error('DoubleSweepViewerGUIPlate:constructor:notADoubleSweep', ...
                    ['DoubleSweepViewerGUIPlate needs a genuine double sweep (numel(angles)>1 AND ' ...
                     'numel(wavVecs)>1); got %d angle(s), %d wavevector(s) -- use ' ...
                     'CombinedViewerGUIPlate directly for an angle-only sweep, or DispersionViewerGUIPlate ' ...
                     'directly for a wavevector-only one.'], ...
                    numel(sweepResult.angles), numel(sweepResult.wavVecs));
            end
            if isempty(sweepResult.sweeps)
                error('DoubleSweepViewerGUIPlate:constructor:notCalculated', ...
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
            obj.fig = uifigure('Name', 'GuliweRR -- Double Sweep Viewer (Plate)', ...
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

            obj.contentParent = uigridlayout(mainLayout, [1 1]);
            obj.contentParent.Layout.Row = 2;
            obj.contentParent.Padding = [0 0 0 0];

            obj.statusLine = GuiCommon.statusLabel(mainLayout, 3);
        end

        function buildAngularContent(obj)
        %%% (Re)build angular-map mode: ModalShapeViewerGUIPlate +
        %%% AngSweepViewerGUIPlate on sweepResult.sweeps{currentWavNo}, wired
        %%% exactly as CombinedViewerGUIPlate wires them.
        %%% >> buildAngularContent(obj)
            delete(obj.contentParent.Children);

            innerLayout = uigridlayout(obj.contentParent, [2 1]);
            innerLayout.RowHeight = {36, '1x'};
            innerLayout.Padding = [0 0 0 0];

            wavRow = uigridlayout(innerLayout, [1 2]);
            wavRow.Layout.Row = 1;
            wavRow.ColumnWidth = {90, 160};
            wavRow.Padding = [0 0 0 0];
            uilabel(wavRow, 'Text', 'Wavevector:');
            wavStrs = arrayfun(@(w) sprintf('%.4g', w), obj.sweepResult.wavVecs, 'UniformOutput', false);
            uidropdown(wavRow, 'Items', wavStrs, 'ItemsData', 1:numel(obj.sweepResult.wavVecs), ...
                'Value', obj.currentWavNo, 'ValueChangedFcn', @(src,~) obj.onWavChanged(src.Value));

            [shapeParent, controlsParent, bottomParent, splitLayout] = ...
                GuiCommon.dividedTopBottom(innerLayout, {'44x', 3, '56x'}, [0 0 0 0]);
            splitLayout.Layout.Row = 2;   % below the wavevector stepper row

            sliceResult = obj.sweepResult.sweeps{obj.currentWavNo};

            obj.angularModalShapeViewer = ModalShapeViewerGUIPlate(sliceResult, [], controlsParent, shapeParent);
            obj.angularSweepViewer      = AngSweepViewerGUIPlate(sliceResult, bottomParent);

            obj.angularSweepViewer.OnModeSelected = ...
                @(angle, modeNo, velSnap, velRaw) obj.onAngularModeSelected(angle, modeNo, velRaw);
            obj.angularModalShapeViewer.OnSelectionChanged = ...
                @(angle, vel) obj.onAngularSelectionChanged(angle, vel);

            GuiCommon.restoreAxisRange(obj.angularSweepViewer, obj.savedAxisRangeAngular);

            obj.switchModeBtn.Text = 'Switch to dispersion view ->';
            obj.mode = 'angular';
            GuiCommon.setStatusOk(obj.statusLine, sprintf( ...
                'Angular map at wavevector = %.4g', obj.sweepResult.wavVecs(obj.currentWavNo)));
        end

        function buildDispersionContent(obj)
        %%% (Re)build dispersion-curve mode: DispersionViewerGUIPlate
        %%% embedded via its own parentContainer support, angleIdx=currentAngNo.
        %%% >> buildDispersionContent(obj)
            delete(obj.contentParent.Children);

            obj.dispersionViewer = DispersionViewerGUIPlate(obj.sweepResult, 'combined', obj.currentAngNo, obj.contentParent);
            %%% DispersionViewerGUIPlate's own constructor always goTo(1,1) --
            %%% override the wavevector back to whatever was last active (own
            %%% currentWavNo, tracked independently of angle) so switching
            %%% modes doesn't silently reset wavevector/branch position;
            %%% branch index has no cross-mode equivalent to restore, so it
            %%% stays at 1.
            obj.dispersionViewer.goTo(obj.currentWavNo, 1);
            GuiCommon.restoreAxisRange(obj.dispersionViewer, obj.savedAxisRangeDispersion);

            obj.switchModeBtn.Text = '<- Switch to angular view';
            obj.mode = 'dispersion';
            GuiCommon.setStatusOk(obj.statusLine, sprintf( ...
                'Dispersion curve at angle = %.4g deg', obj.sweepResult.angles(obj.currentAngNo)));
        end
    end

    %%% ===================================================================
    %%% Public entry point + callbacks
    %%% ===================================================================
    methods
        function onSwitchMode(obj)
        %%% Switch between angular-map and dispersion-curve mode, carrying
        %%% over the wavevector/angle currently active, and saving/restoring
        %%% each mode's own axis-range fields -- see DoubleSweepViewerGUI's
        %%% own onSwitchMode docblock for the drilled-in-view edge case this
        %%% mirrors exactly (AngSweepViewerGUIPlate's own isDrilledIn).
        %%% >> onSwitchMode(obj)
            if strcmp(obj.mode, 'angular')
                if obj.angularSweepViewer.isDrilledIn
                    obj.angularSweepViewer.returnToSweepView();
                end
                obj.savedAxisRangeAngular = GuiCommon.captureAxisRange(obj.angularSweepViewer);
                obj.buildDispersionContent();
            else
                %%% no cMin/cMax here -- DispersionViewerGUIPlate has no colour
                %%% axis to control, which captureAxisRange detects for itself
                obj.savedAxisRangeDispersion = GuiCommon.captureAxisRange(obj.dispersionViewer);
                obj.currentWavNo = obj.dispersionViewer.currentWavNo;
                obj.buildAngularContent();
            end
        end

        function onWavChanged(obj, newWavNo)
        %%% Angular-mode wavevector dropdown callback -- rebuilds against a
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
        %%% AngSweepViewerGUIPlate.OnModeSelected handler -- forwards to the
        %%% modal-shape viewer AND tracks the clicked angle as
        %%% obj.currentAngNo, so a later switch to dispersion mode lands on
        %%% "the selected angle".
        %%% >> onAngularModeSelected(obj, angle, modeNo, velRaw)
            [~, obj.currentAngNo] = min(abs(obj.sweepResult.angles - angle));
            obj.angularModalShapeViewer.goToSelection(angle, modeNo, velRaw);
        end

        function onAngularSelectionChanged(obj, angle, vel)
        %%% ModalShapeViewerGUIPlate.OnSelectionChanged handler -- forwards
        %%% to the map's crosshair AND tracks the angle as obj.currentAngNo.
        %%% >> onAngularSelectionChanged(obj, angle, vel)
            [~, obj.currentAngNo] = min(abs(obj.sweepResult.angles - angle));
            obj.angularSweepViewer.showCrosshairAt(angle, vel);
        end
    end
end
