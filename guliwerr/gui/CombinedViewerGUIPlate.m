classdef CombinedViewerGUIPlate < handle
%%% CombinedViewerGUIPlate -- combination of ModalShapeViewerGUIPlate and
%%% AngSweepViewerGUIPlate, for CalcSinglePlate/
%%% CalcAngSweepPlate results. Plate counterpart of CombinedViewerGUI, same
%%% windowMode toggle, bidirectional wiring and single-angle-result dispatch
%%% -- see CombinedViewerGUI's own docblock for the full design rationale.
%%%
%%% A DELIBERATELY SEPARATE class from CombinedViewerGUI, not a
%%% class-dispatch branch inside it -- same "separate parallel sibling"
%%% convention as every other Plate GUI class.
%%%
%%%   CombinedViewerGUIPlate(sweepResult)              -- ONE window
%%%       (default, 'combined')
%%%   CombinedViewerGUIPlate(sweepResult, 'split')     -- TWO windows
%%%   CombinedViewerGUIPlate(sweepResult, 'fullSplit') -- THREE windows
%%%
%%% obj.fig is only set in 'combined' mode.
%%%
%%% SINGLE-ANGLE RESULTS: sweepResult may also be a CalcSinglePlate (already
%%% calc()'d) -- obj.angSweepViewer then holds a SingleAngleViewerGUIPlate
%%% instance instead (1D velocity curve rather than a 2D angle x velocity
%%% map); obj.modalShapeViewer is unaffected since ModalShapeViewerGUIPlate
%%% already accepts both result families directly.
%%%
%%% GuliweRR v6.0

    properties
        fig
        modalShapeViewer   % ModalShapeViewerGUIPlate instance
        angSweepViewer     % AngSweepViewerGUIPlate instance (angular-sweep sweepResult) OR
                            % SingleAngleViewerGUIPlate instance (single-angle sweepResult)
    end

    properties (SetAccess = private)
        windowMode   % 'combined' (default) | 'split' | 'fullSplit'
    end

    methods
        function obj = CombinedViewerGUIPlate(sweepResult, windowMode)
        %%% Construct the Plate combined viewer
        %%% >> obj = CombinedViewerGUIPlate(sweepResult, windowMode)
        %%%
        %%% INPUT: 1. sweepResult .. CalcAngSweepPlate, already calc()'d --
        %%%                          OR CalcSinglePlate, already calc()'d
        %%%                          (see class docblock's "SINGLE-ANGLE
        %%%                          RESULTS" note)
        %%%  (opt) 2. windowMode  .. 'combined' (default) | 'split' | 'fullSplit'
            if nargin < 2 || isempty(windowMode); windowMode = 'combined'; end
            GuiCommon.validateWindowMode('CombinedViewerGUIPlate', windowMode, ...
                {'combined', 'split', 'fullSplit'});
            obj.windowMode = windowMode;

            switch windowMode
                case 'combined'
                    obj.fig = uifigure('Name', 'GuliweRR -- Combined Viewer (Plate)', ...
                        'Position', [80 60 1100 900], 'Visible', 'off');
                    [shapeParent, controlsParent, bottomParent] = ...
                        GuiCommon.dividedTopBottom(obj.fig, {'44x', 3, '56x'}, [4 4 4 4]);

                    obj.modalShapeViewer = Gui.modalShapeViewerPlate(sweepResult, [], controlsParent, shapeParent);
                    obj.angSweepViewer   = obj.buildBottomViewer(sweepResult, bottomParent);

                case 'split'
                    obj.fig = [];
                    obj.modalShapeViewer = Gui.modalShapeViewerPlate(sweepResult);
                    obj.angSweepViewer   = obj.buildBottomViewer(sweepResult);
                    obj.modalShapeViewer.fig.Position(1:2) = [80, 460];
                    obj.angSweepViewer.fig.Position(1:2)   = [80, 60];

                case 'fullSplit'
                    obj.fig = [];
                    obj.modalShapeViewer = Gui.modalShapeViewerPlate(sweepResult, 'split');
                    obj.angSweepViewer   = obj.buildBottomViewer(sweepResult);
                    obj.modalShapeViewer.fig.Position(1:2)      = [80, 460];
                    obj.modalShapeViewer.shapeFig.Position(1:2) = [720, 460];
                    obj.angSweepViewer.fig.Position(1:2)        = [80, 60];
            end

            obj.angSweepViewer.OnModeSelected = ...
                @(angle, modeNo, velSnap, velRaw) obj.modalShapeViewer.goToSelection(angle, modeNo, velRaw);
            obj.modalShapeViewer.OnSelectionChanged = ...
                @(angle, vel) obj.angSweepViewer.showCrosshairAt(angle, vel);
        end

        function show(obj)
        %%% Make every window visible (constructed invisible so headless/test
        %%% use is safe).
        %%% >> show(obj)
            obj.modalShapeViewer.show();
            obj.angSweepViewer.show();
            if ~isempty(obj.fig)   % 'combined': the one shared window opens maximised
                GuiCommon.showMaximised(obj.fig);
            end
        end

        function setAxisRange(obj, axisName, lims)
        %%% Set one axis range of the bottom map/curve viewer from code
        %%% >> setAxisRange(obj, axisName, lims)    e.g. setAxisRange(obj, 'c', [-3 0])
        %%%
        %%% Forwarded to angSweepViewer.setAxisRange -- see GuiCommon.setAxisRange
            obj.angSweepViewer.setAxisRange(axisName, lims);
        end
    end

    methods (Access = private)
        function viewer = buildBottomViewer(~, sweepResult, parentContainer)
        %%% Dispatch the bottom-window viewer class on sweepResult's own type --
        %%% angular sweep (CalcAngSweepPlate) -> AngSweepViewerGUIPlate, single
        %%% angle (CalcSinglePlate) -> SingleAngleViewerGUIPlate. Both target
        %%% classes implement the same small duck-typed interface (.fig,
        %%% .show(), .OnModeSelected, .showCrosshairAt(angle, vel)), so the
        %%% bidirectional-linking code right after the switch block above
        %%% works unmodified regardless of which one obj.angSweepViewer ends
        %%% up holding.
        %%% >> viewer = buildBottomViewer(obj, sweepResult, parentContainer)
            if nargin < 3; parentContainer = []; end
            if isa(sweepResult, 'CalcSinglePlate')
                viewer = Gui.singleAngleViewerPlate(sweepResult, parentContainer);
            else
                viewer = Gui.angSweepViewerPlate(sweepResult, parentContainer);
            end
        end
    end
end
