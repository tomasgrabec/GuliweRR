classdef CombinedViewerGUI < handle
%%% CombinedViewerGUI -- combination of ModalShapeViewerGUI and
%%% AngSweepViewerGUI, with a windowMode toggle (2nd constructor argument) for
%%% how many physical windows that combination opens:
%%%
%%%   CombinedViewerGUI(sweepResult)              -- ONE window (default,
%%%       'combined'): modal shapes on top with their controls to the right,
%%%       angular-sweep map on the bottom, in a 44%/56% split.
%%%   CombinedViewerGUI(sweepResult, 'split')     -- TWO windows:
%%%       ModalShapeViewerGUI standalone/'combined' (shapes and controls
%%%       together) plus AngSweepViewerGUI standalone.
%%%   CombinedViewerGUI(sweepResult, 'fullSplit') -- THREE windows:
%%%       ModalShapeViewerGUI standalone/'split' (itself two windows, controls
%%%       and shapes separately) plus AngSweepViewerGUI standalone.
%%%
%%% obj.fig is set only in 'combined' mode; there is no single shared figure in
%%% 'split'/'fullSplit', where it stays []. Use show(), which delegates to each
%%% sub-viewer and is correct in every mode, rather than touching obj.fig.
%%%
%%% Both underlying classes stay fully independent in every mode. This class
%%% only supplies 'combined' mode with parent-container arguments and wires
%%% AngSweepViewerGUI.OnModeSelected / ModalShapeViewerGUI.OnSelectionChanged
%%% together — the same wiring a caller would do by hand for two standalone
%%% windows, needing only the two sub-viewer handles, not how many physical
%%% windows they occupy.
%%%
%%% SINGLE-ANGLE RESULTS: sweepResult may also be a CalcSingleHS/
%%% CalcSingleLayered (calcAll()'d, not calc()'d — those have no angle
%%% dimension to sweep). obj.angSweepViewer then holds a SingleAngleViewerGUI
%%% instead, a 1D velocity curve rather than a 2D angle x velocity map, while
%%% obj.modalShapeViewer is unaffected: ModalShapeViewerGUI accepts both result
%%% families directly. The property keeps one name across both cases, external
%%% code needing only the shared duck-typed interface (.fig/.show()/
%%% .OnModeSelected/.OnClassificationChanged/.showCrosshairAt), never the concrete
%%% class.
%%%
%%% RECLASSIFY: the bottom viewer's classification bar reclassifies the result; its
%%% OnClassificationChanged then refreshes the modal-shape viewer
%%% (ModalShapeViewerGUI.refreshClassification), which reads the classification
%%% through the result handle.
%%%
%%% GuliweRR v6.0

    properties
        fig
        modalShapeViewer   % ModalShapeViewerGUI instance
        angSweepViewer     % AngSweepViewerGUI instance (angular-sweep sweepResult) OR
                            % SingleAngleViewerGUI instance (single-angle sweepResult) --
                            % see class docblock's "SINGLE-ANGLE RESULTS" note
    end

    properties (SetAccess = private)
        windowMode   % 'combined' (default) | 'split' | 'fullSplit' -- see class docblock
    end

    methods
        function obj = CombinedViewerGUI(sweepResult, windowMode)
        %%% Construct the combined viewer
        %%% >> obj = CombinedViewerGUI(sweepResult, windowMode)
        %%%
        %%% INPUT: 1. sweepResult .. CalcAngSweepHS or CalcAngSweepLayered,
        %%%                          already calc()'d -- OR CalcSingleHS/
        %%%                          CalcSingleLayered, already calcAll()'d
        %%%                          (see class docblock's "SINGLE-ANGLE
        %%%                          RESULTS" note)
        %%%  (opt) 2. windowMode  .. 'combined' (default) | 'split' | 'fullSplit'
        %%%                          -- see class docblock
            if nargin < 2 || isempty(windowMode); windowMode = 'combined'; end
            GuiCommon.validateWindowMode('CombinedViewerGUI', windowMode, ...
                {'combined', 'split', 'fullSplit'});
            obj.windowMode = windowMode;

            switch windowMode
                case 'combined'
                    obj.fig = uifigure('Name', 'GuliweRR -- Combined Viewer', ...
                        'Position', [80 60 1100 900], 'Visible', 'off');
                    %%% shapes + their controls on top, a thin divider, the
                    %%% map (or a drilled-down single-angle curve) below --
                    %%% 44%/56% top/bottom split
                    [shapeParent, controlsParent, bottomParent] = ...
                        GuiCommon.dividedTopBottom(obj.fig, {'44x', 3, '56x'}, [4 4 4 4]);

                    obj.modalShapeViewer = Gui.modalShapeViewer(sweepResult, [], controlsParent, shapeParent);
                    obj.angSweepViewer   = obj.buildBottomViewer(sweepResult, bottomParent);

                case 'split'
                    obj.fig = [];
                    obj.modalShapeViewer = Gui.modalShapeViewer(sweepResult);   % standalone, own default 'combined' -> one window
                    obj.angSweepViewer   = obj.buildBottomViewer(sweepResult);  % standalone -> one window
                    obj.modalShapeViewer.fig.Position(1:2) = [80, 460];
                    obj.angSweepViewer.fig.Position(1:2)   = [80, 60];

                case 'fullSplit'
                    obj.fig = [];
                    obj.modalShapeViewer = Gui.modalShapeViewer(sweepResult, 'split');   % standalone, two windows
                    obj.angSweepViewer   = obj.buildBottomViewer(sweepResult);           % standalone -> one window
                    obj.modalShapeViewer.fig.Position(1:2)      = [80, 460];
                    obj.modalShapeViewer.shapeFig.Position(1:2) = [720, 460];
                    obj.angSweepViewer.fig.Position(1:2)        = [80, 60];
            end

            %%% automatic bidirectional linking, identical in every mode -- map
            %%% click drives the shape viewer (Ritz-mode snap or raw Green's
            %%% velocity, depending on the shape viewer's OWN current source);
            %%% the shape viewer's own state changes (angle/mode stepping, Go
            %%% To buttons, etc.) mark the map's crosshair in return
            obj.angSweepViewer.OnModeSelected = ...
                @(angle, modeNo, velSnap, velRaw) obj.modalShapeViewer.goToSelection(angle, modeNo, velRaw);
            obj.modalShapeViewer.OnSelectionChanged = ...
                @(angle, vel) obj.angSweepViewer.showCrosshairAt(angle, vel);
            %%% refresh chain of a Reclassify in the bottom viewer's classification
            %%% bar: it redraws itself (and its drilled-in curve), then the shape
            %%% viewer rebuilds its SAW / limiting panel
            obj.angSweepViewer.OnClassificationChanged = ...
                @() obj.modalShapeViewer.refreshClassification();
        end

        function show(obj)
        %%% Make every window visible (constructed invisible so headless/test
        %%% use is safe) -- delegates to each sub-viewer's own show(), which
        %%% already correctly handles however many figures it owns, so this
        %%% one call covers every windowMode uniformly.
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
        %%% angular sweep (CalcAngSweepHS/CalcAngSweepLayered) -> AngSweepViewerGUI
        %%% (2D angle x velocity map), single angle (CalcSingleHS/
        %%% CalcSingleLayered) -> SingleAngleViewerGUI (1D velocity curve). Same
        %%% isa(sweepResult, ...) dispatch, same convention used throughout this
        %%% GUI layer (see docs/reference/architecture/gui_viewers.md) -- ModalShapeViewerGUI's own top-half
        %%% construction needs no equivalent branch here since it already
        %%% accepts both result families directly (see its own docblock).
        %%% Both target classes implement the same small duck-typed interface
        %%% (.fig, .show(), .OnModeSelected, .showCrosshairAt(angle, vel)), so
        %%% the bidirectional-linking code right after the switch block below
        %%% works unmodified regardless of which one obj.angSweepViewer ends up
        %%% holding.
        %%% >> viewer = buildBottomViewer(obj, sweepResult, parentContainer)
            if nargin < 3; parentContainer = []; end
            if isa(sweepResult, 'CalcSingleHS') || isa(sweepResult, 'CalcSingleLayered')
                viewer = Gui.singleAngleViewer(sweepResult, parentContainer);
            else
                viewer = Gui.angSweepViewer(sweepResult, parentContainer);
            end
        end
    end
end
