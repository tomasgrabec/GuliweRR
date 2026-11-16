classdef IsofreqCombinedViewerGUI < handle
%%% IsofreqCombinedViewerGUI -- combination of IsofreqModalShapeViewerGUI
%%% (top) and IsofreqAngSweepViewerGUI/IsofreqSingleLayeredViewerGUI
%%% (bottom, dispatched on the result's class), the isoFreq counterpart of
%%% CombinedViewerGUI. Line-for-line mirror of that class's windowMode
%%% (combined/split/fullSplit) and bidirectional-wiring mechanics, with the
%%% renamed OnPointSelected/goToPoint/OnSelectionChanged/showCrosshairAt
%%% interface (docs/reference/architecture/gui_viewers.md's isoFreq family section --
%%% kParRange is a continuous axis, not a discrete Ritz mode index) in place
%%% of CombinedViewerGUI's OnModeSelected/goToSelection pair.
%%%
%%% GuliweRR v6.0

    properties
        fig
        modalShapeViewer   % IsofreqModalShapeViewerGUI instance
        angSweepViewer     % IsofreqAngSweepViewerGUI instance (angular-sweep sweepResult) OR
                            % IsofreqSingleLayeredViewerGUI instance (single-angle sweepResult)
    end

    properties (SetAccess = private)
        windowMode   % 'combined' (default) | 'split' | 'fullSplit' -- see CombinedViewerGUI's own docblock
    end

    methods
        function obj = IsofreqCombinedViewerGUI(sweepResult, windowMode)
        %%% Construct the isoFreq combined viewer
        %%% >> obj = IsofreqCombinedViewerGUI(sweepResult, windowMode)
        %%%
        %%% INPUT: 1. sweepResult .. CalcIsofreqAngSweepLayered, already
        %%%                          calc()'d -- OR CalcIsofreqSingleLayered,
        %%%                          already calc()'d
        %%%  (opt) 2. windowMode  .. 'combined' (default) | 'split' | 'fullSplit'
        %%%                          -- see CombinedViewerGUI's own docblock
            if nargin < 2 || isempty(windowMode); windowMode = 'combined'; end
            GuiCommon.validateWindowMode('IsofreqCombinedViewerGUI', windowMode, ...
                {'combined', 'split', 'fullSplit'});
            obj.windowMode = windowMode;

            switch windowMode
                case 'combined'
                    obj.fig = uifigure('Name', 'GuliweRR -- isoFreq Combined Viewer', ...
                        'Position', [80 60 1100 900], 'Visible', 'off');
                    %%% 34/66 rather than CombinedViewerGUI's 44/56 -- the
                    %%% isoFreq shape half is a single Green's-response row, not
                    %%% a {3 x nLay} stack, so it needs less of the window
                    [shapeParent, controlsParent, bottomParent] = ...
                        GuiCommon.dividedTopBottom(obj.fig, {'34x', 3, '66x'}, [4 4 4 4]);

                    obj.modalShapeViewer = Gui.isofreqModalShapeViewer(sweepResult, [], controlsParent, shapeParent);
                    obj.angSweepViewer   = obj.buildBottomViewer(sweepResult, bottomParent);

                case 'split'
                    obj.fig = [];
                    obj.modalShapeViewer = Gui.isofreqModalShapeViewer(sweepResult);   % standalone, own default 'combined' -> one window
                    obj.angSweepViewer   = obj.buildBottomViewer(sweepResult);          % standalone -> one window
                    obj.modalShapeViewer.fig.Position(1:2) = [80, 460];
                    obj.angSweepViewer.fig.Position(1:2)   = [80, 60];

                case 'fullSplit'
                    obj.fig = [];
                    obj.modalShapeViewer = Gui.isofreqModalShapeViewer(sweepResult, 'split');   % standalone, two windows
                    obj.angSweepViewer   = obj.buildBottomViewer(sweepResult);                   % standalone -> one window
                    obj.modalShapeViewer.fig.Position(1:2)      = [80, 460];
                    obj.modalShapeViewer.shapeFig.Position(1:2) = [720, 460];
                    obj.angSweepViewer.fig.Position(1:2)        = [80, 60];
            end

            %%% automatic bidirectional linking, identical in every mode --
            %%% map/curve click drives the shape viewer; the shape viewer's
            %%% own state changes (angle/kParRange stepping) mark the
            %%% map/curve crosshair in return. Same mechanism as
            %%% CombinedViewerGUI, renamed interface (class docblock).
            obj.angSweepViewer.OnPointSelected = ...
                @(angle, kParVal) obj.modalShapeViewer.goToPoint(angle, kParVal);
            obj.modalShapeViewer.OnSelectionChanged = ...
                @(angle, kParVal) obj.angSweepViewer.showCrosshairAt(angle, kParVal);
        end

        function show(obj)
        %%% Make every window visible (constructed invisible so headless/test
        %%% use is safe) -- delegates to each sub-viewer's own show().
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
        %%% Dispatch the bottom-window viewer class on sweepResult's own
        %%% type -- angular sweep (CalcIsofreqAngSweepLayered) ->
        %%% IsofreqAngSweepViewerGUI (2D heatmap), single angle
        %%% (CalcIsofreqSingleLayered) -> IsofreqSingleLayeredViewerGUI (1D
        %%% curve). Same isa(sweepResult, ...) dispatch CombinedViewerGUI's
        %%% own buildBottomViewer uses.
        %%% >> viewer = buildBottomViewer(obj, sweepResult, parentContainer)
            if nargin < 3; parentContainer = []; end
            if isa(sweepResult, 'CalcIsofreqSingleLayered')
                viewer = Gui.isofreqSingleLayeredViewer(sweepResult, parentContainer);
            else
                viewer = Gui.isofreqAngSweepViewer(sweepResult, parentContainer);
            end
        end
    end
end
