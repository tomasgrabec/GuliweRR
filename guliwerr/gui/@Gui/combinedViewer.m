function obj = combinedViewer(varargin)
%%% Thin Gui namespace launcher for CombinedViewerGUI -- see
%%% CombinedViewerGUI.m.
%%% >> obj = Gui.combinedViewer(sweepResult, windowMode)

    obj = CombinedViewerGUI(varargin{:});
end
