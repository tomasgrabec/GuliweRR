function obj = combinedViewerPlate(varargin)
%%% Thin Gui namespace launcher for CombinedViewerGUIPlate -- see
%%% CombinedViewerGUIPlate.m.
%%% >> obj = Gui.combinedViewerPlate(sweepResult, windowMode)

    obj = CombinedViewerGUIPlate(varargin{:});
end
