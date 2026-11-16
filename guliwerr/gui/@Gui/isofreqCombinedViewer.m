function obj = isofreqCombinedViewer(varargin)
%%% Thin Gui namespace launcher for IsofreqCombinedViewerGUI -- see
%%% IsofreqCombinedViewerGUI.m.
%%% >> obj = Gui.isofreqCombinedViewer(sweepResult, windowMode)

    obj = IsofreqCombinedViewerGUI(varargin{:});
end
