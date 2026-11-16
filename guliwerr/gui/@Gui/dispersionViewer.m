function obj = dispersionViewer(varargin)
%%% Thin Gui namespace launcher for DispersionViewerGUI -- see
%%% DispersionViewerGUI.m.
%%% >> obj = Gui.dispersionViewer(sweepResult, windowMode)

    obj = DispersionViewerGUI(varargin{:});
end
