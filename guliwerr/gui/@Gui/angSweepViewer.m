function obj = angSweepViewer(varargin)
%%% Thin Gui namespace launcher for AngSweepViewerGUI -- see
%%% AngSweepViewerGUI.m.
%%% >> obj = Gui.angSweepViewer(sweepResult, parentContainer, classifyTarget)

    obj = AngSweepViewerGUI(varargin{:});
end
