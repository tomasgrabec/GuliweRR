function obj = dispersionViewerPlate(varargin)
%%% Thin Gui namespace launcher for DispersionViewerGUIPlate -- see
%%% DispersionViewerGUIPlate.m.
%%% >> obj = Gui.dispersionViewerPlate(sweepResult, windowMode)

    obj = DispersionViewerGUIPlate(varargin{:});
end
