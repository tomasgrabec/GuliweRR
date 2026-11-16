function obj = modalShapeViewerPlate(varargin)
%%% Thin Gui namespace launcher for ModalShapeViewerGUIPlate -- see
%%% ModalShapeViewerGUIPlate.m.
%%% >> obj = Gui.modalShapeViewerPlate(sweepResult, windowMode, controlsParent, shapeParent)

    obj = ModalShapeViewerGUIPlate(varargin{:});
end
