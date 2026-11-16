function obj = modalShapeViewer(varargin)
%%% Thin Gui namespace launcher for ModalShapeViewerGUI -- see
%%% ModalShapeViewerGUI.m.
%%% >> obj = Gui.modalShapeViewer(sweepResult, windowMode, controlsParent, shapeParent)

    obj = ModalShapeViewerGUI(varargin{:});
end
