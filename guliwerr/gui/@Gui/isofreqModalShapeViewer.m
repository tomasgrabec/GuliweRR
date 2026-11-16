function obj = isofreqModalShapeViewer(varargin)
%%% Thin Gui namespace launcher for IsofreqModalShapeViewerGUI -- see
%%% IsofreqModalShapeViewerGUI.m.
%%% >> obj = Gui.isofreqModalShapeViewer(sweepResult, windowMode, controlsParent, shapeParent)

    obj = IsofreqModalShapeViewerGUI(varargin{:});
end
