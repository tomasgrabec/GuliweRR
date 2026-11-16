function obj = structureControl(varargin)
%%% Thin Gui namespace launcher for StructureControlGUI -- see
%%% StructureControlGUI.m.
%%% >> obj = Gui.structureControl(initialStructure)

    obj = StructureControlGUI(varargin{:});
end
