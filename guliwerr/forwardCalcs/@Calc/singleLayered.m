function obj = singleLayered(varargin)
%%% Thin Calc namespace launcher for CalcSingleLayered -- see
%%% CalcSingleLayered.m.
%%% >> obj = Calc.singleLayered(structure, angle, vMax, precompStruct, calcMode, retentionLevel)

    obj = CalcSingleLayered(varargin{:});
end
