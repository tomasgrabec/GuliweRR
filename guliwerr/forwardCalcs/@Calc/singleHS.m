function obj = singleHS(varargin)
%%% Thin Calc namespace launcher for CalcSingleHS -- see CalcSingleHS.m.
%%% >> obj = Calc.singleHS(structure, angle, precompStruct, calcMode, retentionLevel, nVelGreen)

    obj = CalcSingleHS(varargin{:});
end
