function obj = singlePlate(varargin)
%%% Thin Calc namespace launcher for CalcSinglePlate -- see CalcSinglePlate.m.
%%% >> obj = Calc.singlePlate(structure, angle, realThickness, realWavVec, precompStruct, retentionLevel)

    obj = CalcSinglePlate(varargin{:});
end
