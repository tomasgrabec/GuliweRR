function obj = angSweepPlate(varargin)
%%% Thin Calc namespace launcher for CalcAngSweepPlate -- see CalcAngSweepPlate.m.
%%% >> obj = Calc.angSweepPlate(structure, angles, realThickness, realWavVec, retentionLevel, useParallel)

    obj = CalcAngSweepPlate(varargin{:});
end
