function obj = angSweepHS(varargin)
%%% Thin Calc namespace launcher for CalcAngSweepHS -- see CalcAngSweepHS.m.
%%% >> obj = Calc.angSweepHS(structure, angles, calcMode, retentionLevel, nVelGreen, useParallel)

    obj = CalcAngSweepHS(varargin{:});
end
