function obj = angSweepLayered(varargin)
%%% Thin Calc namespace launcher for CalcAngSweepLayered -- see
%%% CalcAngSweepLayered.m.
%%% >> obj = Calc.angSweepLayered(structure, angles, vMax, calcMode, retentionLevel, [], rawByDepthCell, useParallel)

    obj = CalcAngSweepLayered(varargin{:});
end
