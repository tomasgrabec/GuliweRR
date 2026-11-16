function [calcResult, meta] = loadResult(varargin)
%%% Thin Calc namespace launcher for CalcCommon.loadCalc -- load a Calc*
%%% result saved by Calc.saveResult (or by a plain save())
%%% >> [calcResult, meta] = Calc.loadResult(filePath)
%%%
%%% See CalcCommon.loadCalc.

    [calcResult, meta] = CalcCommon.loadCalc(varargin{:});
end
