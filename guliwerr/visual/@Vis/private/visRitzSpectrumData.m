function [absS, velAxis] = visRitzSpectrumData(singleResult, velAxis)
%%% The Ritz spectrum |S(v)| of one angle, rebuilt from the stored modal numerators
%%% >> [absS, velAxis] = visRitzSpectrumData(singleResult)
%%% >> absS = visRitzSpectrumData(singleResult, velAxis)
%%%
%%% @Vis-private helper: the one evaluation of |S(v)| for the Ritz spectrum curve, the
%%% Ritz spectrum map (each angle on a common axis) and the classification curve's
%%% projection, so all three draw the same numbers.
%%%
%%% INPUT:  1. singleResult .. CalcSingleHS or CalcSingleLayered with stored numerators
%%%                            (check with visHasRitzSpectrum first)
%%%  (opt)  2. velAxis      .. velocities [km/s] to evaluate on (default: the result's
%%%                            own ritzScoreVelAxis)
%%% OUTPUT: absS    .. |RitzScore.spectrum| on velAxis, column
%%%         velAxis .. the axis used, column
    num = singleResult.ritzScoreNum;
    if nargin < 2 || isempty(velAxis)
        velAxis = singleResult.ritzScoreVelAxis;
    end
    velAxis = velAxis(:);
    absS = abs(RitzScore.spectrum(singleResult.ritzSolver.vels(1:numel(num)), num, velAxis));
end
