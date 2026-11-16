function [cVals, cLabel] = visClassificationColorVals(colorBy, vels, scoreDb, greenPeakVel)
%%% Continuous colour values for a classification plot's 'score'/'confidence'
%%% colouring
%%% >> [cVals, cLabel] = visClassificationColorVals(colorBy, vels, scoreDb, greenPeakVel)
%%%
%%% @Vis-private helper shared by Vis.plotClassificationDiagnostic and
%%% Vis.plotClassificationCurve, so the same word ('score', 'confidence')
%%% means the same number in both.
%%%
%%% INPUT: 1. colorBy      .. 'score' | 'confidence'
%%%        2. vels         .. per-point bound Ritz mode velocity (the entries'
%%%                           modeVel; read by 'confidence' only)
%%%        3. scoreDb      .. per-point classification score in dB (the entry's own
%%%                           scoreDb, on its method's scale; never recomputed here)
%%%        4. greenPeakVel .. matched Green peak velocity, NaN where the
%%%                           classifier made no Green match
%%% OUTPUT: cVals  .. per-point colour value (NaN allowed -- scatter simply
%%%                   leaves such a point uncoloured)
%%%         cLabel .. colour-bar label
%%%
%%% 'score': the score the classification thresholded, in dB (NaN on limiting
%%% entries). Comparable only within one method -- a plot shows one result's active
%%% classification, so all its points share one.
%%% 'confidence': how far the bound Ritz mode sits from the Green peak it was
%%% matched to, abs(modeVel - greenPeakVel) -- so 0 is a perfect match and NaN
%%% means "never matched", NOT "matched badly". Not the entry's vel: a
%%% 'green'/'refined' entry's vel IS its greenPeakVel.

    if strcmp(colorBy, 'score')
        cVals  = scoreDb;
        cLabel = 'score [dB]';
    else
        cVals  = abs(vels - greenPeakVel);   % NaN where no Green match
        cLabel = '|v_{mode} - v_{Green peak}| [km/s]';
    end
end
