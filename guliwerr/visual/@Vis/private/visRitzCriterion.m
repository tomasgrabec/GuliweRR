function [vals, info] = visRitzCriterion(singleResult, criterion)
%%% Per-mode values of a Ritz-mode criterion, and how to show them
%%% >> [vals, info] = visRitzCriterion(singleResult, criterion)
%%% >> names = visRitzCriterion()
%%%
%%% @Vis-private helper shared by Vis.plotRitzModeMap (every angle of a sweep) and
%%% Vis.plotRitzModeCurve (one angle), so a criterion means the same number in both.
%%% Called without arguments it returns the criterion names, the list both plot
%%% functions validate against.
%%%
%%% INPUT: 1. singleResult .. CalcSingleHS or CalcSingleLayered, calcAll()'d
%%%        2. criterion    .. 'ritzScore' | 'freqDisp' | 'depthInsensitivity'
%%% OUTPUT: vals .. [nCompModes x 1] value per mode: 'ritzScore' in dB
%%%                 (20 x Vis.ritzModeScores, NaN = the mode supports no Ritz
%%%                 candidate); the other two read straight off ritzSolver.diagnostics
%%%         info .. struct: label (axis / colour-bar label), vetoLevel (the freqDisp
%%%                 veto RitzScore.MIN_FREQDISP on a half-space; NaN otherwise --
%%%                 the veto is not used on layered stacks)
    names = {'ritzScore', 'freqDisp', 'depthInsensitivity'};
    if nargin == 0
        vals = names;
        return
    end
    ritz = singleResult.ritzSolver;
    nM   = ritz.nCompModes;
    info = struct('label', criterion, 'vetoLevel', NaN);
    switch criterion
        case 'ritzScore'
            vals = 20 * Vis.ritzModeScores(singleResult);
            info.label = 'Ritz score [dB]';
        case 'freqDisp'
            vals = ritz.diagnostics.freqDisp(1:nM);
            if ~visIsLayered(singleResult)
                info.vetoLevel = RitzScore.MIN_FREQDISP;
            end
        case 'depthInsensitivity'
            vals = ritz.diagnostics.depthInsensitivity(1:nM);
        otherwise
            %%% unreachable through the plot functions, which validate first
            error('Vis:visRitzCriterion:badCriterion', 'Unknown criterion ''%s''.', criterion);
    end
    vals = vals(:);
end
