function scores = ritzModeScores(singleResult)
%%% Per-mode Ritz score: each Ritz mode gets the score of the Ritz candidate it supports
%%% >> scores = Vis.ritzModeScores(singleResult)
%%%
%%% INPUT:  1. singleResult .. CalcSingleHS or CalcSingleLayered, calcAll()'d (either
%%%                           calcMode: the Ritz candidates are stored on both)
%%% OUTPUT: scores .. [nCompModes x 1] the Ritz rule's score (log10; x 20 for dB) of
%%%                   the candidate whose supporting mode (RitzScore.localMode) is
%%%                   this mode -- RitzScore.candidateScore (half-space dScore,
%%%                   layered dScore + polTerm), the value the rule thresholds. The highest one where a mode
%%%                   supports several candidates; NaN where it supports none.
%%%
%%% Every scored candidate counts, accepted or not (a vetoed or gated candidate still
%%% has its score), so the view shows where the rule looked, not only what it kept.
%%% The criterion 'ritzScore' of Vis.plotRitzModeMap / plotRitzModeCurve.
    nM = singleResult.ritzSolver.nCompModes;
    scores = NaN(nM, 1);
    cands = singleResult.ritzScoreCands;
    if isempty(cands); return; end
    sc = RitzScore.candidateScore(cands, visIsLayered(singleResult));   % what the rule thresholds
    for j = 1:numel(cands)
        m = cands(j).modeNo;
        if m < 1 || m > nM; continue; end
        if isnan(scores(m)) || sc(j) > scores(m)
            scores(m) = sc(j);
        end
    end
end
