function tf = visHasRitzSpectrum(result)
%%% Whether a result stores the Ritz modal numerators |S(v)| is rebuilt from
%%% >> tf = visHasRitzSpectrum(result)
%%%
%%% @Vis-private helper: the one "has a Ritz spectrum" test of the Ritz spectrum and
%%% classification plots.
%%%
%%% INPUT:  1. result .. a single-angle result (CalcSingleHS / CalcSingleLayered) or
%%%                      an angular sweep (CalcAngSweepHS / CalcAngSweepLayered)
%%% OUTPUT: tf .. single angle: ritzScoreNum and ritzScoreVelAxis both non-empty;
%%%               sweep: true for every angle (and at least one angle)
    if isprop(result, 'cs')
        tf = ~isempty(result.cs) && all(cellfun(@hasOwn, result.cs));
    else
        tf = hasOwn(result);
    end
end

function tf = hasOwn(single)
    tf = ~isempty(single.ritzScoreNum) && ~isempty(single.ritzScoreVelAxis);
end
