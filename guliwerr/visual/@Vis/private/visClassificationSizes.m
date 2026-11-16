function sizeVals = visClassificationSizes(pts, isLayered, minSz, maxSz, defaultSz, limitingSz)
%%% Marker size of each classified point by its score, clipped to its method's range
%%% >> sizeVals = visClassificationSizes(pts, isLayered, minSz, maxSz, defaultSz, limitingSz)
%%%
%%% @Vis-private helper shared by Vis.plotClassificationDiagnostic and
%%% Vis.plotClassificationCurve.
%%%
%%% INPUT: 1. pts        .. visFlattenSignificantModes output (score, method, type)
%%%        2. isLayered  .. true for a layered result (selects the method's scale)
%%%        3-6. minSz, maxSz, defaultSz, limitingSz .. marker areas: the score range
%%%             maps linearly onto [minSz, maxSz]; defaultSz where a point has no
%%%             score; limitingSz for every 'limiting' point
%%% OUTPUT: sizeVals .. [nPoints x 1]
%%%
%%% The range is Vis.classificationScoreScale(method, isLayered).sizeRange: from the
%%% method's looser ('high') threshold, where the weakest accepted entry sits, to a
%%% decade or more above its stricter one. A score is a log10 amplitude ratio, so the
%%% size grows with the response level in dB.
    sizeVals = defaultSz * ones(numel(pts.vel), 1);
    meths = unique(pts.method(~cellfun(@isempty, pts.method)));
    for k = 1:numel(meths)
        sRange = Vis.classificationScoreScale(meths{k}, isLayered).sizeRange;
        mask = strcmp(pts.method, meths{k})' & ~isnan(pts.score);
        frac = min(max((pts.score(mask) - sRange(1)) / (sRange(2) - sRange(1)), 0), 1);
        sizeVals(mask) = minSz + (maxSz - minSz) * frac;
    end
    sizeVals(strcmp(pts.type, 'limiting')) = limitingSz;
end
