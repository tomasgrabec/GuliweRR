function cutoffIdx = visDecayCutoffIdx(dispMag, thresholdFrac)
%%% First depth index (scanning inward from the surface) at which a
%%% displacement profile has decayed below a fraction of its surface value
%%% >> cutoffIdx = visDecayCutoffIdx(dispMag, thresholdFrac)
%%%
%%% @Vis-private helper shared by Vis.calcSubstrateCutoffDepth's 'decay'
%%% branch and Vis.plotGreenModalShape's own decay crop -- the same scan, on
%%% the same row convention, in both places.
%%%
%%% INPUT: 1. dispMag       .. [nMeshDepth x 1] displacement magnitude,
%%%                            index 1 = deepest, index end = the reference
%%%                            surface (a free surface, or the substrate/film
%%%                            interface)
%%%        2. thresholdFrac .. threshold as a fraction of the value AT that
%%%                            reference surface (e.g. 0.05)
%%% OUTPUT: cutoffIdx .. index to truncate at; 1 (keep everything) if the
%%%                      profile never drops below the threshold
%%%
%%% Scanning CONTIGUOUSLY from the surface inward, rather than looking for
%%% "anywhere below threshold", is what makes an oscillating (overtone) depth
%%% profile behave sensibly: the first crossing on the way down is the edge of
%%% the region that still carries the mode, not an incidental node.

    reference = max(dispMag(end), eps);   % value AT the reference surface
    threshold = thresholdFrac * reference;

    cutoffIdx = 1;   % default: never dropped below -> keep the whole profile
    for k = numel(dispMag):-1:1
        if dispMag(k) < threshold
            cutoffIdx = k;
            break
        end
    end
end
