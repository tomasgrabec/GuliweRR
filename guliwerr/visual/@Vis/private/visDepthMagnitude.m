function dispMag = visDepthMagnitude(fieldsPerDir)
%%% Per-depth-row displacement magnitude of a signed 2D field triplet
%%% >> dispMag = visDepthMagnitude(fieldsPerDir)
%%%
%%% @Vis-private helper: the depth profile a decay cutoff is measured against,
%%% built identically wherever it is needed (Vis.plotModalShape's 'decay'
%%% substrate cutoff, Vis.plotGreenModalShape's own decay crop).
%%%
%%% INPUT: 1. fieldsPerDir .. 3-element cell, each [nMeshDepth x nMeshLateral]
%%%                           signed real field for one displacement direction
%%%                           (as produced by Vis.calcDispField2D /
%%%                           Vis.calcGreenDispField2D)
%%% OUTPUT: dispMag .. [nMeshDepth x 1], row-wise magnitude across the three
%%%                    directions
%%%
%%% Each direction contributes its own lateral MAXIMUM (not an RMS or a
%%% mean): the lateral factor is a harmonic, so a row's peak amplitude is the
%%% physically meaningful measure of "how much this depth still moves",
%%% independent of where in the lateral period the sample happens to sit.

    rowMagSq = zeros(size(fieldsPerDir{1}, 1), 1);
    for iDir = 1:3
        rowMagSq = rowMagSq + max(fieldsPerDir{iDir}, [], 2).^2;
    end
    dispMag = sqrt(rowMagSq);
end
