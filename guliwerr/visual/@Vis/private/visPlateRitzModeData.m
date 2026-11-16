function [vels, types, conf] = visPlateRitzModeData(ritz)
%%% Per-mode velocity, mode type and type confidence from a Lamb Ritz solver,
%%% all as COLUMNS
%%% >> [vels, types, conf] = visPlateRitzModeData(ritz)
%%%
%%% @Vis-private helper shared by Vis.plotRitzModeMapPlate and
%%% Vis.plotRitzModeCurvePlate.
%%%
%%% INPUT: 1. ritz .. RitzSolverLamb (vels/nCompModes/modesType/
%%%                   modeTypeConfidence read directly -- works at any
%%%                   retentionLevel, since pruneEigVecs only ever touches
%%%                   eigVecs_con, never these fields)
%%% OUTPUT: vels, types, conf .. [nCompModes x 1] each
%%%
%%% The reshape is load-bearing, not cosmetic: modesType/modeTypeConfidence
%%% are stored as ROW vectors (RitzSolverLamb's own zeros(1, nCompModes)
%%% convention) while vels is a COLUMN. Logical-indexing those rows with a
%%% column keep-mask silently preserves row orientation, which corrupted
%%% plotRitzModeMapPlate's vertical [allType; ...] accumulation into an
%%% [nAngles x nM] stack instead of a flat column -- caught then by
%%% scatter's own 'C must be m-by-1' error rather than by a wrong answer.
%%% Forcing columns here fixes that at the single point where the data is
%%% read, for every caller at once.

    nM    = ritz.nCompModes;
    vels  = reshape(ritz.vels(1:nM), [], 1);
    types = reshape(ritz.modesType(1:nM), [], 1);
    conf  = reshape(ritz.modeTypeConfidence(1:nM), [], 1);
end
