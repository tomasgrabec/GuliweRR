function dispField = calcDispField2D(ritzSolver, modeNo, meshSpacing)
%%% Signed 2D (depth x lateral) displacement field for one Ritz mode
%%% >> dispField = Vis.calcDispField2D(ritzSolver, modeNo, meshSpacing)
%%%
%%% INPUT:  1. ritzSolver  .. RitzSolver (or subclass) object, eigVecs_con/
%%%                           nullConstrMat populated for the requested mode
%%%                           (retentionLevel='full'; there is no
%%%                           recompute-on-demand for a pruned RitzSolver)
%%%         2. modeNo      .. which mode (column of eigVecs_con) to evaluate
%%%  (opt)  3. meshSpacing .. depth/lateral mesh spacing (default 0.01, matching
%%%                           RitzSolver.calcDisplacement's cached-basis default)
%%% OUTPUT: dispField .. {3 x nLay} cell array, dispField{iDir,iLay} = real signed
%%%                      2D field [nMeshDepth x nMeshLateral], normalised by one
%%%                      global max across every direction and layer together.
%%%                      That max is over POSITIVE values only, not max(abs(...)),
%%%                      so a negative-dominant mode can fall outside [-1,1].
%%%
%%% Built on RitzSolver.calcDisplacement; this function only slices that output
%%% into per-layer depth profiles. The lateral factor, outer product and
%%% normalisation are the @Vis-private visDispField2D helper.

    if nargin < 3 || isempty(meshSpacing)
        meshSpacing = 0.01;
    end

    displDir     = ritzSolver.calcDisplacement(modeNo, meshSpacing);   % [depthPts*nLay x modeNo x 3]
    depthProfile = squeeze(displDir(:, end, :));                       % mode modeNo only, [depthPts*nLay x 3]

    nLay      = ritzSolver.structure.nLay;
    displSize = round(2 / meshSpacing) + 1;

    depthProfiles = cell(1, nLay);
    for ilay = 1:nLay
        depthProfiles{ilay} = depthProfile((ilay-1)*displSize + 1 : ilay*displSize, :);
    end

    dispField = visDispField2D(depthProfiles, meshSpacing);
end
