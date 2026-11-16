function dispField = visDispField2D(depthProfiles, meshSpacing)
%%% Shared depth-profile -> signed 2D (depth x lateral) field assembly
%%% >> dispField = visDispField2D(depthProfiles, meshSpacing)
%%%
%%% @Vis-private helper behind BOTH Vis.calcDispField2D (Ritz eigenmode, nLay
%%% profiles) and Vis.calcGreenDispField2D (Green's response, a single half-space
%%% profile): they differ only in where the complex depth profile comes from, not
%%% in how it becomes a field.
%%%
%%% INPUT: 1. depthProfiles .. {1 x nLay} cell, each [nMeshDepth x 3] COMPLEX
%%%                            depth profile for one layer, row 1 = deepest
%%%                            (ti=-1), row end = shallowest (ti=+1)
%%%        2. meshSpacing   .. lateral mesh spacing (already defaulted by the caller)
%%% OUTPUT: dispField .. {3 x nLay} cell, dispField{iDir,iLay} = real signed 2D
%%%                      field [nMeshDepth x nMeshLateral], normalised by one
%%%                      global max across every direction and layer together
%%%
%%% That max is over POSITIVE values only, not max(abs(...)), so a
%%% negative-dominant mode can fall outside [-1,1], which the callers' fixed
%%% clim([-1,1]) clips.
%%%
%%% lateralFactor must stay genuinely COMPLEX through the multiplication: the
%%% depth profile is complex (Hermitian eigenproblem), so simplifying to a real
%%% cosine would silently drop its imaginary part.

    ti            = -1 : meshSpacing : 1;
    lateralFactor = exp(1i * pi * (ti + 1));   % [1 x nMeshLateral], complex

    nLay      = numel(depthProfiles);
    dispField = cell(3, nLay);
    maxPos    = 0;
    for ilay = 1:nLay
        for iDir = 1:3
            field2D = real(depthProfiles{ilay}(:, iDir) * lateralFactor);   % outer product, real part
            dispField{iDir, ilay} = field2D;
            maxPos = max([maxPos, max(field2D(:))]);
        end
    end

    maxPos = max(maxPos, eps);
    for k = 1:numel(dispField)
        dispField{k} = dispField{k} ./ maxPos;
    end
end
