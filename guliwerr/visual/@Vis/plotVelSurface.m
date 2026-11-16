function [hOut, hAx] = plotVelSurface(mater, surfaceType, nTheta, nPhi, hFig)
%%% 3D directional bulk-wave velocity (or slowness) surfaces
%%% >> [hOut, hAx] = plotVelSurface(mater, surfaceType, nTheta, nPhi, hFig)
%%%
%%% INPUT: 1. mater       .. Material object (elastic constants + density +
%%%                          dirMat orientation) -- a pure bulk-wave/Christoffel
%%%                          quantity, no Ritz/Green/ModeClassification
%%%                          dependency, so this takes a Material directly
%%%                          rather than a Structure/CalcSingle* wrapper
%%%  (opt) 2. surfaceType .. 'velocity' (default) | 'slowness' (= 1/velocity,
%%%                          computed slot-wise, i.e. per wave-type branch,
%%%                          not re-sorted after inversion)
%%%  (opt) 3. nTheta      .. polar-angle mesh resolution (default 61)
%%%  (opt) 4. nPhi        .. azimuthal mesh resolution (default 121)
%%%  (opt) 5. hFig        .. target figure (created fresh if omitted)
%%% OUTPUT: hOut .. 1x3 array of surf handles [qL, qT1, qT2]
%%%         hAx  .. 1x3 array of axes handles
%%%
%%% The tensor is taken in the standard frame used throughout this project
%%% (Ritz/Christoffel/GreenSolver): axis 1 = dirMat's in-plane propagation
%%% reference (dirMat row 2), axis 2 = the remaining in-plane direction (dirMat
%%% row 3), axis 3 = the SURFACE NORMAL (dirMat row 1). A full spherical sweep
%%% covers the whole sphere in whatever frame elTensor is expressed in, but the
%%% axis LABELS below name them correctly so the plot reads physically.
%%%
%%% Group-velocity surfaces are Vis.plotGroupVelSurface's job: there the plotted
%%% point is a genuine vector that need not lie along the wave normal.

    if nargin < 2 || isempty(surfaceType); surfaceType = 'velocity'; end
    visValidateChoice(surfaceType, {'velocity', 'slowness'}, 'plotVelSurface', 'surfaceType');
    if nargin < 3 || isempty(nTheta); nTheta = 61;  end
    if nargin < 4 || isempty(nPhi);   nPhi   = 121; end
    if nargin < 5 || isempty(hFig)
        hFig = visWaveSurfaceFigure();   % wide enough for a 1x3 'axis equal' row
    end

    [elTensor, rho] = visStandardFrameTensor(mater);
    [~, ~, X, Y, Z] = visSphereDirGrid(nTheta, nPhi);

    V = zeros(nPhi, nTheta, 3);
    for it = 1:nTheta
        for ip = 1:nPhi
            dirVec = [X(ip, it); Y(ip, it); Z(ip, it)];
            V(ip, it, :) = Christoffel.phaseVels(elTensor, dirVec, rho);
        end
    end

    if strcmp(surfaceType, 'slowness')
        V = 1 ./ V;
    end

    %%% each sheet is the unit wave normal scaled by that wave type's own
    %%% scalar radius -- unlike plotGroupVelSurface, where the plotted point
    %%% is a genuine vector that need not lie along the wave normal
    [hOut, hAx] = visWaveSurfacePanels(hFig, X .* V, Y .* V, Z .* V, V);
end
