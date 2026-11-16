function [hOut, hAx] = plotGroupVelSurface(mater, nTheta, nPhi, hFig)
%%% 3D group-velocity ("ray") surface -- the direction of energy/signal
%%% propagation, which in general deviates from the wave-normal direction
%%% >> [hOut, hAx] = plotGroupVelSurface(mater, nTheta, nPhi, hFig)
%%%
%%% INPUT: 1. mater  .. Material object (elastic constants + density +
%%%                     dirMat orientation) -- same convention as
%%%                     plotVelSurface
%%%  (opt) 2. nTheta .. polar-angle mesh resolution (default 61)
%%%  (opt) 3. nPhi   .. azimuthal mesh resolution (default 121)
%%%  (opt) 4. hFig   .. target figure (created fresh if omitted)
%%% OUTPUT: hOut .. 1x3 array of surf handles [qL, qT1, qT2]
%%%         hAx  .. 1x3 array of axes handles
%%%
%%% Same (theta,phi) full-sphere sweep and frame convention as plotVelSurface
%%% (axis 3 = surface normal). Unlike plotVelSurface, the surface point at
%%% (theta,phi) is NOT dirVec scaled by a scalar: it is the actual group-velocity
%%% VECTOR Christoffel.groupVelVectorsForDir returns, which can point in a
%%% genuinely different direction from the wave normal that produced it. A
%%% self-intersecting, "multiply-folded" surface is therefore real physics, not a
%%% plotting bug (Every et al. 2013, Wave Motion 50, Sec. 2.3).
%%%
%%% For a 2D planar cut (e.g. the sagittal plane), pass this function's returned
%%% surf handles' XData/YData/ZData to extractPlanarCut.m. Computing group
%%% velocity over the WHOLE sphere and then filtering is necessary here, unlike
%%% plotVelSurface, where restricting the sweep to one plane from the start
%%% already gives the right answer.

    if nargin < 2 || isempty(nTheta); nTheta = 61;  end
    if nargin < 3 || isempty(nPhi);   nPhi   = 121; end
    if nargin < 4 || isempty(hFig)
        hFig = visWaveSurfaceFigure();   % wide enough for a 1x3 'axis equal' row
    end

    [elTensor, rho] = visStandardFrameTensor(mater);
    [~, ~, Nx, Ny, Nz] = visSphereDirGrid(nTheta, nPhi);

    X = zeros(nPhi, nTheta, 3);
    Y = zeros(nPhi, nTheta, 3);
    Z = zeros(nPhi, nTheta, 3);
    for it = 1:nTheta
        for ip = 1:nPhi
            dirVec = [Nx(ip, it); Ny(ip, it); Nz(ip, it)];
            vG = Christoffel.groupVelVectorsForDir(elTensor, dirVec, rho);
            X(ip, it, :) = vG(1, :);
            Y(ip, it, :) = vG(2, :);
            Z(ip, it, :) = vG(3, :);
        end
    end

    %%% colour by group-velocity MAGNITUDE -- the surface radius is not itself
    %%% that magnitude here (the point need not lie along the wave normal)
    [hOut, hAx] = visWaveSurfacePanels(hFig, X, Y, Z, sqrt(X.^2 + Y.^2 + Z.^2));
end
