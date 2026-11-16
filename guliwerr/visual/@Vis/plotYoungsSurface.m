function [hOut, hAx] = plotYoungsSurface(mater, nTheta, nPhi, hFig)
%%% 3D directional Young's-modulus surface
%%% >> [hOut, hAx] = plotYoungsSurface(mater, nTheta, nPhi, hFig)
%%%
%%% INPUT: 1. mater  .. Material object (elastic constants + density + dirMat
%%%                     orientation)
%%%  (opt) 2. nTheta .. polar-angle mesh resolution (default 61)
%%%  (opt) 3. nPhi   .. azimuthal mesh resolution (default 121)
%%%  (opt) 4. hFig   .. target figure (created fresh if omitted)
%%% OUTPUT: hOut .. surf handle
%%%         hAx  .. axes used
%%%
%%% Single-panel analogue of Vis.plotVelSurface: Young's modulus is one scalar
%%% field (not three wave-type sheets), so this does not route through
%%% visWaveSurfacePanels -- the plotted point is the unit direction scaled by
%%% E(direction), same "radius = scalar modulus" convention as plotVelSurface's
%%% velocity sheets.
%%%
%%% The tensor is taken in the standard frame used throughout this project
%%% (Ritz/Christoffel/GreenSolver): axis 1 = dirMat's in-plane propagation
%%% reference (dirMat row 2), axis 2 = the remaining in-plane direction (dirMat
%%% row 3), axis 3 = the SURFACE NORMAL (dirMat row 1).

    if nargin < 2 || isempty(nTheta); nTheta = 61;  end
    if nargin < 3 || isempty(nPhi);   nPhi   = 121; end
    if nargin < 4 || isempty(hFig);   hFig   = figure(); end

    [elTensor, ~] = visStandardFrameTensor(mater);
    [~, ~, X, Y, Z] = visSphereDirGrid(nTheta, nPhi);

    E = zeros(nPhi, nTheta);
    for it = 1:nTheta
        for ip = 1:nPhi
            n1 = [X(ip, it), Y(ip, it), Z(ip, it)];

            %%% cheap per-point orthonormal completion (n1,n2,n3) -- avoids
            %%% ElastConsts.RightHandFrame's warning side effect, which would
            %%% otherwise fire on nearly every one of the nTheta*nPhi points
            ref = [0, 0, 1];
            if abs(dot(n1, ref)) > 0.9   % n1 too close to ref -- switch axis
                ref = [1, 0, 0];
            end
            n2 = cross(ref, n1); n2 = n2 / norm(n2);
            n3 = cross(n1, n2);

            E(ip, it) = ElastConsts.Tensor2Youngs(elTensor, [n1; n2; n3]);
        end
    end

    hAx = axes('Parent', hFig);
    hOut = surf(hAx, X .* E, Y .* E, Z .* E, E, 'FaceColor', 'interp', 'EdgeColor', 'none');
    axis(hAx, 'equal');
    xlabel(hAx, 'x_1 (in-plane ref)');
    ylabel(hAx, 'x_2 (in-plane)');
    zlabel(hAx, 'x_3 (surface normal)');
    title(hAx, 'Young''s modulus [GPa]');
    colorbar(hAx);
end
