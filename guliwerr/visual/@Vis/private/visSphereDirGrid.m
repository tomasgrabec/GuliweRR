function [theta, phi, Nx, Ny, Nz] = visSphereDirGrid(nTheta, nPhi)
%%% Full-sphere (theta,phi) direction grid and its unit direction cosines
%%% >> [theta, phi, Nx, Ny, Nz] = visSphereDirGrid(nTheta, nPhi)
%%%
%%% @Vis-private helper shared by Vis.plotVelSurface and
%%% Vis.plotGroupVelSurface, so both sweep exactly the same directions -- a
%%% precondition for Vis.extractPlanarCut, which takes the SAME theta/phi
%%% vectors back and would silently mis-locate crossings against a grid built
%%% even slightly differently.
%%%
%%% INPUT: 1. nTheta .. polar-angle resolution
%%%        2. nPhi   .. azimuthal resolution
%%% OUTPUT: theta, phi   .. the plain parameter vectors (hand these to
%%%                         Vis.extractPlanarCut unchanged)
%%%         Nx, Ny, Nz   .. [nPhi x nTheta] unit direction cosines

    theta = linspace(0, pi, nTheta);
    phi   = linspace(0, 2*pi, nPhi);
    [Theta, Phi] = meshgrid(theta, phi);   % both [nPhi x nTheta]

    Nx = sin(Theta) .* cos(Phi);
    Ny = sin(Theta) .* sin(Phi);
    Nz = cos(Theta);
end
