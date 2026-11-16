function [uCoord, wCoord, ptsXYZ] = extractPlanarCut(X, Y, Z, theta, phi, planeNormal, inPlaneRef)
%%% Extract the planar cut of a (theta,phi)-parametrised 3D surface
%%% >> [uCoord, wCoord, ptsXYZ] = extractPlanarCut(X, Y, Z, theta, phi, planeNormal, inPlaneRef)
%%%
%%% INPUT: 1-3. X,Y,Z     .. [nPhi x nTheta] surface coordinates for ONE
%%%                          sheet/wave-type (call once per sheet for a
%%%                          multi-sheet surface, e.g. plotVelSurface's or
%%%                          plotGroupVelSurface's 3 wave types -- typically
%%%                          hOut(k).XData/.YData/.ZData directly)
%%%        4-5. theta,phi .. the parameter grids used to build X,Y,Z (plain
%%%                          vectors, e.g. linspace(0,pi,nTheta)/
%%%                          linspace(0,2*pi,nPhi) -- the SAME nTheta/nPhi
%%%                          passed to plotVelSurface/plotGroupVelSurface)
%%%        6. planeNormal .. 3-vector, normal of the cutting plane (through
%%%                          the origin) -- e.g. [0;0;1] for the sagittal
%%%                          plane in the standard convention (axis 3 =
%%%                          surface normal) used throughout this project,
%%%                          including plotVelSurface/plotGroupVelSurface.
%%%  (opt) 7. inPlaneRef  .. 3-vector, one of the two in-plane basis
%%%                          directions to project onto (default: an
%%%                          arbitrary vector perpendicular to planeNormal)
%%% OUTPUT: uCoord, wCoord .. 2D in-plane coordinates of the cut curve --
%%%                           possibly several disconnected segments,
%%%                           NaN-separated (matches MATLAB's own plot/
%%%                           contour convention for broken curves; SAW/
%%%                           pseudoSAW cuts in particular need not be
%%%                           continuous)
%%%         ptsXYZ         .. [N x 3] the same points in original 3D
%%%                           coordinates, same NaN-row separation
%%%
%%% The full-sphere sweep is needed first because a vector-valued surface's
%%% plotted point is not simply the wave normal scaled by a scalar — group
%%% velocity can point off-plane.
%%%
%%% Method: the surface point at (theta,phi) lies in the cutting plane
%%% exactly where f(theta,phi) = dot([X,Y,Z](theta,phi), planeNormal) = 0.
%%% Finds that zero contour with MATLAB's own contourc over the (theta,phi)
%%% parameter domain, then bilinearly interpolates X,Y,Z at each contour
%%% vertex to recover the exact 3D crossing point and projects onto two
%%% orthonormal in-plane basis vectors.

    n = planeNormal(:) / norm(planeNormal(:));
    if nargin < 7 || isempty(inPlaneRef)
        %%% pick whichever standard-basis axis is LEAST aligned with n, to
        %%% avoid a near-degenerate (n, candidate) pair
        candidates = eye(3);
        [~, idx] = min(abs(candidates' * n));
        v = candidates(:, idx);
    else
        v = inPlaneRef(:);
    end
    u = v - (v' * n) * n;
    u = u / norm(u);
    w = cross(n, u);

    F = X .* n(1) + Y .* n(2) + Z .* n(3);   % signed distance from the plane

    if all(F(:) > 0) || all(F(:) < 0)
        %%% surface never crosses the plane -- contourc would error on a
        %%% degenerate level; return empty rather than let that happen
        uCoord = []; wCoord = []; ptsXYZ = zeros(0, 3);
        return
    end

    C = contourc(theta, phi, F, [0 0]);

    uCoord = []; wCoord = []; ptsXYZ = zeros(0, 3);
    idx = 1;
    while idx < size(C, 2)
        nPts = C(2, idx);
        seg  = C(:, idx+1 : idx+nPts);   % row 1 = theta, row 2 = phi
        segTheta = seg(1, :);
        segPhi   = seg(2, :);

        segX = interp2(theta, phi, X, segTheta, segPhi);
        segY = interp2(theta, phi, Y, segTheta, segPhi);
        segZ = interp2(theta, phi, Z, segTheta, segPhi);

        segU = segX * u(1) + segY * u(2) + segZ * u(3);
        segW = segX * w(1) + segY * w(2) + segZ * w(3);

        uCoord = [uCoord, segU, NaN];                                       %#ok<AGROW>
        wCoord = [wCoord, segW, NaN];                                       %#ok<AGROW>
        ptsXYZ = [ptsXYZ; [segX(:), segY(:), segZ(:)]; [NaN, NaN, NaN]];     %#ok<AGROW>

        idx = idx + nPts + 1;
    end
end
