function [hOut, hAx] = plotYoungsCurve(mater, angleDeg, hAx)
%%% Directional Young's modulus vs. in-plane angle, at the fixed surface normal
%%% in mater.dirMat
%%% >> [hOut, hAx] = plotYoungsCurve(mater, angleDeg, hAx)
%%%
%%% INPUT: 1. mater    .. Material object (elastic constants + density + dirMat
%%%                       orientation)
%%%  (opt) 2. angleDeg .. in-plane rotation angles [deg] (default 0:1:360, i.e.
%%%                       a full circle at 1 deg resolution)
%%%  (opt) 3. hAx       .. target axes, created fresh if omitted
%%% OUTPUT: hOut .. line handle
%%%         hAx  .. axes used
%%%
%%% The 2D-curve analogue of Vis.plotYoungsSurface (which sweeps the full
%%% sphere): this sweeps one great circle -- the in-plane rotation of
%%% mater.dirMat about its fixed surface normal (dirMat row 1), the same angle
%%% convention Vis.plotVelCurve/Christoffel.phaseVelScan use for bulk waves.
%%% mater.GetRotMat(angleDeg) already returns the 3x2 [propagation|normal]
%%% format ElastConsts.Tensor2Youngs expects directly.

    if nargin < 2 || isempty(angleDeg); angleDeg = linspace(0, 360, 361); end
    if nargin < 3 || isempty(hAx)
        figure();
        hAx = axes();
    end

    E = zeros(size(angleDeg));
    for ia = 1:numel(angleDeg)
        rotMat = mater.GetRotMat(angleDeg(ia));
        E(ia)  = mater.elConsts.GetYoungsModulus(rotMat);
    end

    hOut = plot(hAx, angleDeg, E);
    xlabel(hAx, 'angle [deg]');
    ylabel(hAx, 'Young''s modulus [GPa]');
end
