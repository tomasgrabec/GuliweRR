function [hOut, hAx] = plotVelCurve(mater, angleDeg, hAx)
%%% Bulk-wave phase velocity vs. in-plane angle, at the fixed surface normal in mater.dirMat
%%% >> [hOut, hAx] = plotVelCurve(mater, angleDeg, hAx)
%%%
%%% INPUT: 1. mater    .. Material object (elastic constants + density + dirMat
%%%                       orientation) -- a pure bulk-wave/Christoffel quantity,
%%%                       no Ritz/Green/ModeClassification dependency
%%%  (opt) 2. angleDeg .. in-plane rotation angles [deg] (default 0:1:360, i.e.
%%%                       a full circle at 1 deg resolution)
%%%  (opt) 3. hAx       .. target axes, created fresh if omitted
%%% OUTPUT: hOut .. 1x3 array of line handles [qL, qT1, qT2]
%%%         hAx  .. axes used
%%%
%%% The 2D-curve analogue of Vis.plotVelSurface: that sweeps the full sphere,
%%% this sweeps one great circle -- the in-plane rotation of mater.dirMat about
%%% its fixed surface normal (dirMat row 1), the same angle convention every
%%% other angular sweep in this project uses (Christoffel.phaseVelScan,
%%% CalcAngSweepHS/Layered, Vis.plotAngularMap).

    if nargin < 2 || isempty(angleDeg); angleDeg = linspace(0, 360, 361); end
    if nargin < 3 || isempty(hAx)
        figure();
        hAx = axes();
    end

    vels = Christoffel.phaseVelScan(mater, angleDeg);   % [vL;vT1;vT2] x nAngles

    hold(hAx, 'on');
    labels = {'qL', 'qT1', 'qT2'};
    hOut = gobjects(1, 3);
    for k = 1:3
        hOut(k) = plot(hAx, angleDeg, vels(k, :), 'DisplayName', labels{k});
    end
    xlabel(hAx, 'angle [deg]');
    ylabel(hAx, 'phase velocity [km/s]');
    legend(hAx);
end
