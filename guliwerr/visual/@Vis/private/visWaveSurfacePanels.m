function [hOut, hAx] = visWaveSurfacePanels(hFig, X, Y, Z, C)
%%% Render the qL/qT1/qT2 row of 3D wave-surface panels
%%% >> [hOut, hAx] = visWaveSurfacePanels(hFig, X, Y, Z, C)
%%%
%%% @Vis-private helper: the shared rendering half of Vis.plotVelSurface and
%%% Vis.plotGroupVelSurface. Only the surface GEOMETRY differs between those
%%% two (a scalar phase velocity/slowness scaling the wave normal vs. the
%%% genuine group-velocity vector, which need not point along it) -- the panel
%%% layout, labels, titles and colour treatment are identical, so they live
%%% here once.
%%%
%%% INPUT: 1. hFig    .. target figure
%%%        2-4. X,Y,Z .. [nPhi x nTheta x 3] surface coordinates, third
%%%                      dimension = wave type (qL, qT1, qT2 in that order)
%%%        5. C       .. [nPhi x nTheta x 3] colour data, same indexing
%%% OUTPUT: hOut .. 1x3 array of surf handles
%%%         hAx  .. 1x3 array of axes handles
%%%
%%% Axis labels name the standard frame's axes (x_3 = surface normal) -- see
%%% visStandardFrameTensor's docblock for the traced derivation.

    titles = {'qL', 'qT1', 'qT2'};
    hOut = gobjects(1, 3);
    hAx  = gobjects(1, 3);
    for k = 1:3
        hAx(k)  = subplot(1, 3, k, 'Parent', hFig);
        hOut(k) = surf(hAx(k), X(:, :, k), Y(:, :, k), Z(:, :, k), C(:, :, k), ...
            'FaceColor', 'interp', 'EdgeColor', 'none');
        axis(hAx(k), 'equal');
        xlabel(hAx(k), 'x_1 (in-plane ref)');
        ylabel(hAx(k), 'x_2 (in-plane)');
        zlabel(hAx(k), 'x_3 (surface normal)');
        title(hAx(k), titles{k});
        colorbar(hAx(k));
    end
end
