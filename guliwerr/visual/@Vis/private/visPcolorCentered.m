function hOut = visPcolorCentered(hAx, x, y, C)
%%% pcolor a map with each sample drawn CENTRED on its own coordinate
%%% >> hOut = visPcolorCentered(hAx, x, y, C)
%%%
%%% @Vis-private helper: the three-line "edges + pcolor + shading flat" idiom
%%% every Cartesian heatmap in visual/ uses (plotGreenMap, plotDispersionMap,
%%% plotDispersionCurves' Green background), written once. The edge
%%% construction itself is the public, separately-tested
%%% Vis.pcolorCenteredGrid -- this only bundles it with the pcolor call that
%%% must accompany it, so no caller can accidentally pcolor the raw sample
%%% coordinates instead (plain pcolor(x,y,C) puts each sample at a cell's
%%% lower-left corner, offsetting the whole map from its own axis ticks --
%%% docs/pitfalls.md).
%%%
%%% INPUT: 1. hAx .. target axes
%%%        2. x   .. [1xN] or [Nx1] sample coordinates (monotonic)
%%%        3. y   .. [1xM] or [Mx1] sample coordinates (monotonic)
%%%        4. C   .. [MxN] colour data
%%% OUTPUT: hOut .. pcolor handle

    [xEdges, yEdges, CPadded] = Vis.pcolorCenteredGrid(x, y, C);
    hOut = pcolor(hAx, xEdges, yEdges, CPadded);
    shading(hAx, 'flat');
end
