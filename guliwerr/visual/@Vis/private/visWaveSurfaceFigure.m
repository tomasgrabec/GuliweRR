function hFig = visWaveSurfaceFigure()
%%% Fresh figure sized for a 1x3 row of 3D wave-surface panels
%%% >> hFig = visWaveSurfaceFigure()
%%%
%%% @Vis-private helper shared by Vis.plotVelSurface and
%%% Vis.plotGroupVelSurface.
%%%
%%% A default figure() is ~4:3 landscape; split 1x3, each panel cell is narrower
%%% than tall, so a 3D 'axis equal' surf shrinks to fit that width while the
%%% colorbar is sized to the cell's full height — a tiny plot beside an oversized
%%% colorbar. Widening the figure per-panel fixes it. Used only when creating a
%%% fresh figure, never to override a caller-supplied hFig's size.

    defPos = get(0, 'defaultfigureposition');
    hFig   = figure('Position', [defPos(1), defPos(2), 3 * 480, 450]);
end
