function hCb = visGreenColourAxis(hAx, scaleMode, climRange, cbLabel)
%%% Apply the shared Green's-magnitude colour treatment to an axes
%%% >> hCb = visGreenColourAxis(hAx, scaleMode, climRange, cbLabel)
%%%
%%% @Vis-private helper shared by plotGreenMap, plotDispersionMap's heatmap cases
%%% and plotDispersionCurves' Green background, so every Green |G| map in visual/
%%% reads on the same colour scale, and a function's Cartesian and polar branches
%%% get identical colour handling — Vis.plotPolarHeatmap does no colour work of
%%% its own.
%%%
%%% INPUT: 1. hAx       .. target axes
%%%        2. scaleMode .. 'log' | 'linear'. clim is applied for 'log' only: a
%%%                        'linear' map's numeric range is material- and
%%%                        units-dependent, so a fixed limit pair would be
%%%                        meaningless and MATLAB's auto-scaling is right
%%%        3. climRange .. 2-element colour-axis limits, log10(normalised) units;
%%%                        the caller's own default, which differs legitimately
%%%                        between scan families
%%%        4. cbLabel   .. colorbar label (Vis.scaleGreenData's second output)
%%% OUTPUT: hCb .. colorbar handle
%%%
%%% ColMapGrayOrange: |G| is a magnitude, so a sequential map with no meaningful
%%% zero to centre on is the right family.

    colormap(hAx, Vis.ColMapGrayOrange());
    if strcmp(scaleMode, 'log')
        clim(hAx, climRange);
    end
    hCb = colorbar(hAx);
    ylabel(hCb, cbLabel);
end
