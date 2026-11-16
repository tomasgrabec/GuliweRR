function [hBar, hAx, hHeat] = plotSensitivitySpectrum(sensitivityAnalysis, hFig)
%%% Plot the sensitivity spectrum: normalised sqrt(eigenvalue) bars (the
%%% error-amplification factor of each transformed combination relative to
%%% the best-determined one) plus the SIGNED transformation matrix as a
%%% labelled heat-map, in the form of Grabec et al. 2024's Tables A.1-A.4.
%%% >> [hBar, hAx, hHeat] = Vis.plotSensitivitySpectrum(sensitivityAnalysis, hFig)
%%%
%%% INPUT:  1. sensitivityAnalysis .. a computed SensitivityAnalysis object
%%%  (opt)  2. hFig                .. figure handle (new figure if omitted)
%%% OUTPUT: hBar  .. bar-chart handle (normalisedSensitivityRoots, log scale)
%%%         hAx   .. [1x2] axes handles [barAxes, heatAxes]
%%%         hHeat .. heat-map image handle (SIGNED transformMatrix)
%%%
%%% Bars at/above sensitivityCutoff (the optimised combinations) are shown in
%%% the default colour; bars below it (frozen) are greyed out.
%%%
%%% The heat-map plots transformMatrix SIGNED, not abs(): the sign is what
%%% tells, e.g. for cubic symmetry, a longitudinal-like "c11+c12" combination
%%% (same sign on both) from a shear-like "c11-c12" one (opposite signs).
%%% Diverging Vis.ColormapRWB, symmetric about zero so white always means zero,
%%% as for other signed fields (e.g. displacement modal shapes).

    if nargin < 2 || isempty(hFig); hFig = figure(); end

    L = sensitivityAnalysis.nIndepConsts;
    roots = sensitivityAnalysis.normalisedSensitivityRoots;
    P = sensitivityAnalysis.nOptimisedConstants;
    T = sensitivityAnalysis.transformMatrix;

    %%% "c11"/"c12"/"c13"/... labels for each independent constant, from its
    %%% [row,col] position in the 6x6 Voigt matrix (1-indexed) -- e.g. cubic
    %%% (symmClassNo=2) has exactly 3: c11, c12, c44.
    posList = ElastConsts.IndepVoigtPositions(sensitivityAnalysis.symmClassNo);
    constLabels = cellfun(@(rc) sprintf('c%d%d', rc(1), rc(2)), posList, 'UniformOutput', false);

    hAx = gobjects(1, 2);

    hAx(1) = subplot(1, 2, 1, 'Parent', hFig);
    barColors = repmat([0 0.4470 0.7410], L, 1);
    if P < L
        barColors(P+1:end, :) = repmat([0.7 0.7 0.7], L-P, 1);
    end
    hBar = bar(hAx(1), 1:L, roots, 'FaceColor', 'flat');
    hBar.CData = barColors;
    hold(hAx(1), 'on');
    yline(hAx(1), sensitivityAnalysis.sensitivityCutoff, '--r', 'cutoff');
    set(hAx(1), 'YScale', 'log');
    xlabel(hAx(1), 'transformed combination p');
    ylabel(hAx(1), '(\surd\psi_p)_{norm} = \surd(\psi_p/\psi_1)');
    title(hAx(1), sprintf('Sensitivity spectrum (%d of %d optimised)', P, L));
    xticks(hAx(1), 1:L);

    hAx(2) = subplot(1, 2, 2, 'Parent', hFig);
    hHeat = imagesc(hAx(2), T);
    colorbar(hAx(2));
    colormap(hAx(2), Vis.ColormapRWB());
    maxAbsT = max(abs(T(:)));
    clim(hAx(2), [-maxAbsT, maxAbsT]);   % symmetric about zero -- white = 0, not the data min
    xlabel(hAx(2), 'transformed combination p');
    ylabel(hAx(2), 'independent Voigt constant');
    title(hAx(2), 'T -- signed composition of each combination');
    xticks(hAx(2), 1:L);
    yticks(hAx(2), 1:L);
    yticklabels(hAx(2), constLabels);
    axis(hAx(2), 'square');

    %%% print each cell's value (legible at L<=~10) -- two entries that both
    %%% render as "strong blue" can still differ in magnitude
    for iRow = 1:L
        for jCol = 1:L
            textColor = 'k';
            if abs(T(iRow, jCol))/maxAbsT > 0.6   % dark background -> white text
                textColor = 'w';
            end
            text(hAx(2), jCol, iRow, sprintf('%.2f', T(iRow, jCol)), ...
                'HorizontalAlignment', 'center', 'Color', textColor, 'FontSize', 8);
        end
    end
end
