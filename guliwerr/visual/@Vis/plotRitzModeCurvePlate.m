function [hOut, hAx] = plotRitzModeCurvePlate(singleResult, criterion, velCutoffFactor, hAx)
%%% Plot every computed Ritz mode for a CalcSinglePlate result as a
%%% velocity x criterion scatter -- Plate analogue of plotRitzModeCurve.m
%%% (1D, single propagation angle, no angle axis to sweep).
%%% >> [hOut, hAx] = Vis.plotRitzModeCurvePlate(singleResult, criterion, velCutoffFactor, hAx)
%%%
%%% Separate from plotRitzModeCurve rather than a class-dispatch branch inside
%%% it: the only Plate mode-quality concept is mode-TYPE classification
%%% (RitzSolverLamb.modesType/modeTypeConfidence), SAW's
%%% diagnostics.depthInsensitivity/freqDisp having no Plate analogue, and there is
%%% no Green's-function curve option at all.
%%%
%%% INPUT: 1. singleResult     .. CalcSinglePlate result object, calc()'d. Reads
%%%                              .ritzSolver directly (vels/nCompModes/
%%%                              modesType/modeTypeConfidence), so it works at
%%%                              any retentionLevel
%%%  (opt) 2. criterion        .. 'modeType' (default) | 'confidence'.
%%%                              'modeType' puts the categorical mode type on
%%%                              the y-axis (Asymm/Symm/TH-Asymm/TH-Symm, ticks
%%%                              labelled accordingly): which family this mode
%%%                              belongs to at this velocity. 'confidence' puts
%%%                              modeTypeConfidence there, in [0,1]: how much
%%%                              that label should be trusted. Marker colour
%%%                              always encodes mode type either way, so a
%%%                              confidence view stays linked to which family
%%%                              each point belongs to
%%%  (opt) 3. velCutoffFactor  .. upper velocity cutoff, as a multiple of the
%%%                              structure's bottom-layer (and, if different,
%%%                              top-layer) Christoffel longitudinal limiting
%%%                              velocity (default 1.10). Computed fresh here,
%%%                              since CalcSinglePlate has no stored longLimVel.
%%%                              Pass Inf to disable ([] = the default)
%%%  (opt) 4. hAx              .. target axes (created fresh if omitted)
%%% OUTPUT: hOut .. scatter handle
%%%         hAx  .. axes used

    if nargin < 2 || isempty(criterion); criterion = 'modeType'; end
    visValidateChoice(criterion, {'modeType', 'confidence'}, 'plotRitzModeCurvePlate', 'criterion');
    if nargin < 3 || isempty(velCutoffFactor); velCutoffFactor = 1.10; end
    if nargin < 4 || isempty(hAx)
        figure();
        hAx = axes();
    end
    hold(hAx, 'on');

    if isfinite(velCutoffFactor)
        %%% [] angle: CalcSinglePlate.structure is ALREADY rotated for this
        %%% result's propagation angle, so no further rotation
        velCutoff = velCutoffFactor * visPlateLongLimVel(singleResult.structure, []);
    else
        velCutoff = Inf;
    end

    %%% all three come back as COLUMNS (see visPlateRitzModeData)
    [vels, types, conf] = visPlateRitzModeData(singleResult.ritzSolver);
    keep    = vels <= velCutoff;

    velKeep  = vels(keep);
    typeKeep = types(keep);
    confKeep = conf(keep);

    [velSort, order] = sort(velKeep, 'ascend');
    typeSort = typeKeep(order);
    confSort = confKeep(order);

    [typeNames, typeColors] = Vis.plateModeTypeStyle();
    colSort = typeColors(typeSort, :);

    if strcmp(criterion, 'modeType')
        ySort = typeSort;
    else
        ySort = confSort;
    end

    hOut = scatter(hAx, velSort, ySort, 36, colSort, 'filled');
    xlabel(hAx, 'velocity [km/s]');
    if strcmp(criterion, 'modeType')
        ylabel(hAx, 'mode type');
        hAx.YTick = 1:4;
        hAx.YTickLabel = typeNames;
        hAx.YLim = [0.5, 4.5];
    else
        ylabel(hAx, 'mode-type confidence');
        hAx.YLim = [0, 1];
    end
    title(hAx, sprintf('Ritz mode curve (angle = %.1f deg) -- %s', singleResult.angle, criterion), ...
        'Interpreter', 'none');
end
