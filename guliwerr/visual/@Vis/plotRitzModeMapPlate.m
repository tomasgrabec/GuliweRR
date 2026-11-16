function [hOut, hAx] = plotRitzModeMapPlate(sweepResult, criterion, velCutoffFactor, hAx)
%%% Plot every computed Ritz mode for a CalcAngSweepPlate result as an
%%% angle x velocity scatter, shaded by mode-type classification -- Plate
%%% analogue of plotRitzModeMap.m.
%%% >> [hOut, hAx] = Vis.plotRitzModeMapPlate(sweepResult, criterion, velCutoffFactor, hAx)
%%%
%%% Separate from plotRitzModeMap rather than a class-dispatch branch inside it:
%%% mode-TYPE classification (RitzSolverLamb.modesType/modeTypeConfidence)
%%% replaces SAW's diagnostics criteria, and there is no Green's-function map
%%% option at all.
%%%
%%% INPUT: 1. sweepResult      .. CalcAngSweepPlate result object, calc()'d.
%%%                              Reads .cs{iAngle}.ritzSolver directly
%%%                              (vels/nCompModes/modesType/modeTypeConfidence),
%%%                              so it works at any retentionLevel
%%%  (opt) 2. criterion        .. 'modeType' (default) | 'confidence'. Both drive
%%%                              marker colour directly, x/y already being
%%%                              angle/velocity: 'modeType' uses the same
%%%                              4-colour map plotDispersionCurvesPlate uses with
%%%                              a 4-entry legend, a discrete classification
%%%                              having no natural colour bar; 'confidence' uses a
%%%                              continuous `gray` colormap plus colour bar, as
%%%                              plotRitzModeMap does
%%%  (opt) 3. velCutoffFactor  .. upper velocity cutoff, as a multiple of the
%%%                              MAXIMUM over every angle of the structure's
%%%                              bottom-layer (and, if different, top-layer)
%%%                              Christoffel longitudinal limiting velocity
%%%                              (default 1.10). Computed fresh here, since
%%%                              CalcAngSweepPlate has no stored longVelMap.
%%%                              Pass Inf to disable ([] = the default)
%%%  (opt) 4. hAx              .. target axes (created fresh if omitted)
%%% OUTPUT: hOut .. scatter handle (the real data; legend/colorbar proxies are
%%%                 separate handles, not returned)
%%%         hAx  .. axes used
%%%
%%% Points are sorted so 'modeType' draws in a fixed, deterministic order and
%%% 'confidence' draws highest-confidence last, on top.

    if nargin < 2 || isempty(criterion); criterion = 'modeType'; end
    visValidateChoice(criterion, {'modeType', 'confidence'}, 'plotRitzModeMapPlate', 'criterion');
    if nargin < 3 || isempty(velCutoffFactor); velCutoffFactor = 1.10; end
    if nargin < 4 || isempty(hAx)
        figure();
        hAx = axes();
    end
    hold(hAx, 'on');

    angles  = sweepResult.angles(:);
    nAngles = numel(angles);

    if isfinite(velCutoffFactor)
        %%% CalcAngSweepPlate.structure is kept UNROTATED, so each angle is
        %%% rotated to inside the shared helper
        longLimByAngle = arrayfun(@(a) visPlateLongLimVel(sweepResult.structure, a), angles);
        velCutoff = velCutoffFactor * max(longLimByAngle);
    else
        velCutoff = Inf;
    end

    allAng = []; allVel = []; allType = []; allConf = [];
    for iAngle = 1:nAngles
        %%% all three come back as COLUMNS -- the orientation matters, see
        %%% visPlateRitzModeData
        [vels, types, conf] = visPlateRitzModeData(sweepResult.cs{iAngle}.ritzSolver);
        keep = vels <= velCutoff;
        allAng  = [allAng;  repmat(angles(iAngle), sum(keep), 1)]; %#ok<AGROW>
        allVel  = [allVel;  vels(keep)];                           %#ok<AGROW>
        allType = [allType; types(keep)];                          %#ok<AGROW>
        allConf = [allConf; conf(keep)];                           %#ok<AGROW>
    end

    [typeNames, typeColors] = Vis.plateModeTypeStyle();

    if strcmp(criterion, 'modeType')
        [typeSort, order] = sort(allType, 'ascend');
        angSort = allAng(order);
        velSort = allVel(order);
        colSort = typeColors(typeSort, :);

        hOut = scatter(hAx, angSort, velSort, 20, colSort, 'filled', 'HandleVisibility', 'off');
        for it = 1:4
            scatter(hAx, NaN, NaN, 36, typeColors(it, :), 'filled', 'DisplayName', typeNames{it});
        end
        legend(hAx, 'Location', 'best');
        hAx.Color = [1, 1, 1];
    else
        [confSort, order] = sort(allConf, 'ascend');   % highest confidence plotted last (on top)
        angSort = allAng(order);
        velSort = allVel(order);

        hOut = scatter(hAx, angSort, velSort, 20, confSort, 'filled');
        colormap(hAx, gray);
        hAx.Color = [0, 0, 0];
        hCb = colorbar(hAx);
        hCb.Label.String = 'mode-type confidence';
        clim(hAx, [0, 1]);
    end

    xlabel(hAx, 'angle [deg]');
    ylabel(hAx, 'velocity [km/s]');
    title(hAx, sprintf('Ritz mode map (Plate) -- %s', criterion), 'Interpreter', 'none');
end
