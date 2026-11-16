function [hOut, hAx] = plotRitzModeMap(sweepResult, criterion, velCutoffFactor, hAx)
%%% Plot every computed Ritz mode -- not just ModeClassification's
%%% significant-mode subset -- as an angle x velocity scatter, shaded by a
%%% mode-recognition criterion.
%%% >> [hOut, hAx] = plotRitzModeMap(sweepResult, criterion, velCutoffFactor, hAx)
%%%
%%% INPUT: 1. sweepResult      .. CalcAngSweepHS or CalcAngSweepLayered result
%%%                              object, calc()'d. Reads .cs{iAngle}.ritzSolver
%%%                              directly (vels/nCompModes/diagnostics), so it
%%%                              works at any retentionLevel: pruneEigVecs only
%%%                              ever touches eigVecs_con
%%%  (opt) 2. criterion        .. 'ritzScore' (default, as in the viewers) |
%%%                              'freqDisp' | 'depthInsensitivity':
%%%                                'ritzScore' -- the Ritz rule's own value: each
%%%                                  mode gets the score (dB) of the Ritz candidate
%%%                                  it supports (Vis.ritzModeScores); modes that
%%%                                  support none are drawn faint
%%%                                'freqDisp' -- the half-space veto quantity; the
%%%                                  colour scale is fixed to [0 1] and shades the
%%%                                  vetoed range below RitzScore.MIN_FREQDISP red
%%%                                  (half-space only: not used on layered stacks,
%%%                                  where the film makes true guided modes
%%%                                  dispersive, so no veto is marked there)
%%%                                'depthInsensitivity' -- how little a mode's
%%%                                  velocity changes with the Ritz domain depth
%%%                                  (log10 1/|dv/dD|): large for surface modes
%%%                              freqDisp / depthInsensitivity are read straight
%%%                              off ritzSolver.diagnostics, never recomputed here
%%%  (opt) 3. velCutoffFactor  .. upper velocity cutoff, as a multiple of the
%%%                              sweep's max longitudinal limiting velocity
%%%                              (default 1.10). Drops the spurious
%%%                              high-frequency Ritz modes that otherwise
%%%                              dominate the velocity axis. Pass Inf to
%%%                              disable ([] = the default)
%%%  (opt) 4. hAx              .. target axes (created fresh if omitted)
%%% OUTPUT: hOut .. scatter handle; 'ritzScore': [scored, faint] (two scatter
%%%                 handles, the faint one empty-but-valid if every mode is scored)
%%%         hAx  .. axes used
%%%
%%% Points are sorted ascending by criterion value before scattering so the
%%% highest-scoring ones draw last, on top: scatter() has no explicit z-order
%%% control, and a handful of bright points would otherwise be buried under a sea
%%% of low-score ones.
%%%
%%% Black axes background with an unflipped `gray` colormap: low values render
%%% dark-on-dark, high values stand out bright white (freqDisp: the vetoed range in
%%% red shades, see visRitzCriterionColour).
%%%
%%% An exploratory diagnostic over the WHOLE Ritz spectrum, not a
%%% results-presentation plot -- plotClassificationDiagnostic is the
%%% classified-modes-only equivalent.

    if nargin < 2 || isempty(criterion); criterion = 'ritzScore'; end
    visValidateChoice(criterion, visRitzCriterion(), 'plotRitzModeMap', 'criterion');
    if nargin < 3 || isempty(velCutoffFactor); velCutoffFactor = 1.10; end
    if nargin < 4 || isempty(hAx)
        figure();
        hAx = axes();
    end

    if isfinite(velCutoffFactor)
        %%% Duck-typed: HS's longVelMap vs Layered's substrate-only subLongVelMap
        longVel = visPropOrEmpty(sweepResult, {'longVelMap', 'subLongVelMap'});
        if isempty(longVel) || all(isnan(longVel))
            %%% 'ritzOnly' sweeps assemble no limiting-velocity map, but every angle
            %%% holds its own (the Ritz score needs it)
            longVel = cellfun(@(c) visScalarOrNaN(visPropOrEmpty(c, {'longLimVel', 'subLongLimVel'})), ...
                sweepResult.cs);
        end
        if isempty(longVel) || all(isnan(longVel))
            error('Vis:plotRitzModeMap:noLongVel', ...
                ['longitudinal limiting velocity is unavailable (calcMode was not ''full'') -- ' ...
                 'cannot apply velCutoffFactor; pass Inf to disable it.']);
        end
        velCutoff = velCutoffFactor * max(longVel);
    else
        velCutoff = Inf;
    end

    angles  = sweepResult.angles(:);
    nAngles = numel(angles);
    allAng = cell(nAngles, 1); allVel = cell(nAngles, 1); allCrit = cell(nAngles, 1);
    for iAngle = 1:nAngles
        ritz = sweepResult.cs{iAngle}.ritzSolver;
        vels = ritz.vels(1:ritz.nCompModes);
        [critVal, info] = visRitzCriterion(sweepResult.cs{iAngle}, criterion);
        keep = vels(:) <= velCutoff;
        allAng{iAngle}  = repmat(angles(iAngle), sum(keep), 1);
        allVel{iAngle}  = vels(keep);
        allCrit{iAngle} = critVal(keep);
    end
    allAng = vertcat(allAng{:}); allVel = vertcat(allVel{:}); allCrit = vertcat(allCrit{:});

    %%% 'ritzScore': modes supporting no candidate (NaN) first, faint and small
    isFaint = isnan(allCrit);
    hFaint = gobjects(0);
    wasHeld = ishold(hAx);
    if strcmp(criterion, 'ritzScore')
        hFaint = scatter(hAx, allAng(isFaint), allVel(isFaint), 6, [0.35 0.35 0.35], 'filled');
        hold(hAx, 'on');
    end
    [critSort, order] = sort(allCrit(~isFaint), 'ascend');   % highest plotted last (on top)
    angKeep = allAng(~isFaint); velKeep = allVel(~isFaint);

    hOut = scatter(hAx, angKeep(order), velKeep(order), 20, critSort, 'filled');
    hOut = [hOut, hFaint];
    if ~wasHeld; hold(hAx, 'off'); end
    hAx.Color = [0, 0, 0];
    hCb = colorbar(hAx);
    hCb.Label.String = info.label;   % makes a criterion switch unmistakable: colour
                                     % auto-scales to each criterion's own min/max,
                                     % so two correlated criteria can look alike
    visRitzCriterionColour(hAx, hCb, info);
    xlabel(hAx, 'angle [deg]');
    ylabel(hAx, 'velocity [km/s]');
    title(hAx, sprintf('Ritz mode map -- %s', info.label), 'Interpreter', 'none');
end

function v = visScalarOrNaN(v)
%%% A scalar limiting velocity, or NaN when the property is empty
    if isempty(v); v = NaN; end
end
