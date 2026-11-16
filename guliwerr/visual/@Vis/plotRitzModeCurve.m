function [hOut, hAx] = plotRitzModeCurve(singleResult, criterion, velCutoffFactor, hAx)
%%% Plot every computed Ritz mode (not just ModeClassification's
%%% significant-mode subset) as a velocity x criterion curve, for a SINGLE
%%% propagation angle -- 1D analogue of plotRitzModeMap.m for CalcSingleHS/
%%% CalcSingleLayered results (no angle axis to sweep).
%%% >> [hOut, hAx] = plotRitzModeCurve(singleResult, criterion, velCutoffFactor, hAx)
%%%
%%% INPUT: 1. singleResult     .. CalcSingleHS or CalcSingleLayered result
%%%                              object, calcAll()'d. Reads .ritzSolver directly
%%%                              (vels/nCompModes/diagnostics), so it works at any
%%%                              retentionLevel
%%%  (opt) 2. criterion        .. 'ritzScore' (default) | 'freqDisp' |
%%%                              'depthInsensitivity', the same choices plotRitzModeMap
%%%                              offers (see its docblock). 'freqDisp' on a
%%%                              half-space also draws the veto level
%%%                              (RitzScore.MIN_FREQDISP) as a dashed red line;
%%%                              'ritzScore' draws the scored modes as a line with
%%%                              small markers and the modes supporting no
%%%                              candidate as faint grey dots on a baseline below
%%%  (opt) 3. velCutoffFactor  .. upper velocity cutoff, as a multiple of
%%%                              singleResult.longLimVel (default 1.10) — as in
%%%                              plotRitzModeMap, but against a scalar limiting
%%%                              velocity instead of a per-angle map. Pass Inf
%%%                              to disable ([] = the default)
%%%  (opt) 4. hAx              .. target axes (created fresh if omitted)
%%% OUTPUT: hOut .. line handle; 'ritzScore': [scored, faint] (two line handles,
%%%                 the faint one markers only)
%%%         hAx  .. axes used
%%%
%%% A line rather than plotRitzModeMap's scatter, whose colour encodes the
%%% criterion over a 2D (angle, velocity) domain. Points are sorted ascending
%%% by velocity first: plot() draws strictly in array order, so unsorted points
%%% would zig-zag between Ritz modes instead of tracing one curve across the
%%% velocity axis.

    if nargin < 2 || isempty(criterion); criterion = 'ritzScore'; end
    visValidateChoice(criterion, visRitzCriterion(), 'plotRitzModeCurve', 'criterion');
    if nargin < 3 || isempty(velCutoffFactor); velCutoffFactor = 1.10; end
    if nargin < 4 || isempty(hAx)
        figure();
        hAx = axes();
    end

    if isfinite(velCutoffFactor)
        %%% Duck-typed: longLimVel (HS) vs subLongLimVel (layered)
        longLimVel = visPropOrEmpty(singleResult, {'longLimVel', 'subLongLimVel'});
        if isempty(longLimVel)
            error('Vis:plotRitzModeCurve:noLongVel', ...
                ['longitudinal limiting velocity is unavailable (calcMode was not ''full'') -- ' ...
                 'cannot apply velCutoffFactor; pass Inf to disable it.']);
        end
        velCutoff = velCutoffFactor * longLimVel;
    else
        velCutoff = Inf;
    end

    ritz    = singleResult.ritzSolver;
    nM      = ritz.nCompModes;
    vels    = ritz.vels(1:nM);
    [critVal, info] = visRitzCriterion(singleResult, criterion);
    keep    = vels(:) <= velCutoff;

    [velSort, order] = sort(vels(keep), 'ascend');
    critValKeep = critVal(keep);
    critSort    = critValKeep(order);

    wasHeld = ishold(hAx);
    if strcmp(criterion, 'ritzScore')
        %%% discrete per-mode values, most of them NaN: the scored modes as a line
        %%% with small markers; the modes supporting no candidate as faint dots just
        %%% below the lowest score, so the velocity axis still shows where every
        %%% mode is
        isFaint = isnan(critSort);
        if all(isFaint); yBase = 0; else; yBase = min(critSort(~isFaint)) - 6; end
        hScored = plot(hAx, velSort(~isFaint), critSort(~isFaint), '-o', 'MarkerSize', 4, ...
            'Color', [0 0.4470 0.7410], 'MarkerFaceColor', [0 0.4470 0.7410]);
        hold(hAx, 'on');
        hFaint = plot(hAx, velSort(isFaint), repmat(yBase, sum(isFaint), 1), '.', ...
            'Color', [0.7 0.7 0.7]);
        hOut = [hScored, hFaint];
    else
        hOut = plot(hAx, velSort, critSort, '-');
        if ~isnan(info.vetoLevel)
            hold(hAx, 'on');
            yline(hAx, info.vetoLevel, '--', 'veto', 'Color', [0.8 0.1 0.1], ...
                'HandleVisibility', 'off');
        end
    end
    if ~wasHeld; hold(hAx, 'off'); end
    xlabel(hAx, 'velocity [km/s]');
    ylabel(hAx, info.label, 'Interpreter', 'none');
    title(hAx, sprintf('Ritz mode curve (angle = %.1f deg) -- %s', singleResult.angle, info.label), ...
        'Interpreter', 'none');
end
