function [hOut, hAx] = plotGreenCurve(singleResult, scaleMode, hAx, component)
%%% Plot the Green's-function velocity curve -- |G| vs. velocity, for a
%%% SINGLE propagation angle -- 1D analogue of plotGreenMap.m for
%%% CalcSingleHS/CalcSingleLayered results (no angle axis to sweep).
%%% >> [hOut, hAx] = plotGreenCurve(singleResult, scaleMode, hAx, component)
%%%
%%% INPUT: 1. singleResult .. CalcSingleHS or CalcSingleLayered result object,
%%%                          calcAll()'d with calcMode='full' — the Green fields
%%%                          stay empty under 'ritzOnly', which errors clearly
%%%                          here. .greenVelAxis and the component's curve are
%%%                          read directly, both classes naming them the same.
%%%                          A sweep's maps are these same curves stacked per angle
%%%  (opt) 2. scaleMode    .. 'log' (default) | 'linear' -- the same
%%%                          normalisation plotGreenMap uses: 'log' plots
%%%                          log10(|G| / max|G|), normalised to this curve's
%%%                          own peak
%%%  (opt) 3. hAx          .. target axes (created fresh if omitted)
%%%  (opt) 4. component    .. 'G13eff' (default, the stored |G13_eff| -- the TGS
%%%                          response the existence score reads) |
%%%                          any 'GIJ' / 'GIJeff', e.g. 'G13' (one-sided G13_+),
%%%                          'G33': recomputed for this call, not stored -- see
%%%                          Vis.greenComponentData
%%% OUTPUT: hOut .. line handle
%%%         hAx  .. axes used
%%%
%%% No climRange/plotType/colorbar parameters: those are plotGreenMap's
%%% heatmap-specific concerns, and a single-angle curve has no second axis to
%%% spend colour on.

    if nargin < 4 || isempty(component); component = 'G13eff'; end
    visValidateGreenComponent(component, 'plotGreenCurve');
    [absCurve, componentLabel, velAxis] = Vis.greenComponentData(singleResult, component);

    if nargin < 2 || isempty(scaleMode); scaleMode = 'log'; end
    visValidateChoice(scaleMode, {'log', 'linear'}, 'plotGreenCurve', 'scaleMode');
    if nargin < 3 || isempty(hAx)
        figure();
        hAx = axes();
    end

    [plotData, yLabel] = Vis.scaleGreenData(absCurve, scaleMode, componentLabel);

    hOut = plot(hAx, velAxis, plotData, '-', 'DisplayName', ['|' componentLabel '|']);
    xlabel(hAx, 'velocity [km/s]');
    ylabel(hAx, yLabel);
    title(hAx, sprintf('Green''s function curve (%s, angle = %.1f deg)', componentLabel, singleResult.angle));
end
