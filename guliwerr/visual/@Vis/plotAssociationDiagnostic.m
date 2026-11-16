function [hOut, hAx] = plotAssociationDiagnostic(association, hAx, overlayMode)
%%% Plot every associated point's propagation angle vs. calculated velocity,
%%% coloured by its achieved MAC value -- a mode-tracking debugging view, since
%%% a collapsing MAC is the observable signature of a binding going bad. Points
%%% with modeCollision=true are ringed in red.
%%% >> [hOut, hAx] = Vis.plotAssociationDiagnostic(association, hAx, overlayMode)
%%%
%%% INPUT:  1. association .. ModeAssociation, after at least one rebind()
%%%                          (e.g. inversion.association after calc())
%%%  (opt)  2. hAx          .. axes handle (new figure if omitted)
%%%  (opt)  3. overlayMode  .. false (default): standalone plot with the
%%%                          Vis.ColormapRYG() colormap, a colorbar, labels and
%%%                          title. true: markers ONLY, coloured by explicit
%%%                          per-point RGB so the axes' colormap/CLim are left
%%%                          untouched -- for overlaying on a background that owns
%%%                          them, e.g. Vis.plotGreenMap, whose angle/velocity
%%%                          axes this matches exactly
%%% OUTPUT: hOut .. struct of the scatter handles actually drawn
%%%         hAx  .. axes handle
%%%
%%% 'limiting'-class points are routed through Christoffel, not MAC-tracked
%%% (sidestepping the mode-identity problem entirely), so their macValue is
%%% always NaN -- and scatter silently DROPS NaN-coloured points. They are drawn
%%% as a separate, labelled, colour-neutral series instead.

    if nargin < 2 || isempty(hAx); hAx = axes(figure()); end %#ok<LAXES>
    if nargin < 3 || isempty(overlayMode); overlayMode = false; end

    pts = association.points;
    bindingClass = {pts.bindingClass};
    macTracked = ismember(bindingClass, {'ritzMode', 'pseudoSAW'});

    hold(hAx, 'on');
    hOut = struct();

    %%% not bound to a Ritz mode, so no MAC to diagnose
    if any(~macTracked)
        pNA = pts(~macTracked);
        hOut.notApplicable = scatter(hAx, [pNA.angleDeg], [pNA.calculatedVelocity], ...
            55, [0.55 0.55 0.55], 'd', 'LineWidth', 1.25, 'DisplayName', ...
            sprintf('limiting (Christoffel-routed, no MAC) (%d)', numel(pNA)));
    end

    %%% MAC-tracked points (ritzMode/pseudoSAW). Overlay mode passes scatter an
    %%% explicit [N x 3] RGB spec, which bypasses the axes colormap, instead of
    %%% indexed CData, so a background image's colormap is left untouched.
    pMAC = pts(macTracked);
    angleDeg = [pMAC.angleDeg];
    velCalc  = [pMAC.calculatedVelocity];
    macValue = [pMAC.macValue];
    collided = [pMAC.modeCollision];

    if overlayMode
        macColour = localConfidenceToRGB(macValue);
        hOut.macTracked = scatter(hAx, angleDeg, velCalc, 45, macColour, 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'DisplayName', 'Ritz-tracked, coloured by MAC');
    else
        hOut.macTracked = scatter(hAx, angleDeg, velCalc, 50, macValue, 'filled', ...
            'DisplayName', 'Ritz-tracked, coloured by MAC');
    end
    if any(collided)
        hOut.collided = scatter(hAx, angleDeg(collided), velCalc(collided), 90, 'r', 'o', ...
            'LineWidth', 1.5, 'DisplayName', 'collision / no match');
    end

    nCollided = sum(collided);

    if overlayMode
        return
    end

    cb = colorbar(hAx);
    cb.Label.String = 'MAC value (Ritz-tracked points only)';
    clim(hAx, [0 1]);
    colormap(hAx, Vis.ColormapRYG());   % 0 = bad/red, 1 = good/green; see ColormapRYG

    xlabel(hAx, 'propagation angle [deg]');
    ylabel(hAx, 'calculated velocity [km/s]');
    if nCollided > 0
        title(hAx, sprintf('Association diagnostic (%d point(s) ringed = collision/no match)', nCollided));
    else
        title(hAx, 'Association diagnostic (no collisions)');
    end
    legend(hAx, 'Location', 'best');
    grid(hAx, 'on');
end

function rgb = localConfidenceToRGB(vals)
%%% Map confidence values in [0,1] to explicit RGB rows via Vis.ColormapRYG(),
%%% for use as a scatter colour spec that does not touch the axes colormap.
%%% NaN (e.g. an unmatched point) maps to a neutral mid-grey, not a colormap
%%% extreme -- NaN is "no confidence value", not "zero confidence".
    cMap = Vis.ColormapRYG();
    n = size(cMap, 1);
    vals = vals(:);
    idx = round(1 + max(0, min(1, vals)) * (n - 1));
    rgb = cMap(idx, :);
    rgb(isnan(vals), :) = 0.6;
end
