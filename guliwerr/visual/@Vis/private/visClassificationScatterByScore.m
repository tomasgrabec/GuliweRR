function hPoints = visClassificationScatterByScore(hAx, xVals, yVals, sizeVals, pts, colorBy, isWeak)
%%% Scatter classified mode points in one series, coloured by score or confidence on
%%% their own grey scale, with their own colour bar
%%% >> hPoints = visClassificationScatterByScore(hAx, xVals, yVals, sizeVals, pts, colorBy, isWeak)
%%%
%%% @Vis-private helper shared by Vis.plotClassificationDiagnostic and
%%% Vis.plotClassificationCurve, the continuous counterpart of
%%% visClassificationScatterByType.
%%%
%%% INPUT: 1. hAx      .. target axes
%%%        2-3. xVals, yVals .. per-point coordinates
%%%        4. sizeVals .. per-point marker size
%%%        5. pts      .. visFlattenSignificantModes output (modeVel, scoreDb, greenPeakVel)
%%%        6. colorBy  .. 'score' | 'confidence' (visClassificationColorVals)
%%%        7. isWeak   .. per-point logical: weak entries are drawn open (face alpha 0)
%%% OUTPUT: hPoints .. the scatter handle. UserData: .isWeak (per point), .cVals (the
%%%                    per-point colour value, NaN where none), .colorbar (the points'
%%%                    colour bar, [] if none could be placed)
%%%
%%% OWN COLOUR SCALE. A MATLAB axes has one colormap and one CLim, and the Green or
%%% Ritz-spectrum background drawn first already owns them (ColMapGrayOrange, log
%%% |G| limits). The points therefore carry explicit RGB colours from `gray` over
%%% their own finite range (black low, white high), and their colour bar sits on a
%%% hidden helper axes placed to the right of the axes' own colour bar, if any; the
%%% background keeps its colormap, limits and colour bar. A point with no value (NaN:
%%% a limiting entry, or 'confidence' with no Green match) is drawn as an open ring.
%%% The helper axes needs a plain container (figure, panel, tab): in a tiledlayout
%%% or uigridlayout the points are still coloured, but no colour bar is drawn (warning
%%% 'Vis:classificationScatterByScore:noColourbar').
%%%
%%% DARK OUTLINE. A black ring is drawn under every point (a second, unpickable
%%% scatter, deleted with hPoints), so white points stay visible on a white or bright
%%% background and weak (open) points keep their colour.
    [cVals, cLabel] = visClassificationColorVals(colorBy, pts.modeVel, pts.scoreDb, pts.greenPeakVel);
    %%% modeVel, not vel: a Green entry's vel IS its greenPeakVel, so the
    %%% confidence |vel - greenPeakVel| would always be 0
    cVals = cVals(:);
    isWeak = isWeak(:);
    hasVal = isfinite(cVals);

    cmap = gray(256);
    if any(hasVal)
        cRange = [min(cVals(hasVal)), max(cVals(hasVal))];
    else
        cRange = [0, 1];
    end
    if cRange(2) == cRange(1); cRange(2) = cRange(1) + 1; end
    idx = round(1 + (cVals - cRange(1)) / diff(cRange) * (size(cmap, 1) - 1));
    rgb = zeros(numel(cVals), 3);                 % no value: black edge, open face
    rgb(hasVal, :) = cmap(min(max(idx(hasVal), 1), size(cmap, 1)), :);

    hHalo = scatter(hAx, xVals, yVals, sizeVals, [0 0 0], 'LineWidth', 3, ...
        'PickableParts', 'none', 'HitTest', 'off', 'HandleVisibility', 'off');
    hPoints = scatter(hAx, xVals, yVals, sizeVals, rgb, 'filled', ...
        'MarkerEdgeColor', 'flat', 'LineWidth', 1.5, 'AlphaData', double(~isWeak & hasVal), ...
        'AlphaDataMapping', 'none', 'MarkerFaceAlpha', 'flat');

    pAx0 = hAx.Position; pc0 = hAx.PositionConstraint;
    [hCb, hCbAx] = localScoreColourbar(hAx, cRange, cLabel);
    hPoints.UserData = struct('isWeak', isWeak, 'cVals', cVals, 'colorbar', hCb);
    hPoints.DeleteFcn = @(~, ~) localDeleteValid({hHalo, hCbAx}, hAx, pAx0, pc0, ~isempty(hCbAx));
end

function [hCb, hCbAx] = localScoreColourbar(hAx, cRange, cLabel)
%%% Colour bar for the points on a hidden helper axes, right of hAx's own colour bar
    hCb = []; hCbAx = [];
    par = hAx.Parent;
    if isa(par, 'matlab.graphics.layout.TiledChartLayout') || isa(par, 'matlab.ui.container.GridLayout')
        warning('Vis:classificationScatterByScore:noColourbar', ...
            'The axes sit in a %s -- the points are coloured, but their colour bar is not drawn.', class(par));
        return
    end
    hasOwnCb = isprop(hAx, 'Colorbar') && ~isempty(hAx.Colorbar) && isvalid(hAx.Colorbar);
    fig = ancestor(hAx, 'figure');
    pAx = hAx.Position;          % already shrunk for hAx's own colour bar, if any
    hCbAx = axes('Parent', par, 'Units', hAx.Units, 'Position', pAx, 'Visible', 'off', ...
        'HitTest', 'off', 'PickableParts', 'none', 'HandleVisibility', 'off', ...
        'Tag', 'visClassificationScoreColourbar');
    colormap(hCbAx, gray(256));
    clim(hCbAx, cRange);
    hCb = colorbar(hCbAx);
    hCb.Label.String = cLabel;
    hCb.Visible = 'on';
    %%% the helper shrank by the width its colour bar needs; hAx gives up the same
    %%% width, so hAx's own colour bar moves left and the helper's takes its old place
    hAx.Position = hCbAx.Position;
    if hasOwnCb
        hCbAx.Position = pAx;
    end
    if ~isempty(fig); fig.CurrentAxes = hAx; end
end

function localDeleteValid(hs, hAx, pAx0, pc0, moved)
    for k = 1:numel(hs)
        if ~isempty(hs{k}) && isgraphics(hs{k}); delete(hs{k}); end
    end
    %%% give back the width the points' colour bar took from hAx, so redraws into the
    %%% same axes do not shrink it step by step
    if moved && isgraphics(hAx) && strcmp(hAx.BeingDeleted, 'off')
        hAx.Position = pAx0; hAx.PositionConstraint = pc0;
    end
end
