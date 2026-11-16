function [hOut, hAxCell] = visRgbShapeStack(dispField, layBottom, layHeightFrac, hFig, hAxCellIn, fnName)
%%% Render a displacement field as one polarisation-RGB panel per layer
%%% >> [hOut, hAxCell] = visRgbShapeStack(dispField, layBottom, layHeightFrac, hFig, hAxCellIn, fnName)
%%%
%%% @Vis-private helper: the dirs = 'rgb' branch of Vis.plotModalShape,
%%% plotGreenModalShape and plotGreenModalShapeLayered. The callers have already
%%% cropped the field and computed the per-layer axis heights, so this only
%%% replaces the "one signed panel per direction" loop.
%%%
%%% INPUT: 1. dispField     .. {3 x nLay} signed real fields, already cropped
%%%        2. layBottom     .. [1 x nLay] figure-fraction bottom of each layer row
%%%        3. layHeightFrac .. [1 x nLay] figure-fraction height of each row
%%%        4. hFig          .. classic figure to draw into (ignored with hAxCellIn)
%%%        5. hAxCellIn     .. [] or a {1 x nLay} cell of axes/uiaxes to render
%%%                            into (same contract as the callers' hAxCellIn:
%%%                            no Position assignment, no key -- a GUI draws
%%%                            Vis.plotPolarisationKey into an axes of its own)
%%%        6. fnName        .. calling function's name, for error ids
%%% OUTPUT: hOut    .. [1 x nLay] image handles
%%%         hAxCell .. {1 x nLay} axes handles
%%%
%%% In a classic figure the panel sits in the left part (about the width of
%%% the three-column view's two left columns) and the colour key on the right.

    nLay = size(dispField, 2);
    [rgbCell, ampCell] = Vis.calcPolarisationRgb(dispField);

    useGivenAxes = ~isempty(hAxCellIn);
    if useGivenAxes
        validateattributes(hAxCellIn, {'cell'}, {'size', [1, nLay]}, fnName, 'hAxCellIn');
    else
        clf(hFig);
    end

    hOut    = gobjects(1, nLay);
    hAxCell = cell(1, nLay);
    for ilay = 1:nLay
        if useGivenAxes
            ax = hAxCellIn{1, ilay};
            cla(ax);
        else
            ax = axes('Parent', hFig, 'Position', [0.04, layBottom(ilay), 0.5, layHeightFrac(ilay)]);
        end
        hOut(ilay)    = polarisationPanel(ax, rgbCell{ilay}, ampCell{ilay}, ilay == nLay);
        hAxCell{ilay} = ax;
    end

    if ~useGivenAxes
        keyAx = axes('Parent', hFig, 'Position', [0.6, 0.35, 0.38, 0.38]);
        Vis.plotPolarisationKey(keyAx);
    end
end

function h = polarisationPanel(ax, rgb, amp, isTop)
%%% One panel: the RGB image on the same normalised 0-1 axes as
%%% visModalShapePanel, thin amplitude contours (as in the paper's figures),
%%% and the black frame.
    h = image(ax, [0, 1], [0, 1], rgb);
    ax.YDir = 'normal';
    hold(ax, 'on');

    if max(amp(:)) - min(amp(:)) > eps
        xNorm = linspace(0, 1, size(amp, 2));
        yNorm = linspace(0, 1, size(amp, 1));
        contour(ax, xNorm, yNorm, amp, 0.2:0.2:0.8, ...
            'LineColor', [0.5, 0.5, 0.5], 'LineWidth', 0.3);
    end

    xlim(ax, [0, 1]);
    ylim(ax, [0, 1]);
    axis(ax, 'off');
    rectangle(ax, 'Position', [0, 0, 1, 1], 'EdgeColor', 'k', 'LineWidth', 1);
    if isTop
        title(ax, 'u (polarisation RGB)');
    end
end
