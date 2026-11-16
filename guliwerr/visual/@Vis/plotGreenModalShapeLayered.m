function [hOut, hAxCell] = plotGreenModalShapeLayered(greenObj, structure, forceDir, substrateCutoff, cutoffParam, dirs, meshSpacing, hFig, hAxCellIn)
%%% Plot a LAYERED Green's-function response's signed 2D (depth x lateral)
%%% displacement field, for one already-selected velocity and force
%%% direction. Layered counterpart of the half-space-only plotGreenModalShape;
%%% its multi-layer axis stacking mirrors plotModalShape.m instead.
%%% >> [hOut, hAxCell] = plotGreenModalShapeLayered(greenObj, structure, forceDir, substrateCutoff, cutoffParam, dirs, meshSpacing, hFig, hAxCellIn)
%%%
%%% INPUT: 1. greenObj        .. GreenSolver object from GreenSolver.calc called
%%%                              with nargout>1 and a SCALAR velocity on a
%%%                              LAYERED structure -- see
%%%                              calcGreenDispField2DLayered
%%%        2. structure       .. the SAME Structure passed to GreenSolver.calc
%%%        3. forceDir        .. which force direction's response to show (1-3)
%%%  (opt) 4. substrateCutoff .. 'none' | 'decay' (default) | 'fixedMultiple' |
%%%                              'boosted' -- see calcSubstrateCutoffDepth.m.
%%%                              Defaults to 'decay', not plotModalShape's
%%%                              'none': GreenSolver's substrate is analytically
%%%                              semi-infinite, so its display extent is not a
%%%                              physical "true scale" the way a Ritz substrate
%%%                              truncation is, and shown in full it is usually
%%%                              a mostly-empty panel
%%%  (opt) 5. cutoffParam     .. threshold fraction (decay, default 0.05) or
%%%                              multiple (fixedMultiple, default 4); ignored
%%%                              for 'none'/'boosted'
%%%  (opt) 6. dirs            .. which displacement directions to draw (default
%%%                              1:3), or 'rgb' for the merged polarisation view
%%%                              -- as plotModalShape's dirs
%%%  (opt) 7. meshSpacing     .. depth/lateral mesh spacing passed to
%%%                              calcGreenDispField2DLayered (default 0.01)
%%%  (opt) 8. hFig            .. target figure (a classic figure() created fresh
%%%                              if omitted). IGNORED when hAxCellIn is given
%%%  (opt) 9. hAxCellIn       .. pre-built {numel(dirs) x nLay} cell of
%%%                              axes/uiaxes to render into -- same contract as
%%%                              plotModalShape's hAxCellIn ({1 x nLay} with
%%%                              dirs='rgb')
%%% OUTPUT: hOut    .. array of pcolor object handles [numel(dirs) x nLay]
%%%                    (dirs='rgb': image handles [1 x nLay])
%%%         hAxCell .. {numel(dirs) x nLay} cell of axes handles ({1 x nLay})
%%%
%%% Layer stacking as in plotModalShape.m: layer 1 (substrate) at the BOTTOM of
%%% the figure, layer nLay (outermost film) at the TOP. Default colormap
%%% ColormapRWB.

    if nargin < 4 || isempty(substrateCutoff); substrateCutoff = 'decay'; end
    if nargin < 5; cutoffParam = []; end
    if nargin < 6 || isempty(dirs); dirs = 1:3; end
    if nargin < 7; meshSpacing = []; end
    if nargin < 9; hAxCellIn = []; end
    isRgb = visIsRgbDirs(dirs, 'plotGreenModalShapeLayered');
    useGivenAxes = ~isempty(hAxCellIn);
    if ~useGivenAxes && (nargin < 8 || isempty(hFig))
        hFig = figure();
    end

    relDepths = [structure.calcParams.relDepth];
    nLay      = numel(relDepths);

    dispField = Vis.calcGreenDispField2DLayered(greenObj, structure, forceDir, meshSpacing);   % {3 x nLay}

    %%% substrate cutoff (axis-height allocation AND data crop), as in
    %%% plotModalShape.m; dispField{iDir,1} starts at the FULL, uncropped
    %%% relDepths(1) extent, just as a Ritz substrate does there
    dispRelDepths = relDepths;
    if strcmp(substrateCutoff, 'decay')
        substrateDispMag = visDepthMagnitude(dispField(:, 1));   % {3 x 1}, this layer's 3 directions
        dispRelDepths(1) = Vis.calcSubstrateCutoffDepth(substrateCutoff, cutoffParam, relDepths, substrateDispMag);
    elseif strcmp(substrateCutoff, 'fixedMultiple')
        dispRelDepths(1) = Vis.calcSubstrateCutoffDepth(substrateCutoff, cutoffParam, relDepths);
    elseif ~strcmp(substrateCutoff, 'none') && ~strcmp(substrateCutoff, 'boosted')
        error('Vis:plotGreenModalShapeLayered:badSubstrateCutoff', ...
            'substrateCutoff must be ''none'', ''decay'', ''fixedMultiple'', or ''boosted'' (got ''%s'').', substrateCutoff);
    end

    %%% crop the substrate's field to the displayed depth, keeping the shallow
    %%% end -- the LAST rows, since row 1 is deepest
    if dispRelDepths(1) < relDepths(1)
        keepFrac = dispRelDepths(1) / relDepths(1);
        for iDir = 1:3
            nRows   = size(dispField{iDir, 1}, 1);
            keepN   = max(1, round(keepFrac * nRows));
            dispField{iDir, 1} = dispField{iDir, 1}(end-keepN+1:end, :);
        end
    end

    %%% axis-height allocation, bottom (substrate) to top (outermost film), as
    %%% in plotModalShape.m
    topMargin    = 0.15;
    bottomMargin = 0.025;
    plotHeight   = 1 - topMargin - bottomMargin;

    if strcmp(substrateCutoff, 'boosted')
        layHeightFrac = Vis.calcModalShapeRowHeights(dispRelDepths) * plotHeight;
    else
        totalDepth    = sum(dispRelDepths);
        layHeightFrac = (dispRelDepths / totalDepth) * plotHeight;
    end
    layBottom = bottomMargin + [0, cumsum(layHeightFrac(1:end-1))];

    if isRgb
        if useGivenAxes; hFig = []; end   % ignored with hAxCellIn, may be unset
        [hOut, hAxCell] = visRgbShapeStack(dispField, layBottom, layHeightFrac, ...
            hFig, hAxCellIn, 'plotGreenModalShapeLayered');
        return;
    end

    nDirs    = numel(dirs);
    colWidth = 1 / nDirs;
    cMap     = Vis.ColormapRWB();

    if useGivenAxes
        validateattributes(hAxCellIn, {'cell'}, {'size', [nDirs, nLay]}, ...
            'plotGreenModalShapeLayered', 'hAxCellIn');
    else
        clf(hFig);
    end
    hOut    = gobjects(nDirs, nLay);
    hAxCell = cell(nDirs, nLay);

    for ilay = 1:nLay
        for ic = 1:nDirs
            iDir = dirs(ic);

            if useGivenAxes
                %%% see plotModalShape's hAxCellIn branch for why this must
                %%% be a manually-positioned uipanel/uiaxes, not a
                %%% uigridlayout, when the caller stacks 2+ panels
                ax = hAxCellIn{ic, ilay};
                cla(ax);
            else
                ax = axes('Parent', hFig, 'Position', ...
                    [(ic-1)*colWidth + 0.02*colWidth, layBottom(ilay), 0.96*colWidth, layHeightFrac(ilay)]);
            end
            %%% direction titles on the TOP row only (the outermost film) --
            %%% one label per column, not one per panel
            if ilay == nLay; titleDir = iDir; else; titleDir = []; end
            h = visModalShapePanel(ax, dispField{iDir, ilay}, cMap, titleDir);

            hOut(ic, ilay)    = h;
            hAxCell{ic, ilay} = ax;
        end
    end
end
