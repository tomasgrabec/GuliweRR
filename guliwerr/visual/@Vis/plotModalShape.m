function [hOut, hAxCell] = plotModalShape(ritzSolver, modeNo, substrateCutoff, cutoffParam, dirs, meshSpacing, hFig, hAxCellIn)
%%% Plot a Ritz mode's signed 2D (depth x lateral) displacement field
%%% >> [hOut, hAxCell] = plotModalShape(ritzSolver, modeNo, substrateCutoff, cutoffParam, dirs, meshSpacing, hFig, hAxCellIn)
%%%
%%% INPUT: 1. ritzSolver      .. RitzSolver (or subclass) object, eigVecs_con/
%%%                              nullConstrMat already populated for modeNo. There
%%%                              is no retention-aware recompute-on-demand: this
%%%                              errors if eigVecs_con lacks the column, rather
%%%                              than silently doing the wrong thing
%%%        2. modeNo          .. which mode to plot
%%%  (opt) 3. substrateCutoff .. 'none' (default, true-scale) | 'decay' |
%%%                              'fixedMultiple' -- see calcSubstrateCutoffDepth.m
%%%                              -- | 'boosted': no cropping at all (the full
%%%                              substrate depth is always rendered), but row
%%%                              heights are NOT plain-proportional to relDepth
%%%                              either -- see calcModalShapeRowHeights.m. Use
%%%                              this when a thin film next to a deep substrate
%%%                              needs to stay visible without losing any of
%%%                              the substrate's own displayed data.
%%%  (opt) 4. cutoffParam     .. threshold fraction (decay, default 0.05) or
%%%                              multiple (fixedMultiple, default 4); ignored
%%%                              for 'none'/'boosted'
%%%  (opt) 5. dirs            .. which displacement directions to draw (default
%%%                              1:3), or 'rgb': one merged panel per layer,
%%%                              coloured by polarisation direction with
%%%                              brightness = |u| (Vis.calcPolarisationRgb), plus
%%%                              the colour key (Vis.plotPolarisationKey) when
%%%                              drawing into a fresh figure
%%%  (opt) 6. meshSpacing     .. depth/lateral mesh spacing passed to
%%%                              calcDispField2D (default 0.01 — 201x201 points
%%%                              spanning the FULL relDepth range of the deepest
%%%                              layer). A substrateCutoff/ylim zoom into a
%%%                              shallow sliver of that range only sees the points
%%%                              falling inside it, so a zoomed-in caller usually
%%%                              wants this finer than the default
%%%  (opt) 7. hFig            .. target figure (a classic figure() created fresh if
%%%                              omitted — classic figure/axes, not uifigure/
%%%                              uiaxes, this being general-purpose plotting
%%%                              output rather than a GUI window). IGNORED when
%%%                              hAxCellIn is given
%%%  (opt) 8. hAxCellIn       .. pre-built {numel(dirs) x nLay} cell of axes/uiaxes
%%%                              to render into. A caller building a combined
%%%                              uifigure needs this: classic axes() cannot be
%%%                              parented inside a uifigure at all. When given,
%%%                              hFig/clf are skipped and each axes is cla()'d and
%%%                              rendered into directly, with no Position
%%%                              assignment — sizing and placement stay the
%%%                              caller's responsibility. With dirs='rgb': a
%%%                              {1 x nLay} cell, and no key is drawn
%%% OUTPUT: hOut    .. array of pcolor object handles [numel(dirs) x nLay]
%%%                    (dirs='rgb': image handles [1 x nLay])
%%%         hAxCell .. {numel(dirs) x nLay} cell of axes handles ({1 x nLay})
%%%
%%% Layer stacking: layer 1 (substrate) at the BOTTOM of the figure (deepest),
%%% layer nLay (outermost film, nearest the free surface) at the TOP — matching
%%% RitzSolver.calcDisplacement's own row order, where ilay loops 1:nLay filling
%%% rows bottom-up and row "end" of the LAST layer is the true free surface.
%%% Columns are the requested directions. Default colormap ColormapRWB.

    if nargin < 3 || isempty(substrateCutoff); substrateCutoff = 'none'; end
    if nargin < 4; cutoffParam = []; end
    if nargin < 5 || isempty(dirs); dirs = 1:3; end
    if nargin < 6; meshSpacing = []; end
    if nargin < 8; hAxCellIn = []; end
    isRgb = visIsRgbDirs(dirs, 'plotModalShape');
    useGivenAxes = ~isempty(hAxCellIn);
    if ~useGivenAxes && (nargin < 7 || isempty(hFig))
        hFig = figure();
    end

    nColsAvail = size(ritzSolver.eigVecs_con, 2);
    if isempty(ritzSolver.eigVecs_con) || nColsAvail < modeNo
        error('Vis:plotModalShape:modeNotAvailable', ...
            ['Mode %d is not available on this RitzSolver (eigVecs_con has %d column(s)) -- ' ...
             'retention-aware recompute-on-demand is not built yet; ' ...
             'construct/calc this RitzSolver with retentionLevel=''full''.'], modeNo, nColsAvail);
    end

    relDepths = [ritzSolver.structure.calcParams.relDepth];
    nLay      = numel(relDepths);

    dispField = Vis.calcDispField2D(ritzSolver, modeNo, meshSpacing);   % {3 x nLay}

    %%% substrate cutoff (axis-height allocation only; 'none'/'boosted' are a
    %%% no-op here -- 'boosted' never crops data, see layHeightFrac below)
    dispRelDepths = relDepths;
    if strcmp(substrateCutoff, 'decay')
        substrateDispMag = visDepthMagnitude(dispField(:, 1));   % {3 x 1}, this layer's 3 directions
        dispRelDepths(1) = Vis.calcSubstrateCutoffDepth(substrateCutoff, cutoffParam, relDepths, substrateDispMag);
    elseif strcmp(substrateCutoff, 'fixedMultiple')
        dispRelDepths(1) = Vis.calcSubstrateCutoffDepth(substrateCutoff, cutoffParam, relDepths);
    elseif ~strcmp(substrateCutoff, 'none') && ~strcmp(substrateCutoff, 'boosted')
        error('Vis:plotModalShape:badSubstrateCutoff', ...
            'substrateCutoff must be ''none'', ''decay'', ''fixedMultiple'', or ''boosted'' (got ''%s'').', substrateCutoff);
    end

    %%% crop the substrate's field to the displayed depth (keep the shallow end,
    %%% i.e. the LAST rows -- see docblock on row ordering)
    if dispRelDepths(1) < relDepths(1)
        keepFrac = dispRelDepths(1) / relDepths(1);
        for iDir = 1:3
            nRows   = size(dispField{iDir, 1}, 1);
            keepN   = max(1, round(keepFrac * nRows));
            dispField{iDir, 1} = dispField{iDir, 1}(end-keepN+1:end, :);
        end
    end

    %%% axis-height allocation, bottom (substrate) to top (outermost film).
    %%% topMargin/bottomMargin reserve figure-fraction space above and below the
    %%% stack for the top row's per-panel title and any caller-added sgtitle.
    %%% 0.15 at the top, not less: MATLAB does not count an sgtitle's own text
    %%% extent when placing it against a plain, non-tiledlayout figure.
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
            hFig, hAxCellIn, 'plotModalShape');
        return;
    end

    nDirs   = numel(dirs);
    colWidth = 1 / nDirs;
    cMap     = Vis.ColormapRWB();

    if useGivenAxes
        validateattributes(hAxCellIn, {'cell'}, {'size', [nDirs, nLay]}, ...
            'plotModalShape', 'hAxCellIn');
    else
        clf(hFig);
    end
    hOut    = gobjects(nDirs, nLay);
    hAxCell = cell(nDirs, nLay);

    for ilay = 1:nLay
        for ic = 1:nDirs
            iDir = dirs(ic);

            if useGivenAxes
                %%% Sizing/placement is the CALLER's responsibility -- no
                %%% Position assignment here. A caller stacking 2+ panels must
                %%% NOT use a uigridlayout: with 2+ children it re-manages each
                %%% child's PositionConstraint/InnerPosition on every drawnow,
                %%% silently reverting any "fill the whole cell, no uiaxes label
                %%% margin" override. A plain uipanel with each uiaxes manually
                %%% positioned is what makes adjacent layers' black frames meet
                %%% with zero gap.
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
