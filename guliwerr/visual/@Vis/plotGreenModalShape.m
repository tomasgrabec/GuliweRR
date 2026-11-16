function [hOut, hAxCell] = plotGreenModalShape(greenObj, forceDir, cutoffParam, dirs, meshSpacing, hFig, hAxCellIn)
%%% Plot a Green's-function response's signed 2D (depth x lateral) displacement
%%% field -- Green's-response analogue of plotModalShape, for one already-
%%% selected velocity and force direction.
%%% >> [hOut, hAxCell] = plotGreenModalShape(greenObj, forceDir, cutoffParam, dirs, meshSpacing, hFig, hAxCellIn)
%%%
%%% INPUT: 1. greenObj    .. GreenSolver object, .dispField populated -- see
%%%                          calcGreenDispField2D's own docblock (call
%%%                          GreenSolver.calc with nargout>1 and a SCALAR
%%%                          velocity first; half-space only)
%%%        2. forceDir    .. which force direction's response to show (1-3)
%%%  (opt) 3. cutoffParam .. decay threshold fraction, default 0.02: the display
%%%                          truncates once the depth-profile magnitude drops
%%%                          below this fraction of its value AT the surface.
%%%                          Always a plain decay cutoff, unlike
%%%                          plotModalShape's substrateCutoff — a Green's
%%%                          response is ALWAYS a half-space, so there is no
%%%                          film-stack choice to make and
%%%                          calcSubstrateCutoffDepth's nLay>=2 guard does not
%%%                          fit. The crop is applied here directly, but the
%%%                          decay SCAN is the same visDecayCutoffIdx helper
%%%  (opt) 4. dirs         .. which displacement directions to draw (default
%%%                           1:3), or 'rgb' for the merged polarisation view
%%%                           -- as plotModalShape's dirs
%%%  (opt) 5. meshSpacing  .. lateral mesh spacing (default 0.01)
%%%  (opt) 6. hFig         .. target figure (a classic figure() created fresh if
%%%                           omitted, as in plotModalShape). IGNORED when
%%%                           hAxCellIn is given
%%%  (opt) 7. hAxCellIn    .. pre-built {numel(dirs) x 1} cell of axes/uiaxes to
%%%                           render into -- same contract as plotModalShape's
%%%                           hAxCellIn: skips hFig/clf/Position entirely, and
%%%                           cla()'s and renders into each given handle.
%%%                           With dirs='rgb': a {1 x 1} cell, no key drawn
%%% OUTPUT: hOut    .. array of pcolor handles [numel(dirs) x 1] (dirs='rgb':
%%%                    one image handle)
%%%         hAxCell .. {numel(dirs) x 1} cell of axes handles ({1 x 1})
%%%
%%% Visual language matches plotModalShape (ColormapRWB, [-1,1] clim, thin black
%%% frame, axis off, direction titles) so a Ritz modal shape and a Green's
%%% response can sit side by side — both render each panel through the same
%%% visModalShapePanel helper. Always exactly one "layer": Green's full output is
%%% half-space only, so there is no multi-layer stacking to allocate.

    if nargin < 3 || isempty(cutoffParam); cutoffParam = 0.02; end
    validateattributes(cutoffParam, {'double'}, {'scalar', '>', 0, '<', 1}, ...
        'plotGreenModalShape', 'cutoffParam');
    if nargin < 4 || isempty(dirs); dirs = 1:3; end
    if nargin < 5; meshSpacing = []; end
    if nargin < 7; hAxCellIn = []; end
    isRgb = visIsRgbDirs(dirs, 'plotGreenModalShape');
    useGivenAxes = ~isempty(hAxCellIn);
    if ~useGivenAxes && (nargin < 6 || isempty(hFig))
        hFig = figure();
    end

    dispField = Vis.calcGreenDispField2D(greenObj, forceDir, meshSpacing);   % {3 x 1}

    %%% decay-based crop: row index END is the free surface
    %%% (calcGreenDispField2D's re-oriented convention), so the scan runs inward
    %%% from the surface.
    keepIdx = visDecayCutoffIdx(visDepthMagnitude(dispField), cutoffParam);
    for iDir = 1:3
        dispField{iDir} = dispField{iDir}(keepIdx:end, :);
    end

    if isRgb
        if useGivenAxes; hFig = []; end   % ignored with hAxCellIn, may be unset
        [hOut, hAxCell] = visRgbShapeStack(dispField, 0.05, 0.85, ...
            hFig, hAxCellIn, 'plotGreenModalShape');
        hOut = hOut(:); hAxCell = hAxCell(:);
        return;
    end

    nDirs    = numel(dirs);
    colWidth = 1 / nDirs;
    cMap     = Vis.ColormapRWB();

    if useGivenAxes
        validateattributes(hAxCellIn, {'cell'}, {'size', [nDirs, 1]}, ...
            'plotGreenModalShape', 'hAxCellIn');
    else
        clf(hFig);
    end
    hOut    = gobjects(nDirs, 1);
    hAxCell = cell(nDirs, 1);

    for ic = 1:nDirs
        iDir = dirs(ic);

        if useGivenAxes
            %%% see plotModalShape's hAxCellIn branch: a uigridlayout-managed
            %%% axes cannot reliably have its label margin overridden, so
            %%% callers stacking cells need a manually-positioned uipanel
            ax = hAxCellIn{ic};
            cla(ax);
        else
            ax = axes('Parent', hFig, 'Position', ...
                [(ic-1)*colWidth + 0.02*colWidth, 0.05, 0.96*colWidth, 0.85]);
        end
        %%% single row, so every panel is titled -- plotModalShape titles its
        %%% top row only
        hOut(ic)    = visModalShapePanel(ax, dispField{iDir}, cMap, iDir);
        hAxCell{ic} = ax;
    end
end
