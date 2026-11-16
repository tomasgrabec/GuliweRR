function h = visModalShapePanel(ax, field2D, cMap, titleDir)
%%% Render one signed displacement-field panel into an existing axes
%%% >> h = visModalShapePanel(ax, field2D, cMap, titleDir)
%%%
%%% @Vis-private helper: the shared visual language of Vis.plotModalShape and
%%% Vis.plotGreenModalShape, so a Ritz modal shape and a Green's response can sit
%%% side by side without a visual-language shift.
%%%
%%% INPUT: 1. ax       .. target axes/uiaxes, already positioned and cla'd by the
%%%                       caller — placement stays the caller's business
%%%        2. field2D  .. [nMeshDepth x nMeshLateral] signed real field
%%%        3. cMap     .. colormap (Vis.ColormapRWB for both current callers)
%%%  (opt) 4. titleDir .. 1-3 to title the panel with that displacement
%%%                       direction's name; [] or omitted for no title
%%% OUTPUT: h .. the pcolor object handle
%%%
%%% Both axes are drawn on a normalised 0-1 range: lateral = x/wavelength, depth =
%%% fraction of this panel's own displayed depth. The callers' fields are globally
%%% normalised, so the fixed clim([-1,1]) is what makes panels comparable; it
%%% clips a negative-dominant mode by design (see visDispField2D).
%%%
%%% The thin black frame is drawn explicitly: axis('off') hides the axes Box along
%%% with the ticks and labels, and ColormapRWB's white zero-crossing needs a
%%% visible boundary against the figure's own white background.

    nMeshDepth   = size(field2D, 1);
    nMeshLateral = size(field2D, 2);
    xNorm = linspace(0, 1, nMeshLateral);   % x / wavelength
    yNorm = linspace(0, 1, nMeshDepth);      % fraction of this panel's displayed depth

    h = pcolor(ax, xNorm, yNorm, field2D);
    shading(ax, 'flat');
    colormap(ax, cMap);
    clim(ax, [-1, 1]);
    axis(ax, 'off');

    hold(ax, 'on');
    rectangle(ax, 'Position', [0, 0, 1, 1], 'EdgeColor', 'k', 'LineWidth', 1);

    if nargin >= 4 && ~isempty(titleDir)
        %%% Standard convention: x1=propagation, x2=in-plane, x3=surface
        %%% normal (matches Material.GetStandardDirMat / GreenSolver).
        dirNames = {'u_x (propag.dir.)', 'u_y (in-plane)', 'u_z (out-of-plane)'};
        title(ax, dirNames{titleDir});
    end
end
