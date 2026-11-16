function [hPatch, hAx] = plotPolarisationKey(hAx, nSub)
%%% Colour key of the polarisation-RGB displacement view: the positive octant of
%%% the unit sphere of displacement directions
%%% >> [hPatch, hAx] = Vis.plotPolarisationKey(hAx, nSub)
%%%
%%% INPUT: (opt) 1. hAx  .. target axes/uiaxes (a new classic figure and axes if
%%%                         omitted or empty); cla()'d first
%%%        (opt) 2. nSub .. subdivisions per octant edge (default 40)
%%% OUTPUT: hPatch .. the patch object of the octant
%%%         hAx    .. the axes drawn into
%%%
%%% Each point of the octant is a unit direction p = (p1, p2, p3), p_i >= 0,
%%% coloured p1*red + p2*blue + p3*green at full amplitude (the base colours of
%%% visPolarisationColours) -- exactly the hue Vis.calcPolarisationRgb gives a
%%% displacement in that direction. The octant
%%% is viewed along (1,1,1), oriented as in the UTGS paper's key: u3
%%% (out-of-plane, green) up, u2 (in-plane transverse, blue) lower left, u1
%%% (propagation, red) lower right. Lower amplitude fades the same hue towards
%%% white; the key shows the saturated end only.

    if nargin < 1 || isempty(hAx)
        figure();
        hAx = axes();
    end
    if nargin < 2 || isempty(nSub); nSub = 40; end
    validateattributes(nSub, {'numeric'}, {'scalar', 'integer', '>=', 2}, 'plotPolarisationKey', 'nSub');
    cla(hAx);

    %%% barycentric grid on the flat triangle e1-e2-e3, pushed out onto the sphere
    [iGrid, jGrid] = meshgrid(0:nSub, 0:nSub);
    keep = iGrid + jGrid <= nSub;
    iV = iGrid(keep); jV = jGrid(keep);
    kV = nSub - iV - jV;
    P  = [iV, jV, kV] / nSub;              % rows = (p1, p2, p3) on the flat triangle
    P  = P ./ sqrt(sum(P.^2, 2));          % onto the unit sphere

    %%% triangles of the grid: index lookup (i,j) -> vertex number
    idx = zeros(nSub + 1);
    idx(sub2ind(size(idx), iV + 1, jV + 1)) = 1:numel(iV);
    faces = zeros(nSub^2, 3);
    nF = 0;
    for i = 0:nSub - 1
        for j = 0:nSub - 1 - i
            nF = nF + 1;
            faces(nF, :) = [idx(i+1, j+1), idx(i+2, j+1), idx(i+1, j+2)];
            if i + j < nSub - 1
                nF = nF + 1;
                faces(nF, :) = [idx(i+2, j+1), idx(i+2, j+2), idx(i+1, j+2)];
            end
        end
    end
    faces = faces(1:nF, :);

    [xy, toScreen] = projectOctant(P);
    colours = P * visPolarisationColours();
    hPatch = patch(hAx, 'Faces', faces, 'Vertices', xy, ...
        'FaceVertexCData', colours, 'FaceColor', 'interp', 'EdgeColor', 'none');
    hold(hAx, 'on');

    %%% octant outline: the three great-circle arcs between the axis directions
    t = linspace(0, pi/2, 60)';
    arcs = {[cos(t), sin(t), 0*t], [0*t, cos(t), sin(t)], [sin(t), 0*t, cos(t)]};
    for k = 1:3
        a = toScreen(arcs{k});
        plot(hAx, a(:, 1), a(:, 2), 'k-', 'LineWidth', 0.75);
    end

    %%% corner labels, pushed outward from the centre
    corners = toScreen(eye(3));
    labels  = {'u_1 (propag.)', 'u_2 (in-plane)', 'u_3 (out-of-plane)'};
    cornerColours = 0.8 * visPolarisationColours();   % label text a shade darker
    vAlign  = {'top', 'top', 'bottom'};
    labelOffset = [0, -0.12; 0, -0.12; 0, 0.06];   % below the lower corners, above the top one
    for k = 1:3
        pos = corners(k, :) + labelOffset(k, :);
        text(hAx, pos(1), pos(2), labels{k}, 'Color', cornerColours(k, :), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', vAlign{k}, 'FontSize', 9);
    end

    axis(hAx, 'equal');
    axis(hAx, 'off');
    xlim(hAx, [-1.45, 1.45]);
    ylim(hAx, [-0.85, 1.15]);
    title(hAx, 'polarisation');
end

function [xy, toScreen] = projectOctant(P)
%%% Orthographic projection along (1,1,1): screen-up = the projection of u3,
%%% screen-right = u1 - u2.
    eRight = [1, -1, 0] / sqrt(2);
    eUp    = [-1, -1, 2] / sqrt(6);
    toScreen = @(Q) [Q * eRight.', Q * eUp.'];
    xy = toScreen(P);
end
