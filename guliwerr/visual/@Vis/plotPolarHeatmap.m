function [hOut, hAx] = plotPolarHeatmap(angles, radius, C, hAx)
%%% Draw a 2D angle x radius heatmap in polar form, on ordinary Cartesian axes
%%% >> [hOut, hAx] = plotPolarHeatmap(angles, radius, C, hAx)
%%%
%%% INPUT: 1. angles .. [nAngles x 1] angles [deg] -- plotted exactly as
%%%                     given, NOT mirrored/extended to a full circle. Pass
%%%                     a 0-360 deg array directly if a full disc is wanted;
%%%                     a 0-180 deg array draws a half-disc wedge.
%%%        2. radius  .. [nRadius x 1] radial axis (e.g. velocity, kParRange)
%%%        3. C       .. [nRadius x nAngles] colour data -- already scaled
%%%                     (log-normalised or not); this function does no
%%%                     scaling itself, it only draws
%%%  (opt) 4. hAx      .. target axes (created fresh if omitted) -- MUST be
%%%                     plain axes(), not polaraxes(): the latter cannot host
%%%                     pcolor/surf in any current MATLAB release (see
%%%                     docs/reference/architecture/visual.md).
%%% OUTPUT: hOut .. pcolor handle
%%%         hAx  .. axes used (axis 'equal' and 'off' -- the hand-drawn grid
%%%                 replaces the default x/y ticks)
%%%
%%% Colour scaling (log/linear), colormap, clim and colorbar are left to the
%%% CALLER: this replaces only the pcolor + axis-decoration half of a Cartesian
%%% call, so both display modes get identical colour treatment.
%%%
%%% Generic over any (angle, radius, colour) map, not Green-function-specific.

    validateattributes(angles, {'double'}, {'vector', 'real', 'finite'}, ...
        'plotPolarHeatmap', 'angles');
    validateattributes(radius, {'double'}, {'vector', 'real', 'finite', 'nonnegative'}, ...
        'plotPolarHeatmap', 'radius');
    angles = angles(:)';
    radius = radius(:);
    validateattributes(C, {'double'}, {'size', [numel(radius), numel(angles)]}, ...
        'plotPolarHeatmap', 'C');

    if nargin < 4 || isempty(hAx)
        figure();
        hAx = axes();
    end
    hold(hAx, 'on');

    [angEdges, radEdges, CPadded] = Vis.pcolorCenteredGrid(angles, radius, C);
    radEdges = max(radEdges, 0);   % an extrapolated inner edge can dip below the
                                    % innermost sample; a negative radius is meaningless
                                    % here and would flip that ring to the opposite
                                    % angle under pol2cart
    [Theta, R] = meshgrid(deg2rad(angEdges), radEdges);
    [X, Y]     = pol2cart(Theta, R);
    hOut = pcolor(hAx, X, Y, CPadded);
    shading(hAx, 'flat');
    axis(hAx, 'equal');
    axis(hAx, 'off');

    drawPolarGrid(hAx, angles, radius);
end

function drawPolarGrid(hAx, angles, radius)
%%% Hand-drawn radial rings and angular spokes, labelled -- the decoration a real
%%% polaraxes() would give for free, needed here because pcolor requires a plain
%%% axes().
    gridColor = [0.55, 0.55, 0.55];
    rMin = min(radius);
    rMax = max(radius);
    if rMax <= rMin
        return   % degenerate single-radius input -- nothing sensible to draw
    end

    angMin = min(angles);
    angMax = max(angles);
    thArc  = linspace(deg2rad(angMin), deg2rad(angMax), 200);   % the data's actual angular
                                                                   % extent, so a 0-180 deg input
                                                                   % draws half-rings and the grid
                                                                   % never implies absent data

    nRings    = 4;
    ringRadii = linspace(rMin, rMax, nRings + 1);
    ringRadii = ringRadii(2:end);   % skip r=0 (a point, no ring to draw)
    for r = ringRadii
        plot(hAx, r * cos(thArc), r * sin(thArc), ':', ...
            'Color', gridColor, 'LineWidth', 0.75, 'HandleVisibility', 'off');
        text(hAx, r * cos(thArc(end)), r * sin(thArc(end)), sprintf(' %.3g', r), ...
            'Color', gridColor, 'FontSize', 8, 'VerticalAlignment', 'bottom');
    end

    angStep  = 30;   % deg
    angTicks = angMin:angStep:angMax;
    if isempty(angTicks) || angTicks(end) ~= angMax
        angTicks = [angTicks, angMax];   %#ok<AGROW>
    end
    for a = angTicks
        plot(hAx, [0, rMax * cosd(a)], [0, rMax * sind(a)], ':', ...
            'Color', gridColor, 'LineWidth', 0.75, 'HandleVisibility', 'off');
        text(hAx, 1.06 * rMax * cosd(a), 1.06 * rMax * sind(a), sprintf('%g\xB0', a), ...
            'Color', gridColor, 'FontSize', 8, 'HorizontalAlignment', 'center');
    end
end
