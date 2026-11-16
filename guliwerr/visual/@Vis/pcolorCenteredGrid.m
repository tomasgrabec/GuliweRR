function [xEdges, yEdges, CPadded] = pcolorCenteredGrid(x, y, C)
%%% Cell-edge coordinates + padded C so pcolor(...)+shading('flat') draws
%%% each C(i,j) cell CENTERED on (x(j), y(i))
%%% >> [xEdges, yEdges, CPadded] = Vis.pcolorCenteredGrid(x, y, C)
%%%
%%% pcolor's native behaviour treats (x(j), y(i)) as a cell's LOWER-LEFT
%%% corner, not its centre, and 'flat' shading additionally drops the last
%%% row/column of C -- so C(:,j) ends up drawn from x(j) to x(j+1) rather
%%% than centred on x(j) (docs/pitfalls.md). This function converts sample
%%% coordinates to cell edges instead, so the caller's own pcolor(xEdges,
%%% yEdges, CPadded) call draws each sample centred where it belongs.
%%%
%%% INPUT: 1. x .. [1xN] or [Nx1] sample coordinates (need not be evenly
%%%                spaced)
%%%        2. y .. [1xM] or [Mx1] sample coordinates
%%%        3. C .. [MxN] colour data (double or single), C(i,j) belongs to
%%%                (x(j), y(i)); CPadded keeps its class
%%% OUTPUT: xEdges  .. [1x(N+1)] cell-edge x-coordinates
%%%         yEdges  .. [(M+1)x1] cell-edge y-coordinates
%%%         CPadded .. [(M+1)x(N+1)] -- last row/column duplicated from the
%%%                    last real row/column; 'flat' shading never draws it,
%%%                    it exists only to satisfy pcolor's size requirement
%%%
%%% Interior edges are midpoints between consecutive samples; the two open ends
%%% are linearly extrapolated by their nearest interior half-gap. Under
%%% non-uniform spacing this brackets each sample by half its gap to each
%%% neighbour, giving unequal margins where the gaps differ — equal margins on
%%% both sides is not achievable without overlapping cells or leaving gaps.
%%%
%%% x/y need only be MONOTONIC, not increasing: CalcIsofreqAngSweepLayered's
%%% kParRange (= omega0 ./ vAxis) is a naturally decreasing axis.

    x = x(:)';
    y = y(:);
    validateattributes(x, {'double'}, {'vector', 'real', 'finite'}, ...
        'pcolorCenteredGrid', 'x');
    validateattributes(y, {'double'}, {'vector', 'real', 'finite'}, ...
        'pcolorCenteredGrid', 'y');
    validateattributes(C, {'double', 'single'}, {'size', [numel(y), numel(x)]}, ...
        'pcolorCenteredGrid', 'C');
    if ~isMonotonic(x)
        error('Vis:pcolorCenteredGrid:xNotMonotonic', 'x must be monotonic (increasing or decreasing; repeated values are fine).');
    end
    if ~isMonotonic(y)
        error('Vis:pcolorCenteredGrid:yNotMonotonic', 'y must be monotonic (increasing or decreasing; repeated values are fine).');
    end

    xEdges = localEdges(x);
    yEdges = localEdges(y.').';
    CPadded = [C, C(:, end); C(end, :), C(end, end)];
end

function e = localEdges(v)
    if isscalar(v)
        e = [v - 0.5, v + 0.5];   % degenerate single-sample axis -- arbitrary unit width
        return
    end
    mid = (v(1:end-1) + v(2:end)) / 2;
    e = [2*v(1) - mid(1), mid, 2*v(end) - mid(end)];
end

function tf = isMonotonic(v)
%%% non-strict on purpose: tolerates repeated-equal values, e.g. two identical
%%% radius samples, and rejects only an actual direction reversal
    d = diff(v);
    tf = isscalar(v) || all(d >= 0) || all(d <= 0);
end
