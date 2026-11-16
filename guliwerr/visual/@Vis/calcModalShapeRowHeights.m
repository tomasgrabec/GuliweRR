function frac = calcModalShapeRowHeights(relDepths, minFrac)
%%% Non-proportional per-layer height fractions for the modal-shape depth stack
%%% >> frac = Vis.calcModalShapeRowHeights(relDepths, minFrac)
%%%
%%% INPUT: 1. relDepths .. [nLay x 1] every layer's relDepth, uncropped --
%%%                        layer 1 = substrate, layer nLay = outermost film
%%%  (opt) 2. minFrac    .. floor fraction of the total displayed height
%%%                         every layer is guaranteed, before the remainder
%%%                         is split proportionally to true relDepth
%%%                         (default 0.3). Auto-clamped down if nLay*minFrac
%%%                         would reach/exceed 1 (many thin films).
%%% OUTPUT: frac .. [nLay x 1], sums to 1
%%%
%%% Every layer is guaranteed a floor share up front so a thin film stays visible
%%% next to a much deeper substrate.

    if nargin < 2 || isempty(minFrac); minFrac = 0.3; end
    validateattributes(relDepths, {'double'}, {'vector', 'positive'}, ...
        'calcModalShapeRowHeights', 'relDepths');
    validateattributes(minFrac, {'double'}, {'scalar', '>', 0, '<', 1}, ...
        'calcModalShapeRowHeights', 'minFrac');
    nLay = numel(relDepths);

    if nLay * minFrac >= 1
        minFrac = 0.9 / nLay;
    end

    remainPool = 1 - nLay * minFrac;
    %%% shape follows the caller's own relDepths orientation -- a ROW in
    %%% plotModalShape, whose layHeightFrac is built via [0, cumsum(...)]
    frac = minFrac + remainPool * (relDepths / sum(relDepths));
end
