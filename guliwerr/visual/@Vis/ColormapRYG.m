function cMap = ColormapRYG()
%%% Red-Yellow-Green sequential colormap, for a confidence/quality scalar in [0,1]
%%% where the VALUE, not just the ordering, has a meaning: 0 = bad, 1 = good.
%%% >> cMap = Vis.ColormapRYG()
%%%
%%% OUTPUT: cMap .. 101x3 RGB colormap, red at 0, yellow at 0.5, green at 1
%%%
%%% Deliberately not a perceptually-uniform sequential map like parula: parula puts
%%% its brightest, most eye-catching colour (yellow) at the HIGH end, which reads as
%%% "flagged"/"warning" for a quality score where high = good, not as reassurance.
%%% Traffic-light semantics instead -- e.g. ModeAssociation's per-point MAC value,
%%% where red should mean "look at this" and green "this tracked cleanly".

    n = 101;
    t = linspace(0, 1, n)';
    lowHalf = t <= 0.5;

    r = zeros(n, 1); g = zeros(n, 1); b = zeros(n, 1);

    %%% red (0.85,0.10,0.10) -> yellow (0.95,0.85,0.10)
    frac1 = 2 * t(lowHalf);
    r(lowHalf) = 0.85 + (0.95-0.85) * frac1;
    g(lowHalf) = 0.10 + (0.85-0.10) * frac1;
    b(lowHalf) = 0.10;

    %%% yellow (0.95,0.85,0.10) -> green (0.10,0.60,0.15)
    frac2 = 2 * (t(~lowHalf) - 0.5);
    r(~lowHalf) = 0.95 + (0.10-0.95) * frac2;
    g(~lowHalf) = 0.85 + (0.60-0.85) * frac2;
    b(~lowHalf) = 0.10 + (0.15-0.10) * frac2;

    cMap = [r, g, b];
end
