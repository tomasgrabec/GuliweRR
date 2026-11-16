function [plotData, dataLabel] = scaleGreenData(absG, scaleMode, componentLabel)
%%% Apply the shared Green's-magnitude display scaling, with its axis label
%%% >> [plotData, dataLabel] = Vis.scaleGreenData(absG, scaleMode, componentLabel)
%%%
%%% INPUT: 1. absG           .. |G| data of any shape (curve, map, or a
%%%                             partially-NaN interpolated map -- see below)
%%%        2. scaleMode      .. 'log' | 'linear'
%%%        3. componentLabel .. TeX label for the quantity absG holds, used
%%%                             verbatim in dataLabel, e.g. 'G_{13}^{eff}' |
%%%                             'G_{33}' | 'S' (the Ritz spectrum); see
%%%                             Vis.greenComponentData. Required: the label must
%%%                             name the data the caller passes
%%% OUTPUT: plotData  .. the values to plot/colour by, double (the Calc classes
%%%                      store |G| in single precision; plotting and interpolating
%%%                      are done in double)
%%%         dataLabel .. matching colorbar/axis label string
%%%
%%% 'log' plots log10(|G| / max|G|) -- normalised to THIS data's own
%%% peak first, so the display range is material- and units-independent: it
%%% shows how sharp the SAW peak is relative to its own background rather
%%% than its absolute magnitude. That is the sensible default everywhere it
%%% is offered: a linear scale buries the Rayleigh peak's actual sharpness
%%% against a much broader background. Zeros are floored to eps before
%%% log10 so CData stays real and finite; out-of-range values are clipped for
%%% DISPLAY only, via the caller's clim()/ylim().
%%%
%%% NaN handling: NaN entries (a real input -- plotDispersionCurves
%%% interpolates per-wavelength Green columns onto one shared velocity axis
%%% and leaves NaN outside each column's own range) are excluded from the
%%% peak, and restored explicitly afterwards. Without that restore,
%%% max(NaN, eps) resolves to eps in MATLAB's two-array max (NaN is treated
%%% as absent whenever the other operand isn't also NaN), which would paint
%%% every out-of-range cell at the floor colour instead of leaving it blank.
%%%
%%% Public rather than @Vis-private because GUI code outside visual/ needs the
%%% same normalisation to place its own overlays on these plots — as does
%%% plotClassificationCurve's y-projection, which must land exactly on the curve
%%% plotGreenCurve drew.

    narginchk(3, 3);
    absG = double(absG);

    if strcmp(scaleMode, 'log')
        normVal   = absG ./ max(absG(:));
        plotData  = log10(max(normVal, eps));
        plotData(isnan(normVal)) = NaN;
        dataLabel = sprintf('log_{10}(|%s| / max|%s|)', componentLabel, componentLabel);
    else
        plotData  = absG;
        dataLabel = sprintf('|%s|', componentLabel);
    end
end
