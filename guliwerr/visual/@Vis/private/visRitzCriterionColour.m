function visRitzCriterionColour(hAx, hCb, info)
%%% Colour scale of a Ritz-mode criterion map: unflipped `gray`, with the freqDisp veto marked
%%% >> visRitzCriterionColour(hAx, hCb, info)
%%%
%%% @Vis-private helper of Vis.plotRitzModeMap.
%%%
%%% INPUT: 1. hAx  .. the map axes
%%%        2. hCb  .. its colour bar
%%%        3. info .. visRitzCriterion's info (label, vetoLevel)
%%%
%%% Without a veto level: `gray` (black = low, white = high), auto limits. With one
%%% (freqDisp on a half-space): limits fixed to [0 1], the freqDisp range, so the
%%% veto level sits at a fixed place on the bar; values below it, which the Ritz rule
%%% vetoes, are shaded dark-to-light red, the rest grey; the bar gets a tick at the
%%% level labelled "veto".
    N = 256;
    if isnan(info.vetoLevel)
        colormap(hAx, gray(N));
        return
    end
    nVeto = round(N * info.vetoLevel);
    ramp  = linspace(0.25, 1, nVeto)';
    vetoRows = [ramp, 0.25 * ramp, 0.25 * ramp];
    keepRows = gray(N - nVeto);
    colormap(hAx, [vetoRows; keepRows]);
    clim(hAx, [0 1]);
    ticks = unique([0, info.vetoLevel, 0.5, 1]);
    labels = arrayfun(@(t) sprintf('%.1f', t), ticks, 'UniformOutput', false);
    labels{ticks == info.vetoLevel} = sprintf('%.1f veto', info.vetoLevel);
    hCb.Ticks = ticks;
    hCb.TickLabels = labels;
end
