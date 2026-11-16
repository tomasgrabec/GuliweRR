function [hOut, hAx] = plotDispersionMap(sweepResult, omegaIdx, scaleMode, climRange, plotType, hAx)
%%% Plot a 2D (or 1D) dispersion map -- angle x wavelength/kParRange/omega0
%%% >> [hOut, hAx] = plotDispersionMap(sweepResult, omegaIdx, scaleMode, climRange, plotType, hAx)
%%%
%%% INPUT: 1. sweepResult .. one of four different result classes, dispatched
%%%                          by isa() since their underlying data shapes are
%%%                          genuinely different, not just differently-named
%%%                          properties on the same shape (unlike
%%%                          plotAngularMap/plotGreenMap's HS-vs-Layered
%%%                          duck-typing):
%%%           - CalcWavevecAngSweep2Layered (isoWavevec 2D: angle x wavelength)
%%%             -> scatter of classified SAW points (angle, wavelength),
%%%                coloured by LOCAL RANK within each point's own cluster.
%%%                sawVelMap2D{iw,ia} is already velocity-sorted ascending, so
%%%                rank is position in that array and no new sort is needed.
%%%                These are ranks, NOT tracked branches — for connected branch
%%%                lines use Vis.plotDispersionCurves. No Green magnitude here,
%%%                so scaleMode/climRange/plotType are ignored: always a
%%%                Cartesian scatter.
%%%           - CalcFreqAngSweep2Layered (isoFreq 2D: angle x omega0, each
%%%             omega0 its own angle x kParRange Green |G| heatmap slice, not
%%%             eagerly aggregated) -> plots ONE selected omega0 slice, see
%%%             omegaIdx
%%%           - CalcIsofreqAngSweepLayered (isoFreq, one omega0: angle x
%%%             kParRange) -> Green |G| heatmap, the same ColMapGrayOrange
%%%             treatment plotGreenMap uses
%%%           - CalcIsofreqSingleLayered (isoFreq, one angle, one omega0:
%%%             kParRange only) -> 1D line, absG vs kParRange. No angle axis at
%%%             all, so plotType is ignored: always a Cartesian line.
%%%  (opt) 2. omegaIdx    .. which omega0 slice to plot -- ONLY used when
%%%                          sweepResult is a CalcFreqAngSweep2Layered
%%%                          (default 1); ignored otherwise
%%%  (opt) 3. scaleMode    .. 'log' (default) | 'linear', as in plotGreenMap.
%%%                          Applies to the three Green-magnitude cases;
%%%                          ignored for CalcWavevecAngSweep2Layered, which has
%%%                          no Green magnitude to scale. 'log' plots
%%%                          log10(|G_33| / max|G_33|), normalised to this
%%%                          plot's own peak so the display range is
%%%                          material- and units-independent; 'linear' gives
%%%                          plain |G_33|, which buries the peak's sharpness
%%%                          against a much broader background
%%%  (opt) 4. climRange    .. 2-element [lo hi]. Colour-axis limits for the 2D
%%%                          heatmap cases (log10-normalised units when
%%%                          scaleMode='log'), default [-5, -2]; ylim for the 1D
%%%                          CalcIsofreqSingleLayered line, default [-5, 0] (the
%%%                          normalised peak sits at 0, so the heatmap default
%%%                          would cut it off). Ignored for 'linear' unless
%%%                          given explicitly
%%%  (opt) 5. plotType     .. 'cartesian' (default, angle x kParRange axes) |
%%%                          'polar' (angle as azimuth, kParRange as radius).
%%%                          Affects the two 2D heatmap cases only, delegating
%%%                          to plotPolarHeatmap.m. sweepResult.angles is
%%%                          plotted exactly as given, not mirrored to a full
%%%                          circle
%%%  (opt) 6. hAx          .. target axes (created fresh if omitted)
%%% OUTPUT: hOut .. scatter/pcolor/line handle, whichever case applied
%%%         hAx  .. axes used
%%%
%%% kParRange is plotted as given, not converted to a velocity axis.
%%%
%%% The climRange default here ([-5,-2]) is tuned independently of
%%% plotGreenMap's ([-4,-1]) — they cover different scan types.

    if nargin < 3 || isempty(scaleMode); scaleMode = 'log'; end
    visValidateChoice(scaleMode, {'log', 'linear'}, 'plotDispersionMap', 'scaleMode');
    climGiven = nargin >= 4 && ~isempty(climRange);
    if ~climGiven; climRange = [-5, -2]; end
    if nargin < 5 || isempty(plotType); plotType = 'cartesian'; end
    visValidateChoice(plotType, {'cartesian', 'polar'}, 'plotDispersionMap', 'plotType');
    if nargin < 6 || isempty(hAx)
        figure();
        hAx = axes();
    end
    hold(hAx, 'on');

    if isa(sweepResult, 'CalcWavevecAngSweep2Layered')
        [nWav, nAngles] = size(sweepResult.sawVelMap2D);
        allAng = []; allWav = []; allRank = [];
        for iw = 1:nWav
            for ia = 1:nAngles
                v = sweepResult.sawVelMap2D{iw, ia};
                n = numel(v);
                if n == 0; continue; end
                allAng  = [allAng;  repmat(sweepResult.angles(ia), n, 1)];      %#ok<AGROW>
                allWav  = [allWav;  repmat(sweepResult.wavelengths(iw), n, 1)]; %#ok<AGROW>
                allRank = [allRank; (1:n)'];                                    %#ok<AGROW>
            end
        end
        hOut = scatter(hAx, allAng, allWav, 36, allRank, 'filled');
        colormap(hAx, parula);
        colorbar(hAx);
        xlabel(hAx, 'angle [deg]');
        ylabel(hAx, 'wavelength');
        title(hAx, 'colour = local branch rank (no branch tracking yet)');

    elseif isa(sweepResult, 'CalcFreqAngSweep2Layered')
        if nargin < 2 || isempty(omegaIdx); omegaIdx = 1; end
        validateattributes(omegaIdx, {'double'}, ...
            {'scalar', 'integer', 'positive', '<=', numel(sweepResult.omega0List)}, ...
            'plotDispersionMap', 'omegaIdx');
        slice = sweepResult.sweeps{omegaIdx};
        hOut = plotIsofreqAngSweepHeatmap(slice, hAx, scaleMode, climRange, plotType);
        title(hAx, sprintf('omega0 = %.4g', sweepResult.omega0List(omegaIdx)));

    elseif isa(sweepResult, 'CalcIsofreqAngSweepLayered')
        hOut = plotIsofreqAngSweepHeatmap(sweepResult, hAx, scaleMode, climRange, plotType);

    elseif isa(sweepResult, 'CalcIsofreqSingleLayered')
        [plotData, yLabelStr] = Vis.scaleGreenData(sweepResult.absG, scaleMode, 'G_{33}');
        hOut = plot(hAx, sweepResult.kParRange, plotData);
        if strcmp(scaleMode, 'log')
            if climGiven; ylim(hAx, climRange); else; ylim(hAx, [-5, 0]); end
        end
        %%% fixed, not auto: MATLAB's auto x-limits skip data clipped by a fixed ylim
        xlim(hAx, [min(sweepResult.kParRange), max(sweepResult.kParRange)]);
        xlabel(hAx, 'kParRange');
        ylabel(hAx, yLabelStr);

    else
        error('Vis:plotDispersionMap:badInput', ...
            ['sweepResult must be a CalcWavevecAngSweep2Layered, CalcFreqAngSweep2Layered, ' ...
             'CalcIsofreqAngSweepLayered, or CalcIsofreqSingleLayered (got %s).'], class(sweepResult));
    end
end

function hOut = plotIsofreqAngSweepHeatmap(slice, hAx, scaleMode, climRange, plotType)
%%% Shared angle x kParRange Green |G| heatmap, used both directly
%%% (CalcIsofreqAngSweepLayered) and for one selected slice of
%%% CalcFreqAngSweep2Layered
    [plotData, cbLabel] = Vis.scaleGreenData(slice.absGMap, scaleMode, 'G_{33}');

    if strcmp(plotType, 'cartesian')
        hOut = visPcolorCentered(hAx, slice.angles, slice.kParRange, plotData);
        xlabel(hAx, 'angle [deg]');
        ylabel(hAx, 'kParRange');
    else
        hOut = Vis.plotPolarHeatmap(slice.angles, slice.kParRange, plotData, hAx);
    end

    visGreenColourAxis(hAx, scaleMode, climRange, cbLabel);
end
