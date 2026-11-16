function [hOut, hAx, hCb] = plotRitzSpectrumMap(sweepResult, scaleMode, climRange, plotType, hAx, showCands)
%%% Angle x velocity map of the Ritz spectrum |S(v)| -- the Ritz-side counterpart of
%%% plotGreenMap's |G13_eff| map, rebuilt from the stored modal numerators (no solve)
%%% >> [hOut, hAx, hCb] = plotRitzSpectrumMap(sweepResult)
%%% >> [hOut, hAx, hCb] = plotRitzSpectrumMap(sweepResult, scaleMode, climRange, plotType, hAx, showCands)
%%%
%%% INPUT: 1. sweepResult .. CalcAngSweepHS or CalcAngSweepLayered, calc()'d (either
%%%                         calcMode)
%%%  (opt) 2. scaleMode   .. 'log' (default) | 'linear', as plotGreenMap
%%%  (opt) 3. climRange   .. colour limits in the log scale (default [-4, -1], as
%%%                         plotGreenMap)
%%%  (opt) 4. plotType    .. 'cartesian' (default and only value; the slot keeps the
%%%                         argument order of plotGreenMap, so hAx is 5th in both)
%%%  (opt) 5. hAx         .. target axes (created fresh if omitted)
%%%  (opt) 6. showCands   .. false (default, as plotRitzSpectrumCurve) | true: mark
%%%                         every angle's Ritz candidates at their supporting modes'
%%%                         velocities
%%% OUTPUT: hOut .. the map's surface handle
%%%         hAx  .. axes used
%%%         hCb  .. colour bar
%%%
%%% Each angle's S(v) (RitzScore.spectrum) is evaluated on one common axis (the union
%%% of the angles' Ritz score axes, RitzScore.VEL_STEP), normalised to the map's maximum
%%% like plotGreenMap. The default background of a 'ritzOnly' result, which has no
%%% Green map.
    if nargin < 2 || isempty(scaleMode); scaleMode = 'log'; end
    visValidateChoice(scaleMode, {'log', 'linear'}, 'plotRitzSpectrumMap', 'scaleMode');
    if nargin < 3 || isempty(climRange); climRange = [-4, -1]; end
    if nargin < 4 || isempty(plotType); plotType = 'cartesian'; end
    visValidateChoice(plotType, {'cartesian'}, 'plotRitzSpectrumMap', 'plotType');
    if nargin < 6 || isempty(showCands); showCands = false; end
    cs = sweepResult.cs;
    if ~visHasRitzSpectrum(sweepResult)
        error('Vis:plotRitzSpectrumMap:noRitzSpectrum', ...
            'This sweep stores no Ritz modal numerators (calc() stores them).');
    end
    if nargin < 5 || isempty(hAx)
        figure();
        hAx = axes();
    end

    vLo = min(cellfun(@(c) c.ritzScoreVelAxis(1), cs));
    vHi = max(cellfun(@(c) c.ritzScoreVelAxis(end), cs));
    velAxis = (vLo:RitzScore.VEL_STEP:vHi)';
    nA = numel(cs);
    absMap = zeros(numel(velAxis), nA);
    for ia = 1:nA
        absMap(:, ia) = visRitzSpectrumData(cs{ia}, velAxis);
    end
    [plotData, cbLabel] = Vis.scaleGreenData(absMap, scaleMode, 'S');

    hOut = visPcolorCentered(hAx, sweepResult.angles, velAxis, plotData);
    hOut.DisplayName = '|S|';
    xlabel(hAx, 'angle [deg]');
    ylabel(hAx, 'velocity [km/s]');
    title(hAx, 'Ritz spectrum map', 'Interpreter', 'none');
    hCb = visGreenColourAxis(hAx, scaleMode, climRange, cbLabel);
    if showCands
        hold(hAx, 'on');
        for ia = 1:nA
            c = cs{ia}.ritzScoreCands;
            if isempty(c); continue; end
            plot(hAx, repmat(sweepResult.angles(ia), 1, numel(c)), [c.modeVel], 'r.', ...
                'HandleVisibility', 'off');
        end
        hold(hAx, 'off');
    end
end
