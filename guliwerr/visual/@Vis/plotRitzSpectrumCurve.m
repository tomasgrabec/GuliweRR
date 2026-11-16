function [hOut, hAx, spec] = plotRitzSpectrumCurve(singleResult, scaleMode, hAx, showCands)
%%% Plot the Ritz spectrum |S(v)| of one propagation angle -- the Ritz-side counterpart
%%% of plotGreenCurve's |G13_eff|, rebuilt from the stored modal numerators (no solve)
%%% >> [hOut, hAx, spec] = plotRitzSpectrumCurve(singleResult)
%%% >> [hOut, hAx, spec] = plotRitzSpectrumCurve(singleResult, scaleMode, hAx, showCands)
%%%
%%% INPUT: 1. singleResult .. CalcSingleHS or CalcSingleLayered, calcAll()'d (either
%%%                          calcMode: the Ritz candidates are stored on both)
%%%  (opt) 2. scaleMode    .. 'log' (default) | 'linear' -- the normalisation of
%%%                          plotGreenCurve (log10, normalised to this curve's peak)
%%%  (opt) 3. hAx          .. target axes (created fresh if omitted)
%%%  (opt) 4. showCands    .. false (default, as plotRitzSpectrumMap) | true: mark the
%%%                          Ritz candidates -- down-triangles at the modal-sum peaks,
%%%                          circles at their supporting modes
%%% OUTPUT: hOut .. line handle
%%%         hAx  .. axes used
%%%         spec .. struct: velAxis, absS (|S| on velAxis), cands (ritzScoreCands), so a
%%%                 caller can reuse the spectrum without rebuilding it
%%%
%%% S(v) = sum_n u3_n conj(u1_n) / (v_n^2 - v^2 (1 + i ETA)) (RitzScore.spectrum), the
%%% G13-like modal sum the Ritz score is built on, on the result's Ritz score axis
%%% (ritzScoreVelAxis). Errors Vis:plotRitzSpectrumCurve:noRitzSpectrum without stored
%%% numerators.
    if nargin < 2 || isempty(scaleMode); scaleMode = 'log'; end
    visValidateChoice(scaleMode, {'log', 'linear'}, 'plotRitzSpectrumCurve', 'scaleMode');
    if nargin < 4 || isempty(showCands); showCands = false; end
    if ~visHasRitzSpectrum(singleResult)
        error('Vis:plotRitzSpectrumCurve:noRitzSpectrum', ...
            'This result stores no Ritz modal numerators (calcAll() stores them).');
    end
    if nargin < 3 || isempty(hAx)
        figure();
        hAx = axes();
    end

    [absS, velAxis] = visRitzSpectrumData(singleResult);
    [plotData, yLabel] = Vis.scaleGreenData(absS, scaleMode, 'S');

    hOut = plot(hAx, velAxis, plotData, '-', 'DisplayName', '|S|');
    cands = singleResult.ritzScoreCands;
    if showCands && ~isempty(cands)
        hold(hAx, 'on');
        yPk = interp1(velAxis, plotData, [cands.vel], 'linear', NaN);
        scatter(hAx, [cands.vel], yPk, 36, 'r', 'v', 'filled', 'HandleVisibility', 'off');
        yMode = interp1(velAxis, plotData, [cands.modeVel], 'linear', NaN);
        scatter(hAx, [cands.modeVel], yMode, 30, 'k', 'o', 'HandleVisibility', 'off');
        hold(hAx, 'off');
    end
    xlabel(hAx, 'velocity [km/s]');
    ylabel(hAx, yLabel);
    title(hAx, sprintf('Ritz spectrum (angle = %.1f deg)', singleResult.angle));
    spec = struct('velAxis', velAxis, 'absS', absS, 'cands', cands);
end
