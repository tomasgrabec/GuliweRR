function [hPoints, hAx, hBg, hLim] = plotClassificationCurve(singleResult, colorBy, showGreenBackground, refLevel, hAx, component)
%%% Comparative mode-classification diagnostic for a SINGLE propagation
%%% angle: classified points marked on the Green's-function curve, with a
%%% Christoffel limiting-velocity reference line -- 1D analogue of
%%% plotClassificationDiagnostic.m for CalcSingleHS/CalcSingleLayered
%%% results (no angle axis to sweep).
%%% >> [hPoints, hAx, hBg, hLim] = plotClassificationCurve(singleResult, colorBy, showGreenBackground, refLevel, hAx, component)
%%%
%%% INPUT: 1. singleResult         .. CalcSingleHS or CalcSingleLayered result
%%%                                  object, calcAll()'d. Reads
%%%                                  .classification.significantModes directly
%%%                                  and never recomputes a parallel score
%%%  (opt) 2. colorBy             .. 'type' (default) | 'score' | 'confidence',
%%%                                  as in plotClassificationDiagnostic
%%%  (opt) 3. showGreenBackground .. draw the response curve first (default true;
%%%                                  skipped with a warning if the chosen curve's
%%%                                  data is missing, e.g. a Green component on a
%%%                                  'ritzOnly' result)
%%%  (opt) 4. refLevel            .. how much classification detail to show
%%%                                  (default 'full') -- the same 3 tiers
%%%                                  plotClassificationDiagnostic documents:
%%%                                  'sawOnly' | 'sawLimit' | 'full'
%%%  (opt) 5. hAx                 .. target axes (created fresh if omitted)
%%%  (opt) 6. component           .. response curve drawn and projected onto:
%%%                                  'G13eff' (|G13_eff| -- the component the Green
%%%                                  score reads, so points sit on the curve they
%%%                                  were decided on; the default when the result
%%%                                  has Green data) | any 'GIJ'/'GIJeff',
%%%                                  recomputed on request (see
%%%                                  Vis.greenComponentData) | 'ritzSpectrum'
%%%                                  (plotRitzSpectrumCurve, |S(v)|; the default on
%%%                                  a 'ritzOnly' result)
%%% OUTPUT: hPoints .. colorBy='type': array of scatter handles, one per type
%%%                    present (legend shown); 'score'/'confidence': a single
%%%                    scatter handle (its own colour bar: UserData.colorbar)
%%%         hAx     .. axes used
%%%         hBg     .. the background curve's line handle, or [] if not drawn
%%%         hLim    .. struct('long', <xline handle or []>, 'shear', <array of
%%%                    xline handles or []>). Only ONE is ever populated: long
%%%                    for HS, shear for layered
%%%
%%% Y-AXIS: velocity occupies the x-axis here (the swept-away dimension is
%%% angle, not velocity), so each classified mode is marked at its own velocity
%%% and projected onto the log-normalised Green curve's y-value there by
%%% interpolation, reading as a conventional peak-picking plot. Falls back to a
%%% fixed y=0 baseline when the chosen curve's data is missing, so classified points
%%% stay visible rather than being suppressed for want of a curve to project onto.
%%%
%%% MARKER SIZE and GRADE: as plotClassificationDiagnostic -- size by score over the
%%% method's range (Vis.classificationScoreScale), 'weak' entries open markers.
%%%
%%% 'score'/'confidence' colouring is `gray` unflipped — black low, white high —
%%% on the points' own range and colour bar, with a dark outline, as in
%%% plotClassificationDiagnostic.

    if nargin < 2 || isempty(colorBy); colorBy = 'type'; end
    visValidateChoice(colorBy, {'type', 'score', 'confidence'}, ...
        'plotClassificationCurve', 'colorBy');
    if nargin < 3 || isempty(showGreenBackground); showGreenBackground = true; end
    if nargin < 4 || isempty(refLevel);            refLevel            = 'full'; end
    visValidateChoice(refLevel, {'sawOnly', 'sawLimit', 'full'}, ...
        'plotClassificationCurve', 'refLevel');

    if nargin < 6 || isempty(component)
        if isempty(singleResult.greenVelAxis) && visHasRitzSpectrum(singleResult)
            component = 'ritzSpectrum';   % 'ritzOnly': the Ritz-side response
        else
            component = 'G13eff';
        end
    end
    if ~strcmp(component, 'ritzSpectrum')
        visValidateGreenComponent(component, 'plotClassificationCurve');
    end

    if nargin < 5 || isempty(hAx)
        figure();
        hAx = axes();
    end
    hold(hAx, 'on');

    %%% the chosen response curve's data present (Green: calcMode 'full')
    switch component
        case 'spot'
            hasGreen = isprop(singleResult, 'spotGreen') && ~isempty(singleResult.spotGreen);
        case 'ritzSpectrum'
            hasGreen = visHasRitzSpectrum(singleResult);
        otherwise
            hasGreen = ~isempty(singleResult.greenVelAxis);
    end

    hBg = [];
    if showGreenBackground
        if ~hasGreen
            warning('Vis:plotClassificationCurve:noGreenBackground', '%s -- skipping the background.', ...
                visMissingBackgroundReason(component, false));
        elseif strcmp(component, 'ritzSpectrum')
            hBg = Vis.plotRitzSpectrumCurve(singleResult, [], hAx, false);
            hold(hAx, 'on');
        else
            hBg = Vis.plotGreenCurve(singleResult, [], hAx, component);
        end
    end

    %%% flatten significantModes into flat per-point arrays, filtered per refLevel
    typeAllowed = visClassificationTypesAllowed(refLevel);
    %%% one "angle" cell entry and no angle axis to place points on, so the shared
    %%% flattener's angle column is unused here (NaN)
    pts = visFlattenSignificantModes({singleResult.classification.significantModes}, ...
        NaN, typeAllowed);

    %%% y-position: interpolate onto the same log-normalised curve plotGreenCurve
    %%% draws, recomputed here rather than read off hBg.YData, which is [] whenever
    %%% showGreenBackground=false. Fixed y=0 baseline when Green data is absent.
    if hasGreen && strcmp(component, 'ritzSpectrum')
        [absCurve, velAxis] = visRitzSpectrumData(singleResult);
        [curveYData, yLabelStr] = Vis.scaleGreenData(absCurve, 'log', 'S');
        allYPos = interp1(velAxis, curveYData, pts.vel, 'linear', 'extrap');
    elseif hasGreen
        [absCurve, componentLabel, velAxis] = Vis.greenComponentData(singleResult, component);
        [curveYData, yLabelStr] = Vis.scaleGreenData(absCurve, 'log', componentLabel);
        allYPos    = interp1(velAxis, curveYData, pts.vel, 'linear', 'extrap');
    else
        allYPos = zeros(size(pts.vel));
    end

    %%% marker size by score, weak entries open -- see docblock. Larger overall than
    %%% plotClassificationDiagnostic's sizes: a single-angle curve carries far fewer
    %%% points and can afford them.
    MIN_SZ = 30; MAX_SZ = 140; DEFAULT_SZ = 60; LIMITING_SZ = 26;
    isLayered = visIsLayered(singleResult);
    sizeVals = visClassificationSizes(pts, isLayered, ...
        MIN_SZ, MAX_SZ, DEFAULT_SZ, LIMITING_SZ);
    isWeak = strcmp(pts.grade, 'weak');

    switch colorBy
        case 'type'
            %%% limiting-velocity points get a smaller FIXED size (LIMITING_SZ), so a
            %%% limiting REFERENCE marker cannot hide a real mode point
            hPoints = visClassificationScatterByType(hAx, pts.vel, allYPos, sizeVals, pts.type, isWeak);

        otherwise   % 'score' or 'confidence': black=low, white=high -- see docblock
            hPoints = visClassificationScatterByScore(hAx, pts.vel, allYPos, sizeVals, pts, colorBy, isWeak);
    end

    %%% ONE reference line only, picked by family: HS gets longitudinal, layered
    %%% gets substrate shear
    hLim = struct('long', [], 'shear', []);
    if ~strcmp(refLevel, 'sawOnly')
        if ~isLayered
            longLimVel = singleResult.longLimVel;
            if ~isempty(longLimVel) && ~isnan(longLimVel)
                hLim.long = xline(hAx, longLimVel, '-', 'Color', [0.9290 0.6940 0.1250], ...
                    'LineWidth', 1.5, 'DisplayName', 'longitudinal limit');
            end
        else
            shearVels = visPropOrEmpty(singleResult, {'subShearLimVels'});
            if ~isempty(shearVels)
                hShear = gobjects(1, numel(shearVels));
                for k = 1:numel(shearVels)
                    %%% one legend entry for the whole family: 'DisplayName' on the
                    %%% first line only, 'HandleVisibility' off on the rest, or each
                    %%% legends separately as a generic 'dataN'
                    if k == 1
                        hShear(k) = xline(hAx, shearVels(k), '--', 'Color', [0.4660 0.6740 0.1880], ...
                            'DisplayName', 'shear limit');
                    else
                        hShear(k) = xline(hAx, shearVels(k), '--', 'Color', [0.4660 0.6740 0.1880], ...
                            'HandleVisibility', 'off');
                    end
                end
                hLim.shear = hShear;
            end
        end
    end

    xlabel(hAx, 'velocity [km/s]');
    if hasGreen
        ylabel(hAx, yLabelStr);
    end
    title(hAx, sprintf('Classification curve (angle = %.1f deg)', singleResult.angle));
end
