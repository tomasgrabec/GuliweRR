function [hPoints, hAx, hBg, hLim] = plotClassificationDiagnostic(sweepResult, colorBy, showGreenBackground, refLevel, hAx, component)
%%% Comparative mode-classification diagnostic: classified points over the
%%% Green map, with a Christoffel limiting-velocity reference line
%%% >> [hPoints, hAx, hBg, hLim] = plotClassificationDiagnostic(sweepResult, colorBy, showGreenBackground, refLevel, hAx, component)
%%%
%%% INPUT: 1. sweepResult         .. CalcAngSweepHS or CalcAngSweepLayered result
%%%                                  object, calc()'d. Reads .significantModesMap
%%%                                  directly — the already-classified struct
%%%                                  array per angle — and never recomputes a
%%%                                  parallel score
%%%  (opt) 2. colorBy             .. 'type' (default, categorical) | 'score'
%%%                                  (continuous: the entry's score in dB, on its
%%%                                  method's scale) | 'confidence' (continuous:
%%%                                  abs(modeVel-greenPeakVel), the bound Ritz
%%%                                  mode's distance from its Green peak; NaN
%%%                                  where no Green match was made)
%%%  (opt) 3. showGreenBackground .. draw the response background first (default
%%%                                  true; skipped with a warning if the chosen
%%%                                  background's data is missing, e.g. a Green
%%%                                  component on a 'ritzOnly' result)
%%%  (opt) 4. refLevel            .. how much classification detail to show
%%%                                  (default 'full'), in 3 cumulative tiers.
%%%                                  'other'-type points are never shown at any
%%%                                  tier — this is a classification read-out, not
%%%                                  a catalogue of every candidate considered.
%%%                                    'sawOnly'  -- SAW-family points only
%%%                                                  (Rayleigh/pseudoRayleigh/SAW),
%%%                                                  no reference line
%%%                                    'sawLimit' -- + ONE Christoffel reference
%%%                                                  line, by family: longitudinal
%%%                                                  for HS (Rayleigh is already
%%%                                                  subsonic w.r.t. shear, so the
%%%                                                  longitudinal bound is the
%%%                                                  informative one pseudoRayleigh
%%%                                                  approaches), substrate shear
%%%                                                  for layered (the bound above
%%%                                                  which the wave leaks into the
%%%                                                  substrate)
%%%                                    'full'     -- + 'limiting'-type points, the
%%%                                                  Ritz modes sitting AT a
%%%                                                  bulk-wave limiting velocity.
%%%                                                  For layered results this is
%%%                                                  what shows the substrate-shear
%%%                                                  boundary as discrete data
%%%  (opt) 5. hAx                 .. target axes (created fresh if omitted)
%%%  (opt) 6. component           .. response background: 'G13eff' (|G13_eff|,
%%%                                  what the Green score read; the default when
%%%                                  the result has Green maps) | any
%%%                                  'GIJ'/'GIJeff', recomputed on request
%%%                                  (Vis.greenComponentData) | 'spot' (the stored
%%%                                  calcSpotGreenMap result) | 'ritzSpectrum'
%%%                                  (plotRitzSpectrumMap, |S(v)| from the stored
%%%                                  modal numerators; the default on a 'ritzOnly'
%%%                                  result, which has no Green map)
%%% OUTPUT: hPoints .. colorBy='type': array of scatter handles, one per type
%%%                    present (legend shown); 'score'/'confidence': a single
%%%                    scatter handle (its own colour bar: UserData.colorbar)
%%%         hAx     .. axes used
%%%         hBg     .. the background's pcolor handle, or [] if not drawn
%%%         hLim    .. struct('long', <handle or []>, 'shear', <handle or []>).
%%%                    Only ONE is ever populated — long for HS, shear for layered
%%%                    — returned as a struct so a caller can toggle either
%%%                    without knowing which family sweepResult belongs to
%%%
%%% MARKER SIZE is scaled by each point's score: linear in the score (a log10
%%% amplitude ratio, so in dB) over its method's range,
%%% Vis.classificationScoreScale(...).sizeRange, clipped at both ends. Points with no
%%% score (none on a current result apart from 'limiting') get a fixed size;
%%% 'limiting' points a smaller fixed one.
%%%
%%% GRADE: entries graded 'weak' (accepted at the active, looser standard but not
%%% at the method's 'low' one) are drawn as open markers, 'strong' ones filled.
%%%
%%% 'score'/'confidence' colouring is plain `gray`, which already maps low data to
%%% black and high to white — NOT flipud(gray) — over the points' own range, with
%%% its own colour bar right of the background's: the background keeps its
%%% ColMapGrayOrange colormap and limits (explicit RGB points,
%%% visClassificationScatterByScore). Every point has a dark outline, so white
%%% points read against a white or bright background too.
%%%
%%% HS vs layered by class (visIsLayered); the per-angle limiting velocities are
%%% read duck-typed -- see docs/reference/architecture/visual.md, "Duck-typing across
%%% HS/Layered result classes".

    if nargin < 2 || isempty(colorBy); colorBy = 'type'; end
    visValidateChoice(colorBy, {'type', 'score', 'confidence'}, ...
        'plotClassificationDiagnostic', 'colorBy');
    if nargin < 3 || isempty(showGreenBackground); showGreenBackground = true; end
    if nargin < 4 || isempty(refLevel);            refLevel            = 'full'; end
    visValidateChoice(refLevel, {'sawOnly', 'sawLimit', 'full'}, ...
        'plotClassificationDiagnostic', 'refLevel');

    if nargin < 6 || isempty(component)
        if isempty(sweepResult.greenVelAxis) && visHasRitzSpectrum(sweepResult)
            component = 'ritzSpectrum';   % 'ritzOnly': the Ritz-side response
        else
            component = 'G13eff';
        end
    end
    if ~strcmp(component, 'ritzSpectrum')
        visValidateGreenComponent(component, 'plotClassificationDiagnostic');
    end

    if nargin < 5 || isempty(hAx)
        figure();
        hAx = axes();
    end
    hold(hAx, 'on');

    hBg = [];
    if showGreenBackground
        switch component
            case 'spot'
                noBg = ~isprop(sweepResult, 'spotGreenMap') || isempty(sweepResult.spotGreenMap);
            case 'ritzSpectrum'
                noBg = ~visHasRitzSpectrum(sweepResult);
            otherwise
                noBg = isempty(sweepResult.greenVelAxis);   % calcMode 'ritzOnly'
        end
        if noBg
            warning('Vis:plotClassificationDiagnostic:noGreenBackground', '%s -- skipping the background.', ...
                visMissingBackgroundReason(component, true));
        elseif strcmp(component, 'ritzSpectrum')
            hBg = Vis.plotRitzSpectrumMap(sweepResult, [], [], [], hAx);
            hold(hAx, 'on');
        else
            hBg = Vis.plotGreenMap(sweepResult, [], [], [], hAx, component);
        end
    end

    %%% flatten significantModesMap into flat per-point arrays, filtered per
    %%% refLevel -- shared rule, see visClassificationTypesAllowed
    typeAllowed = visClassificationTypesAllowed(refLevel);
    angles  = sweepResult.angles(:);
    nAngles = numel(angles);
    pts = visFlattenSignificantModes(sweepResult.significantModesMap, angles, typeAllowed);

    %%% marker size by score, weak entries open -- see docblock
    MIN_SZ = 12; MAX_SZ = 100; DEFAULT_SZ = 36; LIMITING_SZ = 16;
    isLayered = visIsLayered(sweepResult);
    sizeVals = visClassificationSizes(pts, isLayered, ...
        MIN_SZ, MAX_SZ, DEFAULT_SZ, LIMITING_SZ);
    isWeak = strcmp(pts.grade, 'weak');

    switch colorBy
        case 'type'
            %%% colours, draw order and legend come from
            %%% visClassificationScatterByType; the smaller LIMITING_SZ above
            %%% shrinks the halo left when a larger SAW/Rayleigh marker does not
            %%% fully cover a nearly-coincident limiting point
            hPoints = visClassificationScatterByType(hAx, pts.angle, pts.vel, sizeVals, pts.type, isWeak);

        otherwise   % 'score' or 'confidence': black=low, white=high -- see docblock
            hPoints = visClassificationScatterByScore(hAx, pts.angle, pts.vel, sizeVals, pts, colorBy, isWeak);
    end

    %%% ONE reference line only, picked by family: HS gets longitudinal, layered
    %%% gets substrate shear (see refLevel above)
    hLim = struct('long', [], 'shear', []);
    if ~strcmp(refLevel, 'sawOnly')
        if ~isLayered
            longVel = sweepResult.longVelMap;
            if ~isempty(longVel) && ~all(isnan(longVel))
                hLim.long = plot(hAx, angles, longVel, '-', ...
                    'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1.5, 'DisplayName', 'longitudinal limit');
            end
        else
            shearAngles = []; shearVels = [];
            for ia = 1:nAngles
                sv = visPropOrEmpty(sweepResult.cs{ia}, {'subShearLimVels'});
                if isempty(sv); continue; end
                shearAngles = [shearAngles; repmat(angles(ia), numel(sv), 1)]; %#ok<AGROW>
                shearVels   = [shearVels;   sv(:)];                            %#ok<AGROW>
            end
            if ~isempty(shearVels)
                hLim.shear = plot(hAx, shearAngles, shearVels, '+', ...
                    'Color', [0.4660 0.6740 0.1880], 'DisplayName', 'shear limit');
            end
        end
    end

    xlabel(hAx, 'angle [deg]');
    ylabel(hAx, 'velocity [km/s]');
end
