function [hLines, hAx, hChristLines] = plotDispersionCurvesPlate(sweepResult, curveType, angleIdx, colorMode, confidenceThresh, showChristoffelBounds, hAx)
%%% Plot fixed-angle dispersion curves for a CalcWavevecAngSweep2Plate
%%% wavevector sweep -- one connected line per tracked physical branch,
%%% optionally with Christoffel bulk/limiting-velocity reference lines.
%%% >> [hLines, hAx, hChristLines] = Vis.plotDispersionCurvesPlate(sweepResult, curveType, angleIdx, colorMode, confidenceThresh, showChristoffelBounds, hAx)
%%%
%%% Separate from Vis.plotDispersionCurves (the CalcWavevecAngSweep2Layered
%%% analogue) rather than a class-dispatch branch inside it, the two cases
%%% differing in more than cosmetics:
%%%   - the swept axis is CalcWavevecAngSweep2Plate's own `wavVecs` directly, a
%%%     real wavevector, not a wavelength converted to k; wavelength here is
%%%     DERIVED as 2*pi./wavVecs, the reverse of the Layered function.
%%%   - there is no "significant SAW/Sezawa mode" concept to highlight — every
%%%     Plate mode is a genuine physical branch — so highlighting uses
%%%     classifyBranches' mode-type labels, with its confidence-driven fallback
%%%     to plain branch identity, instead of significantModesMap.
%%%   - there is no Green's-function background: no analogue exists for a free
%%%     plate, so there is no showGreenBackground parameter at all.
%%%   - Christoffel bulk/limiting velocities stay meaningful without Green data,
%%%     and appear here as optional reference lines.
%%%
%%% INPUT: 1. sweepResult .. CalcWavevecAngSweep2Plate, already calc()'d, with
%%%                          retentionLevel='full' — classifyBranches' own
%%%                          requirement, whose errors propagate from here
%%%  (opt) 2. curveType    .. 'wavelength-freq' (default) | 'wavelength-vphase'
%%%                          | 'k-freq' | 'k-vphase' | 'freq-vphase' |
%%%                          'freq-vgroup' — the same 6 names
%%%                          plotDispersionCurves uses. 'k' here IS
%%%                          sweepResult.wavVecs, 'wavelength' = 2*pi./wavVecs,
%%%                          and 'freq' is sweepResult.freqsPhysicalMap3D,
%%%                          already in whatever real units realThickness/wavVecs
%%%                          were given in and not re-derived here
%%%  (opt) 3. angleIdx     .. which angle to plot at (default 1)
%%%  (opt) 4. colorMode    .. 'type' (default) | 'branch'. 'type' colours each
%%%                          branch by mode-type (Asymm/Symm/TH-Asymm/TH-Symm)
%%%                          where classifyBranches trusts the label, by its
%%%                          parity class (Vis.plateModeTypeStyle's parityColors,
%%%                          dashed) where only the parity is exact, evaluated
%%%                          at this curve's FIRST wavevector step — a
%%%                          per-branch, not per-point decision, so a branch
%%%                          whose confidence drops later keeps its first-step
%%%                          colour. 'branch' colours by branch identity
%%%                          (MATLAB `lines`, cyclic beyond 7), ignoring
%%%                          classification: the escape hatch for a stack where
%%%                          the labels are not trusted
%%%  (opt) 5. confidenceThresh .. passed to classifyBranches (default 0.3, with
%%%                          the same unvalidated-default caveat)
%%%  (opt) 6. showChristoffelBounds .. default true when curveType has a velocity
%%%                          y-axis, false otherwise — a phase-velocity bound has
%%%                          no meaning on a frequency or group-velocity axis, so
%%%                          it is silently not drawn rather than an error, this
%%%                          being an additive overlay and not an explicit data
%%%                          request. Draws the bottom-layer and, if different,
%%%                          top-layer longitudinal and shear Christoffel
%%%                          limiting velocities as horizontal lines, rotated to
%%%                          sweepResult.angles(angleIdx). Both surfaces are
%%%                          shown because at high kd an asymmetric stack's
%%%                          surface modes decouple toward each surface's OWN
%%%                          Rayleigh-adjacent velocity independently
%%%  (opt) 7. hAx          .. target axes (created fresh if omitted)
%%% OUTPUT: hLines       .. [1 x nCompModes] line handles, one per tracked branch,
%%%                         each tagged UserData.branchIdx for click-to-select
%%%         hAx          .. axes used
%%%         hChristLines .. Christoffel reference-line handles (yline objects),
%%%                         or empty gobjects(1,0) if not drawn

    if nargin < 2 || isempty(curveType); curveType = 'wavelength-freq'; end
    visValidateChoice(curveType, visDispersionCurveTypes(), 'plotDispersionCurvesPlate', 'curveType');
    if nargin < 3 || isempty(angleIdx); angleIdx = 1; end
    validateattributes(angleIdx, {'double'}, ...
        {'scalar', 'integer', 'positive', '<=', numel(sweepResult.angles)}, ...
        'plotDispersionCurvesPlate', 'angleIdx');
    if nargin < 4 || isempty(colorMode); colorMode = 'type'; end
    visValidateChoice(colorMode, {'type', 'branch'}, 'plotDispersionCurvesPlate', 'colorMode');
    if nargin < 5 || isempty(confidenceThresh); confidenceThresh = 0.3; end
    VPHASE_CURVE_TYPES = {'wavelength-vphase', 'k-vphase', 'freq-vphase'};
    if nargin < 6 || isempty(showChristoffelBounds)
        showChristoffelBounds = ismember(curveType, VPHASE_CURVE_TYPES);
    end
    if nargin < 7 || isempty(hAx)
        figure();
        hAx = axes();
    end
    hold(hAx, 'on');

    classification = sweepResult.classifyBranches(angleIdx, confidenceThresh);
    branchModeNoMap = classification.branchModeNoMap;
    [nWav, nCompModes] = size(branchModeNoMap);
    wavVecs = sweepResult.wavVecs(:);

    vel      = NaN(nWav, nCompModes);
    velGroup = NaN(nWav, nCompModes);
    freq     = NaN(nWav, nCompModes);
    for iw = 1:nWav
        rawIdx          = branchModeNoMap(iw, :);
        vel(iw, :)      = sweepResult.velsMap3D(rawIdx, angleIdx, iw).';
        velGroup(iw, :) = sweepResult.velsGroupMap3D(rawIdx, angleIdx, iw).';
        freq(iw, :)     = sweepResult.freqsPhysicalMap3D(rawIdx, angleIdx, iw).';
    end
    curveData = struct( ...
        'wavelength',      2*pi ./ wavVecs, ...   % DERIVED here, unlike the Layered side
        'k',               wavVecs, ...           % the swept wavevector itself
        'vel',             vel, ...
        'velGroup',        velGroup, ...
        'freq',            freq, ...
        'wavelengthLabel', 'wavelength', ...
        'kLabel',          'k', ...
        'freqLabel',       'frequency');
    [xAll, yAll, xLabelStr, yLabelStr] = visDispersionCurveAxes(curveType, curveData);

    [typeNames, typeColors, parityNames, parityColors] = Vis.plateModeTypeStyle();
    linesCmap = lines(max(nCompModes, 7));

    hLines = gobjects(1, nCompModes);
    for im = 1:nCompModes
        isTypeMode = strcmp(colorMode, 'type');
        typeNo     = classification.modesTypeByBranch(1, im);
        lineStyle  = '-';
        if isTypeMode && classification.trustLabelByBranch(1, im)
            col      = typeColors(typeNo, :);
            dispName = typeNames{typeNo};
        elseif isTypeMode && classification.isExactMirrorSymmetric
            %%% parity exact, Lamb/TH not: the family colour, dashed
            col       = parityColors(typeNo, :);
            dispName  = parityNames{typeNo};
            lineStyle = '--';
        else
            col      = linesCmap(mod(im-1, size(linesCmap, 1)) + 1, :);
            dispName = sprintf('branch %d', im);
        end
        hLines(im) = plot(hAx, xAll(:, im), yAll(:, im), 'LineWidth', 1.5, ...
            'Color', col, 'LineStyle', lineStyle, 'DisplayName', dispName);
        hLines(im).UserData.branchIdx = im;
    end
    xlabel(hAx, xLabelStr);
    ylabel(hAx, yLabelStr);

    hChristLines = gobjects(1, 0);
    if showChristoffelBounds
        structure = sweepResult.structure;
        %%% CalcWavevecAngSweep2Plate.structure is kept UNROTATED, so the angle
        %%% is passed through for rotation here
        angleDeg  = sweepResult.angles(angleIdx);
        [longLimBot, shearLimBot] = visPlateLimVels(structure, 1, angleDeg);
        hChristLines(end+1) = yline(hAx, shearLimBot, '--', 'bottom surface shear limit', ...
            'Color', [0.3 0.3 0.3], 'LabelHorizontalAlignment', 'left');
        hChristLines(end+1) = yline(hAx, longLimBot, ':', 'bottom surface long. limit', ...
            'Color', [0.3 0.3 0.3], 'LabelHorizontalAlignment', 'left');
        if structure.nLay > 1
            [longLimTop, shearLimTop] = visPlateLimVels(structure, structure.nLay, angleDeg);
            sameSurface = abs(shearLimTop - shearLimBot) < 1e-9 && abs(longLimTop - longLimBot) < 1e-9;
            if ~sameSurface
                hChristLines(end+1) = yline(hAx, shearLimTop, '--', 'top surface shear limit', ...
                    'Color', [0.6 0.6 0.6], 'LabelHorizontalAlignment', 'right');
                hChristLines(end+1) = yline(hAx, longLimTop, ':', 'top surface long. limit', ...
                    'Color', [0.6 0.6 0.6], 'LabelHorizontalAlignment', 'right');
            end
        end
    end
end
