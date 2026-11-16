function [hLines, hAx, hBg] = plotDispersionCurves(sweepResult, curveType, angleIdx, showRitzBranches, showGreenBackground, scaleMode, climRange, hAx)
%%% Plot fixed-angle dispersion curves for a CalcWavevecAngSweep2Layered
%%% wavelength sweep -- one connected line per tracked physical branch,
%%% optionally over a Green's-|G| background
%%% >> [hLines, hAx, hBg] = Vis.plotDispersionCurves(sweepResult, curveType, angleIdx, showRitzBranches, showGreenBackground, scaleMode, climRange, hAx)
%%%
%%% INPUT: 1. sweepResult .. CalcWavevecAngSweep2Layered, already calc()'d,
%%%                          with retentionLevel='full' --
%%%                          DispersionBranches.track's own
%%%                          requirement (this function calls it directly
%%%                          and lets its errors propagate)
%%%  (opt) 2. curveType    .. 'wavelength-freq' (default) | 'wavelength-vphase'
%%%                          | 'k-freq' | 'k-vphase' | 'freq-vphase' | 'freq-vgroup'
%%%  (opt) 3. angleIdx     .. which angle to plot the dispersion curve at
%%%                          (default 1). sweepResult may carry more than one
%%%                          angle, but this is always a SINGLE fixed-angle
%%%                          curve set, never several at once
%%%  (opt) 4. showRitzBranches .. draw the tracked branch lines (default true).
%%%                          The lines are fully computed and returned in
%%%                          hLines either way — false only sets Visible='off'
%%%                          — so a caller can flip them back on without a
%%%                          redraw, the same convention the GUI viewers'
%%%                          per-type checkboxes use
%%%  (opt) 5. showGreenBackground .. draw a Green's |G13_eff| heatmap behind the
%%%                          lines (default false, unlike
%%%                          plotClassificationDiagnostic's default-true:
%%%                          Green data here needs calcMode='full' on the
%%%                          WHOLE wavelength sweep, and only 2 of the 6
%%%                          curveTypes can show it). Requires curveType
%%%                          'wavelength-vphase' or 'k-vphase' (velocity
%%%                          y-axis) and ERRORS otherwise, as it does when
%%%                          sweepResult was calc()'d with calcMode='ritzOnly'
%%%                          — an explicit opt-in cannot be silently ignored.
%%%                          The 4 frequency-based curveTypes are excluded
%%%                          because freq = vel/wavelength rescales the y-axis
%%%                          per wavelength column, making the background a
%%%                          non-rectangular grid. Both flags false is valid:
%%%                          an empty, labelled axes
%%%  (opt) 6. scaleMode    .. 'log' (default) | 'linear' -- same meaning as
%%%                          plotGreenMap's; used only with the background
%%%  (opt) 7. climRange    .. 2-element [lo hi], default [-4, -1], matching
%%%                          plotGreenMap (the same |G13_eff| quantity, with
%%%                          wavelength standing in for angle) rather than
%%%                          plotDispersionMap's [-5, -2], which is tuned for
%%%                          the isoFreq angle x kParRange family. Used only
%%%                          with the background
%%%  (opt) 8. hAx          .. target axes (created fresh if omitted)
%%% OUTPUT: hLines .. [1 x nCompModes] line handles, one per tracked branch, each
%%%                   tagged UserData.branchIdx for the GUI's click-to-select
%%%                   convention (ButtonDownFcn is set by the caller, not here)
%%%         hAx    .. axes used
%%%         hBg    .. background pcolor handle, or [] if showGreenBackground was false
%%%
%%% Branch continuity comes from DispersionBranches.track.
%%% Separate from plotDispersionMap.m rather than a case dispatched there:
%%% plotDispersionMap's CalcWavevecAngSweep2Layered branch is a multi-angle 2D
%%% scatter, this is a single-fixed-angle curve view — connected lines, a
%%% structurally different rendering.
%%%
%%% Branches classified SAW/limiting at the FIRST wavelength step are drawn
%%% highlighted (thicker, coloured), every other raw Ritz branch thin grey.
%%% Highlighting by physical significance rather than by symmetric/antisymmetric
%%% parity: a layered film-on-substrate structure has no such parity.
%%%
%%% freq = vel/wavelength, k = 2*pi/wavelength.
%%%
%%% GREEN BACKGROUND: |G13_eff| (greenAbsG13EffMap, the component the existence
%%% score reads and the only one a result stores). Each sweeps{iw}
%%% carries its own (greenVelAxis, greenAbsG13EffMap) under calcMode='full'. The axes normally coincide (the
%%% substrate limiting velocities and vMax do not change with wavelength), but
%%% nothing guarantees it, so plotWavevecGreenBackground interpolates every
%%% column onto one shared velocity axis first, leaving NaN outside a
%%% wavelength's own range, which pcolor draws blank.

    if nargin < 2 || isempty(curveType); curveType = 'wavelength-freq'; end
    visValidateChoice(curveType, visDispersionCurveTypes(), 'plotDispersionCurves', 'curveType');
    if nargin < 3 || isempty(angleIdx); angleIdx = 1; end
    validateattributes(angleIdx, {'double'}, ...
        {'scalar', 'integer', 'positive', '<=', numel(sweepResult.angles)}, ...
        'plotDispersionCurves', 'angleIdx');
    if nargin < 4 || isempty(showRitzBranches); showRitzBranches = true; end
    if nargin < 5 || isempty(showGreenBackground); showGreenBackground = false; end
    if nargin < 6 || isempty(scaleMode); scaleMode = 'log'; end
    visValidateChoice(scaleMode, {'log', 'linear'}, 'plotDispersionCurves', 'scaleMode');
    if nargin < 7 || isempty(climRange); climRange = [-4, -1]; end
    VPHASE_CURVE_TYPES = {'wavelength-vphase', 'k-vphase'};
    if showGreenBackground && ~ismember(curveType, VPHASE_CURVE_TYPES)
        error('Vis:plotDispersionCurves:backgroundNeedsVphase', ...
            ['showGreenBackground=true requires curveType ''wavelength-vphase'' or ''k-vphase'' ' ...
             '(velocity y-axis) -- got ''%s''. The 4 frequency-based curve types would need a ' ...
             'non-rectangular background grid (freq = vel/wavelength rescales per column); ' ...
             'not supported.'], curveType);
    end
    if nargin < 8 || isempty(hAx)
        figure();
        hAx = axes();
    end
    hold(hAx, 'on');

    hBg = [];
    if showGreenBackground
        hBg = plotWavevecGreenBackground(sweepResult, curveType, angleIdx, scaleMode, climRange, hAx);
    end

    branchModeNoMap = DispersionBranches.track(sweepResult, angleIdx);
    [nWav, nCompModes] = size(branchModeNoMap);
    wavelengths = sweepResult.wavelengths(:);

    sm1 = sweepResult.sweeps{1}.significantModesMap{angleIdx};
    isSig1 = ismember({sm1.type}, {'SAW', 'limiting'});
    highlightModeNos1 = [sm1(isSig1).modeNo];

    vel      = zeros(nWav, nCompModes);
    velGroup = zeros(nWav, nCompModes);
    for iw = 1:nWav
        ritz      = sweepResult.sweeps{iw}.cs{angleIdx}.ritzSolver;
        modeNoRow = branchModeNoMap(iw, :);
        vel(iw, :)      = ritz.vels(modeNoRow);
        velGroup(iw, :) = ritz.velsGroup(modeNoRow);
    end
    curveData = struct( ...
        'wavelength',      wavelengths, ...
        'k',               2*pi ./ wavelengths, ...   % [nWav x 1], standard angular wavenumber
        'vel',             vel, ...
        'velGroup',        velGroup, ...
        'freq',            vel ./ wavelengths, ...    % [nWav x nCompModes], broadcast over columns
        'wavelengthLabel', 'wavelength', ...
        'kLabel',          'k [2\pi/wavelength]', ...
        'freqLabel',       'frequency (= v_\phi / wavelength)');
    [xAll, yAll, xLabelStr, yLabelStr] = visDispersionCurveAxes(curveType, curveData);

    hLines = gobjects(1, nCompModes);
    for im = 1:nCompModes
        if ismember(branchModeNoMap(1, im), highlightModeNos1)
            hLines(im) = plot(hAx, xAll(:, im), yAll(:, im), 'LineWidth', 2, 'Color', [0.80 0.20 0.02]);
        else
            hLines(im) = plot(hAx, xAll(:, im), yAll(:, im), 'LineWidth', 0.5, 'Color', [0.6 0.6 0.6]);
        end
        hLines(im).UserData.branchIdx = im;
        %%% hidden, not skipped -- XData/YData stay usable (e.g. for a
        %%% crosshair) and a caller can flip Visible back on without a redraw
        if ~showRitzBranches
            hLines(im).Visible = 'off';
        end
    end
    xlabel(hAx, xLabelStr);
    ylabel(hAx, yLabelStr);
end

function hBg = plotWavevecGreenBackground(sweepResult, curveType, angleIdx, scaleMode, climRange, hAx)
%%% Green's |G13_eff| heatmap over (wavelength-or-k, velocity), fixed angle -- the
%%% wavevec-sweep counterpart of plotGreenMap's angle x velocity map. Every
%%% wavelength step's own greenVelAxis/greenAbsG13EffMap column is
%%% interpolated onto one shared velocity axis first; see the file docblock.
%%% >> hBg = plotWavevecGreenBackground(sweepResult, curveType, angleIdx, scaleMode, climRange, hAx)
    nWav = numel(sweepResult.sweeps);
    if isempty(sweepResult.sweeps{1}.greenAbsG13EffMap)
        error('Vis:plotDispersionCurves:noGreenData', ...
            ['greenAbsG13EffMap is empty -- sweepResult must be calc()''d with calcMode=''full'' ' ...
             'to populate Green''s function data (showGreenBackground=true needs it).']);
    end

    nVelGreen = numel(sweepResult.sweeps{1}.greenVelAxis);
    vMin =  Inf;
    vMax = -Inf;
    for iw = 1:nWav
        axLocal = sweepResult.sweeps{iw}.greenVelAxis;
        vMin = min(vMin, min(axLocal));
        vMax = max(vMax, max(axLocal));
    end
    velCommon = linspace(vMin, vMax, nVelGreen)';

    C = NaN(nVelGreen, nWav);
    for iw = 1:nWav
        axLocal = sweepResult.sweeps{iw}.greenVelAxis;
        gLocal  = double(sweepResult.sweeps{iw}.greenAbsG13EffMap(:, angleIdx));
        C(:, iw) = interp1(axLocal, gLocal, velCommon, 'linear', NaN);
    end

    %%% the out-of-range NaN cells above neither participate in the peak nor get
    %%% painted at the floor colour -- Vis.scaleGreenData's NaN contract
    [plotData, cbLabel] = Vis.scaleGreenData(C, scaleMode, 'G_{13}^{eff}');

    if strcmp(curveType, 'wavelength-vphase')
        xAxis = sweepResult.wavelengths(:);
    else   % 'k-vphase'
        xAxis = 2*pi ./ sweepResult.wavelengths(:);   % decreasing for increasing wavelengths;
                                                         % Vis.pcolorCenteredGrid supports either
                                                         % direction
    end

    hBg = visPcolorCentered(hAx, xAxis, velCommon, plotData);
    hBg.DisplayName = '|G_{13}^{eff}|';
    visGreenColourAxis(hAx, scaleMode, climRange, cbLabel);
end
