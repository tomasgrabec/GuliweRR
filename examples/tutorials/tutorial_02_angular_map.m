%% Tutorial: anisotropic half-space -- full angular sweep + visualisation
%
% Cu(001) surface, propagation swept [100] -> [010]. A 'full' sweep (Ritz +
% Christoffel + Green), then the standard plots in the order you would reach
% for them: SAW velocity/slowness vs. angle, the |G_13,eff| Green's-function
% map (cartesian and polar), the classification diagnostic, every Ritz mode
% (exploratory), and a recomputed one-sided |G_13+| map for comparison.
%
% Copy this script as a starting point for an angular sweep. Call sequence:
% docs/quickStart.md; gotchas: docs/pitfalls.md.
%
% Run run('guliwerr/init.m') from the project root first (tutorial_01 says why).

%% Material -- Cu(001), cubic. dirMat is Material's default here, set
% explicitly to show where an orientation goes (tutorial_S1 builds others).
% The Cu constants and density come from the materials library
% (Material.ListLibrary() shows values and sources).
dirMat = [1 0 0; ...         % surface normal
          0 0 1];            % propagation reference (Material's default)
mat    = Material.FromLibrary('Cu', dirMat);

%% Discretisation -- "easy" tier (docs/pitfalls.md MTH1): Cu's Zener
% ratio 2*C44/(C11-C12) = 3.2 < 5 on a principal plane. Use the "complex"
% tier CalcParams(30, 70) for a non-principal cut or a Zener ratio of 5-10.
relDepth = 20;   % depth of the computational domain, in wavelengths
legO     = 50;   % Legendre basis order across that depth; keep >= 2*relDepth
cp = CalcParams(relDepth, legO);

strName   = 'cu_001_100to010';
nMaxModes = 500;   % max number of Ritz modes computed/retained
str = Structure(mat, cp, strName, nMaxModes);

%% Angular sweep, 0:2:90 deg ([100] -> [010]), full mode (Ritz+Christoffel+Green)
angles         = (0:2:90)';
calcMode       = 'full';       % 'full' (Ritz+Christoffel+Green) | 'ritzOnly'
                               % (SAW velocity only, no Green map -- much faster)
retentionLevel = 'saw_only';   % 'full' | 'saw_only' | 'none' -- 'full' only
                               % needed for modal-shape plots (docs/pitfalls.md NUM4)
nVelGreen      = [];     % number of Green-axis points ([] = greenVelStep steps, 6 m/s)
useParallel    = "auto"; % "auto" (thread pool for >=20 angles) | "always" |
                         % "never" (use when debugging -- docs/pitfalls.md NUM5)
sweep = Calc.angSweepHS(str, angles, calcMode, retentionLevel, nVelGreen, useParallel);
sweep = sweep.calc();

%% First look -- SAW velocity vs. angle (also available from a 'ritzOnly'
% sweep). Cartesian (velocity) and polar (slowness).
Vis.plotAngularMap(sweep);
title('Cu (001), [100]\rightarrow[010]: SAW velocity vs. angle');

plotType = 'polar';   % 'cartesian' (default, velocity) | 'polar' (slowness = 1/velocity)
Vis.plotAngularMap(sweep, plotType);
thetalim([0 90])   % cubic (001): 0-90 deg covers every direction
title('Cu (001), [100]\rightarrow[010]: SAW slowness (polar)');

%% |G13_eff| map (log scale by default) -- the +-k-symmetrised TGS response,
% the spectrum the classification read and the only Green component stored
Vis.plotGreenMap(sweep);
clim([-3 -1])
ylim([1 5.5])   % the Green axis starts at 0.3 km/s, far below any surface wave
title('Cu (001), [100]\rightarrow[010]: |G_{13}^{eff}|');

%% Same map, polar -- the angle x velocity heatmaps (plotGreenMap,
% plotDispersionMap) take plotType as their 4th argument.
scaleMode = [];        % [] = default 'log' (log10(|G|/max|G|)) | 'linear' (raw |G|)
climRange = [];        % [] = default [-4 -1] for 'log' (1e-4..1e-1 of the peak); auto for 'linear'
plotType  = 'polar';   % 'cartesian' (default) | 'polar'
Vis.plotGreenMap(sweep, scaleMode, climRange, plotType);
clim([-3 -1])
axis([0 5.5 0 5.5])   % plain axes, velocity as radius: crop to 5.5 km/s
title('Cu (001), [100]\rightarrow[010]: |G_{13}^{eff}| (polar)');

%% Classification diagnostic: classified SAW/limiting points and the
% longitudinal Christoffel reference line over the Green map, to check the
% classified points against the Green peaks.
colorBy              = 'type';   % 'type' (categorical) | 'score' (classification score, dB) |
                                 % 'confidence' (|v_Ritz - v_GreenPeak|)
showGreenBackground  = true;     % draw |G13_eff| underneath (a 'ritzOnly' result draws its Ritz spectrum instead)
refLevel             = 'full';   % 'sawOnly' (SAW points only) | 'sawLimit' (+ longitudinal
                                 % reference line, HS) | 'full' (+ limiting-type points)
Vis.plotClassificationDiagnostic(sweep, colorBy, showGreenBackground, refLevel);
clim([-3 -1])
ylim([1 5.5])
title('Cu (001), [100]\rightarrow[010]: mode classification vs. |G_{13}^{eff}|');
% These are the default 'green' classification at its 'high' standard. Switch
% without recomputing: sweep.reclassify('low') (another standard),
% sweep.classify('ritz') or sweep.classify('refined') (another method; 'refined'
% solves its spot directions once) -- docs/quickStart.md, "Choosing the classification".

%% Exploratory: every computed Ritz mode, not only the classified ones,
% shaded by a criterion -- shows what the classification left out.
criterion       = 'depthInsensitivity';   % 'ritzScore' (default) | 'freqDisp' | 'depthInsensitivity'
velCutoffFactor = 1.10;   % hide modes above 1.1x the max longitudinal limiting velocity ([] = default, Inf = none)
Vis.plotRitzModeMap(sweep, criterion, velCutoffFactor);
clim([1 8])
title('Cu (001), [100]\rightarrow[010]: all Ritz modes (depthInsensitivity)', 'Interpreter', 'tex');
% (plotRitzModeMap sets its own title with Interpreter 'none'; 'tex' restores \rightarrow)

%% Optional: any other Green component, e.g. the out-of-plane |G_33|.
% Not stored: recomputed here, one solve per angle.
component = 'G33';   % 'G13eff' (default, stored) | any 'GIJ' / 'GIJeff', e.g. 'G13', 'G33'
Vis.plotGreenMap(sweep, scaleMode, climRange, 'cartesian', [], component);
title('Cu (001), [100]\rightarrow[010]: |G_{33}| (recomputed)');
clim([-3 -1])
ylim([1 5.5])
