%% Tutorial 02 (code) -- anisotropic half-space: Cu(001) angular sweep + standard plots
% Narrated version: examples/tutorials/tutorial_02_angular_map.m. Run guliwerr/init.m first.

%% Material -- Cu, cubic, (001) surface
dirMat   = [1 0 0; 0 0 1];          % row 1 = surface normal, row 2 = propagation at angle 0
material = Material.FromLibrary('Cu', dirMat);   % constants and density from the materials library

%% Discretisation -- "easy" tier (Zener ratio 2*C44/(C11-C12) = 3.2, principal plane)
relDepth   = 20;
legOrder   = 50;
calcParams = CalcParams(relDepth, legOrder);

%% Structure
structureName = 'cu_001_100to010';
nMaxModes     = 500;
structure     = Structure(material, calcParams, structureName, nMaxModes);

%% Angular sweep -- Ritz + Christoffel + Green
anglesDeg      = (0:2:90)';        % [100] -> [010]
calcMode       = 'full';           % 'full' | 'ritzOnly'
retentionLevel = 'saw_only';       % 'full' | 'saw_only' | 'none'
nVelGreen      = [];               % Green-axis points ([] = greenVelStep steps, 6 m/s)
useParallel    = "auto";           % "auto" | "always" | "never"
sweep = Calc.angSweepHS(structure, anglesDeg, calcMode, retentionLevel, nVelGreen, useParallel);
sweep = sweep.calc();

%% SAW velocity vs. angle -- cartesian (velocity) and polar (slowness)
plotTypeCartesian = 'cartesian';   % 'cartesian' | 'polar'
Vis.plotAngularMap(sweep, plotTypeCartesian);
title('Cu (001): SAW velocity');
plotTypePolar = 'polar';
Vis.plotAngularMap(sweep, plotTypePolar);
thetalim([0 90])
title('Cu (001): SAW slowness');

%% |G_13,eff| Green map (the spectrum the classification read) -- cartesian and polar
scaleMode = 'log';                 % 'log' | 'linear'
climRange = [-3, -1];              % log10(|G|/max|G|) colour limits (default [-4 -1])
Vis.plotGreenMap(sweep, scaleMode, climRange, plotTypeCartesian);
ylim([1 5.5])
title('Cu (001): |G_{13}^{eff}|');
Vis.plotGreenMap(sweep, scaleMode, climRange, plotTypePolar);
axis([0 5.5 0 5.5])
title('Cu (001): |G_{13}^{eff}| (polar)');

%% Classification diagnostic -- classified modes over the Green map
colorBy             = 'type';      % 'type' | 'score' | 'confidence'
showGreenBackground = true;
refLevel            = 'full';      % 'sawOnly' | 'sawLimit' | 'full'
Vis.plotClassificationDiagnostic(sweep, colorBy, showGreenBackground, refLevel);
clim([-3 -1])
ylim([1 5.5])
title('Cu (001): classification vs. |G_{13}^{eff}|');
% switch without recomputing: sweep.reclassify('low'), sweep.classify('ritz' | 'refined')

%% All Ritz modes, shaded by a recognition criterion (diagnostic)
criterion       = 'depthInsensitivity';   % 'ritzScore' (default) | 'freqDisp' | 'depthInsensitivity'
velCutoffFactor = 1.10;                   % x max longitudinal velocity; Inf = no cutoff
Vis.plotRitzModeMap(sweep, criterion, velCutoffFactor);
clim([1 8])
title('Cu (001): all Ritz modes');

%% Any other component, recomputed on request (not stored) -- here |G_33|
targetAx  = [];                    % [] = new figure
component = 'G33';                 % 'G13eff' (default, stored) | 'GIJ' | 'GIJeff', e.g. 'G13', 'G33'
Vis.plotGreenMap(sweep, scaleMode, climRange, plotTypeCartesian, targetAx, component);
ylim([1 5.5])
title('Cu (001): |G_{33}| (recomputed)');
