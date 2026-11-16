%% Tutorial 05 (code) -- layered 'full' angular sweep: Ni film on Cu(001), Green map + classification
% Narrated version: examples/tutorials/tutorial_05_layered_green.m. Run guliwerr/init.m first.

%% Substrate -- Cu(001), cubic (materials library, see below)
dirMatSub    = [1 0 0; 0 0 1];     % surface normal; propagation reference

%% Film -- Ni(001), cubic (materials library, see below)
dirMatFilm   = [1 0 0; 0 0 1];

%% Discretisation -- "easy" tier ((001) principal plane); film relDepth set directly (no ApplyWavelength)
relDepthSub  = 20;                 % substrate depth, wavelengths
relDepthFilm = 0.2;                % film thickness, wavelengths
legOrder     = 50;

matSub  = Material.FromLibrary('Cu', dirMatSub);   % constants + density from the materials library
matFilm = Material.FromLibrary('Ni', dirMatFilm);
cpSub   = CalcParams(relDepthSub, legOrder);
cpFilm  = CalcParams(relDepthFilm, legOrder);

%% Structure -- layer 1 = substrate, AddLayer puts the film on top
structureName = 'cu_sub_ni_film';
nMaxModes     = 500;
structure = Structure(matSub, cpSub, structureName, nMaxModes);
structure = structure.AddLayer(matFilm, cpFilm);

%% Angular sweep -- vMax is optional ([] = Green axis ends at 1.05 x the fastest vL of any layer)
vMax           = [];               % km/s; a larger value only extends the displayed map
anglesDeg      = (0:2:90)';
calcMode       = 'full';           % 'full' | 'ritzOnly'
retentionLevel = 'full';           % 'full' | 'saw_only' | 'none'
nVelGreen      = [];               % unused slot for layered stacks: must be empty (4 m/s axis)
sweep = Calc.angSweepLayered(structure, anglesDeg, vMax, calcMode, retentionLevel, nVelGreen);
sweep = sweep.calc();

%% |G_13,eff| map with a widened colour range (6th argument e.g. 'G33' recomputes |G_33|)
scaleMode = 'log';                 % 'log' | 'linear'
climRange = [-3, -0.5];            % log10(|G|/max|G|) colour limits (default [-4 -1])
Vis.plotGreenMap(sweep, scaleMode, climRange);
title('Cu / Ni film: |G_{13}^{eff}|');

%% Classification diagnostic -- same climRange applied via the returned axes
colorBy             = 'type';      % 'type' | 'score' | 'confidence'
showGreenBackground = true;
refLevel            = 'full';      % 'sawOnly' | 'sawLimit' | 'full'
[~, axClass] = Vis.plotClassificationDiagnostic(sweep, colorBy, showGreenBackground, refLevel);
clim(axClass, climRange);
title(axClass, 'Cu / Ni film: classification vs. |G_{13}^{eff}|');
