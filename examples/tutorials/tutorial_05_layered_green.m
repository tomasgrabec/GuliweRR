%% Tutorial: layered full angular sweep -- Green's function on a film+substrate stack
%
% The layered counterpart of tutorial_02's half-space sweep: a film-on-
% substrate stack in calcMode 'full' (Ritz + Christoffel + Green). The layered
% launchers take one extra argument, vMax; it is optional and only extends
% the Green velocity axis beyond its default end.
%
% Copy this script as a starting point for a layered 'full' calculation.
% Call sequence: docs/quickStart.md.
%
% Run run('guliwerr/init.m') from the project root first (tutorial_01 says why).

%% Substrate -- Cu(001), cubic, from the materials library
% (Material.ListLibrary() shows values and sources).
dirMatSub = [1 0 0; ...      % surface normal
             0 0 1];         % propagation reference (Material's default)

%% Film -- Ni(001), cubic (constants and density from the materials library)
dirMatFilm = [1 0 0; ...      % surface normal
              0 0 1];         % propagation reference (Material's default)

%% Discretisation -- "easy" tier (docs/pitfalls.md MTH1): both layers are
% (001) cuts, a principal plane, with Zener ratios 2*C44/(C11-C12) of 3.2 (Cu)
% and 2.7 (Ni), both in the "easy" tier as a starting point. The film
% relDepth is set directly here; tutorial_06/07 derive it from a real
% thickness with Structure.ApplyWavelength. An entry's vel is its Green peak
% velocity. A leaky entry (above the substrate's shear velocity) is bound to the
% Ritz mode nearest that peak, but the slab has no mode exactly at a leaky
% resonance: the mode's own velocity (modeVel) can sit up to half a mode spacing
% (here a few m/s) off the peak, at any tier.
relDepthSub  = 20;   % substrate depth, in wavelengths
relDepthFilm = 0.2;  % film thickness, in wavelengths
legO         = 50;   % Legendre basis order per layer; keep >= 2*relDepth

matSub = Material.FromLibrary('Cu', dirMatSub);
cpSub  = CalcParams(relDepthSub, legO);

matFilm = Material.FromLibrary('Ni', dirMatFilm);
cpFilm  = CalcParams(relDepthFilm, legO);

strName   = 'cu_sub_ni_film';
nMaxModes = 500;   % max number of Ritz modes computed/retained
str = Structure(matSub, cpSub, strName, nMaxModes);   % layer 1 = substrate
str = str.AddLayer(matFilm, cpFilm);                       % layer 2 = film, on top

%% Green velocity-axis upper end -- optional. [] ends the axis at 1.05 x the
% fastest longitudinal limiting velocity of any layer (the end of the band
% the classification searches); a larger value only extends the displayed
% Green map.
vMax = [];

%% Angular sweep, 0:2:90 deg ([100] -> [010]), full mode (Ritz+Christoffel+Green)
angles         = (0:2:90)';
calcMode       = 'full';
retentionLevel = 'full';   % 'full' | 'saw_only' | 'none' -- 'full' keeps modal shapes for a viewer
nVelGreen      = [];       % unused slot for layered stacks: must be empty (fixed 4 m/s Green axis)
sweep = Calc.angSweepLayered(str, angles, vMax, calcMode, retentionLevel, nVelGreen);
sweep = sweep.calc();

%% |G13_eff| map -- the spectrum the classification read; log scale, colour range
% [-3 -0.5] instead of the [-4 -1] default (a narrower, higher window that
% keeps the strong branches distinct; lower climRange(1) to see weaker ones).
% Pass e.g. 'G33' as the 6th argument (component) to recompute |G33| instead.
scaleMode = 'log';       % 'log' (default) | 'linear'
climRange = [-3, -0.5];  % log10(|G|/max|G|) colour limits
Vis.plotGreenMap(sweep, scaleMode, climRange);
title('Cu substrate / Ni film, [100]\rightarrow[010]: |G_{13}^{eff}|');

%% Classification diagnostic: classified points and the substrate shear
% reference line over the Green map. On a stack every accepted wave has the
% single type 'SAW' -- no Rayleigh/Sezawa split (docs/pitfalls.md PHY2);
% marker size scales with the score. The background uses the default colour
% range, so the returned axes get the same climRange afterwards.
colorBy             = 'type';   % 'type' | 'score' | 'confidence'
showGreenBackground = true;     % draw plotGreenMap underneath
refLevel            = 'full';   % 'sawOnly' | 'sawLimit' (+ substrate shear line, layered) |
                                % 'full' (+ limiting-type points)
[~, axClass] = Vis.plotClassificationDiagnostic(sweep, colorBy, showGreenBackground, refLevel);
clim(axClass, climRange);
title(axClass, 'Cu substrate / Ni film, [100]\rightarrow[010]: mode classification vs. |G_{13}^{eff}|');

