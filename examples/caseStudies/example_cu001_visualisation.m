%% Example: Cu (001) surface, [100] -> [010] angular sweep and its standard plots
%
% G13_eff Green's-function map, SAW/limiting-velocity classification map, and
% SAW modal-shape displacement fields at 0/20/40 deg.
%
% Surface: (001), propagation swept 0:1:90 deg from [100] to [010].
% Material's default dirMat is Cu(001), so no dirMat is set explicitly below.
%
% Run run('guliwerr/init.m') from the project root first (see
% examples/tutorials/tutorial_01_isotropic_halfspace.m for why not from inside it).

%% Material -- Cu-like cubic constants
cuConsts = [170, 122, 74];   % [C11, C12, C44] GPa
cuSymm   = 2;                % cubic
cuRho    = 8.96;             % g/cm^3

%% Mesh resolution -- "easy" tier (docs/pitfalls.md MTH1)
legO     = 50;
relDepth = 20;

ec  = ElastConsts(cuConsts, cuSymm);
mat = Material(ec, cuRho);              % default dirMat = Cu(001)
cp  = CalcParams(relDepth, legO);

strName  = 'cu_001_100to010';
maxModes = 500;
str = Structure(mat, cp, strName, maxModes);

%% Angular sweep, 0:1:90 deg ([100] -> [010]), full mode (Ritz+Christoffel+Green)
angles         = (0:1:90)';
calcMode       = 'full';
retentionLevel = 'full';
nVelGreen      = [];     % number of Green-axis points ([] = greenVelStep steps, 6 m/s)
sweep = Calc.angSweepHS(str, angles, calcMode, retentionLevel, nVelGreen);
sweep = sweep.calc();

%% |G13_eff| map (log-scale default, see plotGreenMap.m) -- what the classification read
Vis.plotGreenMap(sweep);
ylim([1 5])   % the Green axis starts at 0.3 km/s
title('Cu (001), [100]\rightarrow[010]: |G_{13}^{eff}|');

%% Classification diagnostic: classified SAW/limiting points and the
% longitudinal Christoffel reference line over the Green map
colorBy             = 'type';
showGreenBackground = true;
refLevel            = 'full';
Vis.plotClassificationDiagnostic(sweep, colorBy, showGreenBackground, refLevel);
title('Cu (001), [100]\rightarrow[010]: mode classification vs. |G_{13}^{eff}|');

%% SAW modal-shape displacement fields at 0, 20, 40 deg
% substrateCutoff 'none': true-scale depth (the other options are for films)
dispAngles      = [0, 20, 40];
substrateCutoff = 'none';
cutoffParam     = [];
dirs            = 1:3;
meshSpacing     = 0.01 / 2;
for a = dispAngles
    idx = find(sweep.angles == a, 1);
    if isempty(idx) || isempty(sweep.sawModeNoMap{idx})
        warning('example_cu001_visualisation:noSAWMode', ...
            'No SAW mode identified at %g deg -- skipping modal shape.', a);
        continue
    end
    modeNo = sweep.sawModeNoMap{idx}(1);
    vSAW   = sweep.sawVelMap{idx}(1);
    ritzSolver = sweep.cs{idx}.ritzSolver;
    Vis.plotModalShape(ritzSolver, modeNo, substrateCutoff, cutoffParam, dirs, meshSpacing);
    sgtitle(sprintf('Cu (001), \\theta = %g deg: %s mode, v = %.3f km/s', ...
        a, sweep.sawTypeMap{idx}{1}, vSAW));
end
