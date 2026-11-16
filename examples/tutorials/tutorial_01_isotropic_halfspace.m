%% Tutorial: isotropic half-space -- minimal SAW velocity calculation
%
% Shortest path from "a material" to a SAW velocity. Checks the textbook
% approximation v_SAW ~ 0.92*v_T (transverse bulk velocity) for an isotropic
% solid; its SAW velocity does not depend on angle, so one angle suffices.
%
% Copy this script as a starting point. Call sequence: docs/quickStart.md;
% the legO/relDepth rule and other gotchas: docs/pitfalls.md.
%
% Run run('guliwerr/init.m') once per MATLAB session from the project root
% first -- not from inside this script, since run() changes into the
% script's own folder and the relative path then fails.

%% Material -- isotropic Al. The isotropic class takes [c11, c44] only
% (c12 = c11 - 2*c44 is implied; ElastConsts.SymmClassLibrary lists every class).
% The values are the Hill average of the single-crystal Al constants from the
% materials library (Material.ListLibrary()), the usual isotropic estimate.
matAlCrystal   = Material.FromLibrary('Al');
[alC11, alC44] = matAlCrystal.elConsts.GetHillAverages();   % GPa
alConsts = [alC11, alC44];
alSymm   = 1;                     % isotropic
alRho    = matAlCrystal.rho;      % g/cm^3

ec  = ElastConsts(alConsts, alSymm);
mat = Material(ec, alRho);

%% Discretisation -- CalcParams(relDepth, legO), "easy" tier (isotropic)
relDepth = 20;  % depth of the computational domain, in wavelengths
legO     = 50;  % highest order of the Legendre basis spanning that depth --
                % keep legO >= 2*relDepth or the depth profile is
                % under-resolved (docs/pitfalls.md MTH1)
cp = CalcParams(relDepth, legO);

%% Structure -- a single half-space layer
strName   = 'al_isotropic_halfspace';
nMaxModes = 500;   % max number of Ritz modes computed/retained
str = Structure(mat, cp, strName, nMaxModes);

%% One-angle, Ritz-only sweep
angle          = 0;            % deg
calcMode       = 'ritzOnly';   % 'ritzOnly' (Ritz only) | 'full' (adds Christoffel + Green's function)
retentionLevel = 'saw_only';   % Ritz eigenvectors kept: 'full' (every mode) |
                               % 'saw_only' (classified modes) | 'none'
nVelGreen      = [];           % Green-axis points ([] = 6 m/s steps); unused for 'ritzOnly'
sweep = Calc.angSweepHS(str, angle, calcMode, retentionLevel, nVelGreen);
sweep = sweep.calc();

vSAW = sweep.sawVelMap{1}(1);   % km/s, first angle, first SAW
vT   = sqrt(alC44 / alRho);     % transverse bulk velocity, km/s

fprintf('v_SAW = %.4f km/s, v_T = %.4f km/s, ratio = %.3f\n', vSAW, vT, vSAW / vT);
