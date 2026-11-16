%% Tutorial 01 (code) -- isotropic half-space: one SAW velocity, checked against v_SAW ~ 0.92*v_T
% Narrated version: examples/tutorials/tutorial_01_isotropic_halfspace.m. Run guliwerr/init.m first.

%% Material -- isotropic Al (isotropic class takes [C11, C44] only), Hill average of the library's single crystal
matAlCrystal   = Material.FromLibrary('Al');
[alC11, alC44] = matAlCrystal.elConsts.GetHillAverages();   % GPa
alConsts   = [alC11, alC44];
symmNo     = 1;                    % 1 = isotropic (ElastConsts.SymmClassLibrary lists all)
density    = matAlCrystal.rho;     % g/cm^3
elastConsts = ElastConsts(alConsts, symmNo);
material    = Material(elastConsts, density);

%% Discretisation -- "easy" tier, legOrder >= 2*relDepth (docs/pitfalls.md MTH1)
relDepth   = 20;                   % domain depth, in wavelengths
legOrder   = 50;                   % Legendre basis order across that depth
calcParams = CalcParams(relDepth, legOrder);

%% Structure -- one half-space layer
structureName = 'al_isotropic_halfspace';
nMaxModes     = 500;               % max Ritz modes computed/retained
structure     = Structure(material, calcParams, structureName, nMaxModes);

%% Sweep -- one angle, Ritz only
angleDeg       = 0;                % isotropic: any angle gives the same result
calcMode       = 'ritzOnly';       % 'ritzOnly' | 'full' (adds Christoffel + Green)
retentionLevel = 'saw_only';       % 'full' | 'saw_only' | 'none'
nVelGreen      = [];               % Green-axis points ([] = 6 m/s steps); unused for 'ritzOnly'
sweep = Calc.angSweepHS(structure, angleDeg, calcMode, retentionLevel, nVelGreen);
sweep = sweep.calc();

%% Result -- SAW velocity vs. transverse bulk velocity
vSAW = sweep.sawVelMap{1}(1);      % km/s, first angle, first SAW
vT   = sqrt(alC44 / density);      % km/s
fprintf('v_SAW = %.4f km/s, v_T = %.4f km/s, ratio = %.3f\n', vSAW, vT, vSAW / vT);
