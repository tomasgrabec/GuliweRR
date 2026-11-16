%% Tutorial 06 (code) -- Cu film on Al: SAW velocity vs. film relDepth, one tracked branch
% Narrated version: examples/tutorials/tutorial_06_film_on_substrate.m. Run guliwerr/init.m first.

%% Substrate -- isotropic Al (faster than the Cu film)
matAlCrystal    = Material.FromLibrary('Al');   % single crystal; isotropic Al = its Hill average
[alC11, alC44]  = matAlCrystal.elConsts.GetHillAverages();
alConsts        = [alC11, alC44];  % [C11, C44] GPa
isoSymmNo       = 1;
alDensity       = matAlCrystal.rho;   % g/cm^3
relDepthSub     = 20;              % substrate depth, wavelengths

%% Film -- Cu(001), cubic, real thickness in nm
dirMatFilm      = [1 0 0; 0 0 1];
filmThicknessNm = 200;             % ApplyWavelength divides this by each wavelength

%% Discretisation -- legOrder >= 2x the largest relDepth reached (3 here)
legOrder             = 50;
relDepthPlaceholder  = 1;          % film relDepth, overwritten per wavelength

ecSub   = ElastConsts(alConsts, isoSymmNo);
matSub  = Material(ecSub, alDensity);
matFilm = Material.FromLibrary('Cu', dirMatFilm);   % constants + density from the materials library
cpSub   = CalcParams(relDepthSub, legOrder);
cpFilm  = CalcParams(relDepthPlaceholder, legOrder, filmThicknessNm);

structureName = 'al_sub_cu_film';
nMaxModes     = 500;
structure = Structure(matSub, cpSub, structureName, nMaxModes);
structure = structure.AddLayer(matFilm, cpFilm);

%% Wavelength sweep at one angle -- dense steps help branch tracking
filmRelDepths  = [0.01:0.01:1, 1.05:0.05:3];
wavelengthsNm  = filmThicknessNm ./ filmRelDepths;   % descending
angleDeg       = 0;
vMax           = [];               % not needed for 'ritzOnly'
calcMode       = 'ritzOnly';
retentionLevel = 'full';           % REQUIRED by DispersionBranches.track
sweep = Calc.wavevecAngSweep2Layered(structure, angleDeg, wavelengthsNm, vMax, calcMode, retentionLevel);
sweep = sweep.calc();

%% Branch tracking by eigenvector overlap, seeded at the thinnest film
vSAWTracked = DispersionBranches.trackSeededBranch(sweep).';   % seed: slowest SAW at the thinnest film

%% Half-space limits -- each material alone
calcModeHS = 'ritzOnly';
cpHalfSpace = CalcParams(relDepthSub, legOrder);   % each material as its own half-space
structAlHS = Structure(matSub, cpHalfSpace, 'al_hs', nMaxModes);
singleAl   = Calc.singleHS(structAlHS, angleDeg);
singleAl   = singleAl.calcAll(calcModeHS);
smAl       = singleAl.classification.significantModes;
vSAWAl     = min([smAl(strcmp({smAl.type}, 'Rayleigh')).vel]);
structCuHS = Structure(matFilm, cpHalfSpace, 'cu_hs', nMaxModes);
singleCu   = Calc.singleHS(structCuHS, angleDeg);
singleCu   = singleCu.calcAll(calcModeHS);
smCu       = singleCu.classification.significantModes;
vSAWCu     = min([smCu(strcmp({smCu.type}, 'Rayleigh')).vel]);

%% Plot 1 -- tracked branch between the two half-space limits
figure();
semilogx(filmRelDepths, vSAWTracked, 'o-');
yline(vSAWAl, '--', sprintf('Al substrate %.3f km/s', vSAWAl));
yline(vSAWCu, '--', sprintf('Cu film %.3f km/s', vSAWCu));
xlabel('film relDepth (\propto kd)'); ylabel('v_{SAW} [km/s]');

%% Plot 2 -- every tracked Ritz branch, SAW branches highlighted
curveType = 'k-vphase';            % 'wavelength-freq' | 'wavelength-vphase' | 'k-freq' | 'k-vphase' | 'freq-vphase' | 'freq-vgroup'
Vis.plotDispersionCurves(sweep, curveType);
xlim([0 0.03]); ylim([1 4]);

%% Interactive dispersion viewer, preset to k-vphase
gui = Gui.dispersionViewer(sweep);
gui.show();
gui.curveTypeDropdown.Value = curveType;
gui.onCurveTypeChanged();
gui.setAxisRange('x', [0 0.03]); gui.setAxisRange('y', [1 4]);

%% Optional: 'full' mode -- tracked branches over the layered |G13_eff| map
sweepFull = Calc.wavevecAngSweep2Layered(structure, angleDeg, wavelengthsNm, vMax, 'full', retentionLevel);
sweepFull = sweepFull.calc();
angleIdx = 1; showRitzBranches = true; showGreenBackground = true;
scaleMode = 'log'; climRange = [-3, -0.5];
Vis.plotDispersionCurves(sweepFull, curveType, angleIdx, showRitzBranches, showGreenBackground, scaleMode, climRange);
title('Al / Cu film: tracked branches over |G_{13}^{eff}|');
xlim([0 0.03]); ylim([1 4]);

% dispersion viewer = the combined viewer's counterpart for wavelength sweeps
guiFull = Gui.dispersionViewer(sweepFull);
guiFull.show();
guiFull.curveTypeDropdown.Value = curveType;
guiFull.showGreenBgCheckbox.Value = true;
guiFull.onCurveTypeChanged();
guiFull.setAxisRange('x', [0 0.03]); guiFull.setAxisRange('y', [1 4]);
