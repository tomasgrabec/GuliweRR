%% Tutorial 07 (code) -- Cu film on Si(100): angle sweep, wavelength sweep, double sweep + viewers
% Narrated version: examples/tutorials/tutorial_07_dispersion_map_and_viewer.m. Run guliwerr/init.m first.

%% Substrate -- Si(100), cubic (materials library)
dirMatSi     = [1 0 0; 0 1 0];     % (100) surface, [010] reference
relDepthSub  = 20;                 % not swept by ApplyWavelength

%% Film -- sputtered Cu, isotropic via the Hill average of cubic Cu (materials library)
matCuCrystal    = Material.FromLibrary('Cu');   % single crystal
cuDensity       = matCuCrystal.rho;             % g/cm^3
isoSymmNo       = 1;
[cuC11Hill, cuC44Hill] = matCuCrystal.elConsts.GetHillAverages();   % isotropic [M, G] = [C11, C44]
filmThicknessUm = 3;              % um, same unit as the wavelengths

%% Discretisation -- substrate relDepth 20 (film at most 3/1.5 = 2) -> legOrder 50
legOrder            = 50;
relDepthPlaceholder = 1;           % film relDepth, overwritten per wavelength

ecFilm  = ElastConsts([cuC11Hill, cuC44Hill], isoSymmNo);
matSub  = Material.FromLibrary('Si', dirMatSi);
matFilm = Material(ecFilm, cuDensity);
cpSub   = CalcParams(relDepthSub, legOrder);
cpFilm  = CalcParams(relDepthPlaceholder, legOrder, filmThicknessUm);

structureName = 'si100_sub_cu_film';
nMaxModes     = 500;
structure = Structure(matSub, cpSub, structureName, nMaxModes);
structure = structure.AddLayer(matFilm, cpFilm);

%% Shared axes and vMax
anglesDeg     = (0:2:90)';
wavelengthsUm = (1.5:0.25:12)';     % um, evenly spaced; film relDepth 2 .. 0.25
vMax     = [];                     % optional: [] = Green axis ends at 1.05 x the fastest vL of any layer
calcMode = 'full';                 % 'ritzOnly': no Green map, trapped waves only

%% Part 1 -- one wavelength, many angles
wavelengthAngular = wavelengthsUm(end);                       % um, thinnest relative film
structureAtWav    = structure.ApplyWavelength(wavelengthAngular); % returns a NEW Structure
retentionAngular  = 'saw_only';    % viewer recomputes pruned angles on demand
sweepAngular = Calc.angSweepLayered(structureAtWav, anglesDeg, vMax, calcMode, retentionAngular);
sweepAngular = sweepAngular.calc();
velLimsAngular = [2.5 3];          % km/s, figure and viewer alike
Vis.plotAngularMap(sweepAngular);
ylim(velLimsAngular);
title(sprintf('Si(100)/%g \\mum Cu: SAW velocity vs. angle, %.1f \\mum', filmThicknessUm, wavelengthAngular));

guiAngular = Gui.combinedViewer(sweepAngular);
guiAngular.show();
guiAngular.setAxisRange('y', velLimsAngular);
guiAngular.setAxisRange('c', [-2 0]);

%% Part 2 -- one angle, many wavelengths, branch-tracked
angleDeg        = 0;
retentionCurve  = 'full';          % REQUIRED by DispersionBranches.track
sweepCurve = Calc.wavevecAngSweep2Layered(structure, angleDeg, wavelengthsUm, vMax, calcMode, retentionCurve);
sweepCurve = sweepCurve.calc();

curveType           = 'k-vphase';  % 6 options, see Vis.plotDispersionCurves help
angleIdx            = 1;           % index into the sweep's angles
showRitzBranches    = true;
showGreenBackground = true;        % needs calcMode='full' and a velocity y-axis
kAxisLimits         = round(2*pi ./ [max(wavelengthsUm), min(wavelengthsUm)], 3, 'significant');   % rad/um
velAxisLimits       = [1 6];       % km/s

gui = Gui.dispersionViewer(sweepCurve);
gui.show();
gui.curveTypeDropdown.Value   = curveType;
gui.showGreenBgCheckbox.Value = showGreenBackground;   % set both BEFORE onCurveTypeChanged
gui.onCurveTypeChanged();
gui.setAxisRange('x', kAxisLimits); gui.setAxisRange('y', velAxisLimits);

%% Static equivalent -- branches over the Green background
Vis.plotDispersionCurves(sweepCurve, curveType, angleIdx, showRitzBranches, showGreenBackground);
xlim(kAxisLimits); ylim(velAxisLimits);

%% Part 3 -- genuine double sweep: angle x wavelength in one calculation
anglesDouble      = (0:10:90)';                              % coarse: RAM scales with the grid
wavelengthsDouble = (1.5:0.5:12)';  % um
retentionDouble   = 'full';        % 'none' suffices for angular maps only
sweepDouble = Calc.wavevecAngSweep2Layered(structure, anglesDouble, wavelengthsDouble, vMax, calcMode, retentionDouble);
sweepDouble = sweepDouble.calc();

%% Angular |G_13,eff| map at one wavelength -- each sweeps{iWav} is a full angular sweep
scaleMode   = 'log';
climRange   = [-3, -1];
iWavShort   = 1;
Vis.plotGreenMap(sweepDouble.sweeps{iWavShort}, scaleMode, climRange);
ylim(velAxisLimits);
title(sprintf('Si(100)/%g \\mum Cu: |G_{13}^{eff}| at %.1f \\mum', filmThicknessUm, wavelengthsDouble(iWavShort)));

%% Dispersion curve at one angle of the same grid
curveTypeDouble = 'wavelength-vphase';
angleIdxDouble  = find(anglesDouble == 0);
Vis.plotDispersionCurves(sweepDouble, curveTypeDouble, angleIdxDouble, showRitzBranches, showGreenBackground, scaleMode, climRange);
wavAxisLimits = wavelengthsDouble([1 end])';   % um, the swept span
xlim(wavAxisLimits); ylim(velAxisLimits);

%% Interactive double-sweep viewer -- toggle angular map <-> dispersion curve
guiDouble = Gui.doubleSweepViewer(sweepDouble);
guiDouble.show();
guiDouble.setAxisRange('y', velAxisLimits);   % angular map at iWavShort, as its figure
guiDouble.setAxisRange('c', climRange);
guiDouble.onSwitchMode();   % open on the dispersion curve, same ranges as above
guiDouble.dispersionViewer.curveTypeDropdown.Value = curveTypeDouble;
guiDouble.dispersionViewer.showGreenBgCheckbox.Value = showGreenBackground;
guiDouble.dispersionViewer.onCurveTypeChanged();
guiDouble.setAxisRange('x', wavAxisLimits); guiDouble.setAxisRange('y', velAxisLimits);
clim(guiDouble.dispersionViewer.mapAxes, climRange);   % no colour control in the viewer
