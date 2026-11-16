%% Tutorial 08 (code) -- isoFreq scan: Ni film on Cu(001) at fixed real frequency, |G| vs k_parallel
% Narrated version: examples/tutorials/tutorial_08_isofreq_scan.m. Run guliwerr/init.m first.

%% Substrate -- Cu(001); treated analytically, its CalcParams are never read
dirMatSub   = [1 0 0; 0 0 1];
matSub      = Material.FromLibrary('Cu', dirMatSub);   % constants + density from the materials library
cpSub       = CalcParams();        % defaults, kept only so the Structure is well-formed

%% Film -- Ni(001), 100 nm
dirMatFilm      = [1 0 0; 0 0 1];
filmThicknessNm = 100;
matFilm         = Material.FromLibrary('Ni', dirMatFilm);   % materials library
relDepthPlaceholder = 1;           % overwritten by ApplyWavelength below
legOrder            = 50;
cpFilm = CalcParams(relDepthPlaceholder, legOrder, filmThicknessNm);

%% Structure, normalised to a reference wavelength
structureName = 'cu_sub_ni_film';
nMaxModes     = 500;
lambdaRefNm   = 300;               % free normalisation, ~ the expected SAW wavelength
structure = Structure(matSub, cpSub, structureName, nMaxModes);
structure = structure.AddLayer(matFilm, cpFilm);
structure = structure.ApplyWavelength(lambdaRefNm);

%% Fixed frequency -> normalised omega0, and a k_parallel axis from a velocity axis
freqGHz   = 5;
omega0    = Calc.freqToOmega0(freqGHz * 1e9, lambdaRefNm);   % lambdaRefNm must match ApplyWavelength
velMax    = 6;                     % km/s, brackets the expected bulk/SAW velocities
velMin    = 1;                     % km/s
nVel      = 300;
velAxis   = linspace(velMax, velMin, nVel)';   % descending -> kParRange ascending
kParRange = omega0 ./ velAxis;

%% Step 1 -- single angle: 1D |G_33| vs k_parallel
angleDeg   = 0;
singleFreq = Calc.isofreqSingleLayered(structure, angleDeg, omega0, kParRange);
singleFreq = singleFreq.calc();
Vis.plotDispersionMap(singleFreq);
ylim([-2.5 0])                     % normalised log: peak at 0
title(sprintf('Cu(001)/%g nm Ni, 0 deg, %g GHz: isoFreq spectrum', filmThicknessNm, freqGHz));

%% Step 2 -- angle sweep at the same frequency: 2D map
anglesDeg   = (0:2:90)';
useParallel = "auto";              % "auto" | "always" | "never"
sweep = Calc.isofreqAngSweepLayered(structure, anglesDeg, omega0, kParRange, useParallel);
sweep = sweep.calc();
Vis.plotDispersionMap(sweep);
clim([-3.5 -1.5])
axis tight
title(sprintf('Cu(001)/%g nm Ni, %g GHz: isoFreq |G_{33}| map', filmThicknessNm, freqGHz));

%% Step 3 -- several frequencies: one shared kParRange must span all of them
freqListGHz    = [3, 5, 7];
omega0List     = Calc.freqToOmega0(freqListGHz(:) * 1e9, lambdaRefNm);
nK             = 400;
kParRangeMulti = linspace(min(omega0List) / velMax, max(omega0List) / velMin, nK)';
sweepMulti = Calc.freqAngSweep2Layered(structure, anglesDeg, omega0List, kParRangeMulti, useParallel);
sweepMulti = sweepMulti.calc();
for omegaIdx = 1:numel(omega0List)
    Vis.plotDispersionMap(sweepMulti, omegaIdx);   % one frequency slice per figure
    clim([-3.5 -1.5])
    axis tight
    title(sprintf('Cu(001)/%g nm Ni, %g GHz: isoFreq |G_{33}| map', filmThicknessNm, freqListGHz(omegaIdx)));
end

%% Interactive frequency-slice viewer
gui = Gui.freqAngSweep2LayeredViewer(sweepMulti);
gui.setAxisRange('c', [-3.5 -1.5]);
gui.show();
