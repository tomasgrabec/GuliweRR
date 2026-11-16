%% Tutorial: isoFreq scan -- fixed real frequency, Green response vs k_parallel
%
% The isoFreq family asks: "at one fixed real frequency, which in-plane
% wavevectors k_parallel carry a guided mode?" The other layered tutorials
% are isoWavevec: at fixed |k_parallel| = 1, which velocities are supported?
% The two coincide only for a non-dispersive medium (docs/background.md).
%
% Green's function only: Ritz fixes |k_parallel| = 1 by construction, so it
% cannot hold the frequency fixed while k_parallel varies. calcMode and
% retentionLevel do not apply.
%
% Copy this script as a starting point for an isoFreq calculation.
%
% Run run('guliwerr/init.m') from the project root first (tutorial_01 says why).

%% Substrate -- Cu(001), cubic. The Green solver treats it as a semi-infinite
% half-space, so its CalcParams are never read; the defaults just keep the
% Structure valid (and usable for a Ritz calculation later).
% The Cu constants and density come from the materials library
% (Material.ListLibrary() shows values and sources).
dirMatSub = [1 0 0; ...      % surface normal
             0 0 1];         % propagation reference (Material's default)

matSub = Material.FromLibrary('Cu', dirMatSub);
cpSub  = CalcParams();

%% Film -- Ni(001), cubic, thin, on top (constants and density from the
% materials library)
dirMatFilm = [1 0 0; ...      % surface normal
              0 0 1];         % propagation reference (Material's default)

filmThicknessNm = 100;   % real film thickness [nm]
lambdaRefNm     = 300;   % reference wavelength [nm]: a free normalisation
                         % (any positive value gives the same physics); pick it
                         % near the expected SAW wavelength so the derived
                         % relDepth stays well-conditioned (help Calc.freqToOmega0)

matFilm = Material.FromLibrary('Ni', dirMatFilm);

relDepthFilmPlaceholder = 1;   % replaced by ApplyWavelength below
legO = 50;   % not read by the Green solver; a valid value for later Ritz use
cpFilm = CalcParams(relDepthFilmPlaceholder, legO, filmThicknessNm);

strName   = 'cu_sub_ni_film';
nMaxModes = 500;   % max number of Ritz modes (Ritz use only)
str = Structure(matSub, cpSub, strName, nMaxModes);   % layer 1 = substrate
str = str.AddLayer(matFilm, cpFilm);                   % layer 2 = film, on top
str = str.ApplyWavelength(lambdaRefNm);   % film relDepth = thickness/lambdaRefNm

%% Fixed real frequency -> GreenSolver's normalised omega0
freqGHz = 5;
omega0  = Calc.freqToOmega0(freqGHz * 1e9, lambdaRefNm);

%% k_parallel axis, built from a velocity axis bracketing the expected
% SAW/bulk velocities. velAxis descending -> kParRange = omega0./velAxis ascending.
nVel    = 300;
velAxis = linspace(6, 1, nVel)';   % km/s
kParRange = omega0 ./ velAxis;

%% Step 1: single angle -- 1D spectrum, |G_33| vs k_parallel
angleSingle = 0;   % deg
singleFreq = Calc.isofreqSingleLayered(str, angleSingle, omega0, kParRange);
singleFreq = singleFreq.calc();

Vis.plotDispersionMap(singleFreq);
ylim([-2.5 0])   % log10(|G|/max|G|): the peak sits at 0 by construction (default [-5 0])
title(sprintf('Cu(001)/%g nm Ni film, 0 deg, f = %g GHz -- isoFreq spectrum', filmThicknessNm, freqGHz));

%% Step 2: angle sweep at the same fixed frequency -- 2D map, angle x k_parallel
angles = (0:2:90)';
sweep  = Calc.isofreqAngSweepLayered(str, angles, omega0, kParRange);
sweep  = sweep.calc();

Vis.plotDispersionMap(sweep);
clim([-3.5 -1.5])
axis tight   % crop the empty margin around the map
title(sprintf('Cu(001)/%g nm Ni film, f = %g GHz -- isoFreq |G_{33}| map', filmThicknessNm, freqGHz));

%% Step 3: several frequencies at once -- Calc.freqAngSweep2Layered repeats
% Step 2 for every omega0 of a list. ONE kParRange serves every frequency, so
% it must span all of them: from omega0/v at the fastest velocity and lowest
% frequency to the slowest velocity at the highest.
freqListGHz   = [3, 5, 7];
omega0List    = Calc.freqToOmega0(freqListGHz(:) * 1e9, lambdaRefNm);
velMax        = 6;   % km/s
velMin        = 1;   % km/s
kParRangeMulti = linspace(min(omega0List)/velMax, max(omega0List)/velMin, 400)';
useParallel   = "auto";   % "auto" | "always" | "never"
sweepMulti = Calc.freqAngSweep2Layered(str, angles, omega0List, kParRangeMulti, useParallel);
sweepMulti = sweepMulti.calc();

for omegaIdx = 1:numel(omega0List)   % one frequency slice per figure
    Vis.plotDispersionMap(sweepMulti, omegaIdx);
    clim([-3.5 -1.5])
    axis tight
    title(sprintf('Cu(001)/%g nm Ni film, f = %g GHz -- isoFreq |G_{33}| map', ...
        filmThicknessNm, freqListGHz(omegaIdx)));
end

%% Interactive: step through the frequency slices with Green's-response
% fields, no recomputation. Needs an interactive MATLAB session.
gui = Gui.freqAngSweep2LayeredViewer(sweepMulti);
gui.setAxisRange('c', [-3.5 -1.5]);   % colour limits, kept across frequency steps
gui.show();

% Not shown: GreenPeaks.find extracts peak positions from (absG, kParRange);
% calc() does not pick peaks itself.
