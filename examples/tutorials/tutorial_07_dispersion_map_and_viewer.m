%% Tutorial: dispersion curves -- angle sweep and wavelength sweep, two viewers
%
% Three views of one film-on-substrate pair, in real units:
%   Part 1 -- ONE wavelength, MANY angles: SAW velocity vs. angle
%             (Vis.plotAngularMap) and the combined viewer.
%   Part 2 -- ONE angle, MANY wavelengths: branch-tracked dispersion curves
%             in the dispersion viewer.
%   Part 3 -- MANY angles x MANY wavelengths in one double sweep: angular
%             Green maps at fixed wavelengths, a dispersion curve at a fixed
%             angle, and the double-sweep viewer switching between the two.
%
% Si(100) substrate / sputtered Cu film: the Cu is polycrystalline and assumed untextured, so
% isotropic, and all angular variation comes from the cubic Si.
%
% Copy this script as a starting point for an angle sweep (Part 1) or a
% wavelength sweep (Part 2) on a film-on-substrate pair.
%
% Run run('guliwerr/init.m') from the project root first (tutorial_01 says why).

%% Substrate -- Si(100), cubic, from the materials library (constants and
% density with their source: Material.ListLibrary())
dirMatSi = [1 0 0; ...                % surface normal ([100])
            0 1 0];                   % propagation reference ([010]; any
                                       % <100> in-plane axis is equivalent)
relDepthSub = 20;   % substrate depth, in wavelengths; unlike the film's, it
                    % does not change with the swept wavelength

%% Film -- Cu, isotropic (sputtered, so no crystal orientation and no dirMat):
% the Hill average of the single-crystal Cu constants (materials library, as
% in the other tutorials). Hill is the mean of the Voigt (uniform strain,
% upper) and Reuss (uniform stress, lower) bounds, the usual practical choice
% for a polycrystalline film; for Cu it gives a clearly softer shear constant
% (C44 47.0 vs 54.1 GPa Voigt), so it matters here. tutorial_S1 compares the
% three averages.
matCuCrystal = Material.FromLibrary('Cu');
ecCuCubic    = matCuCrystal.elConsts;
cuRho        = matCuCrystal.rho;   % g/cm^3
[cuC11Hill, cuC44Hill] = ecCuCubic.GetHillAverages();   % -> [M, G], isotropic pair
filmThicknessUm = 3;   % real film thickness [um]; divided by each wavelength
                       % (same unit) to get the film relDepth

%% Discretisation -- legO must stay >= 2x the largest relDepth in the sweep.
% The largest is the substrate's 20, so legO = 50 clears the rule; the film
% reaches only relDepth 3/1.5 = 2 at the shortest wavelength (1.5 um).
legO = 50;

matSub = Material.FromLibrary('Si', dirMatSi);
cpSub  = CalcParams(relDepthSub, legO);

ecFilm  = ElastConsts([cuC11Hill, cuC44Hill], 1);   % isotropic class -- no dirMat needed
matFilm = Material(ecFilm, cuRho);
relDepthFilmPlaceholder = 1;   % replaced at every wavelength of the sweeps
cpFilm = CalcParams(relDepthFilmPlaceholder, legO, filmThicknessUm);

strName   = 'si100_sub_cu_film';
nMaxModes = 500;   % max number of Ritz modes computed/retained
str = Structure(matSub, cpSub, strName, nMaxModes);   % layer 1 = substrate
str = str.AddLayer(matFilm, cpFilm);                   % layer 2 = film, on top
% No ApplyWavelength here: Calc.wavevecAngSweep2Layered applies each
% wavelength to its own copy of str.

%% ===================================================================
%% Part 1 -- fixed layer stack, SAW velocity vs. angle (Calc.angSweepLayered)
%% ===================================================================
%
% The layered counterpart of Calc.angSweepHS (tutorial_02): the film
% relDepth is fixed by one wavelength, and only the angle is swept.
% Part 1 uses the longest wavelength of the vector Part 2 sweeps.
angles      = (0:2:90)';                             % deg
wavelengths = (1.5:0.25:12)';   % um, evenly spaced; film relDepth 2 .. 0.25
% At relDepth 2 the film is thick enough that the substrate no longer shifts
% the Rayleigh wave (its velocity is within 1e-4 of a Cu half-space's).

%% vMax = [] (optional): the Green axis ends at 1.05 x the fastest
% longitudinal limiting velocity of any layer; a larger value only extends
% the displayed map. calcMode 'full', not 'ritzOnly': 'ritzOnly' has no Green
% map and reports trapped waves only (no leaky branch above the substrate's
% shear limit), while sweep.classify('ritz') on a 'full' result still gives
% the Ritz-only entries at no extra cost.
vMax     = [];
calcMode = 'full';

wavAngular = wavelengths(end);   % um; longest wavelength -> thinnest film (relDepth 0.25)
strAtWav   = str.ApplyWavelength(wavAngular);   % returns a NEW Structure; str is unchanged
retentionLevelAng = 'saw_only';   % enough: the viewer below recomputes a
                                  % pruned angle's modal shapes on demand
sweepAng = Calc.angSweepLayered(strAtWav, angles, vMax, calcMode, retentionLevelAng);
sweepAng = sweepAng.calc();

velLimsAng = [2.5 3];     % km/s; the same range in the figure and the viewer below
Vis.plotAngularMap(sweepAng);   % creates its own figure/axes when none is passed
ylim(velLimsAng);
title(sprintf('Si(100)/%g \\mum Cu film -- SAW velocity vs. angle, wavelength = %.1f \\mum', filmThicknessUm, wavAngular));

%% Interactive: the combined viewer pairs modal shapes with the same
% angular map; click to inspect. Needs an interactive MATLAB session.
guiAng = Gui.combinedViewer(sweepAng);
guiAng.show();
guiAng.setAxisRange('y', velLimsAng);   % sets the viewer's own range fields too
guiAng.setAxisRange('c', [-2 0]);   % colour: the two branches stand apart from the background

%% ===================================================================
%% Part 2 -- ONE angle from Part 1, followed as branch-tracked curves
%% ===================================================================

angle          = 0;      % deg
retentionLevel = 'full'; % REQUIRED: branch tracking compares every eigenvector
                         % of consecutive wavelength steps
calcModeCurve  = 'full'; % not 'ritzOnly': each wavelength step keeps its Green
                         % map, needed for the Green background below
sweep1D = Calc.wavevecAngSweep2Layered(str, angle, wavelengths, vMax, calcModeCurve, retentionLevel);
sweep1D = sweep1D.calc();

branchModeNoMap = DispersionBranches.track(sweep1D);   % Ritz mode number of every branch at every step

%% The dispersion viewer -- 'combined' (default) is one window; 'split' gives
% two. Click a curve, or step Wavelength/Branch, to see that mode's
% displacement field. Needs an interactive MATLAB session.
gui = Gui.dispersionViewer(sweep1D);
gui.show();

%% k-vphase with a fixed axis range (the viewer opens on 'wavelength-freq',
% autoscaled), over the Green background -- available for the two v_phase
% curve types on a 'full' sweep -- to see which tracked branch follows a
% Green peak. Set both controls before the one onCurveTypeChanged() call:
% each call autoscales mapAxes again.
gui.curveTypeDropdown.Value = 'k-vphase';
gui.showGreenBgCheckbox.Value = true;
gui.onCurveTypeChanged();
kLims   = round(2*pi ./ [max(wavelengths), min(wavelengths)], 3, 'significant');   % k [rad/um], the swept span
velLims = [1 6];                          % km/s
gui.setAxisRange('x', kLims);   % as typed into the viewer's range fields
gui.setAxisRange('y', velLims);

%% The same view as a static plot
Vis.plotDispersionCurves(sweep1D, 'k-vphase', [], true, true);   % angleIdx [] = 1; branches; Green background
xlim(kLims);
ylim(velLims);
title(sprintf('Si(100)/%g \\mum Cu film -- k-dispersion, angle = 0 deg, with Green''s background', filmThicknessUm));

%% "Show Ritz branches" off: the Green background alone (static equivalent:
% Vis.plotDispersionCurves(sweep1D, 'k-vphase', [], false, true))
gui.showRitzBranchesCheckbox.Value = false;
gui.onCurveTypeChanged();
gui.setAxisRange('x', kLims);
gui.setAxisRange('y', velLims);

%% ===================================================================
%% Part 3 -- genuine double sweep: angle x wavelength in one calculation
%% ===================================================================
%
% Calc.wavevecAngSweep2Layered with several angles AND several wavelengths
% fills the whole grid: sweep2D.sweeps{iw} is a complete angular sweep at
% wavelength iw, so every Part 1 plot works on a slice. Under
% retentionLevel 'full' RAM scales with nAngles x nWavelengths -- keep the
% grid coarse (docs/pitfalls.md NUM4); 'full' is needed only for the
% dispersion curves (branch tracking).
anglesDouble      = (0:10:90)';                          % deg, coarse
wavelengthsDouble = (1.5:0.5:12)';   % um, coarse
retentionDouble   = 'full';   % 'none' is enough if only angular maps are wanted
sweep2D = Calc.wavevecAngSweep2Layered(str, anglesDouble, wavelengthsDouble, vMax, ...
    calcMode, retentionDouble);
sweep2D = sweep2D.calc();

%% Angular |G_13,eff| maps at the shortest and longest wavelength (thickest vs.
% thinnest relative film) -- each slice is an ordinary angular sweep result.
% At 1.5 um (film 2 wavelengths thick) the Rayleigh wave is that of a Cu
% half-space, flat in angle; the extra lines are standing waves across the
% film thickness. At 12 um the Si below shows through in the angle dependence.
iWavShort = 1;
iWavLong  = numel(wavelengthsDouble);
climRange = [-3, -1];   % one decade above the [-4 -1] default: the film branches stand out
velLimsGreen = {velLims, [1 Inf]};   % km/s; the long wavelength reaches faster, Si-like velocities
iWavShow = [iWavShort, iWavLong];
for iFig = 1:numel(iWavShow)
    iWav = iWavShow(iFig);
    Vis.plotGreenMap(sweep2D.sweeps{iWav}, 'log', climRange);
    ylim(velLimsGreen{iFig});
    title(sprintf('Si(100)/%g \\mum Cu film -- |G_{13}^{eff}|, wavelength = %.1f \\mum', filmThicknessUm, wavelengthsDouble(iWav)));
end

%% Dispersion curve at one fixed angle of the SAME grid -- angleIdx picks it
angleIdxDouble = find(anglesDouble == 0);   % index into anglesDouble, not an angle in deg
wavLimsDouble  = wavelengthsDouble([1 end])';   % um, the swept span
Vis.plotDispersionCurves(sweep2D, 'wavelength-vphase', angleIdxDouble, true, true, 'log', climRange);
title(sprintf('Si(100)/%g \\mum Cu film -- wavelength dispersion at 0 deg, from the double sweep', filmThicknessUm));
xlim(wavLimsDouble);
ylim(velLims);

%% Interactive: the double-sweep viewer -- a button switches between the
% angular map (at a wavelength you step through) and the dispersion curve (at
% a fixed angle), keeping each side's axis ranges. Needs an interactive session.
guiDouble = Gui.doubleSweepViewer(sweep2D);
guiDouble.show();
% setAxisRange acts on the mode currently shown, and each mode keeps its own
% ranges. The viewer opens on the angular map at the shortest wavelength:
% give it the ranges of that wavelength's Green-map figure above.
guiDouble.setAxisRange('y', velLimsGreen{1});
guiDouble.setAxisRange('c', climRange);
% Then open on the dispersion curve, as in the figure above (the button
% switches back to the angular map). The dispersion viewer has no colour
% control, so its Green background takes climRange directly on the axes
% (a display toggle there redraws it at the default range).
guiDouble.onSwitchMode();
guiDouble.dispersionViewer.curveTypeDropdown.Value = 'wavelength-vphase';
guiDouble.dispersionViewer.showGreenBgCheckbox.Value = true;
guiDouble.dispersionViewer.onCurveTypeChanged();
guiDouble.setAxisRange('x', wavLimsDouble);
guiDouble.setAxisRange('y', velLims);
clim(guiDouble.dispersionViewer.mapAxes, climRange);
