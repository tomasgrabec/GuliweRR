%% Tutorial: free plate -- Lamb-wave dispersion, mode types and ZGV points
%
% A free-standing PLATE instead of a surface wave on a semi-infinite
% substrate: both faces traction-free, finite thickness -- the Lamb-wave
% problem (Grabec, Sedlák, Seiner 2020, Wave Motion 96, 102567). Same Ritz
% method, different boundary conditions (docs/reference/architecture/plate.md):
%   - no Green's function and no calcMode argument: only the Ritz solve runs;
%   - every computed mode is a physical branch, so there is no "which mode is
%     the SAW" question; each mode gets a TYPE instead: Asymm (A0-like), Symm
%     (S0-like), TH-Asymm, TH-Symm (shear-horizontal-like);
%   - retentionLevel is 'full' | 'none' only.
%
% Physical scaling: pass the plate's real thickness and a real wavevector k
% in CONSISTENT units -- here mm and rad/mm. Velocities stay in km/s (= mm/us),
% so frequency f = v*k/(2*pi) comes out in MHz.
%
%   Part 1 -- one angle, one wavevector: the mode list itself
%   Part 2 -- one angle, many wavevectors: dispersion curves, branch
%             classification, zero-group-velocity (ZGV) points
%   Part 3 -- many angles, one wavevector: angular map + combined viewer
%   Part 4 -- a coated (two-layer, NON-mirror-symmetric) plate: when mode-type
%             labels stop being rigorous
%   Part 5 -- a measured reference: the S1 ZGV resonances of a Si(001) wafer,
%             compared with published measurements
%
% Copy this script as a starting point for a plate calculation; every call is
% listed in docs/reference/functionReference.md, "Free plate (Lamb waves)".
%
% Run run('guliwerr/init.m') from the project root first (tutorial_01 says why).

%% Material -- Cu(001), cubic: an anisotropic plate with ZGV points for
% Part 2 to find. Constants and density come from the materials library
% (Material.ListLibrary() shows values and sources).
dirMat = [1 0 0; ...         % plate normal
          0 0 1];            % in-plane propagation reference (Material's default)
mat    = Material.FromLibrary('Cu', dirMat);

%% Discretisation -- for a plate, relDepth is the thickness in wavelengths
% (k*d/(2*pi)), set at every point from realThickness and realWavVec below;
% the value here is a placeholder. legO must resolve the thickest relDepth
% reached, kMax*d/(2*pi) = 20/(2*pi) ~ 3.2; legO = 20 and 30 give the same
% branches below 6 MHz (Part 2), so 20 is ample there.
relDepthPlaceholder = 0.2;
legO = 20;   % Legendre basis order across the thickness
cp   = CalcParams(relDepthPlaceholder, legO);

strName   = 'cu_001_plate';
nMaxModes = 500;
str = Structure(mat, cp, strName, nMaxModes);   % nLay = 1: homogeneous plate

realThickness = 1;   % mm -- plate thickness, FIXED throughout this script

%% ===================================================================
%% Part 1 -- one angle, one wavevector: the mode list
%% ===================================================================
angle      = 0;     % deg
realWavVec = 2;     % rad/mm -> wavelength 2*pi/k ~ 3.1 mm, k*d = 2
retentionLevel = 'full';   % 'full' (keep eigenvectors: modal shapes, tracking) | 'none'

singleResult = Calc.singlePlate(str, angle, realThickness, realWavVec);
singleResult = singleResult.calcAll(retentionLevel);
ritz   = singleResult.ritzSolver;

modeTypeNames = {'Asymm', 'Symm', 'TH-Asymm', 'TH-Symm'};   % modesType codes 1..4
nShow = 8;
fprintf('\nCu(001) plate, d = %g mm, k = %g rad/mm, angle = %g deg\n', realThickness, realWavVec, angle);
fprintf('%5s %10s %10s %10s %10s %11s\n', 'mode', 'v [km/s]', 'vg [km/s]', 'f [MHz]', 'type', 'confidence');
for iMode = 1:nShow
    fprintf('%5d %10.4f %10.4f %10.4f %10s %11.2f\n', iMode, ritz.vels(iMode), ritz.velsGroup(iMode), ...
        ritz.freqsPhysical(iMode), modeTypeNames{ritz.modesType(iMode)}, ritz.modeTypeConfidence(iMode));
end
% The type comes from mid-plane parity: each displacement component is split
% into its even and odd parts about the plate's mid-plane, and the mode goes
% to the class holding most of its energy (Symm: u1 even + u3 odd, Asymm:
% u1 odd + u3 even, TH-Symm: u2 even, TH-Asymm: u2 odd).
% modeTypeConfidence = largest fraction - runner-up, in [0, 1]: near 1 =
% pure (every mode here, along the cubic axis), low = mixed: SH-Lamb
% coupling within one parity class off the axes, parity mixing on an
% asymmetric stack (Part 4)

%% 1D view: every mode at this one point, coloured by type
criterion       = 'modeType';   % 'modeType' (categorical) | 'confidence' (continuous)
velCutoffFactor = 2.5;          % x the longitudinal velocity: the 8 modes of the table above
                                % ([] = default 1.1, a surface-wave habit that hides most plate modes)
Vis.plotRitzModeCurvePlate(singleResult, criterion, velCutoffFactor);
title(sprintf('Cu(001) plate, k = %g rad/mm: Ritz modes by type', realWavVec));

%% ===================================================================
%% Part 2 -- one angle, many wavevectors: dispersion curves and ZGV
%% ===================================================================
% Calc.wavevecAngSweep2Plate sweeps the real wavevector at a fixed
% realThickness. Dense steps help branch tracking, which matches
% CONSECUTIVE steps by eigenvector overlap (MAC). The span matters as much as
% the density: a branch is drawn only over the k range swept, so the low
% branches need k up to ~20 rad/mm to reach 6 MHz in the group-velocity plot.
wavVecs         = (0.1:0.1:20)';   % rad/mm
retentionSweep  = 'full';     % REQUIRED for branch tracking, classifyBranches and locateZGVs
useParallel     = "auto";     % "auto" | "always" | "never"
sweepDisp = Calc.wavevecAngSweep2Plate(str, angle, wavVecs, realThickness, retentionSweep, useParallel);
sweepDisp = sweepDisp.calc();

%% Dispersion curves -- one tracked line per branch. curveType picks the axes:
% 'wavelength-freq' (default) | 'wavelength-vphase' | 'k-freq' | 'k-vphase' |
% 'freq-vphase' | 'freq-vgroup'. 'k' is wavVecs directly.
angleIdx              = 1;        % which angle of the sweep (only one here)
colorMode             = 'type';   % 'type' (mode type where trusted) | 'branch' (branch identity)
confidenceThresh      = 0.3;      % default; below it a label falls back to 'branch N' (unvalidated)
showChristoffelBounds = true;     % bulk shear/longitudinal reference lines (velocity y-axes only)
velLimsDisp           = [0 8];    % km/s, phase-velocity axis (here and in the viewer below)
freqLimsDisp          = [0 6];    % MHz; higher branches are less resolved and crowd the plot
Vis.plotDispersionCurvesPlate(sweepDisp, 'k-vphase', angleIdx, colorMode, confidenceThresh, showChristoffelBounds);
ylim(velLimsDisp);
xlabel('k [rad/mm]');   % the plot cannot know the length unit: units follow realThickness and wavVecs
title('Cu(001) plate: phase velocity vs. k');

Vis.plotDispersionCurvesPlate(sweepDisp, 'k-freq', angleIdx, colorMode, confidenceThresh);
ylim(freqLimsDisp);
xlabel('k [rad/mm]');
ylabel('frequency [MHz]');   % km/s x rad/mm / (2*pi) = MHz
title('Cu(001) plate: frequency vs. k');

Vis.plotDispersionCurvesPlate(sweepDisp, 'freq-vgroup', angleIdx, colorMode, confidenceThresh);
xlim(freqLimsDisp);
xlabel('frequency [MHz]');
yline(0, 'k:');   % a branch crossing zero group velocity at finite k = a ZGV point
title('Cu(001) plate: group velocity vs. frequency');

%% Branch classification -- Cu(001) along [100] is symmetric under both the
% mid-plane mirror (exact parity: Symm/Asymm) and the sagittal mirror (exact
% Lamb/TH split), so every label is rigorous here. Off the symmetry axes
% (e.g. 30 deg) only the parity is exact, and a cut like (111) has neither;
% Part 4 shows a stack without the mid-plane mirror
classResult = sweepDisp.classifyBranches(angleIdx, confidenceThresh);
fprintf('\nMid-plane mirror: %d, sagittal mirror: %d (both: labels rigorous regardless of confidence)\n', ...
    classResult.isExactMirrorSymmetric, classResult.isExactSagittalSymmetric);
fprintf('Branches 1-6 at the first k: %s\n', strjoin(classResult.displayLabel(1, 1:6), ', '));

%% Zero-group-velocity points -- a scan along each tracked branch for a sign
% change of the group velocity, then a Newton refinement that solves the plate
% eigenproblem and the zero-group-velocity condition together (no mode number
% is held fixed, so a branch crossing nearby does not matter). A point is kept
% only where an independent solve confirms that the group velocity vanishes
maxCount = Inf;   % cap on returned points (Inf = all)
zgvList  = sweepDisp.locateZGVs(angleIdx, maxCount);
zgvShow  = find([zgvList.freqZGV] <= freqLimsDisp(2));   % the plotted frequency range
fprintf('\nZGV points found: %d, of them below %g MHz:\n', numel(zgvList), freqLimsDisp(2));
for iZgv = zgvShow
    fprintf('  branch %2d: k = %.4f rad/mm, f = %.4f MHz, true ZGV: %d\n', ...
        zgvList(iZgv).branchNo, zgvList(iZgv).wavVecZGV, zgvList(iZgv).freqZGV, ...
        zgvList(iZgv).isTrueZGV);
end
% Check higher-branch points for convergence (raise legO, refine wavVecs)
% before trusting them: high-order modes are the least resolved. Part 5 checks
% the same search against published measurements on a silicon wafer.
%
% The search runs along k at ONE fixed angle, so it finds where the group
% velocity along k vanishes. Along a symmetry direction of the plate ([100]
% here) the group velocity is parallel to k, and that is a true ZGV point.
% Along any other direction the group velocity also has a component
% perpendicular to k (velGroupTrans, from ritzSolver.velsGroupTrans), which
% need not vanish there; isTrueZGV says whether it does (Part 5).

%% Interactive: dispersion viewer -- click a curve (or step wavevector/branch)
% to see that mode's displacement field through the thickness.
% Needs an interactive MATLAB session.
windowMode = 'combined';   % 'combined' (one window) | 'split' (two)
guiDisp = Gui.dispersionViewerPlate(sweepDisp, windowMode);
guiDisp.show();
% the viewer opens on 'wavelength-freq', autoscaled: switch to the k-vphase
% view of the first figure, then fix its range (each curve-type change
% autoscales again, so the range comes last)
guiDisp.curveTypeDropdown.Value = 'k-vphase';
guiDisp.onCurveTypeChanged();
guiDisp.setAxisRange('y', velLimsDisp);   % as typed into the viewer's range fields

%% ===================================================================
%% Part 3 -- many angles, one wavevector: angular dependence
%% ===================================================================
angles = (0:2:90)';   % deg, [100] -> [010] on the (001) plate
sweepAng = Calc.angSweepPlate(str, angles, realThickness, realWavVec, retentionLevel, useParallel);
sweepAng = sweepAng.calc();
% sweepAng.velsMap / .velsGroupMap / .modesTypeMap / .freqsPhysicalMap are
% dense [nModes x nAngles] arrays -- every angle has the same mode count.

% Mode type and its confidence side by side in one window: pass each plot an
% axes of a tiledlayout instead of letting it open its own figure
velLimsAng = [1 5.5];   % km/s; the same range in both maps and the viewer below
figure();
tiledlayout(1, 2);
axType = nexttile;
Vis.plotRitzModeMapPlate(sweepAng, 'modeType', [], axType);
ylim(axType, velLimsAng);
title(axType, sprintf('Cu(001) plate, k = %g rad/mm: by type', realWavVec));
axConf = nexttile;
Vis.plotRitzModeMapPlate(sweepAng, 'confidence', [], axConf);
ylim(axConf, velLimsAng);
title(axConf, 'mode-type confidence');

%% Interactive: combined viewer -- angular map (bottom) linked to the modal
% shape (top); it also takes a single-angle result (Part 1's singleResult)
guiAng = Gui.combinedViewerPlate(sweepAng);
guiAng.show();
guiAng.setAxisRange('y', velLimsAng);   % sets the viewer's own range fields too

%% ===================================================================
%% Part 4 -- coated plate: a thin isotropic Al layer on ONE face
%% ===================================================================
% A different material on one face removes the mid-plane mirror symmetry:
% symmetric and antisymmetric motion mix, and a mode-type label only names
% the dominant parity class, trusted where confidence >= confidenceThresh
% (a wave bound to one face is half even, half odd: confidence ~0). On four
% multi-layer test stacks the frequencies, group velocities and ZGV points
% agree with GEWtool (Kiefer, open source) to round-off; the confidenceThresh
% default is not validated (docs/reference/architecture/plate.md).
%
% realThickness is the TOTAL thickness; every layer's relDepth is scaled by
% the same factor, so the CalcParams relDepths below set the thickness ratio.
% Isotropic (polycrystalline) Al = the Hill average of the single-crystal
% constants from the materials library (tutorial_S1 section 9 shows the averages).
matAlCrystal   = Material.FromLibrary('Al');
[alC11, alC44] = matAlCrystal.elConsts.GetHillAverages();
alConsts = [alC11, alC44];   % [C11, C44] GPa, isotropic
alSymm   = 1;
alRho    = matAlCrystal.rho; % g/cm^3
matAl    = Material(ElastConsts(alConsts, alSymm), alRho);

coatingFraction = 0.1;    % coating = 10% of the total thickness
cpCu = CalcParams((1 - coatingFraction) * relDepthPlaceholder, legO);
cpAl = CalcParams(coatingFraction * relDepthPlaceholder, legO);
strCoated = Structure(mat, cpCu, 'cu_plate_al_coating', nMaxModes);   % layer 1 = bottom face
strCoated = strCoated.AddLayer(matAl, cpAl);                          % layer 2 = top face

fprintf('\nCoated plate mirror-symmetric: %d\n', strCoated.isMirrorSymmetric());

sweepCoated = Calc.wavevecAngSweep2Plate(strCoated, angle, wavVecs, realThickness, retentionSweep, useParallel);
sweepCoated = sweepCoated.calc();

classCoated = sweepCoated.classifyBranches(angleIdx, confidenceThresh);
fprintf('Coated plate: %.0f%% of (step, branch) labels trusted at confidenceThresh = %g\n', ...
    100 * mean(classCoated.trustLabelByBranch(:)), confidenceThresh);

% 'type' colours by trusted label; 'branch' is the escape hatch when labels are not trusted
Vis.plotDispersionCurvesPlate(sweepCoated, 'k-vphase', angleIdx, 'type', confidenceThresh);
ylim(velLimsDisp);
xlabel('k [rad/mm]');
title('Cu plate + 10% Al coating: phase velocity vs. k (coloured by type)');

Vis.plotDispersionCurvesPlate(sweepCoated, 'k-vphase', angleIdx, 'branch');
ylim(velLimsDisp);
xlabel('k [rad/mm]');
title('Cu plate + 10% Al coating: phase velocity vs. k (coloured by branch)');

%% ===================================================================
%% Part 5 -- a measured reference: S1 ZGV resonances of a Si(001) wafer
%% ===================================================================
% Prada, Clorennec, Murray & Royer (2009, J. Acoust. Soc. Am. 126, 620) and
% Kiefer, Mezil & Prada (2023, Sci. Adv. 9, eadk6846) excited a [001]-cut
% silicon wafer with a pulsed laser and measured the S1 ZGV resonances with an
% interferometer. Kiefer et al. give everything needed to repeat the
% calculation: thickness 524.6 um, C11 = 165.6, C12 = 63.9, C44 = 79.5 GPa,
% density 2330 kg/m^3. The ZGV points sit on the principal directions:
%   ZGV1 on <100>: computed 7.7079 MHz at k = 3.421 rad/mm, measured 7.7072 MHz
%   ZGV2 on <110>: computed 7.7551 MHz at k = 3.142 rad/mm, measured 7.7545 MHz
% (computed with their own semi-analytical code, GEWtool). Prada et al. give
% the same points as f*d = 4.044 and 4.069 MHz*mm (computed with DISPERSE).
% Kiefer et al.'s own Si constants, kept on purpose to repeat their calculation;
% the materials library's Si (McSkimin & Andreatch 1964) differs by about 0.1 %.
siConsts = [165.6, 63.9, 79.5];   % [C11, C12, C44] GPa
siRho    = 2.330;                  % g/cm^3
siDirMat = [0 0 1; ...             % plate normal [001]
            1 0 0];                % angle 0 = [100], angle 45 = [110]
matSi = Material(ElastConsts(siConsts, 2), siRho, siDirMat);
strSi = Structure(matSi, CalcParams(relDepthPlaceholder, legO), 'si_001_wafer', nMaxModes);

siThickness = 0.5246;               % mm
siAngles    = [0; 45; 18.4];        % deg: [100], [110], and 26.6 deg from [110]
siWavVecs   = linspace(1, 6, 51)';  % rad/mm, around the S1 minimum near 3.3 rad/mm
sweepSi = Calc.wavevecAngSweep2Plate(strSi, siAngles, siWavVecs, siThickness, retentionSweep, useParallel);
sweepSi = sweepSi.calc();

publishedComputed = [7.7079, 7.7551, 7.727];   % MHz, Kiefer et al. 2023
publishedMeasured = [7.7072, 7.7545, NaN];     % MHz, Kiefer et al. 2023, Fig. 8
fprintf('\nSi(001) wafer, d = %.4f mm: lowest point of zero group velocity along k, per direction\n', siThickness);
fprintf('%9s %12s %12s %12s %12s %6s %16s\n', 'direction', 'k [rad/mm]', 'f [MHz]', 'published', ...
    'measured', 'type', 'c_gtheta [mm/us]');
dirNames  = {'[100]', '[110]', '18.4 deg'};
typeNames = {'TGV', 'ZGV'};   % isTrueZGV false / true
for iAng = 1:numel(siAngles)
    zgvSi = sweepSi.locateZGVs(iAng, Inf);
    [~, iLow] = min([zgvSi.freqZGV]);   % the S1 point; the A2 point is near 15 MHz
    z = zgvSi(iLow);
    fprintf('%9s %12.4f %12.4f %12.4f %12.4f %6s %16.4f\n', dirNames{iAng}, z.wavVecZGV, ...
        z.freqZGV, publishedComputed(iAng), publishedMeasured(iAng), typeNames{z.isTrueZGV + 1}, ...
        z.velGroupTrans);
end
% GuliweRR reproduces the published ZGV frequencies to every digit given, the
% TGV frequency to 1 kHz, and the measured ones to 0.01 %. Prada et al. used
% 2329 kg/m^3, which raises every frequency by sqrt(2330/2329), about 0.02 %;
% that is how their f*d = 4.044 and 4.069 MHz*mm come out.
%
% The two directions differ by about 50 kHz. ZGV1 is a minimum of the dispersion
% surface frequency(kx, ky); ZGV2 is a saddle point of it: a maximum across the
% angle, a minimum along k. Between the two axes, locateZGVs still finds a
% minimum along k at every angle, but there the group velocity has a component
% perpendicular to k, so the energy does not stay put: isTrueZGV is false.
% Kiefer et al. call these TGV (transverse-group-velocity) waves and measured
% them leaving a line source sideways. At 26.6 deg from [110] they computed
% 7.727 MHz and a transverse group velocity of -0.177 mm/us, and measured
% +-0.170 mm/us. Their angle runs from [110] towards [010], so that direction is
% 71.6 deg in GuliweRR's angle; 18.4 deg, used here, is its mirror image across
% [110]: an equivalent wave with the transverse component reversed. GuliweRR
% gives +0.178 mm/us here and -0.178 mm/us at 71.6 deg. Only the eight points on
% the principal directions are true ZGV resonances, and only they persist after
% about 40 us in the measurement.
