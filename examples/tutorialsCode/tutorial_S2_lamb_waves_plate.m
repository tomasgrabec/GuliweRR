%% Tutorial S2 (code) -- free plate (Lamb waves): modes, dispersion, mode types, ZGV, coated plate, Si(001) ZGV check
% Narrated version: examples/tutorials/tutorial_S2_lamb_waves_plate.m. Run guliwerr/init.m first.
% Units: thickness in mm, wavevector in rad/mm, velocity in km/s -> frequency in MHz.

%% Material -- Cu(001) plate
dirMat   = [1 0 0; 0 0 1];          % plate normal; propagation reference
material = Material.FromLibrary('Cu', dirMat);   % constants + density from the materials library

%% Discretisation -- relDepth is rescaled from realThickness * k at every point
relDepthPlaceholder = 0.2;         % overwritten; only the layer RATIO matters
legOrder   = 20;                   % plate up to ~3.2 wavelengths thick; converged below 6 MHz
calcParams = CalcParams(relDepthPlaceholder, legOrder);
nMaxModes  = 500;
structure  = Structure(material, calcParams, 'cu_001_plate', nMaxModes);
realThickness = 1;                 % mm, total plate thickness

%% Part 1 -- one angle, one wavevector: mode list
angleDeg       = 0;
realWavVec     = 2;                % rad/mm, k*d = 2
retentionLevel = 'full';           % 'full' | 'none' (no 'saw_only' for plates)
singleResult = Calc.singlePlate(structure, angleDeg, realThickness, realWavVec);
singleResult = singleResult.calcAll(retentionLevel);
ritzSolver   = singleResult.ritzSolver;

modeTypeNames = {'Asymm', 'Symm', 'TH-Asymm', 'TH-Symm'};   % modesType codes 1-4
nShow = 8;
for iMode = 1:nShow
    fprintf('%2d: v %.4f, vg %.4f km/s, f %.4f MHz, %-8s (confidence %.2f)\n', iMode, ...
        ritzSolver.vels(iMode), ritzSolver.velsGroup(iMode), ritzSolver.freqsPhysical(iMode), ...
        modeTypeNames{ritzSolver.modesType(iMode)}, ritzSolver.modeTypeConfidence(iMode));
end

criterion       = 'modeType';      % 'modeType' | 'confidence'
velCutoffFactor = 2.5;             % x longitudinal velocity ([] = default 1.1 hides most plate modes)
Vis.plotRitzModeCurvePlate(singleResult, criterion, velCutoffFactor);
title(sprintf('Cu(001) plate, k = %g rad/mm: modes by type', realWavVec));

%% Part 2 -- one angle, many wavevectors: dispersion curves
wavVecs        = (0.1:0.1:20)';    % rad/mm, dense for branch tracking, wide for whole branches
retentionSweep = 'full';           % REQUIRED for tracking, classifyBranches, locateZGVs
useParallel    = "auto";           % "auto" | "always" | "never"
sweepDisp = Calc.wavevecAngSweep2Plate(structure, angleDeg, wavVecs, realThickness, retentionSweep, useParallel);
sweepDisp = sweepDisp.calc();

angleIdx              = 1;
colorMode             = 'type';    % 'type' | 'branch'
confidenceThresh      = 0.3;       % default, unvalidated
showChristoffelBounds = true;      % bulk-velocity reference lines
curveTypeVphase = 'k-vphase';      % 'wavelength-freq' | 'wavelength-vphase' | 'k-freq' | 'k-vphase' | 'freq-vphase' | 'freq-vgroup'
curveTypeFreq   = 'k-freq';
curveTypeGroup  = 'freq-vgroup';
velLimsDisp     = [0 8];           % km/s
freqLimsDisp    = [0 6];           % MHz
Vis.plotDispersionCurvesPlate(sweepDisp, curveTypeVphase, angleIdx, colorMode, confidenceThresh, showChristoffelBounds);
ylim(velLimsDisp);
xlabel('k [rad/mm]');              % units follow realThickness and wavVecs
title('Cu(001) plate: phase velocity');
Vis.plotDispersionCurvesPlate(sweepDisp, curveTypeFreq, angleIdx, colorMode, confidenceThresh);
ylim(freqLimsDisp);
xlabel('k [rad/mm]');
ylabel('frequency [MHz]');
title('Cu(001) plate: frequency');
Vis.plotDispersionCurvesPlate(sweepDisp, curveTypeGroup, angleIdx, colorMode, confidenceThresh);
xlim(freqLimsDisp);
xlabel('frequency [MHz]');
yline(0, 'k:');                    % zero group velocity at finite k = ZGV
title('Cu(001) plate: group velocity');

%% Part 2b -- branch labels and zero-group-velocity points
classResult = sweepDisp.classifyBranches(angleIdx, confidenceThresh);
fprintf('Mid-plane / sagittal mirror: %d / %d; first-k labels: %s\n', ...
    classResult.isExactMirrorSymmetric, classResult.isExactSagittalSymmetric, strjoin(classResult.displayLabel(1, 1:6), ', '));
maxCount = Inf;                    % cap on returned ZGV points
zgvList  = sweepDisp.locateZGVs(angleIdx, maxCount);
zgvShow  = find([zgvList.freqZGV] <= freqLimsDisp(2));   % the plotted frequency range
for iZgv = zgvShow                 % fixed-angle search: isTrueZGV only on a symmetry axis (Part 5)
    fprintf('ZGV branch %2d: k = %.4f rad/mm, f = %.4f MHz, true ZGV: %d\n', ...
        zgvList(iZgv).branchNo, zgvList(iZgv).wavVecZGV, zgvList(iZgv).freqZGV, zgvList(iZgv).isTrueZGV);
end

windowMode = 'combined';           % 'combined' | 'split'
guiDisp = Gui.dispersionViewerPlate(sweepDisp, windowMode);
guiDisp.show();
guiDisp.curveTypeDropdown.Value = curveTypeVphase;   % opens on 'wavelength-freq'
guiDisp.onCurveTypeChanged();      % autoscales, so the range comes after it
guiDisp.setAxisRange('y', velLimsDisp);

%% Part 3 -- many angles, one wavevector
anglesDeg = (0:2:90)';
sweepAng  = Calc.angSweepPlate(structure, anglesDeg, realThickness, realWavVec, retentionLevel, useParallel);
sweepAng  = sweepAng.calc();
velLimsAng = [1 5.5];              % km/s, both maps and the viewer
figure();
tiledlayout(1, 2);                 % type and confidence side by side
axType = nexttile;
Vis.plotRitzModeMapPlate(sweepAng, criterion, [], axType);   % .velsMap etc. are dense [nModes x nAngles]
ylim(axType, velLimsAng);
title(axType, sprintf('Cu(001) plate, k = %g rad/mm: by type', realWavVec));
axConf = nexttile;
Vis.plotRitzModeMapPlate(sweepAng, 'confidence', [], axConf);
ylim(axConf, velLimsAng);
title(axConf, 'mode-type confidence');
guiAng = Gui.combinedViewerPlate(sweepAng);
guiAng.show();
guiAng.setAxisRange('y', velLimsAng);

%% Part 4 -- coated plate (10% Al on one face): not mirror-symmetric, labels become heuristic
matAlCrystal   = Material.FromLibrary('Al');   % isotropic Al = Hill average of the single crystal
[alC11, alC44] = matAlCrystal.elConsts.GetHillAverages();
alConsts  = [alC11, alC44];        % [C11, C44] GPa
isoSymmNo = 1;
alDensity = matAlCrystal.rho;      % g/cm^3
ecAl      = ElastConsts(alConsts, isoSymmNo);
matAl     = Material(ecAl, alDensity);
coatingFraction = 0.1;             % relDepth ratio = thickness ratio
cpCu = CalcParams((1 - coatingFraction) * relDepthPlaceholder, legOrder);
cpAl = CalcParams(coatingFraction * relDepthPlaceholder, legOrder);
structureCoated = Structure(material, cpCu, 'cu_plate_al_coating', nMaxModes);
structureCoated = structureCoated.AddLayer(matAl, cpAl);
fprintf('Coated plate mirror-symmetric: %d\n', structureCoated.isMirrorSymmetric());

sweepCoated = Calc.wavevecAngSweep2Plate(structureCoated, angleDeg, wavVecs, realThickness, retentionSweep, useParallel);
sweepCoated = sweepCoated.calc();
classCoated = sweepCoated.classifyBranches(angleIdx, confidenceThresh);
fprintf('Trusted labels: %.0f%%\n', 100 * mean(classCoated.trustLabelByBranch(:)));
colorModeBranch = 'branch';        % escape hatch when labels are not trusted
Vis.plotDispersionCurvesPlate(sweepCoated, curveTypeVphase, angleIdx, colorMode, confidenceThresh);
ylim(velLimsDisp);
xlabel('k [rad/mm]');
title('Cu + 10% Al: phase velocity, by type');
Vis.plotDispersionCurvesPlate(sweepCoated, curveTypeVphase, angleIdx, colorModeBranch, confidenceThresh);
ylim(velLimsDisp);
xlabel('k [rad/mm]');
title('Cu + 10% Al: phase velocity, by branch');

%% Part 5 -- Si(001) wafer: S1 ZGV vs published measurements (Kiefer, Mezil & Prada 2023, Sci. Adv. 9, eadk6846)
siConsts  = [165.6, 63.9, 79.5];   % [C11, C12, C44] GPa, Kiefer et al. (library Si differs ~0.1 %)
siDensity = 2.330;                 % g/cm^3
siDirMat  = [0 0 1; 1 0 0];        % normal [001]; angle 0 = [100], 45 = [110]
cubicSymmNo = 2;
materialSi  = Material(ElastConsts(siConsts, cubicSymmNo), siDensity, siDirMat);
structureSi = Structure(materialSi, calcParams, 'si_001_wafer', nMaxModes);
siThickness = 0.5246;              % mm
siAnglesDeg = [0; 45; 18.4];      % [100], [110], 26.6 deg from [110] (a TGV wave)
siWavVecs   = linspace(1, 6, 51)'; % rad/mm
sweepSi = Calc.wavevecAngSweep2Plate(structureSi, siAnglesDeg, siWavVecs, siThickness, retentionSweep, useParallel);
sweepSi = sweepSi.calc();
publishedComputed = [7.7079, 7.7551, 7.727];   % MHz (GEWtool)
publishedMeasured = [7.7072, 7.7545, NaN];     % MHz
for iAng = 1:numel(siAnglesDeg)
    zgvSi = sweepSi.locateZGVs(iAng, Inf);
    [~, iLow] = min([zgvSi.freqZGV]);   % S1 point (A2 is near 15 MHz)
    z = zgvSi(iLow);
    fprintf('%4.1f deg: k = %.4f rad/mm, f = %.4f MHz (published %.4f, measured %.4f), true ZGV %d, c_gtheta = %.4f mm/us\n', ...
        siAnglesDeg(iAng), z.wavVecZGV, z.freqZGV, publishedComputed(iAng), publishedMeasured(iAng), ...
        z.isTrueZGV, z.velGroupTrans);   % 18.4 deg: TGV, c_gtheta 0.178 (published -0.177 at the mirror direction, 71.6 deg)
end
