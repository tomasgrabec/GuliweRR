%% Tutorial S1 (code) -- building a Material from literature data: Voigt input, conversions, orientation, averaging
% Narrated version: examples/tutorials/tutorial_S1_material_from_literature.m. Run guliwerr/init.m first.

%% 1. Full 6x6 Voigt matrix input, symmetry given explicitly (SiC, cubic)
voigtSiC = [390 134 134   0   0   0; ...
            134 390 134   0   0   0; ...
            134 134 390   0   0   0; ...
              0   0   0 253   0   0; ...
              0   0   0   0 253   0; ...
              0   0   0   0   0 253];   % GPa
symmCubic = 2;                          % ElastConsts.SymmClassLibrary() lists all classes
ecSiC = ElastConsts(voigtSiC, symmCubic);
fprintf('[C11, C12, C44] = [%.1f, %.1f, %.1f] GPa\n', ecSiC.vector);

%% 2. Symmetry auto-detection -- exists but fragile, prefer an explicit symmNo
symmDetected = ElastConsts.FindSymmClass_input(voigtSiC);
fprintf('Auto-detected class %d (%s)\n', symmDetected, ElastConsts.SymmClassLibrary(symmDetected));

%% 3. Cubic (C11, C12, C44) <-> (B, c', c44) round trip
[bulkShearC44, ~] = ElastConsts.CubicConversion(ecSiC.vector);   % [B, (C11-C12)/2, C44]
standardVector    = ElastConsts.CubicDeconversion(bulkShearC44);
fprintf('(B, c'', c44) = [%.2f, %.2f, %.2f], back: [%.1f, %.1f, %.1f] GPa\n', bulkShearC44, standardVector);

%% 4. Orientation from direction rows: row 1 = surface normal, row 2 = propagation reference
dirMatSiC  = [0 0 1; 1 1 0];            % (001) surface, [110] reference
densitySiC = 3.21;                      % g/cm^3
nameSiC    = 'SiC (001)/<110>';
matSiC     = Material(ecSiC, densitySiC, dirMatSiC, nameSiC);

%% 5. Zener ratio -- selects the CalcParams tier (docs/pitfalls.md MTH1)
anisotropyFactor = matSiC.elConsts.GetAnisotropyRatio();
fprintf('Zener ratio 2*C44/(C11-C12) %.2f\n', anisotropyFactor);

%% 6. Orientation from Bunge-ZXZ Euler angles (EBSD form)
phi1Deg = 30; PhiDeg = 45; phi2Deg = 10;
dirMatEuler = Material.Euler2DirMat(phi1Deg, PhiDeg, phi2Deg);
matSiCEuler = Material(ecSiC, densitySiC, dirMatEuler, 'SiC (Euler 30/45/10)');
[phi1Back, PhiBack, phi2Back] = Material.DirMat2Euler(dirMatEuler);
fprintf('Euler round trip (%.1f, %.1f, %.1f) deg\n', phi1Back, PhiBack, phi2Back);

%% 7. Orientation from Miller indices -- normal = plane normal, reference = direction
latticeParams = struct();               % cubic needs none; Material.Miller2DirMat(symmNo) lists fields
normalHkl     = [1 1 1];
propUvw       = [1 -1 0];
normalDirCos  = Material.Miller2DirMat(normalHkl, symmCubic, latticeParams, 'planeNormal');
propDirCos    = Material.Miller2DirMat(propUvw, symmCubic, latticeParams, 'direction');
dirMatMiller  = [normalDirCos; propDirCos];
matSiC111     = Material(ecSiC, densitySiC, dirMatMiller, 'SiC (111)/[1-10]');

%% 8. Rotate the in-plane reference about the surface normal
inPlaneRotationDeg = 30;
dirMatFull    = Material.CompleteDirMat(dirMatSiC);   % RotDirMatFromAngle needs the 3x3 frame
dirMatRotated = Material.RotDirMatFromAngle(dirMatFull, inPlaneRotationDeg);
matSiCRotated = Material(ecSiC, densitySiC, dirMatRotated, 'SiC (001), 30 deg');

%% 9. Polycrystalline film -- Voigt / Reuss / Hill isotropic averages [C11, C44]
matAlCrystal = Material.FromLibrary('Al');   % Al single crystal from the materials library
ecAlCubic    = matAlCrystal.elConsts;        % [C11, C12, C44] GPa
[alC11Voigt, alC44Voigt] = ecAlCubic.GetVoigtAverages();   % upper bound
[alC11Reuss, alC44Reuss] = ecAlCubic.GetReussAverages();   % lower bound
[alC11Hill,  alC44Hill]  = ecAlCubic.GetHillAverages();    % mean, usual choice
fprintf('Al [C11, C44]: Voigt [%.2f %.2f], Reuss [%.2f %.2f], Hill [%.2f %.2f] GPa\n', ...
    alC11Voigt, alC44Voigt, alC11Reuss, alC44Reuss, alC11Hill, alC44Hill);
symmIsotropic = 1;
ecAlHill   = ElastConsts([alC11Hill, alC44Hill], symmIsotropic);
densityAl  = matAlCrystal.rho;            % g/cm^3
dirMatNone = [];                        % isotropic: orientation irrelevant
matAlFilm  = Material(ecAlHill, densityAl, dirMatNone, 'Al film (Hill)');

%% 10. Two-layer Structure opened in the structure-building GUI
filmThicknessNm = 200;
relDepthSub     = 20;                   % "easy" tier as a starting point
legOrder        = 50;
relDepthPlaceholder = 1;                % film relDepth, set per wavelength
nMaxModes       = 500;
cpSubstrate = CalcParams(relDepthSub, legOrder);
cpFilm      = CalcParams(relDepthPlaceholder, legOrder, filmThicknessNm);
structure   = Structure(matSiC, cpSubstrate, 'SiC_sub_Al_film', nMaxModes);
structure   = structure.AddLayer(matAlFilm, cpFilm);
gui = Gui.structureControl(structure);
gui.show();

%% 11. Standard crystals from the materials library (default row, or by 'T' / 'variant')
Material.ListLibrary();
matCu001 = Material.FromLibrary('Cu', [0 0 1; 1 1 0]);   % Cu(001), reference [110]
matMgO   = Material.FromLibrary('MgO', [0 0 1; 1 0 0]);   % MgO(001), reference [100]
[matTi, entryTi] = Material.FromLibrary('titanium');       % name or alias, any case
fprintf('%s: rho = %.3f g/cm^3, C11 = %.1f GPa (source: %s)\n', ...
    matTi.name, matTi.rho, matTi.elConsts.voigt(1,1), entryTi.source);
