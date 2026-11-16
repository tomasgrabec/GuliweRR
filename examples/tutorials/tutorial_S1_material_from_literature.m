%% Tutorial: building a Material from real literature/lab data
%
% The other tutorials start from GuliweRR's compact constant vectors
% ([C11, C12, C44] for cubic, etc.) at the default orientation. Data from a
% paper rarely look like that: a full Voigt matrix, or another
% parameterisation (e.g. bulk modulus and shear constants), and a sample cut
% that is not the default. This script builds a Material from such data:
% full-Voigt input, symmetry detection, converting parameterisations, and
% the crystal orientation from direction rows, Euler angles, Miller indices
% or an in-plane rotation; then isotropic averaging for a polycrystalline
% film and the structure GUI. The SiC constants below are example input
% values; they make the round trips checkable.
%
% Copy the parts you need to build a Material from your own data.
%
% Run run('guliwerr/init.m') from the project root first (tutorial_01 says why).

%% 1. Full Voigt matrix input, symmetry given explicitly (the safe path).
% ElastConsts takes the 6x6 matrix directly.
% SiC (3C, cubic): C11 = 390, C12 = 134, C44 = 253 GPa.
voigtSiC = [390  134  134    0    0    0; ...
            134  390  134    0    0    0; ...
            134  134  390    0    0    0; ...
              0    0    0  253    0    0; ...
              0    0    0    0  253    0; ...
              0    0    0    0    0  253];   % GPa, full symmetric 6x6

symmCubic = 2;   % see ElastConsts.SymmClassLibrary() for every class's number/constant order
ecSiC = ElastConsts(voigtSiC, symmCubic);

fprintf('Parsed from full Voigt matrix: [C11,C12,C44] = [%.1f, %.1f, %.1f] GPa\n', ecSiC.vector);

%% 2. Symmetry auto-detection (ElastConsts.FindSymmClass_input, used when
% symmNo is omitted) counts the distinct values of the Voigt matrix, so it
% misidentifies the class when two independent constants happen to be equal.
% Supply symmNo explicitly, as in section 1.
symmDetected = ElastConsts.FindSymmClass_input(voigtSiC);
fprintf('Auto-detected symmetry class: %d (%s)\n', symmDetected, ElastConsts.SymmClassLibrary(symmDetected));

%% 3. Converting between parameterisations. Some papers give cubic
% elasticity as (B, c', c44) -- bulk modulus, tetragonal shear constant
% c' = (C11-C12)/2, c44. ElastConsts.CubicConversion goes from (C11, C12, C44)
% to (B, c', c44), CubicDeconversion back; a round trip checks the numbers.
[BcPrimeC44, ~] = ElastConsts.CubicConversion(ecSiC.vector);
fprintf('SiC as (B, c'', c44)         = [%.2f, %.2f, %.2f] GPa\n', BcPrimeC44);

%% A paper reporting (B, c', c44) -- recover the (C11, C12, C44) ElastConsts needs:
recovered = ElastConsts.CubicDeconversion(BcPrimeC44);
fprintf('Recovered [C11,C12,C44]     = [%.1f, %.1f, %.1f] GPa (matches section 1)\n', recovered);
% Tetragonal and hexagonal have the same pair: TetragonalConversion /
% TetragonalDeconversion and HexagonalConversion / HexagonalDeconversion.

%% 4. Crystal orientation -- give Material a dirMat in the crystal's own
% axes; no manual tensor rotation (ElastConsts.RotateTensor) is needed.
%   row 1 = surface normal
%   row 2 = in-plane propagation reference (lab angle = 0)
%   row 3 = optional; a 2x3 dirMat is completed to a right-handed frame
%
% SiC (001) surface, <110>-type in-plane reference direction:
dirMatSiC = [0 0 1; ...   % surface normal  = [001]
             1 1 0];      % propagation ref = [110]

rhoSiC     = 3.21;   % g/cm^3
matNameSiC = 'SiC (001)/<110>';
matSiC = Material(ecSiC, rhoSiC, dirMatSiC, matNameSiC);

%% 5. Zener ratio -- selects the CalcParams (relDepth/legOrder) tier
% (docs/pitfalls.md MTH1). Cubic: the exact Zener ratio; other symmetries:
% (v_shear_max/v_shear_min)^2 from an angular scan.
A = matSiC.elConsts.GetAnisotropyRatio();
fprintf('SiC Zener ratio 2*C44/(C11-C12) = %.2f\n', A);
fprintf(['-> ratio < 5 on (001): start from the "easy" tier\n' ...
         '   (relDepth=20, legO=50); check convergence (docs/pitfalls.md MTH1).\n']);

%% 6. Orientation from Bunge-ZXZ Euler angles (phi1, Phi, phi2), the form
% EBSD/texture data come in. Material.Euler2DirMat returns a ready-to-use
% dirMat (row 1 = surface normal). Convention as in Dryburgh et al.
% 2020/2022; not yet checked against a published EBSD example.
phi1Deg = 30;   % deg
PhiDeg  = 45;   % deg
phi2Deg = 10;   % deg
dirMatEuler = Material.Euler2DirMat(phi1Deg, PhiDeg, phi2Deg);
matSiCEuler = Material(ecSiC, rhoSiC, dirMatEuler, 'SiC (Euler 30/45/10)');
[phi1Back, PhiBack, phi2Back] = Material.DirMat2Euler(dirMatEuler);   % round trip
fprintf('Euler round trip: (%.1f, %.1f, %.1f) deg\n', phi1Back, PhiBack, phi2Back);

%% 7. Orientation from Miller indices -- Material.Miller2DirMat turns an index
% into a direction-cosine row. The surface normal is a PLANE normal (hkl), the
% propagation reference a real-space DIRECTION [uvw]; the two differ for
% every non-cubic lattice. Material.Miller2DirMat(symmNo) alone lists the
% latticeParams fields that symmetry needs (cubic needs none).
latticeParamsCubic = struct();   % cubic: no lattice-shape parameters
normalDirCos = Material.Miller2DirMat([1 1 1], symmCubic, latticeParamsCubic, 'planeNormal');
propDirCos   = Material.Miller2DirMat([1 -1 0], symmCubic, latticeParamsCubic, 'direction');
dirMatMiller = [normalDirCos; propDirCos];   % (111) surface, [1-10] propagation reference
matSiC111    = Material(ecSiC, rhoSiC, dirMatMiller, 'SiC (111)/[1-10]');

%% 8. Same cut, in-plane reference rotated about the surface normal --
% Material.RotDirMatFromAngle needs the completed 3x3 frame, which
% Material.CompleteDirMat builds from the 2-row input.
inPlaneRotationDeg = 30;   % deg
dirMatFull    = Material.CompleteDirMat(dirMatSiC);
dirMatRotated = Material.RotDirMatFromAngle(dirMatFull, inPlaneRotationDeg);
matSiCRotated = Material(ecSiC, rhoSiC, dirMatRotated, 'SiC (001), 30 deg');
% A sweep's `angles` already rotate per angle; rotate dirMat only to move
% where angle = 0 points.

%% 9. Polycrystalline (e.g. sputtered) film: isotropic constants averaged from
% single-crystal ones, each as the isotropic pair [M, G] = [C11, C44]: Voigt
% (upper bound, uniform strain), Reuss (lower bound, uniform stress), Hill
% (their mean -- the usual practical choice).
matAlCrystal = Material.FromLibrary('Al');   % Al single crystal from the library
ecAlCubic    = matAlCrystal.elConsts;        % [C11, C12, C44] = [107.3, 60.8, 28.3] GPa
[alC11Voigt, alC44Voigt] = ecAlCubic.GetVoigtAverages();
[alC11Reuss, alC44Reuss] = ecAlCubic.GetReussAverages();
[alC11Hill,  alC44Hill]  = ecAlCubic.GetHillAverages();
fprintf('Al isotropic [C11, C44]: Voigt [%.2f, %.2f], Reuss [%.2f, %.2f], Hill [%.2f, %.2f] GPa\n', ...
    alC11Voigt, alC44Voigt, alC11Reuss, alC44Reuss, alC11Hill, alC44Hill);
symmIsotropic = 1;
ecAlHill  = ElastConsts([alC11Hill, alC44Hill], symmIsotropic);
rhoAl     = matAlCrystal.rho;   % g/cm^3, the library's single-crystal density
matAlFilm = Material(ecAlHill, rhoAl, [], 'Al film (Hill-averaged)');   % [] dirMat: isotropic

%% 10. Interactive: the structure GUI on a two-layer stack of the Materials
% above. Gui.structureControl(initialStructure) fills every tab from the
% Structure and can launch a calculation. Needs an interactive MATLAB session.
filmThicknessNm = 200;   % real film thickness [nm]
cpSubstrate = CalcParams(20, 50);                    % "easy" tier, see section 5
cpFilm      = CalcParams(1, 50, filmThicknessNm);    % relDepth placeholder, set per wavelength
strDemo = Structure(matSiC, cpSubstrate, 'SiC_sub_Al_film', 500);
strDemo = strDemo.AddLayer(matAlFilm, cpFilm);
gui = Gui.structureControl(strDemo);
gui.show();

%% 11. Standard crystals from the materials library -- no typing of constants.
% Material.ListLibrary() prints every stored data set (constants, density,
% temperature, source, status); Material.FromLibrary builds the Material,
% with the orientation given exactly as in sections 4-8. Every row names its
% primary measurement; a material with several rows (another temperature
% or source) is chosen by 'T' or 'variant', otherwise its default row.
Material.ListLibrary();
matCu001 = Material.FromLibrary('Cu', [0 0 1; 1 1 0]);   % Cu(001), reference [110]
matMgO   = Material.FromLibrary('MgO', [0 0 1; 1 0 0]);   % MgO(001) substrate, reference [100]
[matTi, entryTi] = Material.FromLibrary('titanium');       % aliases work too
fprintf('%s: rho = %.3f g/cm^3, C11 = %.1f GPa (source: %s)\n', ...
    matTi.name, matTi.rho, matTi.elConsts.voigt(1,1), entryTi.source);
% The same list is the "library material" dropdown in the structure GUI (section 10).

%% From here, matSiC goes on as in tutorial_01/02: CalcParams -> Structure -> Calc.*
