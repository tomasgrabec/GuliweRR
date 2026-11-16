%% Example: epitaxial NiTi film on MgO -- austenite vs. martensite Green's maps
%
% Walkthrough of the model system in Grabec et al. 2024 (Ultrasonics 138,
% 107211), "Guided acoustic waves in thin epitaxial films: Experiment and
% inverse problem solution for NiTi" -- a 3 um NiTi shape-memory film,
% epitaxially grown on a single-crystalline MgO(100) substrate,
% TGS-characterised at 21 degC (film in the twinned martensite phase) and
% 111 degC (film in single-crystalline austenite).
%
% Two thin Cr/V buffer layers (nominal 5 nm / 50 nm) sit between film and
% substrate. IGNORED HERE, matching the paper's own inverse procedure: the
% MgO elastic constants used below are the paper's "effective substrate"
% values -- fitted directly from the bare-substrate (film-free) TGS scan --
% which the paper explicitly checked already approximately absorbs the
% buffer layers' influence ("found the difference in velocity of using the
% effective elasticity versus the MgO plus buffer layers ... to be minimal
% ... we continue with the effective elasticity only", Section 4.1).
%
% Follows the paper's own three-step structure:
%   1. MgO substrate alone -- the effective elasticity the paper's inverse
%      procedure found from the bare-substrate scan (Section 4.1, Fig. 4).
%   2. NiTi austenite film (111 degC) on MgO (111 degC) -- Section 4.2, Fig. 5.
%   3. NiTi martensite film (21 degC) on MgO (21 degC) -- Section 4.3, Fig. 7.
% Each phase is measured (and so recalculated here) at its own temperature,
% since the paper independently re-fitted the MgO substrate constants at
% both temperatures (Section 4.1's second paragraph) -- reused as-is, not
% re-derived.
%
% Density note: the paper does not state a NiTi film density. rho = 6.45 g/cm^3
% below is an example input value, not taken from the paper. The paper treats any density change
% between austenite and martensite as insignificant ("the change in density
% was considered insignificant", Section 4.1), so one value is used for both.
%
% Run run('guliwerr/init.m') from the project root first (see
% examples/tutorials/tutorial_01_isotropic_halfspace.m for why not from inside it).

%% Materials -- MgO substrate, both temperatures (Section 4.1)
% Cubic (symmClassNo=2), "effective" substrate elasticity -- already
% inverted FROM the angular TGS data, not literature MgO single-crystal
% constants (the film-inversion step then treats this as the known
% substrate). dirMat is the trivial identity mapping (surface normal =
% local [001], in-plane reference = local [100]) -- cubic symmetry makes
% the actual crystallographic axis assignment irrelevant to the resulting
% Voigt tensor.
symmCubic = 2;

elC_MgO_21C  = [298.3, 105.4, 149.7];   % GPa, at 21 degC (room temp, Fig. 4)
elC_MgO_111C = [290.2, 104.0, 149.1];   % GPa, at 111 degC
rho_MgO      = 3.58;                    % g/cm^3, Archimedes method (both temperatures)
dM_MgO       = [0 0 1; 1 0 0];

mat_MgO_21C  = Material(ElastConsts(elC_MgO_21C,  symmCubic), rho_MgO, dM_MgO, 'MgO (21 degC)');
mat_MgO_111C = Material(ElastConsts(elC_MgO_111C, symmCubic), rho_MgO, dM_MgO, 'MgO (111 degC)');

%% Materials -- NiTi film, both phases (Sections 4.2/4.3, Eqs. (15) and (17))
% Tetragonal type I (symmClassNo=3): [c11,c12,c13,c33,c44,c66], c-axis
% (index 3) perpendicular to the film surface -- the direction singled out
% by the epitaxial strain from the substrate, regardless of which cubic
% <100> axis it originally was.
%
% In-plane orientation is NOT arbitrary: the angle=0 reference direction is
% a LAB direction shared with the substrate, and the two materials' <100>-
% type in-plane axes are NOT aligned with each other. The epitaxial growth
% relation MgO(100)[001] || NiTi B2(100)[011] (Section 3) means the
% austenite/martensite lattice sits 45 deg rotated in-plane relative to MgO
% -- "[100] of the MgO corresponds to [110] of the austenite film" (Fig.
% 3c). dM_MgO supplies the lab angle=0 reference (its own
% <100>-type in-plane axis); dM_NiTi's in-plane row is [1 1 0], not [1 0 0],
% so the film's <110>-type direction lines up with lab angle=0, matching
% the growth relation. Applied identically to both phases -- the paper's
% martensite analysis assumes the effective (twin-averaged) tetragonal axes
% stay aligned the same way the austenite's were (Section 4.3).
symmTetrI = 3;

elC_NiTi_A = [184.7, 151.5, 150.9, 183.3, 33.98, 33.70];   % GPa, austenite @ 111 degC, Eq. (15)
elC_NiTi_M = [192.4, 116.5, 117.1, 201.1, 18.40, 20.40];   % GPa, martensite @ 21 degC,  Eq. (17)
rho_NiTi   = 6.45;                                          % g/cm^3, see density note above
dM_NiTi    = [0 0 1; 1 1 0];                                % 45 deg rotated in-plane vs. MgO -- see above
h_NiTi_nm  = 3060;                                          % SEM-measured thickness, Fig. 3d

mat_NiTi_A = Material(ElastConsts(elC_NiTi_A, symmTetrI), rho_NiTi, dM_NiTi, 'NiTi austenite (111 degC)');
mat_NiTi_M = Material(ElastConsts(elC_NiTi_M, symmTetrI), rho_NiTi, dM_NiTi, 'NiTi martensite (21 degC)');

%% Mesh resolution -- "easy" tier (docs/pitfalls.md MTH1): both MgO(100)
% and the film's c-axis are principal-plane surface normals, and neither
% material's Zener ratio 2*c44/(c11-c12) exceeds ~2.2 (MgO: 1.55;
% NiTi austenite/martensite: reported Zener ratio A ~ 2/2.2, Sections 4.1/
% 4.3), so relDepth = 20 and legOrder = 50.
legO         = 50;
relDepthSubs = 20;
cp_substrate = CalcParams(relDepthSubs, legO);

% Film relDepth is derived from real dimensions (Structure.ApplyWavelength
% below), not set directly -- matches the paper's own TGS excitation
% wavelength lambda = 10.04 +/- 0.01 um (Section 2.1), so the calculation
% below reproduces the same (thickness, wavelength) operating point as
% Figs. 5b/7b, not an arbitrary normalisation choice.
lambdaRefNm = 10040;   % 10.04 um, paper's TGS fringe spacing
cp_film     = CalcParams(1, legO, h_NiTi_nm);   % relDepth placeholder=1, overwritten below

%% Structures -- bottom (substrate) to top (film)
maxModes = 500;

strAustenite = Structure(mat_MgO_111C, cp_substrate, 'NiTi_austenite_on_MgO', maxModes);
strAustenite = strAustenite.AddLayer(mat_NiTi_A, cp_film);
strAustenite = strAustenite.ApplyWavelength(lambdaRefNm);

strMartensite = Structure(mat_MgO_21C, cp_substrate, 'NiTi_martensite_on_MgO', maxModes);
strMartensite = strMartensite.AddLayer(mat_NiTi_M, cp_film);
strMartensite = strMartensite.ApplyWavelength(lambdaRefNm);

%% Angular sweep grid -- shared by all three calculations
% 0-180 deg suffices (velocities are even in propagation direction); a
% 2 deg step still resolves the Rayleigh/Sezawa/Love angular variation of
% the paper's Figs. 4b/5b/7b.
angles = (0:2:180)';

% The layered Green axis ends by default at 1.05 x the fastest longitudinal
% limiting velocity of any layer (here MgO's, ~9-10 km/s), well above the
% substrate shear limits (~6.2-6.5 km/s) where the "above/below substrate
% ST" transition (Section 4.2) sits, so no vMax is passed.

%% 1. MgO substrate alone -- "the resulting elasticity" (Section 4.1, Fig. 4)
strMgO  = Structure(mat_MgO_21C, cp_substrate, 'MgO_substrate_21C', maxModes);
calcMgO = Calc.angSweepHS(strMgO, angles);
calcMgO = calcMgO.calc();

%% 2. NiTi austenite film on MgO (Section 4.2, Fig. 5)
calcAustenite = Calc.angSweepLayered(strAustenite, angles);
calcAustenite = calcAustenite.calc();

%% 3. NiTi martensite film on MgO (Section 4.3, Fig. 7)
calcMartensite = Calc.angSweepLayered(strMartensite, angles);
calcMartensite = calcMartensite.calc();

%% Plot Green's maps -- compare directly to the paper's Figs. 4b, 5b, 7b
% Both the bare-MgO map (1.) and the two LAYERED maps (2./3.) reproduce the
% paper's Rayleigh/Sezawa multi-branch structure (Figs. 4b/5b/7b): Ritz and
% GreenSolver agree to <0.01% on every guided-mode branch.
[~, axMgO]  = Vis.plotGreenMap(calcMgO);
title(axMgO, 'MgO substrate, 21 degC');

[~, axAust] = Vis.plotGreenMap(calcAustenite);
title(axAust, 'NiTi austenite (111 degC) on MgO');

[~, axMart] = Vis.plotGreenMap(calcMartensite);
title(axMart, 'NiTi martensite (21 degC) on MgO');

%% ===================================================================
%% 4. NiTi austenite ONLY, thickness comparison
%% ===================================================================
%
% Same austenite film/substrate materials as step 2 (mat_NiTi_A on
% mat_MgO_111C), now for THREE film thicknesses at angle 0 over several real
% wavelengths -- the fixed-angle wavelength sweep of
% examples/tutorials/tutorial_06_film_on_substrate.m, once per thickness.
% The thickness is part of each film's CalcParams (realThickness), not a
% sweep axis, so each thickness needs its own Structure and its own sweep.
%
% h_NiTi_nm (3060 nm, step 2's SEM-measured value, Fig. 3d) is the "3 micron"
% case, alongside two thinner ones.
niTiThicknessesNm = [1000, 1600, h_NiTi_nm];   % nm -- 1, 1.6, 3.06 um
wavelengthsUm     = 3:2:11;                     % um
wavelengthsNm     = wavelengthsUm * 1000;       % nm

vSAWByThickness = nan(numel(niTiThicknessesNm), numel(wavelengthsNm));
for it = 1:numel(niTiThicknessesNm)
    cp_film_it = CalcParams(1, legO, niTiThicknessesNm(it));
    strAusteniteIt = Structure(mat_MgO_111C, cp_substrate, ...
        sprintf('NiTi_austenite_%dnm_on_MgO', niTiThicknessesNm(it)), maxModes);
    strAusteniteIt = strAusteniteIt.AddLayer(mat_NiTi_A, cp_film_it);

    sweepIt = Calc.wavevecAngSweep2Layered(strAusteniteIt, 0, wavelengthsNm, [], 'ritzOnly', 'full');
    sweepIt = sweepIt.calc();

    % branch seeded at the slowest SAW of the shortest wavelength
    vSAWByThickness(it, :) = DispersionBranches.trackSeededBranch(sweepIt);
end

figure();
hold on;
thicknessColors = lines(numel(niTiThicknessesNm));
for it = 1:numel(niTiThicknessesNm)
    plot(wavelengthsUm, vSAWByThickness(it, :), 'o-', 'Color', thicknessColors(it, :), ...
        'DisplayName', sprintf('%.2f um NiTi', niTiThicknessesNm(it) / 1000));
end
xlabel('wavelength [\mum]');
ylabel('v_{SAW} [km/s]');
legend('show', 'Location', 'best');
title('NiTi austenite on MgO -- tracked branch vs. wavelength, by film thickness');

%% ===================================================================
%% 5. Literature MgO + Cr/V buffer layers (first case only: 3 um NiTi, 10 um wavelength)
%% ===================================================================
%
% Checks the paper's own "buffers are negligible" claim (Section 4.1,
% quoted in this file's header) with GuliweRR directly, rather than only
% repeating it -- adds the Cr/V buffer layers steps 2/3 (and the paper
% itself) deliberately left out, plus swaps the paper's fitted "effective"
% MgO substrate constants for literature single-crystal MgO values (same
% physical substrate cut as dM_MgO above, only the elastic-constant SOURCE
% differs).
%
% Epitaxial relation for Cr and V (Grabec et al. 2024, main text): MgO(100)
% [001] || Cr,V,NiTi (100)[011] -- i.e. Cr and V share NiTi's OWN 45 deg
% in-plane rotation relative to MgO (dM_NiTi above), NOT MgO's own [100]
% in-plane reference (dM_MgO) -- confirmed directly against the paper, not
% assumed by analogy with the NiTi relation alone.
elC_MgO_lit = [295.9, 95.4, 153.9];   % GPa, literature single crystal
mat_MgO_lit = Material(ElastConsts(elC_MgO_lit, symmCubic), rho_MgO, dM_MgO, 'MgO (literature)');

elC_Cr  = [350, 67.8, 100.8];   % GPa, cubic
rho_Cr  = 7.1988;                % g/cm^3
h_Cr_nm = 5;                      % nm, nominal (see the header)
mat_Cr  = Material(ElastConsts(elC_Cr, symmCubic), rho_Cr, dM_NiTi, 'Cr buffer');

elC_V  = [228.7, 119, 43.15];   % GPa, cubic
rho_V  = 6.1;                    % g/cm^3
h_V_nm = 50;                      % nm, nominal
mat_V  = Material(ElastConsts(elC_V, symmCubic), rho_V, dM_NiTi, 'V buffer');

cp_Cr = CalcParams(1, legO, h_Cr_nm);   % relDepth placeholder -- set by ApplyWavelength below
cp_V  = CalcParams(1, legO, h_V_nm);

strAusteniteBuff = Structure(mat_MgO_lit, cp_substrate, 'NiTi_austenite_buffered_on_MgO_lit', maxModes);
strAusteniteBuff = strAusteniteBuff.AddLayer(mat_Cr, cp_Cr);
strAusteniteBuff = strAusteniteBuff.AddLayer(mat_V,  cp_V);
strAusteniteBuff = strAusteniteBuff.AddLayer(mat_NiTi_A, cp_film);   % cp_film = step 2's own, h_NiTi_nm
strAusteniteBuff = strAusteniteBuff.ApplyWavelength(lambdaRefNm);

% The Green axis end follows the fastest longitudinal limit of ALL FOUR layers
% by itself -- adding Cr/V could shift which layer is fastest.
calcAusteniteBuff = Calc.angSweepLayered(strAusteniteBuff, angles);
calcAusteniteBuff = calcAusteniteBuff.calc();

[~, axAustBuff] = Vis.plotGreenMap(calcAusteniteBuff);
title(axAustBuff, 'NiTi austenite on literature MgO + Cr/V buffers');
% Compare with axAust above (step 2: no buffers, effective-fit MgO); both
% use the same colour scale and axes.
