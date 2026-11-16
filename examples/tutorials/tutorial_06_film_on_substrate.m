%% Tutorial: thin film on substrate -- SAW velocity vs. kd (wavevector-thickness product)
%
% A film-on-substrate SAW velocity depends on kd (here the film relDepth =
% realThickness/wavelength): near relDepth = 0 v_SAW sits at the substrate's
% own value; for relDepth >> 1 it approaches the film's own half-space value.
%
% Al substrate / Cu film is a slow film on a fast substrate (Cu is 3.3x
% denser, so its bulk velocities are lower despite its larger moduli). For
% why a per-point pick can jump between branches in a layered sweep, and how
% the branch tracking below avoids it, see docs/pitfalls.md MTH11.
%
% The relDepth sweep runs at one fixed angle through
% Calc.wavevecAngSweep2Layered with a one-element angle vector rather than a
% hand-written loop: the sweep object holds every wavelength step, which
% DispersionBranches.trackSeededBranch needs to follow one physical branch
% by Ritz-eigenvector overlap (modal assurance criterion, MAC).
%
% Copy this script as a starting point. Call sequence: docs/quickStart.md.
%
% Run run('guliwerr/init.m') from the project root first (tutorial_01 says why).

%% Substrate -- isotropic Al, faster than the Cu film in both
% v_L = sqrt(C11/rho) and v_T = sqrt(C44/rho). Isotropic (polycrystalline) Al
% = the Hill average of the single-crystal constants from the materials
% library (tutorial_S1 section 9 shows the averages).
matAlCrystal   = Material.FromLibrary('Al');
[alC11, alC44] = matAlCrystal.elConsts.GetHillAverages();
alConsts = [alC11, alC44];   % [C11, C44] GPa, isotropic class
alSymm   = 1;
alRho    = matAlCrystal.rho; % g/cm^3
relDepthSub = 20;   % substrate depth, in wavelengths

%% Film -- Cu(001), cubic, slower than the Al substrate
% Cu constants and density come from the materials library.
dirMatFilm = [1 0 0; ...     % surface normal
              0 0 1];        % propagation reference (Material's default)
filmThicknessNm = 200;   % real film thickness [nm]; the sweep divides it by
                         % each wavelength to get the film relDepth

%% Discretisation -- legO must stay >= 2x the largest relDepth of any layer
% (substrate 20, film up to 3). Cu's Zener ratio 2*C44/(C11-C12) = 3.2
% on a principal plane: the "easy" tier (docs/pitfalls.md MTH1).
legO      = 50;   % Legendre basis order per layer

ecSub  = ElastConsts(alConsts, alSymm);
matSub = Material(ecSub, alRho);
cpSub  = CalcParams(relDepthSub, legO);

matFilm = Material.FromLibrary('Cu', dirMatFilm);
relDepthFilmPlaceholder = 1;   % replaced at every wavelength of the sweep
cpFilm = CalcParams(relDepthFilmPlaceholder, legO, filmThicknessNm);

strName   = 'al_sub_cu_film';
nMaxModes = 500;   % max number of Ritz modes computed/retained
str = Structure(matSub, cpSub, strName, nMaxModes);   % layer 1 = substrate
str = str.AddLayer(matFilm, cpFilm);                   % layer 2 = film, on top

%% Sweep -- ONE fixed angle, many real wavelengths (relDepth =
% filmThicknessNm/wavelength). Descending wavelength = ascending relDepth, so
% tracking starts from the thinnest film, where the branches are clearest.
% Tracking matches CONSECUTIVE steps by eigenvector overlap, so the steps
% are dense (0.01) up to relDepth = 1 and coarser (0.05) beyond. The sweep
% starts at relDepth = 0.01 (a 200 nm film at a 20 um wavelength), where the
% branch is still close to the Al value but the film is no longer negligible.
filmRelDepths = [0.01:0.01:1, 1.05:0.05:3];
wavelengths   = filmThicknessNm ./ filmRelDepths;   % nm, descending

angle          = 0;      % deg -- fixed propagation angle for this sweep
vMax           = [];     % not used by calcMode='ritzOnly'
calcMode       = 'ritzOnly';
retentionLevel = 'full'; % REQUIRED: branch tracking compares every eigenvector
                         % of consecutive wavelength steps

sweep = Calc.wavevecAngSweep2Layered(str, angle, wavelengths, vMax, calcMode, retentionLevel);
sweep = sweep.calc();

%% Branch tracking -- one physical branch, followed by Ritz-eigenvector
% overlap (MAC). The seed is the slowest classified SAW at the thinnest film:
% the Rayleigh-type wave. (DispersionBranches.track gives every branch.)
vSAW = DispersionBranches.trackSeededBranch(sweep).';

fprintf('film relDepth :'); fprintf(' %6.2f', filmRelDepths); fprintf('\n');
fprintf('v_SAW [km/s]  :'); fprintf(' %6.3f', vSAW);          fprintf('\n');

%% Reference lines -- each material's own half-space SAW velocity at this
% angle, the two limits the tracked branch moves between
strAl = Structure(matSub, CalcParams(relDepthSub, legO), 'al_hs', nMaxModes);
csAl  = Calc.singleHS(strAl, angle);
csAl  = csAl.calcAll('ritzOnly');
smAl    = csAl.classification.significantModes;
vSAW_Al = min([smAl(strcmp({smAl.type}, 'Rayleigh')).vel]);   % the Rayleigh wave

strCu = Structure(matFilm, CalcParams(relDepthSub, legO), 'cu_hs', nMaxModes);
csCu  = Calc.singleHS(strCu, angle);
csCu  = csCu.calcAll('ritzOnly');
smCu    = csCu.classification.significantModes;
vSAW_Cu = min([smCu(strcmp({smCu.type}, 'Rayleigh')).vel]);

fprintf('Al half-space v_SAW: %.4f km/s\n', vSAW_Al);
fprintf('Cu half-space v_SAW: %.4f km/s\n', vSAW_Cu);

%% Plot 1 -- the tracked branch vs. relDepth, with the two half-space limits
figure();
semilogx(filmRelDepths, vSAW, 'o-');
hold on;
yline(vSAW_Al, '--', sprintf('Al (substrate), %.3f km/s', vSAW_Al), ...
    'LabelHorizontalAlignment', 'left', 'Color', [0.4660 0.6740 0.1880]);
yline(vSAW_Cu, '--', sprintf('Cu (film), %.3f km/s', vSAW_Cu), ...
    'LabelHorizontalAlignment', 'left', 'Color', [0.8500 0.3250 0.0980]);
xlabel('film relDepth (\propto kd)');
ylabel('v_{SAW} [km/s]');
title('Al substrate / Cu film -- tracked branch vs. film relDepth');

%% Plot 2 -- every tracked Ritz branch (thin grey), SAW/limiting branches
% highlighted: Plot 1's branch in the context of the whole Ritz spectrum
curveType = 'k-vphase';   % 'wavelength-freq' | 'wavelength-vphase' | 'k-freq' | 'k-vphase' | 'freq-vphase' | 'freq-vgroup'
Vis.plotDispersionCurves(sweep, curveType);
title('Al substrate / Cu film -- all tracked branches (Vis.plotDispersionCurves)');
xlim([0 0.03])   % k [rad/nm]
ylim([1 4])      % km/s

%% Interactive: click a curve (or step Wavelength/Branch) to see that mode's
% displacement field. Needs an interactive MATLAB session.
gui = Gui.dispersionViewer(sweep);
gui.show();

% Same curve type and axis range as Plot 2 (the viewer opens on
% 'wavelength-freq', autoscaled)
gui.curveTypeDropdown.Value = curveType;
gui.onCurveTypeChanged();
gui.setAxisRange('x', [0 0.03]);   % also fills the viewer's range fields
gui.setAxisRange('y', [1 4]);

%% Optional: the same sweep in 'full' mode -- adds a layered Green map
% (|G13_eff|) at every wavelength, the response a measurement would see,
% behind the tracked Ritz branches. A few seconds extra here (about 0.2 s per
% wavelength on this stack).
calcModeFull = 'full';   % Ritz + substrate Christoffel + Green
sweepFull = Calc.wavevecAngSweep2Layered(str, angle, wavelengths, vMax, calcModeFull, retentionLevel);
sweepFull = sweepFull.calc();

angleIdx            = 1;        % the only angle
showRitzBranches    = true;
showGreenBackground = true;     % needs a v_phase curveType ('k-vphase' | 'wavelength-vphase')
scaleMode           = 'log';
climRange           = [-3, -0.5];   % log10(|G|/max|G|) colour limits
Vis.plotDispersionCurves(sweepFull, curveType, angleIdx, showRitzBranches, ...
    showGreenBackground, scaleMode, climRange);
title('Al substrate / Cu film -- tracked branches over |G_{13}^{eff}|');
xlim([0 0.03])   % k [rad/nm]
ylim([1 4])      % km/s

% Interactive, Green background on. Gui.combinedViewer takes angular and
% single-angle results only; for a wavelength sweep the dispersion viewer is
% its counterpart (curves linked to the modal shape).
guiFull = Gui.dispersionViewer(sweepFull);
guiFull.show();
guiFull.curveTypeDropdown.Value = curveType;
guiFull.showGreenBgCheckbox.Value = true;
guiFull.onCurveTypeChanged();
guiFull.setAxisRange('x', [0 0.03]);
guiFull.setAxisRange('y', [1 4]);
