%% Tutorial 09 (code) -- inversion: recover Cu constants from synthetic SAW/pSAW/limiting velocities
% Narrated version: examples/tutorials/tutorial_09_inversion_synthetic.m. Run guliwerr/init.m first.
% EXPERIMENTAL: the inverse module is validated on synthetic data only; the interface may change.

%% Ground truth, deliberately wrong initial guess, and a Cu {001} plane
dirMat        = [1, 0, 0; 0, 1, 0];   % normal [100], reference [010] (cubic: same as the default)
matTrue       = Material.FromLibrary('Cu', dirMat);   % "true" Cu from the materials library
cuConstsTrue  = matTrue.elConsts.vector;   % [C11, C12, C44] GPa
density       = matTrue.rho;               % g/cm^3
cuConstsGuess = cuConstsTrue .* [1.10, 0.90, 1.08];   % +10%, -10%, +8%
symmNo        = 2;                 % cubic

%% Discretisation -- "easy" tier (Zener ratio ~3.2, principal (001) plane)
relDepth   = 20;
legOrder   = 50;
calcParams = CalcParams(relDepth, legOrder);
anglesDeg  = (0:2:180)';
nMaxModes  = 500;

%% Step 1 -- forward sweep of the true material
calcMode       = 'full';           % Green needed for pseudoRayleigh points
retentionLevel = 'full';           % Inv.association needs every eigenvector
strTrue   = Structure(matTrue, calcParams, 'cuTrue', nMaxModes);
sweepTrue = Calc.angSweepHS(strTrue, anglesDeg, calcMode, retentionLevel);
sweepTrue = sweepTrue.calc();

%% Step 2a -- harvest Rayleigh and the most prominent pseudoRayleigh per angle
rayAngle    = zeros(0, 1); rayVel    = zeros(0, 1);
pseudoAngle = zeros(0, 1); pseudoVel = zeros(0, 1);
for iAngle = 1:numel(anglesDeg)
    sigModes = sweepTrue.significantModesMap{iAngle};
    idxRay = find(strcmp({sigModes.type}, 'Rayleigh'), 1);
    if ~isempty(idxRay)
        rayAngle(end+1, 1) = sweepTrue.angles(iAngle);   %#ok<AGROW>
        rayVel(end+1, 1)   = sigModes(idxRay).modeVel;   %#ok<AGROW>
    end
    idxPseudo = find(strcmp({sigModes.type}, 'pseudoRayleigh'));
    if ~isempty(idxPseudo)
        [~, iBest] = max([sigModes(idxPseudo).relPeakHeight]);   % tallest Green peak wins
        pseudoAngle(end+1, 1) = sweepTrue.angles(iAngle);        %#ok<AGROW>
        pseudoVel(end+1, 1)   = sigModes(idxPseudo(iBest)).modeVel;  %#ok<AGROW>
    end
end

%% Step 2b -- Cu(001) Rayleigh-observability gap, longitudinal points, dropout, noise
gapAngleLow  = 36;                 % deg -- no Rayleigh velocity kept here (pseudoRayleigh range)
gapAngleHigh = 54;                 % deg
keepRay  = ~(rayAngle >= gapAngleLow & rayAngle <= gapAngleHigh);
rayAngle = rayAngle(keepRay);
rayVel   = rayVel(keepRay);

longAngleStep = 10;                % deg
longAngle     = anglesDeg(mod(anglesDeg, longAngleStep) == 0);
[~, longIdx]  = ismember(longAngle, sweepTrue.angles);
longVel       = sweepTrue.longVelMap(longIdx);

allAngle = [rayAngle; pseudoAngle; longAngle(:)];
allVel   = [rayVel; pseudoVel; longVel(:)];
allType  = [repmat({'rayleigh'}, numel(rayAngle), 1); ...
            repmat({'pseudoRayleigh'}, numel(pseudoAngle), 1); ...
            repmat({'longitudinal'}, numel(longAngle), 1)];

randomSeed      = 20260911;
dropoutFraction = 0.15;
noiseSigma      = 0.005;           % km/s, assumed noise level
rng(randomSeed, 'twister');
nAll     = numel(allAngle);
keepMask = true(nAll, 1);
keepMask(randperm(nAll, round(dropoutFraction * nAll))) = false;
allAngle = allAngle(keepMask);
allVel   = allVel(keepMask) + noiseSigma * randn(sum(keepMask), 1);
allType  = allType(keepMask);

measured = Inv.measured(allAngle, allVel, allType, noiseSigma);   % waveType list: MeasuredVelocities.WAVE_TYPES
fprintf('%d measured points\n', measured.nPoints);

%% Step 3 -- forward sweep at the guess, association built against IT
ecGuess    = ElastConsts(cuConstsGuess, symmNo);
matGuess   = Material(ecGuess, density, dirMat);
strGuess   = Structure(matGuess, calcParams, 'cuGuess', nMaxModes);
sweepGuess = Calc.angSweepHS(strGuess, anglesDeg, calcMode, retentionLevel);
sweepGuess = sweepGuess.calc();

velocityTolerance = 0.5;           % km/s, widened from the 0.3 default for a far guess
angleTolerance    = 1e-6;          % deg
association = Inv.association(sweepGuess, measured, velocityTolerance, angleTolerance);

%% Step 4 -- sensitivity: which constant combinations the data determine
invertedLayerNo = 1;
sensitivity = Inv.sensitivity(strGuess, association, invertedLayerNo, measured);
fprintf('%d of %d combinations well-determined\n', sensitivity.nOptimisedConstants, sensitivity.nIndepConsts);
Vis.plotSensitivitySpectrum(sensitivity);

%% Step 5 -- inversion (Name-Value options: Stepwise, AnchorPolicy, IncludePseudoSAW, PseudoSAWRefinement)
stepwise = true;                   % true: grow the optimised set one at a time | false: all at once
inversion = Inv.elasticity(strGuess, measured, association, sensitivity, invertedLayerNo, 'Stepwise', stepwise);
inversion = inversion.calc();

constNames = {'c11', 'c12', 'c44'};
for iConst = 1:numel(constNames)
    fprintf('%s: true %.2f, guess %.2f, recovered %.2f +/- %.2f GPa\n', constNames{iConst}, ...
        cuConstsTrue(iConst), cuConstsGuess(iConst), ...
        inversion.result.voigtElasticVec(iConst), inversion.result.voigtElasticVecErrors(iConst));
end
rmse = sqrt(inversion.result.costFunctionValue / inversion.result.nDataPoints);
fprintf('RMSE %.4f km/s (noise %.4f km/s)\n', rmse, noiseSigma);

%% Step 5b -- second pass from the recovered constants (fixes points matched to the wrong branch from the far guess)
constsPass1 = inversion.result.voigtElasticVec(:)';
strPass2    = Structure(Material(ElastConsts(constsPass1, symmNo), density, dirMat), calcParams, 'cuPass2', nMaxModes);
sweepPass2  = Calc.angSweepHS(strPass2, anglesDeg, calcMode, retentionLevel);
sweepPass2  = sweepPass2.calc();
associationPass2 = Inv.association(sweepPass2, measured, velocityTolerance, angleTolerance);
sensitivityPass2 = Inv.sensitivity(strPass2, associationPass2, invertedLayerNo, measured);
inversion = Inv.elasticity(strPass2, measured, associationPass2, sensitivityPass2, invertedLayerNo, 'Stepwise', stepwise);
inversion = inversion.calc();
for iConst = 1:numel(constNames)
    fprintf('%s: true %.2f, pass 1 %.2f, pass 2 %.2f +/- %.2f GPa\n', constNames{iConst}, ...
        cuConstsTrue(iConst), constsPass1(iConst), ...
        inversion.result.voigtElasticVec(iConst), inversion.result.voigtElasticVecErrors(iConst));
end
rmse = sqrt(inversion.result.costFunctionValue / inversion.result.nDataPoints);
fprintf('RMSE %.4f km/s (noise %.4f km/s)\n', rmse, noiseSigma);

%% Step 6 -- fit quality and mode-tracking diagnostic
Vis.plotInversionFit(measured, inversion.association, association);   % association = initial-guess snapshot
Vis.plotAssociationDiagnostic(inversion.association);

%% Step 7 -- Green maps: true vs. recovered, with the association overlaid
vMaxRecovered  = [];               % optional (layered stacks only)
sweepRecovered = inversion.recalcForwardAtResult(anglesDeg, vMaxRecovered, calcMode, retentionLevel);
scaleMode   = [];                  % [] = defaults
climRange   = [];
plotType    = [];
overlayMode = true;                % markers only, keep the Green colormap
figure('Position', [100 100 1100 850]);
sweepsToShow = {sweepTrue, sweepRecovered};
panelTitles  = {'True material', 'Recovered material'};
for iPanel = 1:2
    axTop = subplot(2, 2, iPanel);
    Vis.plotGreenMap(sweepsToShow{iPanel}, scaleMode, climRange, plotType, axTop);
    title(axTop, panelTitles{iPanel});
    axBottom = subplot(2, 2, iPanel + 2);
    Vis.plotGreenMap(sweepsToShow{iPanel}, scaleMode, climRange, plotType, axBottom);
    Vis.plotAssociationDiagnostic(inversion.association, axBottom, overlayMode);
    title(axBottom, [panelTitles{iPanel}, ' + association']);
end
