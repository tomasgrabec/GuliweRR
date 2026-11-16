%% Tutorial: elastic-constant inversion from synthetic SAW/limiting velocities
%
% The INVERSE problem: given measured surface-wave velocities -- here
% synthesised from a known material, so the answer is known -- recover the
% elastic constants that produced them.
%
% EXPERIMENTAL: the inverse module is validated on synthetic data only, and
% its interface may change. Its pseudo-surface-wave tracking is new and not
% yet published.
%
% Copper, cubic, on a {001} plane (surface normal [100], propagation
% reference [010]; for a cubic crystal equivalent to Material's default
% dirMat). The "measurement" is Cu's own Rayleigh, pseudo-Rayleigh and
% longitudinal limiting velocities with noise added; with real TGS data
% (Inv.measured / MeasuredVelocities.fromTable) everything downstream is the same.
%
% Run run('guliwerr/init.m') from the project root first (tutorial_01 says why).

%% Ground truth and initial (deliberately wrong) guess
dirMat     = [1, 0, 0; 0, 1, 0];   % normal [100], reference [010]
matTrue    = Material.FromLibrary('Cu', dirMat);   % the "true" Cu: materials library
cuVecTrue  = matTrue.elConsts.vector;              % [c11 c12 c44] GPa
cuRho      = matTrue.rho;                          % g/cm^3
cuVecGuess = cuVecTrue .* [1.10, 0.90, 1.08];      % +10%, -10%, +8% -- a deliberately wrong starting point
symmNo     = 2;                                    % cubic

relDepth = 20;   % depth of the computational domain [wavelengths]
legOrder = 50;   % Legendre basis order; "easy" tier (docs/pitfalls.md MTH1):
                 % Cu's Zener ratio 2*c44/(c11-c12) ~ 3.2 on a principal plane
calcParams = CalcParams(relDepth, legOrder);
angleGrid  = (0:2:180)';

%% Step 1 -- forward-calculate the "true" material, to synthesise a measurement from
strTrue   = Structure(matTrue, calcParams, 'cuTrue', 500);
calcMode       = 'full';   % Green is needed: pseudoRayleigh points come from Green peaks
retentionLevel = 'full';   % Inv.association needs every Ritz eigenvector
sweepTrue = Calc.angSweepHS(strTrue, angleGrid, calcMode, retentionLevel);
sweepTrue = sweepTrue.calc();

%% Step 2 -- synthesise the "measured" dataset: Rayleigh and pseudo-Rayleigh
% velocities wherever each is classified, longitudinal limiting velocities
% every 10 deg, a Cu(001) observability gap (no Rayleigh velocity kept at
% 36-54 deg, an illustrative range where the pseudo-Rayleigh wave takes over in
% the calculation), a further random 15% dropout, and 0.005 km/s Gaussian
% noise (an assumed value for this example).
%
% The pseudo-Rayleigh points matter: over the 36-54 deg gap (and its mirror,
% 126-144 deg) they are the only surface-wave data, which is why
% IncludePseudoSAW defaults to true in Inv.sensitivity/Inv.elasticity. Where
% one branch ends and the other starts is a threshold on a continuum
% (docs/pitfalls.md PHY1).
% modeVel: the Ritz eigenvelocity of each entry's mode, the quantity the inversion
% fits (an entry's vel is its Green peak velocity).
rayAngle = zeros(0,1); rayVel = zeros(0,1);
pseudoAngle = zeros(0,1); pseudoVel = zeros(0,1);
for ia = 1:numel(angleGrid)
    sm = sweepTrue.significantModesMap{ia};

    idxR = find(strcmp({sm.type}, 'Rayleigh'), 1);
    if ~isempty(idxR)
        rayAngle(end+1,1) = sweepTrue.angles(ia); %#ok<AGROW>
        rayVel(end+1,1)   = sm(idxR).modeVel; %#ok<AGROW>
    end

    % Several pseudoRayleigh entries can coexist at one angle: take the most
    % prominent, as an experimenter reading the map by eye would
    idxP = find(strcmp({sm.type}, 'pseudoRayleigh'));
    if ~isempty(idxP)
        [~, iBest] = max([sm(idxP).relPeakHeight]);
        pseudoAngle(end+1,1) = sweepTrue.angles(ia); %#ok<AGROW>
        pseudoVel(end+1,1)   = sm(idxP(iBest)).modeVel; %#ok<AGROW>
    end
end
gapAngleLow  = 36;   % Cu(001) Rayleigh-observability gap [deg]: no Rayleigh
gapAngleHigh = 54;   % velocity is kept here
keepMask = ~(rayAngle >= gapAngleLow & rayAngle <= gapAngleHigh);
rayAngle = rayAngle(keepMask);
rayVel   = rayVel(keepMask);

longitudinalAngleStep = 10;   % [deg] spacing of the synthetic longitudinal-limiting points
longAngle = angleGrid(mod(angleGrid, longitudinalAngleStep) == 0);
[~, longIdx] = ismember(longAngle, sweepTrue.angles);
longVel = sweepTrue.longVelMap(longIdx);

allAngle = [rayAngle; pseudoAngle; longAngle(:)];
allVel   = [rayVel;   pseudoVel;   longVel(:)];
allType  = [repmat({'rayleigh'}, numel(rayAngle), 1); ...
            repmat({'pseudoRayleigh'}, numel(pseudoAngle), 1); ...
            repmat({'longitudinal'}, numel(longAngle), 1)];
nCandidateBeforeDropout = numel(allAngle);   % after the Cu(001) gap, before the random dropout

randomSeed            = 20260911;   % fixed for a reproducible tutorial run
randomDropoutFraction = 0.15;       % fraction of the remaining points dropped at random
noiseSigma = 0.005;   % [km/s] additive Gaussian noise; also the stored velocityStd
                      % and the noise floor the fit's RMSE is compared with

rng(randomSeed, 'twister');
nAll = numel(allAngle);
dropIdx = randperm(nAll, round(randomDropoutFraction*nAll));
keepMask2 = true(nAll, 1);
keepMask2(dropIdx) = false;
allAngle = allAngle(keepMask2);
allVel   = allVel(keepMask2) + noiseSigma*randn(sum(keepMask2), 1);
allType  = allType(keepMask2);

measured = Inv.measured(allAngle, allVel, allType, noiseSigma);
fprintf('Synthesised %d measured points (%d dropped at random from the %d that survived the Cu(001) gap):\n', ...
    measured.nPoints, nCandidateBeforeDropout - measured.nPoints, nCandidateBeforeDropout);
fprintf('  rayleigh %d, pseudoRayleigh %d, longitudinal %d\n', ...
    sum(strcmp(measured.waveType, 'rayleigh')), sum(strcmp(measured.waveType, 'pseudoRayleigh')), ...
    sum(strcmp(measured.waveType, 'longitudinal')));

%% Step 3 -- forward-calculate at the (wrong) initial guess and build the
% initial mode association against it -- associating against the truth
% would smuggle the answer into the input.
strGuess   = Structure(Material(ElastConsts(cuVecGuess, symmNo), cuRho, dirMat), calcParams, 'cuGuess', 500);
sweepGuess = Calc.angSweepHS(strGuess, angleGrid, calcMode, retentionLevel);
sweepGuess = sweepGuess.calc();

velocityTolerance = 0.5;   % [km/s] max |v_calc - v_meas| for a mode match; wider than the
                           % 0.3 km/s default, which would drop most Rayleigh points this far off
angleTolerance = 1e-6;     % [deg] max angle mismatch -- tight, both use the same angleGrid
association = Inv.association(sweepGuess, measured, velocityTolerance, angleTolerance);

%% Step 4 -- sensitivity analysis: which combinations of elastic constants
% can this dataset actually determine?
invertedLayerNo = 1;   % which layer's constants to recover (1 = the half-space itself)
sensitivity = Inv.sensitivity(strGuess, association, invertedLayerNo, measured);
fprintf('Sensitivity: %d of %d transformed combinations judged well-determined ', ...
    sensitivity.nOptimisedConstants, sensitivity.nIndepConsts);
fprintf('(cutoff = 1/%.0f error amplification).\n', 1/sensitivity.sensitivityCutoff);

Vis.plotSensitivitySpectrum(sensitivity);

%% Step 5 -- run the inversion
% Name-Value options (all optional): 'Stepwise' (true: grow the optimised set one
% combination at a time | false: all at once), 'AnchorPolicy' ('chained' | 'fixed'),
% 'IncludePseudoSAW' (true), 'PseudoSAWRefinement' ('green' | 'none').
inversion = Inv.elasticity(strGuess, measured, association, sensitivity, invertedLayerNo);
inversion = inversion.calc();

fprintf('\n%8s %12s %12s %12s\n', '', 'true', 'guess', 'recovered +/- error');
labels = {'c11', 'c12', 'c44'};
for i = 1:3
    fprintf('%8s %12.3f %12.3f %8.3f +/- %.3f GPa\n', labels{i}, cuVecTrue(i), cuVecGuess(i), ...
        inversion.result.voigtElasticVec(i), inversion.result.voigtElasticVecErrors(i));
end
fprintf('\nFinal cost function value: %.4g (RMSE %.4f km/s, noise floor %.4f km/s)\n', ...
    inversion.result.costFunctionValue, ...
    sqrt(inversion.result.costFunctionValue/inversion.result.nDataPoints), noiseSigma);

%% Step 5b -- second pass, restarted from the recovered constants
% From a far guess, Step 3 matches some measured points to the wrong branch, and
% the fit then stops short of the true minimum (RMSE well above the noise).
% Re-running the forward sweep and the association at the first result, then
% inverting again, corrects those matches. Do the same with real data: repeat
% until the constants stop changing and the RMSE reaches the noise level.
vecPass1 = inversion.result.voigtElasticVec(:)';
strPass2 = Structure(Material(ElastConsts(vecPass1, symmNo), cuRho, dirMat), calcParams, 'cuPass2', 500);
sweepPass2 = Calc.angSweepHS(strPass2, angleGrid, calcMode, retentionLevel);
sweepPass2 = sweepPass2.calc();
associationPass2 = Inv.association(sweepPass2, measured, velocityTolerance, angleTolerance);
sensitivityPass2 = Inv.sensitivity(strPass2, associationPass2, invertedLayerNo, measured);
inversion = Inv.elasticity(strPass2, measured, associationPass2, sensitivityPass2, invertedLayerNo);
inversion = inversion.calc();

fprintf('\nSecond pass:\n%8s %12s %12s %12s\n', '', 'true', 'pass 1', 'recovered +/- error');
for i = 1:3
    fprintf('%8s %12.3f %12.3f %8.3f +/- %.3f GPa\n', labels{i}, cuVecTrue(i), vecPass1(i), ...
        inversion.result.voigtElasticVec(i), inversion.result.voigtElasticVecErrors(i));
end
fprintf('\nFinal cost function value: %.4g (RMSE %.4f km/s, noise floor %.4f km/s)\n', ...
    inversion.result.costFunctionValue, ...
    sqrt(inversion.result.costFunctionValue/inversion.result.nDataPoints), noiseSigma);

% Mode-tracking health per binding class -- check it whenever pseudo-Rayleigh
% data are included: a mis-tracked pseudo-SAW point can look perfectly
% confident on its Ritz MAC alone. greenAgreement* reports the Green's-function
% cross-check that catches this.
trk = inversion.result.tracking;
fprintf('\nMode-tracking health (%s pseudo-SAW refinement):\n', trk.pseudoSAWRefinement);
fprintf('  %-12s %8s %14s %10s\n', 'class', 'nPoints', 'minConfidence', 'nBelowThr');
for cls = ["ritzMode", "pseudoSAW", "limiting"]
    c = char(cls);
    if trk.nPointsPerClass.(c) == 0; continue; end
    fprintf('  %-12s %8d %14.4f %10d\n', c, trk.nPointsPerClass.(c), trk.minConfidence.(c), ...
        trk.nBelowThresholdPerClass.(c));
end
fprintf('  collisions=%d, limitingRootsLost=%d, excludedAtFinalEval=%d\n', ...
    trk.nModeCollisions, trk.nLimitingRootsLost, trk.nExcludedAtFinalEval);
if strcmp(trk.pseudoSAWRefinement, 'green')
    fprintf('  Green rescans=%d, pseudo-SAW/Green agreement: max %.4f, median %.4f\n', ...
        trk.nGreenRescans, trk.greenAgreementMax, trk.greenAgreementMedian);
end

%% Step 6 -- fit quality and mode-tracking diagnostic
% One plot, every stage, dropped points as grey x. `association` (Step 3) is
% a value object, so inversion.calc() leaves it at the initial guess;
% inversion.association holds the converged binding.
Vis.plotInversionFit(measured, inversion.association, association);
title('Measured, initial guess, and recovered fit -- all in one');

Vis.plotAssociationDiagnostic(inversion.association);

%% Step 7 -- Green's-function maps: true vs. recovered, with the tracking overlaid
% Top row: the |G13_eff| map of the true and of the recovered material --
% nearly identical if the inversion converged. Bottom row: the same maps with
% every associated point overlaid (overlay mode: markers only, the map keeps
% its colormap), coloured by tracking confidence; grey diamonds are limiting
% (Christoffel) points. Comparing the two bottom panels shows whether every
% point, pseudo-SAW included, sits on the same ridge in both.
%
% recalcForwardAtResult() costs one extra forward sweep at the inverted
% elasticity; the optimisation itself never runs one.
vMaxRecovered  = [];   % optional (layered stacks only): extends the Green axis
sweepRecovered = inversion.recalcForwardAtResult(angleGrid, vMaxRecovered, calcMode, retentionLevel);

figure('Position', [100 100 1100 850]);

axTrueTop = subplot(2, 2, 1);
Vis.plotGreenMap(sweepTrue, [], [], [], axTrueTop);
title(axTrueTop, 'True material');

axRecTop = subplot(2, 2, 2);
Vis.plotGreenMap(sweepRecovered, [], [], [], axRecTop);
title(axRecTop, 'Recovered material');

axTrueBottom = subplot(2, 2, 3);
Vis.plotGreenMap(sweepTrue, [], [], [], axTrueBottom);
Vis.plotAssociationDiagnostic(inversion.association, axTrueBottom, true);
title(axTrueBottom, 'True material + association overlay');

axRecBottom = subplot(2, 2, 4);
Vis.plotGreenMap(sweepRecovered, [], [], [], axRecBottom);
Vis.plotAssociationDiagnostic(inversion.association, axRecBottom, true);
title(axRecBottom, 'Recovered material + association overlay');

sgtitle('Green''s-function maps: true vs. recovered, with mode tracking overlaid');

%% What to try next
% - Loosen/tighten cuVecGuess to see how far a starting guess can be from the
%   truth before mode association or convergence starts failing.
% - Try a general (non-principal) crystallographic plane: build dirMat from
%   an arbitrary surface normal and CalcParams(30, 70) instead. The
%   recovered crystal-frame constants should agree with this script's,
%   independent of which plane the "measurement" was cut on.
% - Add measured.addPhaseVelocity(...)/measured.setBulkModulus(...) for the
%   two auxiliary data channels, and re-run from Step 3.
