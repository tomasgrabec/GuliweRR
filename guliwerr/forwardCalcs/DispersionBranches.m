classdef DispersionBranches
%%% DispersionBranches — tracks each Ritz mode's physical identity across a fixed-angle
%%% wavelength/wavevector sweep (CalcWavevecAngSweep2Layered, CalcWavevecAngSweep2Plate),
%%% so a dispersion curve can be drawn as continuous branches. A separate utility, not
%%% part of the surface-wave classification. See
%%% docs/reference/architecture/forwardCalcs_modeClassification.md, "DispersionBranches.track".
%%%
%%% Value class with Static methods only.
%%%
%%% GuliweRR v6.0

    methods (Static)
        function branchModeNoMap = track(sweepResult, angleIdx)
        %%% Track each Ritz mode's physical identity across a fixed-angle wavelength
        %%% sweep, by maximum eigenvector overlap (MAC) with the previous step
        %%% >> branchModeNoMap = DispersionBranches.track(sweepResult)
        %%% >> branchModeNoMap = DispersionBranches.track(sweepResult, angleIdx)
        %%%
        %%% INPUT:  1. sweepResult .. CalcWavevecAngSweep2Layered or
        %%%                           CalcWavevecAngSweep2Plate (both expose the same
        %%%                           sweeps{iw}.cs{angleIdx}.ritzSolver shape),
        %%%                           already calc()'d with retentionLevel='full'.
        %%%                           That is the only level keeping every Ritz
        %%%                           eigenvector column rather than the classified
        %%%                           ones alone, which matters here because a branch
        %%%                           can drift through classification wobble between
        %%%                           steps
        %%%  (opt) 2. angleIdx      .. which angle to track at (default 1). This stays
        %%%                           a FIXED-ANGLE dispersion-curve case; angleIdx
        %%%                           only picks which slice of a multi-angle
        %%%                           sweepResult to run on
        %%% OUTPUT: branchModeNoMap .. [nWav x nCompModes] double. Column im is one
        %%%                            tracked physical branch; branchModeNoMap(iw,im)
        %%%                            is the raw Ritz mode index at wavelength-step iw
        %%%                            continuing branch im's identity.
        %%%                            branchModeNoMap(1,:) == 1:nCompModes by
        %%%                            definition, and every row is a permutation of
        %%%                            1:nCompModes (each raw mode on exactly one branch)
        %%%
        %%% Algorithm: seed every branch at its own raw index at step 1, then at each
        %%% later step compute one [nCompModes x nCompModes] MAC matrix between the
        %%% current reference eigenvectors and that step's full eigenvector matrix in
        %%% a single matrix product, assign branches to candidates ONE-TO-ONE by
        %%% maximising the total MAC (matchpairs, a linear assignment), and chain the
        %%% reference forward to the matched columns.
        %%%
        %%% MAC here is MASS-WEIGHTED: |a^H sqrt(M_a M_b) b| on the unconstrained
        %%% Legendre coefficients, M being each step's own block-diagonal mass matrix
        %%% (layer factor relDepth*rho), i.e. the overlap in the physical inner product
        %%% integral(rho u^H u dz) (massWeightedUnit). The plain Euclidean overlap of
        %%% the coefficients weights every layer alike, however thin: for a film of
        %%% relDepth 0.001 its coefficients hold ~98 % of a SAW eigenvector's norm but
        %%% <1 % of its mass norm, unrelated modes then reach MAC ~0.98 with it, and
        %%% the assignment handed the SAW branch to a 13 km/s film mode.
        %%%
        %%% Why one-to-one: an independent per-branch argmax would let several branches
        %%% claim one candidate, and a candidate no branch claimed would fall off every
        %%% branch for the rest of the sweep. Inside a dense, near-degenerate cluster (a
        %%% layered stack's substrate box modes above its bulk shear velocity) that
        %%% happens at almost every step, so a guided mode emerging from the cluster as k
        %%% grows would never be plotted. Where the argmax gives distinct matches the
        %%% assignment is the same. Simultaneous angle-and-wavelength tracking stays out of scope:
        %%% MAC-chaining across a genuine 2D grid needs a real 2D continuity scheme,
        %%% not a choice of which 1D slice to chain along.
            validateattributes(sweepResult, {'CalcWavevecAngSweep2Layered', 'CalcWavevecAngSweep2Plate'}, ...
                {'scalar'}, 'DispersionBranches.track', 'sweepResult');
            if nargin < 2 || isempty(angleIdx); angleIdx = 1; end
            validateattributes(angleIdx, {'double'}, ...
                {'scalar', 'integer', 'positive', '<=', numel(sweepResult.angles)}, ...
                'DispersionBranches.track', 'angleIdx');

            nWav = numel(sweepResult.sweeps);
            ritz1 = sweepResult.sweeps{1}.cs{angleIdx}.ritzSolver;
            nCompModes = ritz1.nCompModes;
            if isempty(ritz1.eigVecs_con)
                error('DispersionBranches:track:missingEigVecs', ...
                    'sweepResult.eigVecs_con is empty -- DispersionBranches.track needs the full Ritz eigenvector matrix at every wavelength step; construct sweepResult with retentionLevel=''full''.');
            end

            branchModeNoMap = zeros(nWav, nCompModes);
            branchModeNoMap(1, :) = 1:nCompModes;

            refUnit = DispersionBranches.massWeightedUnit(ritz1, nCompModes);
            for iw = 2:nWav
                ritzIw = sweepResult.sweeps{iw}.cs{angleIdx}.ritzSolver;
                if isempty(ritzIw.eigVecs_con)
                    error('DispersionBranches:track:missingEigVecs', ...
                        'sweepResult.eigVecs_con is empty at wavelength step %d -- DispersionBranches.track needs the full Ritz eigenvector matrix at every wavelength step; construct sweepResult with retentionLevel=''full''.', iw);
                end
                candUnit = DispersionBranches.massWeightedUnit(ritzIw, nCompModes);

                MAC = abs(refUnit' * candUnit);          % [nCompModes x nCompModes], rows = branches, cols = candidates
                %%% one-to-one: minimise total (1 - MAC); the unmatched cost (above any
                %%% possible 1 - MAC) forces a complete assignment of the square matrix
                pairs = matchpairs(1 - MAC, nCompModes + 1);
                matchIdx = zeros(nCompModes, 1);
                matchIdx(pairs(:, 1)) = pairs(:, 2);
                branchModeNoMap(iw, :) = matchIdx.';

                refUnit = candUnit(:, matchIdx);          % chain forward
            end
        end

        function [vels, modeNos, seedModeNo] = trackSeededBranch(sweepResult, angleIdx)
        %%% Velocity of the branch seeded at the slowest accepted 'SAW' entry of the
        %%% first wavelength step, followed through the sweep by DispersionBranches.track
        %%% >> vels = DispersionBranches.trackSeededBranch(sweepResult)
        %%% >> [vels, modeNos, seedModeNo] = DispersionBranches.trackSeededBranch(sweepResult, angleIdx)
        %%%
        %%% INPUT:  1. sweepResult .. CalcWavevecAngSweep2Layered, calc()'d with
        %%%                           retentionLevel='full' (as for track)
        %%%  (opt) 2. angleIdx      .. which angle of the sweep (default 1)
        %%% OUTPUT: vels       .. [nWav x 1] Ritz velocity of the seeded branch at every
        %%%                       wavelength step [km/s]
        %%%         modeNos    .. [nWav x 1] its raw Ritz mode index at every step
        %%%         seedModeNo .. the seed entry's modeNo at step 1
        %%%
        %%% The seed is the accepted 'SAW' entry with the slowest bound Ritz mode (modeVel) of
        %%% sweepResult.sweeps{1}.significantModesMap{angleIdx}, i.e. the classification
        %%% under the sweep's active method and standard at the first wavelength step (the
        %%% Rayleigh-type wave when the steps run from thin to thick film). Errors if
        %%% that step has no accepted 'SAW' entry.
            validateattributes(sweepResult, {'CalcWavevecAngSweep2Layered'}, {'scalar'}, ...
                'DispersionBranches.trackSeededBranch', 'sweepResult');
            if nargin < 2 || isempty(angleIdx); angleIdx = 1; end
            branchModeNoMap = DispersionBranches.track(sweepResult, angleIdx);

            sm = sweepResult.sweeps{1}.significantModesMap{angleIdx};
            if ~isempty(sm)
                sm = sm(strcmp({sm.type}, 'SAW'));
            end
            if isempty(sm)
                error('DispersionBranches:trackSeededBranch:noSeed', ...
                    'The first wavelength step has no accepted ''SAW'' entry at angle index %d to seed the branch from.', angleIdx);
            end
            [~, iSeed] = min(ModeClassification.modeVels(sm));   % the bound Ritz mode, the quantity tracked (vel is the Green peak)
            seedModeNo = sm(iSeed).modeNo;

            nWav    = numel(sweepResult.sweeps);
            modeNos = branchModeNoMap(:, seedModeNo);
            vels    = zeros(nWav, 1);
            for iw = 1:nWav
                vels(iw) = sweepResult.sweeps{iw}.cs{angleIdx}.ritzSolver.vels(modeNos(iw));
            end
        end
    end

    methods (Static, Access = private)
        function X = massWeightedUnit(ritz, nCompModes)
        %%% First nCompModes eigenvectors in mass-weighted, unit-norm form
        %%% >> X = DispersionBranches.massWeightedUnit(ritz, nCompModes)
        %%%
        %%% X = sqrt(w) .* (N * eigVecs_con), w = each Legendre coefficient's layer
        %%% factor (GammaMatrix.LayerNormFac), i.e. the diagonal of the unconstrained
        %%% M. X'*X is then the overlap in the physical inner product
        %%% integral(rho u^H u dz), as M_con is. Each step uses its own w, so the
        %%% cross-step product X_a' * X_b carries the geometric mean sqrt(M_a M_b):
        %%% symmetric in the two steps and <= 1 without extra normalisation.
            s   = ritz.structure;
            bVL = [s.calcParams.basisVecL];
            w   = repelem(arrayfun(@(il) GammaMatrix.LayerNormFac(s, il), 1:s.nLay), bVL).';
            X   = sqrt(w) .* (ritz.nullConstrMat * ritz.eigVecs_con(:, 1:nCompModes));
            X   = X ./ vecnorm(X);   % already 1 for M-orthonormal vectors; guards round-off
        end
    end
end
