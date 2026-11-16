classdef SensitivityAnalysis
%%% SensitivityAnalysis — the sensitivity/Jacobian matrix G, the Gauss-Newton Hessian
%%%   S = 2*G'*G, its eigendecomposition (transformMatrix T, sensitivityEigenvalues psi),
%%%   and the cutoff deciding how many transformed constants are well-determined enough
%%%   to optimise.
%%%
%%% EXPERIMENTAL: validated on synthetic data only; the interface may change.
%%%
%%%   Row k of G is one data point (an associated Ritz mode, a limiting bulk velocity, a
%%%   pulse-echo phase velocity, or the bulk-modulus pseudo-datum). Column i is one
%%%   independent elastic constant of the inverted layer, in ElastConsts.vector order.
%%%   G_ki = d(v_k)/d(c_i), evaluated at the initial guess.
%%%
%%%   The Ritz-mode rows use the analytic Rayleigh-quotient formula
%%%     G_ki = Re[ alpha_k' (N' dK/dc_i N) alpha_k ] / (8 pi^2 v_k)
%%%   with GammaMatrix.layerBlockDerivs for dK/dc_i -- the same mathematical pattern as
%%%   RitzSolver.calcGroupVelocity's dw2/dD = alpha' (dK/dD) alpha, differentiating w.r.t.
%%%   an elastic constant instead of the domain depth.
%%%
%%%   The limiting/phase-velocity/bulk-modulus rows use a central finite difference on the
%%%   SAME validated forward-model calls (CalcCommon.limitingAndBulk, Christoffel.matrix,
%%%   ElastConsts.Voigt2Bulk) rather than a hand-derived analytic derivative -- a deliberate
%%%   simplification (Christoffel's own Hellmann-Feynman derivative helpers are private
%%%   and tied to a numerical-regularisation scheme not reused here). Cheap: one 3x3
%%%   Christoffel eigenproblem (or small root-find) per constant, negligible next to the
%%%   Ritz eigensolves.
%%%
%%%   JOINT MULTI-DATASET INVERSION (several crystallographic cuts of the SAME crystal,
%%%   sharing one target elasticity): `structure`, `association` and `measured` in calc()
%%%   are accepted as OBJECT ARRAYS, one element per dataset d = 1..nDatasets (length 1 is
%%%   the ordinary single-dataset case). Every dataset must share the SAME invertedLayerNo
%%%   LAYER INDEX, the SAME symmClassNo and the SAME nIndepConsts for that layer (checked);
%%%   everything else per dataset (orientation, other layers, CalcParams/legOrder/relDepth)
%%%   may differ freely. The per-dataset G_d blocks are vertically stacked
%%%   G = [G_1; G_2; ...; G_D] BEFORE the single pooled SVD. Row order: all of dataset 1's
%%%   rows (Ritz-class, then limiting, then phaseVelocities, then bulkModulus), then all of
%%%   dataset 2's rows in the same sub-order, etc. -- pointWeights must be supplied in that
%%%   same concatenated order. The canonical initial guess (initialVoigtVec/
%%%   initialTransformedVec) is taken from DATASET 1's stored elConsts for the inverted
%%%   layer; every other dataset's stored elConsts for that layer is unused by calc() (only
%%%   its dirMat/rho/other layers matter), though physically it should be the same crystal.
%%%
%%% Value class: computed once and then read.
%%%
%%% GuliweRR v6.0

    properties (Constant)
        DEFAULT_SENSITIVITY_CUTOFF = 1/160   % default of SensitivityCutoff (see the property)
    end

    properties (SetAccess = protected)
        nIndepConsts(1,1)         double = 0     % L
        symmClassNo(1,1)          double = 0
        nDatasets(1,1)            double = 1     % D; number of jointly-inverted datasets
                                                   % (crystallographic cuts) pooled into G;
                                                   % 1 for an ordinary single-dataset call

        initialVoigtVec           double = []    % c at the initial guess [L x 1]
        initialTransformedVec     double = []    % T' * initialVoigtVec   [L x 1]

        sensitivityMatrix         double = []    % G   [N_data x L], (km/s)/GPa
        sensitivityEigenvalues    double = []    % psi, descending [L x 1]
        transformMatrix           double = []    % T, orthogonal [L x L]
        normalisedSensitivityRoots double = []   % sqrt(psi/psi(1)) [L x 1]

        sensitivityCutoff(1,1)    double = SensitivityAnalysis.DEFAULT_SENSITIVITY_CUTOFF
                                                   % keep p while normalisedSensitivityRoots(p)
                                                   % >= sensitivityCutoff (error amplification
                                                   % <= 160x); an empirical choice
                                                   % (Grabec et al. 2024), not derived
        nOptimisedConstants(1,1)  double = 0     % automatic cutoff P

        pointWeights              double = []    % [N_data x 1]; default ones = unweighted
        includePseudoSAW(1,1)     logical = true  % pseudo-Rayleigh points contribute rows to G
                                                   % by DEFAULT: a TGS angular map can be carried
                                                   % by the pseudo-SAW branch alone over wide
                                                   % angle ranges, and excluding it throws away
                                                   % the only data there

        nRitzRows(1,1)            double = 0
        nLimitingRows(1,1)        double = 0
        nPhaseRows(1,1)           double = 0
        hasBulkModulusRow(1,1)    logical = false
    end

    methods
        function obj = SensitivityAnalysis()
            %%% Initialise an empty SensitivityAnalysis
            %%% >> obj = SensitivityAnalysis()
        end
    end

    methods (Static)
        function obj = calc(structure, association, invertedLayerNo, measured, varargin)
            %%% Build the sensitivity matrix and its eigendecomposition
            %%% >> obj = SensitivityAnalysis.calc(structure, association, invertedLayerNo, measured, ...)
            %%%
            %%% INPUT: 1. structure       .. Structure at the INITIAL GUESS; scalar for a single
            %%%                             dataset, or a [1 x nDatasets]/[nDatasets x 1] OBJECT
            %%%                             ARRAY for a JOINT multi-dataset inversion (several
            %%%                             crystallographic cuts of the SAME crystal) -- see the
            %%%                             class docblock.
            %%%        2. association    .. ModeAssociation (or array, one per dataset, same
            %%%                             length as structure), built at the same initial guess
            %%%        3. invertedLayerNo .. index of the layer whose elasticity is inverted;
            %%%                             SHARED across every dataset (same layer position)
            %%%        4. measured        .. MeasuredVelocities (or array, one per dataset) --
            %%%                             for phaseVelocities/bulkModulus
            %%% Name-value: PointWeights (default ones), IncludePseudoSAW (default TRUE --
            %%%                 see the includePseudoSAW property; pass false to exclude
            %%%                 pseudo-Rayleigh points from G),
            %%%             UseParallel ("auto"|"always"|"never", default "auto"),
            %%%             SensitivityCutoff (default 1/160)

            p = inputParser;
            addParameter(p, 'PointWeights', []);
            addParameter(p, 'IncludePseudoSAW', true);
            addParameter(p, 'UseParallel', "auto");
            addParameter(p, 'SensitivityCutoff', SensitivityAnalysis.DEFAULT_SENSITIVITY_CUTOFF);
            parse(p, varargin{:});
            opts = p.Results;
            if ~isempty(opts.PointWeights)
                validateattributes(opts.PointWeights, {'double'}, ...
                    {'vector', 'real', 'finite', 'nonnegative'}, 'SensitivityAnalysis.calc', 'PointWeights');
            end
            validateattributes(opts.IncludePseudoSAW, {'logical', 'numeric'}, {'scalar', 'binary'}, ...
                'SensitivityAnalysis.calc', 'IncludePseudoSAW');
            opts.IncludePseudoSAW = logical(opts.IncludePseudoSAW);
            validateattributes(opts.SensitivityCutoff, {'double'}, ...
                {'scalar', 'real', 'finite', 'positive', '<=', 1}, 'SensitivityAnalysis.calc', 'SensitivityCutoff');

            validateattributes(structure, {'Structure'}, {'nonempty'}, ...
                'SensitivityAnalysis.calc', 'structure');
            validateattributes(association, {'ModeAssociation'}, {'nonempty'}, ...
                'SensitivityAnalysis.calc', 'association');
            validateattributes(measured, {'MeasuredVelocities'}, {'nonempty'}, ...
                'SensitivityAnalysis.calc', 'measured');

            [nDatasets, symmNoAll, LAll] = SensitivityAnalysis.checkDatasets(structure, ...
                association, measured, invertedLayerNo, 'SensitivityAnalysis:calc');
            if any(symmNoAll ~= symmNoAll(1)) || any(LAll ~= LAll(1))
                error('SensitivityAnalysis:calc:symmetryMismatch', ...
                    ['All jointly-inverted datasets must share the same symmetry class and ', ...
                     'number of independent constants for the inverted layer (it is the same ', ...
                     'crystal): symmClassNo = [%s], nIndepConsts = [%s].'], ...
                    num2str(symmNoAll(:).'), num2str(LAll(:).'));
            end
            symmNo = symmNoAll(1);
            L      = LAll(1);
            %%% Canonical initial guess: dataset 1's stored elConsts for the inverted layer.
            %%% Every other dataset's own stored value for that layer is not read again below
            %%% (only its dirMat/rho/other-layers matter) -- see the class docblock.
            cVec0  = structure(1).mater(invertedLayerNo).elConsts.vector(:);

            includeSet = ModeAssociation.activeRitzClasses(opts.IncludePseudoSAW);

            %%% ---- Per-dataset G_d blocks, then pooled ----
            Gcell     = cell(nDatasets, 1);
            nRitzAll  = zeros(nDatasets, 1);
            nLimAll   = zeros(nDatasets, 1);
            nPhaseAll = zeros(nDatasets, 1);
            hasBulkAll = false(nDatasets, 1);
            for d = 1:nDatasets
                [Gcell{d}, counts_d] = SensitivityAnalysis.datasetSensitivityRows( ...
                    structure(d), association(d), invertedLayerNo, measured(d), ...
                    cVec0, symmNo, L, includeSet, opts);
                nRitzAll(d)   = counts_d.nRitz;
                nLimAll(d)    = counts_d.nLim;
                nPhaseAll(d)  = counts_d.nPhase;
                hasBulkAll(d) = counts_d.hasBulk;
            end

            G = cell2mat(Gcell);
            nData = size(G, 1);

            if isempty(opts.PointWeights)
                w = ones(nData, 1);
            else
                w = opts.PointWeights(:);
                if numel(w) ~= nData
                    error('SensitivityAnalysis:calc:badWeights', ...
                        'PointWeights has %d entries; expected %d (one per data row).', ...
                        numel(w), nData);
                end
            end
            G = diag(w) * G;

            %%% Eigendecomposition via SVD of G, not eig(G'G) -- forming the normal matrix
            %%% squares the condition number and loses precision precisely in the small
            %%% eigenvalues the cutoff decision depends on.
            [~, Sv, V] = svd(G, 'econ');
            T   = V;
            psi = 2 * diag(Sv).^2;

            %%% pin each column's sign so results are reproducible run to run
            for pCol = 1:size(T, 2)
                [~, iMax] = max(abs(T(:, pCol)));
                if T(iMax, pCol) < 0
                    T(:, pCol) = -T(:, pCol);
                end
            end

            normRoots = sqrt(psi / psi(1));
            nOpt      = sum(normRoots >= opts.SensitivityCutoff);

            obj = SensitivityAnalysis();
            obj.nIndepConsts     = L;
            obj.symmClassNo      = symmNo;
            obj.nDatasets        = nDatasets;
            obj.initialVoigtVec  = cVec0;
            obj.initialTransformedVec = T' * cVec0;
            obj.sensitivityMatrix  = G;
            obj.sensitivityEigenvalues = psi;
            obj.transformMatrix   = T;
            obj.normalisedSensitivityRoots = normRoots;
            obj.sensitivityCutoff = opts.SensitivityCutoff;
            obj.nOptimisedConstants = nOpt;
            obj.pointWeights      = w;
            obj.includePseudoSAW  = opts.IncludePseudoSAW;
            obj.nRitzRows         = sum(nRitzAll);
            obj.nLimitingRows     = sum(nLimAll);
            obj.nPhaseRows        = sum(nPhaseAll);
            obj.hasBulkModulusRow = any(hasBulkAll);
        end
    end

    methods (Static, Access = private)
        function [G_d, counts] = datasetSensitivityRows(structure, association, invertedLayerNo, ...
                measured, cVec0, symmNo, L, includeSet, opts)
            %%% Build ONE dataset's sensitivity-matrix block G_d
            %%% >> [G_d, counts] = SensitivityAnalysis.datasetSensitivityRows(structure, ...
            %%%        association, invertedLayerNo, measured, cVec0, symmNo, L, includeSet, opts)
            %%% OUTPUT: counts .. struct with fields nRitz, nLim, nPhase, hasBulk
            %%%
            %%% Ritz Rayleigh-quotient rows, then limiting/phase/bulk central-FD rows. calc()
            %%% vertically stacks the per-dataset blocks before the single pooled SVD.

            matLay = structure.mater(invertedLayerNo);
            rho    = matLay.rho;

            %%% ---- Ritz-mode rows (analytic) ----
            ritzIdxGlobal = find(ismember({association.points.bindingClass}, includeSet));
            nRitz  = numel(ritzIdxGlobal);
            G_ritz = zeros(nRitz, L);

            if nRitz > 0
                precompStruct = CalcCommon.ritzSAWPrecompute(structure);
                N         = precompStruct.nullConstrMat;
                dKdC_cell = precompStruct.dKdC_cell;

                uniqueAngles = unique([association.points(ritzIdxGlobal).angleDeg]);
                useP = ParallelPolicy.shouldUseParallel(opts.UseParallel, numel(uniqueAngles)) ...
                     && ParallelPolicy.ensureThreadPool();

                rowsCell = cell(numel(uniqueAngles), 1);
                idxCell  = cell(numel(uniqueAngles), 1);
                if useP
                    parfor ia = 1:numel(uniqueAngles)
                        [rowsCell{ia}, idxCell{ia}] = SensitivityAnalysis.ritzRowsAtAngle( ...
                            structure, invertedLayerNo, dKdC_cell, N, uniqueAngles(ia), ...
                            association, ritzIdxGlobal, includeSet); %#ok<PFBNS>
                    end
                else
                    for ia = 1:numel(uniqueAngles)
                        [rowsCell{ia}, idxCell{ia}] = SensitivityAnalysis.ritzRowsAtAngle( ...
                            structure, invertedLayerNo, dKdC_cell, N, uniqueAngles(ia), ...
                            association, ritzIdxGlobal, includeSet);
                    end
                end
                for ia = 1:numel(uniqueAngles)
                    G_ritz(idxCell{ia}, :) = rowsCell{ia};
                end
            end

            %%% ---- Limiting-velocity rows (central finite difference) ----
            %%% Anchored LOCAL refinement from each point's own tracked theta*/branch
            %%% (limitingVelocityTracked), at the wider limiting-channel step. NOT a full
            %%% sagittal rescan per perturbation, and NOT the 1e-6*|c| step: both give WRONG
            %%% derivatives, not merely slow ones, because the coarse-grid bracket and
            %%% fminbnd's default TolX of 1e-4 rad make the velocity non-smooth at the
            %%% perturbation scale. G
            %%% feeds T and the nOptimisedConstants cutoff, so an error here can change how
            %%% many transformed constants get optimised on a limiting-heavy dataset.
            limIdxGlobal = find(strcmp({association.points.bindingClass}, 'limiting'));
            nLim  = numel(limIdxGlobal);
            G_lim = zeros(nLim, L);
            for k = 1:nLim
                pt  = association.points(limIdxGlobal(k));
                ctx = struct('symmNo', symmNo, 'matLay', matLay, 'pt', pt);
                [G_lim(k, :), ok] = SensitivityAnalysis.channelDerivative('limiting', cVec0, L, @(c) c, ctx);
                if ~ok
                    error('SensitivityAnalysis:calc:limitingRootLost', ...
                        ['Limiting point %d (%s at %.6g deg): the tracked Christoffel root vanished ', ...
                         'within the finite-difference step, so its sensitivity row is undefined. ', ...
                         'Remove the point or start from another initial guess.'], ...
                        limIdxGlobal(k), pt.waveType, pt.angleDeg);
                end
            end

            %%% ---- Pulse-echo phase-velocity rows (central finite difference) ----
            nPhase  = numel(measured.phaseVelocities);
            G_phase = zeros(nPhase, L);
            for k = 1:nPhase
                ctx = struct('symmNo', symmNo, 'rho', rho, 'pv', measured.phaseVelocities(k));
                G_phase(k, :) = SensitivityAnalysis.channelDerivative('phase', cVec0, L, @(c) c, ctx);
            end

            %%% ---- Bulk-modulus pseudo-datum row (central finite difference) ----
            hasBulk = ~isnan(measured.bulkModulus);
            if hasBulk
                ctx = struct('symmNo', symmNo, 'rho', rho);
                G_bulk = SensitivityAnalysis.channelDerivative('bulk', cVec0, L, @(c) c, ctx);
            else
                G_bulk = zeros(0, L);
            end

            G_d = [G_ritz; G_lim; G_phase; G_bulk];
            counts = struct('nRitz', nRitz, 'nLim', nLim, 'nPhase', nPhase, 'hasBulk', hasBulk);
        end

        function [rows, localIdx] = ritzRowsAtAngle(structure, invertedLayerNo, dKdC_cell, N, ...
                angleDeg, association, ritzIdxGlobal, includeSet)
            %%% Sensitivity rows for every Ritz-bound point at one angle, in one block
            %%% quadratic form per constant (not per mode)
            %%% >> [rows, localIdx] = SensitivityAnalysis.ritzRowsAtAngle(structure, ...
            %%%        invertedLayerNo, dKdC_cell, N, angleDeg, association, ritzIdxGlobal, includeSet)
            %%%
            %%% includeSet must be the SAME set the caller used to build ritzIdxGlobal:
            %%% otherwise a point selected here but excluded by the caller gets
            %%% localIdx = 0 from ismember() below, and the caller's G_ritz(idxCell{ia},:)
            %%% assignment throws MATLAB:badsubscript.
            globalIdxAtAngle = association.pointsAtAngle(angleDeg, includeSet);
            [~, localIdx]    = ismember(globalIdxAtAngle, ritzIdxGlobal);

            %%% dK/dc is non-zero only on the inverted layer's diagonal block, so the quadratic
            %%% forms A'(N'(dK/dc)N)A = U'(dK/dc)_block U need only U = N(rows,:)*A -- no
            %%% projection of dK/dc into the constrained basis at all.
            [dG_blk, rowsLay] = GammaMatrix.layerBlockDerivs(structure, invertedLayerNo, dKdC_cell, angleDeg);
            L      = size(dG_blk, 3);
            nPts   = numel(globalIdxAtAngle);

            A    = zeros(size(N, 2), nPts);
            vels = zeros(nPts, 1);
            for k = 1:nPts
                A(:, k) = association.points(globalIdxAtAngle(k)).anchorEigVec_con;
                vels(k) = association.points(globalIdxAtAngle(k)).calculatedVelocity;
            end

            U    = N(rowsLay, :) * A;   % full-space eigenvectors, inverted layer's rows
            rows = zeros(nPts, L);
            for i = 1:L
                quadForms  = real(sum(conj(U) .* (dG_blk(:,:,i) * U), 1)).';
                rows(:, i) = SensitivityAnalysis.ritzVelocityDerivative(quadForms, vels);
            end
        end

    end

    methods (Static, Access = {?SensitivityAnalysis, ?ElasticityInversion, ?ModeAssociation})
        %%% Shared forward-model evaluators for the three non-Ritz data channels, also
        %%% used by ElasticityInversion's residual computation -- one formula per channel,
        %%% never duplicated.
        %%%
        %%% NEVER resolve 'fastShear'/'slowShear' as max()/min() of shearLimVels: that
        %%% presumes it holds exactly two entries, one per branch, whereas either branch
        %%% can hold several roots or none, and the slow-shear branch can supply the
        %%% maximum -- max() then silently returns a number from the wrong branch.
        %%% limitingRootsAt returns the full labelled root set instead, and the caller
        %%% tracks ONE root by identity.

        function roots = limitingRootsAt(matLay, symmNo, cVec, angleDeg)
            %%% The full LABELLED Christoffel limiting-root set at one angle
            %%% >> roots = SensitivityAnalysis.limitingRootsAt(matLay, symmNo, cVec, angleDeg)
            %%%
            %%% INPUT:  1. matLay   .. Material of the inverted layer, UNROTATED
            %%%         2. symmNo   .. symmetry class number
            %%%         3. cVec     .. independent elastic constants [L x 1], Voigt order
            %%%         4. angleDeg .. propagation angle [deg]
            %%% OUTPUT: roots .. [1 x nRoots] struct array, one entry per ACCEPTED
            %%%                  sagittal-plane root of all three wave types, fields:
            %%%                    iwt   .. 1 = longitudinal, 2 = fast shear, 3 = slow
            %%%                             shear (LimitingVels.findRoots' descending-
            %%%                             eigenvalue convention == waveTypeIndex's)
            %%%                    vel   .. limiting velocity [km/s]
            %%%                    pol   .. [3x1] polarisation at the root
            %%%                    theta .. sagittal root angle theta* [rad]
            %%%
            %%% Returns every labelled root rather than one max/min reduction, so that a
            %%% caller can re-identify the ONE root it is tracking by polarisation/
            %%% velocity/theta/branch. See ModeAssociation.rebindLimiting.
            matPerturbed = matLay.SetElastConsts(ElastConsts(cVec, symmNo));
            matRot = matPerturbed.SetDirMat(matPerturbed.GetRotDirMat(angleDeg));
            lim = CalcCommon.limitingAndBulk(matRot);

            nShear = numel(lim.shearLimVels);
            roots  = struct('iwt', cell(1, nShear + 1), 'vel', [], 'pol', [], 'theta', []);

            roots(1).iwt   = 1;
            roots(1).vel   = lim.longLimVel;
            roots(1).pol   = lim.longLimPol(:);
            roots(1).theta = lim.longLimTheta;
            for q = 1:nShear
                roots(q+1).iwt   = lim.shearLimIwt(q);
                roots(q+1).vel   = lim.shearLimVels(q);
                roots(q+1).pol   = lim.shearLimPols(:, q);
                roots(q+1).theta = lim.shearLimTheta(q);
            end
        end

        function [v, theta, pol, ok] = limitingVelocityNear(matLay, symmNo, cVec, angleDeg, ...
                thetaAnchor, iwt, halfWidth)
            %%% ONE anchored limiting root, re-found locally instead of re-scanned
            %%% >> [v, theta, pol, ok] = SensitivityAnalysis.limitingVelocityNear(matLay, symmNo, cVec, angleDeg, thetaAnchor, iwt, halfWidth)
            %%%
            %%% INPUT:  1-4. as limitingRootsAt
            %%%         5. thetaAnchor .. tracked root's theta* [rad] from the previous
            %%%                          ACCEPTED iterate -- the bracket centre, which
            %%%                          must NOT depend on cVec (see refineRoot)
            %%%         6. iwt         .. tracked root's branch index
            %%%  (opt)  7. halfWidth   .. bracket half-width [rad] (LimitingVels.refineRoot
            %%%                          default, 0.05)
            %%% OUTPUT: v, theta, pol, ok .. as LimitingVels.refineRoot. ok = false means
            %%%                          the anchored root no longer exists; v is then
            %%%                          meaningless and the caller must rescan.
            %%%
            %%% Much cheaper than the full scan, which matters because the limiting
            %%% channel's finite differences dominate an objective evaluation once several
            %%% constants are being optimised. It is ALSO what makes that difference
            %%% correct at all (see LimitingVels.refineRoot).
            if nargin < 7; halfWidth = []; end
            rotElTensor = SensitivityAnalysis.rotatedTensorFor(matLay, symmNo, cVec, angleDeg);
            [v, theta, pol, ~, ok] = LimitingVels.refineRoot(rotElTensor, matLay.rho, ...
                thetaAnchor, iwt, halfWidth);
        end

        function [v, theta, ok] = limitingVelocityTracked(matLay, symmNo, cVec, angleDeg, pt, halfWidth)
            %%% Anchored local refinement WITH the mandatory full-rescan fallback
            %%% >> [v, theta, ok] = SensitivityAnalysis.limitingVelocityTracked(matLay, symmNo, cVec, angleDeg, pt, halfWidth)
            %%%
            %%% INPUT:  1-4. as limitingRootsAt
            %%%         5. pt        .. ONE ModeAssociation.points entry of bindingClass
            %%%                        "limiting" (reads limRootTheta, limBranchIwt,
            %%%                        limAnchorPol, calculatedVelocity)
            %%%  (opt)  6. halfWidth .. bracket half-width [rad]
            %%% OUTPUT: v, theta .. the tracked root's velocity and refined theta*
            %%%         ok       .. false only if BOTH the refinement and the rescan
            %%%                     failed to find the anchored root at all
            %%%
            %%% The single place the "refine, and if the root vanished rescan" rule
            %%% lives: SensitivityAnalysis.calc and ElasticityInversion.costAndJacobian
            %%% both call THIS, so a perturbed finite-difference evaluation can never
            %%% silently return a meaningless value from a bracket whose root has gone
            %%% (the vanishing-root case described in ModeAssociation.rebindLimiting).
            if nargin < 6; halfWidth = []; end
            [v, theta, ~, ok] = SensitivityAnalysis.limitingVelocityNear(matLay, symmNo, ...
                cVec, angleDeg, pt.limRootTheta, pt.limBranchIwt, halfWidth);
            if ok
                return
            end
            roots = SensitivityAnalysis.limitingRootsAt(matLay, symmNo, cVec, angleDeg);
            if isempty(roots)
                v = NaN; theta = pt.limRootTheta; ok = false;
                return
            end
            score = ModeAssociation.limitingRootScore(pt, roots);
            [bestScore, q] = max(score);
            if bestScore <= 0
                v = NaN; theta = pt.limRootTheta; ok = false;
                return
            end
            v     = roots(q).vel;
            theta = roots(q).theta;
            ok    = true;
        end

        function idx = waveTypeIndex(waveType)
            %%% Map a limiting/phase-velocity wave-type label onto the Christoffel
            %%% branch index iwt (descending eigenvalue order, the same convention
            %%% LimitingVels.findRoots and Christoffel.matrix use)
            %%% >> idx = SensitivityAnalysis.waveTypeIndex(waveType)
            switch waveType
                case 'longitudinal'; idx = 1;
                case 'fastShear';    idx = 2;
                case 'slowShear';    idx = 3;
                otherwise
                    error('SensitivityAnalysis:waveTypeIndex:badWaveType', ...
                        '"%s" has no Christoffel branch index (expected longitudinal, fastShear or slowShear).', ...
                        waveType);
            end
        end

        function v = phaseVelocityAt(cVec, symmNo, dirVec, rho, waveType)
            %%% Bulk phase velocity along dirVec, given in the INVERTED LAYER'S CRYSTAL FRAME
            %%% (the frame ElastConsts.tensor's components are expressed in), not the
            %%% dirMat/lab frame the angular-sweep channel uses
            %%% >> v = SensitivityAnalysis.phaseVelocityAt(cVec, symmNo, dirVec, rho, waveType)
            idx = SensitivityAnalysis.waveTypeIndex(waveType);
            elTensor = ElastConsts.Any2Tensor(cVec, symmNo);
            eigVals  = Christoffel.matrix(elTensor, dirVec, rho);
            v = sqrt(eigVals(idx));
        end

        function v = bulkVelocityAt(cVec, symmNo, rho)
            %%% Bulk-modulus pseudo-velocity sqrt(B/rho)
            %%% >> v = SensitivityAnalysis.bulkVelocityAt(cVec, symmNo, rho)
            elConsts = ElastConsts(cVec, symmNo);
            B = ElastConsts.Voigt2Bulk(elConsts.voigt);
            v = sqrt(B / rho);
        end

        function h = limitingStep(cAbs)
            %%% Central-difference step for the LIMITING channel only [GPa]
            %%% >> h = SensitivityAnalysis.limitingStep(cAbs)
            %%%
            %%% h = max(1e-4*|c|, 1e-3), NOT the 1e-6*|c| the other channels use. With
            %%% the bracket held fixed and fminbnd at TolX 1e-12 (LimitingVels.refineRoot)
            %%% the limiting velocity is smooth over a wide range of h, but 1e-6*|c| sits
            %%% at its noisy end. The pulse-echo and bulk-modulus channels are analytic
            %%% functions of the constants and keep their own near-optimal 1e-6*|c| step.
            h = max(1e-4 * abs(cAbs), 1e-3);
        end

        function h = phaseBulkStep(cAbs)
            %%% Central-difference step of the pulse-echo and bulk-modulus channels [GPa]
            %%% >> h = SensitivityAnalysis.phaseBulkStep(cAbs)
            %%%
            %%% h = max(1e-6*|c|, 1e-6): both channels are analytic functions of the
            %%% constants (see limitingStep for why the limiting channel differs).
            h = max(1e-6 * abs(cAbs), 1e-6);
        end

        function [row, ok] = channelDerivative(kind, x, nCols, toVoigt, ctx)
            %%% Central-difference row dv/dx_i of one non-Ritz data point, i = 1..nCols
            %%% >> [row, ok] = SensitivityAnalysis.channelDerivative(kind, x, nCols, toVoigt, ctx)
            %%%
            %%% INPUT:  1. kind    .. 'limiting' | 'phase' | 'bulk'
            %%%         2. x       .. [L x 1] point of the difference: Voigt constants
            %%%                       (SensitivityAnalysis) or transformed constants
            %%%                       (ElasticityInversion)
            %%%         3. nCols   .. differentiate against x(1:nCols)
            %%%         4. toVoigt .. function handle: a perturbed x -> the Voigt vector the
            %%%                       channel is evaluated at (@(c) c, or @(c) T*c)
            %%%         5. ctx     .. struct: symmNo; rho ('phase'/'bulk'); pv, one
            %%%                       MeasuredVelocities.phaseVelocities entry ('phase');
            %%%                       matLay, pt, one ModeAssociation point ('limiting')
            %%% OUTPUT: row .. [1 x nCols], +dv/dx (the sensitivity sign; the inversion's
            %%%                Jacobian row, of measured - calculated, is -row)
            %%%         ok  .. false if a perturbed limiting evaluation lost its root
            %%%                (that entry of row is then NaN)
            %%%
            %%% The one finite-difference loop of SensitivityAnalysis.calc and
            %%% ElasticityInversion.costAndJacobian. The step is limitingStep for the
            %%% limiting channel and phaseBulkStep for the other two.
            row = zeros(1, nCols);
            ok  = true;
            for i = 1:nCols
                switch kind
                    case 'limiting'
                        h = SensitivityAnalysis.limitingStep(x(i));
                    case {'phase', 'bulk'}
                        h = SensitivityAnalysis.phaseBulkStep(x(i));
                    otherwise
                        error('SensitivityAnalysis:channelDerivative:badKind', ...
                            '"%s" is not a data channel (expected limiting, phase or bulk).', kind);
                end
                xP = x; xP(i) = xP(i) + h;
                xM = x; xM(i) = xM(i) - h;
                [vP, okP] = SensitivityAnalysis.channelValue(kind, toVoigt(xP), ctx);
                [vM, okM] = SensitivityAnalysis.channelValue(kind, toVoigt(xM), ctx);
                row(i) = (vP - vM) / (2*h);
                ok = ok && okP && okM;
            end
        end

        function [v, ok] = channelValue(kind, cVec, ctx)
            %%% One non-Ritz data point's calculated velocity at the Voigt vector cVec
            %%% >> [v, ok] = SensitivityAnalysis.channelValue(kind, cVec, ctx)
            %%%
            %%% INPUT: kind, ctx .. as channelDerivative
            %%% OUTPUT: v .. [km/s]; ok .. false only for a limiting point whose root was
            %%%              lost (limitingVelocityTracked), v is then NaN
            ok = true;
            switch kind
                case 'limiting'
                    [v, ~, ok] = SensitivityAnalysis.limitingVelocityTracked(ctx.matLay, ...
                        ctx.symmNo, cVec, ctx.pt.angleDeg, ctx.pt);
                case 'phase'
                    v = SensitivityAnalysis.phaseVelocityAt(cVec, ctx.symmNo, ...
                        ctx.pv.propagationDir, ctx.rho, ctx.pv.waveType);
                case 'bulk'
                    v = SensitivityAnalysis.bulkVelocityAt(cVec, ctx.symmNo, ctx.rho);
                otherwise
                    error('SensitivityAnalysis:channelValue:badKind', ...
                        '"%s" is not a data channel (expected limiting, phase or bulk).', kind);
            end
        end

        function d = ritzVelocityDerivative(quadForms, vels)
            %%% Ritz-channel dv/dc from the Rayleigh-quotient numerators
            %%% >> d = SensitivityAnalysis.ritzVelocityDerivative(quadForms, vels)
            %%%
            %%% INPUT:  1. quadForms .. Re[alpha' (dK/dc) alpha], one per mode (real)
            %%%         2. vels      .. the modes' phase velocities [km/s], same size
            %%% OUTPUT: d .. quadForms ./ (8 pi^2 v): with v = sqrt(lambda)/(2 pi) and
            %%%              dlambda/dc = alpha' (dK/dc) alpha for an M-orthonormal alpha
            d = quadForms ./ (8*pi^2 .* vels);
        end

        function [nDatasets, symmNoAll, LAll] = checkDatasets(structure, association, ...
                measured, invertedLayerNo, idPrefix)
            %%% Dataset-array checks shared by SensitivityAnalysis.calc and the
            %%% ElasticityInversion constructor
            %%% >> [nDatasets, symmNoAll, LAll] = SensitivityAnalysis.checkDatasets(structure, association, measured, invertedLayerNo, idPrefix)
            %%%
            %%% INPUT:  1-3. structure, association, measured .. the dataset arrays
            %%%         4. invertedLayerNo .. shared inverted-layer index
            %%%         5. idPrefix        .. 'Class:method' of the caller's error ids
            %%% OUTPUT: nDatasets .. numel(structure)
            %%%         symmNoAll, LAll .. [nDatasets x 1] symmClassNo / nIndepConsts of
            %%%                            every dataset's inverted layer (the caller
            %%%                            decides what they must agree with)
            %%%
            %%% Errors <idPrefix>:datasetCountMismatch (arrays of unequal length) and
            %%% <idPrefix>:badLayer (invertedLayerNo out of range for any dataset).
            nDatasets = numel(structure);
            if numel(association) ~= nDatasets || numel(measured) ~= nDatasets
                error([idPrefix ':datasetCountMismatch'], ...
                    ['structure (%d), association (%d) and measured (%d) must all supply the ', ...
                     'same number of datasets.'], nDatasets, numel(association), numel(measured));
            end
            if invertedLayerNo < 1 || any([structure.nLay] < invertedLayerNo)
                error([idPrefix ':badLayer'], ...
                    'invertedLayerNo %d is out of range for at least one dataset (nLay = [%s]).', ...
                    invertedLayerNo, num2str([structure.nLay]));
            end
            symmNoAll = zeros(nDatasets, 1);
            LAll      = zeros(nDatasets, 1);
            for d = 1:nDatasets
                ecLay_d      = structure(d).mater(invertedLayerNo).elConsts;
                symmNoAll(d) = ecLay_d.symmClassNo;
                LAll(d)      = ecLay_d.nIndepConsts;
            end
        end
    end

    methods (Static, Access = private)
        function rotElTensor = rotatedTensorFor(matLay, symmNo, cVec, angleDeg)
            %%% The rotated elastic tensor CalcCommon.limitingAndBulk would build
            %%% internally, for callers that need LimitingVels directly instead
            %%% >> rotElTensor = SensitivityAnalysis.rotatedTensorFor(matLay, symmNo, cVec, angleDeg)
            %%%
            %%% Must stay bit-identical to limitingAndBulk's own chain (propagation
            %%% rotation of the dirMat, then GetStandardElastTensor(0) -- 0 because the
            %%% dirMat already carries the angle; see docs/reference/architecture/forwardCalcs.md,
            %%% "Rotation convention: pass 0 to solvers after construction"), otherwise a refined root and a rescanned
            %%% one would not agree.
            matPerturbed = matLay.SetElastConsts(ElastConsts(cVec, symmNo));
            matRot  = matPerturbed.SetDirMat(matPerturbed.GetRotDirMat(angleDeg));
            rotElTensor = matRot.GetStandardElastTensor(0);
        end
    end
end
