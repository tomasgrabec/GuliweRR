classdef RitzSolverSAW < RitzSolver
% RitzSolverSAW — Ritz solver for the half-space (surface acoustic waves).
%
%%% GuliweRR v6.0

    properties (SetAccess = protected)
        diagnostics      struct = struct('freqDisp', [], 'oopSurfDispl', [], ...
                                        'depthInsensitivity', [])
                                        % per-mode diagnostics (modeDiagnostics): freqDisp,
                                        % oopSurfDispl, depthInsensitivity
    end

    methods
        function obj = RitzSolverSAW(structure, angleDeg)
        %%% Initialise RitzSolverSAW
        %%% >> obj = RitzSolverSAW(structure, angleDeg)
            %%% RitzSolver's constructor guards each argument with ~isempty(),
            %%% so [] is equivalent to omitting the argument.
            if nargin < 1; structure = []; end
            if nargin < 2; angleDeg  = []; end
            obj@RitzSolver(structure, angleDeg);
            obj.diagnostics = struct('freqDisp', [], 'oopSurfDispl', [], 'depthInsensitivity', []);
        end

        function obj = calc(obj, withDiagnostics, precompStruct)
        %%% Calculate SAW modes for half-space domain
        %%% >> obj = calc(obj, withDiagnostics, precompStruct)
        %%%
        %%% INPUT:
        %%%  (opt) 1. withDiagnostics .. logical: also compute the per-mode diagnostics
        %%%                              (modeDiagnostics -> obj.diagnostics; default true)
        %%%  (opt) 2. precompStruct   .. optional pre-computed quantities:
        %%%            .nullConstrMat   — constraint null space
        %%%            .dKdC_cell       — from GammaMatrix.derivsByDepth
        %%%            .d2KdCdD_cell    — from GammaMatrix.derivsByDepth

            if nargin < 2 || isempty(withDiagnostics); withDiagnostics = true; end
            RitzSolver.checkFlag(withDiagnostics, 'RitzSolverSAW:calc:badFlag', 'withDiagnostics');
            if nargin < 3 || isempty(precompStruct);   precompStruct   = struct(); end

            %%% constraint matrix
            obj.nullConstrMat = RitzSolver.precompField(precompStruct, 'nullConstrMat', ...
                @() RitzSolver.buildConstraintMatrix(obj.structure.calcParams, 'saw'));

            %%% dK/dc and d²K/(dc dD) cells -- GammaMatrix.calc builds both when
            %%% both are absent and rejects a half-supplied pair
            dKdC_cell    = RitzSolver.precompField(precompStruct, 'dKdC_cell');
            d2KdCdD_cell = RitzSolver.precompField(precompStruct, 'd2KdCdD_cell');

            %%% build GammaMatrix K and M -- a local, not stored (see projectToConstrained)
            gm = GammaMatrix();
            gm = gm.calc(obj.structure, obj.angle, dKdC_cell, d2KdCdD_cell);

            %%% project into constrained basis
            [K_con, M_con, gauge] = obj.projectToConstrained(gm);

            %%% solve and fill base properties. Layered (nLay>1) only: depthInsensitivity
            %%% from the substrate block of dK/dD alone. See docs/reference/architecture/ritz.md.
            obj = obj.solveEigenvalueProblem(K_con, M_con, gauge);
            obj = obj.calcGroupVelocity(gm, obj.structure.nLay > 1);

            %%% per-mode diagnostics
            if withDiagnostics && obj.nCompModes > 0
                obj = obj.modeDiagnostics();
            end
        end

        function obj = modeDiagnostics(obj)
        %%% Per-mode diagnostics read by the classification and the viewers
        %%% >> obj = modeDiagnostics(obj)
        %%%
        %%% Fills obj.diagnostics with, per computed mode:
        %%%  - freqDisp           .. exp(-10 ((v - vg)/v)^2): 1 for a non-dispersive mode
        %%%                          (a true half-space surface wave), -> 0 for a
        %%%                          dispersive one;
        %%%  - oopSurfDispl       .. |out-of-plane displacement| at the top surface;
        %%%  - depthInsensitivity .. ritz.depthInsensitivity minus its minimum: large for
        %%%                          modes whose velocity barely changes with the
        %%%                          artificial boundary depth (substrate-only for a
        %%%                          layered structure), i.e. modes confined near the
        %%%                          surface.
        %%% They are diagnostics, not a classification: existence is decided by
        %%% ModeClassification (see docs/reference/architecture/ritz.md).
            nM = obj.nCompModes;

            freqDisp = abs(exp(-10 .* ((obj.vels(1:nM) - obj.velsGroup(1:nM)) ./ max(obj.vels(1:nM), eps)).^2));

            %%% only the top-surface u3 row of the displacement profile is needed
            oopSurfDispl = obj.topSurfaceU3(nM);   % x3 = surface normal

            depthInsens = obj.depthInsensitivity(1:nM) - min(obj.depthInsensitivity(1:nM));

            obj.diagnostics.freqDisp           = freqDisp;
            obj.diagnostics.oopSurfDispl       = oopSurfDispl;
            obj.diagnostics.depthInsensitivity = depthInsens;
        end
    end

    methods (Static)
        function r = ritzSensitivities(structure, angleDeg, precompStruct, nModes, vMax)
        %%% EXPERIMENTAL. Exact first derivatives of the lowest Ritz modes of a layered stack
        %%% >> r = RitzSolverSAW.ritzSensitivities(structure, angleDeg, precompStruct, nModes)
        %%% >> r = RitzSolverSAW.ritzSensitivities(structure, angleDeg, precompStruct, nModes, vMax)
        %%%
        %%% Mirrors the Julia function ritzSensitivities (RitzSensitivity.jl): same outputs,
        %%% field names and formulae.
        %%%
        %%% INPUT:  1. structure     .. layered Structure (nLay > 1)
        %%%         2. angleDeg      .. propagation angle [deg]
        %%%  (opt)  3. precompStruct .. RitzSolverSAW precompute (CalcCommon.ritzSAWPrecompute,
        %%%                             e.g. from a CalcCommon.ritzWavelengthSweep at this
        %%%                             wavelength); missing fields are built ([] = none)
        %%%         4. nModes        .. number of lowest modes, >= 1 (keep it small: tens)
        %%%  (opt)  5. vMax          .. keep only modes slower than this [km/s] (default Inf;
        %%%                             e.g. the substrate's slowest shear limit, for bound
        %%%                             modes)
        %%% OUTPUT: r .. struct, nB <= nModes modes (rows), ascending velocity:
        %%%   .vels        [nB x 1] phase velocities; .eigVals [nB x 1] omega^2 (normalised)
        %%%   .eigVecs     [strBVL x nB] full-basis eigenvectors N * eigVecs_con, a^H M a = 1,
        %%%                original (not real-gauged) basis, as for a RitzSolverSAW solve
        %%%   .dlnv_dlnC   {nLay x 1}, each [nB x nIndepConsts]: dln v / dln c_k for the
        %%%                layer's ElastConsts.vector entries
        %%%   .dlnv_dlnRho [nB x nLay]
        %%%   .dlnv_dlnD   [nB x nLay] (layer 1 = substrate: its depth is the truncation
        %%%                depth, so column 1 measures truncation sensitivity, not physics)
        %%%   .dlnv_dTheta [nB x 1] per degree
        %%%   .vgOverVp    [nB x 1] = 1 + sum over films of dlnv_dlnD (= dln f / dln k)
        %%%   .kinShare    [nB x nLay] kinetic-energy share of each layer
        %%%
        %%% One Ritz solve, no finite differences except dc/dtheta. Hellmann-Feynman for
        %%% K a = lambda M a with a^H M a = 1:  d lambda/dp = a^H (dK/dp - lambda dM/dp) a,
        %%% v = sqrt(lambda)/(2 pi), so dln v/dln p = p (d lambda/dp) / (2 lambda). Per parameter:
        %%%  - elastic constants of layer j: K_j = normFac_j sum_i c_red,i dK/dc_i, with
        %%%    c_red = ElastRotReductMat(R(theta), class) * vector, so
        %%%    d lambda/d vector_k = sum_i R(i,k) q_i, q_i = normFac_j a_j^H (dK/dc_i) a_j;
        %%%  - layer depth D_j: the block (K_j, M_j) is normFac_j (K_j(D_j), I) with
        %%%    normFac_j ~ D_j, so d lambda/dln D_j = a_j^H (K_j + D_j dK_j/dD) a_j
        %%%    - lambda a_j^H M_j a_j. The overall normalisation by D_1 rho_1 scales K and M
        %%%    alike and drops out, so this holds for the substrate too;
        %%%  - density rho_j: normFac_j K_j does not depend on rho_j, M_j ~ rho_j, so
        %%%    dln v/dln rho_j = -1/2 (kinetic share of layer j);
        %%%  - angle: dc_red/dtheta by a +-0.01 deg central difference of the constants alone,
        %%%    as angleSlope (no extra solve);
        %%%  - |k| (wavelength): ApplyWavelength scales the films' D_j with k and keeps the
        %%%    substrate's, so dln v/dln k = sum over films of dlnv_dlnD.
        %%% Checks that hold to round-off: sum of dlnv_dlnC = 1/2 and of dlnv_dlnRho = -1/2
        %%% per mode (v^2 is homogeneous of degree 1 in the constants, -1 in the densities).
        %%% Exact for non-degenerate modes; at a crossing the pair's derivatives are those of
        %%% whichever basis the solver returned.
        %%%
        %%% Deliberate divergence from Julia: MATLAB's eig has no partial-range generalised
        %%% solve, so ALL eigenpairs are computed (RitzSolver.solveConstrainedEig, Cholesky
        %%% path, as every Ritz solve) and the lowest nModes kept. On stiff thin-film stacks
        %%% the lowest eigenvalues are then accurate to ~1e-8 relative (dsygvd-like), against
        %%% ~1e-11 for Julia's lowest-mode bisection.
            if nargin < 3; precompStruct = []; end
            if nargin < 4 || isempty(nModes)
                error('RitzSolverSAW:ritzSensitivities:badNModes', 'nModes must be given (>= 1).');
            end
            if nargin < 5 || isempty(vMax); vMax = Inf; end
            if ~(isa(structure, 'Structure') && isscalar(structure))
                error('RitzSolverSAW:ritzSensitivities:badStructure', 'structure must be one Structure object.');
            end
            if ~(isnumeric(angleDeg) && isscalar(angleDeg) && isreal(angleDeg) && isfinite(angleDeg))
                error('RitzSolverSAW:ritzSensitivities:badAngle', 'angleDeg must be a real finite scalar.');
            end
            if ~(isnumeric(nModes) && isscalar(nModes) && nModes == round(nModes) && nModes >= 1)
                error('RitzSolverSAW:ritzSensitivities:badNModes', 'nModes must be an integer >= 1.');
            end
            if ~(isnumeric(vMax) && isscalar(vMax) && isreal(vMax) && ~isnan(vMax))
                error('RitzSolverSAW:ritzSensitivities:badVMax', 'vMax must be a real scalar.');
            end
            nL = structure.nLay;
            if nL < 2
                error('RitzSolverSAW:ritzSensitivities:notLayered', ...
                    'A layered structure (nLay > 1) is expected.');
            end
            if isempty(precompStruct); precompStruct = struct(); end

            %%% precompute: fill-if-missing, as RitzSolverSAW.calc
            N = RitzSolver.precompField(precompStruct, 'nullConstrMat', ...
                @() RitzSolver.buildConstraintMatrix(structure.calcParams, 'saw'));
            dKdC    = RitzSolver.precompField(precompStruct, 'dKdC_cell');
            d2KdCdD = RitzSolver.precompField(precompStruct, 'd2KdCdD_cell');
            if isempty(dKdC) || isempty(d2KdCdD)
                [dKdC, d2KdCdD] = GammaMatrix.derivsByDepth(structure);
            end

            cps  = structure.calcParams;
            bVLs = [cps.basisVecL];
            row0 = [0, cumsum(bVLs)];
            nfs  = zeros(1, nL);
            for j = 1:nL; nfs(j) = GammaMatrix.LayerNormFac(structure, j); end

            %%% K and dK/dD layer blocks only (the dense block-diagonal GammaMatrix is mostly
            %%% zeros); each block is the same product as in GammaMatrix.calc
            redC = GammaMatrix.RedElastCellFor(structure, angleDeg);
            Kb   = cell(1, nL);
            KDb  = cell(1, nL);
            for j = 1:nL
                bVL    = bVLs(j);
                Kb{j}  = complex(reshape(reshape(dKdC{j},    [], 15) * redC{j}(:), bVL, bVL)) .* nfs(j);
                KDb{j} = complex(reshape(reshape(d2KdCdD{j}, [], 15) * redC{j}(:), bVL, bVL)) .* nfs(j);
            end

            %%% N.' K N and N.' M N summed over the layers' row blocks of N, in the real gauge
            %%% when every block passes it (as RitzSolver.realGauge on the full K: the zero
            %%% off-diagonal blocks add nothing to either norm)
            d   = RitzSolver.gaugeVector(cps);
            Kg  = cell(1, nL);
            im2 = 0; all2 = 0;
            for j = 1:nL
                dj    = d(row0(j) + (1:bVLs(j)));
                Kg{j} = conj(dj) .* Kb{j} .* dj.';
                im2   = im2  + norm(imag(Kg{j}), 'fro')^2;
                all2  = all2 + norm(Kg{j}, 'fro')^2;
            end
            useGauge = sqrt(im2) <= RitzSolver.REAL_GAUGE_TOL * sqrt(all2);
            m     = size(N, 2);
            K_con = zeros(m);
            M_con = zeros(m);
            for j = 1:nL
                Nb = N(row0(j) + (1:bVLs(j)), :);
                if useGauge
                    K_con = K_con + Nb.' * (real(Kg{j}) * Nb);
                else
                    K_con = K_con + Nb.' * (Kb{j} * Nb);   % plain .', never ' (docs/pitfalls.md NUM1)
                end
                M_con = M_con + nfs(j) .* (Nb.' * Nb);
            end

            [lamAll, V] = RitzSolver.solveConstrainedEig(K_con, M_con);
            nk   = min(nModes, numel(lamAll));
            vAll = abs(sqrt(lamAll(1:nk))) / (2*pi);
            nB   = min(nk, nnz(vAll < vMax));
            V    = V(:, 1:nB);
            if useGauge
                %%% back to the original basis, as RitzSolver.solveEigenvalueProblem
                V = N.' * (d .* (N * V));
            end
            U   = N * V;                       % full-basis eigenvectors
            lam = lamAll(1:nB);

            %%% per-layer quotients: kinetic share, K and dK/dD, and q_i of every dK/dc_i
            %%% (a^H D a = vec(conj(a) a.').' vec(D) = kron(a, conj(a)).' vec(D))
            kin = zeros(nB, nL); qK = zeros(nB, nL); qD = zeros(nB, nL);
            Q   = cell(1, nL);                 % [15 x nB] per layer
            for j = 1:nL
                Uj = U(row0(j) + (1:bVLs(j)), :);
                kin(:, j) = nfs(j) .* real(sum(conj(Uj) .* Uj, 1)).';
                qK(:, j)  = real(sum(conj(Uj) .* (Kb{j}  * Uj), 1)).';
                qD(:, j)  = real(sum(conj(Uj) .* (KDb{j} * Uj), 1)).';
                P = complex(zeros(bVLs(j)^2, nB));
                for mm = 1:nB
                    P(:, mm) = kron(Uj(:, mm), conj(Uj(:, mm)));
                end
                Q{j} = nfs(j) .* real(reshape(dKdC{j}, [], 15).' * P);
            end
            mtot = sum(kin, 2);
            den  = 2 .* lam .* mtot;

            relDepths   = [cps.relDepth];
            dlnv_dlnD   = (qK + relDepths .* qD - lam .* kin) ./ den;
            dlnv_dlnRho = -0.5 .* kin ./ mtot;
            dlnv_dlnC   = cell(nL, 1);
            dlnv_dTheta = zeros(nB, 1);
            h  = 0.01;                         % deg, as angleSlope
            cP = GammaMatrix.RedElastCellFor(structure, angleDeg + h);
            cM = GammaMatrix.RedElastCellFor(structure, angleDeg - h);
            for j = 1:nL
                ec = structure.mater(j).elConsts;
                Rm = ElastConsts.ElastRotReductMat(GammaMatrix.LayerRotMat(structure, j, angleDeg), ec.symmClassNo);
                dl = Rm.' * Q{j};              % [nIndep x nB]: d lambda / d vector_k
                dlnv_dlnC{j} = (ec.vector(:) .* dl ./ den.').';
                dc = (cP{j}(:) - cM{j}(:)) / (2*h);
                dlnv_dTheta = dlnv_dTheta + (dc.' * Q{j}).' ./ den;
            end
            vgOverVp = 1 + sum(dlnv_dlnD(:, 2:nL), 2);

            r = struct('vels', vAll(1:nB), 'eigVals', lam, 'eigVecs', U, ...
                'dlnv_dlnC', {dlnv_dlnC}, 'dlnv_dlnRho', dlnv_dlnRho, 'dlnv_dlnD', dlnv_dlnD, ...
                'dlnv_dTheta', dlnv_dTheta, 'vgOverVp', vgOverVp, 'kinShare', kin ./ mtot);
        end
    end
end
