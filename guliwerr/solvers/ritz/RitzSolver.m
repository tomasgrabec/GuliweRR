classdef RitzSolver
% RitzSolver — base class for Ritz-Rayleigh forward solvers.
%   Shared eigenvalue solver, group velocity, constraint matrix, displacement.
%   Subclasses: RitzSolverSAW (half-space), RitzSolverLamb (plate).
%
%%% GuliweRR v6.0

    properties (SetAccess = protected)
        structure       Structure   = Structure()     % computational domain
        angle(1,1)      double      = 0               % propagation angle [°] w.r.t. dirMat row 2

        nCompModes(1,1) double      = 0               % number of stored modes

        freqs           double                         % mode frequencies = sqrt(eigVals)/(2π)
        vels            double                         % phase velocities = freqs (relLength=1)
        velsGroup       double                         % group velocities vg = vp + Σ_i d_i·∂vp/∂d_i

        eigVecs_con     double                         % eigenvectors in constrained basis
        eigVals         double                         % eigenvalues = ω² (normalised)

        nullConstrMat   double                         % N — null space of constraint matrix
        basisVecL_con(1,1) double   = 0               % columns of N

        dVel_dDepth     double                         % dvp/dD per mode (Rayleigh quotient deriv.; nLay>1: all layers +1 at once)
        depthInsensitivity double                        % log10(1/|dvp/dD|) — large for surface modes
    end

    properties (Constant)
        REAL_GAUGE_TOL = 1e-13   % realGauge: ||Im(D'KD)||_F / ||D'KD||_F, round-off level only
    end

    methods
        function obj = RitzSolver(structure, angleDeg)
        %%% Initialise RitzSolver base object
        %%% >> obj = RitzSolver(structure, angleDeg)
            if nargin > 0 && ~isempty(structure)
                if ~(isa(structure, 'Structure') && isscalar(structure))
                    error('RitzSolver:constructor:badStructure', 'structure must be one Structure object.');
                end
                obj.structure = structure;
            end
            if nargin > 1 && ~isempty(angleDeg)
                if ~(isnumeric(angleDeg) && isscalar(angleDeg) && isreal(angleDeg) && isfinite(angleDeg))
                    error('RitzSolver:constructor:badAngle', 'angleDeg must be a real finite scalar.');
                end
                obj.angle = angleDeg;
            end
        end

        function displDir = calcDisplacement(obj, lastModeNo, dscrtSize)
        %%% Evaluate the depth-profile displacement for the first lastModeNo modes
        %%% >> displDir = calcDisplacement(obj, lastModeNo, dscrtSize)
        %%%
        %%% OUTPUT: displDir [depthPoints*nLay x lastModeNo x 3] -- COMPLEX, with a
        %%%         fixed phase per mode: the largest in-plane (u1/u2) Legendre
        %%%         coefficient of the mode is made real and positive (the largest
        %%%         coefficient overall for a mode with no in-plane part).
        %%%
        %%% An eigenvector is defined only up to a phase factor e^(i phi), and the
        %%% eigensolver's choice is arbitrary. Fixing it here makes the output
        %%% reproducible (Vis.calcDispField2D takes real() of it, which an arbitrary
        %%% phase would shift laterally). When every layer is invariant under a 2-fold
        %%% rotation about the surface normal (isotropic, cubic (001)/(110), ...), the
        %%% in-plane components then come out real and SIGNED and u3 as i x real
        %%% (ritzHermitianAndRealGauge.md, Sec. 4); otherwise u1, u2, u3 carry
        %%% genuine relative phases. Callers needing magnitude apply abs().
            if nargin < 2 || isempty(lastModeNo); lastModeNo = obj.nCompModes; end
            if nargin < 3 || isempty(dscrtSize);  dscrtSize  = 0.01;           end

            displSize = round(2/dscrtSize) + 1;
            displDir  = zeros(displSize * obj.structure.nLay, lastModeNo, 3);
            eigVecs   = obj.nullConstrMat * obj.eigVecs_con(:, 1:lastModeNo);
            eigVecs   = eigVecs .* RitzSolver.phaseFix(eigVecs, obj.structure.calcParams);
            lLeg      = zeros(obj.structure.nLay + 1, 1);

            for ilay = 1:obj.structure.nLay
                cp = obj.structure.calcParams(ilay);
                if dscrtSize == 0.01
                    pLeg = cp.polyLeg;
                else
                    pLeg = CalcParams.BasisDecomposition2D(cp.legOrder, dscrtSize);
                end
                lLeg(ilay+1) = lLeg(ilay) + size(pLeg, 2);

                bVL_prev = sum([obj.structure.calcParams(1:ilay-1).basisVecL]);
                bVL      = cp.basisVecL;
                legO     = cp.legOrder;

                for iDir = 1:3
                    startD = bVL_prev + (iDir-1)*legO + 1;
                    endD   = bVL_prev +  iDir   *legO;
                    displDir(lLeg(ilay)+1 : lLeg(ilay+1), :, iDir) = ...
                        pLeg.' * eigVecs(startD:endD, 1:lastModeNo);
                end
            end
        end

        function slope = angleSlope(obj, modeNos, dKdC_cell)
        %%% Angular slope dv/dtheta of selected modes (Hellmann-Feynman, no extra solve)
        %%% >> slope = angleSlope(obj, modeNos)
        %%% >> slope = angleSlope(obj, modeNos, dKdC_cell)
        %%%
        %%% INPUT:  1. modeNos   .. mode indices (into vels)
        %%%  (opt)  2. dKdC_cell .. GammaMatrix.derivsByDepth cell, as the solve used
        %%%                         (precompStruct.dKdC_cell); built if omitted
        %%% OUTPUT: slope .. [numel(modeNos) x 1] dv/dtheta [km/s per deg] at this
        %%%                  angle, theta = the in-plane propagation angle
        %%%
        %%% K is linear in each layer's 15 rotated reduced constants c (GammaMatrix.calc:
        %%% K = sum_c c dK/dc), so dK/dtheta = sum_c (dc/dtheta) dK/dc exactly; only
        %%% dc/dtheta is a central difference (+-0.01 deg, on the constants alone, no
        %%% solve). The eigenvalue derivative is then the Rayleigh-quotient one
        %%% (Hellmann-Feynman), as in calcGroupVelocity:
        %%%     dv/dtheta = alpha^H (dK_con/dtheta) alpha / (8 pi^2 v)
        %%% Exact for a non-degenerate mode; at a crossing of two modes it is the
        %%% slope of whichever combination the solver returned.
        %%% Needs the full eigenvector set (before retentionLevel pruning).
            if size(obj.eigVecs_con, 2) < obj.nCompModes
                error('RitzSolver:angleSlope:noEigVecs', ...
                    'The full Ritz eigenvector set is needed (compute before retentionLevel pruning).');
            end
            if nargin < 3 || isempty(dKdC_cell)
                dKdC_cell = GammaMatrix.derivsByDepth(obj.structure);
            end
            modeNos = modeNos(:);
            if isempty(modeNos); slope = zeros(0, 1); return; end
            h  = 0.01;   % deg
            cP = GammaMatrix.RedElastCellFor(obj.structure, obj.angle + h);
            cM = GammaMatrix.RedElastCellFor(obj.structure, obj.angle - h);
            %%% dK/dtheta is block-diagonal in layer order, so alpha^H (N.' dK N) alpha is
            %%% the sum over layers of U_i^H dK_i U_i with U = N*alpha (N real): no
            %%% projection of dK into the constrained basis, no blkdiag
            U     = obj.nullConstrMat * obj.eigVecs_con(:, modeNos);
            acc   = zeros(1, numel(modeNos));
            row0  = 0;
            for ilay = 1:obj.structure.nLay
                bVL = obj.structure.calcParams(ilay).basisVecL;
                dc  = (cP{ilay} - cM{ilay}) / (2 * h);
                dKl = reshape(reshape(dKdC_cell{ilay}, [], 15) * dc(:), bVL, bVL) ...
                      .* GammaMatrix.LayerNormFac(obj.structure, ilay);
                Ub  = U(row0 + (1:bVL), :);
                acc = acc + real(sum(conj(Ub) .* (dKl * Ub), 1));
                row0 = row0 + bVL;
            end
            v     = obj.vels(modeNos);
            slope = acc.' ./ v(:) / (8 * pi^2);
        end

        function obj = pruneEigVecs(obj, keepCols)
        %%% Retention pruning — keep only selected eigenvector columns
        %%% >> obj = pruneEigVecs(obj, keepCols)
        %%%
        %%% INPUT: keepCols .. column indices of eigVecs_con to retain, or
        %%%                    [] to clear eigVecs_con entirely
        %%% OUTPUT: obj with eigVecs_con trimmed/cleared; vels, eigVals and the
        %%%         per-mode diagnostics are untouched — they are cheap
        %%%         scalars, not the large [strBVL x nModes] eigenvector matrix
        %%%
        %%% eigVecs_con is SetAccess=protected, so this is the only route in for
        %%% an external retention policy (e.g. CalcSingleHS's).
            if nargin < 2 || isempty(keepCols)
                obj.eigVecs_con = [];
            else
                validateattributes(keepCols, {'double'}, {'vector', 'integer', 'positive'}, ...
                    'RitzSolver.pruneEigVecs', 'keepCols');
                obj.eigVecs_con = obj.eigVecs_con(:, keepCols);
            end
        end
    end

    methods (Access = protected)
        function [K_con, M_con, gauge] = projectToConstrained(obj, gm)
        %%% Project the assembled K and M into the constrained basis N
        %%% >> [K_con, M_con, gauge] = projectToConstrained(obj, gm)
        %%%
        %%% INPUT: gm .. this solve's GammaMatrix (K, dK/dD, M). It is a local of
        %%%              calc(), not a property: nothing reads it after the solve, and
        %%%              keeping it cost ~3 MB per angle (strBVL 240), more than the
        %%%              full eigenvector set. Rebuild with GammaMatrix.calc if needed.
        %%%
        %%% N'KN with the PLAIN transpose (.'), never the Hermitian ' —
        %%% gammaMatAll is complex, and ' would silently conjugate it
        %%% (docs/pitfalls.md NUM1).
        %%%
        %%% K_con is projected in the real gauge when the structure allows it
        %%% (RitzSolver.realGauge): K_con is then REAL and gauge is the vector d to
        %%% pass on to solveEigenvalueProblem; otherwise gauge = []. M_con stays in the
        %%% original basis (M is unchanged by the gauge).
        %%%
        %%% M_con = [] for a single layer: there normMatAll is exactly the identity
        %%% (normFac = 1) and N has orthonormal columns, so N.'*M*N is the identity to
        %%% round-off and solveConstrainedEig would only test it and take the standard
        %%% path anyway. dK/dD is not projected at all: calcGroupVelocity evaluates its
        %%% Rayleigh quotients on the full-space eigenvectors N*alpha instead.
            N        = obj.nullConstrMat;
            [K, gauge] = RitzSolver.realGauge(gm.gammaMatAll, obj.structure.calcParams);
            K_con    = N.' * K * N;
            if obj.structure.nLay == 1
                M_con = [];
            else
                M_con = N.' * gm.normMatAll * N;
            end
        end

        function obj = solveEigenvalueProblem(obj, K_con, M_con, gauge)
        %%% Solve constrained eigenvalue problem K_con α = ω² M_con α
        %%% >> obj = solveEigenvalueProblem(obj, K_con, M_con)
        %%% >> obj = solveEigenvalueProblem(obj, K_con, M_con, gauge)
        %%% Fills: eigVals, eigVecs_con, freqs, vels, nCompModes, basisVecL_con
        %%%
        %%% (opt) 3. gauge .. d from RitzSolver.realGauge when K_con is the REAL gauged
        %%%                   matrix; the eigenvectors are mapped back to the original
        %%%                   basis, α = N.'(d .* (N β)), so eigVecs_con keeps its meaning
        %%%                   (complex, M-orthonormal), now with u1, u2 real and u3
        %%%                   imaginary. [] (default): K_con is in the original basis.
            if nargin < 4; gauge = []; end

            [eigValSorted, V_sorted] = RitzSolver.solveConstrainedEig(K_con, M_con);

            obj.nCompModes    = min(obj.structure.nMaxModes, length(eigValSorted));
            obj.eigVals       = eigValSorted(1:obj.nCompModes);
            obj.eigVecs_con   = V_sorted(:, 1:obj.nCompModes);
            if ~isempty(gauge)
                %%% range(N) is invariant under D = diag(d) (each constraint acts on one
                %%% displacement component), so N.' recovers the coordinates exactly
                N = obj.nullConstrMat;
                obj.eigVecs_con = N.' * (gauge .* (N * obj.eigVecs_con));
            end
            obj.freqs         = abs(sqrt(obj.eigVals)) / (2*pi);   % ω/(2π)
            obj.vels          = obj.freqs;                           % v = freq (relLength=1, |k|=1)
            obj.basisVecL_con = size(obj.eigVecs_con, 1);
        end

        function obj = calcGroupVelocity(obj, gm, substrateOnlyInsens)
        %%% Compute group velocity from Rayleigh-quotient depth derivative
        %%% >> obj = calcGroupVelocity(obj, gm)
        %%% >> obj = calcGroupVelocity(obj, gm, substrateOnlyInsens)
        %%% vg = vp + Σ_i d_i ∂vp/∂d_i  (= vp + D·dvp/dD for a single layer)
        %%% Fills: dVel_dDepth, velsGroup, depthInsensitivity
        %%%
        %%% INPUT: 1. gm .. this solve's GammaMatrix (its dGammaDMatAll is read)
        %%% (opt) 3. substrateOnlyInsens .. logical (default false): depthInsensitivity from
        %%%          the substrate's (layer 1) block of dK/dD alone instead of all layers.
        %%%          RitzSolverSAW.calc sets it for layered SAW (nLay>1), since
        %%%          depthInsensitivity must track the substrate's truncation depth, not
        %%%          the film's real thickness. See docs/reference/architecture/ritz.md
        %%%          ("Multi-layer stacks").
        %%%
        %%% The quotients alpha^H (N.' dK/dD N) alpha are read per layer block on the
        %%% full-space eigenvectors U = N*alpha (N is real): U^H dK/dD U. dK/dD is
        %%% block-diagonal in layer order, so the quotient is a sum over layers of
        %%% q_i = U_i^H (dK/dD)_i U_i, and every weighting below is a different
        %%% combination of the q_i. This replaces projecting dK/dD (and scaled and
        %%% substrate-only copies of it) into the constrained basis and taking
        %%% diag(V'*A*V), an O(nCon^3) product to read O(nCon) numbers.
            if nargin < 3 || isempty(substrateOnlyInsens); substrateOnlyInsens = false; end

            U  = obj.nullConstrMat * obj.eigVecs_con;
            q  = RitzSolver.layerQuotients(U, gm.dGammaDMatAll, obj.structure.calcParams);

            %%% Rayleigh quotient derivative: dω²/dD = α^H (dK/dD) α
            %%% dv/dD = dω²/dD / (8π² v)    (relLength=1)
            qSum = sum(q, 1).';
            obj.dVel_dDepth  = qSum ./ obj.vels / (8*pi^2);

            %%% vg = dω/dk. v depends on depths only through k·d_i, so
            %%% vg = vp + Σ_i d_i ∂vp/∂d_i — each layer weighted by ITS OWN
            %%% relDepth (the d_i ρ_i normFac terms cancel only in this weighted
            %%% sum). dKdD is blkdiag(∂K_i/∂d_i) in layer order, so weighting the
            %%% quotient q_i of block i by d_i gives Σ_i d_i ∂K/∂d_i.
            %%% For nLay==1 this is exactly relDepth * dVel_dDepth. For nLay>1,
            %%% dVel_dDepth alone (every layer +1 at once) is NOT a group-velocity
            %%% term: vp + D_total·dVel_dDepth would over-weight thin films by
            %%% D_total/d_film. Agrees with a finite-difference dω/dk.
            relDepths = [obj.structure.calcParams.relDepth];
            if isscalar(relDepths)
                obj.velsGroup = obj.vels + relDepths .* obj.dVel_dDepth;
            else
                qWeighted = (relDepths(:).' * q).';
                obj.velsGroup = obj.vels + qWeighted ./ obj.vels / (8*pi^2);
            end

            if substrateOnlyInsens
                qIns = q(1, :).';
            else
                qIns = qSum;
            end
            dVel_dDepth_forInsens = qIns ./ obj.vels / (8*pi^2);
            obj.depthInsensitivity = log10(real(1 ./ max(abs(dVel_dDepth_forInsens), eps)));
        end

        function oop = topSurfaceU3(obj, nM)
        %%% |u3| at the top surface of the first nM modes, without evaluating the profile
        %%% >> oop = topSurfaceU3(obj, nM)
        %%%
        %%% OUTPUT: oop .. [nM x 1] |out-of-plane displacement| on the last depth point of
        %%%                the top layer (x3 = surface normal)
        %%%
        %%% The last row of calcDisplacement(nM)(:,:,3), read directly: the top layer's
        %%% u3 Legendre coefficients (rows of N*alpha) weighted by its polyLeg at the
        %%% surface. calcDisplacement applies a phase fix per mode, which cannot change
        %%% an absolute value, so it is skipped here.
            cps   = obj.structure.calcParams;
            cp    = cps(end);
            legO  = cp.legOrder;
            rows  = sum([cps.basisVecL]) - legO + (1:legO);   % u3 coefficients of the top layer
            coef  = obj.nullConstrMat(rows, :) * obj.eigVecs_con(:, 1:nM);
            oop   = abs((cp.polyLeg(:, end).' * coef).');
        end
    end

    methods (Static, Access = protected)
        function q = layerQuotients(U, K, calcParams)
        %%% Per-layer Rayleigh quotients of a block-diagonal matrix on full-space vectors
        %%% >> q = RitzSolver.layerQuotients(U, K, calcParams)
        %%%
        %%% INPUT:  1. U          .. [strBVL x nModes] vectors in the full basis (N*alpha)
        %%%         2. K          .. [strBVL x strBVL] block-diagonal in layer order
        %%%                          (gammaMatAll, dGammaDMatAll, ...); off-block entries
        %%%                          are not read
        %%%         3. calcParams .. structure.calcParams (per-layer basisVecL)
        %%% OUTPUT: q .. [nLay x nModes] real, q(i, m) = Re(U_i(:,m)' K_i U_i(:,m)) with
        %%%              U_i, K_i the rows/diagonal block of layer i
            nLay = numel(calcParams);
            q    = zeros(nLay, size(U, 2));
            row0 = 0;
            for ilay = 1:nLay
                rows = row0 + (1 : calcParams(ilay).basisVecL);
                row0 = rows(end);
                if nLay == 1
                    Ub = U;
                    Kb = K;
                else
                    Ub = U(rows, :);
                    Kb = K(rows, rows);
                end
                q(ilay, :) = real(sum(conj(Ub) .* (Kb * Ub), 1));
            end
        end

        function ph = phaseFix(eigVecs, calcParams)
        %%% Per-mode phase factor making the largest in-plane coefficient real positive
        %%% >> ph = RitzSolver.phaseFix(eigVecs, calcParams)
        %%%
        %%% INPUT:  1. eigVecs    .. [strBVL x nModes] eigenvectors in the full basis
        %%%         2. calcParams .. structure.calcParams (per-layer legOrder/basisVecL)
        %%% OUTPUT: ph .. [1 x nModes], |ph| = 1; eigVecs .* ph has the fixed phase
        %%%
        %%% The reference is taken on the Legendre coefficients, not on a sampled
        %%% profile, so it does not depend on the evaluation grid.
            isInPlane = false(size(eigVecs, 1), 1);
            row0 = 0;
            for ilay = 1:numel(calcParams)
                isInPlane(row0 + (1 : 2*calcParams(ilay).legOrder)) = true;
                row0 = row0 + calcParams(ilay).basisVecL;
            end
            nM = size(eigVecs, 2);
            [mxIn, iIn]   = max(abs(eigVecs(isInPlane, :)), [], 1);
            [mxAll, iAll] = max(abs(eigVecs), [], 1);
            rowsIn = find(isInPlane);
            refRow = rowsIn(iIn(:)).';
            noInPlane = mxIn <= 1e-8 * mxAll;     % e.g. a purely u3 mode
            refRow(noInPlane) = iAll(noInPlane);
            ref = eigVecs(sub2ind(size(eigVecs), refRow, 1:nM));
            ph  = conj(ref) ./ max(abs(ref), realmin);
        end

        function checkFlag(flag, errId, argName)
        %%% Error errId unless flag is a logical scalar (or numeric 0 / 1)
        %%% >> RitzSolver.checkFlag(flag, errId, argName)
            if ~(isscalar(flag) && (islogical(flag) || (isnumeric(flag) && (flag == 0 || flag == 1))))
                error(errId, '%s must be true or false.', argName);
            end
        end

        function val = precompField(precompStruct, fieldName, defaultFcn)
        %%% Read one field of a precompute struct, or compute a fallback lazily
        %%% >> val = RitzSolver.precompField(precompStruct, fieldName, defaultFcn)
        %%%
        %%% INPUT:  1. precompStruct .. the caller's precompStruct (may be empty
        %%%                             struct(), may lack the field entirely)
        %%%         2. fieldName     .. field to read
        %%%  (opt)  3. defaultFcn    .. zero-argument function handle evaluated
        %%%                             ONLY when the field is absent or empty;
        %%%                             omit to get [] in that case
        %%%
        %%% defaultFcn is a handle, not a value, because the fallbacks are
        %%% expensive (a null-space solve, the reference-plate eigenvectors)
        %%% and must not run when the caller already supplied the field.
            if isfield(precompStruct, fieldName) && ~isempty(precompStruct.(fieldName))
                val = precompStruct.(fieldName);
            elseif nargin < 3 || isempty(defaultFcn)
                val = [];
            else
                val = defaultFcn();
            end
        end
    end

    methods (Static)
        function [eigVals, V] = solveConstrainedEig(K_con, M_con)
        %%% Sorted eigenpairs of K_con α = ω² M_con α on MATLAB's Hermitian path
        %%% >> [eigVals, V] = RitzSolver.solveConstrainedEig(K_con, M_con)
        %%%
        %%% INPUT:  1. K_con .. [n x n] constrained stiffness (Hermitian up to round-off)
        %%%         2. M_con .. [n x n] constrained mass (real symmetric positive
        %%%                     definite up to round-off), or [] for the identity
        %%% OUTPUT: eigVals .. [n x 1] real, ascending
        %%%         V       .. [n x n] eigenvectors in the same order, M-orthonormal
        %%%                    (α^H M_con α = 1); each column's phase is the solver's choice
        %%%
        %%% The one place the Ritz eigen-solve happens: RitzSolver.solveEigenvalueProblem
        %%% and ElasticityInversion's fit loop both call it.
        %%%
        %%% eig tests Hermitian symmetry EXACTLY. N.'*K*N and N.'*M*N miss it by ~1e-16
        %%% round-off, which would send the solve to the general solvers (complex Schur;
        %%% QZ for eig(K,M)): up to 12x slower on stacks, and the eigenvectors are only
        %%% approximately M-orthogonal. Symmetrising first selects the Hermitian (eig(K))
        %%% or Cholesky (eig(K,M)) path; eigenvalues move by <= 1e-8 relative.
            K_con = (K_con + K_con') / 2;

            %%% M_con is often mathematically the identity (e.g. single-layer SAW),
            %%% but routing eye(n) through N.'*normMatAll*N with N = null(constrMat)
            %%% leaves ~1e-15 roundoff on it. The test must therefore be a tolerance
            %%% against the FULL matrix — off-diagonal projection roundoff matters
            %%% too — not exact equality on the diagonal, which would almost always
            %%% fail and force every such solve onto the generalised eigensolver.
            isIdentity = ~isempty(M_con) && ...
                max(max(abs(M_con - eye(size(M_con))))) < 1e-8;
            if isempty(M_con) || isIdentity
                [V, D] = eig(K_con);
            else
                M_con  = real(M_con + M_con') / 2;   % N and M are real
                [V, D] = eig(K_con, M_con);
                %%% renormalise so eigenvectors are M-orthonormal (α^H M α = 1); the
                %%% Cholesky path already returns them so, this only guards a QZ fallback
                %%% (MATLAB uses QZ if M_con is not numerically positive definite)
                aT_M_a = real(sum(conj(V) .* (M_con * V), 1));
                V      = V ./ sqrt(aT_M_a);
            end

            [eigVals, ord] = sort(real(diag(D)));
            V = V(:, ord);
        end

        function [K, d] = realGauge(K, calcParams)
        %%% Real-symmetric gauge u3 -> i u3 of the full stiffness matrix, if valid
        %%% >> [K, d] = RitzSolver.realGauge(K, calcParams)
        %%%
        %%% INPUT:  1. K          .. [strBVL x strBVL] Hermitian stiffness (gammaMatAll)
        %%%         2. calcParams .. structure.calcParams (per-layer legOrder/basisVecL)
        %%% OUTPUT: K .. D^H K D, REAL symmetric, if every layer passes; else K unchanged
        %%%         d .. [strBVL x 1] diagonal of D (1 on the u1, u2 coefficients, i on
        %%%              u3), or [] when the gauge does not apply
        %%%
        %%% If every layer's elastic tensor is invariant under a 2-fold rotation about
        %%% the surface normal (isotropic, cubic (001)/(110), hexagonal, ...), the
        %%% imaginary part of K sits only in the in-plane/u3 blocks and D^H K D is real.
        %%% The real symmetric eigen-solve is 1.6-3x faster than the Hermitian one and
        %%% gives eigenvectors fixed up to a sign. The test needs no symmetry
        %%% bookkeeping, but it must be at round-off level: for cubic (111) the residual
        %%% is ~3e-5 of ||K|| and forcing it real shifts eigenvalues by 9-19 %.
            d  = RitzSolver.gaugeVector(calcParams);
            Kg = conj(d) .* K .* d.';
            if RitzSolver.isGaugeReal(Kg)
                K = real(Kg);
            else
                d = [];
            end
        end

        function d = gaugeVector(calcParams)
        %%% Diagonal of the gauge D = diag(1, 1, i) per Legendre coefficient (see realGauge)
        %%% >> d = RitzSolver.gaugeVector(calcParams)
        %%% OUTPUT: d .. [strBVL x 1]: 1 on every layer's u1, u2 coefficients, i on u3
            d = cell(numel(calcParams), 1);
            for ilay = 1:numel(calcParams)
                L = calcParams(ilay).legOrder;
                d{ilay} = [ones(2*L, 1); 1i * ones(L, 1)];
            end
            d = vertcat(d{:});
        end

        function tf = isGaugeReal(Kg, scale)
        %%% Is a gauged matrix real to round-off? (see realGauge)
        %%% >> tf = RitzSolver.isGaugeReal(Kg)
        %%% >> tf = RitzSolver.isGaugeReal(Kg, scale)   % scale defaults to ||Kg||_F
            if nargin < 2; scale = norm(Kg, 'fro'); end
            tf = norm(imag(Kg), 'fro') <= RitzSolver.REAL_GAUGE_TOL * scale;
        end

        function nullConstrMat = buildConstraintMatrix(calcParams, bcType)
        %%% Build null space of boundary-condition constraint matrix
        %%% >> nullConstrMat = RitzSolver.buildConstraintMatrix(calcParams, bcType)
        %%%
        %%% INPUT:  1. calcParams .. CalcParams array (one per layer)
        %%%         2. bcType     .. 'saw' (rigid bottom) or 'lamb' (free plate)
        %%% OUTPUT: N — null space of constraint matrix   [strBVL x (strBVL - nConstraints)]

            noLays  = length(calcParams);
            legO    = [calcParams.legOrder];
            bVL     = [calcParams.basisVecL];
            delta   = eye(max(legO));

            strBVL  = sum(bVL);
            constrMat = complex(zeros(strBVL));   % rows = constraints, cols = basis; zero rows ignored by null()
            constrNo  = 1;

            %%% SAW: rigid bottom boundary u_i = 0 at the bottom of layer 1
            if strcmp(bcType, 'saw')
                ilay = 1;
                for iiD = 1:3
                    for iiL = 1:bVL(ilay)
                        iL  = mod(iiL-1, legO(ilay));
                        iD  = ceil(iiL / legO(ilay));
                        constrMat(constrNo, iiL) = delta(iD, iiD) * ((-1)^iL) * sqrt((2*iL+1)/2);
                    end
                    constrNo = constrNo + 1;
                end
            end

            %%% Inter-layer displacement compatibility
            if noLays > 1
                for ilay = 1:noLays-1
                    if ilay == 1
                        legIntLow = 1 : bVL(ilay);
                    else
                        legIntLow = sum(bVL(1:ilay-1))+1 : sum(bVL(1:ilay));
                    end
                    legIntHigh = sum(bVL(1:ilay))+1 : sum(bVL(1:ilay+1));

                    for iiD = 1:3
                        %%% top of lower layer (value at +1): sum_m P_m(+1) c_m = sum_m sqrt((2m+1)/2) c_m
                        for iiL = legIntLow
                            iL = mod((iiL - legIntLow(1)+1)-1, legO(ilay));
                            iD = ceil((iiL - legIntLow(1)+1) / legO(ilay));
                            constrMat(constrNo, iiL) = delta(iD,iiD) * sqrt((2*iL+1)/2);
                        end
                        %%% bottom of upper layer (value at -1): sum_m P_m(-1) c_m = sum_m (-1)^m sqrt((2m+1)/2) c_m
                        for iiL = legIntHigh
                            iL = mod((iiL - legIntLow(end))-1, legO(ilay+1));
                            iD = ceil((iiL - legIntLow(end)) / legO(ilay+1));
                            constrMat(constrNo, iiL) = -delta(iD,iiD) * (-1)^iL * sqrt((2*iL+1)/2);
                        end
                        constrNo = constrNo + 1;
                    end
                end
            end

            %%% null space — the vectors orthogonal to all constraints
            if constrNo == 1
                %%% no constraints were added (Lamb single-layer): N = I
                nullConstrMat = eye(strBVL);
            else
                nullConstrMat = null(constrMat);
            end
        end
    end
end
