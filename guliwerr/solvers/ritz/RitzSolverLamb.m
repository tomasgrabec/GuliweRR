classdef RitzSolverLamb < RitzSolver
% RitzSolverLamb — Ritz solver for the plate (Lamb waves).
%
%%% GuliweRR v6.0

    properties (SetAccess = protected)
        modesType       double          % 1=Asymm, 2=Symm, 3=TH-Asymm, 4=TH-Symm: the largest of
                                         % the four mid-plane parity fractions (modeRecognition)
        modeTypeConfidence double       % top parity fraction minus the runner-up, in [0,1]. 1 = a
                                         % pure type; low = a mixed mode: parity mixing (a stack or
                                         % crystal cut without a mid-plane mirror) or, within one
                                         % parity class, SH-Lamb polarisation coupling
        realThickness   double          % physical plate thickness [same units as realWavVec^-1]
        realWavVec(1,1) double  = 0     % physical wavevector k [rad/m] (or consistent unit)
        realWavelength  double          % = 2π / realWavVec
        freqsPhysical   double          % real frequencies = vels * realWavVec / (2π)
        velsGroupTrans  double          % transverse group velocity c_gθ = (1/k) ∂ω/∂θ = ∂vp/∂θ
                                         % [same unit as vels], per mode: the in-plane component
                                         % of the group velocity perpendicular to k, positive
                                         % towards increasing propagation angle. velsGroup is
                                         % the component along k. Zero on a direction with a
                                         % sagittal mirror (e.g. <100>, <110> of cubic (001))
    end

    methods
        function obj = RitzSolverLamb(structure, angleDeg, realThickness, realWavVec)
        %%% Initialise RitzSolverLamb
        %%% >> obj = RitzSolverLamb(structure, angleDeg, realThickness, realWavVec)
        %%%
        %%% INPUT:  (opt) 1. structure     .. Structure (the plate's layers)
        %%%         (opt) 2. angleDeg      .. propagation angle [deg] (default 0)
        %%%         (opt) 3. realThickness .. physical layer thicknesses
        %%%         (opt) 4. realWavVec    .. physical wavevector k; 0 or [] = work in the
        %%%                                   normalised relDepth units. A nonzero k needs
        %%%                                   realThickness (error
        %%%                                   RitzSolverLamb:constructor:needThickness): the
        %%%                                   relDepths are rescaled so that
        %%%                                   k*d = 2*pi*relDepth.
            %%% RitzSolver's constructor guards each argument with ~isempty(),
            %%% so [] is equivalent to omitting the argument.
            if nargin < 1; structure = []; end
            if nargin < 2; angleDeg  = []; end
            obj@RitzSolver(structure, angleDeg);

            if nargin > 2 && ~isempty(realThickness)
                obj.realThickness = realThickness;
            end
            if nargin > 3 && ~isempty(realWavVec) && realWavVec ~= 0
                if isempty(obj.realThickness)
                    error('RitzSolverLamb:constructor:needThickness', ...
                        'A physical wavevector needs the physical thickness (realThickness).');
                end
                obj.realWavVec     = realWavVec;
                obj.realWavelength = 2*pi / realWavVec;

                %%% scale relDepth so that k*d = 2π*relDepth (i.e., relDepth = d/(lambda))
                if obj.structure.nLay > 0
                    totalThick    = sum(obj.realThickness);
                    totalRelDepth = sum([obj.structure.calcParams.relDepth]);
                    scaleFactor   = (totalThick * realWavVec / (2*pi)) / totalRelDepth;
                    if abs(scaleFactor - 1) > 1e-14
                        for ilay = 1:obj.structure.nLay
                            newRD = obj.structure.calcParams(ilay).relDepth * scaleFactor;
                            newCP = obj.structure.calcParams(ilay).SetRelDepth(newRD);
                            obj.structure = obj.structure.ModifyLayer(ilay, [], newCP);
                        end
                    end
                end
            end
        end

        function obj = calc(obj, withModeTypes, precompStruct)
        %%% Calculate Lamb modes for plate domain
        %%% >> obj = calc(obj, withModeTypes, precompStruct)
        %%%
        %%% INPUT:  (opt) 1. withModeTypes .. logical: also classify the mode types
        %%%                                   (modeRecognition -> modesType,
        %%%                                   modeTypeConfidence; default false)
        %%%         (opt) 2. precompStruct .. .nullConstrMat, .rawByDepthCell,
        %%%                                   .gammaByDepthCell (each built if absent;
        %%%                                   rawByDepthCell depends on legOrder only,
        %%%                                   gammaByDepthCell also on angle and material)

            if nargin < 2 || isempty(withModeTypes); withModeTypes = false; end
            RitzSolver.checkFlag(withModeTypes, 'RitzSolverLamb:calc:badFlag', 'withModeTypes');
            if nargin < 3 || isempty(precompStruct);   precompStruct   = struct(); end

            %%% constraint matrix (Lamb: free surfaces, inter-layer compatibility only)
            obj.nullConstrMat = RitzSolver.precompField(precompStruct, 'nullConstrMat', ...
                @() RitzSolver.buildConstraintMatrix(obj.structure.calcParams, 'lamb'));

            %%% Legendre integrals (material-, angle- and depth-independent) and the
            %%% elasticity-contracted pages for the Lamb path. rawByDepthCell is built here
            %%% rather than inside calcLamb because calcTransGroupVelocity needs it too.
            rawByDepthCell   = RitzSolver.precompField(precompStruct, 'rawByDepthCell', ...
                @() GammaMatrix.RawCellFor(obj.structure));
            gammaByDepthCell = RitzSolver.precompField(precompStruct, 'gammaByDepthCell', ...
                @() GammaMatrix.elastByDepth(obj.structure, obj.angle, rawByDepthCell));

            %%% assemble K and M via calcLamb -- a local, not stored
            %%% (see RitzSolver.projectToConstrained)
            gm = GammaMatrix();
            gm = gm.calcLamb(obj.structure, obj.angle, gammaByDepthCell);

            %%% project into constrained basis. Single-layer Lamb has no
            %%% constraints (buildConstraintMatrix returns the identity), so the
            %%% projection is a no-op and skipping it also lets
            %%% solveEigenvalueProblem take its cheaper M_con=[] path.
            if obj.structure.nLay == 1
                [K_con, gauge] = RitzSolver.realGauge(gm.gammaMatAll, obj.structure.calcParams);
                M_con    = [];
            else
                [K_con, M_con, gauge] = obj.projectToConstrained(gm);
            end

            %%% solve
            obj = obj.solveEigenvalueProblem(K_con, M_con, gauge);
            obj = obj.calcGroupVelocity(gm);
            obj = obj.calcTransGroupVelocity(rawByDepthCell);

            %%% physical frequencies
            if obj.realWavVec ~= 0
                obj.freqsPhysical = obj.vels .* obj.realWavVec / (2*pi);
            end

            %%% mode-type classification
            if withModeTypes && obj.nCompModes > 0
                obj = obj.modeRecognition();
            end
        end

        function obj = calcTransGroupVelocity(obj, rawByDepthCell)
        %%% Transverse group velocity c_gθ of every computed mode (Hellmann-Feynman in angle)
        %%% >> obj = calcTransGroupVelocity(obj, rawByDepthCell)
        %%%
        %%% INPUT: (opt) rawByDepthCell .. GammaMatrix.RawCellFor(obj.structure) (built if
        %%%                                omitted)
        %%% Fills: velsGroupTrans [same unit as vels]
        %%%
        %%% c_gθ = (1/k) ∂ω/∂θ = ∂vp/∂θ (θ in rad, k fixed): the group-velocity component
        %%% perpendicular to k, in the plate plane, positive towards increasing angle.
        %%% Same formula as RitzSolver.angleSlope (the Rayleigh-quotient derivative
        %%% dv/dθ = α^H (dK/dθ) α / (8π² v), exact for a non-degenerate mode; dK/dθ exact
        %%% from the ±0.01° central difference of the rotated constants alone), but built
        %%% the Lamb path's way: K is linear in the 15 reduced constants c, so
        %%% elastByDepth contracted with dc/dθ gives the three depth pages of dK/dθ
        %%% directly, one contraction instead of angleSlope's 15 dK/dc pages
        %%% (GammaMatrix.derivsByDepth). Then the same depth factors, density and layer
        %%% normalisation as calcLamb. No extra solve. Agrees with a ±0.01° finite
        %%% difference of the solved velocities and, to four decimals on a Cu(001) plate,
        %%% with GEWtool's transverse energy velocity.
            if nargin < 2 || isempty(rawByDepthCell)
                rawByDepthCell = GammaMatrix.RawCellFor(obj.structure);
            end
            nM = obj.nCompModes;
            if nM == 0
                obj.velsGroupTrans = zeros(0, 1);
                return
            end
            h  = 0.01;   % deg, as RitzSolver.angleSlope
            cP = GammaMatrix.RedElastCellFor(obj.structure, obj.angle + h);
            cM = GammaMatrix.RedElastCellFor(obj.structure, obj.angle - h);
            dc = cell(obj.structure.nLay, 1);
            for ilay = 1:obj.structure.nLay
                dc{ilay} = (cP{ilay} - cM{ilay}) / (2 * h);
            end
            dPages = GammaMatrix.elastByDepth(obj.structure, obj.angle, rawByDepthCell, dc);

            U    = obj.nullConstrMat * obj.eigVecs_con(:, 1:nM);
            acc  = zeros(1, nM);
            row0 = 0;
            for ilay = 1:obj.structure.nLay
                cp  = obj.structure.calcParams(ilay);
                pg  = dPages{ilay};
                dK  = GammaMatrix.ContractDepth(pg(:,:,1), pg(:,:,2), pg(:,:,3), cp.relDepth) ...
                      ./ obj.structure.mater(ilay).rho .* GammaMatrix.LayerNormFac(obj.structure, ilay);
                Ub  = U(row0 + (1:cp.basisVecL), :);
                acc = acc + real(sum(conj(Ub) .* (dK * Ub), 1));
                row0 = row0 + cp.basisVecL;
            end
            slopePerDeg        = acc.' ./ obj.vels(1:nM) / (8 * pi^2);   % dv/dθ per degree
            obj.velsGroupTrans = slopePerDeg * (180 / pi);
        end

        function obj = modeRecognition(obj)
        %%% Classify modes as Asymm / Symm / TH-Asymm / TH-Symm by mid-plane parity
        %%% >> obj = modeRecognition(obj)
        %%%
        %%% Fills modesType (argmax of RitzSolverLamb.parityFractions) and
        %%% modeTypeConfidence (largest fraction minus the second largest).
        %%% Physically defined, so it needs no reference plate: a mode is labelled by
        %%% where its (kinetic-energy-weighted) displacement sits once each component
        %%% is split into its even and odd parts about the plate's mid-plane. See
        %%% docs/reference/architecture/plate.md, "Mode-type classification".
            U = obj.eigVecs_con(:, 1:obj.nCompModes);
            if obj.structure.nLay > 1          % one layer: nullConstrMat = I
                U = obj.nullConstrMat * U;
            end
            fr = RitzSolverLamb.parityFractions(U, obj.structure);
            [frSorted, idx] = sort(fr, 1, 'descend');
            obj.modesType          = idx(1, :);
            obj.modeTypeConfidence = max(frSorted(1, :) - frSorted(2, :), 0);
        end
    end

    methods (Static)
        function fr = parityFractions(U, structure)
        %%% Energy fractions of the four mid-plane parity/polarisation classes per mode
        %%% >> fr = RitzSolverLamb.parityFractions(U, structure)
        %%%
        %%% INPUT:  1. U         .. [strBVL x nModes] eigenvectors in the full Legendre
        %%%                         basis (nullConstrMat * eigVecs_con)
        %%%         2. structure .. the plate's Structure (relDepth, legOrder, rho per layer)
        %%% OUTPUT: fr .. [4 x nModes], rows in modesType order, each column summing to 1:
        %%%               1 Asymm    = E(u1 odd)  + E(u3 even)
        %%%               2 Symm     = E(u1 even) + E(u3 odd)
        %%%               3 TH-Asymm = E(u2 odd)
        %%%               4 TH-Symm  = E(u2 even)
        %%%               E(.) = int rho_s |.|^2 dz / int rho_s |u|^2 dz, even/odd about the
        %%%               plate's mid-plane, rho_s(z) = (rho(z) + rho(D-z))/2
        %%%
        %%% The mirror z -> D-z maps a symmetric (S0-like) field to u1, u2 even and u3
        %%% odd, an antisymmetric (A0-like) one to u1, u2 odd and u3 even. The u2 parts
        %%% are kept as their own two classes (TH-Symm = u2 even, the SH0 family;
        %%% TH-Asymm = u2 odd, SH1), so with SH-Lamb coupling a symmetric mode splits its
        %%% energy between Symm and TH-Symm (never into an antisymmetric class when the
        %%% plate has a mirror plane) and the label names the dominant polarisation.
        %%%
        %%% Weight: the even and odd parts are orthogonal under any weight that is even
        %%% about the mid-plane, so the symmetrised density rho_s makes the four fractions
        %%% sum to 1 exactly; it is the kinetic-energy weight itself whenever the stack is
        %%% mirror-symmetric (always for one layer, where rho cancels).
        %%%
        %%% One layer: the normalised Legendre polynomials are orthonormal on the layer
        %%% and P_n(-x) = (-1)^n P_n(x) about its centre (= the plate's mid-plane), so
        %%% E(even) and E(odd) are sums of |coefficient|^2 over even / odd n -- read off
        %%% the eigenvector, no depth grid. Several layers: the displacement at depth z
        %%% and at its mirror D-z, on a Gauss-Legendre grid over the lower half of the
        %%% plate whose sub-intervals end at every interface and every mirrored
        %%% interface. Within a sub-interval u(z) and u(D-z) are single-layer
        %%% polynomials of degree < max(legOrder), so max(legOrder) nodes per
        %%% sub-interval integrate the products exactly.
            cps  = structure.calcParams;
            nLay = numel(cps);
            nM   = size(U, 2);
            E    = zeros(6, nM);   % rows: u1e u1o u2e u2o u3e u3o
            if nLay == 1
                L = cps(1).legOrder;
                for ic = 1:3
                    rows = (ic-1)*L + (1:L);
                    A2 = abs(U(rows, :)).^2;
                    E(2*ic-1, :) = sum(A2(1:2:end, :), 1);   % n = 0, 2, ...: even
                    E(2*ic,   :) = sum(A2(2:2:end, :), 1);   % n = 1, 3, ...: odd
                end
            else
                [zLow, w] = RitzSolverLamb.mirrorQuadrature(cps);
                zMir = 1 - zLow;
                tB   = [0, cumsum([cps.relDepth])] / sum([cps.relDepth]);
                rhoL = [structure.mater.rho].';
                layLow = RitzSolverLamb.layerOf(zLow, tB);
                layMir = RitzSolverLamb.layerOf(zMir, tB);
                w = w .* (rhoL(layLow) + rhoL(layMir)) / 2;
                row0 = [0, cumsum([cps.basisVecL])];
                uLow = zeros(numel(zLow), nM);
                uMir = zeros(numel(zLow), nM);
                for ic = 1:3
                    uLow(:) = 0;  uMir(:) = 0;
                    for ilay = 1:nLay
                        L    = cps(ilay).legOrder;
                        rows = row0(ilay) + (ic-1)*L + (1:L);
                        q = layLow == ilay;
                        if any(q)
                            uLow(q, :) = RitzSolverLamb.normLegendre(L, ...
                                RitzSolverLamb.localCoord(zLow(q), tB, ilay)) * U(rows, :);
                        end
                        q = layMir == ilay;
                        if any(q)
                            uMir(q, :) = RitzSolverLamb.normLegendre(L, ...
                                RitzSolverLamb.localCoord(zMir(q), tB, ilay)) * U(rows, :);
                        end
                    end
                    E(2*ic-1, :) = w.' * abs(uLow + uMir).^2;
                    E(2*ic,   :) = w.' * abs(uLow - uMir).^2;
                end
            end
            fr = [E(2,:) + E(5,:); E(1,:) + E(6,:); E(4,:); E(3,:)];
            fr = fr ./ max(sum(fr, 1), realmin);
        end

        function [velGroup, objOut] = fastCalcGroupVel(structure, angleDeg, realThickness, modeNoIn, wavVecIn, precompStruct)
        %%% Group velocity (along k) of one raw mode at one wavevector -- a thin wrapper for checks
        %%% >> [velGroup, objOut] = RitzSolverLamb.fastCalcGroupVel(structure, angleDeg, realThickness, modeNoIn, wavVecIn, precompStruct)
            if nargin < 6; precompStruct = struct(); end   % as in locateZGV below
            fc = RitzSolverLamb(structure, angleDeg, realThickness, wavVecIn);
            fc = fc.calc(false, precompStruct);
            velGroup = fc.velsGroup(modeNoIn);
            if nargout > 1; objOut = fc; end
        end

        function [wavVecZGV, freqZGV, objOut] = locateZGV(structure, angleDeg, realThickness, modeNoIn, wavVec0, precompStruct)
        %%% Locate the zero-group-velocity point nearest a seed on one Lamb mode (Newton)
        %%% >> [wavVecZGV, freqZGV, objOut] = RitzSolverLamb.locateZGV(structure, angleDeg, realThickness, modeNoIn, wavVec0, precompStruct)
        %%%
        %%% INPUT:  1. structure     .. Structure (UNROTATED; angleDeg rotates it)
        %%%         2. angleDeg      .. propagation angle [deg]
        %%%         3. realThickness .. physical plate thickness (as the constructor)
        %%%         4. modeNoIn      .. raw Ritz mode number of the seed at wavVec0
        %%%         5. wavVec0       .. seed wavevector (physical, > 0)
        %%%  (opt)  6. precompStruct .. .nullConstrMat, .rawByDepthCell (built if absent)
        %%% OUTPUT: wavVecZGV, freqZGV .. the ZGV point (physical frequency); NaN when the
        %%%                               Newton iteration does not converge or lands
        %%%                               outside [wavVec0/10, wavVec0*10]
        %%%         objOut             .. RitzSolverLamb solved at wavVecZGV (at wavVec0
        %%%                               when not converged)
        %%%
        %%% Newton refinement (ZgvNewton.solve) on the plate operator as a polynomial in k
        %%% (kPolynomial), seeded with mode modeNoIn's frequency and eigenvector at
        %%% wavVec0. It converges onto the nearest point where (k, ω) is a ZGV point of
        %%% ANY branch -- no mode number is held fixed, so a branch crossing between the
        %%% seed and the ZGV does not matter. freqZGV is the Newton frequency; the mode it
        %%% belongs to in objOut need not be modeNoIn.
            if nargin < 6 || isempty(precompStruct); precompStruct = struct(); end
            if ~isfield(precompStruct, 'nullConstrMat') || isempty(precompStruct.nullConstrMat)
                precompStruct.nullConstrMat = RitzSolver.buildConstraintMatrix(structure.calcParams, 'lamb');
            end
            if ~isfield(precompStruct, 'rawByDepthCell') || isempty(precompStruct.rawByDepthCell)
                precompStruct.rawByDepthCell = GammaMatrix.RawCellFor(structure);
            end
            rs0 = RitzSolverLamb(structure, angleDeg, realThickness, wavVec0);
            rs0 = rs0.calc(false, precompStruct);
            poly = RitzSolverLamb.kPolynomial(structure, angleDeg, realThickness, precompStruct);
            [wavVecZGV, omegaZGV, ~, isConverged] = ZgvNewton.solve(poly.L2, poly.L1, poly.L0, poly.M, ...
                wavVec0, 2*pi * rs0.freqsPhysical(modeNoIn));
            if ~isConverged || wavVecZGV < wavVec0/10 || wavVecZGV > wavVec0*10
                wavVecZGV = NaN;
                freqZGV   = NaN;
                objOut    = rs0;
                return
            end
            freqZGV = omegaZGV / (2*pi);
            if nargout > 2
                objOut = RitzSolverLamb(structure, angleDeg, realThickness, wavVecZGV);
                objOut = objOut.calc(false, precompStruct);
            end
        end

        function poly = kPolynomial(structure, angleDeg, realThickness, precompStruct)
        %%% The plate's Ritz eigenproblem as a quadratic polynomial in the wavevector
        %%% >> poly = RitzSolverLamb.kPolynomial(structure, angleDeg, realThickness, precompStruct)
        %%%
        %%% INPUT:  1. structure     .. Structure (UNROTATED; angleDeg rotates it)
        %%%         2. angleDeg      .. propagation angle [deg]
        %%%         3. realThickness .. physical plate thickness (only its sum is used, as
        %%%                             in the constructor; the layers keep the relDepth
        %%%                             ratios of structure.calcParams)
        %%%  (opt)  4. precompStruct .. .nullConstrMat, .rawByDepthCell (built if absent)
        %%% OUTPUT: poly .. struct with L2, L1, L0, M ([n x n], constrained basis N), such that
        %%%                    (k² L2 + k L1 + L0) α = ω² M α
        %%%                 at every physical wavevector k [1/length] with ω = 2π f [rad per
        %%%                 time, f = freqsPhysical]; all Hermitian (L2, L0, M real in every
        %%%                 case checked, L1 imaginary: it carries the i). Also N
        %%%                 (nullConstrMat).
        %%%
        %%% At fixed physical layer thicknesses h_i the relDepths scale with k
        %%% (relDepth_i = h_i k / 2π), and the K = K0 + K1/(relDepth/2) + K2/(relDepth/2)²
        %%% of GammaMatrix.ContractDepth, rescaled from the normalised eigenvalue
        %%% (2π v)² = (2π ω / k)² to ω², splits per layer block into
        %%%     L2 = K0/(2π)²,  L1 = K1/(π h_i),  L0 = 4 K2/h_i²   (each · normFac_i/ρ_i)
        %%% with normFac_i = h_i ρ_i / (h_1 ρ_1) as GammaMatrix.LayerNormFac and M the
        %%% normFac-weighted identity. The same matrices RitzSolverLamb.calc solves, to
        %%% round-off.
            if nargin < 4 || isempty(precompStruct); precompStruct = struct(); end
            N   = RitzSolver.precompField(precompStruct, 'nullConstrMat', ...
                @() RitzSolver.buildConstraintMatrix(structure.calcParams, 'lamb'));
            raw = RitzSolver.precompField(precompStruct, 'rawByDepthCell', ...
                @() GammaMatrix.RawCellFor(structure));
            pages = GammaMatrix.elastByDepth(structure, angleDeg, raw);

            cps  = structure.calcParams;
            nLay = numel(cps);
            relD = [cps.relDepth];
            h    = sum(realThickness) * relD / sum(relD);   % physical layer thicknesses
            rho  = [structure.mater.rho];
            b2 = cell(1, nLay);  b1 = cell(1, nLay);  b0 = cell(1, nLay);  bm = cell(1, nLay);
            for ilay = 1:nLay
                nf = (h(ilay) * rho(ilay)) / (h(1) * rho(1));
                w  = nf / rho(ilay);
                pg = pages{ilay};
                b2{ilay} = pg(:,:,1) * (w / (2*pi)^2);
                b1{ilay} = pg(:,:,2) * (w / (pi * h(ilay)));
                b0{ilay} = pg(:,:,3) * (w * 4 / h(ilay)^2);
                bm{ilay} = eye(cps(ilay).basisVecL) * nf;
            end
            sym  = @(A) (A + A') / 2;
            poly = struct();
            poly.L2 = sym(N.' * blkdiag(b2{:}) * N);   % real in every case checked, kept
            poly.L1 = sym(N.' * blkdiag(b1{:}) * N);   % complex rather than asserted real
            poly.L0 = sym(N.' * blkdiag(b0{:}) * N);
            poly.M  = real(sym(N.' * blkdiag(bm{:}) * N));
            poly.N  = N;
        end
    end

    methods (Static, Access = private)
        function [z, w] = mirrorQuadrature(cps)
        %%% Gauss-Legendre nodes/weights on the lower half [0, 1/2] of the normalised
        %%% plate thickness, sub-intervals ending at every interface and mirrored interface
            t = [0, cumsum([cps.relDepth])] / sum([cps.relDepth]);
            b = sort([min(t, 1 - t), 0.5]);
            b = b([true, diff(b) > 1e-12]);              % merge coincident breakpoints
            n = max([cps.legOrder]);
            k = 1:n-1;
            beta = k ./ sqrt(4*k.^2 - 1);               % Golub-Welsch
            [V, D] = eig(diag(beta, 1) + diag(beta, -1));
            x  = diag(D);
            wx = 2 * V(1, :).'.^2;
            h  = diff(b);
            z  = b(1:end-1) + (x + 1)/2 .* h;            % [n x nSub]
            w  = wx/2 .* h;
            z  = z(:);
            w  = w(:);
        end

        function lay = layerOf(z, tB)
        %%% Layer index of each normalised depth z (tB = layer boundaries, bottom = 0)
            lay = sum(z(:) >= tB(2:end-1), 2) + 1;
        end

        function x = localCoord(z, tB, ilay)
        %%% Layer-local Legendre coordinate in [-1, 1] (bottom of the layer = -1)
            x = 2*(z(:) - tB(ilay)) / (tB(ilay+1) - tB(ilay)) - 1;
        end

        function P = normLegendre(L, x)
        %%% Normalised Legendre polynomials sqrt((2n+1)/2) P_n(x), n = 0..L-1, [numel(x) x L]
        %%% (the Ritz basis, as legendre(n, x, 'norm') of CalcParams.BasisDecomposition2D)
            x = x(:);
            P = ones(numel(x), L);
            if L > 1; P(:, 2) = x; end
            for n = 1:L-2
                P(:, n+2) = ((2*n + 1) * x .* P(:, n+1) - n * P(:, n)) / (n + 1);
            end
            P = P .* sqrt(((0:L-1) + 0.5));
        end
    end
end
