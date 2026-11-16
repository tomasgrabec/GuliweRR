classdef GammaMatrix
% GammaMatrix — stiffness (K) and mass (M) matrices for the Ritz-Rayleigh solver.
%
%   Three-step static pipeline (Stoklasova 2015, Grabec 2024):
%     Step 1  rawByDepth    — [bVL x bVL x 15 x 3] tensor per layer, material-independent
%     Step 2  derivsByDepth — dK/dc and d²K/(dc dD) cells for angle-sweep path
%     Step 3  elastByDepth  — K0,K1,K2 per layer for Lamb frequency-sweep path
%   Instance method  calc      assembles K and M for an angle sweep (SAW path).
%   Instance method  calcLamb  assembles K and M for frequency-sweep (Lamb path).
%
%%% GuliweRR v6.0
%
%%% Note: not Christoffel's GammaMat / acoustic tensor Gamma(i,k)
%%% (guliwerr/solvers/christoffel/@Christoffel/matrix.m) — same name,
%%% different physics object.

    properties (SetAccess = protected)
        gammaMatAll     double      % K matrix (normalised, block-diagonal)
        dGammaDMatAll   double      % dK/dD matrix (for group velocity)
        normMatAll      double      % M matrix (density-normalised identity blocks)
    end

    methods
        function obj = GammaMatrix()
            %%% Initialise empty GammaMatrix object
            %%% >> obj = GammaMatrix()
        end

        function obj = calc(obj, structure, angleDeg, dKdC_cell, d2KdCdD_cell, redElastCell)
        %%% Assemble K and M for an angle sweep (SAW path)
        %%% >> obj = calc(obj, structure, angleDeg, dKdC_cell, d2KdCdD_cell, redElastCell)
        %%%
        %%% INPUT:  1. structure       .. Structure object
        %%%  (opt)  2. angleDeg        .. propagation angle [deg] (default = 0)
        %%%  (opt)  3. dKdC_cell       .. from derivsByDepth (computed if absent)
        %%%  (opt)  4. d2KdCdD_cell    .. from derivsByDepth (computed if absent)
        %%%  (opt)  5. redElastCell    .. reduced elastic constants per layer (computed if absent)
        %%%
        %%% The two derivative cells are supplied together or not at all; a supplied
        %%% cell must have one entry per layer (error GammaMatrix:calc:badPrecompute
        %%% otherwise, rather than a silent rebuild).

            if nargin < 3 || isempty(angleDeg);  angleDeg = 0;  end
            if nargin < 4; dKdC_cell = [];    end
            if nargin < 5; d2KdCdD_cell = []; end
            if nargin < 6; redElastCell = []; end

            %%% compute dKdC and d2KdCdD if not supplied
            if isempty(dKdC_cell) && isempty(d2KdCdD_cell)
                [dKdC_cell, d2KdCdD_cell] = GammaMatrix.derivsByDepth(structure);
            elseif isempty(dKdC_cell) || isempty(d2KdCdD_cell)
                error('GammaMatrix:calc:badPrecompute', ...
                    'dKdC_cell and d2KdCdD_cell must be supplied together (or both left empty).');
            else
                GammaMatrix.checkCellLength(dKdC_cell, structure, 'calc', 'dKdC_cell');
                GammaMatrix.checkCellLength(d2KdCdD_cell, structure, 'calc', 'd2KdCdD_cell');
            end

            %%% compute rotated reduced elastic constants per layer if not supplied
            if isempty(redElastCell)
                redElastCell = GammaMatrix.RedElastCellFor(structure, angleDeg);
            else
                GammaMatrix.checkCellLength(redElastCell, structure, 'calc', 'redElastCell');
            end

            %%% assemble block-diagonal K, dK/dD, M: one block per layer, joined by a
            %%% single blkdiag (growing the matrix layer by layer re-copies it each time)
            nLay      = structure.nLay;
            gBlk      = cell(1, nLay);
            gBlk_dD   = cell(1, nLay);
            nBlk      = cell(1, nLay);

            for ilay = 1:nLay
                bVL      = structure.calcParams(ilay).basisVecL;
                redElast = redElastCell{ilay}(:);
                normFac  = GammaMatrix.LayerNormFac(structure, ilay);

                %%% K = sum_c c dK/dc, as one matrix-vector product over the 15 constants
                gBlk{ilay}    = complex(reshape(reshape(dKdC_cell{ilay},    [], 15) * redElast, bVL, bVL)) .* normFac;
                gBlk_dD{ilay} = complex(reshape(reshape(d2KdCdD_cell{ilay}, [], 15) * redElast, bVL, bVL)) .* normFac;
                nBlk{ilay}    = eye(bVL) .* normFac;
            end

            gMatAll       = blkdiag(gBlk{:});
            gMatAll_dD    = blkdiag(gBlk_dD{:});
            nMatAll       = blkdiag(nBlk{:});

            obj.gammaMatAll       = gMatAll;
            obj.dGammaDMatAll     = gMatAll_dD;
            obj.normMatAll        = nMatAll;
        end

        function obj = calcLamb(obj, structure, angleDeg, gammaByDepthCell)
        %%% Assemble K and M for frequency-sweep (Lamb plate path)
        %%% >> obj = calcLamb(obj, structure, angleDeg, gammaByDepthCell)
        %%%
        %%% INPUT:  1. structure          .. Structure object
        %%%  (opt)  2. angleDeg           .. propagation angle [deg] (default = 0)
        %%%  (opt)  3. gammaByDepthCell   .. from elastByDepth (computed if absent)

            if nargin < 3 || isempty(angleDeg); angleDeg = 0; end

            if nargin < 4 || isempty(gammaByDepthCell)
                gammaByDepthCell = GammaMatrix.elastByDepth(structure, angleDeg);
            end

            nLay    = structure.nLay;
            gBlk    = cell(1, nLay);
            gBlk_dD = cell(1, nLay);
            nBlk    = cell(1, nLay);

            for ilay = 1:nLay
                gByD = gammaByDepthCell{ilay};
                [gamma_woRho, dGamma_dD_woRho] = GammaMatrix.ContractDepth( ...
                    gByD(:,:,1), gByD(:,:,2), gByD(:,:,3), structure.calcParams(ilay).relDepth);

                gamma_lay     = gamma_woRho     ./ structure.mater(ilay).rho;
                dGamma_dD_lay = dGamma_dD_woRho ./ structure.mater(ilay).rho;

                bVL     = structure.calcParams(ilay).basisVecL;
                normFac = GammaMatrix.LayerNormFac(structure, ilay);

                gBlk{ilay}    = gamma_lay     .* normFac;
                gBlk_dD{ilay} = dGamma_dD_lay .* normFac;
                nBlk{ilay}    = eye(bVL)      .* normFac;
            end

            gMatAll    = blkdiag(gBlk{:});
            gMatAll_dD = blkdiag(gBlk_dD{:});
            nMatAll    = blkdiag(nBlk{:});

            obj.gammaMatAll   = gMatAll;
            obj.dGammaDMatAll = gMatAll_dD;
            obj.normMatAll    = nMatAll;
        end
    end

    methods (Static)
        rawTensor = rawByDepth(calcParams)
        [dKdC_cell, d2KdCdD_cell] = derivsByDepth(structure, rawByDepthCell)
        gammaByDepthCell = elastByDepth(structure, angleDeg, rawByDepthCell, redElastCell)
        [dG_str_dummy, dG_lay] = dummyForLayer(structure, noLay, dKdC_cell, angleDeg)
        [dG_blk, rows, dG_lay] = layerBlockDerivs(structure, noLay, dKdC_cell, angleDeg)

        function rawCell = RawCellFor(structure)
        %%% Step-1 Legendre-integral tensor for every layer of a structure
        %%% >> rawCell = GammaMatrix.RawCellFor(structure)
        %%%
        %%% The default-fill used by derivsByDepth and elastByDepth when no
        %%% pre-built cell is handed in. Public so that callers sweeping a
        %%% structure at many wavelengths/angles can build it once and pass it
        %%% back down (CalcWavevecAngSweep2Layered does this).
            rawCell = cell(structure.nLay, 1);
            for ilay = 1:structure.nLay
                rawCell{ilay} = GammaMatrix.rawByDepth(structure.calcParams(ilay));
            end
        end

        function rotMat = LayerRotMat(structure, ilay, angleDeg)
        %%% Rotation matrix taking layer ilay's crystal frame to the standard
        %%% solver frame (x1 = propagation, x2 = in-plane, x3 = depth) at angleDeg
        %%% >> rotMat = GammaMatrix.LayerRotMat(structure, ilay, angleDeg)
        %%%
        %%% GetStandardDirMat applies the propagation angle within the layer's
        %%% own surface plane and reorders to the standard convention, so
        %%% CalcRotMat takes angle 0 — passing angleDeg to BOTH double-rotates.
        %%% See docs/pitfalls.md OTH1.
            stdDirMat = Material.GetStandardDirMat(structure.mater(ilay).dirMat, angleDeg);
            rotMat    = Material.CalcRotMat(stdDirMat, 0);
        end

        function redElastCell = RedElastCellFor(structure, angleDeg)
        %%% Rotated reduced (15-of-21) elastic constants for every layer
        %%% >> redElastCell = GammaMatrix.RedElastCellFor(structure, angleDeg)
        %%%
        %%% The default-fill shared by GammaMatrix.calc and elastByDepth.
            redElastCell = cell(structure.nLay, 1);
            for ilay = 1:structure.nLay
                rotMat = GammaMatrix.LayerRotMat(structure, ilay, angleDeg);
                redElastCell{ilay} = structure.mater(ilay).elConsts.GetReducedElastConsts(rotMat);
            end
        end

        function normFac = LayerNormFac(structure, ilay)
        %%% Per-layer normalisation factor of the block-diagonal K/M assembly
        %%% >> normFac = GammaMatrix.LayerNormFac(structure, ilay)
        %%%
        %%% (relDepth_ilay * rho_ilay) / (relDepth_1 * rho_1) — layer 1 (the
        %%% substrate / bottom layer) is the reference, so its factor is 1.
            normFac = (structure.calcParams(ilay).relDepth * structure.mater(ilay).rho) ...
                    / (structure.calcParams(1).relDepth    * structure.mater(1).rho);
        end

        function [K, dKdD] = ContractDepth(K0, K1, K2, relDepth)
        %%% Apply the depth factors to a (K0, K1, K2) triple
        %%% >> [K, dKdD] = GammaMatrix.ContractDepth(K0, K1, K2, relDepth)
        %%%
        %%%   K     = K0/dd1^2 + K1/(dd1*dd2/2) + K2/(dd2/2)^2
        %%%   dK/dD = -2*K1/(dd1*dd2^2) - 8*K2/dd2^3
        %%% with dd1 = 1 permanently and dd2 = relDepth (Stoklasova 2015,
        %%% Eq. (A1)-(A4)). Used by both the SAW/angle-sweep path
        %%% (derivsByDepth) and the Lamb/frequency-sweep path (calcLamb); they
        %%% differ in what the arrays hold, not in the depth convention.
            K    = K0 ...
                 + K1 ./ (relDepth/2) ...
                 + K2 ./ (relDepth/2)^2;
            dKdD = -2 * K1 ./ relDepth^2 ...
                   -8 * K2 ./ relDepth^3;
        end
    end

    methods (Static, Access = private)
        function checkCellLength(c, structure, methodName, argName)
        %%% Error GammaMatrix:<methodName>:badPrecompute unless c has one entry per layer
        %%% >> GammaMatrix.checkCellLength(c, structure, methodName, argName)
            if ~iscell(c) || numel(c) ~= structure.nLay
                error(['GammaMatrix:' methodName ':badPrecompute'], ...
                    '%s must be a cell with one entry per layer (%d); got %d.', ...
                    argName, structure.nLay, numel(c));
            end
        end
    end
end
