function gammaByDepthCell = elastByDepth(structure, angleDeg, rawByDepthCell, redElastCell)
%%% Step 3 (Lamb path): contract rawByDepth with elasticity, keep depth factors separate
%%% >> gammaByDepthCell = GammaMatrix.elastByDepth(structure, angleDeg, rawByDepthCell, redElastCell)
%%%
%%% OUTPUT: gammaByDepthCell{ilay} .. [bVL x bVL x 3] (K0, K1, K2 with elasticity, without rho/depth)
%%%
%%% A supplied rawByDepthCell / redElastCell must have one entry per layer (error
%%% GammaMatrix:elastByDepth:badPrecompute otherwise); empty = compute.

    if nargin < 2 || isempty(angleDeg);    angleDeg = 0;    end
    if nargin < 3 || isempty(rawByDepthCell)
        rawByDepthCell = GammaMatrix.RawCellFor(structure);
    else
        GammaMatrix.checkCellLength(rawByDepthCell, structure, 'elastByDepth', 'rawByDepthCell');
    end
    if nargin < 4 || isempty(redElastCell)
        redElastCell = GammaMatrix.RedElastCellFor(structure, angleDeg);
    else
        GammaMatrix.checkCellLength(redElastCell, structure, 'elastByDepth', 'redElastCell');
    end

    gammaByDepthCell = cell(structure.nLay, 1);
    for ilay = 1:structure.nLay
        bVL = structure.calcParams(ilay).basisVecL;
        gByD = complex(zeros(bVL, bVL, 3));
        for iC = 1:15
            %%% squeeze() is load-bearing: rawByDepthCell{ilay}(:,:,iC,:) is
            %%% [bVL x bVL x 1 x 3], and adding that straight to gByD
            %%% ([bVL x bVL x 3]) does NOT error — implicit expansion silently
            %%% broadcasts both mismatched trailing dims, growing gByD to
            %%% [bVL x bVL x 3 x 3] with all three K0/K1/K2 pages collapsed onto
            %%% K0. squeeze() drops the singleton so the add is a genuine
            %%% element-wise K0/K1/K2 accumulation.
            gByD = gByD + squeeze(rawByDepthCell{ilay}(:,:,iC,:)) .* redElastCell{ilay}(iC);
        end
        gammaByDepthCell{ilay} = gByD;
    end
end
