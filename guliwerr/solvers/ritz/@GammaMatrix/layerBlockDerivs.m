function [dG_blk, rows, dG_lay] = layerBlockDerivs(structure, noLay, dKdC_cell, angleDeg)
%%% Inverse support: dK/dc of ONE layer as its own diagonal block (no zero padding)
%%% >> [dG_blk, rows, dG_lay] = GammaMatrix.layerBlockDerivs(structure, noLay, dKdC_cell, angleDeg)
%%%
%%% INPUT:  1. structure  .. Structure object
%%%  (opt)  2. noLay      .. layer index (default = top layer)
%%%  (opt)  3. dKdC_cell  .. from derivsByDepth (computed if absent)
%%%  (opt)  4. angleDeg   .. propagation angle (default = 0)
%%% OUTPUT: dG_blk .. [bVL_lay x bVL_lay x nIndep] dK/dc per independent constant, WITH the
%%%                   layer's normalisation factor: exactly dummyForLayer's noLay block
%%%         rows   .. rows/columns of that block in the full-structure basis
%%%         dG_lay .. the same without the normalisation factor
%%%
%%% dK/dc is non-zero only on this block (K is block-diagonal; the layers couple only
%%% through the constraint null space N), so a caller that wants N'(dK/dc)N or a quadratic
%%% form can use N(rows,:) and skip the zero blocks -- dummyForLayer's full-size array is
%%% strBVL^2 * nIndep complex numbers (120 MB for a 600-function stack).

    if nargin < 2 || isempty(noLay);     noLay    = structure.nLay;  end
    if nargin < 3 || isempty(dKdC_cell); [dKdC_cell, ~] = GammaMatrix.derivsByDepth(structure); end
    if nargin < 4 || isempty(angleDeg);  angleDeg = 0;               end

    nIndep       = structure.mater(noLay).elConsts.nIndepConsts;
    bVL_lay      = structure.calcParams(noLay).basisVecL;
    rotMat       = GammaMatrix.LayerRotMat(structure, noLay, angleDeg);
    indep2RedMat = ElastConsts.ElastRotReductMat(rotMat, structure.mater(noLay).elConsts.symmClassNo);
    normFac      = GammaMatrix.LayerNormFac(structure, noLay);

    %%% dK/dc per independent constant for layer noLay
    dG_lay = complex(zeros(bVL_lay, bVL_lay, nIndep));
    for iI = 1:nIndep
        for iR = 1:15
            dG_lay(:,:,iI) = dG_lay(:,:,iI) + dKdC_cell{noLay}(:,:,iR) .* indep2RedMat(iR,iI);
        end
    end

    dG_blk    = dG_lay .* normFac;
    bVL_below = sum([structure.calcParams(1:noLay-1).basisVecL]);
    rows      = bVL_below+1 : bVL_below+bVL_lay;
end
