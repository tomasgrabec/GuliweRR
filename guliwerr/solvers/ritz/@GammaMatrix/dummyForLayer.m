function [dG_str_dummy, dG_lay] = dummyForLayer(structure, noLay, dKdC_cell, angleDeg)
%%% Inverse support: dK/dc with zeros except for layer noLay  (sensitivity matrix input)
%%% >> [dG_str_dummy, dG_lay] = GammaMatrix.dummyForLayer(structure, noLay, dKdC_cell, angleDeg)
%%%
%%% INPUT:  1. structure  .. Structure object
%%%  (opt)  2. noLay      .. layer index (default = top layer)
%%%  (opt)  3. dKdC_cell  .. from derivsByDepth (computed if absent)
%%%  (opt)  4. angleDeg   .. propagation angle (default = 0)
%%% OUTPUT: dG_str_dummy  .. [strBVL x strBVL x nIndep] zeros except noLay block
%%%         dG_lay        .. [bVL_lay x bVL_lay x nIndep] dK/dc for layer noLay

    if nargin < 2 || isempty(noLay);     noLay    = structure.nLay;  end
    if nargin < 3 || isempty(dKdC_cell); [dKdC_cell, ~] = GammaMatrix.derivsByDepth(structure); end
    if nargin < 4 || isempty(angleDeg);  angleDeg = 0;               end

    [dG_blk, rows, dG_lay] = GammaMatrix.layerBlockDerivs(structure, noLay, dKdC_cell, angleDeg);
    nIndep = size(dG_blk, 3);

    %%% build full-structure dummy: zeros except for the noLay block
    dG_str_dummy = complex(zeros(structure.strBasisVecL, structure.strBasisVecL, nIndep));
    dG_str_dummy(rows, rows, :) = dG_blk;
end
