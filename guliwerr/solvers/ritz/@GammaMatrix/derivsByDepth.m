function [dKdC_cell, d2KdCdD_cell] = derivsByDepth(structure, rawByDepthCell)
%%% Step 2: Contract depth factors to get dK/dc and d²K/(dc dD) per layer
%%% >> [dKdC_cell, d2KdCdD_cell] = GammaMatrix.derivsByDepth(structure, rawByDepthCell)
%%%
%%% INPUT:  1. structure        .. Structure object
%%%  (opt)  2. rawByDepthCell   .. cell of rawByDepth tensors per layer
%%% OUTPUT: dKdC_cell{ilay}     .. [bVL x bVL x 15] for K = sum(dK/dc * c)
%%%         d2KdCdD_cell{ilay}  .. [bVL x bVL x 15] for dK/dD = sum(d²K/(dc dD) * c)

    if nargin < 2 || isempty(rawByDepthCell)
        rawByDepthCell = GammaMatrix.RawCellFor(structure);
    end

    dKdC_cell    = cell(structure.nLay, 1);
    d2KdCdD_cell = cell(structure.nLay, 1);

    for ilay = 1:structure.nLay
        raw = rawByDepthCell{ilay};   % [bVL x bVL x 15 x 3]

        %%% depth factors: K0/K1/K2 -> K, dK/dD
        [dKdC, d2KdCdD] = GammaMatrix.ContractDepth( ...
            raw(:,:,:,1), raw(:,:,:,2), raw(:,:,:,3), structure.calcParams(ilay).relDepth);

        dKdC_cell{ilay}    = dKdC    ./ structure.mater(ilay).rho;
        d2KdCdD_cell{ilay} = d2KdCdD ./ structure.mater(ilay).rho;
    end
end
