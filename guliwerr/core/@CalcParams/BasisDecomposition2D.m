function pLeg = BasisDecomposition2D(legO,meshSpacing)
    %%% Sample the Legendre basis functions on a uniform mesh of [-1, 1]
    %%% >> pLeg = BasisDecomposition2D(legO,meshSpacing)
    %%%
    %%% INPUT:  1. legO        .. max order of the Legendre polynomial
    %%%  (opt)  2. meshSpacing .. spacing of the spatial discretisation (default = 0.02)
    %%% OUTPUT: pLeg .. normalised Legendre-function values (legO × nMesh); row k holds
    %%%                 degree k-1. Used for modal-shape reconstruction from RR eigenvectors.

    if nargin < 2 || isempty(meshSpacing)
        meshSpacing = 0.02;
    end

    ti = -1 : meshSpacing : 1;

    pLeg = zeros(legO,length(ti));
    for iLeg = 1:legO
        legVal = legendre(iLeg-1,ti','norm');
        pLeg(iLeg,:) = legVal(1,:);
    end
end
