function [vecOut,transfMat] = HexagonalConversion(vecIn)
    %%% Convert hexagonal elastic constants (c11,c12,c13,c33,c44) to the physical-mode
    %%%  basis (C11,C33,G1,G2,G3): uniaxial elongations C11 (basal plane) and C33 (c-axis),
    %%%  and shear moduli G1 = c44 (axial shear), G2 = (c11-c12)/2 (basal shear) and
    %%%  G3 = (c11+c33-2*c13)/4 (breathing shear, volume-preserving c/a distortion).
    %%%  G1 = G2 = G3 in the isotropic limit.
    %%% >> [vecOut,transfMat] = HexagonalConversion(vecIn)
    %%%
    %%% INPUT:  (opt) vecIn .. hexagonal elastic constants [c11,c12,c13,c33,c44]
    %%%                        (5x1 or 1x5); omit (or pass []) to get the
    %%%                        transformation matrix alone
    %%% OUTPUT: 1. vecOut    .. transformed constants (C11,C33,G1,G2,G3) as a 5x1 column
    %%%         2. transfMat .. transformation matrix, vecOut(5x1) = transfMat(5x5)*vecIn(5x1)

    if nargin < 1; vecIn = []; end
    vecIn = ElastConsts.AsColumnVec(vecIn,5,'HexagonalConversion');

    transfMat = [...
        1    0    0    0   0;   ... % C11
        0    0    0    1   0;   ... % C33
        0    0    0    0   1;   ... % G1 == axial shear
        1/2  -1/2 0    0   0;   ... % G2 == basal shear
        1/4  0   -1/2  1/4 0;   ... % G3 == breathing shear
          ];
    vecOut = transfMat * vecIn;
end
