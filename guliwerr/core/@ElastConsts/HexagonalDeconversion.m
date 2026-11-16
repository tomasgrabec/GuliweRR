function [vecOut,deTransfMat] = HexagonalDeconversion(vecIn)
    %%% De-Convert vector of hexagonal elastic constants (c11,c12,c13,c33,c44) from the
    %%%  physical-mode basis (C11,C33,G1,G2,G3) defined in HexagonalConversion
    %%% >> [vecOut,deTransfMat] = HexagonalDeconversion(vecIn)
    %%%
    %%% INPUT:  (opt) vecIn .. constants (C11,C33,G1,G2,G3) (5x1 or 1x5); omit (or pass [])
    %%%                        to get the de-transformation matrix alone
    %%% OUTPUT: 1. vecOut      .. recovered constants (c11,c12,c13,c33,c44) as a 5x1 column
    %%%         2. deTransfMat .. inverse of HexagonalConversion's transformation matrix

    if nargin < 1; vecIn = []; end
    vecIn = ElastConsts.AsColumnVec(vecIn,5,'HexagonalDeconversion');

    [~,transfMat] = ElastConsts.HexagonalConversion();
    deTransfMat = inv(transfMat); %%% the inverse itself is the second output
    vecOut = transfMat \ vecIn;
end
