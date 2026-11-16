function [vecOut,deTransfMat] = CubicDeconversion(vecIn)
    %%% De-Convert vector of cubic elastic constants (c11,c12,c44) from constants (B,c',c44),
    %%%  where B is bulk modulus, c' = (c11-c12)/2
    %%% >> [vecOut,deTransfMat] = CubicDeconversion(vecIn)
    %%%
    %%% INPUT:  (opt) vecIn .. constants (B,c',c44) (3x1 or 1x3); omit (or pass []) to get
    %%%                        the de-transformation matrix alone
    %%% OUTPUT: 1. vecOut      .. recovered constants (c11,c12,c44) as a 3x1 column
    %%%         2. deTransfMat .. inverse of CubicConversion's transformation matrix

    if nargin < 1; vecIn = []; end
    vecIn = ElastConsts.AsColumnVec(vecIn,3,'CubicDeconversion');

    [~,transfMat] = ElastConsts.CubicConversion();
    deTransfMat = inv(transfMat); %%% the inverse itself is the second output
    vecOut = transfMat \ vecIn;
end
