function [vecOut,deTransfMat] = TetragonalDeconversion(vecIn)
    %%% De-Convert vector of tetragonal elastic constants (c11,c12,c13,c33,c44,c66) from
    %%%  constants (B,c33,c',c'',c44,c66) -- see TetragonalConversion for the definitions
    %%% >> [vecOut,deTransfMat] = TetragonalDeconversion(vecIn)
    %%%
    %%% INPUT:  (opt) vecIn .. constants (B,c33,c',c'',c44,c66) (6x1 or 1x6); omit (or pass
    %%%                        []) to get the de-transformation matrix alone
    %%% OUTPUT: 1. vecOut      .. recovered constants (c11,c12,c13,c33,c44,c66) as a 6x1 column
    %%%         2. deTransfMat .. inverse of TetragonalConversion's transformation matrix

    if nargin < 1; vecIn = []; end
    vecIn = ElastConsts.AsColumnVec(vecIn,6,'TetragonalDeconversion');

    [~,transfMat] = ElastConsts.TetragonalConversion();
    deTransfMat = inv(transfMat); %%% the inverse itself is the second output
    vecOut = transfMat \ vecIn;
end
