function [vecOut,transfMat] = CubicConversion(vecIn)
    %%% Convert vector of cubic elastic constants (c11,c12,c44) to constants (B,c',c44),
    %%%  where B is bulk modulus, c' = (c11-c12)/2
    %%% >> [vecOut,transfMat] = CubicConversion(vecIn)
    %%%
    %%% INPUT:  (opt) vecIn .. cubic elastic constants [c11,c12,c44] (3x1 or 1x3);
    %%%                        omit (or pass []) to get the transformation matrix alone
    %%% OUTPUT: 1. vecOut    .. transformed constants (B,c',c44) as a 3x1 column
    %%%         2. transfMat .. transformation matrix, vecOut(3x1) = transfMat(3x3)*vecIn(3x1)

    if nargin < 1; vecIn = []; end
    vecIn = ElastConsts.AsColumnVec(vecIn,3,'CubicConversion');

    transfMat = [1/3 2/3 0; ... % B
                  .5 -.5 0; ... % c'
                  0   0  1];    % c44
    vecOut = transfMat * vecIn;
end
