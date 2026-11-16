function [vecOut,transfMat] = TetragonalConversion(vecIn)
    %%% Calculate equivalent tetragonal constants B, c33, c', c'', c44, c66
    %%% >> [vecOut,transfMat] = TetragonalConversion(vecIn)
    %%%
    %%% INPUT:  (opt) vecIn .. tetragonal elastic constants (6x1 or 1x6); omit (or pass [])
    %%%                        to get the transformation matrix alone
    %%% OUTPUT: 1. vecOut    .. transformed constants (B,c33,c',c'',c44,c66) as a 6x1 column
    %%%         2. transfMat .. transformation matrix, vecOut(6x1) = transfMat(6x6)*vecIn(6x1)

    if nargin < 1; vecIn = []; end
    vecIn = ElastConsts.AsColumnVec(vecIn,6,'TetragonalConversion');

    transfMat = [...
        2/9  2/9   4/9   1/9     0   0;   ... % B
         0    0     0     1      0   0;   ... % c33
        1/2 -1/2    0     0      0   0;   ... % c'
        1/6  1/6  -4/6   2/6     0   0;   ... % c'' == volume-conserving 'breathing' mode
         0    0     0     0      1   0;   ... % c44
         0    0     0     0      0   1;   ... % c66
          ];
    vecOut = transfMat * vecIn;
end
