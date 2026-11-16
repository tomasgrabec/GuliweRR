function dirMat = CompleteDirMat(inp)
    %%% Build a right-hand orthonormal coordinate system from two direction vectors
    %%% >> dirMat = CompleteDirMat(inp)
    %%%
    %%% INPUT:  inp .. 2×3 or 3×2 matrix containing two row vectors [n1; n2]
    %%%                (a 3×3 is accepted too — only its first two rows are read, so
    %%%                 its third row is REBUILT right-handed, not preserved)
    %%% OUTPUT: dirMat .. 3×3 matrix (row vectors):
    %%%                     row 1 — n1 (first input vector, preserved)
    %%%                     row 2 — n2 (adjusted to be perpendicular to n1; warning
    %%%                             ElastConsts:RightHandFrame:notPerpendicular if it was not)
    %%%                     row 3 — n3 = cross(n1, n2) (right-hand completion)

    dirMat = ElastConsts.RightHandFrame(inp, ...
        'dirMat recalculated based on its first vector (input vectors not perpendicular).');
end
