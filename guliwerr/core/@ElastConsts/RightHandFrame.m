function frame = RightHandFrame(twoVecs,noteMsg)
    %%% Build a right-hand orthonormal 3x3 frame (row vectors) from two direction vectors
    %%% >> frame = ElastConsts.RightHandFrame(twoVecs,noteMsg)
    %%%
    %%% INPUT:  1. twoVecs .. 2x3 or 3x2 matrix of two direction vectors; a 3x3 is also
    %%%                       accepted and only its first two ROWS are read (its third
    %%%                       row is recomputed, so handedness is imposed, not preserved)
    %%%  (opt)  2. noteMsg .. text of the warning ElastConsts:RightHandFrame:notPerpendicular,
    %%%                       raised when the two input vectors are not perpendicular
    %%%                       (|n1.n2| > 1e-3, i.e. when the second one gets adjusted)
    %%% OUTPUT: frame .. 3x3 matrix of row vectors:
    %%%                    row 1 -- n1, the normalised FIRST input vector (preserved)
    %%%                    row 2 -- n2, the second vector made exactly perpendicular to n1
    %%%                    row 3 -- n3 = cross(n1,n2), right-hand completion
    %%%
    %%% ROLE-AGNOSTIC by design: callers disagree on what the two vectors mean
    %%% (Material.CompleteDirMat passes [surface normal; propagation],
    %%% ElastConsts.RotateTensor passes [propagation; surface normal]), so do not add
    %%% role-specific behaviour here.

    if nargin < 2 || isempty(noteMsg)
        noteMsg = 'Coordinate system recalculated based on the first vector (input vectors not perpendicular).';
    end

    %%% accept both 2x3 and 3x2; a 3x3 is left alone (rows 1-2 used)
    if size(twoVecs,1) > size(twoVecs,2) && size(twoVecs,1) == 3
        twoVecs = twoVecs';
    end
    if size(twoVecs,1) < 2 || size(twoVecs,2) ~= 3
        error('ElastConsts:RightHandFrame:badInput', ...
            'Input must hold two direction vectors as rows (2x3, 3x2 or 3x3); got %dx%d.', ...
            size(twoVecs,1), size(twoVecs,2));
    end

    n1 = twoVecs(1,:) ./ norm(twoVecs(1,:));
    n2 = twoVecs(2,:) ./ norm(twoVecs(2,:));

    if abs(dot(n1,n2)) > 1e-3   % inputs not perpendicular: n2 is about to be adjusted
        warning('ElastConsts:RightHandFrame:notPerpendicular', '%s', noteMsg);
    end

    n3 = cross(n1,n2);
    n2 = cross(n3,n1);   % make n2 exactly perpendicular to n1

    frame = [n1; n2; n3];
end
