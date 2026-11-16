function rotMat = CalcRotMat(dirMat,angleDeg)
    %%% Build a rotation matrix (3×2) from a direction matrix and an angle
    %%% >> rotMat = CalcRotMat(dirMat,angleDeg)
    %%%
    %%% INPUT:  1. dirMat   .. 3×3 direction matrix (row vectors; complete right-hand system)
    %%%  (opt)  2. angleDeg .. in-plane rotation angle [degrees] (default = 0)
    %%% OUTPUT: rotMat .. 3×2 matrix for ElastConsts.RotateTensor:
    %%%             column 1 — propagation direction, rotated within the surface plane
    %%%                        (the plane spanned by dirMat rows 2 and 3)
    %%%             column 2 — surface normal (dirMat row 1), the pivot; unchanged

    if nargin < 2 || isempty(angleDeg)
        angleDeg = 0;
    end

    rotDirMat = Material.RotDirMatFromAngle(dirMat,angleDeg);
    rotMat    = [ rotDirMat(2,:)', rotDirMat(1,:)' ];
end
