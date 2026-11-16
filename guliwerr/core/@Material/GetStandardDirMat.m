function stdDirMat = GetStandardDirMat(dirMat, angleDeg)
    %%% Reorder (and optionally rotate) a direction matrix into the standard
    %%% x1=propagation, x2=in-plane, x3=surface-normal convention of the literature
    %%% (Every, Maznev), from the internal row order [plane; direction; completion]
    %%% >> stdDirMat = Material.GetStandardDirMat(dirMat, angleDeg)
    %%%
    %%% INPUT:  1. dirMat   .. 3×3 direction matrix (rows: [plane; direction; completion])
    %%%  (opt)  2. angleDeg .. in-plane rotation angle [degrees] (default = 0)
    %%% OUTPUT: stdDirMat .. 3×3 direction matrix, rows [completion; direction; plane]
    %%%
    %%% The angle is applied here, before the reorder, because RotDirMatFromAngle mixes
    %%% rows 2 & 3 — rotating an already-reordered matrix would swing the propagation
    %%% direction towards the surface normal. So pass angle 0 downstream, e.g.
    %%% Material.CalcRotMat(stdDirMat, 0).

    if nargin < 2 || isempty(angleDeg)
        angleDeg = 0;
    end
    rotDirMat = Material.RotDirMatFromAngle(dirMat, angleDeg);
    stdDirMat = rotDirMat([3,2,1], :);
end
