function rotDirMat = RotDirMatFromAngle(dirMat,angleDeg)
    %%% Rotate a direction matrix in-plane by the given angle
    %%% >> rotDirMat = RotDirMatFromAngle(dirMat,angleDeg)
    %%%
    %%% INPUT:  1. dirMat   .. 3×3 direction matrix (row vectors; complete right-hand system)
    %%%         2. angleDeg .. in-plane rotation angle [degrees]
    %%% OUTPUT: rotDirMat .. 3×3 direction matrix, rows 2 and 3 rotated by angleDeg from row 2
    %%%                      towards row 3 (counter-clockwise seen looking down the surface
    %%%                      normal, row 1, onto the surface); with Euler2DirMat this is from
    %%%                      specimen X towards Y. Row 1 is the pivot and stays unchanged.
    %%%                      Sole definition of the in-plane rotation's sign convention.

    rotDirMat = [ dirMat(1,:); ...
                  dirMat(2,:) .* cosd(angleDeg)    + dirMat(3,:) .* sind(angleDeg); ...
                  dirMat(2,:) .* (-sind(angleDeg)) + dirMat(3,:) .* cosd(angleDeg) ];
end
