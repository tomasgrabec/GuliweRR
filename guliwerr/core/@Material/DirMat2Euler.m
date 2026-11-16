function [phi1,Phi,phi2] = DirMat2Euler(dirMat)
    %%% Extract Bunge (ZXZ) Euler angles from a Material.dirMat orientation
    %%% >> [phi1,Phi,phi2] = Material.DirMat2Euler(dirMat)
    %%%
    %%% INPUT:  dirMat .. 3x3 Convention-A direction matrix (row 1 = surface normal); OR a
    %%%                   2x3/3x2 vector pair [n1;n2], right-hand completed via CompleteDirMat
    %%% OUTPUT: phi1,Phi,phi2 .. Bunge-ZXZ Euler angles [degrees], in [0,360) resp. [0,180],
    %%%                          exact inverse of Material.Euler2DirMat
    %%%
    %%% Right-hand completes the input (a no-op for any genuine Material.dirMat), reorders it
    %%% [2,3,1] onto Bunge's row-3 tilt axis and calls Material.DCM2Euler -- call that
    %%% directly for a dcm already in Bunge's own row order.
    %%%
    %%% Same Bunge convention as MTEX and the Nottingham SRAS orientation output.

    dirMat = Material.CompleteDirMat(dirMat);
    [phi1,Phi,phi2] = Material.DCM2Euler(dirMat([2,3,1],:));
end
