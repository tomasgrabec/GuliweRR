function dirMat = Euler2DirMat(phi1,Phi,phi2)
    %%% Build a Material.dirMat orientation from Bunge (ZXZ) Euler angles
    %%% >> dirMat = Material.Euler2DirMat(phi1,Phi,phi2)
    %%%
    %%% INPUT:  phi1,Phi,phi2 .. Bunge-ZXZ Euler angles [degrees]
    %%% OUTPUT: dirMat .. 3x3 Convention-A direction matrix (row 1 = surface normal, row 2 =
    %%%                   in-plane propagation reference, row 3 = completion), exact inverse
    %%%                   of Material.DirMat2Euler
    %%%
    %%% Reorders Material.Euler2DCM's output [3,1,2], moving Bunge's row-3 tilt axis onto
    %%% row 1 (surface normal) -- call Euler2DCM directly for Bunge's own row order.

    dcm = Material.Euler2DCM(phi1,Phi,phi2);
    dirMat = dcm([3,1,2],:);
end
