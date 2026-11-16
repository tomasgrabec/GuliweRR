function [phi1,Phi,phi2] = DCM2Euler(dcm)
    %%% Extract Bunge (ZXZ) Euler angles from a direction-cosine matrix
    %%% >> [phi1,Phi,phi2] = Material.DCM2Euler(dcm)
    %%%
    %%% INPUT:  dcm .. 3x3 direction-cosine matrix (rows normalised internally, handedness
    %%%                preserved as given); OR a 2x3/3x2 vector pair [n1;n2], right-hand
    %%%                completed via CompleteDirMat
    %%% OUTPUT: phi1,Phi,phi2 .. Bunge-ZXZ Euler angles [degrees], in [0,360) resp. [0,180]
    %%%
    %%% CONVENTION: g = Rz(phi2)*Rx(Phi)*Rz(phi1) transforms SPECIMEN-frame components into
    %%% CRYSTAL-frame components; dcm = g', so dcm's rows are the specimen axes in crystal
    %%% components. Row/column 3 is Bunge's distinguished ("tilt") axis.
    %%%
    %%% Raw Bunge-frame conversion: for a Material.dirMat (input coordinates, row 1 = surface
    %%% normal -- see docs/reference/architecture/core.md) use Material.DirMat2Euler instead.
    %%%
    %%% Do NOT call this on Material.GetStandardDirMat's output: its [3,2,1] reorder is an
    %%% odd permutation, so the result is improper (det = -1) and has no Bunge decomposition.
    %%% Angles still come back, but describe no physical rotation and will not round-trip.
    %%%
    %%% Same convention as MTEX (orientation.byEuler(...).matrix == dcm) and as the Nottingham
    %%% SRAS orientation output.

    if size(dcm,1) == 2 || size(dcm,2) == 2
        dcm = Material.CompleteDirMat(dcm);
    else
        for jj = 1:3
            dcm(jj,:) = dcm(jj,:) ./ norm(dcm(jj,:));
        end
    end

    tol = 1e-9;
    Phi = acosd(dcm(3,3));
    if abs(sind(Phi)) < tol
        %%% gimbal lock (Phi = 0 or 180 deg): only phi1±phi2 is determined, so by
        %%% convention the whole rotation goes into phi1
        phi2 = 0;
        phi1 = atan2d(dcm(2,1), dcm(1,1));
    else
        phi1 = atan2d(dcm(1,3), -dcm(2,3));
        phi2 = atan2d(dcm(3,1),  dcm(3,2));
    end

    phi1 = mod(phi1,360);
    phi2 = mod(phi2,360);
end
