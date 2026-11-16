function dcm = Euler2DCM(phi1,Phi,phi2)
    %%% Build a direction-cosine matrix from Bunge (ZXZ) Euler angles
    %%% >> dcm = Material.Euler2DCM(phi1,Phi,phi2)
    %%%
    %%% INPUT:  phi1,Phi,phi2 .. Bunge-ZXZ Euler angles [degrees]
    %%% OUTPUT: dcm .. 3x3 direction-cosine matrix, exact inverse of Material.DCM2Euler
    %%%                (see there for the convention)
    %%%
    %%% Raw Bunge frame: row/column 3 is Bunge's tilt axis, NOT GULIWERR's surface normal
    %%% (row 1 of a Convention-A dirMat -- see docs/reference/architecture/core.md). For a
    %%% ready-to-use Material.dirMat, use Material.Euler2DirMat.

    c1 = cosd(phi1); s1 = sind(phi1);
    C  = cosd(Phi);  S  = sind(Phi);
    c2 = cosd(phi2); s2 = sind(phi2);

    g = [ c1*c2 - s1*s2*C,   s1*c2 + c1*s2*C,   s2*S; ...
         -c1*s2 - s1*c2*C,  -s1*s2 + c1*c2*C,   c2*S; ...
          s1*S,             -c1*S,              C   ];

    dcm = g';
end
