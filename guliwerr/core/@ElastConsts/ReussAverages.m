function [c11_reuss,c44_reuss] = ReussAverages(elVoigt)
    %%% Reuss averages -- moduli of the equivalent statistically isotropic material under
    %%%  uniform STRESS (Reuss 1929), the lower bound (VoigtAverages gives the upper)
    %%% >> [c11_reuss,c44_reuss] = ReussAverages(elVoigt)
    %%%
    %%% INPUT:  elVoigt .. Voigt matrix of elastic constants (6x6)
    %%% OUTPUT: c11_reuss .. M, P-wave (longitudinal) modulus
    %%%         c44_reuss .. G, S-wave (transverse) modulus
    %%%
    %%% General-symmetry Reuss K and G from the compliance S = inv(elVoigt), combined as
    %%% M = K + 4/3*G; K_reuss equals Voigt2Bulk's bulk modulus.

    s = inv(elVoigt);   %%% compliance matrix; its blocks give the Reuss averages

    K_reuss = 1 / (s(1,1)+s(2,2)+s(3,3) + 2*(s(1,2)+s(1,3)+s(2,3)));
    G_reuss = 15 / (4*(s(1,1)+s(2,2)+s(3,3)) - 4*(s(1,2)+s(1,3)+s(2,3)) ...
        + 3*(s(4,4)+s(5,5)+s(6,6)));

    c11_reuss = K_reuss + 4/3*G_reuss;
    c44_reuss = G_reuss;
end
