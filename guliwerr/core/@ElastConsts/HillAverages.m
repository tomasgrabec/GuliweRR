function [c11_hill,c44_hill] = HillAverages(elVoigt)
    %%% Hill averages -- arithmetic mean of the Voigt and Reuss bounds (Hill 1952), the
    %%%  standard practical estimate of a polycrystal's isotropic elastic moduli
    %%% >> [c11_hill,c44_hill] = HillAverages(elVoigt)
    %%%
    %%% INPUT:  elVoigt .. Voigt matrix of elastic constants (6x6)
    %%% OUTPUT: c11_hill .. M, P-wave (longitudinal) modulus
    %%%         c44_hill .. G, S-wave (transverse) modulus

    [c11_voigt,c44_voigt] = ElastConsts.VoigtAverages(elVoigt);
    [c11_reuss,c44_reuss] = ElastConsts.ReussAverages(elVoigt);

    c11_hill = (c11_voigt + c11_reuss) / 2;
    c44_hill = (c44_voigt + c44_reuss) / 2;
end
