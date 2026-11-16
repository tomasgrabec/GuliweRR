function [c11_voigt,c44_voigt] = VoigtAverages(elVoigt)
    %%% Voigt averages — the moduli of the equivalent statistically isotropic material
    %%% >> [c11_voigt,c44_voigt] = VoigtAverages(elVoigt)
    %%%
    %%% INPUT:  elVoigt .. Voigt matrix of elastic constants (6x6)
    %%% OUTPUT: c11_voigt .. M, P-wave (longitudinal) modulus
    %%%         c44_voigt .. G, S-wave (transverse) modulus
    %%% Kube & Turner 2015, Eq. (11), with its typographical errors corrected.

    c = elVoigt;

    c11_voigt=1/15*(3*c(1,1)+2*c(1,2)+2*c(1,3)+3*c(2,2)+2*c(2,3)+3*c(3,3) ...
        +4*c(4,4)+4*c(5,5)+4*c(6,6));  %c033 in Turner;
    c44_voigt=1/15*(c(1,1)-c(1,2)-c(1,3)+c(2,2)-c(2,3)+c(3,3) ...
        +3*c(4,4)+3*c(5,5)+3*c(6,6));  %c055 in Turner
end
