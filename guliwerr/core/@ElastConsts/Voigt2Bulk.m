function B = Voigt2Bulk(elVoigt)
    %%% Calculate Bulk modulus
    %%% >> B = Voigt2Bulk(elVoigt)
    %%%
    %%% INPUT:  elVoigt .. elastic Voigt matrix (6x6)

    compVoigt = inv(elVoigt); %%% compliance matrix; its 3x3 block sums to 1/B
    B = 1 / sum(compVoigt(1:3,1:3),'all');
end
