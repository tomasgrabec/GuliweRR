function indep2RedElastMat = ElastRotReductMat(rotMat,symmNo)
    %%% Map from independent constants to the 15 reduced constants of the 2D RR system,
    %%% for a given rotation
    %%% >> indep2RedElastMat = ElastRotReductMat(rotMat,symmNo)
    %%%
    %%% INPUT:  1. rotMat .. matrix of direction row vectors (2x3/3x2 or 3x3)
    %%%         2. symmNo .. symmetry class number (see SymmClassLibrary)
    %%% OUTPUT: indep2RedElastMat .. 15 x noIndepConsts
    %%%
    %%% The reduction drops every Voigt pair touching index 2 (the completion axis under the
    %%% standard x1=direction, x2=completion, x3=plane convention): e_22 vanishes for a field
    %%% with no completion-axis dependence, so C_I2/C_2J terms add nothing to the strain energy.

    [indReducedEl,nReductElast] = ElastConsts.GetIndOf2DReduct_15of21();
    [~,noIndepConsts] = ElastConsts.SymmClassLibrary(symmNo);

    indep2RedElastMat = NaN(nReductElast,noIndepConsts);

    for const_no = 1 : noIndepConsts
        dummyConsts = zeros(1,noIndepConsts);
        dummyConsts(const_no) = 1;
        cTen0 = ElastConsts.Any2Tensor(dummyConsts,symmNo);
        cTenRot = ElastConsts.RotateTensor(cTen0,rotMat);
        cVec = ElastConsts.Tensor2Vector(cTenRot,10);

        indep2RedElastMat(:,const_no) = cVec(indReducedEl)';
    end
end
