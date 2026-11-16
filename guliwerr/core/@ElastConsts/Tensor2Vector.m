function vec = Tensor2Vector(inpTen,symmNo)
    %%% Transform tensor of elast consts to vector of indep. consts
    %%% >> vec = Tensor2Vector(inpTen,symmNo)
    %%%
    %%% INPUT:  1. inpTen .. elastic tensor (3x3x3x3)
    %%%  (opt)  2. symmNo .. symmetry class number; inferred if omitted or 0
    %%% OUTPUT: vec .. row vector of independent constants

    if nargin == 1 || symmNo == 0
        symmNo = ElastConsts.FindSymmClass_input(inpTen);
    end

    c = ElastConsts.Tensor2Voigt(inpTen);
    vec = ElastConsts.Voigt2Vector(c,symmNo);
end
