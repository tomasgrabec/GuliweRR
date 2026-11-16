function vec = Voigt2Vector(cij,symmNo)
    %%% Transform Voigt matrix to vector of indep. consts
    %%% >> vec = Voigt2Vector(cij,symmNo)
    %%%
    %%% INPUT:  1. cij    .. Voigt matrix (6x6)
    %%%  (opt)  2. symmNo .. symmetry class number; inferred if omitted or 0
    %%% OUTPUT: vec .. row vector of independent constants, in IndepVoigtPositions order

    if nargin == 1 || symmNo == 0
        symmNo = ElastConsts.FindSymmClass_input(cij);
    end
    if ~(isnumeric(symmNo) && isscalar(symmNo) && any(symmNo == 1:10))
        error('ElastConsts:Voigt2Vector:unknownSymm', ...
            'Unknown symmetry class number: %d (must be 1-10)', symmNo);
    end

    posList = ElastConsts.IndepVoigtPositions(symmNo);
    pos     = vertcat(posList{:});                     % nConsts x 2, [row,col] per entry
    vec     = reshape(cij(sub2ind(size(cij), pos(:,1), pos(:,2))), 1, []);
end
