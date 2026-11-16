function E = Tensor2Youngs(elTen,rotMat)
    %%% Young's modulus along the x1 axis of the (optionally rotated) elastic tensor
    %%% >> E = Tensor2Youngs(elTen,rotMat)
    %%%
    %%% INPUT:  1. elTen  .. elastic tensor (3x3x3x3)
    %%%  (opt)  2. rotMat .. rotation matrix (2x3/3x2, or 3x3 row vectors)

    if nargin == 2 && ~isempty(rotMat)
        elTen = ElastConsts.RotateTensor(elTen,rotMat);
    end
    elVoigt = ElastConsts.Tensor2Voigt(elTen);
    compVoigt = inv(elVoigt); %%% compliance matrix; its (1,1) element is 1/E
    E = 1 / compVoigt(1,1);
end
