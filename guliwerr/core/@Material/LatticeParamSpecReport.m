function spec = LatticeParamSpecReport(fname,symmNo)
    %%% Print (and return) which latticeParams a symmetry class requires
    %%% >> spec = Material.LatticeParamSpecReport(fname,symmNo)
    %%%
    %%% INPUT:  1. fname  .. name of the calling function, for the printed prefix
    %%%         2. symmNo .. symmetry class number (see SymmClassLibrary)
    %%% OUTPUT: spec .. the LatticeParamSpec struct for symmNo
    %%% NOTE:   Private -- backs the requirements-listing mode of Miller2DirMat and
    %%%         DirMat2Miller.

    spec = Material.LatticeParamSpec(symmNo);
    if isempty(spec.required)
        msg = sprintf('%s: no latticeParams needed (angles/ratios are fixed by symmetry).', ...
            spec.className);
    else
        msg = sprintf('%s: latticeParams must include: %s', ...
            spec.className, strjoin(spec.required,', '));
    end
    disp(['!(N) Material.' fname ': ' msg])
end
