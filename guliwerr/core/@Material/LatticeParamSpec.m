function spec = LatticeParamSpec(symmNo)
    %%% Default lattice-shape parameters and which are required, by symmetry class
    %%% >> spec = Material.LatticeParamSpec(symmNo)
    %%%
    %%% INPUT:  symmNo .. symmetry class number (see SymmClassLibrary)
    %%% OUTPUT: spec.defaults  .. struct with all of b_over_a, c_over_a, alpha, beta, gamma
    %%%                           at their symmetry-implied values
    %%%         spec.required  .. field names the caller must supply
    %%%         spec.className .. symmetry class name
    %%%
    %%% CONVENTIONS: trigonal (5/6) in the hexagonal axis setting, not rhombohedral;
    %%% monoclinic (9) unique-axis-b (alpha = gamma = 90, beta free).

    defaults = struct('b_over_a',1,'c_over_a',1,'alpha',90,'beta',90,'gamma',90);
    switch symmNo
        case {1,2}      % isotropic, cubic -- no crystallographic anisotropy in shape
            required = {};
        case {3,4}      % tetragonal, type I/II
            required = {'c_over_a'};
        case {5,6,7}    % trigonal I/II (hex-axis setting), hexagonal
            defaults.gamma = 120;
            required = {'c_over_a'};
        case 8          % orthorhombic
            required = {'b_over_a','c_over_a'};
        case 9          % monoclinic, unique axis b
            required = {'b_over_a','c_over_a','beta'};
        case 10         % triclinic
            required = {'b_over_a','c_over_a','alpha','beta','gamma'};
        otherwise
            error('Material:LatticeParamSpec:badSymmNo', ...
                'Unknown symmetry class number %d (see SymmClassLibrary).', symmNo)
    end

    spec.defaults  = defaults;
    spec.required  = required;
    spec.className = ElastConsts.SymmClassLibrary(symmNo);
end
