function symmNo = FindSymmClass_input(inp)
    %%% Try to identify symmetry class based on input
    %%% >> symmNo = FindSymmClass_input(inp)
    %%%
    %%% INPUT:  inp .. elastic constants (vector / Voigt 6x6 / tensor 3x3x3x3)
    %%% OUTPUT: symmNo .. symmetry class number
    %%%
    %%% GOTCHA: reliable only from a vector of independent constants. Voigt (and tensor,
    %%% which is converted to Voigt first) input counts unique(cij) entries, which
    %%% miscounts whenever two independent constants happen to be numerically equal.
    %%% Prefer supplying symmNo explicitly.

    isVoigt = true;
    switch ElastConsts.FindType(inp)
        case 'tensor'
            cij = ElastConsts.Tensor2Voigt(inp);
        case 'voigt'
            cij = inp;
        case 'vector'
            isVoigt = false;
            NUnique = numel(inp);
        otherwise
            error('ElastConsts:FindSymmClass_input:badInput', ...
                'Input must be a vector, Voigt (6x6), or elastic tensor (3x3x3x3)');
    end

    if isVoigt
        CUnique = unique(cij);
        if CUnique(1) == 0
            CUnique(1) = [];   % the structural zeros of every class but triclinic
        end
        NUnique = length(CUnique);
    end

    %%% symmetry class implied by the number of independent constants
    switch NUnique
        case 2
            symmNo = 1;
            if isVoigt && cij(1,2) ~= (cij(1,1) - 2*cij(4,4))
                error('ElastConsts:FindSymmClass_input:notSymmetric', ...
                    'Symmetry check not successful, symmetries not fulfilled (assumed symmetry class: %d).', symmNo)
            end
        case 3
            symmNo = 2;
        case 5
            symmNo = 7;
        case 6
            if isVoigt
                if cij(1,4) == 0
                    symmNo = 3;       %tetragonal
                else
                    symmNo = 5;       %trigonal
                    if cij(6,6) ~= (cij(1,1) - cij(1,2))/2
                        error('ElastConsts:FindSymmClass_input:notSymmetric', ...
                            'Symmetry check not successful, symmetries not fulfilled (assumed symmetry class: %d).', symmNo)
                    end
                end
            else
                error('ElastConsts:FindSymmClass_input:ambiguousClass', ...
                    ['Symmetry class not given and cannot be determined (number of consts = 6 -> ', ...
                     'either symmClassNo = 3 for tetragonal or symmClassNo = 5 for trigonal).']);
            end
        case 7
            if isVoigt && cij(1,4) == 0
                symmNo = 4;
            elseif isVoigt
                symmNo = 6;
                if cij(6,6) ~= (cij(1,1) - cij(1,2))/2
                    error('ElastConsts:FindSymmClass_input:notSymmetric', ...
                            'Symmetry check not successful, symmetries not fulfilled (assumed symmetry class: %d).', symmNo)
                end
            else
                error('ElastConsts:FindSymmClass_input:ambiguousClass', ...
                    ['Symmetry class not given and cannot be determined (number of consts = 7 -> ', ...
                     'either symmClassNo = 4 for tetragonal or symmClassNo = 6 for trigonal).']);
            end
        case 9
            symmNo = 8;
        case 13
            symmNo = 9;
        case 21
            symmNo = 10;
        otherwise
            error('ElastConsts:FindSymmClass_input:ambiguous', ...
                'Cannot determine symmetry class: %d unique constants found', NUnique);
    end
end
