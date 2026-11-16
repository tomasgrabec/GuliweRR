function type = FindType(input)
    %%% Identify the format of elastic constant input from its dimensions
    %%% >> type = FindType(input)
    %%%
    %%% INPUT:  tensor 3x3x3x3 / Voigt matrix 6x6 / vector (any length)
    %%% OUTPUT: 'tensor' / 'voigt' / 'vector' / 'unknown'

    switch ndims(input)
        case 4
            type = 'tensor';
        case 2
            if any(size(input) == 1)
                type = 'vector';
            elseif all(size(input) == 6)
                type = 'voigt';
            else
                type = 'unknown';
            end
        otherwise
            type = 'unknown';
    end
end
