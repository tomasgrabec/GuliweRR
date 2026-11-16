function outTen = Any2Tensor(inp,symmNo)
    %%% Transform any form of elastic constants (vector/voigt/tensor) to tensor
    %%% >> outTen = Any2Tensor(inp,symmNo)
    %%%
    %%% INPUT:  1. inp    .. elastic constants (vector / Voigt 6x6 / tensor 3x3x3x3)
    %%%  (opt)  2. symmNo .. symmetry class number (see SymmClassLibrary)
    %%% OUTPUT: outTen .. elastic tensor (3x3x3x3)

    if nargin < 2 || isempty(symmNo)
        symmNo = 0;
    end

    outTen = zeros(3,3,3,3);
    switch ElastConsts.FindType(inp)
        case 'tensor'
            outTen = inp;

        case 'voigt'
            p = [1 6 5;           % Voigt index pair -> tensor index pair
                6 2 4;
                5 4 3];

            for i=1:3
                for j=1:3
                    for k=1:3
                        for l=1:3
                            outTen(i,j,k,l) = inp(p(i,j),p(k,l));
                        end
                    end
                end
            end

        case 'vector'
            cij = ElastConsts.Vector2Voigt(inp,symmNo);
            outTen = ElastConsts.Any2Tensor(cij,symmNo);
        otherwise
            error('ElastConsts:Any2Tensor:badInput', ...
                'Input must be a vector, Voigt (6x6), or elastic tensor (3x3x3x3)');
    end
end
