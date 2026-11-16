function posList = IndepVoigtPositions(symmNo)
    %%% Grid positions of each independent Voigt constant, for symmetry class symmNo
    %%% >> posList = IndepVoigtPositions(symmNo)
    %%%
    %%% INPUT:  symmNo   .. symmetry class number (see SymmClassLibrary)
    %%% OUTPUT: posList  .. cell array of [row,col] pairs (upper-triangle 6x6 Voigt
    %%%                     positions), one per obj.vector entry, in the SAME order
    %%%                     Vector2Voigt reads inpVec(k) for that symmNo
    %%%
    %%% Derived from Vector2Voigt, the one table of the symmetry classes: constant k
    %%% sits at the first upper-triangle position, row by row, where the Voigt matrix
    %%% of the unit vector e_k holds a 1 and that of every other e_j a 0 (e.g. cubic
    %%% c11 at [1,1], not at [2,2]; hexagonal c66 = (c11-c12)/2 is no position).
    %%% Cached per class (the ElastConsts constructor calls this through Voigt2Vector).

    persistent cache
    if isempty(cache); cache = cell(1, 10); end

    if ~(isnumeric(symmNo) && isscalar(symmNo) && any(symmNo == 1:10))
        error('ElastConsts:IndepVoigtPositions:unknownSymm', ...
            'Unknown symmetry class number: %d (must be 1-10)', symmNo);
    end

    if isempty(cache{symmNo})
        [~, nConsts] = ElastConsts.SymmClassLibrary(symmNo);
        B = zeros(6, 6, nConsts);                  % Voigt matrix of each unit vector
        for k = 1:nConsts
            e = zeros(1, nConsts);
            e(k) = 1;
            B(:, :, k) = ElastConsts.Vector2Voigt(e, symmNo);
        end
        alone = sum(B ~= 0, 3) == 1;               % positions fed by one constant only
        list  = cell(1, nConsts);
        for k = 1:nConsts
            hit = alone & B(:, :, k) == 1;
            [cols, rows] = find(triu(hit).');      % transposed: row-major scan order
            if isempty(rows)
                error('ElastConsts:IndepVoigtPositions:noPosition', ...
                    'Constant %d of symmetry class %d has no position of its own.', k, symmNo);
            end
            list{k} = [rows(1), cols(1)];
        end
        cache{symmNo} = list;
    end
    posList = cache{symmNo};
end
