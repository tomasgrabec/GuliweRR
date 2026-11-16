function lib = ReadLibraryFile()
    %%% Read the materials library (guliwerr/data/materialsLibrary.csv) into a table
    %%% >> lib = ReadLibraryFile()
    %%%
    %%% OUTPUT: lib .. table, one row per data set; columns as in the file header
    %%%                (name, aliases, formula, structure, symmNo, consts, rho, T, variant,
    %%%                 default, status, reference, source, note); consts is a cell of
    %%%                numeric row vectors, symmNo/rho/default are numeric, the rest strings
    %%%
    %%% The file is found from this file's own location (guliwerr/core/@Material ->
    %%% guliwerr/data), so it does not need to be on the MATLAB path.

    libDir  = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    libPath = fullfile(libDir, 'data', 'materialsLibrary.csv');
    if ~isfile(libPath)
        error('Material:ReadLibraryFile:noFile', 'Materials library not found: %s', libPath);
    end

    %%% plain line split (no readtable): '#' lines are comments, the first other line is
    %%% the header; the Julia port parses the file the same way
    txt = readlines(libPath);
    txt = txt(strlength(strtrim(txt)) > 0 & ~startsWith(strtrim(txt), '#'));
    header = split(txt(1), ',')';
    cells = strings(numel(txt) - 1, numel(header));
    for r = 2:numel(txt)
        fields = split(txt(r), ',')';
        if numel(fields) ~= numel(header)
            error('Material:ReadLibraryFile:badRow', ...
                'Library row %d has %d fields, the header has %d: %s', ...
                r - 1, numel(fields), numel(header), txt(r));
        end
        cells(r - 1, :) = strtrim(fields);
    end
    lib = array2table(cells, 'VariableNames', cellstr(header));

    lib.symmNo  = str2double(lib.symmNo);
    lib.rho     = str2double(lib.rho);
    lib.default = str2double(lib.default);
    lib.consts  = arrayfun(@(s) str2double(split(strtrim(s)))', lib.consts, 'UniformOutput', false);
end
