function [obj,entry] = FromLibrary(name,dirMat,varargin)
    %%% Build a Material from the materials library (guliwerr/data/materialsLibrary.csv)
    %%% >> [obj,entry] = Material.FromLibrary(name,dirMat,'T',T,'variant',variant)
    %%% >> [obj,entry] = Material.FromLibrary(name,'T',T)        % dirMat may be left out
    %%%
    %%% INPUT:  1. name    .. library name or alias, case-insensitive (e.g. 'Cu', 'copper')
    %%%  (opt)  2. dirMat  .. orientation matrix (see Material); [] or omitted: Material default
    %%%  (opt)  'T'        .. temperature: number in K or 'RT'; must match a stored row
    %%%  (opt)  'variant'  .. variant label of the row (Material.ListLibrary shows them)
    %%% OUTPUT: 1. obj   .. Material (name = library name, plus the variant if not default)
    %%%         2. entry .. one-row table with the data set's reference, source and notes
    %%%
    %%% Without 'T' or 'variant' the material's default row is returned. For rows with status
    %%% 'provisional' (secondary source only, or a doubt given in the row's note) a one-line
    %%% note saying why is printed.

    if nargin < 2
        dirMat = [];
    elseif ischar(dirMat) || isstring(dirMat)   % name-value pairs straight after name
        varargin = [{dirMat}, varargin];
        dirMat = [];
    end
    p = inputParser;
    p.addParameter('T', []);
    p.addParameter('variant', "");
    p.parse(varargin{:});

    lib = Material.ReadLibraryFile();
    key = lower(string(name));
    isMat = lower(lib.name) == key | ...
            arrayfun(@(a) any(lower(split(a, '|')) == key), lib.aliases);
    if ~any(isMat)
        error('Material:FromLibrary:unknownName', ...
            'No material "%s" in the library. Available: %s.', name, ...
            strjoin(unique(lib.name, 'stable'), ', '));
    end
    rows = lib(isMat, :);

    sel = true(height(rows), 1);
    if ~isempty(p.Results.T)
        sel = sel & strcmpi(rows.T, string(p.Results.T));
    end
    if strlength(string(p.Results.variant)) > 0
        sel = sel & strcmpi(rows.variant, string(p.Results.variant));
    end
    if isempty(p.Results.T) && strlength(string(p.Results.variant)) == 0
        sel = rows.default == 1;
    end

    available = strjoin(compose("T=%s variant=%s", rows.T, rows.variant), '; ');
    if ~any(sel)
        error('Material:FromLibrary:noMatch', ...
            'No "%s" row matches the request. Available rows: %s.', rows.name(1), available);
    elseif nnz(sel) > 1
        error('Material:FromLibrary:ambiguous', ...
            'Several "%s" rows match; give ''T'' and/or ''variant''. Available rows: %s.', ...
            rows.name(1), available);
    end
    entry = rows(sel, :);

    matName = entry.name;
    if entry.default ~= 1
        matName = sprintf('%s (%s, T=%s)', entry.name, entry.variant, entry.T);
    end
    if entry.status == "provisional"
        %%% why it is provisional: no primary reference yet (secondary source only), or a
        %%% doubt stated in the note
        if entry.reference == "pending"
            why = sprintf('secondary source only: %s', entry.source);
        else
            why = entry.note;
        end
        fprintf('Note: %s -- provisional library values (%s)\n', matName, why);
    end
    obj = Material(ElastConsts(entry.consts{1}, entry.symmNo), entry.rho, dirMat, matName);
end
