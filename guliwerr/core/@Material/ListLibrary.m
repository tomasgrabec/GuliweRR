function lib = ListLibrary()
    %%% List the materials library (guliwerr/data/materialsLibrary.csv)
    %%% >> lib = Material.ListLibrary()
    %%%
    %%% OUTPUT: lib .. table, one row per data set (see Material.ReadLibraryFile);
    %%%                with no output argument a short overview is printed instead

    lib = Material.ReadLibraryFile();
    if nargout == 0
        fprintf('%-13s %-13s %-5s %-44s %-5s %-8s %-8s %-12s %s\n', 'name', 'structure', ...
            'symm', 'consts [GPa]', 'T', 'rho', 'variant', 'status', 'default');
        for r = 1:height(lib)
            dflt = '';
            if lib.default(r) == 1
                dflt = '*';
            end
            fprintf('%-13s %-13s %-5d %-44s %-5s %-8.5g %-8s %-12s %s\n', lib.name(r), ...
                lib.structure(r), lib.symmNo(r), strtrim(sprintf('%g ', lib.consts{r})), lib.T(r), lib.rho(r), ...
                lib.variant(r), lib.status(r), dflt);
        end
        clear lib
    end
end
