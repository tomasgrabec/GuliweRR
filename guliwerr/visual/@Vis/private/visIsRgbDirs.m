function isRgb = visIsRgbDirs(dirs, fnName)
%%% Resolve a modal-shape plot's dirs argument: 1-3 index vector or 'rgb'
%%% >> isRgb = visIsRgbDirs(dirs, fnName)
%%%
%%% @Vis-private helper shared by Vis.plotModalShape, plotGreenModalShape and
%%% plotGreenModalShapeLayered. dirs = 'rgb' asks for the merged polarisation
%%% view (Vis.calcPolarisationRgb) instead of one signed panel per direction.
%%%
%%% INPUT: 1. dirs   .. numeric vector of directions (1-3) or 'rgb'
%%%        2. fnName .. calling function's name, for the error id
%%% OUTPUT: isRgb .. true for 'rgb'; errors 'Vis:<fnName>:badDirs' on anything else
%%%                  that is neither 'rgb' nor a vector of 1-3

    if ischar(dirs) || isstring(dirs)
        if ~strcmpi(dirs, 'rgb')
            error(['Vis:' fnName ':badDirs'], ...
                'dirs must be a vector of directions 1-3 or ''rgb'' (got ''%s'').', char(dirs));
        end
        isRgb = true;
        return;
    end
    if ~isnumeric(dirs) || isempty(dirs) || ~isvector(dirs) || any(~ismember(dirs, 1:3))
        error(['Vis:' fnName ':badDirs'], 'dirs must be a vector of directions 1-3 or ''rgb''.');
    end
    isRgb = false;
end
