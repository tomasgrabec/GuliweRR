function init()
% INIT  Add all GuliweRR subfolders to the MATLAB path.
%
% Usage (run once per MATLAB session, from the project root):
%   run('guliwerr/init.m')
%   -- or, if guliwerr/ is already on the path --
%   init()
%
% -------------------------------------------------------------------------
% TIP: To load GuliweRR automatically at every MATLAB startup, add the
% following line to your personal startup.m file.
% (In MATLAB: edit(fullfile(userpath,'startup.m')) to open or create it.)
%
%   run('/full/path/to/Guliwerr/guliwerr/init.m')
%
% IMPORTANT: Do NOT use addpath(genpath(...)) on the guliwerr/ folder.
% genpath() recurses into @ClassName/ subfolders and adds them directly to
% the path, which breaks MATLAB's class-method lookup for those classes.
% This init.m adds only the parent folders, which is what MATLAB requires.
% -------------------------------------------------------------------------

    base = fileparts(mfilename('fullpath'));

    addpath(base);
    addpath(fullfile(base, 'core'));
    addpath(fullfile(base, 'solvers', 'christoffel'));
    addpath(fullfile(base, 'solvers', 'green'));
    addpath(fullfile(base, 'solvers', 'ritz'));
    addpath(fullfile(base, 'forwardCalcs'));
    addpath(fullfile(base, 'inverseCalcs'));
    addpath(fullfile(base, 'visual'));
    addpath(fullfile(base, 'gui'));

    fprintf('GuliweRR v%s ready.\n', guliwerrVersion());
end
