function [typeNames, typeColors, parityNames, parityColors] = plateModeTypeStyle()
%%% Display names and colours for the four Lamb/plate mode types and their two parity classes
%%% >> [typeNames, typeColors, parityNames, parityColors] = Vis.plateModeTypeStyle()
%%%
%%% OUTPUT: typeNames    .. {1 x 4} cell of names, indexed by
%%%                         RitzSolverLamb.modesType (1-4)
%%%         typeColors   .. [4 x 3] RGB, same indexing
%%%         parityNames  .. {1 x 4} cell, the parity class of each type
%%%                         ('Asymm-class' for 1 and 3, 'Symm-class' for 2 and 4) --
%%%                         the labels classifyBranches shows where only the parity is exact
%%%         parityColors .. [4 x 3] RGB, same indexing: the midpoint of the class's two
%%%                         type colours (Asymm-class teal between blue and green,
%%%                         Symm-class crimson between red and purple), so a parity-only
%%%                         branch reads as "this family, Lamb/TH undecided"
%%%
%%% Indexed BY modesType directly -- typeColors(modesType, :) -- so a caller
%%% never has to re-derive the mapping.
%%%
%%% Public rather than @Vis-private because the mapping is a presentation
%%% convention shared with callers OUTSIDE visual/: a GUI showing "this mode is
%%% Symm" beside a plot drawn here must use the same names and colours.

    typeNames  = {'Asymm', 'Symm', 'TH-Asymm', 'TH-Symm'};
    typeColors = [0.00 0.45 0.74;    % Asymm    -- blue
                  0.85 0.10 0.10;    % Symm     -- red
                  0.30 0.70 0.30;    % TH-Asymm -- green
                  0.60 0.20 0.70];   % TH-Symm  -- purple

    parityNames  = {'Asymm-class', 'Symm-class', 'Asymm-class', 'Symm-class'};
    asymmClass   = (typeColors(1, :) + typeColors(3, :)) / 2;   % teal
    symmClass    = (typeColors(2, :) + typeColors(4, :)) / 2;   % crimson
    parityColors = [asymmClass; symmClass; asymmClass; symmClass];
end
