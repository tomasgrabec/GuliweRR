function [indSet_ReductElast,nReductElast] = GetIndOf2DReduct_15of21()
    %%% Indices of the 15 (of 21) elastic constants surviving the 2D reduction
    %%% >> [indSet_ReductElast,nReductElast] = GetIndOf2DReduct_15of21()
    %%%
    %%% See ElastRotReductMat for which Voigt pairs are dropped and why.
    indSet_ReductElast = [1 3 4 5 6 12 13 14 15 16 17 18 19 20 21];
    nReductElast = numel(indSet_ReductElast);
end
