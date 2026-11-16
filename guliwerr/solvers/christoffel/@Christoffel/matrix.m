function [eigVals, eigVecs, GammaMat] = matrix(elTensor, dirVec, rho, normalizeDir)
%%% Christoffel matrix and its eigendecomposition
%%% >> [eigVals, eigVecs, GammaMat] = Christoffel.matrix(elTensor, dirVec, rho, normalizeDir)
%%%
%%% INPUT:  1. elTensor    .. 3×3×3×3 tensor of elastic constants
%%%         2. dirVec      .. propagation direction (3-vector)
%%%  (opt)  3. rho         .. density (default = 1; pass 1 for density-normalised input)
%%%  (opt)  4. normalizeDir .. normalise dirVec before use (default = true)
%%% OUTPUT: 1. eigVals    .. eigenvalues sorted descending (= rho*v_P^2 for L, T1, T2)
%%%         2. eigVecs    .. corresponding eigenvectors as columns (polarisations)
%%%         3. GammaMat   .. 3×3 Christoffel matrix

%%% check input
if nargin < 4 || isempty(normalizeDir)
    normalizeDir = true;
end
if nargin < 3 || isempty(rho)
    rho = 1;
end

validateattributes(elTensor, {'double'}, {'size', [3 3 3 3]}, 'Christoffel.matrix', 'elTensor');
validateattributes(dirVec,   {'double'}, {'vector', 'numel', 3}, 'Christoffel.matrix', 'dirVec');

%%% normalise direction
n = dirVec(:);
if normalizeDir
    n = n ./ norm(n);
end

%%% Christoffel matrix: Gamma(i,k) = C(i,j,k,l)*n(j)*n(l) / rho
%%% Not the Ritz GammaMatrix class (guliwerr/solvers/ritz/@GammaMatrix) —
%%% same name, different physics object.
GammaMat = Christoffel.GammaFromKron(elTensor, kron(n, n), rho);

%%% eigendecomposition, sorted descending by eigenvalue
[eigVals, eigVecs] = Christoffel.SortedEig(GammaMat);
end
