function [vG, vP, eigVecsOut] = groupVelVectorsForDir(elTensor, dirVec, rho)
%%% Bulk group-velocity vectors for an ARBITRARY 3D direction
%%% >> [vG, vP, eigVecsOut] = Christoffel.groupVelVectorsForDir(elTensor, dirVec, rho)
%%%
%%% INPUT:  1. elTensor .. 3x3x3x3 tensor of elastic constants
%%%         2. dirVec   .. propagation direction (3-vector, need not be normalised)
%%%  (opt)  3. rho      .. density (default = 1)
%%% OUTPUT: vG         .. 3x3, columns are group velocity vectors [qL, qT1, qT2]
%%%         vP         .. [3x1] phase velocities, descending
%%%         eigVecsOut .. 3x3, columns are polarisation eigenvectors
%%%
%%% Royer–Dieulesaint Eq. 4.24 contraction, vG_i = C(i,j,k,l) g(j) g(k) s(l) / rho
%%% (g the polarisation, s = n / vP the slowness), for an arbitrary 3D
%%% direction, e.g. full-sphere ray surfaces. See docs/reference/architecture/christoffel.md.

    if nargin < 3 || isempty(rho); rho = 1; end
    %%% cheap checks (a ray surface calls this per grid direction); the tensor itself is
    %%% validated by Christoffel.matrix
    if ~(isnumeric(rho) && isscalar(rho) && isreal(rho) && rho > 0 && isfinite(rho))
        error('Christoffel:groupVelVectorsForDir:badRho', 'rho must be a positive finite scalar.');
    end
    if ~(isnumeric(dirVec) && numel(dirVec) == 3 && any(dirVec(:) ~= 0))
        error('Christoffel:groupVelVectorsForDir:badDirVec', 'dirVec must be a nonzero 3-vector.');
    end

    normC = elTensor ./ rho;
    n = dirVec(:) ./ norm(dirVec(:));

    normC_perm = reshape(permute(normC, [1,4,2,3]), 9, 9);

    [vG, vP, eigVecsOut] = Christoffel.GroupVelForDir(normC, normC_perm, n);
end
