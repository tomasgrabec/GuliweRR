function vels = phaseVels(elTensor, dirVec, rho)
%%% Phase velocities for given tensor and direction
%%% >> vels = Christoffel.phaseVels(elTensor, dirVec, rho)
%%%
%%% INPUT:  1. elTensor .. 3×3×3×3 tensor of elastic constants
%%%         2. dirVec   .. propagation direction (3-vector)
%%%  (opt)  3. rho      .. density (default = 1)
%%% OUTPUT: [vL; vT1; vT2] phase velocities, sorted descending

if nargin < 3 || isempty(rho)
    rho = 1;
end

eigVals = Christoffel.matrix(elTensor, dirVec, rho);
vels = sqrt(eigVals);
end
