function [vels, dVels, eigVecs] = phaseVelAngleDeriv(elTensor, angleRad, rho)
%%% Phase velocities and their first derivative w.r.t. sagittal-plane angle
%%% >> [vels, dVels, eigVecs] = Christoffel.phaseVelAngleDeriv(elTensor, angleRad, rho)
%%%
%%% INPUT:  1. elTensor  .. 3×3×3×3 tensor of elastic constants
%%%         2. angleRad  .. sagittal-plane angle IN RADIANS (direction
%%%                         [cos; 0; sin]), unlike phaseVelScan's degrees.
%%%                         Angle-derivative formulas are radian-native; mixing
%%%                         the two costs a silent factor of 180/pi.
%%%  (opt)  3. rho        .. density (default = 1)
%%% OUTPUT: 1. vels   .. [3×1] phase velocities [vL; vT1; vT2], descending
%%%         2. dVels  .. [3×1] d(vels)/d(angleRad), same order
%%%         3. eigVecs .. [3×3] polarisations (columns, same order) -- exactly what
%%%                       Christoffel.matrix returns along [cos; 0; sin], from the same solve
%%%
%%% Derivation: docs/reference/architecture/christoffel.md.

if nargin < 3 || isempty(rho)
    rho = 1;
end

validateattributes(elTensor, {'double'}, {'size', [3 3 3 3]}, 'Christoffel.phaseVelAngleDeriv', 'elTensor');

n  = [cos(angleRad); 0; sin(angleRad)];
dn = [-sin(angleRad); 0; cos(angleRad)];

%%% The tensor contraction is a [9x9] matrix times kron(n,n) (see GammaFromKron); the
%%% permute/reshape to that matrix is done once and serves both Gamma and dGamma/dtheta.
%%% The direction is normalised as Christoffel.matrix does (a no-op to round-off, kept
%%% so the values are bit-identical to calling matrix()).
C9 = reshape(permute(elTensor, [1,3,2,4]), 9, 9);
nU = n ./ norm(n);

%%% Christoffel matrix and its eigendecomposition
[eigVals, eigVecs] = Christoffel.SortedEig(reshape(C9 * kron(nU, nU), 3, 3) / rho);
vels = sqrt(eigVals);

%%% dGamma/dtheta: kron(n,n) replaced by its total differential
dGamma = reshape(C9 * (kron(dn, n) + kron(n, dn)), 3, 3) / rho;

dVels = zeros(3, 1);
for iwt = 1:3
    g = eigVecs(:, iwt);
    dVels(iwt) = (g' * dGamma * g) / (2 * vels(iwt));
end
end
