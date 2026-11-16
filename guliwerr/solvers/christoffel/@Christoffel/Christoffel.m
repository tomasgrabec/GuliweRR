classdef Christoffel
%%% Static methods for bulk-wave Christoffel calculations.
%%% All methods are static — do not instantiate this class.
%%%
%%% Usage:  vels = Christoffel.phaseVels(elTensor, dirVec)
%%%         [vG, vP] = Christoffel.groupVelVectorsForDir(elTensor, dirVec)
%%%
%%% Input convention:
%%%   phaseVelScan takes a Material object; the other methods (matrix,
%%%   phaseVels, phaseVelAngleDeriv, groupVelVectorsForDir) take tensor + rho
%%%   directly.
%%%
%%%   Exception: phaseVelAngleDeriv takes its angle in RADIANS, not degrees
%%%   (angle-derivative formulas are radian-native) — see its own docblock.
%%%
%%% GuliweRR v6.0

    methods (Static)

        %%% Christoffel matrix and phase velocities
        [eigVals, eigVecs, GammaMat] = matrix(elTensor, dirVec, rho, normalizeDir)
        vels                         = phaseVels(elTensor, dirVec, rho)
        vels                         = phaseVelScan(mater, angleDeg)
        [vels, dVels, eigVecs]         = phaseVelAngleDeriv(elTensor, angleRad, rho)

        %%% Group velocities
        [vG, vP, eigVecsOut]         = groupVelVectorsForDir(elTensor, dirVec, rho)

    end

    methods (Static, Access = private)

        function GammaMat = GammaFromKron(elTensor, nnKron, rho)
        %%% Christoffel-style contraction Gamma(i,k) = C(i,j,k,l)*nn(j,l) / rho
        %%% >> GammaMat = Christoffel.GammaFromKron(elTensor, nnKron, rho)
        %%%
        %%% INPUT:  1. elTensor .. 3×3×3×3 tensor of elastic constants
        %%%         2. nnKron   .. 9×1 kron(n,n) — or its total differential
        %%%                        kron(dn,n)+kron(n,dn) for the angle derivative
        %%%  (opt)  3. rho      .. density (default = 1)
        %%%
        %%% The [1,3,2,4] permutation is what pairs (j,l) with kron's index order.
            if nargin < 3 || isempty(rho)
                rho = 1;
            end
            GammaMat = reshape(reshape(permute(elTensor, [1,3,2,4]), 9, 9) * nnKron, 3, 3) / rho;
        end

        function [eigVals, eigVecs] = SortedEig(M)
        %%% Eigendecomposition sorted DESCENDING by eigenvalue
        %%% >> [eigVals, eigVecs] = Christoffel.SortedEig(M)
        %%% Descending order is the [L, T1, T2] wave-type convention used
        %%% throughout this class.
        %%%
        %%% M is symmetrised first. A Christoffel matrix is symmetric physically,
        %%% but one built from a ROTATED elastic tensor is not, at round-off level,
        %%% and eig() then takes its general (non-symmetric) solver: for an exactly
        %%% degenerate pair (the two shear waves of an isotropic solid) that can
        %%% return a complex pair (3.7398 +- 1.8e-16i km/s, fused silica at 150 deg),
        %%% i.e. a complex limiting velocity. With exact symmetry eig() uses the
        %%% symmetric solver: real eigenvalues, orthonormal eigenvectors.
            [eigVecs, eigValsD] = eig((M + M.') / 2);
            [eigVals, eigOrder] = sort(diag(eigValsD), 'descend');
            eigVecs = eigVecs(:, eigOrder);
        end

        function [vG, vP, eigVecs] = GroupVelForDir(normC, normCperm, n)
        %%% Phase/group velocities and polarisations for ONE direction
        %%% >> [vG, vP, eigVecs] = Christoffel.GroupVelForDir(normC, normCperm, n)
        %%%
        %%% INPUT:  1. normC     .. density-NORMALISED 3×3×3×3 elastic tensor
        %%%         2. normCperm .. reshape(permute(normC,[1,4,2,3]),9,9), hoisted
        %%%                         by the caller (no direction dependence)
        %%%         3. n         .. propagation direction, 3×1
        %%% OUTPUT: 1. vG      .. 3×3, columns are group velocity vectors [L,T1,T2]
        %%%         2. vP      .. 3×1 phase velocities, descending
        %%%         3. eigVecs .. 3×3, columns are polarisation eigenvectors
        %%%
        %%% vG_i = C(i,j,k,l)*g(j)*g(k)*s(l), with s = n/vP the slowness vector
        %%% (Royer–Dieulesaint, Elastic Waves in Solids I, Eq. 4.24).
        %%%
        %%% NOTE: n is used AS GIVEN for the slowness — callers that need it
        %%% normalised normalise before calling; Christoffel.matrix normalises
        %%% internally for the eigenproblem regardless.
            [eigVals, eigVecs] = Christoffel.matrix(normC, n, 1);
            vP = sqrt(eigVals);

            vG = zeros(3, 3);
            for iwt = 1:3
                g = eigVecs(:, iwt);
                s = n ./ vP(iwt);                       % slowness vector
                P = g * g';                             % 3×3 outer product g⊗g
                T = reshape(normCperm * P(:), 3, 3);    % T(i,l) = sum_jk normC(i,j,k,l)*P(j,k)
                vG(:, iwt) = T * s;
            end
        end
    end

end
