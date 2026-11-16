classdef ZgvNewton
% ZgvNewton — Newton refinement of zero-group-velocity (ZGV) points of a free plate.
%
%   Solves, for real k, real ω and a unit vector u,
%       (k² L2 + k L1 + L0 − ω² M) u = 0            (the plate's Ritz eigenproblem)
%       uᴴ (2k L2 + L1) u            = 0            (group velocity along k is zero)
%       uᴴ u                         = 1
%   by Newton's method on p = [u; k; ω²] with the bordered Jacobian and the complex
%   correction β of Lu & Su, which keeps k and ω² real. Locally quadratic; needs no mode
%   number, so it does not care which branch is which at a crossing.
%   L2, L0, M real symmetric (or Hermitian), L1 Hermitian (it carries the i): the
%   operator of RitzSolverLamb.kPolynomial.
%
%   Adapted from GEWtool's utilities/ZGVNewtonBeta.m (v2.1.1) by Daniel A. Kiefer; the
%   algorithm was designed by Bor Plestenjak. Differences from the original: the operator
%   sign convention (GULIWERR's L1 already carries the i, ω² enters with a minus sign), the
%   unknowns are scaled by the seed (k/k0, ω²/ω0²) so one tolerance fits every plate and
%   unit system, an ill-conditioned or non-finite Newton step stops the iteration as not
%   converged, and the eigenvector seed comes from the eigenproblem at k0 (seed method)
%   rather than from an SVD.
%   Reference: D. A. Kiefer, B. Plestenjak, H. Gravenkamp, C. Prada, "Computing
%   zero-group-velocity points in anisotropic elastic waveguides: Globally and locally
%   convergent methods", J. Acoust. Soc. Am. 153(2), 1386-1398 (2023),
%   doi:10.1121/10.0017252.
%
%   GEWtool licence (applies to the adapted parts of this file):
%   ISC License
%
%   Copyright (c) 2022-2026 Daniel A. Kiefer
%
%   Permission to use, copy, modify, and/or distribute this software for any
%   purpose with or without fee is hereby granted, provided that the above
%   copyright notice and this permission notice appear in all copies.
%
%   THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
%   REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
%   AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
%   INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
%   LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
%   OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
%   PERFORMANCE OF THIS SOFTWARE.
%
%%% GuliweRR v6.0

    properties (Constant)
        MAX_STEPS  = 10      % Newton steps (GEWtool default)
        TOL        = 1e-14   % residual tolerance, relative to ||[L0 L1 L2 M]||_F of the
                             % scaled operator (GEWtool default)
        HERM_TOL   = 1e-10   % largest relative anti-Hermitian part accepted in L2, L1, L0, M
        RCOND_MIN  = 1e-15   % a Jacobian with rcond below this stops the iteration
        KMIN_REL   = 1e-4    % k below KMIN_REL*k0 = a cut-off (k -> 0), not a ZGV point
                             % (GEWtool default)
    end

    methods (Static)
        function [k, omega, u, isConverged, nSteps] = solve(L2, L1, L0, M, k0, omega0, u0)
        %%% Newton iteration for one ZGV point from a seed (k0, omega0, u0)
        %%% >> [k, omega, u, isConverged, nSteps] = ZgvNewton.solve(L2, L1, L0, M, k0, omega0, u0)
        %%%
        %%% INPUT:  1.-4. L2, L1, L0, M .. [n x n] operator matrices (see the class header)
        %%%         5. k0     .. seed wavevector (real, > 0)
        %%%         6. omega0 .. seed angular frequency (real, > 0), same units as the
        %%%                      operator: (k² L2 + k L1 + L0) u = ω² M u
        %%%  (opt)  7. u0     .. seed vector [n x 1]; default (omitted or []): the
        %%%                      eigenvector of ZgvNewton.seedVector, and omega0 is then
        %%%                      replaced by that eigenvector's own frequency
        %%% OUTPUT: k, omega   .. the ZGV point; NaN when not converged
        %%%         u          .. unit vector there (NaN when not converged)
        %%%         isConverged.. residual below TOL, finite steps, Jacobian regular there
        %%%                       (rcond >= RCOND_MIN: an isolated ZGV), k >= KMIN_REL*k0
        %%%         nSteps     .. Newton steps taken
        %%%
        %%% Errors ZgvNewton:solve:nonHermitian when an operator matrix is not Hermitian
        %%% to HERM_TOL: the complex correction β and the group-velocity row assume a
        %%% lossless (Hermitian) problem, and a non-Hermitian one would need GEWtool's
        %%% ZGVNewtonComplex, not ported.
            validateattributes(k0, {'double'}, {'scalar', 'real', 'finite', 'positive'}, 'ZgvNewton.solve', 'k0');
            validateattributes(omega0, {'double'}, {'scalar', 'real', 'finite', 'positive'}, 'ZgvNewton.solve', 'omega0');
            mats = {L2, L1, L0, M};
            for im = 1:4
                A = mats{im};
                if norm(A - A', 'fro') > ZgvNewton.HERM_TOL * max(norm(A, 'fro'), realmin)
                    error('ZgvNewton:solve:nonHermitian', ...
                        'The operator matrices must be Hermitian (lossless plate); matrix %d is not.', im);
                end
            end
            if nargin < 7 || isempty(u0)
                [u0, omega0] = ZgvNewton.seedVector(L2, L1, L0, M, k0, omega0);
            end

            %%% scaled unknowns kap = k/k0, m = ω²/ω0²: A(kap, m) = kap² B2 + kap B1 + B0 − m BM
            mu0 = omega0^2;
            B2  = (k0^2 / mu0) * (L2 + L2') / 2;
            B1  = (k0   / mu0) * (L1 + L1') / 2;
            B0  = (1    / mu0) * (L0 + L0') / 2;
            BM  = (M + M') / 2;
            n   = size(B0, 1);
            relTol = ZgvNewton.TOL * norm([B0, B1, B2, BM], 'fro');

            kap = 1;  m = 1;
            u   = u0(:) / norm(u0);
            ek  = zeros(n + 2, 1);  ek(n + 1) = 1;
            isConverged = false;
            nSteps = 0;
            ws = warning('off', 'MATLAB:nearlySingularMatrix');
            ws2 = warning('off', 'MATLAB:singularMatrix');
            cleanupW = onCleanup(@() warning([ws, ws2]));
            while true
                A  = kap^2 * B2 + kap * B1 + B0 - m * BM;
                Ak = 2 * kap * B2 + B1;                   % dA/dkap
                F  = [A * u; u' * Ak * u; u' * u - 1];
                J  = [A,              Ak * u,            -BM * u;
                      2 * u' * Ak,    2 * (u' * B2 * u),  0;
                      2 * u',         0,                  0];
                %%% one LU factorisation serves the condition estimate and the step
                isRegular = all(isfinite(J(:)));
                if isRegular
                    dJ = decomposition(J, 'lu', 'CheckCondition', false);
                    isRegular = rcond(dJ) >= ZgvNewton.RCOND_MIN;
                end
                if norm(F) < relTol
                    %%% a solution with a singular Jacobian is not an isolated ZGV point
                    %%% (e.g. a band flat in k, or a ZGV coinciding with a crossing)
                    isConverged = isRegular;
                    break
                end
                if nSteps >= ZgvNewton.MAX_STEPS || ~isRegular
                    break
                end
                nSteps = nSteps + 1;
                kor = -(dJ \ [F, ek]);
                q = kor(:, 1);
                s = kor(:, 2);
                beta = -imag(q(n + 1)) / real(s(n + 1));  % Lu & Su: keeps kap real
                if isfinite(beta) && abs(beta) > 1e-10
                    dp = q + 1i * beta * s;
                else
                    dp = q;
                end
                if ~all(isfinite(dp))
                    break
                end
                u   = u + dp(1:n);
                kap = kap + real(dp(n + 1));             % imaginary parts are round-off
                m   = m + real(dp(n + 2));
                u   = u / norm(u);
            end

            k     = kap * k0;
            omega = sqrt(max(m, 0) * mu0);
            if ~isConverged || ~(kap >= ZgvNewton.KMIN_REL) || m <= 0
                isConverged = false;
                k = NaN;  omega = NaN;  u = NaN(n, 1);
            end
        end

        function [u0, omegaSeed] = seedVector(L2, L1, L0, M, k0, omega0)
        %%% Seed vector: the eigenvector at k0 whose eigenvalue is nearest omega0²
        %%% >> [u0, omegaSeed] = ZgvNewton.seedVector(L2, L1, L0, M, k0, omega0)
        %%%
        %%% OUTPUT: u0        .. unit eigenvector
        %%%         omegaSeed .. its angular frequency at k0
        %%%
        %%% The plate's own eigenproblem at k0 (a point on the dispersion surface), not
        %%% GEWtool's right singular vector of A(k0, ω0): both pick the mode nearest the
        %%% seed, but the eigenvector is the exact one of a nearby point of the branch.
            A  = k0^2 * L2 + k0 * L1 + L0;
            A  = (A + A') / 2;
            Ms = real(M + M') / 2;
            [V, D] = eig(A, Ms);
            lam = real(diag(D));
            [~, j] = min(abs(lam - omega0^2));
            u0 = V(:, j) / norm(V(:, j));
            omegaSeed = sqrt(max(lam(j), realmin));
        end
    end
end
