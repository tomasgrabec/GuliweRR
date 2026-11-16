classdef LimitingVels
%%% Static methods for computing limiting bulk-wave velocities.
%%%
%%% Limiting velocities belong to the bulk waves with wavevector in the sagittal
%%% plane whose group velocity there has no surface-normal component (vg2 = 0):
%%% the energy flux runs parallel to the surface, and the limiting velocity is
%%% its along-surface component vg1, equal to the phase velocity along the
%%% surface, v/cos(theta*). The lowest one bounds the existence domain of SAW.
%%%   findRoots  — limiting velocities for one pre-rotated tensor (the
%%%                rotation to the standard frame is the caller's, e.g.
%%%                CalcCommon.limitingAndBulk)
%%%   refineRoot — re-find ONE already-located root in a small bracket
%%%
%%% GuliweRR v6.0

    properties (Constant, Access = private)
        FMIN_OPTS = optimset('TolX', 1e-12)   % root refinement: see findRoots (theta is pinned
                                              % to ~sqrt(eps)|theta|, not to TolX)
    end

    properties (Constant)
        %%% Two roots of DIFFERENT wave types are one degenerate root (an acoustic axis in
        %%% the sagittal plane, e.g. cubic (001) at 0 deg, or any direction of an isotropic
        %%% medium) when their velocities agree to DEGENERATE_REL_TOL (relative) and their
        %%% theta* to DEGENERATE_THETA_TOL [rad]. Both are set from the root finder's own
        %%% precision: fminbnd stops at sqrt(eps)|theta| ~ 1.5e-8 rad, which moves the
        %%% velocity by up to ~3e-8 relative at the steepest crossings seen (the same root in
        %%% MATLAB and Julia, or a degenerate pair, differs by 0 to 1.4e-8). 1e-6 is 30x
        %%% above that jitter and far below any split the eigensolver resolves: at a gap of
        %%% 1e-6 relative the polarisations are defined to ~eps/gap ~ 1e-10.
        DEGENERATE_REL_TOL   = 1e-6
        DEGENERATE_THETA_TOL = 1e-6
    end

    methods (Static)

        function [vLim_sA, eigenVecs_sA, thetaRoots_sA] = findRoots(rotElTensor, rho, alphaPrecision, relTol)
            %%% Limiting velocities for a single pre-rotated elastic tensor
            %%% >> [vLim_sA, eigenVecs_sA, thetaRoots_sA] = LimitingVels.findRoots(rotElTensor, rho, alphaPrecision, relTol)
            %%%
            %%% INPUT:  1. rotElTensor    .. 3×3×3×3 rotated elastic tensor
            %%%  (opt)  2. rho            .. density (default = 1)
            %%%  (opt)  3. alphaPrecision .. sagittal-plane scan step [rad] (default = 1e-2)
            %%%  (opt)  4. relTol         .. accept threshold, relative to the in-plane
            %%%                              group-velocity magnitude at the candidate
            %%%                              root (default = 1e-3)
            %%% OUTPUT: 1. vLim_sA{iwt}      .. limiting velocities for wave type iwt
            %%%         2. eigenVecs_sA{iwt}  .. corresponding eigenvectors
            %%%         3. thetaRoots_sA{iwt} .. sagittal-plane angle theta* [rad] of each
            %%%                                 accepted root, SAME ORDER as vLim_sA{iwt}
            %%%
            %%% Coarse islocalmin scan over |vg2|, then a bounded fminbnd
            %%% refinement of every candidate before the accept/reject test.
            %%% islocalmin (rather than a sign change) is what also catches
            %%% tangential touches — bulk waves grazing near-zero out-of-plane
            %%% group velocity without crossing sign, a real surface-skimming
            %%% feature. Algorithm: docs/reference/architecture/christoffel.md.
            %%%
            %%% fminbnd runs at TolX 1e-12, not its default 1e-4 rad. At a steep
            %%% sign crossing the residual |vg2| left by a theta error is slope x
            %%% error: a strongly anisotropic cubic can cross at |dvg2/dtheta| of tens
            %%% of km/s/rad, where the default TolX leaves a residual above the accept
            %%% tolerance relTol x vg1 and a genuine root would be rejected. The theta
            %%% error would also move the accepted velocity, as vg1 is not stationary
            %%% at the root. Precision actually reached: fminbnd stops at
            %%% sqrt(eps)|theta| + TolX/3, so the relative term (~1.5e-8 rad at theta ~ 1)
            %%% dominates and a TolX below ~1e-9 changes nothing; theta* is good to about
            %%% 1e-8 rad, enough for the accept test and for refineRoot's finite
            %%% differences at steps >= 1e-3 GPa.
            %%%
            %%% thetaRoots_sA lets a caller anchor a root by its theta* and
            %%% re-find it with refineRoot instead of re-running the whole scan.

            if nargin < 4 || isempty(relTol)
                relTol = 1e-3;
            end
            if nargin < 3 || isempty(alphaPrecision)
                alphaPrecision = 1e-2;
            end
            if nargin < 2 || isempty(rho)
                rho = 1;
            end
            LimitingVels.checkTensorRho(rotElTensor, rho, 'findRoots');
            validateattributes(alphaPrecision, {'numeric'}, {'scalar','real','positive','<',1}, ...
                'LimitingVels.findRoots', 'alphaPrecision');
            validateattributes(relTol, {'numeric'}, {'scalar','real','positive'}, ...
                'LimitingVels.findRoots', 'relTol');

            normC = rotElTensor ./ rho;

            %%% coarse scan, radians throughout (phaseVelAngleDeriv is radian-native)
            theta  = -pi/2 + 1e-4 : alphaPrecision : pi/2 - 1e-4;
            nTheta = length(theta);
            vg2    = zeros(3, nTheta);
            for itheta = 1:nTheta
                vg2(:, itheta) = LimitingVels.vg2All(normC, theta(itheta));
            end

            vLim_sA       = cell(3, 1);
            eigenVecs_sA  = cell(3, 1);
            thetaRoots_sA = cell(3, 1);

            for iwt = 1:3
                candidates = find(islocalmin(abs(vg2(iwt, :))));

                vTemp  = [];
                eVtemp = [];
                thTemp = [];
                for icand = 1:length(candidates)
                    ic = candidates(icand);
                    lo = theta(max(ic - 1, 1));
                    hi = theta(min(ic + 1, nTheta));

                    [vg1Root, thetaRoot, ~, ok, polRoot] = LimitingVels.refineAt(normC, lo, hi, iwt, relTol);

                    if ok
                        vTemp(end + 1, :)  = vg1Root; %#ok<AGROW>
                        thTemp(end + 1, :) = thetaRoot; %#ok<AGROW>
                        if nargout > 1
                            eVtemp(end + 1, :) = polRoot'; %#ok<AGROW>
                        end
                    end
                end

                vLim_sA{iwt}       = vTemp;
                thetaRoots_sA{iwt} = thTemp;
                if nargout > 1
                    eigenVecs_sA{iwt} = eVtemp;
                end
            end

            if nargout > 1
                eigenVecs_sA = LimitingVels.canonicaliseDegenerate(normC, vLim_sA, ...
                    eigenVecs_sA, thetaRoots_sA);
            end
        end

        function [vLim, thetaRoot, pol, residual, ok] = refineRoot(rotElTensor, rho, ...
                thetaGuess, iwt, halfWidth, relTol)
            %%% Re-find ONE already-located sagittal-plane root in a small bracket
            %%% >> [vLim, thetaRoot, pol, residual, ok] = LimitingVels.refineRoot(rotElTensor, rho, thetaGuess, iwt, halfWidth, relTol)
            %%%
            %%% INPUT:  1. rotElTensor .. 3×3×3×3 rotated elastic tensor
            %%%         2. rho         .. density
            %%%         3. thetaGuess  .. sagittal-plane angle [rad] of the root being
            %%%                           tracked, from the previous accepted iterate
            %%%         4. iwt         .. wave-type branch index (1 = longitudinal,
            %%%                           2 = fast shear, 3 = slow shear), the same
            %%%                           descending-eigenvalue convention findRoots uses
            %%%  (opt)  5. halfWidth   .. bracket half-width [rad] (default = 0.05, about
            %%%                           2.9 deg, five steps of findRoots' coarse grid)
            %%%  (opt)  6. relTol      .. accept threshold, relative to the in-plane
            %%%                           group-velocity magnitude (default = 1e-3, the
            %%%                           same test findRoots applies)
            %%% OUTPUT: 1. vLim      .. limiting velocity at the refined root
            %%%         2. thetaRoot .. refined theta* [rad]
            %%%         3. pol       .. [3×1] polarisation at the refined root
            %%%         4. residual  .. |vg2| at the refined root
            %%%         5. ok        .. accept test passed. FALSE is the
            %%%                         root-DISAPPEARANCE detector: a limiting root
            %%%                         can vanish outright under a small change in
            %%%                         the elastic constants, and then the minimum of
            %%%                         |vg2| inside the bracket stops being near zero.
            %%%                         A caller MUST check it and fall back to a full
            %%%                         findRoots rescan; vLim is meaningless when false.
            %%%
            %%% Use this, not findRoots, inside a finite-difference derivative.
            %%% findRoots brackets each root from whichever coarse-grid cell islocalmin
            %%% flagged, which makes the returned velocity a non-smooth function of the
            %%% elastic constants at the 1e-4 GPa scale a central finite difference
            %%% perturbs them by — a tiny elasticity change can move the flagged cell
            %%% by one step and jump the bracket. Here the bracket is fixed by the
            %%% CALLER (identical for the +h and -h evaluations), which makes the value
            %%% smooth and the difference valid.
            %%%
            %%% Freezing theta* instead and differentiating only the explicit
            %%% c-dependence is NOT a valid shortcut: v_lim = v_g1(theta*(c), c) and
            %%% v_g1 is not stationary in theta at the root (only v_g2 is), so the
            %%% dtheta*/dc term is not small — it can exceed the frozen-theta term
            %%% several-fold.

            if nargin < 6 || isempty(relTol)
                relTol = 1e-3;
            end
            if nargin < 5 || isempty(halfWidth)
                halfWidth = 0.05;
            end
            LimitingVels.checkTensorRho(rotElTensor, rho, 'refineRoot');
            validateattributes(thetaGuess, {'numeric'}, {'scalar','real','finite'}, ...
                'LimitingVels.refineRoot', 'thetaGuess');
            validateattributes(iwt, {'numeric'}, {'scalar','integer','>=',1,'<=',3}, ...
                'LimitingVels.refineRoot', 'iwt');
            validateattributes(halfWidth, {'numeric'}, {'scalar','real','positive'}, ...
                'LimitingVels.refineRoot', 'halfWidth');
            validateattributes(relTol, {'numeric'}, {'scalar','real','positive'}, ...
                'LimitingVels.refineRoot', 'relTol');

            normC = rotElTensor ./ rho;

            %%% same domain findRoots scans, so a bracket can never leave it
            thetaLo = -pi/2 + 1e-4;
            thetaHi =  pi/2 - 1e-4;
            lo = max(thetaGuess - halfWidth, thetaLo);
            hi = min(thetaGuess + halfWidth, thetaHi);

            [vLim, thetaRoot, residual, ok, polRoot] = LimitingVels.refineAt(normC, lo, hi, iwt, relTol);

            if nargout > 2
                pol = polRoot;
            else
                pol = [];
            end
        end

    end

    methods (Static)

        function B = canonicalBasis(E)
            %%% Basis-independent orthonormal basis of the span of orthonormal columns E
            %%% >> B = LimitingVels.canonicalBasis(E)
            %%%
            %%% From the projector P = E E' (the same matrix for ANY orthonormal basis of
            %%% the subspace, to round-off): take P's longest column (ties within 1e-10
            %%% relative -> the lowest index), normalise it, deflate P by it, repeat. Its
            %%% own entry P(c,c) > 0 fixes the sign. Cu(001) at 0 deg (eigenspace y-z)
            %%% gives (0,1,0) then (0,0,1).
            P = E * E.';
            P = (P + P.') / 2;
            g = size(E, 2);
            B = zeros(3, g);
            for k = 1:g
                nrm = sqrt(sum(P.^2, 1));
                c = find(nrm >= max(nrm) * (1 - 1e-10), 1);
                b = P(:, c) / nrm(c);
                B(:, k) = b;
                P = P - b * b.';
            end
        end

    end

    methods (Static, Access = private)

        function checkTensorRho(rotElTensor, rho, methodName)
            %%% A real 3x3x3x3 tensor and a positive scalar density
            validateattributes(rotElTensor, {'double'}, {'size',[3 3 3 3],'real','finite'}, ...
                ['LimitingVels.' methodName], 'rotElTensor');
            validateattributes(rho, {'numeric'}, {'scalar','real','finite','positive'}, ...
                ['LimitingVels.' methodName], 'rho');
        end

        function [vLim, thetaRoot, residual, ok, pol] = refineAt(normC, lo, hi, iwt, relTol)
            %%% Refine one root of |vg2| in [lo, hi] and apply the accept test
            %%% (shared by findRoots and refineRoot, so the two cannot drift apart)
            %%%
            %%% vLim = vg1 = v cos(theta) - dv/dtheta sin(theta) at the root; accepted when
            %%% |vg2| < relTol x |vg1| and vg1 > 0. pol is the [3x1] polarisation of wave
            %%% type iwt at the root, from the same Christoffel solve as vLim.
            [thetaRoot, residual] = fminbnd(@(th) LimitingVels.vg2Abs(normC, th, iwt), ...
                lo, hi, LimitingVels.FMIN_OPTS);
            [vRoot, dvRoot, eigVecsRoot] = Christoffel.phaseVelAngleDeriv(normC, thetaRoot, 1);
            pol  = eigVecsRoot(:, iwt);
            vLim = vRoot(iwt) * cos(thetaRoot) - dvRoot(iwt) * sin(thetaRoot);
            ok = residual < relTol * abs(vLim) && vLim > 0;
        end

        function eigenVecs_sA = canonicaliseDegenerate(normC, vLim_sA, eigenVecs_sA, thetaRoots_sA)
            %%% Give the roots of a degenerate group a canonical basis of their eigenspace
            %%% >> eigenVecs_sA = LimitingVels.canonicaliseDegenerate(normC, vLim_sA, eigenVecs_sA, thetaRoots_sA)
            %%%
            %%% A degenerate group: roots of different wave types at one velocity and one
            %%% theta* (DEGENERATE_REL_TOL / DEGENERATE_THETA_TOL) whose phase velocities
            %%% coincide there too. Their eigenvalue is multiple, so each root's own
            %%% eigenvector is an arbitrary vector of the 2-D eigenspace, decided by
            %%% round-off: Cu(001) at 0 deg returned the SAME vector (0,1,0) for both
            %%% shear roots, MATLAB and Julia return different ones. Every polarisation
            %%% consumer (the limiting-mode association, ModeAssociation's anchors) then
            %%% inherits the coin flip. The group's polarisations are replaced by
            %%% canonicalBasis of the eigenspace at the first member's theta*, in wave-type
            %%% order: the same vectors in either language, and a basis that spans the
            %%% eigenspace (which the two raw vectors need not). Non-degenerate roots are
            %%% untouched.
            iwtAll = []; rowAll = []; vAll = []; thAll = [];
            for iwt = 1:3
                n = numel(vLim_sA{iwt});
                iwtAll = [iwtAll; iwt * ones(n, 1)];       %#ok<AGROW>
                rowAll = [rowAll; (1:n).'];                %#ok<AGROW>
                vAll   = [vAll; vLim_sA{iwt}(:)];          %#ok<AGROW>
                thAll  = [thAll; thetaRoots_sA{iwt}(:)];   %#ok<AGROW>
            end
            nAll = numel(vAll);
            done = false(nAll, 1);
            for a = 1:nAll
                if done(a); continue; end
                grp = a;
                for b = a+1:nAll
                    if ~done(b) && ~any(iwtAll(grp) == iwtAll(b)) ...
                            && abs(vAll(b) - vAll(a)) <= LimitingVels.DEGENERATE_REL_TOL * abs(vAll(a)) ...
                            && abs(thAll(b) - thAll(a)) <= LimitingVels.DEGENERATE_THETA_TOL
                        grp(end+1) = b; %#ok<AGROW>
                    end
                end
                if numel(grp) < 2; continue; end
                [vPh, ~, E] = Christoffel.phaseVelAngleDeriv(normC, thAll(a), 1);
                iw = iwtAll(grp);
                if max(vPh(iw)) - min(vPh(iw)) > LimitingVels.DEGENERATE_REL_TOL * max(vPh(iw))
                    continue   % same vg1 and theta*, but not a multiple eigenvalue
                end
                done(grp) = true;
                B = LimitingVels.canonicalBasis(E(:, iw));
                for k = 1:numel(grp)
                    eigenVecs_sA{iwtAll(grp(k))}(rowAll(grp(k)), :) = B(:, k).';
                end
            end
        end

        function vg2 = vg2All(normC, thetaRad)
            %%% Out-of-plane (x3) group-velocity component for all three wave
            %%% types at one sagittal-plane angle:
            %%%     vg2(theta) = v*sin(theta) + dv/dtheta * cos(theta)
            %%%
            %%% Both outputs must be captured on their own line — MATLAB
            %%% silently keeps only the first output of a multi-output call
            %%% used inside a nested expression.
            [v, dv] = Christoffel.phaseVelAngleDeriv(normC, thetaRad, 1);
            vg2 = v .* sin(thetaRad) + dv .* cos(thetaRad);
        end

        function y = vg2Abs(normC, thetaRad, iwt)
            %%% fminbnd objective for findRoots: |vg2| for one wave type.
            vg2 = LimitingVels.vg2All(normC, thetaRad);
            y = abs(vg2(iwt));
        end

    end
end
