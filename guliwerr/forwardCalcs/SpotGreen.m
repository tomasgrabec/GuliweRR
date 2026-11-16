classdef SpotGreen
%%% SpotGreen — the half-space Green's function as a finite measurement sees it: the
%%% coherent sum of the complex Green's function over the wavevectors a ProbeKernel
%%% couples to. A picture tool (CalcAngSweepHS.calcSpotGreenMap,
%%% CalcSingleHS.calcSpotGreen); it never changes a classification.
%%% EXPERIMENTAL: tried on Fe3Al(316) only, not validated against measurements; the
%%% interface may change.
%%%
%%% For one displayed direction theta and kernel nodes (offset o_j, |k|/k0 - 1 = delta_m)
%%% with weights w_mj (ProbeKernel.sample):
%%%
%%%     X(v) = sum_j sum_m w_mj R(theta + o_j, v / (1 + delta_m))
%%%
%%% R(phi, v) = the one-sided Green component in phi's own propagation frame at |k| = 1.
%%% Directions are solved on [0, 180) only: phi + 180 is the same solve read through
%%% reciprocity, R(phi + 180) = s_I s_J G_JI(phi) = s_I s_J G_minus(phi), s = (-1, -1, 1)
%%% (a half turn flips x1 and x2), so a +-k symmetric kernel's two lobes give the
%%% Maznev +-k combination: the point kernel with 'G13' gives G13_eff exactly.
%%%
%%% The delta sum resamples each direction's spectrum at v / (1 + delta) by Keys cubic
%%% convolution on the uniform solve axis (local, exact at the nodes, coefficients
%%% shared by every direction): no extra Green solve for the |k| spread.
%%%
%%% Axes: every direction is solved on ONE structure-wide axis VEL_MIN + n x velStep
%%% (from below VEL_MIN / (1 + delta_max) to above vEnd / (1 + delta_min), vEnd =
%%% 1.05 x the largest bulk longitudinal velocity over in-plane directions -- an upper
%%% bound of every direction's longLimVel), and X is reported on GreenScore.gridAxis(vEnd,
%%% velStep). The contributions to one angle are added in a fixed order (sorted
%%% direction, then +k before -k lobe, then node), so a single angle equals the same
%%% angle in any sweep -- bit for bit as far as the BLAS library rounds one column of a
%%% matrix product independently of the others (true on the machines tested; otherwise
%%% to single-precision round-off).
%%%
%%% Cost and precision: the direction solves are cheap (~0.3 s for a 360-direction SRAS
%%% map); the spot sum dominates. Per direction and lobe, the resampled block T is formed
%%% once and multiplied by the weight columns of every (angle, node) pair that reads it --
%%% and the two quadrature checks' columns -- in ONE matrix product, in single precision
%%% (X is a picture stored as single anyway: |X| changes by ~1e-6 of its peak against a
%%% double-precision sum; the direction solves stay double).
%%%
%%% Quadrature check, for free: the directions on the doubled step and the |k| nodes on
%%% the doubled step are subsets of those already used; X is re-summed on each subset
%%% (weights renormalised) and the largest change of |X|, relative to max |X| of that
%%% angle, is reported (info.dirConvergence / relKConvergence). A large dirConvergence
%%% is the "comb": the branch moves further between neighbouring directions than a
%%% line is wide, so the sum shows separate poles instead of a smooth ridge.
%%%
%%% Static methods only. Half-space only.
%%%
%%% GuliweRR v6.0

    properties (Constant)
        VEL_STEP   = 0.004   % km/s -- default solve and output axis step
        END_FACTOR = 1.05    % axis end: END_FACTOR x max bulk longitudinal velocity
        CONV_WARN  = 0.05    % warn if a quadrature check changes |X| by more than this
        ETA_STEPS  = 4       % default |k| step at most Eta / ETA_STEPS: neighbouring
                             % nodes' poles (v x step apart) must overlap (pole full
                             % width ~ Eta x v), or the |k| sum shows them as a comb;
                             % at Eta/4 the half-density check (Eta/2) stays ~0.1 %
        MAX_DIRECTIONS = 20000   % error above this many direction solves
    end

    methods (Static)
        function opts = options(kernel, varargin)
        %%% Name-Value options of calcSpotGreenMap / calcSpotGreen, defaults filled in
        %%% >> opts = SpotGreen.options(kernel, 'DirStep', 0.1, 'RelKStep', 1e-3, 'Eta', 3e-3, 'VelStep', 0.004)
        %%%
        %%% Defaults: DirStep = kernel.dirStep; RelKStep = min(kernel.relKStep, Eta /
        %%% ETA_STEPS); Eta = SawSpectrum.ETA (the damping of the stored Green maps); VelStep
        %%% = VEL_STEP. An explicit RelKStep is used as given.
            ip = inputParser;
            ip.FunctionName = 'SpotGreen.options';
            ip.addParameter('DirStep', kernel.dirStep, ...
                @(x) validateattributes(x, {'double'}, {'scalar', 'real', 'finite', 'positive'}));
            ip.addParameter('RelKStep', [], ...
                @(x) validateattributes(x, {'double'}, {'scalar', 'real', 'finite', 'positive'}));
            ip.addParameter('Eta', SawSpectrum.ETA, ...
                @(x) validateattributes(x, {'double'}, {'scalar', 'real', 'finite', 'positive'}));
            ip.addParameter('VelStep', SpotGreen.VEL_STEP, ...
                @(x) validateattributes(x, {'double'}, {'scalar', 'real', 'finite', 'positive', '<=', GreenScore.MAX_STEP}));
            ip.parse(varargin{:});
            r = ip.Results;
            if isempty(r.RelKStep); r.RelKStep = min(kernel.relKStep, r.Eta / SpotGreen.ETA_STEPS); end
            opts = struct('dirStep', r.DirStep, 'relKStep', r.RelKStep, 'eta', r.Eta, 'velStep', r.VelStep);
        end

        function [absX, vOut, info] = compute(structureBase, angles, kernel, opts, useParallel)
        %%% |X| for a set of displayed directions
        %%% >> [absX, vOut, info] = SpotGreen.compute(structureBase, angles, kernel, opts, useParallel)
        %%%
        %%% INPUT:  1. structureBase .. UNROTATED half-space Structure
        %%%         2. angles        .. displayed directions [deg]
        %%%         3. kernel        .. ProbeKernel
        %%%         4. opts          .. SpotGreen.options(kernel, ...)
        %%%  (opt)  5. useParallel   .. "auto" (default) | "always" | "never" -- thread pool
        %%%                             over the direction solves (the result does not
        %%%                             depend on it)
        %%% OUTPUT: absX .. [nVel x nAngles] single, |X|
        %%%         vOut .. [nVel x 1] the output axis [km/s]
        %%%         info .. struct: kernel, the opts, vEnd, nDirections, nNodes ([nD nO]),
        %%%                 dirConvergence / relKConvergence ([nAngles x 1], NaN where the
        %%%                 kernel has a single node in that coordinate)
            if nargin < 5 || isempty(useParallel); useParallel = "auto"; end
            validateattributes(structureBase, {'Structure'}, {'scalar'}, 'SpotGreen.compute', 'structureBase');
            if structureBase.nLay ~= 1
                error('SpotGreen:compute:layered', ...
                    'SpotGreen is for a half-space (nLay = 1): on a layered stack |k| sets kd.');
            end
            validateattributes(angles, {'double'}, {'vector', 'real', 'finite'}, 'SpotGreen.compute', 'angles');
            validateattributes(kernel, {'ProbeKernel'}, {'scalar'}, 'SpotGreen.compute', 'kernel');
            angles = angles(:);
            nA = numel(angles);

            q = kernel.sample(opts.dirStep, opts.relKStep);
            nO = numel(q.offsets); nD = numel(q.delta);
            vEnd = SpotGreen.velocityEnd(structureBase);
            vOut = GreenScore.gridAxis(vEnd, opts.velStep);
            vS   = SpotGreen.solveAxis(vOut(end), q.delta, opts.velStep);
            M    = SpotGreen.resampleMatrix(vOut, q.delta, vS(1), opts.velStep, numel(vS));

            %%% Node bookkeeping: every (angle, node direction) pair -> folded direction key
            %%% and lobe (phi in [180, 360) mod 360 reads the -k side of the solve)
            phi = angles + q.offsets;                                   % [nA x nO]
            key180 = int64(round(phi * 1e6));
            key = mod(key180, int64(180e6));
            flip = mod(floor(double(key180) / 180e6), 2) == 1;
            [uKey, ~, dirOf] = unique(key(:));
            if numel(uKey) > SpotGreen.MAX_DIRECTIONS
                error('SpotGreen:compute:tooManyDirections', ...
                    ['%d direction solves needed (limit %d): use a coarser DirStep, fewer ' ...
                     'angles, or angles on the DirStep grid so that they share directions.'], ...
                    numel(uKey), SpotGreen.MAX_DIRECTIONS);
            end
            dc = SpotGreen.solveDirections(structureBase, double(uKey) / 1e6, vS, opts.eta, ...
                kernel.component, useParallel);

            [I, J] = SpotGreen.componentIndices(kernel.component);
            s = [-1, -1, 1];
            sIJ = s(I) * s(J);

            %%% Weights of the two quadrature checks (subsets, renormalised)
            evenO = mod(q.dirIndex, 2) == 0;
            evenD = mod(q.kIndex, 2) == 0;
            wDir = q.w .* evenO; wDir = wDir / sum(abs(wDir(:)));
            wK   = q.w .* evenD; wK   = wK / sum(abs(wK(:)));
            checkDir = any(q.dirIndex ~= 0);
            checkK   = any(q.kIndex ~= 0);

            %%% single precision from here on (see the class header): the weights once,
            %%% each T as it is formed
            wS = single(q.w); wDirS = single(wDir); wKS = single(wK);
            X  = complex(zeros(numel(vOut), nA, 'single'));
            XD = X; XK = X;
            [aIdx, jIdx] = ndgrid(1:nA, 1:nO);
            %%% fixed accumulation order: direction, lobe (+k first), angle, node
            ord = sortrows([dirOf(:), flip(:), aIdx(:), jIdx(:)]);
            grpStart = [1; find(any(diff(ord(:, 1:2), 1, 1), 2)) + 1; size(ord, 1) + 1];
            for g = 1:numel(grpStart) - 1
                rows = ord(grpStart(g):grpStart(g+1)-1, :);
                d = rows(1, 1);
                if rows(1, 2)
                    R = sIJ * dc.Gm(:, d);
                else
                    R = dc.Gp(:, d);
                end
                T = single(reshape(M * R, numel(vOut), nD));
                %%% every (angle, node) pair reading T, plus the checks' columns: one product
                as = rows(:, 3); js = rows(:, 4); nR = numel(js);
                dSel = false(nR, 1);
                if checkDir; dSel = evenO(js(:)); end
                W = wS(:, js);
                if checkDir; W = [W, wDirS(:, js(dSel))]; end %#ok<AGROW>
                if checkK;   W = [W, wKS(:, js)];         end %#ok<AGROW>
                P = T * W;
                for r = 1:nR
                    X(:, as(r)) = X(:, as(r)) + P(:, r);
                end
                c = nR;
                if checkDir
                    aD = as(dSel);
                    for r = 1:numel(aD)
                        XD(:, aD(r)) = XD(:, aD(r)) + P(:, c + r);
                    end
                    c = c + numel(aD);
                end
                if checkK
                    for r = 1:nR
                        XK(:, as(r)) = XK(:, as(r)) + P(:, c + r);
                    end
                end
            end

            absX = abs(X);
            peak = double(max(absX, [], 1).');
            dirConv = nan(nA, 1); relKConv = nan(nA, 1);
            if checkDir; dirConv  = double(max(abs(absX - abs(XD)), [], 1).') ./ peak; end
            if checkK;   relKConv = double(max(abs(absX - abs(XK)), [], 1).') ./ peak; end
            info = struct('kernel', kernel, 'dirStep', q.dirStep, 'relKStep', q.relKStep, ...
                'eta', opts.eta, 'velStep', opts.velStep, 'vEnd', vEnd, ...
                'nDirections', numel(uKey), 'nNodes', [nD, nO], ...
                'dirConvergence', dirConv, 'relKConvergence', relKConv);
            worst = max([dirConv; relKConv; 0], [], 'omitnan');
            if worst > SpotGreen.CONV_WARN
                warning('SpotGreen:notConverged', ...
                    ['The quadrature check changes |X| by up to %.0f%% of its peak (direction ' ...
                     'step %g deg: %.0f%%, |k| step %g: %.0f%%) -- refine DirStep/RelKStep.'], ...
                    100 * worst, q.dirStep, 100 * max([dirConv; 0], [], 'omitnan'), ...
                    q.relKStep, 100 * max([relKConv; 0], [], 'omitnan'));
            end
        end

        function vEnd = velocityEnd(structureBase)
        %%% 1.05 x the largest bulk longitudinal velocity over in-plane directions [km/s]
        %%% >> vEnd = SpotGreen.velocityEnd(structureBase)
        %%%
        %%% The longitudinal limiting velocity of a direction never exceeds the bulk
        %%% longitudinal velocity along it (the limiting wave's slowness projection is the
        %%% largest in the sagittal plane), so this bounds every direction's longLimVel.
        %%% One Christoffel eigenproblem per half degree; a structure property, so a
        %%% single angle and any sweep share it.
            mat = structureBase.mater(structureBase.nLay);
            C = mat.GetStandardElastTensor(0);
            phi = 0:0.5:179.5;
            vMax = 0;
            for i = 1:numel(phi)
                ev = Christoffel.matrix(C, [cosd(phi(i)); sind(phi(i)); 0], mat.rho);
                vMax = max(vMax, sqrt(max(real(ev))));
            end
            vEnd = SpotGreen.END_FACTOR * vMax;
        end

        function vS = solveAxis(vOutEnd, delta, velStep)
        %%% The structure-wide solve axis VEL_MIN + n x velStep covering every v / (1 + delta)
        %%% >> vS = SpotGreen.solveAxis(vOutEnd, delta, velStep)
            v0 = SawSpectrum.VEL_MIN;
            nLo = floor((v0 / (1 + max(delta)) - v0) / velStep) - 2;
            nHi = ceil((vOutEnd / (1 + min(delta)) - v0) / velStep) + 3;
            vS = SawSpectrum.latticeAxis(nLo, nHi, velStep);
        end

        function M = resampleMatrix(vOut, delta, vS1, velStep, nS)
        %%% Sparse Keys-cubic resampling: rows (m-1) nOut + i give R(vOut(i) / (1 + delta(m)))
        %%% >> M = SpotGreen.resampleMatrix(vOut, delta, vS1, velStep, nS)
        %%%
        %%% OUTPUT: M .. sparse [nOut nD x nS], four entries per row (one where the point
        %%%              falls on a solve node, as for delta = 0)
            nOut = numel(vOut); nD = numel(delta);
            vq = vOut(:) ./ (1 + delta(:).');                 % [nOut x nD]
            t = (vq - vS1) / velStep;                         % fractional solve index, 0-based
            tr = round(t);
            onNode = abs(t - tr) < 1e-9;
            t(onNode) = tr(onNode);
            i0 = floor(t(:)); f = t(:) - i0;                  % [nOut nD x 1], m-major blocks
            %%% Keys (a = -1/2) weights of the neighbours i0-1, i0, i0+1, i0+2
            c = [(-0.5*f + f.^2 - 0.5*f.^3), (1 - 2.5*f.^2 + 1.5*f.^3), ...
                 (0.5*f + 2*f.^2 - 1.5*f.^3), (-0.5*f.^2 + 0.5*f.^3)];
            cols = i0 + (0:3);                                % 0-based i0-1..i0+2 -> 1-based
            if any(cols(:) < 1) || any(cols(:) > nS)
                error('SpotGreen:resampleMatrix:outOfAxis', 'The solve axis does not cover the resampling.');
            end
            rowsIdx = repmat((1:nOut*nD)', 1, 4);
            keep = c ~= 0;
            M = sparse(rowsIdx(keep), cols(keep), c(keep), nOut * nD, nS);
        end

        function dc = solveDirections(structureBase, dirs, vS, eta, component, useParallel)
        %%% One Green solve per direction: one-sided G_IJ(+k) and G_IJ(-k) = G_JI(+k)
        %%% >> dc = SpotGreen.solveDirections(structureBase, dirs, vS, eta, component, useParallel)
        %%%
        %%% OUTPUT: dc .. struct: dir [nDir x 1] deg, Gp / Gm [nS x nDir] complex, vS, eta,
        %%%               component
        %%%
        %%% Each direction is solved alone (angle 0 on its own rotated copy -- the
        %%% double-rotation rule), so its spectrum does not depend on which other
        %%% directions were solved.
            [I, J] = SpotGreen.componentIndices(component);
            useP = ~ParallelPolicy.isOnWorker() && ...
                CalcCommon.resolveParallel(string(useParallel), numel(dirs));
            sol = CalcCommon.forEachElement(num2cell(dirs(:)), useP, ...
                @(x, ~) SpotGreen.solveOne(structureBase, x, vS, eta, I, J));
            nS = numel(vS); nDir = numel(dirs);
            Gp = complex(zeros(nS, nDir)); Gm = Gp;
            for i = 1:nDir
                Gp(:, i) = sol{i}.Gp; Gm(:, i) = sol{i}.Gm;
            end
            dc = struct('dir', dirs(:), 'Gp', Gp, 'Gm', Gm, 'vS', vS, 'eta', eta, ...
                'component', char(component));
        end

        function s = solveOne(structureBase, x, vS, eta, I, J)
        %%% G_IJ(+k) and G_IJ(-k) of one direction on the solve axis
        %%% >> s = SpotGreen.solveOne(structureBase, x, vS, eta, I, J)
            sr = CalcCommon.rotateStructure(structureBase, x);
            g = GreenSolver.calc(sr, 0, vS(:).', I, J, 'Eta', eta);
            s = struct('Gp', g.G(:), 'Gm', g.G_minus(:));
        end

        function [I, J] = componentIndices(component)
        %%% 'GIJ' -> I, J
            component = ProbeKernel.validateComponent(component);
            I = component(2) - '0';
            J = component(3) - '0';
        end
    end
end
