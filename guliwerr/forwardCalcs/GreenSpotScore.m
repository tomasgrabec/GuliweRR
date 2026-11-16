classdef GreenSpotScore
%%% GreenSpotScore — the refiner of the half-space 'full' classification: it models
%%% the finite probe spot by summing G13_eff coherently over
%%% a Gaussian spread of directions, then decides existence and velocity from the
%%% late-time spectra of that sum. Run on request (classify('refined')), stored
%%% beside the default GreenScore classification, never inside it.
%%%
%%% EXPERIMENTAL: checked against the measured maps published in Grabec et al. 2024 and the
%%% UTGS paper (arXiv:2510.10696); the scores, their scales and the interface may change.
%%%
%%% For one propagation direction theta:
%%%   1. directions theta + o, o = -2.5:0.5:2.5 deg (FWHM 2.5 deg, cut at 2.5 sigma),
%%%      each solved once (GreenSolver, ETA = 7e-3) on its own axis 0.3 km/s + n x
%%%      VEL_STEP up to 1.05 x its own longLimVel + 0.1 km/s; Q(theta + 180) = Q(theta),
%%%      so directions are folded into [0, 180);
%%%   2. X = sum_j w_j G13_eff(theta + o_j), complex (coherent: the spot excites all
%%%      directions at once), on the shortest of the 11 axes;
%%%   3. late-time |spectra| of X with the gate at GATE_DETECT (4) and GATE_POS (2)
%%%      periods of vRef = sqrt(G_Hill/rho) (SawSpectrum.lateTime: pad 4 -> 2 m/s, tail
%%%      taper 0.2), candidates = findpeaks below longLimVel with half-prominence
%%%      width <= SawSpectrum.MAX_REL_WIDTH x velocity, kept at a reported velocity <= vTf (1 +
%%%      SEARCH_HI_TOL) (vTf = fastest shear limit, GreenScore's cap);
%%%   4. 'high': scoreHigh = log10(prom4 / fU1) >= THRESHOLD_HIGH, fU1 = median |X|
%%%      (un-gated) over [BAND_LO vS, BAND_HI vTf];
%%%      'low' : scoreLow = log10(prom4 / ref.fT12) >= THRESHOLD_LOW, fT12 = median over
%%%      REF_ANGLES of each direction's largest gate-4 prominence (reference(),
%%%      memoised per structure);
%%%   5. velocity: the gate-2 peak whose half-prominence half-width contains the gate-4
%%%      peak (nearest if several), else the gate-4 position (velDet). The early,
%%%      shorter gate keeps more of the pole's own lineshape, so it places the peak
%%%      where the measured map does.
%%%
%%% Why (docs/reference/architecture/forwardCalcs_modeClassification.md, section 3): it uses the
%%% spot-averaged velocity, which sits closer to where a measured map puts the peak, and
%%% a modified accept decision (the averaged spectrum rejects peaks the spot dephases).
%%% Cost: twice the Green solves of a 1-deg default map (one direction every half
%%% degree), so it is a refiner, not the default.
%%%
%%% Single angle = the same angle in any sweep, bit for bit: every direction is
%%% solved alone on an axis set by its own longLimVel, so a direction's spectrum does
%%% not depend on which other directions were solved.
%%%
%%% Value class with Constant properties and Static methods only (plus reference()'s
%%% per-session memo, SawSpectrum.memo). Half-space only.
%%%
%%% GuliweRR v6.0

    properties (Constant)
        ETA            = 7e-3    % GreenSolver damping (calibrated with the rest; > SawSpectrum.ETA)
        VEL_STEP       = 0.008   % km/s -- Green axis step of every direction
        FWHM           = 2.5     % deg -- Gaussian direction spread of the probe spot
        FWHM_TO_SIGMA  = 2.3548  % FWHM / sigma (the literal the thresholds were calibrated with)
        TRUNC_SIGMA    = 2.5     % kernel cut, sigmas
        DIR_STEP       = 0.5     % deg -- kernel direction spacing
        GATE_DETECT    = 4       % late-time gate of the detection spectrum, periods of vRef
        GATE_POS       = 2       % late-time gate of the position spectrum, periods of vRef
        BAND_LO        = GreenScore.BAND_LO   % fU1 band start, x vS (as GreenScore)
        BAND_HI        = 1.00    % fU1 band end, x vTf
        SEARCH_HI_TOL  = GreenScore.SEARCH_HI_TOL   % candidates kept at vel <= vTf (1 + SEARCH_HI_TOL)
        THRESHOLD_HIGH = -0.125  % 'high' standard (scoreHigh)
        THRESHOLD_LOW  = -0.550  % 'low' standard (scoreLow, structure-relative)
        REF_ANGLES     = 0:15:165   % deg -- directions of the 'low' reference
        REF_CACHE_SIZE = 32      % reference() results memoised per session
        LABEL_TOL      = 7e-3    % Rayleigh label below minShearSpot x (1 + LABEL_TOL) (= ETA:
                                 % the pole cannot be placed more finely than its width);
                                 % minShearSpot = the spot-weighted minShear (average())
        POLE_WIN_REL   = 0.012   % identify: a spot peak's single-direction pole lies within
        POLE_WIN_ABS   = 0.030   % max(POLE_WIN_REL v, POLE_WIN_ABS) km/s -- the pairing
                                 % tolerance of refinedSawAgreement
    end

    methods (Static)
        function [o, w] = kernel()
        %%% Direction offsets [deg] and Gaussian weights (summing to 1) of the probe spot
        %%% >> [o, w] = GreenSpotScore.kernel()
            sig = GreenSpotScore.FWHM / GreenSpotScore.FWHM_TO_SIGMA;
            m = floor(GreenSpotScore.TRUNC_SIGMA * sig / GreenSpotScore.DIR_STEP + 1e-9);
            o = (-m:m) * GreenSpotScore.DIR_STEP;
            w = exp(-o.^2 / (2 * sig^2)); w = w / sum(w);
        end

        function [x, key] = foldDirection(x)
        %%% Directions folded into [0, 180) deg (Q(theta + 180) = Q(theta)) and their keys
        %%% >> [x, key] = GreenSpotScore.foldDirection(x)
        %%%
        %%% OUTPUT: x   .. folded, rounded to 1e-6 deg
        %%%         key .. int64 round(x x 1e6), the direction cache key
            key = mod(int64(round(x * 1e6)), int64(180e6));
            x = double(key) / 1e6;
        end

        function D = directionSet(angles)
        %%% Every folded direction the kernels of the given angles need, unique, column
        %%% >> D = GreenSpotScore.directionSet(angles)
            o = GreenSpotScore.kernel();
            [~, key] = GreenSpotScore.foldDirection(angles(:) + o);
            D = double(unique(key(:))) / 1e6;
        end

        function d = solveOne(structureBase, x)
        %%% One direction: its own axis and the complex G13_eff on it
        %%% >> d = GreenSpotScore.solveOne(structureBase, x)
        %%%
        %%% INPUT:  structureBase .. UNROTATED half-space Structure; x .. direction [deg]
        %%% OUTPUT: d .. struct G ([n x 1] complex G13_eff at ETA), vL (longLimVel),
        %%%              vS (min(shearLimVels), the direction's minShear), n
        %%%
        %%% The axis is GreenScore.gridAxis(1.05 vL, VEL_STEP): built as 0.3 + n x step,
        %%% so every direction's axis is the leading part of a longer one.
            sr = CalcCommon.rotateStructure(structureBase, x);
            L = CalcCommon.limitingAndBulk(sr.mater(sr.nLay));
            v = GreenScore.gridAxis(SawSpectrum.AXIS_HI_X_VL * L.longLimVel, GreenSpotScore.VEL_STEP);
            g = GreenSolver.calc(sr, 0, v.', 1, 3, 'Eta', GreenSpotScore.ETA);
            d = struct('G', g.G_eff(:), 'vL', L.longLimVel, 'vS', L.vS, 'n', numel(v));
        end

        function dc = solveDirections(structureBase, D, useParallel, dc)
        %%% Direction cache: solve every direction of D not already in dc
        %%% >> dc = GreenSpotScore.solveDirections(structureBase, D)
        %%% >> dc = GreenSpotScore.solveDirections(structureBase, D, useParallel, dc)
        %%%
        %%% INPUT:  1. structureBase .. UNROTATED half-space Structure
        %%%         2. D             .. directions [deg] (folded here)
        %%%  (opt)  3. useParallel   .. "auto" (default) | "always" | "never" (thread pool
        %%%                             over the directions; the result does not depend on it)
        %%%  (opt)  4. dc            .. an existing cache to extend
        %%% OUTPUT: dc .. struct: key (int64 [nD x 1]), dir, vL, vS (minShear), n, G
        %%%               ([nMax x nD] complex, NaN-padded), v (the common axis, [nMax x 1])
            if nargin < 3 || isempty(useParallel); useParallel = "auto"; end
            if nargin < 4 || isempty(dc)
                dc = struct('key', zeros(0, 1, 'int64'), 'dir', zeros(0, 1), 'vL', zeros(0, 1), ...
                    'vS', zeros(0, 1), 'n', zeros(0, 1), 'G', complex(zeros(0, 0)), 'v', zeros(0, 1));
            end
            [x, key] = GreenSpotScore.foldDirection(D(:));
            [key, iu] = unique(key); x = x(iu);
            new = ~ismember(key, dc.key);
            x = x(new); key = key(new);
            if isempty(key); return; end
            useP = ~ParallelPolicy.isOnWorker() && ...
                CalcCommon.resolveParallel(string(useParallel), numel(x));
            sol = CalcCommon.forEachElement(num2cell(x), useP, ...
                @(xi, ~) GreenSpotScore.solveOne(structureBase, xi));
            nNew = cellfun(@(s) s.n, sol);
            nMax = max([size(dc.G, 1); nNew(:)]);
            G = complex(nan(nMax, numel(dc.key) + numel(key)));
            G(1:size(dc.G, 1), 1:numel(dc.key)) = dc.G;
            for i = 1:numel(sol)
                G(1:sol{i}.n, numel(dc.key) + i) = sol{i}.G;
            end
            dc.G   = G;
            dc.key = [dc.key; key];
            dc.dir = [dc.dir; x];
            dc.vL  = [dc.vL; cellfun(@(s) s.vL, sol(:))];
            dc.vS  = [dc.vS; cellfun(@(s) s.vS, sol(:))];
            dc.n   = [dc.n; nNew(:)];
            dc.v   = SawSpectrum.latticeAxis(0, nMax - 1, GreenSpotScore.VEL_STEP);
        end

        function [X, v, minShearSpot] = average(dc, angle)
        %%% The spot-averaged complex G13_eff of one direction, and the spot-weighted minShear
        %%% >> [X, v] = GreenSpotScore.average(dc, angle)
        %%% >> [X, v, minShearSpot] = GreenSpotScore.average(dc, angle)
        %%%
        %%% OUTPUT: X            .. [n x 1] complex, sum_j w_j G13_eff(angle + o_j), on v =
        %%%                         the first n points of the common axis, n = the shortest
        %%%                         kernel axis
        %%%         minShearSpot .. sum_j w_j minShear(angle + o_j) [km/s] -- the boundary the
        %%%                         refined Rayleigh label compares with (identify()): a bound
        %%%                         Rayleigh wave lies below minShear in every direction of the
        %%%                         spot, so its averaged peak lies below the averaged minShear,
        %%%                         while on a curved branch it can lie above minShear(angle)
        %%%
        %%% Summed term by term in kernel order (not as a matrix product), so the result
        %%% does not depend on BLAS blocking.
            [o, w] = GreenSpotScore.kernel();
            [~, key] = GreenSpotScore.foldDirection(angle + o);
            [ok, col] = ismember(key, dc.key);
            if ~all(ok)
                error('GreenSpotScore:average:missingDirection', ...
                    'The direction cache lacks %d of the %d kernel directions of %g deg.', ...
                    nnz(~ok), numel(ok), angle);
            end
            n = min(dc.n(col));
            X = complex(zeros(n, 1));
            for j = 1:numel(col)
                X = X + w(j) * dc.G(1:n, col(j));
            end
            v = dc.v(1:n);
            minShearSpot = 0;
            for j = 1:numel(col)
                minShearSpot = minShearSpot + w(j) * dc.vS(col(j));
            end
        end

        function pk = direction(X, v, vS, vTf, vL, vRef, ref)
        %%% Scored candidates of one direction from its spot-averaged spectrum
        %%% >> pk = GreenSpotScore.direction(X, v, vS, vTf, vL, vRef)
        %%% >> pk = GreenSpotScore.direction(X, v, vS, vTf, vL, vRef, ref)
        %%%
        %%% INPUT:  1-2. X, v      .. average() output
        %%%         3-5. vS, vTf, vL .. slowest / fastest shear and longitudinal limiting
        %%%                         velocity of this direction [km/s]
        %%%         6. vRef        .. sqrt(G_Hill/rho) [km/s]
        %%%  (opt)  7. ref         .. reference(); fills scoreLow
        %%% OUTPUT: pk .. struct array (0..N, one per gate-4 candidate): vel (reported:
        %%%               gate-2-linked), velDet (gate-4), posLinked, prominence, relWidth,
        %%%               height (gate 4), fU1, scoreHigh, scoreLow; identity fields
        %%%               (modeNo, matchScore, confident, type) and accepted empty
        %%%               until identify() / accepted()
            pk = GreenSpotScore.emptyPeaks();
            b = v >= GreenSpotScore.BAND_LO * vS & v <= GreenSpotScore.BAND_HI * vTf;
            fU1 = median(abs(X(b)), 'omitnan');
            %%% gates 4 (detection) and 2 (position) share window, FFT and axis
            [Ag, vO] = SawSpectrum.lateTime(X, v, vRef, [], [], ...
                [GreenSpotScore.GATE_DETECT, GreenSpotScore.GATE_POS]);
            P4 = SawSpectrum.peaksOf(Ag(:, :, 1), vO, vL);
            if isempty(P4.vel); return; end
            P2 = SawSpectrum.peaksOf(Ag(:, :, 2), vO, vL);
            hw = P2.relWidth .* P2.vel / 2;
            fT12 = NaN; if nargin >= 7 && ~isempty(ref); fT12 = ref.fT12; end
            for j = 1:numel(P4.vel)
                vel = P4.vel(j); linked = false;
                dv = abs(P4.vel(j) - P2.vel); ok = find(dv <= hw);
                if ~isempty(ok)
                    [~, mm] = min(dv(ok)); vel = P2.vel(ok(mm)); linked = true;
                end
                pk(j) = struct('vel', vel, 'velDet', P4.vel(j), 'posLinked', linked, ...
                    'prominence', P4.prom(j), 'relWidth', P4.relWidth(j), 'height', P4.height(j), ...
                    'fU1', fU1, 'scoreHigh', log10(P4.prom(j) / fU1), ...
                    'scoreLow', log10(P4.prom(j) / fT12), 'modeNo', NaN, 'matchScore', NaN, ...
                    'confident', false, 'type', '', 'accepted', false);
            end
            %%% Candidates above the fastest shear limit are not surface waves: the same
            %%% cap as GreenScore.candidates. Between vTf and vL the late-time spectrum
            %%% rings at the longitudinal branch point (a surface-skimming bulk wave,
            %%% 0.07-0.2 % below vL); uncapped, such candidates can score close to 'high'
            %%% on strongly anisotropic cuts (computed with GuliweRR). Applied to the
            %%% reported velocity.
            pk = pk([pk.vel] <= vTf * (1 + GreenSpotScore.SEARCH_HI_TOL));
        end

        function pk = identify(pk, ritz, minShearSpot, poleVels)
        %%% Ritz mode identity and Rayleigh/pseudoRayleigh label of every candidate
        %%% >> pk = GreenSpotScore.identify(pk, ritz, minShearSpot)
        %%% >> pk = GreenSpotScore.identify(pk, ritz, minShearSpot, poleVels)
        %%%
        %%% INPUT:  minShearSpot .. the label boundary: the spot-weighted minShear from
        %%%                         average() (CalcSingleHS.calcRefined passes it)
        %%%  (opt)  poleVels     .. the single-direction poles of this angle: the vel of
        %%%                         every GreenScore candidate (greenScorePeaks), accepted
        %%%                         or not; default [] (no pole: the spot rule for every peak)
        %%%
        %%% Mode: the spot peak is mapped to the nearest single-direction pole within
        %%% max(POLE_WIN_REL v, POLE_WIN_ABS), and that pole bound by the 'green' rule
        %%% (ModeClassification.matchRitzToGreenPeak): the wave at this angle is what the
        %%% Ritz basis, solved at this angle only, represents. Why not the spot peak
        %%% itself: the spot average moves it off the single-direction pole by up to
        %%% tens of m/s, as far as the slab-mode spacing, so a rule at the spot velocity
        %%% binds a neighbour of the wave's mode on leaky entries. With no pole in the
        %%% window (rare; mostly on a strongly anisotropic material where the spot
        %%% spectrum forms a comb) the spot rule (ModeClassification.matchSpotPeak)
        %%% decides. Label: 'Rayleigh'
        %%% below minShearSpot x (1 + LABEL_TOL), else 'pseudoRayleigh' -- from the Green
        %%% velocity alone (docs/pitfalls.md MTH2 and MTH3). The spot-weighted boundary,
        %%% not minShear(angle): the spot average lifts a Rayleigh wave hugging minShear on a
        %%% curved branch above minShear(angle).
            if isempty(pk); return; end
            if nargin < 4; poleVels = []; end
            mc = ModeClassification();
            for j = 1:numel(pk)
                v = pk(j).vel;
                [d, i] = min(abs(poleVels(:) - v));
                if ~isempty(d) && d <= max(GreenSpotScore.POLE_WIN_REL * v, GreenSpotScore.POLE_WIN_ABS)
                    [modeNo, score, confident] = ModeClassification.matchRitzToGreenPeak(ritz, ...
                        poleVels(i), mc, false);
                else
                    [modeNo, score, confident] = ModeClassification.matchSpotPeak(ritz, v, mc);
                end
                pk(j).modeNo = modeNo; pk(j).matchScore = score; pk(j).confident = confident;
                if pk(j).vel < minShearSpot * (1 + GreenSpotScore.LABEL_TOL)
                    pk(j).type = 'Rayleigh';
                else
                    pk(j).type = 'pseudoRayleigh';
                end
            end
        end

        function thr = resolveThreshold(t)
        %%% Scale and value of a threshold: 'high' | 'low' | numeric (the 'high' scale)
        %%% >> thr = GreenSpotScore.resolveThreshold(t)
        %%%
        %%% OUTPUT: thr .. struct: scale ('high' -> scoreHigh, 'low' -> scoreLow), value
            [value, scale] = CalcCommon.thresholdValue(t, GreenSpotScore.THRESHOLD_HIGH, ...
                GreenSpotScore.THRESHOLD_LOW, 'GreenSpotScore:resolveThreshold:badName', ...
                'GreenSpotScore.resolveThreshold');
            thr = struct('scale', scale, 'value', value);
        end

        function tf = accepted(pk, threshold)
        %%% Which scored candidates pass a threshold
        %%% >> tf = GreenSpotScore.accepted(pk, threshold)
        %%%
        %%% 'low' is not a subset of 'high': the two divide by different normalisers
        %%% (the direction's own fU1, the structure's fT12).
        %%% A NaN on the compared scale is an error (GreenSpotScore:accepted:noScore):
        %%% scoreLow is NaN when direction() ran without the reference.
            thr = GreenSpotScore.resolveThreshold(threshold);
            if isempty(pk); tf = false(0, 1); return; end
            if strcmp(thr.scale, 'low')
                x = [pk.scoreLow]';
            else
                x = [pk.scoreHigh]';
            end
            if any(isnan(x))
                error('GreenSpotScore:accepted:noScore', ['A candidate has no %s-scale score ' ...
                    '(the ''low'' scale needs direction() with the reference).'], thr.scale);
            end
            tf = x >= thr.value;
        end

        function ref = reference(structure, dc, useParallel)
        %%% The 'low' reference of a structure: fT12, its typical strongest wave
        %%% >> ref = GreenSpotScore.reference(structure)
        %%% >> ref = GreenSpotScore.reference(structure, dc, useParallel)
        %%%
        %%% INPUT:  1. structure   .. UNROTATED half-space Structure
        %%%  (opt)  2. dc          .. a direction cache (solveDirections) of this structure;
        %%%                           directions already in it are not solved again
        %%%                           (default [] = none)
        %%%  (opt)  3. useParallel .. "auto" (default) | "always" | "never" (the solves)
        %%% OUTPUT: ref .. struct: vRef [km/s], fT12 = median over REF_ANGLES (with
        %%%                at least one candidate) of each direction's largest gate-4
        %%%                prominence, pmaxPerAngle, angles
        %%%
        %%% One number per structure, so the 'low' labels never depend on the range
        %%% computed. Up to 12 x 11 = 132 direction solves (~4 s serial, < 1 s on a pool),
        %%% none in a 0:1:180 sweep, whose own directions cover them. Memoised per session
        %%% (keyed by the structure's contents; clearReferenceCache() empties
        %%% it); not memoised inside a pool worker.
            if nargin < 2; dc = []; end
            if nargin < 3 || isempty(useParallel); useParallel = "auto"; end
            onWorker = ParallelPolicy.isOnWorker();
            if ~onWorker
                key = GreenSpotScore.referenceKey(structure);
                ref = GreenSpotScore.referenceCache('get', key);
                if ~isempty(ref); return; end
            end
            angles = GreenSpotScore.REF_ANGLES(:);
            if onWorker; useParallel = "never"; end
            dc = GreenSpotScore.solveDirections(structure, GreenSpotScore.directionSet(angles), ...
                string(useParallel), dc);
            vRef = CalcCommon.surfaceVRef(structure);
            pmax = nan(numel(angles), 1);
            for i = 1:numel(angles)
                sr = CalcCommon.rotateStructure(structure, angles(i));
                L = CalcCommon.limitingAndBulk(sr.mater(sr.nLay));
                [X, v] = GreenSpotScore.average(dc, angles(i));
                [A, vO] = SawSpectrum.lateTime(X, v, vRef, [], [], GreenSpotScore.GATE_DETECT);
                P = SawSpectrum.peaksOf(A, vO, L.longLimVel);
                if ~isempty(P.prom); pmax(i) = max(P.prom); end
            end
            ref = struct('vRef', vRef, ...
                'fT12', median(pmax(isfinite(pmax))), 'pmaxPerAngle', pmax, 'angles', angles);
            if ~onWorker; GreenSpotScore.referenceCache('put', key, ref); end
        end

        function key = referenceKey(structure)
        %%% reference()'s memo key: the structure's physical content (SawSpectrum.structureKey)
        %%% >> key = GreenSpotScore.referenceKey(structure)
            key = SawSpectrum.structureKey(structure);
        end

        function clearReferenceCache()
        %%% Empty reference()'s per-session memo
        %%% >> GreenSpotScore.clearReferenceCache()
            GreenSpotScore.referenceCache('clear');
        end

        function out = referenceCache(action, key, ref)
        %%% reference()'s memo: 'get' (returns [] on a miss) | 'put' | 'clear' | 'size'
        %%% >> ref = GreenSpotScore.referenceCache('get', key)
        %%%
        %%% SawSpectrum.memo under the name 'GreenSpotScore', REF_CACHE_SIZE entries.
            if nargin < 2; key = []; end
            if nargin < 3; ref = []; end
            out = SawSpectrum.memo('GreenSpotScore', GreenSpotScore.REF_CACHE_SIZE, 'GreenSpotScore:referenceCache', ...
                action, key, ref);
        end

        function pk = emptyPeaks()
        %%% The candidate struct with no rows (direction()'s field set)
        %%% >> pk = GreenSpotScore.emptyPeaks()
            pk = struct('vel', {}, 'velDet', {}, 'posLinked', {}, 'prominence', {}, ...
                'relWidth', {}, 'height', {}, 'fU1', {}, 'scoreHigh', {}, 'scoreLow', {}, ...
                'modeNo', {}, 'matchScore', {}, 'confident', {}, 'type', {}, 'accepted', {});
        end
    end
end
