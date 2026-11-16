classdef GreenScore
%%% GreenScore — the default existence score of the half-space 'full' path: each
%%% G13_eff peak is fitted as one damped pole, and its height is corrected
%%% analytically for the late-time gate and for angular-aperture dephasing.
%%%
%%% EXPERIMENTAL: checked against the measured maps published in Grabec et al. 2024 and the
%%% UTGS paper (arXiv:2510.10696); the scores, their scales and the interface may change.
%%%
%%% For one propagation direction, on the band [BAND_LO x vS, BAND_SOLVE_HI x vTf]
%%% (vS / vTf = slowest / fastest shear limiting velocity) of the complex G13_eff
%%% (the +-k symmetrised TGS response at SawSpectrum.ETA):
%%%   1. candidates: local maxima of |G13_eff| at or below vTf(1 + SEARCH_HI_TOL), with
%%%      half-prominence width <= MAX_REL_WIDTH x velocity;
%%%   2. a 3-point complex pole fit (a0 + a1 x)/(x - q) through the peak sample and the
%%%      samples FIT_SPAN either side gives the pole width gamma = -Im q, floored at the
%%%      damping width SawSpectrum.ETA v/2, and the height; a failed fit (no damped pole within
%%%      the stencil) keeps the grid peak height;
%%%   3. score = log10((pole height - base) / fBand)        fBand = median |G13_eff| over
%%%              + log10 exp(-2 pi GATE_PERIODS gamma/vRef)  [BAND_LO vS, vTf]
%%%              + log10 F(gamma, |dv/dtheta| APERTURE_SIGMA)
%%%      The second term is what a late-time gate starting GATE_PERIODS periods of vRef
%%%      in does to a pole of half-width gamma (exact for one pole, no FFT). F is the
%%%      gated height ratio of the same pole when the probe aperture sums a Gaussian
%%%      spread of directions (APERTURE_SIGMA) along a branch of slope dv/dtheta:
%%%      F = rho sqrt(pi/2) erfcx((rho + q)/sqrt(2)) exp(-q^2/2), rho = gamma/dv,
%%%      q = 2 pi GATE_PERIODS dv / vRef. Narrow poles on steep branches lose most.
%%%      The slope comes from two short probe solves at angle +- PROBE_OFFSET, so a
%%%      single angle classifies exactly as inside any sweep, at any angular step. Each
%%%      probe reads the local maximum of |G| nearest the pole (not the stencil's
%%%      arg-max, which jumps to a stronger neighbour's skirt at crossings and cusps);
%%%      the two sides are combined by combineSlopes.
%%%   4. 'high': score >= THRESHOLD_HIGH (everything a good TGS map shows).
%%%      'low' : lowScore = log10(gated, dephased prominence) - reference.logStrengthTop
%%%              >= THRESHOLD_LOW, i.e. at least 10^THRESHOLD_LOW of the structure's
%%%              typical strongest wave (median over N_REF_ANGLES directions of each
%%%              direction's largest; reference(), memoised).
%%%      numeric: score >= value (the 'high' scale).
%%%
%%% Why (details: docs/reference/architecture/forwardCalcs_modeClassification.md): tuned
%%% against measured TGS maps (the published ones are listed in section 11 there), at
%%% a fraction of the cost of the late-time spectra --
%%% the classification needs only the band, any uniform axis step up to MAX_STEP, and
%%% no FFT. The dephasing term is what matches the measured maps on steep branch ends.
%%% The pole fit, the width
%%% floor, the failed-fit rule and the two-sided probe are shared with the layered
%%% score (fitPeaks, probeSlope); they keep strong waves at branch crossings and cusps,
%%% where a probe reading the stencil's arg-max jumps to a stronger neighbour.
%%% Half-space only: layered stacks use GreenScoreLayered (same candidates and probe,
%%% whole-axis band, own-period gate, polarisation term).
%%%
%%% Value class with Constant properties and Static methods only (plus reference()'s
%%% per-session memo, SawSpectrum.memo).
%%%
%%% GuliweRR v6.0

    properties (Constant)
        GATE_PERIODS   = SawSpectrum.GATE_PERIODS   % analytic late-time gate start, periods of
                                 % vRef = sqrt(G_Hill/rho): the same gate as the late-time spectra
                                 % (RitzScore's dephasing term reads it through dephaseFactor)
        APERTURE_SIGMA = 0.4     % deg -- std of the probe's direction spread (dephasing term);
                                 % 0.4 (rather than 0.8) keeps bound Rayleigh waves that hug vS on
                                 % steep branches; the real spread of the spot is not measured
        MAX_REL_WIDTH  = SawSpectrum.MAX_REL_WIDTH   % candidate half-prominence width / velocity ceiling
        BAND_LO        = 0.80    % band start and normaliser start, x vS
        BAND_SOLVE_HI  = 1.03    % band end (candidate search), x vTf
        SEARCH_HI_TOL  = 0.012   % candidates kept at vel <= vTf (1 + SEARCH_HI_TOL)
        FIT_SPAN       = 0.016   % km/s -- pole-fit stencil: samples i0 and i0 +- round(FIT_SPAN/dv)
        PROBE_OFFSET   = 0.5     % deg -- slope probe direction offset (both sides)
        PROBE_STEP     = 0.016   % km/s -- probe stencil spacing
        PROBE_HALF     = 3       % probe stencil half-width, samples (pole +- 48 m/s)
        DEFAULT_STEP   = 0.006   % km/s -- default Green axis step of the Calc classes (pictures)
        MAX_STEP       = 0.016   % km/s -- largest validated axis step (larger steps are capped)
        THRESHOLD_HIGH = 0.16    % 'high' standard (score)
        THRESHOLD_LOW  = -0.79   % 'low' standard (lowScore, structure-relative)
        N_REF_ANGLES   = 12      % directions of the 'low' reference
        REF_CACHE_SIZE = 32      % reference() results memoised per session
    end

    methods (Static)
        function thr = resolveThreshold(t)
        %%% Scale and value of a threshold: 'high' | 'low' | numeric (the 'high' scale)
        %%% >> thr = GreenScore.resolveThreshold(t)
        %%%
        %%% OUTPUT: thr .. struct: scale ('high' -> compared with pk.score, 'low' -> with
        %%%                pk.lowScore), value
            [value, scale] = CalcCommon.thresholdValue(t, GreenScore.THRESHOLD_HIGH, ...
                GreenScore.THRESHOLD_LOW, 'GreenScore:resolveThreshold:badName', 'GreenScore.resolveThreshold');
            thr = struct('scale', scale, 'value', value);
        end

        function tf = accepted(pk, threshold)
        %%% Which scored candidates pass a threshold
        %%% >> tf = GreenScore.accepted(pk, threshold)
        %%%
        %%% INPUT:  1. pk        .. scored candidates (score(), with lowScore() for 'low')
        %%%         2. threshold .. 'high' | 'low' | numeric
        %%% OUTPUT: tf .. logical, one per candidate
        %%%
        %%% A NaN on the compared scale is an error (GreenScore:accepted:noScore), not a
        %%% silent rejection: lowScore is NaN when direction() ran without the reference.
            thr = GreenScore.resolveThreshold(threshold);
            if isempty(pk); tf = false(0, 1); return; end
            if strcmp(thr.scale, 'low')
                x = [pk.lowScore]';
            else
                x = [pk.score]';
            end
            if any(isnan(x))
                error('GreenScore:accepted:noScore', ['A candidate has no %s-scale score ' ...
                    '(the ''low'' scale needs direction() with the reference).'], thr.scale);
            end
            tf = x >= thr.value;
        end

        function step = resolveStep(step, warnId)
        %%% A Green axis step the score is valid for: positive, capped at MAX_STEP
        %%% >> step = GreenScore.resolveStep(step, warnId)
        %%%
        %%% A step above MAX_STEP is replaced by MAX_STEP, with a once-per-session warning
        %%% (warnId; warning('on', warnId) re-enables it).
            validateattributes(step, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
                'GreenScore.resolveStep', 'step');
            if step > GreenScore.MAX_STEP + 1e-12
                if strcmp(warning('query', warnId).state, 'on')
                    warning(warnId, ['Green axis step %.4g km/s is above the largest step the existence ' ...
                        'score (GreenScore) is validated for; using %.4g km/s. (Shown once per session; ' ...
                        'warning(''on'', ''%s'') to re-enable.)'], step, GreenScore.MAX_STEP, warnId);
                    warning('off', warnId);
                end
                step = GreenScore.MAX_STEP;
            end
        end

        function v = gridAxis(vMax, step)
        %%% Calc Green axis 0.3 km/s + n x step up to vMax + SawSpectrum.EDGE_TAPER, column
        %%% >> v = GreenScore.gridAxis(vMax, step)
        %%%
        %%% Built as VEL_MIN + n*step, not with the colon operator (which fills a range
        %%% from both ends), so a short axis is bit-for-bit the leading part of a longer
        %%% one with the same step: a single angle then scores exactly as inside a sweep.
            n = floor((vMax + SawSpectrum.EDGE_TAPER - SawSpectrum.VEL_MIN) / step + 1e-9);
            v = SawSpectrum.latticeAxis(0, n, step);
        end

        function v = bandAxis(vS, vTf, step)
        %%% The band of one direction on the Calc grid 0.3 km/s + n x step, column
        %%% >> v = GreenScore.bandAxis(vS, vTf, step)
            v0 = SawSpectrum.VEL_MIN;
            n1 = floor((GreenScore.BAND_LO * vS - v0) / step);
            n2 = ceil((GreenScore.BAND_SOLVE_HI * vTf - v0) / step);
            v  = SawSpectrum.latticeAxis(max(n1, 0), n2, step);
        end

        function pk = candidates(G, velAxis, vS, vTf)
        %%% Pole-fitted candidates of one direction's complex G13_eff
        %%% >> pk = GreenScore.candidates(G, velAxis, vS, vTf)
        %%%
        %%% INPUT:  1. G       .. [nVel x 1] complex G13_eff at SawSpectrum.ETA
        %%%         2. velAxis .. [nVel x 1] uniform axis [km/s] covering the band (a whole
        %%%                       Calc axis or bandAxis); only the band is read
        %%%         3. vS, vTf .. slowest / fastest shear limiting velocity [km/s]
        %%% OUTPUT: pk .. struct array (0..N), fields vel (log-parabola refined |G| peak,
        %%%               the reported velocity), polVel (fitted pole position), gamma
        %%%               (pole half-width [km/s], >= SawSpectrum.ETA vel/2), poleOk, height, prominence,
        %%%               base, poleHeight (= height when the fit failed), relWidth,
        %%%               relHeight, fBand, rawScore (log10((poleHeight - base)/fBand));
        %%%               probeSlope() and score() add the rest
        %%%
        %%% The search runs over the whole band, past the vTf bound, so that a peak just
        %%% below vTf is a genuine interior maximum (docs/pitfalls.md MTH4). A band of
        %%% fewer than 5 finite samples gives no candidates; if no finite sample lies in
        %%% the normalising band [BAND_LO vS, vTf], fBand is the median over the whole
        %%% band.
            pk = GreenScore.emptyPeaks();
            [found, v, G, s, h, loc, wid, prom] = GreenScore.bandPeaks(G, velAxis, ...
                GreenScore.BAND_LO * vS, GreenScore.BAND_SOLVE_HI * vTf, 'GreenScore:candidates:nonUniform');
            if ~found; return; end
            ok = isfinite(s);
            nb = ok & v >= GreenScore.BAND_LO * vS - 1e-9 & v <= vTf + 1e-9; if ~any(nb); nb = ok; end
            fBand = median(s(nb));
            keep = loc <= vTf * (1 + GreenScore.SEARCH_HI_TOL) & wid ./ loc <= GreenScore.MAX_REL_WIDTH;
            h = h(keep); loc = loc(keep); wid = wid(keep); prom = prom(keep);
            f = GreenScore.fitPeaks(v, G, h, loc, wid);
            hMax = max([h; eps]);
            for j = 1:numel(loc)
                base = h(j) - prom(j);
                pk(j) = struct('vel', f.vel(j), 'polVel', f.polVel(j), 'gamma', f.gamma(j), ...
                    'poleOk', f.poleOk(j), 'height', h(j), 'prominence', prom(j), 'base', base, ...
                    'poleHeight', f.poleHeight(j), 'relWidth', wid(j) / f.vel(j), ...
                    'relHeight', h(j) / hMax, 'fBand', fBand, ...
                    'rawScore', log10(max(f.poleHeight(j) - base, realmin) / fBand), 'slope', NaN, ...
                    'slopeFwd', NaN, 'slopeBwd', NaN, 'gateTerm', NaN, 'dephaseTerm', NaN, ...
                    'score', NaN, 'logStrength', NaN, 'lowScore', NaN);
            end
        end

        function pk = probeSlope(structureBase, angle, pk)
        %%% Branch slope dv/dtheta of every candidate from probes at angle +- PROBE_OFFSET
        %%% >> pk = GreenScore.probeSlope(structureBase, angle, pk)
        %%%
        %%% INPUT:  1. structureBase .. the UNROTATED Structure
        %%%         2. angle         .. this direction [deg]
        %%%         3. pk            .. candidates (fields polVel, poleOk, vel; GreenScore or
        %%%                             GreenScoreLayered)
        %%% OUTPUT: pk with slopeFwd, slopeBwd [km/s per deg] (each NaN if that probe found
        %%%            no interior local maximum) and slope = combineSlopes(slopeFwd,
        %%%            slopeBwd) (NaN if neither did -- score() then charges the stencil
        %%%            bound PROBE_HALF PROBE_STEP / PROBE_OFFSET)
        %%%
        %%% Each probe solves pole +- PROBE_HALF x PROBE_STEP per candidate (shared points
        %%% once), about a dozen points per candidate, and takes the local maximum of |G|
        %%% nearest the stencil centre: the stencil's arg-max would jump to a stronger
        %%% neighbour's skirt next to a crossing or cusp and charge a real wave the
        %%% steepest slope.
            if isempty(pk); return; end
            sF = GreenScore.probeOneSide(structureBase, angle, pk, GreenScore.PROBE_OFFSET);
            sB = GreenScore.probeOneSide(structureBase, angle, pk, -GreenScore.PROBE_OFFSET);
            s2 = GreenScore.combineSlopes(sF, sB);
            for j = 1:numel(pk)
                pk(j).slopeFwd = sF(j); pk(j).slopeBwd = sB(j); pk(j).slope = s2(j);
            end
        end

        function s2 = combineSlopes(sF, sB)
        %%% Two-sided probe slope from the forward and backward probes [km/s per deg]
        %%% >> s2 = GreenScore.combineSlopes(sF, sB)
        %%%
        %%% Same sign and within a factor 2 -> mean; otherwise the gentler one; one side
        %%% NaN -> the other; both NaN -> NaN. (A probe that jumps to another peak only
        %%% ever reads a steeper slope, so the gentler side is the one to trust.)
        %%% Equal and opposite to 1e-9 relative (a mirror direction: the branch is
        %%% stationary there, e.g. Cu(001) 0/45 deg, NiTi/MgO 0 deg) -> the mean |slope|:
        %%% "the gentler side" would pick the sign by round-off. Only |slope| enters
        %%% score(), so this changes no score beyond 1e-9 relative.
            sF = sF(:); sB = sB(:); s2 = sF;
            both  = isfinite(sF) & isfinite(sB);
            agree = both & sign(sF) == sign(sB) & max(abs(sF), abs(sB)) <= 2 * min(abs(sF), abs(sB));
            s2(agree) = (sF(agree) + sB(agree)) / 2;
            mirror = both & ~agree & abs(abs(sF) - abs(sB)) <= 1e-9 * max(abs(sF), abs(sB));
            s2(mirror) = (abs(sF(mirror)) + abs(sB(mirror))) / 2;
            pickB = both & ~agree & ~mirror & abs(sB) < abs(sF); s2(pickB) = sB(pickB);
            onlyB = ~isfinite(sF) & isfinite(sB);      s2(onlyB) = sB(onlyB);
        end

        function pk = score(pk, vRef)
        %%% Gate and dephasing terms, the score and the absolute log strength
        %%% >> pk = GreenScore.score(pk, vRef)
        %%%
        %%% INPUT:  1. pk   .. candidates() after probeSlope()
        %%%         2. vRef .. sqrt(G_Hill/rho) of the surface material [km/s]
        %%% OUTPUT: pk with gateTerm, dephaseTerm, score = rawScore + gateTerm + dephaseTerm
        %%%            and logStrength = score + log10(fBand) (the gated, dephased pole
        %%%            prominence in |G13_eff| units, log10 -- what lowScore compares)
        %%%
        %%% A slope neither probe could determine (no local maximum inside either
        %%% stencil) is set to the stencil's bound: such a candidate is at least that
        %%% steep, and an unknown slope must not escape the dephasing penalty.
            sLB = GreenScore.PROBE_HALF * GreenScore.PROBE_STEP / GreenScore.PROBE_OFFSET;
            for j = 1:numel(pk)
                sl = pk(j).slope; if ~isfinite(sl); sl = sLB; end
                g = pk(j).gamma;
                pk(j).gateTerm = -2 * GreenScore.GATE_PERIODS * pi * g / vRef / log(10);
                F = GreenScore.dephaseFactor(g, abs(sl) * GreenScore.APERTURE_SIGMA, vRef);
                pk(j).dephaseTerm = log10(max(F, realmin));
                pk(j).score = pk(j).rawScore + pk(j).gateTerm + pk(j).dephaseTerm;
                pk(j).logStrength = pk(j).score + log10(pk(j).fBand);
            end
        end

        function pk = lowScore(pk, ref)
        %%% lowScore = logStrength - ref.logStrengthTop (the 'low' scale)
        %%% >> pk = GreenScore.lowScore(pk, ref)
            for j = 1:numel(pk)
                pk(j).lowScore = pk(j).logStrength - ref.logStrengthTop;
            end
        end

        function pk = direction(structureBase, angle, G, velAxis, vS, vTf, vRef, ref)
        %%% All steps for one direction: candidates -> probe slope -> score (-> lowScore)
        %%% >> pk = GreenScore.direction(structureBase, angle, G, velAxis, vS, vTf, vRef)
        %%% >> pk = GreenScore.direction(structureBase, angle, G, velAxis, vS, vTf, vRef, ref)
        %%%
        %%% INPUT:  1. structureBase .. UNROTATED Structure (the probe rotates it itself)
        %%%         2. angle         .. direction [deg]
        %%%         3-4. G, velAxis  .. complex G13_eff at ETA on a uniform axis covering
        %%%                             the band (this direction's rotated structure)
        %%%         5-7. vS, vTf, vRef
        %%%  (opt)  8. ref           .. reference() result; fills lowScore
            pk = GreenScore.candidates(G, velAxis, vS, vTf);
            pk = GreenScore.probeSlope(structureBase, angle, pk);
            pk = GreenScore.score(pk, vRef);
            if nargin >= 8 && ~isempty(ref); pk = GreenScore.lowScore(pk, ref); end
        end

        function ref = reference(structure, nAngles, useParallel)
        %%% The 'low' reference of a structure: its typical strongest wave
        %%% >> ref = GreenScore.reference(structure)
        %%% >> ref = GreenScore.reference(structure, nAngles, useParallel)
        %%%
        %%% INPUT:  1. structure   .. UNROTATED half-space Structure
        %%%  (opt)  2. nAngles     .. equally spaced directions over [0,180) (default N_REF_ANGLES)
        %%%  (opt)  3. useParallel .. "auto" (default) | "always" | "never" (thread pool for
        %%%                           the nAngles directions; the result does not depend on it)
        %%% OUTPUT: ref .. struct: vRef [km/s], logStrengthTop = median over the directions
        %%%                of each direction's largest logStrength, logStrengthPerAngle, angles
        %%%
        %%% One number per structure, shared by every angle of a sweep and by a standalone
        %%% angle, so the 'low' labels never depend on which angles or which range were
        %%% computed. Band solves at DEFAULT_STEP plus the probes (~0.5 s). Memoised per
        %%% session (keyed by the structure's physical content, SawSpectrum.structureKey, and nAngles;
        %%% GreenScore.clearReferenceCache() empties it); not memoised inside a pool worker.
            if nargin < 2 || isempty(nAngles); nAngles = GreenScore.N_REF_ANGLES; end
            if nargin < 3 || isempty(useParallel); useParallel = "auto"; end
            useParallel = string(useParallel);
            if ~ismember(useParallel, ["auto", "always", "never"])
                error('GreenScore:reference:badParallel', ...
                    "useParallel must be 'auto', 'always' or 'never'; got '%s'.", useParallel);
            end
            onWorker = ParallelPolicy.isOnWorker();
            if ~onWorker
                key = {SawSpectrum.structureKey(structure), nAngles};
                ref = GreenScore.referenceCache('get', key);
                if ~isempty(ref); return; end
            end
            angles = (0:nAngles-1)' * 180 / nAngles;
            vRef = CalcCommon.surfaceVRef(structure);
            useP = useParallel ~= "never" && ~onWorker && ParallelPolicy.ensureThreadPool();
            best = CalcCommon.forEachElement(num2cell(angles), useP, ...
                @(a, ~) GreenScore.strongest(structure, a, vRef));
            ref.vRef = vRef;
            ref.logStrengthPerAngle = cell2mat(best(:));
            ref.logStrengthTop = median(ref.logStrengthPerAngle(isfinite(ref.logStrengthPerAngle)));
            ref.angles = angles;
            if ~onWorker; GreenScore.referenceCache('put', key, ref); end
        end

        function clearReferenceCache()
        %%% Empty reference()'s per-session memo
        %%% >> GreenScore.clearReferenceCache()
            GreenScore.referenceCache('clear');
        end

        function out = referenceCache(action, key, ref)
        %%% reference()'s memo: 'get' (returns [] on a miss) | 'put' | 'clear' | 'size'
        %%% >> ref = GreenScore.referenceCache('get', key)
        %%%
        %%% SawSpectrum.memo under the name 'GreenScore', REF_CACHE_SIZE entries.
            if nargin < 2; key = []; end
            if nargin < 3; ref = []; end
            out = SawSpectrum.memo('GreenScore', GreenScore.REF_CACHE_SIZE, 'GreenScore:referenceCache', ...
                action, key, ref);
        end

        function F = dephaseFactor(gamma, dv, vRef)
        %%% Gated height ratio of a pole with a Gaussian velocity spread dv (std) [km/s]
        %%% >> F = GreenScore.dephaseFactor(gamma, dv, vRef)
        %%%
        %%% A damped oscillation exp(-gamma_t t) with Gaussian dephasing exp(-sigma^2 t^2/2),
        %%% counted from the gate start GATE_PERIODS/vRef: its late-time height relative to
        %%% the same pole without spread is
        %%%     F = rho sqrt(pi/2) erfcx((rho + q)/sqrt(2)) exp(-q^2/2)
        %%% with rho = gamma/dv and q = 2 pi GATE_PERIODS dv/vRef (the grating wavelength
        %%% cancels). F = 1 for dv = 0; F -> 0 for a narrow pole on a steep branch.
            dv = abs(dv); F = ones(size(gamma));
            b = isfinite(dv) & dv > 0;
            rho = gamma(b) ./ dv(b); q = 2 * GreenScore.GATE_PERIODS * pi * dv(b) ./ vRef;
            F(b) = rho .* sqrt(pi / 2) .* erfcx((rho + q) / sqrt(2)) .* exp(-q.^2 / 2);
        end
    end

    methods (Static, Access = private)
        function best = strongest(structure, angle, vRef)
        %%% Largest logStrength of one reference direction (-Inf without candidates)
            sr = CalcCommon.rotateStructure(structure, angle);
            L = CalcCommon.limitingAndBulk(sr.mater(sr.nLay));
            vS = L.vS; vTf = L.vTf;
            v = GreenScore.bandAxis(vS, vTf, GreenScore.DEFAULT_STEP);
            g = GreenSolver.calc(sr, 0, v.', 1, 3, 'Eta', SawSpectrum.ETA);
            pk = GreenScore.direction(structure, angle, g.G_eff(:), v, vS, vTf, vRef);
            best = max([pk.logStrength, -Inf]);
        end

        function pk = emptyPeaks()
            pk = struct('vel', {}, 'polVel', {}, 'gamma', {}, 'poleOk', {}, 'height', {}, ...
                'prominence', {}, 'base', {}, 'poleHeight', {}, 'relWidth', {}, 'relHeight', {}, ...
                'fBand', {}, 'rawScore', {}, 'slope', {}, 'slopeFwd', {}, 'slopeBwd', {}, ...
                'gateTerm', {}, 'dephaseTerm', {}, 'score', {}, 'logStrength', {}, 'lowScore', {});
        end

        function sl = probeOneSide(structureBase, angle, pk, offset)
        %%% Slope of every candidate from one probe solve at angle + offset [km/s per deg]
        %%% Stencil pole +- PROBE_HALF x PROBE_STEP; the probe pole is the interior local
        %%% maximum of |G| nearest the stencil centre, refined by the 3-point pole fit (or
        %%% a log-parabola vertex). NaN where the stencil has no interior local maximum.
            n = numel(pk); sl = nan(n, 1);
            vp = [pk.polVel]'; okp = [pk.poleOk]' & isfinite(vp); vp(~okp) = [pk(~okp).vel];
            k  = (-GreenScore.PROBE_HALF:GreenScore.PROBE_HALF)';
            V  = vp' + k * GreenScore.PROBE_STEP;                        % [nk x nPk]
            [vu, ~, iu] = unique(round(V(:) * 1e9) / 1e9);
            sr = CalcCommon.rotateStructure(structureBase, angle + offset);
            g  = GreenSolver.calc(sr, 0, vu.', 1, 3, 'Eta', SawSpectrum.ETA);
            Gs = reshape(g.G_eff(iu), size(V));
            for j = 1:n
                s = abs(Gs(:, j)); m = numel(s);
                lm = find(s(2:m-1) > s(1:m-2) & s(2:m-1) >= s(3:m)) + 1;   % interior local maxima
                if isempty(lm); continue; end
                [~, q] = min(abs(lm - (GreenScore.PROBE_HALF + 1))); i = lm(q);
                [pV, pG] = GreenScore.poleFit3(V(i-1:i+1, j), Gs(i-1:i+1, j));
                if isfinite(pV) && abs(pV - V(i, j)) <= GreenScore.PROBE_STEP && pG > 0
                    v2 = pV;
                else                                                        % log-parabola vertex
                    y = log(s(i-1:i+1)); den = y(1) - 2 * y(2) + y(3); v2 = V(i, j);
                    if den < 0
                        dl = 0.5 * (y(1) - y(3)) / den;
                        if abs(dl) <= 1; v2 = V(i, j) + dl * GreenScore.PROBE_STEP; end
                    end
                end
                sl(j) = (v2 - vp(j)) / offset;
            end
        end
    end

    methods (Static, Access = {?GreenScore, ?GreenScoreLayered})
        %%% Shared with the layered score (GreenScoreLayered), which finds and fits its
        %%% candidates with the same peak search, 3-point pole model and rules.
        function [found, v, G, s, h, loc, wid, prom] = bandPeaks(G, velAxis, lo, hi, errId)
        %%% The band [lo, hi] of a uniform axis and the local maxima of |G| on it
        %%% >> [found, v, G, s, h, loc, wid, prom] = GreenScore.bandPeaks(G, velAxis, lo, hi, errId)
        %%%
        %%% OUTPUT: found .. false if the band holds fewer than 5 (finite) samples (the
        %%%                  other outputs are then not to be used)
        %%%         v, G, s .. the band (+-1.5 steps beyond lo / hi) of velAxis, G and |G|
        %%%         h, loc, wid, prom .. findpeaks of the finite |G| (half-prominence
        %%%                  widths, prominence floor 1e-4 x the largest value), unfiltered
        %%%
        %%% errId: the caller's error for an axis not uniform over the band.
            v = velAxis(:); G = G(:); s = []; h = []; loc = []; wid = []; prom = [];
            found = false;
            if numel(v) < 5; return; end
            dv0 = v(2) - v(1);
            b = v >= lo - 1.5 * dv0 & v <= hi + 1.5 * dv0;
            v = v(b); G = G(b); s = abs(G);
            if numel(v) < 5; return; end
            dv = v(2) - v(1);
            if max(abs(diff(v) - dv)) > 1e-9
                error(errId, 'velAxis must be uniform over the band.');
            end
            m = find(isfinite(s)); if numel(m) < 5; return; end
            [h, loc, wid, prom] = findpeaks(s(m), v(m), 'WidthReference', 'halfprom', ...
                'MinPeakProminence', 1e-4 * max(s(m)));
            found = true;
        end

        function f = fitPeaks(v, G, h, loc, wid)
        %%% Pole fit of every grid peak: reported velocity, pole, half-width, height
        %%% >> f = GreenScore.fitPeaks(v, G, h, loc, wid)
        %%%
        %%% INPUT:  1-2. v, G .. uniform axis [km/s] and complex G13_eff on it (column)
        %%%         3-5. h, loc, wid .. findpeaks' heights, locations (on v) and
        %%%                             half-prominence widths of |G|
        %%% OUTPUT: f .. struct of columns: vel (log-parabola refined |G| peak), polVel,
        %%%              gamma, poleOk, poleHeight
        %%%
        %%% 3-point pole fit through the peak sample and the samples FIT_SPAN either side.
        %%% poleOk: a damped pole (gamma > 0) within the stencil. gamma: the fitted
        %%% half-width, half the half-prominence width when the fit failed, floored at the
        %%% damping width ETA vel/2 (no physical pole is narrower; a weak pole on a strong
        %%% neighbour's skirt can fit narrower). poleHeight: the model's maximum, never
        %%% below the grid peak, and the grid peak itself when the fit failed (a failed
        %%% fit's model maximum can exceed the peak threefold).
            n = numel(loc); dv = v(2) - v(1); s = abs(G);
            hs = max(1, round(GreenScore.FIT_SPAN / dv));
            f = struct('vel', loc(:), 'polVel', nan(n, 1), 'gamma', nan(n, 1), ...
                'poleOk', false(n, 1), 'poleHeight', h(:));
            for j = 1:n
                i0 = find(abs(v - loc(j)) < dv / 4, 1);
                if i0 > 1 && i0 < numel(v) && all(s(i0-1:i0+1) > 0)
                    y = log(s(i0-1:i0+1)); den = y(1) - 2 * y(2) + y(3);
                    if den < 0
                        dl = 0.5 * (y(1) - y(3)) / den;
                        if abs(dl) <= 1; f.vel(j) = loc(j) + dl * dv; end
                    end
                end
                pV = NaN; pG = NaN; hP = NaN;
                if i0 > hs && i0 + hs <= numel(v)
                    ix = [i0 - hs, i0, i0 + hs];
                    [pV, pG, hP] = GreenScore.poleFit3(v(ix), G(ix));
                end
                if ~isfinite(hP) || hP < h(j) || abs(pV - loc(j)) > hs * dv; hP = h(j); end
                poleOk = isfinite(pV) && abs(pV - loc(j)) <= hs * dv && pG > 0;
                if ~poleOk; hP = h(j); end                                  % failed fit: grid height
                gam = pG;
                if ~poleOk || ~isfinite(gam) || gam <= 0; gam = wid(j) / 2; end
                f.polVel(j) = pV; f.poleOk(j) = poleOk; f.poleHeight(j) = hP;
                f.gamma(j) = max(gam, SawSpectrum.ETA * f.vel(j) / 2);         % damping-width floor
            end
        end

        function [pV, pG, hP, r] = poleFit3(v, g)
        %%% 3-point pole model G(x) ~ (a0 + a1 x)/(x - q), x = v - v(2), through 3 complex samples
        %%% pV = Re pole [km/s], pG = -Im pole (half-width, > 0 for a damped resonance), hP =
        %%% max |model| over the 3-point span, r = the residue a0 + a1 q of the pole (the
        %%% model is a1 + r/(x - q)). NaN if ill-conditioned.
            pV = NaN; pG = NaN; hP = NaN; r = NaN;
            v = v(:); g = g(:);
            if numel(v) ~= 3 || any(~isfinite(g)); return; end
            x = v - v(2); sc = max(abs(x)); xs = x / sc;
            M = [ones(3, 1), xs, g]; rhs = g .* xs;
            if rcond(M' * M) < 1e-14; return; end
            c = M \ rhs; q = c(end);
            p = v(2) + q * sc;
            pV = real(p); pG = -imag(p);
            xf = linspace(-1, 1, 81)'; hP = max(abs((c(1) + c(2) * xf) ./ (xf - q)));
            r = (c(1) + c(2) * q) * sc;
        end
    end
end
