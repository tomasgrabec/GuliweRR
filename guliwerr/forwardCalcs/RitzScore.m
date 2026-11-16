classdef RitzScore
%%% RitzScore — the Ritz-only existence score: a gated Ritz modal sum that mimics
%%% G13, deciding which surface waves exist when no Green's function is computed
%%% (calcMode = 'ritzOnly').
%%%
%%% EXPERIMENTAL: checked against the measured maps published in Grabec et al. 2024 and the
%%% UTGS paper (arXiv:2510.10696); the scores, their scales and the interface may change.
%%%
%%% For one propagation direction, with the Ritz modes (v_n) and their surface
%%% displacements u_n (standard frame: x1 = propagation, x3 = surface normal):
%%%     S(v) = sum_n u3_n conj(u1_n) / (v_n^2 - v^2 (1 + i ETA))
%%% i.e. the modal expansion of the surface response to an in-plane force along x1,
%%% read as out-of-plane displacement (the TGS component). S is evaluated on a uniform
%%% velocity axis and scored with the shared late-time machinery (SawSpectrum): the
%%% late-time gate (SawSpectrum.lateTime, 3 periods of vRef), the candidate search and
%%% width filter (SawSpectrum.candidates), and
%%%     score = log10(prominence / floor)
%%% with the floor of THIS direction (ownReference(): the median late-time modal-sum
%%% level over a band -- half-space [0.8 x the slowest, the fastest shear limit], layered
%%% the whole band [0.3 vSmin, 1.05 vLmax]). The own-direction floor costs nothing, so a
%%% single angle needs no extra Ritz solves and classifies exactly as it does inside a
%%% sweep.
%%%
%%% Why it works: for M-orthonormal
%%% eigenvectors the driven surface response of the finite Ritz domain IS this modal
%%% sum; the discrete bulk modes of the finite domain dephase in time like the bulk
%%% continuum of a true half-space, so the late-time gate suppresses them as it
%%% suppresses bulk features of Green's function, while a surface-wave mode keeps
%%% ringing. The G13-like numerator u3 conj(u1) matters more than any other choice
%%% (an |u3|^2 weighting is clearly worse).
%%%
%%% HALF-SPACE RULE (direction(), acceptedHS()). Everything comes from the one
%%% Ritz solve of the angle -- no reference directions, no neighbouring angles:
%%%   1. score   .. the modal-sum score above (own floor), per candidate peak
%%%   2. mode    .. the LOCAL supporting mode: the strongest modal-sum term among the
%%%                 modes inside the peak's half-width (localMode). The global
%%%                 argmax over all modes can name the Rayleigh mode for a pseudo-SAW peak 0.4-0.8 km/s away
%%%   3. dScore  .. score + log10 F, F = GreenScore.dephaseFactor: the late-time height
%%%                 lost when the TGS aperture (GreenScore.APERTURE_SIGMA) averages the
%%%                 pole over neighbouring directions. Needs the supporting mode's angular
%%%                 slope dv/dtheta, exact from the same solve (RitzSolver.angleSlope,
%%%                 Hellmann-Feynman); pole half-width ETA v/2. A branch that moves fast
%%%                 with angle does not ring in the experiment however sharp it is here
%%%   4. veto    .. the supporting mode's freqDisp >= MIN_FREQDISP (non-dispersion)
%%%   5. windowMass (candidates above min shear only, see belowShear below) .. log10 of the summed
%%%                 near-surface weight r E1 fd of the modes within max(1.2 %, 30 m/s)
%%%                 of the peak (windowWeights). A slab artefact rings from one or two
%%%                 bulk-like modes; a real leaky wave is carried by several modes with
%%%                 surface-concentrated energy. Depth-normalised: + log10(relDepth / 30)
%%%                 (HS_WINDOW_DEPTH_REF)
%%%   6. bottom  (candidates above min shear only) .. the supporting mode must not be
%%%                 bottom-trapped (HS_BOTTOM10_MAX, HS_BOTTOM_RATIO_MAX): a U-shaped slab
%%%                 mode of a strongly leaky pseudo-SAW
%%%   'high': veto & dScore >= HS_THRESHOLD_HIGH & (below min shear | windowMass >=
%%%   HS_MIN_WINDOW_MASS & not bottom-trapped). 'low': the SAME gates with dScore >= HS_THRESHOLD_LOW, so
%%%   'low' is always a subset of 'high'.
%%%   Velocity: the supporting mode's eigenvelocity. Label: 'Rayleigh' if it is below
%%%   min shear + slabResolution (the Ritz slab's own discretisation error of the
%%%   slow-shear limit), else 'pseudoRayleigh'.
%%% Validation and limits: docs/reference/architecture/forwardCalcs_modeClassification.md, sections
%%% 11-12. A very leaky pseudo-SAW that shows in TGS can be dropped at 'high' (Fe3Al(316),
%%% 148-170 deg); a weak Rayleigh wave hugging min shear on a steep branch can be absent
%%% at 'low' (Fe3Al(316), 86-116 deg).
%%%
%%% LAYERED RULE (directionLayered(), acceptedLayered()). Built like the half-space
%%% rule and like the layered Green path (GreenScoreLayered): no layer is singled out as
%%% "the substrate".
%%%   band   .. the modal sum's axis, candidate search and floor all cover
%%%             [0.3 vSmin, 1.05 vLmax] (GreenScoreLayered.band: slowest shear / fastest
%%%             longitudinal limit over ALL layers); modes kept up to MODE_VEL_FACTOR x vLmax
%%%   score  .. dScore (own floor, local supporting mode, aperture dephasing, as the
%%%             half-space) + LAY_K_POL (1 - 2 SH), SH = |u2|^2 / |u|^2 of the supporting
%%%             mode at the free surface (0 sagittal .. 1 Love) -- GreenScoreLayered's term
%%%   gate   .. the supporting mode must be TRAPPED: bottom25 (its energy fraction in the
%%%             deepest 25 % of the whole domain, density- and thickness-weighted) <=
%%%             LAY_EBOT_MAX. A guided wave decays into the substrate (bottom25 ~1e-13);
%%%             a slab-continuum mode, and every leaky wave, reaches the rigid bottom
%%%             (>= 0.03). No freqDisp veto: guided waves in a stack are dispersive.
%%%   'low' (default): score >= LAY_THRESHOLD_LOW -- sagittal-type waves; 'high': score
%%%   >= LAY_THRESHOLD_HIGH -- also Love-type; the same gate, so 'low' is a subset of 'high'.
%%%   Velocity: the supporting mode's eigenvelocity; 'SAW' with aboveSubstrateST (label
%%%   only) against the substrate's slowest shear limit.
%%% Tuned on the measured NiTi/MgO TGS map (Grabec et al. 2024), as GreenScoreLayered:
%%% 'low' finds the sagittal-type waves, 'high' also most Love-type waves at the cost of
%%% false positives. Against GreenScoreLayered on four strongly leaky stacks (computed
%%% with GuliweRR): identical on every trapped wave. KNOWN LIMIT: leaky layered waves
%%% (above the substrate shear limit and radiating into it) are never reported -- the gate
%%% removes them with the slab continuum, which no other Ritz measure tested separates
%%% from them; use calcMode 'full' for those.
%%%
%%% Cost beyond the Ritz solve, half-space: Christoffel limiting velocities (~4 ms), the
%%% displacement profiles of the kept modes (~1 ms), sum, gate and candidates (~3 ms),
%%% the slopes of the supporting modes (~4 ms), slab resolution (<1 ms). No other
%%% direction is solved.
%%% See docs/reference/architecture/forwardCalcs_modeClassification.md (the Ritz method).
%%%
%%% Value class with Constant properties and Static methods only.
%%%
%%% GuliweRR v6.0

    properties (Constant)
        VEL_STEP        = 0.002  % km/s -- Ritz score axis: 2 m/s
                                 % (the analytic modal sum costs nothing per point; no padding
                                 % or tail taper -- the thresholds were calibrated this way)
        ETA             = 1e-2   % modal damping, denominator v_n^2 - v^2 (1 + i ETA). Larger
                                 % than SawSpectrum.ETA: a discrete Ritz mode is an exact pole,
                                 % and 1e-2 keeps it resolved on the 2 m/s axis and decayed
                                 % within the FFT window
        MODE_VEL_FACTOR = 1.2    % modes summed: v_n <= MODE_VEL_FACTOR x longLimVel (half-space)
                                 % or x the fastest longitudinal limit of any layer (layered)
        %%% --- Half-space rule (direction, acceptedHS) -- tuned together against measured
        %%% TGS maps: change none of them alone
        HS_THRESHOLD_HIGH  = 0.16    % 'high': dScore threshold
        HS_THRESHOLD_LOW   = 0.70    % 'low': dScore threshold, same other gates as 'high'
        HS_MIN_WINDOW_MASS = -1.55   % log10 window weight a candidate above min shear needs
                                     % (depth-normalised, see HS_WINDOW_DEPTH_REF)
        HS_WINDOW_DEPTH_REF = 30     % relDepth at which windowMass is unshifted: direction()
                                     % adds log10(relDepth / 30), as the summed weight of a
                                     % real leaky wave falls with the slab depth (more modes
                                     % share it). Without it a 60-lambda slab drops a
                                     % real leaky branch (computed with GuliweRR)
        HS_BOTTOM10_MAX     = 0.15   % bottom-trapped gate (pseudo-SAWs only): reject a candidate
        HS_BOTTOM_RATIO_MAX = 0.8    % whose supporting mode has > 15 % of its |u|^2 in the deepest
                                     % 10 % of the slab AND bottom10 / bottom25 > 0.8 (uniform: 0.4).
                                     % A strongly leaky pseudo-SAW the rigid-bottom slab cannot
                                     % radiate splits into "U-shaped" modes, surface AND bottom; the
                                     % lower member wins the peak 55-65 m/s below the real branch
                                     % (computed with GuliweRR). Evidence from ONE strongly
                                     % anisotropic cut; values for relDepth ~30 -- bottom10
                                     % depends on relDepth
        HS_MIN_LEG_PER_DEPTH = 2.3   % basis guard: legOrder >= 2.3 relDepth. Below it the slab
                                     % does not resolve the band up to vL and the rule breaks
                                     % (legOrder < 2 relDepth: hundreds of slab modes near vL
                                     % accepted); checkBasisHS warns
        HS_FLOOR_LO        = 0.8     % floor band start (ownReference), x min shear; the band ends
                                     % at the fastest shear limit
        WINDOW_REL         = 0.012   % window half-width max(WINDOW_REL v, WINDOW_ABS): the
        WINDOW_ABS         = 0.030   % velocity tolerance of a match [-, km/s]
        WINDOW_TOP         = 0.98    % modes above WINDOW_TOP x vL carry no window weight
                                     % (slab modes crowd the longitudinal limit)
        LOCAL_MIN_HALF     = 0.004   % km/s -- smallest half-width localMode searches (2 axis steps)

        %%% --- Layered rule (directionLayered, acceptedLayered) -- tuned together on
        %%% the measured NiTi/MgO map (Grabec et al. 2024): change none alone
        LAY_THRESHOLD_LOW  = 0.65    % 'low' (default): sagittal-type waves. Centre of the window
                                     % 0.5-0.8 that admits no Love-type wave. Lower thresholds
                                     % recover the Rayleigh-cap ends but admit the lower Love
                                     % arm first: the modal sum's pole half-width ETA v/2
                                     % dephases a steep arm far less than Green's does
        LAY_THRESHOLD_HIGH = -0.7    % 'high': also Love-type waves
        LAY_K_POL          = 0.75    % polarisation weight: score = dScore + LAY_K_POL (1 - 2 SH)
                                     % (as GreenScoreLayered.K_POL; 0-1 give the same 'low' result,
                                     % 0.75 orders Love waves before spurious peaks best)
        LAY_EBOT_MAX       = 1e-2    % trapped gate: the supporting mode's bottom25 <= this.
                                     % Guided waves 1e-22-1e-12, slab continuum and leaky waves
                                     % >= 0.03 (four stacks, substrate depths 20 and 30 lambda);
                                     % 3e-3 to 1e-1 classify NiTi/MgO identically
        MIN_FREQDISP    = 0.3    % HALF-SPACE ONLY non-dispersion veto: an accepted
                                 % modal-sum peak whose supporting Ritz mode (localMode) has
                                 % freqDisp below this is rejected (0.3 also best for the
                                 % half-space rule: 0.5 drops real branches on a strongly
                                 % anisotropic cut).
                                 % Physics: a homogeneous half-space has no length scale, so
                                 % a true SAW is non-dispersive (v_g = v_p); a mode of the
                                 % finite-depth Ritz slab is set by the domain depth and is
                                 % dispersive. The cut allows |v_p - v_g|/v_p up to ~35 %.
                                 % Removes bulk-like slab modes that ring coherently near vL
                                 % (e.g. Ni(110) ~5 km/s, v/vL 0.92-0.98).
                                 % NOT applied to layered stacks: guided waves in a layered
                                 % stack are genuinely dispersive -- the film thickness sets a
                                 % length scale, so v_g ~= v_p for true modes -- and freqDisp
                                 % cannot separate them from finite-domain artefacts (the
                                 % layered rule uses the trapped gate LAY_EBOT_MAX instead).
                                 % Cuts 0.2-0.3 similar, 0.3 best.
    end

    methods (Static)
        function thr = resolveThresholdHS(t)
        %%% Half-space dScore threshold from 'high' | 'low' | a number (see acceptedHS)
        %%% >> thr = RitzScore.resolveThresholdHS(t)
            thr = CalcCommon.thresholdValue(t, RitzScore.HS_THRESHOLD_HIGH, RitzScore.HS_THRESHOLD_LOW, ...
                'RitzScore:resolveThresholdHS:badName', 'RitzScore.resolveThresholdHS');
        end

        function sc = candidateScore(pk, isLayered)
        %%% The score the Ritz rule thresholds, per candidate
        %%% >> sc = RitzScore.candidateScore(pk, isLayered)
        %%%
        %%% INPUT:  1. pk        .. candidates from direction() / directionLayered()
        %%%         2. isLayered .. true for a layered stack
        %%% OUTPUT: sc .. [numel(pk) x 1]: half-space dScore; layered dScore + polTerm
        %%%               (the polarisation term). acceptedHS / acceptedLayered compare
        %%%               it with the threshold; Vis.ritzModeScores shows it.
            if isempty(pk); sc = zeros(0, 1); return; end
            sc = [pk.dScore]';
            if isLayered
                sc = sc + [pk.polTerm]';
            end
        end

        function acc = acceptedHS(pk, threshold)
        %%% Half-space acceptance of scored candidates (direction()) at a threshold
        %%% >> acc = RitzScore.acceptedHS(pk, threshold)
        %%%
        %%% INPUT:  1. pk        .. struct array from direction()
        %%%         2. threshold .. 'high' | 'low' | numeric (a dScore threshold with the
        %%%                         'high' gates)
        %%% OUTPUT: acc .. [numel(pk) x 1] logical:
        %%%                ~vetoed & dScore >= thr & (belowShear | windowMass >= HS_MIN_WINDOW_MASS
        %%%                & ~bottomTrapped)
        %%%                bottomTrapped = bottom10 > HS_BOTTOM10_MAX & bottom10 / bottom25 >
        %%%                HS_BOTTOM_RATIO_MAX (a U-shaped supporting mode, see the constants;
        %%%                never true when the fields are NaN)
        %%% The gates are the same at every threshold, so a higher threshold only ever
        %%% removes entries ('low' is a subset of 'high').
            thr = RitzScore.resolveThresholdHS(threshold);
            if isempty(pk); acc = false(0, 1); return; end
            b10 = [pk.bottom10]'; b25 = [pk.bottom25]';
            trapped = b10 > RitzScore.HS_BOTTOM10_MAX & b10 ./ b25 > RitzScore.HS_BOTTOM_RATIO_MAX;
            acc = ~[pk.vetoed]' & RitzScore.candidateScore(pk, false) >= thr & ([pk.belowShear]' | ...
                ([pk.windowMass]' >= RitzScore.HS_MIN_WINDOW_MASS & ~trapped));
        end

        function pk = direction(ritz, num, E1, r, velAxis, ref, vS, vL, dKdC_cell, Ebot)
        %%% Scored half-space candidates of one direction (threshold-independent)
        %%% >> pk = RitzScore.direction(ritz, num, E1, r, velAxis, ref, vS, vL, dKdC_cell, Ebot)
        %%%
        %%% INPUT:  1. ritz       .. RitzSolverSAW, calc()'d, eigenvectors present
        %%%         2. num, E1, r .. [nKeep x 1] from modeProfiles
        %%%         3. velAxis    .. score axis (SawSpectrum.velocityAxis(1.05 vL, VEL_STEP))
        %%%         4. ref        .. ownReference struct (vRef, floor) plus slabRes
        %%%                          (slabResolution, km/s): the Rayleigh side's tolerance
        %%%         5. vS, vL     .. min shear and longitudinal limiting velocity [km/s]
        %%%         6. dKdC_cell  .. the precompute's GammaMatrix.derivsByDepth cell
        %%%                          (built if empty -- on every call, ~0.1 s at CalcParams(30,70);
        %%%                          the Calc classes pass it, ad-hoc callers should too)
        %%%         7. Ebot       .. [nKeep x 2] from modeProfiles (the bottom-trapped gate)
        %%% Both ref.slabRes and Ebot are required (RitzScore:direction:missingInput): each
        %%% changes the result when left out.
        %%% OUTPUT: pk .. struct array in velocity order, one per SawSpectrum.candidates
        %%%               peak, fields as ModeClassification.ritzScorePeaks: vel (modal-sum
        %%%               peak), modeNo / modeVel (local supporting mode, localMode), score,
        %%%               dScore (= score + dephaseTerm), prominence, relWidth, relHeight,
        %%%               share, freqDisp, slope [km/s per deg], windowMass (depth-
        %%%               normalised: + log10(relDepth / HS_WINDOW_DEPTH_REF)), belowShear
        %%%               (supporting mode below vS + ref.slabRes: the Rayleigh side; the
        %%%               peak itself can jitter across vS), vetoed (freqDisp <
        %%%               MIN_FREQDISP), bottom10 / bottom25 (the supporting mode's |u|^2
        %%%               fraction in the deepest 10 % / 25 % of the slab; diagnostics),
        %%%               accepted (false here; ModeClassification sets it)
        %%%
        %%% Everything a reclassification needs, so no eigenvector is read after this.
        %%% Only the supporting modes' slopes are evaluated (one dK/dtheta, one quadratic
        %%% form per distinct mode).
            pk = RitzScore.emptyPeaks();
            nK = numel(num);
            if ~isfield(ref, 'slabRes') || nargin < 10 || ~isequal(size(Ebot), [nK, 2])
                error('RitzScore:direction:missingInput', ['RitzScore.direction needs ref.slabRes ' ...
                    '(slabResolution) and Ebot [nKeep x 2] (modeProfiles).']);
            end
            vn = ritz.vels(1:nK);      vn = vn(:);
            vg = ritz.velsGroup(1:nK); vg = vg(:);
            %%% freqDisp as RitzSolverSAW.modeDiagnostics defines it, computed here so the
            %%% rule does not depend on the diagnostics having run
            fd = exp(-10 * ((vn - vg) ./ max(vn, eps)).^2);
            c  = SawSpectrum.candidates(RitzScore.spectrum(vn, num, velAxis), velAxis, vL, ref, 1, 0);
            nP = numel(c);
            if nP == 0; return; end
            w  = RitzScore.windowWeights(vn, fd, E1, r, vS, vL);
            %%% Which side of min shear a candidate is on is decided by its supporting
            %%% mode, with the slab's resolution as tolerance -- exactly the Rayleigh label
            %%% (ModeClassification). The modal-sum PEAK of a Rayleigh wave hugging vS
            %%% jitters by +-1-3 m/s on the 2 m/s axis and can land above vS: judged by
            %%% the peak it would face the pseudo-SAW window-mass gate, which a deep
            %%% vS-hugging mode fails (synthetic cubic, Zener ratio 20, (112) cut, 120-121 deg:
            %%% dropped).
            vB = vS + ref.slabRes;
            wmShift = log10(ritz.structure.calcParams(1).relDepth / RitzScore.HS_WINDOW_DEPTH_REF);
            md = zeros(nP, 1); sh = zeros(nP, 1); wm = zeros(nP, 1); below = false(nP, 1);
            for j = 1:nP
                [md(j), sh(j)] = RitzScore.localMode(vn, num, c(j).vel, c(j).relWidth);
                below(j) = vn(md(j)) < vB;
                wm(j) = RitzScore.windowMass(vn, w, c(j).vel, vS, ~below(j)) + wmShift;
            end
            um = unique(md);
            slU = ritz.angleSlope(um, dKdC_cell);
            [~, ii] = ismember(md, um);
            sl = slU(ii);
            dT = RitzScore.dephaseTerm(vn(md), sl, ref.vRef);
            for j = 1:nP
                pk(j) = struct('vel', c(j).vel, 'modeNo', md(j), 'modeVel', vn(md(j)), ...
                    'score', c(j).score, 'dScore', c(j).score + dT(j), ...
                    'prominence', c(j).prominence, 'relWidth', c(j).relWidth, ...
                    'relHeight', c(j).relHeight, 'share', sh(j), 'freqDisp', fd(md(j)), ...
                    'slope', sl(j), 'windowMass', wm(j), 'belowShear', below(j), ...
                    'vetoed', fd(md(j)) < RitzScore.MIN_FREQDISP, ...
                    'bottom10', Ebot(md(j), 1), 'bottom25', Ebot(md(j), 2), ...
                    'SH', NaN, 'polTerm', NaN, 'accepted', false);
            end
        end

        function thr = resolveThresholdLayered(t)
        %%% Layered score threshold from 'low' | 'high' | a number (see acceptedLayered)
        %%% >> thr = RitzScore.resolveThresholdLayered(t)
            thr = CalcCommon.thresholdValue(t, RitzScore.LAY_THRESHOLD_HIGH, RitzScore.LAY_THRESHOLD_LOW, ...
                'RitzScore:resolveThresholdLayered:badName', 'RitzScore.resolveThresholdLayered');
        end

        function acc = acceptedLayered(pk, threshold)
        %%% Layered acceptance of scored candidates (directionLayered) at a threshold
        %%% >> acc = RitzScore.acceptedLayered(pk, threshold)
        %%%
        %%% OUTPUT: acc .. [numel(pk) x 1] logical: bottom25 <= LAY_EBOT_MAX (the supporting
        %%%                mode is trapped) & dScore + polTerm >= threshold. The gate is the
        %%%                same at every threshold, so 'low' is a subset of 'high'.
            thr = RitzScore.resolveThresholdLayered(threshold);
            if isempty(pk); acc = false(0, 1); return; end
            acc = [pk.bottom25]' <= RitzScore.LAY_EBOT_MAX & ...
                RitzScore.candidateScore(pk, true) >= thr;
        end

        function pk = directionLayered(ritz, num, Ebot, SH, velAxis, ref, vLo, vHi, vSsub, dKdC_cell)
        %%% Scored layered candidates of one direction (threshold-independent)
        %%% >> pk = RitzScore.directionLayered(ritz, num, Ebot, SH, velAxis, ref, vLo, vHi, vSsub, dKdC_cell)
        %%%
        %%% INPUT:  1. ritz       .. RitzSolverSAW of the layered stack, calc()'d
        %%%         2-4. num, Ebot, SH .. modeProfiles outputs of the kept modes
        %%%         5. velAxis    .. layeredAxis(vLo, vHi)
        %%%         6. ref        .. ownReference over [vLo, vHi] (vRef, floor)
        %%%         7-8. vLo, vHi .. band (GreenScoreLayered.band of this direction)
        %%%         9. vSsub      .. the substrate's slowest shear limit [km/s] -- a label
        %%%                          only (belowShear, aboveSubstrateST), never a score input
        %%%        10. dKdC_cell  .. the precompute's GammaMatrix.derivsByDepth cell (built
        %%%                          if empty -- on every call; the Calc classes pass it)
        %%% OUTPUT: pk .. struct array in velocity order (emptyPeaks fields): vel (modal-sum
        %%%               peak), modeNo / modeVel (local supporting mode), score, dScore
        %%%               (= score + aperture dephasing at the gate's vRef), prominence,
        %%%               relWidth, relHeight, share, slope, SH (supporting mode's surface
        %%%               |u2|^2 / |u|^2), polTerm (LAY_K_POL (1 - 2 SH)), bottom10 /
        %%%               bottom25 (its energy fraction in the deepest 10 % / 25 % of the
        %%%               domain), belowShear (modeVel < vSsub), freqDisp / windowMass NaN,
        %%%               vetoed false, accepted false (ModeClassification sets it)
        %%% Candidates are SawSpectrum.candidates' peaks of the late-time modal sum below vHi,
        %%% kept from vLo up.
            pk = RitzScore.emptyPeaks();
            nK = numel(num);
            vn = ritz.vels(1:nK); vn = vn(:);
            c  = SawSpectrum.candidates(RitzScore.spectrum(vn, num, velAxis), velAxis, vHi, ref, 1, 0);
            if isempty(c); return; end
            c  = c([c.vel] >= vLo);
            nP = numel(c);
            if nP == 0; return; end
            md = zeros(nP, 1); sh = zeros(nP, 1);
            for j = 1:nP
                [md(j), sh(j)] = RitzScore.localMode(vn, num, c(j).vel, c(j).relWidth);
            end
            um = unique(md);
            slU = ritz.angleSlope(um, dKdC_cell);
            [~, ii] = ismember(md, um);
            sl = slU(ii);
            dT = RitzScore.dephaseTerm(vn(md), sl, ref.vRef);
            for j = 1:nP
                m = md(j);
                pk(j) = struct('vel', c(j).vel, 'modeNo', m, 'modeVel', vn(m), ...
                    'score', c(j).score, 'dScore', c(j).score + dT(j), ...
                    'prominence', c(j).prominence, 'relWidth', c(j).relWidth, ...
                    'relHeight', c(j).relHeight, 'share', sh(j), 'freqDisp', NaN, ...
                    'slope', sl(j), 'windowMass', NaN, 'belowShear', vn(m) < vSsub, ...
                    'vetoed', false, 'bottom10', Ebot(m, 1), 'bottom25', Ebot(m, 2), ...
                    'SH', SH(m), 'polTerm', RitzScore.LAY_K_POL * (1 - 2 * SH(m)), ...
                    'accepted', false);
            end
        end

        function v = layeredAxis(vLo, vHi)
        %%% Layered modal-sum axis: VEL_STEP lattice SawSpectrum.VEL_MIN + n x VEL_STEP, column
        %%% >> v = RitzScore.layeredAxis(vLo, vHi)
        %%%
        %%% From min(VEL_MIN, the lattice point at or below vLo) to vHi + 0.1 km/s +
        %%% SawSpectrum.EDGE_TAPER, so the search (below vHi) never reaches the tapered edge
        %%% and a very slow stack (vLo < 0.3 km/s) keeps its whole band. On the lattice, a
        %%% single angle and any sweep share the same points.
            v0 = SawSpectrum.VEL_MIN; dv = RitzScore.VEL_STEP;
            n1 = min(0, floor((vLo - v0) / dv + 1e-9));
            n2 = floor((vHi + 0.1 + SawSpectrum.EDGE_TAPER - v0) / dv + 1e-9);
            v  = SawSpectrum.latticeAxis(n1, n2, dv);
        end

        function pk = emptyPeaks()
        %%% The empty candidate struct array (ModeClassification.ritzScorePeaks field set)
        %%% >> pk = RitzScore.emptyPeaks()
        %%%
        %%% Both rules fill the same fields (half-space: SH, polTerm NaN; layered: freqDisp,
        %%% windowMass NaN), so the two paths concatenate. Also the default of
        %%% ModeClassification.ritzScorePeaks.
            pk = struct('vel', {}, 'modeNo', {}, 'modeVel', {}, 'score', {}, 'dScore', {}, ...
                'prominence', {}, 'relWidth', {}, 'relHeight', {}, 'share', {}, ...
                'freqDisp', {}, 'slope', {}, 'windowMass', {}, 'belowShear', {}, ...
                'vetoed', {}, 'bottom10', {}, 'bottom25', {}, 'SH', {}, 'polTerm', {}, ...
                'accepted', {});
        end

        function [num, E1, r, Ebot, SH] = modeProfiles(ritz, vLong)
        %%% Modal-sum weights and near-surface measures of the kept modes, one displacement call
        %%% >> [num, E1, r] = RitzScore.modeProfiles(ritz, vLong)
        %%% >> [num, E1, r, Ebot, SH] = RitzScore.modeProfiles(ritz, vLong)
        %%%
        %%% INPUT:  as modalNumerators (layered: vLong = the fastest longitudinal limit of
        %%%         any layer)
        %%% OUTPUT: num .. [nKeep x 1] u3 conj(u1) at the free surface (= modalNumerators)
        %%%         E1  .. [nKeep x 1] fraction of the mode's |u|^2 within one wavelength of
        %%%                the surface (depth-integrated on the profile grid)
        %%%         r   .. [nKeep x 1] |u3(surface)|^2 / depth-mean |u|^2 / relDepth: the
        %%%                surface-normal amplitude a TGS probe reads, relative to the mode's
        %%%                spread over the slab
        %%%         Ebot .. [nKeep x 2] fraction of the mode's |u|^2 in the deepest 10 % / 25 %
        %%%                of the slab (0.10 / 0.25 for a uniform profile). A leaky wave
        %%%                the rigid-bottom slab cannot radiate can form "U-shaped" modes,
        %%%                with energy at the surface AND the bottom (computed with
        %%%                GuliweRR: Ebot10 / Ebot25 >= 0.84 on such modes, <= 0.73 on
        %%%                true pseudo-SAWs)
        %%%         SH  .. [nKeep x 1] |u2|^2 / |u|^2 at the free surface (0 sagittal .. 1
        %%%                shear-horizontal)
        %%% Half-space: plain |u|^2 on layer 1's uniform grid. Layered (nLay > 1): every
        %%% layer's profile (calcDisplacement: the same point count per layer, layer 1 the
        %%% deepest) weighted by its density and its own trapezoid depth step -- the layers
        %%% differ in thickness, so unweighted row sums would over-weight a thin film --
        %%% and depths measured over the whole domain: E1 = energy fraction within one
        %%% wavelength of the surface, r = |u3(surface)|^2 / depth-integrated rho |u|^2,
        %%% Ebot = energy fraction in the deepest 10 % / 25 % of the domain (the layered
        %%% rule's trapped gate, LAY_EBOT_MAX).
            if isempty(ritz.eigVecs_con)
                error('RitzScore:modalNumerators:noEigVecs', ...
                    'Ritz eigenvectors are needed (compute before retentionLevel pruning).');
            end
            nKeep = find(ritz.vels(1:ritz.nCompModes) <= RitzScore.MODE_VEL_FACTOR * vLong, 1, 'last');
            if isempty(nKeep)
                num = complex(zeros(0, 1)); E1 = zeros(0, 1); r = zeros(0, 1); Ebot = zeros(0, 2);
                SH = zeros(0, 1); return;
            end
            D   = ritz.calcDisplacement(nKeep);      % [nPts x nKeep x 3], surface = last row
            us  = reshape(D(end, :, :), nKeep, 3);
            num = us(:, 3) .* conj(us(:, 1));
            if nargout > 4
                u2 = abs(us).^2;
                SH = u2(:, 2) ./ max(sum(u2, 2), realmin);
            end
            if ritz.structure.nLay > 1 && nargout > 1
                [E1, r, Ebot] = layeredProfiles(ritz.structure, D, us);
                return
            end
            if nargout > 1
                P   = sum(abs(D).^2, 3);             % [nPts x nKeep]
                relDepth = ritz.structure.calcParams(1).relDepth;
                z   = flipud(linspace(0, relDepth, size(P, 1))');   % depth below the surface
                E1  = sum(P(z <= 1, :), 1).' ./ max(sum(P, 1).', eps);
                r   = abs(us(:, 3)).^2 ./ max(mean(P, 1).', eps) / relDepth;
            end
            if nargout > 3
                Et   = max(sum(P, 1).', eps);
                Ebot = [sum(P(z >= 0.90 * relDepth, :), 1).' ./ Et, ...
                        sum(P(z >= 0.75 * relDepth, :), 1).' ./ Et];
            end
        end

        function [modeNo, share] = localMode(vn, num, vPeak, relWidth)
        %%% The supporting Ritz mode of a modal-sum peak: the strongest term near it
        %%% >> [modeNo, share] = RitzScore.localMode(vn, num, vPeak, relWidth)
        %%%
        %%% OUTPUT: modeNo .. argmax |num_n / (v_n^2 - vPeak^2 (1 + i ETA))| over the modes
        %%%                   with |v_n - vPeak| <= max(relWidth vPeak / 2, LOCAL_MIN_HALF);
        %%%                   the mode nearest vPeak if none is that close
        %%%         share  .. that term's fraction of sum_n |term_n| over ALL modes
        %%% Why local: a strong Rayleigh pole keeps a large term far from its own velocity,
        %%% so the global argmax over all modes can name it for a pseudo-SAW peak hundreds of
        %%% m/s higher, and the veto/slope/velocity would then belong to the wrong mode.
            vn = vn(:);
            t  = abs(num(:) ./ (vn.^2 - vPeak^2 * (1 + 1i * RitzScore.ETA)));
            in = find(abs(vn - vPeak) <= max(relWidth * vPeak / 2, RitzScore.LOCAL_MIN_HALF));
            if isempty(in)
                [~, modeNo] = min(abs(vn - vPeak));
            else
                [~, k] = max(t(in)); modeNo = in(k);
            end
            share = t(modeNo) / max(sum(t), eps);
        end

        function w = windowWeights(vn, fd, E1, r, vS, vL)
        %%% Per-mode near-surface weight used by windowMass
        %%% >> w = RitzScore.windowWeights(vn, fd, E1, r, vS, vL)
        %%%
        %%% w = r fd below vS; r E1 fd at or above vS (a leaky mode must ALSO keep its
        %%% energy near the surface -- a bound mode is surface-confined anyway); 0 above
        %%% WINDOW_TOP x vL.
            vn = vn(:);
            w  = r(:) .* fd(:);
            up = vn >= vS;
            w(up) = w(up) .* E1(up);
            w(vn > RitzScore.WINDOW_TOP * vL) = 0;
        end

        function wm = windowMass(vn, w, vPeak, vS, above)
        %%% log10 of the window weight around a peak, from the modes on one side of vS
        %%% >> wm = RitzScore.windowMass(vn, w, vPeak, vS, above)
        %%%
        %%% Sum of w over the modes within max(WINDOW_REL vPeak, WINDOW_ABS) of vPeak,
        %%% taking only modes at or above vS if above is true, only modes below vS
        %%% otherwise (a bound Rayleigh mode must not lend weight to a pseudo-SAW just
        %%% above vS, nor the reverse). direction() passes the candidate's side
        %%% (~belowShear). Floored at 1e-12 before the log.
            vn  = vn(:);
            sel = abs(vn - vPeak) <= max(RitzScore.WINDOW_REL * vPeak, RitzScore.WINDOW_ABS) ...
                & (vn >= vS) == above;
            wm  = log10(max(sum(w(sel)), 1e-12));
        end

        function dT = dephaseTerm(vMode, slope, vRef)
        %%% log10 aperture dephasing factor of supporting modes
        %%% >> dT = RitzScore.dephaseTerm(vMode, slope, vRef)
        %%%
        %%% log10 GreenScore.dephaseFactor(gamma, dv, vRef), with the modal sum's own pole
        %%% half-width gamma = ETA vMode / 2 and the velocity spread dv = |slope| x
        %%% GreenScore.APERTURE_SIGMA. <= 0; floored at -6.
            F  = GreenScore.dephaseFactor(RitzScore.ETA * vMode(:) / 2, ...
                abs(slope(:)) * GreenScore.APERTURE_SIGMA, vRef);
            dT = log10(max(F, 1e-6));
        end

        function dD = slabResolution(mat, vS, thetaStar, relDepth)
        %%% How far above the true min shear velocity the Ritz slab can place it [km/s]
        %%% >> dD = RitzScore.slabResolution(mat, vS, thetaStar, relDepth)
        %%%
        %%% INPUT:  1. mat       .. the half-space Material of the ROTATED structure
        %%%         2. vS        .. min(shearLimVels) [km/s]
        %%%         3. thetaStar .. its sagittal root angle (CalcCommon.limitingAndBulk
        %%%                         shearLimTheta of that entry)
        %%%         4. relDepth  .. Ritz domain depth [wavelengths]
        %%% OUTPUT: dD .. >= 0 [km/s]
        %%%
        %%% The slow-shear limit is a bulk wave grazing the surface, q = k_z/k -> q0. A
        %%% slab of depth D (rigid bottom) cannot hold a single grazing partial wave: its
        %%% lowest standing pattern pairs two partial waves whose k_z differ by pi/D
        %%% (quarter wave), i.e. a window of width w = 1/(2 D) in q. Its trace velocity is
        %%% at least the larger of the two partial waves' trace velocities
        %%%     v_tr(q) = v_phase(n(q)) sqrt(1 + q^2),  n = [1 0 q]/sqrt(1 + q^2),
        %%% so the slab's minimum is min over window positions [a, a + w] around q0 of
        %%% max(v_tr(a), v_tr(a + w)). dD = that minus vS. Ritz eigenvalues are upper
        %%% bounds, so the error is one-sided: a Rayleigh wave hugging vS can sit up to dD
        %%% above it in the slab. (Synthetic cubic, Zener ratio 20, (112) cut: modes <= 0.18 m/s
        %%% above vS, dD 0.21-0.30 m/s; other cuts computed: dD 0.2-0.6 m/s.)
        %%% thetaStar is turned into q0 by whichever of tan, -tan, cot reproduces vS
        %%% (LimitingVels' angle convention); 11 window positions, ~24 3x3 eigenproblems.
            C9  = reshape(permute(mat.GetStandardElastTensor(0), [1 3 2 4]), 9, 9);   % rows (i,k), cols (j,l)
            rho = mat.rho;
            vtr = @(q) traceVel(C9, rho, q);
            qs  = [tan(thetaStar); -tan(thetaStar); cot(thetaStar)];
            qs  = qs(isfinite(qs));
            [~, b] = min(arrayfun(@(q) abs(vtr(q) - vS), qs));
            q0  = qs(b);
            w   = 1 / (2 * relDepth);
            a   = q0 - w + (0:10)' * w / 10;
            dD  = min(max(arrayfun(vtr, a), arrayfun(vtr, a + w))) - vS;
            dD  = max(dD, 0);
        end

        function num = modalNumerators(ritz, vLong)
        %%% Per-mode modal-sum weights u3 conj(u1) at the free surface
        %%% >> num = RitzScore.modalNumerators(ritz, vLong)
        %%%
        %%% INPUT:  1. ritz  .. RitzSolverSAW, calc()'d, eigenvectors still present
        %%%                     (i.e. BEFORE any retentionLevel pruning)
        %%%         2. vLong .. longitudinal limiting velocity [km/s] of the half-space,
        %%%                     or the fastest of any layer (layered)
        %%% OUTPUT: num .. [nKeep x 1] complex, one entry per mode 1..nKeep, where nKeep is
        %%%                the last mode with vel <= MODE_VEL_FACTOR x vLong. Modes are
        %%%                velocity-sorted, so num(n) belongs to ritz.vels(n)
        %%%
        %%% The surface displacement is the top row of calcDisplacement (top point of the
        %%% top layer = the free surface). The numerators survive pruning of the
        %%% eigenvectors: they are all that reclassification needs later.
            num = RitzScore.modeProfiles(ritz, vLong);   % one output: no profile measures
        end

        function S = spectrum(vn, num, velAxis)
        %%% Complex Ritz modal sum on a velocity axis
        %%% >> S = RitzScore.spectrum(vn, num, velAxis)
        %%%
        %%% INPUT:  1. vn      .. [nKeep x 1] Ritz phase velocities [km/s]
        %%%         2. num     .. [nKeep x 1] RitzScore.modalNumerators
        %%%         3. velAxis .. [nVel x 1] km/s (SawSpectrum.velocityAxis)
        %%% OUTPUT: S .. [nVel x 1] complex
        %%%
        %%% A loop over modes, not an [nKeep x nVel] matrix: the same cost, without a
        %%% ~10-20 MB temporary per angle (which a thread pool would multiply).
            w = velAxis(:).^2 * (1 + 1i * RitzScore.ETA);
            S = complex(zeros(numel(w), 1));
            for n = 1:numel(num)
                S = S + num(n) ./ (vn(n)^2 - w);
            end
        end

        function share = decompose(vn, num, v)
        %%% Each Ritz mode's share of the modal sum S at one velocity
        %%% >> share = RitzScore.decompose(vn, num, v)
        %%%
        %%% INPUT:  1. vn  .. [nKeep x 1] Ritz phase velocities [km/s]
        %%%         2. num .. [nKeep x 1] RitzScore.modalNumerators
        %%%         3. v   .. velocity [km/s]
        %%% OUTPUT: share .. [nKeep x 1] |term_n| / sum_m |term_m|, term_n = num_n /
        %%%                  (v_n^2 - v^2 (1 + i ETA)) -- the contributing modes of S at v
        %%%                  (the measure localMode's share uses)
            t = abs(num(:) ./ (vn(:).^2 - v^2 * (1 + 1i * RitzScore.ETA)));
            share = t / max(sum(t), eps);
        end

        function ref = ownReference(structure, vn, num, velAxis, fLo, fHi)
        %%% Floor reference of ONE direction from its already computed modal sum
        %%% >> ref = RitzScore.ownReference(structure, vn, num, velAxis, 0.8 * vS, vTf)   % half-space
        %%% >> ref = RitzScore.ownReference(structure, vn, num, velAxis, vLo, vHi)        % layered
        %%%
        %%% INPUT:  1. structure .. Structure (any orientation: only the surface layer's
        %%%                         Hill shear modulus and density are read)
        %%%         2. vn, num   .. [nKeep x 1] Ritz velocities and modalNumerators
        %%%         3. velAxis   .. [nVel x 1] the direction's Ritz score axis
        %%%         4-5. fLo, fHi .. floor band [km/s]: half-space [0.8 x min shear, the
        %%%                         fastest shear limit vTf]; layered the whole band
        %%%                         [0.3 vSmin, 1.05 vLmax] (GreenScoreLayered.band)
        %%% OUTPUT: ref .. struct: vRef [km/s] = sqrt(G_Hill/rho) of the surface layer (the
        %%%                late-time gate's and the dephasing's time scale); floor = median
        %%%                late-time |S| over [fLo, fHi]
        %%%
        %%% The classification's floor on the Ritz-only path (the Calc classes call it per
        %%% angle).
        %%% Why the half-space band ends at vTf: between vTf and vL the slab's modes are not
        %%% converged in legOrder at the production bases, so a floor taken there would sit
        %%% 0.15-0.2 dec too high and every score, and the best thresholds, would move with
        %%% legOrder / relDepth. Over [0.8 vS, vTf] the floor scales with relDepth exactly as
        %%% the prominences do, so the score does not depend on the basis. The layered band is GreenScoreLayered's: no
        %%% layer is singled out. Why the surface layer's vRef on layers too (a near-surface
        %%% property): a gate at the slowest shear limit of any layer classifies NiTi/MgO's
        %%% sagittal-type waves the same but ranks Love waves far worse against spurious
        %%% peaks.
            ref.vRef = CalcCommon.surfaceVRef(structure);
            A = SawSpectrum.lateTime(RitzScore.spectrum(vn, num, velAxis), velAxis, ref.vRef, 1, 0);
            b = isfinite(A) & velAxis(:) >= fLo & velAxis(:) <= fHi;
            ref.floor = median(A(b));
        end

        function ok = checkBasisHS(calcParams)
        %%% Whether a half-space Ritz basis resolves the band the Ritz-only rule scores;
        %%% warns (once per session) if not
        %%% >> ok = RitzScore.checkBasisHS(calcParams)
        %%%
        %%% ok = legOrder >= HS_MIN_LEG_PER_DEPTH x relDepth. The rule's thresholds hold
        %%% (nearly unchanged results) for legOrder 2.3-3.5 x relDepth (relDepth 20-40) and
        %%% collapse below 2 x relDepth. The defaults (20, 50) and
        %%% the production bases (20, 50), (30, 70) pass. Warns once per session;
        %%% warning('on', id) re-enables it.
            ok = calcParams.legOrder >= RitzScore.HS_MIN_LEG_PER_DEPTH * calcParams.relDepth;
            id = 'RitzScore:unresolvedBasis';
            if ~ok && strcmp(warning('query', id).state, 'on')
                warning(id, ['Ritz-only half-space classification: legOrder %g < %.1f x relDepth ' ...
                    '(%g). The slab does not resolve the band up to the longitudinal limit and ' ...
                    'the calibrated thresholds do not hold (many slab modes accepted). Use ' ...
                    'legOrder >= %g. (Shown once per session; warning(''on'', ''%s'') to ' ...
                    're-enable.)'], calcParams.legOrder, RitzScore.HS_MIN_LEG_PER_DEPTH, ...
                    calcParams.relDepth, ceil(RitzScore.HS_MIN_LEG_PER_DEPTH * calcParams.relDepth), id);
                warning('off', id);
            end
        end

        function fl = floorAtAngle(structure, angleDeg, precompStruct)
        %%% Median late-time modal-sum level of one direction, recomputed from the structure
        %%% (an independent check of ownReference)
        %%% >> fl = RitzScore.floorAtAngle(structure, angleDeg, precompStruct)
            sr = CalcCommon.rotateStructure(structure, angleDeg);
            rs = RitzSolverSAW(sr, 0);                 % sr already rotated: angle 0
            rs = rs.calc(false, precompStruct);        % no per-mode diagnostics needed here
            if sr.nLay == 1                            % the half-space band (CalcSingleHS)
                L = CalcCommon.limitingAndBulk(sr.mater(1));
                num = RitzScore.modalNumerators(rs, L.longLimVel);
                v = SawSpectrum.velocityAxis(SawSpectrum.AXIS_HI_X_VL * L.longLimVel, RitzScore.VEL_STEP);
                fLo = RitzScore.HS_FLOOR_LO * L.vS; fHi = L.vTf;
            else                                       % the layered band (CalcSingleLayered)
                [vSmin, vLmax] = GreenScoreLayered.bandLimits(sr);
                [fLo, fHi] = GreenScoreLayered.band(vSmin, vLmax);
                num = RitzScore.modalNumerators(rs, vLmax);
                v = RitzScore.layeredAxis(fLo, fHi);
            end
            own = RitzScore.ownReference(structure, rs.vels(1:numel(num)), num, v, fLo, fHi);
            fl = own.floor;
        end
    end
end

function v = traceVel(C9, rho, q)
%%% Trace velocity along x1 of the slowest bulk wave with k_z/k_x = q (slabResolution)
%%% >> v = traceVel(C9, rho, q)
%%% C9 .. elastic tensor reshaped to 9x9 with rows (i,k), columns (j,l), so that
%%% reshape(C9 * kron(n, n), 3, 3) is the Christoffel matrix C_ijkl n_j n_l.
    n = [1; 0; q] / sqrt(1 + q^2);
    G = reshape(C9 * kron(n, n), 3, 3);
    v = sqrt(min(eig((G + G.') / 2)) / rho) * sqrt(1 + q^2);
end

function [E1, r, Ebot] = layeredProfiles(s, D, us)
%%% Energy-weighted near-surface / bottom measures of a layered slab's modes (modeProfiles)
%%% >> [E1, r, Ebot] = layeredProfiles(structure, D, us)
%%% D .. calcDisplacement output [nPts*nLay x nKeep x 3]: nPts points per layer, layer 1
%%% (the deepest) first, each layer from its bottom to its top; us .. surface row [nKeep x 3].
    nP = size(D, 1) / s.nLay;
    z = zeros(nP, s.nLay); rho = z; wz = z;
    zTop = 0;
    for il = s.nLay:-1:1                                  % depth below the surface, top layer first
        h = s.calcParams(il).relDepth;
        z(:, il) = flipud(linspace(zTop, zTop + h, nP)');
        rho(:, il) = s.mater(il).rho;
        dz = h / (nP - 1); w = dz * ones(nP, 1); w([1 end]) = dz / 2;   % trapezoid
        wz(:, il) = w;
        zTop = zTop + h;
    end
    z = z(:); wz = wz(:) .* rho(:); dTot = zTop;
    Pw   = sum(abs(D).^2, 3);                             % [nRows x nKeep]
    Etot = max((wz' * Pw)', realmin);
    E1   = ((wz .* (z <= 1))' * Pw)' ./ Etot;
    r    = abs(us(:, 3)).^2 ./ Etot;
    Ebot = [((wz .* (z >= 0.90 * dTot))' * Pw)', ((wz .* (z >= 0.75 * dTot))' * Pw)'] ./ Etot;
end
