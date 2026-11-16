classdef GreenScoreLayered
%%% GreenScoreLayered — the existence score of the layered 'full' path: each G13_eff
%%% peak is fitted as one damped pole, corrected analytically for a late-time gate of
%%% the wave's own periods and for angular-aperture dephasing, and weighted by its
%%% surface polarisation. Built from the Green's function of the stack alone: no layer
%%% is singled out as "the substrate", and no 12-direction reference is used.
%%%
%%% EXPERIMENTAL: checked against the measured maps published in Grabec et al. 2024 and the
%%% UTGS paper (arXiv:2510.10696); the scores, their scales and the interface may change.
%%%
%%% For one propagation direction, on the band [BAND_LO_X_VSMIN x vSmin,
%%% BAND_HI_X_VLMAX x vLmax] (vSmin / vLmax = slowest shear / fastest longitudinal
%%% limiting velocity over ALL layers -- a range only, never a score input) of the
%%% complex G13_eff at SawSpectrum.ETA:
%%%   1. candidates: local maxima of |G13_eff| inside the band, half-prominence width
%%%      <= GreenScore.MAX_REL_WIDTH x velocity; 3-point complex pole fit as GreenScore
%%%      (pole half-width gamma, floored at the damping width ETA v/2). A failed fit
%%%      (no damped pole within the stencil) keeps the grid peak height.
%%%   2. slope dv/dtheta from two probe solves, at angle +- GreenScore.PROBE_OFFSET, each
%%%      reading the LOCAL maximum of |G| nearest the pole (not the stencil's arg-max,
%%%      which jumps to a stronger neighbour's skirt); both sides combined (mean if
%%%      consistent, else the gentler slope -- a jump only ever steepens). Steps 1-2
%%%      are GreenScore's own (fitPeaks, probeSlope), shared with the half-space.
%%%   3. surface polarisation SH = |r22| / (|r11| + |r22| + |r33|) from the residues of
%%%      the diagonal G11, G22, G33 at the pole (one full-tensor solve, 3 points per
%%%      candidate): ~0 for a sagittal (Rayleigh/Sezawa-type) wave, ~1 for a Love wave.
%%%   4. score = log10(pole height - base)                      the pole prominence,
%%%              + log10 exp(-2 pi gamma tau)                    gated tau = GATE_OWN_PERIODS
%%%              + log10 F(gamma, |dv/dtheta| APERTURE_SIGMA)     periods of the wave itself,
%%%              - N                                             dephased (GreenScore's F),
%%%              + K_POL (1 - 2 SH)                              direction-own normaliser,
%%%                                                              polarisation (SH NaN: -K_POL)
%%%      N = log10 median |G13_eff| over the band of THIS direction.
%%%   5. 'low'  (default): score >= THRESHOLD_LOW  -- sagittal-type waves (Rayleigh, cap,
%%%                        Sezawa families);
%%%      'high': score >= THRESHOLD_HIGH -- also the Love-type waves;
%%%      numeric: score >= value. One scale: 'low' is always a subset of 'high'.
%%%
%%% Why (docs/reference/architecture/forwardCalcs_modeClassification.md, "Layered stacks"):
%%% tuned on one measured TGS map of a film on a substrate (NiTi on MgO, Grabec et al.
%%% 2024): 'low' finds the sagittal-type waves; 'high' also finds most Love-type
%%% waves, at the cost of false positives (sections 11-12 there).
%%% The half-space
%%% GreenScore's band and gate come from the surface material, which in a stack can be
%%% slower than the Rayleigh wave itself; an own-period gate and an all-layer band need
%%% no layer property. A wave's strength alone cannot separate the weak Rayleigh/Love
%%% handover ends from the strongest Love points (both ~0.09 of the Rayleigh line, in
%%% the calculation and the measurement); the polarisation term does. K_POL and the
%%% thresholds rest on that one measured stack.
%%%
%%% Value class with Constant properties and Static methods only.
%%%
%%% GuliweRR v6.0

    properties (Constant)
        THRESHOLD_LOW    = 0.70     % 'low' standard (default): sagittal-type waves
        THRESHOLD_HIGH   = 0.20     % 'high' standard: also Love-type waves
        K_POL            = 0.75     % polarisation weight: score + K_POL (1 - 2 SH)
        GATE_OWN_PERIODS = 10       % late-time gate start, periods of the wave itself
        BAND_LO_X_VSMIN  = 0.3      % band start, x the slowest shear limit of any layer
        BAND_HI_X_VLMAX  = 1.05     % band end, x the fastest longitudinal limit of any layer
    end

    methods (Static)
        function thr = resolveThreshold(t)
        %%% Value of a threshold: 'low' | 'high' | numeric (the score scale)
        %%% >> thr = GreenScoreLayered.resolveThreshold(t)
            thr = CalcCommon.thresholdValue(t, GreenScoreLayered.THRESHOLD_HIGH, ...
                GreenScoreLayered.THRESHOLD_LOW, 'GreenScoreLayered:resolveThreshold:badName', ...
                'GreenScoreLayered.resolveThreshold');
        end

        function tf = accepted(pk, threshold)
        %%% Which scored candidates pass a threshold
        %%% >> tf = GreenScoreLayered.accepted(pk, threshold)
        %%%
        %%% INPUT:  1. pk        .. scored candidates (direction())
        %%%         2. threshold .. 'low' | 'high' | numeric
        %%% OUTPUT: tf .. [numel(pk) x 1] logical, score >= threshold
        %%%
        %%% A NaN score (an empty normalising band) is an error
        %%% (GreenScoreLayered:accepted:noScore), not a silent rejection.
            thr = GreenScoreLayered.resolveThreshold(threshold);
            if isempty(pk); tf = false(0, 1); return; end
            x = [pk.score]';
            if any(isnan(x))
                error('GreenScoreLayered:accepted:noScore', ...
                    'A candidate has no score (normaliser undefined: empty band).');
            end
            tf = x >= thr;
        end

        function [vSmin, vLmax] = bandLimits(structureRotated, lims)
        %%% Slowest shear and fastest longitudinal limiting velocity over ALL layers [km/s]
        %%% >> [vSmin, vLmax] = GreenScoreLayered.bandLimits(structureRotated)
        %%% >> [vSmin, vLmax] = GreenScoreLayered.bandLimits(structureRotated, lims)
        %%%
        %%% INPUT:  1. structureRotated .. the Structure rotated to this direction
        %%%  (opt)  2. lims             .. {1 x nLay} cell, lims{k} = CalcCommon.limitingAndBulk(
        %%%                               structureRotated.mater(k)) already computed by the
        %%%                               caller (CalcSingleLayered.calcLim needs layer 1's for
        %%%                               the substrate fields anyway); default: computed here
        %%% The band's edges only, never a score input: which layer is the "substrate"
        %%% for a wave changes with the wavelength, so no layer is singled out.
            vSmin = Inf; vLmax = -Inf;
            for k = 1:structureRotated.nLay
                if nargin < 2
                    L = CalcCommon.limitingAndBulk(structureRotated.mater(k));
                else
                    L = lims{k};
                end
                vSmin = min(vSmin, min(L.shearLimVels));
                vLmax = max(vLmax, L.longLimVel);
            end
        end

        function [vLo, vHi] = band(vSmin, vLmax)
        %%% Search and normalisation band [km/s] from bandLimits()
        %%% >> [vLo, vHi] = GreenScoreLayered.band(vSmin, vLmax)
            vLo = GreenScoreLayered.BAND_LO_X_VSMIN * vSmin;
            vHi = GreenScoreLayered.BAND_HI_X_VLMAX * vLmax;
        end

        function v = gridAxis(vLo, vMax)
        %%% Uniform Green axis on the lattice SawSpectrum.VEL_MIN + n x SawSpectrum.VEL_STEP, column
        %%% >> v = GreenScoreLayered.gridAxis(vLo, vMax)
        %%%
        %%% INPUT:  1. vLo  .. lowest velocity the axis must reach [km/s] (the band start)
        %%%         2. vMax .. axis end before SawSpectrum.EDGE_TAPER [km/s]
        %%% OUTPUT: v .. from min(VEL_MIN, the lattice point at or below vLo) to vMax +
        %%%              EDGE_TAPER, 4 m/s steps
        %%%
        %%% n may be negative, so a very slow stack (0.3 vSmin below 0.3 km/s) still gets
        %%% its whole band. Every point is VEL_MIN + n*step with the same integer n in a
        %%% single angle and in a sweep, so the band of a sweep's common axis is bit for
        %%% bit the single angle's band (and a single angle's axis is the leading part of
        %%% the sweep's whenever both start at VEL_MIN).
            v0 = SawSpectrum.VEL_MIN; dv = SawSpectrum.VEL_STEP;
            n1 = min(0, floor((vLo - v0) / dv + 1e-9));
            n2 = floor((vMax + SawSpectrum.EDGE_TAPER - v0) / dv + 1e-9);
            v = SawSpectrum.latticeAxis(n1, n2, dv);
        end

        function pk = direction(structureBase, angle, G, velAxis, vSmin, vLmax)
        %%% All steps for one direction: candidates -> probes -> polarisation -> score
        %%% >> pk = GreenScoreLayered.direction(structureBase, angle, G, velAxis, vSmin, vLmax)
        %%%
        %%% INPUT:  1. structureBase .. UNROTATED Structure (probes and the residue solve
        %%%                             rotate it themselves)
        %%%         2. angle         .. direction [deg]
        %%%         3-4. G, velAxis  .. complex G13_eff at SawSpectrum.ETA on a uniform axis
        %%%                             covering the band (this direction's rotated stack)
        %%%         5-6. vSmin, vLmax .. bandLimits() of this direction
        %%% OUTPUT: pk .. scored candidates in velocity order (emptyPeaks() fields),
        %%%               threshold-independent: accepted() decides at any threshold
        %%%
        %%% Cost beyond the axis solve: two probe solves of ~7 points per candidate and one
        %%% full-tensor solve of 3 points per candidate.
            [vLo, vHi] = GreenScoreLayered.band(vSmin, vLmax);
            pk = GreenScoreLayered.candidates(G, velAxis, vLo, vHi);
            if isempty(pk); return; end
            N  = GreenScoreLayered.normaliser(G, velAxis, vLo, vHi);
            pk = GreenScore.probeSlope(structureBase, angle, pk);
            pk = GreenScoreLayered.polarisation(CalcCommon.rotateStructure(structureBase, angle), pk);
            pk = GreenScoreLayered.score(pk, N);
        end

        function pk = candidates(G, velAxis, vLo, vHi)
        %%% Pole-fitted candidates of one direction's complex G13_eff on [vLo, vHi]
        %%% >> pk = GreenScoreLayered.candidates(G, velAxis, vLo, vHi)
        %%%
        %%% INPUT:  1. G       .. [nVel x 1] complex G13_eff at SawSpectrum.ETA
        %%%         2. velAxis .. [nVel x 1] uniform over the band [km/s]
        %%%         3-4. vLo, vHi .. band (band()); peaks must be interior to it
        %%% OUTPUT: pk .. struct array (emptyPeaks() fields; slope, SH and the score parts
        %%%               NaN until GreenScore.probeSlope / polarisation / score)
        %%%
        %%% GreenScore's peak finder (bandPeaks), width filter and pole fit (GreenScore.fitPeaks:
        %%% gamma floored at ETA v/2, a failed fit keeps the grid peak height), on the
        %%% whole band instead of the half-space [0.8 vS, vTf] cap, without fBand.
            pk = GreenScoreLayered.emptyPeaks();
            [found, v, G, ~, h, loc, wid, prom] = GreenScore.bandPeaks(G, velAxis, vLo, vHi, ...
                'GreenScoreLayered:candidates:nonUniform');
            if ~found; return; end
            keep = loc >= vLo & loc <= vHi & wid ./ loc <= GreenScore.MAX_REL_WIDTH;
            h = h(keep); loc = loc(keep); wid = wid(keep); prom = prom(keep);
            f = GreenScore.fitPeaks(v, G, h, loc, wid);
            hMax = max([h; eps]);
            for j = 1:numel(loc)
                pk(j) = struct('vel', f.vel(j), 'polVel', f.polVel(j), 'gamma', f.gamma(j), ...
                    'poleOk', f.poleOk(j), 'height', h(j), 'prominence', prom(j), 'base', h(j) - prom(j), ...
                    'poleHeight', f.poleHeight(j), 'relWidth', wid(j) / f.vel(j), 'relHeight', h(j) / hMax, ...
                    'slope', NaN, 'slopeFwd', NaN, 'slopeBwd', NaN, 'SH', NaN, 'norm', NaN, 'gateTerm', NaN, ...
                    'dephaseTerm', NaN, 'polTerm', NaN, 'logStrength', NaN, 'score', NaN);
            end
        end

        function N = normaliser(G, velAxis, vLo, vHi)
        %%% Direction-own normaliser: log10 of the median |G13_eff| over the band
        %%% >> N = GreenScoreLayered.normaliser(G, velAxis, vLo, vHi)
        %%%
        %%% The band is wide on purpose: keeping part of the region below the slowest
        %%% shear limit makes N track the stack's strongest wave across wavelengths better
        %%% than a half-space-like [0.8 vS, ...] band.
            v = velAxis(:); s = abs(G(:));
            x = s(v >= vLo & v <= vHi & isfinite(s) & s > 0);
            if isempty(x); N = NaN; return; end
            N = log10(median(x));
        end

        function pk = polarisation(structureRotated, pk)
        %%% Surface polarisation SH of every candidate from the diagonal-G residues at its pole
        %%% >> pk = GreenScoreLayered.polarisation(structureRotated, pk)
        %%%
        %%% SH = |r22| / (|r11| + |r22| + |r33|), r_ii the residue of G_ii in the 3-point
        %%% model a1 + r / (x - q) through pole and pole +- FIT_SPAN (the constant absorbs
        %%% a neighbour's skirt). One full-tensor solve, 3 points per candidate. NaN if a
        %%% fit is ill-conditioned (score() then gives the polarisation term -K_POL).
            if isempty(pk); return; end
            n = numel(pk);
            vp = [pk.polVel]'; bad = ~[pk.poleOk]' | ~isfinite(vp); vp(bad) = [pk(bad).vel];
            V = vp' + [-1; 0; 1] * GreenScore.FIT_SPAN;
            g = GreenSolver.calc(structureRotated, 0, V(:).', 1, 3, 'ComputeMinus', false, 'Eta', SawSpectrum.ETA);
            GF = reshape(g.GFull, 3, 3, 3, n);
            for j = 1:n
                pk(j).SH = GreenScoreLayered.surfaceSH(V(:, j), squeeze(GF(1, 1, :, j)), ...
                    squeeze(GF(2, 2, :, j)), squeeze(GF(3, 3, :, j)));
            end
        end

        function SH = surfaceSH(v, g11, g22, g33)
        %%% Shear-horizontal share of a pole from 3 samples of G11, G22, G33 around it
        %%% >> SH = GreenScoreLayered.surfaceSH(v, g11, g22, g33)
        %%%
        %%% INPUT:  v .. [3 x 1] velocities (pole - FIT_SPAN, pole, pole + FIT_SPAN)
        %%%         g11, g22, g33 .. [3 x 1] complex samples
        %%% OUTPUT: SH = |r22| / (|r11| + |r22| + |r33|) in [0, 1] (NaN if a fit is
        %%%         ill-conditioned), r_ii the residue of the 3-point model a1 + r/(x - q)
            [~, ~, ~, r11] = GreenScore.poleFit3(v, g11);
            [~, ~, ~, r22] = GreenScore.poleFit3(v, g22);
            [~, ~, ~, r33] = GreenScore.poleFit3(v, g33);
            r = [r11, r22, r33];
            SH = abs(r(2)) / max(sum(abs(r)), realmin);
        end

        function pk = score(pk, N)
        %%% Gate, dephasing and polarisation terms and the score
        %%% >> pk = GreenScoreLayered.score(pk, N)
        %%%
        %%% INPUT:  1. pk .. candidates() after GreenScore.probeSlope() and polarisation()
        %%%         2. N  .. normaliser() of this direction
        %%% OUTPUT: pk with norm, gateTerm, dephaseTerm, logStrength (log10 of the gated,
        %%%         dephased pole prominence, |G13_eff| units), polTerm, score
        %%%
        %%% Gate tau = GATE_OWN_PERIODS / v (in grating periods per km/s): a pole of
        %%% half-width gamma decays as exp(-2 pi gamma t / lambda). GreenScore.dephaseFactor
        %%% takes the gate as GATE_PERIODS / vRef, so it is called with vRef = GATE_PERIODS
        %%% / tau. A slope NaN (no probe maximum) is set to the probe stencil's bound, as
        %%% in GreenScore.score.
            sLB = GreenScore.PROBE_HALF * GreenScore.PROBE_STEP / GreenScore.PROBE_OFFSET;
            for j = 1:numel(pk)
                tau = GreenScoreLayered.GATE_OWN_PERIODS / pk(j).vel;
                sl = abs(pk(j).slope); if ~isfinite(sl); sl = sLB; end
                g = pk(j).gamma;
                pk(j).gateTerm = -2 * pi * g * tau / log(10);
                F = GreenScore.dephaseFactor(g, sl * GreenScore.APERTURE_SIGMA, GreenScore.GATE_PERIODS / tau);
                pk(j).dephaseTerm = log10(max(F, realmin));
                pk(j).logStrength = log10(max(pk(j).poleHeight - pk(j).base, realmin)) + ...
                    pk(j).gateTerm + pk(j).dephaseTerm;
                pol = GreenScoreLayered.K_POL * (1 - 2 * pk(j).SH);
                if ~isfinite(pol); pol = -GreenScoreLayered.K_POL; end
                pk(j).polTerm = pol;
                pk(j).norm  = N;
                pk(j).score = pk(j).logStrength - N + pol;
            end
        end

        function pk = emptyPeaks()
        %%% The candidate struct, no entries
        %%% >> pk = GreenScoreLayered.emptyPeaks()
        %%%
        %%% Fields: vel (log-parabola |G| peak, the reported velocity), polVel (fitted
        %%% pole), gamma (pole half-width, floored, km/s), poleOk, height, prominence,
        %%% base, poleHeight, relWidth, relHeight (to this direction's tallest), slope /
        %%% slopeFwd / slopeBwd (km/s per deg), SH, norm, gateTerm, dephaseTerm, polTerm,
        %%% logStrength, score.
            pk = struct('vel', {}, 'polVel', {}, 'gamma', {}, 'poleOk', {}, 'height', {}, ...
                'prominence', {}, 'base', {}, 'poleHeight', {}, 'relWidth', {}, 'relHeight', {}, ...
                'slope', {}, 'slopeFwd', {}, 'slopeBwd', {}, 'SH', {}, 'norm', {}, 'gateTerm', {}, ...
                'dephaseTerm', {}, 'polTerm', {}, 'logStrength', {}, 'score', {});
        end
    end
end
