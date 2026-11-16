classdef SawSpectrum
%%% SawSpectrum -- the late-time spectrum machinery the existence scores share: the
%%% damping of the default existence solves (ETA), the axis lattice (VEL_MIN,
%%% latticeAxis, EDGE_TAPER, AXIS_HI_X_VL), the late-time gate (lateTime, GATE_PERIODS)
%%% and the peak search (peaksOf, candidates, MAX_REL_WIDTH). GreenSpotScore and
%%% RitzScore read the gated spectra and the peak search; GreenScore uses the analytic
%%% form of the same gate (GATE_PERIODS) and the width cap.
%%%
%%% EXPERIMENTAL: checked against the measured maps published in Grabec et al. 2024 and the
%%% UTGS paper (arXiv:2510.10696); the scores, their scales and the interface may change.
%%%
%%% These are what make the scores comparable -- the same linewidth, the same axis
%%% start and the same late-time cut -- so they belong to none of the methods.
%%%  - ETA = 3e-3, not the GreenSolver default 3e-4: at 3e-4 a bound Rayleigh pole
%%%    is ~0.5 m/s wide, narrower than the axis step, so its sampled height is grid
%%%    luck and dwarfs a leaky pseudo-SAW of the same residue. A measurement-like
%%%    linewidth conserves the residue and makes heights comparable.
%%%  - Late-time gate (lateTime): surface-wave poles ring for tens of periods;
%%%    bulk-continuum features (limiting-velocity shoulders, Van Hove edges) decay
%%%    within a few. The gate starts GATE_PERIODS periods of vRef = sqrt(G_Hill/rho)
%%%    of the surface layer, so no grating wavelength enters.
%%%  - Zero-padding + tail taper: the gated time signal is zero-padded PAD_FACTOR x,
%%%    which interpolates the late-time spectrum exactly onto a finer axis; a
%%%    raised-cosine fall over the last TAIL_TAPER of the causal window stops a hard
%%%    cut from rippling the spectrum into spurious candidates.
%%%  - candidates: prominence, not height, since most false candidates are bumps on
%%%    the flank of a stronger feature.
%%% See docs/reference/architecture/forwardCalcs_modeClassification.md.
%%%
%%% Value class with Constant properties and Static methods only.
%%%
%%% GuliweRR v6.0

    properties (Constant)
        ETA            = 3e-3    % GreenSolver damping of the default existence solves
        GATE_PERIODS   = 3       % late-time gate start, periods of vRef
        RAMP_PERIODS   = 0.6     % gate turn-on (effectively abrupt, like the measured front cut)
        MAX_REL_WIDTH  = 0.12    % candidate half-prominence width / velocity ceiling
        VEL_MIN        = 0.3     % km/s -- velocity axis start (below any SAW seen)
        VEL_STEP       = 0.004   % km/s -- default uniform axis step (FFT-able)
        PAD_FACTOR     = 4       % zero-padding of the gated time signal = exact band-limited
                                 % interpolation of the late-time spectrum
        TAIL_TAPER     = 0.2     % raised-cosine fall over the last 20 % of the causal time window:
                                 % at 4 m/s the window (1/(2 dv) = 125 ns) ends before long-lived
                                 % ringing has decayed, and a hard cut ripples the spectrum
        EDGE_TAPER     = 0.1     % km/s -- Tukey taper at both ends of the axis
        AXIS_HI_X_VL   = 1.05    % half-space axes run to AXIS_HI_X_VL x longLimVel (+ EDGE_TAPER),
                                 % so a peak at the longitudinal limit is still interior
    end

    methods (Static)
        function velAxis = velocityAxis(vMax, step)
        %%% Uniform velocity axis from VEL_MIN to vMax [km/s], column
        %%% >> velAxis = SawSpectrum.velocityAxis(vMax)
        %%% >> velAxis = SawSpectrum.velocityAxis(vMax, step)
            if nargin < 2 || isempty(step); step = SawSpectrum.VEL_STEP; end
            velAxis = (SawSpectrum.VEL_MIN:step:vMax + SawSpectrum.EDGE_TAPER)';
        end

        function v = latticeAxis(nLo, nHi, step)
        %%% Axis on the lattice VEL_MIN + n x step, n = nLo..nHi, column [km/s]
        %%% >> v = SawSpectrum.latticeAxis(nLo, nHi, step)
        %%%
        %%% The one builder of every Calc / score axis (GreenScore.gridAxis / bandAxis,
        %%% GreenScoreLayered.gridAxis, RitzScore.layeredAxis, SpotGreen.solveAxis, the
        %%% GreenSpotScore direction cache): built as VEL_MIN + n*step with integer n, not
        %%% with the colon operator on velocities (which fills a range from both ends), so a
        %%% point has the same value in a short axis and in a longer one with the same step
        %%% -- what makes a single angle score exactly as inside a sweep. The callers choose
        %%% the n range.
            v = SawSpectrum.VEL_MIN + (nLo:nHi)' * step;
        end

        function [A, vOut] = lateTime(G, velAxis, vRef, pad, tailFrac, gatePeriods)
        %%% |late-time spectrum|: the first GATE_PERIODS periods of vRef removed
        %%% >> [A, vOut] = SawSpectrum.lateTime(G, velAxis, vRef)
        %%% >> [A, vOut] = SawSpectrum.lateTime(G, velAxis, vRef, pad, tailFrac)
        %%% >> [A, vOut] = SawSpectrum.lateTime(G, velAxis, vRef, pad, tailFrac, gatePeriods)
        %%%
        %%% INPUT:  1. G        .. [nVel x nAng] complex spectrum on the uniform velAxis
        %%%         2. velAxis  .. [nVel x 1] km/s
        %%%         3. vRef     .. reference velocity [km/s]
        %%%  (opt)  4. pad      .. zero-padding factor of the gated time signal (default
        %%%                        PAD_FACTOR): A is returned on vOut, pad x finer than
        %%%                        velAxis (exact band-limited interpolation -- the gated
        %%%                        signal has no content past half the window)
        %%%  (opt)  5. tailFrac .. raised-cosine fall over this fraction at the end of the
        %%%                        causal window (default TAIL_TAPER; 0 = hard cut)
        %%%  (opt)  6. gatePeriods .. gate start in periods of vRef (default GATE_PERIODS;
        %%%                        GreenSpotScore reads gates 2 and 4). A vector gives one
        %%%                        spectrum per gate from a single forward FFT (page ig of A)
        %%% The defaults are the Green-path standard. RitzScore passes pad = 1,
        %%% tailFrac = 0: its analytic spectrum is evaluated on its own 2 m/s axis, and its
        %%% thresholds were calibrated without padding or taper.
        %%% OUTPUT: A    .. [pad*nVel x nAng (x numel(gatePeriods))], NaN on the tapered axis ends
        %%%         vOut .. [pad*nVel x 1] km/s, the axis of A
        %%%
        %%% At fixed |k| velocity is frequency: with a formal wavelength of 1 um,
        %%% f [1/ns] = v [km/s] and one period at vRef lasts 1/vRef ns. FFT gives the
        %%% impulse response on one period 1/dv; causal content sits in the first half
        %%% (the retarded Green's function), negative times are discarded. The gate
        %%% (raised-cosine turn-on of RAMP_PERIODS) removes the early, short-lived bulk
        %%% contributions; surface-wave poles survive. Isolated NaNs from GreenSolver
        %%% are interpolated first -- one NaN would poison the whole FFT column.
            nVel = numel(velAxis);
            dv   = velAxis(2) - velAxis(1);
            tAx  = (0:nVel-1)' / (nVel * dv);          % ns at lambda = 1 um
            edgeW = max(2, round(SawSpectrum.EDGE_TAPER / dv));
            w = ones(nVel, 1);
            r = 0.5 * (1 - cos(pi * (0:edgeW-1)' / edgeW));
            w(1:edgeW) = r;  w(end-edgeW+1:end) = flipud(r);
            if any(~isfinite(G(:)))
                G = fillmissing(G, 'linear', 1, 'EndValues', 'nearest');
            end
            S = fft(G .* w, [], 1);
            if nargin < 6 || isempty(gatePeriods); gatePeriods = SawSpectrum.GATE_PERIODS; end
            if nargin < 4 || isempty(pad); pad = SawSpectrum.PAD_FACTOR; end
            if nargin < 5 || isempty(tailFrac); tailFrac = SawSpectrum.TAIL_TAPER; end
            vOut = velAxis(1) + (0:pad * nVel - 1)' * dv / pad;
            wOut = interp1(velAxis(:), w, vOut, 'linear', 0);
            half = tAx < tAx(end) / 2;
            ramp = SawSpectrum.RAMP_PERIODS / vRef;
            nGate = numel(gatePeriods);
            A = zeros(pad * nVel, size(G, 2), nGate);
            for ig = 1:nGate       % one forward FFT serves every gate
                tCut = gatePeriods(ig) / vRef;
                gate = double(tAx >= tCut & half);
                rise = gate > 0 & tAx < tCut + ramp;
                gate(rise) = 0.5 * (1 - cos(pi * (tAx(rise) - tCut) / ramp));
                if tailFrac > 0
                    tEnd = tAx(end) / 2; t0 = tEnd * (1 - tailFrac);
                    fall = gate > 0 & tAx > t0;
                    gate(fall) = gate(fall) .* 0.5 .* (1 + cos(pi * (tAx(fall) - t0) / (tEnd - t0)));
                end
                Sg = S .* gate;
                if pad > 1
                    Sg = [Sg; zeros((pad - 1) * nVel, size(Sg, 2))];
                end
                Ag = abs(ifft(Sg, [], 1)) * pad;
                Ag(wOut < 1, :) = NaN;
                A(:, :, ig) = Ag;
            end
        end

        function P = peaksOf(A, velAxis, vSearchMax)
        %%% Peaks of one late-time spectrum: below vSearchMax, width <= MAX_REL_WIDTH x velocity
        %%% >> P = SawSpectrum.peaksOf(A, velAxis, vSearchMax)
        %%%
        %%% INPUT:  1. A          .. [nVel x 1] |late-time spectrum| (lateTime output, NaN
        %%%                          on the tapered ends)
        %%%         2. velAxis    .. [nVel x 1] its axis [km/s]
        %%%         3. vSearchMax .. search bound [km/s]
        %%% OUTPUT: P .. struct of column vectors vel, prom, height, relWidth
        %%%
        %%% The local maxima (findpeaks, half-prominence widths, prominence floor 1e-4 x
        %%% the largest value) of A below vSearchMax with half-prominence width <=
        %%% MAX_REL_WIDTH x velocity. The search runs to vSearchMax + 0.05 km/s so that a
        %%% peak at the bound is still interior (docs/pitfalls.md MTH4).
            P = struct('vel', zeros(0, 1), 'prom', zeros(0, 1), 'height', zeros(0, 1), ...
                'relWidth', zeros(0, 1));
            m = find(isfinite(A) & velAxis(:) < vSearchMax + 0.05);
            if numel(m) < 3; return; end
            [h, loc, wid, prom] = findpeaks(A(m), velAxis(m), 'WidthReference', 'halfprom', ...
                'MinPeakProminence', 1e-4 * max(A(m)));
            keep = loc < vSearchMax & wid ./ loc <= SawSpectrum.MAX_REL_WIDTH;
            P.vel = loc(keep); P.prom = prom(keep); P.height = h(keep);
            P.relWidth = wid(keep) ./ loc(keep);
        end

        function key = structureKey(structure)
        %%% Memo key of a Structure's physical content (no names, no nMaxModes)
        %%% >> key = SawSpectrum.structureKey(structure)
        %%%
        %%% The byte stream, per layer, of the elastic tensor and its symmetry class, the
        %%% density, dirMat, and the CalcParams relDepth, legOrder and realThickness. The
        %%% score references (GreenScore.reference, GreenSpotScore.reference) are a function
        %%% of these only, so the same stack under another Structure or Material name, or
        %%% with another nMaxModes, hits the memo instead of re-solving 12-132 directions.
        %%% A byte stream, not the numbers: isequal(NaN, NaN) is false, and realThickness is
        %%% NaN when unset.
            c = cell(structure.nLay, 1);
            for il = 1:structure.nLay
                m  = structure.mater(il);
                cp = structure.calcParams(il);
                c{il} = {m.elConsts.tensor, m.elConsts.symmClassNo, m.rho, m.dirMat, ...
                         cp.relDepth, cp.legOrder, cp.realThickness};
            end
            key = getByteStreamFromArray(c);
        end

        function out = memo(cacheName, maxSize, errPrefix, action, key, value)
        %%% Per-session memo shared by the score references: 'get' | 'put' | 'clear' | 'size'
        %%% >> v = SawSpectrum.memo(cacheName, maxSize, errPrefix, 'get', key)   % [] on a miss
        %%% >> SawSpectrum.memo(cacheName, maxSize, errPrefix, 'put', key, value)
        %%%
        %%% INPUT:  1. cacheName .. one memo per name (e.g. 'GreenScore', 'GreenSpotScore')
        %%%         2. maxSize   .. entries kept; the oldest goes first
        %%%         3. errPrefix .. error-id prefix for an unknown action (<errPrefix>:badAction)
        %%%         4. action, key, value
        %%%
        %%% Keys are compared with isequal (the callers pass getByteStreamFromArray of their
        %%% inputs; structures through SawSpectrum.structureKey). A persistent variable: emptied by clear classes, not shared with pool
        %%% workers (the callers do not memoise there).
            persistent store
            if isempty(store); store = struct(); end
            if ~isfield(store, cacheName); store.(cacheName) = struct('keys', {{}}, 'vals', {{}}); end
            c = store.(cacheName);
            out = [];
            switch action
                case 'get'
                    for i = 1:numel(c.keys)
                        if isequal(c.keys{i}, key); out = c.vals{i}; return; end
                    end
                case 'put'
                    c.keys{end+1} = key; c.vals{end+1} = value;
                    if numel(c.keys) > maxSize
                        c.keys(1) = []; c.vals(1) = [];
                    end
                case 'clear'
                    c.keys = {}; c.vals = {};
                case 'size'
                    out = numel(c.keys);
                otherwise
                    error([errPrefix ':badAction'], 'Unknown action ''%s''.', action);
            end
            store.(cacheName) = c;
        end

        function pk = candidates(G, velAxis, vSearchMax, ref, pad, tailFrac)
        %%% Scored candidates (late-time peaks) of one direction's spectrum
        %%% >> pk = SawSpectrum.candidates(G, velAxis, vSearchMax, ref)
        %%% >> pk = SawSpectrum.candidates(G, velAxis, vSearchMax, ref, pad, tailFrac)
        %%%
        %%% INPUT:  1. G          .. [nVel x 1] complex spectrum (e.g. RitzScore's modal sum)
        %%%         2. velAxis    .. [nVel x 1] uniform km/s (SawSpectrum.velocityAxis)
        %%%         3. vSearchMax .. search bound [km/s]
        %%%         4. ref        .. struct: vRef (gate reference velocity, km/s) and floor
        %%%                          (the score's normaliser)
        %%%  (opt)  5-6. pad, tailFrac .. passed to lateTime (defaults: the Green standard,
        %%%                          PAD_FACTOR and TAIL_TAPER; RitzScore passes 1, 0)
        %%% OUTPUT: pk .. struct array (0..N), fields vel, score, prominence, height,
        %%%               relWidth, relHeight (height / tallest candidate here -- for
        %%%               display only, never part of the decision)
        %%%
        %%% Candidates are the peaks of the late-time spectrum (peaksOf).
            if nargin < 5; pad = []; end
            if nargin < 6; tailFrac = []; end
            [A, velAxis] = SawSpectrum.lateTime(G(:), velAxis(:), ref.vRef, pad, tailFrac);
            pk = struct('vel', {}, 'score', {}, 'prominence', {}, 'height', {}, ...
                'relWidth', {}, 'relHeight', {});
            P = SawSpectrum.peaksOf(A, velAxis, vSearchMax);
            hMax = max([P.height; eps]);
            for j = 1:numel(P.vel)
                pk(j) = struct('vel', P.vel(j), 'score', log10(P.prom(j) / ref.floor), ...
                    'prominence', P.prom(j), 'height', P.height(j), 'relWidth', P.relWidth(j), ...
                    'relHeight', P.height(j) / hMax);
            end
        end
    end
end
