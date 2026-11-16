classdef ProbeKernel
%%% ProbeKernel — the wavevectors a measurement couples to: its k-space weight W(k).
%%% EXPERIMENTAL: tried on Fe3Al(316) only, not validated against measurements; the
%%% interface may still change.
%%%
%%% A measurement never excites or detects one wavevector. A finite source pattern and
%%% a finite detection weight couple to a cloud of wavevectors around the nominal +-k0
%%% (k0 = 2 pi / grating period), and the signal is the coherent (complex) sum
%%%
%%%     s(omega) = integral d^2k  W(k) G(k, omega),     W(k) = D^(k) S^(k),
%%%
%%% S^ / D^ the 2D Fourier transforms of the source pattern and of the detection weight
%%% (any phase included, e.g. a detector placed beside the pattern). On a half-space G
%%% has no length scale, G(k, omega) = g(theta, omega/k)/k, so the direction spread
%%% needs one Green solve per direction and the |k| spread only resamples each
%%% direction's spectrum in velocity (SpotGreen). The 1/k cancels the polar area element
%%% k dk dtheta, so the quadrature weight of a node is W x dtheta x d(delta).
%%%
%%% This class only DESCRIBES W; it knows nothing about Green's functions. The
%%% technique enters only through the constructors (ProbeKernel.tgs, ...), so the same
%%% machinery serves any patterned-source / finite-detection measurement.
%%%
%%% Coordinates: W is a function of kappa = k / k0 in the frame of the displayed
%%% direction (kappa_x along it), over the WHOLE plane -- the lobe at kappa = (+1, 0)
%%% and the one at (-1, 0) are both part of the kernel. Each node is evaluated with the
%%% Green component in its own propagation frame (force/displacement axes turn with
%%% k), so a +-k symmetric kernel reproduces the Maznev +-k symmetrised response: a
%%% point kernel with component 'G13' gives exactly G13_eff = (G13_+ - G13_-)/2, the
%%% TGS response every Calc result stores. A kernel that is not +-k symmetric (a
%%% detector on one side) weights the two sides differently, as it should.
%%%
%%% Constructors:
%%%   ProbeKernel.point(component)                       -- no smearing (reference)
%%%   ProbeKernel.gaussian(sigmaRelK, sigmaAngDeg, ...)  -- generic, for fitting
%%%   ProbeKernel.tgs(gratingPeriod, sourceRadii, probeRadii, ...) -- from geometry
%%%   ProbeKernel.sras(gratingPeriod, nFringes, ...)     -- binary mask + detector beside it
%%%   ProbeKernel.fromRealSpace(sourceImg, detectorImg, pixelSize, gratingPeriod, ...)
%%%                                                      -- any pattern, e.g. a camera image
%%%   ProbeKernel(weightFcn, angHalfWidth, relKHalfWidth, ...)     -- any W
%%%
%%% Value class (a copy, not a shared reference, travels with every assignment).
%%% Half-space use only for now (SpotGreen); a layered stack would need extra
%%% solves for the |k| spread, since there |k| sets kd.
%%%
%%% GuliweRR v6.0

    properties (Constant)
        TRUNC_SIGMA    = 3       % Gaussian constructors: nodes out to 3 sigma
        STEPS_PER_SIGMA = 4      % default quadrature: at least 4 nodes per sigma
        MAX_DIR_STEP   = 0.5     % deg -- default direction step ceiling
        NICE_DIR_STEPS = [0.5 0.25 0.2 0.125 0.1 0.05 0.025 0.02 0.0125 0.01 0.005]
                                 % deg -- default direction steps are taken from this
                                 % list (1/step an integer), so integer-degree angles
                                 % share the solved directions in a sweep
        JINC_ZEROS = [3.831705970207512 7.015586669815619 10.173468135062722 ...
                      13.323691936314223 16.470630050877634 19.615858510468243]
                                 % zeros of J1: the rings of a disc's transform
        SRAS_HARMONICS = 7       % sras: mask harmonics |n| <= this summed (their tails
                                 % reach the fundamental lobes)
        SAMPLES_PER_OSC = 8      % sras / fromRealSpace: default quadrature (and stored
                                 % grid) step = the kernel's shortest oscillation / this
    end

    properties (SetAccess = immutable)
        name          char = 'point'   % short label (plots, records)
        weightFcn                      % W(kx, ky): complex weight, kappa = k/k0 in the
                                       % displayed direction's frame, whole plane;
                                       % elementwise on equal-size arrays
        angHalfWidth  double = 0       % deg -- nodes within +-angHalfWidth of 0 and of 180
        relKHalfWidth double = 0       % nodes at |k|/k0 - 1 within +-relKHalfWidth
        dirStep       double = 1       % deg -- default direction step of the quadrature
        relKStep      double = 1       % default relative |k| step of the quadrature
        component     char = 'G13'     % Green component, own propagation frame: 'GIJ',
                                       % I = force, J = displacement (1 = propagation,
                                       % 3 = surface normal)
        params        struct = struct() % the constructor's physical inputs, for the record
    end

    methods
        function obj = ProbeKernel(weightFcn, angHalfWidth, relKHalfWidth, dirStep, relKStep, ...
                component, name, params)
        %%% A kernel from any weight function
        %%% >> k = ProbeKernel(weightFcn, angHalfWidth, relKHalfWidth, dirStep, relKStep)
        %%% >> k = ProbeKernel(weightFcn, angHalfWidth, relKHalfWidth, dirStep, relKStep, component, name, params)
        %%%
        %%% INPUT:  1. weightFcn     .. @(kx, ky) complex W, kappa = k/k0 in the displayed
        %%%                             direction's frame (kx along it), whole plane --
        %%%                             include the -k0 lobe if the measurement has one
        %%%         2. angHalfWidth  .. deg, >= 0: nodes within +-this of 0 and of 180 deg
        %%%         3. relKHalfWidth .. >= 0, < 1: nodes at |kappa| - 1 within +-this
        %%%         4. dirStep       .. deg, > 0: default direction step (ignored if
        %%%                             angHalfWidth = 0)
        %%%         5. relKStep      .. > 0: default |k| step (ignored if relKHalfWidth = 0)
        %%%  (opt)  6. component     .. 'GIJ', own propagation frame (default 'G13')
        %%%  (opt)  7. name          .. label (default 'custom')
        %%%  (opt)  8. params        .. struct of physical inputs, for the record
            if nargin == 0
                obj.weightFcn = @ProbeKernel.pointWeight;
                return
            end
            validateattributes(weightFcn, {'function_handle'}, {'scalar'}, 'ProbeKernel', 'weightFcn');
            validateattributes(angHalfWidth, {'double'}, {'scalar', 'real', 'finite', '>=', 0, '<', 90}, ...
                'ProbeKernel', 'angHalfWidth');
            validateattributes(relKHalfWidth, {'double'}, {'scalar', 'real', 'finite', '>=', 0, '<', 1}, ...
                'ProbeKernel', 'relKHalfWidth');
            validateattributes(dirStep, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
                'ProbeKernel', 'dirStep');
            validateattributes(relKStep, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
                'ProbeKernel', 'relKStep');
            if nargin < 6 || isempty(component); component = 'G13'; end
            if nargin < 7 || isempty(name); name = 'custom'; end
            if nargin < 8 || isempty(params); params = struct(); end
            obj.weightFcn     = weightFcn;
            obj.angHalfWidth  = angHalfWidth;
            obj.relKHalfWidth = relKHalfWidth;
            obj.dirStep       = dirStep;
            obj.relKStep      = relKStep;
            obj.component     = ProbeKernel.validateComponent(component);
            obj.name          = char(name);
            obj.params        = params;
        end

        function q = sample(obj, dirStep, relKStep)
        %%% Quadrature nodes and normalised complex weights of the kernel
        %%% >> q = k.sample()
        %%% >> q = k.sample(dirStep, relKStep)
        %%%
        %%% INPUT:  (opt) 1. dirStep  .. deg (default obj.dirStep)
        %%%         (opt) 2. relKStep .. (default obj.relKStep)
        %%% OUTPUT: q .. struct:
        %%%              offsets  [1 x nO] deg -- node directions relative to the displayed
        %%%                       one: (-m:m) x dirStep around 0, then the same around 180
        %%%              dirIndex [1 x nO] integer -m..m of each offset within its lobe
        %%%              delta    [nD x 1] -- |k|/k0 - 1 = (-n:n)' x relKStep
        %%%              kIndex   [nD x 1] integer -n..n
        %%%              w        [nD x nO] complex, sum(abs(w(:))) = 1
        %%%              dirStep, relKStep -- the steps used
        %%%
        %%% Weight = W(kappa) x dirStep x relKStep at kappa = (1 + delta)(cos o, sin o): the
        %%% polar area element k dk dtheta and G's 1/k cancel (see the class header).
        %%% Normalised by the total |w| so that a ridge flat in angle and |k| keeps its
        %%% height when the kernel has one phase (Gaussian TGS: the result at zero width
        %%% is exactly G_eff); a kernel with a phase ramp is reduced by its own
        %%% incoherence, which is physical. Directions or |k| rows whose weight is below
        %%% 1e-9 of the largest are dropped (they cost a solve each and add nothing).
            if nargin < 2 || isempty(dirStep);  dirStep  = obj.dirStep;  end
            if nargin < 3 || isempty(relKStep); relKStep = obj.relKStep; end
            validateattributes(dirStep, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
                'ProbeKernel.sample', 'dirStep');
            validateattributes(relKStep, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
                'ProbeKernel.sample', 'relKStep');
            m = floor(obj.angHalfWidth / dirStep + 1e-9);
            n = floor(obj.relKHalfWidth / relKStep + 1e-9);
            idx = -m:m;
            offsets  = [idx * dirStep, 180 + idx * dirStep];
            dirIndex = [idx, idx];
            kIndex = (-n:n)';
            delta  = kIndex * relKStep;
            dA = 1;
            if m > 0; dA = dA * dirStep; end
            if n > 0; dA = dA * relKStep; end
            [O, D] = meshgrid(offsets, delta);
            r = 1 + D;
            W = obj.weightFcn(r .* cosd(O), r .* sind(O)) * dA;
            if ~isequal(size(W), size(O)) || any(~isfinite(W(:)))
                error('ProbeKernel:sample:badWeight', ...
                    'weightFcn must return finite values of the input size (got %s for %s).', ...
                    mat2str(size(W)), mat2str(size(O)));
            end
            a = abs(W);
            if max(a(:)) == 0
                error('ProbeKernel:sample:zeroKernel', 'The kernel is zero at every node.');
            end
            keepO = max(a, [], 1) > 1e-9 * max(a(:));
            keepD = max(a, [], 2) > 1e-9 * max(a(:));
            W = W(keepD, keepO);
            q = struct('offsets', offsets(keepO), 'dirIndex', dirIndex(keepO), ...
                'delta', delta(keepD), 'kIndex', kIndex(keepD), ...
                'w', W / sum(abs(W(:))), 'dirStep', dirStep, 'relKStep', relKStep);
        end
    end

    methods (Static)
        function k = point(component)
        %%% No smearing: the nominal +-k0 only (the reference, = the plain Green map)
        %%% >> k = ProbeKernel.point()            % 'G13': gives G13_eff
        %%% >> k = ProbeKernel.point('G33')
            if nargin < 1 || isempty(component); component = 'G13'; end
            k = ProbeKernel(@ProbeKernel.pointWeight, 0, 0, 1, 1, component, 'point', struct());
        end

        function k = gaussian(sigmaRelK, sigmaAngDeg, varargin)
        %%% A +-k symmetric Gaussian kernel in |k| and direction (polar form)
        %%% >> k = ProbeKernel.gaussian(sigmaRelK, sigmaAngDeg)
        %%% >> k = ProbeKernel.gaussian(sigmaRelK, sigmaAngDeg, 'Component', 'G13', 'DirStep', 0.25, 'RelKStep', 1e-3)
        %%%
        %%% INPUT:  1. sigmaRelK   .. std of |k|/k0 (relative velocity smearing), >= 0
        %%%         2. sigmaAngDeg .. std of the direction [deg], >= 0
        %%%  (opt)  Name-Value:
        %%%         'Component'    .. default 'G13'
        %%%         'DirStep', 'RelKStep' .. quadrature steps (default: >= STEPS_PER_SIGMA
        %%%                 nodes per sigma, direction step from NICE_DIR_STEPS; SpotGreen
        %%%                 further caps the |k| step at its damping / 4)
        %%%
        %%% W = exp(-(|kappa| - 1)^2 / (2 sigmaRelK^2) - dphi^2 / (2 sigmaAng^2)) around
        %%% each of the two lobes (dphi from 0 or 180 deg); a zero sigma is a delta in
        %%% that coordinate. Nodes out to TRUNC_SIGMA sigma.
            validateattributes(sigmaRelK, {'double'}, {'scalar', 'real', 'finite', '>=', 0, '<', 0.3}, ...
                'ProbeKernel.gaussian', 'sigmaRelK');
            validateattributes(sigmaAngDeg, {'double'}, {'scalar', 'real', 'finite', '>=', 0, '<', 20}, ...
                'ProbeKernel.gaussian', 'sigmaAngDeg');
            [component, dirStep, relKStep] = ProbeKernel.quadOptions(varargin, 'ProbeKernel.gaussian');
            if isempty(dirStep);  dirStep  = ProbeKernel.defaultDirStep(sigmaAngDeg); end
            if isempty(relKStep); relKStep = ProbeKernel.defaultRelKStep(sigmaRelK); end
            fcn = @(kx, ky) ProbeKernel.gaussPolar(kx, ky, sigmaRelK, sigmaAngDeg);
            p = struct('sigmaRelK', sigmaRelK, 'sigmaAngDeg', sigmaAngDeg);
            k = ProbeKernel(fcn, ProbeKernel.TRUNC_SIGMA * sigmaAngDeg, ...
                ProbeKernel.TRUNC_SIGMA * sigmaRelK, dirStep, relKStep, component, ...
                sprintf('gaussian(%.3g, %.3g deg)', sigmaRelK, sigmaAngDeg), p);
        end

        function k = tgs(gratingPeriod, sourceRadii, probeRadii, varargin)
        %%% Transient-grating kernel from the measurement geometry
        %%% >> k = ProbeKernel.tgs(gratingPeriod, sourceRadii, probeRadii)
        %%% >> k = ProbeKernel.tgs(gratingPeriod, sourceRadii, probeRadii, 'Component', 'G13', 'DirStep', 0.25)
        %%%
        %%% INPUT:  1. gratingPeriod .. Lambda (any length unit)
        %%%         2. sourceRadii   .. 1/e^2 INTENSITY radius of the excitation spot, same
        %%%                             unit: scalar (circular) or [along k, across k]
        %%%         3. probeRadii    .. the same for the detection (probe/reference) spot;
        %%%                             [] = a detection weight much wider than the source
        %%%  (opt)  Name-Value 'Component', 'DirStep', 'RelKStep' .. as ProbeKernel.gaussian
        %%%
        %%% Interference fringes under a Gaussian envelope, I = E (1 + cos k0 x), have the
        %%% spectrum E^(k -+ k0)/2: a Gaussian cloud around each of +-k0, with
        %%% sigma_k = 2 / radius for E = exp(-2 x^2/a^2 - 2 y^2/b^2). The detection weight
        %%% (probe field x reference field, both from the probe spot) multiplies it, so
        %%% 1/sigma^2 adds: in units of k0,
        %%%     sigmaRelK = Lambda / (pi sqrt(a_s^2 + a_p^2))   (|k| spread, along k)
        %%%     sigmaAng  = Lambda / (pi sqrt(b_s^2 + b_p^2))   [rad] (direction spread)
        %%% The exact Cartesian Gaussian is used (not its polar approximation). Both
        %%% lobes have the same weight: the zero-width limit is the G13_eff convention.
        %%% UNCHECKED against a measured TGS map; the detection-weight model (heterodyne,
        %%% reference sharing the probe spot) is an assumption.
            validateattributes(gratingPeriod, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
                'ProbeKernel.tgs', 'gratingPeriod');
            [aS, bS] = ProbeKernel.radiiPair(sourceRadii, 'sourceRadii');
            if nargin < 3 || isempty(probeRadii)
                aP = 0; bP = 0; probeRadii = [];
            else
                [aP, bP] = ProbeKernel.radiiPair(probeRadii, 'probeRadii');
            end
            [component, dirStep, relKStep] = ProbeKernel.quadOptions(varargin, 'ProbeKernel.tgs');
            sx = gratingPeriod / (pi * sqrt(aS^2 + aP^2));    % sigma of kappa_x
            sy = gratingPeriod / (pi * sqrt(bS^2 + bP^2));    % sigma of kappa_y
            if sx >= 0.3 || sy >= 0.3
                error('ProbeKernel:tgs:tooFewFringes', ...
                    'The spot is too small for its grating period (sigma = %.2g, %.2g of k0).', sx, sy);
            end
            sigAngDeg = rad2deg(sy);
            if isempty(dirStep);  dirStep  = ProbeKernel.defaultDirStep(sigAngDeg); end
            if isempty(relKStep); relKStep = ProbeKernel.defaultRelKStep(sx); end
            fcn = @(kx, ky) ProbeKernel.gaussCartesian(kx, ky, sx, sy);
            p = struct('gratingPeriod', gratingPeriod, 'sourceRadii', [aS bS], ...
                'probeRadii', probeRadii, 'sigmaRelK', sx, 'sigmaAngDeg', sigAngDeg);
            k = ProbeKernel(fcn, ProbeKernel.TRUNC_SIGMA * sigAngDeg, ProbeKernel.TRUNC_SIGMA * sx, ...
                dirStep, relKStep, component, ...
                sprintf('tgs(%.3g, %.3g%% / %.3g deg)', gratingPeriod, 100 * sx, sigAngDeg), p);
        end

        function k = sras(gratingPeriod, nFringes, varargin)
        %%% SRAS kernel: binary fringes in a disc, point-like detector beside the pattern
        %%% EXPERIMENTAL: not yet compared with a measured SRAS map; the source,
        %%% detection and beam models are assumptions (see the end of this help).
        %%% >> k = ProbeKernel.sras(gratingPeriod, nFringes)
        %%% >> k = ProbeKernel.sras(gratingPeriod, nFringes, 'DetectorOffset', [x y], 'KnifeEdgeDir', 0, ...)
        %%%
        %%% INPUT:  1. gratingPeriod .. Lambda (any length unit)
        %%%         2. nFringes      .. fringe periods across the disc diameter: R = nFringes
        %%%                             Lambda / 2
        %%%  (opt)  Name-Value (lengths in the unit of gratingPeriod):
        %%%         'DutyCycle'      .. bright fraction of a period, (0, 1) (default 0.5)
        %%%         'DetectorOffset' .. [along k, across k] from the disc centre (default
        %%%                             [1.1 R, 0]: on the centre line, beside the disc)
        %%%         'DetectorRadius' .. 1/e^2 intensity radius of the probe spot (finite, default 0)
        %%%         'KnifeEdgeDir'   .. deg from k: the knife edge measures the surface slope
        %%%                             along this direction (default 0); NaN = the detector
        %%%                             measures the displacement itself (no slope factor)
        %%%         'KFactors'       .. true (default): amplitude factors |k| of the
        %%%                             thermoelastic force (i k T) and of the slope
        %%%                             (i k.n); false: their directions only
        %%%         'EdgeBlur'       .. std of a Gaussian blur of the imaged mask (default 0)
        %%%         'BeamRadius'     .. 1/e^2 intensity radius of a Gaussian illuminating beam,
        %%%                             centred on the disc (default Inf: flat-top)
        %%%         'Rings'          .. Airy rings kept around each lobe (default 2): the
        %%%                             lobe is cut on a zero of the disc's transform
        %%%         'Component'      .. default 'G13'
        %%%         'DirStep', 'RelKStep' .. quadrature steps (default: SAMPLES_PER_OSC nodes
        %%%                             per shortest oscillation of W; SpotGreen caps the |k|
        %%%                             step at its damping / 4)
        %%%
        %%% Source: the mask's Fourier series, sq(x) = D + sum_n sin(n pi D)/(n pi) e^{i n k0 x}
        %%% (bright stripe centred on the disc), times the disc, so each harmonic is an Airy
        %%% lobe A(|k - n k0 x^|), A(q) = 2 J1(qR)/(qR). Harmonics |n| <= SRAS_HARMONICS are
        %%% summed, so the tails of the DC and 3k0 lobes reach the fundamental lobes as they
        %%% should; the 3k0 lobe itself (apparent velocity 3v) and, for D ~= 0.5, the 2k0 one
        %%% (apparent 2v) are outside the kernel (relKHalfWidth < 1) and not modelled.
        %%% A Gaussian beam (BeamRadius w) makes the envelope exp(-2 r^2 / w^2) inside the disc:
        %%% each lobe is then its Hankel transform (tabulated once per kernel, normalised to 1
        %%% at q = 0) instead of the Airy pattern. A soft envelope is what lowers the Airy rings,
        %%% i.e. the side lobes of the measured line; EdgeBlur cannot do it without washing out
        %%% the fringes (it acts at |k| ~ k0, nearly constant over the lobe).
        %%% Detection: a probe spot at x_d gives D^(-k) = e^{i k.x_d} exp(-k^2 r_d^2 / 8), the
        %%% knife edge a factor i k.n. The phase ramp selects the lobe whose packet travels
        %%% towards the detector: the other lobe cancels coherently, so SRAS sees the
        %%% one-sided G(+k) towards the detector, not the symmetrised G13_eff. W is cut to
        %%% zero outside a circle of radius JINC_ZEROS(Rings + 1) / (k0 R) around each lobe
        %%% (sum |W| over Airy rings does not converge, so the cut fixes the normalisation),
        %%% or 4 sigma of the beam's Gaussian lobe, Lambda / (pi w), if that is larger.
        %%% Depends only on nFringes, x_d / Lambda and r_d / Lambda: a half-space has no
        %%% length scale. ASSUMPTIONS (unchecked against a measured SRAS map): thermoelastic
        %%% in-plane force along k (own-frame G13), knife edge = surface slope, |k| amplitude
        %%% factors, Gaussian beam (BeamRadius) centred on the disc.
            ip = inputParser;
            ip.FunctionName = 'ProbeKernel.sras';
            ip.addParameter('DutyCycle', 0.5, @(x) isscalar(x) && isreal(x) && x > 0 && x < 1);
            ip.addParameter('DetectorOffset', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 2 && all(isfinite(x))));
            ip.addParameter('DetectorRadius', 0, @(x) isscalar(x) && isreal(x) && isfinite(x) && x >= 0);
            ip.addParameter('KnifeEdgeDir', 0, @(x) isscalar(x) && isreal(x));
            ip.addParameter('KFactors', true, @(x) isscalar(x) && (islogical(x) || isnumeric(x)));
            ip.addParameter('EdgeBlur', 0, @(x) isscalar(x) && isreal(x) && x >= 0);
            ip.addParameter('BeamRadius', Inf, @(x) isscalar(x) && isreal(x) && x > 0);
            ip.addParameter('Rings', 2, @(x) isscalar(x) && x == round(x) && x >= 0 && ...
                x < numel(ProbeKernel.JINC_ZEROS));
            ip.addParameter('Component', 'G13');
            ip.addParameter('DirStep', []);
            ip.addParameter('RelKStep', []);
            ip.parse(varargin{:});
            r = ip.Results;
            validateattributes(gratingPeriod, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
                'ProbeKernel.sras', 'gratingPeriod');
            validateattributes(nFringes, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
                'ProbeKernel.sras', 'nFringes');
            R = nFringes * gratingPeriod / 2;
            if isempty(r.DetectorOffset); r.DetectorOffset = [1.1 * R, 0]; end
            k0 = 2 * pi / gratingPeriod;
            Rt = k0 * R;                                   % = pi nFringes
            rho = max(ProbeKernel.JINC_ZEROS(r.Rings + 1) / Rt, 4 * gratingPeriod / (pi * r.BeamRadius));
            if rho >= 0.9
                error('ProbeKernel:sras:tooFewFringes', ...
                    ['%g fringes are too few for %d rings (or the beam too narrow): the lobe ' ...
                     'would reach %.2g k0 from k0 (limit 0.9).'], nFringes, r.Rings, rho);
            end
            if abs(r.DutyCycle - 0.5) > 1e-12
                warning('ProbeKernel:sras:evenHarmonics', ...
                    ['Duty cycle %.3g: the mask has even harmonics; the 2k0 lobe (apparent ' ...
                     'velocity 2v) is not part of the kernel.'], r.DutyCycle);
            end
            n = -ProbeKernel.SRAS_HARMONICS:ProbeKernel.SRAS_HARMONICS;
            c = sin(n * pi * r.DutyCycle) ./ (n * pi);
            c(n == 0) = r.DutyCycle;
            if isinf(r.BeamRadius)
                lobe = @ProbeKernel.jinc;
            else
                lobe = ProbeKernel.gaussDiscLobe(R / r.BeamRadius, Rt * (max(abs(n)) + 1 + rho));
            end
            xdt = k0 * r.DetectorOffset(:).';
            p = struct('gratingPeriod', gratingPeriod, 'nFringes', nFringes, 'radius', R, ...
                'dutyCycle', r.DutyCycle, 'detectorOffset', r.DetectorOffset(:).', ...
                'detectorRadius', r.DetectorRadius, 'knifeEdgeDir', r.KnifeEdgeDir, ...
                'kFactors', logical(r.KFactors), 'edgeBlur', r.EdgeBlur, ...
                'beamRadius', r.BeamRadius, 'rings', r.Rings, 'lobeRadius', rho);
            fcn = @(kx, ky) ProbeKernel.srasWeight(kx, ky, Rt, n, c, lobe, xdt, ...
                k0 * r.DetectorRadius, k0 * r.EdgeBlur, r.KnifeEdgeDir, logical(r.KFactors), rho);
            %%% shortest oscillation of W in kappa: ring spacing ~ pi / Rt, phase ramp 2 pi / |xdt|
            dq = min(pi / Rt, 2 * pi / max(norm(xdt), eps)) / ProbeKernel.SAMPLES_PER_OSC;
            [dirStep, relKStep] = ProbeKernel.oscSteps(dq, r.DirStep, r.RelKStep);
            k = ProbeKernel(fcn, asind(rho), rho, dirStep, relKStep, r.Component, ...
                sprintf('sras(%g, %g fringes, det %.2g R)', gratingPeriod, nFringes, ...
                r.DetectorOffset(1) / R), p);
        end

        function k = fromRealSpace(sourceImg, detectorImg, pixelSize, gratingPeriod, varargin)
        %%% Kernel from real-space images of the source pattern and of the detection weight
        %%% EXPERIMENTAL: checked against ProbeKernel.sras only.
        %%% >> k = ProbeKernel.fromRealSpace(sourceImg, detectorImg, pixelSize, gratingPeriod)
        %%% >> k = ProbeKernel.fromRealSpace(..., 'LobeRadius', 0.4, 'KnifeEdgeDir', 0, ...)
        %%%
        %%% INPUT:  1. sourceImg   .. real matrix, source intensity S(y, x) (rows = y, columns
        %%%                           = x; the grating vector k0 along +x)
        %%%         2. detectorImg .. real matrix of the same size, detection weight D(y, x)
        %%%                           (e.g. the probe spot); the two share one pixel grid,
        %%%                           whose origin does not matter
        %%%         3. pixelSize   .. pixel pitch (any length unit)
        %%%         4. gratingPeriod .. Lambda, same unit
        %%%  (opt)  Name-Value:
        %%%         'LobeRadius'   .. W is kept within this distance (in k0) of +-k0 (default
        %%%                           0.4; < 0.9)
        %%%         'KnifeEdgeDir', 'KFactors', 'Component', 'DirStep', 'RelKStep' .. as
        %%%                           ProbeKernel.sras
        %%%
        %%% W(k) = S^(k) D^(-k), with S^(k) = sum S e^{-i k.r} evaluated by a separable direct
        %%% DFT on a Cartesian grid around each lobe (no padding, exact at the grid points),
        %%% then interpolated (cubic) for the quadrature nodes. The grid and default steps
        %%% resolve the shortest oscillation the images can produce, Lambda / L_support (L =
        %%% diagonal of the box holding every pixel above 1e-3 of either image's maximum),
        %%% with SAMPLES_PER_OSC points. Knife-edge and |k| factors are applied analytically
        %%% on top, as in sras. Pixels are treated as point samples.
            ip = inputParser;
            ip.FunctionName = 'ProbeKernel.fromRealSpace';
            ip.addParameter('LobeRadius', 0.4, @(x) isscalar(x) && isreal(x) && x > 0 && x < 0.9);
            ip.addParameter('KnifeEdgeDir', 0, @(x) isscalar(x) && isreal(x));
            ip.addParameter('KFactors', true, @(x) isscalar(x) && (islogical(x) || isnumeric(x)));
            ip.addParameter('Component', 'G13');
            ip.addParameter('DirStep', []);
            ip.addParameter('RelKStep', []);
            ip.parse(varargin{:});
            r = ip.Results;
            validateattributes(sourceImg, {'numeric'}, {'2d', 'real', 'finite', 'nonempty'}, ...
                'ProbeKernel.fromRealSpace', 'sourceImg');
            validateattributes(detectorImg, {'numeric'}, {'2d', 'real', 'finite', 'size', size(sourceImg)}, ...
                'ProbeKernel.fromRealSpace', 'detectorImg');
            validateattributes(pixelSize, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
                'ProbeKernel.fromRealSpace', 'pixelSize');
            validateattributes(gratingPeriod, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
                'ProbeKernel.fromRealSpace', 'gratingPeriod');
            S = double(sourceImg); D = double(detectorImg);
            if max(abs(S(:))) == 0 || max(abs(D(:))) == 0
                error('ProbeKernel:fromRealSpace:emptyImage', 'Both images must be nonzero somewhere.');
            end
            mask = abs(S) > 1e-3 * max(abs(S(:))) | abs(D) > 1e-3 * max(abs(D(:)));
            [iy, ix] = find(mask);
            L = pixelSize * hypot(max(ix) - min(ix) + 1, max(iy) - min(iy) + 1);
            k0 = 2 * pi / gratingPeriod;
            dq = gratingPeriod / L / ProbeKernel.SAMPLES_PER_OSC;
            rho = r.LobeRadius;
            [ny, nx] = size(S);
            x = ((0:nx-1) * pixelSize).'; y = (0:ny-1) * pixelSize;
            ky = (-ceil(rho / dq):ceil(rho / dq)) * dq;
            kxRel = ky;
            lobes = cell(1, 2);
            for s = [1 2]
                kx = (3 - 2 * s) + kxRel;                         % +1 lobe, then -1
                Ey = exp(-1i * k0 * ky(:) * y);                   % [nky x ny]
                Ex = exp(-1i * k0 * x * kx);                      % [nx x nkx]
                W = (Ey * S * Ex) .* (conj(Ey) * D * conj(Ex));   % S^(k) D^(-k), D real
                lobes{s} = griddedInterpolant({ky, kx}, W, 'cubic', 'none');
            end
            fcn = @(kxq, kyq) ProbeKernel.gridWeight(kxq, kyq, lobes, r.KnifeEdgeDir, ...
                logical(r.KFactors), rho);
            [dirStep, relKStep] = ProbeKernel.oscSteps(dq, r.DirStep, r.RelKStep);
            p = struct('pixelSize', pixelSize, 'gratingPeriod', gratingPeriod, ...
                'imageSize', [ny nx], 'supportDiagonal', L, 'gridStep', dq, ...
                'knifeEdgeDir', r.KnifeEdgeDir, 'kFactors', logical(r.KFactors), 'lobeRadius', rho);
            k = ProbeKernel(fcn, asind(rho), rho, dirStep, relKStep, r.Component, ...
                sprintf('fromRealSpace(%g, %dx%d px)', gratingPeriod, ny, nx), p);
        end

        function component = validateComponent(component)
        %%% 'GIJ' with I, J in 1..3 (own propagation frame, one-sided notation)
        %%% >> component = ProbeKernel.validateComponent('G13')
            if ~(ischar(component) || (isstring(component) && isscalar(component))) || ...
                    isempty(regexp(char(component), '^G[123][123]$', 'once'))
                error('ProbeKernel:validateComponent:badComponent', ...
                    ['component must be ''GIJ'' with I, J in 1..3 (own propagation frame; the ' ...
                     'kernel''s two lobes make the +-k combination) -- got ''%s''.'], char(component));
            end
            component = char(component);
        end

        function s = defaultDirStep(sigmaAngDeg)
        %%% The largest NICE_DIR_STEPS entry <= sigma / STEPS_PER_SIGMA (and <= MAX_DIR_STEP)
            s = ProbeKernel.niceDirStep(sigmaAngDeg / ProbeKernel.STEPS_PER_SIGMA);
            if sigmaAngDeg == 0; s = 1; end
        end

        function s = defaultRelKStep(sigmaRelK)
        %%% sigma / STEPS_PER_SIGMA (SpotGreen.options caps it at Eta / ETA_STEPS)
            s = sigmaRelK / ProbeKernel.STEPS_PER_SIGMA;
            if sigmaRelK == 0; s = 1; end
        end
    end

    methods (Static, Access = private)
        function [component, dirStep, relKStep] = quadOptions(args, fname)
        %%% 'Component', 'DirStep', 'RelKStep' Name-Value options of gaussian / tgs
            ip = inputParser;
            ip.FunctionName = fname;
            ip.addParameter('Component', 'G13');
            ip.addParameter('DirStep', []);
            ip.addParameter('RelKStep', []);
            ip.parse(args{:});
            component = ip.Results.Component;
            dirStep   = ip.Results.DirStep;
            relKStep  = ip.Results.RelKStep;
        end

        function W = srasWeight(kx, ky, Rt, n, c, lobe, xdt, rdt, blurt, knifeDir, kFactors, rho)
        %%% sras's W (see there); kappa units, Rt = k0 R, xdt = k0 x_d, rdt = k0 r_d, lobe(qR)
        %%% the envelope's transform (jinc for a flat-top disc)
            S = zeros(size(kx));
            for i = 1:numel(n)
                S = S + c(i) * lobe(Rt * hypot(kx - n(i), ky));
            end
            k2 = kx.^2 + ky.^2;
            W = S .* exp(-k2 * (blurt^2 / 2 + rdt^2 / 8)) .* exp(1i * (kx * xdt(1) + ky * xdt(2)));
            W = ProbeKernel.detectionFactors(W, kx, ky, knifeDir, kFactors, rho);
        end

        function W = gridWeight(kx, ky, lobes, knifeDir, kFactors, rho)
        %%% fromRealSpace's W: the interpolated DFT of the nearer lobe, then the factors
            W = complex(zeros(size(kx)));
            plus = kx >= 0;
            W(plus)  = lobes{1}(ky(plus), kx(plus));
            W(~plus) = lobes{2}(ky(~plus), kx(~plus));
            W(isnan(W)) = 0;
            W = ProbeKernel.detectionFactors(W, kx, ky, knifeDir, kFactors, rho);
        end

        function W = detectionFactors(W, kx, ky, knifeDir, kFactors, rho)
        %%% Knife-edge slope i k.n (NaN direction: none), force |k| (kFactors), lobe cut
            kr = hypot(kx, ky);
            if ~isnan(knifeDir)
                kn = kx * cosd(knifeDir) + ky * sind(knifeDir);
                if ~kFactors; kn = kn ./ max(kr, eps); end
                W = W .* (1i * kn);
            end
            if kFactors; W = W .* kr; end
            W(min(hypot(kx - 1, ky), hypot(kx + 1, ky)) > rho) = 0;
        end

        function f = gaussDiscLobe(aR, xMax)
        %%% Transform of the envelope exp(-2 r^2 / w^2) circ(r / R), as a function of x = qR,
        %%% normalised to 1 at x = 0; aR = R / w. Tabulated (step 0.05, ~16 points per
        %%% oscillation) by adaptive quadrature of the Hankel integral
        %%% int_0^1 exp(-2 aR^2 s^2) J0(x s) s ds, then a spline.
            x = (0:0.05:xMax + 1).';
            g = @(s) exp(-2 * aR^2 * s.^2) .* s;
            E = integral(@(s) g(s) .* besselj(0, x * s), 0, 1, 'ArrayValued', true, ...
                'AbsTol', 1e-12, 'RelTol', 1e-9);
            E = E / integral(g, 0, 1, 'AbsTol', 1e-14, 'RelTol', 1e-12);
            F = griddedInterpolant(x, E, 'spline', 'nearest');
            f = @(q) F(q);
        end

        function j = jinc(x)
        %%% 2 J1(x) / x, = 1 at x = 0
            j = ones(size(x));
            nz = x ~= 0;
            j(nz) = 2 * besselj(1, x(nz)) ./ x(nz);
        end

        function [dirStep, relKStep] = oscSteps(dq, dirStep, relKStep)
        %%% Default steps from the shortest oscillation step dq (kappa): the direction step
        %%% the largest NICE_DIR_STEPS entry <= dq [rad] (and <= MAX_DIR_STEP)
            if isempty(dirStep)
                dirStep = ProbeKernel.niceDirStep(rad2deg(dq));
            end
            if isempty(relKStep); relKStep = dq; end
        end

        function s = niceDirStep(target)
        %%% The largest NICE_DIR_STEPS entry <= min(target, MAX_DIR_STEP) [deg]; the smallest
        %%% entry if none is
            ok = ProbeKernel.NICE_DIR_STEPS(ProbeKernel.NICE_DIR_STEPS <= ...
                min(ProbeKernel.MAX_DIR_STEP, target) + 1e-12);
            if isempty(ok); s = ProbeKernel.NICE_DIR_STEPS(end); else; s = ok(1); end
        end

        function W = pointWeight(kx, ky)
        %%% 1 at kappa = (+-1, 0), 0 elsewhere
            W = double(abs(abs(kx) - 1) < 1e-12 & abs(ky) < 1e-12);
        end

        function W = gaussPolar(kx, ky, sR, sA)
        %%% Polar Gaussian around both lobes (dphi from the nearer of 0 / 180 deg)
            r = hypot(kx, ky);
            phi = atan2d(ky, kx);                       % (-180, 180]
            dphi = min(abs(phi), 180 - abs(phi));       % distance to the nearer lobe axis
            W = ProbeKernel.gauss1(r - 1, sR) .* ProbeKernel.gauss1(dphi, sA);
        end

        function W = gaussCartesian(kx, ky, sx, sy)
        %%% Cartesian Gaussian clouds at kappa = (+1, 0) and (-1, 0)
            W = (ProbeKernel.gauss1(kx - 1, sx) + ProbeKernel.gauss1(kx + 1, sx)) .* ...
                ProbeKernel.gauss1(ky, sy);
        end

        function g = gauss1(x, s)
        %%% exp(-x^2/(2 s^2)); a delta (1 at x = 0) for s = 0
            if s == 0
                g = double(abs(x) < 1e-9);
            else
                g = exp(-x.^2 / (2 * s^2));
            end
        end

        function [a, b] = radiiPair(r, argName)
            validateattributes(r, {'double'}, {'vector', 'real', 'finite', 'positive'}, ...
                'ProbeKernel.tgs', argName);
            if isscalar(r); a = r; b = r;
            elseif numel(r) == 2; a = r(1); b = r(2);
            else
                error('ProbeKernel:tgs:badRadii', '%s must be a scalar or [along, across].', argName);
            end
        end
    end
end
