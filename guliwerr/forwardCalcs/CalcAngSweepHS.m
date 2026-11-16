classdef CalcAngSweepHS < handle & matlab.mixin.Copyable
%%% CalcAngSweepHS — angular sweep of CalcSingleHS objects for the homogeneous
%%% half-space forward calculation. Drives Ritz -> Christoffel -> Green across a
%%% propagation-angle array so that one common Green velocity axis serves every
%%% angle (one Green solve per angle) — not achievable by simply looping
%%% CalcSingleHS.calcAll(), whose axis depends on each angle's own longLimVel. See
%%% docs/reference/architecture/forwardCalcs.md, "Angle-sweep orchestrators" and "The
%%% 0-180 degree angle grid, and Green's function on it".
%%%
%%% GuliweRR v6.0

    properties (SetAccess = protected)
        structure       Structure   = Structure()      % UNROTATED base structure -- each
                                                          % cs{iAngle} rotates its own internal copy
        angles          double      = (0:2:180)'          % [nAngles x 1] propagation angles [deg].
                                                          % 0-180 only, not a full circle: v(theta)
                                                          % is exactly even in propagation direction,
                                                          % Christoffel being quadratic in n

        %%% calcFlags
        calcMode        char        = 'full'             % 'full' (default) | 'ritzOnly'
        retentionLevel  char        = 'full'             % 'full' (default) | 'saw_only' | 'none'
        nVelGreen(1,1)  double      = NaN                % number of points of the common Green
                                                          % axis (sets the step, capped at
                                                          % GreenScore.MAX_STEP); NaN (default) =
                                                          % use greenVelStep
        useParallel(1,1) string {mustBeMember(useParallel, ["auto","always","never"])} = "auto"
                                                          % "auto" (default, parallel once nAngles >=
                                                          % ParallelPolicy.AUTO_THRESHOLD) | "always"
                                                          % | "never" -- applies uniformly to every
                                                          % phase loop; thread-pool based, not
                                                          % process-based

        %%% cs -- per-angle CalcSingleHS array
        cs              cell        = {}                 % {nAngles x 1} of CalcSingleHS objects

        %%% maps -- greenVelAxis + per-angle assembled distributions, filled by
        %%% assembleMaps() at the end of calc(), not by a logic-bearing post-loop
        greenVelAxis    double      = []                 % [nVel x 1] -- ONE common Green axis
                                                          % (full mode only): 0.3 km/s to 1.05 x
                                                          % max over angles of longLimVel + 0.1
                                                          % km/s, in greenVelStep steps. Each angle's
                                                          % own axis is
                                                          % its leading part
        longVelMap            double = []   % [nAngles x 1] -- full mode only
        bulkVelMap             double = []   % [3 x nAngles] -- full mode only, fixed count
        %%% Green map: full mode only, single precision, at SawSpectrum.ETA -- assembled
        %%% from each angle's one calcGreen() solve. Only |G13_eff| is kept (the TGS
        %%% response, which the existence score reads); any other component is
        %%% recomputed on request (greenComponentMap())
        greenAbsG13EffMap      single = []   % [nVel x nAngles] |G13_eff| = |(G13_+-G13_-)/2|
                                              %   = |(G13_+-G31_+)/2| on greenVelAxis
        sawVelMap              cell   = {}   % {nAngles x 1} of variable-length double vectors
        sawModeNoMap           cell   = {}   % {nAngles x 1} of variable-length double vectors
        sawTypeMap             cell   = {}   % {nAngles x 1} of cellstr ('Rayleigh'/'pseudoRayleigh')
        significantModesMap    cell   = {}   % {nAngles x 1} of significantModes struct arrays,
                                              %   copied whole from cs{iAngle}.classification
    end

    properties (SetAccess = protected, Transient)
        %%% Angle-independent Ritz precompute, built ONCE from the unrotated structure
        %%% and shared by every cs{iAngle} -- the same fields CalcSingleHS would build
        %%% for itself, hoisted here so it happens once per sweep, not once per angle.
        %%% calc() hands it to every angle's Ritz solve explicitly: a thread-pool worker
        %%% receives each cs{iAngle} without its own (Transient) copy, and rebuilding it
        %%% there costs more than the Ritz solve. Transient: not written by save(),
        %%% empty after load() -- calc() rebuilds it then, once.
        precompStructForRitz struct = struct()
    end

    properties (SetAccess = protected)
        %%% 'low' reference of the base structure,
        %%% computed once in calc() and shared by every angle -- see GreenScore.reference
        greenScoreRef = []
        %%% ('ritzOnly' mode needs no shared reference: each angle scores against its own
        %%% floor, RitzScore.ownReference -- cs{i}.ritzScoreRef)

        %%% Refined classification: the GreenSpotScore reference of the base structure, set
        %%% when classify('refined') first solves the spot directions (each cs{i} holds its
        %%% candidates in refinedSaw). Cleared by calc().
        refinedSawRef        = []

        %%% The active classification of every angle (classify())
        activeMethod  char = ''   % 'green' | 'refined' | 'ritz'
        activeThreshold    = []   % 'high' | 'low' | numeric

        %%% Spot-averaged Green map -- calcSpotGreenMap() result (SpotGreen, on request):
        %%% |G| as a measurement with a finite pattern sees it (ProbeKernel). A picture
        %%% only; never read by any classification. Needs only the structure and the
        %%% angles, so calc() keeps it.
        spotGreenMap      single = []   % [nVel x nAngles] |X| on spotGreenVelAxis
        spotGreenVelAxis  double = []   % [nVel x 1] km/s
        spotGreenInfo     struct = struct([])   % SpotGreen.compute info: kernel, steps,
                                                % eta, quadrature checks per angle
    end

    properties
        %%% Classification method calc() runs on every angle: '' (default: follows each
        %%% run's calcMode -- 'green' on 'full', 'ritz' on 'ritzOnly'; the method in effect
        %%% is activeMethod) | 'green' | 'refined' | 'ritz' (see CalcSingleHS; on
        %%% 'ritzOnly' only 'ritz'). Set before calc(); classify() switches the active
        %%% classification without changing it.
        sawClassMethod char = ''
        %%% Its standard: 'high' | 'low' | numeric on the method's own scale (needs
        %%% sawClassMethod set explicitly); [] (default) = the method's default ('high'
        %%% for 'green' / 'refined', 'low' for 'ritz').
        sawThreshold = []
        %%% Green axis step [km/s]: the picture resolution; validated up to
        %%% GreenScore.MAX_STEP (16 m/s), larger steps are capped. Set before calc().
        greenVelStep(1,1) double = GreenScore.DEFAULT_STEP
    end

    methods
        function obj = CalcAngSweepHS(structure, angles, calcMode, retentionLevel, nVelGreen, useParallel)
        %%% Construct CalcAngSweepHS for an angular sweep
        %%% >> obj = CalcAngSweepHS(structure, angles, calcMode, retentionLevel, nVelGreen, useParallel)
        %%%
        %%% INPUT:  1. structure       .. Structure object, UNROTATED (each cs{iAngle}
        %%%                              rotates its own internal copy -- callers never
        %%%                              rotate by hand, same convention as CalcSingleHS)
        %%%  (opt)  2. angles          .. propagation angles [deg] (default = (0:2:180)',
        %%%                              91 angles; see the angles property for why not
        %%%                              a full circle)
        %%%  (opt)  3. calcMode        .. 'full' (default) | 'ritzOnly'
        %%%  (opt)  4. retentionLevel  .. 'full' (default) | 'saw_only' | 'none'
        %%%  (opt)  5. nVelGreen       .. number of points of the common Green axis (sets
        %%%                              the step, capped at GreenScore.MAX_STEP); default:
        %%%                              the greenVelStep property (6 m/s)
        %%%  (opt)  6. useParallel     .. "auto" (default) | "always" | "never"

            if nargin == 0
                return
            end

            validateattributes(structure, {'Structure'}, {'scalar'}, 'CalcAngSweepHS', 'structure');
            if nargin < 2 || isempty(angles); angles = (0:2:180)'; end
            validateattributes(angles, {'double'}, {'vector', 'real', 'finite'}, 'CalcAngSweepHS', 'angles');
            if nargin < 3; calcMode = []; end
            if nargin < 4; retentionLevel = []; end
            CalcCommon.validateConstructorFlags(calcMode, retentionLevel, 'CalcAngSweepHS:constructor', {'full', 'saw_only', 'none'});

            obj.structure = structure;      % UNROTATED -- kept as-is on this object
            obj.angles    = angles(:);

            if nargin >= 3 && ~isempty(calcMode);       obj.calcMode       = calcMode;       end
            if nargin >= 4 && ~isempty(retentionLevel); obj.retentionLevel = retentionLevel; end
            if nargin >= 5 && ~isempty(nVelGreen)
                validateattributes(nVelGreen, {'double'}, {'scalar', 'real', '>=', 2}, ...
                    'CalcAngSweepHS', 'nVelGreen');
                obj.nVelGreen = nVelGreen;
            end
            if nargin >= 6 && ~isempty(useParallel);    obj.useParallel    = useParallel;    end

            %%% Build the angle-independent Ritz precompute ONCE, from the UNROTATED
            %%% structure. Passing no partial struct means every field is built here.
            obj.precompStructForRitz = CalcCommon.ritzSAWPrecompute(structure);

            nAngles = numel(obj.angles);
            obj.cs  = cell(nAngles, 1);
            for iAngle = 1:nAngles
                obj.cs{iAngle} = Calc.singleHS(structure, obj.angles(iAngle), ...
                    obj.precompStructForRitz, obj.calcMode, obj.retentionLevel);
                %%% CalcSingleHS's constructor sees nullConstrMat/dKdC_cell/d2KdCdD_cell
                %%% already populated and skips rebuilding them -- this is what makes
                %%% the precompute genuinely shared.
            end
        end

        function obj = calc(obj, mode, retentionLevel)
        %%% Run the sweep: dispatch to the three-phase 'full' algorithm or a
        %%% per-angle 'ritzOnly' loop, then assemble the convenience maps.
        %%% >> obj = calc(obj, mode, retentionLevel)
        %%%
        %%% INPUT:  (opt) 1. mode           .. 'full' (default) | 'ritzOnly'
        %%%         (opt) 2. retentionLevel .. 'full' (default) | 'saw_only' | 'none'

            if nargin < 2 || isempty(mode);           mode           = obj.calcMode;       end
            if nargin < 3 || isempty(retentionLevel); retentionLevel = obj.retentionLevel;  end

            CalcCommon.validateCalcMode(mode, 'CalcAngSweepHS:calc');
            CalcCommon.validateRetentionLevel(retentionLevel, 'CalcAngSweepHS:calc', ...
                {'full', 'saw_only', 'none'});
            method = CalcCommon.resolveSawClassMethod(obj.sawClassMethod, obj.sawThreshold, ...
                mode, false, 'CalcAngSweepHS:calc');   % unset: 'green' / 'ritz' per mode, every run

            obj.calcMode       = mode;
            obj.retentionLevel = retentionLevel;
            obj = obj.clearRefinedSaw();   % a new calculation: refine again on request
            for ia = 1:numel(obj.cs)     % every angle carries the run's method and standard
                obj.cs{ia}.sawClassMethod = method;
                obj.cs{ia}.sawThreshold   = obj.sawThreshold;
            end

            %%% The ONE shared Ritz precompute, passed to every angle's Ritz solve (see
            %%% the property).
            precompLocal = obj.ritzPrecompute();

            if strcmp(mode, 'ritzOnly')
                %%% calcAll runs calcRitz -> calcLim -> modal-sum weights -> identifyModes
                %%% -> prune per angle (ritzOnlyStep). No classification input is shared
                %%% between angles: no velocity axis (no Green) and no floor reference
                %%% (each angle's own, RitzScore.ownReference), so a standalone angle
                %%% classifies identically.
                useP = CalcCommon.resolveParallel(obj.useParallel, numel(obj.cs));
                thrLocal = obj.sawThreshold;
                obj.cs = CalcCommon.forEachElement(obj.cs, useP, ...
                    @(cs, ~) CalcCommon.ritzOnlyStep(cs, retentionLevel, thrLocal, precompLocal));
            else
                obj = obj.calcFullPhases(precompLocal, method);   % classifies 'green'
                if ~strcmp(method, 'green')
                    obj.activeMethod = 'green';
                    obj = obj.classify(method, obj.sawThreshold);
                    return
                end
            end

            obj = obj.assembleMaps();
        end

        function precomp = ritzPrecompute(obj)
        %%% The shared Ritz precompute of the sweep, rebuilt first if missing (e.g. after load())
        %%% >> precomp = ritzPrecompute(obj)
        %%%
        %%% OUTPUT: precomp .. precompStructForRitz with every field populated
        %%%
        %%% Fill-if-missing (CalcCommon.ritzSAWPrecompute), so the constructor's copy is
        %%% returned at no cost; a rebuilt one is stored back on obj, so it is built at
        %%% most once per sweep.
            obj.precompStructForRitz = CalcCommon.ritzSAWPrecompute(obj.structure, obj.precompStructForRitz);
            precomp = obj.precompStructForRitz;
        end

        function obj = classify(obj, method, threshold)
        %%% Set the active classification of every angle (see CalcSingleHS.classify)
        %%% >> obj = classify(obj)
        %%% >> obj = classify(obj, method)
        %%% >> obj = classify(obj, method, threshold)
        %%%
        %%% INPUT: (opt) method    .. 'green' | 'refined' | 'ritz' (default: the active one)
        %%%        (opt) threshold .. 'high' | 'low' | numeric on the method's own scale
        %%%                           (default: the method's default standard)
        %%% OUTPUT: every cs{ia}.classification and the saw*/significantModes maps,
        %%%         refreshed; activeMethod, activeThreshold
        %%%
        %%% No Green/Ritz recomputation: every method's candidates are stored per angle.
        %%% 'refined' solves, on first use, every direction the angles' spot kernels need,
        %%% once (0:1:180 -> 360 directions on the half-degree grid, also covering the
        %%% reference; about twice the cost of the default 1-deg map), in parallel per
        %%% useParallel, and prints the cost; each angle then equals a standalone
        %%% CalcSingleHS of it, bit for bit. With retentionLevel 'saw_only', modes newly
        %%% accepted at a looser standard have no stored eigenvector. A stale result is
        %%% refused: recompute it.
            if nargin < 3; threshold = []; end
            if nargin < 2 || isempty(method); method = obj.activeMethod; end
            if isempty(method); method = obj.sawClassMethod; end
            if isempty(method)   % neither classified nor set: calcMode's own
                method = CalcCommon.resolveSawClassMethod('', threshold, obj.calcMode, ...
                    false, 'CalcAngSweepHS:classify');
            end
            method = CalcCommon.validateSawClassMethod(method, false, 'CalcAngSweepHS:classify');
            if isempty(threshold)
                threshold = CalcCommon.defaultThreshold(method, false);
            end
            CalcCommon.validateSawThreshold(threshold, 'CalcAngSweepHS:classify');
            CalcCommon.refuseStale(obj, 'CalcAngSweepHS:classify');
            if isempty(obj.cs) || isempty(obj.cs{1}.classification)
                error('CalcAngSweepHS:classify:notCalculated', 'calc() must run first.');
            end
            if ~strcmp(method, 'ritz') && ~strcmp(obj.calcMode, 'full')
                error('CalcAngSweepHS:classify:needsFull', ...
                    "Method '%s' needs calcMode='full'; this sweep is '%s'.", method, obj.calcMode);
            end
            if strcmp(method, 'refined') && ~all(cellfun(@(c) ~isempty(c.refinedSaw), obj.cs))
                t0 = tic;
                D = GreenSpotScore.directionSet(obj.angles);
                if isempty(GreenSpotScore.referenceCache('get', GreenSpotScore.referenceKey(obj.structure)))
                    D = [D; GreenSpotScore.directionSet(GreenSpotScore.REF_ANGLES)];   % solved once, shared
                end
                dc  = GreenSpotScore.solveDirections(obj.structure, D, obj.useParallel);
                ref = GreenSpotScore.reference(obj.structure, dc, obj.useParallel);
                obj.refinedSawRef = ref;
                for ia = 1:numel(obj.cs)     % milliseconds per angle
                    obj.cs{ia} = obj.cs{ia}.calcRefined(dc, ref);
                end
                fprintf('refined classification: %d spot directions solved in %.1f s\n', numel(D), toc(t0));
            end
            useP = CalcCommon.resolveParallel(obj.useParallel, numel(obj.cs));
            if strcmp(method, 'refined'); useP = false; end   % milliseconds per angle
            obj.cs = CalcCommon.forEachElement(obj.cs, useP, @(cs, ~) cs.classify(method, threshold));
            obj = obj.assembleMaps();
        end

        function obj = reclassify(obj, threshold)
        %%% Re-threshold the active classification of every angle (no recomputation)
        %%% >> obj = reclassify(obj, threshold)
        %%%
        %%% = classify(obj, activeMethod, threshold).
            obj = obj.classify(obj.activeMethod, threshold);
        end

        function tf = isStale(obj)
        %%% Whether any angle's classification carries another stamp than the current one
        %%% >> tf = isStale(obj)
            tf = CalcCommon.isStale(obj);
        end

        function T = refinedSawAgreement(obj)
        %%% Default vs refined entries of every angle, paired one to one (see CalcSingleHS)
        %%% >> T = refinedSawAgreement(obj)
        %%%
        %%% OUTPUT: T .. the per-angle tables of CalcSingleHS.refinedSawAgreement, stacked
            if isempty(obj.cs) || any(cellfun(@(c) isempty(c.refinedSaw), obj.cs))
                error('CalcAngSweepHS:refinedSawAgreement:notRefined', ...
                    'classify(''refined'') must run first.');
            end
            parts = cellfun(@(c) c.refinedSawAgreement(), obj.cs, 'UniformOutput', false);
            T = vertcat(parts{:});
        end

        function absGMap = greenComponentMap(obj, component)
        %%% Any Green component over the sweep on greenVelAxis, recomputed now (not stored)
        %%% >> absGMap = greenComponentMap(obj, 'G33')
        %%% >> absGMap = greenComponentMap(obj, 'G13')      % one-sided G13_+
        %%%
        %%% INPUT:  component .. 'GIJ' | 'GIJeff' (see CalcSingleHS.greenComponent);
        %%%                      'G13eff' returns the stored greenAbsG13EffMap
        %%% OUTPUT: absGMap   .. [nVel x nAngles] single
        %%%
        %%% A sweep stores only |G13_eff|; any other component costs one Green solve
        %%% per angle (about the sweep's original Green phase), on the same common
        %%% axis and damping. Requires calc('full').
            absGMap = CalcCommon.greenComponentMapOf(obj, component);
        end

        function obj = calcSpotGreenMap(obj, kernel, varargin)
        %%% |G| map as a finite-pattern measurement sees it: coherent sum over a ProbeKernel
        %%% Experimental: tried on Fe3Al(316) only, not validated against measurements; the
        %%% interface may change.
        %%% >> obj = calcSpotGreenMap(obj, kernel)
        %%% >> obj = calcSpotGreenMap(obj, kernel, 'DirStep', 0.1, 'RelKStep', 1e-3, 'Eta', 3e-3, 'VelStep', 0.004)
        %%%
        %%% INPUT:  1. kernel .. ProbeKernel (e.g. ProbeKernel.tgs(Lambda, spotRadii, probeRadii),
        %%%                      ProbeKernel.gaussian(sigmaRelK, sigmaAngDeg))
        %%%  (opt)  Name-Value: 'DirStep', 'RelKStep', 'Eta', 'VelStep' -- see
        %%%                      CalcSingleHS.calcSpotGreen
        %%% OUTPUT: obj with spotGreenMap, spotGreenVelAxis, spotGreenInfo
        %%%
        %%% Every displayed angle is the coherent sum of the complex Green's function over
        %%% the kernel's directions and |k| (SpotGreen.compute): a steep branch broadens and
        %%% drops, a flat one keeps its height. One Green solve per distinct direction,
        %%% shared between angles (with integer-degree angles and the default direction
        %%% steps, a 0:1:180 map needs 180 / DirStep solves: ~1 min serial at 0.125 deg,
        %%% the thread pool divides it). Only |X| is kept, same size as greenAbsG13EffMap
        %%% on the 4 m/s axis. Any calcMode; half-space only. Plot with
        %%% Vis.plotGreenMap(obj, [], [], [], [], 'spot'). A warning (SpotGreen:notConverged)
        %%% reports a quadrature step too coarse for the branch slopes (the "comb").
            if obj.structure.nLay ~= 1
                error('CalcAngSweepHS:calcSpotGreenMap:layered', ...
                    'calcSpotGreenMap is for a half-space (nLay = 1) only.');
            end
            opts = SpotGreen.options(kernel, varargin{:});
            [obj.spotGreenMap, obj.spotGreenVelAxis, obj.spotGreenInfo] = ...
                SpotGreen.compute(obj.structure, obj.angles, kernel, opts, obj.useParallel);
        end
    end

    methods (Access = private)
        function obj = clearRefinedSaw(obj)
        %%% Empty the sweep-level refined results (each cs{i} clears its own in calcRitz)
            obj.refinedSawRef = [];
        end

        function obj = calcFullPhases(obj, precomp, method)
        %%% Three-phase 'full' sweep: Ritz (all angles) -> Christoffel (all
        %%% angles) -> common Green axis preamble -> Green (all angles, one solve
        %%% each, identify + prune immediately per angle).
        %%% >> obj = calcFullPhases(obj, precomp, method)
        %%%
        %%% INPUT: precomp .. the shared Ritz precompute (calc() completes it)
        %%%        method  .. the run's resolved method (its standard reaches 'green' only
        %%%                   when it is 'green')
            %%% Resolved ONCE and reused by all three phase loops below.
            useP = CalcCommon.resolveParallel(obj.useParallel, numel(obj.cs));

            %%% Phase 1 -- Ritz, all angles. No identification/pruning yet: full-mode
            %%% identifyModes() needs Christoffel and Green too. The precompute goes in
            %%% explicitly: each worker's cs{iAngle} arrives without its own copy.
            obj.cs = CalcCommon.forEachElement(obj.cs, useP, @(cs, ~) cs.calcRitz(precomp));

            %%% Phase 2 -- Christoffel, all angles.
            obj.cs = CalcCommon.forEachElement(obj.cs, useP, @(cs, ~) cs.calcLim());

            %%% Common Green axis -- preamble to Phase 3: 0.3 km/s to 1.05 x the
            %%% largest longLimVel + 0.1 km/s, in ONE step for every angle (greenVelStep,
            %%% or the step giving nVelGreen points). Every angle's own axis
            %%% (CalcSingleHS.greenAxis) is a leading part of it, so each angle solves
            %%% once on the common axis (a genuine [nVel x nAngles] map) and is scored
            %%% exactly as a standalone angle with the same step is (GreenScore reads
            %%% the band).
            longVels = cellfun(@(c) c.longLimVel, obj.cs);
            vBound = SawSpectrum.AXIS_HI_X_VL * max(longVels);
            step = obj.greenVelStep;
            if isfinite(obj.nVelGreen)
                step = (vBound + SawSpectrum.EDGE_TAPER - SawSpectrum.VEL_MIN) / (obj.nVelGreen - 1);
            end
            step = GreenScore.resolveStep(step, 'CalcAngSweepHS:greenStepCapped');
            obj.greenVelAxis = GreenScore.gridAxis(vBound, step);
            for iAngle = 1:numel(obj.cs)     % handle objects: set in place
                obj.cs{iAngle}.greenVelStep   = step;
                obj.cs{iAngle}.nVelGreen      = NaN;
            end

            %%% Phase 3 -- Green, all angles; identify + prune per angle IMMEDIATELY,
            %%% not in a separate post-loop, which keeps peak RAM at one angle's
            %%% eigenvector matrix.
            %%%
            %%% calcGreen() gives one-sided G_+ (what the score reads) plus |G13_eff|
            %%% from the same solve: G_-(theta) is G_+(theta) transposed (reciprocity),
            %%% so the 0-180 grid needs no companion angle for the symmetrised map.
            %%% Score reference -- ONE per structure, shared by every angle, so a sweep
            %%% and a single angle classify identically (coarse Green-only sweep):
            %%% GreenScore's 'low' reference.
            obj.greenScoreRef = GreenScore.reference(obj.structure, [], obj.useParallel);
            sawRefLocal      = obj.greenScoreRef;
            greenVelAxisLocal   = obj.greenVelAxis;    % captured read-only by the step below
            retentionLevelLocal = obj.retentionLevel;  % captured read-only by the step below
            sawThrLocal         = [];                   % 'green' at its own default unless chosen
            if strcmp(method, 'green'); sawThrLocal = obj.sawThreshold; end
            obj.cs = CalcCommon.forEachElement(obj.cs, useP, ...
                @(cs, ~) CalcCommon.greenStep(cs, greenVelAxisLocal, retentionLevelLocal, ...
                sawRefLocal, sawThrLocal, precomp));
        end

        function obj = assembleMaps(obj)
        %%% Extract the sweep-level convenience maps from cs -- no new logic.
        %%% >> obj = assembleMaps(obj)
            %%% 'pseudoRayleigh' is the half-space-only half of this filter; the
            %%% layered sweep passes {'SAW'} to the same helper.
            [obj.sawVelMap, obj.sawModeNoMap, obj.sawTypeMap, obj.significantModesMap] = ...
                CalcCommon.sawMapsFrom(obj.cs, {'Rayleigh', 'pseudoRayleigh'});
            if ~isempty(obj.cs) && ~isempty(obj.cs{1}.activeMethod)
                obj.activeMethod    = obj.cs{1}.activeMethod;
                obj.activeThreshold = obj.cs{1}.activeThreshold;
            end

            if strcmp(obj.calcMode, 'full')
                %%% The Green fields are [] in 'ritzOnly' mode, and the cellfuns below
                %%% error on that -- hence the full-mode guard.
                obj.longVelMap = cellfun(@(c) c.longLimVel, obj.cs);
                obj.bulkVelMap = cell2mat(cellfun(@(c) c.bulkVels(:), obj.cs(:)', 'UniformOutput', false));
                obj.greenAbsG13EffMap = CalcCommon.greenG13EffMapFrom(obj.cs);
            end
        end
    end

    methods (Access = protected)
        function cp = copyElement(obj)
        %%% Deep-copy override (matlab.mixin.Copyable)
        %%% >> cp = copyElement(obj)
        %%% obj.cs holds handle-class CalcSingleHS objects -- the default
        %%% shallow copy@matlab.mixin.Copyable would copy the CELL ARRAY but
        %%% alias every CalcSingleHS handle inside it, so mutating cp.cs{i}
        %%% would also mutate obj.cs{i}. Deep-copy each element explicitly.
            cp = copyElement@matlab.mixin.Copyable(obj);
            cp.cs = cellfun(@copy, obj.cs, 'UniformOutput', false);
        end
    end
end
