classdef CalcAngSweepLayered < handle & matlab.mixin.Copyable
%%% CalcAngSweepLayered — angular sweep of CalcSingleLayered objects for the
%%% layered (film(s)-on-substrate) forward calculation at a FIXED layer stack;
%%% sweeping relDepth/wavelength is CalcWavevecAngSweep2Layered's job. Same
%%% relationship to CalcSingleLayered that CalcAngSweepHS has to CalcSingleHS:
%%% it exists so one identical Green velocity axis serves every angle. See
%%% docs/reference/architecture/forwardCalcs.md, "Angle-sweep orchestrators" and "The
%%% 0-180 degree angle grid, and Green's function on it".
%%%
%%% ONE STRUCTURAL DIFFERENCE from CalcAngSweepHS: Phase 2 is SUBSTRATE-ONLY
%%% Christoffel — no system-wide limiting-velocity scan exists for a layered
%%% stack — so the common Green axis ends at 1.05 x the fastest longitudinal
%%% limiting velocity of any layer (the GreenScoreLayered band end), or at the
%%% optional vMax (CalcSingleLayered.vMax) when that lies beyond it.
%%%
%%% GuliweRR v6.0

    properties (SetAccess = protected)
        structure       Structure   = Structure()      % UNROTATED base structure, nLay >= 2
                                                           % (enforced by CalcSingleLayered's own
                                                           % constructor) -- each cs{iAngle}
                                                           % rotates its own internal copy
        angles          double      = (0:2:180)'          % [nAngles x 1] propagation angles [deg].
                                                           % 0-180 default; the evenness argument
                                                           % behind it has not been independently
                                                           % re-verified for a layered assembly

        %%% calcFlags
        calcMode        char        = 'full'             % 'full' (default) | 'ritzOnly'
        retentionLevel  char        = 'full'             % 'full' (default) | 'saw_only' | 'none'
        useParallel(1,1) string {mustBeMember(useParallel, ["auto","always","never"])} = "auto"
                                                           % "auto" (default, parallel once nAngles
                                                           % >= ParallelPolicy.AUTO_THRESHOLD) |
                                                           % "always" | "never" -- applies uniformly
                                                           % to every phase loop; thread-pool based,
                                                           % not process-based
    end

    properties
        %%% Mirrors CalcSingleLayered.vMax: optional, only extends the Green axis
        %%% beyond the band end (NaN = the band end).
        vMax(1,1)       double      = NaN
    end

    properties (SetAccess = protected)
        %%% cs -- per-angle CalcSingleLayered array
        cs              cell        = {}                 % {nAngles x 1} of CalcSingleLayered objects

        %%% maps / christoffel -- substrate Christoffel per-angle maps, convenience
        %%% extractions from cs, filled by assembleMaps() at the end of calc()
        subLongVelMap          double = []   % [nAngles x 1] -- full mode only; SUBSTRATE-ONLY
                                              %   reference, not a system-wide bound
        subBulkVelMap          double = []   % [3 x nAngles] -- full mode only, substrate-only

        %%% maps / greenSolver + modeClassification -- common Green axis (built as the
        %%% Phase 3 preamble) plus per-angle Green/SAW maps. Green map: full mode only,
        %%% single precision, at SawSpectrum.ETA, from each angle's one solve; only
        %%% |G13_eff| is kept, other components via greenComponentMap()
        greenVelAxis    double      = []                 % [nVel x 1] commonGreenAxis(): the
                                                           % 4 m/s lattice spanning every angle's
                                                           % own axis (CalcSingleLayered.greenAxis)
        greenAbsG13EffMap       single = []   % [nVel x nAngles] |G13_eff| = |(G13_+-G13_-)/2|
                                              %   = |(G13_+-G31_+)/2| on greenVelAxis
        sawVelMap               cell  = {}   % {nAngles x 1} of variable-length double vectors --
                                              %   every 'SAW'-type entry (a single type, no
                                              %   Rayleigh/Sezawa split), contrasting
                                              %   CalcAngSweepHS's 'Rayleigh'/'pseudoRayleigh'
        sawModeNoMap            cell  = {}   % {nAngles x 1} of variable-length double vectors
        sawTypeMap               cell  = {}   % {nAngles x 1} of cellstr, every entry 'SAW'
        significantModesMap     cell  = {}   % {nAngles x 1} of significantModes struct arrays,
                                              %   copied whole from cs{iAngle}.classification --
                                              %   carries aboveSubstrateST per entry
    end

    properties (SetAccess = protected, Transient)
        %%% Angle-independent Ritz precompute, built ONCE from the unrotated structure
        %%% and shared by every cs{iAngle}: GammaMatrix.derivsByDepth contracts
        %%% calcParams(ilay).relDepth and mater(ilay).rho, never dirMat/rotation.
        %%% calc() hands it to every angle's Ritz solve explicitly (a thread-pool worker
        %%% receives each cs{iAngle} without its own Transient copy -- see
        %%% CalcAngSweepHS). Transient: not written by save(), empty after load() --
        %%% ritzPrecompute() rebuilds it then, once.
        precompStructForRitz struct = struct()
    end

    properties (SetAccess = protected)
        %%% No shared classification reference in either mode: 'full' scores each angle
        %%% against its own band median (GreenScoreLayered), 'ritzOnly' against its own
        %%% floor (RitzScore.ownReference -- cs{i}.ritzScoreRef)

        %%% The active classification of every angle (classify())
        activeMethod  char = ''   % 'green' | 'ritz'
        activeThreshold    = []   % 'high' | 'low' | numeric
    end

    properties
        %%% Classification method calc() runs on every angle: '' (default: follows each
        %%% run's calcMode -- 'green' on 'full', 'ritz' on 'ritzOnly'; the method in effect
        %%% is activeMethod) | 'green' | 'ritz' (see CalcSingleLayered; on 'ritzOnly' only
        %%% 'ritz'). Set before calc(); classify() switches the active classification
        %%% without changing it.
        sawClassMethod char = ''
        %%% Its standard: 'low' (sagittal-type waves) | 'high' (also Love-type waves) |
        %%% numeric on the method's own scale (needs sawClassMethod set explicitly);
        %%% [] (default) = 'low'.
        sawThreshold = []
    end

    methods
        function obj = CalcAngSweepLayered(structure, angles, vMax, calcMode, retentionLevel, unused6, rawByDepthCell, useParallel)
        %%% Construct CalcAngSweepLayered for an angular sweep
        %%% >> obj = CalcAngSweepLayered(structure, angles, vMax, calcMode, retentionLevel, [], rawByDepthCell, useParallel)
        %%%
        %%% INPUT:  1. structure       .. Structure object, UNROTATED (each cs{iAngle}
        %%%                              rotates its own internal copy -- callers never
        %%%                              rotate by hand, same convention as CalcSingleLayered);
        %%%                              must have nLay >= 2 (at least one film plus a substrate)
        %%%  (opt)  2. angles          .. propagation angles [deg] (default = (0:2:180)',
        %%%                              91 angles)
        %%%  (opt)  3. vMax            .. Green velocity-axis upper end (the axis runs to
        %%%                              max(vMax, 1.05 x the fastest longitudinal limiting
        %%%                              velocity of any layer)), propagated to every
        %%%                              cs{iAngle} at construction time ([] or NaN
        %%%                              (default) = the band end)
        %%%  (opt)  4. calcMode        .. 'full' (default) | 'ritzOnly'
        %%%  (opt)  5. retentionLevel  .. 'full' (default) | 'saw_only' | 'none'
        %%%  (opt)  6. (unused)        .. must be empty: the Green axis is each angle's
        %%%                              own (CalcSingleLayered.greenAxis)
        %%%  (opt)  7. rawByDepthCell  .. pre-built GammaMatrix.rawByDepth cell (one entry per
        %%%                              layer). Lets CalcWavevecAngSweep2Layered share the
        %%%                              expensive, relDepth-independent Legendre-integral
        %%%                              tensor across every kd value in its sweep instead of
        %%%                              rebuilding it inside each CalcAngSweepLayered.
        %%%                              Omitted, this class builds its own.
        %%%  (opt)  8. useParallel     .. "auto" (default) | "always" | "never"

            if nargin == 0
                return
            end

            validateattributes(structure, {'Structure'}, {'scalar'}, 'CalcAngSweepLayered', 'structure');
            if structure.nLay < 2
                error('CalcAngSweepLayered:constructor:needsFilmAndSubstrate', ...
                    'CalcAngSweepLayered needs at least one film plus a substrate (nLay >= 2); use CalcAngSweepHS for a homogeneous half-space.');
            end
            if nargin < 2 || isempty(angles); angles = (0:2:180)'; end
            validateattributes(angles, {'double'}, {'vector', 'real', 'finite'}, 'CalcAngSweepLayered', 'angles');
            if nargin < 4; calcMode = []; end
            if nargin < 5; retentionLevel = []; end
            CalcCommon.validateConstructorFlags(calcMode, retentionLevel, 'CalcAngSweepLayered:constructor', {'full', 'saw_only', 'none'});

            obj.structure = structure;      % UNROTATED -- kept as-is on this object
            obj.angles    = angles(:);

            if nargin >= 3 && ~isempty(vMax);           obj.vMax           = vMax;           end
            if nargin >= 4 && ~isempty(calcMode);       obj.calcMode       = calcMode;       end
            if nargin >= 5 && ~isempty(retentionLevel); obj.retentionLevel = retentionLevel; end
            if nargin >= 6 && ~isempty(unused6)
                error('CalcAngSweepLayered:constructor:unusedArgument', ...
                    'The 6th argument is unused (the Green axis is each angle''s own); pass [].');
            end
            if nargin >= 8 && ~isempty(useParallel);    obj.useParallel    = useParallel;    end

            %%% Build the angle-independent Ritz precompute ONCE, from the unrotated
            %%% structure. When rawByDepthCell is omitted, derivsByDepth builds its own,
            %%% so a bare CalcAngSweepLayered(structure, angles) stays standalone.
            if nargin < 7; rawByDepthCell = {}; end
            obj.precompStructForRitz = CalcCommon.ritzSAWPrecompute(structure, struct(), rawByDepthCell);

            nAngles = numel(obj.angles);
            obj.cs  = cell(nAngles, 1);
            for iAngle = 1:nAngles
                obj.cs{iAngle} = Calc.singleLayered(structure, obj.angles(iAngle), obj.vMax, ...
                    obj.precompStructForRitz, obj.calcMode, obj.retentionLevel);
                %%% CalcSingleLayered's constructor sees nullConstrMat/dKdC_cell/
                %%% d2KdCdD_cell already populated and skips rebuilding them.
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

            CalcCommon.validateCalcMode(mode, 'CalcAngSweepLayered:calc');
            CalcCommon.validateRetentionLevel(retentionLevel, 'CalcAngSweepLayered:calc', ...
                {'full', 'saw_only', 'none'});
            method = CalcCommon.resolveSawClassMethod(obj.sawClassMethod, obj.sawThreshold, ...
                mode, true, 'CalcAngSweepLayered:calc');   % unset: 'green' / 'ritz' per mode, every run

            obj.calcMode       = mode;
            obj.retentionLevel = retentionLevel;
            for ia = 1:numel(obj.cs)     % every angle carries the run's method and standard
                obj.cs{ia}.sawClassMethod = method;
                obj.cs{ia}.sawThreshold   = obj.sawThreshold;
            end

            if strcmp(mode, 'ritzOnly')
                %%% Per angle: calcRitz -> substrate calcLim -> modal-sum weights ->
                %%% identifyModes -> prune (ritzOnlyStep); nothing is shared between angles
                %%% (each scores against its own floor, see CalcAngSweepHS.calc).
                useP = CalcCommon.resolveParallel(obj.useParallel, numel(obj.cs));
                thrLocal  = obj.sawThreshold;
                precompLocal = obj.ritzPrecompute();   % shared, passed to every angle
                obj.cs = CalcCommon.forEachElement(obj.cs, useP, ...
                    @(cs, ~) CalcCommon.ritzOnlyStep(cs, retentionLevel, thrLocal, precompLocal));
            else
                obj = obj.calcFullPhases(method);      % classifies 'green'
                if ~strcmp(method, 'green')
                    obj.activeMethod = 'green';
                    obj = obj.classify(method, obj.sawThreshold);
                    return
                end
            end

            obj = obj.assembleMaps();
        end

        function obj = setElement(obj, iAngle, csElem)
        %%% Write a single per-angle result back into cs -- cs is SetAccess=protected,
        %%% so this is the narrow accessor CalcWavevecAngSweep2Layered's flat-dispatch
        %%% path uses to scatter a flat-parfor result back in.
        %%% >> obj = setElement(obj, iAngle, csElem)
            obj.cs{iAngle} = csElem;
        end

        function obj = setGreenVelAxis(obj, axis)
        %%% Write the common Green velocity axis. Normally built inside
        %%% calcFullPhases(), but the flat-dispatch path computes it per-wavelength
        %%% itself and needs to write it back in.
        %%% >> obj = setGreenVelAxis(obj, axis)
            obj.greenVelAxis = axis;
        end

        function velAxis = commonGreenAxis(obj)
        %%% The sweep's common Green axis, from every angle's calcLim() result
        %%% >> velAxis = commonGreenAxis(obj)
        %%%
        %%% GreenScoreLayered.gridAxis from the lowest band start (0.3 x bandVSmin) to the
        %%% largest bound max(vMax, 1.05 x bandVLmax) over the angles: the 4 m/s lattice
        %%% spanning every angle's own axis (CalcSingleLayered.greenAxis), so each angle's
        %%% band is bit for bit the one it scores on alone. Read from the cs{} objects
        %%% themselves (calcLim() results). Public so that the flat dispatch of
        %%% CalcWavevecAngSweep2Layered builds the same axis.
            lo = cellfun(@(c) GreenScoreLayered.BAND_LO_X_VSMIN * c.bandVSmin, obj.cs);
            hi = cellfun(@(c) max(c.vMax, GreenScoreLayered.BAND_HI_X_VLMAX * c.bandVLmax), obj.cs);
            velAxis = GreenScoreLayered.gridAxis(min(lo), max(hi));
        end

        function obj = finalizeMaps(obj)
        %%% Public wrapper around the private assembleMaps() -- lets
        %%% CalcWavevecAngSweep2Layered's flat-dispatch path trigger map assembly after
        %%% scattering results back via setElement.
        %%% >> obj = finalizeMaps(obj)
            obj = obj.assembleMaps();
        end

        function precomp = ritzPrecompute(obj)
        %%% The shared Ritz precompute of the sweep, rebuilt first if missing (e.g. after load())
        %%% >> precomp = ritzPrecompute(obj)
        %%%
        %%% OUTPUT: precomp .. precompStructForRitz with every field populated
        %%%
        %%% Same contract as CalcAngSweepHS.ritzPrecompute. A rebuild after load() has
        %%% no rawByDepthCell to reuse, so derivsByDepth builds its own (once).
            obj.precompStructForRitz = CalcCommon.ritzSAWPrecompute(obj.structure, obj.precompStructForRitz);
            precomp = obj.precompStructForRitz;
        end

        function obj = classify(obj, method, threshold)
        %%% Set the active classification of every angle (see CalcSingleLayered.classify)
        %%% >> obj = classify(obj)
        %%% >> obj = classify(obj, method)
        %%% >> obj = classify(obj, method, threshold)
        %%%
        %%% INPUT: (opt) method    .. 'green' | 'ritz' (default: the active one)
        %%%        (opt) threshold .. 'low' | 'high' | numeric on the method's own scale
        %%%                           (default 'low')
        %%% OUTPUT: every cs{ia}.classification and the saw*/significantModes maps,
        %%%         refreshed; activeMethod, activeThreshold. No Green/Ritz recomputation:
        %%%         both methods' candidates are stored per angle. With retentionLevel
        %%%         'saw_only', modes newly accepted at a looser standard have no stored
        %%%         eigenvector. A stale result is refused: recompute it.
            if nargin < 3; threshold = []; end
            if nargin < 2 || isempty(method); method = obj.activeMethod; end
            if isempty(method); method = obj.sawClassMethod; end
            if isempty(method)   % neither classified nor set: calcMode's own
                method = CalcCommon.resolveSawClassMethod('', threshold, obj.calcMode, ...
                    true, 'CalcAngSweepLayered:classify');
            end
            method = CalcCommon.validateSawClassMethod(method, true, 'CalcAngSweepLayered:classify');
            if isempty(threshold)
                threshold = CalcCommon.defaultThreshold(method, true);
            end
            CalcCommon.validateSawThreshold(threshold, 'CalcAngSweepLayered:classify');
            CalcCommon.refuseStale(obj, 'CalcAngSweepLayered:classify');
            if isempty(obj.cs) || isempty(obj.cs{1}.classification)
                error('CalcAngSweepLayered:classify:notCalculated', 'calc() must run first.');
            end
            if strcmp(method, 'green') && ~strcmp(obj.calcMode, 'full')
                error('CalcAngSweepLayered:classify:needsFull', ...
                    "Method 'green' needs calcMode='full'; this sweep is '%s'.", obj.calcMode);
            end
            useP = CalcCommon.resolveParallel(obj.useParallel, numel(obj.cs));
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
    end

    methods (Access = private)
        function obj = calcFullPhases(obj, method)
        %%% Three-phase 'full' sweep: Ritz (all angles) -> substrate-only
        %%% Christoffel (all angles) -> common Green axis preamble -> Green (all
        %%% angles, one solve each, identify + prune immediately per angle).
        %%%
        %%% Every phase body is either a one-line method call on the per-angle object
        %%% or CalcCommon.greenStep -- exactly what CalcWavevecAngSweep2Layered's
        %%% flat-dispatch path runs on its flattened (wavelength, angle) grid, so the
        %%% nested and flat paths cannot silently diverge.
        %%% >> obj = calcFullPhases(obj, method)      % method: the run's resolved one
            %%% Resolved ONCE and reused by all three phase loops.
            useP = CalcCommon.resolveParallel(obj.useParallel, numel(obj.cs));

            %%% obj.vMax is public and may have been changed after construction:
            %%% hand the current value to every angle (cs{} are handles).
            for ia = 1:numel(obj.cs)
                obj.cs{ia}.vMax = obj.vMax;
            end

            %%% Phase 1 -- Ritz, all angles. No identification/pruning yet: full-mode
            %%% identifyModes() needs Christoffel and Green too. The shared precompute
            %%% goes in explicitly: each worker's cs{iAngle} arrives without its own copy.
            precomp = obj.ritzPrecompute();
            obj.cs = CalcCommon.forEachElement(obj.cs, useP, @(cs, ~) cs.calcRitz(precomp));

            %%% Phase 2 -- SUBSTRATE-ONLY Christoffel, all angles (each angle's
            %%% aboveSubstrateST reference), plus the all-layer band limits.
            obj.cs = CalcCommon.forEachElement(obj.cs, useP, @(cs, ~) cs.calcLim());

            %%% Common Green axis -- preamble to Phase 3: spans every angle's own axis,
            %%% up to max(vMax, 1.05 x the largest all-layer longitudinal limit). The
            %%% optional vMax only extends the axis. Each angle scores on its own band of it.
            obj.greenVelAxis = obj.commonGreenAxis();

            %%% Phase 3 -- Green, all angles; identify + prune per angle IMMEDIATELY,
            %%% not in a separate post-loop, keeping peak RAM at one angle's eigenvectors.
            %%% No reference: every angle is normalised by its own band (GreenScoreLayered),
            %%% so a sweep and a single angle classify identically.
            greenVelAxisLocal   = obj.greenVelAxis;    % captured read-only by the step below
            retentionLevelLocal = obj.retentionLevel;  % captured read-only by the step below
            sawThrLocal         = [];                   % 'green' at its own default unless chosen
            if strcmp(method, 'green'); sawThrLocal = obj.sawThreshold; end
            obj.cs = CalcCommon.forEachElement(obj.cs, useP, ...
                @(cs, ~) CalcCommon.greenStep(cs, greenVelAxisLocal, retentionLevelLocal, ...
                [], sawThrLocal, precomp));
        end

        function obj = assembleMaps(obj)
        %%% Extract the sweep-level convenience maps from cs (type 'SAW' on both the
        %%% Green and the Ritz-only path)
        %%% >> obj = assembleMaps(obj)
            [obj.sawVelMap, obj.sawModeNoMap, obj.sawTypeMap, obj.significantModesMap] = ...
                CalcCommon.sawMapsFrom(obj.cs, {'SAW'});
            if ~isempty(obj.cs) && ~isempty(obj.cs{1}.activeMethod)
                obj.activeMethod    = obj.cs{1}.activeMethod;
                obj.activeThreshold = obj.cs{1}.activeThreshold;
            end

            if strcmp(obj.calcMode, 'full')
                %%% The Green fields are [] in 'ritzOnly' mode, and the cellfuns below
                %%% error on that -- hence the guard.
                obj.subLongVelMap = cellfun(@(c) c.subLongLimVel, obj.cs);
                obj.subBulkVelMap = cell2mat(cellfun(@(c) c.subBulkVels(:), obj.cs(:)', 'UniformOutput', false));
                obj.greenAbsG13EffMap = CalcCommon.greenG13EffMapFrom(obj.cs);
            end
        end
    end

    methods (Access = protected)
        function cp = copyElement(obj)
        %%% Deep-copy override (matlab.mixin.Copyable)
        %%% >> cp = copyElement(obj)
        %%% obj.cs holds handle-class CalcSingleLayered objects -- the default shallow
        %%% copy@matlab.mixin.Copyable would copy the CELL ARRAY but alias every
        %%% CalcSingleLayered handle inside it, so mutating cp.cs{i} would also mutate
        %%% obj.cs{i}. Deep-copy each element explicitly.
            cp = copyElement@matlab.mixin.Copyable(obj);
            cp.cs = cellfun(@copy, obj.cs, 'UniformOutput', false);
        end
    end
end
