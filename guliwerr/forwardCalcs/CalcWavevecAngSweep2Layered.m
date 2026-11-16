classdef CalcWavevecAngSweep2Layered < handle & matlab.mixin.Copyable
%%% CalcWavevecAngSweep2Layered — 2D sweep (propagation angle AND real
%%% wavelength) for the layered forward calculation. Built on
%%% CalcAngSweepLayered: constructs one per wavelength (via
%%% structure.ApplyWavelength) and drives it through its own calc(). See
%%% docs/reference/architecture/forwardCalcs.md, "2D sweep orchestrators", "Real-unit
%%% awareness lives here, not below" and "`retentionLevel='none'` default, and
%%% no top-level aggregation".
%%%
%%% WARNING: legOrder must stay constant across the whole wavelength sweep --
%%% rawByDepthCell is built once and shared across every wavelength; varying
%%% legOrder silently invalidates that cache (nothing enforces this
%%% automatically).
%%%
%%% WARNING FOR CALLERS: the substrate (layer 1) is not swept -- choose its
%%% relDepth conservatively deep relative to the SMALLEST value in
%%% `wavelengths`; not enforced programmatically (no way to know a priori
%%% how deep is "enough" without a convergence study).
%%%
%%% GuliweRR v6.0

    properties (SetAccess = protected)
        %%% Inputs
        structure       Structure   = Structure()      % REFERENCE structure, UNROTATED. Every FILM
                                                           % layer's CalcParams.realThickness must
                                                           % already be set (real units, e.g.
                                                           % microns); the constructor validates
                                                           % this up front. The substrate (layer 1)
                                                           % realThickness is NOT read -- its
                                                           % relDepth/legOrder stay as the reference
                                                           % structure has them, unaffected by the
                                                           % wavelength sweep
        angles          double      = (0:2:180)'          % [nAngles x 1] propagation angles [deg]
        wavelengths     double      = []                 % [nWav x 1] positive real wavelengths,
                                                           % same units as every film layer's
                                                           % realThickness -- REQUIRED, no default
        vMax(1,1)       double      = NaN                % optional Green velocity-axis upper
                                                           % end, propagated to every sweeps{iWav}
                                                           % at calc() time; NaN = the band end
        calcMode        char        = 'full'             % 'full' (default) | 'ritzOnly'
        retentionLevel  char        = 'none'              % 'full' | 'saw_only' | 'none' (default) --
                                                           % differs from the 1D classes' 'full'
                                                           % default because a 2D sweep multiplies
                                                           % eigenvector RAM across the outer loop
        useParallel(1,1) string {mustBeMember(useParallel, ["auto","always","never"])} = "auto"
                                                           % "auto" (default) | "always" | "never" --
                                                           % forwarded UNCHANGED to every
                                                           % sweeps{iWav}, which makes its own
                                                           % parallel/serial decision at the angle
                                                           % level; this class's outer iWav loop
                                                           % stays serial, MATLAB forbidding nested
                                                           % parfor

        %%% Per-wavelength containers -- one CalcAngSweepLayered per wavelength
        sweeps          cell        = {}                 % {nWav x 1} of CalcAngSweepLayered objects
    end

    properties
        %%% flattenSweep -- settable any time before calc(). false (default): nested
        %%% behaviour, outer iWav loop serial, each sweeps{iWav} parallelising its own
        %%% angle loop once nAngles crosses ParallelPolicy.AUTO_THRESHOLD. true:
        %%% dispatch across the flattened (wavelength, angle) grid instead, decided
        %%% from nWav*nAngles rather than nAngles alone -- useful when neither axis
        %%% alone is large enough but their product is.
        %%%
        %%% WARNING: peak RAM is higher than the nested path, since ALL nWav*nAngles
        %%% elements are in flight through Phases 1 and 2 together rather than one
        %%% wavelength's worth at a time. A nested 2D sweep already multiplies
        %%% eigenvector RAM across the outer loop (why this class defaults
        %%% retentionLevel to 'none'); flattenSweep=true multiplies that peak by nWav
        %%% again in 'full' mode with retentionLevel='full', which is what
        %%% DispersionBranches.track requires.
        flattenSweep(1,1) logical = false

        %%% Classification method ('' default: follows each run's calcMode, 'green' on
        %%% 'full', 'ritz' on 'ritzOnly' | 'green' | 'ritz'; on 'ritzOnly' only 'ritz') and
        %%% standard ([] = 'low'; a numeric one needs the method set), forwarded to every
        %%% sweeps{iWav} (see CalcAngSweepLayered); classify() switches the active
        %%% classification without changing them. No reference is computed: every grid
        %%% point is scored against its own band / floor.
        sawClassMethod char = ''
        sawThreshold = []
    end

    properties (SetAccess = protected, Transient)
        %%% Built ONCE for the whole 2D sweep, in the constructor, before any
        %%% wavelength-dependent work. Shared across every wavelength AND every angle
        %%% -- rawByDepth depends on neither. Transient: not written by save(),
        %%% empty after load(); calc() rebuilds it if a loaded object is re-run.
        rawByDepthCell  cell        = {}
    end

    properties (SetAccess = protected)
        %%% Convenience 2D extraction -- filled by assembleMaps() at the end of calc(),
        %%% re-indexing each sweeps{iWav}'s already-computed 1D maps
        sawVelMap2D              cell = {}   % {nWav x nAngles} of variable-length double vectors
        sawModeNoMap2D           cell = {}   % {nWav x nAngles}
        sawTypeMap2D             cell = {}   % {nWav x nAngles} of cellstr
        significantModesMap2D    cell = {}   % {nWav x nAngles} of significantModes struct arrays
    end

    methods
        function obj = CalcWavevecAngSweep2Layered(structure, angles, wavelengths, vMax, calcMode, retentionLevel, unused7, useParallel)
        %%% Construct CalcWavevecAngSweep2Layered for a 2D (angle x wavelength) sweep
        %%% >> obj = CalcWavevecAngSweep2Layered(structure, angles, wavelengths, vMax, calcMode, retentionLevel, [], useParallel)
        %%%
        %%% INPUT:  1. structure       .. Structure object, UNROTATED, nLay >= 2 (at least one
        %%%                              film plus a substrate). Every FILM layer's
        %%%                              CalcParams.realThickness must already be set --
        %%%                              validated here, before any work starts (fail fast).
        %%%  (opt)  2. angles          .. propagation angles [deg] (default = (0:2:180)')
        %%%         3. wavelengths     .. [nWav x 1] positive real wavelengths, same units as
        %%%                              realThickness -- REQUIRED, no default
        %%%  (opt)  4. vMax            .. Green velocity-axis upper end, only extends the
        %%%                              axis beyond 1.05 x the fastest longitudinal
        %%%                              limit of any layer ([] or NaN (default) = that
        %%%                              band end)
        %%%  (opt)  5. calcMode        .. 'full' (default) | 'ritzOnly'
        %%%  (opt)  6. retentionLevel  .. 'full' | 'saw_only' | 'none' (default = 'none',
        %%%                              differing from CalcAngSweepLayered's 'full')
        %%%  (opt)  7. (unused)        .. must be empty: the Green axis is each angle's
        %%%                              own (CalcSingleLayered.greenAxis)
        %%%  (opt)  8. useParallel     .. "auto" (default) | "always" | "never" -- forwarded
        %%%                              unchanged to every sweeps{iWav}

            if nargin == 0
                return
            end

            validateattributes(structure, {'Structure'}, {'scalar'}, 'CalcWavevecAngSweep2Layered', 'structure');
            if structure.nLay < 2
                error('CalcWavevecAngSweep2Layered:constructor:needsFilmAndSubstrate', ...
                    'CalcWavevecAngSweep2Layered needs at least one film plus a substrate (nLay >= 2); use CalcAngSweepHS for a homogeneous half-space.');
            end

            if nargin < 2 || isempty(angles); angles = (0:2:180)'; end
            validateattributes(angles, {'double'}, {'vector', 'real', 'finite'}, ...
                'CalcWavevecAngSweep2Layered', 'angles');
            if nargin < 5; calcMode = []; end
            if nargin < 6; retentionLevel = []; end
            CalcCommon.validateConstructorFlags(calcMode, retentionLevel, ...
                'CalcWavevecAngSweep2Layered:constructor', {'full', 'saw_only', 'none'});

            %%% wavelengths is REQUIRED -- caught here so the caller gets a
            %%% project-standard error ID, not a bare "not enough input arguments".
            if nargin < 3 || isempty(wavelengths)
                error('CalcWavevecAngSweep2Layered:constructor:missingWavelengths', ...
                    'wavelengths is required (no default) -- pass a vector of positive real wavelengths.');
            end
            validateattributes(wavelengths, {'double'}, {'vector', 'real', 'finite', 'positive'}, ...
                'CalcWavevecAngSweep2Layered', 'wavelengths');

            %%% Fail fast: every FILM layer needs realThickness set BEFORE any work
            %%% starts, rather than surfacing partway through calc() via
            %%% ApplyWavelength's own lazy error at whichever wavelength runs first.
            for ilay = 2:structure.nLay
                if isnan(structure.calcParams(ilay).realThickness)
                    error('CalcWavevecAngSweep2Layered:constructor:missingRealThickness', ...
                        'Layer %d has no realThickness set -- cannot sweep wavelength.', ilay);
                end
            end

            obj.structure    = structure;      % UNROTATED -- kept as-is on this object
            obj.angles       = angles(:);
            obj.wavelengths  = wavelengths(:);

            if nargin >= 4 && ~isempty(vMax);           obj.vMax           = vMax;           end
            if nargin >= 5 && ~isempty(calcMode);       obj.calcMode       = calcMode;       end
            if nargin >= 6 && ~isempty(retentionLevel); obj.retentionLevel = retentionLevel; end
            if nargin >= 7 && ~isempty(unused7)
                error('CalcWavevecAngSweep2Layered:constructor:unusedArgument', ...
                    'The 7th argument is unused (the Green axis is each angle''s own); pass [].');
            end
            if nargin >= 8 && ~isempty(useParallel);    obj.useParallel    = useParallel;    end

            %%% Build rawByDepthCell ONCE here for the whole 2D sweep -- the same work
            %%% GammaMatrix.derivsByDepth would do internally, lifted out so it runs
            %%% nLay times total rather than once per wavelength.
            obj = obj.buildRawByDepthCell();
        end

        function obj = calc(obj, mode, retentionLevel)
        %%% Run the 2D sweep: one CalcAngSweepLayered per wavelength (sharing
        %%% rawByDepthCell), then assemble the 2D convenience maps.
        %%% >> obj = calc(obj, mode, retentionLevel)
        %%%
        %%% INPUT:  (opt) 1. mode           .. 'full' (default) | 'ritzOnly'
        %%%         (opt) 2. retentionLevel .. 'full' | 'saw_only' | 'none'
        %%%                                    (default = obj.retentionLevel, 'none'
        %%%                                    unless overridden at construction)
        %%%
        %%% Each iWav iteration runs CalcAngSweepLayered's own three-phase loop
        %%% (Ritz -> substrate Christoffel -> Green) across every angle; this class
        %%% adds one outer loop over wavelength and nothing new at the physics level.
        %%% It never constructs RitzSolverSAW/GreenSolver/CalcSingleLayered directly.
        %%% The outer iWav loop stays serial (MATLAB forbids nested parfor);
        %%% obj.useParallel is forwarded to each sweeps{iWav} instead.

            if nargin < 2 || isempty(mode);           mode           = obj.calcMode;       end
            if nargin < 3 || isempty(retentionLevel); retentionLevel = obj.retentionLevel;  end

            CalcCommon.validateCalcMode(mode, 'CalcWavevecAngSweep2Layered:calc');
            CalcCommon.validateRetentionLevel(retentionLevel, 'CalcWavevecAngSweep2Layered:calc', ...
                {'full', 'saw_only', 'none'});
            method = CalcCommon.resolveSawClassMethod(obj.sawClassMethod, obj.sawThreshold, ...
                mode, true, 'CalcWavevecAngSweep2Layered:calc');   % unset: 'green' / 'ritz' per mode, every run

            obj.calcMode       = mode;
            obj.retentionLevel = retentionLevel;

            if isempty(obj.rawByDepthCell)   % Transient -- empty on a load()'d object
                obj = obj.buildRawByDepthCell();
            end

            if obj.flattenSweep
                obj = obj.calcFlat(mode, retentionLevel, method);
                return
            end

            nWav = numel(obj.wavelengths);
            obj.sweeps = cell(nWav, 1);
            for iWav = 1:nWav
                structAtWav = obj.structure.ApplyWavelength(obj.wavelengths(iWav));
                obj.sweeps{iWav} = Calc.angSweepLayered(structAtWav, obj.angles, obj.vMax, ...
                    mode, retentionLevel, [], obj.rawByDepthCell, obj.useParallel);
                %%% rawByDepthCell (7th arg) passes straight through to
                %%% GammaMatrix.derivsByDepth, skipping the expensive tensor rebuild.
                %%% Only the cheap per-relDepth rescaling (dKdC_cell/d2KdCdD_cell) is
                %%% redone at each iWav.
                obj.sweeps{iWav}.sawClassMethod = obj.sawClassMethod;
                obj.sweeps{iWav}.sawThreshold   = obj.sawThreshold;
                obj.sweeps{iWav} = obj.sweeps{iWav}.calc(mode, retentionLevel);
            end

            obj = obj.assembleMaps();
        end

        function obj = classify(obj, method, threshold)
        %%% Set the active classification of every wavelength's sweep (CalcAngSweepLayered.classify)
        %%% >> obj = classify(obj, method)
        %%% >> obj = classify(obj, method, threshold)
        %%%
        %%% INPUT: (opt) method    .. 'green' | 'ritz' (default: the active one)
        %%%        (opt) threshold .. 'low' | 'high' | numeric (default 'low')
        %%% OUTPUT: every sweeps{iWav} reclassified, the 2D maps reassembled
            if nargin < 2; method = []; end
            if nargin < 3; threshold = []; end
            if isempty(obj.sweeps)
                error('CalcWavevecAngSweep2Layered:classify:notCalculated', 'calc() must run first.');
            end
            for iWav = 1:numel(obj.sweeps)
                obj.sweeps{iWav} = obj.sweeps{iWav}.classify(method, threshold);
            end
            obj = obj.assembleMaps();
        end

        function obj = reclassify(obj, threshold)
        %%% Re-threshold the active classification of every wavelength (no recomputation)
        %%% >> obj = reclassify(obj, threshold)
            obj = obj.classify([], threshold);
        end
    end

    methods (Access = private)
        function obj = buildRawByDepthCell(obj)
        %%% Build the shared GammaMatrix.rawByDepth cell, one entry per layer
        %%% >> obj = buildRawByDepthCell(obj)
            obj.rawByDepthCell = cell(obj.structure.nLay, 1);
            for ilay = 1:obj.structure.nLay
                obj.rawByDepthCell{ilay} = GammaMatrix.rawByDepth(obj.structure.calcParams(ilay));
            end
        end

        function obj = calcFlat(obj, mode, retentionLevel, method)
        %%% Flat-dispatch variant of calc(), used when obj.flattenSweep is true
        %%% >> obj = calcFlat(obj, mode, retentionLevel, method)   % method: the run's resolved one
        %%%
        %%% Constructs every sweeps{iWav} serially (object construction only), then
        %%% dispatches phase-by-phase across the FLATTENED (wavelength, angle) grid, one
        %%% parfor per phase. It runs the same per-point steps
        %%% CalcAngSweepLayered.calcFullPhases() runs, so the two paths cannot diverge.
        %%%
        %%% 'full' mode has a genuine per-wavelength BARRIER: Phase 3's common
        %%% greenVelAxis is built from Phase 2's result AT THAT WAVELENGTH. So Phases 1
        %%% and 2 flatten fully across (wavelength, angle), the barrier is a cheap serial
        %%% per-wavelength reduction, and Phase 3 flattens again with each element
        %%% reading its own wavelength's axis. 'ritzOnly' has no barrier and flattens in
        %%% a single pass.

            nWav    = numel(obj.wavelengths);
            nAngles = numel(obj.angles);

            obj.sweeps = cell(nWav, 1);
            for iWav = 1:nWav
                structAtWav = obj.structure.ApplyWavelength(obj.wavelengths(iWav));
                obj.sweeps{iWav} = Calc.angSweepLayered(structAtWav, obj.angles, obj.vMax, ...
                    mode, retentionLevel, [], obj.rawByDepthCell, "never");
                %%% useParallel = "never" -- the flat dispatch below owns parallelism;
                %%% sweeps{iWav}.calc() is never called in this path.
                obj.sweeps{iWav}.sawClassMethod = obj.sawClassMethod;
                obj.sweeps{iWav}.sawThreshold   = obj.sawThreshold;
                for ia = 1:nAngles       % every grid point carries the run's method and standard
                    c = obj.sweeps{iWav}.cs{ia};          % a handle: set in place
                    c.sawClassMethod = method;
                    c.sawThreshold   = obj.sawThreshold;
                end
            end

            [csFlat, wavIdx, angIdx] = CalcCommon.flattenSweeps(obj.sweeps, nAngles);

            %%% Decided from nWav*nAngles, not nAngles alone -- the point of the flat path.
            useP = CalcCommon.resolveParallel(obj.useParallel, numel(csFlat));

            %%% Each wavelength's shared Ritz precompute (relDepth changes with the
            %%% wavelength), handed to every grid point's Ritz solve by flat index: a
            %%% worker's csFlat{ell} arrives without its own Transient copy.
            precompPerWav = cellfun(@(sw) sw.ritzPrecompute(), obj.sweeps, 'UniformOutput', false);

            if strcmp(mode, 'ritzOnly')
                %%% Every grid point scores against its own floor (RitzScore.ownReference):
                %%% no per-wavelength reference is needed.
                thrLocal  = obj.sawThreshold;
                retentionLocal = retentionLevel;
                csFlat = CalcCommon.forEachElement(csFlat, useP, ...
                    @(cs, ell) CalcCommon.ritzOnlyStep(cs, retentionLocal, thrLocal, ...
                    precompPerWav{wavIdx(ell)}));
            else
                %%% Phase 1 -- Ritz, flattened across the whole grid.
                csFlat = CalcCommon.forEachElement(csFlat, useP, ...
                    @(cs, ell) cs.calcRitz(precompPerWav{wavIdx(ell)}));

                %%% Phase 2 -- substrate Christoffel, flattened across the whole grid.
                csFlat = CalcCommon.forEachElement(csFlat, useP, @(cs, ~) cs.calcLim());

                %%% Barrier -- scatter back so each wavelength's common greenVelAxis is
                %%% built from THAT wavelength's own angles only, exactly as
                %%% CalcAngSweepLayered.calcFullPhases builds it (commonGreenAxis).
                obj.sweeps = CalcCommon.scatterToSweeps(obj.sweeps, csFlat, wavIdx, angIdx);
                greenAxisPerWav = cell(nWav, 1);
                for iWav = 1:nWav
                    axisLocal = obj.sweeps{iWav}.commonGreenAxis();
                    obj.sweeps{iWav} = obj.sweeps{iWav}.setGreenVelAxis(axisLocal);
                    greenAxisPerWav{iWav} = axisLocal;
                end

                %%% Phase 3 -- Green, flattened across the whole grid again; each element
                %%% uses ITS OWN wavelength's axis, looked up by flat index through wavIdx.
                retentionLocal = retentionLevel;
                greenThrLocal  = [];                   % 'green' at its own default unless chosen
                if strcmp(method, 'green'); greenThrLocal = obj.sawThreshold; end
                csFlat = CalcCommon.forEachElement(csFlat, useP, ...
                    @(cs, ell) CalcCommon.greenStep(cs, greenAxisPerWav{wavIdx(ell)}, retentionLocal, ...
                    [], greenThrLocal, precompPerWav{wavIdx(ell)}));
            end

            obj.sweeps = CalcCommon.scatterToSweeps(obj.sweeps, csFlat, wavIdx, angIdx);
            for iWav = 1:nWav
                obj.sweeps{iWav} = obj.sweeps{iWav}.finalizeMaps();
                if strcmp(mode, 'full') && strcmp(method, 'ritz')
                    obj.sweeps{iWav} = obj.sweeps{iWav}.classify('ritz', obj.sawThreshold);
                end
            end

            obj = obj.assembleMaps();
        end

        function obj = assembleMaps(obj)
        %%% 2D re-indexing of each wavelength-slice's already-computed 1D maps
        %%% (sweeps{iWav}.sawVelMap etc.) -- no new classification logic.
        %%% >> obj = assembleMaps(obj)
            nWav    = numel(obj.sweeps);
            nAngles = numel(obj.angles);
            obj.sawVelMap2D           = cell(nWav, nAngles);
            obj.sawModeNoMap2D        = cell(nWav, nAngles);
            obj.sawTypeMap2D          = cell(nWav, nAngles);
            obj.significantModesMap2D = cell(nWav, nAngles);

            for iWav = 1:nWav
                obj.sawVelMap2D(iWav, :)           = obj.sweeps{iWav}.sawVelMap(:)';
                obj.sawModeNoMap2D(iWav, :)        = obj.sweeps{iWav}.sawModeNoMap(:)';
                obj.sawTypeMap2D(iWav, :)          = obj.sweeps{iWav}.sawTypeMap(:)';
                obj.significantModesMap2D(iWav, :) = obj.sweeps{iWav}.significantModesMap(:)';
            end
        end
    end

    methods (Access = protected)
        function cp = copyElement(obj)
        %%% Deep-copy override (matlab.mixin.Copyable)
        %%% >> cp = copyElement(obj)
        %%% obj.sweeps holds handle-class CalcAngSweepLayered objects, which the
        %%% default shallow copy would alias. One level of cellfun(@copy, ...) here is
        %%% enough for a fully independent tree, since CalcAngSweepLayered's own
        %%% copyElement deep-copies its cs{} array in turn.
            cp = copyElement@matlab.mixin.Copyable(obj);
            cp.sweeps = cellfun(@copy, obj.sweeps, 'UniformOutput', false);
        end
    end
end
