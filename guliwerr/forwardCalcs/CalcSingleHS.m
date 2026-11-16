classdef CalcSingleHS < handle & matlab.mixin.Copyable
%%% CalcSingleHS — single-angle container for the homogeneous half-space
%%% forward calculation (Christoffel + GreenSolver + RitzSolverSAW).
%%%
%%% Given one Structure (half-space) and one propagation angle, runs the three
%%% solvers for that direction and holds the assembled results together. Usable
%%% standalone, and the per-angle building block CalcAngSweepHS drives across an
%%% angle array.
%%%
%%% Handle class: holds potentially large per-object data (Ritz eigenvectors,
%%% Green spectrum) and the GUI viewers need write-back.
%%%
%%% identifyModes() delegates to ModeClassification — Green-primary,
%%% Ritz-assisted classification when calcGreen() has run ('full' mode, existence
%%% by GreenScore), the gated Ritz modal sum otherwise ('ritzOnly' mode, existence by RitzScore).
%%%
%%% GuliweRR v6.0

    properties (SetAccess = protected)
        structure       Structure     = Structure()   % rotated structure for this angle;
                                                        % CalcParams lives inside it already
                                                        % (structure.calcParams, one per layer)
        angle(1,1)      double        = 0             % propagation angle [deg]
        structureBase   Structure     = Structure()   % the UNROTATED structure passed in:
                                                        % GreenScore.reference is built from it,
                                                        % as a sweep builds its own

        %%% calcFlags
        calcMode        char          = 'full'         % 'full' (default) | 'ritzOnly'

        %%% ritzSolver — calcRitz() result
        ritzSolver      RitzSolverSAW = RitzSolverSAW()

        %%% christoffel — calcLim() result
        longLimVel      double  = []   % scalar — longitudinal limiting velocity
        longLimPol      double  = []   % [3x1]  — polarisation at longLimVel
        shearLimVels    double  = []   % [nShear x 1] — variable count; wave types 2 and 3 concatenated
        shearLimPols    double  = []   % [3 x nShear]
        shearLimThetas  double  = []   % [nShear x 1] sagittal root angles theta* of
                                        % shearLimVels (CalcCommon.limitingAndBulk); the
                                        % Ritz-only Rayleigh label reads the min-shear one
        bulkVels        double  = []   % [3x1] = [vL; vT1; vT2], descending — cheap direct
        bulkPols        double  = []   % [3x3] columns = corresponding polarisations

        %%% greenSolver — calcGreen() result: ONE solve at SawSpectrum.ETA = 3e-3. Only
        %%% |G13_eff| is kept -- the TGS response (in-plane force, out-of-plane
        %%% displacement, +-k symmetrised), which the existence score also reads; any
        %%% other component is recomputed on request (greenComponent())
        greenVelAxis    double  = []   % [nVel x 1] — axis the solve ran on: this angle's
                                        % own axis (greenAxis(): 0.3 km/s to 1.05 x
                                        % longLimVel + 0.1 km/s in greenVelStep steps), or
                                        % a sweep's common (longer) one
        greenAbsG13Eff  single  = []   % [nVel x 1] |(G13_+ - G13_-)/2| = |(G13_+ - G31_+)/2|
                                        % on greenVelAxis

        %%% GreenScore existence criterion ('full' mode) --
        %%% calcGreen() keeps the complex G13_eff on the band, calcGreenScore() adds the
        %%% probe slope, the scores and the 'low' reference (both part of calcAll('full'))
        greenScoreVelAxis double = []   % [nBand x 1] band of greenVelAxis (GreenScore.candidates)
        greenScoreG13Eff  double = []   % [nBand x 1] complex G13_eff on greenScoreVelAxis
        greenScorePeaks   struct = struct([])   % scored candidates (GreenScore.direction);
                                        % survive eigenvector pruning, so reclassify() works
        greenScoreRef            = []   % GreenScore.reference struct (vRef, logStrengthTop) --
                                        % shared across a sweep, or memoised per structure

        %%% Ritz-only existence criterion -- calcRitzScore() result (part of
        %%% calcAll('ritzOnly'))
        ritzScoreNum     double = []   % [nKeep x 1] complex modal-sum weights u3 conj(u1)
                                        % of modes 1..nKeep (RitzScore.modalNumerators);
                                        % survives eigenvector pruning, so reclassify() works
        ritzScoreVelAxis double = []   % [nSaw x 1] uniform axis (SawSpectrum.velocityAxis)
        ritzScoreRef            = []   % RitzScore.ownReference struct (vRef, floor) of
                                        % this direction, + slabRes [km/s]
                                        % (RitzScore.slabResolution, the Rayleigh label's
                                        % tolerance above min shear)
        ritzScoreCands   struct = struct([])   % scored candidates (RitzScore.direction:
                                        % local mode, dScore, slope, windowMass, ...),
                                        % threshold-free; survive eigenvector pruning, so
                                        % reclassify() works

        %%% modeClassification — identifyModes() result
        classification          = []   % ModeClassification once identifyModes() has run
        classificationVersion(1,1) double = 0   % ModeClassification.CLASSIFICATION_VERSION
                                        % this result was classified with (0: unclassified,
                                        % or saved before the stamp existed); see isStale()
        activeMethod     char   = ''   % the method the classification holds: 'green' |
                                        % 'refined' | 'ritz' (classify())
        activeThreshold         = []   % its standard: 'high' | 'low' | numeric
        limitingEntries  struct = struct([])   % the limiting-velocity entries, associated
                                        % once per calculation (they need eigenvectors) and
                                        % carried by every 'green' / 'refined' classification
        limitingDone     logical = false   % limitingEntries hold this calculation's

        %%% Refined candidates -- calcRefined() result (GreenSpotScore; run by
        %%% classify('refined') on first use; 'full' half-space only). struct([]) until
        %%% then; cleared by calcRitz()/calcGreen() (so by calcAll() and a sweep's calc()),
        %%% kept by classify(). Fields: threshold, ref (GreenSpotScore.reference), peaks
        %%% (every candidate: vel, velDet, posLinked, prominence, relWidth, height, fU1,
        %%% scoreHigh, scoreLow, modeNo, matchScore, confident, type, accepted), info (nVel, velEnd, minShearSpot: the
        %%% spot-weighted minShear the Rayleigh label compares with)
        refinedSaw       struct = struct([])

        %%% Spot-averaged Green's function -- calcSpotGreen() result (SpotGreen, on
        %%% request): |G| as a measurement with a finite pattern sees it (ProbeKernel). A
        %%% picture only; never read by any classification. Needs only the structure, so
        %%% it is kept by calcAll()/reclassify().
        spotGreen        single = []   % [nVel x 1] |X| on spotGreenVelAxis
        spotGreenVelAxis double = []   % [nVel x 1] km/s (GreenScore.gridAxis of SpotGreen's vEnd)
        spotGreenInfo    struct = struct([])   % SpotGreen.compute info: kernel, steps, eta,
                                               % quadrature checks
    end

    properties
        %%% Public SetAccess so CalcAngSweepHS can call applyRetention() directly
        %%% between phases. See docs/reference/architecture/forwardCalcs.md, "retentionLevel".
        retentionLevel  char    = 'full'   % 'full' (default) | 'saw_only' | 'none'

        %%% Classification method calcAll() runs: '' (default: follows each run's calcMode
        %%% -- 'green' on 'full', 'ritz' on 'ritzOnly'; the method in effect is
        %%% activeMethod) | 'green' (GreenScore on the Green's function, needs calcMode
        %%% 'full') | 'refined' (GreenSpotScore on the spot-averaged Green's function,
        %%% 'full' only; computed on first use) | 'ritz' (the gated Ritz modal sum; the only
        %%% method on 'ritzOnly'). Set before calcAll(); classify() switches the active
        %%% classification without changing it.
        sawClassMethod  char    = ''
        %%% Its standard: 'high' (sensitive -- everything a good TGS map shows) | 'low'
        %%% (distinct modes only) | numeric on the method's own scale (needs
        %%% sawClassMethod set explicitly). [] (default) = the method's default ('high'
        %%% for 'green' / 'refined', 'low' for 'ritz').
        sawThreshold            = []

        %%% Green axis step [km/s] of 'full' mode: it sets the picture resolution (the
        %%% existence score, GreenScore, is validated up to GreenScore.MAX_STEP = 16 m/s;
        %%% larger steps are capped). Default GreenScore.DEFAULT_STEP (6 m/s).
        greenVelStep(1,1) double = GreenScore.DEFAULT_STEP

        %%% Alternative to greenVelStep (e.g. positional constructor calls): the number
        %%% of points of this angle's Green axis; the step follows from it (capped at
        %%% GreenScore.MAX_STEP). NaN (default) = use greenVelStep.
        nVelGreen(1,1)  double  = NaN
    end

    properties (SetAccess = protected, Transient)
        %%% Angle-independent Ritz precompute (nullConstrMat, dKdC_cell, d2KdCdD_cell),
        %%% built once at construction. These depend only on calcParams/density, not
        %%% on angle or dirMat, so a sweep can build them once for all angles.
        %%% Transient: NOT written by save() -- in memory a sweep's copies share one
        %%% array, but a .mat file would store one full copy per angle (most of the
        %%% file) and load() would bring them back unshared. Empty after load();
        %%% read it through ritzPrecompute(), which rebuilds it on first use.
        %%% Transient also keeps it out of every thread-pool transfer (a parfor
        %%% serialises the object), so it is EMPTY on a pool worker too: a sweep hands
        %%% its shared copy to calcRitz()/calcAll() explicitly instead.
        precompStructForRitz struct = struct()
    end

    methods
        function obj = CalcSingleHS(structure, angle, precompStruct, calcMode, retentionLevel, nVelGreen)
        %%% Construct CalcSingleHS for one propagation angle
        %%% >> obj = CalcSingleHS(structure, angle, precompStruct, calcMode, retentionLevel, nVelGreen)
        %%%
        %%% INPUT:  1. structure       .. Structure object, UNROTATED (the constructor
        %%%                              rotates it internally — callers never rotate
        %%%                              by hand)
        %%%  (opt)  2. angle           .. propagation angle [deg] (default = 0)
        %%%  (opt)  3. precompStruct   .. angle-independent Ritz precompute fields
        %%%                              (nullConstrMat, dKdC_cell, d2KdCdD_cell) — built
        %%%                              from the rotated structure if omitted, so a bare
        %%%                              CalcSingleHS(structure, angle) is fully usable
        %%%                              standalone
        %%%  (opt)  4. calcMode        .. 'full' (default) | 'ritzOnly'
        %%%  (opt)  5. retentionLevel  .. 'full' (default) | 'saw_only' | 'none'
        %%%  (opt)  6. nVelGreen       .. number of Green axis points (sets the step,
        %%%                              capped at GreenScore.MAX_STEP); default: the
        %%%                              greenVelStep property (6 m/s)

            if nargin == 0
                return
            end

            validateattributes(structure, {'Structure'}, {'scalar'}, 'CalcSingleHS', 'structure');
            if nargin < 2 || isempty(angle); angle = 0; end
            validateattributes(angle, {'double'}, {'scalar', 'real', 'finite'}, 'CalcSingleHS', 'angle');
            if nargin < 4; calcMode = []; end
            if nargin < 5; retentionLevel = []; end
            CalcCommon.validateConstructorFlags(calcMode, retentionLevel, 'CalcSingleHS:constructor', {'full', 'saw_only', 'none'});

            obj.angle     = angle;
            obj.structureBase = structure;
            obj.structure = CalcCommon.rotateStructure(structure, angle);

            if nargin < 3; precompStruct = struct(); end
            %%% Fill-if-missing: a sweep hands its shared precompute down and only
            %%% the gaps are built; a bare CalcSingleHS(structure, angle) builds
            %%% everything itself and stays standalone.
            obj.precompStructForRitz = CalcCommon.ritzSAWPrecompute(obj.structure, precompStruct);

            if nargin >= 4 && ~isempty(calcMode)
                obj.calcMode = calcMode;
            end
            if nargin >= 5 && ~isempty(retentionLevel)
                obj.retentionLevel = retentionLevel;
            end
            if nargin >= 6 && ~isempty(nVelGreen)
                validateattributes(nVelGreen, {'double'}, {'scalar', 'real', '>=', 2}, ...
                    'CalcSingleHS', 'nVelGreen');
                obj.nVelGreen = nVelGreen;
            end
        end

        function obj = calcRitz(obj, precomp)
        %%% Run RitzSolverSAW for this angle
        %%% >> obj = calcRitz(obj)
        %%% >> obj = calcRitz(obj, precomp)
        %%%
        %%% INPUT: (opt) precomp .. the angle-independent Ritz precompute
        %%%                         (CalcCommon.ritzSAWPrecompute), as a sweep holds it.
        %%%                         Default: this object's own (ritzPrecompute()).
        %%%
        %%% A sweep always passes its shared copy: on a thread-pool worker this object's
        %%% own copy is empty (Transient properties do not cross to a worker), and
        %%% rebuilding it there costs more than the Ritz solve, once per angle. A passed
        %%% precompute is completed if partial (free when it is complete) and is not
        %%% stored on obj.
        %%%
        %%% withDiagnostics = true stores the per-mode diagnostics (freqDisp, oopSurfDispl,
        %%% depthInsensitivity) on ritzSolver.diagnostics, which ModeClassification then
        %%% reads rather than recomputing.
        %%%
        %%% obj.structure is already rotated, so pass 0 here, not obj.angle, or
        %%% RitzSolverSAW rotates it again — see docs/reference/architecture/forwardCalcs.md,
        %%% "Rotation convention: pass 0 to solvers after construction".
            if nargin < 2 || isempty(precomp)
                precomp = obj.ritzPrecompute();
            else
                precomp = CalcCommon.completeRitzPrecompute(obj.structure, precomp, class(obj));
            end
            rs = RitzSolverSAW(obj.structure, 0);
            obj.ritzSolver = rs.calc(true, precomp);
            obj.refinedSaw = struct([]);   % its mode identity refers to the previous Ritz modes
            obj.limitingEntries = struct([]); obj.limitingDone = false;   % new eigenvectors
        end

        function precomp = ritzPrecompute(obj)
        %%% The Ritz precompute, rebuilt first if missing (e.g. after load())
        %%% >> precomp = ritzPrecompute(obj)
        %%%
        %%% OUTPUT: precomp .. precompStructForRitz with every field populated
        %%%
        %%% Fill-if-missing (CalcCommon.ritzSAWPrecompute), so a populated
        %%% precompute is returned unchanged at no cost. A rebuilt one is stored
        %%% back on obj, so it is built at most once per object. A rebuild on a pool
        %%% worker warns (CalcCommon.warnRitzPrecomputeBuiltOnWorker): the sweep
        %%% should have passed its copy to calcRitz().
            %%% No own "is it complete?" test: ritzSAWPrecompute checks every field
            %%% itself, so a field added there later is covered here automatically.
            obj.precompStructForRitz = CalcCommon.completeRitzPrecompute(obj.structure, ...
                obj.precompStructForRitz, class(obj));
            precomp = obj.precompStructForRitz;
        end

        function obj = calcLim(obj)
        %%% Christoffel limiting velocities and cheap direct bulk velocities
        %%% >> obj = calcLim(obj)
        %%%
        %%% Fills longLimVel/longLimPol/shearLimVels/shearLimPols (LimitingVels
        %%% scan) and bulkVels/bulkPols (one Christoffel eigendecomposition). The
        %%% structure's dirMat already carries the propagation-angle rotation, so
        %%% the elastic tensor is taken with zero further rotation.

            surfIdx = obj.structure.nLay;   %%% surface layer -- for a homogeneous
                                            %%% half-space the whole medium IS this
                                            %%% one material, so these are genuine
                                            %%% SYSTEM-WIDE limiting velocities. NOT 1,
                                            %%% which is CalcSingleLayered's substrate
                                            %%% convention and wrong here
            lim = CalcCommon.limitingAndBulk(obj.structure.mater(surfIdx));

            obj.longLimVel   = lim.longLimVel;
            obj.longLimPol   = lim.longLimPol;
            obj.shearLimVels = lim.shearLimVels;
            obj.shearLimPols = lim.shearLimPols;
            obj.shearLimThetas = lim.shearLimTheta;
            obj.bulkVels     = lim.bulkVels;
            obj.bulkPols     = lim.bulkPols;
        end

        function velAxis = greenAxis(obj)
        %%% This angle's own Green axis: 0.3 km/s to 1.05 x longLimVel + 0.1 km/s (needs calcLim)
        %%% >> velAxis = greenAxis(obj)
        %%%
        %%% Step: greenVelStep, or the step that gives nVelGreen points if that is set,
        %%% capped at GreenScore.MAX_STEP. A sweep's common axis has the same origin and step, so this axis is its
        %%% leading part.
            if isempty(obj.longLimVel)
                error('CalcSingleHS:greenAxis:noLim', ...
                    'calcLim() must run first: the axis runs to 1.05 x longLimVel.');
            end
            vBound = obj.scoreAxisBound();
            step = obj.greenVelStep;
            if isfinite(obj.nVelGreen)
                step = (vBound + SawSpectrum.EDGE_TAPER - SawSpectrum.VEL_MIN) / (obj.nVelGreen - 1);
            end
            step = GreenScore.resolveStep(step, 'CalcSingleHS:greenStepCapped');
            velAxis = GreenScore.gridAxis(vBound, step);
        end

        function obj = calcGreen(obj, velAxis)
        %%% The one Green solve of this angle: |G13_eff| + the existence-score input
        %%% >> obj = calcGreen(obj)
        %%% >> obj = calcGreen(obj, velAxis)
        %%%
        %%% INPUT: (opt) velAxis .. phase-velocity axis [km/s]. Default: this angle's
        %%%                         own axis, greenAxis() (needs calcLim() first). A
        %%%                         sweep passes its common axis, which runs further.
        %%%
        %%% One GreenSolver call at SawSpectrum.ETA (CalcCommon.greenScoreSpectrum).
        %%% Stores greenVelAxis and |G13_eff| (greenAbsG13Eff, single) and the complex
        %%% G13_eff the existence score reads: the band [0.8 vS, 1.03 vTf] of a uniform
        %%% velAxis (greenScoreVelAxis/greenScoreG13Eff); the band of a sweep's common axis
        %%% equals a standalone angle's, so a sweep angle scores as a standalone one.
        %%% On any other axis (or before calcLim) the input is left empty and
        %%% calcGreenScore() solves it separately.
        %%%
        %%% Same double-rotation trap as calcRitz(); CalcCommon.greenScoreSpectrum
        %%% is where the angle-0 argument is enforced.
            if nargin < 2 || isempty(velAxis)
                velAxis = obj.greenAxis();
            end
            g = CalcCommon.greenScoreSpectrum(obj.structure, velAxis, 'CalcSingleHS.calcGreen');

            obj.refinedSaw     = struct([]);   % a new calculation: refine again on request
            obj.greenVelAxis   = g.velAxis;
            obj.greenAbsG13Eff = g.absG13Eff;
            obj.greenScoreVelAxis = []; obj.greenScoreG13Eff = [];
            if isempty(obj.longLimVel); return; end

            v = g.velAxis; dv = v(2) - v(1);
            [vS, vTf] = obj.shearRange();
            b = v >= GreenScore.BAND_LO * vS - 1.5 * dv & v <= GreenScore.BAND_SOLVE_HI * vTf + 1.5 * dv;
            if nnz(b) >= 5 && max(abs(diff(v(b)) - dv)) < 1e-9
                obj.greenScoreVelAxis = v(b);
                obj.greenScoreG13Eff  = g.G13Eff(b);
            end
        end

        function obj = calcGreenScore(obj, ref)
        %%% Score this angle's candidates (GreenScore): probe slope, scores, 'low' reference
        %%% >> obj = calcGreenScore(obj)
        %%% >> obj = calcGreenScore(obj, ref)
        %%%
        %%% INPUT: (opt) ref .. GreenScore.reference struct shared by a sweep; if omitted or
        %%%                     empty, that of the UNROTATED structure this object was built
        %%%                     from (structureBase) -- the one a sweep over it uses, and
        %%%                     memoised, so further angles of the structure pay nothing
        %%%
        %%% Requires calcLim(). Cost beyond calcGreen(): two probe solves of ~7 points per
        %%% candidate at angle +- GreenScore.PROBE_OFFSET. If calcGreen() left no band (not
        %%% run, or a custom non-uniform axis), the band is solved here.
            if isempty(obj.longLimVel)
                error('CalcSingleHS:calcGreenScore:noLim', 'calcLim() must run first.');
            end
            if nargin < 2 || isempty(ref)
                ref = GreenScore.reference(obj.structureBase);   % = the sweep's; memoised
            end
            obj.greenScoreRef = ref;
            [vS, vTf] = obj.shearRange();
            if isempty(obj.greenScoreG13Eff)
                v = GreenScore.bandAxis(vS, vTf, GreenScore.resolveStep(obj.greenVelStep, ...
                    'CalcSingleHS:greenStepCapped'));
                g = CalcCommon.greenScoreSpectrum(obj.structure, v, 'CalcSingleHS.calcGreenScore');
                obj.greenScoreVelAxis = g.velAxis;
                obj.greenScoreG13Eff  = g.G13Eff;
            end
            obj.greenScorePeaks = GreenScore.direction(obj.structureBase, obj.angle, obj.greenScoreG13Eff, ...
                obj.greenScoreVelAxis, vS, vTf, ref.vRef, ref);
        end

        function absG = greenComponent(obj, component)
        %%% |G| of any Green component on greenVelAxis, recomputed now (not stored)
        %%% >> absG = greenComponent(obj, 'G33')
        %%% >> absG = greenComponent(obj, 'G13')      % one-sided G13_+
        %%%
        %%% INPUT:  component .. 'GIJ' (one-sided G_IJ(+k); I = force, J = displacement
        %%%                      direction, 1 = propagation, 3 = surface normal) or
        %%%                      'GIJeff' (+-k symmetrised); 'G13eff' returns the
        %%%                      stored greenAbsG13Eff without a solve
        %%% OUTPUT: absG      .. [nVel x 1] single on greenVelAxis
        %%%
        %%% A calc result stores only |G13_eff|; every other component is one Green
        %%% solve on the same axis and damping (CalcCommon.greenComponentSpectrum).
        %%% Requires calcGreen() (or calcAll('full')) to have run.
            absG = CalcCommon.greenComponentOf(obj, component);
        end

        function obj = calcSpotGreen(obj, kernel, varargin)
        %%% |G| as a finite-pattern measurement sees it: coherent sum over a ProbeKernel
        %%% Experimental: tried on Fe3Al(316) only, not validated against measurements; the
        %%% interface may change.
        %%% >> obj = calcSpotGreen(obj, kernel)
        %%% >> obj = calcSpotGreen(obj, kernel, 'DirStep', 0.1, 'RelKStep', 1e-3, 'Eta', 3e-3, 'VelStep', 0.004)
        %%%
        %%% INPUT:  1. kernel .. ProbeKernel (e.g. ProbeKernel.tgs(Lambda, spotRadii, probeRadii),
        %%%                      ProbeKernel.gaussian(sigmaRelK, sigmaAngDeg))
        %%%  (opt)  Name-Value:
        %%%           'DirStep'  .. direction step [deg] (default kernel.dirStep)
        %%%           'RelKStep' .. relative |k| step (default kernel.relKStep)
        %%%           'Eta'      .. Green damping (default SawSpectrum.ETA, as the stored map)
        %%%           'VelStep'  .. solve/output axis step [km/s] (default 0.004)
        %%% OUTPUT: obj with spotGreen, spotGreenVelAxis, spotGreenInfo
        %%%
        %%% One Green solve per kernel direction (SpotGreen.compute); the |k| spread costs
        %%% no extra solve. Equal, bit for bit, to this angle's column of a sweep's
        %%% calcSpotGreenMap with the same kernel and options. Plot it with
        %%% Vis.plotGreenCurve(obj, [], [], 'spot'). Any calcMode; half-space only.
            if obj.structure.nLay ~= 1
                error('CalcSingleHS:calcSpotGreen:layered', 'calcSpotGreen is for a half-space (nLay = 1) only.');
            end
            opts = SpotGreen.options(kernel, varargin{:});
            [obj.spotGreen, obj.spotGreenVelAxis, obj.spotGreenInfo] = ...
                SpotGreen.compute(obj.structureBase, obj.angle, kernel, opts, "auto");
        end

        function obj = calcRitzScore(obj, precomp)
        %%% Ritz-only existence input: the modal-sum weights of this angle's Ritz modes,
        %%% the uniform score axis, this direction's own floor reference, and the scored
        %%% half-space candidates
        %%% >> obj = calcRitzScore(obj)
        %%% >> obj = calcRitzScore(obj, precomp)
        %%%
        %%% INPUT: (opt) precomp .. the Ritz precompute calcRitz() used (its dKdC_cell gives
        %%%                         the modes' angular slopes); default: ritzPrecompute()
        %%%
        %%% Requires calcRitz() (eigenvectors present, i.e. before applyRetention) and
        %%% calcLim() (the axis and the kept modes run to the longitudinal limiting
        %%% velocity). Fills ritzScoreNum/VelAxis/Ref and ritzScoreCands
        %%% (RitzScore.direction: candidates, local supporting mode, dephased score,
        %%% window mass; everything a reclassify() needs). Cost beyond that: ~10 ms, no
        %%% other direction is solved.
            if isempty(obj.longLimVel)
                error('CalcSingleHS:calcRitzScore:noLim', 'calcLim() must run first.');
            end
            if nargin < 2 || isempty(precomp)
                precomp = obj.ritzPrecompute();
            end
            dKdC = [];   % a partial precompute: angleSlope builds it (as calcRitz would)
            if isfield(precomp, 'dKdC_cell'); dKdC = precomp.dKdC_cell; end
            [vS, vTf] = obj.shearRange(); vL = obj.longLimVel;
            RitzScore.checkBasisHS(obj.structure.calcParams(1));
            [num, E1, r, Ebot] = RitzScore.modeProfiles(obj.ritzSolver, vL);
            obj.ritzScoreNum     = num;
            obj.ritzScoreVelAxis = SawSpectrum.velocityAxis(SawSpectrum.AXIS_HI_X_VL * vL, RitzScore.VEL_STEP);
            ref = RitzScore.ownReference(obj.structure, obj.ritzSolver.vels(1:numel(num)), ...
                num, obj.ritzScoreVelAxis, RitzScore.HS_FLOOR_LO * vS, vTf);
            [~, iS] = min(obj.shearLimVels);
            ref.slabRes = RitzScore.slabResolution(obj.structure.mater(1), vS, ...
                obj.shearLimThetas(iS), obj.structure.calcParams(1).relDepth);
            obj.ritzScoreRef   = ref;
            obj.ritzScoreCands = RitzScore.direction(obj.ritzSolver, num, E1, r, ...
                obj.ritzScoreVelAxis, ref, vS, vL, dKdC, Ebot);
        end

        function obj = classify(obj, method, threshold)
        %%% Set the active classification: method and standard (re-thresholds stored candidates)
        %%% >> obj = classify(obj)                      % the active method, its default standard
        %%% >> obj = classify(obj, method)
        %%% >> obj = classify(obj, method, threshold)
        %%%
        %%% INPUT: (opt) method    .. 'green' | 'refined' | 'ritz' (default: the active one)
        %%%        (opt) threshold .. 'high' | 'low' | numeric on the method's own scale
        %%%                           (default: the method's default standard,
        %%%                           CalcCommon.defaultThreshold)
        %%% OUTPUT: obj with classification (significantModes), activeMethod,
        %%%         activeThreshold and classificationVersion set
        %%%
        %%% The candidates of every method are stored side by side, so switching only
        %%% re-thresholds them. 'refined' solves its spot directions on first use (11
        %%% Green solves, ~0.4 s, plus the memoised reference on a structure's first
        %%% use) and prints the cost. 'green' and 'refined' need calcMode 'full'. A stale
        %%% result (another classification stamp, see isStale) is refused: recompute it.
            if nargin < 3; threshold = []; end
            if nargin < 2 || isempty(method); method = obj.activeMethod; end
            if isempty(method); method = obj.sawClassMethod; end
            if isempty(method)   % neither classified nor set: calcMode's own
                method = CalcCommon.resolveSawClassMethod('', threshold, obj.calcMode, ...
                    false, 'CalcSingleHS:classify');
            end
            method = CalcCommon.validateSawClassMethod(method, false, 'CalcSingleHS:classify');
            if isempty(threshold)
                threshold = CalcCommon.defaultThreshold(method, false);
            end
            CalcCommon.validateSawThreshold(threshold, 'CalcSingleHS:classify');
            CalcCommon.refuseStale(obj, 'CalcSingleHS:classify');
            if ~strcmp(method, 'ritz') && ~strcmp(obj.calcMode, 'full')
                error('CalcSingleHS:classify:needsFull', ...
                    "Method '%s' needs calcMode='full'; this result is '%s'.", method, obj.calcMode);
            end
            if strcmp(method, 'refined') && isempty(obj.refinedSaw)
                t0 = tic;
                obj = obj.calcRefined();
                fprintf('refined classification: spot directions solved in %.1f s\n', toc(t0));
            end
            obj = obj.identifyModes(method, threshold);
        end

        function obj = reclassify(obj, threshold)
        %%% Re-threshold the active classification (no recomputation)
        %%% >> obj = reclassify(obj, threshold)
        %%%
        %%% = classify(obj, activeMethod, threshold).
            obj = obj.classify(obj.activeMethod, threshold);
        end

        function tf = isStale(obj)
        %%% Whether this result's classification carries another stamp than the current one
        %%% >> tf = isStale(obj)
            tf = CalcCommon.isStale(obj);
        end

        function obj = calcRefined(obj, dirCache, ref)
        %%% Candidates of the refined classification (GreenSpotScore) in refinedSaw
        %%% >> obj = calcRefined(obj)
        %%% >> obj = calcRefined(obj, dirCache, ref)
        %%%
        %%% INPUT:  (opt) 1. dirCache .. GreenSpotScore.solveDirections cache holding this angle's
        %%%                              kernel directions (a sweep passes its own); default:
        %%%                              the 11 directions are solved here
        %%%         (opt) 2. ref      .. GreenSpotScore.reference of structureBase (a sweep passes
        %%%                              its own); default: GreenSpotScore.reference(structureBase),
        %%%                              memoised -- the same one a sweep uses
        %%% OUTPUT: obj.refinedSaw: threshold, ref, peaks (every candidate, with
        %%%         its identity and label), info. The classification is set by
        %%%         classify('refined').
        %%%
        %%% Needs calcAll('full') on a half-space (each spot peak is identified via its
        %%% nearest calcGreenScore pole); errors CalcSingleHS:calcRefined:notFull / :layered /
        %%% :notCalculated / :noGreenScore. A new calculation
        %%% (calcRitz/calcGreen/calcAll) clears it. Cost: 11 Green solves (~0.4 s), plus the
        %%% reference on a structure's first use (up to 132 solves, memoised). Identical,
        %%% bit for bit, to the same angle in any sweep.
            if ~strcmp(obj.calcMode, 'full')
                error('CalcSingleHS:calcRefined:notFull', ...
                    "The refined classification needs calcMode 'full'; this result is '%s'.", obj.calcMode);
            end
            if obj.structure.nLay ~= 1
                error('CalcSingleHS:calcRefined:layered', ...
                    'The refined classification is for a half-space (nLay = 1) only.');
            end
            if isempty(obj.longLimVel) || isempty(obj.ritzSolver.vels)
                error('CalcSingleHS:calcRefined:notCalculated', ...
                    'calcAll(''full'') must run first.');
            end
            if isempty(obj.greenScoreVelAxis) || isempty(obj.greenScoreRef)
                error('CalcSingleHS:calcRefined:noGreenScore', ...
                    'calcGreen() and calcGreenScore() must run first (calcAll(''full'') runs both): the spot peaks are identified via the single-direction poles.');
            end
            if nargin < 3 || isempty(ref)
                ref = GreenSpotScore.reference(obj.structureBase);
            end
            if nargin < 2 || isempty(dirCache)
                o = GreenSpotScore.kernel();
                dirCache = GreenSpotScore.solveDirections(obj.structureBase, obj.angle + o, "never");
            end
            [X, v, minShearSpot] = GreenSpotScore.average(dirCache, obj.angle);
            [vS, vTf] = obj.shearRange();
            pk = GreenSpotScore.direction(X, v, vS, vTf, ...
                obj.longLimVel, ref.vRef, ref);
            poleVels = [];   % the single-direction poles (GreenScore candidates) of this angle
            if ~isempty(obj.greenScorePeaks); poleVels = [obj.greenScorePeaks.vel]; end
            pk = GreenSpotScore.identify(pk, obj.ritzSolver, minShearSpot, poleVels);
            acc = GreenSpotScore.accepted(pk, 'high');
            for j = 1:numel(acc)
                pk(j).accepted = acc(j);
            end
            obj.refinedSaw = struct('threshold', 'high', ...
                'ref', ref, 'peaks', pk, 'info', struct('nVel', numel(v), 'velEnd', v(end), ...
                'minShearSpot', minShearSpot));
        end

        function T = refinedSawAgreement(obj)
        %%% Default vs refined surface-wave entries of this angle, paired one to one
        %%% >> T = refinedSawAgreement(obj)
        %%%
        %%% OUTPUT: T .. table, one row per pair or unpaired entry: angle, defaultVel
        %%%              (greenPeakVel of a Rayleigh/pseudoRayleigh significantModes
        %%%              entry), refinedVel (an accepted refinedSaw peak), dv = refined -
        %%%              default [km/s], defaultType, refinedType, status ('both' |
        %%%              'defaultOnly' | 'refinedOnly')
        %%%
        %%% Pairs: nearest first, within max(1.2 % x v, 30 m/s) (the velocity
        %%% tolerance of a match). The 'green' entries are those at the standard the refined
        %%% candidates were last accepted at (a numeric one: 'green''s default). Needs the
        %%% refined candidates (classify('refined') or calcRefined()) first.
            if isempty(obj.refinedSaw)
                error('CalcSingleHS:refinedSawAgreement:notRefined', ...
                    'classify(''refined'') must run first.');
            end
            thr = obj.refinedSaw.threshold;
            if ~ischar(thr); thr = []; end
            cp = copy(obj);
            cp = cp.classify('green', thr);
            sm = cp.classification.significantModes;
            sm = sm(ismember({sm.type}, {'Rayleigh', 'pseudoRayleigh'}));
            pk = obj.refinedSaw.peaks([obj.refinedSaw.peaks.accepted]);
            T = CalcSingleHS.pairEntries(obj.angle, [sm.greenPeakVel], {sm.type}, [pk.vel], {pk.type});
        end

        function obj = identifyModes(obj, method, thr)
        %%% Classify with a method at a standard (ModeClassification)
        %%% >> obj = identifyModes(obj)                 % sawClassMethod at sawThreshold
        %%% >> obj = identifyModes(obj, method, thr)
        %%%
        %%% 'green'  .. GreenScore candidates (needs calcGreenScore());
        %%% 'refined'.. GreenSpotScore candidates (needs calcRefined());
        %%% 'ritz'   .. the gated Ritz modal sum (needs calcLim() and calcRitzScore()).
        %%% method empty: sawClassMethod, or calcMode's own; thr empty = the method's
        %%% default standard. The limiting-velocity
        %%% entries are associated once, on the first Green-path classification after a
        %%% calculation (they need the eigenvectors), stored in limitingEntries and
        %%% carried over by every later 'green' / 'refined' classification; the 'ritz'
        %%% classification has none, so it equals a 'ritzOnly' result's.
        %%% longLimPol is forwarded so the longitudinal limiting velocity is
        %%% cross-linked to a Ritz mode as the shear branches are.
            if nargin < 3; thr = obj.sawThreshold; end
            if nargin < 2 || isempty(method); method = obj.sawClassMethod; end
            if isempty(method)
                method = CalcCommon.resolveSawClassMethod('', thr, obj.calcMode, false, ...
                    'CalcSingleHS:identifyModes');
            end
            method = CalcCommon.validateSawClassMethod(method, false, 'CalcSingleHS:identifyModes');
            if isempty(thr); thr = CalcCommon.defaultThreshold(method, false); end
            switch method
                case 'ritz'
                    if isempty(obj.ritzScoreVelAxis) || isempty(obj.longLimVel)
                        error('CalcSingleHS:identifyModes:noRitzScore', ['Ritz classification ' ...
                            'needs calcLim() and calcRitzScore() first (calcAll runs both).']);
                    end
                    mc = ModeClassification('threshold', thr);
                    obj.classification = mc.calc(obj.ritzSolver, obj.structure, [], [], ...
                        struct('cands', obj.ritzScoreCands));
                otherwise   % 'green' | 'refined'
                    christoffelResult = struct('longLimVel', obj.longLimVel, 'longLimPol', obj.longLimPol, ...
                        'shearLimVels', obj.shearLimVels, 'shearLimPols', obj.shearLimPols);
                    if strcmp(method, 'green')
                        if isempty(obj.greenScoreRef)
                            error('CalcSingleHS:identifyModes:noScore', ...
                                'calcGreenScore() must run first (calcAll(''full'') runs it).');
                        end
                        greenResult = struct('greenScore', struct('peaks', obj.greenScorePeaks));
                    else
                        if isempty(obj.refinedSaw)
                            error('CalcSingleHS:identifyModes:noRefined', ...
                                'calcRefined() must run first (classify(''refined'') runs it).');
                        end
                        greenResult = struct('refined', struct('peaks', obj.refinedSaw.peaks));
                    end
                    [obj.classification, obj.limitingEntries] = CalcCommon.classifyGreenPath(obj, ...
                        thr, christoffelResult, greenResult);
                    obj.limitingDone = true;
                    if strcmp(method, 'refined')
                        acc = GreenSpotScore.accepted(obj.refinedSaw.peaks, thr);
                        for j = 1:numel(acc)
                            obj.refinedSaw.peaks(j).accepted = acc(j);
                        end
                        obj.refinedSaw.threshold = thr;
                    end
            end
            obj.activeMethod          = method;
            obj.activeThreshold       = thr;
            obj.classificationVersion = ModeClassification.CLASSIFICATION_VERSION;
        end

        function obj = calcAll(obj, mode, retentionLevel, precomp)
        %%% Run the full (or Ritz-only) pipeline for this angle, then prune
        %%% >> obj = calcAll(obj, mode, retentionLevel)
        %%% >> obj = calcAll(obj, mode, retentionLevel, precomp)
        %%%
        %%% INPUT:  (opt) 1. mode           .. 'full' (default) | 'ritzOnly'
        %%%         (opt) 2. retentionLevel .. 'full' (default) | 'saw_only' | 'none'
        %%%         (opt) 3. precomp        .. shared Ritz precompute, passed to calcRitz()
        %%%                                    (default: this object's own)
        %%%
        %%% 'full':     calcRitz() -> calcLim() -> calcGreen() -> calcGreenScore() ->
        %%%             calcRitzScore() -> identifyModes() ('green', then sawClassMethod
        %%%             if another) -> prune
        %%%             ONE Green solve: calcGreen() on this angle's axis (0.3 km/s to
        %%%             1.05 x longLimVel + 0.1 km/s, greenVelStep steps) gives both the
        %%%             display magnitudes and the score input; calcGreenScore() adds the
        %%%             probe solves (~7 points per candidate each) and the 'low' reference
        %%%             (memoised per structure); calcRitzScore() the Ritz candidates
        %%%             (no solve). 'refined' adds its spot solves (classify).
        %%% 'ritzOnly': calcRitz() -> calcLim() -> calcRitzScore() -> identifyModes()
        %%%             ('ritz') -> prune (no Green; calcLim is milliseconds).
        %%% An unset sawClassMethod runs the mode's own method ('green' / 'ritz'), worked
        %%% out afresh on every run (CalcCommon.resolveSawClassMethod: a numeric
        %%% sawThreshold needs it set; on 'ritzOnly' any method but 'ritz' is an error).
        %%%
        %%% Pruning always happens AFTER identifyModes() returns, never before —
        %%% classification (and the modal-sum weights) need the full eigenvector set.

            if nargin < 2 || isempty(mode);           mode           = obj.calcMode;       end
            if nargin < 3 || isempty(retentionLevel);  retentionLevel = obj.retentionLevel;  end

            CalcCommon.validateCalcMode(mode, 'CalcSingleHS:calcAll');
            CalcCommon.validateRetentionLevel(retentionLevel, 'CalcSingleHS:calcAll', ...
                {'full', 'saw_only', 'none'});

            target = CalcCommon.resolveSawClassMethod(obj.sawClassMethod, obj.sawThreshold, ...
                mode, false, 'CalcSingleHS:calcAll');      % unset: 'green' / 'ritz' per mode
            targetThr = obj.sawThreshold;
            obj.calcMode = mode;

            if nargin < 4; precomp = []; end
            obj = obj.calcRitz(precomp);

            obj.classificationVersion = 0;
            if strcmp(mode, 'full')
                obj = obj.calcLim();
                obj = obj.calcGreen();                         % the one Green solve
                obj = obj.calcGreenScore(obj.greenScoreRef);     % + probe solves
                obj = obj.calcRitzScore(precomp);           % Ritz candidates (no solve)
                %%% 'green' first: it associates the limiting entries while the
                %%% eigenvectors are in memory; then the chosen method, if another
                greenThr = [];
                if strcmp(target, 'green'); greenThr = targetThr; end
                obj = obj.identifyModes('green', greenThr);
                if ~strcmp(target, 'green')
                    obj = obj.classify(target, targetThr);
                end
            else
                obj = obj.calcLim();
                obj = obj.calcRitzScore(precomp);
                obj = obj.identifyModes(target, targetThr);    % 'ritz', the only method here
            end

            obj = obj.applyRetention(retentionLevel);
            obj.retentionLevel = retentionLevel;
        end

        function obj = applyRetention(obj, retentionLevel)
        %%% Apply eigenvector retention policy to ritzSolver (RAM management)
        %%% >> obj = applyRetention(obj, retentionLevel)
        %%%
        %%% 'full':     keep ritzSolver.eigVecs_con untouched.
        %%% 'saw_only': reduce to the columns identifyModes() classified as
        %%%             significant — multi-column in general (e.g. a pseudoRayleigh
        %%%             cluster or limiting modes), not a single sawModeNo.
        %%% 'none':     clear ritzSolver.eigVecs_con entirely (velocities/scores stay).
        %%%
        %%% Public so CalcAngSweepHS can prune each angle's eigenvectors directly.
            obj.ritzSolver = CalcCommon.pruneRitzTo(obj.ritzSolver, obj.classification, ...
                retentionLevel, 'CalcSingleHS:applyRetention');
        end
    end

    methods (Hidden)
        function obj = setRunMode(obj, mode)
        %%% Stamp the calcMode of a sweep's run (CalcCommon.greenStep: its 'full' phases
        %%% call the steps, not calcAll)
        %%% >> obj = setRunMode(obj, mode)
            CalcCommon.validateCalcMode(mode, 'CalcSingleHS:setRunMode');
            obj.calcMode = mode;
        end
    end

    methods (Static, Hidden)
        function T = pairEntries(angle, vA, tA, vB, tB)
        %%% One-to-one nearest pairing of two velocity lists (refinedSawAgreement)
        %%% >> T = CalcSingleHS.pairEntries(angle, vA, tA, vB, tB)
            vA = vA(:); vB = vB(:); tA = tA(:); tB = tB(:);
            D = abs(vB' - vA); tol = max(GreenSpotScore.POLE_WIN_REL * vA, GreenSpotScore.POLE_WIN_ABS);
            D(D > tol) = Inf;
            iA = zeros(0, 1); iB = zeros(0, 1);
            while any(isfinite(D(:)))
                [~, k] = min(D(:)); [a, b] = ind2sub(size(D), k);
                iA(end+1, 1) = a; iB(end+1, 1) = b; %#ok<AGROW>
                D(a, :) = Inf; D(:, b) = Inf;
            end
            oA = setdiff((1:numel(vA))', iA); oB = setdiff((1:numel(vB))', iB);
            n = numel(iA) + numel(oA) + numel(oB);
            dV = [vA(iA); vA(oA); nan(numel(oB), 1)];
            rV = [vB(iB); nan(numel(oA), 1); vB(oB)];
            dT = [tA(iA); tA(oA); repmat({''}, numel(oB), 1)];
            rT = [tB(iB); repmat({''}, numel(oA), 1); tB(oB)];
            st = [repmat({'both'}, numel(iA), 1); repmat({'defaultOnly'}, numel(oA), 1); ...
                repmat({'refinedOnly'}, numel(oB), 1)];
            T = table(repmat(angle, n, 1), dV, rV, rV - dV, dT, rT, st, 'VariableNames', ...
                {'angle', 'defaultVel', 'refinedVel', 'dv', 'defaultType', 'refinedType', 'status'});
        end
    end

    methods (Access = protected)
        function [vS, vTf] = shearRange(obj)
        %%% Slowest and fastest shear limiting velocity of this angle (needs calcLim)
        %%% >> [vS, vTf] = shearRange(obj)
        %%%
        %%% Plain min / max of shearLimVels, whatever the branch: vS is the Rayleigh
        %%% boundary, vTf the candidate cap (CalcCommon.limitingAndBulk's vS / vTf).
            vS = min(obj.shearLimVels); vTf = max(obj.shearLimVels);
        end

        function vBound = scoreAxisBound(obj)
        %%% Bound of this angle's Green axis: 1.05 x longLimVel (needs calcLim)
        %%% >> vBound = scoreAxisBound(obj)
        %%%
        %%% The axis runs from 0.3 km/s to vBound + 0.1 km/s; the candidate search stops at
        %%% longLimVel, and the 5 % beyond keeps a peak at the bound interior and clears the
        %%% end taper.
            vBound = SawSpectrum.AXIS_HI_X_VL * obj.longLimVel;
        end

        function cp = copyElement(obj)
        %%% Deep-copy override (matlab.mixin.Copyable)
        %%% >> cp = copyElement(obj)
        %%% All current properties are value types (Structure, RitzSolverSAW,
        %%% double, char, struct), so this is equivalent to the default shallow
        %%% copy; kept as an explicit override so any future handle-typed property
        %%% is not silently shallow-copied.
            cp = copyElement@matlab.mixin.Copyable(obj);
        end
    end
end
