classdef CalcSingleLayered < handle & matlab.mixin.Copyable
%%% CalcSingleLayered — single-angle container for the layered forward
%%% calculation (film(s)-on-substrate: RitzSolverSAW + substrate-only
%%% Christoffel + GreenSolver).
%%%
%%% Given one Structure (nLay >= 2: at least one film plus a substrate) and
%%% one propagation angle, runs Ritz and Green for that direction, plus a
%%% SUBSTRATE-ONLY Christoffel reference (not a system-wide limiting-velocity
%%% treatment — that concept is undefined for a layered stack, see calcLim()),
%%% and holds the assembled results together. Usable standalone, and the
%%% per-angle building block CalcAngSweepLayered drives across an angle array.
%%%
%%% A SEPARATE class from CalcSingleHS, not a generalisation of it — layered
%%% systems have an arbitrary number of Sezawa modes and no system-wide
%%% Christoffel bound. See docs/reference/architecture/forwardCalcs.md ("Why half-space
%%% and layered are separate class families").
%%%
%%% Handle class: holds potentially large per-object data (Ritz eigenvectors,
%%% Green spectrum) and the GUI viewers need write-back.
%%%
%%% identifyModes() delegates to ModeClassification's recognizeLayered path
%%% ('full' mode, existence by GreenScoreLayered) or its Ritz-only path ('ritzOnly'
%%% mode): a single 'SAW' type per accepted peak, no Rayleigh/Sezawa split, with the
%%% aboveSubstrateST flag reported per entry.
%%%
%%% GuliweRR v6.0

    properties (SetAccess = protected)
        structure       Structure     = Structure()   % rotated structure for this angle;
                                                         % CalcParams lives inside it already
                                                         % (structure.calcParams, one per layer)
        angle(1,1)      double        = 0             % propagation angle [deg]
        structureBase   Structure     = Structure()   % the UNROTATED structure passed in:
                                                        % GreenScoreLayered's probes rotate it
                                                        % to angle +- PROBE_OFFSET

        %%% calcFlags
        calcMode        char          = 'full'         % 'full' (default) | 'ritzOnly'

        %%% ritzSolver -- calcRitz() result
        ritzSolver      RitzSolverSAW = RitzSolverSAW()

        %%% christoffel -- calcLim() result -- SUBSTRATE-ONLY Christoffel, NOT a
        %%% system-wide limiting-velocity treatment (undefined for a layered stack;
        %%% see calcLim()). Same field shapes as CalcSingleHS for reuse in
        %%% ModeClassification, but computed from structure.mater(1) (substrate),
        %%% not mater(nLay) (surface) -- prefixed "sub" throughout so a bare
        %%% longLimVel-style name never invites treating it as a system-wide bound.
        subLongLimVel   double  = []   % scalar -- substrate longitudinal limiting velocity
        subLongLimPol   double  = []   % [3x1]
        subShearLimVels double  = []   % [nShear x 1] -- variable count; wave types 2 and 3 concatenated
        subShearLimPols double  = []   % [3 x nShear]
        subBulkVels     double  = []   % [3x1] = [vL; vT1; vT2], descending -- substrate only
        subBulkPols     double  = []   % [3x3] columns = corresponding polarisations

        %%% Band limits over ALL layers -- calcLim() result, from GreenScoreLayered.bandLimits:
        %%% the edges of the existence score's band [0.3 bandVSmin, 1.05 bandVLmax] and of
        %%% the Green axis. A range only, never a score input.
        bandVSmin       double  = []   % scalar -- slowest shear limiting velocity of any layer
        bandVLmax       double  = []   % scalar -- fastest longitudinal limiting velocity of any layer

        %%% greenSolver -- calcGreen() result, as CalcSingleHS: ONE solve at
        %%% SawSpectrum.ETA = 3e-3, only |G13_eff| kept (single); any other component is
        %%% recomputed on request (greenComponent())
        greenVelAxis    double  = []   % [nVel x 1] -- this angle's own axis (greenAxis()), or
                                        % a sweep's common (wider) one
        greenAbsG13Eff  single  = []   % [nVel x 1] |(G13_+ - G13_-)/2| = |(G13_+ - G31_+)/2|
                                        % on greenVelAxis

        %%% GreenScoreLayered existence criterion -- calcGreen() keeps the complex G13_eff
        %%% on the band, calcGreenScore() scores it (both part of calcAll('full'))
        greenScoreVelAxis double = []   % [nBand x 1] band [0.3 bandVSmin, 1.05 bandVLmax]
                                        % (+- 1.5 steps) of greenVelAxis
        greenScoreG13Eff  double = []   % [nBand x 1] complex G13_eff on greenScoreVelAxis
        greenScorePeaks   struct = struct([])   % scored candidates (GreenScoreLayered.direction,
                                                % threshold-independent); struct([]) until
                                                % calcGreenScore() has run

        %%% Ritz-only existence criterion -- calcRitzScore() result (part of
        %%% calcAll('ritzOnly')); as CalcSingleHS
        ritzScoreNum     double = []   % [nKeep x 1] complex modal-sum weights u3 conj(u1),
                                        % modes up to 1.2 x bandVLmax
        ritzScoreVelAxis double = []   % [nSaw x 1] uniform axis (RitzScore.layeredAxis
                                        % over the band [0.3 bandVSmin, 1.05 bandVLmax])
        ritzScoreRef            = []   % RitzScore.ownReference struct (vRef, floor) of
                                        % this direction
        ritzScoreCands   struct = struct([])   % scored candidates (RitzScore.directionLayered,
                                                % threshold-independent); struct([]) until
                                                % calcRitzScore() has run

        %%% modeClassification -- identifyModes() result
        classification          = []   % ModeClassification once identifyModes() has run
        classificationVersion(1,1) double = 0   % ModeClassification.CLASSIFICATION_VERSION
                                        % this result was classified with (0: unclassified,
                                        % or saved before the stamp existed); see isStale()
        activeMethod     char   = ''   % the method the classification holds: 'green' |
                                        % 'ritz' (classify())
        activeThreshold         = []   % its standard: 'high' | 'low' | numeric
        limitingEntries  struct = struct([])   % the limiting-velocity entries, associated
                                        % once per calculation (they need eigenvectors) and
                                        % carried by every 'green' classification
        limitingDone     logical = false   % limitingEntries hold this calculation's
    end

    properties
        %%% See docs/reference/architecture/forwardCalcs.md's HS-vs-Layered property table.
        retentionLevel  char    = 'full'   % 'full' (default) | 'saw_only' | 'none'

        %%% Classification method calcAll() runs: '' (default: follows each run's calcMode
        %%% -- 'green' on 'full', 'ritz' on 'ritzOnly'; the method in effect is
        %%% activeMethod) | 'green' (GreenScoreLayered on the Green's function, needs
        %%% calcMode 'full') | 'ritz' (the gated Ritz modal sum; the only method on
        %%% 'ritzOnly'). No 'refined' on layered stacks. Set before calcAll(); classify()
        %%% switches the active classification without changing it.
        sawClassMethod  char    = ''
        %%% Its standard: 'low' (sagittal-type waves: Rayleigh, Sezawa families) | 'high'
        %%% (also Love-type waves) | numeric on the method's own scale (needs
        %%% sawClassMethod set explicitly). [] (default) = 'low' for both methods.
        sawThreshold            = []

        %%% Optional upper end of the Green velocity axis: the axis runs to
        %%% max(vMax, 1.05 x bandVLmax) + 0.1 km/s, so a vMax above the band end only
        %%% extends the displayed |G13_eff|; NaN (default) = the band end.
        vMax(1,1)       double  = NaN
    end

    properties (SetAccess = protected, Transient)
        %%% Angle-independent Ritz precompute (nullConstrMat, dKdC_cell, d2KdCdD_cell),
        %%% built once at construction; a sweep can supply a shared one instead.
        %%% Transient: not written by save() (see CalcSingleHS) -- empty after
        %%% load(); read it through ritzPrecompute(), which rebuilds it on first use.
        %%% Empty on a thread-pool worker too, so a sweep passes its shared copy to
        %%% calcRitz()/calcAll() explicitly (see CalcSingleHS.calcRitz).
        precompStructForRitz struct = struct()
    end

    methods
        function obj = CalcSingleLayered(structure, angle, vMax, precompStruct, calcMode, retentionLevel)
        %%% Construct CalcSingleLayered for one propagation angle
        %%% >> obj = CalcSingleLayered(structure, angle, vMax, precompStruct, calcMode, retentionLevel)
        %%%
        %%% INPUT:  1. structure       .. Structure object, UNROTATED (the constructor
        %%%                              rotates it internally — callers never rotate
        %%%                              by hand); must have nLay >= 2 (at least one
        %%%                              film plus a substrate)
        %%%  (opt)  2. angle           .. propagation angle [deg] (default = 0)
        %%%  (opt)  3. vMax            .. Green velocity-axis upper end (the axis runs to
        %%%                              max(vMax, 1.05 x bandVLmax)); [] or NaN
        %%%                              (default) = 1.05 x bandVLmax
        %%%  (opt)  4. precompStruct   .. angle-independent Ritz precompute fields
        %%%                              (nullConstrMat, dKdC_cell, d2KdCdD_cell) — built
        %%%                              from the rotated structure if omitted, so a bare
        %%%                              CalcSingleLayered(structure, angle) is fully usable
        %%%                              standalone
        %%%  (opt)  5. calcMode        .. 'full' (default) | 'ritzOnly'
        %%%  (opt)  6. retentionLevel  .. 'full' (default) | 'saw_only' | 'none'

            if nargin == 0
                return
            end

            validateattributes(structure, {'Structure'}, {'scalar'}, 'CalcSingleLayered', 'structure');
            if structure.nLay < 2
                error('CalcSingleLayered:constructor:needsFilmAndSubstrate', ...
                    'CalcSingleLayered needs at least one film plus a substrate (nLay >= 2); use CalcSingleHS for a homogeneous half-space.');
            end
            if nargin < 2 || isempty(angle); angle = 0; end
            validateattributes(angle, {'double'}, {'scalar', 'real', 'finite'}, 'CalcSingleLayered', 'angle');
            if nargin < 5; calcMode = []; end
            if nargin < 6; retentionLevel = []; end
            CalcCommon.validateConstructorFlags(calcMode, retentionLevel, 'CalcSingleLayered:constructor', {'full', 'saw_only', 'none'});

            obj.angle     = angle;
            obj.structureBase = structure;
            obj.structure = CalcCommon.rotateStructure(structure, angle);

            if nargin >= 3 && ~isempty(vMax)
                obj.vMax = vMax;
            end

            if nargin < 4; precompStruct = struct(); end
            %%% Fill-if-missing, so a shared sweep precompute passes straight
            %%% through and a bare CalcSingleLayered(structure, angle) builds its own.
            obj.precompStructForRitz = CalcCommon.ritzSAWPrecompute(obj.structure, precompStruct);

            if nargin >= 5 && ~isempty(calcMode)
                obj.calcMode = calcMode;
            end
            if nargin >= 6 && ~isempty(retentionLevel)
                obj.retentionLevel = retentionLevel;
            end
        end

        function obj = calcRitz(obj, precomp)
        %%% Run RitzSolverSAW for this angle
        %%% >> obj = calcRitz(obj)
        %%% >> obj = calcRitz(obj, precomp)
        %%%
        %%% INPUT: (opt) precomp .. the sweep's shared Ritz precompute; default: this
        %%%                         object's own (ritzPrecompute()). Same contract as
        %%%                         CalcSingleHS.calcRitz -- a sweep always passes it.
        %%%
        %%% withDiagnostics = true stores the per-mode diagnostics on ritzSolver.diagnostics
        %%% (oopSurfDispl, freqDisp, depthInsensitivity), which identifyModes() reads;
        %%% freqDisp carries no identity information on a layered stack.
        %%%
        %%% Pass 0, not obj.angle — obj.structure is already rotated.
            if nargin < 2 || isempty(precomp)
                precomp = obj.ritzPrecompute();
            else
                precomp = CalcCommon.completeRitzPrecompute(obj.structure, precomp, class(obj));
            end
            rs = RitzSolverSAW(obj.structure, 0);
            obj.ritzSolver = rs.calc(true, precomp);
            obj.limitingEntries = struct([]); obj.limitingDone = false;   % new eigenvectors
        end

        function precomp = ritzPrecompute(obj)
        %%% The Ritz precompute, rebuilt first if missing (e.g. after load())
        %%% >> precomp = ritzPrecompute(obj)
        %%%
        %%% OUTPUT: precomp .. precompStructForRitz with every field populated
        %%%
        %%% Same fill-if-missing contract as CalcSingleHS.ritzPrecompute.
            %%% No own "is it complete?" test: ritzSAWPrecompute checks every field
            %%% itself, so a field added there later is covered here automatically.
            obj.precompStructForRitz = CalcCommon.completeRitzPrecompute(obj.structure, ...
                obj.precompStructForRitz, class(obj));
            precomp = obj.precompStructForRitz;
        end

        function obj = calcLim(obj)
        %%% Substrate-only Christoffel limiting velocities and bulk velocities
        %%% >> obj = calcLim(obj)
        %%%
        %%% There is no layered analogue of CalcSingleHS.calcLim()'s
        %%% longLimVel/shearLimVels. What IS well-defined is the SUBSTRATE's own
        %%% limiting velocities, computed exactly as CalcSingleHS.calcLim() does
        %%% but anchored to structure.mater(1) (the substrate), NOT mater(nLay).
        %%%
        %%% Downstream this is a SOFT reference velocity: the substrate shear
        %%% limit is the threshold above which Sezawa modes are flagged leaky
        %%% (aboveSubstrateST), not a hard bound on vMax. Leaky/pseudo-Sezawa
        %%% modes do exist above it.
        %%%
        %%% Also fills bandVSmin / bandVLmax: the slowest shear and fastest longitudinal
        %%% limiting velocity over ALL layers, the edges of the GreenScoreLayered band and
        %%% of the Green axis (no layer singled out there).
        %%%
        %%% CAVEAT, unresolved — do not extend silently: the substrate-only
        %%% reference is well-motivated for the classic single-film-on-substrate
        %%% case. For a stack with more than one film layer it is UNTESTED whether
        %%% the substrate alone remains the right reference, whether an
        %%% intermediate film's shear velocity matters too, or whether the concept
        %%% needs rethinking.

            subIdx = 1;   %%% substrate = bottom of the stack (Structure.m convention),
                          %%% NOT nLay -- that is CalcSingleHS's convention, correct
                          %%% there because the whole half-space is one material
            %%% every layer's limits once: layer subIdx fills the substrate fields below, all
            %%% of them give the band edges (bandLimits would otherwise redo the substrate's)
            lims = cell(1, obj.structure.nLay);
            for k = 1:obj.structure.nLay
                lims{k} = CalcCommon.limitingAndBulk(obj.structure.mater(k));
            end
            lim = lims{subIdx};

            obj.subLongLimVel   = lim.longLimVel;
            obj.subLongLimPol   = lim.longLimPol;
            obj.subShearLimVels = lim.shearLimVels;
            obj.subShearLimPols = lim.shearLimPols;
            obj.subBulkVels     = lim.bulkVels;
            obj.subBulkPols     = lim.bulkPols;

            %%% The existence score's band uses every layer, not the substrate alone
            [obj.bandVSmin, obj.bandVLmax] = GreenScoreLayered.bandLimits(obj.structure, lims);
        end

        function velAxis = greenAxis(obj)
        %%% This angle's own Green axis (needs calcLim)
        %%% >> velAxis = greenAxis(obj)
        %%%
        %%% GreenScoreLayered.gridAxis(0.3 x bandVSmin, max(vMax, 1.05 x bandVLmax)): 4 m/s
        %%% steps on the lattice SawSpectrum.VEL_MIN + n x 4 m/s, from 0.3 km/s (or lower, for a
        %%% stack whose slowest shear limit is below 1 km/s) to the bound + 0.1 km/s. A
        %%% sweep's common axis lies on the same lattice and covers every angle's band, so
        %%% the band a sweep angle scores on is bit for bit this one.
            if isempty(obj.bandVSmin)
                error('CalcSingleLayered:greenAxis:noLim', ...
                    'calcLim() must run first: the axis runs to 1.05 x the fastest longitudinal limit of any layer.');
            end
            [vLo, ~] = GreenScoreLayered.band(obj.bandVSmin, obj.bandVLmax);
            velAxis = GreenScoreLayered.gridAxis(vLo, obj.scoreAxisBound());
        end

        function obj = calcGreen(obj, velAxis)
        %%% The one Green solve of this angle: display magnitudes + the existence-score input
        %%% >> obj = calcGreen(obj)
        %%% >> obj = calcGreen(obj, velAxis)
        %%%
        %%% INPUT: (opt) velAxis .. phase-velocity axis [km/s]. Default: this angle's own
        %%%                         axis, greenAxis() (needs calcLim() first). A sweep
        %%%                         passes its common axis.
        %%%
        %%% One CalcCommon.greenScoreSpectrum call (SawSpectrum.ETA, |G13_eff| kept). When
        %%% calcLim() has run and velAxis is uniform over the whole GreenScoreLayered band,
        %%% the complex G13_eff on the band is kept too (greenScoreVelAxis/greenScoreG13Eff);
        %%% otherwise calcGreenScore() solves the band itself. Clears any previous scores.
            if nargin < 2 || isempty(velAxis)
                velAxis = obj.greenAxis();
            end
            g = CalcCommon.greenScoreSpectrum(obj.structure, velAxis, 'CalcSingleLayered.calcGreen');

            obj.greenVelAxis   = g.velAxis;
            obj.greenAbsG13Eff = g.absG13Eff;
            obj.greenScoreVelAxis = []; obj.greenScoreG13Eff = []; obj.greenScorePeaks = struct([]);
            if isempty(obj.bandVSmin); return; end

            [vLo, vHi] = GreenScoreLayered.band(obj.bandVSmin, obj.bandVLmax);
            v = g.velAxis; dv = v(2) - v(1);
            b = v >= vLo - 1.5 * dv & v <= vHi + 1.5 * dv;
            covers = v(1) <= vLo && v(end) >= vHi;
            if covers && nnz(b) >= 5 && max(abs(diff(v(b)) - dv)) < 1e-9
                obj.greenScoreVelAxis = v(b);
                obj.greenScoreG13Eff  = g.G13Eff(b);
            end
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

        function obj = calcGreenScore(obj)
        %%% Score this angle's candidates (GreenScoreLayered): probes, polarisation, scores
        %%% >> obj = calcGreenScore(obj)
        %%%
        %%% Requires calcLim(). Cost beyond calcGreen(): two probe solves of ~7 points per
        %%% candidate (angle +- GreenScore.PROBE_OFFSET) and one full-tensor solve of 3
        %%% points per candidate. No reference: every score is normalised by this
        %%% direction's own band median, so a single angle and the same angle in any sweep
        %%% score identically. If calcGreen() left no band (not run, or a custom axis not
        %%% covering the band), the band is solved here.
            if isempty(obj.bandVSmin)
                error('CalcSingleLayered:calcGreenScore:noLim', 'calcLim() must run first.');
            end
            if isempty(obj.greenScoreG13Eff)
                [vLo, vHi] = GreenScoreLayered.band(obj.bandVSmin, obj.bandVLmax);
                v = GreenScoreLayered.gridAxis(vLo, vHi); dv = v(2) - v(1);
                v = v(v >= vLo - 1.5 * dv & v <= vHi + 1.5 * dv);
                g = CalcCommon.greenScoreSpectrum(obj.structure, v, 'CalcSingleLayered.calcGreenScore');
                obj.greenScoreVelAxis = g.velAxis;
                obj.greenScoreG13Eff  = g.G13Eff;
            end
            obj.greenScorePeaks = GreenScoreLayered.direction(obj.structureBase, obj.angle, ...
                obj.greenScoreG13Eff, obj.greenScoreVelAxis, obj.bandVSmin, obj.bandVLmax);
        end

        function obj = calcRitzScore(obj, precomp)
        %%% Ritz-only existence input: modal-sum weights, score axis, this direction's own
        %%% floor reference and the scored layered candidates (RitzScore layered rule)
        %%% >> obj = calcRitzScore(obj)
        %%% >> obj = calcRitzScore(obj, precomp)
        %%%
        %%% INPUT: (opt) precomp .. the Ritz precompute calcRitz() used (its dKdC_cell gives
        %%%                         the modes' angular slopes); default: ritzPrecompute()
        %%%
        %%% Requires calcRitz() (eigenvectors present, i.e. before applyRetention) and
        %%% calcLim() (bandVSmin / bandVLmax). No vMax needed. The axis, the kept modes, the
        %%% candidate search and the floor all span GreenScoreLayered's band [0.3 bandVSmin,
        %%% 1.05 bandVLmax] -- no layer singled out; the substrate's shear limit is used
        %%% only for the aboveSubstrateST label. Fills ritzScoreNum/VelAxis/Ref and
        %%% ritzScoreCands (RitzScore.directionLayered: everything a reclassify() needs).
            if isempty(obj.bandVLmax)
                error('CalcSingleLayered:calcRitzScore:noLim', 'calcLim() must run first.');
            end
            if nargin < 2 || isempty(precomp)
                precomp = obj.ritzPrecompute();
            end
            dKdC = [];   % a partial precompute: angleSlope builds it (as calcRitz would)
            if isfield(precomp, 'dKdC_cell'); dKdC = precomp.dKdC_cell; end
            [vLo, vHi] = GreenScoreLayered.band(obj.bandVSmin, obj.bandVLmax);
            [num, ~, ~, Ebot, SH] = RitzScore.modeProfiles(obj.ritzSolver, obj.bandVLmax);
            obj.ritzScoreNum     = num;
            obj.ritzScoreVelAxis = RitzScore.layeredAxis(vLo, vHi);
            obj.ritzScoreRef     = RitzScore.ownReference(obj.structure, ...
                obj.ritzSolver.vels(1:numel(num)), num, obj.ritzScoreVelAxis, vLo, vHi);
            obj.ritzScoreCands   = RitzScore.directionLayered(obj.ritzSolver, num, Ebot, SH, ...
                obj.ritzScoreVelAxis, obj.ritzScoreRef, vLo, vHi, min(obj.subShearLimVels), dKdC);
        end

        function obj = classify(obj, method, threshold)
        %%% Set the active classification: method and standard (re-thresholds stored candidates)
        %%% >> obj = classify(obj)
        %%% >> obj = classify(obj, method)
        %%% >> obj = classify(obj, method, threshold)
        %%%
        %%% INPUT: (opt) method    .. 'green' | 'ritz' (default: the active one); 'refined'
        %%%                           is an error (no spot score for layered stacks)
        %%%        (opt) threshold .. 'low' | 'high' | numeric on the method's own scale
        %%%                           (default 'low')
        %%% OUTPUT: obj with classification, activeMethod, activeThreshold and
        %%%         classificationVersion set
        %%%
        %%% As CalcSingleHS.classify: both methods' candidates are stored, so switching only
        %%% re-thresholds; 'green' needs calcMode 'full'; a stale result is refused.
            if nargin < 3; threshold = []; end
            if nargin < 2 || isempty(method); method = obj.activeMethod; end
            if isempty(method); method = obj.sawClassMethod; end
            if isempty(method)   % neither classified nor set: calcMode's own
                method = CalcCommon.resolveSawClassMethod('', threshold, obj.calcMode, ...
                    true, 'CalcSingleLayered:classify');
            end
            method = CalcCommon.validateSawClassMethod(method, true, 'CalcSingleLayered:classify');
            if isempty(threshold)
                threshold = CalcCommon.defaultThreshold(method, true);
            end
            CalcCommon.validateSawThreshold(threshold, 'CalcSingleLayered:classify');
            CalcCommon.refuseStale(obj, 'CalcSingleLayered:classify');
            if strcmp(method, 'green') && ~strcmp(obj.calcMode, 'full')
                error('CalcSingleLayered:classify:needsFull', ...
                    "Method 'green' needs calcMode='full'; this result is '%s'.", obj.calcMode);
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

        function obj = identifyModes(obj, method, thr)
        %%% Classify with a method at a standard (ModeClassification, layered paths)
        %%% >> obj = identifyModes(obj)                 % sawClassMethod at sawThreshold
        %%% >> obj = identifyModes(obj, method, thr)
        %%%
        %%% method empty: sawClassMethod, or calcMode's own; thr empty = 'low'.
        %%% 'green' .. GreenScoreLayered candidates (recognizeLayered; needs calcGreenScore());
        %%% 'ritz'  .. the gated Ritz modal sum (needs calcLim() and calcRitzScore()).
        %%% The limiting entries are associated once per calculation (they need the
        %%% eigenvectors) and carried by every later 'green' classification; 'ritz' has
        %%% none, so it equals a 'ritzOnly' result's. subLongLimPol is forwarded so
        %%% recognizeLayered can cross-link the substrate longitudinal limiting velocity
        %%% to a Ritz mode.
            if nargin < 3; thr = obj.sawThreshold; end
            if nargin < 2 || isempty(method); method = obj.sawClassMethod; end
            if isempty(method)
                method = CalcCommon.resolveSawClassMethod('', thr, obj.calcMode, true, ...
                    'CalcSingleLayered:identifyModes');
            end
            method = CalcCommon.validateSawClassMethod(method, true, 'CalcSingleLayered:identifyModes');
            if isempty(thr); thr = CalcCommon.defaultThreshold(method, true); end
            if strcmp(method, 'ritz')
                if isempty(obj.ritzScoreVelAxis) || isempty(obj.subLongLimVel)
                    error('CalcSingleLayered:identifyModes:noRitzScore', ['Ritz classification ' ...
                        'needs calcLim() and calcRitzScore() first (calcAll runs both).']);
                end
                mc = ModeClassification('threshold', thr);
                obj.classification = mc.calc(obj.ritzSolver, obj.structure, [], [], ...
                    struct('cands', obj.ritzScoreCands));
            else
                if ~isfield(obj.greenScorePeaks, 'score')
                    error('CalcSingleLayered:identifyModes:noScore', ...
                        'calcGreenScore() must run first (calcAll(''full'') runs it).');
                end
                christoffelResult = struct('subLongLimVel', obj.subLongLimVel, 'subLongLimPol', obj.subLongLimPol, ...
                    'subShearLimVels', obj.subShearLimVels, 'subShearLimPols', obj.subShearLimPols);
                greenResult = struct('greenScoreLayered', struct('peaks', obj.greenScorePeaks));
                [obj.classification, obj.limitingEntries] = CalcCommon.classifyGreenPath(obj, ...
                    thr, christoffelResult, greenResult);
                obj.limitingDone = true;
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
        %%%             calcRitzScore() -> identifyModes() ('green', then 'ritz' if chosen) -> prune
        %%%             ONE axis Green solve, on this angle's own axis (greenAxis(): 0.3 km/s,
        %%%             or 0.3 x bandVSmin if lower, to max(vMax, 1.05 x bandVLmax) +
        %%%             0.1 km/s, 4 m/s steps); calcGreenScore() adds the probe and residue
        %%%             solves (a few points per candidate).
        %%% 'ritzOnly': calcRitz() -> calcLim() -> calcRitzScore() -> identifyModes() ('ritz')
        %%%             -> prune (substrate-only Christoffel, no Green, no vMax).
        %%%
        %%% An unset sawClassMethod runs the mode's own method ('green' / 'ritz'), worked
        %%% out afresh on every run (CalcCommon.resolveSawClassMethod).
        %%%
        %%% Pruning always happens AFTER identifyModes() returns, never before.

            if nargin < 2 || isempty(mode);           mode           = obj.calcMode;       end
            if nargin < 3 || isempty(retentionLevel);  retentionLevel = obj.retentionLevel;  end

            CalcCommon.validateCalcMode(mode, 'CalcSingleLayered:calcAll');
            CalcCommon.validateRetentionLevel(retentionLevel, 'CalcSingleLayered:calcAll', ...
                {'full', 'saw_only', 'none'});
            target = CalcCommon.resolveSawClassMethod(obj.sawClassMethod, obj.sawThreshold, ...
                mode, true, 'CalcSingleLayered:calcAll');  % unset: 'green' / 'ritz' per mode
            targetThr = obj.sawThreshold;
            obj.calcMode = mode;

            if nargin < 4; precomp = []; end
            obj = obj.calcRitz(precomp);

            obj.classificationVersion = 0;
            if strcmp(mode, 'full')
                obj = obj.calcLim();
                obj = obj.calcGreen();                     % the one axis Green solve
                obj = obj.calcGreenScore();                 % probes + residues, no reference
                obj = obj.calcRitzScore(precomp);       % Ritz candidates (no solve)
                %%% 'green' first (associates the limiting entries), then 'ritz' if chosen
                greenThr = [];
                if strcmp(target, 'green'); greenThr = targetThr; end
                obj = obj.identifyModes('green', greenThr);
                if ~strcmp(target, 'green')
                    obj = obj.classify(target, targetThr);
                end
            else
                obj = obj.calcLim();
                obj = obj.calcRitzScore(precomp);
                obj = obj.identifyModes(target, targetThr);   % 'ritz', the only method here
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
        %%%             significant -- multi-column in general (any number of
        %%%             'SAW' entries plus limiting entries).
        %%% 'none':     clear ritzSolver.eigVecs_con entirely (velocities/scores stay).
        %%%
        %%% Public so CalcAngSweepLayered can prune each angle's eigenvectors
        %%% directly between phases.
            obj.ritzSolver = CalcCommon.pruneRitzTo(obj.ritzSolver, obj.classification, ...
                retentionLevel, 'CalcSingleLayered:applyRetention');
        end
    end

    methods (Hidden)
        function obj = setRunMode(obj, mode)
        %%% Stamp the calcMode of a sweep's run (CalcCommon.greenStep: its 'full' phases
        %%% call the steps, not calcAll)
        %%% >> obj = setRunMode(obj, mode)
            CalcCommon.validateCalcMode(mode, 'CalcSingleLayered:setRunMode');
            obj.calcMode = mode;
        end
    end

    methods (Access = protected)
        function vBound = scoreAxisBound(obj)
        %%% Bound of this angle's Green axis: max(vMax, 1.05 x bandVLmax)
        %%% >> vBound = scoreAxisBound(obj)
        %%%
        %%% The candidate search stops at the band end, 1.05 x the fastest longitudinal
        %%% limiting velocity of any layer; vMax only extends the (display) axis beyond
        %%% it. max() ignores an unset (NaN) vMax.
            vBound = max(obj.vMax, GreenScoreLayered.BAND_HI_X_VLMAX * obj.bandVLmax);
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
