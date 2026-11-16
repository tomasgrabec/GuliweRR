classdef CalcWavevecAngSweep2Plate < handle & matlab.mixin.Copyable
%%% CalcWavevecAngSweep2Plate — 2D sweep (propagation angle AND real
%%% wavevector) for the free-plate (Lamb-wave) forward calculation. Built on
%%% top of CalcAngSweepPlate: constructs one CalcAngSweepPlate per wavevector
%%% value (at one FIXED realThickness) and drives it through its own calc().
%%% Mirrors CalcWavevecAngSweep2Layered's per-wavelength composition.
%%%
%%% Simpler than the Layered case in two ways, both following from simpler
%%% physics (see CalcSinglePlate.m):
%%%   1. No Structure.ApplyWavelength step -- RitzSolverLamb accepts
%%%      (realThickness, realWavVec) directly and rescales relDepth internally
%%%      per point. The wavevector axis here is RitzSolverLamb's own
%%%      realWavVec, swept at one fixed realThickness.
%%%   2. No Green/Christoffel phase barrier -- every (angle, wavevector) point
%%%      is an independent Ritz solve, so there is no flattenSweep/calcFlat
%%%      counterpart to the Layered class. The outer wavevector loop stays
%%%      serial and each CalcAngSweepPlate parallelises its own angle loop.
%%%
%%% rawByDepthCell (the Legendre integrals, legOrder only) is built once per
%%% wavevector by each CalcAngSweepPlate and shared by its angles (precompute
%%% field of RitzSolverLamb.calc); sharing it across the wavevector axis too
%%% would save nWav-1 more builds (the Julia port does, its grid being flat).
%%%
%%% GuliweRR v6.0

    properties (SetAccess = protected)
        structure       Structure   = Structure()      % REFERENCE structure, UNROTATED.
        angles          double      = (0:2:180)'          % [nAngles x 1] propagation angles [deg]
        wavVecs         double      = []                 % [nWav x 1] positive real wavevectors --
                                                           % REQUIRED, no default
        realThickness   double      = []                 % physical plate thickness, FIXED across
                                                           % the whole sweep -- REQUIRED. Without it
                                                           % realWavVec has no effect on relDepth,
                                                           % making a wavevector sweep meaningless
        retentionLevel  char        = 'none'              % 'full' | 'none' (default = 'none', RAM
                                                           % multiplying across both axes). Pass
                                                           % 'full' when eigenvectors are needed,
                                                           % e.g. for
                                                           % DispersionBranches.track
        useParallel(1,1) string {mustBeMember(useParallel, ["auto","always","never"])} = "auto"
                                                           % forwarded unchanged to every
                                                           % sweeps{iWav}, which decides at the angle
                                                           % level; this class's outer iWav loop
                                                           % stays serial

        %%% Per-wavevector containers -- one CalcAngSweepPlate per wavevector
        sweeps          cell        = {}                 % {nWav x 1} of CalcAngSweepPlate objects

        %%% Convenience 3D extraction -- filled by assembleMaps() at the end of calc().
        %%% Dense, not ragged/cell as CalcWavevecAngSweep2Layered's maps are, because
        %%% nCompModes is fixed.
        velsMap3D               double = []   % [nModes x nAngles x nWav]
        velsGroupMap3D          double = []   % [nModes x nAngles x nWav] group velocity along k
        velsGroupTransMap3D     double = []   % [nModes x nAngles x nWav] transverse group velocity
                                               % (perpendicular to k), RitzSolverLamb.velsGroupTrans
        modesTypeMap3D          double = []   % [nModes x nAngles x nWav]
        modeTypeConfidenceMap3D double = []   % [nModes x nAngles x nWav]
        freqsPhysicalMap3D      double = []   % [nModes x nAngles x nWav] -- always populated here,
                                                % wavVecs being required positive
    end

    properties (Constant)
        %%% locateZGVs acceptance. ZGV_VG_REL_TOL: largest |vg|/vp accepted as a ZGV, checked
        %%% on an independent Ritz solve at the Newton point (a genuine root sits at
        %%% |vg|/vp ~1e-12 or below); the same bound, on |c_gtheta|/vp, decides isTrueZGV
        %%% (a TGV wave off the symmetry directions sits at ~1e-2, Si(001) at 18.4 deg).
        %%% Near a symmetry direction |c_gtheta|/vp grows linearly with the angle offset, so a
        %%% point a little off it still counts as a true ZGV: on Si(001) the slope is 1e-3 per
        %%% degree off [100] and 7e-4 off [110], i.e. within 0.10 deg / 0.14 deg.
        %%% ZGV_DUP_REL_TOL: two Newton results closer than this in k and in f (relative)
        %%% are one point (two seeds converging onto the same ZGV).
        %%% ZGV_FREQ_MATCH_REL_TOL: largest relative frequency mismatch between the Newton
        %%% point and the mode of the independent solve it is matched to.
        ZGV_VG_REL_TOL         = 1e-4
        ZGV_DUP_REL_TOL        = 1e-8
        ZGV_FREQ_MATCH_REL_TOL = 1e-6
        %%% ZGV_DIP_REL_TOL: a local minimum of |vg| along a tracked branch below this
        %%% fraction of vp, with no sign change next to it, also seeds Newton (two ZGV
        %%% points between two grid steps; see zgvCandidates). Chosen from one case, a
        %%% three-layer stack whose dip sits at |vg|/vp = 1.1e-3; revisit when another case
        %%% misses a pair of points.
        ZGV_DIP_REL_TOL        = 0.01
    end

    methods
        function obj = CalcWavevecAngSweep2Plate(structure, angles, wavVecs, realThickness, retentionLevel, useParallel)
        %%% Construct CalcWavevecAngSweep2Plate for a 2D (angle x wavevector) sweep
        %%% >> obj = CalcWavevecAngSweep2Plate(structure, angles, wavVecs, realThickness, retentionLevel, useParallel)
        %%%
        %%% INPUT:  1. structure       .. Structure object, UNROTATED, nLay >= 1
        %%%  (opt)  2. angles          .. propagation angles [deg] (default = (0:2:180)')
        %%%         3. wavVecs         .. [nWav x 1] positive real wavevectors, strictly
        %%%                              increasing -- REQUIRED, no default
        %%%         4. realThickness   .. physical plate thickness, FIXED across the whole
        %%%                              sweep -- REQUIRED
        %%%  (opt)  5. retentionLevel  .. 'full' | 'none' (default = 'none')
        %%%  (opt)  6. useParallel     .. "auto" (default) | "always" | "never"

            if nargin == 0
                return
            end

            validateattributes(structure, {'Structure'}, {'scalar'}, 'CalcWavevecAngSweep2Plate', 'structure');

            if nargin < 2 || isempty(angles); angles = (0:2:180)'; end
            validateattributes(angles, {'double'}, {'vector', 'real', 'finite'}, ...
                'CalcWavevecAngSweep2Plate', 'angles');
            if nargin >= 5 && ~isempty(retentionLevel)
                CalcCommon.validateRetentionLevel(retentionLevel, 'CalcWavevecAngSweep2Plate:constructor', {'full', 'none'});
            end

            if nargin < 3 || isempty(wavVecs)
                error('CalcWavevecAngSweep2Plate:constructor:missingWavVecs', ...
                    'wavVecs is required (no default) -- pass a vector of positive real wavevectors.');
            end
            validateattributes(wavVecs, {'double'}, {'vector', 'real', 'finite', 'positive'}, ...
                'CalcWavevecAngSweep2Plate', 'wavVecs');
            if any(diff(wavVecs) <= 0)   % locateZGVs brackets each ZGV between neighbouring grid points
                error('CalcWavevecAngSweep2Plate:constructor:wavVecsNotIncreasing', ...
                    'wavVecs must be strictly increasing (sort them and drop duplicates).');
            end

            if nargin < 4 || isempty(realThickness)
                error('CalcWavevecAngSweep2Plate:constructor:missingRealThickness', ...
                    'realThickness is required (no default) -- without it, wavVecs has no effect on relDepth (see RitzSolverLamb''s own constructor), making a wavevector sweep meaningless.');
            end

            obj.structure      = structure;      % UNROTATED -- kept as-is on this object
            obj.angles         = angles(:);
            obj.wavVecs        = wavVecs(:);
            obj.realThickness  = realThickness;

            if nargin >= 5 && ~isempty(retentionLevel); obj.retentionLevel = retentionLevel; end
            if nargin >= 6 && ~isempty(useParallel);    obj.useParallel    = useParallel;    end
        end

        function obj = calc(obj, retentionLevel)
        %%% Run the 2D sweep: one CalcAngSweepPlate per wavevector, then
        %%% assemble the 3D convenience maps.
        %%% >> obj = calc(obj, retentionLevel)
        %%%
        %%% INPUT: (opt) retentionLevel .. 'full' | 'none' (default = obj.retentionLevel,
        %%%                                'none' unless overridden at construction)
        %%%
        %%% Each iWav iteration runs CalcAngSweepPlate's own Ritz-only loop across every
        %%% angle; this class adds one outer loop over wavevector and nothing new at the
        %%% physics level. The outer loop stays serial (MATLAB forbids nested parfor);
        %%% obj.useParallel is forwarded to each sweeps{iWav} instead.

            if nargin < 2 || isempty(retentionLevel); retentionLevel = obj.retentionLevel; end
            CalcCommon.validateRetentionLevel(retentionLevel, 'CalcWavevecAngSweep2Plate:calc', ...
                {'full', 'none'});   %%% no 'saw_only' tier -- see CalcSinglePlate.m
            obj.retentionLevel = retentionLevel;

            nWav = numel(obj.wavVecs);
            obj.sweeps = cell(nWav, 1);
            for iWav = 1:nWav
                obj.sweeps{iWav} = Calc.angSweepPlate(obj.structure, obj.angles, ...
                    obj.realThickness, obj.wavVecs(iWav), retentionLevel, obj.useParallel);
                obj.sweeps{iWav} = obj.sweeps{iWav}.calc(retentionLevel);
            end

            obj = obj.assembleMaps();
        end

        function zgvList = locateZGVs(obj, angleIdx, maxCount)
        %%% Locate zero-group-velocity (ZGV) points along the wavevector axis at one
        %%% fixed angle: candidate scan, Newton refinement, ZGV/TGV flag.
        %%% >> zgvList = locateZGVs(obj, angleIdx, maxCount)
        %%%
        %%% A located point has zero group velocity ALONG k (velsGroup). Whether it is a
        %%% true ZGV point (a minimum or saddle of the dispersion surface omega(kx, ky),
        %%% the whole group-velocity vector zero) depends on the component perpendicular
        %%% to k, c_gtheta = velsGroupTrans: it vanishes on a direction with a sagittal
        %%% mirror (e.g. <100>, <110> of a cubic (001) plate) and generally not elsewhere,
        %%% where the point is a TGV wave that carries energy sideways (Kiefer, Mezil &
        %%% Prada 2023, Sci. Adv. 9, eadk6846). Both kinds are returned; isTrueZGV tells
        %%% them apart. Validated on a Si(001) wafer: tutorial_S2, Part 5.
        %%%
        %%% INPUT:  (opt) 1. angleIdx .. which angle to search at (default 1)
        %%%         (opt) 2. maxCount .. cap on the number of ZGV points
        %%%                             returned (default = Inf, no cap)
        %%% OUTPUT: zgvList .. struct array, one entry per located point:
        %%%                    branchNo      -- tracked-branch index (column of
        %%%                                     DispersionBranches.track's output, NOT a
        %%%                                     raw Ritz mode number)
        %%%                    wavVecZGV     -- refined wavevector
        %%%                    freqZGV       -- refined physical frequency
        %%%                    angle         -- obj.angles(angleIdx) [deg]
        %%%                    velGroupTrans -- c_gtheta there [same unit as the phase
        %%%                                     velocity], positive towards increasing
        %%%                                     angle
        %%%                    isTrueZGV     -- |velGroupTrans| <= ZGV_VG_REL_TOL * vp;
        %%%                                     false = a TGV wave
        %%%
        %%% STAGE 1 (candidate scan): reuses velsGroupMap3D, already a byproduct of the
        %%% Ritz solve at every point -- no extra forward solves. Branch identity across
        %%% the wavevector axis comes from DispersionBranches.track, which needs
        %%% retentionLevel='full' and errors otherwise. A candidate is any wavevector step
        %%% where a tracked branch's group velocity changes sign, and any dip of |vg| without
        %%% a sign change (two ZGV points between two grid steps; zgvCandidates).
        %%%
        %%% STAGE 2 (Newton refinement, refineZGV): the first seed is the linear interpolation
        %%% of the sign change (k0 and f0 along the tracked branch); when Newton does not
        %%% converge from it, or converges outside the bracket widened by one step on each
        %%% side (onto a neighbouring point), the two bracket ends are tried in turn, each
        %%% from the sweep's own eigenvector of the branch's raw mode there. Where two modes
        %%% nearly coincide in frequency the interpolated seed can pick the neighbour's
        %%% eigenvector; the stored one cannot (a three-layer stack: 3 of 82 points were lost
        %%% without the retries). ZgvNewton.solve then
        %%% solves the plate eigenproblem written as a polynomial in k
        %%% (RitzSolverLamb.kPolynomial, built once per call) together with the
        %%% zero-group-velocity condition, for k, omega and the mode shape at once
        %%% (adapted from GEWtool's ZGVNewtonBeta, Kiefer et al. 2023, JASA 153, 1386).
        %%% No raw mode number is held fixed, so a branch crossing inside the bracket does
        %%% not matter. About 3-6 Newton steps, each one linear solve of the size of the
        %%% Ritz basis, plus one Ritz solve to check the point.
        %%%
        %%% ACCEPTANCE: Newton converged, the point lies inside the swept wavevector range,
        %%% and an independent RitzSolverLamb solve there has a mode at that frequency
        %%% (ZGV_FREQ_MATCH_REL_TOL) with |vg| <= ZGV_VG_REL_TOL * vp. Two candidates that
        %%% converge onto the same point are reported once. branchNo is the candidate's
        %%% branch when the point lies inside its bracket, else the tracked branch nearest
        %%% in frequency in the step containing the point. maxCount counts accepted points,
        %%% in candidate order (branch, then wavevector step).
        %%%
        %%% A candidate converging nowhere is dropped (no point is invented). Off-axis
        %%% stationary points of omega(kx, ky) (a true ZGV off the symmetry directions of a
        %%% low-symmetry plate) need a 2-D search in (k, angle), not done here.

            if nargin < 2 || isempty(angleIdx); angleIdx = 1; end
            if nargin < 3 || isempty(maxCount); maxCount = Inf; end
            validateattributes(angleIdx, {'double'}, ...
                {'scalar', 'integer', 'positive', '<=', numel(obj.angles)}, ...
                'CalcWavevecAngSweep2Plate.locateZGVs', 'angleIdx');

            branchModeNoMap = DispersionBranches.track(obj, angleIdx);
            nCompModes = size(branchModeNoMap, 2);

            velsGroupByBranch = CalcWavevecAngSweep2Plate.reindexByBranch( ...
                obj.velsGroupMap3D, angleIdx, branchModeNoMap);
            freqByBranch = CalcWavevecAngSweep2Plate.reindexByBranch( ...
                obj.freqsPhysicalMap3D, angleIdx, branchModeNoMap);
            velsByBranch = CalcWavevecAngSweep2Plate.reindexByBranch( ...
                obj.velsMap3D, angleIdx, branchModeNoMap);

            zgvList = struct('branchNo', {}, 'wavVecZGV', {}, 'freqZGV', {}, 'angle', {}, ...
                'velGroupTrans', {}, 'isTrueZGV', {});
            if maxCount <= 0
                return
            end

            %%% via cs{1}.ritzPrecompute(), not sweeps{1}.precompStructForRitz: that
            %%% one is Transient and empty on a load()'d result. nullConstrMat and
            %%% rawByDepthCell depend on calcParams (legOrder) only.
            precomp0      = obj.sweeps{1}.cs{1}.ritzPrecompute();
            precompStruct = struct('nullConstrMat', precomp0.nullConstrMat);
            if isfield(precomp0, 'rawByDepthCell') && ~isempty(precomp0.rawByDepthCell)
                precompStruct.rawByDepthCell = precomp0.rawByDepthCell;
            else
                precompStruct.rawByDepthCell = GammaMatrix.RawCellFor(obj.structure);
            end
            poly = RitzSolverLamb.kPolynomial(obj.structure, obj.angles(angleIdx), ...
                obj.realThickness, precompStruct);

            k = obj.wavVecs;
            for ib = 1:nCompModes
                cands = CalcWavevecAngSweep2Plate.zgvCandidates(k, velsGroupByBranch(:, ib), ...
                    freqByBranch(:, ib), velsByBranch(:, ib), branchModeNoMap(:, ib));
                for ic = 1:numel(cands)
                    if numel(zgvList) >= maxCount
                        break
                    end
                    c = cands(ic);
                    [wavVecZGV, freqZGV, isRoot, velGroupTrans, velPhase] = ...
                        obj.refineZGV(angleIdx, poly, c.seeds, c.kWin, precompStruct);
                    if ~isRoot
                        continue
                    end
                    isDup = false;
                    for iz = 1:numel(zgvList)
                        isDup = isDup || ...
                            (abs(zgvList(iz).wavVecZGV - wavVecZGV) <= CalcWavevecAngSweep2Plate.ZGV_DUP_REL_TOL * wavVecZGV && ...
                             abs(zgvList(iz).freqZGV - freqZGV) <= CalcWavevecAngSweep2Plate.ZGV_DUP_REL_TOL * freqZGV);
                    end
                    if isDup
                        continue
                    end
                    if wavVecZGV >= c.kBracket(1) && wavVecZGV <= c.kBracket(2)
                        branchNo = ib;
                    else
                        branchNo = CalcWavevecAngSweep2Plate.nearestBranch(k, freqByBranch, wavVecZGV, freqZGV);
                    end
                    zgvList(end+1) = struct('branchNo', branchNo, 'wavVecZGV', wavVecZGV, ...
                        'freqZGV', freqZGV, 'angle', obj.angles(angleIdx), ...
                        'velGroupTrans', velGroupTrans, ...
                        'isTrueZGV', abs(velGroupTrans) <= CalcWavevecAngSweep2Plate.ZGV_VG_REL_TOL * abs(velPhase)); %#ok<AGROW>
                end
            end
        end

        function result = classifyBranches(obj, angleIdx, confidenceThresh)
        %%% Multi-layer-aware mode-type classification with a confidence-driven
        %%% fallback. A per-point confidence score rather than a static geometric
        %%% threshold (e.g. a relDepth-ratio cutoff), because how much a multi-layer
        %%% stack's asymmetry matters depends on kd (wavevector x thickness), not on
        %%% geometry alone, and a fixed cutoff cannot track that.
        %%% >> result = classifyBranches(obj, angleIdx, confidenceThresh)
        %%%
        %%% INPUT:  (opt) 1. angleIdx          .. which angle to classify at (default 1)
        %%%         (opt) 2. confidenceThresh  .. RitzSolverLamb.modeTypeConfidence cutoff
        %%%                                      below which a branch's physical label is
        %%%                                      not trusted (default = 0.3, UNVALIDATED --
        %%%                                      see the note below)
        %%% OUTPUT: result .. struct:
        %%%           isExactMirrorSymmetric -- obj.structure.isMirrorSymmetric(), a
        %%%                                     WHOLE-STRUCTURE flag: every mode's
        %%%                                     mid-plane parity is exact
        %%%           isExactSagittalSymmetric -- obj.structure.isSagittalMirrorSymmetric
        %%%                                     at this angle: the Lamb/TH split is exact.
        %%%                                     Both true: every label is rigorous
        %%%                                     regardless of confidence
        %%%           confidenceThresh       -- echoes the input
        %%%           branchModeNoMap        -- [nWav x nCompModes], from
        %%%                                     DispersionBranches.track
        %%%           modesTypeByBranch      -- [nWav x nCompModes], raw modesType
        %%%                                     (1=Asymm, 2=Symm, 3=TH-Asymm, 4=TH-Symm)
        %%%                                     re-indexed into branch order
        %%%           confidenceByBranch     -- [nWav x nCompModes], likewise for
        %%%                                     modeTypeConfidence
        %%%           trustLabelByBranch     -- [nWav x nCompModes] logical: (both exact
        %%%                                     flags) OR confidence >= confidenceThresh,
        %%%                                     evaluated PER (step, branch) since
        %%%                                     confidence varies along a branch, e.g.
        %%%                                     near a crossing
        %%%           displayLabel           -- {nWav x nCompModes} cellstr: the type name
        %%%                                     where trustLabelByBranch is true; else
        %%%                                     'Asymm-class'/'Symm-class' where only the
        %%%                                     parity is exact (isExactMirrorSymmetric
        %%%                                     alone); else 'branch N'. Feed straight to
        %%%                                     a colour-by-string categorical plot
        %%%
        %%% Requires retentionLevel='full', inherited from DispersionBranches.track.
        %%%
        %%% Default threshold (0.3): confidence = top parity fraction minus the runner-up
        %%% (RitzSolverLamb.modeTypeConfidence), so 0 is exactly tied and 1 is a pure
        %%% type. NOT tuned against measured multi-layer data or the literature. Treat it
        %%% as a starting point, not a derived constant.

            if nargin < 2 || isempty(angleIdx); angleIdx = 1; end
            if nargin < 3 || isempty(confidenceThresh); confidenceThresh = 0.3; end

            branchModeNoMap = DispersionBranches.track(obj, angleIdx);
            [nWav, nCompModes] = size(branchModeNoMap);

            isExactSymmetric = obj.structure.isMirrorSymmetric();
            isExactSagittal  = obj.structure.isSagittalMirrorSymmetric(obj.angles(angleIdx));

            typeNames   = {'Asymm', 'Symm', 'TH-Asymm', 'TH-Symm'};
            parityNames = {'Asymm-class', 'Symm-class', 'Asymm-class', 'Symm-class'};

            modesTypeByBranch  = CalcWavevecAngSweep2Plate.reindexByBranch( ...
                obj.modesTypeMap3D, angleIdx, branchModeNoMap);
            confidenceByBranch = CalcWavevecAngSweep2Plate.reindexByBranch( ...
                obj.modeTypeConfidenceMap3D, angleIdx, branchModeNoMap);

            %%% Per (step, branch), not once globally -- confidence varies along a
            %%% branch. Mid-plane AND sagittal mirrors make every label rigorous; the
            %%% mid-plane mirror alone makes only the parity rigorous.
            trustLabelByBranch = (isExactSymmetric && isExactSagittal) | ...
                (confidenceByBranch >= confidenceThresh);

            displayLabel = cell(nWav, nCompModes);
            for iw = 1:nWav
                for ib = 1:nCompModes
                    if trustLabelByBranch(iw, ib)
                        displayLabel{iw, ib} = typeNames{modesTypeByBranch(iw, ib)};
                    elseif isExactSymmetric
                        displayLabel{iw, ib} = parityNames{modesTypeByBranch(iw, ib)};
                    else
                        displayLabel{iw, ib} = sprintf('branch %d', ib);
                    end
                end
            end

            result = struct();
            result.isExactMirrorSymmetric   = isExactSymmetric;
            result.isExactSagittalSymmetric = isExactSagittal;
            result.confidenceThresh         = confidenceThresh;
            result.branchModeNoMap          = branchModeNoMap;
            result.modesTypeByBranch        = modesTypeByBranch;
            result.confidenceByBranch       = confidenceByBranch;
            result.trustLabelByBranch       = trustLabelByBranch;
            result.displayLabel             = displayLabel;
        end
    end

    methods (Static, Access = private)
        function byBranch = reindexByBranch(map3D, angleIdx, branchModeNoMap)
        %%% Re-index one [nModes x nAngles x nWav] map from RAW Ritz mode order
        %%% into TRACKED BRANCH order, at one fixed angle
        %%% >> byBranch = CalcWavevecAngSweep2Plate.reindexByBranch(map3D, angleIdx, branchModeNoMap)
        %%%
        %%% INPUT:  1. map3D           .. any of this class's *Map3D properties
        %%%         2. angleIdx        .. which angle slice to take
        %%%         3. branchModeNoMap .. [nWav x nCompModes] from
        %%%                              DispersionBranches.track
        %%% OUTPUT: byBranch .. [nWav x nCompModes], byBranch(iw,ib) = that branch's
        %%%                     value at wavevector step iw
        %%%
        %%% The indirection is not optional: the RAW Ritz mode index at a given
        %%% wavevector step is a numerical ordering by eigenvalue, not a physical
        %%% identity, and it shifts between steps as branches approach and cross.
        %%% branchModeNoMap is the MAC-tracked answer to "which raw index continues
        %%% this branch at this step", so reading any per-mode quantity along a branch
        %%% has to go through it.
            [nWav, nBranch] = size(branchModeNoMap);
            byBranch = NaN(nWav, nBranch);
            for iw = 1:nWav
                rawSlice = map3D(:, angleIdx, iw);   % [nCompModes x 1]
                byBranch(iw, :) = rawSlice(branchModeNoMap(iw, :)).';
            end
        end
    end

    methods (Static, Access = private)
        function cands = zgvCandidates(k, vg, f, vp, rawMode)
        %%% Newton seeds for the ZGV points of one tracked branch (see locateZGVs)
        %%% >> cands = CalcWavevecAngSweep2Plate.zgvCandidates(k, vg, f, vp, rawMode)
        %%%
        %%% INPUT:  k, vg, f, vp .. [nWav x 1] wavevectors and the branch's group velocity,
        %%%                        frequency and phase velocity there
        %%%         rawMode      .. [nWav x 1] the branch's raw Ritz mode number per step
        %%% OUTPUT: cands .. struct array in wavevector order, fields
        %%%                  seeds    -- [nSeed x 4] rows (k0, f0, iw, rawMode), tried in
        %%%                              order; iw = 0: no stored eigenvector (interpolated
        %%%                              seed), else the seed vector is the sweep's own
        %%%                              eigenvector of rawMode at step iw
        %%%                  kWin     -- a converged point inside it is taken at once
        %%%                  kBracket -- inside it the point belongs to this branch
        %%%
        %%% (1) A sign change of vg between steps iw and iw+1: the linear interpolation of
        %%%     the zero, then the two ends; window = the bracket widened by one step.
        %%% (2) A dip: an interior step where |vg| is below half its value at both
        %%%     neighbours and below ZGV_DIP_REL_TOL * vp, with no sign change next to it. Two ZGV points of one
        %%%     branch between two grid steps leave vg with the same sign at both ends
        %%%     (seen on a three-layer stack); |vg| then dips at the nearer end. Three
        %%%     single-seed candidates, from the dip step and from its two neighbours, each
        %%%     converging onto the ZGV point nearest it (duplicates are merged later).
            nW = numel(k);
            cands = struct('iw', {}, 'seeds', {}, 'kWin', {}, 'kBracket', {});
            sc = sign(vg(1:end-1)) .* sign(vg(2:end)) < 0;
            for iw0 = find(sc(:))'
                t  = vg(iw0) / (vg(iw0) - vg(iw0 + 1));
                k0 = k(iw0) + t * (k(iw0 + 1) - k(iw0));
                f0 = f(iw0) + t * (f(iw0 + 1) - f(iw0));
                cands(end+1) = struct('iw', iw0, ...
                    'seeds', [k0, f0, 0, 0; k(iw0), f(iw0), iw0, rawMode(iw0); ...
                              k(iw0 + 1), f(iw0 + 1), iw0 + 1, rawMode(iw0 + 1)], ...
                    'kWin', k([max(iw0 - 1, 1), min(iw0 + 2, nW)]).', ...
                    'kBracket', k([iw0, iw0 + 1]).'); %#ok<AGROW>
            end
            a = abs(vg);
            for iw = 2:nW-1
                isDip = a(iw) < 0.5 * min(a(iw - 1), a(iw + 1)) && ~sc(iw - 1) && ~sc(iw) ...
                    && a(iw) <= CalcWavevecAngSweep2Plate.ZGV_DIP_REL_TOL * abs(vp(iw));
                if ~isDip
                    continue
                end
                for j = [iw - 1, iw, iw + 1]
                    cands(end+1) = struct('iw', iw + 0.1 * (j - iw), ...
                        'seeds', [k(j), f(j), j, rawMode(j)], ...
                        'kWin', k([iw - 1, iw + 1]).', 'kBracket', k([iw - 1, iw + 1]).'); %#ok<AGROW>
                end
            end
            [~, o] = sort([cands.iw]);
            cands = cands(o);
        end

        function ib = nearestBranch(wavVecs, freqByBranch, wavVecZGV, freqZGV)
        %%% Tracked branch nearest a point in frequency, in the wavevector step containing it
        %%% >> ib = CalcWavevecAngSweep2Plate.nearestBranch(wavVecs, freqByBranch, wavVecZGV, freqZGV)
        %%%
        %%% Linear interpolation of every branch's frequency to wavVecZGV inside its step
        %%% (the end steps when wavVecZGV sits on the range boundary).
            iw = find(wavVecs <= wavVecZGV, 1, 'last');
            if isempty(iw); iw = 1; end
            iw = min(iw, numel(wavVecs) - 1);
            t  = (wavVecZGV - wavVecs(iw)) / (wavVecs(iw + 1) - wavVecs(iw));
            fI = freqByBranch(iw, :) + t * (freqByBranch(iw + 1, :) - freqByBranch(iw, :));
            [~, ib] = min(abs(fI - freqZGV));
        end
    end

    methods (Access = private)
        function [wavVecZGV, freqZGV, isRoot, velGroupTrans, velPhase] = refineZGV(obj, angleIdx, poly, seeds, kWin, precompStruct)
        %%% One Newton ZGV refinement from up to three seeds, then its check (see locateZGVs)
        %%% >> [wavVecZGV, freqZGV, isRoot, velGroupTrans, velPhase] = refineZGV(obj, angleIdx, poly, seeds, kWin, precompStruct)
        %%%
        %%% INPUT:  seeds .. [nSeed x 4] rows (k0, f0, iw, rawMode), tried in order
        %%%                  (zgvCandidates); iw > 0: the seed vector is the sweep's stored
        %%%                  eigenvector of rawMode at step iw (no eigen-solve)
        %%%         kWin  .. [kLo, kHi]: a converged point inside it is taken at once; else
        %%%                  the first converged point inside the swept range
        %%% OUTPUT: wavVecZGV, freqZGV .. the Newton point (NaN when none converged)
        %%%         isRoot             .. converged, inside the swept k range, and the
        %%%                               independent solve's matching mode has
        %%%                               |vg| <= ZGV_VG_REL_TOL * vp
        %%%         velGroupTrans      .. that mode's c_gtheta (NaN when not isRoot)
        %%%         velPhase           .. that mode's phase velocity (NaN when not isRoot)
            isRoot = false;  velGroupTrans = NaN;  velPhase = NaN;
            wavVecZGV = NaN;  freqZGV = NaN;
            kTol = 1e-12 * max(obj.wavVecs);
            for is = 1:size(seeds, 1)
                u0 = [];
                if seeds(is, 3) > 0
                    u0 = obj.sweeps{seeds(is, 3)}.cs{angleIdx}.ritzSolver.eigVecs_con(:, seeds(is, 4));
                end
                [kN, omegaN, ~, isConverged] = ZgvNewton.solve(poly.L2, poly.L1, poly.L0, poly.M, ...
                    seeds(is, 1), 2*pi * seeds(is, 2), u0);
                if ~isConverged || kN < min(obj.wavVecs) - kTol || kN > max(obj.wavVecs) + kTol
                    continue
                end
                if isnan(wavVecZGV)                       % first converged point in range
                    wavVecZGV = kN;  freqZGV = omegaN / (2*pi);
                end
                if kN >= kWin(1) && kN <= kWin(2)          % local: take it
                    wavVecZGV = kN;  freqZGV = omegaN / (2*pi);
                    break
                end
            end
            if isnan(wavVecZGV)
                return
            end
            rs = RitzSolverLamb(obj.structure, obj.angles(angleIdx), obj.realThickness, wavVecZGV);
            rs = rs.calc(false, precompStruct);
            [dev, m] = min(abs(rs.freqsPhysical - freqZGV));
            if dev > CalcWavevecAngSweep2Plate.ZGV_FREQ_MATCH_REL_TOL * freqZGV
                return
            end
            isRoot = abs(rs.velsGroup(m)) <= ...
                CalcWavevecAngSweep2Plate.ZGV_VG_REL_TOL * abs(rs.vels(m));
            if isRoot
                velGroupTrans = rs.velsGroupTrans(m);
                velPhase      = rs.vels(m);
            end
        end

        function obj = assembleMaps(obj)
        %%% 3D re-indexing of each wavevector-slice's already-computed 2D maps
        %%% (sweeps{iWav}.velsMap etc.) -- no new classification logic.
        %%% >> obj = assembleMaps(obj)
            nWav    = numel(obj.sweeps);
            nAngles = numel(obj.angles);
            nModes  = size(obj.sweeps{1}.velsMap, 1);

            obj.velsMap3D               = NaN(nModes, nAngles, nWav);
            obj.velsGroupMap3D          = NaN(nModes, nAngles, nWav);
            obj.velsGroupTransMap3D     = NaN(nModes, nAngles, nWav);
            obj.modesTypeMap3D          = NaN(nModes, nAngles, nWav);
            obj.modeTypeConfidenceMap3D = NaN(nModes, nAngles, nWav);
            obj.freqsPhysicalMap3D      = NaN(nModes, nAngles, nWav);

            for iWav = 1:nWav
                obj.velsMap3D(:, :, iWav)               = obj.sweeps{iWav}.velsMap;
                obj.velsGroupMap3D(:, :, iWav)           = obj.sweeps{iWav}.velsGroupMap;
                obj.velsGroupTransMap3D(:, :, iWav)      = obj.sweeps{iWav}.velsGroupTransMap;
                obj.modesTypeMap3D(:, :, iWav)           = obj.sweeps{iWav}.modesTypeMap;
                obj.modeTypeConfidenceMap3D(:, :, iWav)  = obj.sweeps{iWav}.modeTypeConfidenceMap;
                obj.freqsPhysicalMap3D(:, :, iWav)       = obj.sweeps{iWav}.freqsPhysicalMap;
            end
        end
    end

    methods (Access = protected)
        function cp = copyElement(obj)
        %%% Deep-copy override (matlab.mixin.Copyable)
        %%% >> cp = copyElement(obj)
        %%% obj.sweeps holds handle-class CalcAngSweepPlate objects, which the default
        %%% shallow copy would alias. Each deep-copies its own cs{} array in turn.
            cp = copyElement@matlab.mixin.Copyable(obj);
            cp.sweeps = cellfun(@copy, obj.sweeps, 'UniformOutput', false);
        end
    end
end
