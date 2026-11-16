classdef ElasticityInversion < handle & matlab.mixin.Copyable
%%% ElasticityInversion — the optimisation orchestrator. Owns the Structure, the
%%%   inverted-layer index, the MeasuredVelocities, the live ModeAssociation, the
%%%   SensitivityAnalysis, the optimiser configuration and the results.
%%%
%%% EXPERIMENTAL: validated on synthetic data only; the interface may change.
%%%
%%%   calc() runs an lsqnonlin schedule over the leading P transformed elastic
%%%   constants, holding the remaining ones frozen at their initial-guess values.
%%%   Two modes (obj.stepwise): the default GROWS P from optimisedConstantRange(1) to
%%%   (2), one lsqnonlin stage per value, each warm-started from the previous stage's
%%%   result -- safer when the initial guess is distant/uncertain, at the cost of extra
%%%   stages; stepwise=false instead runs ONE lsqnonlin call directly at the full
%%%   P = optimisedConstantRange(2), which is usually enough in practice, leaving
%%%   step-wise growth as the fallback for a run that fails to converge cleanly. Each
%%%   objective evaluation (costAndJacobian) runs three passes per dataset:
%%%     PASS 1, the Ritz/pseudo-SAW channel, per unique associated angle -- assemble the
%%%       constrained stiffness by LINEAR SUPERPOSITION of precomputed blocks (K is
%%%       exactly linear in the elastic constants, so no forward sweep is ever re-run),
%%%       solve the eigenproblem, M-orthonormalise, and RE-IDENTIFY each associated mode
%%%       by M-weighted MAC against its anchor eigenvector (optionally STEERED by a live
%%%       Green's-function peak for pseudo-SAW points -- see pseudoSAWRefinement).
%%%       Residual and analytic Jacobian come from the same Rayleigh-quotient pattern as
%%%       the sensitivity matrix, evaluated at the CURRENT eigenvectors.
%%%     PASS 2, the limiting channel, per angle carrying limiting points -- compute the
%%%       labelled Christoffel root set, RE-IDENTIFY each anchored root by
%%%       polarisation x velocity-proximity x theta-proximity x branch continuity
%%%       (ModeAssociation.rebindLimiting), take the residual from the tracked root, and
%%%       build the Jacobian column by central difference using ANCHORED LOCAL
%%%       REFINEMENT from the tracked theta* rather than 2P further full sagittal scans.
%%%     PASS 3, the pulse-echo phase-velocity and bulk-modulus channels -- direct
%%%       evaluation plus a central difference, no association, no tracking.
%%%   Then it promotes the rebound association (and, under the "chained" anchor policy,
%%%   the updated anchors) to the live state ONLY when the pooled cost improved, so a
%%%   rejected trial step cannot corrupt the tracking.
%%%
%%%   Passes 1 and 2 both honour the exclusion rule: a point whose re-identification
%%%   confidence collapsed keeps its ROW but has both its residual and its Jacobian row
%%%   ZEROED (never deleted -- lsqnonlin requires a constant residual length). Zeroing is
%%%   equivalent to exclusion for a least-squares solve and keeps the residual and the
%%%   Jacobian row consistent with each other.
%%%
%%% JOINT MULTI-DATASET INVERSION (several crystallographic cuts of the SAME crystal,
%%%   sharing one target elasticity): `structure`, `measured` and `association` are OBJECT
%%%   ARRAYS, one element per dataset d = 1..nDatasets (length 1 is the ordinary
%%%   single-dataset case). All datasets share the SAME invertedLayerNo LAYER INDEX and
%%%   the same `sensitivity` (built by SensitivityAnalysis.calc from the SAME dataset
%%%   arrays, pooling every dataset's G_d block before ONE eigendecomposition). Each
%%%   dataset keeps its OWN Structure (own dirMat/other-layers/CalcParams) and its OWN live
%%%   ModeAssociation, re-bound independently every objective evaluation (mode tracking is
%%%   strictly per dataset). Every dataset's residual/Jacobian rows are concatenated into
%%%   ONE r/J handed to a single lsqnonlin call: one shared optimisation, not D separate ones.
%%%   The cost-improvement gate is the TOTAL POOLED cost, and when it improves every
%%%   dataset's rebound association is promoted TOGETHER. At the end of the run calc()
%%%   writes the SAME recovered elasticity into every dataset's
%%%   `structure(d).mater(invertedLayerNo).elConsts`; each dataset's OTHER layers/
%%%   orientation/CalcParams are left untouched.
%%%   NOT supported (same-crystal assumption): datasets disagreeing on symmClassNo/
%%%   nIndepConsts for the inverted layer (the constructor and SensitivityAnalysis.calc
%%%   both error). Different LAYER COUNTS elsewhere in the stack are fine, since only the
%%%   shared layer INDEX must match. recalcForwardAtResult() covers DATASET 1 only.
%%%
%%% Handle class: carries the iteration state, mutated in place across the objective
%%%   evaluations of a run.
%%%
%%% GuliweRR v6.0

    properties (SetAccess = protected)
        structure            Structure           = Structure()   % current elasticity (updated at
                                                                    % the end of calc(), not per-eval);
                                                                    % [1 x nDatasets] array for a
                                                                    % joint multi-dataset inversion
        invertedLayerNo(1,1) double = 0                          % SAME layer index in every dataset
        nDatasets(1,1)       double = 1                          % D; numel(structure)
        measured             MeasuredVelocities  = MeasuredVelocities()   % array, one per dataset
        association          ModeAssociation     = ModeAssociation()   % LIVE association(s), array
                                                                          % one per dataset, promoted
                                                                          % TOGETHER on improving evals
                                                                          % (pooled cost gate)
        sensitivity           SensitivityAnalysis = SensitivityAnalysis()   % ONE pooled analysis,
                                                                              % shared by every dataset

        optimisedConstantRange double = []   % [first last], or a single N meaning [N N] (run at
                                               % exactly N constants, one-shot-equivalent); empty =>
                                               % [1, sensitivity.nOptimisedConstants]
        stepwise(1,1)         logical = true % true (default): grow P from range(1) to range(2),
                                               % one lsqnonlin stage per value, warm-started from
                                               % the previous stage's result (the safer choice
                                               % for a distant/uncertain initial guess). false:
                                               % ONE lsqnonlin call directly at P = range(2), all
                                               % free constants optimised at once from the
                                               % initial guess. range(1) is unused when
                                               % stepwise=false.
        anchorPolicy(1,:)     char    = 'chained'   % 'fixed' | 'chained'
        rebindVelocityWindow(1,1) double  = ModeAssociation.DEFAULT_REBIND_WINDOW   % [km/s]
        includePseudoSAW(1,1) logical = true        % pseudo-Rayleigh points contribute residual
                                                      % rows by DEFAULT: a TGS angular map can be
                                                      % carried by the pseudo-SAW branch alone over
                                                      % wide angle ranges, and excluding it
                                                      % discards the only data there. Taken from
                                                      % the SensitivityAnalysis (an explicit value
                                                      % must equal it), or G's rows and the
                                                      % residual's rows would describe different
                                                      % datasets.
        pointWeights          double  = []          % [] => unweighted; one per data row, in
                                                      % the sensitivity's row order (default: the
                                                      % sensitivity's own weights)
        toleranceLevel(1,1)   double  = 2
        lowerBound            double  = []          % on the transformed constants, empty = none
        upperBound            double  = []

        %%% ---- re-identification thresholds and observability ----
        macWarnThreshold(1,1) double  = ModeAssociation.DEFAULT_TRACK_OPTS.macWarnThreshold
                                                % below this MAC a Ritz-class match is low-confidence:
                                                % warn once per (point, stage), zero the row, and do
                                                % NOT let it update the anchor
        limMatchScoreThresh(1,1) double = ModeAssociation.DEFAULT_TRACK_OPTS.limMatchScoreThresh
                                                % below this combined score a limiting-root match is
                                                % rejected as limRootLost -- see
                                                % ModeAssociation.defaultTrackOpts
        limitingRescanCadence(1,1) double = 1 % run the FULL Christoffel root scan every n-th
                                                % objective evaluation, using anchored local
                                                % refinement in between. 1 (default) = always
                                                % rescan; higher values are safe only because a
                                                % failed refinement (LimitingVels.refineRoot's
                                                % ok = false) forces an immediate rescan anyway

        %%% ---- Green's-function pseudo-SAW refinement ----
        pseudoSAWRefinement(1,:) char = 'green'  % 'none' | 'green'. Default 'green': plain
                                                   % MAC tracking can follow a pseudo-SAW point
                                                   % onto a neighbouring member of the Ritz
                                                   % cluster that represents the leaky wave,
                                                   % with a MAC of 1 (see ModeAssociation.rebind).
                                                   % Costs nothing on a dataset with no
                                                   % pseudoSAW points: no such point, no Green
                                                   % solve.
        pseudoSAWGreenCadence(1,1) double = 5    % full Green rescan every n-th PROMOTED evaluation
                                                   % (never on a rejected trial -- a rejected trial
                                                   % must not burn Green solves), plus an immediate
                                                   % single-angle rescan wherever steering failed
                                                   % (a pSAW point found no candidate inside the
                                                   % Green filter: modeCollision)

        costFunctionValue(1,1) double = NaN
        result                struct = struct()
    end

    properties (Access = protected)
        iterState struct = struct()   % precomputed once per calc(): uniqueAngles, M_con,
                                        % K_con_fixed{}, dKdTransformed_con{} (in the real
                                        % gauge where gaugeOn(ia)), gaugeBack
        evalCount(1,1)         double = 0      % objective evaluations this calc(), all trials
        promotedEvalCount(1,1) double = 0      % cost-IMPROVING evaluations this calc()
        nGreenRescans(1,1)     double = 0      % how many single-angle Green solves ran, counting
                                               % one per (angle, rescan); for reporting only
        warnedKeys             cell   = {}     % "warn once per (point, stage)" registry
        initialStructure       Structure       = Structure()         % the constructor's inputs:
        initialAssociation     ModeAssociation = ModeAssociation()   % calc() restarts from them
    end

    methods
        function obj = ElasticityInversion(structure, measured, association, sensitivity, ...
                invertedLayerNo, varargin)
        %%% Construct the inversion orchestrator
        %%% >> obj = ElasticityInversion(structure, measured, association, sensitivity, ...
        %%%              invertedLayerNo, Name, Value, ...)
        %%%
        %%% INPUT: 1. structure       .. Structure at the INITIAL GUESS; scalar for a single
        %%%                             dataset, or a [1 x nDatasets] OBJECT ARRAY for a JOINT
        %%%                             multi-dataset inversion (several crystallographic cuts
        %%%                             of the SAME crystal) -- see the class docblock.
        %%%        2. measured        .. MeasuredVelocities (or array, one per dataset)
        %%%        3. association     .. ModeAssociation (or array, one per dataset), built at
        %%%                             the same initial guess
        %%%        4. sensitivity     .. SensitivityAnalysis, built at the same initial guess
        %%%                             from the SAME dataset arrays (always scalar: one pooled
        %%%                             analysis shared by every dataset)
        %%%        5. invertedLayerNo .. index of the layer whose elasticity is inverted;
        %%%                             SHARED across every dataset (same layer position)
        %%% Name-Value: AnchorPolicy ('fixed'|'chained', default 'chained'),
        %%%             RebindVelocityWindow (default 1.0 km/s),
        %%%             IncludePseudoSAW (default: the sensitivity's; an explicit value
        %%%                 must equal it -- see the property),
        %%%             PointWeights (default: the sensitivity's, i.e. unweighted unless it
        %%%                 was built with weights; one per data row of the sensitivity, in
        %%%                 its row order -- see SensitivityAnalysis),
        %%%             ToleranceLevel (default 2),
        %%%             LowerBound/UpperBound (on the transformed constants, one per
        %%%                 independent constant, default [] = none),
        %%%             OptimisedConstantRange ([first last], or a single N meaning [N N];
        %%%                 default [1, sensitivity.nOptimisedConstants]),
        %%%             Stepwise (default true: grow P stage by stage; false: one
        %%%                 lsqnonlin call at P = range(2) -- see the stepwise property),
        %%%             MacWarnThreshold (default 0.7),
        %%%             LimMatchScoreThresh (default 0.5),
        %%%             LimitingRescanCadence (default 1),
        %%%             PseudoSAWRefinement ('none'|'green', default 'green'),
        %%%             PseudoSAWGreenCadence (default 5)

            if nargin == 0
                return
            end

            validateattributes(structure, {'Structure'}, {'nonempty'}, ...
                'ElasticityInversion', 'structure');
            validateattributes(measured, {'MeasuredVelocities'}, {'nonempty'}, ...
                'ElasticityInversion', 'measured');
            validateattributes(association, {'ModeAssociation'}, {'nonempty'}, ...
                'ElasticityInversion', 'association');
            validateattributes(sensitivity, {'SensitivityAnalysis'}, {'scalar'}, ...
                'ElasticityInversion', 'sensitivity');

            [nDatasets, symmNoAll, LAll] = SensitivityAnalysis.checkDatasets(structure, ...
                association, measured, invertedLayerNo, 'ElasticityInversion:constructor');
            if any(symmNoAll ~= sensitivity.symmClassNo) || any(LAll ~= sensitivity.nIndepConsts)
                error('ElasticityInversion:constructor:symmetryMismatch', ...
                    ['All jointly-inverted datasets'' inverted layer must share the symmetry ', ...
                     'class and independent-constant count the supplied sensitivity was built ', ...
                     'with (it is the same crystal): symmClassNo = [%s] (sensitivity: %d), ', ...
                     'nIndepConsts = [%s] (sensitivity: %d).'], num2str(symmNoAll(:).'), ...
                    sensitivity.symmClassNo, num2str(LAll(:).'), sensitivity.nIndepConsts);
            end

            obj.structure       = structure;
            obj.measured        = measured;
            obj.association     = association;
            obj.sensitivity     = sensitivity;
            obj.invertedLayerNo = invertedLayerNo;
            obj.nDatasets       = nDatasets;

            %%% Every option's default is its property's declared default: obj holds
            %%% those until the assignments below, so they are defined once.
            ip = inputParser;
            addParameter(ip, 'AnchorPolicy',           obj.anchorPolicy);
            addParameter(ip, 'RebindVelocityWindow',   obj.rebindVelocityWindow);
            addParameter(ip, 'IncludePseudoSAW',       obj.includePseudoSAW);
            addParameter(ip, 'PointWeights',           obj.pointWeights);
            addParameter(ip, 'ToleranceLevel',         obj.toleranceLevel);
            addParameter(ip, 'LowerBound',             obj.lowerBound);
            addParameter(ip, 'UpperBound',             obj.upperBound);
            addParameter(ip, 'OptimisedConstantRange', obj.optimisedConstantRange);
            addParameter(ip, 'Stepwise',               obj.stepwise);
            addParameter(ip, 'MacWarnThreshold',       obj.macWarnThreshold);
            addParameter(ip, 'LimMatchScoreThresh',    obj.limMatchScoreThresh);
            addParameter(ip, 'LimitingRescanCadence',  obj.limitingRescanCadence);
            addParameter(ip, 'PseudoSAWRefinement',    obj.pseudoSAWRefinement);
            addParameter(ip, 'PseudoSAWGreenCadence',  obj.pseudoSAWGreenCadence);
            parse(ip, varargin{:});
            opts = ip.Results;

            refinement = char(opts.PseudoSAWRefinement);
            if ~ismember(refinement, {'none', 'green'})
                error('ElasticityInversion:constructor:badPseudoSAWRefinement', ...
                    'PseudoSAWRefinement must be ''none'' or ''green''; got ''%s''.', refinement);
            end
            anchorPolicyIn = char(opts.AnchorPolicy);
            if ~ismember(anchorPolicyIn, {'fixed', 'chained'})
                error('ElasticityInversion:constructor:badAnchorPolicy', ...
                    'AnchorPolicy must be ''fixed'' or ''chained''; got ''%s''.', anchorPolicyIn);
            end
            validateattributes(opts.RebindVelocityWindow, {'double'}, ...
                {'scalar', 'real', 'finite', 'positive'}, 'ElasticityInversion', 'RebindVelocityWindow');
            validateattributes(opts.ToleranceLevel, {'double'}, ...
                {'scalar', 'real', 'finite'}, 'ElasticityInversion', 'ToleranceLevel');
            validateattributes(opts.IncludePseudoSAW, {'logical', 'numeric'}, ...
                {'scalar', 'binary'}, 'ElasticityInversion', 'IncludePseudoSAW');
            validateattributes(opts.Stepwise, {'logical', 'numeric'}, ...
                {'scalar', 'binary'}, 'ElasticityInversion', 'Stepwise');
            L = sensitivity.nIndepConsts;
            for bnd = {'LowerBound', 'UpperBound'}
                b = opts.(bnd{1});
                if isempty(b); continue; end
                validateattributes(b, {'double'}, {'vector', 'real', 'nonnan'}, ...
                    'ElasticityInversion', bnd{1});
                if numel(b) ~= L
                    error('ElasticityInversion:constructor:badBounds', ...
                        '%s has %d entries; expected %d (one per transformed constant).', ...
                        bnd{1}, numel(b), L);
                end
            end

            %%% IncludePseudoSAW and PointWeights decide which rows G has and how they are
            %%% weighted, so the residual must use the same: both default to the
            %%% sensitivity's own, and a contradicting value is an error.
            if ismember('IncludePseudoSAW', ip.UsingDefaults)
                opts.IncludePseudoSAW = sensitivity.includePseudoSAW;
            elseif logical(opts.IncludePseudoSAW) ~= sensitivity.includePseudoSAW
                error('ElasticityInversion:constructor:includePseudoSAWMismatch', ...
                    ['IncludePseudoSAW = %d contradicts the sensitivity, which was built with ', ...
                     'IncludePseudoSAW = %d: G''s rows and the residual''s rows would describe ', ...
                     'different datasets.'], logical(opts.IncludePseudoSAW), sensitivity.includePseudoSAW);
            end
            nDataRows = size(sensitivity.sensitivityMatrix, 1);
            if ismember('PointWeights', ip.UsingDefaults)
                if any(sensitivity.pointWeights ~= 1)
                    opts.PointWeights = sensitivity.pointWeights;   % weights G was built with
                end
            elseif ~isempty(opts.PointWeights)
                validateattributes(opts.PointWeights, {'double'}, ...
                    {'vector', 'real', 'finite', 'nonnegative'}, 'ElasticityInversion', 'PointWeights');
                if numel(opts.PointWeights) ~= nDataRows
                    error('ElasticityInversion:constructor:badWeights', ...
                        'PointWeights has %d entries; expected %d (one per data row of the sensitivity).', ...
                        numel(opts.PointWeights), nDataRows);
                end
            end
            validateattributes(opts.MacWarnThreshold, {'double'}, ...
                {'scalar', 'real', 'finite', '>=', 0, '<=', 1}, ...
                'ElasticityInversion', 'MacWarnThreshold');
            validateattributes(opts.LimMatchScoreThresh, {'double'}, ...
                {'scalar', 'real', 'finite', '>=', 0, '<=', 1}, ...
                'ElasticityInversion', 'LimMatchScoreThresh');
            validateattributes(opts.LimitingRescanCadence, {'double'}, ...
                {'scalar', 'integer', '>=', 1}, 'ElasticityInversion', 'LimitingRescanCadence');
            validateattributes(opts.PseudoSAWGreenCadence, {'double'}, ...
                {'scalar', 'integer', '>=', 1}, 'ElasticityInversion', 'PseudoSAWGreenCadence');

            obj.anchorPolicy           = anchorPolicyIn;
            obj.rebindVelocityWindow   = opts.RebindVelocityWindow;
            obj.includePseudoSAW       = logical(opts.IncludePseudoSAW);
            obj.pointWeights           = opts.PointWeights(:);
            obj.toleranceLevel         = opts.ToleranceLevel;
            obj.lowerBound             = opts.LowerBound;
            obj.upperBound             = opts.UpperBound;
            obj.optimisedConstantRange = opts.OptimisedConstantRange;
            obj.stepwise               = logical(opts.Stepwise);
            obj.macWarnThreshold       = opts.MacWarnThreshold;
            obj.limMatchScoreThresh    = opts.LimMatchScoreThresh;
            obj.limitingRescanCadence  = opts.LimitingRescanCadence;
            obj.pseudoSAWRefinement    = refinement;
            obj.pseudoSAWGreenCadence  = opts.PseudoSAWGreenCadence;

            obj.initialStructure   = structure;
            obj.initialAssociation = association;
            obj.costFunctionValue  = NaN;
        end

        function obj = calc(obj)
        %%% Run the lsqnonlin schedule, step-wise (obj.stepwise = true) or in one go
        %%% >> obj = calc(obj)
        %%% Populates obj.result and updates obj.structure to the inverted elasticity.
        %%% obj is a handle -- the return value is the same object, for chaining only.
        %%% Repeatable: every call starts from the constructor's structure and
        %%% association (the initial guess and its binding), not from a previous run's
        %%% promoted anchors, so a second call reproduces the first.

            L = obj.sensitivity.nIndepConsts;
            range = obj.optimisedConstantRange;
            if isempty(range)
                range = [1, obj.sensitivity.nOptimisedConstants];
            elseif isscalar(range)
                range = [range, range];   % N => [N N], i.e. run at exactly N constants
            end
            if numel(range) ~= 2 || range(1) < 1 || range(2) > L || range(1) > range(2)
                error('ElasticityInversion:calc:badRange', ...
                    'optimisedConstantRange must be a scalar N or a [first last] pair with 1 <= first <= last <= %d (L = %d independent constants); got %s.', ...
                    L, L, mat2str(obj.optimisedConstantRange));
            end

            obj.structure   = obj.initialStructure;
            obj.association = obj.initialAssociation;
            obj.prepareIteration();
            obj.costFunctionValue = NaN;   % reset at the top of calc(), never carried from a prior run
            obj.evalCount          = 0;
            obj.promotedEvalCount  = 0;
            obj.nGreenRescans      = 0;

            xPrev = [];
            convergenceLog = struct('P', {}, 'exitflag', {}, 'output', {}, 'resnorm', {});

            optOpts = optimoptions('lsqnonlin', ...
                'Algorithm', 'trust-region-reflective', ...
                'SpecifyObjectiveGradient', true, ...
                'Display', 'iter-detailed', ...
                'StepTolerance',       10^-(3 + obj.toleranceLevel), ...
                'FunctionTolerance',   10^-(5 + obj.toleranceLevel), ...
                'OptimalityTolerance', 10^-(5 + obj.toleranceLevel));

            if obj.stepwise
                Plist = range(1):range(2);
            else
                Plist = range(2);   % ONE stage, all free constants at once; range(1) unused
            end

            for P = Plist
                obj.costFunctionValue = NaN;   % reset at the top of EVERY stage: otherwise a new
                                                 % stage's early trials are gated against the
                                                 % previous stage's optimum, get rejected, and
                                                 % leave the association state stale
                obj.warnedKeys = {};             % "once per (point, stage)": a point whose tracking
                                                 % collapsed again in a HIGHER-dimensional stage is
                                                 % new information and is reported again

                if P == Plist(1)
                    x0 = obj.sensitivity.initialTransformedVec(1:P);
                else
                    x0 = [xPrev; obj.sensitivity.initialTransformedVec(P)];
                end

                lb = []; ub = [];
                if ~isempty(obj.lowerBound); lb = obj.lowerBound(1:P); end
                if ~isempty(obj.upperBound); ub = obj.upperBound(1:P); end

                [xSol, resnorm, ~, exitflag, output] = lsqnonlin( ...
                    @(x) obj.costAndJacobian(x, P), x0, lb, ub, optOpts);

                xPrev = xSol;
                convergenceLog(end+1) = struct('P', P, 'exitflag', exitflag, ...
                    'output', output, 'resnorm', resnorm); %#ok<AGROW>
            end

            obj.finalizeResult(xPrev, range(2), convergenceLog);
        end

        function sweep = recalcForwardAtResult(obj, angles, vMax, calcMode, retentionLevel)
        %%% Fresh forward sweep at the inverted elasticity, for display/checking
        %%% >> sweep = recalcForwardAtResult(obj, angles, vMax, calcMode, retentionLevel)
        %%% Deliberately separate and opt-in -- costs a full sweep, never called inside calc().
        %%% For a joint inversion this recalculates DATASET 1 only; for any other dataset,
        %%% call Calc.angSweepHS/angSweepLayered directly on obj.structure(d)/obj.association(d).
            if nargin < 2 || isempty(angles);         angles         = obj.association(1).uniqueRitzAngles(); end
            if nargin < 4 || isempty(calcMode);       calcMode       = 'full'; end
            if nargin < 5 || isempty(retentionLevel); retentionLevel = 'full'; end

            structure1 = obj.structure(1);
            if structure1.nLay == 1
                sweep = Calc.angSweepHS(structure1, angles, calcMode, retentionLevel);
            else
                if nargin < 3; vMax = []; end   % optional Green-axis extension (layered)
                sweep = Calc.angSweepLayered(structure1, angles, vMax, calcMode, retentionLevel);
            end
            sweep = sweep.calc();
        end
    end

    methods (Access = private)
        function prepareIteration(obj)
        %%% Precompute the iterate-independent Ritz blocks, once per calc(), per dataset
        %%% >> prepareIteration(obj)
        %%%
        %%% For EVERY dataset independently (each has its own mesh/orientation/other
        %%% layers): the angle-independent M_con, and per unique associated angle the fixed
        %%% (other-layers) stiffness and the per-transformed-constant stiffness derivative.
        %%% K is linear in the elastic constants, so nothing here depends on the current
        %%% iterate -- which is what lets the inner loop skip the forward sweep entirely.
        %%% obj.iterState becomes a [1 x nDatasets] struct array.

            L = obj.sensitivity.nIndepConsts;
            T = obj.sensitivity.transformMatrix;
            invLay = obj.invertedLayerNo;

            iterStateArr = struct('uniqueAngles', {}, 'M_con', {}, ...
                'K_con_fixed', {}, 'dKdTransformed_con', {}, 'gaugeOn', {}, 'gaugeBack', {});

            for d = 1:obj.nDatasets
                structure_d = obj.structure(d);

                precompStruct = CalcCommon.ritzSAWPrecompute(structure_d);
                N = precompStruct.nullConstrMat;
                nCon = size(N, 2);

                uniqueAngles = obj.association(d).uniqueRitzAngles();
                nA = numel(uniqueAngles);

                %%% M_con: elasticity- and angle-independent (GammaMatrix.LayerNormFac depends
                %%% only on relDepth*rho ratios) -- one angle suffices to compute it.
                gm0 = GammaMatrix();
                gm0 = gm0.calc(structure_d, 0, precompStruct.dKdC_cell, precompStruct.d2KdCdD_cell);
                M_con = N.' * gm0.normMatAll * N;

                bVL_below = sum([structure_d.calcParams(1:invLay-1).basisVecL]);
                bVL_lay   = structure_d.calcParams(invLay).basisVecL;
                rowsLay   = bVL_below+1 : bVL_below+bVL_lay;
                Nb        = N(rowsLay, :);
                bVLAll    = [structure_d.calcParams.basisVecL];
                layerEnd  = cumsum(bVLAll);
                layerRowsAll = arrayfun(@(b) layerEnd(b)-bVLAll(b)+1 : layerEnd(b), ...
                    1:numel(bVLAll), 'UniformOutput', false);
                otherLayers = setdiff(1:numel(bVLAll), invLay);

                K_con_fixed        = cell(nA, 1);
                dKdTransformed_con = cell(nA, 1);
                gaugeOn            = false(nA, 1);

                %%% Real gauge u3 -> i u3 (RitzSolver.realGauge), decided per angle: when the
                %%% fixed part and every per-parameter derivative are real in the gauge, so is
                %%% K at every iterate (a real combination of them), and the fit loop solves
                %%% the real symmetric problem; the eigenvectors are mapped back with gaugeBack.
                dg        = RitzSolver.gaugeVector(structure_d.calcParams);
                gaugeBack = N.' * (dg .* N);

                for ia = 1:nA
                    angleDeg = uniqueAngles(ia);

                    gmA = GammaMatrix();
                    gmA = gmA.calc(structure_d, angleDeg, precompStruct.dKdC_cell, precompStruct.d2KdCdD_cell);
                    K_full = gmA.gammaMatAll;
                    K_full(rowsLay, rowsLay) = 0;   % remove the inverted layer's own contribution --
                                                      % it is rebuilt every evaluation from the L-term sum

                    %%% dK/dc is non-zero only on the inverted layer's diagonal block, and K_full
                    %%% is block-diagonal (the layers couple only through N), so everything
                    %%% below is formed and projected per block: N(rows,:)' B N(rows,:).
                    dG_blk = GammaMatrix.layerBlockDerivs(structure_d, invLay, precompStruct.dKdC_cell, angleDeg);
                    bVL    = size(dG_blk, 1);
                    dGT    = reshape(reshape(dG_blk, [], L) * T, bVL, bVL, L);   % dK/d(transformed p)

                    Kg   = conj(dg) .* K_full .* dg.';
                    dgL  = dg(rowsLay);
                    dGTg = conj(dgL) .* dGT .* dgL.';
                    scale = norm(gmA.gammaMatAll, 'fro');   % = ||Kg of the unzeroed K||_F, as |dg| = 1
                    gaugeOn(ia) = RitzSolver.isGaugeReal(Kg, scale);
                    for pIdx = 1:L
                        gaugeOn(ia) = gaugeOn(ia) && RitzSolver.isGaugeReal(dGTg(:,:,pIdx));
                    end
                    if gaugeOn(ia)
                        K_full = real(Kg);
                        dGT    = real(dGTg);
                    end

                    %%% K_con_fixed: the other layers' blocks only (the inverted layer's is zero)
                    K_con_fixed{ia} = zeros(nCon, nCon, 'like', K_full);
                    for ib = otherLayers
                        r = layerRowsAll{ib};
                        Nr = N(r, :);
                        K_con_fixed{ia} = K_con_fixed{ia} + Nr.' * K_full(r, r) * Nr;
                    end
                    dKdT = zeros(nCon, nCon, L, 'like', dGT);
                    for pIdx = 1:L
                        dKdT(:,:,pIdx) = Nb.' * dGT(:,:,pIdx) * Nb;
                    end
                    dKdTransformed_con{ia} = dKdT;
                end

                iterStateArr(d).uniqueAngles       = uniqueAngles;
                iterStateArr(d).M_con              = M_con;
                iterStateArr(d).K_con_fixed         = K_con_fixed;
                iterStateArr(d).dKdTransformed_con  = dKdTransformed_con;
                iterStateArr(d).gaugeOn            = gaugeOn;
                iterStateArr(d).gaugeBack          = gaugeBack;
            end

            obj.iterState = iterStateArr;
        end

        function [r, J] = costAndJacobian(obj, x, P)
        %%% lsqnonlin objective: residuals (measured - calculated) and analytic Jacobian,
        %%% POOLED over every dataset
        %%% >> [r, J] = costAndJacobian(obj, x, P)
        %%%
        %%% One shared candidate elasticity voigtVec, mapped through each dataset's own
        %%% K_con_fixed/dKdTransformed_con. Mutates obj (a handle): promotes EVERY dataset's
        %%% rebound association to the live state TOGETHER, and updates
        %%% obj.costFunctionValue, ONLY when the TOTAL POOLED cost improves on the best seen
        %%% so far in the current stage.

            sens = obj.sensitivity;
            L = sens.nIndepConsts;
            T = sens.transformMatrix;

            chatFull = sens.initialTransformedVec;
            chatFull(1:P) = x(:);
            voigtVec = T * chatFull;   % SAME candidate elasticity for every dataset

            obj.evalCount = obj.evalCount + 1;
            trackOpts = obj.currentTrackOpts();
            updateAnchor = strcmp(obj.anchorPolicy, 'chained');

            rParts = {};
            JParts = {};
            wRow   = zeros(0, 1);   % each row's position in the sensitivity's row order
            rowOffset = 0;          % rows of the datasets before d
            activeCls = ModeAssociation.activeRitzClasses(obj.includePseudoSAW);

            scratchAssocArr = obj.association;   % value-class array; copied on assignment
            greenEnabled    = trackOpts.useGreenSteering;

            for d = 1:obj.nDatasets
                scratchAssoc = scratchAssocArr(d);
                structureGreen = [];   % built lazily, only if a Green rescan is actually needed
                %%% G orders a dataset's Ritz rows by point index; pass 1 below appends them
                %%% by angle. ritzRank maps a point to its G row, so PointWeights reach the
                %%% right rows whatever order the points are in.
                isRitzRow = ismember({scratchAssoc.points.bindingClass}, activeCls);
                ritzRank  = cumsum(isRitzRow);
                nRitzRows = sum(isRitzRow);

                %%% ================= PASS 1: Ritz-mode / pseudo-SAW channel =================
                uniqueAngles = obj.iterState(d).uniqueAngles;
                M_con        = obj.iterState(d).M_con;

                for ia = 1:numel(uniqueAngles)
                    angleDeg = uniqueAngles(ia);
                    K_con = obj.iterState(d).K_con_fixed{ia};
                    dKdT  = obj.iterState(d).dKdTransformed_con{ia};
                    for pIdx = 1:L
                        K_con = K_con + chatFull(pIdx) .* dKdT(:,:,pIdx);
                    end

                    %%% Vs: in the basis K_con and dKdT are stored in (the real gauge when
                    %%% gaugeOn), for the Jacobian's Rayleigh quotients; V: the original basis,
                    %%% as rebind's anchors are. gaugeBack is unitary, so the quotients agree.
                    [eigValSorted, Vs] = RitzSolver.solveConstrainedEig(K_con, M_con);
                    if obj.iterState(d).gaugeOn(ia)
                        V = obj.iterState(d).gaugeBack * Vs;
                    else
                        V = Vs;
                    end
                    vels = abs(sqrt(eigValSorted)) / (2*pi);

                    assocBeforeRebind = scratchAssoc;
                    scratchAssoc = scratchAssoc.rebind(angleDeg, V, M_con, vels, ...
                        obj.rebindVelocityWindow, updateAnchor, trackOpts);

                    %%% IMMEDIATE, out-of-cadence Green rescan at THIS ANGLE ONLY, when
                    %%% steering there just failed. The stored peak was measured at the last
                    %%% PROMOTED iterate; after a large trial step the resonance can move out
                    %%% from under it, rebind's Green filter then empties the candidate set
                    %%% and raises modeCollision, and the zeroed row would make the residual
                    %%% vector's content flicker between evaluations -- which lsqnonlin's
                    %%% trust region handles badly. So re-measure the peak here and redo the
                    %%% match from the PRE-rebind state, so the failed attempt does not
                    %%% double-count against the point's low-confidence tally. This is the ONE
                    %%% sanctioned exception to "no Green solves on the trial path": it is
                    %%% conditional, it is one angle, and it happens exactly where obeying a
                    %%% stale peak would otherwise do damage.
                    if greenEnabled
                        idxPS = scratchAssoc.pointsAtAngle(angleDeg, {'pseudoSAW'});
                        needsRescan = false;
                        for kk = 1:numel(idxPS)
                            if scratchAssoc.points(idxPS(kk)).modeCollision
                                needsRescan = true; break
                            end
                        end
                        if needsRescan
                            if isempty(structureGreen)
                                structureGreen = obj.structureAtElasticity(d, voigtVec);
                            end
                            assocBeforeRebind = obj.refreshGreenPeaksAtAngle( ...
                                assocBeforeRebind, structureGreen, angleDeg);
                            scratchAssoc = assocBeforeRebind.rebind(angleDeg, V, M_con, vels, ...
                                obj.rebindVelocityWindow, updateAnchor, trackOpts);
                        end
                    end

                    idxPts = scratchAssoc.pointsAtAngle(angleDeg, ModeAssociation.RITZ_CLASSES);
                    for kk = 1:numel(idxPts)
                        pt = scratchAssoc.points(idxPts(kk));
                        if strcmp(pt.bindingClass, 'pseudoSAW') && ~obj.includePseudoSAW
                            continue
                        end
                        wRow(end+1, 1) = rowOffset + ritzRank(idxPts(kk)); %#ok<AGROW>
                        if ~pt.activeThisEvaluation
                            rParts{end+1} = 0; %#ok<AGROW>
                            JParts{end+1} = zeros(1, P); %#ok<AGROW>
                            continue
                        end
                        alpha = Vs(:, pt.ritzModeNo);
                        rParts{end+1} = pt.measuredVelocity - pt.calculatedVelocity; %#ok<AGROW>
                        Jrow = zeros(1, P);
                        for pIdx = 1:P
                            Jrow(pIdx) = -SensitivityAnalysis.ritzVelocityDerivative( ...
                                real(alpha' * dKdT(:,:,pIdx) * alpha), pt.calculatedVelocity);
                        end
                        JParts{end+1} = Jrow; %#ok<AGROW>
                    end
                end

                %%% ================= PASS 2: the limiting channel =================
                %%% THIS dataset's own Structure (rho, dirMat), at the SAME shared voigtVec.
                matLay = obj.structure(d).mater(obj.invertedLayerNo);
                symmNo = sens.symmClassNo;
                rho    = matLay.rho;

                %%% Read from and write back into scratchAssoc, NOT obj.association(d)
                %%% directly: a rejected trial step must not be able to move an anchor or
                %%% advance a tracked root angle.
                scratchAssoc = obj.trackLimitingRoots(scratchAssoc, matLay, symmNo, voigtVec, ...
                    updateAnchor, trackOpts);

                limIdx = find(strcmp({scratchAssoc.points.bindingClass}, 'limiting'));
                nLim   = numel(limIdx);
                for kk = 1:nLim
                    pt = scratchAssoc.points(limIdx(kk));
                    wRow(end+1, 1) = rowOffset + nRitzRows + kk; %#ok<AGROW>
                    if ~pt.activeThisEvaluation
                        rParts{end+1} = 0; %#ok<AGROW>
                        JParts{end+1} = zeros(1, P); %#ok<AGROW>
                        continue
                    end
                    %%% Anchored LOCAL refinement from this point's own tracked theta*, at
                    %%% the wider limiting-channel step: cheaper than a full rescan and,
                    %%% more important, a correct derivative. A full sagittal scan per
                    %%% perturbation is non-smooth at the perturbation scale (fminbnd's
                    %%% default TolX, a grid-dependent bracket), which can give the
                    %%% derivative the wrong sign.
                    ctx = struct('symmNo', symmNo, 'matLay', matLay, 'pt', pt);
                    [row, okFD] = SensitivityAnalysis.channelDerivative('limiting', ...
                        chatFull, P, @(c) T*c, ctx);
                    if okFD
                        rParts{end+1} = pt.measuredVelocity - pt.calculatedVelocity; %#ok<AGROW>
                        JParts{end+1} = -row; %#ok<AGROW>
                    else
                        %%% A perturbed evaluation lost the tracked root: the Jacobian row
                        %%% is undefined, so the point is excluded from this evaluation
                        %%% exactly like a root rebindLimiting found gone (warned on
                        %%% promotion as limRootLost).
                        scratchAssoc = scratchAssoc.setLimitingLost(limIdx(kk));
                        rParts{end+1} = 0; %#ok<AGROW>
                        JParts{end+1} = zeros(1, P); %#ok<AGROW>
                    end
                end

                %%% ================= PASS 3: pulse-echo and bulk-modulus =================
                %%% No association, no tracking: a direction and a wave type fully determine
                %%% the answer, so there is no identity problem to solve here.
                nPhase = numel(obj.measured(d).phaseVelocities);
                for kk = 1:nPhase
                    pv = obj.measured(d).phaseVelocities(kk);
                    wRow(end+1, 1) = rowOffset + nRitzRows + nLim + kk; %#ok<AGROW>
                    vCalc = SensitivityAnalysis.phaseVelocityAt(voigtVec, symmNo, pv.propagationDir, rho, pv.waveType);
                    rParts{end+1} = pv.velocity - vCalc; %#ok<AGROW>
                    ctx = struct('symmNo', symmNo, 'rho', rho, 'pv', pv);
                    JParts{end+1} = -SensitivityAnalysis.channelDerivative('phase', ...
                        chatFull, P, @(c) T*c, ctx); %#ok<AGROW>
                end

                hasBulk = ~isnan(obj.measured(d).bulkModulus);
                if hasBulk
                    wRow(end+1, 1) = rowOffset + nRitzRows + nLim + nPhase + 1; %#ok<AGROW>
                    vCalc = SensitivityAnalysis.bulkVelocityAt(voigtVec, symmNo, rho);
                    vMeas = sqrt(obj.measured(d).bulkModulus / rho);
                    rParts{end+1} = vMeas - vCalc; %#ok<AGROW>
                    ctx = struct('symmNo', symmNo, 'rho', rho);
                    JParts{end+1} = -SensitivityAnalysis.channelDerivative('bulk', ...
                        chatFull, P, @(c) T*c, ctx); %#ok<AGROW>
                end

                scratchAssocArr(d) = scratchAssoc;
                rowOffset = rowOffset + nRitzRows + nLim + nPhase + double(hasBulk);
            end

            r = cell2mat(rParts(:));
            J = cell2mat(JParts(:));

            if ~isempty(obj.pointWeights)
                if numel(obj.pointWeights) ~= numel(r)
                    error('ElasticityInversion:costAndJacobian:badWeights', ...
                        ['PointWeights has %d entries but the residual has %d rows: the ', ...
                         'association differs from the one the sensitivity was built from.'], ...
                        numel(obj.pointWeights), numel(r));
                end
                w = obj.pointWeights(wRow);
                r = w .* r;
                J = w .* J;
            end

            %%% Promotion: the gate is the TOTAL pooled cost, and on improvement EVERY
            %%% dataset's rebound association is promoted together. All datasets are
            %%% evaluated at the same candidate x inside one trust-region step, so there is
            %%% one accept/reject decision; promoting per dataset would let one dataset's
            %%% live association drift ahead of another's, with no single iterate the
            %%% associations jointly correspond to.
            cost = sum(r.^2);
            if isnan(obj.costFunctionValue) || cost < obj.costFunctionValue
                obj.costFunctionValue = cost;
                obj.association        = scratchAssocArr;   % ALL datasets promoted together
                obj.promotedEvalCount  = obj.promotedEvalCount + 1;
                obj.emitTrackingWarnings(P);
                obj.maybeRefreshGreenPeaks(voigtVec);
            end
        end

        function opts = currentTrackOpts(obj)
        %%% This object's re-identification thresholds, as the struct ModeAssociation's
        %%% rebind/rebindLimiting consume
        %%% >> opts = currentTrackOpts(obj)
            opts = ModeAssociation.defaultTrackOpts();
            opts.macWarnThreshold    = obj.macWarnThreshold;
            opts.limMatchScoreThresh = obj.limMatchScoreThresh;
            opts.useGreenSteering    = strcmp(obj.pseudoSAWRefinement, 'green');
        end

        function assoc = trackLimitingRoots(obj, assoc, matLay, symmNo, voigtVec, updateAnchor, trackOptsIn)
        %%% Re-identify every tracked limiting root at the current iterate
        %%% >> assoc = trackLimitingRoots(obj, assoc, matLay, symmNo, voigtVec, updateAnchor, trackOptsIn)
        %%%
        %%% The full labelled root set is what detects roots APPEARING and VANISHING, so
        %%% by default (limitingRescanCadence = 1) it is recomputed every evaluation. At a
        %%% higher cadence the anchored local refinement is used in between, which is safe
        %%% only because a refinement that fails (LimitingVels.refineRoot's ok = false)
        %%% forces an immediate full rescan of that angle -- and "the anchored root no
        %%% longer exists" is exactly what ok = false means.
            limAngles = assoc.uniqueLimitingAngles();
            if isempty(limAngles)
                return
            end
            fullScanDue = obj.limitingRescanCadence <= 1 ...
                || mod(obj.evalCount - 1, obj.limitingRescanCadence) == 0;

            for ia = 1:numel(limAngles)
                angleDeg   = limAngles(ia);
                idxAtAngle = assoc.pointsAtAngle(angleDeg, {'limiting'});
                doFullScan = fullScanDue;

                if ~doFullScan
                    vRef  = zeros(numel(idxAtAngle), 1);
                    thRef = zeros(numel(idxAtAngle), 1);
                    for kk = 1:numel(idxAtAngle)
                        pt = assoc.points(idxAtAngle(kk));
                        [vRef(kk), thRef(kk), ~, okk] = SensitivityAnalysis.limitingVelocityNear( ...
                            matLay, symmNo, voigtVec, angleDeg, pt.limRootTheta, pt.limBranchIwt);
                        if ~okk
                            doFullScan = true;
                            break
                        end
                    end
                    if ~doFullScan
                        for kk = 1:numel(idxAtAngle)
                            assoc = assoc.setLimitingRefined(idxAtAngle(kk), vRef(kk), thRef(kk));
                        end
                    end
                end

                if doFullScan
                    roots = SensitivityAnalysis.limitingRootsAt(matLay, symmNo, voigtVec, angleDeg);
                    assoc = assoc.rebindLimiting(angleDeg, roots, updateAnchor, trackOptsIn);
                end
            end
        end

        function maybeRefreshGreenPeaks(obj, voigtVec)
        %%% Re-measure every pseudo-SAW point's Green peak, on cadence, at the iterate
        %%% that was just PROMOTED
        %%% >> maybeRefreshGreenPeaks(obj, voigtVec)
        %%%
        %%% Called only from the promotion branch of costAndJacobian, which is what makes
        %%% "every n-th PROMOTED evaluation" true: a rejected trial step never burns a
        %%% Green solve on this path. The first promotion always rescans (the peaks seeded
        %%% by ModeAssociation.fromSweep are from the initial guess and go stale as soon as
        %%% the optimiser moves), then every pseudoSAWGreenCadence-th. The complementary
        %%% path is the immediate, single-angle rescan inside costAndJacobian's pass 1,
        %%% which fires only where steering has actually failed.
        %%%
        %%% This and refreshGreenPeaksAtAngle are the ONLY places the
        %%% "costAndJacobian never touches a Structure" rule is relaxed, and deliberately
        %%% as narrowly as possible: the ElastConsts -> SetElastConsts -> ModifyLayer chain
        %%% runs once per dataset per rescan, and the rotation once per pSAW-carrying
        %%% ANGLE, never universally.
            if ~strcmp(obj.pseudoSAWRefinement, 'green')
                return
            end
            due = obj.promotedEvalCount == 1 ...
                || mod(obj.promotedEvalCount, obj.pseudoSAWGreenCadence) == 0;
            if ~due
                return
            end

            for d = 1:obj.nDatasets
                assoc = obj.association(d);
                idxP  = find(strcmp({assoc.points.bindingClass}, 'pseudoSAW'));
                if isempty(idxP)
                    continue
                end
                structureD = obj.structureAtElasticity(d, voigtVec);
                anglesP = unique([assoc.points(idxP).angleDeg]);
                for ia = 1:numel(anglesP)
                    assoc = obj.refreshGreenPeaksAtAngle(assoc, structureD, anglesP(ia));
                end
                obj.association(d) = assoc;
            end
        end

        function structureD = structureAtElasticity(obj, d, voigtVec)
        %%% Dataset d's Structure with the inverted layer set to a candidate elasticity
        %%% >> structureD = structureAtElasticity(obj, d, voigtVec)
        %%% UNROTATED -- the caller rotates per angle via CalcCommon.rotateStructure.
            elC = ElastConsts(voigtVec, obj.sensitivity.symmClassNo);
            matLayNew = obj.structure(d).mater(obj.invertedLayerNo).SetElastConsts(elC);
            structureD = obj.structure(d).ModifyLayer(obj.invertedLayerNo, matLayNew, []);
        end

        function assoc = refreshGreenPeaksAtAngle(obj, assoc, structureD, angleDeg)
        %%% One Green solve at ONE angle; store the tracked pseudo-SAW peak on every
        %%% pseudo-SAW point bound there
        %%% >> assoc = refreshGreenPeaksAtAngle(obj, assoc, structureD, angleDeg)
        %%%
        %%% INPUT: 1. assoc      .. ModeAssociation to update (a scratch copy on the trial
        %%%                        path, obj.association(d) on the promotion path)
        %%%        2. structureD .. UNROTATED Structure at the current candidate elasticity
        %%%        3. angleDeg   .. propagation angle [deg]
        %%%
        %%% The peaks are those of the forward classifier (pseudoSAWGreenPeaks), so the
        %%% seed fromSweep stored and every rescan measure the same observable. The peak
        %%% is tracked by continuity with the point's OWN PREVIOUS PEAK, not with its
        %%% current calculatedVelocity. Selecting the peak nearest the tracked velocity
        %%% makes the "independent" evidence a mirror of the Ritz tracking: a point that
        %%% has drifted picks the peak nearest where it drifted to, agrees with it
        %%% perfectly, and reports high confidence while being wrong. The peak sequence is
        %%% seeded at construction from the classifier's own matched peak
        %%% (anchorGreenPeakVel) and walks on its own from there.
            idxP = assoc.pointsAtAngle(angleDeg, {'pseudoSAW'});
            if isempty(idxP)
                return
            end

            peaks = ElasticityInversion.pseudoSAWGreenPeaks(structureD, angleDeg);
            obj.nGreenRescans = obj.nGreenRescans + 1;

            for kk = 1:numel(idxP)
                if isempty(peaks)
                    %%% NaN disables steering for this point until a peak reappears.
                    %%% Falling back to plain MAC is the safe failure mode; silently
                    %%% steering towards a peak that is not there is not.
                    assoc = assoc.setGreenPeak(idxP(kk), NaN);
                    continue
                end
                p = assoc.points(idxP(kk));
                if isfinite(p.greenPeakVel) && p.greenPeakVel > 0
                    vRef = p.greenPeakVel;            % continuity in the PEAK's own sequence
                elseif isfinite(p.anchorGreenPeakVel) && p.anchorGreenPeakVel > 0
                    vRef = p.anchorGreenPeakVel;      % the construction-time seed
                else
                    vRef = p.calculatedVelocity;      % last resort, never after a first rescan
                end
                [~, iNear] = min(abs([peaks.vel] - vRef));
                assoc = assoc.setGreenPeak(idxP(kk), peaks(iNear).vel);
            end
        end


        function emitTrackingWarnings(obj, P)
        %%% Warn ONCE per (dataset, point, stage) about a collapsed re-identification
        %%% >> emitTrackingWarnings(obj, P)
        %%%
        %%% A point whose tracking has failed must never drop out of the fit silently.
        %%% Identifiers are ElasticityInversion:tracking:* because this class, not
        %%% ModeAssociation, owns the once-per-stage registry and emits them -- suppress
        %%% with warning('off', 'ElasticityInversion:tracking:lowMAC') etc.
            opts = obj.currentTrackOpts();
            for d = 1:obj.nDatasets
                pts = obj.association(d).points;
                for k = 1:numel(pts)
                    p = pts(k);
                    if strcmp(p.bindingClass, 'limiting')
                        if ~p.limRootLost
                            continue
                        end
                        id  = 'ElasticityInversion:tracking:limRootLost';
                        msg = sprintf(['Dataset %d, point %d (%s at %.6g deg, stage P = %d): the ', ...
                            'tracked Christoffel root was lost at this iterate (no root cleared ', ...
                            'limMatchScoreThresh %.3g, last score %.3g, or a finite-difference ', ...
                            'step lost it); the point''s residual and Jacobian row are zeroed%s.'], ...
                            d, k, p.waveType, p.angleDeg, P, obj.limMatchScoreThresh, ...
                            p.limMatchScore, ElasticityInversion.pinnedSuffix(p));
                    elseif p.modeCollision
                        id  = 'ElasticityInversion:tracking:modeCollision';
                        msg = sprintf(['Dataset %d, point %d (%s at %.6g deg, stage P = %d): the ', ...
                            'one-to-one assignment gave this point no Ritz mode at all. Its ', ...
                            'residual and Jacobian row are zeroed%s.'], ...
                            d, k, p.waveType, p.angleDeg, P, ElasticityInversion.pinnedSuffix(p));
                    elseif p.macValue < obj.macWarnThreshold ...
                            && ~ElasticityInversion.isGreenCorroborated(p, opts)
                        %%% A Green-corroborated pseudo-SAW match is NOT low-confidence
                        %%% just because its MAC is low: steering deliberately moves it off
                        %%% the mode a drifted anchor preferred, and the MAC measures
                        %%% agreement with that stale anchor. greenAgreement is this point's
                        %%% confidence signal, and it is already gated by the veto inside
                        %%% rebind.
                        id  = 'ElasticityInversion:tracking:lowMAC';
                        msg = sprintf(['Dataset %d, point %d (%s at %.6g deg, stage P = %d): MAC ', ...
                            '%.3g is below macWarnThreshold (%.3g). Its residual and Jacobian row ', ...
                            'are zeroed and its anchor was not updated%s.'], ...
                            d, k, p.waveType, p.angleDeg, P, p.macValue, obj.macWarnThreshold, ...
                            ElasticityInversion.pinnedSuffix(p));
                    else
                        continue
                    end

                    key = sprintf('%s|%d|%d|%d', id, d, k, P);
                    if any(strcmp(obj.warnedKeys, key))
                        continue
                    end
                    obj.warnedKeys{end+1} = key;
                    warning(id, '%s', msg);
                end
            end
        end

        function finalizeResult(obj, transformedVecFinal, P, convergenceLog)
        %%% Reverse-transform to Voigt constants, compute errors in both spaces, and
        %%% update every obj.structure(d) and obj.result
        %%% >> finalizeResult(obj, transformedVecFinal, P, convergenceLog)
        %%%
        %%% One crystal, one recovered elasticity: the SAME elConsts is written into every
        %%% dataset's inverted layer; each dataset's own orientation/other-layers/CalcParams
        %%% are left untouched.

            sens = obj.sensitivity;
            L = sens.nIndepConsts;
            T = sens.transformMatrix;

            chatFull = sens.initialTransformedVec;
            chatFull(1:P) = transformedVecFinal(:);
            voigtVec = T * chatFull;

            elConsts = ElastConsts(voigtVec, sens.symmClassNo);
            for d = 1:obj.nDatasets
                matFinal_d = obj.structure(d).mater(obj.invertedLayerNo).SetElastConsts(elConsts);
                obj.structure(d) = obj.structure(d).ModifyLayer(obj.invertedLayerNo, matFinal_d, []);
            end

            Fmin = convergenceLog(end).resnorm;

            activeRitzClasses = ModeAssociation.activeRitzClasses(obj.includePseudoSAW);
            nData = 0;
            for d = 1:obj.nDatasets
                assoc_d = obj.association(d);
                nRitzActive_d = sum(ismember({assoc_d.points.bindingClass}, activeRitzClasses));
                nLim_d   = sum(strcmp({assoc_d.points.bindingClass}, 'limiting'));
                nPhase_d = numel(obj.measured(d).phaseVelocities);
                hasBulk_d = ~isnan(obj.measured(d).bulkModulus);
                nData = nData + nRitzActive_d + nLim_d + nPhase_d + double(hasBulk_d);
            end

            %%% delta_p = sqrt(2*Fmin/(N*psi_p)), Grabec et al. 2024, Eq.(13); computed for
            %%% ALL L combinations (including frozen ones), not just the P optimised, and all
            %%% of them enter Eq.(14) below. For a frozen combination it is a conservative
            %%% "how far could this move unnoticed" proxy, not an estimation error.
            deltaHat = sqrt(2 * Fmin ./ (nData .* sens.sensitivityEigenvalues));
            frozen = false(L, 1);
            frozen(P+1:end) = true;

            %%% sigma_i = sqrt(sum_p T_ip^2 delta_p^2), Grabec et al. 2024, Eq.(14)
            voigtErr = sqrt((T.^2) * (deltaHat.^2));

            obj.result = struct( ...
                'transformedElasticVec',       chatFull, ...
                'transformedElasticVecErrors', deltaHat, ...
                'sensitivityEigenvalues',      sens.sensitivityEigenvalues, ...
                'transformMatrix',             T, ...
                'frozenCombinations',          frozen, ...
                'voigtElasticVec',             voigtVec, ...
                'voigtElasticVecErrors',       voigtErr, ...
                'elConsts',                    elConsts, ...
                'costFunctionValue',           Fmin, ...
                'nDataPoints',                 nData, ...
                'convergence',                 convergenceLog, ...
                'tracking',                    obj.buildTrackingReport());
        end

        function rep = buildTrackingReport(obj)
        %%% Per-class re-identification health, for obj.result.tracking
        %%% >> rep = buildTrackingReport(obj)
        %%%
        %%% A point that spent most of the run below its confidence threshold contributed
        %%% almost nothing to the answer, and the user should be told that here rather
        %%% than left to infer it from a suspiciously small nDataPoints.
            classes = {'ritzMode', 'pseudoSAW', 'limiting'};
            minConf  = struct('ritzMode', NaN, 'pseudoSAW', NaN, 'limiting', NaN);
            nPerCls  = struct('ritzMode', 0,   'pseudoSAW', 0,   'limiting', 0);
            nLowCls  = struct('ritzMode', 0,   'pseudoSAW', 0,   'limiting', 0);

            nCollision = 0; nRootLost = 0; nPinned = 0; nExcluded = 0;
            greenAgree = [];
            lowCounts  = cell(obj.nDatasets, 1);
            opts = obj.currentTrackOpts();

            for d = 1:obj.nDatasets
                assoc = obj.association(d);
                pts   = assoc.points;
                if isempty(pts)
                    lowCounts{d} = zeros(0, 1);
                    continue
                end
                conf = assoc.confidence();
                lowCounts{d} = [pts.lowConfidenceCount].';
                %%% a Green-corroborated pseudo-SAW match is not counted as
                %%% below-threshold on its MAC -- see isGreenCorroborated
                corroborated = arrayfun(@(p) ElasticityInversion.isGreenCorroborated(p, opts), pts).';

                for c = 1:numel(classes)
                    m = strcmp({pts.bindingClass}, classes{c});
                    if ~any(m); continue; end
                    nPerCls.(classes{c}) = nPerCls.(classes{c}) + sum(m);
                    cVals = conf(m);
                    cVals = cVals(~isnan(cVals));
                    if ~isempty(cVals)
                        minConf.(classes{c}) = min([minConf.(classes{c}), min(cVals)]);
                    end
                    if strcmp(classes{c}, 'limiting')
                        thr = obj.limMatchScoreThresh;
                        gated = conf(m) < thr;
                    else
                        thr = obj.macWarnThreshold;
                        gated = conf(m) < thr & ~corroborated(m);
                    end
                    nLowCls.(classes{c}) = nLowCls.(classes{c}) + sum(gated(~isnan(conf(m))));
                end

                nCollision = nCollision + sum([pts.modeCollision]);
                nRootLost  = nRootLost  + sum([pts.limRootLost]);
                nPinned    = nPinned    + sum([pts.pinned]);
                nExcluded  = nExcluded  + sum(~[pts.activeThisEvaluation]);
                ga = [pts(strcmp({pts.bindingClass}, 'pseudoSAW')).greenAgreement];
                greenAgree = [greenAgree, ga(isfinite(ga))]; %#ok<AGROW>
            end

            rep = struct( ...
                'minConfidence',            minConf, ...
                'nPointsPerClass',          nPerCls, ...
                'nBelowThresholdPerClass',  nLowCls, ...
                'nModeCollisions',          nCollision, ...
                'nLimitingRootsLost',       nRootLost, ...
                'nPinnedPoints',            nPinned, ...
                'nExcludedAtFinalEval',     nExcluded, ...
                'lowConfidenceCounts',      {lowCounts}, ...
                'pseudoSAWRefinement',      obj.pseudoSAWRefinement, ...
                'nGreenRescans',            obj.nGreenRescans, ...
                'greenAgreementMax',        maxOrNaN(greenAgree), ...
                'greenAgreementMedian',     medianOrNaN(greenAgree), ...
                'nObjectiveEvaluations',    obj.evalCount, ...
                'nPromotedEvaluations',     obj.promotedEvalCount);

            function v = maxOrNaN(x)
                if isempty(x); v = NaN; else; v = max(x); end
            end
            function v = medianOrNaN(x)
                if isempty(x); v = NaN; else; v = median(x); end
            end
        end
    end

    methods (Static, Access = private)
        function peaks = pseudoSAWGreenPeaks(structureD, angleDeg)
        %%% The pseudo-SAW peaks the forward classifier would report at one angle
        %%% >> peaks = ElasticityInversion.pseudoSAWGreenPeaks(structureD, angleDeg)
        %%%
        %%% INPUT:  1. structureD .. UNROTATED half-space Structure at the candidate
        %%%                         elasticity
        %%%         2. angleDeg   .. propagation angle [deg]
        %%% OUTPUT: peaks .. GreenScore candidates (field vel = the peak velocity) that
        %%%                  pass the default 'green' standard and lie on the
        %%%                  pseudoRayleigh side of minShear
        %%%
        %%% The same observable and the same rules as the 'green' classification that
        %%% seeded the points (CalcSingleHS.calcGreen / calcGreenScore and
        %%% ModeClassification.recognizeFull): G13_eff at SawSpectrum.ETA on the score's
        %%% band at the default axis step, GreenScore.direction (pole-fitted candidates,
        %%% probe slope, score), GreenScore.accepted at the default threshold, and the
        %%% Rayleigh/pseudoRayleigh boundary minShear x (1 + RAYLEIGH_LABEL_TOL). The
        %%% structure is rotated once (CalcCommon.rotateStructure) and the solver gets
        %%% angle 0 -- the rotation convention of docs/reference/architecture/forwardCalcs.md.
            structRot = CalcCommon.rotateStructure(structureD, angleDeg);
            lim = CalcCommon.limitingAndBulk(structRot.mater(structRot.nLay));
            v = GreenScore.bandAxis(lim.vS, lim.vTf, GreenScore.DEFAULT_STEP);
            g = CalcCommon.greenScoreSpectrum(structRot, v, 'ElasticityInversion.pseudoSAWGreenPeaks');
            pk = GreenScore.direction(structureD, angleDeg, g.G13Eff, g.velAxis, lim.vS, lim.vTf, ...
                CalcCommon.surfaceVRef(structureD));
            if isempty(pk)
                peaks = pk;
                return
            end
            keep = GreenScore.accepted(pk, CalcCommon.defaultThreshold('green', false)) ...
                & [pk.vel].' >= lim.vS * (1 + ModeClassification.RAYLEIGH_LABEL_TOL);
            peaks = pk(keep);
        end

        function tf = isGreenCorroborated(p, opts)
        %%% Is this a pseudo-SAW point whose match an INDEPENDENT live Green peak agrees
        %%% with, within the veto distance?
        %%% >> tf = ElasticityInversion.isGreenCorroborated(p, opts)
            tf = opts.useGreenSteering ...
                && strcmp(p.bindingClass, 'pseudoSAW') ...
                && isfinite(p.greenPeakVel) && p.greenPeakVel > 0 ...
                && isfinite(p.greenAgreement) && p.greenAgreement <= opts.greenFilterRel;
        end

        function s = pinnedSuffix(p)
        %%% " (pinned -- kept in the fit anyway)" when a warned point is human-confirmed
        %%% >> s = ElasticityInversion.pinnedSuffix(p)
        %%%
        %%% A pinned point is exempt from exclusion but NOT from warning: a human said
        %%% this binding is right and the system says it is uncertain, and both facts
        %%% should be visible. Excluding it would be the system overriding the human by
        %%% a different route.
            if p.pinned
                s = ' (pinned -- kept in the fit anyway)';
            else
                s = '';
            end
        end
    end

    methods (Access = protected)
        function cp = copyElement(obj)
        %%% Deep-copy override (matlab.mixin.Copyable)
        %%% >> cp = copyElement(obj)
        %%% No property is handle-typed, so the default shallow copy already duplicates
        %%% them correctly; implemented explicitly per house convention.
            cp = copyElement@matlab.mixin.Copyable(obj);
        end
    end
end
