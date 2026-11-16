classdef ModeAssociation
%%% ModeAssociation — the binding table: which measured velocity point corresponds to
%%%   which calculated object. For "ritzMode"/"pseudoSAW" points that object is one
%%%   eigenvector of the Ritz eigenproblem; for "limiting" points it is one ROOT of the
%%%   Christoffel sagittal-plane limiting-velocity problem. One entry per
%%%   (measured point -> calculated counterpart) pair.
%%%
%%% EXPERIMENTAL: validated on synthetic data only; the interface may change.
%%%
%%%   Built initially from a forward sweep (fromSweep) at the initial elasticity guess;
%%%   then re-established at every optimisation iteration (rebind for the Ritz classes,
%%%   rebindLimiting for the limiting class) as elasticity -- and therefore mode
%%%   character -- changes. Three words are used precisely throughout, because
%%%   conflating them is what makes this subject hard to discuss:
%%%
%%%     BINDING          the one-time act, at construction, of deciding that measured
%%%                      point k corresponds to a particular calculated object.
%%%     ANCHOR           the stored, elasticity-INDEPENDENT description of identity that
%%%                      lets that correspondence be re-established later. For a Ritz
%%%                      mode it is the eigenvector (anchorEigVec_con); for a limiting
%%%                      root it is the branch index plus the polarisation 3-vector
%%%                      (limBranchIwt + limAnchorPol). An anchor is compared against,
%%%                      never searched for.
%%%     RE-IDENTIFICATION (the code calls it REBINDING) the act, repeated at every
%%%                      objective evaluation, of finding which of this iterate's
%%%                      calculated objects the anchor now refers to. A mode INDEX is
%%%                      worthless across iterations because the spectrum reorders;
%%%                      only the anchor survives.
%%%
%%%   CONFIDENCE is the score the winning candidate achieved at re-identification:
%%%   macValue for a Ritz-tracked point, limMatchScore for a limiting one. A collapsing
%%%   confidence is the observable signature of the tracking failure this whole
%%%   machinery exists to prevent -- see confidence() for the one accessor that knows
%%%   which signal belongs to which class.
%%%
%%%   Also owns the hand-build and pinning API: addPoint and setLimitingAnchor build an
%%%   association point by point (tests, scripts), setPinned/pinnedPoints mark a binding
%%%   as human-confirmed (see setPinned).
%%%
%%% Value class: rebind/rebindLimiting/addPoint/setPinned return the modified object. That is what lets ElasticityInversion re-bind into a SCRATCH
%%%   copy every objective evaluation and promote it only when the cost improved, so a
%%%   rejected trial step can never move an anchor.
%%%
%%% JOINT MULTI-DATASET INVERSION: one object per dataset (mode tracking is strictly a
%%%   per-dataset concern). For a joint inversion of several crystallographic cuts of the
%%%   SAME crystal, pass an ARRAY of ModeAssociation (one per dataset, each built against
%%%   its OWN forward sweep at the SAME shared initial-guess elasticity) to
%%%   SensitivityAnalysis.calc/ElasticityInversion -- see those classes' docblocks.
%%%
%%% GuliweRR v6.0

    properties (Constant)
        %%% Shape constants of the limiting-root re-identification score (see
        %%% limitingRootScore for the formula and the reasons)
        LIM_VELOCITY_WINDOW = 0.15   % Gaussian width for velocity proximity, RELATIVE to
                                       % the anchored root's own velocity
        LIM_THETA_WINDOW    = 0.15   % Gaussian width for sagittal-root-angle proximity [rad]
        LIM_BRANCH_PENALTY  = 0.5    % multiplier when the candidate root's branch index
                                       % differs from the anchor's (a soft preference)

        %%% The binding classes tracked through the Ritz eigenproblem (rebind), as
        %%% opposed to 'limiting' (rebindLimiting). Every per-class selection in
        %%% ModeAssociation/SensitivityAnalysis/ElasticityInversion reads this set (or
        %%% activeRitzClasses): the sensitivity rows and the residual rows must select
        %%% the same points.
        RITZ_CLASSES = {'ritzMode', 'pseudoSAW'}

        %%% Default re-identification thresholds (defaultTrackOpts) and rebind window,
        %%% the one definition ElasticityInversion's property defaults also read
        DEFAULT_TRACK_OPTS = struct( ...
            'macWarnThreshold',    0.7, ...
            'limMatchScoreThresh', 0.5, ...
            'useGreenSteering',    false, ...
            'greenWindowRel',      0.02, ...
            'greenFilterRel',      0.04)
        DEFAULT_REBIND_WINDOW = 1.0   % [km/s]
    end

    properties (SetAccess = protected)
        points struct = struct( ...
            'angleDeg',             {}, ...   % propagation angle [deg]
            'measuredVelocity',     {}, ...   % [km/s]
            'waveType',             {}, ...   % MeasuredVelocities.WAVE_TYPES member
            'ritzModeNo',           {}, ...   % column index into eigVecs_con; NaN for "limiting"
            'calculatedVelocity',   {}, ...   % [km/s], current best match
            'anchorEigVec_con',     {}, ...   % mode-identity eigenvector; [] for "limiting"
            'bindingClass',         {}, ...   % "ritzMode" | "limiting" | "pseudoSAW"
            'macValue',             {}, ...   % achieved MAC at last rebind; NaN until rebound
            'modeCollision',        {}, ...   % true if the last rebind could not uniquely match
            ...
            ... %%% --- limiting-root anchor and live state (all NaN/[] off the limiting class) ---
            'limBranchIwt',         {}, ...   % 1 = longitudinal, 2 = fast shear, 3 = slow shear;
                                        ...   %   LimitingVels.findRoots' branch index. PART OF
                                        ...   %   THE ANCHOR.
            'limAnchorPol',         {}, ...   % [3x1] polarisation at the tracked root. PART OF
                                        ...   %   THE ANCHOR; governed by anchorPolicy.
            'limRootTheta',         {}, ...   % sagittal root angle theta* [rad]. NOT part of the
                                        ...   %   anchor -- live state, refreshed on every
                                        ...   %   promoted evaluation, because it is the bracket
                                        ...   %   centre for the next local refinement.
            'limMatchScore',        {}, ...   % confidence at the last re-identification, [0,1];
                                        ...   %   the limiting channel's counterpart to macValue
            'limRootLost',          {}, ...   % no candidate root cleared the threshold -- the
                                        ...   %   tracked feature has ceased to exist. The
                                        ...   %   counterpart to modeCollision.
            ...
            ... %%% --- Green's-function-derived confidence (Ritz classes) ---
            'anchorGreenPeakVel',   {}, ...   % Green peak velocity this point was bound against
                                        ...   %   at CONSTRUCTION [km/s]; NaN if unavailable
            'greenPeakVel',         {}, ...   % LIVE Green peak velocity at the last rescan
                                        ...   %   [km/s]; NaN disables Green steering for this
                                        ...   %   point. Seeded from anchorGreenPeakVel.
            'greenAgreement',       {}, ...   % |v_matched - greenPeakVel| / greenPeakVel at the
                                        ...   %   last re-identification; NaN if no live peak
            ...
            ... %%% --- human confirmation and per-evaluation bookkeeping ---
            'pinned',               {}, ...   % human-confirmed ANCHOR: never replace it, give
                                        ...   %   this point first claim when candidates are
                                        ...   %   contested, and never exclude it for low
                                        ...   %   confidence. Pinning fixes IDENTITY, not values.
            'lowConfidenceCount',   {}, ...   % how many PROMOTED evaluations this point spent
                                        ...   %   below its class's confidence threshold
            'activeThisEvaluation', {})       % false => this evaluation's residual AND Jacobian
                                              %   row are both zeroed (never deleted: lsqnonlin
                                              %   requires a constant residual length)
    end

    methods
        function obj = ModeAssociation()
            %%% Initialise an empty ModeAssociation
            %%% >> obj = ModeAssociation()
        end

        function n = nPoints(obj)
            %%% Number of bound points
            %%% >> n = nPoints(obj)
            n = numel(obj.points);
        end

        function angles = uniqueRitzAngles(obj)
            %%% Unique propagation angles among Ritz-mode-bound points ("ritzMode"/"pseudoSAW"),
            %%% i.e. the angles at which the inversion loop must solve an eigenproblem.
            %%% >> angles = uniqueRitzAngles(obj)
            mask = ismember({obj.points.bindingClass}, ModeAssociation.RITZ_CLASSES);
            angles = unique([obj.points(mask).angleDeg]);
        end

        function angles = uniqueLimitingAngles(obj)
            %%% Unique propagation angles among "limiting"-bound points, i.e. the angles
            %%% at which the inversion loop must compute a Christoffel root set.
            %%% >> angles = uniqueLimitingAngles(obj)
            if isempty(obj.points)
                angles = [];
                return
            end
            mask = strcmp({obj.points.bindingClass}, 'limiting');
            angles = unique([obj.points(mask).angleDeg]);
        end

        function idx = pointsAtAngle(obj, angleDeg, bindingClasses)
            %%% Indices of points at a given angle, optionally filtered by bindingClass
            %%% >> idx = pointsAtAngle(obj, angleDeg, bindingClasses)
            %%% INPUT: (opt) bindingClasses .. cellstr subset of {"ritzMode","limiting",
            %%%                                "pseudoSAW"}; [] = no filter
            if isempty(obj.points)
                idx = [];
                return
            end
            if nargin < 3 || isempty(bindingClasses)
                mask = true(1, numel(obj.points));
            else
                mask = ismember({obj.points.bindingClass}, bindingClasses);
            end
            idx = find(abs([obj.points.angleDeg] - angleDeg) < 1e-9 & mask);
        end

        function c = confidence(obj)
            %%% The class-appropriate confidence signal for every point, as ONE vector
            %%% >> c = confidence(obj)
            %%% OUTPUT: c .. [nPoints x 1]; macValue for "ritzMode"/"pseudoSAW",
            %%%              limMatchScore for "limiting", NaN where never re-identified
            %%%
            %%% Exists so that reporting and plotting code does not have to know the
            %%% class taxonomy: reading macValue directly shows nothing at all for
            %%% limiting points, which are not MAC-tracked and carry macValue = NaN
            %%% forever.
            n = numel(obj.points);
            c = NaN(n, 1);
            for k = 1:n
                if strcmp(obj.points(k).bindingClass, 'limiting')
                    c(k) = obj.points(k).limMatchScore;
                else
                    c(k) = obj.points(k).macValue;
                end
            end
        end

        function idx = pinnedPoints(obj)
            %%% Indices of human-confirmed (pinned) points, for reporting
            %%% >> idx = pinnedPoints(obj)
            if isempty(obj.points)
                idx = [];
                return
            end
            idx = find([obj.points.pinned]);
        end

        function obj = setPinned(obj, idx, tf)
            %%% Mark one or more points as human-confirmed (or un-confirm them)
            %%% >> obj = setPinned(obj, idx, tf)
            %%%
            %%% INPUT: 1. idx .. point index/indices
            %%%  (opt) 2. tf  .. logical (default true)
            %%%
            %%% pinned = true means: THIS POINT'S ANCHOR IS HUMAN-CONFIRMED. Never
            %%% replace it, and give the point priority when candidates are contested.
            %%% Re-identify it every iteration exactly like any other point -- that is
            %%% what keeps it on the mode the human chose. It does NOT mean "freeze
            %%% ritzModeNo": the spectrum reorders at every iterate, so mode #37 at
            %%% iteration 1 and mode #37 at iteration 40 are different physical modes,
            %%% and freezing the index would guarantee the drift this whole scheme
            %%% exists to prevent. Live state (ritzModeNo, calculatedVelocity,
            %%% macValue/limMatchScore, limRootTheta) is updated on a pinned point
            %%% exactly as on any other.
            if nargin < 3 || isempty(tf); tf = true; end
            validateattributes(idx, {'numeric'}, {'vector', 'positive', 'integer', ...
                '<=', numel(obj.points)}, 'ModeAssociation.setPinned', 'idx');
            for k = idx(:).'
                obj.points(k).pinned = logical(tf);
            end
        end

        function obj = setLimitingAnchor(obj, idx, iwt, pol, thetaRoot)
            %%% Set one limiting point's ROOT ANCHOR by hand
            %%% >> obj = setLimitingAnchor(obj, idx, iwt, pol, thetaRoot)
            %%%
            %%% INPUT: 1. idx       .. point index (must be bindingClass "limiting")
            %%%        2. iwt       .. branch index (1 longitudinal, 2 fast shear, 3 slow shear)
            %%%        3. pol       .. [3x1] polarisation at the root
            %%%        4. thetaRoot .. sagittal root angle theta* [rad]
            %%%
            %%% The limiting-class counterpart of addPoint's anchorEigVec_con argument,
            %%% kept as a separate setter rather than four more positional arguments on
            %%% addPoint. fromSweep fills these itself; this is the hand-build path.
            validateattributes(idx, {'numeric'}, {'scalar', 'positive', 'integer', ...
                '<=', numel(obj.points)}, 'ModeAssociation.setLimitingAnchor', 'idx');
            validateattributes(iwt, {'numeric'}, {'scalar', 'integer', '>=', 1, '<=', 3}, ...
                'ModeAssociation.setLimitingAnchor', 'iwt');
            validateattributes(pol, {'double'}, {'vector', 'numel', 3, 'real', 'finite'}, ...
                'ModeAssociation.setLimitingAnchor', 'pol');
            validateattributes(thetaRoot, {'double'}, {'scalar', 'real', 'finite'}, ...
                'ModeAssociation.setLimitingAnchor', 'thetaRoot');
            if ~strcmp(obj.points(idx).bindingClass, 'limiting')
                error('ModeAssociation:setLimitingAnchor:wrongClass', ...
                    'Point %d has bindingClass "%s"; a root anchor only means anything on "limiting".', ...
                    idx, obj.points(idx).bindingClass);
            end
            obj.points(idx).limBranchIwt = iwt;
            obj.points(idx).limAnchorPol = pol(:);
            obj.points(idx).limRootTheta = thetaRoot;
            obj.points(idx).limRootLost  = false;
        end

        function obj = setLimitingRefined(obj, idx, velocity, thetaRoot)
            %%% Advance ONE limiting point from a successful anchored local refinement,
            %%% without a full re-identification
            %%% >> obj = setLimitingRefined(obj, idx, velocity, thetaRoot)
            %%%
            %%% The limitingRescanCadence > 1 shortcut ONLY (see
            %%% ElasticityInversion.trackLimitingRoots). It is legitimate to skip
            %%% rebindLimiting here precisely because LimitingVels.refineRoot already
            %%% applied findRoots' own accept test inside the anchored bracket and
            %%% returned ok = true: the tracked root demonstrably still exists, at a
            %%% theta* within a few degrees of where it was. limMatchScore is left at its
            %%% last re-identified value rather than being invented, since no scoring
            %%% happened on this evaluation.
            obj.points(idx).calculatedVelocity   = velocity;
            obj.points(idx).limRootTheta         = thetaRoot;
            obj.points(idx).limRootLost          = false;
            obj.points(idx).activeThisEvaluation = true;
        end

        function obj = setLimitingLost(obj, idx)
            %%% Flag ONE limiting point's tracked root as lost at this evaluation
            %%% >> obj = setLimitingLost(obj, idx)
            %%%
            %%% Written by ElasticityInversion when a perturbed evaluation of the point's
            %%% finite-difference Jacobian row lost the root
            %%% (SensitivityAnalysis.channelDerivative's ok = false): the row cannot be
            %%% formed, so the point is excluded from this evaluation exactly as
            %%% rebindLimiting excludes a root that has vanished.
            obj.points(idx).limRootLost          = true;
            obj.points(idx).activeThisEvaluation = false;
            obj.points(idx).lowConfidenceCount   = obj.points(idx).lowConfidenceCount + 1;
        end

        function obj = setGreenPeak(obj, idx, peakVel)
            %%% Store a freshly measured live Green peak velocity on one point
            %%% >> obj = setGreenPeak(obj, idx, peakVel)
            %%%
            %%% INPUT: 1. idx     .. point index
            %%%        2. peakVel .. Green peak velocity at the CURRENT iterate [km/s],
            %%%                     or NaN to disable Green steering for this point
            %%%
            %%% Written only by ElasticityInversion's Green rescan. It NEVER sets
            %%% calculatedVelocity: a Green peak is a velocity, not an eigenvector, so it
            %%% cannot be an anchor and is not the quantity being fitted. It only steers
            %%% which Ritz mode the MAC match selects, which keeps the residual on the Ritz
            %%% eigenvalue and therefore keeps the analytic Rayleigh-quotient Jacobian
            %%% exactly valid.
            obj.points(idx).greenPeakVel = peakVel;
        end

        function obj = addPoint(obj, angleDeg, ritzModeNo, measuredVelocity, calculatedVelocity, ...
                anchorEigVec_con, bindingClass, waveType, pinned)
            %%% Add one measured-point <-> calculated-object binding
            %%% >> obj = addPoint(obj, angleDeg, ritzModeNo, measuredVelocity, ...
            %%%                   calculatedVelocity, anchorEigVec_con, bindingClass, ...
            %%%                   waveType, pinned)
            %%% INPUT: (opt) bindingClass .. default "ritzMode"
            %%%        (opt) waveType     .. default "unknown"
            %%%        (opt) pinned       .. human-confirmed anchor (default false); true
            %%%                              when a human hand-picks the mode, e.g. off an
            %%%                              angular map.
            if nargin < 7 || isempty(bindingClass); bindingClass = 'ritzMode'; end
            if nargin < 8 || isempty(waveType);     waveType     = 'unknown';  end
            if nargin < 9 || isempty(pinned);       pinned       = false;      end
            bindingClass = char(bindingClass);
            waveType     = char(waveType);
            if ~ismember(bindingClass, [ModeAssociation.RITZ_CLASSES, {'limiting'}])
                error('ModeAssociation:addPoint:badBindingClass', ...
                    'bindingClass must be ''ritzMode'', ''pseudoSAW'' or ''limiting''; got ''%s''.', ...
                    bindingClass);
            end
            if ~ismember(waveType, MeasuredVelocities.WAVE_TYPES)
                error('ModeAssociation:addPoint:badWaveType', '"%s" is not one of: %s.', ...
                    waveType, strjoin(MeasuredVelocities.WAVE_TYPES, ', '));
            end

            newPt = ModeAssociation.emptyPoint();
            newPt.angleDeg           = angleDeg;
            newPt.measuredVelocity   = measuredVelocity;
            newPt.waveType           = waveType;
            newPt.ritzModeNo         = ritzModeNo;
            newPt.calculatedVelocity = calculatedVelocity;
            newPt.anchorEigVec_con   = anchorEigVec_con;
            newPt.bindingClass       = bindingClass;
            newPt.pinned             = logical(pinned);
            obj.points(end+1) = newPt;
        end

        function obj = removeDuplicatePoints(obj)
            %%% Keep at most one binding per (angle, ritzModeNo) for "ritzMode"/"pseudoSAW"
            %%% points, and at most one per (angle, waveType) for "limiting" points -- the
            %%% closest-matching duplicate wins, the rest are dropped with a warning.
            %%% >> obj = removeDuplicatePoints(obj)
            if numel(obj.points) < 2
                return
            end

            keyOf = cell(numel(obj.points), 1);
            for k = 1:numel(obj.points)
                p = obj.points(k);
                if strcmp(p.bindingClass, 'limiting')
                    keyOf{k} = sprintf('%.6f|%s', p.angleDeg, p.waveType);
                else
                    keyOf{k} = sprintf('%.6f|%d', p.angleDeg, p.ritzModeNo);
                end
            end
            [uniqueKeys, ~, groupId] = unique(keyOf);

            keep = true(numel(obj.points), 1);
            nDropped = 0;
            for g = 1:numel(uniqueKeys)
                members = find(groupId == g);
                if numel(members) < 2
                    continue
                end
                errs = abs([obj.points(members).calculatedVelocity] - [obj.points(members).measuredVelocity]);
                [~, iBest] = min(errs);
                dropIdx = members(members ~= members(iBest));
                keep(dropIdx) = false;
                nDropped = nDropped + numel(dropIdx);
            end

            if nDropped > 0
                warning('ModeAssociation:removeDuplicatePoints:dropped', ...
                    '%d duplicate point(s) removed (same angle+mode bound more than once).', nDropped);
            end
            obj.points = obj.points(keep);
        end

        function obj = rebind(obj, angleDeg, eigVecs_con, M_con, vels, rebindVelocityWindow, ...
                updateAnchor, trackOpts)
            %%% Re-establish the RITZ-class mode binding at one angle from a freshly
            %%% solved eigenproblem
            %%% >> obj = rebind(obj, angleDeg, eigVecs_con, M_con, vels, rebindVelocityWindow, updateAnchor, trackOpts)
            %%%
            %%% INPUT: 1. angleDeg     .. propagation angle [deg]
            %%%        2. eigVecs_con  .. [nCon x nModes] M-orthonormal eigenvectors, current iterate
            %%%        3. M_con        .. [nCon x nCon] mass matrix, constrained basis
            %%%        4. vels         .. [nModes x 1] phase velocities of eigVecs_con's columns
            %%%  (opt)  5. rebindVelocityWindow .. only candidates within this window [km/s] of
            %%%                          each point's window CENTRE are considered (default 1.0)
            %%%  (opt)  6. updateAnchor .. if true, replace anchorEigVec_con with the newly
            %%%                          matched eigenvector ("chained" anchor policy); if
            %%%                          false, anchorEigVec_con is left untouched ("fixed")
            %%%                          (default false)
            %%%  (opt)  7. trackOpts    .. struct, see ModeAssociation.defaultTrackOpts
            %%%
            %%% M-weighted MAC: MAC_m = |alpha_m^H M_con alpha_anchor|. The inner product is
            %%% M-weighted, not Euclidean: because the eigenvectors are M-orthonormal,
            %%% alpha_m^H M_con alpha_anchor is exactly the expansion coefficient of the
            %%% anchor in the current eigenbasis and sum_m |coef|^2 = 1, so the MAC value is
            %%% a genuine fraction-of-identity. (DispersionBranches.track uses the same
            %%% mass weighting, with each wavelength step's own M.)
            %%%
            %%% GREEN STEERING (pseudoSAW points only, trackOpts.useGreenSteering). A leaky
            %%% pseudo-surface wave radiates real shear energy into the half-space, and
            %%% Ritz's finite, non-radiating domain represents that resonance as a CLUSTER
            %%% of near-degenerate modes whose eigenvector COMPOSITION is very sensitive to
            %%% small elastic-constant changes. The eigenVALUE is precise; the mode SHAPE,
            %%% and hence the MAC, is what is fragile. When a live Green peak is available
            %%% the velocity window is re-centred on IT rather than on the previous accepted
            %%% calculatedVelocity (which drags the window along with a mis-track), the MAC
            %%% is restricted to candidates within trackOpts.greenFilterRel of that peak,
            %%% with a Gaussian of width trackOpts.greenWindowRel breaking ties inside it.
            %%% Green never sets calculatedVelocity: it only decides WHICH mode is matched,
            %%% so the residual stays on the Ritz eigenvalue and the analytic Jacobian stays
            %%% exactly valid. The MAC alone cannot detect this failure: a chained anchor
            %%% follows a mis-track onto a neighbouring cluster member and then agrees with
            %%% it perfectly (macValue = 1).
            %%%
            %%% PINNED points take their own argmax BEFORE matchpairs runs on the rest, so
            %%% an unpinned point with a marginally higher score cannot displace a human's
            %%% confirmed binding onto its second choice.

            if nargin < 6 || isempty(rebindVelocityWindow)
                rebindVelocityWindow = ModeAssociation.DEFAULT_REBIND_WINDOW;
            end
            if nargin < 7 || isempty(updateAnchor);         updateAnchor         = false; end
            if nargin < 8; trackOpts = []; end
            trackOpts = ModeAssociation.mergeTrackOpts(trackOpts);

            idxPts = obj.pointsAtAngle(angleDeg, ModeAssociation.RITZ_CLASSES);
            if isempty(idxPts)
                return
            end

            nModes = numel(vels);
            nPts   = numel(idxPts);
            MAC    = zeros(nModes, nPts);
            score  = zeros(nModes, nPts);
            steered = false(1, nPts);
            for k = 1:nPts
                p = obj.points(idxPts(k));
                useGreen = trackOpts.useGreenSteering && strcmp(p.bindingClass, 'pseudoSAW') ...
                    && isfinite(p.greenPeakVel) && p.greenPeakVel > 0;
                steered(k) = useGreen;
                if useGreen
                    centre = p.greenPeakVel;
                else
                    centre = p.calculatedVelocity;
                end
                inWindow = abs(vels - centre) < rebindVelocityWindow;
                if useGreen
                    %%% HARD candidate filter at the GROSS-ERROR scale, with the MAC left
                    %%% free to choose inside it. Both halves are load-bearing:
                    %%%
                    %%%  - A filter, not a weight, because the mechanism being fought is
                    %%%    HYSTERESIS. One evaluation whose match lands on a neighbouring
                    %%%    cluster member re-seats the chained anchor onto it, after which
                    %%%    that member scores MAC 1.0 forever and the correct mode scores
                    %%%    almost nothing. No realistic Gaussian WEIGHT can pull that back.
                    %%%  - But the filter must NOT be tight enough to re-select the mode on
                    %%%    the fine scale. Which member of a near-degenerate cluster "is" a
                    %%%    pSAW measurement is genuinely ambiguous at the mode spacing, and
                    %%%    the answer depends on the selection RULE: the forward classifier
                    %%%    picks by its identity score (ModeClassification.matchRitzToGreenPeak)
                    %%%    while this filter can only see velocity, and re-selecting on the
                    %%%    fine scale biases every pSAW residual the same way. The MAC chain,
                    %%%    seeded at construction from the classifier's pick, carries that
                    %%%    choice forward unbiased; Green's job is only to stop it wandering
                    %%%    off the resonance altogether.
                    %%%
                    %%% Emptying the candidate set is SAFE, and intended: it raises
                    %%% modeCollision, which ElasticityInversion answers with an immediate
                    %%% single-angle Green rescan and one retry, so a stale peak is corrected
                    %%% rather than obeyed.
                    inWindow = inWindow & ...
                        abs(vels - p.greenPeakVel) <= trackOpts.greenFilterRel * p.greenPeakVel;
                end
                if ~any(inWindow)
                    continue
                end
                macCol = abs(eigVecs_con(:, inWindow)' * M_con * p.anchorEigVec_con);
                MAC(inWindow, k) = macCol;
                if useGreen
                    %%% inside the filter, the Gaussian only breaks ties
                    w = trackOpts.greenWindowRel * p.greenPeakVel;
                    score(inWindow, k) = macCol .* exp(-((vels(inWindow) - p.greenPeakVel) ./ w).^2);
                else
                    score(inWindow, k) = macCol;
                end
            end

            assignedMode = ModeAssociation.assignWithPinnedPriority(score, [obj.points(idxPts).pinned]);

            for k = 1:nPts
                gi = idxPts(k);
                p  = obj.points(gi);
                m  = assignedMode(k);

                matched = ~isnan(m) && score(m, k) > 0;

                if ~matched
                    p.macValue      = 0;
                    p.modeCollision = true;
                    %%% ritzModeNo/calculatedVelocity are deliberately LEFT at their previous
                    %%% values (there is no better estimate), but the row is zeroed so that a
                    %%% stale velocity can never be paired with a current-eigenvector Jacobian
                    %%% row.
                    p.lowConfidenceCount   = p.lowConfidenceCount + 1;
                    p.activeThisEvaluation = p.pinned;   % a human said this binding is right
                else
                    p.ritzModeNo         = m;
                    p.calculatedVelocity = vels(m);
                    p.macValue           = MAC(m, k);
                    p.modeCollision      = false;
                    if isfinite(p.greenPeakVel) && p.greenPeakVel > 0
                        p.greenAgreement = abs(vels(m) - p.greenPeakVel) / p.greenPeakVel;
                    end
                    if steered(k)
                        %%% For a STEERED point the MAC is a diagnostic, not the gate.
                        %%% greenAgreement is the confidence signal here -- and it is already
                        %%% bounded by greenFilterRel, since the filter above admitted no
                        %%% candidate further away. When Green moves the match OFF the mode
                        %%% the drifted anchor preferred -- the whole purpose of steering --
                        %%% the MAC against that stale anchor legitimately drops. Gating on
                        %%% it would zero exactly the rows the correction just rescued and
                        %%% freeze the anchor on the wrong cluster member permanently.
                        lowConf = false;
                    else
                        lowConf = p.macValue < trackOpts.macWarnThreshold;
                    end
                    if lowConf
                        p.lowConfidenceCount   = p.lowConfidenceCount + 1;
                        p.activeThisEvaluation = p.pinned;
                    else
                        p.activeThisEvaluation = true;
                    end
                    %%% A CONFIDENT match is what earns the right to redefine identity: never
                    %%% let a low-confidence one update the anchor, even on a cost-improving
                    %%% evaluation, and never update a pinned point's anchor at all.
                    if updateAnchor && ~p.pinned && ~lowConf
                        p.anchorEigVec_con = eigVecs_con(:, m);
                    end
                end
                obj.points(gi) = p;
            end
        end

        function obj = rebindLimiting(obj, angleDeg, roots, updateAnchor, trackOpts)
            %%% Re-establish the LIMITING-class root binding at one angle from a freshly
            %%% computed Christoffel root set
            %%% >> obj = rebindLimiting(obj, angleDeg, roots, updateAnchor, trackOpts)
            %%%
            %%% INPUT: 1. angleDeg  .. propagation angle [deg]
            %%%        2. roots     .. struct array from SensitivityAnalysis.limitingRootsAt
            %%%                       (fields iwt, vel, pol, theta), current iterate
            %%%  (opt) 3. updateAnchor .. replace limAnchorPol with the newly matched root's
            %%%                       polarisation ("chained"); limBranchIwt and limRootTheta
            %%%                       are updated regardless (default false)
            %%%  (opt) 4. trackOpts .. struct, see ModeAssociation.defaultTrackOpts
            %%%
            %%% Deliberately a sibling of rebind: same shape, same matchpairs step, same
            %%% pinned-priority pass, same promote-discipline contract. What differs is only
            %%% the score -- there is no eigenvector here, so identity is carried by the
            %%% polarisation 3-vector plus the branch index and the sagittal root angle
            %%% (limitingRootScore).
            %%%
            %%% Limiting roots can VANISH, not merely cross: over a small change of the
            %%% elastic constants a shear root can cease to exist, and an untracked
            %%% max()/min() over the shear roots then makes the calculated velocity jump
            %%% to another root with no warning, with a finite-difference Jacobian taken
            %%% straight across the discontinuity. A remedy designed only against crossings
            %%% would not catch it; the polarisation overlap does, because it collapses at
            %%% exactly that step.

            if nargin < 4 || isempty(updateAnchor); updateAnchor = false; end
            if nargin < 5; trackOpts = []; end
            trackOpts = ModeAssociation.mergeTrackOpts(trackOpts);

            idxPts = obj.pointsAtAngle(angleDeg, {'limiting'});
            if isempty(idxPts)
                return
            end

            nRoots = numel(roots);
            nPts   = numel(idxPts);
            score  = zeros(nRoots, nPts);
            for k = 1:nPts
                if nRoots > 0
                    score(:, k) = ModeAssociation.limitingRootScore(obj.points(idxPts(k)), roots);
                end
            end

            assignedRoot = ModeAssociation.assignWithPinnedPriority(score, [obj.points(idxPts).pinned]);

            for k = 1:nPts
                gi = idxPts(k);
                p  = obj.points(gi);
                q  = assignedRoot(k);

                if isnan(q) || score(q, k) <= trackOpts.limMatchScoreThresh
                    if isnan(q)
                        p.limMatchScore = 0;
                    else
                        p.limMatchScore = score(q, k);
                    end
                    p.limRootLost          = true;
                    p.lowConfidenceCount   = p.lowConfidenceCount + 1;
                    p.activeThisEvaluation = p.pinned;
                    %%% calculatedVelocity/limRootTheta/limBranchIwt deliberately left at
                    %%% their previous values -- the row is zeroed instead, which keeps the
                    %%% residual and the Jacobian row consistent with each other. A point
                    %%% whose tracked root has vanished must NOT silently take the best
                    %%% remaining root: that is exactly the jump described above.
                else
                    p.calculatedVelocity   = roots(q).vel;
                    p.limMatchScore        = score(q, k);
                    p.limRootTheta         = roots(q).theta;
                    p.limBranchIwt         = roots(q).iwt;
                    p.limRootLost          = false;
                    p.activeThisEvaluation = true;
                    if updateAnchor && ~p.pinned
                        p.limAnchorPol = roots(q).pol(:);
                    end
                end
                obj.points(gi) = p;
            end
        end
    end

    methods (Static)
        function obj = fromSweep(sweep, measured, velocityTolerance, angleTolerance, invertedLayerNo)
            %%% Build the initial mode association from a forward sweep at the initial guess
            %%% >> obj = ModeAssociation.fromSweep(sweep, measured, velocityTolerance, angleTolerance, invertedLayerNo)
            %%%
            %%% INPUT: 1. sweep             .. CalcAngSweepHS/CalcAngSweepLayered result,
            %%%                                retentionLevel='full', calcMode='full'
            %%%        2. measured          .. MeasuredVelocities
            %%%  (opt) 3. velocityTolerance .. max |v_calc - v_meas| to accept a candidate
            %%%                                [km/s] (default 0.3)
            %%%  (opt) 4. angleTolerance    .. max |angle_meas - angle_sweep| to consider the
            %%%                                same angle [deg] (default 1e-6); errors if a
            %%%                                measured angle has no sweep angle within it
            %%%  (opt) 5. invertedLayerNo   .. which layer's material the LIMITING channel
            %%%                                refers to (default 1). For a half-space there
            %%%                                is only one layer, so the default is always
            %%%                                right there; it exists so the limiting binding
            %%%                                below uses the SAME material
            %%%                                SensitivityAnalysis/ElasticityInversion
            %%%                                evaluate the limiting channel with.
            %%%
            %%% Wave-type routing: see ModeAssociation.classifyWaveType.
            %%%
            %%% The two classes bind DIFFERENTLY:
            %%%   Ritz classes -- filter sweep.significantModesMap by the accepted
            %%%     ModeClassification types, take the entry minimising |vel - measured|,
            %%%     reject beyond velocityTolerance. Also stores the classifier's Green
            %%%     peak velocity as the seed of the Green steering (anchorGreenPeakVel).
            %%%   Limiting -- compute the Christoffel root set at that angle directly, filter
            %%%     by iwt == waveTypeIndex(waveType), take the nearest by velocity, reject
            %%%     beyond velocityTolerance. Routing these through significantModesMap
            %%%     instead would (a) let 'longitudinal'/'fastShear'/'slowShear' all compete
            %%%     against EVERY 'limiting' entry at that angle with no branch filter, and
            %%%     (b) apply velocityTolerance to a nearby RITZ MODE'S velocity rather than
            %%%     to the Christoffel root that is actually fitted downstream.

                if nargin < 3 || isempty(velocityTolerance); velocityTolerance = 0.3;  end
                if nargin < 4 || isempty(angleTolerance);    angleTolerance    = 1e-6; end
                if nargin < 5 || isempty(invertedLayerNo);   invertedLayerNo   = 1;    end

                obj = ModeAssociation();
                nDropped = 0;
                nNoGreen = 0;   % pseudoSAW points bound without a Green peak

                matLay = sweep.structure.mater(invertedLayerNo);
                symmNo = matLay.elConsts.symmClassNo;
                cVec0  = matLay.elConsts.vector(:);
                limRootsCache = cell(numel(sweep.angles), 1);

                for i = 1:measured.nPoints
                    angleDeg = measured.angleDeg(i);
                    [dAngle, iAngle] = min(abs(sweep.angles - angleDeg));
                    if dAngle > angleTolerance
                        error('ModeAssociation:fromSweep:angleNotFound', ...
                            ['Measured angle %.6g deg has no sweep angle within tolerance ', ...
                             '%.2g deg (closest sweep angle: %.6g deg).'], ...
                            angleDeg, angleTolerance, sweep.angles(iAngle));
                    end

                    wtMeas = measured.waveType{i};
                    [acceptedTypes, bindingClass] = ModeAssociation.classifyWaveType(wtMeas);

                    newPt = ModeAssociation.emptyPoint();
                    newPt.angleDeg         = sweep.angles(iAngle);
                    newPt.measuredVelocity = measured.velocity(i);
                    newPt.waveType         = wtMeas;
                    newPt.bindingClass     = bindingClass;

                    if strcmp(bindingClass, 'limiting')
                        if isempty(limRootsCache{iAngle})
                            limRootsCache{iAngle} = SensitivityAnalysis.limitingRootsAt( ...
                                matLay, symmNo, cVec0, sweep.angles(iAngle));
                        end
                        roots = limRootsCache{iAngle};
                        iwtWanted = SensitivityAnalysis.waveTypeIndex(wtMeas);
                        cand = roots([roots.iwt] == iwtWanted);
                        if isempty(cand)
                            %%% The branch can have no root at this angle at all; never
                            %%% fall back to a root of the OTHER branch.
                            warning('ModeAssociation:fromSweep:noLimitingRoot', ...
                                ['No %s limiting root exists at %.6g deg; the point was ', ...
                                 'dropped rather than bound to another branch.'], wtMeas, angleDeg);
                            nDropped = nDropped + 1;
                            continue
                        end
                        [dv, iBest] = min(abs([cand.vel] - measured.velocity(i)));
                        if dv > velocityTolerance
                            nDropped = nDropped + 1;
                            continue
                        end
                        newPt.calculatedVelocity = cand(iBest).vel;
                        newPt.ritzModeNo         = NaN;
                        newPt.anchorEigVec_con   = [];
                        newPt.limBranchIwt       = cand(iBest).iwt;
                        newPt.limAnchorPol       = cand(iBest).pol(:);
                        newPt.limRootTheta       = cand(iBest).theta;
                        newPt.limMatchScore      = NaN;
                        newPt.limRootLost        = false;
                    else
                        candidates = sweep.significantModesMap{iAngle};
                        if ~isempty(acceptedTypes)
                            candidates = candidates(ismember({candidates.type}, acceptedTypes));
                        end
                        if isempty(candidates)
                            nDropped = nDropped + 1;
                            continue
                        end
                        %%% modeVel: the bound Ritz mode's eigenvelocity, what the inversion
                        %%% fits (an entry's vel is the Green peak)
                        candVels = ModeClassification.modeVels(candidates);
                        [dv, iBest] = min(abs(candVels - measured.velocity(i)));
                        if dv > velocityTolerance
                            nDropped = nDropped + 1;
                            continue
                        end
                        best = candidates(iBest);

                        newPt.calculatedVelocity = candVels(iBest);
                        newPt.ritzModeNo         = best.modeNo;
                        newPt.anchorEigVec_con   = sweep.cs{iAngle}.ritzSolver.eigVecs_con(:, best.modeNo);

                        %%% An INDEPENDENT, Green-derived estimate of where the physical
                        %%% resonance sits, obtained without a modal decomposition and
                        %%% therefore not subject to the near-degenerate-cluster problem the
                        %%% MAC has.
                        newPt.anchorGreenPeakVel   = best.greenPeakVel;
                        newPt.greenPeakVel         = best.greenPeakVel;
                        if isfinite(best.greenPeakVel) && best.greenPeakVel > 0
                            newPt.greenAgreement = abs(candVels(iBest) - best.greenPeakVel) / best.greenPeakVel;
                        elseif strcmp(bindingClass, 'pseudoSAW')
                            nNoGreen = nNoGreen + 1;
                        end
                    end

                    obj.points(end+1) = newPt;
                end

                if nDropped > 0
                    warning('ModeAssociation:fromSweep:pointsDropped', ...
                        ['%d of %d measured point(s) had no matching mode within tolerance ', ...
                         'and were dropped.'], nDropped, measured.nPoints);
                end
                if nNoGreen > 0
                    %%% A sweep classified without Green (the 'ritz' method) has no peak to
                    %%% seed the steering with: ElasticityInversion's first rescan would then
                    %%% start from the point's own Ritz velocity, the mirror of the Ritz
                    %%% tracking that steering exists to avoid.
                    warning('ModeAssociation:fromSweep:noGreenPeak', ...
                        ['%d pseudoSAW point(s) were bound to entries without a Green peak ', ...
                         '(a sweep classified without Green, e.g. the ''ritz'' method): Green ', ...
                         'steering has no seed for them. Classify the sweep with ''green'' ', ...
                         'for a pseudo-SAW inversion.'], nNoGreen);
                end

                obj = obj.removeDuplicatePoints();
            end

        function opts = defaultTrackOpts()
            %%% Default re-identification thresholds shared by rebind/rebindLimiting
            %%% >> opts = ModeAssociation.defaultTrackOpts()
            %%%
            %%% OUTPUT: opts .. struct with fields
            %%%   macWarnThreshold    .. below this MAC a Ritz-class match is treated as
            %%%                          low-confidence: the row is zeroed and the anchor is
            %%%                          not updated (0.7)
            %%%   limMatchScoreThresh .. below this combined score a limiting-root match is
            %%%                          rejected outright as limRootLost (0.5): a vanished
            %%%                          root collapses the score far below it, a tracked one
            %%%                          stays near 1
            %%%   useGreenSteering    .. let a live Green peak steer pseudoSAW matches (false
            %%%                          here; ElasticityInversion sets it from
            %%%                          pseudoSAWRefinement)
            %%%   greenWindowRel      .. width of the tie-breaking Gaussian on relative
            %%%                          velocity distance from the live Green peak (0.02 --
            %%%                          ModeClassification.matchRitzToGreenPeak's own
            %%%                          velMatchTol, and the same exp(-(relDist/w)^2) idiom)
            %%%   greenFilterRel      .. HOW FAR from the live Green peak a candidate Ritz
            %%%                          mode may sit and still be eligible at all (0.04). A
            %%%                          HARD filter, set at the GROSS-ERROR scale: between a
            %%%                          correctly matched pSAW mode (well under 1% from the
            %%%                          peak) and a mis-track onto another cluster member.
            %%%                          Do NOT tighten it towards greenWindowRel: see rebind.
            opts = ModeAssociation.DEFAULT_TRACK_OPTS;
        end

        function score = limitingRootScore(pt, roots)
            %%% Re-identification score of every candidate root against ONE point's anchor
            %%% >> score = ModeAssociation.limitingRootScore(pt, roots)
            %%%
            %%% INPUT: 1. pt    .. one points entry of bindingClass "limiting"
            %%%        2. roots .. struct array from SensitivityAnalysis.limitingRootsAt
            %%% OUTPUT: score .. [nRoots x 1] in [0,1]
            %%%
            %%%   score = |pol_a' pol_k| / (|pol_a| |pol_k|)          polarisation overlap
            %%%         * exp(-((v_k - v_a) / (LIM_VELOCITY_WINDOW*v_a))^2)
            %%%         * exp(-((theta_k - theta_a) / LIM_THETA_WINDOW)^2)
            %%%         * (iwt_k == iwt_a ? 1 : LIM_BRANCH_PENALTY)
            %%%
            %%% The polarisation overlap is the dominant term and the two Gaussians are
            %%% tie-breakers -- the same multiplicative projection x Gaussian x Gaussian
            %%% idiom as ModeClassification.matchRitzToGreenPeak/associateLimitingMode: a
            %%% linear penalty is too weak, and a bare projection with only an "> eps"
            %%% floor picks whatever is nearest when nothing is genuinely near. The branch
            %%% penalty is a SOFT preference, not a hard filter: near an acoustic axis the
            %%% T1/T2 labels are interchangeable (the descending-eigenvalue ordering that
            %%% defines them is arbitrary where the branches touch), and a hard filter
            %%% would reject the correct root exactly there.
            %%% abs() on the overlap because an eigenvector's overall sign is
            %%% arbitrary -- which means a polarisation and its exact reverse score 1.0.
            %%% For a limiting root that is physically correct (the same standing feature),
            %%% but it would be WRONG if this score were ever reused for something
            %%% directional such as group-velocity orientation. Do not reuse it carelessly.
            n = numel(roots);
            score = zeros(n, 1);
            polA = pt.limAnchorPol(:);
            nA   = norm(polA);
            if nA == 0 || isempty(polA)
                return
            end
            vA   = pt.calculatedVelocity;
            thA  = pt.limRootTheta;
            iwtA = pt.limBranchIwt;
            wV   = ModeAssociation.LIM_VELOCITY_WINDOW * abs(vA);
            for q = 1:n
                polK = roots(q).pol(:);
                nK   = norm(polK);
                if nK == 0
                    continue
                end
                overlap = abs(polA' * polK) / (nA * nK);
                fVel    = exp(-((roots(q).vel - vA) / wV)^2);
                fTheta  = exp(-((roots(q).theta - thA) / ModeAssociation.LIM_THETA_WINDOW)^2);
                if roots(q).iwt == iwtA
                    fBranch = 1;
                else
                    fBranch = ModeAssociation.LIM_BRANCH_PENALTY;
                end
                score(q) = overlap * fVel * fTheta * fBranch;
            end
        end
    end

    methods (Static, Access = {?ModeAssociation, ?SensitivityAnalysis, ?ElasticityInversion})
        function classes = activeRitzClasses(includePseudoSAW)
            %%% The Ritz binding classes that contribute data rows
            %%% >> classes = ModeAssociation.activeRitzClasses(includePseudoSAW)
            %%%
            %%% OUTPUT: classes .. {'ritzMode'}, plus 'pseudoSAW' when includePseudoSAW
            %%%                    (the IncludePseudoSAW option of SensitivityAnalysis.calc
            %%%                    and ElasticityInversion)
            if includePseudoSAW
                classes = ModeAssociation.RITZ_CLASSES;
            else
                classes = ModeAssociation.RITZ_CLASSES(1);
            end
        end
    end

    methods (Static, Access = private)
        function pt = emptyPoint()
            %%% One fully-populated points entry with every field at its neutral default
            %%% >> pt = ModeAssociation.emptyPoint()
            %%%
            %%% EVERY field of a MATLAB struct array must be present on EVERY entry --
            %%% concatenation does NOT fill a missing one with [] -- so every construction
            %%% path (addPoint, fromSweep) must go through here rather than assembling its
            %%% own subset.
            pt = struct( ...
                'angleDeg',             NaN, ...
                'measuredVelocity',     NaN, ...
                'waveType',             'unknown', ...
                'ritzModeNo',           NaN, ...
                'calculatedVelocity',   NaN, ...
                'anchorEigVec_con',     [], ...
                'bindingClass',         'ritzMode', ...
                'macValue',             NaN, ...
                'modeCollision',        false, ...
                'limBranchIwt',         NaN, ...
                'limAnchorPol',         [], ...
                'limRootTheta',         NaN, ...
                'limMatchScore',        NaN, ...
                'limRootLost',          false, ...
                'anchorGreenPeakVel',   NaN, ...
                'greenPeakVel',         NaN, ...
                'greenAgreement',       NaN, ...
                'pinned',               false, ...
                'lowConfidenceCount',   0, ...
                'activeThisEvaluation', true);
        end

        function opts = mergeTrackOpts(userOpts)
            %%% Fill any field the caller omitted from defaultTrackOpts
            %%% >> opts = ModeAssociation.mergeTrackOpts(userOpts)
            opts = ModeAssociation.defaultTrackOpts();
            if isempty(userOpts)
                return
            end
            fn = fieldnames(userOpts);
            for k = 1:numel(fn)
                if ~isfield(opts, fn{k})
                    error('ModeAssociation:mergeTrackOpts:unknownField', ...
                        '"%s" is not a recognised tracking option.', fn{k});
                end
                opts.(fn{k}) = userOpts.(fn{k});
            end
        end

        function assigned = assignWithPinnedPriority(score, pinnedMask)
            %%% One-to-one candidate assignment, pinned points claiming first
            %%% >> assigned = ModeAssociation.assignWithPinnedPriority(score, pinnedMask)
            %%%
            %%% INPUT: 1. score      .. [nCandidates x nPoints], >= 0, 0 = ineligible
            %%%        2. pinnedMask .. [1 x nPoints] logical
            %%% OUTPUT: assigned .. [1 x nPoints], candidate row index or NaN
            %%%
            %%% matchpairs (Optimization Toolbox) with a non-assignment cost of 0 maximises
            %%% the total score while forbidding two points at one angle from claiming the
            %%% same candidate, and leaves an all-zero column unmatched. Pinned points are
            %%% removed from that problem first, by their own argmax over the still-free
            %%% candidates: without this, an unpinned point with a marginally higher score
            %%% could displace a human's confirmed binding onto its second choice, which is
            %%% precisely the automatic system silently overriding a human confirmation.
            [nCand, nPts] = size(score);
            assigned = NaN(1, nPts);
            if nCand == 0 || nPts == 0
                return
            end
            freeRows = true(nCand, 1);

            for k = find(pinnedMask)
                col = score(:, k);
                col(~freeRows) = 0;
                [bestScore, m] = max(col);
                if bestScore > 0
                    assigned(k) = m;
                    freeRows(m) = false;
                end
            end

            remaining = find(~pinnedMask);
            if isempty(remaining) || ~any(freeRows)
                return
            end
            rowMap   = find(freeRows);
            subScore = score(freeRows, remaining);
            assignment = matchpairs(-subScore, 0);   % maximise total score <=> minimise -score
            for a = 1:size(assignment, 1)
                assigned(remaining(assignment(a, 2))) = rowMap(assignment(a, 1));
            end
        end

        function [acceptedTypes, bindingClass] = classifyWaveType(waveType)
            %%% Map a MeasuredVelocities wave-type label onto ModeClassification `type`
            %%% names and a ModeAssociation bindingClass
            %%% >> [acceptedTypes, bindingClass] = ModeAssociation.classifyWaveType(waveType)
            %%%
            %%% acceptedTypes is meaningful for the Ritz classes only: the limiting class does
            %%% not route through significantModesMap at all (see fromSweep), so its
            %%% {'limiting'} entry is documentation, not a filter.
            switch waveType
                case 'rayleigh'
                    acceptedTypes = {'Rayleigh', 'SAW'};
                    bindingClass  = 'ritzMode';
                case {'sezawa', 'love'}
                    acceptedTypes = {'SAW'};
                    bindingClass  = 'ritzMode';
                case 'pseudoRayleigh'
                    acceptedTypes = {'pseudoRayleigh'};
                    bindingClass  = 'pseudoSAW';
                case {'longitudinal', 'fastShear', 'slowShear'}
                    acceptedTypes = {'limiting'};
                    bindingClass  = 'limiting';
                case 'unknown'
                    acceptedTypes = {};   % no type filter -- accept any candidate
                    bindingClass  = 'ritzMode';
                otherwise
                    error('ModeAssociation:classifyWaveType:badWaveType', ...
                        '"%s" is not a recognised wave type.', waveType);
            end
        end
    end
end
