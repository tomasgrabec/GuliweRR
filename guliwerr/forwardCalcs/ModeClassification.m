classdef ModeClassification
%%% ModeClassification — cross-links Ritz (+ optionally Christoffel, Green)
%%% results into a classified list of significant surface-wave modes, for
%%% either a homogeneous half-space (recognizeFull) or a layered
%%% film-on-substrate structure (recognizeLayered); calc() dispatches on
%%% structure.nLay > 1. Without Green (calcMode = 'ritzOnly') the Ritz-only
%%% path (recognizeRitzOnlyScore) serves both.
%%%
%%% EXPERIMENTAL: checked against the measured maps published in Grabec et al. 2024 and the
%%% UTGS paper (arXiv:2510.10696); the scores, their scales and the interface may change.
%%%
%%% Two separate questions. EXISTENCE (is there a surface wave an experiment would
%%% see here?) is decided by one score per peak: GreenScore (pole-fitted G13_eff,
%%% analytic gate + aperture dephasing) on the half-space Green path,
%%% GreenScoreLayered (the same pole fit with an own-period gate, a whole-stack band
%%% and a polarisation term) on the layered Green path, GreenSpotScore on the refined
%%% path, RitzScore (a late-time Ritz modal sum) on the Ritz-only path. IDENTITY
%%% (which Ritz mode represents it) comes from Ritz. The Rayleigh/pseudoRayleigh label (or, layered, 'SAW' plus
%%% aboveSubstrateST) is applied last, from the peak velocity. See
%%% docs/reference/architecture/forwardCalcs_modeClassification.md.
%%%
%%% Value class, default-constructed then calc() is called and the result
%%% reassigned. Distinct from RitzSolverSAW.diagnostics, which holds the solver's
%%% own raw per-mode scores; this class CONSUMES those (plus Christoffel/Green
%%% when available) and produces the cross-linked, classified result.
%%%
%%% GuliweRR v6.0

    properties (Constant)
        %%% Relative tolerance of the half-space Rayleigh/pseudoRayleigh label: a Green
        %%% peak below minShear x (1 + RAYLEIGH_LABEL_TOL) is 'Rayleigh'. Equal to the
        %%% Green damping eta of the existence solve (= SawSpectrum.ETA):
        %%% eta broadens every pole to a full width of about eta x v, and the branch
        %%% point at minShear pulls the |G| maximum as well, so a Green peak cannot be
        %%% placed on either side of minShear more finely than that. A Rayleigh wave
        %%% hugging minShear (synthetic cubic, Zener ratio 20, (112) cut) put its peak 0 to
        %%% 0.77 m/s ABOVE minShear at 18 angles (tolerance there ~2.2 m/s); the
        %%% nearest genuine pseudoRayleigh above minShear in the Ni, Fe3Al, and Cu
        %%% cuts computed with GuliweRR is 39.8 m/s (Cu(112), 84 deg; tolerance ~5.7 m/s).
        RAYLEIGH_LABEL_TOL = 3e-3

        %%% Relative tie window of matchRitzToGreenPeak: Ritz modes within
        %%% GREEN_TIE_REL x vGreen of the nearest mode's distance count as equally near,
        %%% and identityWeight decides between them. The size of the 'green' peak's own
        %%% placement error against a fine (0.1 m/s, eta 3e-4) G13_eff peak: median
        %%% 0.1-1 m/s, p95 0.4-5 m/s on the cuts and stacks computed with GuliweRR
        %%% (~1.4 m/s here at 2.8 km/s). Closer than that, which mode is nearer is noise.
        GREEN_TIE_REL = 5e-4

        %%% Orthogonality tolerance of degenerateLimitGroups (|cos| of the angle between
        %%% two limiting polarisations): LimitingVels.findRoots makes a degenerate group's
        %%% polarisations orthonormal to round-off (~1e-16); a mirror twin's are not.
        LIMIT_ORTH_TOL = 1e-6

        %%% Classification stamp: bumped whenever a classification rule or a stored
        %%% candidate field changes. A Calc result records the value it was classified
        %%% with (classificationVersion); a result with another value is stale: it still
        %%% loads and shows its maps and modes, but cannot be reclassified -- recompute it.
        CLASSIFICATION_VERSION = 3
    end

    properties (SetAccess = protected)
        mode(1,:)        char   = 'full'   % 'ritzOnly' | 'full' | 'layered' — which path produced this result
        method(1,:)      char   = 'green'  % 'green' | 'refined' | 'ritz' — the classification method

        significantModes struct = ModeClassification.emptySignificantModes()
            % type in {'Rayleigh','pseudoRayleigh','limiting','other'} for the
            % 'ritzOnly'/'full' (homogeneous HS) paths, or {'SAW','limiting','other'}
            % for the 'layered' path — no Love. See
            % docs/reference/architecture/forwardCalcs_modeClassification.md, "The
            % significantModes entries".
            % vel -- the velocity of the method that decided existence: greenPeakVel
            % ('green': the pole; 'refined': the spot peak), the supporting mode's
            % eigenvelocity ('ritz'), the associated mode's eigenvelocity (limiting).
            % modeVel -- the eigenvelocity of the Ritz mode modeNo, on every entry.
            % greenPeakVel/limVelAssoc are NaN when not applicable to an entry.
            % aboveSubstrateST is true/false only on 'layered' SAW entries (velocity
            % above/below the substrate shear limiting velocity — informational, not
            % a type split) and NaN everywhere else. Every path must emit the SAME
            % field set: MATLAB struct-array concatenation ([sm1, sm2], sm(end+1)=...)
            % requires identical fields and does NOT fill a missing one with [].
            % relPeakHeight -- a display-only dominance measure in [0,1] (1 = the
            % tallest feature at this angle), never a gate: the entry's peak height /
            % the tallest scored candidate at this angle (G13_eff on the Green paths,
            % the late-time Ritz modal sum on 'ritzOnly').
            % Existence does not establish DOMINANCE: a Rayleigh peak can clear the
            % existence threshold for several degrees after a much taller
            % pseudoRayleigh peak has become the observable surface response.
            % relPeakHeight lets a caller tell those apart without recomputing anything.
            % NaN on 'limiting'/'other'.
            % method -- the classification that produced the entry: 'green' | 'refined'
            % | 'ritz' ('' on limiting entries).
            % score -- the value the existence criterion thresholds, on that method's
            % own scale: GreenScore / GreenScoreLayered score ('green'), GreenSpotScore
            % scoreHigh ('refined'), RitzScore dScore (+ polTerm, layered) ('ritz').
            % Scores are log10 amplitude ratios; scoreDb = 20 x score. Compare a score
            % only with its own method's thresholds. NaN on limiting/'other' entries.
            % grade -- 'strong' if the candidate passes the method's 'low' standard,
            % 'weak' if it passes only the active (looser) one; '' on limiting entries.
            % Default 'ritzOnly' entries: vel = the supporting mode's eigenvelocity
            % (RitzScore.localMode); matchScore = its freqDisp (half-space) or its
            % share of the modal sum (layered), both in [0,1]; score = the thresholded
            % value (half-space dScore, layered dScore + polTerm); greenPeakVel = NaN
            % (the modal-sum peak velocity and the rest are in ritzScorePeaks).

        %%% Empirical thresholds, exposed as properties so they can be tuned without
        %%% touching algorithm code.
        velMatchTol      (1,1) double = 0.02   % relative tolerance for "close" velocity match (2%): the Gaussian width of matchScoreOf (matchSpotPeak's ranking, every entry's matchScore) and associateLimitingMode
        matchScoreThresh (1,1) double = 0.05   % minimum matchScore for a Green / spot peak's Ritz match to count as confident
            % For Green-peak matches this only sets `confident`; it never decides
            % existence, which the Green score alone settles. A weak Ritz-basis representation must not
            % veto a Green-confirmed peak: a fundamental mode whose surface
            % displacement is smeared across several near-degenerate Ritz modes has
            % no single strongly-scoring mode, yet is real. A non-confident match
            % is still reported, with its plain argmax mode. In
            % associateLimitingMode it IS an existence floor, since no Green peak
            % backs a limiting-mode association. See
            % docs/reference/architecture/forwardCalcs_modeClassification.md, "Green decides
            % existence, Ritz only assists identity".

        %%% Existence threshold of whichever path runs: a candidate is a surface wave if
        %%% its score clears it. 'high' (sensitive) | 'low' (distinct only) | numeric, read
        %%% on the scale of the score the caller supplies: greenResult.greenScore ->
        %%% GreenScore, greenResult.refined -> GreenSpotScore (half-space);
        %%% greenResult.greenScoreLayered -> GreenScoreLayered (layered); the ritzScore
        %%% input -> RitzScore.acceptedHS / acceptedLayered (numeric: a dScore, or
        %%% layered dScore + polTerm, threshold with the 'high' gates). A Green input
        %%% without a score is an error. The Calc classes pass the method's default.
        threshold             = 'high'
        %%% Re-thresholding support (CalcSingleHS/CalcSingleLayered.reclassify): the
        %%% limiting-velocity association does not depend on the threshold but needs the
        %%% Ritz eigenvectors, which retentionLevel may have pruned -- so a reclassify
        %%% skips it (associateLimiting = false) and carries the previous 'limiting'
        %%% entries over unchanged (carryEntries, appended at the end).
        associateLimiting logical = true
        carryEntries      struct  = struct([])

        %%% Ritz-only existence criterion (no Green input): every peak of the gated Ritz
        %%% modal sum whose score clears the threshold is a surface wave (RitzScore;
        %%% 0..N entries, labelled Rayleigh/pseudoRayleigh, or 'SAW' + aboveSubstrateST
        %%% for layered). Needs the ritzScore input to calc(); without it, an error.
        %%% OUTPUT of the Ritz-only path: every scored modal-sum candidate at this angle,
        %%% accepted or not (field set RitzScore.emptyPeaks): vel (the modal-sum PEAK
        %%% velocity), modeNo / modeVel (the supporting mode; its velocity is the entry's
        %%% vel), score (own-floor modal-sum score), dScore (score + aperture dephasing;
        %%% the half-space value thresholded), prominence, relWidth, relHeight, share,
        %%% freqDisp (half-space: of the supporting mode; NaN layered), slope (its
        %%% dv/dtheta, km/s per deg), windowMass (half-space; NaN layered), belowShear
        %%% (half-space: supporting mode below min shear + slab resolution, i.e. the
        %%% Rayleigh label; layered: below the substrate's min shear, the aboveSubstrateST
        %%% label), vetoed (half-space only: freqDisp < RitzScore.MIN_FREQDISP),
        %%% bottom10 / bottom25 (the supporting mode's |u|^2 fraction in the deepest 10 % /
        %%% 25 % of the slab; layered: energy-weighted over the whole domain, the trapped
        %%% gate), SH / polTerm (layered: surface shear-horizontal share and the
        %%% polarisation term; the value thresholded is dScore + polTerm; NaN half-space),
        %%% accepted (RitzScore.acceptedHS / acceptedLayered). Empty on every other path.
        ritzScorePeaks struct = RitzScore.emptyPeaks()
        %%% OUTPUT of the Green paths: every candidate at this angle, accepted or not.
        %%% Half-space (GreenScore): vel (the |G13_eff| peak velocity, = the entry's
        %%% greenPeakVel), score (the 'high' scale, = the entry's score), lowScore (the
        %%% structure-relative 'low' scale), its parts rawScore + gateTerm + dephaseTerm,
        %%% gamma (pole half-width, km/s), slope (km/s per deg, NaN = neither probe found
        %%% the branch), relWidth, relHeight, accepted. Layered (GreenScoreLayered): the
        %%% GreenScoreLayered.emptyPeaks fields (vel = the entry's greenPeakVel, score = the
        %%% entry's score, its parts logStrength, norm, polTerm, gateTerm, dephaseTerm; SH
        %%% = surface polarisation, 0 sagittal .. 1 shear-horizontal; gamma, slope, ...)
        %%% plus accepted. Empty on the Ritz-only and refined paths.
        greenScorePeaks struct = struct([])
    end

    methods
        function obj = ModeClassification(varargin)
        %%% Default-constructed, or with name-value overrides of the settings
        %%% >> obj = ModeClassification()
        %%% >> obj = ModeClassification('threshold', 'low')
        %%%
        %%% Settable: threshold, velMatchTol, matchScoreThresh, associateLimiting,
        %%% carryEntries. The outputs (mode, method, significantModes, ritzScorePeaks,
        %%% greenScorePeaks) are set by calc() only; any other name is an error
        %%% (ModeClassification:constructor:notSettable). Lets a caller
        %%% (CalcSingleHS/CalcSingleLayered.identifyModes, a script) set a non-default
        %%% threshold without a friend-class copy-and-override; the dynamic obj.(name)
        %%% assignment is permitted despite SetAccess=protected because the constructor
        %%% runs as the class itself.
            settable = {'threshold', 'velMatchTol', 'matchScoreThresh', 'associateLimiting', ...
                'carryEntries'};
            for k = 1:2:numel(varargin)
                if ~(ischar(varargin{k}) || isstring(varargin{k})) || ~ismember(char(varargin{k}), settable)
                    error('ModeClassification:constructor:notSettable', ...
                        'ModeClassification settings are %s.', CalcCommon.orList(settable));
                end
                obj.(char(varargin{k})) = varargin{k+1};
            end
        end

        function obj = calc(obj, ritzSolver, structure, christoffelResult, greenResult, ritzScore)
        %%% Classify significant modes from Ritz (+ optionally Christoffel/Green)
        %%% >> obj = calc(obj, ritzSolver, structure, christoffelResult, greenResult)
        %%% >> obj = calc(obj, ritzSolver, structure, [], [], ritzScore)
        %%%
        %%% INPUT:  1. ritzSolver        .. RitzSolverSAW result, already calc()'d
        %%%         2. structure         .. Structure for this angle. structure.nLay > 1
        %%%                                is the dispatch signal between the
        %%%                                homogeneous (recognizeFull) and layered
        %%%                                (recognizeLayered) paths
        %%%  (opt)  3. christoffelResult .. Green paths: struct with fields longLimVel,
        %%%                                longLimPol, shearLimVels, shearLimPols
        %%%                                (homogeneous HS) or subLongLimVel,
        %%%                                subLongLimPol, subShearLimVels,
        %%%                                subShearLimPols (layered, a substrate-only
        %%%                                reference). [] on the Ritz-only path, which
        %%%                                does not read it
        %%%  (opt)  4. greenResult       .. Green paths: struct with the score input:
        %%%                                half-space greenScore (struct peaks,
        %%%                                matchGreenScorePeaks); layered greenScoreLayered
        %%%                                (struct peaks, matchGreenScoreLayeredPeaks);
        %%%                                half-space refined (struct peaks: GreenSpotScore
        %%%                                candidates after GreenSpotScore.identify,
        %%%                                recognizeRefined). Without a score input, an
        %%%                                error. [] for the Ritz-only path
        %%%  (opt)  5. ritzScore         .. Ritz-only path only: struct with field cands
        %%%                                (the scored candidates: RitzScore.direction /
        %%%                                directionLayered). Without it, or without
        %%%                                cands, an error
        %%% OUTPUT: obj with mode/significantModes filled (and ritzScorePeaks on the
        %%%         Ritz-only path, greenScorePeaks on the Green paths)

            validateattributes(ritzSolver, {'RitzSolverSAW'}, {'scalar'}, ...
                'ModeClassification.calc', 'ritzSolver');
            if nargin >= 3 && ~isempty(structure)
                validateattributes(structure, {'Structure'}, {'scalar'}, ...
                    'ModeClassification.calc', 'structure');
            end
            if nargin < 6; ritzScore = []; end

            if nargin < 4 || isempty(christoffelResult) || nargin < 5 || isempty(greenResult)
                obj.mode = 'ritzOnly';
                if ~ModeClassification.usesRitzScore(ritzScore)
                    error('ModeClassification:ritzOnly:noScore', ...
                        'The Ritz-only path needs the ritzScore input.');
                end
                if ~isfield(ritzScore, 'cands')
                    error('ModeClassification:ritzOnly:noCandidates', ['ritzScore needs the ' ...
                        'field cands (RitzScore.direction / directionLayered output).']);
                end
                isLayered = nargin >= 3 && ~isempty(structure) && structure.nLay > 1;
                obj.method = 'ritz';
                [obj.significantModes, obj.ritzScorePeaks] = ModeClassification.recognizeRitzOnlyScore( ...
                    obj, ritzSolver, ritzScore, isLayered);
            elseif isstruct(greenResult) && isfield(greenResult, 'refined')
                obj.mode   = 'full';
                obj.method = 'refined';
                obj.significantModes = ModeClassification.recognizeRefined(obj, ritzSolver, ...
                    christoffelResult, greenResult.refined);
            elseif nargin >= 3 && ~isempty(structure) && structure.nLay > 1
                obj.mode             = 'layered';
                obj.method           = 'green';
                [obj.significantModes, obj.greenScorePeaks] = ModeClassification.recognizeLayered( ...
                    obj, ritzSolver, christoffelResult, greenResult);
            else
                obj.mode             = 'full';
                obj.method           = 'green';
                [obj.significantModes, obj.greenScorePeaks] = ModeClassification.recognizeFull(obj, ...
                    ritzSolver, christoffelResult, greenResult);
            end
        end
    end

    methods (Static, Hidden)
        function sm = emptySignificantModes()
        %%% The empty significantModes struct array (the property default)
        %%% >> sm = ModeClassification.emptySignificantModes()
            sm = emptyEntries();
        end

        function v = modeVels(sm)
        %%% The eigenvelocity of each entry's Ritz mode (modeVel), as a row
        %%% >> v = ModeClassification.modeVels(sm)
        %%% Entries saved before CLASSIFICATION_VERSION 3 have no modeVel; their vel
        %%% was that eigenvelocity.
            if isfield(sm, 'modeVel'); v = [sm.modeVel]; else; v = [sm.vel]; end
        end
    end

    methods (Static, Access = {?ModeClassification, ?test_ModeClassification, ?GreenSpotScore})
        %%% Access-restricted to this class and its dedicated test class (the MATLAB
        %%% "friend class" idiom): the helpers below are the internal algorithm, not
        %%% public API — only calc() and the properties above are. GreenSpotScore (the
        %%% refiner) uses matchRitzToGreenPeak and matchSpotPeak for its mode identity.
        %%%
        %%% CALLING CONVENTION: these are Static methods, and MATLAB resolves
        %%% obj.staticMethod(args) by the CLASS of obj only — it does NOT prepend obj
        %%% as an implicit first argument the way it does for instance methods. Any
        %%% helper needing obj (to read matchScoreThresh / velMatchTol) therefore
        %%% takes it as an explicit parameter, passed explicitly at every call site.

        function tf = usesGreenScore(greenResult)
        %%% Whether the GreenScore existence criterion applies (half-space 'full' default)
        %%% >> tf = ModeClassification.usesGreenScore(greenResult)
        %%%
        %%% Chosen by the caller supplying greenResult.greenScore (CalcSingleHS with
        %%% sawClassMethod = 'green').
            tf = isstruct(greenResult) && ...
                isfield(greenResult, 'greenScore') && ~isempty(greenResult.greenScore);
        end

        function tf = usesGreenScoreLayered(greenResult)
        %%% Whether the GreenScoreLayered existence criterion applies (layered 'full')
        %%% >> tf = ModeClassification.usesGreenScoreLayered(greenResult)
        %%%
        %%% Chosen by the caller supplying greenResult.greenScoreLayered (CalcSingleLayered).
            tf = isstruct(greenResult) && ...
                isfield(greenResult, 'greenScoreLayered') && ~isempty(greenResult.greenScoreLayered);
        end

        function tf = usesRitzScore(ritzScore)
        %%% Whether the Ritz-only input is present
        %%% >> tf = ModeClassification.usesRitzScore(ritzScore)
            tf = isstruct(ritzScore) && ~isempty(ritzScore);
        end

        function [matched, peaks] = matchGreenScorePeaks(obj, ritz, ps)
        %%% GreenScore candidates above threshold, each matched to its best Ritz mode
        %%% >> [matched, peaks] = ModeClassification.matchGreenScorePeaks(obj, ritz, ps)
        %%%
        %%% INPUT:  1. obj  .. ModeClassification instance (threshold, read as a
        %%%                    GreenScore threshold: 'high' | 'low' | numeric on the 'high'
        %%%                    scale; matchScoreThresh, velMatchTol)
        %%%         2. ritz .. RitzSolverSAW result, already calc()'d
        %%%         3. ps   .. struct: peaks (GreenScore.direction output, lowScore filled
        %%%                    for 'low')
        %%% OUTPUT: matched .. matchGreenPeaks' struct array (existScore = the 'high'-scale
        %%%                    score; strong = passes 'low')
        %%%         peaks   .. every candidate, the greenScorePeaks shape, with accepted
        %%%
        %%% Existence by the score alone; the Ritz match (always succeeding) only supplies
        %%% the mode identity. The candidates are
        %%% already restricted to vel <= vTf (1 + GreenScore.SEARCH_HI_TOL).
            pk = ps.peaks;
            acc = GreenScore.accepted(pk, obj.threshold);
            strong = GreenScore.accepted(pk, 'low');
            peaks = struct('vel', {}, 'score', {}, 'lowScore', {}, 'rawScore', {}, ...
                'gateTerm', {}, 'dephaseTerm', {}, 'gamma', {}, 'slope', {}, 'relWidth', {}, ...
                'relHeight', {}, 'accepted', {});
            for j = 1:numel(pk)
                peaks(j) = struct('vel', pk(j).vel, 'score', pk(j).score, 'lowScore', pk(j).lowScore, ...
                    'rawScore', pk(j).rawScore, 'gateTerm', pk(j).gateTerm, 'dephaseTerm', pk(j).dephaseTerm, ...
                    'gamma', pk(j).gamma, 'slope', pk(j).slope, 'relWidth', pk(j).relWidth, ...
                    'relHeight', pk(j).relHeight, 'accepted', acc(j));
            end
            matched = ModeClassification.matchGreenPeaks(obj, ritz, pk, acc, strong, false);
        end

        function [matched, peaks] = matchGreenScoreLayeredPeaks(obj, ritz, psl)
        %%% GreenScoreLayered candidates above threshold, each matched to its best Ritz mode
        %%% >> [matched, peaks] = ModeClassification.matchGreenScoreLayeredPeaks(obj, ritz, psl)
        %%%
        %%% INPUT:  1. obj  .. ModeClassification instance (threshold, read as a
        %%%                    GreenScoreLayered threshold: 'low' | 'high' | numeric;
        %%%                    matchScoreThresh, velMatchTol)
        %%%         2. ritz .. RitzSolverSAW result, already calc()'d
        %%%         3. psl  .. struct: peaks (GreenScoreLayered.direction output)
        %%% OUTPUT: matched .. as matchGreenScorePeaks (score = the GreenScoreLayered score)
        %%%         peaks   .. every candidate, the greenScorePeaks shape, with accepted
        %%%
        %%% Existence by the score alone, as matchGreenScorePeaks; identity by the layered
        %%% match (matchRitzToGreenPeak with isLayered = true). The candidates
        %%% already span the whole band [0.3 vSmin, 1.05 vLmax] of all layers.
            pk = psl.peaks;
            acc = GreenScoreLayered.accepted(pk, obj.threshold);
            strong = GreenScoreLayered.accepted(pk, 'low');
            peaks = GreenScoreLayered.emptyPeaks();
            peaks(1).accepted = [];
            peaks = peaks([]);                 % the candidate fields plus accepted, no entries
            for j = 1:numel(pk)
                e = pk(j); e.accepted = acc(j);
                peaks(j) = e;
            end
            matched = ModeClassification.matchGreenPeaks(obj, ritz, pk, acc, strong, true);
        end

        function matched = matchGreenPeaks(obj, ritz, pk, acc, strong, isLayered)
        %%% The accepted Green candidates, each matched to its best Ritz mode
        %%% >> matched = ModeClassification.matchGreenPeaks(obj, ritz, pk, acc, strong, isLayered)
        %%%
        %%% INPUT:  1. obj        .. ModeClassification instance (matchScoreThresh, velMatchTol)
        %%%         2. ritz       .. RitzSolverSAW result, already calc()'d
        %%%         3. pk         .. scored candidates (GreenScore / GreenScoreLayered fields
        %%%                          vel, height, prominence, relWidth, relHeight, score)
        %%%         4-5. acc, strong .. logical per candidate: passes the active / the 'low'
        %%%                          standard
        %%%         6. isLayered  .. the identity weight of matchRitzToGreenPeak's score
        %%% OUTPUT: matched .. struct array, one per accepted candidate in candidate order:
        %%%                    vel (the matched Ritz mode's velocity), modeNo, score,
        %%%                    confident (matchRitzToGreenPeak), greenPeakVel, peakHeight,
        %%%                    peakProminence, peakWidth, relPeakHeight, existScore (the
        %%%                    candidate's score), strong
        %%%
        %%% The core of matchGreenScorePeaks / matchGreenScoreLayeredPeaks.
            matched = struct('vel', {}, 'modeNo', {}, 'score', {}, 'confident', {}, ...
                'greenPeakVel', {}, 'peakHeight', {}, 'peakProminence', {}, 'peakWidth', {}, ...
                'relPeakHeight', {}, 'existScore', {}, 'strong', {});
            for j = 1:numel(pk)
                if ~acc(j); continue; end
                [modeNo, score, confident] = ModeClassification.matchRitzToGreenPeak(ritz, pk(j).vel, obj, isLayered);
                matched(end+1) = struct('vel', ritz.vels(modeNo), 'modeNo', modeNo, ...
                    'score', score, 'confident', confident, 'greenPeakVel', pk(j).vel, ...
                    'peakHeight', pk(j).height, 'peakProminence', pk(j).prominence, ...
                    'peakWidth', pk(j).relWidth * pk(j).vel, ...
                    'relPeakHeight', pk(j).relHeight, 'existScore', pk(j).score, 'strong', strong(j)); %#ok<AGROW>
            end
        end

        function [sm, peaks] = recognizeRitzOnlyScore(obj, ritz, rs, isLayered)
        %%% Ritz-only significant-mode list from the gated Ritz modal sum (RitzScore)
        %%% >> [sm, peaks] = ModeClassification.recognizeRitzOnlyScore(obj, ritz, rs, isLayered)
        %%%
        %%% INPUT:  1. obj       .. ModeClassification instance (threshold)
        %%%         2. ritz      .. RitzSolverSAW result (eigenvectors NOT needed -- the
        %%%                         stored candidates already hold what is used)
        %%%         3. rs        .. struct: cands (see calc())
        %%%         4. isLayered .. true for a layered structure
        %%% OUTPUT: sm    .. significantModes struct array, 0..N entries
        %%%         peaks .. every scored candidate (see the ritzScorePeaks property)
        %%%
        %%% The candidates were scored once, threshold-free, when the angle was computed
        %%% (rs.cands: RitzScore.direction on a half-space, RitzScore.directionLayered
        %%% on a layered stack); here they are only thresholded and labelled.
        %%% HALF-SPACE: RitzScore.acceptedHS (veto, dScore, window mass, bottom-trapped
        %%% gate); 'Rayleigh' if the candidate's belowShear is set (RitzScore.direction:
        %%% the supporting mode's eigenvelocity below minShear + RitzScore.slabResolution,
        %%% the Ritz slab's own discretisation error of the slow-shear limit; Ritz
        %%% eigenvalues are upper bounds, so a Rayleigh mode hugging minShear can sit that
        %%% much above it), else 'pseudoRayleigh'. matchScore = freqDisp, score = dScore.
        %%% LAYERED: RitzScore.acceptedLayered (trapped gate, dScore + polarisation
        %%% term); 'SAW' with aboveSubstrateST = the supporting mode above the substrate's
        %%% slowest shear limit (a label only). matchScore = share (the supporting mode's
        %%% fraction of the modal sum at the peak), score = dScore + polTerm.
        %%% No limiting-velocity entries (neither Ritz-only path emits them).
        %%% greenPeakVel stays NaN: there is no Green peak.
            sm = emptyEntries();
            peaks = rs.cands;
            if isLayered
                acc = RitzScore.acceptedLayered(peaks, obj.threshold);
                strong = RitzScore.acceptedLayered(peaks, 'low');
            else
                acc = RitzScore.acceptedHS(peaks, obj.threshold);
                strong = RitzScore.acceptedHS(peaks, 'low');
            end
            for j = 1:numel(peaks)
                peaks(j).accepted = acc(j);
                if ~acc(j); continue; end
                if isLayered
                    sm(end+1) = makeSAWEntry('SAW', ritz, peaks(j).modeNo, peaks(j).share, NaN, ...
                        peaks(j).relHeight, ~peaks(j).belowShear, ...
                        peaks(j).dScore + peaks(j).polTerm, 'ritz', strong(j)); %#ok<AGROW>
                    continue
                end
                %%% belowShear = modeVel < minShear + slab resolution, set by
                %%% RitzScore.direction: the same side decision that chose the candidate's gates
                if peaks(j).belowShear
                    type = 'Rayleigh';
                else
                    type = 'pseudoRayleigh';
                end
                sm(end+1) = makeSAWEntry(type, ritz, peaks(j).modeNo, peaks(j).freqDisp, NaN, ...
                    peaks(j).relHeight, NaN, peaks(j).dScore, 'ritz', strong(j)); %#ok<AGROW>
            end
            %%% sm follows the candidates' peak-velocity order (SawSpectrum.candidates).
        end

        function [sm, polePeaks] = recognizeFull(obj, ritz, christoffelResult, greenResult)
        %%% Green-primary, Ritz-assisted significant-mode list
        %%% >> sm = ModeClassification.recognizeFull(obj, ritz, christoffelResult, greenResult)
        %%% >> [sm, polePeaks] = ModeClassification.recognizeFull(...)
        %%%
        %%% INPUT:  1. obj               .. ModeClassification instance (thresholds)
        %%%         2. ritz              .. RitzSolverSAW result, already calc()'d
        %%%         3. christoffelResult .. struct: longLimVel, longLimPol, shearLimVels,
        %%%                                shearLimPols
        %%%         4. greenResult       .. struct: greenScore (struct peaks); without it, an
        %%%                                error
        %%% OUTPUT: sm        .. significantModes struct array, 0..N entries
        %%%         polePeaks .. every candidate (greenScorePeaks shape)
        %%%
        %%% Green decides existence; Ritz only assists identity. Every Green peak
        %%% that passes the existence test gets a significantModes entry, and no
        %%% later step removes one:
        %%%   1+2. Existence (greenResult.greenScore): every pole-fitted G13_eff
        %%%      candidate whose GreenScore clears the threshold (matchGreenScorePeaks;
        %%%      candidates up to vTf(1 + SEARCH_HI_TOL)). Each accepted peak is matched
        %%%      to its best-supporting Ritz mode (always succeeds); a low-confidence
        %%%      match is still reported -- the match never decides existence.
        %%%   3. Only then label each entry Rayleigh or pseudoRayleigh, from its
        %%%      Green peak velocity against minShear x (1 + RAYLEIGH_LABEL_TOL) -- a
        %%%      name only, with no effect on which candidates survive.
        %%%
        %%% A single search mask, with no velocity-based restriction on which Ritz
        %%% modes may match a peak, stops a peak being rejected from both sides of
        %%% a Ritz/Christoffel discretisation mismatch at minShear. See
        %%% docs/reference/architecture/forwardCalcs_modeClassification.md, "Green decides
        %%% existence, Ritz only assists identity" and "One search range, one
        %%% deciding quantity".

            %%% Plain min(), regardless of which shear branch it came from, IS the
            %%% correct Rayleigh/pseudoRayleigh boundary here -- do NOT "fix" this into
            %%% something shearLimIwt-aware. Limiting velocities are group velocities
            %%% from a non-convex sagittal-plane shear scan, so >2 roots is normal, and
            %%% the slowest is the true radiation threshold whatever its branch; branch
            %%% identity only matters for the inverse module's mode tracking. See
            %%% docs/reference/architecture/forwardCalcs_modeClassification.md,
            %%% "minShear = min(shearLimVels)".
            minShear  = min(christoffelResult.shearLimVels);
            sm = emptyEntries();

            if ~ModeClassification.usesGreenScore(greenResult)
                error('ModeClassification:recognizeFull:noScore', ['greenResult carries no ' ...
                    'greenScore input: the Green path needs the scored candidates.']);
            end
            %%% 1+2. Existence: every pole-fitted G13_eff candidate whose GreenScore clears
            %%% the threshold, each matched to its best-supporting Ritz mode.
            [matched, polePeaks] = ModeClassification.matchGreenScorePeaks(obj, ritz, greenResult.greenScore);

            %%% 3. Label Rayleigh/pseudoRayleigh LAST, purely from each candidate's
            %%% Green peak velocity against minShear -- a name only, decided after
            %%% existence and mode identity are settled, so it influences neither.
            %%% The boundary carries Green's own resolution, RAYLEIGH_LABEL_TOL: a
            %%% peak within it above minShear is a Rayleigh wave hugging minShear.
            vLabel = minShear * (1 + ModeClassification.RAYLEIGH_LABEL_TOL);
            for k = 1:numel(matched)
                if matched(k).greenPeakVel < vLabel
                    type = 'Rayleigh';
                else
                    type = 'pseudoRayleigh';
                end
                sm(end+1) = makeSAWEntry(type, ritz, matched(k).modeNo, matched(k).score, ...
                    matched(k).greenPeakVel, matched(k).relPeakHeight, NaN, ...
                    matched(k).existScore, 'green', matched(k).strong); %#ok<AGROW>
            end

            %%% 4. Limiting-velocity entries (associated, or carried over on a reclassify)
            sm = ModeClassification.appendLimiting(obj, ritz, christoffelResult, sm);

            %%% 5. everything else: 'other' — not pre-populated, materialise on demand
        end

        function sm = appendLimiting(obj, ritz, lim, sm)
        %%% Append the limiting-velocity entries (or the carried ones)
        %%% >> sm = ModeClassification.appendLimiting(obj, ritz, lim, sm)
        %%%
        %%% lim .. struct longLimVel, longLimPol, shearLimVels, shearLimPols: the half-space's
        %%%        own limits, or the substrate's on a layered stack (recognizeLayered passes
        %%%        its sub* fields under these names).
        %%%
        %%% With associateLimiting = false (a reclassify) the previous limiting entries
        %%% (carryEntries) are appended unchanged: the association needs eigenvectors
        %%% that retentionLevel may have pruned, and does not depend on the threshold.
            if ~obj.associateLimiting
                sm = [sm, obj.carryEntries];
                return
            end
            %%% The displacement profiles of the first modes are evaluated ONCE here, up to the
            %%% highest mode any limiting velocity's +-5 window reaches, and shared by every
            %%% association (calcDisplacement returns modes 1..N independently per mode, so
            %%% the window's columns are the same as from a call per limiting velocity).
            vLims   = [lim.shearLimVels(:); lim.longLimVel];
            lastMode = 0;
            for v = vLims.'
                [~, closestIdx] = min(abs(ritz.vels - v));
                lastMode = max(lastMode, min(numel(ritz.vels), closestIdx + 5));
            end
            displDir = ritz.calcDisplacement(lastMode);

            %%% 4. Limiting velocity modes (each shear limit, then 4b. the longitudinal
            %%% one) — depthInsensitivity x polarisation projection, NOT oopSurfDispl
            %%% (different physical criterion). The longitudinal branch always yields
            %%% exactly one root; it associates the grazing-incidence mode at the
            %%% Christoffel-computed threshold only, exactly as the shear case does — NOT a
            %%% general search for longitudinal bulk modes across the spectrum.
            %%%
            %%% A degenerate group (degenerateLimitGroups: one velocity, orthogonal
            %%% polarisations -- the two shear roots at an acoustic axis) is associated
            %%% JOINTLY by associateLimitingGroup: its polarisations are only a basis of the
            %%% eigenspace, so pairing each with a mode by its own vector would be decided
            %%% by that basis. Every other limit keeps associateLimitingMode.
            pols   = [lim.shearLimPols, lim.longLimPol(:)];
            grpId  = ModeClassification.degenerateLimitGroups(vLims, pols);
            modeOf = cell(numel(vLims), 1);
            for g = unique(grpId).'
                members = find(grpId == g);
                if isscalar(members)
                    [modeOf{members}, ~] = ModeClassification.associateLimitingMode(ritz, ...
                        vLims(members), pols(:, members), obj, displDir);
                else
                    modeNos = ModeClassification.associateLimitingGroup(ritz, ...
                        vLims(members), pols(:, members), obj, displDir);
                    for k = 1:numel(modeNos)
                        modeOf{members(k)} = modeNos(k);
                    end
                end
            end
            for i = 1:numel(vLims)
                if ~isempty(modeOf{i})
                    sm(end+1) = makeLimitingEntry(ritz, modeOf{i}, vLims(i)); %#ok<AGROW>
                end
            end
        end

        function grpId = degenerateLimitGroups(vLims, pols)
        %%% Group the limiting velocities of one angle into degenerate sets
        %%% >> grpId = ModeClassification.degenerateLimitGroups(vLims, pols)
        %%%
        %%% INPUT:  1. vLims .. [n x 1] limiting velocities (shear limits, then longLimVel)
        %%%         2. pols  .. [3 x n] their polarisations
        %%% OUTPUT: grpId .. [n x 1] group index; a group of one is a non-degenerate limit
        %%%
        %%% Same group: velocities within LimitingVels.DEGENERATE_REL_TOL (relative) AND
        %%% polarisations orthogonal to LIMIT_ORTH_TOL. Velocity alone is not enough: a
        %%% mirror-symmetric sagittal plane gives every root a twin at -theta* with the same
        %%% velocity (Cu(001) at 15-45 deg) and a mirrored, NOT orthogonal, polarisation --
        %%% one physical root, associated on its own. LimitingVels.findRoots hands a
        %%% degenerate group an orthonormal canonical basis of its eigenspace, so the test
        %%% separates the two (a mirror twin is orthogonal only if its out-of-plane
        %%% component is exactly 1/sqrt(2)).
            n = numel(vLims);
            grpId = zeros(n, 1);
            nrm = sqrt(sum(pols.^2, 1));
            g = 0;
            for a = 1:n
                if grpId(a) > 0; continue; end
                g = g + 1;
                grpId(a) = g;
                members = a;
                for b = a+1:n
                    if grpId(b) == 0 ...
                            && abs(vLims(b) - vLims(a)) <= LimitingVels.DEGENERATE_REL_TOL * abs(vLims(a)) ...
                            && all(abs(pols(:, members).' * pols(:, b)) ...
                                   <= ModeClassification.LIMIT_ORTH_TOL * nrm(members).' * nrm(b))
                        grpId(b) = g;
                        members(end+1) = b; %#ok<AGROW>
                    end
                end
            end
        end

        function [modeNos, scores] = associateLimitingGroup(ritz, vGroup, polGroup, obj, displDir)
        %%% Associate a degenerate group of limiting velocities with distinct Ritz modes
        %%% >> [modeNos, scores] = ModeClassification.associateLimitingGroup(ritz, vGroup, polGroup, obj, displDir)
        %%%
        %%% INPUT:  1. ritz     .. RitzSolverSAW result, already calc()'d
        %%%         2. vGroup   .. [g x 1] the group's velocities (equal to
        %%%                        LimitingVels.DEGENERATE_REL_TOL); the first sets the
        %%%                        window and the velocity proximity
        %%%         3. polGroup .. [3 x g] polarisations spanning the degenerate eigenspace
        %%%         4. obj      .. ModeClassification instance
        %%%  (opt)  5. displDir .. as associateLimitingMode
        %%% OUTPUT: modeNos .. [1 x m] distinct modes, m <= g, ASCENDING; member k of the
        %%%                    group takes modeNos(k), the members past m get none
        %%%         scores  .. [1 x m] their scores
        %%%
        %%% associateLimitingMode's score with the projection onto the group's SPAN,
        %%% ||Q' displNorm|| for an orthonormal basis Q of polGroup, in place of one
        %%% vector's |pol . displNorm| (the two agree for g = 1). The span projection
        %%% does not depend on which basis of the eigenspace the root finder returned. The
        %%% g best-scoring modes above matchScoreThresh are taken (ties to the lower mode)
        %%% and dealt out in mode order: g modes for g roots, never one mode twice.
            sensN = (ritz.depthInsensitivity - min(ritz.depthInsensitivity)) ...
                  ./ max(max(ritz.depthInsensitivity) - min(ritz.depthInsensitivity), eps);
            vLim = vGroup(1);
            [~, closestIdx] = min(abs(ritz.vels - vLim));
            window = max(1, closestIdx-5) : min(numel(ritz.vels), closestIdx+5);
            if nargin < 5 || isempty(displDir)
                displDir = ritz.calcDisplacement(max(window));
            end
            displSum  = squeeze(sum(abs(displDir(:, window, :)), 1));
            displNorm = displSum ./ max(sum(displSum, 2), eps);
            [Q, ~] = qr(polGroup, 0);
            projection = sqrt(sum((displNorm * Q).^2, 2));
            relDist      = abs(ritz.vels(window) - vLim) / vLim;
            velProximity = exp(-(relDist / obj.velMatchTol).^2);
            s = sensN(window) .* projection(:) .* velProximity(:);
            [sSorted, ord] = sort(s, 'descend');   % stable: ties keep the lower mode first
            ord = ord(sSorted > obj.matchScoreThresh);
            ord = ord(1:min(numel(vGroup), numel(ord)));
            [modeNos, iSort] = sort(window(ord));
            scores = s(ord(iSort)).';
        end

        function sm = recognizeRefined(obj, ritz, christoffelResult, rf)
        %%% Significant-mode list from the refined (spot-averaged) Green score, half-space
        %%% >> sm = ModeClassification.recognizeRefined(obj, ritz, christoffelResult, rf)
        %%%
        %%% INPUT:  1. obj               .. ModeClassification instance (threshold, read on
        %%%                                GreenSpotScore's scale)
        %%%         2. ritz              .. RitzSolverSAW result, already calc()'d
        %%%         3. christoffelResult .. as recognizeFull (limiting entries)
        %%%         4. rf                .. struct: peaks -- every GreenSpotScore candidate after
        %%%                                GreenSpotScore.identify (modeNo, matchScore, type from
        %%%                                the spot-weighted minShear)
        %%% OUTPUT: sm .. significantModes struct array: one entry per accepted candidate
        %%%               (greenPeakVel = the spot peak velocity, score = scoreHigh), then the
        %%%               limiting entries, as recognizeFull
        %%%
        %%% Existence by GreenSpotScore.accepted; identity and label were fixed by
        %%% GreenSpotScore.identify (the single-direction pole, spot-weighted minShear), so the
        %%% entries are plottable and can be drilled into like the green ones.
            sm = emptyEntries();
            pk = rf.peaks;
            if ~isempty(pk)
                acc = GreenSpotScore.accepted(pk, obj.threshold);
                strong = GreenSpotScore.accepted(pk, 'low');
                hMax = max([pk.height]);
                for j = 1:numel(pk)
                    if ~acc(j); continue; end
                    sm(end+1) = makeSAWEntry(pk(j).type, ritz, pk(j).modeNo, pk(j).matchScore, ...
                        pk(j).vel, pk(j).height / hMax, NaN, pk(j).scoreHigh, 'refined', ...
                        strong(j)); %#ok<AGROW>
                end
            end
            sm = ModeClassification.appendLimiting(obj, ritz, christoffelResult, sm);
        end

        function [sm, polePeaks] = recognizeLayered(obj, ritz, subChristoffelResult, greenResult)
        %%% Green-primary, Ritz-assisted significant-mode list for a layered structure
        %%% >> sm = ModeClassification.recognizeLayered(obj, ritz, subChristoffelResult, greenResult)
        %%% >> [sm, polePeaks] = ModeClassification.recognizeLayered(...)
        %%%
        %%% INPUT:  1. obj                  .. ModeClassification instance (thresholds)
        %%%         2. ritz                 .. RitzSolverSAW result, already calc()'d
        %%%         3. subChristoffelResult .. struct: subLongLimVel, subLongLimPol,
        %%%                                    subShearLimVels, subShearLimPols — a
        %%%                                    SUBSTRATE-ONLY reference from
        %%%                                    CalcSingleLayered.calcLim(), NOT a
        %%%                                    system-wide limiting velocity
        %%%         4. greenResult          .. struct: greenScoreLayered (struct peaks); without
        %%%                                    it, an error
        %%% OUTPUT: sm        .. significantModes struct array, 0..N entries; type in
        %%%                      {'SAW','limiting','other'} — no 'Love' type (a Love-type
        %%%                      wave is a 'SAW' entry when the threshold admits it)
        %%%         polePeaks .. every candidate (greenScorePeaks shape)
        %%%
        %%% No below/above-shear-limit split, unlike recognizeFull's Rayleigh/
        %%% pseudoRayleigh labelling: a layered stack has no single system-wide
        %%% shear limit to split on. Every accepted Green peak across the whole
        %%% velocity axis is a SAW-like candidate, and the substrate shear limit
        %%% only informs the aboveSubstrateST flag afterwards, never which peaks
        %%% are searched for.
        %%%
        %%% There is deliberately NO Rayleigh/Sezawa split by velocity rank. Real
        %%% angular sweeps show that convention breaking down under anisotropy: the
        %%% lowest-velocity entry can match a Green peak at 15-30% of the spectrum's
        %%% dominant peak height while the actually-observable dominant peak sits
        %%% higher. Re-ranking by dominance instead is no better — "lowest velocity"
        %%% and "tallest peak" decouple off a principal plane in ways neither simple
        %%% rule captures; the convention only holds cleanly in the isotropic case.
        %%% Every matched candidate is a single 'SAW' type, with vel, matchScore,
        %%% relPeakHeight and aboveSubstrateST all reported per entry so a caller can
        %%% judge branch order, dominance and leakiness from real numbers.
        %%%
        %%% Existence as in recognizeFull, by one score: GreenScoreLayered
        %%% (greenResult.greenScoreLayered; candidates over the band [0.3 vSmin, 1.05 vLmax]
        %%% of ALL layers, no substrate; 'low' = sagittal-type waves, 'high' = also Love-
        %%% type; tuned on the measured NiTi/MgO map, Grabec et al. 2024).
        %%% Mode identity: matchRitzToGreenPeak (the Ritz mode nearest the Green peak);
        %%% its reported matchScore uses the layered weight (oopSurfDispl, without the
        %%% half-space freqDisp factor, since every true guided mode of a layered stack
        %%% is dispersive).

            sm = emptyEntries();

            if ~ModeClassification.usesGreenScoreLayered(greenResult)
                error('ModeClassification:recognizeLayered:noScore', ['greenResult carries ' ...
                    'no greenScoreLayered input: the layered Green path needs the scored candidates.']);
            end
            %%% 1. Existence: every pole-fitted G13_eff candidate whose GreenScoreLayered
            %%% score clears the threshold, each matched to its best-supporting Ritz
            %%% mode (layered identity score).
            [matched, polePeaks] = ModeClassification.matchGreenScoreLayeredPeaks(obj, ritz, ...
                greenResult.greenScoreLayered);

            %%% 2. Every matched candidate is a single 'SAW' type. Sorted by velocity
            %%% for readability only, NOT to derive a Rayleigh/Sezawa label;
            %%% aboveSubstrateST is informational, not a type branch.
            [~, order] = sort([matched.greenPeakVel]);   % = the entries' vel
            matched     = matched(order);
            minSubShear = min(subChristoffelResult.subShearLimVels);
            for k = 1:numel(matched)
                leaky = matched(k).greenPeakVel > minSubShear;
                sm(end+1) = makeSAWEntry('SAW', ritz, matched(k).modeNo, matched(k).score, ...
                    matched(k).greenPeakVel, matched(k).relPeakHeight, ...
                    leaky, matched(k).existScore, 'green', matched(k).strong); %#ok<AGROW>
            end

            %%% 3. Limiting velocity modes (or the carried ones) -- the same association
            %%% as the half-space, fed substrate data instead of system-wide data: it has
            %%% no half-space assumption, taking a bare velocity + polarisation pair.
            lim = struct('longLimVel', subChristoffelResult.subLongLimVel, ...
                'longLimPol', subChristoffelResult.subLongLimPol, ...
                'shearLimVels', subChristoffelResult.subShearLimVels, ...
                'shearLimPols', subChristoffelResult.subShearLimPols);
            sm = ModeClassification.appendLimiting(obj, ritz, lim, sm);

            %%% 4. everything else: 'other' — not pre-populated, same
            %%% convention as recognizeFull
        end

        function [modeNo, score, confident] = matchRitzToGreenPeak(ritz, vGreen, obj, isLayered)
        %%% Match a single-direction Green peak velocity to the nearest Ritz mode
        %%% >> [modeNo, score, confident] = ModeClassification.matchRitzToGreenPeak(ritz, vGreen, obj)
        %%% >> [modeNo, score, confident] = ModeClassification.matchRitzToGreenPeak(ritz, vGreen, obj, isLayered)
        %%%
        %%% INPUT:  1. ritz      .. RitzSolverSAW result, already calc()'d
        %%%         2. vGreen    .. the 'green' peak velocity (GreenScore /
        %%%                         GreenScoreLayered: the pole-fitted G13_eff peak of
        %%%                         this one direction)
        %%%         3. obj       .. ModeClassification instance (velMatchTol,
        %%%                         matchScoreThresh)
        %%%  (opt)  4. isLayered .. true for a layered stack (recognizeLayered's calls);
        %%%                         default false (half-space). Passed by the caller
        %%%                         rather than read from ritz, so plain-struct ritz
        %%%                         stand-ins work
        %%% OUTPUT: modeNo    .. the Ritz mode nearest vGreen; among modes within
        %%%                      GREEN_TIE_REL x vGreen of the nearest distance, the
        %%%                      one with the largest identityWeight -- always a real
        %%%                      mode number, never []
        %%%         score     .. matchScoreOf that mode
        %%%                      (Gaussian velProximity x identityWeight): how well it
        %%%                      matches, not how it was chosen
        %%%         confident .. score > matchScoreThresh. false never rejects the
        %%%                      peak
        %%%
        %%% Why nearest, with the weight only as a tie-break. The 'green' peak is the pole of the same
        %%% direction the Ritz basis was solved for, placed to a few m/s at most, far finer than
        %%% the slab-mode spacing around a leaky (above-shear) wave (1-2 % of v on the slabs
        %%% computed; it shrinks as the slab deepens). A leaky wave is carried by several
        %%% neighbouring slab modes, so a rule that
        %%% weights proximity loosely against the surface weight can bind a neighbour tens of
        %%% m/s away. A trapped wave has a Ritz mode at its peak, so nearest is
        %%% exact there. Within the peak's own placement error (GREEN_TIE_REL) the
        %%% nearer of two modes is chance, so the surface weight decides. The mode that
        %%% dominates the Ritz modal sum at the peak is not used: a weak (Love-type)
        %%% layered peak loses it to a strong neighbour far away in velocity.
        %%% A 'refined' spot peak is not a single-direction pole (the spot average moves
        %%% it by up to tens of m/s): GreenSpotScore.identify first maps it to the nearest
        %%% single-direction pole and binds that pole here.
            if nargin < 4; isLayered = false; end
            d    = abs(ritz.vels(:) - vGreen);
            near = find(d <= min(d) + ModeClassification.GREEN_TIE_REL * vGreen);
            [~, k] = max(ModeClassification.identityWeight(ritz, near, isLayered));
            modeNo = near(k);
            score     = ModeClassification.matchScoreOf(ritz, vGreen, obj, isLayered, modeNo);
            confident = score > obj.matchScoreThresh;
        end

        function [modeNo, score, confident] = matchSpotPeak(ritz, vSpot, obj)
        %%% Match a 'refined' (spot-averaged) half-space peak with no single-direction pole to a Ritz mode
        %%% >> [modeNo, score, confident] = ModeClassification.matchSpotPeak(ritz, vSpot, obj)
        %%%
        %%% INPUT:  1. ritz  .. RitzSolverSAW result, already calc()'d
        %%%         2. vSpot .. GreenSpotScore peak velocity (spot-averaged)
        %%%         3. obj   .. ModeClassification instance (velMatchTol, matchScoreThresh)
        %%% OUTPUT: modeNo    .. argmax mode of the combined score, always a real mode
        %%%         score     .. its combined score
        %%%         confident .. score > matchScoreThresh (never rejects the peak)
        %%%
        %%% Score = Gaussian velProximity x identityWeight: oopSurfDispl x freqDisp
        %%% (half-space; read from ritz.diagnostics, not recomputed).
        %%% The fallback of GreenSpotScore.identify, which normally maps a spot peak to
        %%% its single-direction pole and binds that with matchRitzToGreenPeak; used when
        %%% no pole lies within its window (rare: 3 of 1204 computed 'refined' entries,
        %%% the kernel comb on a synthetic cubic of Zener ratio 20). Why not nearest (as
        %%% matchRitzToGreenPeak): the spot average moves a peak off the single-direction
        %%% pole on a steep branch (synthetic cubic, Zener ratio 20: >= 40 m/s), into the
        %%% dense bulk-like slab modes beside a Rayleigh wave hugging minShear. There the nearest mode is a bulk mode
        %%% (half-space 'refined' trapped entries: 94.5 % nearest vs 100 % this rule),
        %%% and the surface weight is what finds the wave's mode.
        %%% freqDisp separates a true (non-dispersive) SAW from bulk-like Ritz modes;
        %%% it (rather than depthInsensitivity) gave the better argmax accuracy on
        %%% half-space checks, and covers depthInsensitivity's role (an artefact of
        %%% the finite domain depth would show spurious dispersion) -- an assumption
        %%% not verified as rigorously as the accuracy comparison. The score only
        %%% ranks modes for a peak the spot score has already accepted; it is never
        %%% an existence floor.
        %%%
        %%% velProximity must be GAUSSIAN in relative velocity distance, scaled by
        %%% obj.velMatchTol. A linear penalty (1 - relDist) is far too weak: a mode
        %%% with an outsized oopSurfDispl x freqDisp product then still wins the
        %%% argmax against a peak ~10% away in velocity, producing spurious
        %%% duplicate matches — the same dominant Rayleigh mode claimed again for an
        %%% unrelated pseudoRayleigh-range peak — which do not go away under mesh
        %%% refinement.
        %%%
            score = ModeClassification.matchScoreOf(ritz, vSpot, obj, false, 1:numel(ritz.vels));
            [bestScore, modeNo] = max(score);
            score     = bestScore;
            confident = bestScore > obj.matchScoreThresh;
        end

        function s = matchScoreOf(ritz, v, obj, isLayered, modeNos)
        %%% Combined match score of the given modes against a peak velocity
        %%% >> s = ModeClassification.matchScoreOf(ritz, v, obj, isLayered, modeNos)
        %%%
        %%% OUTPUT: s .. [numel(modeNos) x 1]: exp(-(|v_n - v| / v / velMatchTol)^2) x
        %%%              identityWeight -- matchSpotPeak's ranking score, and the
        %%%              matchScore every 'green' and 'refined' entry reports
            vn = ritz.vels(modeNos);
            relDist      = abs(vn(:) - v) / v;
            velProximity = exp(-(relDist / obj.velMatchTol).^2);
            s = velProximity .* ModeClassification.identityWeight(ritz, modeNos, isLayered);
        end

        function w = identityWeight(ritz, modeNos, isLayered)
        %%% Per-mode weight of the Green-peak identity score (velocity proximity aside)
        %%% >> w = ModeClassification.identityWeight(ritz, modeNos, isLayered)
        %%%
        %%% OUTPUT: w .. [numel(modeNos) x 1]: oopSurfDispl x freqDisp (half-space) or
        %%%              oopSurfDispl alone (layered: every true guided mode of a stack
        %%%              is dispersive, because the film sets a length scale, so
        %%%              freqDisp carries no identity information there)
        %%% Never gated on oopSurfDispl's absolute magnitude: a layered structure's
        %%% oopSurfDispl scale is not known to be comparable to a half-space's
        %%% (energy splits across film + substrate).
            w = ritz.diagnostics.oopSurfDispl(modeNos);
            w = w(:);
            if ~isLayered
                fd = ritz.diagnostics.freqDisp(modeNos);
                w  = w .* fd(:);
            end
        end

        function [modeNo, score] = associateLimitingMode(ritz, vLim, polLim, obj, displDir)
        %%% Associate a Christoffel limiting velocity (shear OR longitudinal)
        %%% with a Ritz mode
        %%% >> [modeNo, score] = ModeClassification.associateLimitingMode(ritz, vLim, polLim, obj)
        %%% >> [modeNo, score] = ModeClassification.associateLimitingMode(ritz, vLim, polLim, obj, displDir)
        %%%
        %%% INPUT:  1. ritz   .. RitzSolverSAW result, already calc()'d
        %%%         2. vLim   .. one Christoffel limiting velocity — a shear branch
        %%%                      (shearLimVels(i)) or the longitudinal branch
        %%%                      (longLimVel); the algorithm treats them alike
        %%%         3. polLim .. [3x1] polarisation at that limiting velocity
        %%%         4. obj    .. ModeClassification instance (velMatchTol,
        %%%                      matchScoreThresh)
        %%%  (opt)  5. displDir .. ritz.calcDisplacement(N) with N >= the last mode of this
        %%%                      velocity's window (appendLimiting evaluates it once for all
        %%%                      limiting velocities of the angle); computed here if omitted
        %%% OUTPUT: modeNo .. best-matching Ritz mode index within a +-5 window
        %%%                   of vLim, or [] if no mode clears matchScoreThresh
        %%%         score  .. combined score (NaN if modeNo is [])
        %%%
        %%% IMPORTANT: uses ritz.depthInsensitivity (the raw base-RitzSolver property,
        %%% log10 form, all nCompModes), NOT ritz.diagnostics.depthInsensitivity — the
        %%% latter is already offset by its own min over the whole mode range, and
        %%% re-offsetting it within a small window here would
        %%% double-process it.
        %%%
        %%% The Gaussian velProximity term and the matchScoreThresh floor are both
        %%% load-bearing: a bare sensN.*projection score with only a >eps floor picks
        %%% the dominant Rayleigh mode itself as the "closest available" match when no
        %%% Ritz mode is genuinely near vLim, hiding the real Rayleigh entry behind a
        %%% spurious 'limiting' one.

            sensN = (ritz.depthInsensitivity - min(ritz.depthInsensitivity)) ...
                  ./ max(max(ritz.depthInsensitivity) - min(ritz.depthInsensitivity), eps);   % [0,1], all modes

            [~, closestIdx] = min(abs(ritz.vels - vLim));
            window = max(1, closestIdx-5) : min(numel(ritz.vels), closestIdx+5);

            if nargin < 5 || isempty(displDir)
                displDir = ritz.calcDisplacement(max(window));        % [depthPts x nModes x 3]
            end
            displSum  = squeeze(sum(abs(displDir(:, window, :)), 1));   % [numel(window) x 3]
            displNorm = displSum ./ max(sum(displSum, 2), eps);          % normalise per mode

            projection = abs(displNorm * polLim(:));   % polarisation projection, [numel(window) x 1]
            %%% Every operand below must stay a COLUMN. A row .* column would
            %%% silently broadcast into a [numel(window) x numel(window)] matrix
            %%% instead of an elementwise product.
            relDist      = abs(ritz.vels(window) - vLim) / vLim;
            velProximity = exp(-(relDist / obj.velMatchTol).^2);
            scores = sensN(window) .* projection(:) .* velProximity(:);
            [bestScore, k] = max(scores);

            if bestScore <= obj.matchScoreThresh
                modeNo = [];
                score  = NaN;
            else
                modeNo = window(k);
                score  = bestScore;
            end
        end
    end
end

function sm = emptyEntries()
%%% Build the canonical EMPTY significantModes struct array
%%% >> sm = emptyEntries()
%%%
%%% One place that knows the field set. Every recognize* path seeds its output
%%% with this, and the significantModes property default reads it (through
%%% emptySignificantModes), because MATLAB struct-array concatenation (sm(end+1) = ...,
%%% [sm1, sm2]) requires IDENTICAL field sets across all three paths and does
%%% NOT fill a missing field with []. A field added here must also be added to
%%% makeSAWEntry/makeLimitingEntry.
    sm = struct('type', {}, 'vel', {}, 'modeNo', {}, 'modeVel', {}, 'matchScore', {}, 'greenPeakVel', {}, ...
                'oopDisplScore', {}, 'dispersionScore', {}, 'limVelAssoc', {}, ...
                'aboveSubstrateST', {}, 'relPeakHeight', {}, ...
                'method', {}, 'score', {}, 'scoreDb', {}, 'grade', {});
end

function s = makeSAWEntry(typeStr, ritz, modeNo, matchScore, greenVel, relPeakHeight, aboveSubstrateST, score, method, strong)
%%% Build one surface-wave significantModes entry (any recognize* path)
%%% >> s = makeSAWEntry(typeStr, ritz, modeNo, matchScore, greenVel, relPeakHeight, aboveSubstrateST, score, method, strong)
%%%
%%% A local function after the classdef end, visible to every method above.
%%%
%%% INPUT:  1. typeStr          .. 'Rayleigh' | 'pseudoRayleigh' (half-space) | 'SAW'
%%%                                (layered)
%%%         2. ritz             .. RitzSolverSAW result, already calc()'d
%%%         3. modeNo           .. matched Ritz mode index
%%%         4. matchScore       .. matchScoreOf the matched mode (Green and refined
%%%                                paths) or the supporting mode's freqDisp / modal-sum
%%%                                share (Ritz-only path)
%%%         5. greenVel         .. the Green (or spot) peak velocity that was matched; NaN
%%%                                on the Ritz-only path
%%%         6. relPeakHeight    .. peak height / tallest candidate at that angle, in [0,1]
%%%         7. aboveSubstrateST .. layered only: velocity above the substrate shear limit
%%%                                (informational); NaN on a half-space
%%%         8. score            .. the existence score on the method's own scale
%%%         9. method           .. 'green' | 'refined' | 'ritz'
%%%        10. strong           .. the candidate passes the method's 'low' standard
%%% OUTPUT: s .. one significantModes-shaped struct entry
%%%
%%% vel is the velocity of the method that decided existence: greenVel on the Green paths
%%% (the 'green' pole, the 'refined' spot peak), the mode's eigenvelocity on the Ritz-only
%%% path (greenVel NaN). modeVel is always the bound mode's eigenvelocity. Why: a leaky
%%% wave has no slab mode at its resonance, and the bound mode sits up to half a mode
%%% spacing off it (p95 15 m/s half-space, 26 m/s CuFilm/Si) where the Green peak is within
%%% ~5 m/s; one source per method keeps a branch continuous across the shear limit.
    if strong; grade = 'strong'; else; grade = 'weak'; end
    vel = greenVel;
    if ~isfinite(vel); vel = ritz.vels(modeNo); end
    s = struct('type', typeStr, 'vel', vel, 'modeNo', modeNo, 'modeVel', ritz.vels(modeNo), ...
        'matchScore', matchScore, 'greenPeakVel', greenVel, ...
        'oopDisplScore', ritz.diagnostics.oopSurfDispl(modeNo), ...
        'dispersionScore', ritz.diagnostics.depthInsensitivity(modeNo), 'limVelAssoc', NaN, ...
        'aboveSubstrateST', aboveSubstrateST, 'relPeakHeight', relPeakHeight, ...
        'method', method, 'score', score, ...
        'scoreDb', 20 * score, 'grade', grade);
end

function s = makeLimitingEntry(ritz, modeNo, vLim)
%%% Build one 'limiting' significantModes entry
%%% >> s = makeLimitingEntry(ritz, modeNo, vLim)
%%%
%%% INPUT:  1. ritz   .. RitzSolverSAW result, already calc()'d
%%%         2. modeNo .. Ritz mode associateLimitingMode picked
%%%         3. vLim   .. the Christoffel limiting velocity it was associated with —
%%%                      a shear OR longitudinal branch, system-wide (recognizeFull)
%%%                      or substrate-only (recognizeLayered); the entry makes no
%%%                      distinction, as associateLimitingMode makes none
%%% OUTPUT: s .. one significantModes-shaped struct entry
%%%
%%% NOT a peak-matched entry, so matchScore/greenPeakVel/relPeakHeight/score are NaN
%%% and method/grade empty. oopDisplScore is NaN too: the limiting association scores
%%% on depthInsensitivity x polarisation projection, a different physical criterion. aboveSubstrateST is NaN, being meaningful only on layered
%%% SAW entries.
    s = struct('type', 'limiting', 'vel', ritz.vels(modeNo), 'modeNo', modeNo, 'modeVel', ritz.vels(modeNo), ...
        'matchScore', NaN, 'greenPeakVel', NaN, 'oopDisplScore', NaN, ...
        'dispersionScore', ritz.diagnostics.depthInsensitivity(modeNo), ...
        'limVelAssoc', vLim, 'aboveSubstrateST', NaN, 'relPeakHeight', NaN, ...
        'method', '', 'score', NaN, 'scoreDb', NaN, 'grade', '');
end
