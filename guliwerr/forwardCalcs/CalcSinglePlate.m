classdef CalcSinglePlate < handle & matlab.mixin.Copyable
%%% CalcSinglePlate — single-angle container for the free-plate (Lamb-wave)
%%% forward calculation (RitzSolverLamb only).
%%%
%%% Given one Structure (free plate: rigid-bottom BC does not apply, both
%%% surfaces are free) and one propagation angle, runs the Lamb solver for that
%%% direction and holds the result. Usable standalone, and the per-angle building
%%% block CalcAngSweepPlate drives across an angle array.
%%%
%%% LEANER than CalcSingleHS/CalcSingleLayered -- no calcLim(), no calcGreen(),
%%% no identifyModes()/classification property. Both reasons are physical:
%%%   1. There is no Green's-function analogue for a free plate in this codebase
%%%      (GreenSolver models a half-space/layered-substrate surface response to a
%%%      point force). So there is no 'full' calcMode here, only the Ritz solve.
%%%   2. A free plate's boundary conditions are the TRUE physical ones, not an
%%%      artificial truncation of an infinite domain the way the half-space's
%%%      rigid-bottom BC is. Every converged Ritz mode is therefore a genuine
%%%      Lamb branch, so there is no "which numerical mode is the real one"
%%%      question for ModeClassification to answer. Mode TYPE
%%%      (Asymm/Symm/TH-Asymm/TH-Symm) is a separate question, answered inline by
%%%      RitzSolverLamb.modeRecognition (ritzSolver.modesType,
%%%      .modeTypeConfidence).
%%%
%%% Handle class: Ritz eigenvectors can be large and the GUI viewers need
%%% write-back.
%%%
%%% GuliweRR v6.0

    properties (SetAccess = protected)
        structure       Structure       = Structure()   % rotated structure for this angle;
                                                           % CalcParams lives inside it already
        angle(1,1)      double          = 0             % propagation angle [deg]

        %%% plate-specific physical scaling (see RitzSolverLamb's own docblock) --
        %%% both optional; realWavVec == 0 means "work in normalised relDepth
        %%% units only, no physical-frequency conversion", matching
        %%% RitzSolverLamb's own default
        realThickness   double          = []            % physical plate thickness, one value per
                                                           % layer or scalar (uniform); same units as
                                                           % realWavVec^-1
        realWavVec(1,1) double          = 0              % physical wavevector k [rad/unit]

        %%% ritzSolver -- calcRitz() result; modeRecognition runs inline
        %%% (withModeTypes = true), so modesType/modeTypeConfidence/freqsPhysical are
        %%% populated as soon as calcRitz() returns
        ritzSolver      RitzSolverLamb  = RitzSolverLamb()
    end

    properties
        %%% Public so CalcAngSweepPlate can prune each angle's eigenvectors
        %%% directly between phases.
        retentionLevel  char    = 'full'   % 'full' (default) | 'none' -- no 'saw_only'
                                             % tier here: every Lamb mode is physically
                                             % real, so there is no significant subset
                                             % to prune to, only "keep everything" or
                                             % "discard eigenvectors, keep the cheap
                                             % scalars"
    end

    properties (SetAccess = protected, Transient)
        %%% Transient: not written by save() (see CalcSingleHS) -- empty after
        %%% load(); read it through ritzPrecompute(), which rebuilds it on first use.
        %%%
        %%% ritzSolver precompute, shared across a sweep when supplied.
        %%% ONLY nullConstrMat belongs here unconditionally -- it depends on
        %%% calcParams alone, not on angle or wavevector. A sweep also passes
        %%% rawByDepthCell (legOrder only, built by RitzSolverLamb.calc when absent).
        %%% gammaByDepthCell is deliberately NOT auto-built: GammaMatrix.elastByDepth
        %%% bakes in the ROTATED elastic tensor at a specific angle, so it is reusable
        %%% only across REPEATED WAVEVECTOR evaluations AT ONE FIXED ANGLE, never
        %%% across angles. A caller in that position may still pass one in; this
        %%% constructor cannot invent one, having no angle to invent it for.
        precompStructForRitz struct = struct()
    end

    methods
        function obj = CalcSinglePlate(structure, angle, realThickness, realWavVec, precompStruct, retentionLevel)
        %%% Construct CalcSinglePlate for one propagation angle
        %%% >> obj = CalcSinglePlate(structure, angle, realThickness, realWavVec, precompStruct, retentionLevel)
        %%%
        %%% INPUT:  1. structure       .. Structure object, UNROTATED (the constructor
        %%%                              rotates it internally -- callers never rotate
        %%%                              by hand). nLay >= 1 -- a free plate may be a
        %%%                              single homogeneous layer or a genuine
        %%%                              multi-layer stack (unlike CalcSingleLayered,
        %%%                              no nLay >= 2 requirement).
        %%%  (opt)  2. angle           .. propagation angle [deg] (default = 0)
        %%%  (opt)  3. realThickness   .. physical plate thickness (default = [],
        %%%                              i.e. normalised relDepth units only)
        %%%  (opt)  4. realWavVec      .. physical wavevector k (default = 0, same
        %%%                              meaning)
        %%%  (opt)  5. precompStruct   .. angle-independent Ritz precompute field
        %%%                              (nullConstrMat) -- built from the rotated
        %%%                              structure if omitted, so a bare
        %%%                              CalcSinglePlate(structure, angle) is fully
        %%%                              usable standalone. May also carry
        %%%                              gammaByDepthCell (see precompStructForRitz's
        %%%                              own comment) -- passed through as-is, never
        %%%                              auto-built here.
        %%%  (opt)  6. retentionLevel  .. 'full' (default) | 'none'

            if nargin == 0
                return
            end

            validateattributes(structure, {'Structure'}, {'scalar'}, 'CalcSinglePlate', 'structure');
            if nargin < 2 || isempty(angle); angle = 0; end
            validateattributes(angle, {'double'}, {'scalar', 'real', 'finite'}, 'CalcSinglePlate', 'angle');
            if nargin >= 6 && ~isempty(retentionLevel)
                CalcCommon.validateRetentionLevel(retentionLevel, 'CalcSinglePlate:constructor', {'full', 'none'});
            end

            obj.angle     = angle;
            obj.structure = CalcCommon.rotateStructure(structure, angle);

            if nargin >= 3 && ~isempty(realThickness)
                obj.realThickness = realThickness;
            end
            if nargin >= 4 && ~isempty(realWavVec)
                obj.realWavVec = realWavVec;
            end

            if nargin < 5 || isempty(precompStruct)
                precompStruct = struct();
            end
            if ~isfield(precompStruct, 'nullConstrMat') || isempty(precompStruct.nullConstrMat)
                precompStruct.nullConstrMat = RitzSolver.buildConstraintMatrix(obj.structure.calcParams, 'lamb');
            end
            obj.precompStructForRitz = precompStruct;

            if nargin >= 6 && ~isempty(retentionLevel)
                obj.retentionLevel = retentionLevel;
            end
        end

        function obj = calcRitz(obj, precomp)
        %%% Run RitzSolverLamb for this angle
        %%% >> obj = calcRitz(obj)
        %%% >> obj = calcRitz(obj, precomp)
        %%%
        %%% INPUT: (opt) precomp .. the sweep's shared Ritz precompute (nullConstrMat);
        %%%                         default: this object's own (ritzPrecompute()). A
        %%%                         sweep always passes it: this object's own copy is
        %%%                         Transient and so empty on a thread-pool worker (see
        %%%                         CalcSingleHS.calcRitz). Not stored on obj.
        %%%
        %%% withModeTypes = true -- mode-type classification happens inline.
        %%%
        %%% obj.structure is already rotated by the constructor, so pass 0 here,
        %%% not obj.angle, or GammaMatrix.elastByDepth rotates it again.
            if nargin < 2 || isempty(precomp)
                precomp = obj.ritzPrecompute();
            elseif ~isfield(precomp, 'nullConstrMat') || isempty(precomp.nullConstrMat)
                precomp.nullConstrMat = RitzSolver.buildConstraintMatrix(obj.structure.calcParams, 'lamb');
                CalcCommon.warnRitzPrecomputeBuiltOnWorker(class(obj));
            end
            rs = RitzSolverLamb(obj.structure, 0, obj.realThickness, obj.realWavVec);
            obj.ritzSolver = rs.calc(true, precomp);
        end

        function precomp = ritzPrecompute(obj)
        %%% The Ritz precompute, nullConstrMat rebuilt first if missing (e.g. after load())
        %%% >> precomp = ritzPrecompute(obj)
        %%%
        %%% OUTPUT: precomp .. precompStructForRitz with nullConstrMat populated
        %%%
        %%% Only nullConstrMat is rebuilt -- the one field the constructor itself
        %%% guarantees (a caller-supplied gammaByDepthCell is angle-specific and
        %%% optional, so a missing one is simply left missing). A rebuilt field
        %%% is stored back on obj.
            if ~isfield(obj.precompStructForRitz, 'nullConstrMat') ...
                    || isempty(obj.precompStructForRitz.nullConstrMat)
                obj.precompStructForRitz.nullConstrMat = ...
                    RitzSolver.buildConstraintMatrix(obj.structure.calcParams, 'lamb');
                CalcCommon.warnRitzPrecomputeBuiltOnWorker(class(obj));
            end
            precomp = obj.precompStructForRitz;
        end

        function obj = calcAll(obj, retentionLevel, precomp)
        %%% Run the (only) pipeline for this angle, then prune
        %%% >> obj = calcAll(obj, retentionLevel)
        %%% >> obj = calcAll(obj, retentionLevel, precomp)
        %%%
        %%% INPUT: (opt) 1. retentionLevel .. 'full' (default) | 'none'
        %%%        (opt) 2. precomp        .. shared Ritz precompute, passed to calcRitz()
        %%%                                   (default: this object's own)
        %%%
        %%% calcRitz() -> applyRetention. No mode argument, unlike
        %%% CalcSingleHS/CalcSingleLayered's calcAll(mode, retentionLevel): there
        %%% is only one pipeline here.

            if nargin < 2 || isempty(retentionLevel); retentionLevel = obj.retentionLevel; end
            CalcCommon.validateRetentionLevel(retentionLevel, 'CalcSinglePlate:calcAll', ...
                {'full', 'none'});   %%% no 'saw_only' tier -- see the property comment

            if nargin < 3; precomp = []; end
            obj = obj.calcRitz(precomp);
            obj = obj.applyRetention(retentionLevel);
            obj.retentionLevel = retentionLevel;
        end

        function obj = applyRetention(obj, retentionLevel)
        %%% Apply eigenvector retention policy to ritzSolver (RAM management)
        %%% >> obj = applyRetention(obj, retentionLevel)
        %%%
        %%% 'full': keep ritzSolver.eigVecs_con untouched.
        %%% 'none': clear ritzSolver.eigVecs_con entirely (vels, modesType,
        %%%         modeTypeConfidence, freqsPhysical, velsGroup all stay --
        %%%         cheap scalars, not the large [strBVL x nModes] matrix).
        %%%
        %%% Public so CalcAngSweepPlate can prune each angle's eigenvectors
        %%% directly between phases.
            switch retentionLevel
                case 'full'
                    %%% nothing to do
                case 'none'
                    obj.ritzSolver = obj.ritzSolver.pruneEigVecs([]);
                otherwise
                    %%% NOT routed through CalcCommon.pruneRitzTo: that helper
                    %%% ACCEPTS 'saw_only', and a plate must REJECT it.
                    CalcCommon.validateRetentionLevel(retentionLevel, ...
                        'CalcSinglePlate:applyRetention', {'full', 'none'});
            end
        end
    end

    methods (Access = protected)
        function cp = copyElement(obj)
        %%% Deep-copy override (matlab.mixin.Copyable)
        %%% >> cp = copyElement(obj)
        %%% All current properties are value types (Structure, RitzSolverLamb,
        %%% double, char, struct), so this is equivalent to the default shallow
        %%% copy; kept as an explicit override so any future handle-typed property
        %%% is not silently shallow-copied.
            cp = copyElement@matlab.mixin.Copyable(obj);
        end
    end
end
