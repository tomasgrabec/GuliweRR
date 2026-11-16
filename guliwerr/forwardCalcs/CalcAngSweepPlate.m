classdef CalcAngSweepPlate < handle & matlab.mixin.Copyable
%%% CalcAngSweepPlate — angular sweep of CalcSinglePlate objects for the
%%% free-plate (Lamb-wave) forward calculation, at one FIXED thickness/
%%% wavevector. Mirrors CalcAngSweepHS/CalcAngSweepLayered's per-angle
%%% composition, simplified the same way CalcSinglePlate simplifies
%%% CalcSingleHS/CalcSingleLayered: no Christoffel/Green, and no "which mode is
%%% physically real" question.
%%%
%%% Also the per-wavevector building block CalcWavevecAngSweep2Plate composes
%%% across a wavevector array, as CalcAngSweepLayered is
%%% CalcWavevecAngSweep2Layered's per-wavelength slice. Standalone use is a
%%% single-thickness angular map of Lamb-mode velocities.
%%%
%%% Every angle shares the SAME number of computed modes — nCompModes depends
%%% only on calcParams (legOrder/layer count), never on angle or material
%%% rotation — unlike CalcAngSweepHS's ragged per-angle SAW-mode count. The
%%% convenience maps below are therefore dense [nCompModes x nAngles] matrices,
%%% not cell arrays of variable-length vectors.
%%%
%%% GuliweRR v6.0

    properties (SetAccess = protected)
        structure       Structure   = Structure()      % UNROTATED base structure -- each
                                                          % cs{iAngle} rotates its own internal copy
        angles          double      = (0:2:180)'          % [nAngles x 1] propagation angles [deg].
                                                          % 0-180 only, on the same
                                                          % Christoffel-is-quadratic-in-n symmetry
                                                          % argument as CalcAngSweepHS

        %%% plate-specific physical scaling, FIXED across the whole sweep --
        %%% same meaning as CalcSinglePlate's own properties
        realThickness   double      = []
        realWavVec(1,1) double      = 0

        retentionLevel  char        = 'full'   % 'full' (default) | 'none' -- see
                                                 % CalcSinglePlate.m for why no 'saw_only' tier
        useParallel(1,1) string {mustBeMember(useParallel, ["auto","always","never"])} = "auto"

        %%% cs -- per-angle CalcSinglePlate array
        cs              cell        = {}                 % {nAngles x 1} of CalcSinglePlate objects

        %%% maps -- filled by assembleMaps() at the end of calc(). Dense
        %%% [nCompModes x nAngles], unlike CalcAngSweepHS's ragged cell maps.
        velsMap                 double = []   % [nCompModes x nAngles] phase velocities
        velsGroupMap            double = []   % [nCompModes x nAngles] group velocities (along k)
        velsGroupTransMap       double = []   % [nCompModes x nAngles] transverse group velocities
                                               % (perpendicular to k), RitzSolverLamb.velsGroupTrans
        modesTypeMap            double = []   % [nCompModes x nAngles] 1=Asymm,2=Symm,3=TH-Asymm,4=TH-Symm
        modeTypeConfidenceMap   double = []   % [nCompModes x nAngles] see RitzSolverLamb.modeTypeConfidence
        freqsPhysicalMap        double = []   % [nCompModes x nAngles] -- only when realWavVec ~= 0
    end

    properties (SetAccess = protected, Transient)
        %%% Transient: not written by save(), empty after load() -- calc() rebuilds it
        %%% then, once. calc() hands it to every angle explicitly: a thread-pool worker
        %%% receives each cs{iAngle} without its own Transient copy (see CalcAngSweepHS).
        %%%
        %%% Angle-independent Ritz precompute, built ONCE and shared by every cs{iAngle}:
        %%% nullConstrMat and rawByDepthCell (GammaMatrix.RawCellFor, the Legendre
        %%% integrals; both depend on legOrder alone, not on angle, material or depth).
        %%% gammaByDepthCell from GammaMatrix.elastByDepth bakes in the ROTATED elastic
        %%% tensor at a specific angle and is not shareable across angles; each
        %%% cs{iAngle} builds its own.
        precompStructForRitz struct = struct()
    end

    methods
        function obj = CalcAngSweepPlate(structure, angles, realThickness, realWavVec, retentionLevel, useParallel)
        %%% Construct CalcAngSweepPlate for an angular sweep at fixed thickness/wavevector
        %%% >> obj = CalcAngSweepPlate(structure, angles, realThickness, realWavVec, retentionLevel, useParallel)
        %%%
        %%% INPUT:  1. structure       .. Structure object, UNROTATED (each cs{iAngle}
        %%%                              rotates its own internal copy). nLay >= 1.
        %%%  (opt)  2. angles          .. propagation angles [deg] (default = (0:2:180)')
        %%%  (opt)  3. realThickness   .. physical plate thickness (default = [])
        %%%  (opt)  4. realWavVec      .. physical wavevector k (default = 0)
        %%%  (opt)  5. retentionLevel  .. 'full' (default) | 'none'
        %%%  (opt)  6. useParallel     .. "auto" (default) | "always" | "never"

            if nargin == 0
                return
            end

            validateattributes(structure, {'Structure'}, {'scalar'}, 'CalcAngSweepPlate', 'structure');
            if nargin < 2 || isempty(angles); angles = (0:2:180)'; end
            validateattributes(angles, {'double'}, {'vector', 'real', 'finite'}, 'CalcAngSweepPlate', 'angles');
            if nargin >= 5 && ~isempty(retentionLevel)
                CalcCommon.validateRetentionLevel(retentionLevel, 'CalcAngSweepPlate:constructor', {'full', 'none'});
            end

            obj.structure = structure;      % UNROTATED -- kept as-is on this object
            obj.angles    = angles(:);

            if nargin >= 3 && ~isempty(realThickness); obj.realThickness = realThickness; end
            if nargin >= 4 && ~isempty(realWavVec);    obj.realWavVec    = realWavVec;    end
            if nargin >= 5 && ~isempty(retentionLevel); obj.retentionLevel = retentionLevel; end
            if nargin >= 6 && ~isempty(useParallel);    obj.useParallel    = useParallel;    end

            %%% Build the angle-independent Ritz precompute ONCE.
            obj.precompStructForRitz.nullConstrMat = ...
                RitzSolver.buildConstraintMatrix(structure.calcParams, 'lamb');
            obj.precompStructForRitz.rawByDepthCell = GammaMatrix.RawCellFor(structure);

            nAngles = numel(obj.angles);
            obj.cs  = cell(nAngles, 1);
            for iAngle = 1:nAngles
                obj.cs{iAngle} = Calc.singlePlate(structure, obj.angles(iAngle), ...
                    obj.realThickness, obj.realWavVec, obj.precompStructForRitz, obj.retentionLevel);
                %%% CalcSinglePlate's constructor sees nullConstrMat already populated
                %%% and skips rebuilding it -- genuinely shared, not rebuilt per angle.
            end
        end

        function obj = calc(obj, retentionLevel)
        %%% Run the sweep (Ritz only, every angle), then assemble convenience maps
        %%% >> obj = calc(obj, retentionLevel)
        %%%
        %%% INPUT: (opt) retentionLevel .. 'full' (default) | 'none'

            if nargin < 2 || isempty(retentionLevel); retentionLevel = obj.retentionLevel; end
            CalcCommon.validateRetentionLevel(retentionLevel, 'CalcAngSweepPlate:calc', ...
                {'full', 'none'});   %%% no 'saw_only' tier -- see CalcSinglePlate.m
            obj.retentionLevel = retentionLevel;

            %%% Single flat loop, no phase split: every angle's Lamb solve is fully
            %%% independent (no Green/Christoffel shared-axis barrier exists here),
            %%% sharing only the precompute, which goes in explicitly.
            if ~isfield(obj.precompStructForRitz, 'nullConstrMat') || isempty(obj.precompStructForRitz.nullConstrMat)
                obj.precompStructForRitz.nullConstrMat = ...   % after load()
                    RitzSolver.buildConstraintMatrix(obj.structure.calcParams, 'lamb');
            end
            if ~isfield(obj.precompStructForRitz, 'rawByDepthCell') || isempty(obj.precompStructForRitz.rawByDepthCell)
                obj.precompStructForRitz.rawByDepthCell = GammaMatrix.RawCellFor(obj.structure);
            end
            precompLocal = obj.precompStructForRitz;
            useP   = CalcCommon.resolveParallel(obj.useParallel, numel(obj.cs));
            obj.cs = CalcCommon.forEachElement(obj.cs, useP, @(cs, ~) cs.calcAll(retentionLevel, precompLocal));

            obj = obj.assembleMaps();
        end
    end

    methods (Access = private)
        function obj = assembleMaps(obj)
        %%% Extract the sweep-level convenience maps from cs -- no new logic.
        %%% >> obj = assembleMaps(obj)
            nAngles = numel(obj.cs);
            nModes  = obj.cs{1}.ritzSolver.nCompModes;

            obj.velsMap               = NaN(nModes, nAngles);
            obj.velsGroupMap          = NaN(nModes, nAngles);
            obj.velsGroupTransMap     = NaN(nModes, nAngles);
            obj.modesTypeMap          = NaN(nModes, nAngles);
            obj.modeTypeConfidenceMap = NaN(nModes, nAngles);
            hasPhysicalFreqs = obj.realWavVec ~= 0;
            if hasPhysicalFreqs
                obj.freqsPhysicalMap = NaN(nModes, nAngles);
            end

            for iAngle = 1:nAngles
                rs = obj.cs{iAngle}.ritzSolver;
                obj.velsMap(:, iAngle)               = rs.vels(1:nModes);
                obj.velsGroupMap(:, iAngle)           = rs.velsGroup(1:nModes);
                obj.velsGroupTransMap(:, iAngle)      = rs.velsGroupTrans(1:nModes);
                obj.modesTypeMap(:, iAngle)           = rs.modesType(1:nModes);
                obj.modeTypeConfidenceMap(:, iAngle)  = rs.modeTypeConfidence(1:nModes);
                if hasPhysicalFreqs
                    obj.freqsPhysicalMap(:, iAngle) = rs.freqsPhysical(1:nModes);
                end
            end
        end
    end

    methods (Access = protected)
        function cp = copyElement(obj)
        %%% Deep-copy override (matlab.mixin.Copyable)
        %%% >> cp = copyElement(obj)
        %%% obj.cs holds handle-class CalcSinglePlate objects, which the default
        %%% shallow copy would alias -- deep-copy each element explicitly.
            cp = copyElement@matlab.mixin.Copyable(obj);
            cp.cs = cellfun(@copy, obj.cs, 'UniformOutput', false);
        end
    end
end
