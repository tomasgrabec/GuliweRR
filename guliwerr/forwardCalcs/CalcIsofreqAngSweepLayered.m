classdef CalcIsofreqAngSweepLayered < handle & matlab.mixin.Copyable
%%% CalcIsofreqAngSweepLayered — angular sweep of CalcIsofreqSingleLayered
%%% objects at FIXED reference frequency omega0. Rotates elasticity per angle
%%% and runs one CalcIsofreqSingleLayered (fixed omega0, internal sweep over
%%% k∥) at each, delivering raw per-angle isoFreq spectra — NOT peak-tracked
%%% into connected contours. No Ritz and no phase split, so calc() is a single
%%% flat loop rather than CalcAngSweepLayered's three-phase structure. See
%%% docs/reference/architecture/forwardCalcs.md, "isoWavevec vs. isoFreq".
%%%
%%% GuliweRR v6.0

    properties (SetAccess = protected)
        %%% Inputs
        structure   Structure   = Structure()      % UNROTATED base structure -- each
                                                       % cs{iAngle} rotates its own internal copy
        angles      double      = (0:2:180)'          % [nAngles x 1] propagation angles [deg]
        omega0(1,1) double      = NaN                 % fixed reference frequency -- REQUIRED,
                                                       % no default
        kParRange   double      = []                  % [nK x 1] wavevector-magnitude axis --
                                                       % REQUIRED, shared across every angle, a
                                                       % direct user input rather than derived
                                                       % from any phase's results
        useParallel(1,1) string {mustBeMember(useParallel, ["auto","always","never"])} = "auto"
                                                       % "auto" (default, parallel once nAngles >=
                                                       % ParallelPolicy.AUTO_THRESHOLD) | "always"
                                                       % | "never" -- thread-pool based, not
                                                       % process-based

        cs          cell        = {}                  % {nAngles x 1} of CalcIsofreqSingleLayered

        %%% Convenience extraction -- filled by assembleMaps() at the end of calc()
        absGMap     double      = []                  % [nK x nAngles]
    end

    methods
        function obj = CalcIsofreqAngSweepLayered(structure, angles, omega0, kParRange, useParallel)
        %%% Construct CalcIsofreqAngSweepLayered for an angular sweep at fixed omega0
        %%% >> obj = CalcIsofreqAngSweepLayered(structure, angles, omega0, kParRange, useParallel)
        %%%
        %%% INPUT:  1. structure  .. Structure object, UNROTATED (each cs{iAngle} rotates
        %%%                         its own internal copy -- callers never rotate by
        %%%                         hand, same convention as CalcIsofreqSingleLayered);
        %%%                         any nLay >= 1 (half-space or layered --
        %%%                         GreenSolver.calc is layered-generic)
        %%%  (opt)  2. angles     .. propagation angles [deg] (default = (0:2:180)')
        %%%         3. omega0     .. fixed reference frequency -- REQUIRED, no default
        %%%         4. kParRange  .. [nK x 1] wavevector-magnitude axis to sweep --
        %%%                         REQUIRED, shared across every angle
        %%%  (opt)  5. useParallel .. "auto" (default) | "always" | "never"

            if nargin == 0
                return
            end

            validateattributes(structure, {'Structure'}, {'scalar'}, 'CalcIsofreqAngSweepLayered', 'structure');
            if nargin < 2 || isempty(angles); angles = (0:2:180)'; end
            validateattributes(angles, {'double'}, {'vector', 'real', 'finite'}, 'CalcIsofreqAngSweepLayered', 'angles');

            if nargin < 3 || isempty(omega0)
                error('CalcIsofreqAngSweepLayered:constructor:missingOmega0', ...
                    'omega0 is required (no default) -- a CalcIsofreqAngSweepLayered with no frequency reference is meaningless.');
            end
            validateattributes(omega0, {'double'}, {'scalar', 'real', 'finite', 'positive'}, 'CalcIsofreqAngSweepLayered', 'omega0');

            if nargin < 4 || isempty(kParRange)
                error('CalcIsofreqAngSweepLayered:constructor:missingKParRange', ...
                    'kParRange is required (no default) -- the wavevector-magnitude axis to sweep, shared across every angle.');
            end
            validateattributes(kParRange, {'double'}, {'vector', 'real', 'finite', 'positive'}, 'CalcIsofreqAngSweepLayered', 'kParRange');

            obj.structure = structure;      % UNROTATED -- kept as-is on this object
            obj.angles    = angles(:);
            obj.omega0    = omega0;
            obj.kParRange = kParRange(:);
            if nargin >= 5 && ~isempty(useParallel); obj.useParallel = useParallel; end

            nAngles = numel(obj.angles);
            obj.cs  = cell(nAngles, 1);
            for iAngle = 1:nAngles
                obj.cs{iAngle} = Calc.isofreqSingleLayered(structure, obj.angles(iAngle), obj.omega0, obj.kParRange);
                %%% Nothing to precompute and share across angles here: with no Ritz
                %%% and no GammaMatrix, there is no analogue of
                %%% nullConstrMat/rawByDepthCell. Each construction is independent, and
                %%% CalcIsofreqSingleLayered's own constructor does the per-angle
                %%% rotation -- this class never rotates anything itself.
            end
        end

        function obj = calc(obj)
        %%% Run the isoFreq sweep: one CalcIsofreqSingleLayered.calc() per angle, then
        %%% assemble the convenience maps.
        %%% >> obj = calc(obj)
        %%%
        %%% No phase split, no shared-axis preamble, no retention policy -- every
        %%% angle's calculation is fully independent, so this is a single flat loop
        %%% rather than CalcAngSweepLayered's three-phase structure. cs holds
        %%% VALUE-class objects, so there is no handle-aliasing concern either.
            useP   = CalcCommon.resolveParallel(obj.useParallel, numel(obj.cs));
            obj.cs = CalcCommon.forEachElement(obj.cs, useP, @(cs, ~) cs.calc());
            obj = obj.assembleMaps();
        end

        function obj = setElement(obj, iAngle, csElem)
        %%% Write a single per-angle result back into cs -- cs is SetAccess=protected,
        %%% so this is the narrow accessor CalcFreqAngSweep2Layered's flat-dispatch
        %%% path uses to scatter a flat-parfor result back in.
        %%% >> obj = setElement(obj, iAngle, csElem)
            obj.cs{iAngle} = csElem;
        end

        function obj = finalizeMaps(obj)
        %%% Public wrapper around the private assembleMaps() -- lets
        %%% CalcFreqAngSweep2Layered's flat-dispatch path trigger map assembly after
        %%% scattering results back via setElement.
        %%% >> obj = finalizeMaps(obj)
            obj = obj.assembleMaps();
        end
    end

    methods (Access = private)
        function obj = assembleMaps(obj)
        %%% Extract the sweep-level convenience map from cs.
        %%% >> obj = assembleMaps(obj)
            obj.absGMap = cell2mat(cellfun(@(c) c.absG(:), obj.cs(:)', 'UniformOutput', false));
        end
    end
end
