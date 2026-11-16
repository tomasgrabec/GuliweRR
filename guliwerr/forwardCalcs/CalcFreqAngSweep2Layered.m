classdef CalcFreqAngSweep2Layered < handle & matlab.mixin.Copyable
%%% CalcFreqAngSweep2Layered — 2D sweep (propagation angle AND real reference
%%% frequency omega0) for the layered forward calculation, isoFreq side. The
%%% counterpart to CalcWavevecAngSweep2Layered: constructs one
%%% CalcIsofreqAngSweepLayered per omega0 value and drives it through its own
%%% calc(). See docs/reference/architecture/forwardCalcs.md, "2D sweep orchestrators",
%%% "Precompute sharing differs between the two" and "`retentionLevel='none'`
%%% default, and no top-level aggregation".
%%%
%%% Handle class; implements copyElement() (sweeps holds handle-class
%%% CalcIsofreqAngSweepLayered elements that would otherwise alias on copy).
%%%
%%% GuliweRR v6.0

    properties (SetAccess = protected)
        %%% Inputs
        structure   Structure   = Structure()      % UNROTATED base structure -- each
                                                       % sweeps{iOmega} rotates its own
                                                       % internal copies per angle, callers
                                                       % never rotate by hand
        angles      double      = (0:2:180)'          % [nAngles x 1] propagation angles [deg]
        omega0List  double      = []                  % [nOmega x 1] fixed reference frequencies --
                                                       % REQUIRED, no default
        kParRange   double      = []                  % [nK x 1] wavevector-magnitude axis, shared
                                                       % across every (angle, omega0) point
        useParallel(1,1) string {mustBeMember(useParallel, ["auto","always","never"])} = "auto"
                                                       % "auto" (default) | "always" | "never" --
                                                       % forwarded UNCHANGED to every sweeps{iOmega},
                                                       % which decides at the angle level; this
                                                       % class's outer iOmega loop stays serial,
                                                       % MATLAB forbidding nested parfor

        %%% Per-omega0 containers -- one CalcIsofreqAngSweepLayered per omega0
        sweeps      cell        = {}                  % {nOmega x 1} of CalcIsofreqAngSweepLayered
    end

    properties
        %%% flattenSweep -- settable any time before calc(). false (default): nested
        %%% behaviour, outer iOmega loop serial, each sweeps{iOmega} parallelising its
        %%% own angle loop once nAngles crosses ParallelPolicy.AUTO_THRESHOLD. true:
        %%% dispatch across the flattened (omega0, angle) grid in ONE parfor, decided
        %%% from nOmega*nAngles rather than nAngles alone -- useful when neither axis
        %%% alone is large enough but their product is. This family has NO cross-angle
        %%% dependency (no phase split, no shared axis), so flattening is exact here.
        flattenSweep(1,1) logical = false
    end

    methods
        function obj = CalcFreqAngSweep2Layered(structure, angles, omega0List, kParRange, useParallel)
        %%% Construct CalcFreqAngSweep2Layered for a 2D (angle x omega0) sweep
        %%% >> obj = CalcFreqAngSweep2Layered(structure, angles, omega0List, kParRange, useParallel)
        %%%
        %%% INPUT:  1. structure   .. Structure object, UNROTATED (each sweeps{iOmega}'s
        %%%                          own CalcIsofreqAngSweepLayered/CalcIsofreqSingleLayered
        %%%                          rotates internally -- callers never rotate by hand);
        %%%                          any nLay >= 1 (half-space or layered --
        %%%                          GreenSolver.calc is layered-generic)
        %%%  (opt)  2. angles      .. propagation angles [deg] (default = (0:2:180)')
        %%%         3. omega0List  .. [nOmega x 1] fixed reference frequencies -- REQUIRED,
        %%%                          no default
        %%%         4. kParRange   .. [nK x 1] wavevector-magnitude axis to sweep --
        %%%                          REQUIRED, shared across every (angle, omega0) point
        %%%  (opt)  5. useParallel .. "auto" (default) | "always" | "never" -- forwarded
        %%%                          unchanged to every sweeps{iOmega}

            if nargin == 0
                return
            end

            validateattributes(structure, {'Structure'}, {'scalar'}, 'CalcFreqAngSweep2Layered', 'structure');

            if nargin < 2 || isempty(angles); angles = (0:2:180)'; end
            validateattributes(angles, {'double'}, {'vector', 'real', 'finite'}, ...
                'CalcFreqAngSweep2Layered', 'angles');

            %%% omega0List is REQUIRED -- caught here so the caller gets a
            %%% project-standard error ID, not validateattributes' generic message.
            if nargin < 3 || isempty(omega0List)
                error('CalcFreqAngSweep2Layered:constructor:missingOmega0List', ...
                    'omega0List is required (no default) -- a CalcFreqAngSweep2Layered with no frequency reference is meaningless.');
            end
            validateattributes(omega0List, {'double'}, {'vector', 'real', 'finite', 'positive'}, ...
                'CalcFreqAngSweep2Layered', 'omega0List');

            if nargin < 4 || isempty(kParRange)
                error('CalcFreqAngSweep2Layered:constructor:missingKParRange', ...
                    'kParRange is required (no default) -- the wavevector-magnitude axis to sweep, shared across every (angle, omega0) point.');
            end
            validateattributes(kParRange, {'double'}, {'vector', 'real', 'finite', 'positive'}, ...
                'CalcFreqAngSweep2Layered', 'kParRange');

            obj.structure  = structure;      % UNROTATED -- kept as-is on this object
            obj.angles     = angles(:);
            obj.omega0List = omega0List(:);
            obj.kParRange  = kParRange(:);
            if nargin >= 5 && ~isempty(useParallel); obj.useParallel = useParallel; end
        end

        function obj = calc(obj)
        %%% Run the 2D sweep: one CalcIsofreqAngSweepLayered per omega0, leaving the
        %%% per-omega0 spectra distributed with no top-level aggregation.
        %%% >> obj = calc(obj)
        %%%
        %%% Plain nested loop, no shared precompute -- each sweeps{iOmega} is
        %%% constructed and driven through its own calc(). It never constructs
        %%% CalcIsofreqSingleLayered/GreenSolver directly and never rotates a Structure
        %%% itself. The outer iOmega loop stays serial (MATLAB forbids nested parfor);
        %%% obj.useParallel is forwarded to each sweeps{iOmega} instead.
        %%%
        %%% With obj.flattenSweep true, calcFlat() dispatches across the flattened
        %%% (omega0, angle) grid in one parfor instead.

            if obj.flattenSweep
                obj = obj.calcFlat();
                return
            end

            nOmega = numel(obj.omega0List);
            obj.sweeps = cell(nOmega, 1);
            for iOmega = 1:nOmega
                obj.sweeps{iOmega} = Calc.isofreqAngSweepLayered(obj.structure, obj.angles, ...
                    obj.omega0List(iOmega), obj.kParRange, obj.useParallel);
                obj.sweeps{iOmega} = obj.sweeps{iOmega}.calc();
            end
        end
    end

    methods (Access = private)
        function obj = calcFlat(obj)
        %%% Flat-dispatch variant of calc(), used when obj.flattenSweep is true.
        %%% Constructs every sweeps{iOmega} serially, flattens every
        %%% sweeps{iOmega}.cs{iAngle} into one [nOmega*nAngles x 1] list and runs ONE
        %%% parfor over the whole grid -- this family has no cross-angle dependency, so
        %%% unlike CalcAngSweepLayered there is no barrier to respect. Results scatter
        %%% back via setElement and finalizeMaps.
        %%% >> obj = calcFlat(obj)
            nOmega  = numel(obj.omega0List);
            nAngles = numel(obj.angles);
            obj.sweeps = cell(nOmega, 1);
            for iOmega = 1:nOmega
                obj.sweeps{iOmega} = Calc.isofreqAngSweepLayered(obj.structure, obj.angles, ...
                    obj.omega0List(iOmega), obj.kParRange, "never");
                %%% useParallel = "never" -- the flat dispatch below owns parallelism;
                %%% sweeps{iOmega}.calc() is never called in this path.
            end

            [csFlat, omegaIdx, angleIdx] = CalcCommon.flattenSweeps(obj.sweeps, nAngles);

            %%% Decided from nOmega*nAngles, not nAngles alone -- the point of the flat
            %%% path. Still goes through CalcIsofreqSingleLayered's own calc(), as the
            %%% nested path does.
            useP   = CalcCommon.resolveParallel(obj.useParallel, numel(csFlat));
            csFlat = CalcCommon.forEachElement(csFlat, useP, @(cs, ~) cs.calc());

            obj.sweeps = CalcCommon.scatterToSweeps(obj.sweeps, csFlat, omegaIdx, angleIdx);
            for iOmega = 1:nOmega
                obj.sweeps{iOmega} = obj.sweeps{iOmega}.finalizeMaps();
            end
        end
    end

    methods (Access = protected)
        function cp = copyElement(obj)
        %%% Deep-copy override (matlab.mixin.Copyable)
        %%% >> cp = copyElement(obj)
        %%% obj.sweeps holds handle-class CalcIsofreqAngSweepLayered objects, which the
        %%% default shallow copy would alias. One level of cellfun(@copy, ...) is
        %%% enough for a fully independent tree: their own cs{} arrays hold value-class
        %%% CalcIsofreqSingleLayered elements, which the default copy already handles.
            cp = copyElement@matlab.mixin.Copyable(obj);
            cp.sweeps = cellfun(@copy, obj.sweeps, 'UniformOutput', false);
        end
    end
end
