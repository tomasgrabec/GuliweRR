classdef Inv
%%% Namespace for GuliweRR's inverse-characterisation constructors.
%%% All methods are static -- do not instantiate this class.
%%%
%%% EXPERIMENTAL: validated on synthetic data only; the interface may change.
%%%
%%% Usage:  measured    = Inv.measured(anglesDeg, velocities, waveTypes, velocityStd)
%%%         association = Inv.association(sweep, measured, velocityTolerance, ...
%%%                           angleTolerance, invertedLayerNo)
%%%         sensitivity = Inv.sensitivity(structure, association, invertedLayerNo, measured, ...)
%%%         inversion   = Inv.elasticity(structure, measured, association, sensitivity, ...
%%%                           invertedLayerNo, ...)
%%%
%%% Thin varargin passthrough only -- no defaults live here, and none should be added:
%%% a default duplicated here can silently drift from the real constructor's own default.
%%% Keep this file in sync with MeasuredVelocities/ModeAssociation.fromSweep/
%%% SensitivityAnalysis.calc/ElasticityInversion whenever any of their signatures change.
%%%
%%% Argument order: the two calls do NOT share it. Inv.sensitivity takes
%%% (structure, association, invertedLayerNo, measured), Inv.elasticity takes
%%% (structure, measured, association, sensitivity, invertedLayerNo) -- measured and
%%% association swap places, and invertedLayerNo moves from third to last.
%%%
%%% JOINT MULTI-DATASET INVERSION: Inv.sensitivity/Inv.elasticity accept `structure`/
%%% `association`/`measured` as either scalars (single dataset) or object arrays (several
%%% crystallographic cuts of the SAME crystal, sharing one target elasticity); see
%%% SensitivityAnalysis's and ElasticityInversion's own docblocks.
%%%
%%% GuliweRR v6.0

    methods (Static)
        function obj = measured(varargin)
            %%% >> obj = Inv.measured(anglesDeg, velocities, waveTypes, velocityStd)
            obj = MeasuredVelocities(varargin{:});
        end

        function obj = association(varargin)
            %%% >> obj = Inv.association(sweep, measured, velocityTolerance, angleTolerance, invertedLayerNo)
            obj = ModeAssociation.fromSweep(varargin{:});
        end

        function obj = sensitivity(varargin)
            %%% >> obj = Inv.sensitivity(structure, association, invertedLayerNo, measured, ...)
            obj = SensitivityAnalysis.calc(varargin{:});
        end

        function obj = elasticity(varargin)
            %%% >> obj = Inv.elasticity(structure, measured, association, sensitivity, ...
            %%%              invertedLayerNo, ...)
            obj = ElasticityInversion(varargin{:});
        end
    end
end
