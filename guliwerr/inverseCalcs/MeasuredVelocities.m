classdef MeasuredVelocities
%%% MeasuredVelocities — the measured dataset from one TGS angular sweep, and nothing else.
%%%   Immutable, calculation-free: what was measured, with what uncertainty, in what units.
%%%   No mode number, no eigenvector, no calculated velocity -- that binding lives in
%%%   ModeAssociation. Optionally carries independently measured pulse-echo bulk phase
%%%   velocities and a single soft-fixed bulk modulus.
%%%
%%% EXPERIMENTAL: validated on synthetic data only; the interface may change.
%%%
%%% Value class: addPhaseVelocity/setBulkModulus return the modified object.
%%%
%%% JOINT MULTI-DATASET INVERSION: one object per dataset. For a joint inversion of several
%%%   crystallographic cuts of the SAME crystal, pass an ARRAY of MeasuredVelocities (one per
%%%   dataset) to SensitivityAnalysis.calc/ElasticityInversion -- see those classes' docblocks.
%%%
%%% GuliweRR v6.0

    properties (Constant)
        %%% Fixed wave-type vocabulary, matching ModeClassification's `type` names
        %%% (plus 'sezawa'/'love' for layered systems and 'unknown' as an escape hatch).
        WAVE_TYPES = {'longitudinal','fastShear','slowShear','rayleigh', ...
                      'pseudoRayleigh','sezawa','love','unknown'}
    end

    properties (SetAccess = protected)
        angleDeg        double = zeros(0,1)   % [nPoints x 1] propagation angle [deg]
        velocity        double = zeros(0,1)   % [nPoints x 1] measured phase velocity [km/s]
        waveType        cell   = {}           % [nPoints x 1] cellstr, from WAVE_TYPES
        velocityStd     double = zeros(0,1)   % [nPoints x 1] measurement std [km/s]; NaN = unknown.
                                              % Metadata only: no calculation reads it (weights
                                              % are the PointWeights option of
                                              % SensitivityAnalysis.calc/ElasticityInversion)

        phaseVelocities struct = struct('propagationDir', {}, 'waveType', {}, 'velocity', {})
                                        % pulse-echo bulk phase velocities (optional channel)
        bulkModulus(1,1) double = NaN        % soft-fixed bulk modulus [GPa]; NaN = not set

        nPoints(1,1)    double = 0
    end

    methods
        function obj = MeasuredVelocities(angleDeg, velocity, waveType, velocityStd)
        %%% Construct the measured-velocity dataset
        %%% >> obj = MeasuredVelocities(angleDeg, velocity, waveType, velocityStd)
        %%%
        %%% INPUT:  1. angleDeg    .. [N x 1] or [1 x N] propagation angles [deg]
        %%%         2. velocity    .. [N x 1] or [1 x N] measured phase velocities [km/s]
        %%%  (opt)  3. waveType    .. cellstr [N x 1], or a single char/string applied to
        %%%                          every point; member of MeasuredVelocities.WAVE_TYPES.
        %%%                          Default = 'unknown' for every point.
        %%%  (opt)  4. velocityStd .. [N x 1], or a scalar applied to every point.
        %%%                          Default = NaN (unknown) for every point. Stored as
        %%%                          metadata only; it does not weight the fit.

            if nargin == 0
                return
            end

            validateattributes(angleDeg, {'double'}, {'vector', 'real', 'finite'}, ...
                'MeasuredVelocities', 'angleDeg');
            validateattributes(velocity, {'double'}, {'vector', 'real', 'finite', 'positive'}, ...
                'MeasuredVelocities', 'velocity');

            angleDeg = angleDeg(:);
            velocity = velocity(:);
            n = numel(angleDeg);
            if numel(velocity) ~= n
                error('MeasuredVelocities:constructor:sizeMismatch', ...
                    'angleDeg (%d) and velocity (%d) must have the same length.', n, numel(velocity));
            end

            if nargin < 3 || isempty(waveType)
                waveType = repmat({'unknown'}, n, 1);
            elseif ischar(waveType) || (isstring(waveType) && isscalar(waveType))
                waveType = repmat({char(waveType)}, n, 1);
            else
                waveType = cellstr(waveType(:));
                if numel(waveType) ~= n
                    error('MeasuredVelocities:constructor:sizeMismatch', ...
                        'waveType (%d) must have the same length as angleDeg (%d).', numel(waveType), n);
                end
            end
            badType = ~ismember(waveType, MeasuredVelocities.WAVE_TYPES);
            if any(badType)
                error('MeasuredVelocities:constructor:badWaveType', ...
                    'waveType entry %d ("%s") is not one of: %s.', ...
                    find(badType,1), waveType{find(badType,1)}, strjoin(MeasuredVelocities.WAVE_TYPES, ', '));
            end

            if nargin < 4 || isempty(velocityStd)
                velocityStd = NaN(n, 1);
            elseif isscalar(velocityStd)
                velocityStd = repmat(velocityStd, n, 1);
            else
                velocityStd = velocityStd(:);
                if numel(velocityStd) ~= n
                    error('MeasuredVelocities:constructor:sizeMismatch', ...
                        'velocityStd (%d) must have length 1 or n (%d).', numel(velocityStd), n);
                end
            end

            obj.angleDeg    = angleDeg;
            obj.velocity    = velocity;
            obj.waveType    = waveType;
            obj.velocityStd = velocityStd;
            obj.nPoints     = n;
        end

        function obj = addPhaseVelocity(obj, propagationDir, waveType, velocity)
        %%% Append one independently measured pulse-echo bulk phase velocity
        %%% >> obj = addPhaseVelocity(obj, propagationDir, waveType, velocity)
        %%%
        %%% INPUT: 1. propagationDir .. [1x3] or [3x1] propagation direction in the inverted
        %%%                             layer's CRYSTAL frame, not the dirMat/lab frame (need
        %%%                             not be unit length; normalised on use)
        %%%        2. waveType       .. 'longitudinal' | 'fastShear' | 'slowShear' (a bulk
        %%%                             wave: the Christoffel branch the channel evaluates)
        %%%        3. velocity       .. measured bulk phase velocity [km/s]

            validateattributes(propagationDir, {'double'}, {'vector', 'numel', 3, 'real', 'finite'}, ...
                'MeasuredVelocities.addPhaseVelocity', 'propagationDir');
            waveType = char(waveType);
            bulkTypes = {'longitudinal', 'fastShear', 'slowShear'};
            if ~ismember(waveType, bulkTypes)
                error('MeasuredVelocities:addPhaseVelocity:badWaveType', ...
                    '"%s" is not a bulk wave type (one of: %s).', waveType, strjoin(bulkTypes, ', '));
            end
            validateattributes(velocity, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
                'MeasuredVelocities.addPhaseVelocity', 'velocity');

            newEntry.propagationDir = propagationDir(:).' / norm(propagationDir);
            newEntry.waveType       = waveType;
            newEntry.velocity       = velocity;
            obj.phaseVelocities     = [obj.phaseVelocities, newEntry];
        end

        function obj = setBulkModulus(obj, bulkModulusGPa)
        %%% Set the soft-fixed bulk modulus datum
        %%% >> obj = setBulkModulus(obj, bulkModulusGPa)
            validateattributes(bulkModulusGPa, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
                'MeasuredVelocities.setBulkModulus', 'bulkModulusGPa');
            obj.bulkModulus = bulkModulusGPa;
        end
    end

    methods (Static)
        function obj = fromTable(tbl)
        %%% Build a MeasuredVelocities from a plain MATLAB table (e.g. an importer's export)
        %%% >> obj = MeasuredVelocities.fromTable(tbl)
        %%%
        %%% INPUT: tbl .. table with columns angleDeg, velocity, and optionally
        %%%               waveType, velocityStd

            validateattributes(tbl, {'table'}, {}, 'MeasuredVelocities.fromTable', 'tbl');
            if ~all(ismember({'angleDeg', 'velocity'}, tbl.Properties.VariableNames))
                error('MeasuredVelocities:fromTable:missingColumns', ...
                    'tbl must have columns "angleDeg" and "velocity" at minimum.');
            end

            waveType = [];
            if ismember('waveType', tbl.Properties.VariableNames)
                waveType = tbl.waveType;
            end
            velocityStd = [];
            if ismember('velocityStd', tbl.Properties.VariableNames)
                velocityStd = tbl.velocityStd;
            end

            obj = MeasuredVelocities(tbl.angleDeg, tbl.velocity, waveType, velocityStd);
        end
    end
end
