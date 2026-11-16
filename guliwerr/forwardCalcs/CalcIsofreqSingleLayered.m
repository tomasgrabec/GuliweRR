classdef CalcIsofreqSingleLayered
%%% CalcIsofreqSingleLayered — single-angle, single-frequency container for the
%%% isoFreq Green's function scan: fixes one angle and one real reference
%%% frequency omega0, then sweeps wavevector magnitude k∥ via
%%% GreenSolver.calc's scanMode='isoFreq'. No Ritz and no phase split — Ritz
%%% fixes |k∥|=1 by construction. See docs/reference/architecture/forwardCalcs.md,
%%% "isoWavevec vs. isoFreq".
%%%
%%% GuliweRR v6.0

    properties (SetAccess = protected)
        %%% Inputs
        structure   Structure   = Structure()   % rotated structure for this angle;
                                                  % caller never rotates by hand,
                                                  % same convention as CalcSingleLayered
        angle(1,1)  double      = 0             % propagation angle [deg]
        omega0(1,1) double      = NaN           % fixed reference frequency -- REQUIRED, no
                                                  % default; same velocity-like normalised
                                                  % convention as GreenSolver.calc's omega0
        kParRange   double      = []            % [nK x 1] wavevector-magnitude axis to sweep,
                                                  % the physically meaningful swept quantity
                                                  % here. Must be positive; need not be sorted,
                                                  % but SHOULD be ascending if
                                                  % GreenPeaks.find will be
                                                  % used on the result, since findpeaks(Y,X)
                                                  % requires increasing X

        %%% calc() result
        absG        double      = []            % [nK x 1] abs(G33) at each kParRange(i) --
                                                  % indexed by kParRange, NOT velRangeForGreen
    end

    methods
        function obj = CalcIsofreqSingleLayered(structure, angle, omega0, kParRange)
        %%% Construct CalcIsofreqSingleLayered for one propagation angle and
        %%% one fixed reference frequency
        %%% >> obj = CalcIsofreqSingleLayered(structure, angle, omega0, kParRange)
        %%%
        %%% INPUT:  1. structure  .. Structure object, UNROTATED (the constructor
        %%%                         rotates it internally -- callers never rotate
        %%%                         by hand); any nLay >= 1 (half-space or layered
        %%%                         -- GreenSolver.calc is layered-generic)
        %%%  (opt)  2. angle      .. propagation angle [deg] (default = 0)
        %%%         3. omega0     .. fixed reference frequency -- REQUIRED, no
        %%%                         default (a CalcIsofreqSingleLayered with no
        %%%                         frequency reference is meaningless)
        %%%         4. kParRange  .. [nK x 1] wavevector-magnitude axis to sweep
        %%%                         -- REQUIRED, no default

            if nargin == 0
                return
            end

            validateattributes(structure, {'Structure'}, {'scalar'}, 'CalcIsofreqSingleLayered', 'structure');
            if nargin < 2 || isempty(angle); angle = 0; end
            validateattributes(angle, {'double'}, {'scalar', 'real', 'finite'}, 'CalcIsofreqSingleLayered', 'angle');

            if nargin < 3 || isempty(omega0)
                error('CalcIsofreqSingleLayered:constructor:missingOmega0', ...
                    'omega0 is required (no default) -- a CalcIsofreqSingleLayered with no frequency reference is meaningless.');
            end
            validateattributes(omega0, {'double'}, {'scalar', 'real', 'finite', 'positive'}, 'CalcIsofreqSingleLayered', 'omega0');

            if nargin < 4 || isempty(kParRange)
                error('CalcIsofreqSingleLayered:constructor:missingKParRange', ...
                    'kParRange is required (no default) -- the wavevector-magnitude axis to sweep.');
            end
            validateattributes(kParRange, {'double'}, {'vector', 'real', 'finite', 'positive'}, 'CalcIsofreqSingleLayered', 'kParRange');

            obj.angle     = angle;
            obj.structure = CalcCommon.rotateStructure(structure, angle);
            obj.omega0    = omega0;
            obj.kParRange = kParRange(:);
        end

        function obj = calc(obj)
        %%% Run the isoFreq Green's-function scan over kParRange
        %%% >> obj = calc(obj)
        %%%
        %%% Converts kParRange into the velRange argument GreenSolver.calc's
        %%% isoFreq scan mode expects: |k∥| = omega0/v, so v = omega0/|k∥|.
        %%% GreenSolver.calc returns G in the same element order as the velRange
        %%% it was given, so absG comes back already indexed by kParRange's order.

            velRangeForGreen = obj.omega0 ./ obj.kParRange;

            greenResult = GreenSolver.calc(obj.structure, 0, velRangeForGreen, 3, 3, ...
                'ScanMode', 'isoFreq', 'Omega0', obj.omega0, 'ComputeMinus', false);
                %%% 0, not obj.angle -- obj.structure is already rotated
                %%% (docs/reference/architecture/forwardCalcs.md, "Rotation convention: pass 0 to
                %%% solvers after construction").
                %%% forceDir=displDir=3 (G33) with computeMinus=false: a diagonal
                %%% pair has G_- = G_+ (reciprocity), so G_eff = G_+ adds nothing.

            obj.absG = abs(squeeze(greenResult.GFull(3, 3, :)));
        end
    end
end
