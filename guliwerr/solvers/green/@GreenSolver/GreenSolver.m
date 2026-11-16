classdef GreenSolver
%%% Surface dynamic Green's function G_IJ(k∥, ω) for anisotropic solids.
%%% Half-space and layered geometry (Every et al. 1997, Every 2002, Maznev et al. 1999).
%%%
%%% Value class — copy-on-modify semantics.
%%%
%%% Usage:
%%%   obj = GreenSolver.calc(structure, angle, velRange)
%%%   obj = GreenSolver.calc(structure, angle, velRange, 3, 3)   % G_33 (surface Brillouin scattering)
%%%   obj = GreenSolver.calc(structure, angle, velRange, 1, 3)   % G_13 (TGS)
%%%   obj = GreenSolver.calc(structure, angle, velRange, 3, 3, 'ComputeMinus', false)
%%%       % leave G_minus/G_eff (and their GFull_ forms) empty. Costs the same
%%%       % either way -- G_− is read off G_+ by reciprocity, not solved.
%%%       % obj.GFull(I,J,:) gives every component's G_+.
%%%   obj = GreenSolver.calc(structure, 0, velAxis, 1, 3, 'Eta', 3e-3)
%%%       % larger limiting-absorption damping, emulating a finite measurement
%%%       % linewidth (default 3e-4). CalcCommon.greenScoreSpectrum uses this call
%%%       % (with eta = SawSpectrum.ETA) and reads obj.G_eff (G13_eff).
%%%   obj = GreenSolver.calc(structure, 0, velAxis, 3, 3, 'ScanMode', 'isoFreq', 'Omega0', w0)
%%%       % fixed frequency, swept |k∥| (layered dispersion)
%%%
%%% Output:
%%%   obj.G       — G_IJ(+k∥)  (complex), for the requested (forceDir,displDir)
%%%   obj.G_minus — G_IJ(−k∥)  (complex) = G_JI(+k∥), by elastodynamic
%%%                 reciprocity (exact for any elastic tensor, cut and layer
%%%                 stack); empty if computeMinus=false
%%%   obj.G_eff   — (G_+ + (-1)^n * G_−) / 2, n = number of (I,J) equal to the
%%%                 propagation index (Maznev 1999, parity-signed), i.e.
%%%                 (G_IJ + (-1)^n * G_JI) / 2 at +k∥ -- e.g. G13_eff =
%%%                 (G_13 - G_31)/2; empty if computeMinus=false. Diagonal
%%%                 pairs (I=J, e.g. G_33): G_eff = G_+ exactly. Off-diagonal
%%%                 pairs (e.g. G_13, G_23): G_eff = G_+ only where k1->-k1 is a
%%%                 true symmetry of that direction. A two-fold axis along the
%%%                 surface normal (every cubic (001) and (110) face, isotropic,
%%%                 hexagonal with c along the normal) makes this true at ALL
%%%                 in-plane angles: the rotation R = diag(-1,-1,1) maps k to -k
%%%                 and gives G_IJ(-k) = s_I s_J G_IJ(k), s = (-1,-1,+1), so with
%%%                 reciprocity G_13 = -G_31 and G13_eff = G_13 exactly, while
%%%                 G_12 = G_21 and G_23 = -G_32 give G12_eff = G23_eff = 0
%%%                 identically (G_12, G_23 themselves are not zero). See
%%%                 docs/reference/architecture/green.md, "±k symmetrisation".
%%%   obj.GFull, obj.GFull_minus, obj.GFull_eff — the SAME three quantities for
%%%                 ALL nine (forceDir,displDir) pairs at once, [3×3×nVel];
%%%                 GFull_minus = GFull transposed in (I,J). GFull_minus/
%%%                 GFull_eff follow the same computeMinus=false emptiness rule
%%%                 as G_minus/G_eff.
%%%   obj.computeMinus — whether the G_−/G_eff fields were filled
%%%
%%% Optional full output ([obj, ~] = GreenSolver.calc(...), intended for a SCALAR
%%% velocity): asking for a second output switches it on, even when it is
%%% discarded with ~ (that is the intended call), and it is large (the half-space
%%% dispField is 6003 x 3 x 3 x nVel complex). Half-space -> obj.dispField/
%%% .depthAxis; layered (structure.nLay>1) -> obj.k3OutLayered/
%%% .polarizationULayered/.partialAmplitudeLayered/.hScaleOut instead, with no
%%% dispField grid (Vis.calcGreenDispField2DLayered reconstructs the depth
%%% profile from these plus the caller's Structure).
%%%
%%% Coordinate system:
%%%   Standard convention throughout GULIWERR (Material.GetStandardDirMat, applied
%%%   inside calc.m): I/J = 1 → propagation direction, 3 → surface normal.
%%%   Default I=1, J=3 gives G_13 — the TGS component (horizontal force, vertical disp).
%%%
%%% GuliweRR v6.0

    %% Properties
    properties (SetAccess = protected)
        forceDir(1,1)   double = 1      % I: force direction (column of inv(B))
        displDir(1,1)   double = 3      % J: displacement direction (row of U)

        velRange(1,:)   double = []     % phase-velocity axis used in the calculation
        G(1,:)          double = []     % G_IJ(+k∥) — complex
        G_minus(1,:)    double = []     % G_IJ(−k∥) = G_JI(+k∥) (reciprocity) — complex; empty if computeMinus=false
        G_eff(1,:)      double = []     % (G_+ + (-1)^n * G_−) / 2, parity-signed (see header); empty if computeMinus=false

        computeMinus(1,1) logical = true  % whether G_minus/G_eff/GFull_minus/GFull_eff were filled (calc()'s 'ComputeMinus')

        %%% Full (forceDir,displDir) tensor, [3×3×nVel] — all nine G_IJ pairs from
        %%% the same solve. GFull always populated; GFull_minus/GFull_eff only if
        %%% computeMinus=true.
        GFull                    % G_IJ(+k∥) for all (I,J) — complex [3×3×nVel]
        GFull_minus              % G_IJ(−k∥) = GFull transposed in (I,J) — complex [3×3×nVel]; empty if computeMinus=false
        GFull_eff                % parity-signed ±k mean per (I,J), as G_eff — complex [3×3×nVel]; empty if computeMinus=false

        %%% Optional full output — populated only when nargout > 1 in calc (half-space only)
        dispField                       % displacement field profile u_ij(depth, component, force, velocity)
        depthAxis                       % depth axis for dispField, normalised units (2*pi = one
                                         % wavelength) — [1 × nxx3] row vector, 0 to 30 wavelengths.
                                         % Use this rather than re-deriving the axis in the caller.

        %%% Optional full output for a LAYERED structure — populated only when
        %%% nargout > 1 in calc AND structure.nLay > 1 (the half-space fields
        %%% above stay empty then). {nLay×1} cells, STRUCTURE order (index 1 =
        %%% substrate ... index nLay = outermost film), ragged: each cell entry
        %%% is [6 × nVel] (film) or [3 × nVel] (substrate) for k3OutLayered,
        %%% [3 × 6/3 × nVel] for polarizationULayered, [6/3 × 3 × nVel] for
        %%% partialAmplitudeLayered (rows = partial-wave index, cols = forceDir
        %%% 1:3).
        %%% No dispField grid: each layer's depth range needs Structure's
        %%% relDepth, which this class does not keep, so the depth profile is
        %%% reconstructed downstream by Vis.calcGreenDispField2DLayered.
        %%% See docs/reference/architecture/green.md, "Layered extension".
        k3OutLayered                    % {nLay×1} cell, per-layer partial-wave k3
        polarizationULayered            % {nLay×1} cell, per-layer polarisations
        partialAmplitudeLayered         % {nLay×1} cell, per-layer amplitudes A

        %%% hScale used to build k3OutLayered/... above (1 for 'isoWavevec';
        %%% omega0/v for 'isoFreq'). Scalar, not per-velocity: like the
        %%% half-space full output, it is intended for a single velocity.
        hScaleOut                       % scalar
    end

    %% Constructor
    methods
        function obj = GreenSolver(forceDir, displDir)
            %%% Initialise GreenSolver; optionally set force and displacement directions.
            %%% >> obj = GreenSolver(forceDir, displDir)
            %%%
            %%% INPUT (opt): 1. forceDir .. force direction I (default = 1, TGS)
            %%%              2. displDir .. displacement direction J (default = 3, TGS)
            if nargin >= 1 && ~isempty(forceDir)
                validateattributes(forceDir, {'numeric'}, {'scalar','integer','>=',1,'<=',3}, ...
                    'GreenSolver', 'forceDir');
                obj.forceDir = forceDir;
            end
            if nargin >= 2 && ~isempty(displDir)
                validateattributes(displDir, {'numeric'}, {'scalar','integer','>=',1,'<=',3}, ...
                    'GreenSolver', 'displDir');
                obj.displDir = displDir;
            end
        end
    end

    %% Static methods
    methods (Static)
        [obj, dispFieldOut] = calc(structure, angle, velRange, forceDir, displDir, varargin)
    end

end
