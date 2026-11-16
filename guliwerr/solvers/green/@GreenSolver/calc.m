function [obj, dispFieldOut] = calc(structure, angle, velRange, forceDir, displDir, varargin)
%%% Calculate the G_IJ surface Green's function for a Structure object.
%%% >> obj = GreenSolver.calc(structure, angle, velRange, forceDir, displDir)
%%% >> obj = GreenSolver.calc(..., 'ScanMode', 'isoFreq', 'Omega0', w0, 'ComputeMinus', false, 'Eta', 3e-3)
%%%
%%% INPUT:  1. structure   .. Structure object; mater(1)=substrate, mater(nLay)=surface
%%%         2. angle       .. propagation angle [degrees] — rotates dirMat in-plane
%%%  (opt)  3. velRange    .. phase-velocity scan range [km/s]; default linspace(0.1,10,5000)
%%%  (opt)  4. forceDir    .. force direction I (default = 1 = propagation, for TGS)
%%%  (opt)  5. displDir    .. displacement direction J (default = 3 = surface normal, for TGS)
%%%  (opt)  Name-Value ([] = default):
%%%         'ScanMode'     .. 'isoWavevec' — fix |k∥|=1, sweep ω=v (default; matches TGS
%%%                           at fixed q); 'isoFreq' — fix ω=Omega0, sweep |k∥|=Omega0/v
%%%                           (layered dispersion)
%%%         'Omega0'       .. reference velocity for 'isoFreq' [km/s]; default = max(velRange)
%%%         'ComputeMinus' .. logical: also fill G_minus/G_eff/GFull_minus/GFull_eff
%%%                           (default true). Free either way: G_− is read off the one
%%%                           +k∥ solve by reciprocity (see ±k PARITY below), no second
%%%                           solve. They stay empty when false.
%%%         'Eta'          .. limiting-absorption damping (default = 3e-4, see DAMPING
%%%                           below). Larger values emulate a finite measurement
%%%                           linewidth: omega -> omega*sqrt(1+i*eta), so every
%%%                           mode's field decays as exp(-omega*eta/2*t).
%%%
%%% OUTPUT: GreenSolver object with .G, .G_minus, .G_eff, .velRange, plus .GFull/
%%%         .GFull_minus/.GFull_eff — the same quantities for all nine
%%%         (forceDir,displDir) pairs at once, [3×3×nVel].
%%%         [obj, dispFieldOut] = GreenSolver.calc(...) — any second output,
%%%         also ~ — switches on the depth-resolved full output (intended for a
%%%         scalar velocity; large): half-space (nLay == 1) fills obj.dispField/
%%%         .depthAxis, and dispFieldOut = obj.dispField. For nLay > 1 it
%%%         instead fills obj.k3OutLayered/.polarizationULayered/
%%%         .partialAmplitudeLayered/.hScaleOut (per-layer, Structure order —
%%%         see GreenSolver.m and docs/reference/architecture/green.md, "Layered
%%%         extension") and dispFieldOut = [].
%%%
%%% COORDINATE SYSTEM: x1=propagation, x2=in-plane, x3=surface normal/depth
%%%   (Every/Maznev convention). I=1, J=3 gives G_13 (TGS).
%%%
%%% ±k PARITY: G_eff = (G_plus + (-1)^n * G_minus) / 2 (Maznev 1999, parity-signed).
%%%   Let n be how many of (forceDir,displDir) equal 1, the propagation index.
%%%   Elastodynamic reciprocity (C_ijkl = C_klij plus translation invariance
%%%   along the surface) gives G_IJ(−k∥) = G_JI(+k∥) EXACTLY, for any elastic
%%%   tensor, cut and layer stack, damped or not -- so GFull_minus is GFull
%%%   transposed in (I,J) and only the k1=+1 branch is ever solved. Hence
%%%   G_eff = (G_IJ + (-1)^n * G_JI) / 2 at +k∥, e.g. G13_eff = (G_13 - G_31)/2.
%%%   Diagonal pairs (I=J, e.g. G_33, G_11): G_plus = G_minus, G_eff = G_plus.
%%%   Off-diagonal pairs (I≠J, e.g. G_13, G_23): G_plus = (-1)^n * G_minus,
%%%   i.e. G_eff = G_plus, ONLY where k1->-k1 is a true symmetry of that
%%%   material/direction; elsewhere G_IJ ≠ ±G_JI and the two sides differ.
%%%   Derivation and numerical checks: docs/reference/architecture/green.md.
%%%
%%% DAMPING: the limiting-absorption complex frequency w = omega^2*(1+i*eta) is
%%%   carried through the WHOLE computation — every layer's k3/U/B/G, not a side
%%%   channel. This is the standard retarded-Green's-function prescription
%%%   omega -> omega+i0+ (Every et al. 1997, p. 1347): it broadens the pole and makes
%%%   substrate classification unambiguous (Im(k3)>0, sextic splits exactly 3/3).
%%%   The nargout>1 full output reports that same damped solve.
%%%
%%% NEAR-ISOTROPIC NOTICE: a layer with computed shear anisotropy ratio below
%%%   ~1.01 raises 'GreenSolver:calc:nearIsotropic' — informational only:
%%%   Calc_U_allVel then uses a closed-form SH/SV basis for the degenerate shear
%%%   branch instead of eig()'s discontinuous basis choice. Fires once per thread,
%%%   then self-silences; warning('on', ...) re-enables it. See
%%%   docs/reference/architecture/green.md, "Near-isotropic materials".
%%%
%%% CONDITIONING: B can be genuinely near-singular at specific (angle, velocity)
%%%   points, e.g. a thick film at high |k∥|. MATLAB's own per-point
%%%   'MATLAB:nearlySingularMatrix'/'MATLAB:singularMatrix'/
%%%   'MATLAB:illConditionedMatrix' warnings are suppressed here (a sweep can
%%%   trigger hundreds) and replaced by rcond(B) against a size-scaled
%%%   threshold (nWaves*eps; a non-finite B, RCOND = NaN, counts as singular); one consolidated
%%%   'GreenSolver:calc:nearSingular' warning reports the affected velocity range
%%%   and worst RCOND. See docs/reference/architecture/green.md, "Conditioning".
%%%
%%% Algorithm, ±k symmetrisation, relDepth's 2*pi convention and scan modes:
%%% docs/reference/architecture/green.md.
%%%
%%% References: Every et al. 1997 JASA 102; Every 2015 JASA 138;
%%%             Maznev et al. 1999 JAP 86; Every 2002 MST 13.

%%% ------- numerical constants -----------------------------------------------
OMEGA_TOLERANCE    = 5e-3;  %%% frequency match for degenerate k3 detection in Calc_U
ZERO_TOLERANCE     = 1e-4;  %%% k3 degeneracy gap threshold in CheckIdenticalRecords

%%% ------- limiting absorption (see DAMPING in the docblock) ----------------
%%% eta must sit far above roots()' accuracy on a DOUBLE root (~sqrt(eps) ~
%%% 1.5e-8) and far below the spacing between physically distinct roots (O(1)
%%% in these normalised units). 1e-5 to 1e-3 all satisfy that, so substrate
%%% classification alone does not pin the value down.
%%%
%%% What does: eta is threaded through the whole solve, so it also sets the
%%% SAW/Rayleigh peak's lineshape in G. A bound (non-leaky) Rayleigh peak is a
%%% true pole of the undamped problem, so at very small eta its raw peak HEIGHT
%%% does not converge under any finite velocity-axis sampling. Raising eta
%%% gives the pole a finite, sampling-independent width (Lorentzian-type
%%% broadening, peak height x FWHM conserved) at the cost of compressing the
%%% peak and broadening genuinely leaky (pseudo-SAW) modes: the damping adds
%%% eta*v/2 to a pole's own half-width, so a narrow leaky line can be dominated
%%% by it. 3e-4 is the chosen compromise: residual sampling noise on a bound
%%% peak of a few per cent (docs/reference/architecture/green.md).
LIMITING_ABSORPTION_ETA = 3e-4;

%%% ------- conditioning: suppress MATLAB's per-point singular-matrix warnings,
%%% replaced at the end of this call by one consolidated warning. onCleanup
%%% restores them so the suppression never outlives this call.
%%% 'MATLAB:illConditionedMatrix' is the RCOND = NaN variant (a non-finite B).
warnStateSing1 = warning('off', 'MATLAB:nearlySingularMatrix');
warnStateSing2 = warning('off', 'MATLAB:singularMatrix');
warnStateSing3 = warning('off', 'MATLAB:illConditionedMatrix');
cleanupSingWarn = onCleanup(@() [warning(warnStateSing1), warning(warnStateSing2), ...
                                 warning(warnStateSing3)]);
nearSingCount   = 0;      %%% velocity points flagged
nearSingVMin    = Inf;    %%% smallest velocity among flagged points
nearSingVMax    = -Inf;   %%% largest velocity among flagged points
nearSingRcondWorst = Inf; %%% worst (smallest) rcond(B) seen
laFallbackCount    = 0;   %%% times the limiting-absorption selection did not give
                          %%% exactly 3 waves and the Im(k3) tiebreak was used

%%% ------- defaults and validation -------------------------------------------
if nargin < 3 || isempty(velRange)
    velRange = linspace(0.1, 10, 5000);
end
if nargin < 4 || isempty(forceDir)
    forceDir = 1;
end
if nargin < 5 || isempty(displDir)
    displDir = 3;
end
%%% Name-Value options: a plain loop rather than inputParser, since a sweep calls
%%% calc() once per point
opts = struct('ScanMode', 'isoWavevec', 'Omega0', [], 'ComputeMinus', true, ...
    'Eta', LIMITING_ABSORPTION_ETA);
optNames = fieldnames(opts);
if mod(numel(varargin), 2) ~= 0
    error('GreenSolver:calc:badOption', 'Options must be Name-Value pairs.');
end
for iOpt = 1:2:numel(varargin)
    name = varargin{iOpt};
    hit  = [];
    if ischar(name) || (isstring(name) && isscalar(name))
        hit = find(strcmpi(char(name), optNames));
    end
    if isempty(hit)
        error('GreenSolver:calc:badOption', ...
            'Unknown option (expected ''ScanMode'', ''Omega0'', ''ComputeMinus'' or ''Eta'').');
    end
    if ~isempty(varargin{iOpt + 1})
        opts.(optNames{hit}) = varargin{iOpt + 1};
    end
end
scanMode     = opts.ScanMode;
omega0       = opts.Omega0;
computeMinus = opts.ComputeMinus;
eta          = opts.Eta;
if isempty(omega0)
    omega0 = max(velRange);
end

validateattributes(structure, {'Structure'}, {'scalar'}, 'GreenSolver.calc', 'structure');
validateattributes(angle, {'double'}, {'scalar', 'real', 'finite'}, 'GreenSolver.calc', 'angle');
validateattributes(velRange, {'double'}, {'vector', 'real', 'finite', 'positive'}, 'GreenSolver.calc', 'velRange');
validateattributes(forceDir, {'double'}, {'scalar', 'integer', '>=', 1, '<=', 3}, 'GreenSolver.calc', 'forceDir');
validateattributes(displDir, {'double'}, {'scalar', 'integer', '>=', 1, '<=', 3}, 'GreenSolver.calc', 'displDir');
if ~(isscalar(computeMinus) && (islogical(computeMinus) || ...
        (isnumeric(computeMinus) && (computeMinus == 0 || computeMinus == 1))))
    error('GreenSolver:calc:badComputeMinus', 'ComputeMinus must be true or false.');
end
computeMinus = logical(computeMinus);
validateattributes(eta, {'double'}, {'scalar', 'real', 'finite', 'positive'}, 'GreenSolver.calc', 'Eta');
validateattributes(omega0, {'double'}, {'scalar', 'real', 'finite', 'positive'}, 'GreenSolver.calc', 'Omega0');
if ~(ischar(scanMode) || (isstring(scanMode) && isscalar(scanMode))) || ...
        ~ismember(char(scanMode), {'isoWavevec', 'isoFreq'})
    error('GreenSolver:calc:badScanMode', ...
        'ScanMode must be ''isoWavevec'' or ''isoFreq''; got ''%s''', string(scanMode));
end
if structure.nLay < 1
    error('GreenSolver:calc:emptyStructure', 'Structure must contain at least one layer (substrate).');
end

%%% ------- layer setup -------------------------------------------------------
%%% Internal ordering: iLay=1 = surface (top), iLay=nLay = substrate.
%%% Structure ordering: mater(1) = substrate (bottom), mater(nLay) = surface (top).
%%% The internal reversal matches the Every 2002 convention (free surface first).
nLay = structure.nLay;

elTen   = cell(nLay, 1);   %%% rotated elastic tensors, internal order
density = zeros(nLay, 1);  %%% densities [g/cm³], internal order

for iLay = 1:nLay
    vIdx = nLay - iLay + 1;  %%% index into structure.mater (reversed)
    mat0 = structure.mater(vIdx);

    %%% Standard (Every/Maznev) coordinate convention: x1=propagation, x3=surface
    %%% normal/depth, rotated within the native surface plane to the propagation angle.
    elTen{iLay}   = mat0.GetStandardElastTensor(angle);
    density(iLay) = mat0.rho;
end

%%% ------- near-isotropic material flag/warning -----------------------------
%%% The check is on the material's actual computed anisotropy ratio, not on
%%% symmClassNo. See docs/reference/architecture/green.md ("Near-isotropic materials").
NEAR_ISO_ANGLE_STEP  = 15;    %%% coarse angular scan step [deg] for the screen
NEAR_ISO_WARN_THRESH = 1.01;  %%% (fastest shear/slowest shear)^2 threshold
isIsotropic = false(nLay, 1);  %%% per layer: use Calc_U_allVel's closed-form
                                %%% SH/SV degenerate branch instead of eig()

%%% Dedup via MATLAB's own warning on/off state: a sweep calls calc() once per
%%% point, so without it this notice fires at every point. The state is
%%% per-thread, so a parallel sweep warns at most once per pool worker, and any
%%% serial use warns exactly once. TRADE-OFF: after the first material trips it,
%%% a second, genuinely different near-isotropic material on the same thread
%%% stays silent too — acceptable for a notice about an internal numerical
%%% branch, not a correctness warning.
WARN_ID = 'GreenSolver:calc:nearIsotropic';

for iLay = 1:nLay
    %%% The ratio is rotation-invariant: ElastConsts stores it exactly for isotropic
    %%% and cubic constants (shearAnisotropyRatio), so those need no scan here. NaN
    %%% (lower symmetry): scan this layer's rotated tensor.
    anisoRatio = structure.mater(nLay - iLay + 1).elConsts.shearAnisotropyRatio;
    if ~isfinite(anisoRatio)
        anisoRatio = ElastConsts.CalcAnisotropyRatio(elTen{iLay}, NEAR_ISO_ANGLE_STEP);
    end
    isIsotropic(iLay) = anisoRatio < NEAR_ISO_WARN_THRESH;
    if isIsotropic(iLay) && strcmp(warning('query', WARN_ID).state, 'on')
        warning(WARN_ID, ...
            ['Layer %d/%d (internal order, surface-first) has shear anisotropy ' ...
             'ratio %.6f (< %.3f threshold) -- close enough to isotropic that ' ...
             'GreenSolver uses a closed-form SH/SV basis for the degenerate shear ' ...
             'branch instead of eig() (exact for isotropic, an approximation that ' ...
             'improves as the ratio -> 1 for near-isotropic materials). See ' ...
             'docs/reference/architecture/green.md, "Near-isotropic materials", for the full finding. ' ...
             '(Shown once per session -- further occurrences, including from OTHER near-' ...
             'isotropic materials, are suppressed; warning(''on'', ''%s'') to re-enable.)'], ...
            iLay, nLay, anisoRatio, NEAR_ISO_WARN_THRESH, WARN_ID);
        warning('off', WARN_ID);
    end
end

%%% Layer thicknesses in internal order (surface-first).
%%% relDepth = h/lambda is converted to the k_par=1 normalised units used
%%% here via h_norm = 2*pi*relDepth, matching GammaMatrix.rawByDepth's own
%%% convention. See docs/reference/architecture/green.md ("relDepth -> physical phase
%%% factor") for the derivation.
nFilm = nLay - 1;  %%% number of film layers (0 for half-space)
if nFilm > 0
    h = zeros(nFilm, 1);
    for iLay = 1:nFilm
        vIdx    = nLay - iLay + 1;  %%% internal film iLay corresponds to Structure layer vIdx
        h(iLay) = 2*pi * structure.calcParams(vIdx).relDepth;
    end
    %%% cumDepthFromBelow(iLay) = sum of film thicknesses from film iLay to the bottom film.
    %%% This is the depth of the top of film iLay below the substrate interface (x3=0).
    H_total          = sum(h);
    cumDepthFromBelow = H_total - [0; cumsum(h(1:end-1))];  %%% nFilm × 1
end

%%% ------- scan setup -------------------------------------------------------
velRange = velRange(:)';    %%% ensure row vector
nVel     = length(velRange);
k2       = 0;               %%% wave vector component 2 (homogeneous direction); always 0

wantFullOutput = nargout > 1;

%%% Only the k1=+1 branch is solved: the k1=-1 one follows from it by
%%% reciprocity (see ±k PARITY in the docblock and the output assembly below).
k1 = 1;   %%% normalised wave vector component (|k∥|=1, k2=0)

%%% pre-allocate outputs
G_plus      = complex(NaN(1, nVel));
%%% Full (forceDir,displDir) tensor -- cheap byproduct of the same solve
GFull_plus  = complex(NaN(3, 3, nVel));

if wantFullOutput && nFilm == 0
    %%% displacement field profile: depth × 3 components × 3 force dirs × nVel
    %%% depth axis defined below after xx3 is set
    nxx3      = 6003;
    dispField = complex(NaN(nxx3, 3, 3, nVel));
    xx3       = linspace(0, 2*pi*30, nxx3);  %%% depth axis in normalised units
end

if wantFullOutput && nFilm > 0
    %%% Per-layer partial-wave diagnostics (k3, U, A), one cell per layer in
    %%% INTERNAL order (iLay=1 = surface film ... iLay=nLay = substrate);
    %%% reversed into Structure order after the velocity loop. Ragged: 6
    %%% partial waves per film, 3 for the substrate. No dispField grid (see
    %%% GreenSolver.m; Vis.calcGreenDispField2DLayered reconstructs it).
    k3StoreLay = cell(nLay, 1);
    UStoreLay  = cell(nLay, 1);
    AStoreLay  = cell(nLay, 1);
    for iLay = 1:nLay
        if iLay <= nFilm
            nWavesLay = 6;
        else
            nWavesLay = 3;
        end
        k3StoreLay{iLay} = complex(NaN(nWavesLay, nVel));
        UStoreLay{iLay}  = complex(NaN(3, nWavesLay, nVel));
        AStoreLay{iLay}  = complex(NaN(nWavesLay, 3, nVel));
    end
    hScaleOutVal = NaN;
end

%%% ------- conditioning threshold ---------------------------------------
%%% B's size is fixed for the whole call (depends only on nFilm, not v).
%%% nWaves*eps is the standard "numerically singular" convention.
nWavesForRcond  = 3 + 6*nFilm;
rcondWarnThresh = nWavesForRcond * eps;

%%% ------- Christoffel Det polynomial, whole velocity axis -----------------
%%% The characteristic sextic of the partial-wave method (Farnell 1970; Every
%%% et al. 1997, Eq. 10), solved as a polynomial in k3 -- not as the 6x6
%%% eigenvalue problem of the Stroh (1958) / Barnett-Lothe (1973) formalism,
%%% which holds the same roots.
%%% DetW [7 x 4] holds the velocity-INDEPENDENT part of the sextic k3-finder,
%%% built once per layer -- see Calc_DetW_byChristoffel. Since
%%% w = omega^2 = v^2 is known for the whole axis up front, every point's
%%% degree-6 coefficient vector comes from one [7x4]*[4xnVel] product:
%%% DetAllD_lay{...}(:,iv) is ready to hand straight to roots().
%%% Cost: 7*nVel complex doubles per layer (~0.6 MB at nVel=5000).
%%%
%%% Only the DAMPED polynomial is built, for all layers. velRangeD is the
%%% complex damped velocity axis whose square is wAllD; it is passed to
%%% Calc_U_allVel so the polarisation solve's own rho*omega^2 term is damped
%%% consistently with the k3 root-find.
wAll      = velRange(:).' .^ 2;                             %%% 1 x nVel (real)
wAllD     = wAll * (1 + 1i*eta);         %%% 1 x nVel (complex)
WPowD     = [ones(1, nVel); wAllD; wAllD.^2; wAllD.^3];      %%% 4 x nVel
velRangeD = sqrt(wAllD);                                     %%% complex damped velocity axis
DetAllD_lay = cell(nLay, 1);
for iLay = 1:nLay
    DetW = Calc_DetW_byChristoffel(elTen{iLay}, density(iLay), k1, k2);
    DetAllD_lay{iLay} = DetW * WPowD;     %%% 7 x nVel, damped
end

%%% ------- per-layer traction/Gamma contraction matrices --------------------
%%% Both depend only on the layer's elastic tensor, so they are built once per
%%% layer here: C3 for calcTractionCols, elTenPerm for Calc_U_allVel's Gamma.
C3_lay        = cell(nLay, 1);
elTenPerm_lay = cell(nLay, 1);
for iLay = 1:nLay
    C3_lay{iLay}        = reshape(squeeze(elTen{iLay}(3,:,:,:)), [3,9]);  %%% rows=l, cols=(q-1)*3+p
    elTenPerm_lay{iLay} = reshape(permute(elTen{iLay}, [1,3,2,4]), 9, 9);
end

%%% ------- whole-axis k3 / polarisation pre-pass ----------------------------
%%% The root-find and the polarisation solve depend only on (layer, velocity)
%%% -- never on the boundary conditions, the layer stack geometry or
%%% the scan mode -- so they are computed here for the whole axis rather than
%%% inside the velocity loop. That lets both per-point solves run as batched
%%% page operations instead of one LAPACK call per point (per-call dispatch,
%%% not arithmetic, is what costs):
%%%   - the sextic roots: one pageeig over all companion matrices, built
%%%     exactly as roots.m builds its own (first row -c(2:end)/c(1), ones on the
%%%     subdiagonal). pageeig runs the same eig per page, so the roots AND their
%%%     order are bit-identical to a roots() loop (checked on half-space,
%%%     low-symmetry and layered fixtures); the degeneracy/partial-wave-selection
%%%     logic below depends on that exact k3 ordering. roots() itself only
%%%     strips leading/trailing zero coefficients before building the companion
%%%     matrix, so any point whose leading or trailing coefficient is zero (or
%%%     non-finite) is handed to roots() and keeps its exact behaviour,
%%%     including the rootCount error.
%%%   - the 3x3 null vectors: row cross products (see Calc_U_allVel).
%%% Memory stays bounded by chunking (see Calc_U_allVel); the companion pages
%%% are 6x6x nVel complex (~0.6 MB per 1000 points), built once per layer.
k3_pre  = cell(nLay, 1);   %%% each [6 x nVel], DAMPED roots
U_pre   = cell(nLay, 1);   %%% each [3 x 6 x nVel], DAMPED polarisations
adm_pre = [];              %%% [6 x nVel] logical: substrate admissibility
for iLay = 1:nLay
    D    = DetAllD_lay{iLay};
    k3M  = complex(zeros(6, nVel));
    %%% roots.m's own test for "companion matrix of the full degree-6 polynomial":
    %%% finite coefficients, non-zero leading and trailing coefficient, and finite
    %%% normalised coefficients (roots() drops a leading coefficient whose ratio
    %%% overflows). Anything else goes through roots() so its stripping/error
    %%% behaviour is untouched.
    dNorm   = D(2:7, :) ./ D(1, :);
    okBatch = all(isfinite(D), 1) & (D(1, :) ~= 0) & (D(7, :) ~= 0) & all(isfinite(dNorm), 1);
    idxB    = find(okBatch);
    if ~isempty(idxB)
        nB    = numel(idxB);
        compA = repmat(diag(ones(1, 5), -1), 1, 1, nB);   %%% ones on the subdiagonal
        %%% -(c/c1) exactly as roots.m writes it: (-c)/c1 gives the same values but
        %%% the opposite sign of zero in the imaginary parts, and eig's rounding
        %%% then differs in the last bits.
        compA(1, :, :) = reshape(-dNorm(:, idxB), 1, 6, nB);
        k3M(:, idxB) = 1i * reshape(pageeig(compA), 6, nB);   %%% sextic-equation convention: P = -i k3
    end
    for iv = reshape(find(~okBatch), 1, [])
        r = roots(D(:, iv));
        if numel(r) ~= 6
            %%% everything downstream assumes exactly 6 roots
            error('GreenSolver:calc:rootCount', ...
                ['Christoffel polynomial gave %d roots (expected 6) at v = %.4f km/s, ' ...
                 'layer %d -- leading/trailing coefficient vanished.'], numel(r), velRange(iv), iLay);
        end
        k3M(:, iv) = 1i * r;   %%% sextic-equation convention: P = -i k3
    end
    %%% Substrate only: admissibility is read directly off Im(k3)>0 on the
    %%% same roots already computed for the physics. With Im(omega) > 0 no root
    %%% is real, so the damped sextic splits exactly 3/3 (docs/reference/architecture/green.md,
    %%% step 3): no threshold and no group-velocity evaluation are needed.
    if iLay == nLay
        adm_pre = imag(k3M) > 0;
    end

    %%% degenerate k3 pairs, whole axis at once
    identM = CheckIdenticalRecords_allVel(k3M, ZERO_TOLERANCE);
    k3_pre{iLay} = k3M;
    U_pre{iLay}  = Calc_U_allVel(elTenPerm_lay{iLay}, density(iLay), ...
        k1, k2, k3M, velRangeD, identM, OMEGA_TOLERANCE, ZERO_TOLERANCE, isIsotropic(iLay));
end

%%% ------- velocity loop ---------------------------------------------------
for iv = 1:nVel
    v     = velRange(iv);
    omegaD = velRangeD(iv);   %%% normalised |k∥|=1, so omega = v; damped: v*sqrt(1+i*eta),
                              %%% matching the damped k3/U

    %%% scan-mode scale for effective layer depth (changes at each v in 'isoFreq' mode)
    if strcmp(scanMode, 'isoFreq')
        hScale = omega0 / v;  %%% |k∥| = ω0/v; effective h_norm = h × |k∥|
    else
        hScale = 1;           %%% isoWavevec: |k∥|=1, h_norm = h (constant)
    end

    pointRcondWorst = Inf;  %%% worst rcond(B) at this v

    %%% per-layer k3 and polarisation vectors
    k3_lay = cell(nLay, 1);
    U_lay  = cell(nLay, 1);

    for iLay = 1:nLay
        %%% read off the whole-axis pre-pass above
        k3_all6 = k3_pre{iLay}(:, iv);
        U_all6  = U_pre{iLay}(:, :, iv);

        if iLay == nLay
            %%% substrate: keep only the 3 physically admissible partial
            %%% waves (decaying into depth). Im(k3)>0 on the damped roots
            %%% gives exactly 3 by construction (the 3/3 split of the
            %%% decaying/outgoing choice, Every et al. 1997, pp. 1346-1347). The
            %%% fallback stays inside the same criterion rather than switching
            %%% methods: take the 3 waves with the largest Im(k3), i.e. the most
            %%% strongly decaying.
            admCol = adm_pre(:, iv);
            if sum(admCol) == 3
                inds = find(admCol);
            else
                laFallbackCount = laFallbackCount + 1;
                [~, ord] = sort(imag(k3_all6), 'descend');
                inds = ord(1:3);
            end
            indsS        = sort(inds);
            k3_lay{iLay} = k3_all6(indsS);
            U_lay{iLay}  = U_all6(:, indsS);
        else
            %%% film: retain all 6 partial waves
            k3_lay{iLay} = k3_all6;
            U_lay{iLay}  = U_all6;
        end
    end

    %%% build boundary condition matrix and compute G_IJ(v)
    if nFilm == 0
        %%% half-space: 3×3 boundary condition matrix (Every 2015 Eq. 5 and 6)
        B    = calcTractionCols(C3_lay{1}, U_lay{1}, k1, k2, k3_lay{1}, omegaD);
        rc = rcond(B);
        if isnan(rc), rc = 0; end  %%% non-finite B: counted as singular (CONDITIONING)
        pointRcondWorst = min(pointRcondWorst, rc);  %%% see CONDITIONING, docblock
        invB = inv(B);  %%% full inverse needed: all elements of invB are used in G_IJ sum

        %%% Full 3x3 G_IJ tensor for this v/k-sign; the requested
        %%% (forceDir,displDir) is one entry of it. omegaD's imaginary part
        %%% is what regularises G at the pole.
        GvalFull = calcGFullFromPartialWaves(invB, U_lay{1}, 1, omegaD);
        Gval     = GvalFull(forceDir, displDir);

        %%% optional full output (k1=+1, half-space only) -- the depth profile
        %%% of the damped solve, i.e. the same computation .G comes from
        if wantFullOutput
            A     = (1i/omegaD) * invB;  %#ok<MINV>  %%% scaled full inverse, every element used
            u_ij  = complex(zeros(nxx3, 3, 3));
            for ci = 1:3
                for cj = 1:3
                    for imn = 1:3
                        u_ij(:,ci,cj) = u_ij(:,ci,cj) + ...
                            A(imn,cj) * U_lay{1}(ci,imn) .* exp(1i * k3_lay{1}(imn) .* xx3');
                    end
                end
            end
            dispField(:,:,:,iv) = u_ij;
        end

    else
        %%% layered: (6*nFilm+3)×(6*nFilm+3) boundary condition matrix (Every 2002).
        %%% phaseVec [6x1] is film 1's free-surface phase factor, built as part
        %%% of B's own free-surface block and reused for G below.
        [B, phaseVec] = buildBCMatrix(C3_lay, U_lay, k3_lay, k1, k2, omegaD, nFilm, cumDepthFromBelow, hScale);
        %%% Only the first 3 columns of inv(B) are ever read (the force direction ranges
        %%% over 1:3: G_IJ and the partial-wave amplitudes), so B is solved against
        %%% eye(n,3) instead of forming the full inverse. Same algebra, rounding-level
        %%% differences, and one LU instead of an LU plus n triangular solves.
        E3 = eye(size(B, 1), 3);
        %%% EQUILIBRATION: a thick film (relDepth >> 1) makes B's column norms
        %%% span an enormous range -- the exp(-Im(k3)*h) factors for
        %%% h = O(2*pi*relDepth) put some columns near double precision's
        %%% overflow end and others near its underflow end (ratio up to ~1e54
        %%% at relDepth=20). rcond(B) reports that as catastrophically
        %%% singular, but it is a SCALING problem, not a rank deficiency.
        %%% equilibrate() is exact algebra -- invB = C*inv(R*P*B*C)*R*P is the
        %%% identical true inverse -- and fixes the reported conditioning.
        %%%
        %%% Gated on the rcond(B) check already needed for CONDITIONING
        %%% tracking: equilibrate() does real per-call work, and running it on
        %%% every well-conditioned point measures ~18% slower for no benefit.
        rc = rcond(B);
        if isnan(rc), rc = 0; end  %%% non-finite B: counted as singular (CONDITIONING)
        if rc < rcondWarnThresh && all(isfinite(B(:)))
            %%% At extreme relDepth (~O(100), already discouraged by the
            %%% legO >= 2*relDepth rule in docs/pitfalls.md) exp(Im(k3)*h) can
            %%% overflow to Inf while B is being formed -- the data is lost at
            %%% source and no post-hoc scaling recovers it. equilibrate()
            %%% errors on non-finite input, so the isfinite guard falls back
            %%% to the raw rcond/inv(B) path instead of crashing.
            [P, R, C] = equilibrate(B);
            Beq   = R * P * B * C;
            rcEq  = rcond(Beq);
            if rcEq > rc
                rc   = rcEq;
                invB = C * (Beq \ (R * P * E3));   %%% columns 1:3 of C*inv(Beq)*R*P
            else
                invB = B \ E3;   %%% equilibration didn't help here; fall back
            end
        else
            invB = B \ E3;
        end
        pointRcondWorst = min(pointRcondWorst, rc);  %%% see CONDITIONING, docblock

        %%% optional full output (k1=+1, layered) -- every layer's
        %%% partial-wave amplitudes are already in this SAME invB, not just
        %%% film 1's 6 rows. Stored in INTERNAL (surface-first) order;
        %%% reversed into Structure order after the velocity loop.
        if wantFullOutput
            hScaleOutVal = hScale;
            for iLay = 1:nLay
                if iLay <= nFilm
                    rowsLay = 6*(iLay-1)+1 : 6*iLay;
                else
                    rowsLay = 6*nFilm+1 : 6*nFilm+3;
                end
                %%% Rows=partial-wave index, cols=forceDir 1:3.
                AStoreLay{iLay}(:,:,iv)  = (1i/omegaD) * invB(rowsLay, 1:3);
                k3StoreLay{iLay}(:,iv)   = k3_lay{iLay};
                UStoreLay{iLay}(:,:,iv)  = U_lay{iLay};
            end
        end

        %%% G_IJ from layer-1 waves (columns 1..6 of B → rows 1..6 of invB),
        %%% with phaseVec's free-surface factor at x3 = −H_total
        GvalFull = calcGFullFromPartialWaves(invB, U_lay{1}, phaseVec, omegaD);
        Gval     = GvalFull(forceDir, displDir);
    end

    G_plus(iv)         = Gval;
    GFull_plus(:,:,iv) = GvalFull;


    if pointRcondWorst < rcondWarnThresh
        nearSingCount      = nearSingCount + 1;
        nearSingVMin       = min(nearSingVMin, v);
        nearSingVMax       = max(nearSingVMax, v);
        nearSingRcondWorst = min(nearSingRcondWorst, pointRcondWorst);
    end
end  %%% velocity loop

if laFallbackCount > 0
    warning('GreenSolver:calc:absorptionFallback', ...
        ['%d/%d substrate selections did not yield exactly 3 admissible partial ' ...
         'waves under the limiting-absorption test and fell back to a largest-' ...
         'Im(k3) tiebreak.'], laFallbackCount, nVel);
end

if nearSingCount > 0
    %%% `angle` is deliberately left out of the message: several callers
    %%% (e.g. CalcIsofreqSingleLayered.calc) pass angle=0 here having already
    %%% rotated the structure's Materials upstream, so it is not reliably the
    %%% physical propagation direction. One warning per calc() call still gives
    %%% a per-angle report; the caller's own loop identifies which angle.
    warning('GreenSolver:calc:nearSingular', ...
        ['%d/%d velocity points (v = %.4g-%.4g km/s) had a near-singular ' ...
         'boundary-condition matrix (worst RCOND = %.3e < threshold %.3e) ' ...
         '-- G values there may be unreliable. See ' ...
         'docs/reference/architecture/green.md, "Conditioning".'], ...
        nearSingCount, nVel, nearSingVMin, nearSingVMax, ...
        nearSingRcondWorst, rcondWarnThresh);
end

%%% ------- assemble output object -------------------------------------------
obj              = GreenSolver(forceDir, displDir);
obj.velRange     = velRange;
obj.computeMinus = computeMinus;
obj.G            = G_plus;
obj.GFull        = GFull_plus;
if computeMinus
    %%% G_IJ(−k∥) = G_JI(+k∥) by reciprocity (see ±k PARITY in the docblock):
    %%% the k1=-1 tensor is the k1=+1 one transposed in (I,J), no second solve.
    %%% It agrees with an independent k1=-1 solve to ~1e-12 relative on half-space
    %%% and layered stacks, both scan modes.
    GFull_minus     = permute(GFull_plus, [2 1 3]);
    obj.GFull_minus = GFull_minus;
    obj.G_minus     = reshape(GFull_minus(forceDir, displDir, :), 1, []);
    %%% ±k symmetrisation (Maznev 1999), parity-signed: G_eff = (G_+ + (-1)^n * G_−) / 2,
    %%% n = how many of (forceDir,displDir) equal 1. reflSign(I)*reflSign(J) = (-1)^n.
    reflSign        = @(k) 1 - 2*(k == 1);   %%% +1 for k=2,3 ; -1 for k=1
    obj.G_eff       = (G_plus + reflSign(forceDir)*reflSign(displDir) * obj.G_minus) / 2;
    signMat         = reflSign([1 2 3])' * reflSign([1 2 3]);   %%% [3x3], S(I,J) = (-1)^n
    obj.GFull_eff   = (GFull_plus + signMat .* GFull_minus) / 2;
end
%%% computeMinus=false: G_minus/G_eff/GFull_minus/GFull_eff stay at their
%%% class-default empty [], the convention for "not computed".

%%% ------- optional full output ----------------------------------------------
if wantFullOutput
    if nFilm == 0
        obj.dispField = dispField;
        obj.depthAxis = xx3;
        dispFieldOut  = dispField;
    else
        %%% Reverse INTERNAL (surface-first) order into STRUCTURE order
        %%% (index 1 = substrate ... index nLay = outermost film), the
        %%% project-wide layer convention.
        obj.k3OutLayered            = k3StoreLay(end:-1:1);
        obj.polarizationULayered    = UStoreLay(end:-1:1);
        obj.partialAmplitudeLayered = AStoreLay(end:-1:1);
        obj.hScaleOut               = hScaleOutVal;
        %%% dispField/depthAxis and the direct output stay half-space only.
        dispFieldOut = [];
    end
end

end  %%% end of calc


%% ============================================================
%%  LOCAL FUNCTIONS
%% ============================================================

function DetW = Calc_DetW_byChristoffel(elTen, rho, k1, k2)
%%% Velocity-INDEPENDENT part of the sextic k3-finder.
%%% >> DetW = Calc_DetW_byChristoffel(elTen, rho, k1, k2)
%%%
%%% The k3-finder needs the determinant of a 3x3 matrix of quadratics in
%%% P (= -i k3)
%%%
%%%     d(ll,n)(P) = c_ln*P^2 + b_ln*P + (a_ln + rho*delta_ln*w),   w = omega^2
%%%
%%% a_ln, b_ln, c_ln contain no omega, so the ONLY velocity dependence is the
%%% additive rho*delta_ln*w on the three DIAGONAL entries. A determinant is
%%% multilinear in its rows and each row contributes at most one factor of w, so
%%%
%%%     Det(P; w) = D0(P) + D1(P)*w + D2(P)*w^2 + D3(P)*w^3
%%%
%%% exactly — a cubic in w with polynomial-in-P coefficients, which is what
%%% lets the whole velocity axis be evaluated in one matrix product.
%%%
%%% OUTPUT: DetW [7 x 4] — row i = coefficient of P^(7-i) (descending, roots()'
%%% convention), column j = coefficient of w^(j-1). The caller recovers the
%%% per-velocity polynomial as Det = DetW * [1; w; w^2; w^3].
%%%
%%% D0..D3 come out exactly (no fit, no Symbolic Math Toolbox) by writing each
%%% d entry as a BIVARIATE polynomial in (P, w) — a 3x2 coefficient matrix,
%%% rows = P^2/P^1/P^0, columns = w^0/w^1 — so each polynomial multiplication
%%% is a conv2(). Convolution is index-additive, so the mixed descending-P /
%%% ascending-w ordering is harmless as long as it stays consistent.
%%%
%%% NOTE: k1 enters a quadratically and b linearly, so DetW depends on the
%%% k-sign; calc() only ever builds k1=+1 (k1=-1 follows by reciprocity). k2 is
%%% architecturally 0 in this file but is kept general here.

d2 = cell(3, 3);   %%% d2{ll,n} = [3 x 2] bivariate coefficient matrix in (P, w)
for ll = 1:3
    for n = 1:3
        a_ln = -elTen(ll,1,n,1)*k1*k1 ...
               - (elTen(ll,1,n,2)+elTen(ll,2,n,1))*k1*k2 ...
               - elTen(ll,2,n,2)*k2*k2;
        b_ln = -1i*(elTen(ll,1,n,3)+elTen(ll,3,n,1))*k1 ...
               - 1i*(elTen(ll,2,n,3)+elTen(ll,3,n,2))*k2;
        c_ln =  elTen(ll,3,n,3);
        M = complex(zeros(3, 2));
        M(:,1) = [c_ln; b_ln; a_ln];       %%% w^0 column: the original quadratic
        M(3,2) = rho * (ll == n);           %%% w^1 column: only the P^0 diagonal term
        d2{ll,n} = M;
    end
end

tp = @(x, y, z) conv2(x, conv2(y, z));   % triple bivariate polynomial product

DetW = tp(d2{3,3}, d2{1,1}, d2{2,2}) ...
     + tp(d2{3,1}, d2{1,2}, d2{2,3}) ...
     + tp(d2{2,1}, d2{1,3}, d2{3,2}) ...
     - tp(d2{3,1}, d2{1,3}, d2{2,2}) ...
     - tp(d2{3,3}, d2{1,2}, d2{2,1}) ...
     - tp(d2{1,1}, d2{2,3}, d2{3,2});
end


function identM = CheckIdenticalRecords_allVel(k3M, zeroTolerance)
%%% Detect degenerate k3 pairs (two k3 values nearly equal) — occurs at high-symmetry
%%% directions — for the WHOLE velocity axis at once.
%%% >> identM = CheckIdenticalRecords_allVel(k3M, zeroTolerance)
%%%
%%% INPUT:  k3M    .. [6 x nVel] Christoffel roots, one column per velocity
%%% OUTPUT: identM .. [6 x nVel]; non-zero entries are INDICES into that
%%%         column's k3 that form a degenerate pair
%%%
%%% The sequential `counts` bookkeeping matters: the five gap slots are visited
%%% in order ii=1..5 and capped at `counts <= 4`, so a column records at most
%%% two pairs. sort() is stable in MATLAB, so a column-wise sort breaks ties the
%%% same way a per-column sort of a vector would.

nVel = size(k3M, 2);

%%% mean(abs(real(k3))) > tol, per column -- sum/6 avoids mean()'s dispatch
useReal = sum(abs(real(k3M)), 1) / 6 > zeroTolerance;      %%% 1 x nVel

[~, posR] = sort(real(k3M), 1);
[~, posI] = sort(imag(k3M), 1);
position  = posR;
position(:, ~useReal) = posI(:, ~useReal);                  %%% 6 x nVel

%%% k3(position) per column, then the consecutive gaps
linIdx = position + (0:nVel-1) * 6;
gap    = diff(k3M(linIdx), 1, 1);                           %%% 5 x nVel
small  = abs(gap) < zeroTolerance;                          %%% 5 x nVel

identM = zeros(6, nVel);
counts = ones(1, nVel);
for ii = 1:5
    act = small(ii, :) & counts <= 4;
    cc  = find(act);
    if isempty(cc), continue; end
    identM(counts(cc)   + (cc-1)*6) = position(ii,   cc);
    identM(counts(cc)+1 + (cc-1)*6) = position(ii+1, cc);
    counts(cc) = counts(cc) + 2;
end
end


function UAll = Calc_U_allVel(elTenPerm, rho, k1, k2, k3M, velRange, identM, omegaTolerance, zeroTolerance, isIsotropic)
%%% Polarisation vectors (null space of the Christoffel matrix) for all 6 k3
%%% values at EVERY velocity point of one (layer, k-sign), returned as
%%% [3 x 6 x nVel].
%%%
%%% >> UAll = Calc_U_allVel(elTenPerm, rho, k1, k2, k3M, velRange, identM, ...)
%%%
%%% INPUT: k3M    .. [6 x nVel] Christoffel roots, one column per velocity
%%%        identM .. [6 x nVel] CheckIdenticalRecords output, one column per
%%%                  velocity (non-zero entries are INDICES forming degenerate pairs)
%%%
%%% Branch logic:
%%%   non-degenerate                      -> largest row cross product of Gamma - rho*omega^2*I (null vector)
%%%   degenerate, isotropic/near-iso      -> closed-form SH/SV basis
%%%   degenerate, generically anisotropic -> eig() with explicit orthogonalisation
%%%
%%% Gamma(l,j) = sum_mn elTen(l,m,j,n)*kVec(m)*kVec(n), with m varying fastest
%%% in the K(m,n) layout. NOTE: no ' anywhere -- k3 is complex and the Hermitian
%%% transpose would silently conjugate it (docs/pitfalls.md NUM1).
%%%
%%% ORDERING: the generic-anisotropic degenerate branch READS its partner's
%%% already-assigned U(:,pair2), so it must run as a scalar loop over the
%%% degenerate indices in ASCENDING order. Pre-filling the non-degenerate
%%% columns via the batched null-vector step cannot disturb that: every index in
%%% identicalRecord is by construction a degenerate index, so pair2 always lies
%%% inside the degenerate set and its column stays zero until its own turn --
%%% which is exactly what the `any(U(:,pair2) ~= 0)` test below probes.
%%%
%%% MEMORY: the velocity axis is processed in chunks of CHUNK points, so peak
%%% working set is O(CHUNK) rather than O(nVel) — a few MB for any velRange.

nVel = size(k3M, 2);
UAll = complex(zeros(3, 6, nVel));

%%% Boolean form of identM: degenMask(ii,iv) == any(identM(:,iv) == ii)
degenMask = false(6, nVel);
nz = identM > 0;
if any(nz(:))
    [~, ccNz] = find(nz);
    degenMask(sub2ind([6, nVel], identM(nz), ccNz(:))) = true;
end

CHUNK = 4096;
for c0 = 1:CHUNK:nVel
    c1   = min(c0 + CHUNK - 1, nVel);
    idxV = c0:c1;
    nc   = numel(idxV);
    Mc   = 6 * nc;

    %%% Column-major linearisation throughout: partial-wave index varies
    %%% FASTEST, velocity slowest -- so page (jj-1)*6+ii is wave ii of the
    %%% jj-th velocity in this chunk.
    kVecs    = kVecCols(k1, k2, k3M(:, idxV));
    Kall     = outerCols(kVecs, kVecs);
    GammaAll = reshape(elTenPerm * Kall, 3, 3, Mc);

    %%% rho*omega^2*I, with omega = v repeated once per partial wave
    wRep     = reshape(repmat(rho * velRange(idxV).^2, 6, 1), 1, 1, Mc);
    DetPages = GammaAll - wRep .* eye(3);

    dm    = degenMask(:, idxV);
    idxND = find(~dm(:));
    Uc    = complex(zeros(3, Mc));
    if ~isempty(idxND)
        %%% non-degenerate: Gamma - rho*omega^2*I has rank 2 (a simple eigenvalue),
        %%% so the cross product (bilinear, NO conjugation) of any two of its rows
        %%% is orthogonal to all three rows -- it spans the null space. Take, per
        %%% page, the largest-norm one of the three row-pair cross products (this
        %%% avoids a nearly parallel pair) and normalise it. Same ray as the last
        %%% right singular vector of an SVD, several times cheaper; the overall
        %%% phase is arbitrary (as it is for the SVD) and cancels in every G.
        %%% Tie-break order matches the Julia port (GreenSextic.jl _null_vector).
        P   = DetPages(:, :, idxND);
        a11 = reshape(P(1,1,:), 1, []); a12 = reshape(P(1,2,:), 1, []); a13 = reshape(P(1,3,:), 1, []);
        a21 = reshape(P(2,1,:), 1, []); a22 = reshape(P(2,2,:), 1, []); a23 = reshape(P(2,3,:), 1, []);
        a31 = reshape(P(3,1,:), 1, []); a32 = reshape(P(3,2,:), 1, []); a33 = reshape(P(3,3,:), 1, []);
        %%% row1 x row2, row1 x row3, row2 x row3 (each 3 x nND)
        X1 = [a12.*a23 - a13.*a22; a13.*a21 - a11.*a23; a11.*a22 - a12.*a21];
        X2 = [a12.*a33 - a13.*a32; a13.*a31 - a11.*a33; a11.*a32 - a12.*a31];
        X3 = [a22.*a33 - a23.*a32; a23.*a31 - a21.*a33; a21.*a32 - a22.*a31];
        n1 = sum(abs(X1).^2, 1);  n2 = sum(abs(X2).^2, 1);  n3 = sum(abs(X3).^2, 1);
        use1 = n1 >= n2 & n1 >= n3;
        use2 = ~use1 & n2 >= n3;
        Xb   = X3;  nb = n3;
        Xb(:, use2) = X2(:, use2);  nb(use2) = n2(use2);
        Xb(:, use1) = X1(:, use1);  nb(use1) = n1(use1);
        Uc(:, idxND) = Xb ./ sqrt(nb);
        %%% Numerically rank-1 page (largest cross product < 1e-6 of the matrix's own
        %%% scale, i.e. sigma_2/sigma_1 ~< 1e-6): a double eigenvalue that
        %%% CheckIdenticalRecords did not pair -- e.g. both shear waves at k3 ~ 0 on a
        %%% symmetry axis with eta ~ 1e-12. The null space is then (nearly) 2-D, the
        %%% cross product is noise, and the two waves of the pair would get the SAME
        %%% vector (singular B). Fall back to the SVD's null vector for those pages;
        %%% default eta (3e-4) keeps sigma_2/sigma_1 of order eta even at such a point, so the
        %%% fallback is never taken in ordinary use.
        frob2 = reshape(sum(sum(abs(P).^2, 1), 2), 1, []);
        bad   = find(~(nb > 1e-12 * frob2.^2) | ~isfinite(nb));
        if ~isempty(bad)
            [~, ~, VVbad]     = pagesvd(P(:, :, bad));
            Uc(:, idxND(bad)) = reshape(VVbad(:, 3, :), 3, numel(bad));
        end
    end
    Uc = reshape(Uc, 3, 6, nc);

    %%% Degenerate points: per-point, ascending-index scalar logic. NOTE the
    %%% chunk-local velocity index is `jp`, NOT `jj` -- the generic-anisotropic
    %%% branch below has its own `for jj = 1:3` eigenvalue scan, which would
    %%% silently clobber an outer `jj`.
    for jp = reshape(find(any(dm, 1)), 1, [])
        iv              = idxV(jp);
        omega           = velRange(iv);
        identicalRecord = reshape(identM(:, iv), 1, 6);
        U               = Uc(:, :, jp);
        for ii = reshape(find(dm(:, jp)), 1, [])
            kVec  = [k1, k2, k3M(ii, iv)];
            Gamma = GammaAll(:, :, (jp-1)*6 + ii);

            if isIsotropic
                %%% Isotropic/near-isotropic degenerate shear branch: closed form.
                %%% For isotropic elasticity Gamma(k)*v = mu*(k.k)*v + (lambda+mu)*k*(k.v)
                %%% (bilinear, no conjugation), so ANY v with k.v=0 is an eigenvector of
                %%% the degenerate shear eigenvalue — where eig()'s basis choice for a
                %%% repeated eigenvalue is arbitrary and discontinuous. This file's
                %%% Christoffel solve always has k2=0 (the angle is folded into elTen's
                %%% rotation upstream, not into k1/k2), so the sagittal plane is always
                %%% x1-x3 and P-SV-SH collapses to two fixed, orthogonal, continuous-in-k
                %%% vectors:
                eSH = [0; 1; 0];             %%% pure horizontal shear: k.eSH = k1*0+k2*1+k3*0 = k2 = 0
                eSV = cross(kVec(:), eSH);   %%% sagittal-plane vector, orthogonal to both k and eSH
                pairPos = find(identicalRecord == ii, 1);
                if mod(pairPos, 2) == 1
                    U(:,ii) = eSH;
                else
                    U(:,ii) = eSV / norm(eSV);
                end
            else
                %%% degenerate, generically anisotropic material: eig() gives two orthogonal
                %%% eigenvectors for the repeated eigenvalue
                [VV, DD] = eig(Gamma);
                freq = sqrt(abs(diag(DD)) / rho);

                %%% find which eigenvalues are at omega. freq is a REAL magnitude
                %%% (sqrt(abs(.))), so it must be compared against real(omega), not
                %%% omega itself: omega is the complex, limiting-absorption-damped
                %%% frequency (see DAMPING in this file's docblock), and its
                %%% imaginary part -- v*eta/2, which exceeds omegaTolerance once
                %%% eta is O(1e-3) or larger at the top of a typical velRange --
                %%% would otherwise make a genuine match fail this test purely
                %%% because of damping, not because the eigenvalue is actually off:
                %%% at larger eta both eigenvalues of a real degenerate pair would be
                %%% missed, and a boundary-condition matrix built from two
                %%% uncoordinated fallback eigenvectors can be numerically singular
                %%% (NaN/Inf G).
                degenLocation = zeros(1, 3);
                count = 1;
                for jj = 1:3
                    if abs(freq(jj) - real(omega)) < omegaTolerance && count <= 2
                        degenLocation(count) = jj;
                        count = count + 1;
                    end
                end

                %%% partner index in identicalRecord -- needed both for the normal
                %%% orthogonal pick below and for the degenLocation-not-found
                %%% fallback, so it is resolved once, up front.
                pairPos = find(identicalRecord == ii, 1);
                if isempty(pairPos)
                    pair2 = 0;
                elseif mod(pairPos, 2) == 1
                    pair2 = identicalRecord(pairPos + 1);
                else
                    pair2 = identicalRecord(pairPos - 1);
                end
                partnerAssigned = pair2 > 0 && any(U(:,pair2) ~= 0);

                if degenLocation(1) == 0 || degenLocation(2) == 0
                    warning('GreenSolver:calc:degenEig', ...
                        'Cannot find correct degenerate eigenvalue for k3(%d); result may be inaccurate.', ii);
                    if partnerAssigned
                        %%% Still honour the orthogonal-to-partner requirement that
                        %%% would otherwise be skipped here: an uncoordinated pick
                        %%% (e.g. always VV(:,1)) can come out parallel to the
                        %%% partner's own column, which makes the boundary-condition
                        %%% matrix B built from both columns exactly singular (see
                        %%% note above) rather than merely "may be inaccurate".
                        [~, best] = min(abs(VV' * U(:,pair2)));
                        U(:,ii) = VV(:, best);
                    else
                        U(:,ii) = VV(:,1);
                    end
                    continue
                end

                if isempty(pairPos)
                    U(:,ii) = VV(:, degenLocation(1));
                    continue
                end

                %%% pick the eigenvector orthogonal to the already-assigned partner
                if partnerAssigned
                    if abs(U(:,pair2)' * VV(:,degenLocation(1))) < zeroTolerance
                        U(:,ii) = VV(:, degenLocation(1));
                    elseif abs(U(:,pair2)' * VV(:,degenLocation(2))) < zeroTolerance
                        U(:,ii) = VV(:, degenLocation(2));
                    else
                        %%% fallback: pick whichever is not equal to the partner
                        if ~isequal(U(:,pair2), VV(:,degenLocation(1)))
                            U(:,ii) = VV(:, degenLocation(1));
                        else
                            U(:,ii) = VV(:, degenLocation(2));
                        end
                    end
                else
                    %%% partner not yet assigned: assign first eigenvector
                    U(:,ii) = VV(:, degenLocation(1));
                end
            end
        end  %%% partial-wave (degenerate indices only)
        Uc(:, :, jp) = U;
    end  %%% degenerate velocity points in this chunk

    UAll(:, :, idxV) = Uc;
end  %%% velocity chunk
end


function GvalFull = calcGFullFromPartialWaves(invB, U, phase, omega)
%%% Build the full 3x3 G_IJ tensor (all nine force/displacement pairs) at one
%%% velocity point and one k-sign, from already-solved partial-wave data.
%%% >> GvalFull = calcGFullFromPartialWaves(invB, U, phase, omega)
%%%
%%% An O(nWaves*9) contraction, negligible next to the root-find and matrix
%%% inversion that produced invB and U, so all nine pairs are computed every
%%% time rather than only the requested one.
%%%
%%% INPUT:  1. invB    .. [nWaves x nWaves] inverse boundary-condition matrix
%%%                       (nWaves=3 half-space, 6 layered); only columns 1:3
%%%                       are ever read (forceDir/displDir range over 1:3)
%%%         2. U       .. [3 x nWaves] surface-layer (film 1) polarisation vectors
%%%         3. phase   .. scalar 1 (half-space, no phase factor) or [nWaves x 1]
%%%                       (layered: per-partial-wave interface phase factor)
%%%         4. omega   .. scalar angular frequency; COMPLEX (damped) in this file
%%%                       -- its imaginary part is what regularises G at the pole
%%% OUTPUT: GvalFull   .. [3 x 3], GvalFull(I,J) = G_IJ for this v/k-sign

nWaves  = size(U, 2);
invBSub = invB(1:nWaves, 1:3);                                            %%% [nWaves x 3], (n,I)
UT      = U.';                                                             %%% [nWaves x 3], (n,J)
psiFull = reshape(invBSub, [nWaves,3,1]) .* reshape(UT, [nWaves,1,3]);     %%% [nWaves x 3 x 3], (n,I,J)
if isscalar(phase)
    psiFull = psiFull * phase;                                            %%% half-space: no per-wave phase
else
    psiFull = psiFull .* reshape(phase(:), [nWaves,1,1]);                 %%% layered: per-partial-wave phase
end
GvalFull = squeeze((1i/omega) * sum(psiFull, 1));                          %%% [3 x 3]
end


function T = calcTractionCols(C3, U, k1, k2, k3, omega)
%%% Traction columns for ALL partial waves of one layer at once.
%%% >> T = calcTractionCols(C3, U, k1, k2, k3, omega)
%%%
%%% T(l,n) = Σ_{p,q} C_{3,l,p,q} * U_p(n) * (k_q(n)/ω)
%%%
%%% C3 = reshape(squeeze(elTen(3,:,:,:)),[3,9]) is passed in, precomputed once
%%% per layer in calc's setup — it has no velocity, k-sign or partial-wave
%%% dependence. Its column index is (q-1)*3+p, p varying fastest, which is
%%% exactly what outerCols produces.
S  = kVecCols(k1, k2, k3) / omega;   %%% 3 x nW slowness
KU = outerCols(U, S);                %%% 9 x nW, (q-1)*3+p
T  = C3 * KU;                        %%% 3 x nW
end


function K = kVecCols(k1, k2, k3)
%%% Wave-vector (or, divided by omega, slowness) columns for a set of waves.
%%% >> K = kVecCols(k1, k2, k3)
%%%
%%% INPUT:  k1, k2 .. scalar in-plane components, shared by every wave
%%%         k3     .. [nW x 1] or [1 x nW] depth components, one per wave
%%% OUTPUT: K      .. [3 x nW], column n = [k1; k2; k3(n)]
n = numel(k3);
K = [repmat(k1, 1, n); repmat(k2, 1, n); reshape(k3, 1, n)];
end


function P = outerCols(A, B)
%%% Column-wise 3x3 outer products, flattened column-major.
%%% >> P = outerCols(A, B)
%%%
%%% INPUT:  A, B .. [3 x nW]
%%% OUTPUT: P    .. [9 x nW], with P((q-1)*3+p, n) = A(p,n) * B(q,n)
%%%
%%% i.e. the column-major linearisation of A(:,n)*B(:,n).'. That p-fastest
%%% layout is what kron(x,y) produces and what the 9x9 contraction matrices in
%%% this file (C3, elTenPerm) all index by — getting it wrong fails silently.
%%%
%%% NOTE: no ' anywhere -- the inputs are complex (U, k3, slowness) and the
%%% Hermitian transpose would silently conjugate them (docs/pitfalls.md NUM1).
%%% Callers that DO want a conjugate pass conj() in explicitly.
n = size(A, 2);
P = reshape(reshape(A, [3,1,n]) .* reshape(B, [1,3,n]), 9, n);
end


function [B, surfPhase] = buildBCMatrix(C3_lay, U_lay, k3_lay, k1, k2, omega, nFilm, cumDepthFromBelow, hScale)
%%% Build (6*nFilm+3)×(6*nFilm+3) boundary condition matrix for N films on a substrate.
%%% Following Every 2002: substrate interface at x3=0, free surface at x3=−H_total.
%%%
%%% 2nd output surfPhase [6x1]: film 1's per-partial-wave phase factor at the
%%% free surface, exp(-i*k3*H_total*hScale). Built here for B's own free-surface
%%% block, and returned because it is also exactly the factor calc() needs to
%%% carry G_IJ from the substrate-interface reference plane out to the surface.
%%%
%%% Column assignment:
%%%   cols 1..6         — film 1 (surface) partial waves
%%%   cols 7..12        — film 2 partial waves
%%%   ...
%%%   cols 6*(nFilm-1)+1..6*nFilm — film nFilm (bottom film) partial waves
%%%   cols 6*nFilm+1..6*nFilm+3  — substrate partial waves
%%%
%%% Row assignment:
%%%   rows 1..3                         — traction-free at free surface (film 1, with phase)
%%%   rows 3+6*(i-1)+1 .. 3+6*i        — displacement + traction continuity at interface i
%%%     (interface i is between film i and film i+1, or between film nFilm and substrate)
%%%
%%% C3_lay:   cell{1..nFilm} = film C3 matrices, cell{nFilm+1} = substrate C3
%%% U_lay:    cell{1..nFilm} = film polarisations (3×6), cell{nFilm+1} = substrate (3×3)
%%% k3_lay:   cell{1..nFilm} = film k3 values (6×1), cell{nFilm+1} = substrate (3×1)
%%% cumDepthFromBelow(iLay) = total depth from substrate interface to top of film iLay.
%%% hScale:   depth scaling factor (=1 for 'isoWavevec'; =omega0/v for 'isoFreq').

nWaves = 6*nFilm + 3;
B      = complex(zeros(nWaves, nWaves));

%%% per-layer traction blocks, computed once (depth/phase enters only as a
%%% per-column scale factor applied below)
T_lay = cell(nFilm+1, 1);
for iLay = 1:nFilm+1
    T_lay{iLay} = calcTractionCols(C3_lay{iLay}, U_lay{iLay}, k1, k2, k3_lay{iLay}, omega);
end

%%% Block 1: traction-free at free surface (rows 1-3, film 1 columns 1-6)
hEff1     = cumDepthFromBelow(1) * hScale;              %%% effective total depth to free surface
surfPhase = exp(-1i * k3_lay{1} * hEff1);               %%% [6x1], one factor per partial wave
ph1       = reshape(surfPhase, 1, 6);                   %%% row [1x6]; .' semantics, no conjugation
B(1:3, 1:6) = T_lay{1} .* ph1;

%%% Interface blocks (rows 3+6*(iLay-1)+1 .. 3+6*iLay) for iLay = 1..nFilm
for iLay = 1:nFilm
    rD = 3 + 6*(iLay-1) + 1;  %%% first displacement-continuity row
    rT = rD + 3;                %%% first traction-continuity row

    %%% depth of the lower interface of film iLay
    if iLay < nFilm
        dIface = cumDepthFromBelow(iLay+1) * hScale;  %%% inter-film interface
    else
        dIface = 0;  %%% film nFilm / substrate interface is always at x3=0 (reference)
    end

    cS = 6*(iLay-1) + 1;  %%% column start for this film's partial waves

    %%% this film's partial waves at its lower boundary
    ph = reshape(exp(-1i * k3_lay{iLay} * dIface), 1, 6);
    B(rD:rD+2, cS:cS+5) =  U_lay{iLay} .* ph;
    B(rT:rT+2, cS:cS+5) =  T_lay{iLay} .* ph;

    if iLay < nFilm
        %%% next film's partial waves at the same interface (upper boundary of film iLay+1)
        cN  = 6*iLay + 1;
        phN = reshape(exp(-1i * k3_lay{iLay+1} * dIface), 1, 6);
        B(rD:rD+2, cN:cN+5) = -U_lay{iLay+1} .* phN;
        B(rT:rT+2, cN:cN+5) = -T_lay{iLay+1} .* phN;
    else
        %%% substrate partial waves at x3=0 (phase=1; substrate is cell{nFilm+1})
        cSub = 6*nFilm + 1;
        B(rD:rD+2, cSub:cSub+2) = -U_lay{nFilm+1};
        B(rT:rT+2, cSub:cSub+2) = -T_lay{nFilm+1};
    end
end
end
