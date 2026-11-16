classdef CalcCommon
%%% CalcCommon — shared internals for the forwardCalcs Calc* family.
%%%
%%% Mechanical steps only — nothing here decides WHICH combination of solves to
%%% run (each Calc* class's own job) and nothing here does wave physics.
%%% Stateless static class: never instantiate it. Not fronted by Calc.*, with one
%%% exception: saveCalc/loadCalc are user-facing (a finished result to/from a
%%% .mat file) and are launched as Calc.saveResult/Calc.loadResult.
%%%
%%% Static helpers rather than a shared base class. The half-space/layered/plate
%%% families genuinely differ in mode taxonomy, in which Christoffel reference is
%%% meaningful, and in whether a system-wide limiting velocity exists at all;
%%% static calls leave each class's physics decisions visible at its own call
%%% site, where inheritance would hide them in an override chain. See
%%% docs/reference/architecture/forwardCalcs.md, "Why half-space and layered are separate
%%% class families".
%%%
%%% Error-ID convention: helpers that raise errors on a caller's behalf take an
%%% `errPrefix` argument (e.g. 'CalcSingleHS:calcAll') and append their own final
%%% component, so the thrown identifier stays the caller's own
%%% (CalcSingleHS:calcAll:badMode) and greppable from the calling class.
%%%
%%% GuliweRR v6.0

    methods (Static)
        function structRot = rotateStructure(structure, angleDeg)
        %%% Rotate every layer's Material in-plane to a propagation angle
        %%% >> structRot = CalcCommon.rotateStructure(structure, angleDeg)
        %%%
        %%% INPUT:  1. structure .. Structure object, UNROTATED
        %%%         2. angleDeg  .. propagation angle [deg]
        %%% OUTPUT: structRot .. same Structure with every layer's dirMat rotated
        %%%
        %%% CalcParams is preserved per layer (ModifyLayer with calcParams0 omitted).
        %%%
        %%% Callers rotate ONCE, in their own constructor, and must then pass 0
        %%% (never obj.angle) to any solver called afterwards — the solvers
        %%% rotate internally too. See docs/reference/architecture/forwardCalcs.md, "Rotation
        %%% convention: pass 0 to solvers after construction".
            structRot = structure;
            for iLay = 1:structure.nLay
                mat0 = structure.mater(iLay);
                mat0 = mat0.SetDirMat(mat0.GetRotDirMat(angleDeg));
                structRot = structRot.ModifyLayer(iLay, mat0);
            end
        end

        function useP = resolveParallel(useParallelSetting, nItems)
        %%% Resolve a useParallel setting into an actual parallel/serial decision
        %%% >> useP = CalcCommon.resolveParallel(useParallelSetting, nItems)
        %%%
        %%% INPUT:  1. useParallelSetting .. "auto" | "always" | "never"
        %%%         2. nItems             .. number of independent iterations
        %%% OUTPUT: useP .. true only if the policy says yes AND a thread pool
        %%%                 is actually available
        %%%
        %%% The two ParallelPolicy calls must travel together: asking the policy
        %%% without confirming a pool exists runs `parfor` with no pool — silently
        %%% serial, plus a startup attempt per loop. Callers resolve this ONCE per
        %%% calc() and reuse the answer for every phase loop.
            useP = ParallelPolicy.shouldUseParallel(useParallelSetting, nItems) ...
                && ParallelPolicy.ensureThreadPool();
        end

        function outCell = forEachElement(inCell, useP, stepFcn)
        %%% Map one phase step over a cell of per-point calc objects, parallel or serial
        %%% >> outCell = CalcCommon.forEachElement(inCell, useP, stepFcn)
        %%%
        %%% INPUT:  1. inCell  .. {n x 1} cell of per-point Calc* objects
        %%%         2. useP    .. from CalcCommon.resolveParallel
        %%%         3. stepFcn .. function handle, outCell{i} = stepFcn(inCell{i}, i)
        %%% OUTPUT: outCell .. {n x 1}, same order as inCell
        %%%
        %%% The element index is passed as stepFcn's 2nd argument so a phase
        %%% needing per-element data (e.g. the flat 2D path's per-wavelength
        %%% Green velocity axis) can index into it from inside the handle,
        %%% rather than needing a second sliced loop variable.
        %%%
        %%% outCell is seeded from inCell and then both read and written at the
        %%% loop index — a plain sliced variable, which is what `parfor` requires.
            nItems  = numel(inCell);
            outCell = inCell;
            if useP
                parfor iItem = 1:nItems
                    outCell{iItem} = stepFcn(outCell{iItem}, iItem); %#ok<PFBNS>
                    %%% PFBNS: stepFcn is a broadcast variable. Harmless --
                    %%% ensureThreadPool guarantees a THREAD pool, whose workers share
                    %%% this process's memory, so a broadcast handle (and anything it
                    %%% captures) is not copied per worker.
                end
            else
                for iItem = 1:nItems
                    outCell{iItem} = stepFcn(outCell{iItem}, iItem);
                end
            end
        end

        function validateCalcMode(mode, errPrefix)
        %%% Validate a calcMode argument ('full' | 'ritzOnly')
        %%% >> CalcCommon.validateCalcMode(mode, errPrefix)
        %%%
        %%% INPUT:  1. mode      .. value to check
        %%%         2. errPrefix .. caller's error-ID prefix, e.g. 'CalcSingleHS:calcAll'
        %%%
        %%% Throws <errPrefix>:badMode.
            if ~ismember(mode, {'full', 'ritzOnly'})
                error([errPrefix ':badMode'], ...
                    "mode must be %s; got '%s'.", CalcCommon.orList({'full', 'ritzOnly'}), mode);
            end
        end

        function validateRetentionLevel(retentionLevel, errPrefix, allowed)
        %%% Validate a retentionLevel argument against the caller's own allowed set
        %%% >> CalcCommon.validateRetentionLevel(retentionLevel, errPrefix, allowed)
        %%%
        %%% INPUT:  1. retentionLevel .. value to check
        %%%         2. errPrefix      .. caller's error-ID prefix, e.g. 'CalcAngSweepHS:calc'
        %%%         3. allowed        .. cellstr of permitted values
        %%%
        %%% `allowed` is explicit rather than a fixed list because the families
        %%% differ: the SAW side has {'full','saw_only','none'}, the Plate side
        %%% only {'full','none'} — every converged Lamb mode is physically real,
        %%% so there is no "significant subset" to prune to.
        %%%
        %%% Throws <errPrefix>:badRetentionLevel.
            if ~ismember(retentionLevel, allowed)
                error([errPrefix ':badRetentionLevel'], ...
                    "retentionLevel must be %s; got '%s'.", CalcCommon.orList(allowed), retentionLevel);
            end
        end

        function validateConstructorFlags(calcMode, retentionLevel, errPrefix, allowedRetention)
        %%% Validate a constructor's calcMode / retentionLevel arguments (each skipped if empty)
        %%% >> CalcCommon.validateConstructorFlags(calcMode, retentionLevel, errPrefix, allowedRetention)
        %%%
        %%% So a bad flag fails at construction, not at the first calc(). Throws
        %%% <errPrefix>:badMode / <errPrefix>:badRetentionLevel (validateCalcMode,
        %%% validateRetentionLevel); errPrefix is e.g. 'CalcSingleHS:constructor'.
            if ~isempty(calcMode); CalcCommon.validateCalcMode(calcMode, errPrefix); end
            if ~isempty(retentionLevel)
                CalcCommon.validateRetentionLevel(retentionLevel, errPrefix, allowedRetention);
            end
        end

        function validateSawThreshold(t, errPrefix)
        %%% Validate a sawThreshold: [] (the method's default) | 'high' | 'low' | a finite real number
        %%% >> CalcCommon.validateSawThreshold(t, errPrefix)
        %%%
        %%% Throws <errPrefix>:badSawThreshold. Checked where a classification starts
        %%% (classify, calcAll, a sweep's calc), so one id names the caller whatever score
        %%% class would read the threshold.
            if isempty(t); return; end
            if ischar(t) || isstring(t)
                ok = isscalar(string(t)) && ismember(lower(char(t)), {'high', 'low'});
            else
                ok = isa(t, 'double') && isscalar(t) && isreal(t) && isfinite(t);
            end
            if ~ok
                error([errPrefix ':badSawThreshold'], ...
                    "sawThreshold must be 'high', 'low', a finite real number or [] (the method's default).");
            end
        end

        function method = resolveSawClassMethod(method, threshold, mode, isLayered, errPrefix)
        %%% The classification a calculation runs: an unset sawClassMethod ('') follows calcMode
        %%% >> method = CalcCommon.resolveSawClassMethod(method, threshold, mode, isLayered, errPrefix)
        %%%
        %%% INPUT:  1. method    .. sawClassMethod: '' (unset) | 'green' | 'refined' | 'ritz'
        %%%         2. threshold .. sawThreshold: [] | 'high' | 'low' | numeric
        %%%         3. mode      .. calcMode: 'full' | 'ritzOnly'
        %%%         4. isLayered .. layered stack (no 'refined')
        %%%         5. errPrefix .. e.g. 'CalcSingleHS:calcAll'
        %%% OUTPUT: method .. as char: unset -> modeSawClassMethod(mode) ('green' on 'full',
        %%%                  'ritz' on 'ritzOnly')
        %%%
        %%% Validates both inputs. Throws <errPrefix>:numericThresholdNeedsMethod (a number
        %%% is read on one method's own scale, so it needs the method set explicitly),
        %%% <errPrefix>:ritzOnlyNeedsRitz (any method but 'ritz' on 'ritzOnly'), besides
        %%% validateSawThreshold's / validateSawClassMethod's own errors.
            CalcCommon.validateSawThreshold(threshold, errPrefix);
            if isstring(method); method = char(method); end
            if isempty(method)
                if isnumeric(threshold) && ~isempty(threshold)
                    error([errPrefix ':numericThresholdNeedsMethod'], ...
                        ['A numeric sawThreshold needs sawClassMethod set explicitly (each method ' ...
                         'has its own scale); or use ''high'', ''low'' or [] (the default).']);
                end
                method = CalcCommon.modeSawClassMethod(mode);
            end
            method = CalcCommon.validateSawClassMethod(method, isLayered, errPrefix);
            if strcmp(mode, 'ritzOnly') && ~strcmp(method, 'ritz')
                error([errPrefix ':ritzOnlyNeedsRitz'], ...
                    ['calcMode ''ritzOnly'' classifies with ''ritz'' only; sawClassMethod is ''%s'' ' ...
                     '(set it to ''ritz'', or leave it unset).'], method);
            end
        end

        function method = modeSawClassMethod(mode)
        %%% The classification method of a calcMode when sawClassMethod is unset
        %%% >> method = CalcCommon.modeSawClassMethod(mode)
        %%%
        %%% OUTPUT: method .. 'ritz' for 'ritzOnly', 'green' otherwise ('full')
            if strcmp(mode, 'ritzOnly'); method = 'ritz'; else; method = 'green'; end
        end

        function [ids, names] = classificationMethods(isLayered)
        %%% The classification methods of a geometry: ids and display names
        %%% >> [ids, names] = CalcCommon.classificationMethods(isLayered)
        %%%
        %%% OUTPUT: ids   .. {'green', 'refined', 'ritz'}; layered {'green', 'ritz'}
        %%%                  (there is no spot score for layered stacks)
        %%%         names .. the matching display names ('Green', 'Refined', 'Ritz')
        %%% The one list the validation, the viewers' classification bar and the
        %%% StructureControlGUI method dropdown read.
            ids   = {'green', 'refined', 'ritz'};
            names = {'Green', 'Refined', 'Ritz'};
            if isLayered
                keep = ~strcmp(ids, 'refined');
                ids = ids(keep); names = names(keep);
            end
        end

        function method = validateSawClassMethod(method, isLayered, errPrefix)
        %%% Validate a classification method: 'green' | 'refined' | 'ritz' (layered: no 'refined')
        %%% >> method = CalcCommon.validateSawClassMethod(method, isLayered, errPrefix)
        %%%
        %%% OUTPUT: method .. as char. Throws <errPrefix>:badSawClassMethod, or
        %%%                  <errPrefix>:noRefinedLayered for 'refined' on a layered stack
        %%%                  (there is no spot score for layered stacks).
            if isstring(method); method = char(method); end
            if isLayered && ischar(method) && strcmp(method, 'refined')
                error([errPrefix ':noRefinedLayered'], ...
                    "Method 'refined' is not available for layered stacks (no spot score for layered stacks).");
            end
            allowed = CalcCommon.classificationMethods(isLayered);
            if ~ischar(method) || ~ismember(method, allowed)
                error([errPrefix ':badSawClassMethod'], "sawClassMethod must be %s.", ...
                    CalcCommon.orList(allowed));
            end
        end

        function [value, scale] = thresholdValue(t, high, low, errId, fname)
        %%% Value and scale of a standard: 'high' | 'low' | a number (read on the 'high' scale)
        %%% >> [value, scale] = CalcCommon.thresholdValue(t, high, low, errId, fname)
        %%%
        %%% INPUT:  1. t           .. 'high' | 'low' (char or string) | real finite scalar
        %%%         2-3. high, low .. the score class's values of the two standards
        %%%         4. errId       .. error id for an unknown name (the caller's own)
        %%%         5. fname       .. caller name for validateattributes
        %%% OUTPUT: value .. the threshold; scale .. 'high' | 'low' (a number: 'high')
        %%%
        %%% The one resolver behind GreenScore / GreenSpotScore / GreenScoreLayered /
        %%% RitzScore .resolveThreshold*: each passes its own constants and error id.
            if ischar(t) || isstring(t)
                switch lower(char(t))
                    case 'high'; value = high; scale = 'high';
                    case 'low';  value = low;  scale = 'low';
                    otherwise
                        error(errId, 'threshold must be ''high'', ''low'' or numeric; got ''%s''', char(t));
                end
            else
                validateattributes(t, {'double'}, {'scalar', 'real', 'finite'}, fname, 'threshold');
                value = t; scale = 'high';
            end
        end

        function thr = defaultThreshold(method, isLayered)
        %%% The default standard of a classification method
        %%% >> thr = CalcCommon.defaultThreshold(method, isLayered)
        %%%
        %%% Half-space: 'green' 'high', 'refined' 'high', 'ritz' 'low' (what Ritz-only is
        %%% reliable for). Layered: 'low' for both (sagittal-type waves; 'high' adds the
        %%% Love-type ones).
            if isLayered || strcmp(method, 'ritz')
                thr = 'low';
            else
                thr = 'high';
            end
        end

        function [precompStruct, built] = ritzSAWPrecompute(structure, precompStruct, rawByDepthCell)
        %%% Build any missing field of the angle-independent RitzSolverSAW precompute
        %%% >> [precompStruct, built] = CalcCommon.ritzSAWPrecompute(structure, precompStruct, rawByDepthCell)
        %%%
        %%% INPUT:  1. structure      .. Structure to build from
        %%%  (opt)  2. precompStruct  .. partially-filled precompute (default: empty
        %%%                             struct, i.e. build everything)
        %%%  (opt)  3. rawByDepthCell .. pre-built GammaMatrix.rawByDepth cell to hand
        %%%                             to derivsByDepth (default: {}, let
        %%%                             derivsByDepth build its own)
        %%% OUTPUT: 1. precompStruct .. with nullConstrMat, dKdC_cell and d2KdCdD_cell
        %%%                             populated
        %%%         2. built         .. true if any field had to be built (false: the
        %%%                             input was complete and is returned unchanged,
        %%%                             at no cost)
        %%%
        %%% Every field depends only on calcParams (legOrder/relDepth) and density
        %%% — never on dirMat or propagation angle. That is what makes it correct,
        %%% not merely convenient, for a sweep class to build this ONCE from the
        %%% unrotated structure and share it across every angle.
        %%%
        %%% Fill-if-missing throughout, so one call serves both cases: a sweep
        %%% passes nothing and gets everything built; a single-point class passes
        %%% whatever its sweep handed it and only the gaps are filled.
        %%%
        %%% Lamb-side precompute is deliberately NOT handled here: it needs
        %%% RitzSolver.buildConstraintMatrix(..., 'lamb') (a different boundary
        %%% condition) and must NOT carry dKdC_cell at all.
        %%%
        %%% >> precompStruct = CalcCommon.ritzSAWPrecompute(structure, sweep)
        %%% With a CalcCommon.ritzWavelengthSweep struct as 2nd argument: the full precompute
        %%% of a structure that differs from the sweep's only in its layer depths
        %%% (wavelength), elastic constants or densities; see ritzWavelengthSweep.
            if nargin < 2 || isempty(precompStruct); precompStruct = struct(); end
            if nargin < 3; rawByDepthCell = {}; end
            if CalcCommon.isRitzWavelengthSweep(precompStruct)
                precompStruct = CalcCommon.ritzSAWPrecomputeFromSweep(structure, precompStruct);
                built = true;
                return
            end
            built = false;

            if ~isfield(precompStruct, 'nullConstrMat') || isempty(precompStruct.nullConstrMat)
                built = true;
                precompStruct.nullConstrMat = ...
                    RitzSolver.buildConstraintMatrix(structure.calcParams, 'saw');
            end
            if ~isfield(precompStruct, 'dKdC_cell')    || isempty(precompStruct.dKdC_cell) ...
            || ~isfield(precompStruct, 'd2KdCdD_cell') || isempty(precompStruct.d2KdCdD_cell)
                built = true;
                if isempty(rawByDepthCell)
                    [precompStruct.dKdC_cell, precompStruct.d2KdCdD_cell] = ...
                        GammaMatrix.derivsByDepth(structure);
                else
                    %%% derivsByDepth accepts a pre-built rawByDepthCell as its own 2nd
                    %%% argument — skips rebuilding the expensive, relDepth-independent
                    %%% Legendre-integral tensor once per wavelength in a 2D sweep.
                    [precompStruct.dKdC_cell, precompStruct.d2KdCdD_cell] = ...
                        GammaMatrix.derivsByDepth(structure, rawByDepthCell);
                end
            end
        end

        function sweep = ritzWavelengthSweep(structure)
        %%% The wavelength-independent part of the RitzSolverSAW precompute of a structure
        %%% >> sweep = CalcCommon.ritzWavelengthSweep(structure)
        %%%
        %%% INPUT:  structure .. Structure to build from
        %%% OUTPUT: sweep .. struct: legOrders [1 x nLay] (the Legendre order of every layer,
        %%%                  which is all the other fields depend on), rawByDepthCell
        %%%                  (GammaMatrix.RawCellFor) and nullConstrMat
        %%%
        %%% Mirrors the Julia ritzWavelengthSweep / RitzWavelengthSweep. ApplyWavelength
        %%% changes only the films' relDepth, while the raw Legendre tensor and the constraint
        %%% null space N depend on the Legendre orders alone, so a sweep over wavelengths
        %%% (|k|) builds them once; at each wavelength
        %%% CalcCommon.ritzSAWPrecompute(structureAtWavelength, sweep) then only contracts
        %%% the depth factors. Elastic constants and densities may change between those
        %%% structures too (neither the raw tensor nor N carries them); the Legendre orders
        %%% may not. Julia's sweep also caches the real-gauge back-map and the Gram blocks of
        %%% M_con; MATLAB keeps neither (measured < 1 % of a solve here).
            validateattributes(structure, {'Structure'}, {'scalar'}, ...
                'CalcCommon.ritzWavelengthSweep', 'structure');
            cps   = structure.calcParams;
            sweep = struct('legOrders', [cps.legOrder], ...
                'rawByDepthCell', {GammaMatrix.RawCellFor(structure)}, ...
                'nullConstrMat', RitzSolver.buildConstraintMatrix(cps, 'saw'));
        end

        function warnRitzPrecomputeBuiltOnWorker(className)
        %%% Warn when a Calc* object had to build its Ritz precompute on a pool worker
        %%% >> CalcCommon.warnRitzPrecomputeBuiltOnWorker(className)
        %%%
        %%% INPUT: className .. class of the object, e.g. 'CalcSingleHS'
        %%%
        %%% Called by ritzPrecompute() after it had to build something; silent on the
        %%% client. The precompute is Transient, and a thread-pool parfor serialises each
        %%% object it hands to a worker, so a worker's copy arrives without it. Rebuilding
        %%% it there costs more than the Ritz solve itself, once per angle -- the sweeps
        %%% therefore pass their shared copy into calcRitz()/calcAll() explicitly. This
        %%% warning means some parallel step does not, i.e. a lost speed-up, not a wrong
        %%% result. Not switched off after the first time: it should stay loud. (A
        %%% worker's warning is printed on the client but does not set lastwarn; the
        %%% tests capture it with evalc.)
            if ParallelPolicy.isOnWorker()
                warning([className ':ritzPrecomputeBuiltOnWorker'], ...
                    ['%s: the Ritz precompute was rebuilt on a parallel worker (the sweep did ' ...
                     'not pass its shared copy) -- correct, but this costs more than the Ritz ' ...
                     'solve. Pass the precompute to calcRitz()/calcAll().'], className);
            end
        end

        function vRef = surfaceVRef(structure)
        %%% sqrt(G_Hill / rho) of the surface layer [km/s]: the score gates' time scale
        %%% >> vRef = CalcCommon.surfaceVRef(structure)
        %%%
        %%% INPUT: structure .. any orientation (the Hill average is rotation-invariant);
        %%%                     the surface layer is mater(nLay)
            top = structure.mater(structure.nLay);
            [~, gHill] = top.elConsts.GetHillAverages();
            vRef = sqrt(gHill / top.rho);
        end

        function lim = limitingAndBulk(mat)
        %%% Christoffel limiting velocities + cheap direct bulk velocities for ONE material
        %%% >> lim = CalcCommon.limitingAndBulk(mat)
        %%%
        %%% INPUT:  mat .. Material, taken from an ALREADY-ROTATED Structure
        %%% OUTPUT: lim .. struct with fields longLimVel [1x1], longLimPol [3x1],
        %%%                longLimTheta [1x1], shearLimVels [nShear x 1],
        %%%                shearLimPols [3 x nShear], shearLimIwt [nShear x 1],
        %%%                shearLimTheta [nShear x 1], vS / vTf [1x1] (slowest /
        %%%                fastest shear limiting velocity, whatever the branch: the
        %%%                Rayleigh boundary and the half-space candidate cap),
        %%%                bulkVels [3x1] = [vL; vT1; vT2] descending,
        %%%                bulkPols [3x3] (columns = corresponding polarisations)
        %%%
        %%% A LimitingVels.findRoots sagittal-plane scan (longitudinal branch
        %%% iwt=1, shear branches iwt=2,3 concatenated) plus one Christoffel
        %%% eigendecomposition along [1;0;0] for the bulk velocities.
        %%%
        %%% shearLimVels CONCATENATES two physically distinct shear branches, each
        %%% with a variable number of roots, so the branch label is carried
        %%% separately in shearLimIwt (2 or 3 per entry). Never max()/min()
        %%% shearLimVels to recover fast- vs slow-shear: the slow branch can supply
        %%% the maximum, and either branch can have no root at a given angle.
        %%% shearLimTheta/longLimTheta are the sagittal root angles theta*, which
        %%% LimitingVels.refineRoot needs as its bracket centre.
        %%%
        %%% WHICH material to pass is the caller's physics decision. CalcSingleHS
        %%% passes mater(nLay), the surface layer — for a homogeneous half-space
        %%% that is the whole medium, so the result is a genuine system-wide
        %%% limiting velocity. CalcSingleLayered passes mater(1), the substrate,
        %%% and stores the result under sub*-prefixed names: a SOFT reference
        %%% velocity (the threshold above which a mode is flagged leaky), NOT a
        %%% system-wide bound, since a layered stack has none. Swapping the two
        %%% indices is a physics error in either direction.
        %%%
        %%% No further rotation here (GetRotElastTensor(0)): the caller's
        %%% Structure already carries the propagation-angle rotation.
            rotElTensor = mat.GetStandardElastTensor(0);
            rho         = mat.rho;

            [vLim_sA, eigVec_sA, thetaRoots_sA] = LimitingVels.findRoots(rotElTensor, rho);

            %%% Longitudinal: physically exactly one root. unique() rather than a bare
            %%% (1) so that, should the islocalmin scan ever report two distinct iwt=1
            %%% roots, the lowest velocity is taken rather than whichever came first
            %%% in theta.
            [velL, iL] = unique(vLim_sA{1});
            lim.longLimVel   = velL(1);
            lim.longLimPol   = eigVec_sA{1}(iL(1), :).';
            lim.longLimTheta = thetaRoots_sA{1}(iL(1));

            lim.shearLimVels  = [vLim_sA{2}(:);    vLim_sA{3}(:)];
            lim.shearLimPols  = [eigVec_sA{2}; eigVec_sA{3}].';
            lim.shearLimIwt   = [2*ones(numel(vLim_sA{2}), 1); 3*ones(numel(vLim_sA{3}), 1)];
            lim.shearLimTheta = [thetaRoots_sA{2}(:); thetaRoots_sA{3}(:)];
            lim.vS  = min(lim.shearLimVels);
            lim.vTf = max(lim.shearLimVels);

            [eigVals, eigVecs] = Christoffel.matrix(rotElTensor, [1; 0; 0], rho);
            lim.bulkVels = sqrt(eigVals);
            lim.bulkPols = eigVecs;
        end

        function g = greenScoreSpectrum(structure, velAxis, fname)
        %%% The one Green solve of a 'full' angle: G13_eff at SawSpectrum.ETA
        %%% >> g = CalcCommon.greenScoreSpectrum(structure, velAxis, fname)
        %%%
        %%% INPUT:  1. structure .. ALREADY-ROTATED Structure for this angle
        %%%         2. velAxis   .. phase-velocity axis [km/s]; in the Calc pipeline a
        %%%                        uniform axis from SawSpectrum.VEL_MIN
        %%%         3. fname     .. caller name for validateattributes messages,
        %%%                        e.g. 'CalcSingleHS.calcGreen'
        %%% OUTPUT: g .. struct with
        %%%              velAxis   [nVel x 1] double
        %%%              G13Eff    [nVel x 1] complex double: G13_eff = (G13_+ - G13_-)/2,
        %%%                        the existence-score input
        %%%              absG13Eff [nVel x 1] single: abs(G13Eff), what the viewers draw
        %%%
        %%% ONE solve serves display and classification, and G13_eff is all a Calc
        %%% result keeps: it is the TGS response (in-plane force, out-of-plane
        %%% displacement, +-k symmetrised) and costs no more than any one-sided
        %%% component -- reciprocity gives G13_- = G31_+ from the same +k solve
        %%% (docs/reference/architecture/green.md, "+-k symmetrisation"). Any other component
        %%% is recomputed on request (greenComponentSpectrum), never stored. The
        %%% damping is SawSpectrum.ETA (3e-3), not GreenSolver's default 3e-4, because the
        %%% score is calibrated at it; every displayed peak therefore has the
        %%% linewidth the classifier saw.
        %%%
        %%% The magnitude is single precision (display needs no more); the complex
        %%% G13Eff stays double for the score's late-time FFT.
        %%%
        %%% The angle argument is 0, NEVER the caller's obj.angle: the caller's
        %%% constructor has already rotated the structure and GreenSolver.calc
        %%% rotates internally too, so passing the angle rotates twice. See
        %%% docs/reference/architecture/forwardCalcs.md, "Rotation convention: pass 0 to
        %%% solvers after construction".
        %%%
        %%% isoWavevec side only. CalcIsofreqSingleLayered keeps its own
        %%% GreenSolver.calc call: scanMode='isoFreq' with a real omega0 is a
        %%% physically different measurement (fixed frequency, swept |k∥|).
            validateattributes(velAxis, {'double'}, {'vector', 'real', 'finite', 'positive'}, ...
                fname, 'velAxis');

            greenResult = GreenSolver.calc(structure, 0, velAxis(:).', 1, 3, 'Eta', SawSpectrum.ETA);

            g.velAxis   = velAxis(:);
            g.G13Eff    = greenResult.G_eff(:);               % (G13_+ - G31_+)/2
            g.absG13Eff = single(abs(g.G13Eff));
        end

        function absG = greenComponentSpectrum(structure, velAxis, component, fname)
        %%% |G| of any one Green component, recomputed on request (never stored)
        %%% >> absG = CalcCommon.greenComponentSpectrum(structure, velAxis, component, fname)
        %%%
        %%% INPUT:  1. structure .. ALREADY-ROTATED Structure for this angle
        %%%         2. velAxis   .. phase-velocity axis [km/s], normally the result's
        %%%                        own greenVelAxis
        %%%         3. component .. 'GIJ' (one-sided G_IJ(+k), I = force direction,
        %%%                        J = displacement direction, 1 = propagation, 3 =
        %%%                        surface normal -- e.g. 'G13' = G13_+, 'G33') or
        %%%                        'GIJeff' (+-k symmetrised, e.g. 'G13eff')
        %%%         4. fname     .. caller name for error messages
        %%% OUTPUT: absG .. [nVel x 1] single, matches velAxis
        %%%
        %%% Same damping (SawSpectrum.ETA) and the same single solve as
        %%% greenScoreSpectrum, so any component compares like for like with the
        %%% stored |G13_eff|. Angle argument 0 (the double-rotation trap).
            validateattributes(velAxis, {'double'}, {'vector', 'real', 'finite', 'positive'}, ...
                fname, 'velAxis');
            [I, J, isEff] = CalcCommon.parseGreenComponent(component);
            greenResult = GreenSolver.calc(structure, 0, velAxis(:).', I, J, 'ComputeMinus', isEff, 'Eta', SawSpectrum.ETA);
            if isEff
                absG = single(abs(greenResult.G_eff(:)));
            else
                absG = single(abs(greenResult.G(:)));
            end
        end

        function [I, J, isEff] = parseGreenComponent(component)
        %%% 'GIJ' / 'GIJeff' -> force index I, displacement index J, symmetrised flag
        %%% >> [I, J, isEff] = CalcCommon.parseGreenComponent('G13eff')   % 1, 3, true
            tok = regexp(char(component), '^G([123])([123])(eff)?$', 'tokens', 'once');
            if isempty(tok)
                error('CalcCommon:parseGreenComponent:badComponent', ...
                    ['Green component must be ''GIJ'' or ''GIJeff'' with I, J in 1..3 ' ...
                     '(e.g. ''G13eff'', ''G13'', ''G33''); got ''%s''.'], char(component));
            end
            I = str2double(tok{1});
            J = str2double(tok{2});
            isEff = numel(tok) >= 3 && ~isempty(tok{3});
        end

        function g = greenSpectrum(structure, velAxis, fname)
        %%% One isoWavevec |G33| solve at GreenSolver's default damping (not the Calc pipeline)
        %%% >> g = CalcCommon.greenSpectrum(structure, velAxis, fname)
        %%%
        %%% INPUT:  1. structure .. ALREADY-ROTATED Structure for this angle
        %%%         2. velAxis   .. phase-velocity axis [km/s]
        %%%         3. fname     .. caller name for validateattributes messages
        %%% OUTPUT: g .. struct with velAxis [nVel x 1], absGFull [3 x 3 x nVel],
        %%%              absG33 [nVel x 1] (double)
        %%%
        %%% GreenSolver's default eta (3e-4), any axis: the isoWavevec |G33| reference
        %%% the isoFreq tests and validations compare against. The Calc classes' 'full'
        %%% mode uses greenScoreSpectrum instead.
        %%%
        %%% Angle argument 0, never the caller's angle (double-rotation trap).
        %%% computeMinus=false: for G33 (a diagonal pair) G_- = G_+ exactly.
            validateattributes(velAxis, {'double'}, {'vector', 'real', 'finite', 'positive'}, ...
                fname, 'velAxis');

            greenResult = GreenSolver.calc(structure, 0, velAxis, 3, 3, 'ComputeMinus', false);

            g.velAxis  = velAxis(:);
            g.absGFull = abs(greenResult.GFull);
            g.absG33   = squeeze(g.absGFull(3, 3, :));   % [nVel x 1], matches velAxis
        end

        function velAxis = buildAdaptiveVelAxis(candidateVels, velRange, denseHalfWidth, nDensePerPeak, nSparse)
        %%% Non-uniform Green velocity axis: dense near candidate velocities, sparse elsewhere
        %%% >> velAxis = CalcCommon.buildAdaptiveVelAxis(candidateVels, velRange, denseHalfWidth, nDensePerPeak, nSparse)
        %%%
        %%% Standalone helper, not used by the Calc classes: their axis is uniform,
        %%% because the existence scores need uniform sampling.
        %%%
        %%% INPUT:  1. candidateVels  .. velocities worth resolving finely (e.g. Ritz mode
        %%%                             velocities and the Christoffel limiting velocities)
        %%%         2. velRange       .. [vMin vMax], the range to cover
        %%%  (opt)  3. denseHalfWidth .. km/s half-width of each dense window (default 0.03)
        %%%  (opt)  4. nDensePerPeak  .. points per dense window (default 150)
        %%%  (opt)  5. nSparse        .. background-fill points across the whole range
        %%%                             (default 150)
        %%% OUTPUT: velAxis .. sorted, de-duplicated column vector: a dense linspace
        %%%                    window around every candidate (clamped to velRange),
        %%%                    unioned with one sparse linspace fill
        %%%
        %%% The defaults are chosen, not derived, and were calibrated on Fe3Al only --
        %%% retune per material class if needed. A bound (non-leaky) Rayleigh peak's
        %%% HEIGHT stays sampling-sensitive even at this density (it sits on a
        %%% near-singular pole); only its existence and location are well-defined,
        %%% so this axis is for finding WHERE features are, not for converged peak
        %%% heights.
        %%%
        %%% Empty candidateVels degrades to a plain uniform axis (the sparse fill
        %%% only). GreenSolver and GreenPeaks.find both accept a non-uniform axis.
            if nargin < 3 || isempty(denseHalfWidth); denseHalfWidth = 0.03; end
            if nargin < 4 || isempty(nDensePerPeak);  nDensePerPeak  = 150;  end
            if nargin < 5 || isempty(nSparse);        nSparse        = 150;  end

            vMin = velRange(1);
            vMax = velRange(2);

            candidateVels = candidateVels(:).';
            candidateVels = candidateVels(candidateVels >= vMin & candidateVels <= vMax);

            denseParts = cell(1, numel(candidateVels));
            for k = 1:numel(candidateVels)
                v0 = candidateVels(k);
                denseParts{k} = linspace(max(vMin, v0 - denseHalfWidth), ...
                                          min(vMax, v0 + denseHalfWidth), nDensePerPeak);
            end
            velAxis = unique([denseParts{:}, linspace(vMin, vMax, nSparse)]).';
        end

        function ritzSolver = pruneRitzTo(ritzSolver, classification, retentionLevel, errPrefix)
        %%% Apply the SAW-family eigenvector retention policy to a RitzSolverSAW result
        %%% >> ritzSolver = CalcCommon.pruneRitzTo(ritzSolver, classification, retentionLevel, errPrefix)
        %%%
        %%% INPUT:  1. ritzSolver      .. RitzSolverSAW result to prune
        %%%         2. classification  .. that point's ModeClassification result
        %%%         3. retentionLevel  .. 'full' | 'saw_only' | 'none'
        %%%         4. errPrefix       .. caller's error-ID prefix, e.g.
        %%%                              'CalcSingleHS:applyRetention'
        %%% OUTPUT: ritzSolver .. pruned (velocities/scores always survive)
        %%%
        %%% 'full':     keep eigVecs_con untouched.
        %%% 'saw_only': reduce to the columns identifyModes() classified as
        %%%             significant — multi-column in general (a pseudoRayleigh
        %%%             or SAW cluster plus limiting-mode entries), not a single
        %%%             sawModeNo.
        %%% 'none':     clear eigVecs_con entirely.
        %%%
        %%% SAW family only (CalcSingleHS/CalcSingleLayered). CalcSinglePlate keeps
        %%% its own two-level version: it has no 'saw_only' tier, and routing it
        %%% through here would make a plate silently ACCEPT 'saw_only'.
            switch retentionLevel
                case 'full'
                    %%% nothing to do
                case 'saw_only'
                    if ~isempty(ritzSolver.eigVecs_con)
                        keepCols = [classification.significantModes.modeNo];
                        ritzSolver = ritzSolver.pruneEigVecs(keepCols);
                        %%% keepCols is [] when no significant modes were found --
                        %%% pruneEigVecs([]) clears eigVecs_con entirely, the right
                        %%% behaviour when nothing was classified as significant
                    end
                case 'none'
                    ritzSolver = ritzSolver.pruneEigVecs([]);
                otherwise
                    error([errPrefix ':badRetentionLevel'], ...
                        "retentionLevel must be %s; got '%s'.", ...
                        CalcCommon.orList({'full', 'saw_only', 'none'}), retentionLevel);
            end
        end

        function tf = isStale(cs)
        %%% Whether a classified Calc result (single angle, or any angle of a sweep) carries
        %%% another classification stamp than the current one
        %%% >> tf = CalcCommon.isStale(result)
        %%%
        %%% INPUT:  cs .. CalcSingleHS/CalcSingleLayered, or a sweep (its cs{:} are checked)
        %%% OUTPUT: tf .. true if classified (classification not empty) with a
        %%%               classificationVersion other than ModeClassification.CLASSIFICATION_VERSION
            if isprop(cs, 'cs')
                tf = any(cellfun(@CalcCommon.isStale, cs.cs));
                return
            end
            tf = ~isempty(cs.classification) && ...
                cs.classificationVersion ~= ModeClassification.CLASSIFICATION_VERSION;
        end

        function refuseStale(cs, errPrefix)
        %%% Error <errPrefix>:stale if the result's classification is stale (see isStale)
        %%% >> CalcCommon.refuseStale(result, errPrefix)
            if CalcCommon.isStale(cs)
                error([errPrefix ':stale'], ['This result was classified by an older rule set ' ...
                    '(classification stamp other than %d) and cannot be reclassified: recompute it.'], ...
                    ModeClassification.CLASSIFICATION_VERSION);
            end
        end

        function precomp = completeRitzPrecompute(structure, precomp, className)
        %%% The Ritz precompute with every field built, warning if built on a pool worker
        %%% >> precomp = CalcCommon.completeRitzPrecompute(structure, precomp, className)
        %%%
        %%% ritzSAWPrecompute (fill-if-missing, free when complete) plus
        %%% warnRitzPrecomputeBuiltOnWorker(className) when anything had to be built: the
        %%% common body of the Calc classes' ritzPrecompute() and calcRitz(precomp).
            [precomp, built] = CalcCommon.ritzSAWPrecompute(structure, precomp);
            if built; CalcCommon.warnRitzPrecomputeBuiltOnWorker(className); end
        end

        function absG = greenComponentOf(cs, component)
        %%% |G| of any Green component on a single-angle result's greenVelAxis (not stored)
        %%% >> absG = CalcCommon.greenComponentOf(cs, component)
        %%%
        %%% The body of CalcSingleHS/CalcSingleLayered.greenComponent: 'G13eff' returns the
        %%% stored greenAbsG13Eff, any other component is one solve on the same axis and
        %%% damping (greenComponentSpectrum). Errors <class>:greenComponent:noGreenAxis
        %%% before calcGreen().
            if isempty(cs.greenVelAxis)
                error([class(cs) ':greenComponent:noGreenAxis'], ...
                    'calcGreen() (or calcAll(''full'')) must run first -- greenVelAxis is empty.');
            end
            if strcmp(component, 'G13eff') && ~isempty(cs.greenAbsG13Eff)
                absG = cs.greenAbsG13Eff;
                return
            end
            absG = CalcCommon.greenComponentSpectrum(cs.structure, cs.greenVelAxis, component, ...
                [class(cs) '.greenComponent']);
        end

        function absGMap = greenComponentMapOf(sweep, component)
        %%% Any Green component over a sweep on its greenVelAxis (not stored)
        %%% >> absGMap = CalcCommon.greenComponentMapOf(sweep, component)
        %%%
        %%% The body of CalcAngSweepHS/CalcAngSweepLayered.greenComponentMap: 'G13eff'
        %%% returns the stored map, any other component one solve per angle
        %%% (greenComponentMapFrom). Errors <class>:greenComponentMap:notCalculated before
        %%% calc('full').
            if isempty(sweep.greenVelAxis)
                error([class(sweep) ':greenComponentMap:notCalculated'], ...
                    'calc() must be run with calcMode=''full'' first.');
            end
            if strcmp(component, 'G13eff') && ~isempty(sweep.greenAbsG13EffMap)
                absGMap = sweep.greenAbsG13EffMap;
                return
            end
            absGMap = CalcCommon.greenComponentMapFrom(sweep.cs, component, sweep.useParallel);
        end

        function [classification, limitingEntries] = classifyGreenPath(cs, thr, christoffelResult, greenResult)
        %%% Classify one angle on a Green path, associating or carrying the limiting entries
        %%% >> [classification, limitingEntries] = CalcCommon.classifyGreenPath(cs, thr, christoffelResult, greenResult)
        %%%
        %%% INPUT:  cs .. CalcSingleHS/CalcSingleLayered (ritzSolver, structure, limitingDone,
        %%%               limitingEntries read); thr .. the standard; christoffelResult,
        %%%               greenResult .. ModeClassification.calc inputs
        %%% OUTPUT: classification   .. ModeClassification result
        %%%         limitingEntries  .. the limiting entries this calculation carries (the
        %%%                             new association on the first classification after a
        %%%                             calculation, else cs.limitingEntries unchanged)
        %%%
        %%% The limiting association needs the eigenvectors, which retention may prune, and
        %%% does not depend on the threshold: it runs once per calculation and every later
        %%% classification carries its entries (associateLimiting = false).
            if cs.limitingDone
                mc = ModeClassification('threshold', thr, 'associateLimiting', false, ...
                    'carryEntries', cs.limitingEntries);
            else
                mc = ModeClassification('threshold', thr);
            end
            classification = mc.calc(cs.ritzSolver, cs.structure, christoffelResult, greenResult);
            limitingEntries = cs.limitingEntries;
            if ~cs.limitingDone
                sm = classification.significantModes;
                limitingEntries = sm(strcmp({sm.type}, 'limiting'));
            end
        end

        function cs1 = greenStep(cs0, greenVelAxis, retentionLevel, greenScoreRef, greenThreshold, precomp)
        %%% Phase 3 body for one angle: the Green solve -> existence score -> classify -> prune
        %%% >> cs1 = CalcCommon.greenStep(cs0, greenVelAxis, retentionLevel, greenScoreRef, greenThreshold, precomp)
        %%%
        %%% INPUT:  1. cs0            .. one CalcSingleHS/CalcSingleLayered, Phases 1-2 done
        %%%         2. greenVelAxis   .. the sweep-wide common axis (every angle's own axis
        %%%                              is a leading part of it)
        %%%         3. retentionLevel .. 'full' | 'saw_only' | 'none'
        %%%  (opt)  4. greenScoreRef    .. the reference shared by the whole sweep (half-space
        %%%                              only): GreenScore.reference (computed per object if
        %%%                              omitted). Ignored for a layered angle:
        %%%                              GreenScoreLayered needs no reference
        %%%  (opt)  5. greenThreshold .. the 'green' standard: 'high' | 'low' | numeric
        %%%                              (default: the method's default)
        %%%  (opt)  6. precomp        .. the shared Ritz precompute (the Ritz candidates'
        %%%                              angular slopes; default: cs0's own, empty on a
        %%%                              pool worker)
        %%% OUTPUT: cs1 .. same object with Green, classification and retention applied
        %%%
        %%% One Green solve per angle: calcGreen stores the display magnitudes on the
        %%% common axis and the complex score input; calcGreenScore adds the probe solves
        %%% and the scores (half-space GreenScore, or GreenScoreLayered for a layered
        %%% angle); calcRitzScore stores the Ritz candidates from the eigenvectors
        %%% still in memory, so 'ritz' can be selected later at no cost. The angle is
        %%% stamped calcMode 'full' and classified 'green' (the sweep switches the method
        %%% afterwards if asked).
        %%%
        %%% ORDERING IS LOAD-BEARING: prune runs AFTER identifyModes, never
        %%% before — classification needs the full eigenvector set to know which
        %%% columns are significant. Pruning per angle inside the loop, rather
        %%% than in a post-loop, keeps peak RAM at ONE angle's eigenvector matrix
        %%% instead of every angle's at once.
            if nargin < 4; greenScoreRef = []; end
            if nargin < 5; greenThreshold = []; end
            if nargin < 6; precomp = []; end
            cs0 = cs0.setRunMode('full');
            cs1 = cs0.calcGreen(greenVelAxis);
            if isa(cs1, 'CalcSingleLayered')
                cs1 = cs1.calcGreenScore();
            else
                cs1 = cs1.calcGreenScore(greenScoreRef);
            end
            cs1 = cs1.calcRitzScore(precomp);   % the Ritz candidates too (no solve)
            cs1 = cs1.identifyModes('green', greenThreshold);
            cs1 = cs1.applyRetention(retentionLevel);
            cs1.retentionLevel = retentionLevel;
        end

        function cs1 = ritzOnlyStep(cs0, retentionLevel, ritzThreshold, precomp)
        %%% 'ritzOnly' body for one angle: Ritz -> (Christoffel + modal sum) -> classify -> prune
        %%% >> cs1 = CalcCommon.ritzOnlyStep(cs0, retentionLevel, ritzThreshold, precomp)
        %%%
        %%% INPUT:  1. cs0                .. one CalcSingleHS/CalcSingleLayered
        %%%         2. retentionLevel     .. 'full' | 'saw_only' | 'none'
        %%%         3. ritzThreshold      .. 'high' | 'low' | numeric ([] = the default)
        %%%  (opt)  4. precomp            .. the sweep's shared Ritz precompute, handed to
        %%%                                  calcRitz (default: cs0's own -- empty on a
        %%%                                  pool worker, see warnRitzPrecomputeBuiltOnWorker)
        %%% OUTPUT: cs1 .. same object, calcAll('ritzOnly') done
        %%%
        %%% The Ritz-only counterpart of greenStep. No classification input is shared
        %%% between angles: each angle scores against its own floor, so a sweep and a
        %%% single angle classify identically. calcAll keeps the load-bearing order
        %%% (modal-sum weights and classification BEFORE pruning).
            if nargin < 4; precomp = []; end
            cs0.sawClassMethod = 'ritz';           % the angle carries the run's method
            cs0.sawThreshold   = ritzThreshold;
            cs1 = cs0.calcAll('ritzOnly', retentionLevel, precomp);
        end

        function [sawVelMap, sawModeNoMap, sawTypeMap, significantModesMap] = ...
                sawMapsFrom(csCell, sawTypes)
        %%% Per-angle SAW convenience maps, extracted from a cell of single-point objects
        %%% >> [sawVelMap, sawModeNoMap, sawTypeMap, significantModesMap] = ...
        %%%        CalcCommon.sawMapsFrom(csCell, sawTypes)
        %%%
        %%% INPUT:  1. csCell   .. {nAngles x 1} of CalcSingleHS/CalcSingleLayered,
        %%%                       each with classification already filled
        %%%         2. sawTypes .. cellstr of significantModes types counting as SAW
        %%% OUTPUT: four {nAngles x 1} cell maps (ragged per angle by nature —
        %%%         the number of significant modes varies with angle)
        %%%
        %%% Extraction only, no classification logic. `sawTypes` is explicit
        %%% because the families filter on different type names and the difference
        %%% is load-bearing: CalcAngSweepHS uses {'Rayleigh','pseudoRayleigh'},
        %%% CalcAngSweepLayered {'SAW'} (both its Green and Ritz-only paths emit 'SAW').
            nAngles             = numel(csCell);
            sawVelMap           = cell(nAngles, 1);
            sawModeNoMap        = cell(nAngles, 1);
            sawTypeMap          = cell(nAngles, 1);
            significantModesMap = cell(nAngles, 1);

            for iAngle = 1:nAngles
                sm = csCell{iAngle}.classification.significantModes;
                significantModesMap{iAngle} = sm;
                isSAW = ismember({sm.type}, sawTypes);
                sawVelMap{iAngle}    = [sm(isSAW).vel];
                sawModeNoMap{iAngle} = [sm(isSAW).modeNo];
                sawTypeMap{iAngle}   = {sm(isSAW).type};
                %%% sm can be 0x0 (struct with named fields, no elements) if no significant
                %%% mode was found at this angle -- {sm.type} is then a valid 1x0 cellstr and
                %%% every line above degrades gracefully to an empty result, not an error.
            end
        end

        function greenAbsG13EffMap = greenG13EffMapFrom(csCell)
        %%% The sweep's |G13_eff| map, assembled from a cell of single-point objects
        %%% >> greenAbsG13EffMap = CalcCommon.greenG13EffMapFrom(csCell)
        %%%
        %%% INPUT:  csCell .. {nAngles x 1} of CalcSingleHS/CalcSingleLayered, each
        %%%                   with greenAbsG13Eff filled on ONE common axis (calcGreen())
        %%% OUTPUT: greenAbsG13EffMap .. [nVel x nAngles] single; [] if any element
        %%%                   lacks greenAbsG13Eff
        %%%
        %%% Assembly only — no Green solve here. FULL MODE ONLY: callers guard this
        %%% behind a calcMode check.
            if ~all(cellfun(@(c) ~isempty(c.greenAbsG13Eff), csCell(:)))
                greenAbsG13EffMap = single([]);
                return
            end
            greenAbsG13EffMap = cell2mat(cellfun(@(c) c.greenAbsG13Eff(:), csCell(:)', ...
                'UniformOutput', false));
        end

        function absGMap = greenComponentMapFrom(csCell, component, useParallel)
        %%% Any Green component over a sweep, recomputed per angle on request
        %%% >> absGMap = CalcCommon.greenComponentMapFrom(csCell, component, useParallel)
        %%%
        %%% INPUT:  1. csCell      .. {nAngles x 1} of CalcSingleHS/CalcSingleLayered
        %%%                           with calcGreen() run on ONE common axis
        %%%         2. component   .. 'GIJ' | 'GIJeff' (CalcCommon.parseGreenComponent)
        %%%         3. useParallel .. as the sweep's own useParallel
        %%% OUTPUT: absGMap .. [nVel x nAngles] single, nothing stored
        %%%
        %%% One Green solve per angle (each element's greenComponent()), so a whole
        %%% sweep costs about as much as its original Green phase.
            useP  = CalcCommon.resolveParallel(useParallel, numel(csCell));
            cols  = CalcCommon.forEachElement(csCell(:), useP, @(cs, ~) cs.greenComponent(component));
            absGMap = cell2mat(cellfun(@(c) c(:), cols(:)', 'UniformOutput', false));
        end

        function [csFlat, outerIdx, angleIdx] = flattenSweeps(sweeps, nAngles)
        %%% Flatten a {nOuter x 1} cell of angle-sweep objects into one linear list
        %%% >> [csFlat, outerIdx, angleIdx] = CalcCommon.flattenSweeps(sweeps, nAngles)
        %%%
        %%% INPUT:  1. sweeps  .. {nOuter x 1} of CalcAngSweepLayered or
        %%%                      CalcIsofreqAngSweepLayered (anything exposing .cs)
        %%%         2. nAngles .. angles per sweep (identical for every element)
        %%% OUTPUT: csFlat   .. {nOuter*nAngles x 1} of the per-point objects
        %%%         outerIdx .. [N x 1] which sweeps{} each element came from
        %%%         angleIdx .. [N x 1] which cs{} slot within that sweep
        %%%
        %%% The flat-dispatch (flattenSweep=true) counterpart to the nested path:
        %%% it lets ParallelPolicy see the PRODUCT of both axes rather than
        %%% nAngles alone, which matters when neither axis alone crosses
        %%% AUTO_THRESHOLD but their product does. Pair with scatterToSweeps.
            nOuter   = numel(sweeps);
            N        = nOuter * nAngles;
            csFlat   = cell(N, 1);
            outerIdx = zeros(N, 1);
            angleIdx = zeros(N, 1);
            ell = 0;
            for iOuter = 1:nOuter
                for iAngle = 1:nAngles
                    ell = ell + 1;
                    csFlat{ell}   = sweeps{iOuter}.cs{iAngle};
                    outerIdx(ell) = iOuter;
                    angleIdx(ell) = iAngle;
                end
            end
        end

        function sweeps = scatterToSweeps(sweeps, csFlat, outerIdx, angleIdx)
        %%% Scatter a flattened result list back into the sweep objects it came from
        %%% >> sweeps = CalcCommon.scatterToSweeps(sweeps, csFlat, outerIdx, angleIdx)
        %%%
        %%% INPUT:  1. sweeps                    .. {nOuter x 1} of angle-sweep objects
        %%%         2. csFlat, outerIdx, angleIdx .. as returned by flattenSweeps
        %%% OUTPUT: sweeps .. same cell, each element updated
        %%%
        %%% Writes through each sweep class's own narrow setElement() accessor
        %%% rather than relaxing SetAccess on its protected cs property.
            for ell = 1:numel(csFlat)
                sweeps{outerIdx(ell)} = sweeps{outerIdx(ell)}.setElement(angleIdx(ell), csFlat{ell});
            end
        end

        function names = resultClasses()
        %%% Every Calc* result class saveCalc/loadCalc (and Gui.viewResult) accept
        %%% >> names = CalcCommon.resultClasses()
        %%%
        %%% OUTPUT: names .. {1 x n} cellstr of class names
            names = {'CalcSingleHS', 'CalcAngSweepHS', ...
                     'CalcSingleLayered', 'CalcAngSweepLayered', 'CalcWavevecAngSweep2Layered', ...
                     'CalcIsofreqSingleLayered', 'CalcIsofreqAngSweepLayered', 'CalcFreqAngSweep2Layered', ...
                     'CalcSinglePlate', 'CalcAngSweepPlate', 'CalcWavevecAngSweep2Plate'};
        end

        function filePath = saveCalc(calcResult, filePath, note)
        %%% Save a finished Calc* result, with provenance metadata, to a .mat file
        %%% >> filePath = CalcCommon.saveCalc(calcResult, filePath, note)
        %%%
        %%% INPUT:  1. calcResult .. any CalcCommon.resultClasses() object, normally
        %%%                          already calc()'d / calcAll()'d
        %%%  (opt)  2. filePath   .. target .mat file (default = '': an auto-named
        %%%                          file in <project>/results/, see below);
        %%%                          '.mat' is appended if missing
        %%%  (opt)  3. note       .. free-text note stored in meta.note (default = '')
        %%% OUTPUT: filePath .. absolute path of the file written
        %%%
        %%% The file holds two variables: calcResult and meta (className, savedAt,
        %%% guliwerrVersion, gitCommit, matlabVersion, note). Load it back with
        %%% CalcCommon.loadCalc or straight into a viewer with Gui.viewSaved.
        %%%
        %%% Auto name: <project>/results/<class>_<structure name>_<yyyymmdd-HHMMSS>.mat,
        %%% <project> being the folder that holds guliwerr/ -- fixed, so it does not
        %%% depend on the current folder of whichever session saved it.
        %%%
        %%% Always -v7.3 (HDF5-based): no 2 GB per-variable limit, and measured no
        %%% larger or slower than -v7 for these objects. The Ritz precompute is
        %%% Transient on every Calc* class and therefore not saved (see
        %%% CalcSingleHS.precompStructForRitz).
        %%%
        %%% guliwerrVersion is guliwerrVersion() at save time -- the provenance that
        %%% survives a download without .git. gitCommit is HEAD of the repository
        %%% holding guliwerr/ at save time; uncommitted edits are not captured.
            if nargin < 2; filePath = ''; end
            if nargin < 3 || isempty(note); note = ''; end
            CalcCommon.validateResultObject(calcResult, 'CalcCommon:saveCalc');
            if ~(ischar(note) || (isstring(note) && isscalar(note)))
                error('CalcCommon:saveCalc:badNote', 'note must be a character vector or string scalar.');
            end

            if isempty(filePath)
                filePath = CalcCommon.defaultResultPath(calcResult);
            else
                filePath = char(filePath);
                [~, ~, ext] = fileparts(filePath);
                if ~strcmpi(ext, '.mat'); filePath = [filePath '.mat']; end
            end
            targetDir = fileparts(filePath);
            if ~isempty(targetDir) && ~isfolder(targetDir)
                mkdir(targetDir);
            end

            meta = struct( ...
                'className',     class(calcResult), ...
                'savedAt',       char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), ...
                'guliwerrVersion', guliwerrVersion(), ...
                'gitCommit',     CalcCommon.gitCommit(), ...
                'matlabVersion', version, ...
                'note',          char(note));
            save(filePath, 'calcResult', 'meta', '-v7.3');

            %%% absolute path for the caller to hand on (e.g. to Gui.viewSaved)
            fileInfo = dir(filePath);
            filePath = fullfile(fileInfo.folder, fileInfo.name);
        end

        function filePath = defaultResultPath(calcResult)
        %%% saveCalc's auto-generated file name for a result
        %%% >> filePath = CalcCommon.defaultResultPath(calcResult)
        %%%
        %%% OUTPUT: filePath .. <project>/results/<class>_<structure name>_<yyyymmdd-HHMMSS>.mat
        %%%                     (structure name with anything but letters, digits,
        %%%                     '_' and '-' replaced by '_')
            stamp    = char(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
            safeName = regexprep(calcResult.structure.name, '[^\w-]', '_');
            filePath = fullfile(CalcCommon.projectRoot(), 'results', ...
                sprintf('%s_%s_%s.mat', class(calcResult), safeName, stamp));
        end

        function [calcResult, meta] = loadCalc(filePath)
        %%% Load a Calc* result saved by saveCalc (or by a plain save())
        %%% >> [calcResult, meta] = CalcCommon.loadCalc(filePath)
        %%%
        %%% INPUT:  1. filePath .. .mat file
        %%% OUTPUT: 1. calcResult .. the saved Calc* object
        %%%         2. meta       .. saveCalc's metadata struct (fields empty for a
        %%%                          file written by a plain save(); guliwerrVersion
        %%%                          empty for a file saved before it was recorded)
        %%%
        %%% A plain save() file is accepted if it holds exactly one Calc* result
        %%% variable. Warns (CalcCommon:loadCalc:commitMismatch) when the file was
        %%% saved at a different git commit than the current code -- harmless unless
        %%% the Calc* classes changed in between. Errors if MATLAB could not rebuild
        %%% the object and fell back to a struct, which is what load() does when
        %%% the class definition has changed incompatibly since saving.
            if ~(ischar(filePath) || (isstring(filePath) && isscalar(filePath)))
                error('CalcCommon:loadCalc:badPath', 'filePath must be a character vector or string scalar.');
            end
            filePath = char(filePath);
            if ~isfile(filePath)
                error('CalcCommon:loadCalc:fileNotFound', 'File not found: %s', filePath);
            end

            S = load(filePath);
            meta = struct('className', '', 'savedAt', '', 'guliwerrVersion', '', ...
                'gitCommit', '', 'matlabVersion', '', 'note', '');
            if isfield(S, 'calcResult')
                calcResult = S.calcResult;
                if isfield(S, 'meta'); meta = S.meta; end
                if ~isfield(meta, 'guliwerrVersion'); meta.guliwerrVersion = ''; end
            else
                %%% plain save(): take the one variable that is a Calc* result
                varNames = fieldnames(S);
                isResult = cellfun(@(v) CalcCommon.isResultObject(S.(v)), varNames);
                if nnz(isResult) ~= 1
                    error('CalcCommon:loadCalc:noUniqueResult', ...
                        ['%s holds %d Calc* result variables (expected exactly one, ' ...
                         'or a ''calcResult'' variable as written by saveCalc).'], ...
                        filePath, nnz(isResult));
                end
                calcResult = S.(varNames{isResult});
                meta.className = class(calcResult);
            end

            if ~CalcCommon.isResultObject(calcResult)
                error('CalcCommon:loadCalc:notAnObject', ...
                    ['%s: the saved result came back as a %s, not a Calc* object -- the ' ...
                     'class definition has most likely changed incompatibly since it was ' ...
                     'saved (meta.gitCommit = ''%s''). Check out that commit to load it, ' ...
                     'or recompute.'], filePath, class(calcResult), meta.gitCommit);
            end

            nowCommit = CalcCommon.gitCommit();
            if ~isempty(meta.gitCommit) && ~isempty(nowCommit) && ~strcmp(meta.gitCommit, nowCommit)
                warning('CalcCommon:loadCalc:commitMismatch', ...
                    ['%s was saved at commit %s; the code is now at %s. Fine unless the ' ...
                     'Calc* classes changed in between.'], filePath, ...
                    meta.gitCommit(1:min(7, end)), nowCommit(1:min(7, end)));
            end
        end

        function s = orList(items)
        %%% Format a cellstr as a quoted "'a', 'b', or 'c'" fragment for error messages
        %%% >> s = CalcCommon.orList(items)
        %%%
        %%% Lets the validators above build the message from the allowed-value
        %%% list itself, so the two cannot drift apart.
            quoted = cellfun(@(c) sprintf('''%s''', c), items, 'UniformOutput', false);
            if isscalar(quoted)
                s = quoted{1};
            elseif numel(quoted) == 2
                s = sprintf('%s or %s', quoted{1}, quoted{2});
            else
                s = sprintf('%s, or %s', strjoin(quoted(1:end-1), ', '), quoted{end});
            end
        end

        function tf = isResultObject(x)
        %%% True if x is a scalar object of one of CalcCommon.resultClasses()
        %%% >> tf = CalcCommon.isResultObject(x)
            tf = isscalar(x) && isobject(x) && ismember(class(x), CalcCommon.resultClasses());
        end

        function validateResultObject(x, errPrefix)
        %%% Throw <errPrefix>:badResult unless x is a Calc* result object
        %%% >> CalcCommon.validateResultObject(x, errPrefix)
            if ~CalcCommon.isResultObject(x)
                error([errPrefix ':badResult'], ...
                    'Expected a scalar Calc* result object (%s); got a %s.', ...
                    CalcCommon.orList(CalcCommon.resultClasses()), class(x));
            end
        end

        function root = projectRoot()
        %%% The folder holding guliwerr/ (this file is guliwerr/forwardCalcs/CalcCommon.m)
        %%% >> root = CalcCommon.projectRoot()
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
        end

        function sha = gitCommit()
        %%% HEAD commit hash of the repository holding guliwerr/, '' if not found
        %%% >> sha = CalcCommon.gitCommit()
        %%%
        %%% Read straight from .git/HEAD (and the ref file or packed-refs it points
        %%% to) rather than by calling git through system(): no dependency on a git
        %%% binary on the path, and no shell call. In a linked worktree .git is a
        %%% file "gitdir: <dir>": HEAD is read from <dir>, and the refs from <dir>
        %%% and from the shared repository named by <dir>/commondir.
            sha    = '';
            root   = CalcCommon.projectRoot();
            gitDir = fullfile(root, '.git');
            if isfile(gitDir)
                tok = regexp(fileread(gitDir), '(?m)^gitdir:\s*(.+?)\s*$', 'tokens', 'once');
                if isempty(tok); return; end
                gitDir = CalcCommon.resolveGitPath(root, tok{1});
            end
            headFile = fullfile(gitDir, 'HEAD');
            if ~isfile(headFile); return; end
            head = strtrim(fileread(headFile));
            if ~startsWith(head, 'ref: ')
                sha = head;   % detached HEAD holds the hash itself
                return
            end
            refName = strtrim(head(6:end));
            refDirs = {gitDir};
            commonFile = fullfile(gitDir, 'commondir');
            if isfile(commonFile)
                refDirs{end+1} = CalcCommon.resolveGitPath(gitDir, strtrim(fileread(commonFile)));
            end
            for iDir = 1:numel(refDirs)
                refFile = fullfile(refDirs{iDir}, refName);
                if isfile(refFile)
                    sha = strtrim(fileread(refFile));
                    return
                end
            end
            for iDir = 1:numel(refDirs)
                packedFile = fullfile(refDirs{iDir}, 'packed-refs');
                if ~isfile(packedFile); continue; end
                %%% packed-refs lines: "<sha> <refName>"
                tok = regexp(fileread(packedFile), ...
                    ['(?m)^([0-9a-f]{40}) ' regexptranslate('escape', refName) '\s*$'], 'tokens', 'once');
                if ~isempty(tok); sha = tok{1}; return; end
            end
        end
    end

    methods (Static, Access = private)
        function p = resolveGitPath(baseDir, p)
        %%% A path from a .git file or a commondir file: absolute as given, else
        %%% relative to baseDir
        %%% >> p = CalcCommon.resolveGitPath(baseDir, p)
            isAbs = startsWith(p, '/') || ~isempty(regexp(p, '^[A-Za-z]:[\\/]', 'once'));
            if ~isAbs
                p = fullfile(baseDir, p);
            end
        end

        function tf = isRitzWavelengthSweep(s)
        %%% Is s a CalcCommon.ritzWavelengthSweep struct (rather than a precompute)?
        %%% >> tf = CalcCommon.isRitzWavelengthSweep(s)
            tf = isstruct(s) && isscalar(s) && isfield(s, 'legOrders') && isfield(s, 'rawByDepthCell');
        end

        function precomp = ritzSAWPrecomputeFromSweep(structure, sweep)
        %%% The full RitzSolverSAW precompute of a structure from a wavelength sweep
        %%% >> precomp = CalcCommon.ritzSAWPrecomputeFromSweep(structure, sweep)
        %%%
        %%% dKdC_cell/d2KdCdD_cell are identical (isequal) to the from-scratch build, and
        %%% nullConstrMat is the sweep's own matrix (the same values; shared in RAM by
        %%% copy-on-write). Error CalcCommon:ritzSAWPrecompute:sweepMismatch when the
        %%% Legendre orders differ from the sweep's.
            cps = structure.calcParams;
            if ~isequal([cps.legOrder], sweep.legOrders)
                error('CalcCommon:ritzSAWPrecompute:sweepMismatch', ...
                    ['The structure''s Legendre orders differ from those the wavelength sweep ' ...
                     'was built for.']);
            end
            [dKdC_cell, d2KdCdD_cell] = GammaMatrix.derivsByDepth(structure, sweep.rawByDepthCell);
            precomp = struct('nullConstrMat', sweep.nullConstrMat, ...
                'dKdC_cell', {dKdC_cell}, 'd2KdCdD_cell', {d2KdCdD_cell});
        end
    end
end
