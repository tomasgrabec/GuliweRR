classdef Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The GuliweRR computational domain — a stack of layers, each described by a Material
%  and a CalcParams object.
%
% Value class: AddLayer, RemoveLayer and ModifyLayer return the modified object.
%
%%% GuliweRR v6.0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = protected)
        name            char          = 'str 1'   % structure name
        nMaxModes       double        = 500       % maximum number of computed modes

        mater           Material                  % layer array: Material(1)=bottom … Material(nLay)=top
        calcParams      CalcParams                % layer array: CalcParams(1)=bottom … CalcParams(nLay)=top

        nLay            double        = 0         % number of layers in the structure
        strBasisVecL    double        = 0         % sum of basisVecL across all layers
    end

    %% Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function obj = Structure(materIn,calcParamsIn,name,nMaxModes)
            %%% Initialise Structure object, optionally adding a first layer
            %%% >> obj = Structure(materIn,calcParamsIn,name,nMaxModes)
            %%%
            %%% INPUT (all opt): 1. materIn      .. Material for the first layer; [] adds no layer
            %%%                  2. calcParamsIn .. CalcParams for it (default if empty)
            %%%                  3. name         .. structure name [string]
            %%%                  4. nMaxModes    .. max number of modes to compute (default = 500)
            %%%
            %%% An empty materIn is what suppresses AddLayer, not the argument count, so
            %%% Structure([],[],'strRef',500) sets name/nMaxModes on an empty structure.

            if nargin > 0 && ~isempty(materIn)
                if nargin < 2 || isempty(calcParamsIn)
                    calcParamsIn = CalcParams();
                end
                obj = obj.AddLayer(materIn,calcParamsIn);
            end
            if nargin > 2 && ~isempty(name)
                obj.name = name;
            end
            if nargin > 3 && ~isempty(nMaxModes)
                obj.nMaxModes = nMaxModes;
            end
        end

        function obj = AddLayer(obj,mater0,calcParams0)
            %%% Append a layer to the top of the structure
            %%% >> obj = AddLayer(obj,mater0,calcParams0)
            %%%
            %%% INPUT:  1. mater0      .. Material for the new layer (default if empty)
            %%%  (opt)  2. calcParams0 .. CalcParams for the new layer (default if empty)

            if nargin < 2 || isempty(mater0)
                mater0 = Material();
            end
            if nargin < 3 || isempty(calcParams0)
                calcParams0 = CalcParams();
            end
            Structure.checkLayerInputs(mater0, calcParams0, 'AddLayer');

            obj.nLay                   = obj.nLay + 1;
            obj.mater(obj.nLay)        = mater0;
            obj.calcParams(obj.nLay)   = calcParams0;
            obj.strBasisVecL           = obj.strBasisVecL + calcParams0.basisVecL;
        end

        function obj = RemoveLayer(obj,ind)
            %%% Remove a layer from the structure by index
            %%% >> obj = RemoveLayer(obj,ind)
            %%%
            %%% INPUT: ind .. index of the layer to remove (1 = bottom, nLay = top)

            if ~isnumeric(ind) || ~isscalar(ind)
                error('Structure:RemoveLayer:badInput', ...
                    'Layer index must be a numeric scalar');
            end
            if ind < 1 || ind > obj.nLay
                error('Structure:RemoveLayer:outOfRange', ...
                    'Layer index %d is out of range (structure has %d layers)', ind, obj.nLay);
            end

            %%% read basisVecL before the deletion, not after
            removedBasisVecL = obj.calcParams(ind).basisVecL;

            obj.nLay             = obj.nLay - 1;
            obj.mater(ind)       = [];
            obj.calcParams(ind)  = [];
            obj.strBasisVecL     = obj.strBasisVecL - removedBasisVecL;
        end

        function obj = ModifyLayer(obj,indLay,mater0,calcParams0)
            %%% Replace the material and/or calcParams of an existing layer
            %%% >> obj = ModifyLayer(obj,indLay,mater0,calcParams0)
            %%%
            %%% INPUT:  1. indLay      .. layer index to modify (empty → top-most layer)
            %%%         2. mater0      .. new Material (empty → keep existing)
            %%%  (opt)  3. calcParams0 .. new CalcParams (empty → keep existing)

            if isempty(indLay)
                indLay = obj.nLay;
            end
            if ~isnumeric(indLay) || ~isscalar(indLay)
                error('Structure:ModifyLayer:badInput', ...
                    'Layer index must be a numeric scalar');
            end
            if indLay < 1 || indLay > obj.nLay
                error('Structure:ModifyLayer:outOfRange', ...
                    'Layer index %d is out of range (structure has %d layers)', indLay, obj.nLay);
            end

            if isempty(mater0)
                mater0 = obj.mater(indLay);
            end
            if nargin < 4 || isempty(calcParams0)
                calcParams0 = obj.calcParams(indLay);
            end
            Structure.checkLayerInputs(mater0, calcParams0, 'ModifyLayer');

            oldBasisVecL = obj.calcParams(indLay).basisVecL;

            obj.mater(indLay)      = mater0;
            obj.calcParams(indLay) = calcParams0;
            obj.strBasisVecL       = obj.strBasisVecL - oldBasisVecL + calcParams0.basisVecL;
        end

        function obj = ApplyWavelength(obj,wavelengthReal)
            %%% Derive relDepth = realThickness / wavelengthReal for every FILM layer
            %%% >> obj = ApplyWavelength(obj,wavelengthReal)
            %%%
            %%% INPUT: wavelengthReal .. real wavelength, same units as realThickness, which
            %%%                          every film layer (2..nLay) must have set
            %%%
            %%% The substrate (layer 1) is skipped: legOrder stays fixed across a wavelength
            %%% sweep, so a derived substrate relDepth would grow arbitrarily large at a
            %%% long-wavelength point and undermesh the substrate's Legendre basis. Set it
            %%% deliberately with SetRelDepth if needed.

            validateattributes(wavelengthReal, {'double'}, {'scalar','real','finite','positive'}, ...
                'Structure.ApplyWavelength', 'wavelengthReal');

            for ilay = 2:obj.nLay   % layer 1 = substrate, skipped
                cp0 = obj.calcParams(ilay);
                if isnan(cp0.realThickness)
                    error('Structure:ApplyWavelength:missingRealThickness', ...
                        'Layer %d has no realThickness set -- cannot derive relDepth.', ilay);
                end
                cp0 = cp0.SetRelDepth(cp0.realThickness / wavelengthReal);
                obj = obj.ModifyLayer(ilay, [], cp0);
            end
        end

        function tf = isMirrorSymmetric(obj, tol)
            %%% Is the plate invariant under the mid-plane mirror z -> D-z?
            %%% >> tf = isMirrorSymmetric(obj, tol)
            %%%
            %%% INPUT (opt): tol .. relative tolerance for the comparisons (default = 1e-10)
            %%% OUTPUT: tf .. true iff, for every i and j = nLay+1-i (the middle layer of an
            %%%               odd stack paired with itself), relDepth and density match and
            %%%               layer j's standard-frame tensor equals the MIRROR IMAGE of
            %%%               layer i's: C_j = M C_i M, M = diag(1,1,-1)
            %%%
            %%% This is the exact condition for the free-plate operator to commute with the
            %%% mirror, i.e. for every mode to have a definite mid-plane parity (u1, u2 even
            %%% and u3 odd, or the reverse). Pairing a layer with its mirror image, not with
            %%% an identical copy, makes one test cover both requirements: the stack must be
            %%% symmetric AND each crystal must look the same from either face -- a single
            %%% cubic (111) layer fails (3-fold axis along the normal, C14-type entries),
            %%% (001)/(110) and isotropic layers pass. It also accepts a stack of two
            %%% mirror-twin layers, and identical layers given by equivalent dirMats.
            %%%
            %%% The x3 mirror commutes with every in-plane rotation, so the result does not
            %%% depend on the propagation angle (checked at angle 0). Conservative: a plate
            %%% that is symmetric only after merging identical adjacent layers (e.g. one
            %%% material split 0.1/0.3) returns false; classifyBranches then falls back to
            %%% confidence, which is safe.
            %%%
            %%% Exact parity does NOT make the Lamb/TH half of the four mode-type labels
            %%% exact: that needs the sagittal mirror too (isSagittalMirrorSymmetric).
            if nargin < 2 || isempty(tol); tol = 1e-10; end

            signZ = Structure.mirrorSign(3);
            tf = true;
            for i = 1:ceil(obj.nLay / 2)
                j = obj.nLay + 1 - i;
                Ci = obj.mater(i).GetStandardElastTensor(0);
                Cj = obj.mater(j).GetStandardElastTensor(0);
                dI = obj.calcParams(i).relDepth;  dJ = obj.calcParams(j).relDepth;
                rI = obj.mater(i).rho;            rJ = obj.mater(j).rho;
                if abs(dI - dJ) > tol * max(abs(dI), abs(dJ)) ...
                        || abs(rI - rJ) > tol * max(abs(rI), abs(rJ)) ...
                        || max(abs(Cj - signZ .* Ci), [], 'all') > tol * max(abs(Ci), [], 'all')
                    tf = false;
                    return
                end
            end
        end

        function tf = isSagittalMirrorSymmetric(obj, angleDeg, tol)
            %%% Is the plate invariant under the sagittal mirror x2 -> -x2 at this angle?
            %%% >> tf = isSagittalMirrorSymmetric(obj, angleDeg, tol)
            %%%
            %%% INPUT: 1. angleDeg .. propagation angle [deg], as the sweeps take it
            %%%  (opt) 2. tol      .. relative tolerance (default = 1e-10)
            %%% OUTPUT: tf .. true iff every layer's standard-frame tensor at angleDeg is
            %%%               unchanged by M = diag(1,-1,1): C = M C M
            %%%
            %%% Where this holds, the SH (u2) motion decouples from the sagittal (u1, u3)
            %%% motion, so the Lamb/TH split of the mode labels is exact. Together with
            %%% isMirrorSymmetric the tensor is orthotropic in the propagation frame and all
            %%% four labels are exact. Depends on the angle (e.g. Cu(001): true along <100>
            %%% and <110>, false at 30 deg); layer order and thickness play no role.
            if nargin < 3 || isempty(tol); tol = 1e-10; end

            signY = Structure.mirrorSign(2);
            tf = true;
            for i = 1:obj.nLay
                C = obj.mater(i).GetStandardElastTensor(angleDeg);
                if max(abs(C - signY .* C), [], 'all') > tol * max(abs(C), [], 'all')
                    tf = false;
                    return
                end
            end
        end
    end

    methods (Static, Access = private)
        function s = mirrorSign(ax)
            %%% Sign pattern of a mirror normal to axis ax on a 4th-order tensor
            %%% >> s = Structure.mirrorSign(ax)
            %%%
            %%% M = diag with -1 at ax, so (M C M)_ijkl = s_ijkl C_ijkl with
            %%% s_ijkl = (-1)^(how many of i,j,k,l equal ax)
            v = ones(3, 1);  v(ax) = -1;
            s = reshape(kron(kron(v, v), kron(v, v)), 3, 3, 3, 3);
        end

        function checkLayerInputs(mater0, calcParams0, methodName)
            %%% Error Structure:<methodName>:badInput unless the layer inputs are one
            %%% Material and one CalcParams
            if ~(isa(mater0, 'Material') && isscalar(mater0))
                error(['Structure:' methodName ':badInput'], 'mater0 must be one Material object.');
            end
            if ~(isa(calcParams0, 'CalcParams') && isscalar(calcParams0))
                error(['Structure:' methodName ':badInput'], 'calcParams0 must be one CalcParams object.');
            end
        end
    end
end
