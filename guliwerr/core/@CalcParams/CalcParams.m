classdef CalcParams
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation parameters for the GuliweRR Ritz-Rayleigh solver: Legendre integration
%  matrices and the sampled Legendre basis.
%
% Value class. Set*/Get* methods take the object; other methods are static.
%
%%% GuliweRR v6.0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = protected)
        relDepth        double      = 20       % domain depth relative to the wavelength, which is
                                                 % fixed at 1 throughout the solvers (|k_par| = 1)
        legOrder        double      = 50       % max order of Legendre polynomials

        realThickness   double      = NaN      % real thickness, same units as the wavelength
                                                 % Structure.ApplyWavelength is called with;
                                                 % NaN = unset

        legDiag         double                 % Legendre integration matrix — diagonal part   (legOrder × legOrder)
        legNondiag      double                 % Legendre integration matrix — non-diagonal     (legOrder × legOrder)

        basisVecL       double      = 3*50     % length of the full basis vector = 3 × legOrder

        polyLeg         double                 % Legendre-function values at spatial mesh  (legOrder × nMesh)
    end

    %% Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function obj = CalcParams(relDepth,legOrder,realThickness)
            %%% Initialise CalcParams object
            %%% >> obj = CalcParams(relDepth,legOrder,realThickness)
            %%%
            %%% INPUT (opt): 1. relDepth      .. relative domain depth (default = 20)
            %%%              2. legOrder      .. max Legendre order (default = 50)
            %%%              3. realThickness .. real thickness, in the units ApplyWavelength
            %%%                                  will use (default = NaN, unset)

            if nargin > 0 && ~isempty(relDepth)
                obj = obj.SetRelDepth(relDepth);
            end
            if nargin > 1 && ~isempty(legOrder)
                CalcParams.checkLegOrder(legOrder);
                obj.legOrder = legOrder;
            end
            if nargin > 2 && ~isempty(realThickness)
                obj = obj.SetRealThickness(realThickness);
            end

            obj = obj.RecomputeLegBasis();
        end

        function obj = SetRelDepth(obj,newRelDepth)
            %%% Set relative depth of the domain (positive scalar; NaN = a film's placeholder
            %%% until Structure.ApplyWavelength derives it)
            %%% >> obj = SetRelDepth(obj,newRelDepth)
            CalcParams.checkPositiveOrNaN(newRelDepth, 'SetRelDepth', 'relDepth');
            obj.relDepth = newRelDepth;
        end
        function obj = SetRealThickness(obj,newRealThickness)
            %%% Set real physical thickness (same units ApplyWavelength will use; NaN = unset)
            %%% >> obj = SetRealThickness(obj,newRealThickness)
            CalcParams.checkPositiveOrNaN(newRealThickness, 'SetRealThickness', 'realThickness');
            obj.realThickness = newRealThickness;
        end
        function obj = SetLegOrder(obj,newLegOrder)
            %%% Set max Legendre polynomial order (positive integer) and recompute all
            %%% dependent matrices
            %%% >> obj = SetLegOrder(obj,newLegOrder)
            CalcParams.checkLegOrder(newLegOrder);
            obj.legOrder = newLegOrder;
            obj = obj.RecomputeLegBasis();
        end
    end

    %% Private methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = private)
        function obj = RecomputeLegBasis(obj)
            %%% Recompute every field derived from legOrder — sole owner of that
            %%% derivation, so no object can hold a mesh stale for its own legOrder
            %%% >> obj = RecomputeLegBasis(obj)
            [obj.legDiag,obj.legNondiag] = CalcParams.LegendreMatrixDiagNondiag(obj.legOrder);
            obj.basisVecL                = 3 * obj.legOrder;
            obj.polyLeg                  = CalcParams.BasisDecomposition2D(obj.legOrder, 0.01);
        end
    end

    methods (Static, Access = private)
        function checkLegOrder(legOrder)
            %%% legOrder must be a positive integer scalar
            validateattributes(legOrder, {'numeric'}, {'scalar','integer','positive'}, ...
                'CalcParams', 'legOrder');
        end
        function checkPositiveOrNaN(x, methodName, argName)
            %%% Error CalcParams:<methodName>:badInput unless x is a positive finite scalar or NaN
            if ~(isnumeric(x) && isscalar(x) && isreal(x) && (isnan(x) || (isfinite(x) && x > 0)))
                error(['CalcParams:' methodName ':badInput'], ...
                    '%s must be a positive finite scalar, or NaN (unset).', argName);
            end
        end
    end

    %% Static methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Static)
        [d,n] = LegendreMatrixDiagNondiag(legO)
        pLeg = BasisDecomposition2D(legO,meshSpacing)
    end
end
