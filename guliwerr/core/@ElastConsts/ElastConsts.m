classdef ElastConsts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Elastic constants in the three formats GuliweRR uses: independent-constant vector,
%  Voigt 6x6 matrix, and full 3x3x3x3 tensor, kept in sync with a symmetry class.
%
% Value class. Set*/Get* methods take the object; other methods are static.
%
%%% GuliweRR v6.0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = protected)
        vector              double      = NaN           % short vector format of independent constants
        voigt(6,6)          double      = NaN(6,6)      % voigt matrix   (6x6)
        tensor(3,3,3,3)     double      = NaN(3,3,3,3)  % tensor of elastic constants (3^4)
        symmClassNo(1,1)    double      = 0             % symmetry class number -- see SymmClassLibrary
        nIndepConsts(1,1)   double      = 0             % number of independent constants
        shearAnisotropyRatio(1,1) double = NaN          % (v_shear,max / v_shear,min)^2 over ALL
                                                        % directions, >= 1, rotation-invariant;
                                                        % stored only where a closed form exists:
                                                        % isotropic 1, cubic max(A_Z, 1/A_Z). NaN
                                                        % otherwise (no closed form;
                                                        % CalcAnisotropyRatio scans)
    end

    properties (Constant, Access = private)
        DEFAULT_ROT_MAT = [1 0 0; 0 0 1]'   % columns [propagation | surface normal]:
                                            % propagation [100], normal [001]
    end

    %% Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        %% Initialization method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = ElastConsts(consts,symmNo)
            %%% Initialise elastic constants object
            %%% >> obj = ElastConsts(consts,symmNo)
            %%%
            %%% INPUT: 1. consts .. elastic constants (vector / Voigt 6x6 / tensor 3x3x3x3)
            %%%  (opt) 2. symmNo .. symmetry class number (see SymmClassLibrary); if omitted,
            %%%                     inferred via FindSymmClass_input — see its caveats

            if nargin < 1  % no input: blank triclinic object of NaNs
                consts = NaN(6,6);
                symmNo = 10;
            elseif nargin < 2
                symmNo = ElastConsts.FindSymmClass_input(consts);
            end

            [~,expectedNumberOfConsts] = ElastConsts.SymmClassLibrary(symmNo);

            %%% for vector input, validate count before Vector2Voigt silently truncates
            if isvector(consts) && numel(consts) > 1 ...
                    && numel(consts) ~= expectedNumberOfConsts
                error('ElastConsts:constructor:wrongCount', ...
                    'Number of input constants does not match the symmetry class.')
            end

            obj.tensor = ElastConsts.Any2Tensor(consts,symmNo);
            obj.voigt = ElastConsts.Tensor2Voigt(obj.tensor);
            obj.vector = ElastConsts.Tensor2Vector(obj.tensor,symmNo);

            obj.nIndepConsts = length(obj.vector);

            if obj.nIndepConsts ~= expectedNumberOfConsts
                error('ElastConsts:constructor:wrongCount', ...
                    'Number of input constants does not match the symmetry class.')
            else
                obj.symmClassNo = symmNo;
            end

            %%% Shear anisotropy ratio, rotation-invariant, so stored once here instead of
            %%% being rescanned by every GreenSolver call (twice per angle of a sweep).
            %%% Only the free closed forms: a direction scan (~0.8 ms) would cost more
            %%% than this whole constructor, which finite-difference loops call many
            %%% times (SensitivityAnalysis). Cubic: the shear extremes are c44 and
            %%% (c11-c12)/2, so (vmax/vmin)^2 = max(A_Z, 1/A_Z) --
            %%% NOT the Zener ratio itself, which is < 1 for 2 c44 < c11 - c12.
            if symmNo == 1
                obj.shearAnisotropyRatio = 1;
            elseif symmNo == 2
                zener = ElastConsts.ZenerRatio(obj.vector);
                if isfinite(zener) && zener > 0
                    obj.shearAnisotropyRatio = max(zener, 1/zener);
                end
            end
        end

        %% Methods to modify object %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = SetNewSymmClass(obj,newSymmNo)
            %%% Rebuild the object under a different symmetry class (from its Voigt matrix)
            %%% >> obj = SetNewSymmClass(obj,newSymmNo)
            voigtMat = obj.voigt;
            obj = ElastConsts(voigtMat,newSymmNo);
        end
        %%  Object-facing wrappers of the static methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function redElast = GetReducedElastConsts(obj,rotMat)
            %%% Get the 15 reduced elastic constants (2D RR system) for a given rotation
            %%% >> redElast = GetReducedElastConsts(obj,rotMat)
            %%%
            %%% INPUT (opt): rotMat .. rotation matrix, columns [propagation | surface normal]
            %%%                       (3x2; default propagation [100], normal [001])
            if nargin < 2 || isempty(rotMat)
                rotMat = ElastConsts.DEFAULT_ROT_MAT;
            end
            indep2RedElastMat = ElastConsts.ElastRotReductMat(rotMat,obj.symmClassNo);
            redElast = indep2RedElastMat * obj.vector';
        end
        function E = GetYoungsModulus(obj,rotMat)
            %%% Young's modulus along the rotated x1 direction
            %%% >> E = GetYoungsModulus(obj,rotMat)
            %%%
            %%% INPUT (opt): rotMat .. rotation matrix, columns [propagation | surface normal]
            %%%                       (3x2; default propagation [100], normal [001])
            if nargin < 2 || isempty(rotMat)
                rotMat = ElastConsts.DEFAULT_ROT_MAT;
            end
            E = ElastConsts.Tensor2Youngs(obj.tensor,rotMat);
        end
        function [c11_voigt,c44_voigt] = GetVoigtAverages(obj)
            %%% Voigt averages: c11_voigt == M (longitudinal), c44_voigt == G (shear)
            %%% >> [c11_voigt,c44_voigt] = GetVoigtAverages(obj)
            [c11_voigt,c44_voigt] = ElastConsts.VoigtAverages(obj.voigt);
        end
        function [c11_reuss,c44_reuss] = GetReussAverages(obj)
            %%% Reuss averages: c11_reuss == M (longitudinal), c44_reuss == G (shear)
            %%% >> [c11_reuss,c44_reuss] = GetReussAverages(obj)
            [c11_reuss,c44_reuss] = ElastConsts.ReussAverages(obj.voigt);
        end
        function [c11_hill,c44_hill] = GetHillAverages(obj)
            %%% Hill averages: c11_hill == M (longitudinal), c44_hill == G (shear)
            %%% >> [c11_hill,c44_hill] = GetHillAverages(obj)
            [c11_hill,c44_hill] = ElastConsts.HillAverages(obj.voigt);
        end
        function A = GetAnisotropyRatio(obj,angleStep)
            %%% Anisotropy ratio: exact Zener 2*c44/(c11-c12) for cubic, otherwise the
            %%% generalised form (v_shear_max/v_shear_min)^2 from an angular scan
            %%% >> A = GetAnisotropyRatio(obj,angleStep)
            %%%
            %%% INPUT (opt): angleStep .. scan step [degrees] (default = 4); ignored for
            %%%                           isotropic/cubic, which need no scan
            if nargin < 2 || isempty(angleStep)
                angleStep = 4;
            end
            if obj.symmClassNo == 1
                A = 1;
            elseif obj.symmClassNo == 2
                A = ElastConsts.ZenerRatio(obj.vector);
            else
                A = ElastConsts.CalcAnisotropyRatio(obj.tensor,angleStep);
            end
        end

        function rotElTen = GetRotElastTensor(obj,rotMat)
            %%% Elastic tensor (3x3x3x3) rotated by the given rotation matrix
            %%% >> rotElTen = GetRotElastTensor(obj,rotMat)
            %%%
            %%% INPUT (opt): rotMat .. rotation matrix, columns [propagation | surface normal]
            %%%                        (3x2; default propagation [100], normal [001])
            if nargin < 2 || isempty(rotMat)
                rotMat = ElastConsts.DEFAULT_ROT_MAT;
            end
            rotElTen = ElastConsts.RotateTensor(obj.tensor,rotMat);
        end
        function rotElVoigt = GetRotElastVoigt(obj,rotMat)
            %%% Voigt matrix (6x6) of the elastic tensor rotated by the given rotation matrix
            %%% >> rotElVoigt = GetRotElastVoigt(obj,rotMat)
            %%%
            %%% INPUT (opt): rotMat .. rotation matrix, columns [propagation | surface normal]
            %%%                        (3x2; default propagation [100], normal [001])
            if nargin < 2 || isempty(rotMat)
                rotMat = ElastConsts.DEFAULT_ROT_MAT;
            end
            rotElTen = ElastConsts.RotateTensor(obj.tensor,rotMat);
            rotElVoigt = ElastConsts.Tensor2Voigt(rotElTen);
        end
        %% Symmetry-specific methods
        function [vecOut,transfMatrix] = GetCubicConversion(obj)
            %%% Equivalent cubic constants [B, c', c44] (cubic symmetry only)
            %%% >> [vecOut,transfMatrix] = GetCubicConversion(obj)
            %%%
            %%% OUTPUT: vecOut(3x1) = transfMatrix(3x3) * [c11;c12;c44]
            if obj.symmClassNo == 2
                [vecOut,transfMatrix] = ElastConsts.CubicConversion(obj.vector);
            else
                error('ElastConsts:GetCubicConversion:notCubic', ...
                    'Not possible for non-cubic elasticities.')
            end
        end
        function [vecOut,transfMatrix] = GetTetragonalConversion(obj)
            %%% Equivalent tetragonal constants [B, c33, c', c'', c44, c66] (tetragonal type I only)
            %%% >> [vecOut,transfMatrix] = GetTetragonalConversion(obj)
            %%%
            %%% OUTPUT: vecOut(6x1) = transfMatrix(6x6) * vecInput(6x1)
            if obj.symmClassNo == 3
                [vecOut,transfMatrix] = ElastConsts.TetragonalConversion(obj.vector);
            else
                error('ElastConsts:GetTetragonalConversion:notTetragonal', ...
                    'Not possible for non-tetragonal elasticity.')
            end
        end
        function [vecOut,transfMatrix] = GetHexagonalConversion(obj)
            %%% Equivalent hexagonal physical-mode constants [C11, C33, G1, G2, G3] (hexagonal only)
            %%% >> [vecOut,transfMatrix] = GetHexagonalConversion(obj)
            %%%
            %%% OUTPUT: vecOut(5x1) = transfMatrix(5x5) * [c11;c12;c13;c33;c44]
            if obj.symmClassNo == 7
                [vecOut,transfMatrix] = ElastConsts.HexagonalConversion(obj.vector);
            else
                error('ElastConsts:GetHexagonalConversion:notHexagonal', ...
                    'Not possible for non-hexagonal elasticity.')
            end
        end
    end
    %% Static methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Static)
        %% Symmetry class related
        symmNo = FindSymmClass_input(inp)
        [symmC,nConsts] = SymmClassLibrary(symmNo)
        %% Conversion methods
        outTen = Any2Tensor(inp,symmNo)
        cij = Vector2Voigt(inpVec,symmNo)
        cij = Tensor2Voigt(inpTen)
        vec = Voigt2Vector(cij,symmNo)
        vec = Tensor2Vector(inpTen,symmNo)
        posList = IndepVoigtPositions(symmNo)
        %% Derived moduli
        E = Tensor2Youngs(elTen,rotMat)
        B = Voigt2Bulk(elVoigt)
        [c11_voigt,c44_voigt] = VoigtAverages(elVoigt)
        [c11_reuss,c44_reuss] = ReussAverages(elVoigt)
        [c11_hill,c44_hill] = HillAverages(elVoigt)
        %% Symmetry-specific
        [vecOut,transfMat] = CubicConversion(vecIn)
        [vecOut,deTransfMat] = CubicDeconversion(vecIn)
        [vecOut,transfMat] = TetragonalConversion(vecIn)
        [vecOut,deTransfMat] = TetragonalDeconversion(vecIn)
        [vecOut,transfMat] = HexagonalConversion(vecIn)
        [vecOut,deTransfMat] = HexagonalDeconversion(vecIn)
        %% Rotation
        frame = RightHandFrame(twoVecs,noteMsg)
        rotTen = RotateTensor(inpTen,rotMat)
        %% Reduced constants
        [indSet_ReductElast,nReductElast] = GetIndOf2DReduct_15of21()
        indep2RedElastMat = ElastRotReductMat(rotMat,symmNo)
        %% Anisotropy
        A = CalcAnisotropyRatio(normElTen,angleStep)
    end

    methods (Static, Access = private)
        type = FindType(input)
        vecOut = AsColumnVec(vecIn,nExpected,fname)

        function A = ZenerRatio(cubicVec)
            %%% Zener ratio 2*c44/(c11-c12) of a cubic vector [c11, c12, c44]
            A = 2*cubicVec(3) / (cubicVec(1) - cubicVec(2));
        end
    end
end
