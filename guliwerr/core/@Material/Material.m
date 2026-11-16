classdef Material
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Material parameters for the GuliweRR Ritz-Rayleigh calculation: elastic constants,
%  density, and orientation (direction matrix).
%
% Value class. Set*/Get* methods take the object; other methods are static.
%
%%% GuliweRR v6.0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = protected)
        name         string       = 'mater'                     % material name

        elConsts     ElastConsts  = ElastConsts([100,30],1)     % elastic constants (ElastConsts object)
        rho(1,1)     double       = 1                           % density  [g·cm⁻³]

                                                                % orientation matrix, row vectors
                                                                % (right-hand: row3 = cross(row1,row2)):
        dirMat(3,3)  double       = [1  0  0; ...               %   row 1 — surface normal (plane)
                                     0  0  1; ...               %   row 2 — propagation direction
                                     0 -1  0];                  %   row 3 — in-plane completion
                                                                % Two frames derive from it, and they differ:
                                                                % CalcRotMat gives x1=propagation, x2=plane,
                                                                % x3=completion; GetStandardDirMat gives
                                                                % x1=propagation, x2=completion, x3=plane —
                                                                % the literature convention the solvers use.
    end

    %% Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function obj = Material(elConsts,rho,dirMat,name)
            %%% Initialise Material object
            %%% >> obj = Material(elConsts,rho,dirMat,name)
            %%%
            %%% INPUT (all opt): 1. elConsts .. ElastConsts object, numeric vector, or
            %%%                                 cell {consts} / {consts, symmNo}
            %%%                  2. rho      .. density [g·cm⁻³]
            %%%                  3. dirMat   .. orientation matrix (2×3, 3×2, or 3×3 row vectors)
            %%%                  4. name     .. material name [string]

            if nargin > 0 && ~isempty(elConsts)
                obj = obj.SetElastConsts(elConsts);
            end
            if nargin > 1 && ~isempty(rho)
                obj = obj.SetDensity(rho);
            end
            if nargin > 2 && ~isempty(dirMat)
                obj = obj.SetDirMat(dirMat);
            end
            if nargin > 3 && ~isempty(name)
                obj.name = name;
            end
        end

        %% Set methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = SetElastConsts(obj,elConsts)
            %%% Set elastic constants of the material
            %%% >> obj = SetElastConsts(obj,elConsts)
            %%%
            %%% INPUT: elConsts .. ElastConsts object, numeric vector, or cell {consts[,symmNo]}
            if isa(elConsts,'ElastConsts')
                obj.elConsts = elConsts;
            elseif iscell(elConsts)
                if length(elConsts) > 1
                    obj.elConsts = ElastConsts(elConsts{1}, elConsts{2});
                else
                    obj.elConsts = ElastConsts(elConsts{1});
                end
            elseif isnumeric(elConsts)
                obj.elConsts = ElastConsts(elConsts);
            else
                error('Material:SetElastConsts:badInput', ...
                    'elConsts must be an ElastConsts object, numeric vector, or cell array');
            end
        end
        function obj = SetDirMat(obj,matr)
            %%% Set orientation matrix; the result is always right-hand orthonormal
            %%% >> obj = SetDirMat(obj,matr)
            %%%
            %%% INPUT: matr .. 2×3, 3×2, or 3×3 matrix of direction row vectors
            obj.dirMat = Material.CompleteDirMat(matr);
        end
        function obj = SetDensity(obj,rho)
            %%% Set density of the material [g·cm⁻³] (positive scalar)
            %%% >> obj = SetDensity(obj,rho)
            validateattributes(rho, {'numeric'}, {'scalar','real','finite','positive'}, ...
                'Material.SetDensity', 'rho');
            obj.rho = rho;
        end
        function obj = SetName(obj,name)
            %%% Set material name
            %%% >> obj = SetName(obj,name)
            obj.name = name;
        end

        %% Rotation methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function rotMat = GetRotMat(obj,angleDeg,relativeToggle)
            %%% Get rotation matrix for a given in-plane rotation angle
            %%% >> rotMat = GetRotMat(obj,angleDeg,relativeToggle)
            %%%
            %%% INPUT:  1. angleDeg       .. rotation angle [degrees]
            %%%  (opt)  2. relativeToggle .. logical: false — rotate within obj.dirMat (default)
            %%%                             true — rotate in the absolute [100]/[010] frame
            %%% OUTPUT: rotMat .. 3×2, columns [rotated propagation | surface normal]
            if nargin < 3 || isempty(relativeToggle)
                relativeToggle = false;
            end
            if ~(isscalar(relativeToggle) && (islogical(relativeToggle) || ...
                    (isnumeric(relativeToggle) && any(relativeToggle == [0 1]))))
                error('Material:GetRotMat:badInput', ...
                    'relativeToggle must be false (relative to dirMat) or true (absolute frame)');
            end
            if relativeToggle
                dMat = Material.CompleteDirMat([1 0 0; 0 1 0]);
            else
                dMat = obj.dirMat;
            end
            rotMat = Material.CalcRotMat(dMat,angleDeg);
        end
        function rotElTensor = GetRotElastTensor(obj,angleDeg,relativeToggle)
            %%% Elastic tensor (3×3×3×3) rotated by the given in-plane angle
            %%% >> rotElTensor = GetRotElastTensor(obj,angleDeg,relativeToggle)
            %%%
            %%% INPUT:  1. angleDeg       .. rotation angle [degrees], required (0 = the
            %%%                             dirMat frame; the unrotated crystal tensor is
            %%%                             obj.elConsts.tensor)
            %%%  (opt)  2. relativeToggle .. false — relative to dirMat (default); true — absolute frame
            if nargin < 3
                relativeToggle = [];   % GetRotMat owns the default and its validation
            end
            if nargin < 2 || isempty(angleDeg)
                error('Material:GetRotElastTensor:noAngle', ...
                    'An angle is required (0 for the dirMat frame; obj.elConsts.tensor is the unrotated tensor).');
            end
            rotMat = obj.GetRotMat(angleDeg,relativeToggle);
            rotElTensor = obj.elConsts.GetRotElastTensor(rotMat);
        end
        function rotElVoigt = GetRotElastVoigt(obj,angleDeg,relativeToggle)
            %%% Voigt matrix (6×6) of the elastic tensor rotated by the given in-plane angle
            %%% >> rotElVoigt = GetRotElastVoigt(obj,angleDeg,relativeToggle)
            %%%
            %%% INPUT:  1. angleDeg       .. rotation angle [degrees], required (0 = the
            %%%                             dirMat frame; the unrotated crystal matrix is
            %%%                             obj.elConsts.voigt)
            %%%  (opt)  2. relativeToggle .. false — relative to dirMat (default); true — absolute frame
            if nargin < 3
                relativeToggle = [];   % GetRotMat owns the default and its validation
            end
            if nargin < 2 || isempty(angleDeg)
                error('Material:GetRotElastVoigt:noAngle', ...
                    'An angle is required (0 for the dirMat frame; obj.elConsts.voigt is the unrotated matrix).');
            end
            rotMat = obj.GetRotMat(angleDeg,relativeToggle);
            rotElVoigt = obj.elConsts.GetRotElastVoigt(rotMat);
        end
        function elTensor = GetStandardElastTensor(obj,angleDeg)
            %%% Elastic tensor (3×3×3×3) in the standard solver frame: x1 = propagation,
            %%% x2 = in-plane completion, x3 = surface normal (Every/Maznev convention)
            %%% >> elTensor = GetStandardElastTensor(obj,angleDeg)
            %%%
            %%% INPUT (opt): angleDeg .. propagation angle [degrees] within the surface
            %%%                          plane, from dirMat row 2 (default = 0)
            %%%
            %%% The angle is applied once, by GetStandardDirMat; the tensor rotation then
            %%% takes angle 0 (passing the angle to both would rotate twice;
            %%% docs/pitfalls.md OTH1). The one place the solvers, the limiting
            %%% velocities and the plots get their standard-frame tensor from.
            if nargin < 2 || isempty(angleDeg)
                angleDeg = 0;
            end
            matStd   = obj.SetDirMat(Material.GetStandardDirMat(obj.dirMat, angleDeg));
            elTensor = matStd.GetRotElastTensor(0);
        end
        function rotDirMat = GetRotDirMat(obj,angleDeg)
            %%% Direction matrix rotated in-plane by the given angle (row 1, the surface
            %%% normal, is the pivot and stays unchanged)
            %%% >> rotDirMat = GetRotDirMat(obj,angleDeg)
            rotDirMat = Material.RotDirMatFromAngle(obj.dirMat,angleDeg);
        end
    end

    %% Static methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Static)
        dirMat = CompleteDirMat(inp)
        rotMat = CalcRotMat(dirMat,angleDeg)
        rotDirMat = RotDirMatFromAngle(dirMat,angleDeg)
        stdDirMat = GetStandardDirMat(dirMat,angleDeg)
        [phi1,Phi,phi2] = DirMat2Euler(dirMat)
        dirMat = Euler2DirMat(phi1,Phi,phi2)
        [phi1,Phi,phi2] = DCM2Euler(dcm)
        dcm = Euler2DCM(phi1,Phi,phi2)
        %% Miller index <-> direction cosine
        dcm = Miller2DirMat(hkl,symmNo,latticeParams,indexType)
        [hkl,residualDeg] = DirMat2Miller(dirCos,symmNo,latticeParams,indexType,maxIndex)
        %% Materials library (guliwerr/data/materialsLibrary.csv)
        [obj,entry] = FromLibrary(name,dirMat,varargin)
        lib = ListLibrary()
    end

    methods (Static, Access = private)
        spec = LatticeParamSpec(symmNo)
        spec = LatticeParamSpecReport(fname,symmNo)
        lib = ReadLibraryFile()
    end
end
