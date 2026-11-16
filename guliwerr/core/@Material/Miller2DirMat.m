function dcm = Miller2DirMat(hkl,symmNo,latticeParams,indexType)
    %%% Convert a Miller (or Miller-Bravais) index to a direction-cosine row vector
    %%% >> dcm = Material.Miller2DirMat(hkl,symmNo,latticeParams,indexType)
    %%% >> spec = Material.Miller2DirMat(symmNo)          -- list required latticeParams fields
    %%% >> spec = Material.Miller2DirMat([],symmNo)       -- same, explicit form
    %%%
    %%% INPUT:  1. hkl           .. Miller index [h k l]; OR Miller-Bravais [h k i l]
    %%%                             (i = -(h+k); trigonal/hexagonal PLANE NORMALS only)
    %%%         2. symmNo        .. symmetry class number (see SymmClassLibrary)
    %%%  (opt)  3. latticeParams .. struct of lattice-shape parameters -- ratios and angles
    %%%                             only, absolute scale never affects a direction cosine.
    %%%                             Call with hkl=[] (or just symmNo) to list the fields
    %%%                             symmNo requires; unrequired fields are ignored:
    %%%                               b_over_a, c_over_a .. axis-length ratios
    %%%                               alpha,beta,gamma   .. interaxial angles [degrees]
    %%%                                 (alpha = angle(b,c), beta = angle(a,c), gamma = angle(a,b))
    %%%  (opt)  4. indexType     .. 'direction' (real-space [uvw], DEFAULT) or 'planeNormal'
    %%%                             (reciprocal-space (hkl)). NOT parallel except for cubic
    %%%                             symmetry -- a surface normal is a planeNormal, a
    %%%                             propagation direction is a direction.
    %%% OUTPUT: dcm .. 1x3 unit direction-cosine row vector in the Cartesian crystallophysical
    %%%                frame (a along X, b in the XY-plane, c completing the cell); OR, in
    %%%                requirements-listing mode, the LatticeParamSpec struct

    %%% requirements-listing mode
    if nargin < 2
        if nargin < 1 || ~isscalar(hkl)
            error('Material:Miller2DirMat:needSymmNo', ...
                'Give a symmetry class number: Material.Miller2DirMat(symmNo) or Miller2DirMat([],symmNo).')
        end
        symmNo = hkl;
        hkl    = [];
    end
    if isempty(hkl)
        dcm = Material.LatticeParamSpecReport('Miller2DirMat',symmNo);
        return
    end

    if nargin < 3 || isempty(latticeParams)
        latticeParams = struct();
    end
    if nargin < 4 || isempty(indexType)
        indexType = 'direction';
    end

    spec = Material.LatticeParamSpec(symmNo);
    for k = 1:numel(spec.required)
        f = spec.required{k};
        if ~isfield(latticeParams,f) || isempty(latticeParams.(f))
            error('Material:Miller2DirMat:missingParam', ...
                '%s requires latticeParams.%s (missing). Full requirement list: %s', ...
                spec.className, f, strjoin(spec.required,', '))
        end
    end
    p = spec.defaults;
    for k = 1:numel(spec.required)
        f = spec.required{k};
        p.(f) = latticeParams.(f);
    end

    %%% Miller-Bravais 4-index [h k i l] -> 3-index [h k l]; plane normals only, since the
    %%% real-space zone-axis form [UVTW] is a different transform and is not implemented
    if numel(hkl) == 4
        if strcmpi(indexType,'direction')
            error('Material:Miller2DirMat:zoneAxisNotSupported', ...
                ['4-index Miller-Bravais [h k i l] is a PLANE-NORMAL notation; the equivalent ' ...
                 '4-index zone-axis notation for a DIRECTION ([UVTW]) uses a different, ' ...
                 'non-trivial transform not implemented here. Supply the 3-index direction ' ...
                 '[u v w] instead.'])
        end
        if abs(hkl(3) + hkl(1) + hkl(2)) > 1e-9
            error('Material:Miller2DirMat:badMillerBravais', ...
                'Miller-Bravais index must satisfy i = -(h+k); got h=%g, k=%g, i=%g.', ...
                hkl(1), hkl(2), hkl(3))
        end
        hkl = hkl([1,2,4]);
    end
    hkl = hkl(:)';

    %%% direct lattice vectors; a fixed to 1, as absolute scale cannot affect a direction cosine
    a = 1; b = p.b_over_a; c = p.c_over_a;
    al = deg2rad(p.alpha); be = deg2rad(p.beta); ga = deg2rad(p.gamma);

    aVec = [a 0 0];
    bVec = [b*cos(ga), b*sin(ga), 0];
    cx = c*cos(be);
    cy = c*(cos(al) - cos(be)*cos(ga)) / sin(ga);
    cz2 = c^2 - cx^2 - cy^2;
    if cz2 < -1e-9
        error('Material:Miller2DirMat:badLatticeParams', ...
            'latticeParams do not describe a valid unit cell (negative cell volume^2).')
    end
    cVec = [cx, cy, sqrt(max(cz2,0))];

    switch lower(indexType)
        case 'direction'
            vCart = hkl(1)*aVec + hkl(2)*bVec + hkl(3)*cVec;
        case 'planenormal'
            V = dot(aVec, cross(bVec,cVec));
            aStar = cross(bVec,cVec) / V;
            bStar = cross(cVec,aVec) / V;
            cStar = cross(aVec,bVec) / V;
            vCart = hkl(1)*aStar + hkl(2)*bStar + hkl(3)*cStar;
        otherwise
            error('Material:Miller2DirMat:badIndexType', ...
                'indexType must be ''direction'' or ''planeNormal''.')
    end

    nrm = norm(vCart);
    if nrm < 1e-12
        error('Material:Miller2DirMat:zeroVector', 'hkl maps to the zero vector.')
    end
    dcm = vCart / nrm;
end
