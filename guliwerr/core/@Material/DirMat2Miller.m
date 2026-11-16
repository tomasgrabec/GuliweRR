function [hkl,residualDeg] = DirMat2Miller(dirCos,symmNo,latticeParams,indexType,maxIndex)
    %%% Find the smallest-integer Miller index whose direction best matches a direction cosine
    %%% >> [hkl,residualDeg] = Material.DirMat2Miller(dirCos,symmNo,latticeParams,indexType,maxIndex)
    %%% >> spec = Material.DirMat2Miller(symmNo)          -- list required latticeParams fields
    %%% >> spec = Material.DirMat2Miller([],symmNo)       -- same, explicit form
    %%%
    %%% Indexing has no exact general solution -- an arbitrary direction need not lie along
    %%% any low-index crystallographic direction. ALWAYS check residualDeg.
    %%%
    %%% INPUT:  1. dirCos        .. 1x3 (or 3x1) direction-cosine vector (need not be
    %%%                             pre-normalised)
    %%%         2. symmNo        .. symmetry class number (see SymmClassLibrary)
    %%%  (opt)  3. latticeParams .. same as Miller2DirMat (call with dirCos=[] to see requirements)
    %%%  (opt)  4. indexType     .. 'direction' or 'planeNormal' (see Miller2DirMat)
    %%%                             [default 'direction']
    %%%  (opt)  5. maxIndex      .. search |h|,|k|,|l| up to this integer [default 12]
    %%% OUTPUT: hkl         .. best-fit small-integer Miller index -- [h k l], or Miller-Bravais
    %%%                        [h k i l] (i=-(h+k)) for trigonal/hexagonal with
    %%%                        indexType 'planeNormal'
    %%%         residualDeg .. angular misfit [degrees] between hkl's direction and dirCos;
    %%%                        0 = exact match
    %%%
    %%% Brute-force over reduced integer triples scored against Miller2DirMat: exact for all
    %%% 7 crystal systems, but O(maxIndex^3) -- not for a hot loop.

    %%% requirements-listing mode
    if nargin < 2
        if nargin < 1 || ~isscalar(dirCos)
            error('Material:DirMat2Miller:needSymmNo', ...
                'Give a symmetry class number: Material.DirMat2Miller(symmNo) or DirMat2Miller([],symmNo).')
        end
        symmNo = dirCos;
        dirCos = [];
    end
    if isempty(dirCos)
        hkl = Material.LatticeParamSpecReport('DirMat2Miller',symmNo);
        residualDeg = [];
        return
    end

    if nargin < 3 || isempty(latticeParams); latticeParams = struct(); end
    if nargin < 4 || isempty(indexType);     indexType = 'direction';   end
    if nargin < 5 || isempty(maxIndex);      maxIndex = 12;             end

    dirCos = dirCos(:)' / norm(dirCos(:));

    bestDot = -inf;
    bestH = [1 0 0];
    for h = -maxIndex:maxIndex
        for k = -maxIndex:maxIndex
            for l = -maxIndex:maxIndex
                if h==0 && k==0 && l==0
                    continue
                end
                if gcd(gcd(abs(h),abs(k)),abs(l)) > 1
                    continue   % skip non-reduced triples (e.g. [2 2 0] duplicates [1 1 0])
                end
                cand = Material.Miller2DirMat([h k l],symmNo,latticeParams,indexType);
                d = dot(cand,dirCos);
                if d > bestDot
                    bestDot = d;
                    bestH = [h k l];
                end
            end
        end
    end

    residualDeg = acosd(min(bestDot,1));
    hkl = bestH;
    if ismember(symmNo,[5,6,7]) && strcmpi(indexType,'planeNormal')
        hkl = [hkl(1), hkl(2), -(hkl(1)+hkl(2)), hkl(3)];
    end
end
