function A = CalcAnisotropyRatio(normElTen,angleStep)
    %%% Calculate generalised anisotropy ratio A
    %%%  as (fastest shear / slowest shear)^2 over all propagation directions
    %%% >> A = CalcAnisotropyRatio(normElTen,angleStep)
    %%%
    %%% INPUT: 1. normElTen  .. elastic tensor (3x3x3x3), density-normalised
    %%%  (opt) 2. angleStep  .. angular step in degrees (default = 4)
    %%% OUTPUT: generalised anisotropy ratio (>= 1; = 1 for isotropic)
    %%% NOTE: For cubic symmetry prefer GetAnisotropyRatio (exact Zener formula).

    if ~isequal(size(normElTen),[3 3 3 3])
        error('ElastConsts:CalcAnisotropyRatio:badInput', ...
            'normElTen must be a 3x3x3x3 tensor');
    end
    if nargin < 2 || isempty(angleStep)
        angleStep = 4;
    end
    if ~isnumeric(angleStep) || ~isscalar(angleStep) || angleStep <= 0
        error('ElastConsts:CalcAnisotropyRatio:badInput', ...
            'angleStep must be a positive scalar (degrees)');
    end

    alpha    = -45 : angleStep : 45;
    tanAlpha = tand(alpha);
    nAlpha   = length(alpha);

    phaseVels  = zeros(3,nAlpha,nAlpha,3);
    phaseVels2 = zeros(2,nAlpha,3*nAlpha);
    dirV       = zeros(3,3);

    %%% scan all directions — three orthogonal projections of the ±45° tangent grid
    for ii = 1:nAlpha
        for jj = 1:nAlpha
            dirV(1,:) = [tanAlpha(ii), tanAlpha(jj), 1           ];
            dirV(2,:) = [tanAlpha(ii),           1 , tanAlpha(jj)];
            dirV(3,:) = [          1 , tanAlpha(ii), tanAlpha(jj)];
            for idir = 1:3
                dirV(idir,:) = dirV(idir,:) ./ norm(dirV(idir,:));
                phaseVels(idir,ii,jj,:) = Christoffel.phaseVels(normElTen,dirV(idir,:));
            end
        end
    end

    %%% shear modes only (mode 1 is longitudinal)
    for iwt = 2:3
        phaseVels2(iwt-1,:,:) = [squeeze(phaseVels(1,:,:,iwt)) ...
                                  squeeze(phaseVels(2,:,:,iwt)) ...
                                  squeeze(phaseVels(3,:,:,iwt))];
    end

    %%% velocity^2 ratio -> elastic-stiffness units
    A = (max(phaseVels2,[],'all') / min(phaseVels2,[],'all'))^2;
end
