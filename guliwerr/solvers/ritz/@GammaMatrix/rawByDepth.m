function rawTensor = rawByDepth(calcParams)
%%% Step 1: Legendre-integral tensor (bVL x bVL x 15 x 3) — material-independent
%%% >> rawTensor = GammaMatrix.rawByDepth(calcParams)
%%% Ref: Stoklasova 2015, Eq. (A1)-(A4); Grabec 2024
%%%
%%% dim 3: index over 15 reduced elastic constants
%%% dim 4: depth-factor group:
%%%   (:,:,:,1) = dK0/dc  — factor 1/dd1^2  (= 1 since dd1=1)
%%%   (:,:,:,2) = dK1/dc  — factor 1/(dd1*dd2/2)
%%%   (:,:,:,3) = dK2/dc  — factor 1/(dd2/2)^2

    legO  = calcParams.legOrder;
    bVL   = calcParams.basisVecL;   % = 3*legO
    legD  = calcParams.legDiag;     % D(m,n) = integral(dPm/dx * dPn/dx)
    legN  = calcParams.legNondiag;  % N(m,n) = integral(dPm/dx * Pn)
    delta = eye(legO);

    [indRedEl, nRedEl] = ElastConsts.GetIndOf2DReduct_15of21();
    rawTensor = zeros(bVL, bVL, nRedEl, 3);

    for iconst = 1:nRedEl
        consts           = zeros(1, 21);
        consts(indRedEl(iconst)) = 1;
        cTen = ElastConsts.Any2Tensor(consts, 10);   % triclinic dummy (symmNo=10)

        triu_K0  = zeros(bVL, bVL);
        triu_K1a = zeros(bVL, bVL);
        triu_K1b = zeros(bVL, bVL);
        triu_K2  = zeros(bVL, bVL);

        for ind1 = 1:bVL
            iL1 = mod(ind1-1, legO);      % 0-indexed Legendre order of row
            iD1 = ceil(ind1 / legO);      % direction index 1-3
            for ind2 = ind1:bVL           % upper triangle only
                iL2 = mod(ind2-1, legO);
                iD2 = ceil(ind2 / legO);

                triu_K0 (ind1,ind2) = cTen(iD1,1,iD2,1) * (2*pi)^2 * delta(iL1+1, iL2+1);
                triu_K1a(ind1,ind2) = cTen(iD1,1,iD2,3) * ( 2*pi)  * legN(iL2+1, iL1+1);
                triu_K1b(ind1,ind2) = cTen(iD1,3,iD2,1) * (-2*pi)  * legN(iL1+1, iL2+1);
                triu_K2 (ind1,ind2) = cTen(iD1,3,iD2,3)             * legD(iL1+1, iL2+1);
            end
        end

        triu_K1 = triu_K1a + triu_K1b;

        %%% symmetrise: K0 and K2 are real-symmetric; K1 is Hermitian (1i x antisymmetric)
        rawTensor(:,:,iconst,1) = triu_K0  + triu(triu_K0, 1).';
        rawTensor(:,:,iconst,2) = 1i * (triu_K1 - triu(triu_K1, 1).');
        rawTensor(:,:,iconst,3) = triu_K2  + triu(triu_K2, 1).';
    end
end
