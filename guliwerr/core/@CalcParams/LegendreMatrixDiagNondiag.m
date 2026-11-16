function [d,n] = LegendreMatrixDiagNondiag(legO)
    %%% Generate diagonal and non-diagonal Legendre integration matrices
    %%% >> [d,n] = LegendreMatrixDiagNondiag(legO)
    %%%
    %%% INPUT:  legO .. max order of the Legendre polynomials
    %%% OUTPUT: d    .. D_mn = integral( dP_m/dx * dP_n/dx )  (legO × legO, symmetric)
    %%%         n    .. N_mn = integral( dP_m/dx * P_n )      (legO × legO, lower-triangular)
    %%%   Polynomials are L2-normalised over [-1, 1]; both matrices are checkerboard-sparse
    %%%   (N nonzero only where m+n is odd, D only where m+n is even).

    %%% N_mn = sqrt((2m-1)(2n-1)) where m > n and m+n is odd, 0 elsewhere
    v = (1:legO);
    A = sqrt((2*v-1)' * (2*v-1)) / 2;
    [M, N_idx] = ndgrid(v, v);
    n = tril(A .* (1 - (-1).^(M + N_idx)), -1);

    %%% D = N*N': the N_mk are the coefficients of P_m' in the orthonormal basis {P_k},
    %%% so orthonormality collapses the D integral to sum_k N_mk N_nk
    d = n * n';
end
