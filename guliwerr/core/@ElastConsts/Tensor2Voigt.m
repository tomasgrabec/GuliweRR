function cij = Tensor2Voigt(inpTen)
    %%% Transform elastic tensor (3x3x3x3) to Voigt matrix (6x6)
    %%% >> cij = Tensor2Voigt(inpTen)

    cij= [inpTen(1,1,1,1) inpTen(1,1,2,2) inpTen(1,1,3,3) inpTen(1,1,2,3) inpTen(1,1,1,3) inpTen(1,1,1,2);
        inpTen(1,1,2,2) inpTen(2,2,2,2) inpTen(2,2,3,3) inpTen(2,2,2,3) inpTen(2,2,1,3) inpTen(2,2,1,2);
        inpTen(1,1,3,3) inpTen(2,2,3,3) inpTen(3,3,3,3) inpTen(3,3,2,3) inpTen(3,3,1,3) inpTen(3,3,1,2);
        inpTen(1,1,2,3) inpTen(2,2,2,3) inpTen(3,3,2,3) inpTen(2,3,2,3) inpTen(2,3,1,3) inpTen(2,3,1,2);
        inpTen(1,1,1,3) inpTen(2,2,1,3) inpTen(3,3,1,3) inpTen(2,3,1,3) inpTen(1,3,1,3) inpTen(1,3,1,2);
        inpTen(1,1,1,2) inpTen(2,2,1,2) inpTen(3,3,1,2) inpTen(2,3,1,2) inpTen(1,3,1,2) inpTen(1,2,1,2)];
end
