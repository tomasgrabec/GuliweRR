function cij = Vector2Voigt(inpVec,symmNo)
    %%% Transform vector of independent elast. consts into Voigt matrix
    %%% >> cij = Vector2Voigt(inpVec,symmNo)
    %%%
    %%% INPUT:  1. inpVec .. vector of independent elastic constants
    %%%  (opt)  2. symmNo .. symmetry class number; inferred if omitted or 0
    %%% OUTPUT: cij .. Voigt matrix (6x6)
    %%%
    %%% The one table of the symmetry classes: IndepVoigtPositions is derived from it.

    if nargin == 1 || symmNo == 0
        symmNo = ElastConsts.FindSymmClass_input(inpVec);
    end

    switch symmNo
        case 1 %iso
            c11 = inpVec(1);
            c44 = inpVec(2);
            c12 = c11-2*c44;
            cij=[c11 c12 c12  0   0   0
                c12 c11 c12  0   0   0
                c12 c12 c11  0   0   0
                0   0   0  c44  0   0
                0   0   0   0  c44  0
                0   0   0   0   0  c44];

        case  2 %cub
            c11 = inpVec(1);
            c12 = inpVec(2);
            c44 = inpVec(3);
            cij=[c11 c12 c12  0   0   0
                c12 c11 c12  0   0   0
                c12 c12 c11  0   0   0
                0   0   0 c44   0   0
                0   0   0   0  c44  0
                0   0   0   0   0  c44];

        case  3 %tetr_I
            c11 = inpVec(1);
            c12 = inpVec(2);
            c13 = inpVec(3);
            c33 = inpVec(4);
            c44 = inpVec(5);
            c66 = inpVec(6);
            cij=[c11 c12 c13  0   0   0
                c12 c11 c13  0   0   0
                c13 c13 c33  0   0   0
                0   0   0 c44   0   0
                0   0   0   0  c44  0
                0   0   0   0   0  c66];

        case  4 %tetr_II
            c11 = inpVec(1);
            c12 = inpVec(2);
            c13 = inpVec(3);
            c16 = inpVec(4);
            c33 = inpVec(5);
            c44 = inpVec(6);
            c66 = inpVec(7);
            cij=[c11  c12 c13  0   0   c16
                c12  c11 c13  0   0  -c16
                c13  c13 c33  0   0    0
                0    0   0 c44   0    0
                0    0   0   0  c44   0
                c16 -c16 0   0   0   c66];

        case  5 %trig_I
            c11 = inpVec(1);
            c12 = inpVec(2);
            c13 = inpVec(3);
            c14 = inpVec(4);
            c33 = inpVec(5);
            c44 = inpVec(6);
            c66 = (c11-c12)/2;
            cij=[c11 c12 c13  c14 0   0
                c12 c11 c13 -c14 0   0
                c13 c13 c33  0   0   0
                c14 -c14 0  c44  0   0
                0   0   0   0  c44 c14
                0   0   0   0  c14  c66];

        case  6 %trig_II
            c11 = inpVec(1);
            c12 = inpVec(2);
            c13 = inpVec(3);
            c14 = inpVec(4);
            c15 = inpVec(5);
            c33 = inpVec(6);
            c44 = inpVec(7);
            c66 = (c11-c12)/2;
            cij=[c11 c12 c13  c14  c15  0
                c12 c11 c13 -c14 -c15  0
                c13 c13 c33   0    0   0
                c14 -c14 0   c44   0  -c15
                c15 -c15 0    0   c44  c14
                0   0   0  -c15  c14  c66];

        case  7 %hex
            c11 = inpVec(1);
            c12 = inpVec(2);
            c13 = inpVec(3);
            c33 = inpVec(4);
            c44 = inpVec(5);
            c66 = (c11-c12)/2;
            cij=[c11 c12 c13  0   0   0
                c12 c11 c13  0   0   0
                c13 c13 c33  0   0   0
                0   0   0 c44   0   0
                0   0   0   0  c44  0
                0   0   0   0   0  c66];

        case 8 %orto
            c11=inpVec(1);
            c12=inpVec(2);
            c13=inpVec(3);
            c22=inpVec(4);
            c23=inpVec(5);
            c33=inpVec(6);
            c44=inpVec(7);
            c55=inpVec(8);
            c66=inpVec(9);
            cij=[c11 c12 c13  0  0  0
                c12 c22 c23  0  0  0
                c13 c23 c33  0  0  0
                0   0   0 c44  0  0
                0   0   0   0 c55 0
                0   0   0   0  0 c66];

        case 9 %mono
            c11=inpVec(1);
            c12=inpVec(2);
            c13=inpVec(3);
            c15=inpVec(4);
            c22=inpVec(5);
            c23=inpVec(6);
            c25=inpVec(7);
            c33=inpVec(8);
            c35=inpVec(9);
            c44=inpVec(10);
            c46=inpVec(11);
            c55=inpVec(12);
            c66=inpVec(13);
            cij=[c11  c12  c13  0  c15  0
                c12  c22  c23  0  c25  0
                c13  c23  c33  0  c35  0
                0    0     0  c44  0  c46
                c15  c25  c35  0  c55  0
                0    0     0  c46  0  c66];

        case 10 %tric
            c11=inpVec(1);
            c12=inpVec(2);
            c13=inpVec(3);
            c14=inpVec(4);
            c15=inpVec(5);
            c16=inpVec(6);
            c22=inpVec(7);
            c23=inpVec(8);
            c24=inpVec(9);
            c25=inpVec(10);
            c26=inpVec(11);
            c33=inpVec(12);
            c34=inpVec(13);
            c35=inpVec(14);
            c36=inpVec(15);
            c44=inpVec(16);
            c45=inpVec(17);
            c46=inpVec(18);
            c55=inpVec(19);
            c56=inpVec(20);
            c66=inpVec(21);
            cij=[c11 c12 c13  c14  c15  c16
                c12 c22 c23  c24  c25  c26
                c13 c23 c33  c34  c35  c36
                c14 c24 c34  c44  c45  c46
                c15 c25 c35  c45  c55  c56
                c16 c26 c36  c46  c56 c66];
        otherwise
            error('ElastConsts:Vector2Voigt:unknownSymm', ...
                'Unknown symmetry class number: %d (must be 1-10)', symmNo);
    end
end
