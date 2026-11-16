function [symmC,nConsts] = SymmClassLibrary(symmNo)
    %%% Library of symmetry classes and number of independent constants
    %%% >> [symmC,nConsts] = SymmClassLibrary(symmNo)
    %%%
    %%% INPUT:  symmNo .. symmetry class number; omit to print and return all ten
    %%% OUTPUT: 1. symmC   .. symmetry class name (or all names, for no input)
    %%%         2. nConsts .. number of independent constants (or all counts)

    classNames = { ...
        'isotropic: [c11,c44]'; ...
        'cubic: [c11,c12,c44]'; ...
        'tetragonal, type I: [c11,c12,c13,c33,c44,c66]'; ...
        'tetragonal, type II  (c_16 ~= 0): [c11,c12c,c13,c16,c33,c44,c66]'; ...
        'trigonal, type I: [c11,c12,c13,c14,c33,c44]'; ...
        'trigonal, type II  (c_15 ~= 0): [c11,c12,c13,c14,c15,c33,c44]'; ...
        'hexagonal: [c11,c12,c13,c33,c44]'; ...
        'orthogonal: [c11,c12,c13,c22,c23,c33,c44,c55,c66]'; ...
        'monoclinic: [c11,c12,c13,c15,c22,c23,c25,c33,c35,c44,c46,c55,c66]'; ...
        'triclinic: [c11,c12,...,c16,c22,c23,...,c26,c33,...]'};
    constCounts = [2,3,6,7,6,7,5,9,13,21];

    if nargin == 0
        symmC   = classNames;
        nConsts = constCounts;

        scLines = cell(11,1);
        scLines{1} = ['Class No.','   ','No.el.consts.','   ','Class name'];
        for ii = 1:10
            scLines{ii+1} = [' ', num2str(ii,'%02d'), '         ', ...
                             num2str(nConsts(ii),'%02d'), '               ', symmC{ii}];
        end
        disp(strjoin(scLines, newline))

    elseif isnumeric(symmNo) && isscalar(symmNo) && any(symmNo == 1:10)
        symmC   = classNames{symmNo};
        nConsts = constCounts(symmNo);
    else
        error('ElastConsts:SymmClassLibrary:unknownSymm', ...
            'Unknown symmetry class number (must be 1-10).');
    end
end
