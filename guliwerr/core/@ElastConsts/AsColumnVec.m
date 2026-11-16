function vecOut = AsColumnVec(vecIn,nExpected,fname)
    %%% Normalise an optional elastic-constant vector input to an nExpected x 1 column
    %%% >> vecOut = ElastConsts.AsColumnVec(vecIn,nExpected,fname)
    %%%
    %%% INPUT:  1. vecIn     .. nExpected-element vector (row or column), or [] / omitted
    %%%         2. nExpected .. required number of elements
    %%%         3. fname     .. calling function's name, for the error message
    %%% OUTPUT: vecOut .. nExpected x 1 column vector; all-zeros for empty input, which is
    %%%                   how the *Conversion functions are called when only their
    %%%                   transformation matrix (their second output) is wanted

    if isempty(vecIn)
        vecOut = zeros(nExpected,1);
        return
    end
    if ~isnumeric(vecIn) || ~isvector(vecIn) || numel(vecIn) ~= nExpected
        error('ElastConsts:AsColumnVec:wrongSize', ...
            'ElastConsts.%s: wrong input size -- expected a %d-element vector.', ...
            fname, nExpected)
    end
    vecOut = vecIn(:);
end
