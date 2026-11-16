function rotTen = RotateTensor(inpTen,rotMat)
    %%% Rotate tensor based on the rotation matrix
    %%% >> rotTen = RotateTensor(inpTen,rotMat)
    %%% INPUT:  1. tensor 3x3x3x3
    %%%         2. rotation matrix (2x3/3x2, or 3x3 row vectors)
    %%% OUTPUT: rotated tensor 3x3x3x3

    switch numel(rotMat)
        case 6      % 2x3 or 3x2: two given vectors, third completed right-hand
            R = ElastConsts.RightHandFrame(rotMat, ...
                'Rotation matrix recalculated based on its first vector (input vectors not perpendicular).');

        case 9      % 3x3: used as given, once checked orthogonal
            if all(size(rotMat) == 3)
                R = rotMat;
                if any(abs(R'*R - eye(3)) > 1e-6, 'all')
                    error('ElastConsts:RotateTensor:notOrthogonal', ...
                        'Rotation matrix not orthogonal.')
                end
                for jj = 1:3
                    R(jj,:)=R(jj,:)./norm(R(jj,:));
                end
            else
                error('ElastConsts:RotateTensor:badFormat', ...
                'Unknown rotation matrix format (set a 2x3 or 3x3 matrix).');
            end
        otherwise
            error('ElastConsts:RotateTensor:badFormat', ...
                'Unknown rotation matrix format (set a 2x3 or 3x3 matrix).');
    end

    %%% C'_{i1i2i3i4} = sum_{j1j2j3j4} R_{i1j1} R_{i2j2} R_{i3j3} R_{i4j4} C_{j1j2j3j4},
    %%% evaluated as four sequential mode-n matrix products
    rotTen = inpTen;
    for d = 1:4
        perm        = [d, 1:d-1, d+1:4];   % bring dimension d to front
        iperm       = zeros(1,4);
        iperm(perm) = 1:4;
        rotTen = permute(rotTen, perm);
        rotTen = reshape(R * reshape(rotTen, 3, []), [3 3 3 3]);
        rotTen = permute(rotTen, iperm);
    end
end
