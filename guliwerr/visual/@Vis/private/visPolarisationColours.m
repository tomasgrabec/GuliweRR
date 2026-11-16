function baseColours = visPolarisationColours()
%%% Base colour of each displacement direction in the polarisation-RGB view
%%% >> baseColours = visPolarisationColours()
%%%
%%% @Vis-private helper shared by Vis.calcPolarisationRgb and
%%% Vis.plotPolarisationKey, so the panels and the key cannot drift apart.
%%%
%%% OUTPUT: baseColours .. [3 x 3], row i = the RGB colour of a displacement
%%%                        purely along direction i (1 = propagation, 2 =
%%%                        in-plane transverse, 3 = surface normal)
%%%
%%% Red and blue are the display primaries; green is darkened to 0.65.
%%% Pure (0,1,0) has about 3.4x the luminance of pure red and 10x that of
%%% blue, so on the white background a faint green region was barely visible
%%% and u3 read differently from u1/u2 at the same amplitude. 0.65 brings
%%% green's luminance (~0.27) close to red's (0.21). It also matches the
%%% UTGS paper's figures (H. Seiner), whose greens sit largely at G = 120-240
%%% on a 0-255 scale, not at the full primary.

    baseColours = [1, 0,    0;     % u1 -- propagation direction: red
                   0, 0,    1;     % u2 -- in-plane transverse: blue
                   0, 0.65, 0];    % u3 -- out-of-plane: darkened green
end
