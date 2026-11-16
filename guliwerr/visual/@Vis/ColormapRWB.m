function cMap = ColormapRWB()
%%% Red-White-Blue diverging colormap, for signed data (e.g. displacement fields)
%%% >> cMap = Vis.ColormapRWB()
%%%
%%% OUTPUT: cMap .. 201x3 RGB colormap, red at low end, blue at high end, true
%%%                 white at the centre
%%%
%%% Formula courtesy of IA Veres. The centre entry is 1, giving a true white
%%% zero-crossing — the physically meaningful choice for signed modal-shape
%%% fields. A 0 there renders pure black instead.

    t = (0:0.01:0.99);
    cMap = [ones(1,length(t)) 1 fliplr(t); ...
            t                1 fliplr(t); ...
            t                1 ones(1,length(t))]';
end
