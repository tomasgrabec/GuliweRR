function [rgbCell, ampCell] = calcPolarisationRgb(dispField)
%%% Merge a displacement-field triplet into one true-colour polarisation image
%%% >> [rgbCell, ampCell] = Vis.calcPolarisationRgb(dispField)
%%%
%%% INPUT: 1. dispField .. {3 x nLay} cell of signed real 2D fields
%%%                        [nMeshDepth x nMeshLateral], one per displacement
%%%                        direction (standard frame: 1 = propagation, 2 =
%%%                        in-plane transverse, 3 = surface normal), as from
%%%                        Vis.calcDispField2D / calcGreenDispField2D(Layered).
%%%                        The layers may differ in row count (cropped substrate)
%%% OUTPUT: rgbCell .. {1 x nLay}, each [nMeshDepth x nMeshLateral x 3] RGB
%%%                    image in [0,1]
%%%         ampCell .. {1 x nLay}, each [nMeshDepth x nMeshLateral] amplitude
%%%                    |u| / max|u| in [0,1], the same normaliser for every layer
%%%
%%% The colour scheme of Seiner's displacement-field figures in the UTGS paper
%%% (Grabec et al., rev2, Figs 3(g) and 4(b)):
%%%   - hue = the DIRECTION of the displacement at that point: the unit vector
%%%     p = (|u1|, |u2|, |u3|)/|u| weights one base colour per direction,
%%%     hue = p1*red + p2*blue + p3*green (visPolarisationColours: red and blue
%%%     are the primaries, green is darkened to 0.65 to balance its luminance).
%%%     Every direction maps to one point of the positive octant of the unit
%%%     sphere, which is what Vis.plotPolarisationKey draws;
%%%   - brightness = the amplitude |u|, white at zero, the saturated hue at the
%%%     field's global maximum: pixel = (1 - a) + a * hue.
%%% Signs are dropped (|u_i|), so each component shows two maxima per
%%% wavelength. The field is one snapshot (the real part at one phase), so two
%%% components that are 90 deg out of phase -- u1 and u3 of a Rayleigh wave --
%%% appear as separate stripes a quarter wavelength apart, not as one mixed
%%% colour.
%%%
%%% One normaliser across all layers (max of |u|, not the positive maximum
%%% used by the signed component panels), so nothing is clipped and layers stay
%%% comparable.

    validateattributes(dispField, {'cell'}, {'nrows', 3}, 'calcPolarisationRgb', 'dispField');
    nLay = size(dispField, 2);

    baseColours = visPolarisationColours();   % [3 x 3], row = direction's colour

    ampCell = cell(1, nLay);
    absCell = cell(1, nLay);
    ampMax  = 0;
    for ilay = 1:nLay
        absU = cat(3, abs(dispField{1, ilay}), abs(dispField{2, ilay}), abs(dispField{3, ilay}));
        amp  = sqrt(sum(absU.^2, 3));
        absCell{ilay} = absU;
        ampCell{ilay} = amp;
        ampMax = max(ampMax, max(amp(:)));
    end
    ampMax = max(ampMax, eps);

    rgbCell = cell(1, nLay);
    for ilay = 1:nLay
        amp = ampCell{ilay};
        p   = absCell{ilay} ./ max(amp, eps);   % unit direction (|u1|,|u2|,|u3|); 0 where amp = 0
        [nR, nC, ~] = size(p);
        hue = reshape(reshape(p, [], 3) * baseColours, nR, nC, 3);
        a   = amp ./ ampMax;
        rgbCell{ilay} = min(max(1 - a + a .* hue, 0), 1);
        ampCell{ilay} = a;
    end
end
