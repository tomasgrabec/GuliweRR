function [absG, componentLabel, velAxis] = greenComponentData(result, component)
%%% |G| data, TeX label and velocity axis of one Green component of a Calc result
%%% >> [absG, componentLabel] = Vis.greenComponentData(result)            % |G13_eff|
%%% >> [absG, componentLabel] = Vis.greenComponentData(result, 'G33')     % recomputed
%%% >> [absG, componentLabel, velAxis] = Vis.greenComponentData(result, 'spot')
%%%
%%% INPUT: 1. result    .. CalcSingleHS/CalcSingleLayered (a curve) or
%%%                        CalcAngSweepHS/CalcAngSweepLayered (a map) result,
%%%                        calc()'d with calcMode='full'
%%%  (opt) 2. component .. 'G13eff' (default) | 'GIJ' | 'GIJeff', I = force and
%%%                        J = displacement direction (1 = propagation, 3 = normal)
%%%                          'G13eff' -- |(G13_+ - G13_-)/2| = |(G13_+ - G31_+)/2|,
%%%                                      the +-k symmetrised TGS response: the
%%%                                      ONLY component a result stores, and the
%%%                                      one the existence score reads
%%%                          'GIJ'    -- one-sided |G_IJ(+k)|, e.g. 'G13' (G13_+),
%%%                                      'G33'; 'GIJeff' -- +-k symmetrised
%%%                        Anything but 'G13eff' is RECOMPUTED on the result's own
%%%                        axis and damping (greenComponent()/greenComponentMap():
%%%                        one Green solve per angle) and not stored.
%%%                          'spot'   -- the stored spot-averaged |G| (a half-space
%%%                                      result's calcSpotGreen()/calcSpotGreenMap()
%%%                                      run with a ProbeKernel), on its own axis
%%% OUTPUT: absG           .. [nVel x 1] (single result) or [nVel x nAngles] (sweep),
%%%                           on velAxis, single precision
%%%         componentLabel .. TeX label: 'G_{13}^{eff}', 'G_{13}^{+}', 'G_{33}',
%%%                           'G_{13}^{spot}', ...
%%%         velAxis        .. [nVel x 1] km/s: result.greenVelAxis, or
%%%                           result.spotGreenVelAxis for 'spot'
%%%
%%% Errors clearly when the data is missing: calcMode was not 'full'
%%% (Vis:greenComponentData:noGreenData), or the stored G13_eff is empty
%%% (Vis:greenComponentData:noG13EffData); 'spot' without the stored spot data
%%% (:noSpotGreenData) or on a layered result (:noSpotSupport).
%%%
%%% Public rather than @Vis-private because the viewers place their own overlays
%%% on these curves/maps and must read the same component the plot shows.

    if nargin < 2 || isempty(component); component = 'G13eff'; end
    visValidateGreenComponent(component, 'greenComponentData');
    component = char(component);

    if strcmp(component, 'spot')
        [absG, componentLabel, velAxis] = spotData(result);
        return
    end
    velAxis = result.greenVelAxis;

    isSweep = isprop(result, 'greenAbsG13EffMap');
    if isSweep
        effField  = 'greenAbsG13EffMap';
    else
        effField  = 'greenAbsG13Eff';
    end
    [I, J, isEff] = CalcCommon.parseGreenComponent(component);
    if I == J   % diagonal: G_+ = G_- (reciprocity), so no side to label
        componentLabel = sprintf('G_{%d%d}', I, J);
    else
        componentLabel = sprintf('G_{%d%d}^{%s}', I, J, ternarySide(isEff));
    end

    if isempty(result.greenVelAxis)
        error('Vis:greenComponentData:noGreenData', ...
            ['greenVelAxis is empty -- the result must be calculated with calcMode=''full'' ' ...
             'to populate Green''s function data.']);
    end
    if strcmp(component, 'G13eff')
        absG = result.(effField);
        if isempty(absG)
            error('Vis:greenComponentData:noG13EffData', ...
                ['%s is empty -- the result has no Green data (calcMode ''ritzOnly'' ' ...
                 'or not calculated).'], effField);
        end
    elseif isSweep
        absG = result.greenComponentMap(component);
    else
        absG = result.greenComponent(component);
    end
end

function [absG, componentLabel, velAxis] = spotData(result)
%%% The stored spot-averaged map/curve (calcSpotGreenMap / calcSpotGreen)
    if isprop(result, 'spotGreenMap')
        field = 'spotGreenMap'; how = 'sweepResult.calcSpotGreenMap(kernel)';
    elseif isprop(result, 'spotGreen')
        field = 'spotGreen'; how = 'singleResult.calcSpotGreen(kernel)';
    else
        error('Vis:greenComponentData:noSpotSupport', ...
            '''spot'' needs a half-space result (CalcAngSweepHS / CalcSingleHS); got %s.', class(result));
    end
    absG = result.(field);
    if isempty(absG)
        error('Vis:greenComponentData:noSpotGreenData', ...
            '%s is empty -- run %s first (kernel: a ProbeKernel).', field, how);
    end
    velAxis = result.spotGreenVelAxis;
    componentLabel = sprintf('G_{%s}^{spot}', result.spotGreenInfo.kernel.component(2:3));
end

function s = ternarySide(isEff)
    if isEff; s = 'eff'; else; s = '+'; end
end
