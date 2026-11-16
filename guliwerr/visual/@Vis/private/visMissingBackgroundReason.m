function msg = visMissingBackgroundReason(component, isSweep)
%%% Why a classification plot's response background is missing, per component
%%% >> msg = visMissingBackgroundReason(component, isSweep)
%%%
%%% @Vis-private helper of Vis.plotClassificationDiagnostic / plotClassificationCurve
%%% (their noGreenBackground warning): each background has its own source.
%%%
%%% INPUT:  1. component .. 'spot' | 'ritzSpectrum' | a Green component ('G13eff', ...)
%%%         2. isSweep   .. true for the diagnostic (a sweep), false for the curve
%%% OUTPUT: msg .. e.g. "No G13eff map (calcMode was not 'full')"
    if isSweep; noun = 'map'; else; noun = 'curve'; end
    switch component
        case 'spot'
            if isSweep; how = 'calcSpotGreenMap'; else; how = 'calcSpotGreen'; end
            msg = sprintf('No spot-averaged Green %s (run %s first)', noun, how);
        case 'ritzSpectrum'
            msg = 'No Ritz spectrum (this result stores no Ritz modal numerators)';
        otherwise
            msg = sprintf('No %s %s (calcMode was not ''full'')', component, noun);
    end
end
