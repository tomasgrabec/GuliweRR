function [hSaw, hAx] = plotAngularMap(sweepResult, plotType, hAx)
%%% Plot SAW velocity vs. angle -- Cartesian (velocity) or polar (slowness)
%%% >> [hSaw, hAx] = plotAngularMap(sweepResult, plotType, hAx)
%%%
%%% INPUT: 1. sweepResult .. CalcAngSweepHS or CalcAngSweepLayered result object,
%%%                          already calc()'d -- .angles/.sawVelMap/.sawTypeMap
%%%                          read directly (populated under both 'full' and
%%%                          'ritzOnly' calcMode -- mode classification itself
%%%                          doesn't need Green)
%%%  (opt) 2. plotType    .. 'cartesian' (default, velocity vs angle) | 'polar'
%%%                          (SLOWNESS = 1/velocity vs angle, NOT velocity — the
%%%                          standard slowness-surface representation, and a real
%%%                          behavioural difference between the two modes, not
%%%                          just a coordinate-system swap)
%%%  (opt) 3. hAx         .. target axes, created fresh if omitted — axes() for
%%%                          cartesian, polaraxes() for polar. A caller passing
%%%                          one in must supply the matching TYPE
%%% OUTPUT: hSaw .. plot object handle for the SAW-velocity series: every SAW-type
%%%                 point across every angle, flattened from the ragged
%%%                 sawVelMap{iAngle} cell into one series
%%%         hAx  .. the axes used
%%%
%%% Each angle's sawVelMap entry may be empty or hold several velocities, so the
%%% flattening below never assumes one point per angle.
%%%
%%% Limiting-velocity reference lines are deliberately NOT drawn here —
%%% plotClassificationDiagnostic is where they belong.

    if nargin < 2 || isempty(plotType); plotType = 'cartesian'; end
    visValidateChoice(plotType, {'cartesian', 'polar'}, 'plotAngularMap', 'plotType');

    angles  = sweepResult.angles(:);
    nAngles = numel(angles);

    allAngles = [];
    allVels   = [];
    for ia = 1:nAngles
        v = sweepResult.sawVelMap{ia}(:);
        allAngles = [allAngles; repmat(angles(ia), numel(v), 1)]; %#ok<AGROW>
        allVels   = [allVels;   v];                                %#ok<AGROW>
    end

    if nargin < 3 || isempty(hAx)
        figure();
        if strcmp(plotType, 'cartesian')
            hAx = axes();
        else
            hAx = polaraxes();
        end
    end

    hold(hAx, 'on');
    if strcmp(plotType, 'cartesian')
        hSaw = plot(hAx, allAngles, allVels, '+');
        xlabel(hAx, 'angle [deg]');
        ylabel(hAx, 'SAW velocity [km/s]');
    else
        hSaw = polarplot(hAx, deg2rad(allAngles), 1 ./ allVels, '+');
    end
    hSaw.Color = 'r';
end
