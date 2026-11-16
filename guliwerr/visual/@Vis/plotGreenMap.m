function [hOut, hAx, hCb] = plotGreenMap(sweepResult, scaleMode, climRange, plotType, hAx, component)
%%% Plot the Green's-function angular map -- |G| heatmap, angle x velocity
%%% >> [hOut, hAx, hCb] = plotGreenMap(sweepResult, scaleMode, climRange, plotType, hAx, component)
%%%
%%% INPUT: 1. sweepResult .. CalcAngSweepHS or CalcAngSweepLayered result object,
%%%                          calc()'d with calcMode='full' — the Green maps stay
%%%                          empty under 'ritzOnly', which errors clearly here.
%%%                          .angles/.greenVelAxis and the component's map are
%%%                          read directly, both classes naming them the same
%%%  (opt) 2. scaleMode    .. 'log' (default) | 'linear'. 'log' plots
%%%                          log10(|G| / max|G|), normalised to the map's own
%%%                          peak so the display range is material- and
%%%                          units-independent and shows how sharp the SAW peak
%%%                          is relative to the background rather than its
%%%                          absolute magnitude
%%%  (opt) 3. climRange    .. 2-element [lo hi] colour-axis limits. For 'log' this
%%%                          is in log10(normalised) units, default [-4, -1], i.e.
%%%                          1e-4x to 1e-1x of the peak; auto for 'linear' unless
%%%                          given explicitly
%%%  (opt) 4. plotType     .. 'cartesian' (default, angle x velocity axes) |
%%%                          'polar' (angle as azimuth, velocity as radius,
%%%                          delegating to plotPolarHeatmap.m).
%%%                          sweepResult.angles is plotted exactly as given, not
%%%                          mirrored to a full circle — a 0-180 deg sweep draws
%%%                          a half-disc wedge
%%%  (opt) 5. hAx          .. target axes (created fresh if omitted). ALWAYS a
%%%                          plain axes(), never polaraxes(), for both plotTypes
%%%  (opt) 6. component    .. 'G13eff' (default, the stored |G13_eff| map -- the
%%%                          TGS response the existence score reads)
%%%                          | any 'GIJ' / 'GIJeff', e.g. 'G13' (one-sided G13_+),
%%%                          'G33': recomputed for this call (one Green solve per
%%%                          angle), not stored -- see Vis.greenComponentData
%%%                          | 'spot': the stored spot-averaged map (half-space result,
%%%                          sweepResult.calcSpotGreenMap(kernel) first), drawn on
%%%                          its own spotGreenVelAxis
%%% OUTPUT: hOut .. pcolor object handle
%%%         hAx  .. the axes used
%%%         hCb  .. colorbar handle
%%%
%%% Every component is drawn on the same angles x greenVelAxis grid at the same
%%% damping (SawSpectrum.ETA = 3e-3): the stored |G13_eff|, and any other
%%% component recomputed on request (Vis.greenComponentData).
%%% The axis starts at 0.3 km/s (SawSpectrum.VEL_MIN), so the lowest part is
%%% usually empty; set ylim, or use the viewers' range fields.
%%%
%%% The Cartesian branch draws each sample CENTRED on its (angle, velocity) value
%%% via Vis.pcolorCenteredGrid; plain pcolor(angles, vel, C) would draw angle(j)'s
%%% column from angle(j) to angle(j+1) instead — see docs/pitfalls.md.
%%%
%%% Default colormap ColMapGrayOrange. CData holds real, finite numbers (zeros
%%% floored to eps before log10); out-of-range values are clipped for DISPLAY
%%% only, via clim().

    if nargin < 6 || isempty(component); component = 'G13eff'; end
    visValidateGreenComponent(component, 'plotGreenMap');
    [absMap, componentLabel, velAxis] = Vis.greenComponentData(sweepResult, component);

    if nargin < 2 || isempty(scaleMode); scaleMode = 'log'; end
    visValidateChoice(scaleMode, {'log', 'linear'}, 'plotGreenMap', 'scaleMode');
    if nargin < 3 || isempty(climRange); climRange = [-4, -1]; end
    if nargin < 4 || isempty(plotType); plotType = 'cartesian'; end
    visValidateChoice(plotType, {'cartesian', 'polar'}, 'plotGreenMap', 'plotType');
    if nargin < 5 || isempty(hAx)
        figure();
        hAx = axes();
    end

    [plotData, cbLabel] = Vis.scaleGreenData(absMap, scaleMode, componentLabel);

    if strcmp(plotType, 'cartesian')
        hOut = visPcolorCentered(hAx, sweepResult.angles, velAxis, plotData);
        xlabel(hAx, 'angle [deg]');
        ylabel(hAx, 'velocity [km/s]');
    else
        hOut = Vis.plotPolarHeatmap(sweepResult.angles, velAxis, plotData, hAx);
    end

    hOut.DisplayName = ['|' componentLabel '|'];   % otherwise a generic 'dataN' if
                                                     % this axes' legend is later shown
    hCb = visGreenColourAxis(hAx, scaleMode, climRange, cbLabel);
end
