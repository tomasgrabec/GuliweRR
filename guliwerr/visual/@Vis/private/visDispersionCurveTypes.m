function validTypes = visDispersionCurveTypes()
%%% The six curveType names both dispersion-curve plotters accept
%%% >> validTypes = visDispersionCurveTypes()
%%%
%%% @Vis-private helper: Vis.plotDispersionCurves (Layered) and
%%% Vis.plotDispersionCurvesPlate stay separate functions but offer the SAME six
%%% curve types, for familiarity. One list keeps that promise literal.

    validTypes = {'wavelength-freq', 'wavelength-vphase', 'k-freq', ...
                  'k-vphase', 'freq-vphase', 'freq-vgroup'};
end
