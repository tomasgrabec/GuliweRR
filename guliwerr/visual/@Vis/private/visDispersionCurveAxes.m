function [xAll, yAll, xLabelStr, yLabelStr] = visDispersionCurveAxes(curveType, data)
%%% Map a curveType onto the x/y data and axis labels of a dispersion-curve plot
%%% >> [xAll, yAll, xLabelStr, yLabelStr] = visDispersionCurveAxes(curveType, data)
%%%
%%% @Vis-private helper shared by Vis.plotDispersionCurves and
%%% Vis.plotDispersionCurvesPlate. The two derive their quantities very
%%% differently -- the Layered side sweeps real wavelengths and computes
%%% freq = vel/wavelength, the Plate side sweeps real wavevectors and reads
%%% freqsPhysicalMap3D straight out of RitzSolverLamb -- so each caller does
%%% its OWN physics and hands the finished arrays here. Only the six-way
%%% "which of these goes on which axis" mapping is shared, since that part
%%% genuinely is the same decision in both.
%%%
%%% INPUT: 1. curveType .. one of visDispersionCurveTypes()
%%%        2. data      .. struct with fields:
%%%             .wavelength      [nWav x 1]
%%%             .k               [nWav x 1]
%%%             .vel, .velGroup, .freq   [nWav x nCompModes]
%%%             .wavelengthLabel, .kLabel, .freqLabel  -- axis-label strings,
%%%                per-caller because the same axis means something slightly
%%%                different on each side (the Layered side's frequency is
%%%                derived as v_phi/wavelength and says so; the Plate side's
%%%                is a physical frequency in the units realThickness/wavVecs
%%%                were given in, and its 'k' IS the swept wavevector rather
%%%                than 2*pi/wavelength)
%%% OUTPUT: xAll, yAll .. [nWav x nCompModes] plot data
%%%         xLabelStr, yLabelStr .. axis labels
%%%
%%% Velocity labels are NOT caller-supplied: phase and group velocity mean
%%% exactly the same thing, in the same units, on both sides.

    nCompModes = size(data.vel, 2);
    VPHASE_LABEL = 'v_\phi [km/s]';
    VGROUP_LABEL = 'v_g [km/s]';

    switch curveType
        case 'wavelength-freq'
            xAll = repmat(data.wavelength, 1, nCompModes); yAll = data.freq;
            xLabelStr = data.wavelengthLabel; yLabelStr = data.freqLabel;
        case 'wavelength-vphase'
            xAll = repmat(data.wavelength, 1, nCompModes); yAll = data.vel;
            xLabelStr = data.wavelengthLabel; yLabelStr = VPHASE_LABEL;
        case 'k-freq'
            xAll = repmat(data.k, 1, nCompModes); yAll = data.freq;
            xLabelStr = data.kLabel; yLabelStr = data.freqLabel;
        case 'k-vphase'
            xAll = repmat(data.k, 1, nCompModes); yAll = data.vel;
            xLabelStr = data.kLabel; yLabelStr = VPHASE_LABEL;
        case 'freq-vphase'
            xAll = data.freq; yAll = data.vel;
            xLabelStr = data.freqLabel; yLabelStr = VPHASE_LABEL;
        case 'freq-vgroup'
            xAll = data.freq; yAll = data.velGroup;
            xLabelStr = data.freqLabel; yLabelStr = VGROUP_LABEL;
    end
end
