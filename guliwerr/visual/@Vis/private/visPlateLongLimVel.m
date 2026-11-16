function longLimVel = visPlateLongLimVel(structure, angleDeg)
%%% Longitudinal limiting velocity of a Plate structure, taken over BOTH free
%%% surfaces
%%% >> longLimVel = visPlateLongLimVel(structure, angleDeg)
%%%
%%% @Vis-private helper shared by Vis.plotRitzModeMapPlate and
%%% Vis.plotRitzModeCurvePlate, which both need a velocity-axis cutoff for a
%%% structure that (unlike the HS/Layered result classes) stores no
%%% longVelMap/longLimVel property of its own.
%%%
%%% INPUT: 1. structure .. Plate Structure
%%%        2. angleDeg  .. propagation angle [deg] to rotate to, or [] if the
%%%                        structure is already rotated -- see visPlateLimVels
%%% OUTPUT: longLimVel .. max over the bottom (layer 1) and, when the stack has
%%%                       more than one layer, the top (layer nLay) surface
%%%
%%% Max over both surfaces so a genuinely asymmetric stack's cutoff is never
%%% set by the slower surface alone (which would silently crop real modes
%%% belonging to the faster one).

    longLimVel = visPlateLimVels(structure, 1, angleDeg);
    if structure.nLay > 1
        longLimVel = max(longLimVel, visPlateLimVels(structure, structure.nLay, angleDeg));
    end
end
