function [longLimVel, shearLimVel] = visPlateLimVels(structure, layerIdx, angleDeg)
%%% Christoffel longitudinal + slowest-shear limiting velocities for one
%%% layer (surface) of a Plate structure
%%% >> [longLimVel, shearLimVel] = visPlateLimVels(structure, layerIdx, angleDeg)
%%%
%%% @Vis-private helper: the same LimitingVels.findRoots machinery
%%% CalcSingleHS.calcLim uses, applied to a single named layer. Shared by
%%% Vis.plotDispersionCurvesPlate (reference lines), Vis.plotRitzModeMapPlate and
%%% Vis.plotRitzModeCurvePlate (velocity cutoff).
%%%
%%% INPUT: 1. structure .. Plate Structure
%%%        2. layerIdx  .. which layer's material (1 = bottom surface,
%%%                        structure.nLay = top surface)
%%%        3. angleDeg  .. propagation angle [deg] to rotate the material to,
%%%                        or [] if the structure's dirMat ALREADY carries the
%%%                        propagation-angle rotation
%%% OUTPUT: longLimVel  .. longitudinal limiting velocity
%%%         shearLimVel .. slowest (Rayleigh-adjacent) shear limiting velocity
%%%
%%% The [] vs. angleDeg distinction is LOAD-BEARING, not a convenience:
%%% CalcWavevecAngSweep2Plate's and CalcAngSweepPlate's .structure are kept
%%% UNROTATED, so those callers must rotate here, whereas
%%% CalcSinglePlate.structure is already rotated for its single propagation angle
%%% and rotating again would double-apply it.

    mater = structure.mater(layerIdx);
    if ~isempty(angleDeg)
        mater = mater.SetDirMat(mater.GetRotDirMat(angleDeg));
    end

    [elTensor, rho] = visStandardFrameTensor(mater);
    [vLim_sA, ~]    = LimitingVels.findRoots(elTensor, rho);

    longLimVel  = vLim_sA{1}(1);
    shearLimVel = min([vLim_sA{2}(:); vLim_sA{3}(:)]);
end
