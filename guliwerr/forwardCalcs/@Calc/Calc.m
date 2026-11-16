classdef Calc
%%% Namespace for GuliweRR's forward-calculation orchestrators.
%%% Static methods are thin launchers over the Calc* handle classes in
%%% guliwerr/forwardCalcs/ -- same @ClassName/ pattern as
%%% Christoffel/GammaMatrix/GreenSolver.
%%% Do not instantiate this class.
%%%
%%% Usage:  sweep = Calc.angSweepHS(structure, angles, calcMode, retentionLevel, nVelGreen, useParallel)
%%%         sweep = sweep.calc();
%%%         path  = Calc.saveResult(sweep);   % view later: Gui.viewSaved(path)

    methods (Static)
        %%% Homogeneous half-space
        obj = singleHS(structure, angle, precompStruct, calcMode, retentionLevel, nVelGreen)
        obj = angSweepHS(structure, angles, calcMode, retentionLevel, nVelGreen, useParallel)

        %%% Single film on a substrate (layered)
        obj = singleLayered(structure, angle, vMax, precompStruct, calcMode, retentionLevel)
        obj = angSweepLayered(structure, angles, vMax, calcMode, retentionLevel, unused6, rawByDepthCell, useParallel)
        obj = wavevecAngSweep2Layered(structure, angles, wavelengths, vMax, calcMode, retentionLevel, unused7, useParallel)

        %%% Iso-frequency (fixed omega0, swept k-parallel) sweeps
        obj = isofreqSingleLayered(structure, angle, omega0, kParRange)
        obj = isofreqAngSweepLayered(structure, angles, omega0, kParRange, useParallel)
        obj = freqAngSweep2Layered(structure, angles, omega0List, kParRange, useParallel)

        %%% Free plate (Lamb waves)
        obj = singlePlate(structure, angle, realThickness, realWavVec, precompStruct, retentionLevel)
        obj = angSweepPlate(structure, angles, realThickness, realWavVec, retentionLevel, useParallel)
        obj = wavevecAngSweep2Plate(structure, angles, wavVecs, realThickness, retentionLevel, useParallel)

        %%% Mode classification and helpers
        obj    = modeClassification()
        omega0 = freqToOmega0(freqHz, lambdaRefNm)

        %%% Saving / loading a finished result (see CalcCommon.saveCalc/loadCalc)
        filePath           = saveResult(calcResult, filePath, note)
        [calcResult, meta] = loadResult(filePath)
    end
end
