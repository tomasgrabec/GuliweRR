function [elTensor, rho] = visStandardFrameTensor(mater)
%%% A Material's elastic tensor in the project's standard frame, plus density
%%% >> [elTensor, rho] = visStandardFrameTensor(mater)
%%%
%%% @Vis-private helper: the one incantation every Christoffel-consuming
%%% plotting function needs before it can call Christoffel.phaseVels /
%%% Christoffel.groupVelVectorsForDir / LimitingVels.findRoots, written once
%%% instead of once per caller.
%%%
%%% INPUT: 1. mater .. Material whose dirMat ALREADY carries whatever
%%%                    orientation the caller wants evaluated (a propagation
%%%                    angle rotation, if any, is the caller's business -- see
%%%                    Vis.christoffelLimVels' docblock for why that split
%%%                    matters)
%%% OUTPUT: elTensor .. rotated elastic tensor, standard frame
%%%         rho      .. density, passed straight through
%%%
%%% Standard frame (Material.GetStandardElastTensor: GetStandardDirMat ->
%%% Material.CalcRotMat -> ElastConsts.RotateTensor): axis 1 = dirMat's
%%% in-plane propagation reference (dirMat row 2), axis 2 = the remaining
%%% in-plane direction (right-hand completion, dirMat row 3), axis 3 = the
%%% SURFACE NORMAL (dirMat row 1) -- the same frame Ritz/Christoffel/
%%% GreenSolver use throughout the project.

    elTensor = mater.GetStandardElastTensor(0);
    rho      = mater.rho;
end
