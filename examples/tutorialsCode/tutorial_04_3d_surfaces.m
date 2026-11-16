%% Tutorial 04 (code) -- directional bulk quantities: velocity/group-velocity/Young's surfaces and curves
% Narrated version: examples/tutorials/tutorial_04_3d_surfaces.m. Run guliwerr/init.m first.

%% Material -- Cu(001); no Structure or sweep needed
dirMat   = [1 0 0; 0 0 1];          % surface normal; propagation reference
material = Material.FromLibrary('Cu', dirMat);   % constants and density from the materials library

%% Mesh -- finer than the 61 x 121 default
nTheta = 91;                       % polar-angle points
nPhi   = 181;                      % azimuthal points

%% Bulk phase-velocity surfaces (qL, qT1, qT2)
surfaceType = 'velocity';          % 'velocity' | 'slowness'
Vis.plotVelSurface(material, surfaceType, nTheta, nPhi);
sgtitle('Cu (001): bulk phase velocity');

%% Bulk group-velocity (ray) surfaces -- folds are real physics
Vis.plotGroupVelSurface(material, nTheta, nPhi);
sgtitle('Cu (001): bulk group velocity');

%% Young's modulus surface
Vis.plotYoungsSurface(material, nTheta, nPhi);
title('Cu (001): Young''s modulus [GPa]');

%% 2D in-plane cuts at the dirMat surface normal
inPlaneAnglesDeg = 0:1:360;
Vis.plotVelCurve(material, inPlaneAnglesDeg);
title('Cu (001): bulk phase velocity vs. in-plane angle');
Vis.plotYoungsCurve(material, inPlaneAnglesDeg);
title('Cu (001): Young''s modulus vs. in-plane angle');
