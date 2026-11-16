%% Tutorial: 3D directional surfaces -- bulk waves and group velocity
%
% Direction-dependent bulk quantities of one anisotropic Material. Each plot
% sweeps direction itself, so only a Material is needed -- no Structure,
% CalcParams or sweep. Exploratory plots from the elastic tensor alone (no Ritz/Green):
% - Vis.plotVelSurface       -- bulk-wave (qL/qT1/qT2) phase velocity or
%                                slowness vs. wave-normal direction
% - Vis.plotGroupVelSurface  -- bulk-wave group velocity (the energy direction,
%                                which deviates from the wave normal); folds
%                                in strongly anisotropic materials are real
% - Vis.plotVelCurve / Vis.plotYoungsCurve -- the cheap 2D counterparts on the
%                                in-plane great circle of the dirMat surface
% - Vis.plotYoungsSurface    -- directional Young's modulus
%
% The SAW velocity has no such 3D surface: it exists only for a chosen cut
% and varies with the in-plane direction on it -- see tutorial_02
% (Vis.plotAngularMap) for one cut.
%
% Axes: the standard frame used by Ritz/Christoffel/GreenSolver (axis 3 =
% surface normal = dirMat row 1, axis 1 = dirMat row 2).
%
% Run run('guliwerr/init.m') from the project root first (tutorial_01 says why).

%% Material -- Cu(001), cubic (dirMat = Material's default, set explicitly)
% The Cu constants and density come from the materials library
% (Material.ListLibrary() shows values and sources).
dirMat = [1 0 0; ...         % surface normal
          0 0 1];            % propagation reference (Material's default)
mat    = Material.FromLibrary('Cu', dirMat);

%% Bulk-wave phase-velocity surfaces (qL/qT1/qT2), on a mesh finer than the
% 61 x 121 default
surfaceType = 'velocity';   % 'velocity' (default) | 'slowness' (= 1/velocity)
nTheta      = 91;           % polar-angle mesh points (default 61)
nPhi        = 181;          % azimuthal mesh points (default 121)
Vis.plotVelSurface(mat, surfaceType, nTheta, nPhi);
sgtitle('Cu (001) -- bulk phase-velocity surfaces');

%% Bulk-wave group-velocity ("ray") surfaces -- their shape can differ
% from the phase-velocity surfaces above
Vis.plotGroupVelSurface(mat, nTheta, nPhi);
sgtitle('Cu (001) -- bulk group-velocity surfaces');

%% 2D cut: bulk phase velocities vs. in-plane angle on the (001) surface --
% the bulk counterpart of the SAW angular map (tutorial_02), same angle axis
inPlaneAngles = 0:1:360;   % deg (default: full 1-deg circle)
Vis.plotVelCurve(mat, inPlaneAngles);
title('Cu (001) -- bulk phase velocities vs. in-plane angle');

%% Directional Young's modulus -- 3D surface and the same 2D in-plane cut
Vis.plotYoungsSurface(mat, nTheta, nPhi);
title('Cu (001) -- directional Young''s modulus [GPa]');

Vis.plotYoungsCurve(mat, inPlaneAngles);
title('Cu (001) -- Young''s modulus vs. in-plane angle');
