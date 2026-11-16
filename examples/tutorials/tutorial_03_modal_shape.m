%% Tutorial: modal-shape displacement fields -- Ritz mode + Green's response
%
% The SAW's displacement field (depth x lateral), two ways:
% - Vis.plotModalShape      -- the Ritz eigenmode (needs retentionLevel='full')
% - Vis.plotGreenModalShape -- the Green's-function response at the same
%                              velocity (half-space only)
% plus a |G| vs. velocity spectrum from a direct GreenSolver.calc call (own
% velocity axis and component), and the interactive viewer on the same result.
%
% Run run('guliwerr/init.m') from the project root first (tutorial_01 says why).

%% Material -- Cu(001), cubic (dirMat = Material's default, set explicitly)
% The Cu constants and density come from the materials library
% (Material.ListLibrary() shows values and sources).
dirMat = [1 0 0; ...         % surface normal
          0 0 1];            % propagation reference (Material's default)
mat    = Material.FromLibrary('Cu', dirMat);

%% Discretisation -- "easy" tier (Cu Zener ratio 2*C44/(C11-C12) = 3.2, principal (001) plane)
relDepth = 20;  % depth of the computational domain, in wavelengths
legO     = 50;  % Legendre basis order across that depth; keep >= 2*relDepth
cp = CalcParams(relDepth, legO);

strName   = 'cu_001_20deg';
nMaxModes = 500;   % max number of Ritz modes computed/retained
str = Structure(mat, cp, strName, nMaxModes);

%% One-angle 'full' sweep
angle          = 20;
calcMode       = 'full';   % Ritz+Christoffel+Green together
retentionLevel = 'full';   % plotModalShape needs the full eigenvector matrix
nVelGreen      = [];       % Green-axis points ([] = greenVelStep steps, 6 m/s)
sweep = Calc.angSweepHS(str, angle, calcMode, retentionLevel, nVelGreen);
sweep = sweep.calc();

if isempty(sweep.sawModeNoMap{1})
    error('tutorial_03_modal_shape:noSAWMode', ...
        'No SAW mode identified at %g deg for this material/angle.', angle);
end
modeNo = sweep.sawModeNoMap{1}(1);
vSAW   = sweep.sawVelMap{1}(1);

%% Ritz modal shape, true-scale depth
ritzSolver      = sweep.cs{1}.ritzSolver;
substrateCutoff = 'none';   % 'none' (true-scale depth) | 'decay' | 'fixedMultiple' |
                            % 'boosted' -- the last three crop or rescale a deep substrate under a film
cutoffParam     = [];       % decay fraction ('decay', default 0.05) or depth multiple
                            % ('fixedMultiple', default 4); ignored for 'none'/'boosted'
dirs            = 1:3;      % draw all three displacement directions
meshSpacing     = 0.01;     % depth/lateral mesh spacing passed to calcDispField2D
Vis.plotModalShape(ritzSolver, modeNo, substrateCutoff, cutoffParam, dirs, meshSpacing);
sgtitle(sprintf('Cu (001), \\theta = %g deg: Ritz %s mode, v = %.3f km/s', ...
    angle, sweep.sawTypeMap{1}{1}, vSAW));

%% The same mode as ONE polarisation-coloured panel (dirs = 'rgb')
% Colour = direction of the displacement (red u_x = propagation, blue u_y =
% in-plane transverse, green u_z = out-of-plane; key on the right), brightness
% = |u| (white where nothing moves). It is one snapshot: u_x and u_z of a
% Rayleigh wave are 90 deg out of phase, so they show as red and green stripes
% a quarter wavelength apart. 'rgb' also works for Vis.plotGreenModalShape /
% plotGreenModalShapeLayered and in the viewers (Display dropdown).
Vis.plotModalShape(ritzSolver, modeNo, substrateCutoff, cutoffParam, 'rgb', meshSpacing);
sgtitle(sprintf('Cu (001), \\theta = %g deg: polarisation view, v = %.3f km/s', angle, vSAW));

%% Green's-response field at that same SAW velocity
% A SCALAR velocity plus nargout>1 makes GreenSolver fill .dispField (half-space only).
[greenObj, ~] = GreenSolver.calc(str, angle, vSAW);
forceDir = 3;   % 1 = propagation, 2 = in-plane transverse, 3 = surface normal
Vis.plotGreenModalShape(greenObj, forceDir);
sgtitle(sprintf('Cu (001), \\theta = %g deg: Green response, v = %.3f km/s', ...
    angle, vSAW));

%% Green's-function spectrum |G| vs. velocity at this one angle, called
% directly to choose the velocity axis and component (here G_33). Peaks mark
% surface waves: a sharp pole = SAW, a sharp finite peak = pseudo-SAW.
longLimVel   = sweep.cs{1}.longLimVel;              % longitudinal limiting velocity, km/s
velMinGreen  = 1.5;                                 % km/s, lower end of the scan
velAxisGreen = linspace(velMinGreen, 1.05*longLimVel, 1000);
forceDirG    = 3;   % force direction I: 1 | 2 | 3 (x1 = propagation, x3 = surface normal)
displDirG    = 3;   % displacement direction J: 1 | 2 | 3 -- (1,3) gives G_13 (TGS)
greenSpec    = GreenSolver.calc(str, angle, velAxisGreen, forceDirG, displDirG);
absG         = abs(greenSpec.G);

figure();
semilogy(velAxisGreen, absG, 'LineWidth', 1.5);
hold on;
xline(vSAW, '--', sprintf('Ritz SAW %.3f km/s', vSAW));
xlabel('velocity [km/s]'); ylabel(sprintf('|G_{%d%d}| [a.u.]', forceDirG, displDirG));
title(sprintf('Cu (001), \\theta = %g deg: Green''s function vs. velocity', angle));

%% Interactive: the same result in the combined viewer -- map or curve
% (bottom) linked to the modal shape (top); click a point to see its shape.
% Needs an interactive MATLAB session (a -batch run has no display).
gui = Gui.combinedViewer(sweep);
gui.show();
gui.angSweepViewer.drillDownToAngle(1);   % open on the single-angle curve
% (same as the 'View single angle ->' button, top right; it toggles back to the map)
