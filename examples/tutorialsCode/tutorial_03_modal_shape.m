%% Tutorial 03 (code) -- modal shapes: Ritz mode, Green's response, Green spectrum, viewer
% Narrated version: examples/tutorials/tutorial_03_modal_shape.m. Run guliwerr/init.m first.

%% Material -- Cu(001)
dirMat   = [1 0 0; 0 0 1];          % surface normal; propagation reference
material = Material.FromLibrary('Cu', dirMat);   % constants and density from the materials library

%% Discretisation + structure -- "easy" tier
relDepth      = 20;
legOrder      = 50;
calcParams    = CalcParams(relDepth, legOrder);
structureName = 'cu_001_20deg';
nMaxModes     = 500;
structure     = Structure(material, calcParams, structureName, nMaxModes);

%% Single-angle full sweep, keeping every eigenvector
angleDeg       = 20;
calcMode       = 'full';
retentionLevel = 'full';           % plotModalShape needs the full eigenvector matrix
nVelGreen      = [];               % Green-axis points ([] = 6 m/s steps)
sweep = Calc.angSweepHS(structure, angleDeg, calcMode, retentionLevel, nVelGreen);
sweep = sweep.calc();

modeNo = sweep.sawModeNoMap{1}(1); % Ritz index of the first classified SAW
vSAW   = sweep.sawVelMap{1}(1);    % km/s

%% Ritz modal shape
ritzSolver      = sweep.cs{1}.ritzSolver;
substrateCutoff = 'none';          % 'none' | 'decay' | 'fixedMultiple' | 'boosted'
cutoffParam     = [];              % decay fraction or depth multiple; unused for 'none'
dispDirs        = 1:3;             % displacement components to draw
meshSpacing     = 0.01;            % depth/lateral mesh spacing
Vis.plotModalShape(ritzSolver, modeNo, substrateCutoff, cutoffParam, dispDirs, meshSpacing);
sgtitle(sprintf('Cu (001), %g deg: Ritz mode, v = %.3f km/s', angleDeg, vSAW));

%% Same mode, one polarisation-RGB panel (red u_x, blue u_y, green u_z; brightness = |u|)
rgbDirs = 'rgb';                   % merged view instead of one panel per component
Vis.plotModalShape(ritzSolver, modeNo, substrateCutoff, cutoffParam, rgbDirs, meshSpacing);
sgtitle(sprintf('Cu (001), %g deg: polarisation view, v = %.3f km/s', angleDeg, vSAW));

%% Green's-response displacement field at the same velocity (half-space only)
[greenField, ~] = GreenSolver.calc(structure, angleDeg, vSAW);   % scalar velocity + 2 outputs fills .dispField
forceDir = 3;                      % 1 = propagation, 2 = in-plane transverse, 3 = normal
Vis.plotGreenModalShape(greenField, forceDir);
sgtitle(sprintf('Cu (001), %g deg: Green response, v = %.3f km/s', angleDeg, vSAW));

%% Green spectrum |G_IJ| vs. velocity, own axis and component
longLimVel   = sweep.cs{1}.longLimVel;                     % km/s
velMin       = 1.5;                                        % km/s
nVel         = 1000;
velAxisGreen = linspace(velMin, 1.05 * longLimVel, nVel);
forceDirI    = 3;                  % force direction I (1-3)
displDirJ    = 3;                  % displacement direction J (1-3); (1,3) = G_13 (TGS)
greenSpectrum = GreenSolver.calc(structure, angleDeg, velAxisGreen, forceDirI, displDirJ);
figure();
semilogy(velAxisGreen, abs(greenSpectrum.G), 'LineWidth', 1.5);
xline(vSAW, '--', 'Ritz SAW');
xlabel('velocity [km/s]'); ylabel(sprintf('|G_{%d%d}|', forceDirI, displDirJ));

%% Interactive viewer -- Green/Ritz curve linked to the modal shape
gui = Gui.combinedViewer(sweep);
gui.show();
gui.angSweepViewer.drillDownToAngle(1);   % open on the single-angle curve
% (same as the 'View single angle ->' button, top right; it toggles back to the map)
