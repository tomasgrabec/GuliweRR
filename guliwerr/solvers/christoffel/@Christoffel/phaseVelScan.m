function vels = phaseVelScan(mater, angleDeg)
%%% Phase velocities scanned over rotation angles for a Material object
%%% >> vels = Christoffel.phaseVelScan(mater, angleDeg)
%%%
%%% INPUT:  1. mater    .. Material object
%%%  (opt)  2. angleDeg .. rotation angles in degrees (default = 0)
%%% OUTPUT: vels(wavetype, angle) — [vL; vT1; vT2] at each angle

if nargin < 2 || isempty(angleDeg)
    angleDeg = 0;
end

if ~isa(mater, 'Material')
    error('Christoffel:phaseVelScan:badInput', 'mater must be a Material object.');
end

vels = zeros(3, length(angleDeg));

for iangle = 1:length(angleDeg)
    rotDirMat = mater.GetRotDirMat(angleDeg(iangle));
    dirVec    = rotDirMat(2,:);    % propagation direction for this angle
    vels(:, iangle) = Christoffel.phaseVels(mater.elConsts.tensor, dirVec, mater.rho);
end
end
