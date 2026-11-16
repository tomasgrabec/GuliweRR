function dispField = calcGreenDispField2DLayered(greenObj, structure, forceDir, meshSpacing)
%%% Signed 2D (depth x lateral) displacement field for a LAYERED Green's-
%%% function response (one velocity, one force direction) -- the layered
%%% counterpart of calcGreenDispField2D, built from GreenSolver.calc's
%%% k3OutLayered/polarizationULayered/partialAmplitudeLayered/hScaleOut
%%% instead of its half-space-only .dispField.
%%% >> dispField = Vis.calcGreenDispField2DLayered(greenObj, structure, forceDir, meshSpacing)
%%%
%%% INPUT:  1. greenObj    .. GreenSolver object from GreenSolver.calc(structure,
%%%                           angle, velSelected, ...) called with nargout>1 on
%%%                           a LAYERED structure (nLay>=2). velSelected should
%%%                           be a SCALAR velocity: as in calcGreenDispField2D,
%%%                           only the FIRST stored velocity point is read
%%%         2. structure   .. the SAME Structure passed to GreenSolver.calc,
%%%                           which keeps no copy of it. Only
%%%                           .calcParams(:).relDepth is read
%%%         3. forceDir    .. which force direction's response to show (1-3)
%%%  (opt)  4. meshSpacing .. depth/lateral mesh spacing (default 0.01)
%%% OUTPUT: dispField .. {3 x nLay} cell, dispField{iDir,iLay} = real signed
%%%                      2D field [nMeshDepth x nMeshLateral], Structure order
%%%                      (1 = substrate ... nLay = outermost film), normalised
%%%                      by ONE global max across every direction and layer --
%%%                      the same convention as calcDispField2D, so
%%%                      plotGreenModalShapeLayered can mirror plotModalShape.
%%%                      The substrate comes back at its FULL relDepth(1)
%%%                      extent; substrateCutoff cropping is the caller's job
%%%
%%% RECONSTRUCTION: every layer's solved amplitude vector is referenced to the
%%% SAME global x3 as the substrate (x3=0 at the substrate interface,
%%% x3=-H_total at the free surface; Every 2002 convention) -- buildBCMatrix
%%% never re-references a film's phase to its own local boundaries, only to
%%% whichever interface plane a boundary-condition row enforces. Per layer
%%% (nWaves = 6 film, 3 substrate):
%%%
%%%   u(x3,ci,cj) = sum_n A(n,cj) * U(ci,n) * exp(i*k3(n)*x3)
%%%
%%% ci = displacement component, cj = force direction. x3 maps linearly per
%%% layer onto RitzSolver's ti in [-1,1] (ti=-1 deep, ti=+1 shallow), the
%%% row-1-deepest ordering visDispField2D expects. Derivation:
%%% docs/reference/architecture/green.md, "Layered extension". At the free surface the
%%% reconstruction reproduces the solver's surface G_IJ, and a film of the
%%% substrate's own material gives the plain half-space field to round-off.
%%%
%%% SUBSTRATE DISPLAY EXTENT: GreenSolver's substrate is analytically
%%% semi-infinite, so structure.calcParams(1).relDepth plays no part in its
%%% physics. It is reused here purely as a finite display extent -- the same
%%% substrate truncation a companion RitzSolver on this Structure uses, so at
%%% matching meshSpacing both solvers share an identical ti grid per layer,
%%% index for index, and a direct Ritz-vs-Green shape comparison is meaningful.
%%%
%%% EVANESCENT-OVERFLOW LIMITATION: exp(i*k3*x3) for a single partial wave can
%%% overflow to Inf for a thick film, a deep substrate extent, or a layer far
%%% from x3=0, even though the true summed field is finite (Inf times a tiny
%%% amplitude is not the small residual it should be). Not fixed here -- see
%%% docs/reference/architecture/green.md, "Conditioning".

    if isempty(greenObj.k3OutLayered)
        error('Vis:calcGreenDispField2DLayered:noFullOutput', ...
            ['greenObj.k3OutLayered is empty -- call GreenSolver.calc(...) with nargout>1 ' ...
             'on a LAYERED structure (nLay>=2) first.']);
    end
    validateattributes(structure, {'Structure'}, {'scalar'}, 'calcGreenDispField2DLayered', 'structure');
    if structure.nLay ~= numel(greenObj.k3OutLayered)
        error('Vis:calcGreenDispField2DLayered:structureMismatch', ...
            ['structure.nLay (%d) does not match greenObj.k3OutLayered (%d layers) -- ' ...
             'pass the SAME Structure used to build greenObj.'], structure.nLay, numel(greenObj.k3OutLayered));
    end
    validateattributes(forceDir, {'double'}, {'scalar', 'integer', '>=', 1, '<=', 3}, ...
        'calcGreenDispField2DLayered', 'forceDir');
    if nargin < 4 || isempty(meshSpacing)
        meshSpacing = 0.01;
    end

    nLay   = structure.nLay;
    hScale = greenObj.hScaleOut;
    ti     = -1 : meshSpacing : 1;
    nMesh  = numel(ti);

    %%% film thickness per layer (index 1 = substrate, unused), in calc.m's
    %%% h = 2*pi*relDepth convention -- docs/reference/architecture/green.md, "relDepth
    %%% -> physical phase factor"
    filmThickness = zeros(nLay, 1);
    for vIdx = 2:nLay
        filmThickness(vIdx) = 2*pi * structure.calcParams(vIdx).relDepth;
    end

    depthProfiles = cell(1, nLay);
    for vIdx = 1:nLay
        k3 = greenObj.k3OutLayered{vIdx}(:, 1);             % [nWaves x 1], first (only) velocity
        U  = greenObj.polarizationULayered{vIdx}(:, :, 1);  % [3 x nWaves]
        A  = greenObj.partialAmplitudeLayered{vIdx}(:, :, 1); % [nWaves x 3]
        nWaves = numel(k3);

        if vIdx == 1
            %%% substrate (see SUBSTRATE DISPLAY EXTENT): ti=-1 (deep) -> x3 =
            %%% full extent, ti=+1 (shallow) -> x3 = 0, the first-film interface
            substrateExtent = 2*pi * structure.calcParams(1).relDepth * hScale;
            x3 = substrateExtent * (1 - ti) / 2;
        else
            %%% film: x3 of its lower (deep, ti=-1) and upper (shallow, ti=+1)
            %%% boundaries, offset by the films strictly BELOW it (index 2 = the
            %%% bottom film, touching the substrate)
            dBelow = sum(filmThickness(2:vIdx-1)) * hScale;   % 0 when vIdx==2
            dThis  = filmThickness(vIdx) * hScale;
            x3     = -dBelow - (ti + 1)/2 * dThis;
        end

        profile = complex(zeros(nMesh, 3));
        for ci = 1:3
            for n = 1:nWaves
                profile(:, ci) = profile(:, ci) + A(n, forceDir) * U(ci, n) .* exp(1i * k3(n) .* x3(:));
            end
        end
        depthProfiles{vIdx} = profile;
    end

    dispField = visDispField2D(depthProfiles, meshSpacing);
end
