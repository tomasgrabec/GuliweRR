function dispField = calcGreenDispField2D(greenObj, forceDir, meshSpacing)
%%% Signed 2D (depth x lateral) displacement field for one Green's-function
%%% response (one velocity, one force direction) -- the Green's-response
%%% analogue of calcDispField2D, feeding ModalShapeViewerGUI's Green's-response
%%% source.
%%% >> dispField = Vis.calcGreenDispField2D(greenObj, forceDir, meshSpacing)
%%%
%%% INPUT:  1. greenObj    .. GreenSolver object with .dispField populated: call
%%%                           GreenSolver.calc(structure, angle, velSelected, ...)
%%%                           with nargout>1 first. Full output is half-space
%%%                           only. velSelected should be a SCALAR velocity (a
%%%                           1-element velRange) — dispField is built over every
%%%                           point of whatever velRange is given, so a single
%%%                           value keeps the expensive full-partial-wave solve
%%%                           to the one velocity being displayed
%%%         2. forceDir    .. which force direction's response to show (1-3)
%%%  (opt)  3. meshSpacing .. lateral mesh spacing (default 0.01)
%%% OUTPUT: dispField .. {3 x 1} cell (half-space only, so nLay=1 always),
%%%                      dispField{iDir} = real signed 2D field
%%%                      [nMeshDepth x nMeshLateral], normalised by the global max
%%%                      across all 3 directions — the same shape and
%%%                      normalisation convention as calcDispField2D's {3 x nLay}
%%%                      output, so plotGreenModalShape mirrors plotModalShape.
%%%                      Row 1 = deepest, row end = the free surface, RE-ORIENTED
%%%                      from GreenSolver.calc's own storage (index 1 = surface)
%%%                      to match calcDispField2D's row convention
%%%
%%% No physical depth-axis reconstruction is attempted: GreenSolver.calc never
%%% stores its xx3 axis, and plotModalShape only uses a normalised row-index
%%% fraction for its y-axis anyway.

    if isempty(greenObj.dispField)
        error('Vis:calcGreenDispField2D:noFullOutput', ...
            ['greenObj.dispField is empty -- call GreenSolver.calc(...) with nargout>1 ' ...
             '(full output) first; full output is available for half-space structures only.']);
    end
    validateattributes(forceDir, {'double'}, {'scalar', 'integer', '>=', 1, '<=', 3}, ...
        'calcGreenDispField2D', 'forceDir');
    if nargin < 3 || isempty(meshSpacing)
        meshSpacing = 0.01;
    end

    depthProfile = squeeze(greenObj.dispField(:, :, forceDir, 1));   % [nxx3 x 3], complex
    depthProfile = flipud(depthProfile);                              % row 1 = deepest now

    %%% one "layer" only (half-space), so the shared helper's global
    %%% normalisation reduces to "across the 3 directions"
    dispField = visDispField2D({depthProfile}, meshSpacing);
end
