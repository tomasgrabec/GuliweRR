function cutoffRelDepth = calcSubstrateCutoffDepth(substrateCutoff, cutoffParam, relDepths, substrateDispMag)
%%% Effective substrate relDepth to DISPLAY for plotModalShape's axis layout
%%% >> cutoffRelDepth = Vis.calcSubstrateCutoffDepth(substrateCutoff, cutoffParam, relDepths, substrateDispMag)
%%%
%%% INPUT: 1. substrateCutoff  .. 'none' | 'decay' | 'fixedMultiple'
%%%        2. cutoffParam      .. threshold fraction (decay, e.g. 0.05) or multiple
%%%                               (fixedMultiple, e.g. 4); ignored for 'none'
%%%        3. relDepths        .. [nLay x 1] every layer's relDepth, layer 1 = substrate,
%%%                               layer nLay = outermost film (nearest the free surface)
%%%  (opt) 4. substrateDispMag .. required only for 'decay': [nMeshDepth x 1]
%%%                               displacement magnitude across the SUBSTRATE's own
%%%                               depth mesh, index 1 = deepest (ti=-1), index end =
%%%                               the substrate/first-film interface (ti=+1) -- or the
%%%                               free surface itself if nLay==1
%%% OUTPUT: cutoffRelDepth .. effective substrate relDepth to use for axis-height
%%%                          allocation (<= relDepths(1))

    if nargin < 2; cutoffParam = []; end
    validateattributes(substrateCutoff, {'char', 'string'}, {'scalartext'}, ...
        'calcSubstrateCutoffDepth', 'substrateCutoff');
    validateattributes(relDepths, {'double'}, {'vector', 'positive'}, ...
        'calcSubstrateCutoffDepth', 'relDepths');
    nLay = numel(relDepths);

    switch substrateCutoff
        case 'none'
            cutoffRelDepth = relDepths(1);

        case 'fixedMultiple'
            if nLay < 2
                error('Vis:calcSubstrateCutoffDepth:noFilmLayers', ...
                    '''fixedMultiple'' needs at least one film layer to reference -- use ''none'' for a plain half-space.');
            end
            if isempty(cutoffParam); cutoffParam = 4; end
            validateattributes(cutoffParam, {'double'}, {'scalar', 'positive'}, ...
                'calcSubstrateCutoffDepth', 'cutoffParam');
            filmStackDepth = sum(relDepths(2:end));
            cutoffRelDepth = min(relDepths(1), cutoffParam * filmStackDepth);

        case 'decay'
            if nLay < 2
                error('Vis:calcSubstrateCutoffDepth:noFilmLayers', ...
                    '''decay'' needs at least one film layer to reference -- use ''none'' for a plain half-space.');
            end
            if nargin < 4 || isempty(substrateDispMag)
                error('Vis:calcSubstrateCutoffDepth:missingDispMag', ...
                    '''decay'' requires substrateDispMag (the computed mode''s displacement magnitude across the substrate depth mesh).');
            end
            if isempty(cutoffParam); cutoffParam = 0.05; end
            validateattributes(cutoffParam, {'double'}, {'scalar', '>', 0, '<', 1}, ...
                'calcSubstrateCutoffDepth', 'cutoffParam');
            validateattributes(substrateDispMag, {'double'}, {'vector', 'nonnegative'}, ...
                'calcSubstrateCutoffDepth', 'substrateDispMag');

            %%% scan from the interface (index end, ti=+1) down toward the
            %%% deepest point (index 1, ti=-1)
            nMesh     = numel(substrateDispMag);
            cutoffIdx = visDecayCutoffIdx(substrateDispMag, cutoffParam);

            ti = linspace(-1, 1, nMesh);
            cutoffRelDepth = relDepths(1) * (1 - ti(cutoffIdx)) / 2;

        otherwise
            error('Vis:calcSubstrateCutoffDepth:badMode', ...
                'substrateCutoff must be ''none'', ''decay'', or ''fixedMultiple'' (got ''%s'').', substrateCutoff);
    end
end
