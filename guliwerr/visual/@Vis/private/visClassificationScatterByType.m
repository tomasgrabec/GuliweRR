function hPoints = visClassificationScatterByType(hAx, xVals, yVals, sizeVals, allTypes, isWeak)
%%% Scatter classified mode points, one series per type, in the fixed draw
%%% order, with a legend
%%% >> hPoints = visClassificationScatterByType(hAx, xVals, yVals, sizeVals, allTypes, isWeak)
%%%
%%% @Vis-private helper shared by Vis.plotClassificationDiagnostic (x=angle,
%%% y=velocity) and Vis.plotClassificationCurve (x=velocity, y=Green-curve
%%% projection). The two plot different coordinates, but their type colouring,
%%% draw order and legend behaviour must match: a reader comparing a swept map
%%% with a single-angle curve should not have to re-learn what blue means.
%%%
%%% INPUT: 1. hAx      .. target axes
%%%        2-3. xVals, yVals .. per-point coordinates (whatever the caller's own
%%%                             axes mean)
%%%        4. sizeVals .. per-point marker size
%%%        5. allTypes .. cell of per-point type names
%%%        6. isWeak   .. per-point logical: grade 'weak' (accepted only at the
%%%                       active, looser standard) -- drawn as an OPEN marker
%%% OUTPUT: hPoints .. array of scatter handles, one per type actually present.
%%%                    Each keeps its points' isWeak in UserData.isWeak, so an
%%%                    overlay restyle (GuiCommon.styleClassificationOverlay) can
%%%                    tell the weak points apart again
%%%
%%% Strong points are filled, weak points open (edge only), within ONE scatter per
%%% type: per-point face transparency (AlphaData 1 / 0, MarkerFaceAlpha 'flat') with
%%% an always-opaque edge in the type colour. One handle per type keeps the
%%% per-type checkboxes and the legend unchanged.
%%%
%%% Draw order is FIXED here, NOT unique()'s alphabetical order: limiting-velocity
%%% points are a Christoffel REFERENCE marker, not a classified mode, so they draw
%%% early and a genuine SAW/Rayleigh/pseudoRayleigh point at nearly the same spot
%%% draws on top. Alphabetical order puts 'limiting' after 'SAW' ('S' < 'l' in
%%% ASCII), hiding real SAW points behind limiting markers exactly where it
%%% matters physically — where a SAW branch approaches a bulk-wave limiting
%%% velocity. Any type not in the fixed order is appended in stable order.
%%%
%%% 'SAW' is the layered path's single type, sharing pseudoRayleigh's colour since
%%% the two never co-occur in one result.

    typeColors = struct( ...
        'Rayleigh',       [0      0.4470 0.7410], ...
        'pseudoRayleigh', [0.3010 0.7450 0.9330], ...
        'SAW',            [0.3010 0.7450 0.9330], ...
        'limiting',       [0.4660 0.6740 0.1880], ...
        'other',          [0.5    0.5    0.5   ]);
    typeDrawOrder = {'other', 'limiting', 'pseudoRayleigh', 'Rayleigh', 'SAW'};

    uTypesPresent = unique(allTypes);
    uTypes = [typeDrawOrder(ismember(typeDrawOrder, uTypesPresent)), ...
              setdiff(uTypesPresent, typeDrawOrder, 'stable')];

    hPoints = gobjects(1, numel(uTypes));
    for k = 1:numel(uTypes)
        mask = strcmp(allTypes, uTypes{k});
        weak = isWeak(mask);
        hPoints(k) = scatter(hAx, xVals(mask), yVals(mask), sizeVals(mask), ...
            typeColors.(uTypes{k}), 'filled', 'DisplayName', uTypes{k}, ...
            'MarkerEdgeColor', typeColors.(uTypes{k}), 'LineWidth', 1, ...
            'AlphaData', double(~weak(:)), 'AlphaDataMapping', 'none', ...
            'MarkerFaceAlpha', 'flat', 'MarkerEdgeAlpha', 1);
        hPoints(k).UserData = struct('isWeak', weak(:));
    end
    if ~isempty(hPoints); legend(hAx, 'show'); end
end
