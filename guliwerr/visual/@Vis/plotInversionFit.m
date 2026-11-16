function [hOut, hAx] = plotInversionFit(measured, associationFinal, associationInitial, hAx)
%%% Plot every stage of the inversion pipeline against propagation angle, on
%%% one axes: which measured points were associated vs. dropped, their
%%% initial-guess calculated velocity (optional), and their final/recovered
%%% calculated velocity.
%%% >> [hOut, hAx] = Vis.plotInversionFit(measured, associationFinal, associationInitial, hAx)
%%%
%%% INPUT:  1. measured           .. MeasuredVelocities, the FULL dataset -- including
%%%                                 any point ModeAssociation.fromSweep dropped (the
%%%                                 "N of M measured point(s) had no matching mode
%%%                                 within tolerance" warning)
%%%         2. associationFinal   .. ModeAssociation to plot as "recovered", e.g.
%%%                                 inversion.association AFTER inversion.calc()
%%%  (opt)  3. associationInitial .. an EARLIER ModeAssociation snapshot of the same
%%%                                 points (e.g. the variable Inv.association(...)
%%%                                 returned BEFORE calc() ran) to overlay as the
%%%                                 initial-guess fit. [] or omitted to skip that series.
%%%                                 Works because ModeAssociation is a value class: the
%%%                                 variable you built the ElasticityInversion FROM is
%%%                                 never touched by calc() -- only the instance's own
%%%                                 copy (inversion.association) is
%%%  (opt)  4. hAx                .. axes handle (new figure if omitted)
%%% OUTPUT: hOut .. struct of the scatter handles actually drawn (fields depend on
%%%                 which binding classes and series are present)
%%%         hAx  .. axes handle
%%%
%%% A measured point counts as "associated" if associationFinal contains a point at
%%% the same angle with the same measuredVelocity (both are copied verbatim from
%%% `measured` by ModeAssociation.fromSweep and never altered by rebind, so this match
%%% is exact, not a tolerance guess). Points failing that match were dropped -- either
%%% no candidate of the right wave type existed, or none was within velocityTolerance.

    if nargin < 3; associationInitial = []; end
    if nargin < 4 || isempty(hAx); hAx = axes(figure()); end %#ok<LAXES>

    pts = associationFinal.points;
    assocAngle = [pts.angleDeg];
    assocVelMeas = [pts.measuredVelocity];

    isAssociated = false(measured.nPoints, 1);
    for i = 1:measured.nPoints
        onAngle = abs(assocAngle - measured.angleDeg(i)) < 1e-6;
        isAssociated(i) = any(onAngle) && any(abs(assocVelMeas(onAngle) - measured.velocity(i)) < 1e-9);
    end

    markerOf = containers.Map({'ritzMode', 'limiting', 'pseudoSAW'}, {'o', 's', '^'});

    hold(hAx, 'on');
    hOut = struct();

    %%% Dropped measured points -- no calculated counterpart exists for these at any stage
    if any(~isAssociated)
        hOut.measuredDropped = scatter(hAx, measured.angleDeg(~isAssociated), measured.velocity(~isAssociated), ...
            50, [0.6 0.6 0.6], 'x', 'LineWidth', 1.5, 'DisplayName', ...
            sprintf('measured, dropped (%d)', sum(~isAssociated)));
    end

    %%% Measured / initial-guess / recovered, per binding class, for every ASSOCIATED point
    bindingClass = {pts.bindingClass};
    classes = unique(bindingClass);
    for k = 1:numel(classes)
        mask   = strcmp(bindingClass, classes{k});
        marker = markerOf(classes{k});
        angleK = assocAngle(mask);
        velMeasK = assocVelMeas(mask);

        hOut.(['measured_' classes{k}]) = scatter(hAx, angleK, velMeasK, ...
            40, 'k', marker, 'DisplayName', sprintf('measured (%s)', classes{k}));

        if ~isempty(associationInitial)
            velGuess = localMatchCalculatedVelocity(associationInitial.points, angleK, velMeasK);
            hOut.(['guess_' classes{k}]) = scatter(hAx, angleK, velGuess, ...
                35, [0.8500 0.3250 0.0980], marker, 'DisplayName', sprintf('initial guess (%s)', classes{k}));
        end

        velRecovered = [pts(mask).calculatedVelocity];
        hOut.(['recovered_' classes{k}]) = scatter(hAx, angleK, velRecovered, ...
            30, 'r', marker, 'filled', 'DisplayName', sprintf('recovered (%s)', classes{k}));
    end

    xlabel(hAx, 'propagation angle [deg]');
    ylabel(hAx, 'velocity [km/s]');
    title(hAx, sprintf('Inversion fit -- %d of %d measured point(s) associated', ...
        sum(isAssociated), measured.nPoints));
    legend(hAx, 'Location', 'best');
    grid(hAx, 'on');
end

function velOut = localMatchCalculatedVelocity(pts, angleList, velMeasList)
%%% For each (angle, measuredVelocity) pair, find the matching entry in a DIFFERENT
%%% ModeAssociation snapshot of the same points and return ITS calculatedVelocity.
    velOut = nan(size(angleList));
    ptsAngle  = [pts.angleDeg];
    ptsVelMeas = [pts.measuredVelocity];
    ptsVelCalc = [pts.calculatedVelocity];
    for k = 1:numel(angleList)
        onAngle = abs(ptsAngle - angleList(k)) < 1e-6;
        match = onAngle & abs(ptsVelMeas - velMeasList(k)) < 1e-9;
        if any(match)
            iMatch = find(match, 1);
            velOut(k) = ptsVelCalc(iMatch);
        end
    end
end
