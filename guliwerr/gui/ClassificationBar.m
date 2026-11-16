classdef ClassificationBar < handle
%%% ClassificationBar -- the viewers' classification row: method, standard, Reclassify
%%%
%%% One row of controls, built into a row of a viewer's own layout (the second
%%% top-bar row of SingleAngleViewerGUI, AngSweepViewerGUI and DispersionViewerGUI):
%%%   Classification method .. Green | Refined | Ritz (layered: Green | Ritz)
%%%   Standard   .. high | low | custom; custom enables a numeric field on the
%%%                 method's own score scale (the scale 'high' is compared on)
%%%   dB label   .. the selected standard's value in dB (20 x value)
%%%   Reclassify .. classifyTarget.classify(method, threshold): re-thresholds the
%%%                 stored candidates, no recomputation (Refined solves its spot
%%%                 directions on first use and prints the cost)
%%%   status     .. the active method and standard, or the stale banner
%%% After a successful Reclassify the bar refreshes itself and calls OnReclassified,
%%% which the owning viewer uses to redraw and to pass the refresh on (nested
%%% single-angle viewer, linked modal-shape viewer).
%%%
%%% UNAVAILABLE METHODS. uidropdown has no per-item Enable, so:
%%%   - a method that does not exist for the geometry is left out of the list
%%%     (Refined on a layered stack: there is no spot score for stacks);
%%%   - a method that exists but needs a recompute stays in the list, greyed
%%%     (italic grey item style), and selecting it reverts to the previous method
%%%     with the reason on the status line (Green and Refined on a 'ritzOnly'
%%%     result: they need calcMode 'full').
%%% Keeping such items visible says that the method exists and why it is off, which
%%% removing them would hide. The dropdown tooltip lists every method's
%%% availability: the reason, or Refined's first-use cost.
%%%
%%% STANDARDS. 'high' and 'low' are the method's own named standards
%%% (Vis.classificationScoreScale reads the values off the score classes). For
%%% 'green' and 'refined' on a half-space, 'low' is compared with a separate
%%% structure-relative score; the label says so. Changing the method resets the
%%% standard to that method's default (CalcCommon.defaultThreshold) and the custom
%%% value to its 'high' value, since the scales differ between methods.
%%%
%%% STALE RESULTS. A result classified by an older rule set (CalcCommon.isStale)
%%% shows a "stale -- recompute" banner and Reclassify is disabled; classify would
%%% refuse it anyway. A result from before the active method was stored (empty
%%% activeMethod / activeThreshold, only ever a stale one) shows its sawClassMethod
%%% (if unset, calcMode's own) and that method's default standard in the controls: a
%%% starting point for the recompute, not a statement about how it was classified (the
%%% banner says so).
%%%
%%% The bar reads its state (calcMode, activeMethod, activeThreshold, staleness)
%%% from `result` and calls classify on `classifyTarget`, by default the same
%%% object. A layered wavevector sweep (CalcWavevecAngSweep2Layered) passes one of
%%% its CalcAngSweepLayered slices as `result` and itself as `classifyTarget`, so
%%% every wavelength is reclassified and its 2D maps reassembled. Both are fixed at
%%% construction; the owner sets its reclassify handler with setReclassifyHandler.
%%% The methods offered come from CalcCommon.classificationMethods.
%%%
%%% GuliweRR v6.0

    properties
        grid                 % the bar's own uigridlayout
        methodDropdown
        standardDropdown
        customField          % numeric uieditfield, enabled for 'custom' only
        dbLabel
        reclassifyBtn
        statusLabel
    end

    properties (SetAccess = private)
        classifyTarget       % object whose classify(method, threshold) Reclassify calls
        OnReclassified = []  % @() called after a successful Reclassify (setReclassifyHandler)
        result               % CalcSingleHS/CalcSingleLayered or CalcAngSweepHS/CalcAngSweepLayered
        isLayered            % true for a layered stack
        unavailable = struct()   % field = method, value = the reason it is greyed
        lastMethod = ''      % method dropdown value before the last change (revert target)
    end

    methods
        function obj = ClassificationBar(parent, result, rowNo, classifyTarget)
        %%% Build the classification bar into one row of a parent layout
        %%% >> obj = ClassificationBar(parent, result, rowNo)
        %%% >> obj = ClassificationBar(parent, result, rowNo, classifyTarget)
        %%%
        %%% INPUT: 1. parent         .. uigridlayout to place the bar in
        %%%        2. result         .. CalcSingleHS/CalcSingleLayered or
        %%%                             CalcAngSweepHS/CalcAngSweepLayered, classified
        %%%                             (state source, see class docblock)
        %%%        3. rowNo          .. Layout.Row of the bar in parent
        %%%  (opt) 4. classifyTarget .. object whose classify() Reclassify calls
        %%%                             (default: result)
            validateattributes(result, {'CalcSingleHS', 'CalcSingleLayered', ...
                'CalcAngSweepHS', 'CalcAngSweepLayered'}, {'scalar'}, 'ClassificationBar', 'result');
            if nargin < 4 || isempty(classifyTarget); classifyTarget = result; end
            obj.result         = result;
            obj.classifyTarget = classifyTarget;
            obj.isLayered      = isa(result, 'CalcSingleLayered') || isa(result, 'CalcAngSweepLayered');
            obj.build(parent, rowNo);
            obj.refresh();
        end

        function setReclassifyHandler(obj, fcn)
        %%% Set the function called after a successful Reclassify -- the owning
        %%% viewer's redraw, which passes the refresh on
        %%% >> setReclassifyHandler(obj, fcn)
        %%%
        %%% INPUT: 1. fcn .. @() function handle, or [] for none
            obj.OnReclassified = fcn;
        end

        function refresh(obj)
        %%% Sync every control to the result's active classification and staleness
        %%% >> refresh(obj)
            [method, thr] = obj.activeState();
            obj.methodDropdown.Value = method;
            obj.lastMethod = method;
            if ischar(thr) || isstring(thr)
                obj.standardDropdown.Value = char(thr);
                obj.customField.Value = Vis.classificationScoreScale(method, obj.isLayered).high;
            else
                obj.standardDropdown.Value = 'custom';
                obj.customField.Value = thr;
            end
            obj.syncStandardControls();
            obj.showActiveStatus();
        end

        function [method, threshold] = selection(obj)
        %%% The method and standard the controls currently select
        %%% >> [method, threshold] = selection(obj)
        %%%
        %%% OUTPUT: method    .. 'green' | 'refined' | 'ritz'
        %%%         threshold .. 'high' | 'low' | numeric (custom)
            method = obj.methodDropdown.Value;
            threshold = obj.standardDropdown.Value;
            if strcmp(threshold, 'custom')
                threshold = obj.customField.Value;
            end
        end

        function onMethodChanged(obj)
        %%% Method dropdown callback: revert from a greyed method (reason on the
        %%% status line), else reset the standard to the method's default
        %%% >> onMethodChanged(obj)
            method = obj.methodDropdown.Value;
            if isfield(obj.unavailable, method)
                obj.methodDropdown.Value = obj.lastMethod;
                GuiCommon.setStatusError(obj.statusLabel, obj.unavailable.(method));
                return
            end
            obj.lastMethod = method;
            obj.standardDropdown.Value = CalcCommon.defaultThreshold(method, obj.isLayered);
            obj.customField.Value = Vis.classificationScoreScale(method, obj.isLayered).high;
            obj.syncStandardControls();
            obj.showActiveStatus();
        end

        function onStandardChanged(obj)
        %%% Standard dropdown / custom field callback: enable the field, relabel dB
        %%% >> onStandardChanged(obj)
            obj.syncStandardControls();
        end

        function onReclassify(obj)
        %%% Reclassify button callback: classifyTarget.classify(method, threshold),
        %%% then refresh and call OnReclassified. Errors go to the status line.
        %%% >> onReclassify(obj)
            [method, threshold] = obj.selection();
            try
                obj.classifyTarget.classify(method, threshold);
            catch err
                GuiCommon.setStatusError(obj.statusLabel, err.message);
                return
            end
            obj.refresh();
            if ~isempty(obj.OnReclassified)
                obj.OnReclassified();
            end
        end
    end

    methods (Access = private)
        function build(obj, parent, rowNo)
        %%% Build the controls (see class docblock)
        %%% >> build(obj, parent, rowNo)
            obj.grid = uigridlayout(parent, [1 8]);
            obj.grid.Layout.Row = rowNo;
            obj.grid.ColumnWidth = {135, 90, 65, 75, 70, 190, 90, '1x'};
            obj.grid.Padding = [0 0 0 0];

            uilabel(obj.grid, 'Text', 'Classification method:');
            [ids, names] = CalcCommon.classificationMethods(obj.isLayered);
            obj.methodDropdown = uidropdown(obj.grid, 'Items', names, 'ItemsData', ids, ...
                'ValueChangedFcn', @(~,~) obj.onMethodChanged());

            uilabel(obj.grid, 'Text', 'Standard:');
            obj.standardDropdown = uidropdown(obj.grid, 'Items', {'high', 'low', 'custom'}, ...
                'Tooltip', ['high: everything a good TGS map shows; low: the stronger waves ' ...
                            '(a subset of high); custom: a value on the method''s score scale'], ...
                'ValueChangedFcn', @(~,~) obj.onStandardChanged());
            obj.customField = uieditfield(obj.grid, 'numeric', 'Enable', 'off', ...
                'Tooltip', 'Custom standard on the method''s score scale (log10 amplitude ratio)', ...
                'ValueChangedFcn', @(~,~) obj.onStandardChanged());
            obj.dbLabel = uilabel(obj.grid, 'Text', '');

            obj.reclassifyBtn = uibutton(obj.grid, 'Text', 'Reclassify', ...
                'Tooltip', 'Apply the method and standard (no recomputation)', ...
                'ButtonPushedFcn', @(~,~) obj.onReclassify());
            obj.statusLabel = GuiCommon.statusLabel(obj.grid);

            obj.markUnavailable(ids, names);
        end

        function markUnavailable(obj, ids, names)
        %%% Grey the methods this result cannot use, and build the dropdown tooltip
        %%% >> markUnavailable(obj, ids, names)
            isRitzOnly = strcmp(obj.result.calcMode, 'ritzOnly');
            tip = cell(1, numel(ids));
            for k = 1:numel(ids)
                if isRitzOnly && ~strcmp(ids{k}, 'ritz')
                    obj.unavailable.(ids{k}) = sprintf(['%s needs calcMode ''full'' (a Green''s ' ...
                        'function); this result is ''ritzOnly'' -- recompute it with ''full''.'], names{k});
                    tip{k} = sprintf('%s: unavailable -- needs calcMode ''full''', names{k});
                    GuiCommon.greyDropdownItem(obj.methodDropdown, ids{k});
                elseif strcmp(ids{k}, 'refined')
                    tip{k} = ['Refined: ' obj.refinedCostText()];
                else
                    tip{k} = sprintf('%s: available, no recomputation', names{k});
                end
            end
            if obj.isLayered
                tip{end+1} = 'Refined: not available for layered stacks (no spot score)';
            end
            obj.methodDropdown.Tooltip = strjoin(tip, newline);
        end

        function txt = refinedCostText(obj)
        %%% Refined's first-use cost (or that it is already solved)
        %%% >> txt = refinedCostText(obj)
            if isprop(obj.result, 'cs')
                solved = all(cellfun(@(c) ~isempty(c.refinedSaw), obj.result.cs));
                cost = ['first use solves the spot directions of every angle (cost: ' ...
                        'docs/reference/architecture/forwardCalcs_modeClassification.md)'];
            else
                solved = ~isempty(obj.result.refinedSaw);
                cost = 'first use solves its spot directions (cost: CalcSingleHS.calcRefined)';
            end
            if solved
                txt = 'spot directions solved; switching costs nothing';
            else
                txt = cost;
            end
        end

        function [method, thr] = activeState(obj)
        %%% The result's active method and standard; for a result that stores none
        %%% (stale, see class docblock) its sawClassMethod and that method's default
        %%% >> [method, thr] = activeState(obj)
            method = obj.result.activeMethod;
            if isempty(method); method = obj.result.sawClassMethod; end
            if isempty(method); method = CalcCommon.modeSawClassMethod(obj.result.calcMode); end
            thr = obj.result.activeThreshold;
            if isempty(thr); thr = CalcCommon.defaultThreshold(method, obj.isLayered); end
        end

        function syncStandardControls(obj)
        %%% Enable the custom field for 'custom' only; show the selected standard in dB
        %%% >> syncStandardControls(obj)
            method = obj.methodDropdown.Value;
            stdName = obj.standardDropdown.Value;
            obj.customField.Enable = strcmp(stdName, 'custom');
            obj.dbLabel.Text = ClassificationBar.standardText(method, stdName, ...
                obj.customField.Value, obj.isLayered);
        end

        function showActiveStatus(obj)
        %%% Status: the active method and standard, or the stale banner
        %%% >> showActiveStatus(obj)
            if CalcCommon.isStale(obj.result)
                obj.reclassifyBtn.Enable = 'off';
                GuiCommon.setStatusError(obj.statusLabel, ...
                    'Stale classification (older rule set) -- recompute to reclassify');
                return
            end
            obj.reclassifyBtn.Enable = 'on';
            [method, thr] = obj.activeState();
            if ischar(thr) || isstring(thr)
                stdName = char(thr); val = NaN;
            else
                stdName = 'custom'; val = thr;
            end
            GuiCommon.setStatusOk(obj.statusLabel, sprintf('Active: %s, %s', ...
                Vis.classificationScoreScale(method, obj.isLayered).name, ...
                ClassificationBar.standardText(method, stdName, val, obj.isLayered)));
        end
    end

    methods (Static)
        function txt = standardText(method, standard, customValue, isLayered)
        %%% "<standard> = <value> (<dB> dB)" of a method's standard
        %%% >> txt = ClassificationBar.standardText(method, standard, customValue, isLayered)
        %%%
        %%% INPUT: 1. method      .. 'green' | 'refined' | 'ritz'
        %%%        2. standard    .. 'high' | 'low' | 'custom'
        %%%        3. customValue .. the custom value (used for 'custom' only)
        %%%        4. isLayered   .. true for a layered stack
        %%% OUTPUT: txt .. e.g. "high = 0.16 (3.2 dB)"; a 'low' standard on a
        %%%                structure-relative scale says so
            s = Vis.classificationScoreScale(method, isLayered);
            switch standard
                case 'high';   val = s.high;
                case 'low';    val = s.low;
                case 'custom'; val = customValue;
                otherwise
                    error('ClassificationBar:standardText:badStandard', ...
                        "standard must be 'high', 'low' or 'custom'; got '%s'.", standard);
            end
            txt = sprintf('%s = %.3g (%.1f dB)', standard, val, 20 * val);
            if strcmp(standard, 'low') && strcmp(s.lowScale, 'structure-relative')
                txt = [txt ', structure-relative'];
            end
        end
    end
end
