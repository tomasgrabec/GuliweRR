classdef StructureControlGUI < handle
%%% StructureControlGUI -- combined structure-builder + control/launcher window
%%%
%%% One uifigure window: a tab per layer (Substrate + film layers, added/removed
%%% dynamically via uitabgroup) for Material/CalcParams input, plus a control panel
%%% that runs the forward-calc classes on the resulting Structure.
%%%
%%% See docs/reference/architecture/gui.md for the full design rationale.
%%%
%%% KNOWN LIMITATIONS:
%%%   - Running a film layer standalone through the HS branch uses a hardcoded
%%%     relDepth=20 truncation depth rather than an editable field.
%%%
%%% Constructing with an existing Structure (initialStructure) reverse-populates
%%% every tab (dirMat, density, symmetry, Cij, relDepth/realThickness, legOrder)
%%% from structIn.mater(k)/calcParams(k) -- see loadFromStructure/
%%% populatePanelFromLayer.
%%%
%%% "Open viewer after Run" checkbox: on a successful Run, opens
%%% Gui.combinedViewer for any Single or AngSweep result,
%%% Gui.dispersionViewer for an isoWavevec 2-axis sweep,
%%% Gui.isofreqCombinedViewer for an isoFreq Single or AngSweep result, or
%%% Gui.freqAngSweep2LayeredViewer for an isoFreq 2-axis sweep. The
%%% Combined/IsofreqCombined viewers dispatch Single vs. AngSweep internally in
%%% their own bottom halves. Text and enabled state track the (layerMode, family,
%%% calcType) toggles live, via refreshOpenViewerOption.
%%%
%%% GuliweRR v6.0

    properties
        fig
        tabGroup
        panels          cell = {}      % {1,k}: panels{k} is a struct with .tab
                                        % .isSubstrate .dirMatFields(2x3 cell)
                                        % .dirMat(3x3) .libraryDropdown .hillToggle
                                        % .hillSource .rhoField .symmetryDropdown
                                        % .cijFields(6x6 cell, lower triangle empty)
                                        % .depthField .legOrderField. A cell, NOT a
                                        % "(1,:) struct"-typed property: a size/class-
                                        % constrained struct property starts as a
                                        % 0-field struct with no default, which then
                                        % rejects the first panels(1)=structWithFields
                                        % assignment as "dissimilar structures"
        statusLine
        addLayerBtn, removeLayerBtn, buildVarNameField, buildBtn

        layerModeGroup, layerModeHSBtn, layerModeFullBtn, layerModePlateBtn
        layerSelectDropdown
        familyGroup, familyIsoWavevecBtn, familyIsoFreqBtn
        calcTypeGroup, calcTypeSingleBtn, calcTypeAngSweepBtn, calcType2AxisBtn
        calcFieldLabels   % 1x3 uilabel
        calcFieldInputs   % 1x3 uieditfield (text)
        advancedToggleBtn
        advancedPanel
        advCalcMode, advRetentionLevel
        advMethod          % classification method -> the Calc objects' sawClassMethod
        advStandard        % its standard: 'method default' | 'high' | 'low' | 'custom' -> sawThreshold
        advCustomThreshold % numeric custom standard on the method's score scale
        runVarNameField, runBtn
        openViewerCheckbox
        lastResult                     % whichever Calc* object [Run] last produced
        lastViewer                     % viewer window [Run] last opened (openViewerCheckbox), or []

        maxGuiLayers    double    = 5   % soft usability cap on FILM layers (+1 for Substrate)
    end

    properties (Access = private)
        userMethod    = 'green'   % the classification method the user last picked; restored
                                  % when calcMode leaves 'ritzOnly' or the layer mode
                                  % offers it again (refreshClassificationOptions)
        lastSetMethod = 'green'   % the method refreshClassificationOptions last showed --
                                  % a different dropdown value means the user picked one
    end

    properties (SetObservable)
        structureModel  Structure = Structure()   % SetObservable drives onStructureChanged
                                                    % via the PostSet listener registered in
                                                    % the constructor -- own block because
                                                    % attributes like SetObservable apply to
                                                    % the whole properties block, not inline
    end

    methods
        function obj = StructureControlGUI(initialStructure)
        %%% Construct the combined structure-builder + control window
        %%% >> obj = StructureControlGUI(initialStructure)
        %%%
        %%% INPUT: (opt) initialStructure .. Structure to reverse-populate the tabs
        %%%                                  from; omit/empty for a blank start
            obj.buildFigure();
            addlistener(obj, 'structureModel', 'PostSet', @(~,~) obj.onStructureChanged());
            obj.onStructureChanged();

            if nargin > 0 && ~isempty(initialStructure) && initialStructure.nLay > 0
                obj.loadFromStructure(initialStructure);
            end
        end

        function show(obj)
        %%% Make the window visible (constructed invisible so headless/test use is safe)
        %%% >> show(obj)
            obj.fig.Visible = 'on';
        end
    end

    %%% ===================================================================
    %%% Figure construction
    %%% ===================================================================
    methods (Access = private)
        function buildFigure(obj)
            % buildControlPanel's gl needs >=596px for its 8 fixed rows alone (508px
            % summed + 8x10px RowSpacing + 8px padding), before row 5's own 'fit'
            % height (content-driven, not resolvable headlessly) or mainLayout's own
            % padding. At 680 the Run button falls below the visible
            % window; 780 leaves real margin rather than chasing the exact total.
            obj.fig = uifigure('Name', 'GuliweRR -- Structure & Control', ...
                'Position', [100 100 1000 780], 'Visible', 'off');

            mainLayout = uigridlayout(obj.fig, [1 2]);
            mainLayout.ColumnWidth = {'2x', '1x'};
            mainLayout.Padding = [4 4 4 4];

            leftLayout = uigridlayout(mainLayout, [3 1]);
            leftLayout.Layout.Row = 1; leftLayout.Layout.Column = 1;
            leftLayout.RowHeight = {'1x', 32, 24};
            leftLayout.Padding = [4 4 4 4];

            obj.tabGroup = uitabgroup(leftLayout);
            obj.tabGroup.Layout.Row = 1;

            btnRow = uigridlayout(leftLayout, [1 4]);
            btnRow.Layout.Row = 2;
            btnRow.Padding = [0 0 0 0];
            btnRow.ColumnWidth = {'1x', '1x', 110, '1x'};
            obj.addLayerBtn = uibutton(btnRow, 'Text', '+ Add layer', ...
                'ButtonPushedFcn', @(~,~) obj.onAddLayer());
            obj.removeLayerBtn = uibutton(btnRow, 'Text', '- Remove active layer', ...
                'Tooltip', 'Removes whichever tab is currently open (Substrate is never removable)', ...
                'ButtonPushedFcn', @(~,~) obj.onRemoveLayer());
            obj.buildVarNameField = uieditfield(btnRow, 'text', 'Value', 'guliStruct', ...
                'Tooltip', 'Base-workspace variable name that Build Structure assigns to');
            obj.buildBtn = uibutton(btnRow, 'Text', 'Build Structure', ...
                'ButtonPushedFcn', @(~,~) obj.onBuild());

            obj.statusLine = GuiCommon.statusLabel(leftLayout, 3);

            obj.buildControlPanel(mainLayout);

            obj.addPanel(1, true);   % Substrate tab, always present
        end

        function buildControlPanel(obj, parent)
            gl = uigridlayout(parent, [9 1]);
            gl.Layout.Row = 1; gl.Layout.Column = 2;
            % Row heights for the two- and three-radiobutton uibuttongroups (rows 1, 3,
            % 4) must clear both the group's title bar AND the full button stack below
            % it, or the title text overlaps the top button.
            gl.RowHeight = {104, 30, 78, 104, 'fit', 26, 200, 22, 34};
            gl.Padding = [4 4 4 4];
            %%% Scrollable safety net -- row 5 ('fit') and every fixed-pixel row
            %%% here are a best guess at natural content height, which can still be
            %%% wrong on an unfamiliar platform/font-metrics combination. If gl's
            %%% real total height ever exceeds what mainLayout's column 2 gives it,
            %%% this turns silent clipping/overlap into a visible scrollbar.
            gl.Scrollable = 'on';

            %%% layer-mode toggle -- Single layer (HS) / Full stack (Layered) /
            %%% Free plate (Lamb); same button-stack spacing as the 3-button
            %%% calcTypeGroup on row 4
            obj.layerModeGroup = uibuttongroup(gl, 'Title', 'Calculation target', ...
                'SelectionChangedFcn', @(~,~) obj.onLayerModeToggle());
            obj.layerModeGroup.Layout.Row = 1;
            obj.layerModeHSBtn = uiradiobutton(obj.layerModeGroup, ...
                'Text', 'Single layer (HS)', 'Position', [10 50 200 22]);
            obj.layerModeFullBtn = uiradiobutton(obj.layerModeGroup, ...
                'Text', 'Full stack (Layered)', 'Position', [10 27 200 22]);
            obj.layerModePlateBtn = uiradiobutton(obj.layerModeGroup, ...
                'Text', 'Free plate (Lamb)', 'Position', [10 4 200 22]);
            obj.layerModeHSBtn.Value = true;

            %%% HS-branch layer selector
            obj.layerSelectDropdown = uidropdown(gl, 'Items', {'Substrate'}, 'ItemsData', 1, ...
                'ValueChangedFcn', @(~,~) obj.onLayerSelected());
            obj.layerSelectDropdown.Layout.Row = 2;

            %%% Layered-branch family toggle -- isoWavevec / isoFreq
            obj.familyGroup = uibuttongroup(gl, 'Title', 'Family', 'Visible', 'off', ...
                'SelectionChangedFcn', @(~,~) obj.onFamilyToggle());
            obj.familyGroup.Layout.Row = 3;
            obj.familyIsoWavevecBtn = uiradiobutton(obj.familyGroup, ...
                'Text', 'isoWavevec', 'Position', [10 24 140 22]);
            obj.familyIsoFreqBtn = uiradiobutton(obj.familyGroup, ...
                'Text', 'isoFreq', 'Position', [10 2 140 22]);
            obj.familyIsoWavevecBtn.Value = true;

            %%% calc-type toggle -- Single / AngSweep / 2-axis sweep
            obj.calcTypeGroup = uibuttongroup(gl, 'Title', 'Calc type', ...
                'SelectionChangedFcn', @(~,~) obj.onCalcTypeToggle());
            obj.calcTypeGroup.Layout.Row = 4;
            obj.calcTypeSingleBtn = uiradiobutton(obj.calcTypeGroup, ...
                'Text', 'Single', 'Position', [10 50 150 22]);
            obj.calcTypeAngSweepBtn = uiradiobutton(obj.calcTypeGroup, ...
                'Text', 'AngSweep', 'Position', [10 27 150 22]);
            obj.calcType2AxisBtn = uiradiobutton(obj.calcTypeGroup, ...
                'Text', '2-axis sweep', 'Position', [10 4 150 22]);
            obj.calcTypeSingleBtn.Value = true;

            %%% generic calc fields, relabelled per active combination. The label is
            %%% stacked ABOVE its field rather than beside it in a fixed 110px column:
            %%% long labels like "angle range, e.g. 0:2:180" get clipped there, and
            %%% stacking gives them the full row width.
            fieldsGl = uigridlayout(gl, [6 1]);
            fieldsGl.Layout.Row = 5;
            %%% 'fit' (auto-size to rendered content) rather than fixed pixel
            %%% heights: a label wrapping to 2 lines on some platform's font
            %%% metrics grows its row instead of overflowing into the one below
            %%% (gl.Scrollable above is the matching outer-grid safety net).
            fieldsGl.RowHeight = repmat({'fit'}, 1, 6);
            fieldsGl.Padding = [0 0 0 0];
            obj.calcFieldLabels = gobjects(1, 3);
            obj.calcFieldInputs = gobjects(1, 3);
            for i = 1:3
                lbl = uilabel(fieldsGl, 'Text', '');
                lbl.Layout.Row = 2*i - 1; lbl.Layout.Column = 1;
                ef = uieditfield(fieldsGl, 'text', 'Value', '');
                ef.Layout.Row = 2*i; ef.Layout.Column = 1;
                obj.calcFieldLabels(i) = lbl;
                obj.calcFieldInputs(i) = ef;
            end

            %%% advanced options -- collapsed by default
            % char(9656)/char(9662) built at runtime, not written as a literal escape:
            % '\x25b8' in a single-quoted char array is NOT an escape sequence in
            % MATLAB (only double-quoted strings support \x) and prints literally.
            obj.advancedToggleBtn = uibutton(gl, 'Text', ['Advanced ' char(9656)], ...
                'ButtonPushedFcn', @(~,~) obj.onAdvancedToggle());
            obj.advancedToggleBtn.Layout.Row = 6;
            obj.advancedPanel = uipanel(gl, 'Title', 'Advanced options', 'Visible', 'off');
            obj.advancedPanel.Layout.Row = 7;
            advGl = uigridlayout(obj.advancedPanel, [5 2]);
            advGl.ColumnWidth = {110, '1x'};
            uilabel(advGl, 'Text', 'calcMode');
            obj.advCalcMode = uidropdown(advGl, 'Items', {'full', 'ritzOnly'}, 'Value', 'full', ...
                'ValueChangedFcn', @(~,~) obj.refreshClassificationOptions());
            uilabel(advGl, 'Text', 'retentionLevel');
            obj.advRetentionLevel = uidropdown(advGl, ...
                'Items', {'full', 'saw_only', 'none'}, 'Value', 'full');
            %%% classification method and standard, next to calcMode (they depend on
            %%% it): passed to the Calc objects as sawClassMethod / sawThreshold.
            %%% Available methods follow calcMode and the layer mode -- see
            %%% refreshClassificationOptions. Switchable later without recomputation
            %%% in the viewers' classification bar.
            uilabel(advGl, 'Text', 'Classification');
            [methodIds, methodNames] = CalcCommon.classificationMethods(false);
            obj.advMethod = uidropdown(advGl, 'Items', methodNames, ...
                'ItemsData', methodIds, 'Value', 'green', ...
                'Tooltip', ['Green: the default (needs calcMode ''full''); Refined: the ' ...
                            'spot-averaged Green score (half-space, about twice the cost); ' ...
                            'Ritz: the Ritz-only score (the only one on ''ritzOnly'')'], ...
                'ValueChangedFcn', @(src, ~) obj.onMethodPicked(src.Value));
            uilabel(advGl, 'Text', 'Standard');
            %%% not 'default': a graphics Value 'default' means "reset to the default"
            obj.advStandard = uidropdown(advGl, 'Items', {'method default', 'high', 'low', 'custom'}, ...
                'Value', 'method default', ...
                'Tooltip', StructureControlGUI.methodDefaultTooltip(methodIds, methodNames), ...
                'ValueChangedFcn', @(~,~) obj.refreshClassificationOptions());
            uilabel(advGl, 'Text', 'Custom value');
            obj.advCustomThreshold = uieditfield(advGl, 'numeric', 'Value', 0, 'Enable', 'off', ...
                'Tooltip', 'Custom standard on the method''s score scale (log10 amplitude ratio)');
            %%% No nVelGreen field: the Green axis step is the Calc classes' greenVelStep
            %%% (default GreenScore.DEFAULT_STEP).

            %%% open-viewer checkbox -- text/enabled state kept in sync with the
            %%% (layerMode, family, calcType) combination by refreshOpenViewerOption,
            %%% called from refreshCalcFields (already wired to all three toggles)
            obj.openViewerCheckbox = uicheckbox(gl, 'Text', 'Open viewer after Run', 'Value', true);
            obj.openViewerCheckbox.Layout.Row = 8;

            %%% run -- result variable name field alongside the button, same pattern as
            %%% the structure-side Build Structure / buildVarNameField
            runRowGl = uigridlayout(gl, [1 2]);
            runRowGl.Layout.Row = 9;
            runRowGl.ColumnWidth = {'1x', 90};
            runRowGl.Padding = [0 0 0 0];
            obj.runVarNameField = uieditfield(runRowGl, 'text', 'Value', 'guliResult', ...
                'Tooltip', 'Base-workspace variable name Run assigns the result to');
            obj.runBtn = uibutton(runRowGl, 'Text', 'Run', 'ButtonPushedFcn', @(~,~) obj.onRunClicked());

            %%% onLayerModeToggle() also calls refreshCalcFields() -- setting
            %%% obj.layerModeHSBtn.Value=true directly above (not via a selection-
            %%% change event) does not fire SelectionChangedFcn, so the initial
            %%% "2-axis sweep disabled in HS mode" state needs an explicit call here
            obj.onLayerModeToggle();
        end

        function addPanel(obj, k, isSubstrate)
        %%% Create tab k's widgets (Substrate or a film Layer), append to obj.panels,
        %%% and seed obj.structureModel with its default-valued layer.
        %%% >> addPanel(obj, k, isSubstrate)
            if isSubstrate
                tabTitle = 'Substrate';
            else
                tabTitle = sprintf('Layer %d', k - 1);
            end
            tab = uitab(obj.tabGroup, 'Title', tabTitle);
            gl = uigridlayout(tab, [9 1]);
            gl.RowHeight = {22, 90, 30, 30, 30, 22, '1x', 30, 30};
            gl.Padding = [4 4 4 4];

            %%% direction matrix (2x3) -- row 1 = surface normal (n1), row 2 = first/
            %%% propagation direction (n2); row order matches Material.m's dirMat
            %%% docblock and CompleteDirMat.m ([n1;n2] in, n1 preserved as row 1).
            dirMatTitleLbl = uilabel(gl, 'Text', 'Direction Matrix', 'FontWeight', 'bold');
            dirMatTitleLbl.Layout.Row = 1;
            dirGl = uigridlayout(gl, [2 5]);
            dirGl.Layout.Row = 2;
            dirGl.ColumnWidth = {110, 70, 70, 70, '1x'};
            dirRowLabels = {'Surface normal', 'First direction'};
            defaultDir = [1 0 0; 0 0 1];
            dirMatFields = cell(2, 3);
            for r = 1:2
                rowLbl = uilabel(dirGl, 'Text', dirRowLabels{r});
                rowLbl.Layout.Row = r; rowLbl.Layout.Column = 1;
                for c = 1:3
                    ef = uieditfield(dirGl, 'numeric', 'Value', defaultDir(r, c));
                    ef.Layout.Row = r; ef.Layout.Column = c + 1;
                    ef.ValueChangedFcn = @(~,~) obj.onDirMatEdited(k);
                    dirMatFields{r, c} = ef;
                end
            end

            %%% materials library -- picking an entry fills density, symmetry and
            %%% Cij (orientation is left as it is); editing any of those by hand
            %%% sets it back to "custom". ItemsData = library row, 0 = custom.
            %%% While a library entry is selected the symmetry dropdown is locked;
            %%% the "Hill avg" toggle beside it is the only symmetry switch then
            %%% (see updateMaterialLock).
            libTooltip = ['Fill density, symmetry and Cij from the materials library ' ...
                '(Material.ListLibrary). "custom" = values typed in by hand.'];
            libGl = obj.labelledFieldRow(gl, 3, 'library material', libTooltip);
            libGl.ColumnWidth = {130, 140, 70, '1x'};
            [libItems, libData] = obj.libraryItems();
            libDropdown = uidropdown(libGl, 'Items', libItems, 'ItemsData', libData, ...
                'Value', 0, 'Tooltip', libTooltip, ...
                'ValueChangedFcn', @(~,~) obj.onLibrarySelected(k));
            libDropdown.Layout.Column = 2;
            hillTooltip = ['Replace the anisotropic Cij by the isotropic Hill averages ' ...
                '(c11 = M, c44 = G; Hill 1952). Off restores the anisotropic constants.'];
            hillToggle = uibutton(libGl, 'state', 'Text', 'Hill avg', 'Value', false, ...
                'Tooltip', hillTooltip, 'ValueChangedFcn', @(~,~) obj.onHillToggled(k));
            hillToggle.Layout.Column = 3;

            %%% density -- the field column is 140 wide so the symmetry dropdown
            %%% below, which shares it, fits its longest item ("tetragonal II")
            rhoGl = obj.labelledFieldRow(gl, 4, 'density [g/cm3]');
            rhoField = uieditfield(rhoGl, 'numeric', 'Value', 1, ...
                'ValueChangedFcn', @(~,~) obj.onMaterialFieldEdited(k, @() obj.syncPanelToStructure(k)));

            %%% symmetry -- labelledFieldRow gives it the same {label, field,
            %%% filler} column pattern as the density/relDepth/legOrder rows,
            %%% so the dropdown box lines up with them
            symGl = obj.labelledFieldRow(gl, 5, 'symmetry');
            symmItems = {'isotropic', 'cubic', 'tetragonal I', 'tetragonal II', ...
                'trigonal I', 'trigonal II', 'hexagonal', 'orthorhombic', 'monoclinic', 'triclinic'};
            symDropdown = uidropdown(symGl, 'Items', symmItems, 'ItemsData', 1:10, 'Value', 1, ...
                'ValueChangedFcn', @(~,~) obj.onMaterialFieldEdited(k, @() obj.onSymmetryChanged(k)));

            %%% Cij grid (upper triangle of 6x6 only)
            cijTitleLbl = uilabel(gl, 'Text', 'Elasticity: Cij (Voigt matrix) [GPa]', 'FontWeight', 'bold');
            cijTitleLbl.Layout.Row = 6;
            cijGl = uigridlayout(gl, [6 6]);
            cijGl.Layout.Row = 7; cijGl.Padding = [0 0 0 0];
            cijFields = cell(6, 6);
            for r = 1:6
                for c = r:6
                    ef = uieditfield(cijGl, 'numeric', 'Value', 0);
                    ef.Layout.Row = r; ef.Layout.Column = c;
                    ef.ValueChangedFcn = @(~,~) obj.onMaterialFieldEdited(k, @() obj.recomputeAndRefreshCij(k));
                    cijFields{r, c} = ef;
                end
            end

            %%% relDepth (Substrate, direct) OR realThickness (film Layer, real units)
            if isSubstrate
                depthGl = obj.labelledFieldRow(gl, 8, 'relDepth');
                depthField = uieditfield(depthGl, 'numeric', 'Value', 20);
            else
                %%% realThickness (CalcParams.m docblock): real physical thickness, same
                %%% units as the wavelength ApplyWavelength is later called with -- not
                %%% fixed to microns by the code, but that is the intended working unit.
                depthTooltip = ['Real physical thickness (realThickness). Converted to ' ...
                    'relDepth via Structure.ApplyWavelength at calc time -- micrometres, ' ...
                    'matching the wavelength field''s own unit.'];
                depthGl = obj.labelledFieldRow(gl, 8, 'thickness [um]', depthTooltip);
                depthField = uieditfield(depthGl, 'numeric', 'Value', 1, 'Tooltip', depthTooltip);
            end
            depthField.ValueChangedFcn = @(~,~) obj.syncPanelToStructure(k);

            %%% legOrder
            legGl = obj.labelledFieldRow(gl, 9, 'legOrder');
            legField = uieditfield(legGl, 'numeric', 'Value', 50, ...
                'ValueChangedFcn', @(~,~) obj.syncPanelToStructure(k));

            obj.panels{k} = struct('tab', tab, 'isSubstrate', isSubstrate, ...
                'dirMatFields', {dirMatFields}, 'dirMat', eye(3), 'libraryDropdown', libDropdown, ...
                'hillToggle', hillToggle, 'hillSource', [], 'rhoField', rhoField, ...
                'symmetryDropdown', symDropdown, 'cijFields', {cijFields}, ...
                'depthField', depthField, 'legOrderField', legField);

            obj.onDirMatEdited(k);
            obj.onSymmetryChanged(k);   % also enables/disables Cij fields + syncs
        end

        function loadFromStructure(obj, structIn)
        %%% Reverse-populate every tab from an existing (non-empty) Structure --
        %%% one extra film-layer tab per structIn.nLay beyond the Substrate that
        %%% buildFigure already created, then every tab's fields.
        %%% >> loadFromStructure(obj, structIn)
        %%%
        %%% No maxGuiLayers cap here (unlike onAddLayer) -- a Structure built
        %%% programmatically is not subject to the GUI's own soft usability cap.
            for k = 2:structIn.nLay
                obj.addPanel(k, false);
            end
            for k = 1:structIn.nLay
                obj.populatePanelFromLayer(k, structIn.mater(k), structIn.calcParams(k));
            end
            GuiCommon.setStatusOk(obj.statusLine, ...
                sprintf('Loaded %d-layer structure "%s"', structIn.nLay, structIn.name));
        end

        function populatePanelFromLayer(obj, k, mater0, calcParams0)
        %%% Write one Material/CalcParams pair's values into panel k's widgets,
        %%% then recompute/sync exactly as a manual edit would (onSymmetryChanged
        %%% and onDirMatEdited each end with syncPanelToStructure).
        %%% >> populatePanelFromLayer(obj, k, mater0, calcParams0)
            p = obj.panels{k};
            for c = 1:3
                p.dirMatFields{1, c}.Value = mater0.dirMat(1, c);
                p.dirMatFields{2, c}.Value = mater0.dirMat(2, c);
            end
            p.libraryDropdown.Value = 0;
            p.hillToggle.Value = false;
            obj.panels{k}.hillSource = [];
            p.rhoField.Value = mater0.rho;
            p.symmetryDropdown.Value = mater0.elConsts.symmClassNo;
            for r = 1:6
                for cc = r:6
                    p.cijFields{r, cc}.Value = mater0.elConsts.voigt(r, cc);
                end
            end
            if p.isSubstrate
                p.depthField.Value = calcParams0.relDepth;
            else
                p.depthField.Value = calcParams0.realThickness;
            end
            p.legOrderField.Value = calcParams0.legOrder;

            obj.onSymmetryChanged(k);   % enables/disables Cij cells, recomputes, syncs
            obj.onDirMatEdited(k);      % completes dirMat, syncs
        end
    end

    %%% ===================================================================
    %%% Structure-building callbacks
    %%% ===================================================================
    methods
        function onLibrarySelected(obj, k)
        %%% Fill panel k's density, symmetry and Cij from the chosen library row,
        %%% Hill-averaged if the Hill avg toggle is on
        %%% >> onLibrarySelected(obj, k)
            p = obj.panels{k};
            row = p.libraryDropdown.Value;
            if row == 0
                % "custom": keep whatever is in the fields, as plain hand-typed values
                p.hillToggle.Value = false;
                obj.panels{k}.hillSource = [];
                obj.updateMaterialLock(k);
                obj.syncPanelToStructure(k);   % the layer loses its library name
                return
            end
            lib = Material.ListLibrary();
            entry = lib(row, :);
            mater0 = Material(ElastConsts(entry.consts{1}, entry.symmNo), entry.rho);
            p.rhoField.Value = mater0.rho;
            p.symmetryDropdown.Value = entry.symmNo;
            for r = 1:6
                for cc = r:6
                    p.cijFields{r, cc}.Value = mater0.elConsts.voigt(r, cc);
                end
            end
            if p.hillToggle.Value
                obj.applyHillAverage(k);
            end
            obj.onSymmetryChanged(k);   % enables/disables Cij cells, recomputes, syncs
            msg = sprintf('%s loaded from the materials library', p.libraryDropdown.Items{ ...
                p.libraryDropdown.ItemsData == row});
            if p.hillToggle.Value
                msg = [msg ' (Hill average)'];
            end
            if entry.status == "provisional"
                msg = sprintf('%s -- provisional values (%s)', msg, entry.source);
            end
            GuiCommon.setStatusOk(obj.statusLine, msg);
        end

        function onHillToggled(obj, k)
        %%% Hill avg on: replace panel k's Cij by its isotropic Hill averages
        %%% (anisotropic constants kept in hillSource); off: restore them
        %%% >> onHillToggled(obj, k)
            p = obj.panels{k};
            if p.hillToggle.Value
                obj.applyHillAverage(k);
            elseif ~isempty(p.hillSource)
                p.symmetryDropdown.Value = p.hillSource.symmNo;
                for r = 1:6
                    for c = r:6
                        p.cijFields{r, c}.Value = p.hillSource.voigt(r, c);
                    end
                end
                obj.panels{k}.hillSource = [];
            end
            obj.onSymmetryChanged(k);   % recomputes, updates the lock, syncs
        end

        function applyHillAverage(obj, k)
        %%% Store panel k's current (anisotropic) Cij in hillSource, then write the
        %%% isotropic Hill averages c11 = M, c44 = G into the grid
        %%% >> applyHillAverage(obj, k)
        %%%
        %%% The caller runs onSymmetryChanged afterwards (dependent cells, lock, sync).
            p = obj.panels{k};
            [vec, symmNo] = obj.readPanelCijVector(k);
            elC = ElastConsts(vec, symmNo);
            obj.panels{k}.hillSource = struct('symmNo', symmNo, 'voigt', elC.voigt);
            [mHill, gHill] = elC.GetHillAverages();
            p.symmetryDropdown.Value = 1;
            indepPos = ElastConsts.IndepVoigtPositions(1);   % isotropic: {c11, c44}
            hillVec = [mHill, gHill];
            for i = 1:numel(indepPos)
                p.cijFields{indepPos{i}(1), indepPos{i}(2)}.Value = hillVec(i);
            end
        end

        function onMaterialFieldEdited(obj, k, updateFcn)
        %%% A hand edit of density/symmetry/Cij: the panel no longer holds a library
        %%% entry or a Hill average, so set the library dropdown back to "custom" and
        %%% the Hill avg toggle off (the edited values stay), then run the update
        %%% >> onMaterialFieldEdited(obj, k, updateFcn)
            obj.panels{k}.libraryDropdown.Value = 0;
            obj.panels{k}.hillToggle.Value = false;
            obj.panels{k}.hillSource = [];
            updateFcn();
            obj.updateMaterialLock(k);
        end

        function updateMaterialLock(obj, k)
        %%% Symmetry dropdown locked while a library entry is selected or the Hill
        %%% average is on; Hill avg usable for any non-isotropic Cij (and to switch off)
        %%% >> updateMaterialLock(obj, k)
            p = obj.panels{k};
            locked = p.libraryDropdown.Value ~= 0 || p.hillToggle.Value;
            p.symmetryDropdown.Enable = matlab.lang.OnOffSwitchState(~locked);
            p.hillToggle.Enable = matlab.lang.OnOffSwitchState( ...
                p.hillToggle.Value || p.symmetryDropdown.Value ~= 1);
        end

        function onSymmetryChanged(obj, k)
        %%% Enable/disable Cij grid cells per symmetry class, then refresh values
        %%% >> onSymmetryChanged(obj, k)
            symmNo = obj.panels{k}.symmetryDropdown.Value;
            indepPos = ElastConsts.IndepVoigtPositions(symmNo);
            isIndep = false(6, 6);
            for i = 1:numel(indepPos)
                isIndep(indepPos{i}(1), indepPos{i}(2)) = true;
            end
            for r = 1:6
                for c = r:6
                    f = obj.panels{k}.cijFields{r, c};
                    if isIndep(r, c)
                        f.Enable = 'on';
                    else
                        f.Enable = 'off';
                    end
                end
            end
            obj.updateMaterialLock(k);
            obj.recomputeAndRefreshCij(k);
        end

        function recomputeAndRefreshCij(obj, k)
        %%% Rebuild ElastConsts from the enabled Cij cells, write dependent cells back
        %%% >> recomputeAndRefreshCij(obj, k)
            [vec, symmNo] = obj.readPanelCijVector(k);
            try
                elC = ElastConsts(vec, symmNo);
                for r = 1:6
                    for c = r:6
                        obj.panels{k}.cijFields{r, c}.Value = elC.voigt(r, c);
                    end
                end
                GuiCommon.setStatusOk(obj.statusLine, '');
            catch err
                GuiCommon.setStatusError(obj.statusLine, err.message);
            end
            obj.syncPanelToStructure(k);
        end

        function onDirMatEdited(obj, k)
        %%% Complete the direction matrix, surface a status message if adjusted
        %%% >> onDirMatEdited(obj, k)
        %%%
        %%% Material.CompleteDirMat has no wasAdjusted output, so the same check
        %%% is recomputed here. See docs/reference/architecture/gui.md.
            df = obj.panels{k}.dirMatFields;
            inputRows = [df{1,1}.Value, df{1,2}.Value, df{1,3}.Value; ...
                         df{2,1}.Value, df{2,2}.Value, df{2,3}.Value];
            try
                n1 = inputRows(1, :) ./ norm(inputRows(1, :));
                n2 = inputRows(2, :) ./ norm(inputRows(2, :));
                wasAdjusted = dot(n1, n2) > 1e-3;   % same tolerance CompleteDirMat uses internally
                obj.panels{k}.dirMat = Material.CompleteDirMat(inputRows);
                if wasAdjusted
                    GuiCommon.setStatusOk(obj.statusLine, 'First direction adjusted for orthogonality');
                else
                    GuiCommon.setStatusOk(obj.statusLine, '');
                end
            catch err
                GuiCommon.setStatusError(obj.statusLine, err.message);
            end
            obj.syncPanelToStructure(k);
        end

        function syncPanelToStructure(obj, k)
        %%% Push panel k's current field values into obj.structureModel
        %%% >> syncPanelToStructure(obj, k)
            p = obj.panels{k};
            [vec, symmNo] = obj.readPanelCijVector(k);
            try
                elC = ElastConsts(vec, symmNo);
            catch
                return   % invalid Cij state -- don't push a broken material
            end
            mater0 = Material(elC, p.rhoField.Value, p.dirMat, obj.panelMaterialName(k));
            if p.isSubstrate
                calcParams0 = CalcParams(p.depthField.Value, p.legOrderField.Value);
            else
                calcParams0 = CalcParams([], p.legOrderField.Value, p.depthField.Value);
            end
            if k > obj.structureModel.nLay
                obj.structureModel = obj.structureModel.AddLayer(mater0, calcParams0);
            else
                obj.structureModel = obj.structureModel.ModifyLayer(k, mater0, calcParams0);
            end
        end

        function onAddLayer(obj)
        %%% Create a new film-layer tab on top of the stack
        %%% >> onAddLayer(obj)
            if numel(obj.panels) >= obj.maxGuiLayers + 1
                GuiCommon.setStatusError(obj.statusLine, ...
                    'At GUI layer cap -- build additional layers programmatically (Structure.AddLayer)');
                return
            end
            k = numel(obj.panels) + 1;
            obj.addPanel(k, false);
            obj.tabGroup.SelectedTab = obj.panels{k}.tab;
        end

        function onRemoveLayer(obj)
        %%% Delete whichever tab is currently open; Substrate (index 1) is never
        %%% removable
        %%% >> onRemoveLayer(obj)
            k = obj.activeTabIndex();
            if isempty(k) || k == 1
                return
            end
            delete(obj.panels{k}.tab);
            obj.panels(k) = [];
            obj.structureModel = obj.structureModel.RemoveLayer(k);
            obj.renumberLayerTabs();
            if k > numel(obj.panels)
                obj.tabGroup.SelectedTab = obj.panels{end}.tab;
            else
                obj.tabGroup.SelectedTab = obj.panels{k}.tab;   % tab that now occupies slot k
            end
        end

        function structureOut = onBuild(obj)
        %%% Assign the current Structure to a base-workspace variable named by the
        %%% buildVarNameField text box. structureModel is already kept live by every
        %%% field edit (syncPanelToStructure), so this button's only job is handing a
        %%% copy out to the command line.
        %%% >> structureOut = onBuild(obj)
            structureOut = obj.structureModel;
            varName = strtrim(obj.buildVarNameField.Value);
            if isempty(varName) || ~isvarname(varName)
                GuiCommon.setStatusError(obj.statusLine, ...
                    sprintf('Build Structure: "%s" is not a valid variable name.', varName));
                return
            end
            assignin('base', varName, structureOut);
            GuiCommon.setStatusOk(obj.statusLine, ...
                sprintf('Structure built: %d layer(s) -- assigned to base workspace as ''%s''', ...
                obj.structureModel.nLay, varName));
        end
    end

    %%% ===================================================================
    %%% Control-panel callbacks
    %%% ===================================================================
    methods
        function onStructureChanged(obj)
        %%% PostSet listener target -- refresh the HS layer dropdown and the
        %%% "Full stack" option's availability whenever structureModel changes.
        %%% >> onStructureChanged(obj)
            n = obj.structureModel.nLay;
            if n == 0
                items = {'Substrate'}; itemsData = 1;
            else
                items = cell(1, n);
                items{1} = 'Substrate';
                for i = 2:n
                    items{i} = sprintf('Layer %d', i - 1);
                end
                itemsData = 1:n;
            end
            prevVal = obj.layerSelectDropdown.Value;
            obj.layerSelectDropdown.Items = items;
            obj.layerSelectDropdown.ItemsData = itemsData;
            if isempty(prevVal) || prevVal > numel(itemsData)
                obj.layerSelectDropdown.Value = itemsData(1);
            else
                obj.layerSelectDropdown.Value = prevVal;
            end

            if n < 2
                obj.layerModeFullBtn.Enable = 'off';
                if isequal(obj.layerModeGroup.SelectedObject, obj.layerModeFullBtn)
                    obj.layerModeGroup.SelectedObject = obj.layerModeHSBtn;
                    obj.onLayerModeToggle();
                end
            else
                obj.layerModeFullBtn.Enable = 'on';
            end
        end

        function sel = selectionState(obj)
        %%% The active (layerMode, family, calcType) toggle combination, read
        %%% once as a struct of booleans
        %%% >> sel = selectionState(obj)
        %%%
        %%% OUTPUT: sel .. struct with fields isHS, isPlate, isIsoWavevec,
        %%%                isSingle, isAngSweep, is2Axis
        %%%
        %%% Public (not Access=private): same convention as every other
        %%% callback here, so a test can read the dispatch state directly.
        %%%
        %%% One reader for all four methods that branch on this combination
        %%% (onLayerModeToggle, refreshCalcFields, describeViewerForCurrentSelection,
        %%% onRunClicked), which MUST agree -- four hand-copied sets of isequal()
        %%% calls is exactly how that agreement lapses. isIsoWavevec is read
        %%% unconditionally even though familyGroup is hidden for HS/Plate: the
        %%% button group still holds a valid selection there, and only the
        %%% Layered branch ever consults it.
            sel.isHS         = isequal(obj.layerModeGroup.SelectedObject, obj.layerModeHSBtn);
            sel.isPlate      = isequal(obj.layerModeGroup.SelectedObject, obj.layerModePlateBtn);
            sel.isIsoWavevec = isequal(obj.familyGroup.SelectedObject, obj.familyIsoWavevecBtn);
            sel.isSingle     = isequal(obj.calcTypeGroup.SelectedObject, obj.calcTypeSingleBtn);
            sel.isAngSweep   = isequal(obj.calcTypeGroup.SelectedObject, obj.calcTypeAngSweepBtn);
            sel.is2Axis      = isequal(obj.calcTypeGroup.SelectedObject, obj.calcType2AxisBtn);
        end

        function onLayerModeToggle(obj)
        %%% Switch visible controls between the HS / Layered / Plate branches
        %%% >> onLayerModeToggle(obj)
            sel = obj.selectionState();
            isHS = sel.isHS;
            obj.layerSelectDropdown.Visible = isHS;
            %%% familyGroup (isoWavevec/isoFreq) is a Layered-only concept -- Plate
            %%% has no isoFreq family at all (no Green's function to scan --
            %%% docs/reference/architecture/plate.md), so it stays hidden for Plate as for HS.
            obj.familyGroup.Visible = ~isHS && ~sel.isPlate;

            if isHS
                %%% no 2-axis sweep class exists for HS (no CalcWavevecAngSweep2HS)
                obj.calcType2AxisBtn.Enable = 'off';
                if isequal(obj.calcTypeGroup.SelectedObject, obj.calcType2AxisBtn)
                    obj.calcTypeGroup.SelectedObject = obj.calcTypeSingleBtn;
                end
            else
                %%% Layered AND Plate both have a genuine 2-axis sweep class
                %%% (CalcWavevecAngSweep2Layered / CalcWavevecAngSweep2Plate).
                obj.calcType2AxisBtn.Enable = 'on';
            end
            obj.refreshCalcFields();
        end

        function onLayerSelected(obj)
        %%% HS-branch layer dropdown -- nothing to precompute, read at Run time
        %%% >> onLayerSelected(obj)
        end

        function onFamilyToggle(obj)
        %%% isoWavevec/isoFreq toggle -- relabel the calc fields
        %%% >> onFamilyToggle(obj)
            obj.refreshCalcFields();
        end

        function onCalcTypeToggle(obj)
        %%% Single/AngSweep/2-axis-sweep toggle -- relabel the calc fields
        %%% >> onCalcTypeToggle(obj)
            obj.refreshCalcFields();
        end

        function onAdvancedToggle(obj)
        %%% Collapse/expand the Advanced options panel, flip the disclosure arrow
        %%% >> onAdvancedToggle(obj)
            obj.advancedPanel.Visible = ~obj.advancedPanel.Visible;
            if obj.advancedPanel.Visible
                obj.advancedToggleBtn.Text = ['Advanced ' char(9662)];   % expanded
            else
                obj.advancedToggleBtn.Text = ['Advanced ' char(9656)];   % collapsed
            end
        end

        function refreshCalcFields(obj)
        %%% Relabel/show-hide the three generic calc-input fields for the active
        %%% (layerMode, family, calcType) combination -- see docs/reference/architecture/gui.md,
        %%% "The three generic calc fields".
        %%% >> refreshCalcFields(obj)
            sel = obj.selectionState();
            isSingle   = sel.isSingle;
            isAngSweep = sel.isAngSweep;

            labels = {'', '', ''};
            visible = [false false false];

            if sel.isHS
                if isSingle
                    labels{1} = 'angle [deg] (opt, 0)';
                else
                    labels{1} = 'angle range, e.g. 0:2:180';
                end
                visible(1) = true;
            elseif sel.isPlate
                %%% field1=thickness/wavevector axis, field2=angle axis, field3=extra
                %%% param -- same role convention as the isoWavevec
                %%% Layered branch below (field1=wavelength, field2=angle,
                %%% field3=vMax), substituting Plate's own parameters
                %%% (realWavVec/realThickness -- both optional for Single/AngSweep,
                %%% see CalcSinglePlate.m/CalcAngSweepPlate.m; wavVecs/realThickness
                %%% both REQUIRED for the 2-axis sweep, CalcWavevecAngSweep2Plate.m).
                if isSingle
                    labels = {'realWavVec (opt)', 'angle [deg] (opt, 0)', 'realThickness (opt)'};
                elseif isAngSweep
                    labels = {'realWavVec (opt)', 'angle range, e.g. 0:2:180', 'realThickness (opt)'};
                else
                    labels = {'wavVecs, e.g. 1:0.5:5', 'angle range, e.g. 0:2:180', 'realThickness'};
                end
                visible = [true true true];
            else
                if sel.isIsoWavevec
                    if isSingle
                        labels = {'wavelength [um]', 'angle [deg] (opt, 0)', 'vMax (opt)'};
                    elseif isAngSweep
                        labels = {'wavelength [um]', 'angle range', 'vMax (opt)'};
                    else
                        labels = {'wavelength range [um]', 'angle range', 'vMax (opt)'};
                    end
                else
                    if isSingle
                        labels = {'omega0 (normalised)', 'angle [deg] (opt, 0)', 'kParRange'};
                    elseif isAngSweep
                        labels = {'omega0 (normalised)', 'angle range', 'kParRange'};
                    else
                        labels = {'omega0 range', 'angle range', 'kParRange'};
                    end
                end
                visible = [true true true];
            end

            for i = 1:3
                obj.calcFieldLabels(i).Text = labels{i};
                obj.calcFieldLabels(i).Visible = visible(i);
                obj.calcFieldInputs(i).Visible = visible(i);
            end

            obj.refreshOpenViewerOption();
            obj.refreshClassificationOptions();
        end

        function refreshClassificationOptions(obj)
        %%% Sync the classification method / standard controls with calcMode and the
        %%% layer mode
        %%% >> refreshClassificationOptions(obj)
        %%%
        %%% 'ritzOnly': Ritz only (the dropdown shows it, disabled). Layered: Green |
        %%% Ritz (no spot score for stacks). Plate and isoFreq: no SAW classification,
        %%% disabled. The custom value is enabled for the 'custom' standard only.
        %%% The user's own pick (userMethod) is restored as soon as the controls offer
        %%% it again, so leaving 'ritzOnly' or a layered stack does not silently change
        %%% the method of the next run.
            sel = obj.selectionState();
            hasClassification = ~sel.isPlate && (sel.isHS || sel.isIsoWavevec);
            [ids, names] = CalcCommon.classificationMethods(~sel.isHS);
            if ~strcmp(obj.advMethod.Value, obj.lastSetMethod)
                obj.userMethod = obj.advMethod.Value;   % set since the last refresh
            end
            obj.advMethod.Items = names;
            obj.advMethod.ItemsData = ids;
            isRitzOnly = strcmp(obj.advCalcMode.Value, 'ritzOnly');
            if isRitzOnly
                obj.advMethod.Value = 'ritz';
            elseif ismember(obj.userMethod, ids)
                obj.advMethod.Value = obj.userMethod;
            else
                obj.advMethod.Value = 'green';
            end
            obj.lastSetMethod = obj.advMethod.Value;
            obj.advMethod.Enable = hasClassification && ~isRitzOnly;
            obj.advStandard.Enable = hasClassification;
            obj.advCustomThreshold.Enable = hasClassification && strcmp(obj.advStandard.Value, 'custom');
        end

        function onMethodPicked(obj, method)
        %%% Method dropdown callback: remember the user's pick (see
        %%% refreshClassificationOptions)
        %%% >> onMethodPicked(obj, method)
            obj.userMethod = method;
            obj.lastSetMethod = method;
        end

        function [method, threshold] = classificationSettings(obj)
        %%% The method and standard Run passes to the Calc objects
        %%% >> [method, threshold] = classificationSettings(obj)
        %%%
        %%% OUTPUT: method    .. 'green' | 'refined' | 'ritz' ('ritz' whenever calcMode
        %%%                      is 'ritzOnly', the only method there)
        %%%         threshold .. [] (the method's default) | 'high' | 'low' | numeric
            method = obj.advMethod.Value;
            if strcmp(obj.advCalcMode.Value, 'ritzOnly'); method = 'ritz'; end
            switch obj.advStandard.Value
                case 'method default'; threshold = [];
                case 'custom';         threshold = obj.advCustomThreshold.Value;
                otherwise;             threshold = obj.advStandard.Value;
            end
        end

        function refreshOpenViewerOption(obj)
        %%% Sync openViewerCheckbox's text/enabled state with the active
        %%% (layerMode, family, calcType) combination (see
        %%% describeViewerForCurrentSelection); a combination with no matching
        %%% viewer is disabled and force-unchecked so onRunClicked never tries to
        %%% open one.
        %%% >> refreshOpenViewerOption(obj)
            [available, label] = obj.describeViewerForCurrentSelection();
            if available
                obj.openViewerCheckbox.Enable = 'on';
                obj.openViewerCheckbox.Value = true;
                obj.openViewerCheckbox.Text = ['Open ' label ' after Run'];
            else
                obj.openViewerCheckbox.Enable = 'off';
                obj.openViewerCheckbox.Value = false;
                obj.openViewerCheckbox.Text = 'Open viewer after Run (n/a)';
            end
        end

        function [available, label] = describeViewerForCurrentSelection(obj)
        %%% Which viewer (if any) matches the CURRENTLY SELECTED (not yet run)
        %%% dispatch state -- mirrors onRunClicked's own isHS/isSingle/isAngSweep/
        %%% isIsoWavevec branching so the checkbox stays truthful about what Run is
        %%% about to produce.
        %%%
        %%% Single (HS or isoWavevec Layered) routes to the SAME combined viewer as
        %%% AngSweep: CombinedViewerGUI's bottom half dispatches on the result's own
        %%% class (CalcSingleHS/CalcSingleLayered -> SingleAngleViewerGUI's 1D curve,
        %%% CalcAngSweepHS/CalcAngSweepLayered -> AngSweepViewerGUI's 2D map), so
        %%% this checkbox does not need to distinguish them.
        %%% >> [available, label] = describeViewerForCurrentSelection(obj)
            sel = obj.selectionState();
            isSingle   = sel.isSingle;
            isAngSweep = sel.isAngSweep;
            is2Axis    = sel.is2Axis;

            available = false;
            label = '';
            if sel.isHS
                if isSingle || isAngSweep
                    available = true;
                    label = 'combined viewer';
                end
            elseif sel.isPlate
                %%% Single/AngSweep -> CombinedViewerGUIPlate (mirrors the HS/
                %%% isoWavevec Layered 'combined viewer' routing above); 2-axis
                %%% -> DispersionViewerGUIPlate. A genuine double sweep
                %%% (numel(angles)>1 AND numel(wavVecs)>1) still only routes to
                %%% the 2-axis dispersion viewer here, not DoubleSweepViewerGUIPlate
                %%% -- consistent with Layered, where openResultViewer never
                %%% auto-opens the double-sweep viewer either (reachable only via
                %%% Gui.doubleSweepViewer/Gui.doubleSweepViewerPlate directly).
                if isSingle || isAngSweep
                    available = true;
                    label = 'combined viewer (plate)';
                elseif is2Axis
                    available = true;
                    label = 'plate dispersion viewer';
                end
            else
                if sel.isIsoWavevec
                    if isSingle || isAngSweep
                        available = true;
                        label = 'combined viewer';
                    elseif is2Axis
                        available = true;
                        label = 'dispersion viewer';
                    end
                else
                    %%% isoFreq family: Single/AngSweep -> Gui.isofreqCombinedViewer
                    %%% (mirrors the isoWavevec routing above: one combined-viewer
                    %%% entry point that dispatches internally, see
                    %%% IsofreqCombinedViewerGUI.buildBottomViewer); 2-axis ->
                    %%% Gui.freqAngSweep2LayeredViewer. docs/reference/architecture/gui_viewers.md.
                    if isSingle || isAngSweep
                        available = true;
                        label = 'isoFreq combined viewer';
                    elseif is2Axis
                        available = true;
                        label = 'isoFreq frequency sweep viewer';
                    end
                end
            end
        end

        function onRunClicked(obj)
        %%% Dispatch to the matching Calc* class for the active toggle state, run it,
        %%% store the result on obj.lastResult -- see docs/reference/architecture/gui.md,
        %%% "Running a calculation and handing off to a viewer".
        %%% >> onRunClicked(obj)
            try
                calcMode = obj.advCalcMode.Value;
                retentionLevel = obj.advRetentionLevel.Value;
                [sawMethod, sawThr] = obj.classificationSettings();
                setClass = @(c) StructureControlGUI.applyClassification(c, sawMethod, sawThr);

                sel = obj.selectionState();
                isSingle   = sel.isSingle;
                isAngSweep = sel.isAngSweep;

                if sel.isHS
                    structureHS = obj.buildHSStructure();
                    if isSingle
                        angle = obj.parseFieldOrDefault(1, 0);
                        cs = setClass(Calc.singleHS(structureHS, angle, [], calcMode, retentionLevel));
                        obj.lastResult = cs.calcAll();
                    else
                        angles = obj.parseFieldOrDefault(1, (0:2:180)');
                        cs = setClass(Calc.angSweepHS(structureHS, angles, calcMode, retentionLevel));
                        obj.lastResult = cs.calc();
                    end
                elseif sel.isPlate
                    %%% No calcMode for Plate (no Green's function --
                    %%% docs/reference/architecture/plate.md) -- advCalcMode is
                    %%% simply not read here. retentionLevel is passed through
                    %%% UNCHECKED ('full'/'saw_only'/'none' are all selectable in the
                    %%% Advanced panel, but Plate only accepts 'full'/'none') --
                    %%% Calc.*Plate's own validation errors clearly if 'saw_only' is
                    %%% selected, caught by this method's try/catch like any other
                    %%% bad input. A Plate-aware retentionLevel dropdown is a known gap.
                    if isSingle
                        realWavVec = obj.parseFieldOrDefault(1, []);
                        angle = obj.parseFieldOrDefault(2, 0);
                        realThickness = obj.parseFieldOrDefault(3, []);
                        cs = Calc.singlePlate(obj.structureModel, angle, realThickness, realWavVec, [], retentionLevel);
                        obj.lastResult = cs.calcAll();
                    elseif isAngSweep
                        realWavVec = obj.parseFieldOrDefault(1, []);
                        angles = obj.parseFieldOrDefault(2, (0:2:180)');
                        realThickness = obj.parseFieldOrDefault(3, []);
                        cs = Calc.angSweepPlate(obj.structureModel, angles, realThickness, realWavVec, retentionLevel);
                        obj.lastResult = cs.calc();
                    else
                        wavVecs = obj.parseField(1);
                        angles = obj.parseFieldOrDefault(2, (0:2:180)');
                        realThickness = obj.parseField(3);
                        cs = Calc.wavevecAngSweep2Plate(obj.structureModel, angles, wavVecs, realThickness, retentionLevel);
                        obj.lastResult = cs.calc();
                    end
                else
                    if sel.isIsoWavevec
                        %%% vMax is optional: it only extends the Green axis beyond the
                        %%% layered band end ([] = the band end)
                        vMax = obj.parseFieldOrDefault(3, []);
                        if isSingle
                            wavelength = obj.parseField(1);
                            angle = obj.parseFieldOrDefault(2, 0);
                            appliedStruct = obj.structureModel.ApplyWavelength(wavelength);
                            cs = setClass(Calc.singleLayered(appliedStruct, angle, vMax, [], calcMode, retentionLevel));
                            obj.lastResult = cs.calcAll();
                        elseif isAngSweep
                            wavelength = obj.parseField(1);
                            angles = obj.parseFieldOrDefault(2, (0:2:180)');
                            appliedStruct = obj.structureModel.ApplyWavelength(wavelength);
                            cs = setClass(Calc.angSweepLayered(appliedStruct, angles, vMax, calcMode, retentionLevel));
                            obj.lastResult = cs.calc();
                        else
                            wavelengths = obj.parseField(1);
                            angles = obj.parseFieldOrDefault(2, (0:2:180)');
                            cs = setClass(Calc.wavevecAngSweep2Layered(obj.structureModel, angles, wavelengths, ...
                                vMax, calcMode, retentionLevel));
                            obj.lastResult = cs.calc();
                        end
                    else
                        kParRange = obj.parseField(3);
                        if isSingle
                            omega0 = obj.parseField(1);
                            angle = obj.parseFieldOrDefault(2, 0);
                            cs = Calc.isofreqSingleLayered(obj.structureModel, angle, omega0, kParRange);
                            obj.lastResult = cs.calc();
                        elseif isAngSweep
                            omega0 = obj.parseField(1);
                            angles = obj.parseFieldOrDefault(2, (0:2:180)');
                            cs = Calc.isofreqAngSweepLayered(obj.structureModel, angles, omega0, kParRange);
                            obj.lastResult = cs.calc();
                        else
                            omega0List = obj.parseField(1);
                            angles = obj.parseFieldOrDefault(2, (0:2:180)');
                            cs = Calc.freqAngSweep2Layered(obj.structureModel, angles, omega0List, kParRange);
                            obj.lastResult = cs.calc();
                        end
                    end
                end
                %%% assign the result out, as Build Structure does
                varName = strtrim(obj.runVarNameField.Value);
                if isempty(varName) || ~isvarname(varName)
                    GuiCommon.setStatusError(obj.statusLine, sprintf(...
                        'Done: %s (result variable name "%s" is invalid -- not assigned)', ...
                        class(obj.lastResult), varName));
                else
                    assignin('base', varName, obj.lastResult);
                    GuiCommon.setStatusOk(obj.statusLine, ...
                        sprintf('Done: %s -- assigned to base workspace as ''%s''', ...
                        class(obj.lastResult), varName));
                end
                if obj.openViewerCheckbox.Value
                    obj.openResultViewer();
                end
            catch err
                GuiCommon.setStatusError(obj.statusLine, err.message);
            end
        end

        function openResultViewer(obj)
        %%% Launch the viewer window matching obj.lastResult's class, store it on
        %%% obj.lastViewer and show() it. Gated on openViewerCheckbox by
        %%% onRunClicked, so this is reached only when
        %%% describeViewerForCurrentSelection found a match at Run time. A
        %%% previously opened viewer is left open, not closed.
        %%% >> openResultViewer(obj)
            if isa(obj.lastResult, 'CalcAngSweepHS') || isa(obj.lastResult, 'CalcAngSweepLayered') || ...
                    isa(obj.lastResult, 'CalcSingleHS') || isa(obj.lastResult, 'CalcSingleLayered')
                obj.lastViewer = Gui.combinedViewer(obj.lastResult);
            elseif isa(obj.lastResult, 'CalcAngSweepPlate') || isa(obj.lastResult, 'CalcSinglePlate')
                obj.lastViewer = Gui.combinedViewerPlate(obj.lastResult);
            elseif isa(obj.lastResult, 'CalcWavevecAngSweep2Layered')
                obj.lastViewer = Gui.dispersionViewer(obj.lastResult);
            elseif isa(obj.lastResult, 'CalcWavevecAngSweep2Plate')
                obj.lastViewer = Gui.dispersionViewerPlate(obj.lastResult);
            elseif isa(obj.lastResult, 'CalcIsofreqSingleLayered') || isa(obj.lastResult, 'CalcIsofreqAngSweepLayered')
                obj.lastViewer = Gui.isofreqCombinedViewer(obj.lastResult);
            elseif isa(obj.lastResult, 'CalcFreqAngSweep2Layered')
                obj.lastViewer = Gui.freqAngSweep2LayeredViewer(obj.lastResult);
            else
                return
            end
            obj.lastViewer.show();
        end
    end

    %%% ===================================================================
    %%% Run-time helpers
    %%% ===================================================================
    methods (Access = private)
        function [items, data] = libraryItems(~)
        %%% Dropdown items for the materials library: "custom" (0), then one item
        %%% per library row -- the default row by its name, others with variant/T
        %%% >> [items, data] = libraryItems(obj)
            lib = Material.ListLibrary();
            items = {'custom'};
            data = 0;
            for r = 1:height(lib)
                if lib.default(r) == 1
                    items{end+1} = char(lib.name(r)); %#ok<AGROW>
                else
                    items{end+1} = sprintf('%s (%s, T=%s)', lib.name(r), lib.variant(r), lib.T(r)); %#ok<AGROW>
                end
                data(end+1) = r; %#ok<AGROW>
            end
        end

        function name = panelMaterialName(obj, k)
        %%% Panel k's material name: the library item while one is selected (with
        %%% " (Hill avg)" when averaged), else '' (Material then keeps its default name)
        %%% >> name = panelMaterialName(obj, k)
            dd = obj.panels{k}.libraryDropdown;
            name = '';
            if dd.Value ~= 0
                name = dd.Items{dd.ItemsData == dd.Value};
                if obj.panels{k}.hillToggle.Value
                    name = [name ' (Hill avg)'];
                end
            end
        end

        function row = labelledFieldRow(~, parent, rowNo, labelText, tooltip)
        %%% One "<label> | <field> | filler" row of a layer tab, ready for the
        %%% caller to drop its own uieditfield/uidropdown into column 2
        %%% >> row = labelledFieldRow(obj, parent, rowNo, labelText, tooltip)
        %%%
        %%% The five rows that use this (library, density, symmetry, relDepth/
        %%% thickness, legOrder) must share one column layout for their fields to line up
        %%% -- the 140px field column exists so the symmetry dropdown's longest
        %%% item ("tetragonal II") isn't truncated, and the others follow it.
            row = uigridlayout(parent, [1 3]);
            row.Layout.Row = rowNo;
            row.ColumnWidth = {130, 140, '1x'};
            row.Padding = [0 0 0 0];
            lbl = uilabel(row, 'Text', labelText);
            if nargin > 4 && ~isempty(tooltip)
                lbl.Tooltip = tooltip;
            end
        end

        function [vec, symmNo] = readPanelCijVector(obj, k)
        %%% Panel k's independent Voigt entries, in the exact order
        %%% ElastConsts(vec, symmNo) expects
        %%% >> [vec, symmNo] = readPanelCijVector(obj, k)
        %%%
        %%% Which cells are independent is a function of the symmetry class
        %%% alone (ElastConsts.IndepVoigtPositions), so this read and the
        %%% symmetry it was taken under always travel together -- returned as a
        %%% pair for that reason. Deliberately does NOT construct the
        %%% ElastConsts itself: each of the three callers wraps that
        %%% construction in its own error handling (status line / silent
        %%% return / propagate to onRunClicked's try), and folding it in here
        %%% would move those boundaries.
            p = obj.panels{k};
            symmNo = p.symmetryDropdown.Value;
            indepPos = ElastConsts.IndepVoigtPositions(symmNo);
            vec = zeros(1, numel(indepPos));
            for i = 1:numel(indepPos)
                vec(i) = p.cijFields{indepPos{i}(1), indepPos{i}(2)}.Value;
            end
        end

        function k = activeTabIndex(obj)
        %%% Index into obj.panels of whichever tab is currently open in the tab group
        %%% >> k = activeTabIndex(obj)
            k = [];
            for i = 1:numel(obj.panels)
                if isequal(obj.panels{i}.tab, obj.tabGroup.SelectedTab)
                    k = i;
                    return
                end
            end
        end

        function renumberLayerTabs(obj)
        %%% Keep film-layer tab titles ("Layer N") contiguous with their position in
        %%% obj.panels after a removal shifts everything above it down by one
        %%% >> renumberLayerTabs(obj)
            for i = 2:numel(obj.panels)
                obj.panels{i}.tab.Title = sprintf('Layer %d', i - 1);
            end
        end

        function structureHS = buildHSStructure(obj)
        %%% Extract the selected layer's material alone as a fresh 1-layer Structure
        %%% >> structureHS = buildHSStructure(obj)
        %%%
        %%% Film-layer-as-HS uses a hardcoded relDepth=20 -- see class docblock's
        %%% KNOWN LIMITATIONS.
            k = obj.layerSelectDropdown.Value;
            p = obj.panels{k};
            [vec, symmNo] = obj.readPanelCijVector(k);
            elC = ElastConsts(vec, symmNo);
            mater0 = Material(elC, p.rhoField.Value, p.dirMat, obj.panelMaterialName(k));
            if p.isSubstrate
                relDepth_HS = p.depthField.Value;
            else
                relDepth_HS = 20;
            end
            structureHS = Structure(mater0, CalcParams(relDepth_HS, p.legOrderField.Value), 'HS', 500);
        end

        function val = parseField(obj, i)
        %%% Parse calc-field i as a required MATLAB numeric expression (scalar/range/list)
        %%% >> val = parseField(obj, i)
            str = obj.calcFieldInputs(i).Value;
            if isempty(strtrim(str))
                error('StructureControlGUI:onRunClicked:missingField', ...
                    '%s is required.', obj.calcFieldLabels(i).Text);
            end
            val = str2num(str); %#ok<ST2NM> -- deliberate: supports scalar/colon-range/list syntax
            if isempty(val)
                error('StructureControlGUI:onRunClicked:badField', ...
                    'Could not parse "%s" as a numeric value/range.', str);
            end
            val = val(:);
            if isscalar(val)
                val = val(1);
            end
        end

        function val = parseFieldOrDefault(obj, i, defaultVal)
        %%% Parse calc-field i, or return defaultVal if left blank
        %%% >> val = parseFieldOrDefault(obj, i, defaultVal)
            str = obj.calcFieldInputs(i).Value;
            if isempty(strtrim(str))
                val = defaultVal;
            else
                val = obj.parseField(i);
            end
        end
    end

    methods (Static, Access = private)
        function cs = applyClassification(cs, method, threshold)
        %%% Set a Calc object's classification method and standard before it runs
        %%% >> cs = StructureControlGUI.applyClassification(cs, method, threshold)
        %%%
        %%% INPUT: 1. cs        .. a Calc object with sawClassMethod / sawThreshold
        %%%        2. method    .. 'green' | 'refined' | 'ritz'
        %%%        3. threshold .. [] | 'high' | 'low' | numeric
            cs.sawClassMethod = method;
            cs.sawThreshold   = threshold;
        end

        function tip = methodDefaultTooltip(ids, names)
        %%% Tooltip of the Standard dropdown: each method's default standard, read off
        %%% CalcCommon.defaultThreshold
        %%% >> tip = StructureControlGUI.methodDefaultTooltip(ids, names)
            parts = cell(1, numel(ids));
            for k = 1:numel(ids)
                parts{k} = sprintf('%s ''%s''', names{k}, CalcCommon.defaultThreshold(ids{k}, false));
            end
            layeredIds = CalcCommon.classificationMethods(true);
            layeredStd = unique(cellfun(@(m) CalcCommon.defaultThreshold(m, true), layeredIds, ...
                'UniformOutput', false));
            tip = sprintf('method default: %s; layered stacks: ''%s''', strjoin(parts, ', '), ...
                strjoin(layeredStd, ''', '''));
        end
    end
end
