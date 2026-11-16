classdef GuiCommon
%%% GuiCommon -- shared internals for the GuliweRR viewer-window family.
%%%
%%% Mechanical GUI steps only — status-line colouring, the crosshair marker,
%%% click-point extraction, the shared layout scaffolding, the min/max
%%% range-field pairs. Nothing here decides what a window shows or how a result
%%% is interpreted (each ...GUI class's own job) and nothing here does physics.
%%% Stateless static class: never instantiate it, never call it through Gui.*.
%%%
%%% Static helpers rather than a shared base class. The HS/Layered vs. Plate vs.
%%% isoFreq viewer families genuinely differ in which controls exist at all, in
%%% what a click means, and in what the y-axis is; static calls leave each
%%% class's decisions visible at its own call site, where inheritance would hide
%%% them in an override chain.
%%%
%%% Error-ID convention: the constructor-argument validators raise errors on the
%%% calling constructor's behalf. They take a `className` argument and throw
%%% <className>:constructor:<condition> (e.g.
%%% ModalShapeViewerGUI:constructor:invalidWindowMode), the three-part form of every
%%% other error, so the identifier stays the caller's own and greppable from it.
%%%
%%% Duck-typed helpers: captureAxisRange/restoreAxisRange read a viewer's
%%% xMinField/xMaxField/yMinField/yMaxField (plus cMinField/cMaxField where the
%%% viewer has them) by NAME rather than taking six handle arguments — the same
%%% small interface CombinedViewerGUI/DoubleSweepViewerGUI rely on.
%%%
%%% GuliweRR v6.0

    properties (Constant)
        STATUS_OK_COLOR    = [0    0.45 0   ]   % dark green -- confirmations/info
        STATUS_ERROR_COLOR = [0.6  0    0   ]   % dark red -- errors and blocked actions
        CROSSHAIR_COLOR    = [0.72 0.20 0.95]   % selection crosshair -- bright purple, stands
                                   % out on grey/white axes and busy maps; not orange,
                                   % which would clash with the ColMapGrayOrange colormap
        DIVIDER_COLOR      = [0.6  0.6  0.6 ]   % thin between-halves divider panel
        GREY_ITEM_COLOR    = [0.55 0.55 0.55]   % font of a greyed dropdown item (greyDropdownItem)
        OVERLAY_WEAK_ALPHA = 0.35               % ring opacity of a 'weak' entry in a
                                   % classification overlay (strong rings: 1)
        SIDE_PANEL_WIDTH   = 230                % px, plotWithSidePanel's control column:
                                   % fixed, so a wide window widens the plot, not the
                                   % range fields and checkboxes
        CONTROLS_WIDTH     = 330                % px, shapeControlsSplit's control column:
                                   % fixed, wide enough for the Angle / Mode (with its
                                   % velocity readout) / Source + Display rows
    end

    %%% ===================================================================
    %%% Status line
    %%% ===================================================================
    methods (Static)
        function lbl = statusLabel(parent, rowNo)
        %%% Build the standard bottom status label (empty, OK colour)
        %%% >> lbl = GuiCommon.statusLabel(parent, rowNo)
        %%%
        %%% INPUT: 1. parent .. uigridlayout to place it in
        %%%  (opt) 2. rowNo  .. Layout.Row to assign; omitted/[] leaves the
        %%%                     layout's own default placement
            lbl = uilabel(parent, 'Text', '', 'FontColor', GuiCommon.STATUS_OK_COLOR);
            if nargin > 1 && ~isempty(rowNo)
                lbl.Layout.Row = rowNo;
            end
        end

        function setStatusOk(statusLabel, msg)
        %%% Show a non-error status message (confirmation/info) in STATUS_OK_COLOR
        %%% >> GuiCommon.setStatusOk(statusLabel, msg)
            statusLabel.FontColor = GuiCommon.STATUS_OK_COLOR;
            statusLabel.Text = msg;
        end

        function setStatusError(statusLabel, msg)
        %%% Show an error/blocked-action status message in STATUS_ERROR_COLOR
        %%% >> GuiCommon.setStatusError(statusLabel, msg)
            statusLabel.FontColor = GuiCommon.STATUS_ERROR_COLOR;
            statusLabel.Text = msg;
        end
    end

    %%% ===================================================================
    %%% Figure / layout scaffolding
    %%% ===================================================================
    methods (Static)
        function [fig, mainLayout] = hostLayout(parentContainer, figName, figPosition, gridSize)
        %%% Host a viewer's content either in its own uifigure or in a
        %%% caller-supplied region
        %%% >> [fig, mainLayout] = GuiCommon.hostLayout(parentContainer, figName, figPosition, gridSize)
        %%%
        %%% INPUT: 1. parentContainer .. [] for standalone (own uifigure), or a
        %%%                              uigridlayout region to build into
        %%%        2. figName         .. uifigure Name, standalone case only
        %%%        3. figPosition     .. uifigure Position, standalone case only
        %%%        4. gridSize        .. [nRows nCols] for the returned mainLayout
        %%% OUTPUT: fig        .. the owned uifigure, or the ancestor figure of
        %%%                       parentContainer (so .fig/.show() work either way)
        %%%         mainLayout .. the top-level uigridlayout, Padding already [4 4 4 4]
        %%%
        %%% RowHeight/ColumnWidth stay the caller's job -- that is where the
        %%% families genuinely differ. Constructed 'Visible','off' so headless/
        %%% test use is safe, the convention every viewer already follows.
            if isempty(parentContainer)
                fig = uifigure('Name', figName, 'Position', figPosition, 'Visible', 'off');
                mainLayout = uigridlayout(fig, gridSize);
            else
                fig = ancestor(parentContainer, 'figure');
                mainLayout = uigridlayout(parentContainer, gridSize);
            end
            mainLayout.Padding = [4 4 4 4];
        end

        function showMaximised(fig)
        %%% Show a one-window composite viewer maximised (fills the screen; not the
        %%% separate macOS full-screen space): the combined viewers are the main
        %%% working windows, and their default Position is only the restored size
        %%% >> GuiCommon.showMaximised(fig)
        %%%
        %%% drawnow between the two: maximised in the same step as it is shown, the
        %%% window's first layout pass skips the modal-shape axes, which shapeAxesGrid
        %%% places by hand inside a plain uipanel, and they keep their pre-maximise
        %%% size in a corner of the larger panel
            fig.Visible = 'on';
            drawnow;
            fig.WindowState = 'maximized';
        end

        function [shapeParent, controlsParent, topLayout] = shapeControlsSplit(parent, padding)
        %%% The shared "shapes left, controls right (fixed CONTROLS_WIDTH)" two-column split
        %%% >> [shapeParent, controlsParent, topLayout] = GuiCommon.shapeControlsSplit(parent, padding)
        %%%
        %%% Each region is its own dedicated 1x1 cell: the modal-shape and control
        %%% classes fill whatever container they are given completely, so every
        %%% region must be pre-positioned within the shared layout rather than
        %%% handed the split itself. topLayout is returned so a caller inside a
        %%% multi-row parent can set its Layout.Row.
            topLayout = uigridlayout(parent, [1 2]);
            topLayout.ColumnWidth = {'1x', GuiCommon.CONTROLS_WIDTH};
            topLayout.Padding = padding;

            shapeParent = uigridlayout(topLayout, [1 1]);
            shapeParent.Layout.Column = 1; shapeParent.Padding = [0 0 0 0];
            controlsParent = uigridlayout(topLayout, [1 1]);
            controlsParent.Layout.Column = 2; controlsParent.Padding = [0 0 0 0];
        end

        function [shapeParent, controlsParent, bottomParent, splitLayout] = dividedTopBottom(parent, rowHeight, padding)
        %%% The shared combined-viewer body: shapes+controls on top, a thin
        %%% divider, a map/curve region below
        %%% >> [shapeParent, controlsParent, bottomParent, splitLayout] = GuiCommon.dividedTopBottom(parent, rowHeight, padding)
        %%%
        %%% INPUT: 1. parent    .. uifigure or uigridlayout region to build into
        %%%        2. rowHeight .. 3-element RowHeight, e.g. {'44x', 3, '56x'}
        %%%                        (the middle entry is the divider's own pixel
        %%%                        height) -- the ratio differs per family, which
        %%%                        is exactly why it is a parameter
        %%%        3. padding   .. Padding for the outer 3-row layout
        %%% OUTPUT: 4. splitLayout .. the outer 3-row layout itself, so a caller
        %%%                           placing it inside a multi-row parent can set
        %%%                           its Layout.Row
        %%%
        %%% The divider exists because the two halves otherwise sit flush against
        %%% each other with no visible boundary.
            splitLayout = uigridlayout(parent, [3 1]);
            splitLayout.RowHeight = rowHeight;
            splitLayout.Padding = padding;

            [shapeParent, controlsParent, topLayout] = GuiCommon.shapeControlsSplit(splitLayout, [0 0 0 0]);
            topLayout.Layout.Row = 1;

            dividerPanel = uipanel(splitLayout, 'BackgroundColor', GuiCommon.DIVIDER_COLOR, ...
                'BorderType', 'none');
            dividerPanel.Layout.Row = 2;

            bottomParent = uigridlayout(splitLayout, [1 1]);
            bottomParent.Layout.Row = 3; bottomParent.Padding = [0 0 0 0];
        end

        function [ax, sideGl, content] = plotWithSidePanel(parent, sideRowHeight, onClick)
        %%% The shared map/curve viewer body: one plot axes with a control side
        %%% panel of fixed width (SIDE_PANEL_WIDTH) beside it
        %%% >> [ax, sideGl, content] = GuiCommon.plotWithSidePanel(parent, sideRowHeight, onClick)
        %%%
        %%% INPUT: 1. parent        .. uifigure or uigridlayout region to build into
        %%%        2. sideRowHeight .. the side panel's RowHeight; its length sets
        %%%                            the row count too, so the two cannot
        %%%                            disagree. Every viewer's own row budget
        %%%                            differs (how many range-field pairs, whether
        %%%                            there is a type-checkbox panel), which is why
        %%%                            this is a parameter and the column split is not.
        %%%        3. onClick       .. ButtonDownFcn for the axes -- the AXES-LEVEL
        %%%                            hit (clicking empty background). Each viewer
        %%%                            additionally sets the same handler on the
        %%%                            plotted objects after every draw, since a hit
        %%%                            on a scatter/line never reaches the axes.
        %%% OUTPUT: ax      .. the uiaxes, already placed at row 1 / column 1
        %%%         sideGl  .. the side panel, ready to have controls appended in
        %%%                    row order (rangeFieldPair consumes two rows each)
        %%%         content .. the 1x2 split itself, so a caller inside a
        %%%                    multi-row parent can set its Layout.Row
            content = uigridlayout(parent, [1 2]);
            content.ColumnWidth = {'1x', GuiCommon.SIDE_PANEL_WIDTH};
            content.Padding = [0 0 0 0];

            ax = uiaxes(content);
            ax.Layout.Row = 1; ax.Layout.Column = 1;
            ax.ButtonDownFcn = onClick;

            sideGl = uigridlayout(content, [numel(sideRowHeight) 1]);
            sideGl.Layout.Row = 1; sideGl.Layout.Column = 2;
            sideGl.RowHeight = sideRowHeight;
            sideGl.Padding = [4 4 4 4];
        end

        function [shapeAxesCell, rgbAxesCell, keyAx] = shapeAxesGrid(shapeParent, relDepths)
        %%% Pre-build the {3 x nLay} uiaxes grid Vis.plotModalShape/
        %%% plotGreenModalShape render into via their hAxCellIn parameter, plus
        %%% the hidden axes of the polarisation-RGB display (dirs = 'rgb')
        %%% >> [shapeAxesCell, rgbAxesCell, keyAx] = GuiCommon.shapeAxesGrid(shapeParent, relDepths)
        %%%
        %%% INPUT: 1. shapeParent .. uigridlayout region to build the grid in
        %%%        2. relDepths   .. every layer's relDepth (layer 1 = substrate,
        %%%                          drawn at the BOTTOM, matching plotModalShape's
        %%%                          own row convention). Pass 1 for a
        %%%                          single-row grid (Green's response is always
        %%%                          exactly one "layer").
        %%% OUTPUT: shapeAxesCell .. {3 x numel(relDepths)} uiaxes
        %%%         rgbAxesCell   .. {1 x numel(relDepths)} uiaxes for the merged
        %%%                          RGB panel, same rows, spanning the left
        %%%                          ~half; built Visible='off'
        %%%         keyAx         .. uiaxes for Vis.plotPolarisationKey, right
        %%%                          half; built Visible='off'
        %%%
        %%% The two sets share one uipanel and are switched by
        %%% GuiCommon.setShapeDisplayMode; building both once keeps the fixed,
        %%% hand-placed geometry below in one place.
        %%%
        %%% Row heights weighted by Vis.calcModalShapeRowHeights — the SAME boosted
        %%% weights every caller's Vis.plotModalShape('boosted',...) call uses.
        %%% They are fixed once at construction rather than rebuilt every redraw,
        %%% so they must already match what the redraw will render with. Plain
        %%% relDepth-proportional weights make a real film-on-substrate case
        %%% (relDepth 20 vs 0.305) render the film as a ~1% sliver.
        %%%
        %%% Manually-positioned uiaxes inside a plain uipanel, NOT a uigridlayout:
        %%% a uigridlayout with 2+ children re-manages each child's
        %%% PositionConstraint on every drawnow, silently reverting any attempt to
        %%% remove a uiaxes' default label/tick margin and reintroducing a gap
        %%% between adjacent layers. RowSpacing does not help — it controls spacing
        %%% BETWEEN cells, while the gap is INSIDE each cell, between its boundary
        %%% and where content is drawn. A plain uipanel has no such layout manager,
        %%% so Position/InnerPosition set once here survive cla()/redraw/drawnow.
        %%%
        %%% topMargin/bottomMargin/layBottom/layHeightFrac below mirror
        %%% plotModalShape.m's classic-figure math, kept in sync by hand since a
        %%% caller-supplied axes grid cannot reuse that function's internal
        %%% Position computation. Both reserve room for the top row's title ONCE,
        %%% outside every axes' box, rather than per cell.
            nLay = numel(relDepths);
            rowWeights = Vis.calcModalShapeRowHeights(relDepths);

            %%% the panel re-pins its axes whenever it changes size (pinShapeAxes):
            %%% a resize -- e.g. the window maximised as it first opens -- can leave
            %%% them at their old pixel size in a corner of the larger panel
            pnl = uipanel(shapeParent, 'BorderType', 'none', 'AutoResizeChildren', 'off', ...
                'SizeChangedFcn', @(src, ~) GuiCommon.pinShapeAxes(src));
            pnl.Layout.Row = 1; pnl.Layout.Column = 1;

            topMargin     = 0.15;
            bottomMargin  = 0.025;
            plotHeight    = 1 - topMargin - bottomMargin;
            layHeightFrac = rowWeights * plotHeight;
            layBottom     = bottomMargin + [0, cumsum(layHeightFrac(1:end-1))];
            colWidth      = 1/3;

            shapeAxesCell = cell(3, nLay);
            rgbAxesCell   = cell(1, nLay);
            for ilay = 1:nLay
                for ic = 1:3
                    pos = [(ic-1)*colWidth + 0.02*colWidth, layBottom(ilay), ...
                           0.96*colWidth, layHeightFrac(ilay)];
                    shapeAxesCell{ic, ilay} = GuiCommon.fixedUiaxes(pnl, pos, 'on', 'shapeComponent');
                end
                %%% the RGB panel: as wide as one and a half component columns,
                %%% the key to its right (same split as visRgbShapeStack's classic figure)
                pos = [0.04, layBottom(ilay), 0.5, layHeightFrac(ilay)];
                rgbAxesCell{ilay} = GuiCommon.fixedUiaxes(pnl, pos, 'off', 'shapeRgb');
            end
            keyAx = GuiCommon.fixedUiaxes(pnl, [0.6, 0.35, 0.38, 0.38], 'off', 'shapeKey');
        end

        function ax = fixedUiaxes(pnl, pos, visible, tag)
        %%% One manually-positioned uiaxes of shapeAxesGrid (see its docblock
        %%% for why the position must be pinned this way). The Tag
        %%% ('shapeComponent' | 'shapeRgb' | 'shapeKey') tells the sets apart
        %%% for anything searching the panel, e.g. findobj(..., 'Tag', 'shapeComponent')
        %%% >> ax = GuiCommon.fixedUiaxes(pnl, pos, visible, tag)
            ax = uiaxes(pnl, 'Units', 'normalized', 'Position', pos, 'Visible', visible, 'Tag', tag);
            ax.PositionConstraint = 'innerposition';
            ax.InnerPosition = pos;
            setappdata(ax, 'pinnedPosition', pos);   % for pinShapeAxes
        end

        function pinShapeAxes(pnl)
        %%% shapeAxesGrid's panel SizeChangedFcn: put every axes built by
        %%% fixedUiaxes back to its own normalised inner position
        %%% >> GuiCommon.pinShapeAxes(pnl)
            axList = findobj(pnl.Children, 'flat', 'Type', 'axes');
            for k = 1:numel(axList)
                if isappdata(axList(k), 'pinnedPosition')
                    axList(k).Units = 'normalized';
                    axList(k).PositionConstraint = 'innerposition';
                    axList(k).InnerPosition = getappdata(axList(k), 'pinnedPosition');
                end
            end
        end

        function dd = shapeDisplayDropdown(parent, onChange)
        %%% The "Display:" choice of every modal-shape viewer: signed component
        %%% panels, or the merged polarisation-RGB panel
        %%% >> dd = GuiCommon.shapeDisplayDropdown(parent, onChange)
        %%%
        %%% INPUT: 1. parent   .. container (a uigridlayout cell)
        %%%        2. onChange .. ValueChangedFcn callback
        %%% OUTPUT: dd .. uidropdown, Value 'components' (default) | 'rgb' --
        %%%               the value is passed straight on as the plot functions'
        %%%               dirs argument via GuiCommon.shapeDirs
            dd = uidropdown(parent, ...
                'Items', {'Components (u_x, u_y, u_z)', 'Polarisation (RGB)'}, ...
                'ItemsData', {'components', 'rgb'}, 'Value', 'components', ...
                'Tooltip', ['Components: one signed panel per displacement direction. ' ...
                            'Polarisation (RGB): one panel, colour = displacement direction ' ...
                            '(red u_x, blue u_y, green u_z), brightness = |u|.'], ...
                'ValueChangedFcn', onChange);
        end

        function dirs = shapeDirs(displayDropdown)
        %%% The plot functions' dirs argument for the current display choice
        %%% >> dirs = GuiCommon.shapeDirs(displayDropdown)
            if strcmp(displayDropdown.Value, 'rgb')
                dirs = 'rgb';
            else
                dirs = 1:3;
            end
        end

        function setShapeDisplayMode(shapeAxesCell, rgbAxesCell, keyAx, isRgb)
        %%% Switch a shapeAxesGrid between its component and RGB axes sets
        %%% >> GuiCommon.setShapeDisplayMode(shapeAxesCell, rgbAxesCell, keyAx, isRgb)
        %%%
        %%% The set not in use is cleared and hidden: hiding an axes does not hide
        %%% its children, so cla() is what actually removes the old picture, and
        %%% its title (not a child) has to be blanked separately. The
        %%% key is drawn here (it does not depend on the mode shown); the panels
        %%% themselves are drawn by the caller's next redraw.
            if isRgb
                hideSet = shapeAxesCell(:);
            else
                hideSet = [rgbAxesCell(:); {keyAx}];
            end
            for k = 1:numel(hideSet)
                cla(hideSet{k});
                title(hideSet{k}, '');
                hideSet{k}.Visible = 'off';
            end
            if isRgb
                Vis.plotPolarisationKey(keyAx);
            end
        end
    end

    %%% ===================================================================
    %%% Stepper rows (label + selector + << >> buttons)
    %%% ===================================================================
    methods (Static)
        function [dropdown, prevBtn, nextBtn] = dropdownStepperRow(parent, rowNo, labelText, itemStrs, colWidths, onSelect, onPrev, onNext)
        %%% Label + index-valued dropdown + << / >> buttons, on one row
        %%% >> [dropdown, prevBtn, nextBtn] = GuiCommon.dropdownStepperRow(parent, rowNo, labelText, itemStrs, colWidths, onSelect, onPrev, onNext)
        %%%
        %%% ItemsData is always 1:numel(itemStrs), so .Value is the INDEX into
        %%% whatever axis the caller is stepping (angles/wavelengths/wavVecs/
        %%% omega0List), never the displayed string.
            row = uigridlayout(parent, [1 4]);
            row.Layout.Row = rowNo;
            row.ColumnWidth = colWidths;
            row.Padding = [0 0 0 0];

            uilabel(row, 'Text', labelText);
            dropdown = uidropdown(row, 'Items', itemStrs, 'ItemsData', 1:numel(itemStrs), ...
                'ValueChangedFcn', onSelect);
            prevBtn = uibutton(row, 'Text', '<<', 'ButtonPushedFcn', onPrev);
            nextBtn = uibutton(row, 'Text', '>>', 'ButtonPushedFcn', onNext);
        end

        function [field, prevBtn, nextBtn, trailLabel, rowLabel] = indexStepperRow(parent, rowNo, labelText, colWidths, onChange, onPrev, onNext)
        %%% Label + positive-integer field + << / >> buttons + a trailing
        %%% free-text label, on one row
        %%% >> [field, prevBtn, nextBtn, trailLabel, rowLabel] = GuiCommon.indexStepperRow(parent, rowNo, labelText, colWidths, onChange, onPrev, onNext)
        %%%
        %%% Limits [1 Inf] / RoundFractionalValues 'on' is the shared "this is a
        %%% 1-based index" setup. ModalShapeViewerGUI relaxes both at run time when
        %%% it re-purposes its field as a continuous velocity entry
        %%% (onSourceChanged) — a per-class override, not a second shape this
        %%% helper builds. rowLabel is returned for the same reason: a caller
        %%% that re-purposes the field also needs to re-word its label.
            row = uigridlayout(parent, [1 5]);
            row.Layout.Row = rowNo;
            row.ColumnWidth = colWidths;
            row.Padding = [0 0 0 0];

            rowLabel = uilabel(row, 'Text', labelText);
            field = uieditfield(row, 'numeric', 'Value', 1, 'Limits', [1 Inf], ...
                'RoundFractionalValues', 'on', 'ValueChangedFcn', onChange);
            prevBtn = uibutton(row, 'Text', '<<', 'ButtonPushedFcn', onPrev);
            nextBtn = uibutton(row, 'Text', '>>', 'ButtonPushedFcn', onNext);
            trailLabel = uilabel(row, 'Text', '');
        end
    end

    %%% ===================================================================
    %%% Axis-range controls
    %%% ===================================================================
    methods (Static)
        function [minField, maxField, label] = rangeFieldPair(parent, labelText, minVal, maxVal, callback)
        %%% A "<quantity> min/max" label over a side-by-side pair of text fields
        %%% >> [minField, maxField, label] = GuiCommon.rangeFieldPair(parent, labelText, minVal, maxVal, callback)
        %%%
        %%% INPUT: 1. parent    .. the side uigridlayout; consumes TWO of its rows
        %%%                        (label, then the field pair)
        %%%        2. labelText .. e.g. 'Angle [deg] min/max'
        %%%        3. minVal    .. initial min text ('' = blank/auto)
        %%%        4. maxVal    .. initial max text
        %%%        5. callback  .. ValueChangedFcn for BOTH fields (always the
        %%%                        owner's applyAxisRange)
        %%%
        %%% 'text', not 'numeric', deliberately: a blank field has to mean
        %%% "leave this axis alone", which a numeric field cannot express --
        %%% applyFieldLimit below reads that back as NaN and skips the axis.
            label = uilabel(parent, 'Text', labelText);
            gl = uigridlayout(parent, [1 2]);
            gl.Padding = [0 0 0 0];
            minField = uieditfield(gl, 'text', 'Value', minVal, 'ValueChangedFcn', callback);
            maxField = uieditfield(gl, 'text', 'Value', maxVal, 'ValueChangedFcn', callback);
        end

        function applyFieldLimit(ax, axisName, minField, maxField)
        %%% Apply one min/max field pair to an axes limit, ignoring blank/
        %%% unparseable/inverted input
        %%% >> GuiCommon.applyFieldLimit(ax, axisName, minField, maxField)
        %%%
        %%% INPUT: 1. ax       .. target axes
        %%%        2. axisName .. 'x' | 'y' | 'c' (colour)
        %%%        3/4. the field pair rangeFieldPair returned
        %%%
        %%% Silently doing nothing (rather than erroring or clamping) is the
        %%% established behaviour: these fields are edited character by
        %%% character, so a half-typed value must not disturb the plot.
            lo = str2double(minField.Value);
            hi = str2double(maxField.Value);
            if isnan(lo) || isnan(hi) || hi <= lo
                return
            end
            switch axisName
                case 'x'
                    xlim(ax, [lo, hi]);
                case 'y'
                    ylim(ax, [lo, hi]);
                case 'c'
                    clim(ax, [lo, hi]);
            end
        end

        function s = captureAxisRange(viewer)
        %%% Snapshot a viewer's axis-range field VALUES (not the handles)
        %%% >> s = GuiCommon.captureAxisRange(viewer)
        %%%
        %%% OUTPUT: s .. struct(xMin, xMax, yMin, yMax [, cMin, cMax])
        %%%
        %%% Values, not handles, because every caller uses this across a teardown:
        %%% the fields themselves are destroyed by the rebuild that follows, so
        %%% only the text survives. cMin/cMax appear only for viewers that have a
        %%% colour axis to control.
            s = struct( ...
                'xMin', viewer.xMinField.Value, 'xMax', viewer.xMaxField.Value, ...
                'yMin', viewer.yMinField.Value, 'yMax', viewer.yMaxField.Value);
            if isprop(viewer, 'cMinField')
                s.cMin = viewer.cMinField.Value;
                s.cMax = viewer.cMaxField.Value;
            end
        end

        function restoreAxisRange(viewer, s)
        %%% Write a captureAxisRange snapshot back into a viewer and apply it
        %%% >> GuiCommon.restoreAxisRange(viewer, s)
        %%%
        %%% A no-op for an empty snapshot (nothing was ever captured), so
        %%% callers need no guard of their own. Calls the viewer's own public
        %%% applyAxisRange() at the end -- restoring the text without applying
        %%% it would leave the fields disagreeing with the axes.
            if isempty(s)
                return
            end
            viewer.xMinField.Value = s.xMin; viewer.xMaxField.Value = s.xMax;
            viewer.yMinField.Value = s.yMin; viewer.yMaxField.Value = s.yMax;
            if isfield(s, 'cMin') && isprop(viewer, 'cMinField')
                viewer.cMinField.Value = s.cMin; viewer.cMaxField.Value = s.cMax;
            end
            viewer.applyAxisRange();
        end

        function setAxisRange(viewer, ax, axisName, lims)
        %%% Set one axis range of a viewer from code, as if typed into its fields
        %%% >> GuiCommon.setAxisRange(viewer, ax, axisName, lims)
        %%%
        %%% INPUT: 1. viewer   .. a viewer with x/y (and optionally c) Min/MaxField
        %%%                       pairs and a public applyAxisRange()
        %%%        2. ax       .. the axes those fields control (fills a kept
        %%%                       limit whose field is blank, see below)
        %%%        3. axisName .. 'x' | 'y' | 'c' (colour)
        %%%        4. lims     .. [lo hi]; a NaN entry keeps that limit, [] keeps both
        %%%
        %%% Writes the FIELDS, then calls the viewer's own applyAxisRange(), so the
        %%% fields never disagree with the axes and every later save/restore
        %%% (drill-down, mode switch) carries the value like a typed one. A kept
        %%% limit is the field's current text, or the axes' current limit where
        %%% that field is blank (DispersionViewerGUI starts blank). Unlike typing,
        %%% a call from code errors on bad input instead of ignoring it; error IDs
        %%% are <viewer class>:setAxisRange:<condition>.
            className = class(viewer);
            axisName = char(axisName);
            if ~ismember(axisName, {'x', 'y', 'c'})
                error([className ':setAxisRange:invalidAxis'], ...
                    'axisName must be ''x'', ''y'' or ''c''; got ''%s''.', axisName);
            end
            if isempty(lims)
                return
            end
            if ~isnumeric(lims) || numel(lims) ~= 2
                error([className ':setAxisRange:invalidLimits'], ...
                    'lims must be [lo hi] (NaN keeps that limit) or [].');
            end
            if isprop(viewer, 'isDrilledIn') && viewer.isDrilledIn
                error([className ':setAxisRange:drilledIn'], ...
                    ['The map''s range fields do not exist while drilled into a single angle ' ...
                     '-- return to the sweep view first.']);
            end
            if axisName == 'c' && (~isprop(viewer, 'cMinField') || strcmp(viewer.cMinField.Enable, 'off'))
                error([className ':setAxisRange:noColourAxis'], ...
                    'This view has no colour range to set.');
            end

            minField = viewer.([axisName 'MinField']);
            maxField = viewer.([axisName 'MaxField']);
            cur = [str2double(minField.Value), str2double(maxField.Value)];
            axLim = ax.([upper(axisName) 'Lim']);
            cur(isnan(cur)) = axLim(isnan(cur));
            new = double(lims(:).');
            new(isnan(new)) = cur(isnan(new));
            if new(2) <= new(1)
                error([className ':setAxisRange:invalidLimits'], ...
                    'Resulting %s range [%g %g] is empty or inverted.', axisName, new(1), new(2));
            end
            minField.Value = sprintf('%.10g', new(1));
            maxField.Value = sprintf('%.10g', new(2));
            viewer.applyAxisRange();
        end
    end

    %%% ===================================================================
    %%% Click handling / crosshair
    %%% ===================================================================
    methods (Static)
        function pt = clickPoint(ax, event)
        %%% Data-space (x, y) of a click, from either a real HitEvent or the axes
        %%% >> pt = GuiCommon.clickPoint(ax, event)
        %%%
        %%% event.IntersectionPoint exists on a real object-level hit (clicking a
        %%% plotted scatter/pcolor/line); an axes-level hit on empty background
        %%% carries no such field. isfield and isprop are checked separately
        %%% because isprop errors outright on a plain struct, which is what tests
        %%% pass here in place of a real HitEvent object.
            if isstruct(event)
                hasIntersection = isfield(event, 'IntersectionPoint') && ~isempty(event.IntersectionPoint);
            else
                hasIntersection = isprop(event, 'IntersectionPoint') && ~isempty(event.IntersectionPoint);
            end
            if hasIntersection
                pt = event.IntersectionPoint(1:2);
            else
                cp = ax.CurrentPoint;
                pt = cp(1, 1:2);
            end
        end

        function h = markCrosshair(ax, h, xVal, yVal)
        %%% Move (or create) the selection crosshair at (xVal, yVal)
        %%% >> h = GuiCommon.markCrosshair(ax, h, xVal, yVal)
        %%%
        %%% INPUT: 2. h .. the existing crosshair handle, or gobjects(0)/[] --
        %%%                a stale handle from an axes that has since been
        %%%                destroyed also counts as "none" (~isgraphics)
        %%% OUTPUT: h .. the live crosshair handle, to store back
        %%%
        %%% Kept out of the legend. HitTest 'off' so the marker is never itself the click target --
        %%% otherwise clicking an already-marked point would hit the crosshair
        %%% instead of the data underneath it.
            if isempty(h) || ~isgraphics(h)
                hold(ax, 'on');
                %%% compact but heavy-stroked, so it stays visible on busy maps
                h = plot(ax, xVal, yVal, '+', 'Color', GuiCommon.CROSSHAIR_COLOR, ...
                    'MarkerSize', 10, 'LineWidth', 3.2, 'HitTest', 'off');
                %%% no legend entry (it showed as 'data1' on a map with a legend);
                %%% not HandleVisibility 'off', which cla() would then leave behind
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
                hold(ax, 'off');
            else
                set(h, 'XData', xVal, 'YData', yVal);
            end
        end

        function [bestIdx, bestLine] = nearestLinePoint(ax, hLines, pt, defaultIdx, defaultLine)
        %%% Nearest plotted point across a set of line handles, in
        %%% axis-range-normalised distance
        %%% >> [bestIdx, bestLine] = GuiCommon.nearestLinePoint(ax, hLines, pt, defaultIdx, defaultLine)
        %%%
        %%% INPUT: 3. pt .. the clicked (x, y), from clickPoint above
        %%%        4/5. what to return when nothing is hittable (the caller's
        %%%             own current selection, so a miss changes nothing)
        %%% OUTPUT: bestIdx  .. index ALONG the winning line's data
        %%%         bestLine .. which line (index into hLines) won
        %%%
        %%% Normalised by the axes' own displayed spans, so wildly different x/y
        %%% scales (wavelength in microns vs. velocity in km/s) do not bias the
        %%% match toward whichever axis spans more raw units.
        %%%
        %%% Skips HIDDEN lines, not just missing ones: a de-selected series
        %%% (Vis.plotDispersionCurves' showRitzBranches=false) leaves real line
        %%% objects behind with Visible='off', and a click near where a hidden
        %%% branch would be must not snap to it — it cannot be seen to aim at.
            xRange = max(diff(xlim(ax)), eps);
            yRange = max(diff(ylim(ax)), eps);

            bestDist = Inf; bestIdx = defaultIdx; bestLine = defaultLine;
            for im = 1:numel(hLines)
                if ~isgraphics(hLines(im)); continue; end
                if strcmp(hLines(im).Visible, 'off'); continue; end
                dx = (hLines(im).XData - pt(1)) / xRange;
                dy = (hLines(im).YData - pt(2)) / yRange;
                [dMin, iMin] = min(dx.^2 + dy.^2);
                if dMin < bestDist
                    bestDist = dMin; bestIdx = iMin; bestLine = im;
                end
            end
        end
    end

    %%% ===================================================================
    %%% Per-type visibility checkboxes
    %%% ===================================================================
    methods (Static)
        function checkboxes = rebuildTypeCheckboxes(panel, hPoints, visState, defaultOn, onToggle, showHeader)
        %%% (Re)build a 'Show types:' checkbox list from a per-type handle array
        %%% >> checkboxes = GuiCommon.rebuildTypeCheckboxes(panel, hPoints, visState, defaultOn, onToggle, showHeader)
        %%%
        %%% INPUT: 1. panel     .. uipanel whose children are replaced wholesale
        %%%        2. hPoints   .. the per-type scatter handles a classification
        %%%                        plot returned under colorBy='type' (one handle
        %%%                        per type actually present, DisplayName=type
        %%%                        name), so the checkbox list always matches
        %%%                        what's really on the plot rather than a
        %%%                        hardcoded type list. Pass gobjects(0) to clear
        %%%                        the panel (after a failed draw, or no
        %%%                        classified points to show).
        %%%  (opt) 3. visState  .. struct of remembered per-type on/off states,
        %%%                        field = type name (default struct()). A type
        %%%                        with no field gets defaultOn
        %%%  (opt) 4. defaultOn .. initial state of a type not in visState
        %%%                        (default true)
        %%%  (opt) 5. onToggle  .. @(typeName, isChecked), called after every
        %%%                        user toggle so the viewer can store the new
        %%%                        state in its own visState (default [])
        %%%  (opt) 6. showHeader .. true (default): a bold 'Show types:' first row;
        %%%                        false: the caller shows its own (typesHeaderRow)
        %%%
        %%% Each checkbox toggles its own handle's Visible directly -- no
        %%% redraw needed, and the legend (built off the same handles) is
        %%% unaffected by a Visible toggle. The initial state is applied to the
        %%% handles here too, so checkbox and plot never disagree. Each row
        %%% carries a colour swatch of its series' own marker colour, since an
        %%% overlay (styleClassificationOverlay) has no legend to explain it.
            if nargin < 3 || isempty(visState);  visState  = struct(); end
            if nargin < 4 || isempty(defaultOn); defaultOn = true;     end
            if nargin < 5; onToggle = []; end
            if nargin < 6 || isempty(showHeader); showHeader = true; end
            delete(panel.Children);
            checkboxes = gobjects(0);
            if isempty(hPoints) || ~all(isgraphics(hPoints))
                return
            end
            n = numel(hPoints);
            nHdr = double(showHeader);
            gl = uigridlayout(panel, [n + nHdr, 2]);
            gl.RowHeight = repmat({18}, 1, n + nHdr);
            gl.ColumnWidth = {12, '1x'};
            gl.Padding = [0 0 0 0];
            gl.RowSpacing = 2;
            gl.Scrollable = 'on';   % never clip a row off a short side panel

            if showHeader
                hdr = uilabel(gl, 'Text', 'Show types:', 'FontWeight', 'bold');
                hdr.Layout.Row = 1; hdr.Layout.Column = [1 2];
            end

            checkboxes = gobjects(1, n);
            for k = 1:n
                typeName = hPoints(k).DisplayName;
                if isfield(visState, typeName)
                    isOn = logical(visState.(typeName));
                else
                    isOn = logical(defaultOn);
                end
                GuiCommon.setSeriesVisible(hPoints(k), isOn);

                swatchColour = hPoints(k).CData;
                if ~(isnumeric(swatchColour) && isequal(size(swatchColour), [1 3]))
                    swatchColour = [1 1 1];   % per-point colours -- no single swatch
                end
                swatch = uilabel(gl, 'Text', '', 'BackgroundColor', swatchColour);
                swatch.Layout.Row = k + nHdr; swatch.Layout.Column = 1;

                cb = uicheckbox(gl, 'Text', typeName, 'Value', isOn, ...
                    'ValueChangedFcn', @(src, ~) GuiCommon.onTypeCheckboxChanged( ...
                        hPoints(k), typeName, src.Value, onToggle));
                cb.Layout.Row = k + nHdr; cb.Layout.Column = 2;
                checkboxes(k) = cb;
            end
        end

        function onTypeCheckboxChanged(hPoint, typeName, isChecked, onToggle)
        %%% One type checkbox's ValueChangedFcn: show/hide the series, then let
        %%% the owning viewer remember the new state
        %%% >> GuiCommon.onTypeCheckboxChanged(hPoint, typeName, isChecked, onToggle)
            GuiCommon.setSeriesVisible(hPoint, isChecked);
            if ~isempty(onToggle)
                onToggle(typeName, isChecked);
            end
        end

        function styleClassificationOverlay(ax, hPoints, hLim, savedLabels)
        %%% Turn a just-drawn classification plot (plotClassificationDiagnostic/
        %%% plotClassificationCurve, showGreenBackground=false) into an overlay
        %%% on top of another graph type already on ax
        %%% >> GuiCommon.styleClassificationOverlay(ax, hPoints, hLim, savedLabels)
        %%%
        %%% INPUT: 1. ax          .. the axes both plots share
        %%%        2. hPoints     .. the classification plot's per-type scatter handles
        %%%        3. hLim        .. its limiting-velocity reference-line struct
        %%%                          ('long'/'shear' handles or []), deleted here:
        %%%                          only the type series have checkboxes
        %%%        4. savedLabels .. the underlying plot's labels, from
        %%%                          GuiCommon.captureAxesLabels(ax) BEFORE the
        %%%                          classification plot overwrote them
        %%%
        %%% Markers become hollow rings (edge in the type colour, no fill) so the
        %%% underlying Ritz shading / Green heatmap stays readable inside them, and
        %%% non-pickable so a click passes through to the underlying plot's own
        %%% ButtonDownFcn. The legend is switched off: the underlying plot has no
        %%% legend entry of its own, and the checkbox swatches name the colours.
        %%% Weak entries (UserData.isWeak, see the Vis classification plots), already
        %%% open there, stay apart here as FAINT rings (edge alpha
        %%% OVERLAY_WEAK_ALPHA), since every overlay ring is open.
            set(hPoints, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'flat', ...
                'LineWidth', 1.5, 'PickableParts', 'none');
            for k = 1:numel(hPoints)
                ud = hPoints(k).UserData;
                if isstruct(ud) && isfield(ud, 'isWeak') && any(ud.isWeak)
                    hPoints(k).AlphaData = 1 - (1 - GuiCommon.OVERLAY_WEAK_ALPHA) * double(ud.isWeak(:));
                    hPoints(k).AlphaDataMapping = 'none';
                    hPoints(k).MarkerEdgeAlpha = 'flat';
                end
            end
            limHandles = [hLim.long(:); hLim.shear(:)];
            delete(limHandles(isgraphics(limHandles)));
            legend(ax, 'off');
            textObjs = {ax.XLabel, ax.YLabel, ax.Title};
            for k = 1:3
                textObjs{k}.Interpreter = savedLabels(k).Interpreter;
                textObjs{k}.String      = savedLabels(k).String;
            end
        end

        function savedLabels = captureAxesLabels(ax)
        %%% x-label, y-label and title String+Interpreter of ax, in that order --
        %%% for styleClassificationOverlay to restore after an overlay draw
        %%% >> savedLabels = GuiCommon.captureAxesLabels(ax)
            textObjs = {ax.XLabel, ax.YLabel, ax.Title};
            savedLabels = struct('String', cell(1, 3), 'Interpreter', cell(1, 3));
            for k = 1:3
                savedLabels(k).String      = textObjs{k}.String;
                savedLabels(k).Interpreter = textObjs{k}.Interpreter;
            end
        end

        function setSeriesVisible(hPoint, isChecked)
        %%% Show/hide one plotted series without redrawing anything else
        %%% >> GuiCommon.setSeriesVisible(hPoint, isChecked)
            if isChecked
                hPoint.Visible = 'on';
            else
                hPoint.Visible = 'off';
            end
        end

        function [hdrLabel, exportBtn] = typesHeaderRow(sideGl, rowNo, onExport)
        %%% The side panel's 'Show types:' header with an 'Export figure' button to
        %%% its right -- one fixed row above the type-checkbox panel, which the
        %%% viewer then fills with rebuildTypeCheckboxes(..., showHeader=false)
        %%% >> [hdrLabel, exportBtn] = GuiCommon.typesHeaderRow(sideGl, rowNo, onExport)
            row = uigridlayout(sideGl, [1 2]);
            row.Layout.Row = rowNo;
            row.ColumnWidth = {'1x', 100};
            row.Padding = [0 0 0 0];
            hdrLabel = uilabel(row, 'Text', 'Show types:', 'FontWeight', 'bold');
            exportBtn = uibutton(row, 'Text', 'Export figure', ...
                'Tooltip', 'Open the shown plot in a separate figure window, to save in any format', ...
                'ButtonPushedFcn', @(~,~) onExport());
        end

        function checkboxes = classificationTypeCheckboxes(panel, isClassificationView, hOut, ...
                drawOverlay, diagnosticState, overlayState, onToggle)
        %%% The 'Show types:' checkboxes after a map/curve draw: the classification
        %%% view's own per-type series, or the classification overlay on any other view
        %%% >> checkboxes = GuiCommon.classificationTypeCheckboxes(panel, isClassificationView, hOut, drawOverlay, diagnosticState, overlayState, onToggle)
        %%%
        %%% INPUT: 1. panel                .. the viewer's type-checkbox uipanel
        %%%        2. isClassificationView .. true: hOut is the classification plot's
        %%%                                  per-type series (all on by default)
        %%%        3. hOut                 .. what the view's plot function returned
        %%%        4. drawOverlay          .. @() returning the overlay's per-type handles,
        %%%                                  called only on another view with a plot
        %%%                                  drawn (overlay off by default)
        %%%        5-6. diagnosticState, overlayState .. the viewer's remembered
        %%%                                  per-type states of the two contexts
        %%%        7. onToggle             .. @(stateProp, typeName, isOn), stateProp
        %%%                                  'diagnosticTypeVisible' | 'overlayTypeVisible'
        %%% OUTPUT: checkboxes .. rebuildTypeCheckboxes' output (the header is the
        %%%                       viewer's own typesHeaderRow)
            if isClassificationView
                checkboxes = GuiCommon.rebuildTypeCheckboxes(panel, hOut, diagnosticState, true, ...
                    @(typeName, isOn) onToggle('diagnosticTypeVisible', typeName, isOn), false);
            elseif ~isempty(hOut)
                hOverlay = drawOverlay();
                checkboxes = GuiCommon.rebuildTypeCheckboxes(panel, hOverlay, overlayState, false, ...
                    @(typeName, isOn) onToggle('overlayTypeVisible', typeName, isOn), false);
            else
                checkboxes = GuiCommon.rebuildTypeCheckboxes(panel, gobjects(0));
            end
        end

        function clearAxes(ax)
        %%% Clear a viewer's plot axes before the next graph type is drawn: cla(),
        %%% plus what cla() leaves behind -- objects drawn with HandleVisibility
        %%% 'off' (the freqDisp veto line) and the labels' interpreter (a plot that
        %%% set Interpreter='none', the Ritz-mode curve/map's criterion label, would
        %%% leave the next graph type's TeX label printed raw)
        %%% >> GuiCommon.clearAxes(ax)
            cla(ax);
            delete(allchild(ax));
            ax.XLabel.Interpreter = 'tex';
            ax.YLabel.Interpreter = 'tex';
            ax.Title.Interpreter  = 'tex';
        end

        function fig = exportAxes(ax)
        %%% Copy a viewer's plot into a new ordinary figure window, as shown: every
        %%% plotted object (hidden types stay hidden), limits, scales, colours,
        %%% colormap, labels, title, legend and colour bar. The copy is static (no
        %%% click callbacks) and is saved from its own File menu or with
        %%% exportgraphics/saveas, in any format.
        %%% >> fig = GuiCommon.exportAxes(ax)
        %%%
        %%% A uiaxes cannot itself be parented to an ordinary figure, so its
        %%% children are copied into a fresh axes and the axes-level settings are
        %%% copied one by one.
            figName = ax.Title.String;
            if iscell(figName); figName = strjoin(figName, ' '); end
            if isempty(figName); figName = 'plot'; end
            fig = figure('Color', 'w', 'Name', ['GuliweRR -- ' char(figName)]);
            newAx = axes(fig);
            %%% allchild, not Children: lines drawn with HandleVisibility 'off' (a
            %%% veto level, a connecting line) are part of the plot too. copyobj of
            %%% the whole list keeps the stacking order (first = on top)
            kids = allchild(ax);
            newKids = copyobj(kids, newAx);
            set(findall(newAx, '-property', 'ButtonDownFcn'), 'ButtonDownFcn', '');

            props = {'XScale', 'YScale', 'XDir', 'YDir', 'Color', 'Box', 'XGrid', 'YGrid', ...
                'Layer', 'View'};
            for k = 1:numel(props)
                newAx.(props{k}) = ax.(props{k});
            end
            %%% limits only where the viewer set them (range fields, a plot's own
            %%% clim); auto limits are recomputed by the copy itself -- a uiaxes not
            %%% yet rendered can still report the previous plot's auto limits
            limNames = {'XLim', 'YLim', 'ZLim', 'CLim'};
            for k = 1:numel(limNames)
                if strcmp(ax.([limNames{k} 'Mode']), 'manual')
                    newAx.(limNames{k}) = ax.(limNames{k});
                end
            end
            colormap(newAx, colormap(ax));
            textNames = {'XLabel', 'YLabel', 'Title'};
            for k = 1:3
                newAx.(textNames{k}).String      = ax.(textNames{k}).String;
                newAx.(textNames{k}).Interpreter = ax.(textNames{k}).Interpreter;
            end

            lg = ax.Legend;
            if ~isempty(lg) && isvalid(lg) && strcmp(lg.Visible, 'on')
                [~, idx] = ismember(lg.PlotChildren, kids);
                idx = idx(idx > 0);
                if ~isempty(idx)
                    legend(newAx, newKids(idx), 'Location', lg.Location, ...
                        'Interpreter', lg.Interpreter);
                end
            end
            cb = ax.Colorbar;
            if ~isempty(cb) && isvalid(cb) && strcmp(cb.Visible, 'on')
                newCb = colorbar(newAx, 'Location', cb.Location);
                newCb.Label.String      = cb.Label.String;
                newCb.Label.Interpreter = cb.Label.Interpreter;
                if strcmp(cb.TicksMode, 'manual'); newCb.Ticks = cb.Ticks; end
                if strcmp(cb.TickLabelsMode, 'manual'); newCb.TickLabels = cb.TickLabels; end
            end
        end
    end

    %%% ===================================================================
    %%% Result lookups shared across viewer families
    %%% ===================================================================
    methods (Static)
        function typeStr = significantModeType(significantModes, modeNo)
        %%% Classification type of a Ritz mode, if it appears in a
        %%% significantModes entry list
        %%% >> typeStr = GuiCommon.significantModeType(significantModes, modeNo)
        %%%
        %%% INPUT: 1. significantModes .. one ModeClassification entry array --
        %%%                               sweepResult.significantModesMap{angNo},
        %%%                               or singleResult.classification.significantModes
        %%% OUTPUT: typeStr .. the entry's own .type, or 'unclassified' when the
        %%%                    mode is not a significant one at all
        %%%
        %%% 'unclassified' rather than an error: every viewer here lets the user
        %%% click/step onto ANY computed Ritz mode, most of which are legitimately
        %%% not significant modes.
            idx = find([significantModes.modeNo] == modeNo, 1);
            if isempty(idx)
                typeStr = 'unclassified';
            else
                typeStr = significantModes(idx).type;
            end
        end

        function dd = ritzCriterionDropdown(parent, isLayered, onChange)
        %%% The Ritz-mode criterion dropdown of the HS/layered viewers
        %%% >> dd = GuiCommon.ritzCriterionDropdown(parent, isLayered, onChange)
        %%%
        %%% INPUT: 1. parent    .. container
        %%%        2. isLayered .. true for a layered result: freqDisp is greyed
        %%%                        (still selectable, with ritzCriterionNote's note)
        %%%        3. onChange  .. ValueChangedFcn
        %%% OUTPUT: dd .. uidropdown, Items 'Ritz score' | 'freqDisp' |
        %%%               'depthInsensitivity', ItemsData = the Vis criterion names,
        %%%               built hidden (shown with the Ritz-mode graph type)
            dd = uidropdown(parent, ...
                'Items',     {'Ritz score', 'freqDisp', 'depthInsensitivity'}, ...
                'ItemsData', {'ritzScore', 'freqDisp', 'depthInsensitivity'}, ...
                'Visible', 'off', 'ValueChangedFcn', onChange, ...
                'Tooltip', ['Ritz score: the score of the Ritz candidate each mode supports; ' ...
                            'freqDisp: the half-space veto quantity; depthInsensitivity: ' ...
                            'how little the mode''s velocity changes with the Ritz domain ' ...
                            'depth (large for surface modes)']);
            if isLayered
                GuiCommon.greyDropdownItem(dd, 'freqDisp');
            end
        end

        function greyDropdownItem(dd, itemData)
        %%% Show one dropdown item greyed (italic grey uistyle) -- uidropdown has no
        %%% per-item Enable, so an item that stays selectable but is not in use is
        %%% marked this way
        %%% >> GuiCommon.greyDropdownItem(dd, itemData)
        %%%
        %%% INPUT: 1. dd       .. uidropdown with ItemsData
        %%%        2. itemData .. the item's ItemsData value (its position is looked up,
        %%%                       so a reordered list cannot grey the wrong item)
            idx = find(strcmp(dd.ItemsData, itemData), 1);
            addStyle(dd, uistyle('FontColor', GuiCommon.GREY_ITEM_COLOR, 'FontAngle', 'italic'), ...
                'item', idx);
        end

        function note = ritzCriterionNote(criterion, isLayered)
        %%% Status note for a Ritz-mode criterion: freqDisp on a layered result
        %%% >> note = GuiCommon.ritzCriterionNote(criterion, isLayered)
        %%%
        %%% OUTPUT: note .. '' unless criterion is 'freqDisp' on a layered result
            note = '';
            if isLayered && strcmp(criterion, 'freqDisp')
                note = ['Note: freqDisp is not used on layered stacks (the film makes ' ...
                        'true guided modes dispersive).'];
            end
        end

        function lines = selectionLines(angleVal, modeNo, velSel, typeStr)
        %%% The four lines of a click's selection readout: angle, mode #, velocity, type
        %%% >> lines = GuiCommon.selectionLines(angleVal, modeNo, velSel, typeStr)
        %%%
        %%% OUTPUT: lines .. 1x4 cell; the viewers append selectionExtras' lines
            lines = { ...
                sprintf('angle = %.1f deg', angleVal), ...
                sprintf('mode # = %d', modeNo), ...
                sprintf('velocity = %.4f km/s', velSel), ...
                sprintf('type = %s', typeStr)};
        end

        function msg = statusText(selectionInfo, showSelection, criterion, isLayered)
        %%% Status-line text of a map/curve viewer after a draw or click
        %%% >> msg = GuiCommon.statusText(selectionInfo, showSelection, criterion, isLayered)
        %%%
        %%% INPUT: 1. selectionInfo .. selectionLines + selectionExtras of the last
        %%%                            click, {} before any
        %%%        2. showSelection .. true: all of it (standalone viewer); false: only
        %%%                            the lines after the first four (classification
        %%%                            and criterion), since an embedding host shows the
        %%%                            selection itself
        %%%        3. criterion     .. the Ritz-mode criterion shown, or ''
        %%%        4. isLayered     .. true for a layered result (ritzCriterionNote)
        %%% OUTPUT: msg .. the parts joined by ',  ', '' when there is nothing to say
            parts = {};
            if showSelection
                parts = selectionInfo;
            elseif numel(selectionInfo) > 4
                parts = selectionInfo(5:end);
            end
            note = GuiCommon.ritzCriterionNote(criterion, isLayered);
            if ~isempty(note); parts{end+1} = note; end
            msg = strjoin(parts, ',  ');
        end

        function lines = selectionExtras(singleResult, significantModes, modeNo, criterion)
        %%% Extra readout lines of a click: classification readout and criterion value
        %%% >> lines = GuiCommon.selectionExtras(singleResult, significantModes, modeNo, criterion)
        %%%
        %%% INPUT: 1. singleResult     .. CalcSingleHS/CalcSingleLayered of the angle
        %%%        2. significantModes .. its entries
        %%%        3. modeNo           .. the selected Ritz mode
        %%%        4. criterion        .. the Ritz-mode criterion shown, or '' when
        %%%                               another graph type is shown
        %%% OUTPUT: lines .. row cell: classificationReadout (if the mode is a scored
        %%%                  entry), then the criterion value -- 'Ritz score x dB' /
        %%%                  'freqDisp = x', with "(vetoed)" below RitzScore.MIN_FREQDISP
        %%%                  on a half-space / 'depthInsensitivity = x'
            lines = {};
            txt = GuiCommon.classificationReadout(significantModes, modeNo);
            if ~isempty(txt); lines{end+1} = txt; end
            ritz = singleResult.ritzSolver;
            switch criterion
                case 'ritzScore'
                    sc = Vis.ritzModeScores(singleResult);
                    if isnan(sc(modeNo))
                        lines{end+1} = 'Ritz score: supports no Ritz candidate';
                    else
                        lines{end+1} = sprintf('Ritz score %.1f dB', 20 * sc(modeNo));
                    end
                case 'freqDisp'
                    fd = ritz.diagnostics.freqDisp(modeNo);
                    lines{end+1} = sprintf('freqDisp = %.3f', fd);
                    if ~isa(singleResult, 'CalcSingleLayered') && fd < RitzScore.MIN_FREQDISP
                        lines{end} = sprintf('%s (vetoed: < %.2g)', lines{end}, RitzScore.MIN_FREQDISP);
                    end
                case 'depthInsensitivity'
                    lines{end+1} = sprintf('depthInsensitivity = %.3f', ritz.diagnostics.depthInsensitivity(modeNo));
                otherwise
                    %%% '' -- another graph type is shown: no criterion line
            end
        end

        function txt = classificationReadout(significantModes, modeNo)
        %%% "<method> score <x.x> dB (<grade>), supporting mode <n>" of a classified
        %%% Ritz mode -- the crosshair readout of the viewers
        %%% >> txt = GuiCommon.classificationReadout(significantModes, modeNo)
        %%%
        %%% INPUT: 1. significantModes .. one ModeClassification entry array
        %%%        2. modeNo           .. the selected Ritz mode
        %%% OUTPUT: txt .. '' when the mode is no scored entry (unclassified, or a
        %%%                'limiting' entry, which has no score)
        %%%
        %%% The supporting mode is the entry's own modeNo: the Ritz mode the
        %%% classification attributes the wave to, which a click snaps to.
            idx = find([significantModes.modeNo] == modeNo, 1);
            txt = '';
            if isempty(idx) || isnan(significantModes(idx).score)
                return
            end
            e = significantModes(idx);
            txt = sprintf('%s score %.1f dB (%s), supporting mode %d', ...
                e.method, e.scoreDb, e.grade, e.modeNo);
        end

        function txt = plateModeTypeText(ritz, modeNo)
        %%% "<type name> (confidence <x.xx>)" for one Lamb/plate mode
        %%% >> txt = GuiCommon.plateModeTypeText(ritz, modeNo)
        %%%
        %%% INPUT: 1. ritz   .. a calc()'d RitzSolverLamb
        %%%        2. modeNo .. mode index into its modesType/modeTypeConfidence
        %%%
        %%% Names come from Vis.plateModeTypeStyle, the one source for the plate
        %%% mode-type visual vocabulary, so a viewer's own label always agrees
        %%% with the plot it sits next to.
            typeNames = Vis.plateModeTypeStyle();
            txt = sprintf('%s (confidence %.2f)', ...
                typeNames{ritz.modesType(modeNo)}, ritz.modeTypeConfidence(modeNo));
        end
    end

    %%% ===================================================================
    %%% Constructor-argument validation
    %%% ===================================================================
    methods (Static)
        function validateWindowMode(className, windowMode, allowedModes)
        %%% Reject a windowMode outside the caller's own allowed set
        %%% >> GuiCommon.validateWindowMode(className, windowMode, allowedModes)
        %%%
        %%% Throws <className>:constructor:invalidWindowMode, built from className so the
        %%% identifier stays greppable at the call site. The message comes from
        %%% allowedModes itself (CalcCommon.orList), so the two cannot drift apart.
            if ~ismember(windowMode, allowedModes)
                error([className ':constructor:invalidWindowMode'], ...
                    'windowMode must be %s (got ''%s'').', ...
                    CalcCommon.orList(allowedModes), windowMode);
            end
        end

        function requireStandaloneForSplit(className, windowMode, parentArg, parentArgNames)
        %%% Reject windowMode='split' together with a caller-supplied parent
        %%% >> GuiCommon.requireStandaloneForSplit(className, windowMode, parentArg, parentArgNames)
        %%%
        %%% Splitting only makes sense when the class owns its own window(s);
        %%% once a caller hosts it, there is no second window to split into.
        %%% Throws <className>:constructor:splitRequiresStandalone.
            if ~isempty(parentArg) && strcmp(windowMode, 'split')
                error([className ':constructor:splitRequiresStandalone'], ...
                    '''split'' windowMode only applies to standalone construction (omit %s).', ...
                    parentArgNames);
            end
        end

        function validateEmbedParents(className, windowMode, controlsParent, shapeParent)
        %%% Validate the two-region embedding contract of the modal-shape family
        %%% >> GuiCommon.validateEmbedParents(className, windowMode, controlsParent, shapeParent)
        %%%
        %%% controlsParent and shapeParent are both-or-neither (a control panel
        %%% with nowhere to draw shapes, or the reverse, is never meaningful),
        %%% and embedding excludes 'split'. Throws
        %%% <className>:constructor:parentArgsMismatch / :splitRequiresStandalone.
            if xor(isempty(controlsParent), isempty(shapeParent))
                error([className ':constructor:parentArgsMismatch'], ...
                    'controlsParent and shapeParent must be given together, or both omitted.');
            end
            GuiCommon.requireStandaloneForSplit(className, windowMode, controlsParent, ...
                'controlsParent/shapeParent');
        end
    end
end
