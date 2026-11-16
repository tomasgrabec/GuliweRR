function [win, calcResult, meta] = viewSaved(filePath, windowMode, showNow)
%%% Load a saved Calc* result and open its default viewer -- no recalculation
%%% >> [win, calcResult, meta] = Gui.viewSaved(filePath, windowMode, showNow)
%%%
%%% INPUT:  (opt) 1. filePath   .. .mat file written by Calc.saveResult (or a
%%%                                plain save() of one Calc* result); omitted or
%%%                                empty opens a file dialog in <project>/results/
%%%         (opt) 2. windowMode .. as Gui.viewResult (default = [])
%%%         (opt) 3. showNow    .. as Gui.viewResult (default = true)
%%% OUTPUT: 1. win        .. the viewer object ([] if the dialog was cancelled)
%%%         2. calcResult .. the loaded result
%%%         3. meta       .. its save metadata (see CalcCommon.saveCalc)
%%%
%%% Typical use: a batch run (no display) ends with Calc.saveResult(...); the
%%% file is then opened here from an interactive MATLAB session.

    if nargin < 2; windowMode = []; end
    if nargin < 3; showNow = []; end
    win = []; calcResult = []; meta = [];

    if nargin < 1 || isempty(filePath)
        startDir = fullfile(CalcCommon.projectRoot(), 'results');
        if ~isfolder(startDir); startDir = pwd; end
        [fileName, folder] = uigetfile(fullfile(startDir, '*.mat'), 'Open a saved GuliweRR result');
        if isequal(fileName, 0)
            return   % cancelled
        end
        filePath = fullfile(folder, fileName);
    end

    [calcResult, meta] = CalcCommon.loadCalc(filePath);
    if ~isempty(meta.note)
        fprintf('%s (%s, saved %s): %s\n', meta.className, filePath, meta.savedAt, meta.note);
    end
    win = Gui.viewResult(calcResult, windowMode, showNow);
end
