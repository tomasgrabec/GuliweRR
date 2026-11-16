function filePath = saveResult(varargin)
%%% Thin Calc namespace launcher for CalcCommon.saveCalc -- save a finished
%%% Calc* result (plus provenance metadata) to a .mat file
%%% >> filePath = Calc.saveResult(calcResult, filePath, note)
%%%
%%% filePath defaults to an auto-named file in <project>/results/. View the
%%% file later, in any MATLAB session, with Gui.viewSaved(filePath).
%%% See CalcCommon.saveCalc.

    filePath = CalcCommon.saveCalc(varargin{:});
end
