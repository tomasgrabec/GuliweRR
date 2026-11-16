function tf = visIsLayered(result)
%%% Whether a result is a layered stack (CalcSingleLayered / CalcAngSweepLayered)
%%% >> tf = visIsLayered(result)
%%%
%%% @Vis-private helper: the one layered test of the classification and Ritz-mode
%%% plots, which pick a method's score scale, the Ritz score's polarisation term and
%%% the limiting-velocity reference line by it.
%%%
%%% INPUT:  1. result .. CalcSingleHS / CalcSingleLayered / CalcAngSweepHS /
%%%                      CalcAngSweepLayered
%%% OUTPUT: tf .. true for CalcSingleLayered and CalcAngSweepLayered
    tf = isa(result, 'CalcSingleLayered') || isa(result, 'CalcAngSweepLayered');
end
