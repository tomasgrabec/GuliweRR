function typeAllowed = visClassificationTypesAllowed(refLevel)
%%% Which classified mode types a refLevel tier displays
%%% >> typeAllowed = visClassificationTypesAllowed(refLevel)
%%%
%%% @Vis-private helper: one source for the refLevel filtering rule, shared by
%%% Vis.plotClassificationDiagnostic (swept) and Vis.plotClassificationCurve
%%% (single angle). These two are meant to show the same modes at the same
%%% tier -- keeping the rule in one place is what actually guarantees that,
%%% rather than two lists that have to be edited together by hand.
%%%
%%% INPUT: 1. refLevel .. 'sawOnly' | 'sawLimit' | 'full'
%%% OUTPUT: typeAllowed .. cell of type names to keep
%%%
%%% 'other'-type points (unclassified/rejected candidates) are never shown at
%%% any tier -- these plots are a classification READ-OUT, not a catalogue of
%%% every candidate considered. 'limiting' points (Ritz modes sitting AT a
%%% bulk-wave limiting velocity) appear at 'full' only; 'sawOnly' and
%%% 'sawLimit' differ from each other in their reference LINE, not in which
%%% points they draw, so they share a list here.

    if strcmp(refLevel, 'full')
        typeAllowed = {'Rayleigh', 'pseudoRayleigh', 'SAW', 'limiting'};
    else
        typeAllowed = {'Rayleigh', 'pseudoRayleigh', 'SAW'};
    end
end
