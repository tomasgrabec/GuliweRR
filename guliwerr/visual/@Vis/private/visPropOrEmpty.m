function val = visPropOrEmpty(obj, names)
%%% Read the first property in a preference list that the object actually has
%%% >> val = visPropOrEmpty(obj, names)
%%%
%%% @Vis-private helper for the duck-typing convention documented in
%%% docs/reference/architecture/visual.md ("Duck-typing across HS/Layered result
%%% classes"): several result properties carry the same physical meaning under
%%% different names on the HS and Layered classes (longVelMap vs
%%% subLongVelMap, longLimVel vs subLongLimVel), and one may simply be absent.
%%% Each caller had its own isprop/elseif/else ladder saying so.
%%%
%%% INPUT: 1. obj   .. the result object to read from
%%%        2. names .. cell of property names, most-preferred first
%%% OUTPUT: val .. the first present property's value, or [] if none exists
%%%
%%% NOTE the two families' names mean genuinely different physical quantities
%%% (HS's system-wide longitudinal limiting velocity vs. Layered's
%%% SUBSTRATE-only one) -- they are interchangeable for the display purposes
%%% that use this helper (a velocity-axis cutoff), not in general. A caller
%%% that needs to KNOW which family it has should branch on isprop itself
%%% rather than call this.

    val = [];
    for k = 1:numel(names)
        if isprop(obj, names{k})
            val = obj.(names{k});
            return
        end
    end
end
