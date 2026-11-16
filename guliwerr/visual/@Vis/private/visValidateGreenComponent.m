function visValidateGreenComponent(component, fnName)
%%% Reject a Green-component argument that isn't 'GIJ' / 'GIJeff' (I, J in 1..3) or 'spot'
%%% >> visValidateGreenComponent(component, fnName)
%%%
%%% @Vis-private helper, the Green-component counterpart of visValidateChoice:
%%% the allowed set is a pattern (CalcCommon.parseGreenComponent), not a list.
%%% Error identifier 'Vis:<fnName>:badComponent', as visValidateChoice would give.
    if ~(ischar(component) || (isstring(component) && isscalar(component))) || ...
            isempty(regexp(char(component), '^(G[123][123](eff)?|spot)$', 'once'))
        error(['Vis:' fnName ':badComponent'], ...
            ['component must be ''GIJ'' or ''GIJeff'' with I, J in 1..3 -- e.g. ''G13eff'' ' ...
             '(stored, default), ''G13'' (one-sided G13_+), ''G33'', or ''spot'' (the stored ' ...
             'spot-averaged map) (got ''%s'').'], char(component));
    end
end
