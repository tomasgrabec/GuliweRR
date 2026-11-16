function visValidateChoice(value, allowed, fnName, paramName)
%%% Reject a string-option argument that isn't one of the allowed values
%%% >> visValidateChoice(value, allowed, fnName, paramName)
%%%
%%% @Vis-private helper: nearly every plotting function in visual/ takes one
%%% or more string options ('log'/'linear', 'cartesian'/'polar', 'type'/
%%% 'score'/'confidence', ...) and each had grown its own hand-written
%%% four-line check, differing in message punctuation but not in behaviour.
%%% This is that check, once.
%%%
%%% INPUT: 1. value     .. the caller-supplied option value
%%%        2. allowed   .. cell array of permitted values
%%%        3. fnName    .. the calling function's name, for the error id
%%%        4. paramName .. the parameter's name, for the error id and message
%%%
%%% The error identifier is 'Vis:<fnName>:bad<ParamName>' (e.g.
%%% 'Vis:plotGreenMap:badScaleMode'), the Vis:function:condition form of every
%%% error in visual/.

    if ~ismember(value, allowed)
        quoted = cellfun(@(s) ['''' s ''''], allowed, 'UniformOutput', false);
        if isscalar(quoted)
            allowedStr = quoted{1};
        else
            allowedStr = [strjoin(quoted(1:end-1), ', ') ' or ' quoted{end}];
        end
        error(['Vis:' fnName ':bad' upper(paramName(1)) paramName(2:end)], ...
            '%s must be %s (got ''%s'').', paramName, allowedStr, char(value));
    end
end
