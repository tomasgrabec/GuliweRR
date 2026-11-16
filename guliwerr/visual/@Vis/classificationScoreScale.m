function s = classificationScoreScale(method, isLayered)
%%% The score scale of a classification method: its standards and the display range
%%% >> s = Vis.classificationScoreScale(method, isLayered)
%%%
%%% INPUT: 1. method    .. 'green' | 'refined' | 'ritz' (layered: 'green' | 'ritz')
%%%        2. isLayered .. true for a layered stack
%%% OUTPUT: s .. struct:
%%%           high, low  .. the 'high' / 'low' threshold values (the classes' own
%%%                         constants, read here, never copied)
%%%           lowScale   .. 'score' if 'low' is compared with the entry score, like
%%%                         'high' and a numeric standard; 'structure-relative' if it
%%%                         is compared with a separate structure-relative score
%%%                         ('green' / 'refined' half-space: lowScore / scoreLow)
%%%           sizeRange  .. [lo hi] score range the marker size spans; scores
%%%                         outside it are clipped
%%%           name       .. display name ('Green' | 'Refined' | 'Ritz')
%%%
%%% Scores are log10 amplitude ratios (dB = 20 x score), comparable only within one
%%% method. sizeRange starts at the looser ('high') threshold, where the weakest
%%% accepted entry sits, and ends one decade (20 dB) above the stricter 'low'
%%% threshold where both are on the score scale, else 1.5 decades (30 dB) above
%%% 'high'. So the strong entries of a typical map spread over the size range
%%% instead of saturating it.
%%%
%%% Shared by Vis.plotClassificationDiagnostic / plotClassificationCurve (marker
%%% size) and the viewers' classification bar (standard values, dB labels).
    if isstring(method); method = char(method); end
    switch method
        case 'green'
            if isLayered
                s = struct('high', GreenScoreLayered.THRESHOLD_HIGH, ...
                    'low', GreenScoreLayered.THRESHOLD_LOW, 'lowScale', 'score');
            else
                s = struct('high', GreenScore.THRESHOLD_HIGH, ...
                    'low', GreenScore.THRESHOLD_LOW, 'lowScale', 'structure-relative');
            end
        case 'refined'
            if isLayered
                error('Vis:classificationScoreScale:noRefinedLayered', ...
                    "Method 'refined' does not exist for layered stacks.");
            end
            s = struct('high', GreenSpotScore.THRESHOLD_HIGH, ...
                'low', GreenSpotScore.THRESHOLD_LOW, 'lowScale', 'structure-relative');
        case 'ritz'
            if isLayered
                s = struct('high', RitzScore.LAY_THRESHOLD_HIGH, ...
                    'low', RitzScore.LAY_THRESHOLD_LOW, 'lowScale', 'score');
            else
                s = struct('high', RitzScore.HS_THRESHOLD_HIGH, ...
                    'low', RitzScore.HS_THRESHOLD_LOW, 'lowScale', 'score');
            end
        otherwise
            error('Vis:classificationScoreScale:badMethod', ...
                "method must be 'green', 'refined' or 'ritz'; got '%s'.", method);
    end
    [ids, names] = CalcCommon.classificationMethods(false);
    s.name = names{strcmp(ids, method)};
    if strcmp(s.lowScale, 'score')
        s.sizeRange = [s.high, max(s.high, s.low) + 1.0];
    else
        s.sizeRange = [s.high, s.high + 1.5];
    end
end
