function v = guliwerrVersion()
%%% GuliweRR version string -- the single place the version number is set
%%% >> v = guliwerrVersion()
%%%
%%% OUTPUT: v .. 'MAJOR.MINOR.PATCH', char
%%%
%%% MAJOR marks a large restructuring, MINOR new or changed functionality (it
%%% may change or remove existing functions/options -- listed under "Breaking"
%%% in docs/ChangelogPublic.md), PATCH fixes only. Each part is an integer, so
%%% 6.9.0 is followed by 6.10.0. A '-dev' suffix marks code between releases
%%% (e.g. '6.1.0-dev' is work towards 6.1.0, not 6.1.0 itself).
%%%
%%% Recorded in every result saved by Calc.saveResult (meta.guliwerrVersion).
    v = '6.0.0';
end
