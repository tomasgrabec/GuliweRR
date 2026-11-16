function pts = visFlattenSignificantModes(smCell, angles, typeAllowed)
%%% Flatten per-angle significantModes struct arrays into flat point arrays
%%% >> pts = visFlattenSignificantModes(smCell, angles, typeAllowed)
%%%
%%% @Vis-private helper shared by Vis.plotClassificationDiagnostic (a whole
%%% angular sweep) and Vis.plotClassificationCurve (a single angle, passed as
%%% a one-element cell). Both plots read the ALREADY-classified struct array
%%% and never recompute a parallel score -- this helper just reshapes it, and
%%% having one copy of the reshape is what keeps the two plots showing the
%%% same per-point quantities.
%%%
%%% INPUT: 1. smCell      .. {1 x nAngles} cell of significantModes struct
%%%                          arrays (ModeClassification's own output shape)
%%%        2. angles       .. [nAngles x 1] angle per cell entry; pass NaN for
%%%                          a single-angle caller with no angle axis to draw
%%%        3. typeAllowed  .. cell of type names to keep
%%%                           (visClassificationTypesAllowed)
%%% OUTPUT: pts .. struct of flat column arrays, one entry per kept point:
%%%           .angle .vel .modeNo .modeVel .greenPeakVel .score .scoreDb   -- and ROW cells
%%%           .type, .method, .grade of per-point names (row because
%%%           strcmp/ismember over a cell of names returns a row mask, which the
%%%           callers index with).
%%%
%%% Every field is present on every significantModes entry regardless of which
%%% recognize* path produced it (ModeClassification's own property
%%% declaration fixes the field set), so no per-field existence check is
%%% needed here.

    pts = struct('angle', [], 'vel', [], 'modeNo', [], 'modeVel', [], 'type', {{}}, 'method', {{}}, ...
                 'grade', {{}}, 'greenPeakVel', [], 'score', [], 'scoreDb', []);

    for ia = 1:numel(smCell)
        sm = smCell{ia};
        sm = sm(ismember({sm.type}, typeAllowed));
        n  = numel(sm);
        if n == 0; continue; end
        pts.angle        = [pts.angle;        repmat(angles(ia), n, 1)];
        pts.vel          = [pts.vel;          [sm.vel]'];
        pts.modeNo       = [pts.modeNo;       [sm.modeNo]'];
        pts.modeVel      = [pts.modeVel;      ModeClassification.modeVels(sm)'];
        pts.type         = [pts.type,         {sm.type}];
        pts.method       = [pts.method,       {sm.method}];
        pts.grade        = [pts.grade,        {sm.grade}];
        pts.greenPeakVel = [pts.greenPeakVel; [sm.greenPeakVel]'];
        pts.score        = [pts.score;        [sm.score]'];
        pts.scoreDb      = [pts.scoreDb;      [sm.scoreDb]'];
    end
end
