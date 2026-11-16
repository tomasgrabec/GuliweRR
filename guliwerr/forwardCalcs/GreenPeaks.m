classdef GreenPeaks
%%% GreenPeaks — generic peak finder for a |G| curve: the local maxima inside a mask
%%% whose height exceeds a fraction of the curve's global maximum.
%%%
%%% Not part of the surface-wave classification (the existence scores GreenScore,
%%% GreenScoreLayered, GreenSpotScore and RitzScore have their own candidate searches).
%%% Used for isoFreq spectra (kParRange standing in for the velocity axis). See
%%% docs/reference/architecture/forwardCalcs_modeClassification.md, "GreenPeaks.find".
%%%
%%% Value class with a Constant property and a Static method only.
%%%
%%% GuliweRR v6.0

    properties (Constant)
        %%% Default height floor (fraction of the global max of absG). Chosen empirically
        %%% against the computed Cu(001) [100]->[010] cut (tutorial_02): the true Rayleigh
        %%% peak's height crosses below ~15 % of the local max between ~32-34 deg there,
        %%% while still exceeding 20 % through 30 deg where it remains dominant. A lower floor
        %%% keeps Rayleigh remnants long past the point where a much taller pseudoRayleigh
        %%% peak has taken over the spectrum. ElasticityInversion's pseudo-SAW search is
        %%% calibrated against this globalMax-relative floor.
        HEIGHT_FRAC = 0.15
    end

    methods (Static)
        function peaks = find(absG, velAxis, mask, heightFrac)
        %%% Peak-find abs(G) within a velocity-axis mask
        %%% >> peaks = GreenPeaks.find(absG, velAxis, mask)
        %%% >> peaks = GreenPeaks.find(absG, velAxis, mask, heightFrac)
        %%%
        %%% INPUT:  1. absG       .. [nVel x 1] abs(G) on velAxis
        %%%         2. velAxis    .. [nVel x 1] axis (must be increasing -- required by
        %%%                          findpeaks(Y,X) when X is supplied)
        %%%         3. mask       .. logical [nVel x 1], the sub-range to search
        %%%  (opt)  4. heightFrac .. height floor as a fraction of max(absG) over the WHOLE
        %%%                          array, not just the mask (default HEIGHT_FRAC)
        %%% OUTPUT: peaks .. struct array with fields vel, height, prominence, width
        %%%                  (0..N entries)
        %%%
        %%% Uses MATLAB's findpeaks (Signal Processing Toolbox). prominence (height above
        %%% the higher of a peak's two bounding valleys) and width (half-prominence width,
        %%% 'WidthReference'='halfprom', in the units of velAxis) are reported for every
        %%% accepted peak but NOT filtered on: the height floor is the only gate.
            if nargin < 4 || isempty(heightFrac); heightFrac = GreenPeaks.HEIGHT_FRAC; end
            if ~any(mask)
                peaks = struct('vel', {}, 'height', {}, 'prominence', {}, 'width', {});
                return
            end
            globalMax = max(absG);
            %%% "no data point clears MinPeakHeight" is an expected outcome here --
            %%% e.g. no pseudoRayleigh feature above shear at this angle -- so
            %%% suppress that one findpeaks warning and restore the state after
            warnState = warning('off', 'signal:findpeaks:largeMinPeakHeight');
            cleanupObj = onCleanup(@() warning(warnState));
            [pkHeights, pkLocs, pkWidths, pkProms] = findpeaks(absG(mask), velAxis(mask), ...
                'MinPeakHeight', heightFrac * globalMax, 'WidthReference', 'halfprom');
            peaks = struct('vel', num2cell(pkLocs(:)'), 'height', num2cell(pkHeights(:)'), ...
                           'prominence', num2cell(pkProms(:)'), 'width', num2cell(pkWidths(:)'));
        end
    end
end
