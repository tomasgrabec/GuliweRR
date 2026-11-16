function omega0 = freqToOmega0(freqHz, lambdaRefNm)
%%% Convert a real frequency to GreenSolver's normalised isoFreq omega0
%%% >> omega0 = Calc.freqToOmega0(freqHz, lambdaRefNm)
%%%
%%% INPUT: 1. freqHz      .. real frequency [Hz] (scalar or array)
%%%        2. lambdaRefNm .. reference wavelength [nm] -- MUST be the same
%%%                          wavelength passed to Structure.ApplyWavelength when
%%%                          building relDepth for the isoFreq structure. That
%%%                          call is what fixes what "|k∥|=1" means in real
%%%                          units; use a different reference here and omega0
%%%                          and relDepth disagree about the length scale
%%% OUTPUT: omega0 .. GreenSolver.calc's isoFreq reference velocity [km/s]
%%%
%%% omega0 is the phase velocity that would give |k∥|=1 at the chosen
%%% lambdaRefNm — the standard v = f*lambda relation at that one reference
%%% wavelength, converted from [Hz]x[nm] to [km/s] by 1e-12.
%%%
%%% lambdaRefNm is a free normalisation choice: any positive value gives the
%%% same physical answer, affecting only the internal relDepth/kParRange
%%% numbers. Pick something of the order of the expected SAW wavelength at the
%%% frequency of interest, to keep relDepth*hScale in a well-conditioned range
%%% (see GreenSolver.calc's CONDITIONING note).

    validateattributes(freqHz, {'double'}, {'real', 'finite', 'positive'}, ...
        'freqToOmega0', 'freqHz');
    validateattributes(lambdaRefNm, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
        'freqToOmega0', 'lambdaRefNm');

    omega0 = freqHz .* lambdaRefNm * 1e-12;   % [Hz]*[nm] -> [km/s]
end
