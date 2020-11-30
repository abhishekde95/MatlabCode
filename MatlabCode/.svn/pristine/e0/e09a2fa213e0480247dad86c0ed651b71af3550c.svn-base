function parse_header(p, codes)
global gl

gl.bkgndrgb = GetValsIfPossible([codes.BKGNDRGBCD codes.BKGNDRGBCD], p, 'double');
spectra = reshape(GetValsIfPossible([codes.MONSPDCD codes.MONSPDCD], p, 'double'), [], 3);
fundamentals = reshape(GetValsIfPossible([codes.FUNDSCD codes.FUNDSCD], p, 'double'), [], 3);
P_device = SplineSpd(linspace(380,780,size(spectra,1))', spectra, ...
    linspace(380,780,size(fundamentals,1))');
gl.M = fundamentals'*P_device;

message_REX('PARSEDHEADER');
