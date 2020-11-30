function [in_gamut,scalars] = gamut_extent(ccs)
global gl
[in_gamut,scalars] = gamutCheck([ccs zeros(size(ccs,1), 3-size(ccs,2))]', gl.bkgndrgb, gl.M, 'both');
