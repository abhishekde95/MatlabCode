% This function (hopefully) eliminates roundoff error
% JPW   12/23/11    Created

function [out] = rndofferr(input, dec)


out = (round(input.*10^dec))./10^dec;
%out = (fix(input.*10^dec))./10^dec;