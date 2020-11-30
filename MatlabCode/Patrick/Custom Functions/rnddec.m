% This code is for rounding to a specified decimal
% JPW   10/28/11    Created

function [out] = rnddec(input, dec)

out = (round(input .* 10^dec)) ./ 10^dec;
