function [numOfRows vector] = TransformToColumnVector(vector)
%function [numOfRows vector] = TransformToColumnVector(vector)
%
% A helper utlity used as part of the Expo Matlab routines
%
% Used in GetEvents, GetAnalog, GetSpikeTimes, GetWaveforms etc.
%
%   Author:      Julian Brown
%   Last updated:  2004-12-28
%   E-mail:      julianb@stanford.edu


if size(vector, 1) == 1 && size(vector, 2) > 1
        % turn row vector into column vector
        vector = vector';
    end

    numOfRows = size(vector, 1);
return
