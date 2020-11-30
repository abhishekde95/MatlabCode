function [sampleTimes, TemporalWindow] = trapezoidalTemporalWindowCreate( temporalParams )
% Author : Abhishek De, 08/16
% Creates a trapezoidal Window

sampleTimes = 0:temporalParams.stimulusSamplingIntervalInSeconds:temporalParams.stimulusDurationInSeconds;
L = numel(sampleTimes);
TemporalWindow = ones(L,1);
TemporalWindow(1:ceil(0.25*L)) = linspace(0,1,ceil(0.25*L));
TemporalWindow(ceil(0.75*L):L) = fliplr(linspace(0,1,L-ceil(0.75*L)+1));
end

