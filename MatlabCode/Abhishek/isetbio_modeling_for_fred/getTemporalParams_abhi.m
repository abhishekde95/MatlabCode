function temporalParams = getTemporalParams_abhi()
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
temporalParams.frameRate = 60;
temporalParams.windowTauInSeconds = 0.1;
temporalParams.stimulusDurationInSeconds = 0.5; % 500 ms stimulus presentation
temporalParams.stimulusSamplingIntervalInSeconds = 1/temporalParams.frameRate;
temporalParams.eyesDoNotMove = false; 
temporalParams.millisecondsToInclude = 500;
temporalParams.millisecondsToIncludeOffset = 35;

end

