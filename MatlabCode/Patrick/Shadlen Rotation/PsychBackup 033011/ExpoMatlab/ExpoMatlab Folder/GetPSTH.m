function [m, sem, t] = GetPSTH(spikeTimes, startTime, binWidth, duration, timeUnit, spikeUnit)
% function [m, sem, t] = GetPSTH(spikeTimes, startTime, binWidth, duration, timeUnit, spikeUnit)
%
% spikeTimes: a cell-array of vectors such as would be returned by the call GetSpikeTimes(..., startTime, duration, 1, ...)
% startTime, binWidth, duration: define a histogram (post-stimulus-time histogram, PSTH) which will be computed
% from the spikes and returned.  
% [optional:]
% timeUnit: the Units of all of the above arguments.  Must be 'sec' or 'msec' [default 'msec'].
%     Only used if spikeUnit is 'impulses/sec' (below).
% spikeUnit: makes the y-axis be scaled in 'impulses/sec' or 'impulses/bin' [default: 'impulses/bin'].
%
% Note: Perhaps a later version will allow startTime and duration to have the more general form that can
% be supplied to GetSpikeTimes.
%
% See also ReadExpoXML, GetSlots, GetPasses, GetEvents, GetSpikeTimes, GetAnalog,
% GetSpikeCounts, PlotPSTH, GetStartTimes, GetEndTimes, GetDuration,
% GetConversionFactor.
%
%   Author:      Jim Muller
%   Last updated:  2005-01-11
%   E-mail:      jim@monkeybiz.stanford.edu



if ~exist('timeUnit'), timeUnit = 'msec'; end % NB: Its default is actually unused/irrelevant.
if ~exist('spikeUnit'), spikeUnit = 'impulses/bin'; end

if strcmp(spikeUnit, 'impulses/sec')
    if strcmp(timeUnit, 'sec')
        scale = 1 ./ binWidth;
    elseif strcmp(timeUnit,'msec')
        scale = 1000 ./ binWidth;
    else
        error(['GetPSTH: ' timeUnit ' is not a valid value for timeUnit']); 
    end
elseif strcmp(spikeUnit, 'impulses/bin')
    scale = 1;
else
     error(['GetPSTH: ' spikeUnit ' is not a valid value for spikeUnit']);
end

        
        
t = startTime:binWidth:(startTime+duration);
ntrials = length(spikeTimes);
% Accumulate raw counts here for every trial for every bin.  This allows s.e.m. to be computed.
pstmatrix = zeros(ntrials, length(t)); 
for idx = 1:ntrials
    if ~isempty(spikeTimes{idx})
        pstmatrix(idx, :) = histc(spikeTimes{idx}, t);
    end
end

m = mean(pstmatrix, 1);
sem = std(pstmatrix, 1) ./ sqrt(ntrials);

% Get rid of the last element, which is just a sentinel for histc.
m = scale * m(1:end-1);
sem = scale * sem(1:end-1);
t = t(1:end-1); 