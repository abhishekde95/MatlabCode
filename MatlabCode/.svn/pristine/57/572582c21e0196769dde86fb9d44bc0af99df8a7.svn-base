function [counts, m, sem] = GetSpikeCounts(spikeTimes, startTimes, durations, timeUnit, spikeUnit)
% function [counts, m, sem] = GetSpikeCounts(spikeTimes, startTimes, durations, timeUnit, spikeUnit)
%
% Compute spike counts (or spike rates) for each imtrial in a raster.  Also return mean m+-sem of those counts.
%
% spikeTimes: a cell-array of vectors such as would be returned by the call GetSpikeTimes(..., startTime, duration, 1, ...)
%   This is also known as a set of spike rasters.
% startTimes: start of interval to count over within each raster
% durations: duration of interval to count over within each raster
% [optional:]
% timeUnit: the Units of all three of the above arguments.  Must be 'sec' or 'msec' [default 'msec'].
% spikeUnit: makes the y-axis be scaled in 'impulses/sec' or 'impulses' [default: 'impulses'].
%
% See also ReadExpoXML, GetSlots, GetPasses, GetEvents, GetSpikeTimes, GetAnalog,
% GetPSTH, PlotPSTH, GetStartTimes, GetEndTimes, GetDuration,
% GetConversionFactor.
%
%   Author:      Jim Muller
%   Last updated:  2005-01-11
%   E-mail:      jim@monkeybiz.stanford.edu


if ~exist('timeUnit'), timeUnit = 'msec'; end % NB: Its default is actually unused/irrelevant.
if ~exist('spikeUnit'), spikeUnit = 'impulses'; end

ntrials = length(spikeTimes);

% Regularize start/dur vectors
[numOfStartTimes startTimes] = TransformToColumnVector(startTimes);
[numOfDurations durations] = TransformToColumnVector(durations);
if numOfStartTimes == 1, startTimes = startTimes * ones(ntrials, 1); end % Vectorize scalars.
if numOfDurations == 1, durations = durations * ones(ntrials, 1); end

if numOfStartTimes > 1  && ntrials ~= numOfStartTimes
    error(sprintf('startTimes must be either a scalar or a vector of same length as spikeTimes. startTimes contains %d values whereas spikeTimes contains %d.', size(startTimes), ntrials));
end

if numOfDurations > 1  && ntrials ~= numOfDurations
    error(sprintf('durations must be either a scalar or a vector of same length as spikeTimes.  durations contains %d values whereas spikeTimes contains %d.', size(durations), ntrials));
end

% Scale as per units.  NB scale is a row-vector, so we can return as a row, as is conventional.        
if strcmp(spikeUnit, 'impulses/sec')
    if strcmp(timeUnit, 'sec')
        scale = 1 ./ durations';
    elseif strcmp(timeUnit,'msec')
        scale = 1000 ./ durations';
    else
        error(['GetPSTH: ' timeUnit ' is not a valid value for timeUnit']); 
    end
elseif strcmp(spikeUnit, 'impulses')
    scale = 1;
else
     error(['GetPSTH: ' spikeUnit ' is not a valid value for spikeUnit']);
end

            
% Accumulate raw counts here for every trial for every bin.  This allows s.e.m. to be computed.
counts = zeros(1, ntrials); 
for idx = 1:ntrials
    counts(idx) = length(find( (startTimes(idx) <= spikeTimes{idx}) & (spikeTimes{idx} < startTimes(idx) + durations(idx)) ));
end

counts = scale .* counts;
if ntrials > 0
    m = mean(counts);
    sem = std(counts) ./ sqrt(ntrials);
else
    m = NaN;
    sem = NaN;
end
