function [analogData firstSampleTimes sampleInterval] ...
    = GetAnalog(expoDataSet, channelID, passIDs, startTimes, endTimes, isDuration, timeUnit, analogUnit, forceCellArray)
%function [analogData firstSampleTimes sampleInterval] ...
%   = GetAnalog(expoDataSet, channelID, passIDs, startTimes, endTimes, isDuration, timeUnit, analogUnit, forceCellArray)
%
% This Expo utility retrieves analog data recorded during a particular set of passes.
% It returns either a matrix or cell array in which each row or cell contains analog data belonging to a time window
% defined by one of the passIDs.  The choice between a matrix and cell array depends on whether the pass windows are all
% of equal length.  If they are the function returns a matrix unless the optional parameter forceCellArray is set to 1.
%
% The function also returns firstSampleTimes, a vector of times relative to startTimes for the first requested
% sample in each pass, and the scalar value sampleinterval.
% All times are in 0.1mS units unless a differnt timeUnit is specified e.g. 'sec' or 'msec'.
% The analog output is in normalized volts (+1 to -1) unless a different analogUnit is specified
% e.g. 'deg'.
%
% The exact placement and size of the window can be modified using the startTimes and endTimes parameters which
% can be either vectors of the same size as the passIDs vector or scalars.
%
% startTimes and endTimes are always relative to the pass time boundaries.
% Setting startTime=0 and endTime=0 ensures the period used is the full width of each pass time window
% i.e. the entire period during which each pass was executed.
% Setting any other scalar value adds a fixed relative displacement to the beginning and end of each window
% If startTimes and endTimes are vectors the size and placement of each
% window is individually adjusted by the appropriate amounts.
%
% Setting isDuration=1 causes endTimes to be interpreted as a list of durations of each window.
% Otherwise with isDuration=0, endTimes is interpreted as a list of times relative to the endo fhte pass.
%
% channelID refers to the particular channel for the analog data
%
% See also ReadExpoXML, GetSlots, GetPasses, GetEvents, GetSpikeTimes,
% GetPSTH, PlotPSTH, GetWaveforms, GetStartTimes, GetEndTimes, GetDuration,
% GetConversionFactor.
%
%   Author:      Julian Brown
%	Version:     1.1
%   Last updated:  2005-03-28
%   E-mail:      julianb@stanford.edu

    analogData = [];
    firstSampleTimes = [];

    matlabImportVersion = '1.1';
    CheckExpoVersion(expoDataSet, matlabImportVersion);

    sampleInterval = expoDataSet.analog.SampleInterval;
    numOfChannels = expoDataSet.analog.NumOfChannels;

    if ~exist('forceCellArray') , forceCellArray = 0; end

    units = expoDataSet.environment.Conversion.units;

    if exist('timeUnit') && length(timeUnit) > 0
        timeUnitNum = GetUnitNumFromName(units, timeUnit, units.U_BASETIME);
    else
        timeUnitNum = units.U_BASETIME;
    end

    if exist('analogUnit') && length(analogUnit) > 0
        analogUnitNum = GetUnitNumFromName(units, analogUnit, units.U_NORMVOLT);
    else
        analogUnitNum = units.U_NORMVOLT;
    end

    [numOfPasses passIDs] = TransformToColumnVector(passIDs);

    if ~isnumeric(passIDs)
        error('passIDs should be a vector of numbers.  It appears to be non numeric.');
    end

    [numOfStartTimes startTimes] = TransformToColumnVector(startTimes);
    [numOfEndTimes endTimes] = TransformToColumnVector(endTimes);

    if timeUnitNum ~= units.U_BASETIME;
        startTimes = ConvertNum(units, startTimes, timeUnitNum, units.U_BASETIME);
        endTimes = ConvertNum(units, endTimes, timeUnitNum, units.U_BASETIME);
    end

    if ~isnumeric(channelID)
        error('channelID is non numeric.');
    elseif size(channelID, 2) ~= 1
        error('channelID should be a scalar value.');
    end

    if channelID + 1 > numOfChannels
        error(sprintf('channelID %d is out of range. There are only %d channels starting at channel 0.', channelID, numOfChannels));
    elseif channelID<0
        error(sprintf('%d is an invalid channelID value. It must be 0 or more.', channelID));
    end


    if numOfStartTimes > 1  && numOfPasses ~= numOfStartTimes
        error(sprintf('startTimes must be either a scalar or a vector of same length as passIDs.  startTimes contains %d values whereas passIDs contains %d.', size(startTimes), numOfPasses));
    end

    if numOfEndTimes > 1  && numOfPasses ~= numOfEndTimes
        error(sprintf('endTimes must be either a scalar or a vector of same length as passIDs.  endTimes contains %d values whereas passIDs contains %d.', size(endTimes), numOfPasses));
    end

    if max(expoDataSet.passes.IDs) < max(passIDs)
       error(sprintf('there appear to be illegal values in the passIDs vector becuase the maximum value within passIDs %d exceeds the maximum value within expoDataSet %d. ', max(passIDs), max(expoDataSet.passes(:, 1))))
    end

    if isDuration~=0 && isDuration ~=1
        error('isDuration should be either 0 or 1.');
    elseif isDuration == 1 && max(endTimes) <= 0
        error('isDuration is set to 1 but the durations are all negative or 0.');
    end

    % allocate space for the analogData
    analogData = cell(numOfPasses, 1);

    conversionFactorForNormalizedVolts = 1/32767;
    conversionFactor = ConvertNum(units, conversionFactorForNormalizedVolts, units.U_NORMVOLT, analogUnitNum);

    isFirstPass = 1;
    windowDurationsVary = 0;

    % loop thru each passID
    for i = 1:numOfPasses

        % get the passNumber (i.e. position within the expoDataSet.passes array) for this passID
        passNum = find (expoDataSet.passes.IDs == passIDs(i));

        if size(passNum, 2) < 1
            warning('Warning: no entry for passID %d found in expoDataSet', i);
            continue
        end


        passStartTime = expoDataSet.passes.StartTimes(passNum);
        passEndTime = expoDataSet.passes.EndTimes(passNum);

        % determine the time window for which we want to extract spikes
        if numOfStartTimes == 1
            adjustedPassStartTime = passStartTime + startTimes;
        else
            adjustedPassStartTime = passStartTime + startTimes(i);
        end

        if numOfEndTimes == 1
            if isDuration == 1
                adjustedPassEndTime = adjustedPassStartTime + endTimes;
            else
                adjustedPassEndTime = passEndTime + endTimes;
            end
        elseif isDuration == 1
            adjustedPassEndTime = adjustedPassStartTime + endTimes(i);
        else
            adjustedPassEndTime = passEndTime + endTimes(i);
        end

        % get window duration
        windowDuration = adjustedPassEndTime - adjustedPassStartTime;

        % check whether durations vary - will be used to decide whether we
        % can turn results into a matrix at the end
        if isFirstPass == 1
            previousWindowDuration = windowDuration;
            isFirstPass = 0;
        elseif windowDurationsVary == 0 && windowDuration ~= previousWindowDuration
            windowDurationsVary = 1;
        end

        if windowDuration<=0, continue, end

        numOfSamplesInWindow = int32(double(windowDuration)/double(sampleInterval));

        % find the appropriate segment
        segmentNum = max(find (expoDataSet.analog.Segments.StartTimes <= adjustedPassStartTime));

        if segmentNum > 0
            % get the start and finish times of this segment
            startTimeOfSegment = expoDataSet.analog.Segments.StartTimes(segmentNum);

            % obtain time difference between the start of analogData and the window onset
            timeDiffForStart = adjustedPassStartTime - startTimeOfSegment;

            numOfSamplesInSegment = expoDataSet.analog.Segments.SampleCounts(segmentNum);

            % calculate boundary conditions
            firstSampleNum = ceil(double(timeDiffForStart) / double(sampleInterval));
            startSamplePos =  firstSampleNum * 2 * numOfChannels + 1;
            endSamplePos = startSamplePos + fix(numOfSamplesInWindow) * 2 * numOfChannels -1;

            firstSampleTimes(i) = firstSampleNum * sampleInterval + startTimeOfSegment - adjustedPassStartTime;


            % check whether the window went beyond the end of the segment
            if startSamplePos > numOfSamplesInSegment * 2 %each sample is two bytes
                analogData{i} = [NaN];
                windowDurationsVary = 1;
                continue
            elseif endSamplePos > numOfSamplesInSegment * 2
                % pick last sample in segment
                lastSamplePos = 2 * (numOfSamplesInSegment - numOfChannels);
            else
                lastSamplePos = endSamplePos;
            end

            % get data from this segment
            segmentData = expoDataSet.analog.Segments.Data{segmentNum};

            % isolate the window we want
            data = double(segmentData(startSamplePos:lastSamplePos));

            % the data is stored as a char array held as an array of shorts
            % i.e. int16s with different channel data interleaved
            % so to extract what we want, first turn into a MxN matrix where M = 2*num of channels
            data2 = reshape(data, 2 * numOfChannels, size(data, 2)/(2 * numOfChannels));

            % now select the channel we want
            rowNum = 2 * channelID +1;
            data3 = data2(rowNum:rowNum + 1, :)';

            % combine values in each row to form a signed int16
            analogData{i} = double(256*((data3(:, 1) > 127)*(-256) + data3(:, 1)) + data3(:,2)) * conversionFactor;


            if endSamplePos > numOfSamplesInSegment * 2
                numOfMissingValues = (endSamplePos/2 - numOfSamplesInSegment)/numOfChannels + 1;
                % tack on appropriate number of NaNs for time beyond the segment
                % generate an array of NaNs using 0/0
                % but only do this for reasonable numbers i.e. not
                % exceeding 1000 otherwise just add one NaN
                if numOfMissingValues<1000
                    warning('off','MATLAB:divideByZero')
                    analogData{i} = cat(1, analogData{i}, zeros(numOfMissingValues, 1)/0);
                    warning('on','MATLAB:divideByZero')
                else
                    analogData{i} = cat(1, analogData{i}, [NaN]);
                    windowDurationsVary = 1;
                end
            end
        else
            analogData{i} = [NaN];
            firstSampleTimes(i) = [NaN];
            windowDurationsVary = 1;
        end


    end

    if forceCellArray ~= 1 && isFirstPass == 0
       % check whether we can convert to a matrix
       if windowDurationsVary == 0
           analogDataMatrix = zeros(numOfPasses, numOfSamplesInWindow, 'double');
           for i=1:numOfPasses
               analogDataMatrix(i, :) = analogData{i}';
           end
           analogData = analogDataMatrix;
       end
    end

    if timeUnitNum ~= units.U_BASETIME;
        firstSampleTimes = ConvertNum(units, firstSampleTimes, units.U_BASETIME, timeUnitNum);
        sampleInterval = ConvertNum(units, sampleInterval, units.U_BASETIME, timeUnitNum);
    end

return

