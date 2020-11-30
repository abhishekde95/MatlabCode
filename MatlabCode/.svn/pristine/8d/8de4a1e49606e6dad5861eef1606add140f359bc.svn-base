function [waveforms firstSampleTimes sampleInterval] = GetWaveforms(expoDataSet, channelID, passIDs, startTimes, ...
                                                                    endTimes, isDuration, timeUnit, waveformUnit, forceCellArray)
%function [waveforms firstSampleTimes sampleInterval] = GetWaveforms(expoDataSet, channelID, passIDs, startTimes, endTimes, isDuration, timeUnit, waveformUnit, forceCellArray)
%
% This Expo utility retrieves waveform data recorded during a particular set of passes.
% It returns either a matrix or cell array in which each row or cell contains waveform (spike) data belonging to a time window
% defined by one of the passIDs.  The choice between a matrix and cell array depends on whether the pass windows are all of  
% equal length.  If they are the function returns a matrix unless the optional parameter forceCellArray is set to 1. 
%
% The function also returns firstSampleTimes, a vector of times relative to startTimes for the first requested 
% sample in each pass, and the scalar value sampleinterval.
% All times are in 0.1mS units unless a differnt timeUnit is specified e.g. 'Seconds'. The waveform output is in 
% normalized volts (+1 to -1) unless a different waveformUnit is specified e.g. 'Volts'.
%
% The exact placement and size of the window can be modified using the startTimes and endTimes parameters which
% can be either vectors of the same size as the passIDs vector or scalars.
%
% startTimes and endTimes are relative to the pass time boundaries.
% Setting startTime=0 and endTime=0 ensures the period used is the full width of each pass time window
% i.e. the entire period during which each pass was executed.
% Setting any other scalar value adds a fixed relative displacement to the beginning and end of each window
% If startTimes and endTimes are vectors the size and placement of each window is individually adjusted by the appropriate amounts.
%
% Setting isDuration=1 causes endTimes to be interpreted as a list of durations for each window.
% Otherwise with isDuration=0, endTimes is interpreted as a list of times relative to the end of the window.
%
% channelID refers to the particular channel for the waveform data
%
% See also ReadExpoXML, GetSlots, GetPasses, GetEvents, GetSpikeTimes, GetAnalog, 
% GetPSTH, PlotPSTH, GetStartTimes, GetEndTimes, GetDuration,
% GetConversionFactor.
%
%   Author:      Julian Brown
%	Version:     1.1
%   Last updated:  2005-03-30
%   E-mail:      julianb@stanford.edu


    waveforms = [];
    firstSampleTimes = [];
    
    matlabImportVersion = '1.1';
    CheckExpoVersion(expoDataSet, matlabImportVersion);
    
    if size(expoDataSet.waveforms, 1) == 0
        error('There is no waveform data.')
    end
    
    numOfPackets = expoDataSet.waveforms.NumOfBuffers;
    numOfChannels = expoDataSet.waveforms.NumOfChannels;
    numOfSamplesPerPacket = expoDataSet.waveforms.FramesPerBuf;
    packetDuration = 10000 * double(expoDataSet.waveforms.FramesPerBuf) / double(expoDataSet.waveforms.SampleRate);        
    sampleInterval = double(packetDuration)/double(numOfSamplesPerPacket);

    [numOfPasses passIDs] = TransformToColumnVector(passIDs);
    
    if ~exist('forceCellArray') , forceCellArray = 0; end
    
    units = expoDataSet.environment.Conversion.units;
   
    if exist('timeUnit') && length(timeUnit) > 0
        timeUnitNum = GetUnitNumFromName(units, timeUnit, units.U_BASETIME); 
    else
        timeUnitNum = units.U_BASETIME; 
    end
    
    if exist('waveformUnit') && length(waveformUnit) > 0 
        waveformUnitNum = GetUnitNumFromName(units, waveformUnit, units.U_NORMVOLT); 
    else
        waveformUnitNum = units.U_NORMVOLT; 
    end

    if ~isnumeric(passIDs)
        error('passIDs should be a vector of numbers. It appears to be non numeric.');
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
    
    if isDuration~=0 && isDuration ~=1
        error('isDuration should be either 0 or 1.');
    elseif isDuration == 1 && max(endTimes) <= 0
        error('isDuration is set to 1 but the durations are all negative or 0.');
    end
    
    [numOfStartTimes startTimes] = TransformToColumnVector(startTimes);
    [numOfEndTimes endTimes] = TransformToColumnVector(endTimes);
    
    if timeUnitNum ~= units.U_BASETIME;
        startTimes = ConvertNum(units, startTimes, timeUnitNum, units.U_BASETIME);
        endTimes = ConvertNum(units, endTimes, timeUnitNum, units.U_BASETIME);
    end
    
    compressionScale = expoDataSet.waveforms.Channels.Scales(channelID+1);
    compressionOffset = expoDataSet.waveforms.Channels.Offsets(channelID+1);


    if numOfStartTimes > 1  && numOfPasses ~= numOfStartTimes
        error(sprintf('startTimes must be either a scalar or a vector of same length as passIDs.  startTimes contains %d values whereas passIDs contains %d.', size(startTimes), numOfPasses));
    end

    if numOfEndTimes > 1  && numOfPasses ~= numOfEndTimes
        error(sprintf('endTimes must be either a scalar or a vector of same length as passIDs.  endTimes contains %d values whereas passIDs contains %d.', size(endTimes), numOfPasses));
    end

    if max(expoDataSet.passes.IDs) < max(passIDs)
       error(sprintf('There appear to be out of range values in the passIDs vector becuase the maximum value within passIDs %d exceeds the maximum value within expoDataSet %d. ', max(passIDs), max(expoDataSet.passes(:, 1))))
    end
        
    % allocate space for the waveforms
    waveforms = cell(numOfPasses, 1);
    
    conversionFactorForNormalizedVolts = 1;
    conversionFactor = ConvertNum(units, conversionFactorForNormalizedVolts, units.U_NORMVOLT, waveformUnitNum);
    
    isFirstPass = 1;
    windowDurationsVary = 0;
    
    % loop thru each passID
    for i = 1:numOfPasses

        % get the passNumber (i.e. position within the expoDataSet.passes array) for this passID
        passNum = find (expoDataSet.passes.IDs == passIDs(i));

        if size(passNum, 2) < 1
            warning('No entry for passID %d found in expoDataSet', i);
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
        
        % get window duration for this pass
        windowDuration = adjustedPassEndTime - adjustedPassStartTime;
        
        % check whether durations vary - will be used to decide whether we
        % can turn results into a matrix at the end
        if isFirstPass == 1
            previousWindowDuration = windowDuration;
            isFirstPass = 0;
        elseif windowDurationsVary == 0 && windowDuration ~= previousWindowDuration 
            windowDurationsVary = 1;
        end
        
        % check whether the window is zero in size
        if windowDuration <= 0
            waveforms{i} = [];
            firstSampleTimes(i) = [NaN];
            continue
        end
        
        totalNumOfSamples = ceil(double(windowDuration)/sampleInterval);
        samplesRemaining = totalNumOfSamples;
        
        % find the first wave packet that spans the onset start
        packetNum = find (expoDataSet.waveforms.Times <= adjustedPassStartTime & expoDataSet.waveforms.Times > adjustedPassStartTime - packetDuration);
        
        isFirstPacket = 1;
        
        if packetNum == 0
            %could not find a packet that spans onset so look for next best
            %thing - one that starts after the onset but before the window end
            packetNum = min(find (expoDataSet.waveforms.Times >= adjustedPassStartTime & expoDataSet.waveforms.Times < adjustedPassEndTime));
            
            % calculate how big the gap is
            startGap = expoDataSet.waveforms.Times(packetNum) - adjustedPassStartTime;
            
            numOfMissingValues = fix(double(startGap)/sampleInterval);
            samplesRemaining = samplesRemaining - numOfMissingValues;
                    
            % generate an array of NaNs using zeros()/0                 
            warning('off','MATLAB:divideByZero')
            waveforms{i} = zeros(1, numOfMissingValues)/0;
            warning('on','MATLAB:divideByZero')
            
            isFirstPacket = 0;
        end
        
        if packetNum > 0
            % start loop that will gather all of the packets required to reconstruct the full waveform data for this pass    
            while 1
                % get the start and finish times of this packet
                startTimeOfPacket = expoDataSet.waveforms.Times(packetNum);
                endTimeOfPacket = startTimeOfPacket + packetDuration;

                % obtain time difference between the start of waveforms and the window onset
                timeDiffForStart = adjustedPassStartTime - startTimeOfPacket;

                % calculate boundary conditions within packet
                if timeDiffForStart <= 0
                    % the pass starts before this packet so take data from the beginning of the packet
                    firstSampleNum = 0;
                    samplesRemainingInPacket = numOfSamplesPerPacket;
                else
                    firstSampleNum = ceil(double(timeDiffForStart) / double(sampleInterval));
                    samplesRemainingInPacket = numOfSamplesPerPacket - firstSampleNum;
                end                        

                startSamplePos =  firstSampleNum * numOfChannels + channelID + 1;

                if samplesRemaining >= samplesRemainingInPacket
                    lastSampleNum = numOfSamplesPerPacket;
                    samplesRemaining = samplesRemaining - samplesRemainingInPacket;
                else
                    lastSampleNum = samplesRemaining + firstSampleNum;
                    samplesRemaining = 0;
                end

                endSamplePos = (lastSampleNum - 1) * numOfChannels + channelID + 1;

                if isFirstPacket 
                    firstSampleTimes(i) = firstSampleNum * sampleInterval + startTimeOfPacket - adjustedPassStartTime;
                end
                
                % get data from this packet
                packetData = expoDataSet.waveforms.Data{packetNum};

                % extract the data we want 
                % the data is stored as a byte array with different channel data interleaved                       
                data = double(packetData(startSamplePos:numOfChannels:endSamplePos));

                % process the values using the appropriate compression and offset parameters
                if isFirstPacket == 1
                    waveforms{i} = conversionFactor*(compressionOffset + data/compressionScale);
                    isFirstPacket = 0;
                else
                    waveforms{i} = cat(2, waveforms{i}, conversionFactor*(compressionOffset + data/compressionScale));
                end

                if samplesRemaining == 0 %adjustedPassEndTime <= endTimeOfPacket
                    % we don't need any more packets for this pass
                    break
                end
                
                % check next packet
                packetNum = packetNum + 1;
                if packetNum > numOfPackets 
                    % we ran out of packets - so tack on a NaN for goodwill and finish
                    waveforms{i} = cat(2, waveforms{i}, [NaN]);
                    windowDurationsVary = 1;
                    break;
                end
                    
                startTimeOfNextPacket = expoDataSet.waveforms.Times(packetNum);
                
                if startTimeOfNextPacket-endTimeOfPacket > sampleInterval
                    % uh oh - like the fossil record there are gaps in the waveform records
                    % do we have any more data for this pass?
                    if startTimeOfNextPacket >= adjustedPassEndTime
                        % we ran out of data - tack on a NaN and finish
                        waveforms{i} = cat(2, waveforms{i}, [NaN]);
                        windowDurationsVary = 1;
                        break;
                    end
                    
                    %we have more data but we'll need to fill the gap with NaNs
                    numOfMissingValues = fix(double(startTimeOfNextPacket - endTimeOfPacket)/sampleInterval);
                    
                    % generate an array of NaNs using zeros()/0                 
                    warning('off','MATLAB:divideByZero')
                    waveforms{i} = cat(2, waveforms{i}, zeros(1, numOfMissingValues)/0);
                    warning('on','MATLAB:divideByZero')
                end
            end
        else
            waveforms{i} = [NaN];
            firstSampleTimes(i) = [NaN];
            windowDurationsVary = 1;
        end        
    end
    
    if forceCellArray ~= 1 && isFirstPass == 0
       % check whether we can convert to a matrix
       if windowDurationsVary == 0
           waveformsMatrix = zeros(totalNumOfSamples, numOfPasses, 'double');
           for i=1:numOfPasses
               waveformsMatrix(:, i) = waveforms{i};
           end
           waveforms = waveformsMatrix;
       end
    end
    
    if timeUnitNum ~= units.U_BASETIME;
        firstSampleTimes = ConvertNum(units, firstSampleTimes, units.U_BASETIME, timeUnitNum);
        sampleInterval = ConvertNum(units, sampleInterval, units.U_BASETIME, timeUnitNum);
    end
   
return

