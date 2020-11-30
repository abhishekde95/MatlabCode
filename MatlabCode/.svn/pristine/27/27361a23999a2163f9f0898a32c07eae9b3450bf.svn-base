function spikeTimes = GetSpikeTimes(expoDataSet, spikeNum, passIDs, ...
                                    startTimes, endTimes, isDuration, timeUnit)
% function spikeTimes = GetSpikeTimes(expoDataSet, spikeNum, passIDs, ...
%                                   startTimes, endTimes, isDuration, timeUnit)
%
% This function can be used to obtain spikes that occurred during a particular set of passes.
% spikeNum is the numbered spike template (aka spike id) from Expo
% It returns a cell array each cell of which contains a vector of spike times that fall within a window
% defined by one of the passIDs. This is also known as a set of spike rasters.  The returned values of 
% spike times are relative to each value specified in startTimes.
%
% If, for example, the third entry in the supplied passIDs vector is 22, then the vector of spike times
% appearing in the third cell will be those that fall within a time window belonging to passID 22.
%
% All times are in 0.1mS units unless a differnt timeUnit is specified e.g. 'sec' or 'msec'.
%
% The exact placement and size of the window can be modified using the startTimes and endTimes parameters which
% can be either vectors of the same size as the passIDs vector or scalars.
%
% startTimes and endTimes are always relative to the pass time boundaries.
% Setting startTime=0 and endTime=0 ensures the period used is the full width of each pass time window
% i.e. the entire period during which each pass was executed.
% Setting any other scalar value adds a fixed relative displacement to the beginning and end of each window
% If startTimes and endTimes are vectors the size and placement of each window is individually adjusted 
% by the appropriate amounts. Note that the returned spikeTimes are such that a spike occurrring at 
% startTime(s) will have the value 0.
%
% Setting isDuration=1 causes endTimes to be interprested as a list of durations of each window.
% Otherwise with isDuration=0, endTimes is interpreted as a list of end times relative to the end of the pass.
%
%
% See also ReadExpoXML, GetSlots, GetPasses, GetEvents, GetAnalog, 
% GetPSTH, PlotPSTH, GetWaveforms, GetStartTimes, GetEndTimes, GetDuration,
% GetConversionFactor.
%
%   Author:      Julian Brown
%	Version:     1.1
%   Last updated:  2005-03-28
%   E-mail:      julian@monkeybiz.stanford.edu


    spikeTimes = [];
    
    matlabImportVersion = '1.1';
    CheckExpoVersion(expoDataSet, matlabImportVersion);
    
    units = expoDataSet.environment.Conversion.units;
    
    if exist('timeUnit') && length(timeUnit) > 0
        timeUnitNum = GetUnitNumFromName(units, timeUnit, units.U_BASETIME); 
    else
        timeUnitNum = units.U_BASETIME; 
    end

    [numOfPasses passIDs] = TransformToColumnVector(passIDs);

    if ~isnumeric(passIDs)
        error('passIDs should be a vector of numbers. It appears to be non numeric.');
    end

    [numOfStartTimes startTimes] = TransformToColumnVector(startTimes);
    [numOfEndTimes endTimes] = TransformToColumnVector(endTimes);
    
    if timeUnitNum ~= units.U_BASETIME;
        startTimes = ConvertNum(units, startTimes, timeUnitNum, units.U_BASETIME);
        endTimes = ConvertNum(units, endTimes, timeUnitNum, units.U_BASETIME);
    end

    if ~isnumeric(spikeNum)
        error('spikeNum is non numeric.');
    elseif size(spikeNum, 2) ~= 1
        error('spikeNum should be a scalar value.');
    end

    spikeTypeNum = find(expoDataSet.spiketimes.IDs == spikeNum);

    if isempty(spikeTypeNum)
        disp(sprintf('Error: there is no spikeNum that matches the value %d.', spikeNum));
        if isempty(expoDataSet.spiketimes.IDs)
            disp('No spikes were collected.')
        else
            disp('Valid values are: ');
            disp(expoDataSet.spiketimes.IDs);
        end
        error('Aborting.');
    end


    if numOfStartTimes > 1  && numOfPasses ~= numOfStartTimes
        error(sprintf('startTimes must be either a scalar or a vector of same length as passIDs. startTimes contains %d values whereas passIDs contains %d.', size(startTimes), numOfPasses));
    end

    if numOfEndTimes > 1  && numOfPasses ~= numOfEndTimes
        error(sprintf('endTimes must be either a scalar or a vector of same length as passIDs.  endTimes contains %d values whereas passIDs contains %d.', size(endTimes), numOfPasses));
    end

    if max(expoDataSet.passes.IDs) < max(passIDs)
       error(sprintf('there appear to be illegal values in the passIDs vector becuase the maximum value within passIDs %d exceeds the maximum value within expoDataSet %d. ', max(passIDs), max(expoDataSet.passes(:, 1))))
    end
    
    if isDuration~=0 & isDuration ~=1
        error('isDuration should be either 0 or 1.');
    end

    % allocate space for the spikeTimes
    spikeTimes = cell(numOfPasses, 1);


    % loop thru each passID
    for i = 1:numOfPasses

        % get the passNumber (i.e. position within the expoDataSet.passes array) for this passID
        passNum = find (expoDataSet.passes.IDs == passIDs(i));

        if size(passNum, 2) < 1
            warning('no entry for passID %d found in expoDataSet', i);
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

        % now extract from spikeTimes data the ones that fall in the time window
        spikeTimesVector = expoDataSet.spiketimes.Times{spikeTypeNum};
        
        spikeTimesVector = spikeTimesVector(find(spikeTimesVector >= adjustedPassStartTime & spikeTimesVector <= adjustedPassEndTime)) - double(adjustedPassStartTime) ;
    
        spikeTimes{i} = ConvertNum(units, spikeTimesVector, units.U_BASETIME, timeUnitNum);
    end
return

