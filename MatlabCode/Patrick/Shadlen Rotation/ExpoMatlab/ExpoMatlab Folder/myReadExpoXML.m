function expoDataSet = myReadExpoXML(filename, doSpikeWaveforms)
% function ReadExpoXML(filename, doSpikeWaveforms)
%
% Expo utilty to import exported data from Expo into Matlab.
% This function requires Matlab v7.0 (aka Release 14) or greater
%
% Filename: Reads data from ExpoX data xml file into matrices and cell arrays and saves as mat file
% doSpikeWaveforms: 1 if you want to process waveform data, else 0. [default 0].
%
% One then loads the mat file to get the following variables:
%
% (Program:)
% blocks - a structure containing information about each block (state) in the program
%   members are
%   (1) vector of BlockIDs
%   (2) vector of Names
%   (3) cell array of Routines, one cell for each block
%
% routines - a structure within blocks containing info about the routines included in each cell of all the blocks
%   members are:
%   (1) vector of routine IDs indexed by blockID
%   (2) vector of routine Names indexed by blockID
%   (3) vector of routine Labels indexed by blockID
%   (4) cell array of parameters
%   (5) routineMap
%   (6) routineEntries
%
% routineMap: a structure within routines useful as a lookup table to see which routines are used in which blocks
%    (1) RoutineID
%    (2) Vector of RountineEntries IDs (see details for Routine Entries table)
%    (3) Routine Name
%    -- Col 4: Routine Label
%
% routineEntries - a structure enumerating every routine entry in all of the blocks
%    (1) RoutineID
%    (2) Position of routine in the block
%    (3) BlockID
%    (4) Position of block in list of sequence of blocks
%
% parameters - each cell within the routines param array contains
%   (1) vector of param Names
%   (2) vector of param VarNames
%   (3) vector of param Ops
%   (4) vector of param Units
%   (5) vector of param Evalues
%
%
% slots: each row describes a slot in the program's flow-control diagram: ID, Name, ID of Block it executes.
%
% (Data:)
% passes: a structure containing information about each pass in the execution of the program.
%         The strucure contains arrays of pass ids, slot ids, block ids, starttimes, endtimes
%         (times are always in .1 msec units) and events.
%         The block ID is included in the pass as well as in the slot because if the pass
%         belongs to Expo's Matrix object the block id is not derivable from the slot.
%
% The events object within passes (expoDataSet.passes.events) contains a list of event structures, one for each pass
% Each event structure corresponds to the run routines executed during a pass and contains
%       (1) Array of RoutineIDs
%       (2) Times stamps at which they executed or finished
%       (3) Cellular array of data for the routine parameters.  Each cell corresponds to the data for each routine
%
% analogSampleInterval: a time, e.g. 29 means 2.9 msec, i.e. 1/344 sec.
% analogNumofChannels: total number of analog channels (e.g. eyetraces etc) collected.
% analogSegments: for each segment there are:
%    (1) Number of Samples (2) Start Time (3) character array containing the samples.
%    Within (3) the channels are interleaved.
%    Since there are AnalogNumofChannels, the samples for channel 1 are chars at positions
%       1, 1+AnalogNumofChannels, 1+2*AnalogNumofChannels, ...
%       Use double(analogSegments.Data{1}(m:x:n)) to turn segment 1 into array of doubles for
%       channel m where x is num of channels and n is the last sample in range
%
% spiketimes - a structure containing:
%       (1) an array of spike IDs - each id referring to a spike template or type
%       (2) an array of channel IDs - one entry for each spike id
%       (3) a cell array of spike times; each cell contains a vector array of times for each spike ID.
%
% waveforms [if raw waveforms are loaded]: structure containing
%   Channels - structure containing vector arrays describing compression parameters applied to data received from each electrode
%   Times - vector of start times for each data packet
%   Data - cell array of data packets containing waveform data - different channels are interleaved
%   NumOfBuffers - num of data packets or buffers
%   FramesPerBuf - number of samples per buffer
%   NumOfChannels
%   SampleRate in Hz
%
% environment structure containing:
%   Version number of Expo used to export data
%   Filename
%   Conversion parameters
%
% See also GetSlots, GetPasses, GetEvents, GetAnalog, GetSpikeTimes,
% GetSpikeCounts, GetPSTH, PlotPSTH, GetStartTimes, GetEndTimes, GetDuration,
% GetConversionFactor, MergeExpoDataSets, GetTransitionProbabilities.
%
%   Author:      Julian Brown
%	Version:     1.1
%   Last updated:  2005-03-28
%   E-mail:      julian@monkeybiz.stanford.edu

    clear global

    global linenum blocknum slotnum passnum analognum spikenum waveChannelNum waveformnum matlabImportVersion
    global expoDataSet matrix blocks slots passes analog spiketimes environment waveforms
    global currentLine doWaveforms routinesHash routineEntries routineEntryNum routinesMap routineMapNum

    % version info % expoVersion = '1.1';
    % the expoVersion # is stored in the CheckExpoVersion routine
    % we have a separate matlabImport version # so that changes in the Matlab
    % code can be released independently of changes in Expo
    matlabImportVersion = '1.1';

    expoDataSet.MatlabImportVersion = matlabImportVersion;

    disp(sprintf('ReadExpoXML Version %s', matlabImportVersion));

    if ~exist('doSpikeWaveforms'), doSpikeWaveforms = 0;, end

    doWaveforms = doSpikeWaveforms;

    % count items so we can create space for arrays
    if (CountItems(filename)==0), return, end

    % preallocate arrays
    InitializeArrays()

    % extract the data
    ParseXML(filename)

    % save the data as mat file
    basefilename = strtok(filename, '.'); % removes file extension (eg '.xml')

    blocks.routinesMap = routinesMap;
    blocks.routinesEntries = routineEntries;

    FileName = expoDataSet.FileName;
    ExpoVersion = expoDataSet.ExpoVersion;
    MatlabImportVersion = expoDataSet.MatlabImportVersion;

    if isfield(expoDataSet, 'Comment')
		Comment = expoDataSet.Comment;
	else
		expoDataSet.Comment='';
		Comment = '';
	end

    if isfield(expoDataSet, 'Notes')
        Notes = expoDataSet.Notes;
    else
        expoDataSet.Notes = '';
        Notes = '';
    end

    % put data into the structure returned by this function
    expoDataSet.matrix = matrix;
    expoDataSet.blocks = blocks;
    expoDataSet.slots = slots;
    expoDataSet.passes = passes;
    expoDataSet.spiketimes = spiketimes;
    expoDataSet.waveforms = waveforms;
    expoDataSet.analog = analog;
    expoDataSet.environment = environment;

    environment.Conversion.units = InitializeUnitConversionMatrices(expoDataSet);
    expoDataSet.environment = environment;

    save(basefilename, '-struct', 'expoDataSet');

    % show a little of the data
    DisplaySampleData()
return

function ProcessMatrix(fid)
    global linenum currentLine matrix
    matrix = [];
    matrix = ParseNode(matrix, 0);
    
    dimensionNum = 0;
    
    
    while 1
        [token eof] = GetNextToken(fid);
        if eof==1, break, end 
        
        dimension = [];
        
        switch token
            case 'Dimension'
                dimensionNum = dimensionNum + 1;
                dimension = ParseNode(dimension, 0);
                dimensions{dimensionNum} = dimension;
            case '/Matrix'
                matrix.Dimensions = dimensions;
                break;
        end
    end
                         
    
return


function ProcessBlocksAttributes()
    global linenum currentLine blocks matrix
    blockInfo = [];
    blockInfo = ParseNode(blockInfo, 0);

    % this is for v1.0 Expo which put the matrixBasedID in the
    % block data
    if isfield(blockInfo, 'MatrixBaseID')
        matrix.MatrixBaseID = blockInfo.MatrixBaseID;
    end
return;


function ProcessBlock(fid)
    global linenum blocknum blocks currentLine routinesMap routineEntries routineEntryNum routineMapNum routinesHash

    blockRegex = '\<Block ID=\"([\w.]*)\" Name=\"([^\"]*)\" NumOfRoutines=\"([\w.]*)\"';
    routineRegex = '\<Routine ID=\"([\w.]*)\" Name=\"([^\"]*)\" Label=\"([^\"]*)\" NumOfParams=\"([\w.]*)\"';
    paramRegex = '\<Param Name=\"([^\"]*)\" VarName=\"([^\"]*)\" Op=\"([\S]*)\" EUnit=\"([\w.-"]*)\" CUnit=\"([\w.-]*)\" Evalue=\"([\w.-]*)\"';

    blocknum = blocknum+1;
    blockAttributes = GetAttributes(currentLine, 'block', blockRegex, 3);
    blockID = sscanf(blockAttributes{1},'%i'); %ID
    blocks.IDs(blocknum) = blockID;
    blocks.Names{blocknum} = blockAttributes{2}; %Name
    numOfRoutines = sscanf(blockAttributes{3},'%i'); %NumOfRoutines
    routines = cell(numOfRoutines,3);
    routineNum = 0;

    % process routines and params for this block

    while 1
        [token eof] = GetNextToken(fid);
        if eof==1, break, end
        switch token
            case 'Routine'
                routineNum = routineNum + 1; % the routine number within a block
                routineEntryNum = routineEntryNum + 1; % the routine number with a global list of routines
                paramNum = 0;

                routineAttributes = GetAttributes(currentLine, 'routine', routineRegex, 4);
                numOfParams = sscanf(routineAttributes{4},'%i');

                routineID = sscanf(routineAttributes{1},'%i'); % routine ID
                routines.IDs(routineNum)= routineID;
                routineName = routineAttributes{2}; % Name
                routines.Names{routineNum} =routineAttributes{2}; % Name
                routineLabel = routineAttributes{3}; % Label
                routines.Labels{routineNum} = routineLabel; % Label

                % add an entry into a global list of routine entries
                % this list keeps pointers to blocks and routine position
                % and is used in conjunction with routinesMap for finding
                % blockIDs and routine entries by routineID
                routineEntries.RoutineIDs(routineEntryNum) = routineID;
                routineEntries.RoutineNums(routineEntryNum) = routineNum;
                routineEntries.BlockIDs(routineEntryNum) = blockID;
                routineEntries.BlockNums(routineEntryNum) = blocknum;

                % is this a new routine?
%                if routinesHash.containsKey(routineID)
%                    % not new so add entry to routineMap in an existing row
%                    routineMapPos = routinesHash.get(routineID);
%                    routineMapCellSize = size(routinesMap.RoutineNums{routineMapPos});
%                    routinesMap.RoutineNums{routineMapPos}(routineMapCellSize(2) +1) =  routineEntryNum;
%                else
%                    % we have a new routine
%                    routineMapNum = routineMapNum + 1;
%                    % add routine details to the Map
%                    routinesMap.RoutineIDs(routineMapNum) =  routineID;
%                    routinesMap.RoutineNums{routineMapNum}(1) =  routineEntryNum;
%                    routinesMap.RoutineNames{routineMapNum} = routineName;
%                    routinesMap.RoutineLabels{routineMapNum} = routineLabel;
%
%                    % add the routineID and map pointer to hashtable
%                    routinesHash.put(routineID, routineMapNum);
%                end

                %params = cell(numOfParams, 4); % Params cell
                continue

            case '/Routine'
                routines.Params{routineNum} = params;

            case 'Param'

                paramNum = paramNum + 1;
                paramAttributes = GetAttributes(currentLine, 'Param', paramRegex, 6);
                params.Names{paramNum} = paramAttributes{1}; % Name
                params.VarNames{paramNum} = paramAttributes{2}; % VarName
                params.Ops{paramNum} = paramAttributes{3}; % Op
                params.EUnits(paramNum) = sscanf(paramAttributes{4}, '%i'); % EUnit
                params.CUnits(paramNum) = sscanf(paramAttributes{5}, '%i'); % CUnit
                params.Evalues{paramNum} = sscanf(paramAttributes{6},'%f'); % Evalue
                continue

            case '/Block'
                blocks.routines{blocknum} = routines;
                break;
        end
    end
return


function ProcessAnalogData(fid)
    global currentLine linenum analognum analog

    analogRegex = '\<AnalogData NumOfSegments=\"([\w.]*)\" NumOfChannels=\"([\w]*)\" SampleInterval=\"([\w.]*)\"';
    segmentRegex = '\<AnalogSegment NumOfSamples=\"([\w.]*)\" StartTime=\"([\w]*)\">([^\<]*)';

    analogAttributes = GetAttributes(currentLine, 'analog', analogRegex, 3);
    analog.NumOfChannels  =  sscanf(analogAttributes{2},'%i'); %NumOfChannels
    analog.SampleInterval  =  sscanf(analogAttributes{3},'%f'); %SampleInterval

    % process analogSegments
    while 1
        [token eof] = GetNextToken(fid);
        if eof==1, break, end

        switch token
            case 'AnalogSegment'
                analognum = analognum + 1;
                segmentAttributes = GetAttributes(currentLine, 'segment', segmentRegex, 3);
                analog.Segments.SampleCounts(analognum) = sscanf(segmentAttributes{1},'%i'); % NumOfSamples
                analog.Segments.StartTimes(analognum) = sscanf(segmentAttributes{2},'%i'); % StartTime
                analog.Segments.Data{analognum} = base64decode(segmentAttributes{3}); % Data

            case '/AnalogData'
                break
        end
    end
return

function ProcessSlotRecord(fid)
    global currentLine linenum slotnum slots
    slotnum = slotnum+1;
    
    slotRegex = '\<Slot ID=\"([\w.]*)\" Label=\"([^\"]*)\" BlockID=\"([^\"]*)\"';
    slotAttributes = GetAttributes(currentLine, 'slot', slotRegex, 3);
    slots.IDs(slotnum) = sscanf(slotAttributes{1},'%i'); %IDs
    slots.Labels{slotnum} = slotAttributes{2}; %Labels
    slots.BlockIDs(slotnum) = sscanf(slotAttributes{3},'%i'); %BlockIDs
    
    % for backward compatibility with Expo v1.0 check to see whether the
    % node ended here
    if length(regexp(currentLine, '<Slot.*(?=/>)')) > 0 
        return % no more info so return
    end

    while 1
        [token eof] = GetNextToken(fid);
        if eof==1, break, end
        
        info = [];
        
        switch token
            case '/Slot'
                break;
            otherwise
                slots.(token)(slotnum) = ParseNode(info, 0);
        end
    end
    
return

function ProcessPassRecord(fid)
    global currentLine linenum passnum passes
    passnum = passnum+1;
    passRegex = '\<Pass ID=\"([\w.]*)\" SlotID=\"([\w.]*)\" BlockID=\"([\w.]*)\" StartTime=\"([\w.]*)\" EndTime=\"([\w.]*)\"';
    passAttributes = GetAttributes(currentLine, 'pass', passRegex, 5);
    passes.IDs(passnum) = sscanf(passAttributes{1},'%i'); %IDs
    passes.SlotIDs(passnum) = sscanf(passAttributes{2},'%i'); %SlotIDs
    passes.BlockIDs(passnum) = sscanf(passAttributes{3},'%i'); %BlockIDs
    passes.StartTimes(passnum) = sscanf(passAttributes{4},'%i'); %StartTimes
    passes.EndTimes(passnum) = sscanf(passAttributes{5},'%i'); %EndTimes
    ProcessEventRecords(fid);

    % show progress
    if mod(passnum,1000) == 0
        disp(sprintf('%i passes processed', passnum))
    end
return


function ProcessEventRecords(fid)
    global currentLine linenum passnum passes

    % Event elements contain RID,Time,DataParam1;DataParam2
    % where there may be 0 to many DataParams eg.
    % 1,0,11620,1;90;65.479353;49.148398;0;0;0;0;0;270
    eventRegex =  '\<Event RID=\"([\w.]*)\" Time=\"([\w.]*)\" Data=\"([^"]*)\"';
    eventnum = 0;

    while 1
        currentLine = fgets(fid);
        linenum = linenum + 1;
        if findstr('</Pass>', currentLine) > 0, break, end
        eventnum = eventnum+1;
        eventAttributes = GetAttributes(currentLine, 'event', eventRegex, 3);
        events.RoutineIDs(eventnum) =  sscanf(eventAttributes{1},'%i'); %RoutineIDs
        events.Times(eventnum) = sscanf(eventAttributes{2},'%i'); %Times

        % could use the more intuitive
        %events.Data{i} = str2num(eventAttributeValue); but str2num is slower than:
        events.Data{eventnum} = eval(['[' eventAttributes{3} ']']); %Data
    end

    if eventnum==0
        events.RoutineIDs = [];
        events.Times = [];
        events.Data{1} = [];
    end

    passes.events{passnum} = events;
return


function ProcessExpoXData()
    global linenum currentLine expoDataSet matlabImportVersion
    expoDataSet = ParseNode(expoDataSet, 0);
    if ~strcmp(expoDataSet.ExpoVersion, matlabImportVersion)
        warning('The xml was generated from version %s of Expo whereas this Matlab code was written for version %s.',expoDataSet.ExpoVersion, matlabImportVersion)
    end
    disp(sprintf('File name %s', expoDataSet.FileName));
return


function ProcessCommentOrNote(fid, elementName)
    global linenum currentLine expoDataSet
    pos = strfind(currentLine, ['<', elementName, '>']);
    pos = pos + length(elementName) + 2;
    comment = currentLine(pos:length(currentLine));

    while 1
        currentLine = fgets(fid);
        if (currentLine == -1)
            return
        end
        pos = strfind(currentLine, ['</', elementName, '>']);
        if pos>0
            if pos>1
                comment = [comment, currentLine(1:pos-1)];
            end
            expoDataSet.Comment = comment;
            return
        end

        comment = [comment, currentLine];
    end

return


function ProcessEnvironmentRecord(fid)
    global linenum currentLine environment
    currentLine = fgets(fid);
    linenum = linenum + 1;
    environment = ParseNode(environment, 1);
return;


function object = ParseNode(object, useElementName)
    global linenum currentLine
    % generic regular expression that extracts element and arbitrary number of attributes
    genRegex = '\<([\w]*) ([\w]*)=\"([^"]*)\"';


    [tok] = regexp(currentLine, genRegex, 'tokens');
    numOfItems = size(tok, 2);
    elementName = tok{1}(1);

    for i=1:numOfItems
        attributeName = tok{i}(2);
        attributeValue = tok{i}(3);

        if strmatch(attributeName, 'Version')
            forceString = 1; %special handling for the Version tag
            attributeName{1} = 'ExpoVersion';
        elseif strmatch(attributeName, 'Name')
            forceString = 1;
        else
            forceString = 0;
        end

        % dynamically create new members of the object structure
        if useElementName == 1
            if length(str2num(attributeValue{1})) == 1 && ~forceString
                object.(elementName{1}).(attributeName{1}) = str2num(attributeValue{1});
            else
                object.(elementName{1}).(attributeName{1}) = attributeValue{1};
            end
        else
            if length(str2num(attributeValue{1})) == 1 && ~forceString
                object.(attributeName{1}) = str2num(attributeValue{1});
            else
                object.(attributeName{1}) = attributeValue{1};
            end
        end
    end
return


function ProcessSpikeTimes(fid)
    global currentLine linenum spikenum spiketimes
    spikenum = spikenum +1;
    spikeRegex = '\<Spike ID=\"([\w.]*)\" Channel=\"([\w]*)\" Times=\"([\w.,]*)\"';
    spikeAttributes = GetAttributes(currentLine, 'spike times', spikeRegex, 3);
    spiketimes.IDs(spikenum) = sscanf(spikeAttributes{1},'%i'); %ID
    spiketimes.Channels(spikenum) = sscanf(spikeAttributes{2},'%i'); %Channel
    spiketimes.Times{spikenum} = eval(['[' spikeAttributes{3} ']']); %Times
return

function ProcessSpikeWaveforms(fid)
    global currentLine linenum waveChannelNum waveformnum waveforms

    waveformRegex = '\<WaveformRecord Time=\"([\w.]*)\" Buffer=\"([\S]*)\"';
    waveformParamsRegex = '\<SpikeWaveforms NumBuffers=\"([\w.]*)\" FramesPerBuf=\"([\w+=/]*)\" NumChannels=\"([\w]*)\" SampleRate=\"([\w.]*)\"';
    waveformParamsAttributes = GetAttributes(currentLine, 'spike waveform params', waveformParamsRegex, 4);

    waveforms.NumOfBuffers = sscanf(waveformParamsAttributes{1},'%i');
    waveforms.FramesPerBuf = sscanf(waveformParamsAttributes{2},'%i');
    waveforms.NumOfChannels = sscanf(waveformParamsAttributes{3},'%i');
    waveforms.SampleRate = sscanf(waveformParamsAttributes{4},'%f');

    while 1
        [token eof] = GetNextToken(fid);
        if eof==1, break, end

        switch token
            case 'Channel'
                waveChannelNum = waveChannelNum + 1;
                waveformCompressionParamsRegex = '\<Channel Num=\"([\w.]*)\" Scale=\"([\w.]*)\" Offset=\"([\w.-]*)\"';
                waveformCompressionParamAttributes = GetAttributes(currentLine, 'spike waveform channel params', waveformCompressionParamsRegex, 3);

                waveforms.Channels.Scales(waveChannelNum) = sscanf(waveformCompressionParamAttributes{2},'%f'); %Scale
                waveforms.Channels.Offsets(waveChannelNum) = sscanf(waveformCompressionParamAttributes{3},'%f'); %Offset

            case 'WaveformRecord'
                waveformnum = waveformnum +1;
                waveformAttributes = GetAttributes(currentLine, 'waveform', waveformRegex, 2);
                waveforms.Times(waveformnum) = sscanf(waveformAttributes{1},'%f'); %Time
                waveforms.Data{waveformnum} = base64decode(waveformAttributes{2}); %Waveform data

                % show progress
                if mod(waveformnum,1000) == 0
                    disp(sprintf('%i waveform records processed. %%%i complete', waveformnum, int16(double(waveformnum)*100/double(length(waveforms.Times)))))
                end

            case '/SpikeWaveforms'
                break
        end


    end
return

function [token eof] = GetNextToken(fid)
    global currentLine linenum

    %use regex to extract token ELEMENTNAME from tag like <ELEMENTNAME BLAH="123" />
     % or <ELEMENTNAME BLAH="123"> or \ELEMENTNAME from </ELEMENTNAME>
    regex = '<(/?\w+).*?/?>';

    while 1
        currentLine = fgets(fid);
        if (currentLine == -1)
            eof = 1;
            token ='';
            return
        else
            eof = 0;
        end

        linenum = linenum + 1;
        [tok] = regexp(currentLine, regex, 'tokens');

        %check we found a token otherwise ignore and continue
        if size(tok,2) == 0, continue, end

        token = tok{1,1}{1,1};
        break
    end
return


function success = CountItems(filename)
    global linenum blocknum slotnum passnum analognum spikenum spiketimes waveChannelNum waveformnum

    % open the file for reading
    fid = fopen(filename,'rt');   % 'rt' means read text
    if (fid < 0)
        error('could not open file');% just abort if error
        success = 0;
        return
    end

    linenum = 0;
    blocknum = 0;
    slotnum = 0;
    passnum = 0;
    analognum = 0;
    spikenum = 0;
    waveChannelNum = 0;
    waveformnum = 0;

    while 1
        [token eof] = GetNextToken(fid);
        if eof==1 , break, end

        if mod(linenum, 10000) == 0
            disp(sprintf('%i lines counted', linenum))
        end

        switch token
            case 'Block'
                blocknum = blocknum + 1;
            case 'Slot'
                slotnum = slotnum + 1;
            case 'Pass'
                passnum = passnum + 1;
            case 'AnalogSegment'
                analognum = analognum + 1;
            case 'Spike'
                spikenum = spikenum + 1;
            case 'Channel'
                waveChannelNum = waveChannelNum + 1;
            case 'WaveformRecord'
                waveformnum = waveformnum + 1;
        end
    end

    fclose(fid);

    display(sprintf('Number of blocks %d', blocknum));
    display(sprintf('Number of slots %d', slotnum));
    display(sprintf('Number of passes %d', passnum));
    display(sprintf('Number of analogs segments %d', analognum));
    display(sprintf('Number of spike templates %d', spikenum));
    display(sprintf('Number of waveform channels %d', waveChannelNum));
    display(sprintf('Number of spike waveform records %d', waveformnum));

    success = 1;
return

function InitializeArrays()
    global linenum blocknum slotnum passnum analognum spikenum waveChannelNum waveformnum doWaveforms
    global routineEntries routineEntryNum routineMapNum routinesMap environment
    global blocks slots passes analog spiketimes waveforms routinesHash

    blocks.IDs = zeros(1, blocknum, 'int16');
    blocks.Names = cell(1, blocknum);
    blocks.routines = cell(1, blocknum);

    slots.IDs = zeros(1, slotnum, 'int16');
    slots.Labels = cell(1, slotnum);
    slots.BlockIDs = zeros(1, slotnum, 'int16');

    passes.IDs =  zeros(1, passnum, 'int32');
    passes.SlotIDs =  zeros(1, passnum, 'int16');
    passes.BlockIDs =  zeros(1, passnum, 'int16');
    passes.StartTimes =  zeros(1, passnum, 'int32');
    passes.EndTimes =  zeros(1, passnum, 'int32');
    passes.events =  cell(1, passnum);

    analog.Segments.SampleCounts = zeros(1, analognum, 'int32');
    analog.Segments.StartTimes = zeros(1, analognum, 'int32');
    analog.Segments.Data = cell(1, analognum);

    spiketimes.IDs = zeros(1, spikenum, 'int16');
    spiketimes.Channels = zeros(1, spikenum, 'int16');
    spiketimes.Times = cell(1, spikenum);

    doWaveforms = doWaveforms && (waveformnum>0);

    if doWaveforms
        waveforms.Channels.Scales = zeros(waveChannelNum, 1, 'double');
        waveforms.Channels.Offsets = zeros(waveChannelNum, 1, 'double');
        waveforms.Times = zeros(waveformnum, 1, 'double');
        waveforms.Data = cell(waveformnum, 2);
    end

    routinesMap = cell(1,4);
    routineEntries = zeros(1, 4, 'int32');

%    routinesHash = java.util.HashMap();

    linenum = 0;
    blocknum = 0;
    slotnum = 0;
    passnum = 0;
    analognum = 0;
    spikenum = 0;
    waveChannelNum = 0;
    waveformnum = 0;
    routineEntryNum = 0;
    routineMapNum = 0;
return


function attributes = GetAttributes(currentLine, elementName, regex, numOfTokens)
    global linenum

    [tok] = regexp(currentLine, regex, 'tokens');
    numOfItems = int32(size(tok{1,1},2));

    if (numOfItems~=numOfTokens)
        error(sprintf('ERROR! Incorrect number of attributes found for %s at line %d', elementName, linenum));
    end

    attributes = tok{1,1};
return


function ParseXML(filename)
    global linenum blocknum slotnum spikenum passnum analognum routineEntryNum
    global doWaveforms

    %open file for reading again
    fid = fopen(filename, 'rt');

    while 1
        [token eof] = GetNextToken(fid);
        if eof==1 , break, end

        switch token
            case 'Matrix'
                ProcessMatrix(fid);
            case 'Blocks'
                ProcessBlocksAttributes();
            case 'Block'
                ProcessBlock(fid);
            case 'Slot'
                ProcessSlotRecord(fid);
            case 'Pass'
                ProcessPassRecord(fid);
            case 'Comment'
                ProcessCommentOrNote(fid, 'Comment')
            case 'Note'
                ProcessCommentOrNote(fid, 'Note')
            case 'Environment'
                ProcessEnvironmentRecord(fid);
            case 'Spike'
                ProcessSpikeTimes(fid);
            case 'SpikeWaveforms'
                if doWaveforms, ProcessSpikeWaveforms(fid), end;
            case 'AnalogData'
                ProcessAnalogData(fid);
            case 'ExpoXData'
                ProcessExpoXData();
        end
    end

    % close the file
    fclose(fid);
return


function DisplaySampleData()
    global blocks slots passes routinesMap
    global blocknum slotnum analognum analog spikenum passnum spiketimes
    global waveforms waveformnum waveChannelNum doWaveforms

    disp(' ')

    %display sample of block data
    if blocknum >0
        for i=1:min(2,blocknum)
            disp(sprintf('Block ID=%d\tName=%s',blocks.IDs(i), blocks.Names{i}));
            disp(blocks.routines{i})
        end
    end


    %display slot data
    if slotnum>0
        for i=1:min(3,slotnum)
            disp(sprintf('Slot ID=%d\t Label=%s \tBlockID=%d',slots.IDs(i), slots.Labels{i}, slots.BlockIDs(i)));
        end

        disp(' ')
    end

    disp('Routine Map')
    disp(routinesMap);

    disp(' ')

    %display sample of pass data
    if passnum > 1
        for i=1:2
            disp(sprintf('Pass ID=%d\t SlotID=%d \t BlockID=%d \tStartTime=%d \tEndTime=%d',passes.IDs(i), passes.SlotIDs(i), passes.BlockIDs(i), passes.StartTimes(i), passes.EndTimes(i)));

            %display the event data for this pass
            for j=1:length(passes.events{i}.RoutineIDs)
                disp(sprintf('Event %d\tRoutineID=%d\tTime=%d\tData=%s',j, passes.events{i}.RoutineIDs(j), passes.events{i}.Times(j), mat2str(passes.events{i}.Data{j}')));
            end
            disp(' ')
        end
    end

    %display a sample of the analog data
    if analognum > 0
        disp(sprintf('Analog NumOfChannels=%d\tSampleInterval=%d',analog.NumOfChannels, analog.SampleInterval));

        %for i=1:1
            disp(sprintf('Segment 1 NumOfSamples=%d\tStartTime=%d\t',analog.Segments.SampleCounts(1), analog.Segments.StartTimes(1)));
            disp(int2str(analog.Segments.Data{1}(1:16)))
        %end

        disp(' ')
    end

    %display a sample of spike times data
    if spikenum > 0
        j = 1; %pick first template
        disp(sprintf('Spike\tTemplateID=%d\tChannelID=%d', spiketimes.IDs(j), spiketimes.Channels(j)));
        if length(spiketimes.Times{j, 1}) >= 5
            for i=1:5
                disp(sprintf('Time=%d', spiketimes.Times{j, 1}(i)));
            end
        end
    end

    if ~doWaveforms, return, end;
    % display a spike waveform record
    if waveformnum>0
        disp(sprintf('Waveforms NumBuffers=%d\tFramesPerBuf=%d NumChannels=%d SampleRate=%d',waveforms.NumOfBuffers, waveforms.FramesPerBuf, waveforms.NumOfChannels, waveforms.SampleRate));

        for i=1:waveChannelNum
            disp(sprintf('Channel Num=%d Scale=%d Offset=%d', i, waveforms.Channels.Scales(i), waveforms.Channels.Offsets(i)));
        end

        y = int16(waveforms.Data{1});
        x = 1:waveforms.NumOfChannels:length(y);
        figure;
        plot(x, waveforms.Channels.Offsets(1) + double(y(x))/waveforms.Channels.Scales(1));
        title('Spike waveform record 1')
    end
return