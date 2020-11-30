%NEX2STRO   Convert NeuroEXplorer files to a STRO structure.
%   STRO = NEX2STRO() prompts the user for the location of a .nex file and
%   returns a MATLAB structure with fields 'sum', 'trial', 'ras', and 'other'.
%   These fields are populated with information from the .nex file and contain
%   data summaries, trial-by-trial parameters, rasters, and any other data. The
%   STRO structure is the standard MATLAB-based representation of all
%   electrophysiological and behavioral data in the Horwitz lab.
%
%   STRO = NEX2STRO(FILENAME) is the same as the above but instead uses the
%   .nex file specified by the absolute path string FILENAME.

% cah       1/21/08 started development...
% gdlh/cah  2/14/08 replaced global variable architecture with nested
%                   functions.
% zalb     10/09/13 removing all non-local variable scopes; boosting perf

function stro = nex2stro(file_name)
stro = [];
if nargin == 0 || isempty(file_name)
    [fname, pathname] = uigetfile([nexfilepath,filesep,'*.nex'], 'Select a NeuroEXplorer file');
    if isequal(fname,0) || isequal(pathname,0)
        return
    end
    file_name = [pathname fname];
elseif exist(file_name, 'dir')
    partial_path = file_name;
    [fname, pathname] = uigetfile('*.nex', 'Select a NeuroEXplorer file', partial_path);
    if isequal(fname,0) || isequal(pathname,0)
        return
    end
    file_name = [pathname fname];
elseif ~exist(file_name, 'file')
    error('Couldn''t find file %s\n', file_name);
end

[stro, spikes, anlg, waves, events] = setupStructures();
stro.sum.fileName = file_name;

[ecodes, spikes, anlg, waves, events] = unpackNexFile(stro.sum.fileName, spikes, anlg, waves, events);
[anlg, events] = cleanupAnalogSignals(anlg, events); % determine if the start times are identical for each channel

% determine which paradigm file to use, then eval the appropriate header
% file. Be sure to return all the appropriate variables (for scoping
% purposes). the feval statement will crash if one of the outputs is
% undefined in paradigmFile
paradigmID = dat2num(ecodes, 8999, 4000, 'int', 0);
stro.sum.paradigmID = paradigmID{1}(1);
paradigmFile = paradigmLibrary(stro.sum.paradigmID);


[C, trialParamsHeader, exptParamsHeader, additionalRasterFields, badLogic, otherInstructions] = feval(paradigmFile);

% Fill up some of the stro.sum contents
stro = getExperimentalParams(stro, ecodes, C.VALOFFSET, exptParamsHeader); % iterates over the exptParamsHeader reconstructing those elements
stro.sum.rasterCells = makeRasterHeader(spikes, anlg, waves, additionalRasterFields); % defined in the 'paradigmFile'.
stro.sum.trialFields = trialParamsHeader; % defined in the 'paradigmFile'
stro.sum.otherFields = otherInstructions;
stro.sum.analog.sigid = anlg.names;
stro.sum.analog.storeRates = anlg.ADFrequency;
stro.sum.analog.ADtoMV = anlg.ADtoMV;
stro.sum.waves.wf_id = waves.names;
stro.sum.waves.storeRates = waves.WFrequency;

% iterate by trials filling up the .raster and .trial fields of stro.
nTrials = length(events.start);
stro.ras = cell(nTrials, length(stro.sum.rasterCells));
stro.trial = nan(nTrials, size(trialParamsHeader, 2));
stro.sum.absTrialNum = nan(nTrials,1);
goodTrial = 0; % the counter of good trials

for trial_num = 1:nTrials
    trialCodes = checkTrial(ecodes, events.start(trial_num), events.stop(trial_num), C, badLogic);
    if ~isempty(trialCodes)
        goodTrial = goodTrial + 1;
        stro = rasterByTrial(stro, spikes, anlg, waves, C.VALOFFSET, ...
            events.start(trial_num), events.stop(trial_num), goodTrial, ...
            trial_num, additionalRasterFields, trialCodes);
        stro = indexByTrial(stro, trialCodes, trialParamsHeader, goodTrial, C.VALOFFSET);
        if ~isempty(otherInstructions)
            stro = otherByTrial(stro, trialCodes, otherInstructions, goodTrial, C);
        end
        stro.sum.absTrialNum(goodTrial) = trial_num; % a mapping b/w actual and good trial nums
    end
end

% now cleanup the bad trials (empty rows) of stro.trial, stro.ras, and
% stro.sum.absTrialNum
stro.ras(goodTrial+1:end, :) = [];
stro.trial(goodTrial+1:end, :) = [];
stro.sum.absTrialNum(goodTrial+1:end) = [];

function [stro, spikes, anlg, waves, events] = setupStructures()
% this function is necessary because later on we'll test the length of
% spikes.names and anlg.names to iterate over neuronal and anlg data. If we
% decide not to record from an anlg channel, setting the field anlg.name to an
% empty vector (i.e., []) will allow this program to function properly. ditto
% for spike channels. If future technology allows us to record from more spike
% or analog channels, change the constants below to reflect that.

SPIKE_CHANNEL_MAX = 64;
ANLG_CHANNEL_MAX = 64;

spikes = struct( ...
            'names', {cell(1,SPIKE_CHANNEL_MAX)}, ...
            'timestamps', {cell(1,SPIKE_CHANNEL_MAX)} ...
         );

anlg = struct( ...
          'names', {cell(1,ANLG_CHANNEL_MAX)}, ...
          'ADFrequency', {cell(1,ANLG_CHANNEL_MAX)}, ...
          'ADtoMV', {cell(1,ANLG_CHANNEL_MAX)}, ...
          'fragStarts', {cell(1,ANLG_CHANNEL_MAX)}, ...
          'nPtsPerFrag', {cell(1,ANLG_CHANNEL_MAX)}, ...
          'data', {cell(1,ANLG_CHANNEL_MAX)} ...
       );

waves = struct( ...
           'names', {{}}, ...
           'WFrequency', {{}}, ...
           'timestamps', {{}}, ...
           'waveforms', {{}} ...
        );

events = struct('start', [], 'stop', []);

stro = struct('sum', struct( ...
                        'fileName', [], ...
                        'date', date, ...
                        'paradigmID', [], ...
                        'waves', struct( ...
                                    'wf_id', [], ...
                                    'storeRates', [] ...
                                 ), ...
                        'analog', struct( ...
                                     'sigid', [], ...
                                     'storeRates', [], ...
                                     'ADtoMV', [] ...
                                  ), ...
                        'trialFields', {{}}, ...
                        'rasterCells', {{}}, ...
                        'otherFields', {{}}, ...
                        'exptParams', struct(), ...
                        'absTrialNum', [], ...
                        'concat', [] ...
                     ), ...
              'trial', {{}}, ...
              'ras', {{}}, ...
              'other', {{}} ...
       );

% This function is based off of readNexFile.m from http://www.neuroexplorer.com/code.html
function [ecodes, spikes, anlg, waves, events] = unpackNexFile(fileName, spikes, anlg, waves, events)
% enforce a little endian ('l') and ASCII character encoding
fid = fopen(fileName, 'r', 'l', 'US-ASCII');

if fid == -1
    error('Unable to open file');
end

% closes the file when `close_obj` goes out of scope (when this function aborts
% or returns normally)
close_obj = onCleanup(@() fclose(fid));

% parse the nex file's header:
% struct NexFileHeader
% {
%     int MagicNumber;
%     int NexFileVersion;
%     char Comment[256];
%     double Frequency;
%     int Beg;
%     int End;
%     int NumVars;
%     int NextFileHeader;
%     char Padding[256];
% };

magic = fread(fid, 1, 'int32');
if magic ~= 827868494
    error('The file is not a valid .nex file');
end

fseek(fid, 260, 'cof'); % skip over NexFileVersion, Comment
freq = fread(fid, 1, 'double');
fseek(fid, 8, 'cof'); % skip over Beg, End
nvar = fread(fid, 1, 'int32');
fseek(fid, 260, 'cof'); % skip over NextFileHeader, Padding

% now read the file putting the appropriate 'type' of Nex info in either
% the ecodes, spikes, or analog vectors. For now ignore the other
% 'types'

neuronCount = 0;
contSigCount = 0;
waveCount = 0;

% here's the layout for variables:
% struct NexVarHeader
% {
%     int Type;
%     int Version;
%     char Name[64];
%     int DataOffset;
%     int Count;
%     int WireNumber;
%     int UnitNumber;
%     int Gain;
%     int Filter;
%     double XPos;
%     double YPos;
%     double WFrequency;
%     double ADtoMV;
%     int NPointsWave;
%     int NMarkers;
%     int MarkerLength;
%     double MVOffset;
%     double PrethresholdTimeInSeconds;
%     char Padding[52];
% };

while nvar
    type = fread(fid, 1, 'int32');
    fseek(fid, 4, 'cof');
    name = deblank(fread(fid, 64, '*char')');
    data_location = fread(fid, 1, 'int32');
    n = fread(fid, 1, 'int32');
    fseek(fid, 32, 'cof');
    
    % only read in what's necessary; fseek over those bytes otherwise
    if type == 3 || type == 5
        WFrequency = fread(fid, 1, 'double');
        ADtoMV = fread(fid, 1, 'double');
        NPointsWave = fread(fid, 1, 'int32');
        fseek(fid, 8, 'cof');
        MVOffset = fread(fid, 1, 'double');
    elseif type == 6
        fseek(fid, 20, 'cof');
        NMarkers = fread(fid, 1, 'int32');
        MarkerLength = fread(fid, 1, 'int32');
        fseek(fid, 8, 'cof');
    else
        fseek(fid, 36, 'cof');
    end
    
    fseek(fid, 60, 'cof');
    next_header = ftell(fid);
    fseek(fid, data_location, 'bof');
    
    if type == 0 % neuron
        neuronCount = neuronCount+1;
        spikes.names{neuronCount} = name;
        spikes.timestamps{neuronCount} = fread(fid, n, 'int32') / freq;
    elseif type == 1 % event (start/stop)
        if ~strcmpi(name, 'event001') % don't care about event001s
            timeStamps = fread(fid, n, 'int32') / freq;
            events.(lower(name)) = timeStamps;
        end
    elseif type == 2 % interval
        % ignore this type for now
        fprintf('Ignoring the ''interval'' type\n');
    elseif type == 3 % waveform
        waveCount = waveCount+1;
        waves.names{waveCount} = name;
        waves.WFrequency{waveCount} = WFrequency;
        waves.timestamps{waveCount} = fread(fid, n, 'int32') / freq;
        waves.waveforms{waveCount} = (fread(fid, [NPointsWave n], 'int16') * ADtoMV + MVOffset)'; % transpose to make trials go down columns and time across rows
    elseif type == 4 % population vector
        % ignore this type for now
        fprintf('Ignoring the ''population vector'' type\n');
    elseif type == 5 % continuous variable (i.e. eye signals)
        contSigCount = contSigCount+1;
        anlg.names{contSigCount} = name;
        anlg.ADtoMV{contSigCount} = ADtoMV;
        anlg.ADFrequency{contSigCount} = WFrequency;
        
        % get the start times of each analog signal snippet
        anlg.fragStarts{contSigCount} = fread(fid, n, 'int32') / freq;
        
        % get the number of points sammpled during each of the snippets
        nFrag = [fread(fid, n, 'int32'); 0];
        nFrag(n+1) = NPointsWave;
        anlg.nPtsPerFrag{contSigCount} = nFrag;
        % now bring in the AtoD data for the entire recording
        anlg.data{contSigCount} = fread(fid, [NPointsWave 1], 'int16') * ADtoMV + MVOffset;
    elseif type == 6 % marker (i.e. ecodes)
        timeStamps = fread(fid, n, 'int32') / freq;
        
        % check to make sure that there are actually markers to retreive before
        % trying
        if NMarkers ~= 1
            error('bad number of markers in file');
        end
        
        % check to make sure that you're about to read in strobed codes
        mname = fread(fid, 64, '*char');
        if ~strncmp('DIO', mname, 3)
            error('unknown marker name');
        end
        
        % now collect the marker (ecode) data. convert them to numeric
        % representations and return them along with timestamps
        markers = deblank(fread(fid, [MarkerLength n], '*char')');
        markers = markers - 48; % convert to numeric
        powersOfTen = 10.^((MarkerLength-2):-1:0)';
        ecodes = [timeStamps markers * powersOfTen];
    else
        fprintf('unknown variable type <%s>\n', type);
    end
    
    fseek(fid, next_header, 'bof'); % go back to the next variable's header
    nvar = nvar - 1;
end

spikes.names = spikes.names(1:neuronCount);
spikes.timestamps = spikes.timestamps(1:neuronCount);
anlg.names = anlg.names(1:contSigCount);
anlg.ADFrequency = anlg.ADFrequency(1:contSigCount);
anlg.ADtoMV = anlg.ADtoMV(1:contSigCount);
anlg.fragStarts = anlg.fragStarts(1:contSigCount);
anlg.nPtsPerFrag = anlg.nPtsPerFrag(1:contSigCount);
anlg.data = anlg.data(1:contSigCount);

function trialCodes = checkTrial(ecodes, trialStart, trialStop, C, badLogic) %#ok<INUSL> % need access to C for `eval`
% pull out the relevant ecodes
ind = ecodes(:,1) >= trialStart & ecodes(:,1) <= trialStop;
trialCodes = ecodes(ind,:);

% test the ecodes to determine if the trial is good or bad. 'badLogic'
% should be defined in the paradigm header file.
bad = eval(badLogic);

% return an empty vector if it's a bad trial
if bad
    trialCodes = [];
end

function stro = rasterByTrial(stro, spikes, anlg, waves, VALOFFSET, start, stop, gtCounter, totTrialCounter, additionalRasterFields, trialCodes)
% fill up the rows of stro.ras trial by trial. treat spike data, anlg sigs, and
% other independently. Format is:
% <neuron 1> <neuron 2> ... <waveForm1>... <anlg 1> <anlg 2> ... <atime> <other>

% each neuron's spikes will go into it's own column:
if ~isempty(spikes.names)
    for ind = 1:length(spikes.names)
        stro.ras{gtCounter, ind} = spikes.timestamps{ind}(spikes.timestamps{ind} >= start & spikes.timestamps{ind} <= stop);
    end
else
    ind = 0;
end

% iterate over the waveform data (if available):
for a = 1:length(waves.names)
    ind = ind + 1;
    whichSpikes = (waves.timestamps{a} >= start & waves.timestamps{a} <= stop);
    stro.ras{gtCounter, ind} = waves.waveforms{a}(whichSpikes,:);
end

% iterate over anlg channels
for a = 1:length(anlg.names)
    ind = ind + 1;
    startInd = anlg.nPtsPerFrag{a}(totTrialCounter) + 1;
    stopInd = anlg.nPtsPerFrag{a}(totTrialCounter+1);
    stro.ras{gtCounter, ind} = anlg.data{a}(startInd:stopInd);
end

% explicitly add the anlg start time. All the analog channels should have the
% same number of start times, so indexing this way is okay
if ~isempty(anlg.fragStarts)
    ind = ind + 1;
    stro.ras{gtCounter, ind} = anlg.fragStarts{1}(totTrialCounter);
end

% now look for additional elements as specified by additionalRasterFields (which
% is defined in the paradigm's 'codes' file)
for a = 1:size(additionalRasterFields, 2)
    ind = ind + 1;
    nums = dat2num(trialCodes(:,2), [additionalRasterFields{4, a}], VALOFFSET, additionalRasterFields{2, a}, additionalRasterFields{3, a});
    stro.ras{gtCounter, ind} = nums{1};
end

function stro = indexByTrial(stro, trialCodes, trialParamsHeader, goodTrial, VALOFFSET)
for T = {'int' 'double' 'float' 'long' 'char' 'short' 'ushort' 'ulong' 'uint'}
    typeInd = strcmp(T{:}, trialParamsHeader(2,:));
    if any(typeInd)
        nums = dat2num(trialCodes(:,2), [trialParamsHeader{4, typeInd}], VALOFFSET, T{:}, 0);
        stro.trial(goodTrial, typeInd) = [nums{:}];
    end
end

% now deal with the 'time' variables
timeInd = find(strcmp('time', trialParamsHeader(2,:)));
timeInd_code = [trialParamsHeader{4, timeInd}];
for a = 1:length(timeInd)
    time = trialCodes(trialCodes(:,2) == timeInd_code(a), 1); % there will be an error if numel(time) > 1
    if ~isempty(time)
        if numel(time) > 1, warning('More than one time code found: using the last!!'); end
        stro.trial(goodTrial, timeInd(a)) = time(end);
    end
end

function stro = otherByTrial(stro, trialCodes, otherInstructions, gtCounter, codeStruct)
for a = 1:length(otherInstructions)
    stro.other{gtCounter, a} = feval(otherInstructions{2,a}, trialCodes, codeStruct);
end

function rasterCells = makeRasterHeader(spikes, anlg, waves, additionalRasterFields)
% start with the spike data then add the anlg data. this means that
% stro.raster will be composed of columns:
% <neuron 1> <neuron 2> ... <neuron 3> <anlg 1> <anlg 2> ... <anlg 3> <atime> <other>
if ~isempty(anlg.fragStarts)
    anlgStartField = 'anlgStartTime';
else
    anlgStartField = [];
end

if ~isempty(additionalRasterFields)
    rasterCells = [spikes.names waves.names anlg.names anlgStartField additionalRasterFields(1,:)];
else
    rasterCells = [spikes.names waves.names anlg.names anlgStartField];
end

function stro = getExperimentalParams(stro, ecodes, VALOFFSET, exptParamsHeader)
% deal with the entries into the stro.sum.exptParams one by one.
for a = 1:size(exptParamsHeader, 2)
    nums = dat2num(ecodes, [exptParamsHeader{4, a}], VALOFFSET, ...
        exptParamsHeader{2, a}, exptParamsHeader{3, a});
    stro.sum.exptParams.(exptParamsHeader{1,a}) = nums{1};
end

function [anlg, events] = cleanupAnalogSignals(anlg, events)
% check to make sure that the number of analog fragments is the same for
% each channel. ditto for the number of AtoD points per fragment

% Ensure an equal number of starts and stops (i.e. check to make sure
% that data acquisition was not arrested part way through a trial)
excessTrial = length(events.start) - length(events.stop);
if excessTrial == 1
    events.start(end) = [];
    fprintf(' Unequal number of STARTs and STOPs. Deleted the terminal START.\n');
elseif excessTrial == -1 && events.stop(1) < events.start(1)
    events.stop(1) = [];
    fprintf(' Unequal number of STARTs and STOPs. Deleted the initial STOP.\n');
elseif excessTrial ~= 0
    error(' Unequal number of STARTS and STOPS.');
end

% now determine if there are the same number of trial starts/stops as
% there are anlg starts/stops. delete any anlg data that occurs before
% the first trial.
if isempty(anlg.fragStarts)
    return
end

err = abs(events.start(1) - anlg.fragStarts{1});
startIdx = find(err == min(err), 1);

if startIdx > 1
    for a = 1:length(anlg.names)
        anlg.fragStarts{a}(1:startIdx-1) = [];
        anlg.nPtsPerFrag{a}(1:startIdx-1) = [];
    end
end

if numel(anlg.names) > 1
    if ~isequal(anlg.fragStarts{:})
        error('frag starts are unequal');
    end
    
    if ~isequal(anlg.nPtsPerFrag{:})
        error('points per frag are unequal');
    end
end
