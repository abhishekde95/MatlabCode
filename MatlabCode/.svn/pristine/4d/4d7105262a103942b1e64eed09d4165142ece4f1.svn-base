function expoDataSet = MergeExpoDataSets(filename, expoDataSets)
% function expoDataSet = MergeExpoDataSets(filename, expoDataSets)
%
% Expo utility to combine the data sets in the cell array expoDataSets, save them in filename, and
% return the combined expoDataSet structure, whose filename field is filename.  
% The input expoDataSets must come from the same Expo program, using the same Environment settings. 
% Each of the expoDataSets can be specified as and expoDataSet structure or as the name of a file 
% containing an expoDataSet.
% The Comment and Notes fields of the combined structure are not merged -- rather, they contain the
% Comment and Notes from the first of the expoDataSets.
%
% If filename is [] then no file is saved.
% The oldstyle call e = MergeExpoDataSets(e1, e2) takes exactly two
% expoDataSets, which must be structures, merges them, and returns them without saving.
%
% See also ReadExpoXML, GetSlots, GetPasses, GetEvents, GetSpikeTimes, GetAnalog,
% GetPSTH, PlotPSTH, GetStartTimes, GetEndTimes, GetDuration,
% GetConversionFactor.
%
%   Author:      Julian Brown
%	Version:     1.1
%   Last updated:  2005-3-28
%   E-mail:      julian@monkeybiz.stanford.edu

if isempty(filename)
    filename = ''; 
end 
if isstruct(filename) & isstruct(expoDataSets)% Old-style call.
    expoDataSet = Merge2ExpoDataSets(filename, expoDataSets); 
    filename = '';
else % Cell array
    for idx = 1:length(expoDataSets)
        if ~isstruct(expoDataSets{idx})
            expoDataSets{idx} = load(expoDataSets{idx});
        end
    end
    while length(expoDataSets) > 2
        expoDataSets = { Merge2ExpoDataSets(expoDataSets{1}, expoDataSets{2}) expoDataSets{3:end}};
    end
    expoDataSet = Merge2ExpoDataSets(expoDataSets{1}, expoDataSets{2});
end
% Save it
if ~strcmp('',filename)
    expoDataSet.FileName = filename;
    save(filename,'-struct', 'expoDataSet')
end




function ds = Merge2ExpoDataSets(ds1, ds2)    

%MatlabImportVersion
if ~strcmp(ds1.MatlabImportVersion,ds2.MatlabImportVersion)
    error('The two datasets have different MatlabImportVersion numbers.  Cannot merge')
end
ds.MatlabImportVersion = ds1.MatlabImportVersion;

%ExpoVersion
if ~strcmp(ds1.ExpoVersion,ds2.ExpoVersion)
    error('The two datasets have different ExpoVersion numbers.  Cannot merge')
end
ds.ExpoVersion = ds1.ExpoVersion;

matlabImportVersion = '1.1';
CheckExpoVersion(ds1, matlabImportVersion);


%FileName
if ~strcmp(ds1.FileName,ds2.FileName)
    error('The two datasets have different FileNames.  Cannot merge')
end
ds.FileName = ds1.FileName;

%Comment
if ~strcmp(ds1.Comment,ds2.Comment)
    warning('The two datasets have different Comment fields.  The new dataset will have the Comment of the first dataset')
end
ds.Comment = ds1.Comment;

%Notes
if ~strcmp(ds1.Notes,ds2.Notes)
    warning('The two datasets have different Notes fields.  The new dataset will have the Notes of the first dataset')
end
ds.Notes = ds1.Notes;

%check the datasets have similar block structures
if length(ds1.blocks.IDs) ~= length(ds2.blocks.IDs)
    error('The two datasets have different numbers of blocks.  Cannot merge')
end
ds.blocks = ds1.blocks;


%slots
if length(ds1.slots.IDs) ~= length(ds2.slots.IDs)
    error('The two datasets have different numbers of slots.  Cannot merge')
end
ds.slots = ds1.slots;

% passes
numOfPassesInds1 = length(ds1.passes.IDs);
startTimeForSecondDataSet = ds1.passes.EndTimes(numOfPassesInds1) + 1;

passFieldNames = fieldnames(ds1.passes);

for i=1:length(passFieldNames)
    field = passFieldNames{i};

    switch field
        case 'IDs'
            ds.passes.IDs = [ds1.passes.IDs ds2.passes.IDs + numOfPassesInds1];
        case 'StartTimes'
            ds.passes.StartTimes = [ds1.passes.StartTimes ds2.passes.StartTimes + startTimeForSecondDataSet];
        case 'EndTimes'
            ds.passes.EndTimes = [ds1.passes.EndTimes ds2.passes.EndTimes + startTimeForSecondDataSet];
        otherwise
            ds.passes.(field) = [ds1.passes.(field) ds2.passes.(field)];
    end
end

%spiketimes
numOfSpikeIDs = length(ds1.spiketimes.IDs);

if numOfSpikeIDs ~= length(ds1.spiketimes.IDs)
    error('The two datasets have different numbers of spike IDs.  Cannot merge')
end

ds.spiketimes.IDs = ds1.spiketimes.IDs;
ds.spiketimes.Channels = ds1.spiketimes.Channels;

for i=1:numOfSpikeIDs
    ds.spiketimes.Times{i} = [ds1.spiketimes.Times{i} ds2.spiketimes.Times{i} + double(startTimeForSecondDataSet)];
end

%waveforms
if isfield(ds1.waveforms, 'NumOfChannels') && isfield(ds2.waveforms, 'NumOfChannels')
    waveformFieldNames = fieldnames(ds1.waveforms);

    for i=1:length(waveformFieldNames)
        field = waveformFieldNames{i};

        switch field
            case 'NumOfChannels'
                if length(ds1.waveforms.NumOfChannels) ~= length(ds2.waveforms.NumOfChannels)
                    error('The two datasets have different numbers of waveform channels.  Cannot merge')
                end
                ds.waveforms.NumOfChannels = ds1.waveforms.NumOfChannels;
            case 'FramesPerBuf'
                if length(ds1.waveforms.FramesPerBuf) ~= length(ds2.waveforms.FramesPerBuf)
                    error('The two datasets differ in the FramesPerBuf for the waveforms.  Cannot merge')
                end
                ds.waveforms.FramesPerBuf = ds1.waveforms.FramesPerBuf;
            case 'SampleRate'
                if length(ds1.waveforms.SampleRate) ~= length(ds2.waveforms.SampleRate)
                    error('The two datasets differ in the SampleRate for the waveforms.  Cannot merge')
                end
                ds.waveforms.SampleRate = ds1.waveforms.SampleRate;
            case 'NumOfBuffers'
                ds.waveforms.NumOfBuffers = ds1.waveforms.NumOfBuffers + ds2.waveforms.NumOfBuffers;
            case 'Times'
                ds.waveforms.Times = [ds1.waveforms.Times + ds2.waveforms.Times];
            case 'Data'
                for j=1:ds1.waveforms.NumOfChannels
                    o1 = ds1.waveforms.Channels.Offsets(j);
                    o2 = ds2.waveforms.Channels.Offsets(j);
                    s1 = ds1.waveforms.Channels.Scales(j);
                    s2 = ds2.waveforms.Channels.Scales(j);
                    boundaries = [o1 o1 + 255/s1 o2 o2 + 255/s2];

                    o3 = min(boundaries);
                    s3 = 255/(max(boundaries) - o3);

                    ds.waveforms.Channels.Offsets(j) = o3;
                    ds.waveforms.Channels.Scales(j) = s3;

					arrayEndPoint = size(ds1.waveforms.Data,1);
                    for k =1:arrayEndPoint
						bytes = floor(s3 * (o1 - o3 + double(ds1.waveforms.Data{k,j}/s1)));
						if isempty(bytes)
							ds.waveforms.Data{k,j} = [];
						else
							if min(bytes)<0 || max(bytes)>255
								error('Eek!  An unexpected error occurred in combining waveform data. Byte data went out of range');
							end
							ds.waveforms.Data{k,j} = char(bytes);
						end
					end

                    for k =1:size(ds2.waveforms.Data,1)
                        bytes = floor(s3 * (o2 - o3 + double(ds2.waveforms.Data{k,j}/s2)));
                        if isempty(bytes)
                            ds.waveforms.Data{arrayEndPoint + k,j} = [];
                        else
                            if min(bytes)<0 || max(bytes)>255
                                error('Eek!  An unexpected error occurred in combining waveform data. Byte data went out of range');
                            end
                            ds.waveforms.Data{arrayEndPoint + k,j} = char(bytes);
                        end
                    end
                end
        end
    end
elseif isfield(ds1.waveforms, 'NumOfChannels') || isfield(ds2.waveforms, 'NumOfChannels')
    error('One of the datasets has waveform data while the other does not.  Cannot merge')
else
 ds.waveforms = [];
end

%analog data
if length(ds1.analog.NumOfChannels) ~= length(ds2.analog.NumOfChannels)
    error('The two datasets have different numbers of analog channels.  Cannot merge')
end

if length(ds1.analog.SampleInterval) ~= length(ds2.analog.SampleInterval)
    error('The two datasets have different analog sample intervals.  Cannot merge')
end

ds.analog.NumOfChannels = ds1.analog.NumOfChannels;
ds.analog.SampleInterval = ds1.analog.SampleInterval;


ds.analog.Segments.SampleCounts = [ds1.analog.Segments.SampleCounts ds2.analog.Segments.SampleCounts];
ds.analog.Segments.StartTimes = [ds1.analog.Segments.StartTimes ds2.analog.Segments.StartTimes + startTimeForSecondDataSet];
ds.analog.Segments.Data = [ds1.analog.Segments.Data ds2.analog.Segments.Data];

%environment
envFieldNames = fieldnames(ds1.environment.Conversion);

for i=1:length(envFieldNames)
    field = envFieldNames{i};

    switch field
        case 'units'
            ds.environment.Conversion.units = ds1.environment.Conversion.units;
        otherwise
            if ds1.environment.Conversion.(field) ~= ds2.environment.Conversion.(field)
                error(sprintf('The two datasets have different values for the conversion property %s.  Cannot merge', field))
            end
            ds.environment.Conversion.(field) = ds1.environment.Conversion.(field);
    end
end
