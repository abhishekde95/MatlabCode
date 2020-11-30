function [fnames, spikenums, spikecodes, uniqueneuronIDs] = fnamesFromTxt(tablename, varargin)
% function [fnames, spikenums, <spikecodes>, <neuron_id>] = fnamesFromTxt(tablename, varargin)
%
% A replacement function for fnamesFromTxt2.m, this function will extract
% filenames and spike numbers from the SQL database and return them in the
% same format that fnamesFromTxt2.m does (a cell array of filenames and a
% vector of indicies into stro.ras.
%
%   EXAMPLE: [fnames, spikeIdx] = fnamesFromTxt('WhiteNoiseLGN_forIS');
%
% GDLH 11/28/17

sel_crit_str =[' WHERE '];
if nargin>1
    if ~rem(nargin,2)
        error('Input argument error: fnamesFromTxt(table name, <col1>, <val1>, <col2>, <val2>...)');
    end
    
    for i = 1:2:(nargin-1)
        whichcol = varargin{i};
        whichvals = varargin{i+1};
        nvals = length(whichvals);
        if nvals > 1
            sel_crit_str = [sel_crit_str,'( '];
        end
        for j = 1:nvals
            sel_crit_str = [sel_crit_str,whichcol,' = ''',whichvals{j},''' '];
            if j < nvals
                sel_crit_str = [sel_crit_str,'OR '];
            end
        end
        if nvals > 1
            sel_crit_str = [sel_crit_str,') '];
        end
        sel_crit_str = [sel_crit_str,'AND ']; % Always need to add this because of quality = 1
    end
end
sel_crit_str = [sel_crit_str,'quality = 1']

conn = database('Nex_Paradigm_Sort','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
str = sprintf('SELECT COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS WHERE TABLE_SCHEMA = ''Nex_Paradigm_Sort'' AND TABLE_NAME = ''%s''',tablename);
column_names = fetch(conn, str, 'DataReturnFormat','cellarray');
allfnames = fetch(conn, ['SELECT fileID FROM ',tablename, sel_crit_str],'DataReturnFormat','cellarray');

if any(strcmp(column_names,'neuron'))
    neuronIDnum = cell2mat(fetch(conn, ['SELECT neuron FROM ',tablename, sel_crit_str], 'DataReturnFormat','cellarray'));
else
    neuronIDnum = [];
end
if any(strcmp(column_names,'spikeCode'))
    allspikecodes = cell2mat(fetch(conn, ['SELECT spikeCode FROM ',tablename, sel_crit_str], 'DataReturnFormat','cellarray'));
else
    allspikecodes = [];
end
close(conn)

if all(isnan(neuronIDnum)) | isempty(neuronIDnum)
    neuronIDnum = [1:length(neuronIDnum)]';
end
uniqueneuronIDs = unique(neuronIDnum);

% Building "filenames" and "spikenums"
fnames = {};
spikenums = [];
if ~isempty(neuronIDnum)
    for neuronID = uniqueneuronIDs'
        L = neuronIDnum == neuronID;
        idx = size(fnames,1)+1;
        fnames = cat(1, fnames, {allfnames(L)'});
        if ~isempty(allspikecodes)
            spikecodes(idx,:) = unique(allspikecodes(L,:),'rows'); % This should error if >1 spike codes are present, which shouldn't happen
            spikenums(idx,:) = abs(spikecodes(idx,end))-96;
        end
    end
else
    fnames = allfnames;
    spikenums = nan(size(fnames));
end
    
