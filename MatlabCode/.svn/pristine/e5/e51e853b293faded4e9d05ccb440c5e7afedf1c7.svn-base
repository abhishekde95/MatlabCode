% Code for pulling out the files where I haven't stored the waveforms
% Author - Abhishek De, 10/19
close all; clearvars;

spikename_options = ['sig001a'; 'sig001b'];
[filename_Lum, spikeIdx_Lum] = fnamesFromTxt2('Lum.txt');
[filename_ColorOpponent, spikeIdx_ColorOpponent] = fnamesFromTxt2('ColorOpponent.txt');

% Loading the Neurothresh files
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
tmp_filename = fetch(conn,'SELECT filename FROM WNthresh');
NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
close(conn);
tmp_filename = tmp_filename(strcmp('subunit',NTmode));
spikeidx_NT = spikeidx_NT(strcmp('subunit',NTmode));
filename_NT = cell(size(tmp_filename));
for kk = 1:size(tmp_filename)
    filename_NT(kk) = {tmp_filename(kk)};
end

% Loading the just the subunit files 
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
tmp_filename = fetch(conn,'SELECT filename FROM WNSubunit');
subunit_mode = fetch(conn,'SELECT mode FROM WNSubunit');
spikeidx_subunit = cell2mat(fetch(conn,'SELECT spikeidx FROM WNSubunit'));
close(conn);
tmp_filename = tmp_filename(strcmp('STA',subunit_mode));
spikeidx_subunit = spikeidx_subunit(strcmp('STA',subunit_mode));
filename_subunit = cell(size(tmp_filename));
for kk = 1:size(tmp_filename)
    filename_subunit(kk) = {tmp_filename(kk)};
end


Input_List = [filename_Lum; filename_ColorOpponent; filename_NT; filename_subunit];
spikeIdx = [spikeIdx_Lum; spikeIdx_ColorOpponent; spikeidx_NT; spikeidx_subunit];


files_not_working = [];
files_without_waveform = [];
count = 1;
load Output_List_waveform.mat
N = size(Output_List,1);

for ii = 1:N
    flag = 0;
    stro = {};
    for jj = 1:size(Input_List{ii},2)
        try 
            tmpstro = nex2stro(findfile(char(Input_List{ii}(jj))));
        catch 
            files_not_working = [files_not_working; Input_List{ii}];
            flag = 1;
            break;
        end
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
    end
    if flag 
        continue;
    end
    
    if ~any(strcmp(stro.sum.rasterCells,'sig001a_wf'))
       files_without_waveform = [files_without_waveform; Input_List{ii}]; 
    end
end

