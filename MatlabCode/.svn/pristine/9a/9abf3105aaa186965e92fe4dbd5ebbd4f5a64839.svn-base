%Written for Greg - takes local file list and organizes files by "cells"
%collected within a day. For each cell, separates them into isosamp or
%white noise files, and after requesting relevant information from the
%user, sends these files to their respective tables in the database.
%Additionally, labels each cell with a number so cells across tables can be
%connected
function addDataToSQLByMultiCell()
%test file: U111617003
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
prompt = {'Filename Root (e.g. GDDMMYY00 or GDDMMYY or GDDMMYY001):', 'File Sub-specifiers (e.g. 1,2,5)'};
dlg_title = 'Multi Cell File Info';
dlg_size = [1 50; 1 50;];
data = inputdlg(prompt,dlg_title, dlg_size);
f_base = data{1};
if size(f_base,2)>9
    f_base = f_base(1:8);
elseif size(f_base,2)==7
    f_base = f_base + '0';
end
f_subspec = data{2};
f_subspec = strsplit(f_subspec, ',');
splitNames = {};
for i = 1:length(f_subspec)
    curr_fs = f_subspec{i};
    curr_fs = strtrim(curr_fs);
    if size(curr_fs, 2) < 2
        curr_fs = strcat('0', curr_fs);
    end
    splitNames{i} = strcat(f_base, curr_fs, '.nex');
end
ras_id = {};
ras_label = 'ABCDEFGHIJKLMNOP';
for j = 1:length(splitNames)
    nf = splitNames{i};
    if ~strcmp(nf(end-3:end), '.nex')
        nf = [nf, '.nex'];
    end
    %paradigmID
    fp = findfile(nf);
    try
        stro = nex2stro(fp);
    catch
        warning = warndlg(['Error entering file ', nf , ': stro'], 'Oops!');
        uiwait(warning);
    end
    if (~isfield(stro,'sum'))
        warning = warndlg(['Error entering file ', nf , ': stro'], 'Oops!');
        uiwait(warning);
        keyboard;
    else
        ssr = stro.sum.rasterCells;
        curr_ras = [];
        for r = 1:size(ssr, 2)
            ras_idx = strncmp(ssr{r}, 'sig0', 4);
            if ras_idx
                curr_ras = [curr_ras; ssr(r)]; %storing stro, spike code
            end
        end
        %curr_ras = strcat(curr_ras, ras_label(j)); %labeling spike code with file specific identifier
        ras_id{j} = {stro, curr_ras, length(curr_ras), zeros(length(curr_ras),1), nf};
    end
end
num_cells = cell2mat(cellfun(@(x) x(3), ras_id));
cellsToAdd = 0;
%cell_comp_mtxs = cellfun(@(x) x(k,:), (cellfun(@(x) x(2), ras_id)), 'UniformOutput',false);
cell_comp_mtxs = cellfun(@(x) x(2), ras_id);
cell_list = [];
for k = 1:length(cell_comp_mtxs)
    cell_list = [cell_list; cell2mat(cell_comp_mtxs{k})];
end
while ~isempty(cell_list)
    cell_wid = cell_list(1,:);
    cell_nid = cell_wid(5:end);
    cellsToAdd = cellsToAdd + 1;
    for k = 1:length(cell_comp_mtxs)
        match_idx = ~cellfun(@isempty, strfind(cell_comp_mtxs{k}, cell_nid)); %find matching suffixes in each cell list, or matching cells
        cell_id = ras_id{k}{4};
        if cell_id(match_idx) == 0
            cell_id(match_idx) = cellsToAdd;
            ras_id{k}{4} = cell_id;
        end
    end
    cell_list(1,:) = [];
end
keyboard;
max_cell_nums = 0;
for l = 1:length(ras_id)
    local_max_cell_num = max(ras_id{l}{4});
    if local_max_cell_num > max_cell_nums
        max_cell_nums = local_max_cell_num;
    end
end
for n = 1:length(ras_id)
    curr_cell_info = ras_id{n};
    curr_stro = curr_cell_info{1};
    nf = curr_cell_info{5};
    PID = curr_stro.sum.paradigmID;
    if PID==107
        tableName = 'IsoSamp_LGN';
    else
        tableName = 'WhiteNoiseLGN_forIS';
    end
    [tempRFX, tempRFY, bf] = dbRFSorter(PID, stro);
    if bf
        warning = warndlg(['Error entering file ', nf , ': rfx,y'], 'Oops!');
        uiwait(warning);
    else
        rfX = tempRFX;
        rfY = tempRFY;
    end
    tempdate = nf(isstrprop(nf, 'digit'));
    if ge(length(tempdate), 6)
        tempdate = tempdate(1:6);
        subjectID = nf(1:strfind(nf, tempdate)-1);
        recordingDate = ['20',tempdate(end-1:end), '-', tempdate(1:2), '-', tempdate(3:4)];
    else
        warning = warndlg('Error entering file date/subjID', 'Oops!');
        uiwait(warning);
    end
    cellNums = fetch(conn, 'SELECT MAX(neuron) from IsoSamp_LGN UNION ALL SELECT MAX(neuron) from WhiteNoiseLGN_forIS');
    cellNums(cellfun(@(x) any(isnan(x)),cellNums)) = [];
    if isempty(cellNums)
        cell = 1;
    else
        cellNums = cell2mat(cellNums);
        cell = max(cellNums)+1;
    end
    for m = 1:max_cell_nums
        spikeContainer = curr_cell_info{4};
        spikeCode = spikeContainer{m};
        cell = cell+m;
        if spikeCode
            dlg_title = strcat(nf, m);
            prompt = {'Notes'; 'Quality (1 or 0)'; 'Cell Class'};
            dlgsize = [1 40; 1 20; 1 20];
            data = inputdlg(prompt,dlg_title, dlgsize);
            notes = data{1};
            quality = str2num(data{2});
            if quality > 1 || quality < 0
                quality = 0;
            end
            cellClass = data{3};
            SQLInsertQuery = sprintf( 'INSERT INTO %s VALUES (NULL, ''%s'', ''%s'', ''%s'', %s, %s, %d, ''%s'', %s, ''%s'', %d, %d);', tableName, nf, recordingDate, subjectID, ...
                rfX, rfY, quality, notes, spikeCode, cellClass, cell, PID);
            addToTable = exec(conn, SQLInsertQuery);
            if addToTable.Message
                warning = warndlg(addToTable.Message,'Error with SQL query');
                uiwait(warning);
                keyboard;
            else
                h = msgbox({['Inserted ', nf, ' to ', tableName, ' for subject ', subjectID]; ['RFs: ', rfX, ',', rfY]; ['Recording Date: ', recordingDate];['Notes: ', notes];
                    ['Spike Number: ', int2str(spikeCode)]; ['Cell Class: ', int2str(cellClass)]; ['Quality: ', int2str(quality)]}, 'Successfully Inserted');
                uiwait(h);
            end
        end
    end
end
close(conn);
end