%Written for Greg a while back for his IsoSamp and WhiteNoise uploading
%needs, now replaced by AddDataToSQLByMultiCell.
function addDataToSQLByCell()
done = 0;
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
nf='';
while ~done
    another = 0;
    cellNums = fetch(conn, 'SELECT MAX(neuron) from IsoSamp_LGN UNION ALL SELECT MAX(neuron) from WhiteNoiseLGN_forIS');
    cellNums(cellfun(@(x) any(isnan(x)),cellNums)) = [];
    if isempty(cellNums)
        cell = 1;
    else
        cellNums = cell2mat(cellNums);
        cell = max(cellNums)+1;
    end
    %get cell file info from user
    data = dialogGen(1, nf);
    f_base = data{1};
    if size(f_base,2)>9
        f_base = f_base(1:8);
    end
    if size(f_base,2)==7
        f_base = f_base + '0';
    end
    f_subspec = data{2};
    f_subspec = strsplit(f_subspec, ',');
    spikeCode = data{3};
    if isempty(spikeCode)
        spikeCode = 0;
    elseif spikeCode == '1'
        spikeCode = 'sig001a';
    elseif spikeCode == '2'
        spikeCode == 'sig001b';
    end
    cellClass = data{4};
    if isempty(cellClass)
        cellClass = 'n/a';
    end
    splitNames = {};
    for i = 1:length(f_subspec)
        %check if any files are duplicates/already in table
        if size(f_subspec{i}, 2) < 2
            f_subspec{i} = strcat('0', f_subspec{i});
        end
        file = strcat(f_base, f_subspec{i}, '.nex');
        dupQuery_IS = sprintf('SELECT * FROM IsoSamp_LGN WHERE fileID = ''%s'' AND spikeCode = ''%s'';', file, spikeCode);
        dupCheck_IS = fetch(conn, dupQuery_IS);
        if ~isempty(dupCheck_IS)
            continue;
        else
            dupQuery_WN = sprintf('SELECT * FROM WhiteNoiseLGN_forIS WHERE fileID = ''%s'' AND spikeCode = ''%s'';', file, spikeCode);
            dupCheck_WN = fetch(conn, dupQuery_WN);
            if ~isempty(dupCheck_WN)
                continue;
            else
                if any(strcmp(splitNames,file))
                    continue;
                else
                    splitNames{i} = file;
                end
            end
        end
    end
    for i = 1:length(splitNames)
        PID = 0;
        rfX = 0;
        rfY = 0;
        tempdate = 0;
        subjectID = '';
        nf = splitNames{i};
        if isempty(nf)
            continue;
        else
            if ~strcmp(nf(end-3:end), '.nex')
                nf = [nf, '.nex'];
            end
            %paradigmID
            fp = findfile(nf);
            try
                PID = getparadigmID(fp);
            catch
                PID = warndlg(['Error entering file ', nf , ': PID'], 'Oops!');
            end
            if (size(PID,1)>1 || isnan(PID))
                PID = warndlg(['Error entering file ', nf , ': PID'], 'Oops!');
            else
                try
                    stro = nex2stro(fp);
                catch
                    PID = warndlg(['Error entering file ', nf , ': stro'], 'Oops!');
                end
                if (~isfield(stro,'sum'))
                    PID = warndlg(['Error entering file ', nf , ': stro'], 'Oops!');
                else
                    [tempRFX, tempRFY, bf] = dbRFSorter(PID, stro);
                    if bf
                        PID = warndlg(['Error entering file ', nf , ': rfx,y'], 'Oops!');
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
                        PID = warndlg(['Error entering file ', nf , ': date/subjID'], 'Oops!');
                    end
                end
            end
            data = dialogGen(2, nf);
            notes = data{1};
            quality = str2num(data{2});
            if quality > 1 || quality < 0
                quality = 0;
            end
            if PID==107
                tableName = 'IsoSamp_LGN';
            else
                tableName = 'WhiteNoiseLGN_forIS';
            end
            SQLInsertQuery = sprintf( 'INSERT INTO %s VALUES (NULL, ''%s'', ''%s'', ''%s'', %s, %s, %d, ''%s'', ''%s'', ''%s'', %d, %d);', tableName, nf, recordingDate, subjectID, ...
                rfX, rfY, quality, notes, spikeCode, cellClass, cell, PID);
            addToTable = exec(conn, SQLInsertQuery);
            if addToTable.Message
                warning = warndlg(addToTable.Message,'Error with SQL query');
                uiwait(warning);
                keyboard;
            else
                h = msgbox({['Inserted ', nf, ' to ', tableName, ' for subject ', subjectID]; ['RFs: ', rfX, ',', rfY]; ['Recording Date: ', recordingDate];['Notes: ', notes];
                    ['Spike Number: ', spikeCode]; ['Cell Class: ', int2str(cellClass)]; ['Quality: ', int2str(quality)]}, 'Successfully Inserted');
                uiwait(h);
            end
        end
    end
    another = questdlg('Do you have another cell?', 'Cell 2?', 'Yes', 'No', 'No');
    if strcmp(another, 'No')
        done = 1;
    end
end
close(conn);
end

function data = dialogGen(option_num, nf)
switch(option_num)
    case 1
        prompt = {'Filename Root (e.g. GDDMMYY00 or GDDMMYY or GDDMMYY001):', 'File Sub-specifiers (e.g. 1,2,5)', 'Enter cell''s spike number:', 'Enter cell classification:'};
        dlg_title = 'First Cell';
        size = [1 50; 1 50; 1 50; 1 50];
    case 2
        dlg_title = nf;
        prompt = {'Notes'; 'Quality (1 or 0)'};
        size = [1 40; 1 20];
    otherwise
        %something went wrong
        keyboard;
end
data = inputdlg(prompt,dlg_title, size);
end