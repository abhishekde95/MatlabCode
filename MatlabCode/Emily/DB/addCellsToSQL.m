%deprecated version of adddatatosqlbymulticell
function addCellsToSQL()
done = 0;
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
while ~done
    another = 0;
    cellnums = fetch(conn, 'SELECT MAX(neuron) from IsoSamp_LGN UNION ALL SELECT MAX(neuron) from WhiteNoiseLGN_forIS');
    %check for "null" values
    if cellnums{1} == cellnums{2}
        cell = cellnums{1} + 1;
    else
        cellnums = cell2mat(cellnums);
        cell = max(cellnums)+1;
    end
    data = dialogGen(1);
    newFilenames = data{1}; %check that this is what you want
    spikeNum = data{2};
    cellClass = data{3};
    %sanitize all manual inputs
    for i = 1:length(newFilenames)
        nf = newFilenames{i};
        try
            PID = getparadigmID(nf);
        catch
            PID = dialogGen(2);
            if ~isnumeric(PID)
                %re-prompt infinitely until correct
            end
        end
        if (size(PID,1)>1 || isnan(PID))
            PIDname = dialogGen(2);
            try
                PID = getParadigmID(PIDname);
            catch
                dialogGen(8);
            end
        else
            try
                stro = nex2stro(nf);
            catch
                data = dialogGen(3); %FIX
            end
            if (~isfield(stro,'sum'))
                data = dialogGen(4); %FIX
            else
                [tempRFX, tempRFY, bf] = dbRFSorter(PID, stro);
                if bf
                    dialogGen(5);
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
                    dialogGen(7);
                end
            end
        end
        data = dialogGen(9);
        notes = data{1};
        quality = data{2};
        if PID==107
            tableName = 'IsoSamp_LGN';
        else
            tableName = 'WhiteNoiseLGN_forIS';
        end
        SQLInsertQuery = ['INSERT INTO ' tableName ' VALUES(NULL, ' nf,''',''', recordingDate, ''',''',subjectID,''',', rfX,', ', rfY, quality, notes, spikeNum, cellClass, cell, int2str(PID), ');'];
        addToTable = exec(conn, SQLInsertQuery);
        if addToTable.messsage()
            warning = warndlg({'Error with SQL query', addToTable.message()},'Oops!');
            %what do
        end
    end
    another = questdlg('Do you have another cell?', 'Cell 2?', 'Yes', 'No', 'No');
    if strcmp(another, 'No')
        done = 1;
    end
end
close(conn);
end

function data = dialogGen(option_num)
switch(option_num)
    case 1
        prompt = {'Enter semicolon separated filenames'; 'Enter cell''s spike number'; 'Enter cell classification'};
        dlg_title = {'Enter information for first cell'};
        size = [1 50; 1 10; 1 10];
    case 2
        prompt = {'Enter Paradigm Name Manually (e.g. WhiteNoise or IsoSamp)'};
        dlg_title = {'header issue - paradigm name could not be extracted from header '+ nf};
        size = [1 20];
    case 3
        prompt = {'Enter Paradigm Name Manually (e.g. WhiteNoise or IsoSamp)'};
        dlg_title = {'header issue - paradigm ID was stored improperly in header '+ nf};
        size = [1 20];
    case 4
        dlg_title = {'stro issue - file could not be opened using nex2stro '+ nf};
        prompt = {'rfX'; 'rfy', 'rec date (YYYY-MM-DD)', 'subj id'};
        size = [1 10; 1 10; 1 10; 1 10];
    case 5
        dlg_title = {'stro issue - sum field could not be found '+ nf};
        prompt = {'rfX'; 'rfy', 'rec date (YYYY-MM-DD)', 'subj id'};
        size = [1 10; 1 10; 1 10; 1 10];
    case 6
        dlg_title = {'receptive field issue - RF could not be determined for file '+nf};
        prompt = {'rfX'; 'rfY'};
        size = [1 10; 1 10];
    case 7
        dlg_title = {'recording date and/or subject ID could not be determined for file '+nf};
        prompt = {'rec date (YYYY-MM-DD)'; 'subj id'};
        size = [1 10; 1 10];
    case 8
        dlg_title = {'Error with Paradigm name entered.'};
        prompt = {'Enter paradigm name manually(e.g. WhiteNoise or IsoSamp)'};
        size = [1 20];
    case 9
        dlg_title = {'Enter Notes and Quality for file ' + nf};
        prompt = {'Notes'; 'Quality (1 or 0)'};
        size = [1 40; 1 2];
    otherwise
        %something went wrong
end
data = inputdlg(prompt,dlg_title, size);
end