function updateDBSingleFile(newest, origPath, newPath)
global gl
%pull subjID, RFX, RFY, recDate and check if in DB
% example call:
% updateDBSingleFile('G010117001.nex','D:\PlexonData','N:\NexFiles\Emily')
conn = database('Nex_Paradigm_Sort','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
checkDBforFileQuery = sprintf('SELECT fileID FROM LMTF WHERE fileID= ''%s'';', newest);
isFileInDB = fetch(conn, checkDBforFileQuery);
if ~isempty(isFileInDB)
    warnDB(4, origPath, newPath, newest); %FILE ALREADY EXISTS IN DB
else
    try
        stro = nex2stro(fullfile(origPath, newest));
    catch
        warnDB(5, origPath, newPath, newest); %HEADER BROKEN
        return;
    end
    %pull relevant information from file
%    date = num2str(sscanf(newest, '%*c%6d*')); %pulls 6 digits from filename
% GDLH 1/1/17 above libe has problems if first character is a zero
    date = newest(end-12:end-7);
    subjID = sprintf('''%s''', sscanf(newest, '%c*')); %pulls just whatever is before the recDate
    recDate = sprintf('''20%s-%s-%s''', date(5:6), date(1:2), date(3:4));
    rfX = abs(stro.sum.exptParams.stim_x);
    rfY = stro.sum.exptParams.stim_y;
    eyewin_w = stro.sum.exptParams.eyewin_w;
    eyewin_h = stro.sum.exptParams.eyewin_h;
    monID = 'ProPixx';
    if strcmp(subjID, ['A', 'U', 'M', 'P'])
        nhp = 1;
    else
        nhp = 0;
    end
    %check if notes in DB (using a quality of 1)
    checkDBforNotes = sprintf('SELECT DISTINCT notes FROM LMTF WHERE rfX = %d AND rfY = %d AND subjID = %s ORDER BY LMTF.recDate DESC', rfX, rfY, subjID);
    totalNotes = fetch(conn, checkDBforNotes);
    %open figure 1, with filename as an editable
    f = figure('Name', 'Update DB with new file', 'Position', [500, 400, 530, 270]);
    finalButton = uicontrol(f, 'callback', 'uiresume(gcbf)', 'HorizontalAlignment', 'center', 'Position', [180, 30, 100, 40]);
    finalLabel= '<HTML><center><FONT color="black"><i>Submit File to LMTF Database</i></Font>';
    set(finalButton, 'string', finalLabel);
    cancelButton = uicontrol(f, 'callback', 'delete(gcf)', 'String', 'Cancel', 'HorizontalAlignment', 'center', 'Position', [280, 30, 50, 40]);
    
    fileNameHandle = uicontrol(f, 'Style', 'Text', 'String', 'Filename: ', 'HorizontalAlignment', 'left', 'position', [130 200 80 40]);
    fileNameTextbox = uicontrol(f, 'Style', 'edit', 'String', newest, 'HorizontalAlignment', 'center', 'position', [180 210 100 40]);
    %option to move file to another directory, with a generic dir path as editable
    moveFileHandle = uicontrol(f, 'Style', 'Text', 'String', 'Move file to new directory?', 'HorizontalAlignment', 'left', 'position', [40 150 140 40]);
    moveFileTextbox = uicontrol(f, 'Style', 'edit', 'String', newPath, 'HorizontalAlignment', 'center', 'position', [180 160 150 40]);
    moveFileConfirm = uicontrol(f, 'String', 'Move', 'HorizontalAlignment', 'left', 'position', [350 160 80 40], 'callback', @fileMover);
    %toggle button to change quality to 0 or back to 1 - CHANGE request to quality = 0 if this is toggled
    qualityToggle = uicontrol(f, 'Style', 'ToggleButton', 'String', 'Quality = 1', 'HorizontalAlignment', 'left', 'position', [20 100 100 40], 'callback', @qualCheck);
    notesText = uicontrol(f, 'Style', 'Text', 'String', 'Notes: ', 'HorizontalAlignment', 'left', 'position', [140 90 100 40]);
    notesListBox = uicontrol(f, 'Style', 'listbox', 'String', totalNotes, 'HorizontalAlignment', 'center', 'position', [180 100 150 40]);
    notesEditBox = uicontrol(f, 'Style', 'edit', 'String', 'Enter Notes Here', 'HorizontalAlignment', 'center', 'position', [180 100 200 40], 'visible', 'off');
    notesListOverride = uicontrol(f, 'String', 'Add New Note', 'HorizontalAlignment', 'left', 'position', [400 100 100 40], 'callback', @noteEditOverride);
    %send info to DB
    uiwait(gcf);
    if ~ishandle(f)
        close(gcf);
        return;
    else
        opposite = get(qualityToggle, 'value');
        quality = ~opposite;
        ifEdit = get(notesEditBox, 'String');
        ifList = get(notesListBox, 'value');
        if strcmp(ifEdit, 'Enter Notes Here') && isempty(totalNotes)
            notes = 'NULL';
        elseif strcmp(ifEdit, 'Enter Notes Here') && ~isempty(totalNotes)
            notes = totalNotes{ifList};
            if strcmp(notes, 'null')
                notes = 'NULL';
            else
                notes = sprintf('''%s''', notes);
            end
        else
            notes = sprintf('''%s''', ifEdit);
        end
        fileName = sprintf('''%s''', get(fileNameTextbox, 'string'));
        sqlInsertQuery = sprintf('INSERT INTO LMTF VALUES(NULL, %s, %s, %s, %d, %d, %d, %s, NULL, NULL, NULL, 157, %d, %d, %d, ''%s'', %d)', fileName, recDate, subjID, rfX, rfY, quality, notes, nhp, eyewin_w, eyewin_h, monID, gl.std_inputs);
        sqlInsert = exec(conn, sqlInsertQuery);
        if sqlInsert.Message
            errordlg(sqlInsert.Message, 'SQL error');
            keyboard;
            close(gcf);
        else
            close(gcf);
        end
    end
end
%callback function definitions
    function qualCheck(hcbo, eventStruct)
        q = get(qualityToggle, 'value');
        if q==1
            guidata(qualityToggle, 'value');
            set(qualityToggle, 'String', 'Quality = 0');
        else
            guidata(qualityToggle, 'value');
            set(qualityToggle, 'String', 'Quality = 1');
        end
    end

    function fileMover(hcbo, eventstruct)
        file = get(fileNameTextbox, 'String');
        NP = get(moveFileTextbox, 'String');
        try
            movefile(file, NP);
        catch
            warndlg('Cannot move the file onto itself. Change new directory and try again.');
        end
    end
    function noteListOverride(hcbo, eventStruct)
        set(notesEditBox, 'visible', 'off');
        set(notesListOverride, 'String', 'Add New Note');
        set(notesListOverride, 'callback', @noteEditOverride);
        set(notesListBox, 'visible', 'on');
    end
    function noteEditOverride(hcbo, eventStruct)
        set(notesEditBox, 'visible', 'on');
        set(notesListBox, 'visible', 'off');
        set(notesListOverride, 'String', 'Select Note');
        set(notesListOverride, 'callback', @noteListOverride);
    end
close(conn);
end