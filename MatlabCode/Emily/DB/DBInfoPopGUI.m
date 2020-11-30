%used by LMTFOnline2 to select LMTF files before each run
function [paradigmName, RFX, RFY, subjectID, notes] = DBInfoPopGUI()
%% initializing
conn = database('Nex_Paradigm_Sort','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
f = figure('Name', 'Requesting Files');
set(f, 'Position', [600, 400, 400, 300]);
subjLMTF = fetch(conn, 'SELECT DISTINCT subjID FROM LMTF');
subjIso = fetch(conn, 'SELECT DISTINCT subjID FROM IsoSamp_old');
cancel = uicontrol(f, 'String', 'Cancel', 'Callback', 'close(gcf)', 'Position', [340, 230, 50, 40]);
standard_inputs = uicontrol(f, 'String', 'Standard Inputs', 'Callback', @std_ipts, 'Position', [300, 275, 80, 25], 'visible', 'off');
si = false;
OK = uicontrol(f, 'String', 'OK', 'Callback', 'uiresume(gcbf)', 'Position', [280, 230, 50, 40], 'visible', 'off');
pTextHandle = uicontrol(f, 'Style', 'Text', 'String', 'Paradigm Name: ', 'HorizontalAlignment', 'left', 'Position', [20 245 100 40]);
pList = {'LMTF', 'IsoSamp'};
pListboxHandle = uicontrol(f, 'Style', 'popup','String', pList, 'Max', 1, 'Position', [120, 250, 150, 40], 'callback', @clicktochoose);
guidata(pListboxHandle, 'value');
sTextHandle = uicontrol(f, 'Style', 'Text', 'String', 'Subject ID: ', 'HorizontalAlignment', 'left','Position', [20 205 100 40]);
%turn off "popupmenu won't be rendered until parameters are valid" warning
%w = warning('query', 'last');
%[warningID] = w.identifier;
%warning('off', warningID);
sListboxHandle = uicontrol(f, 'Style', 'popup', 'String', 'Select a paradigm','HorizontalAlignment', 'left','Position', [120, 210, 150, 40]);
    function clicktochoose(hcbo, eventstruct)
        PID = get(pListboxHandle, 'value');
        set(OK, 'visible', 'on');
        if PID == 1 %eventually make this better, honestly
            set(sListboxHandle, 'String', subjLMTF);
            set(standard_inputs, 'visible', 'on');
            sList = subjLMTF;
            guidata(pListboxHandle, 'value');
        else
            set(sListboxHandle, 'String', subjIso);
            sList = subjIso;
            guidata(pListboxHandle, 'value');
        end
    end

uiwait(gcf);
if ~ishandle(f)
    close(gcf);
    if si
        paradigmName = 'std_inputs';
    else
        paradigmName = 'cancel';
    end
    RFX = 'cancel';
    RFY = 'cancel';
    subjectID = 'cancel';
    notes = 'cancel';
    return;
else
    changeSubj = uicontrol(f, 'String', 'Change Subject', 'Callback', @swapSubj, 'Position', [160, 30, 90, 40]);
    paradigmName = pList{get(pListboxHandle, 'value')};
    subjectID = sList{get(sListboxHandle, 'value')};
    set(OK, 'visible', 'off');
    set(cancel, 'Position', [250, 30, 140, 40]);
    cancelLabel = '<HTML><center><FONT color="black"><i>Start New Experiment</i></Font>';
    set(cancel, 'String', cancelLabel);
    finalButton = uicontrol(f, 'callback', 'uiresume(gcbf)', 'HorizontalAlignment', 'center', 'Position', [10, 30, 150, 40], 'enable', 'off');
    finalLabel= '<HTML><center><FONT color="#989898"><i>Request Filenames from Database</i></Font>';
    set(finalButton, 'string', finalLabel);
    %% Section 2 RF and Dates
    rfQuery = sprintf('SELECT DISTINCT rfX, rfY FROM %s WHERE subjID = ''%s'' ORDER BY %s.recDate DESC', paradigmName, subjectID, paradigmName);
    fetchRF = flipud(fetch(conn, rfQuery));
    for i = 1:size(fetchRF,1)
        RFXY{i} = sprintf('(%d, %d)', fetchRF{i,1}, fetchRF{i,2});
    end
    RFTextHandle = uicontrol(f, 'Style', 'Text', 'String', 'Receptive Field Location: ', 'HorizontalAlignment', 'left', 'Position', [20 150 100 60]);
    RFpopup = uicontrol(f, 'Style', 'popup', 'Max', 1, 'Position', [120, 170, 150, 40], 'callback', @clickForNotes);
    notesTextHandle = uicontrol(f, 'Style', 'Text', 'String', 'Select notes subset: ', 'HorizontalAlignment', 'left', 'Position', [20 100 100 60]);
    notesListBox = uicontrol(f, 'Style', 'listbox', 'Max', 4, 'HorizontalAlignment', 'center', 'Position', [120, 80, 150, 80]);
    set(RFpopup, 'String', RFXY);
    %% Section 3 collecting final data to send off
    uiwait(gcf);
    if ~ishandle(f)
        close(gcf);
        if si
            paradigmName = 'std_inputs';
        else
            paradigmName = 'cancel';
        end
        RFX = 'cancel';
        RFY = 'cancel';
        subjectID = 'cancel';
        notes = 'cancel';
        return;
    else
        try
            rfIDX = get(RFpopup, 'value');
            RFX = fetchRF{rfIDX,1}; RFY = fetchRF{rfIDX,2};
            notesIDX = get(notesListBox, 'value');
            notesOptions = get(notesListBox, 'String');
            notes = notesOptions(notesIDX);
        catch
            RFX = fetchRF{1, 1}; RFY = fetchRF{1,2};
            notes = '*';
        end
        close(f); close(conn);
    end
end
    function clickForNotes(hcbo, eventStruct)
        rfIDX = get(RFpopup, 'value');
        rfVal = RFXY{rfIDX};
        rfSplit = strsplit(rfVal, {'(', ',', ')'});
        RFX = rfSplit{2}; RFY = rfSplit{3};
        notesQuery = sprintf('SELECT DISTINCT notes FROM %s WHERE subjID = ''%s'' AND rfX = %s AND rfY = %s', paradigmName, subjectID, RFX, RFY);
        fetchNotes = flipud(fetch(conn, notesQuery));
        set(notesListBox, 'Value', 1);
        set(notesListBox, 'String', fetchNotes);
        set(finalButton, 'enable', 'on');
        finalLabel2= '<HTML><center><FONT color="black"><i>Request Filenames from Database</i></Font>';
        set(finalButton, 'String', finalLabel2);
        guidata(notesListBox, 'String');
    end

    function swapSubj(hcbo, eventStruct)
        subjIDX = get(sListboxHandle, 'value');
        subjectID = sList{subjIDX};
        rfQuery = sprintf('SELECT DISTINCT rfX, rfY FROM %s WHERE subjID = ''%s''', paradigmName, subjectID);
        fetchRF = fetch(conn, rfQuery);
        for k = 1:size(fetchRF,1)
            new_RFXY{k} = sprintf('(%d, %d)', fetchRF{k,1}, fetchRF{k,2});
        end
        RFXY = new_RFXY;
        set(RFpopup, 'String', RFXY);
    end
    function std_ipts(hcbo, eventstruct)
        si = true;
        close(gcf);
    end
end