%also deprecated practice.
function [paradigmName, RFX, RFY, subjectID, recordingDates, quality, notes] = secondGUI()
%% initializing
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
f = figure('Name', 'Requesting Files');
pTextHandle = uicontrol(f, 'Style', 'Text', 'String', 'Paradigm Name: ', 'HorizontalAlignment', 'left', 'Position', [20 360 100 40]);
pList = {'LMTF', 'IsoSamp', 'Whatever'};
pListboxHandle = uicontrol(f, 'Style', 'listbox','String', pList, 'Max', 1, 'Position', [120, 380, 150, 40]);
OK = uicontrol(f, 'String', 'OK', 'Callback', 'uiresume(gcbf)', 'Position', [280, 380, 50, 40]);
uiwait(gcf);
pIDX = get(pListboxHandle, 'value');
paradigmName = pList{pIDX};
sList = fetch(conn, sprintf('SELECT DISTINCT subjID FROM %s', paradigmName));
RFXlist = fetch(conn, sprintf('SELECT DISTINCT rfX FROM %s', paradigmName));
RFYlist = fetch(conn, sprintf('SELECT DISTINCT rfY FROM %s', paradigmName));
finalButton = uicontrol(f, 'callback', 'uiresume(gcbf)', 'HorizontalAlignment', 'center', 'Position', [200, 40, 150, 40]);
finalLabel= '<HTML><center><FONT color="black"><i>Request Filenames from Database</i></Font>';
set(finalButton, 'string', finalLabel);
%% Section 2: Data collection
%rf
RFTextHandle1 = uicontrol(f, 'Style', 'Text', 'String', 'Receptive Fields: ','HorizontalAlignment', 'left', 'Position', [20 330 100 40]);
RFTextHandle2 = uicontrol(f, 'Style', 'Text', 'String', 'X - ', 'HorizontalAlignment', 'left','Position', [120 330 50 40]);
RFListboxHandleX = uicontrol(f, 'Style', 'listbox', 'String', RFXlist, 'Max', 1, 'Position', [160 330 100 40]);
RFTextHandle3 = uicontrol(f, 'Style', 'Text', 'String', 'Y - ', 'HorizontalAlignment', 'left','Position', [270 330 50 40]);
RFListboxHandleY = uicontrol(f, 'Style', 'listbox', 'String', RFYlist, 'Max', 1,'Position', [310 330 100 40]);
%subjectID
sTextHandle1 = uicontrol(f, 'Style', 'Text', 'String', 'Subject ID: ', 'HorizontalAlignment', 'left','Position', [20 280 100 40]);
sListboxHandle = uicontrol(f, 'Style', 'listbox','String', sList, 'Max', 1, 'HorizontalAlignment', 'left','Position', [120, 280, 150, 40]);
%recordingDates
dateTextHandle = uicontrol(f, 'Style', 'Text', 'String', 'Recording Dates: ', 'Horizontalalignment', 'left', 'Position', [20 230 100 40]);
dateEditBoxHandle = uicontrol(f, 'Style', 'Text','String', '*', 'HorizontalAlignment', 'left', 'BackgroundColor', 'w', 'Position', [120 230 200 40]);
calendarButtonHandle = uicontrol(f, 'String', 'Select dates', 'callback', @pushbutton_cb,'Horizontalalignment', 'center', 'Position', [330 230 100 40]);
    function pushbutton_cb(hcbo, eventStruct)
        uicalendar('Weekend', [0 0 0 0 0 1 1], 'SelectionType', 0, 'DestinationUI', dateEditBoxHandle);
    end
%Other Options (notes, quality)
qList = {'*', '0', '1'};
qHandle = uicontrol(f, 'Style','text', 'String', 'Quality: ', 'HorizontalAlignment', 'left','Position', [20 170 100 40], 'visible', 'off');
qListHandle = uicontrol(f, 'Style', 'listbox', 'String', qList, 'Max', 1, 'HorizontalAlignment', 'left','Position', [80 180 80 40], 'visible', 'off');
nHandle = uicontrol(f, 'Style','text', 'String', 'Notes: ', 'HorizontalAlignment', 'left','Position', [180 170 100 40], 'visible', 'off');
nEditHandle = uicontrol(f, 'Style', 'Edit', 'String', '*', 'HorizontalAlignment', 'left', 'BackgroundColor', 'w', 'Position', [220 180 200 40], 'visible', 'off');
oButtonHandle = uicontrol(f, 'String', 'Other Options', 'Callback', @other,'HorizontalAlignment', 'left', 'Position', [20 180 100 40]);
    function other(hcbo, eventStruct)
        set(oButtonHandle, 'visible', 'off');
        set(qHandle, 'visible', 'on');
        set(qListHandle, 'visible', 'on');
        set(nHandle, 'visible', 'on');
        set(nEditHandle, 'visible', 'on');
    end
drawnow;
%% final section
uiwait(gcf);
RFIDX_x = get(RFListboxHandleX, 'value');
RFX = int2str(RFXlist{RFIDX_x});
RFIDX_y = get(RFListboxHandleY, 'value');
RFY = int2str(RFYlist{RFIDX_y});
sIDX = get(sListboxHandle, 'value');
subjectID = strcat('''', sList{sIDX}, '''');
recordingDates = get(dateEditBoxHandle, 'string');
notes = get(nEditHandle, 'String');
qIDX = get(qListHandle, 'value');
quality = qList{qIDX};
close(f); close(conn);
end