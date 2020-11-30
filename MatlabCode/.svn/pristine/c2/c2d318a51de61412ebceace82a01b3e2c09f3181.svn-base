%Necessary! Part of LMTFOnline2 uploading file to database after each lmtf
%file is collected. Called by fileGrabber
function warnDB(type, origPath, newPath, filename)
wd = dialog('Name', 'Warning', 'Position', [500, 600, 500, 200]);
choice1 = uicontrol('Parent', wd, 'Style', 'pushbutton', 'String', 'Change Original Directory', 'callback', @changeOD, 'Position',[20 80 150 40]);
choice2 = uicontrol('Parent', wd, 'Style', 'pushbutton', 'String', 'Check For Recent Files Again', 'callback', @checkAgain, 'Position',[180 80 150 40]);
choice3 = uicontrol('Parent', wd, 'Style', 'pushbutton', 'String', 'Cancel', 'callback', 'close(gcf)', 'Position',[380 130 80 40]);
choice4 = uicontrol('Parent', wd, 'Style', 'pushbutton', 'String', 'Pull Second Most Recent File', 'callback', @lastFile, 'Position',[340 80 150 40]);
newOD = uicontrol('Parent', wd, 'Style', 'edit', 'String', origPath, 'visible', 'off', 'Position',[20 40 210 40]);
ok = uicontrol('Parent', wd, 'Style', 'pushbutton', 'String', 'Done', 'visible', 'off', 'callback', @continueChange, 'Position',[240 40 50 40]);
switch type
    case 3
        warningString = sprintf('Recent File ''%s%s'' is still in ''.plx'' format. Convert and pull again.', origPath, filename);
    case 1
        warningString = sprintf('No recent files in ''.nex'' or ''.plx'' format found in original directory ''%s''.', origPath);
    case 2
        warningString = sprintf('No files found in original directory ''%s''!', origPath);
        set(choice4, 'visible', 'off');
    case 4
        warningString = sprintf('File ''%s%s'' already exists in DB.', origPath, filename);
    case 5
        warningString = sprintf('File ''%s%s'' is broken and cannot be added to the DB.', origPath, filename);
    otherwise
        warningString = 'Something really broke!';
end
warning = uicontrol('Parent', wd, 'Style', 'text', 'String', warningString, 'Position',[30 120 350 40]);
uiwait(wd);
    function changeOD(hcbo, eventStruct)
        set(newOD, 'visible', 'on');
        set(ok, 'visible', 'on');
        set(choice2, 'visible', 'off');
    end
    function continueChange(hcbo, eventStruct)
        newOrigPath = get(newOD, 'String');
        close(gcf);
        fileGrabber(newOrigPath, newPath, 1);
    end
    function checkAgain(hcbo, eventStruct)
        close(gcf);
        fileGrabber(origPath, newPath, 1);
    end
    function lastFile(hcbo, eventStruct)
        close(gcf);
        fileGrabber(origPath, newPath, 0);
    end
end