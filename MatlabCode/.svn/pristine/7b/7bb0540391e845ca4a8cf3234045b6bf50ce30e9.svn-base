%called by DBPostPopGui - part of LMTFOnline 2 for adding newly collected
%LMTF file to database.
function fileGrabber(origPath, newPath, mostRecent)
cd(origPath); d = dir;
d = d(~[d.isdir]);
if ~isempty(d)
    [dummy, dx]=sort([d.datenum]);
    if mostRecent
        newest = d(dx(end)).name;
    else
        newest = d(dx(end-1)).name;
    end
else
    newest = [];
end
if ~isempty(newest)
    ifNex = strsplit(newest, '.');
    if strcmp(ifNex{2}, 'nex');
        updateDBSingleFile(newest, origPath, newPath);
    elseif strcmp(ifNex{2}, 'plx')
        warnDB(3, origPath, newPath, newest); %CONVERT FILE TO NEX
    else
        warnDB(1, origPath, newPath, newest); %NO NEX OR PLX FILES FOUND
    end
else
    warnDB(2, origPath, newPath, newest);%NO RECENT FILE IN DIRECTORY
end
end