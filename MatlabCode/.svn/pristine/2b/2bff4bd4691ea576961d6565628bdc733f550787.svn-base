%simple function showing how to pull specific lmtf filenames from database
%and display contents. Connect to synology first.
function filenames = pullFilenamesFromDB(paradigm, subject, x, y)
conn = database('FullParadigmSortedDB','admin','adminVector123','Vendor','MySql','Server','localhost'); %connect to database
pullFilenamesRequest = ['Select fileID FROM ', paradigm, ' WHERE subjID = ''', subject, ''' AND (RFX = ', int2str(x), ' OR RFX = -', int2str(x), ...
        ') AND (RFY= ', int2str(y), ' OR RFY = -', int2str(y), ');'];
pullFilenamesCursor = exec(conn, pullFilenamesRequest);
fileFetch = fetch(pullFilenamesCursor);
filenames = fileFetch.Data';
close(conn);
end

%TO CONNECT TO SYNOLOGY:
%conn = database('testing','','','Vendor','MySql','Server','128.95.153.12')