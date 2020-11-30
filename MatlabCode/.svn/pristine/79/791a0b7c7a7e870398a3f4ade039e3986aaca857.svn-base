%Pretty useful for LMTF files. Scrapes entirety of nex files folder, collects all
%LMTF files, organizes data within each file, and sends to database. May
%not currently work as LMTF table in nex paradigm sort database has changed
%a lot since this was written, but can be easily updated.
startDir = 'N:\NexFiles';
brokenfiles = 0;
TablesCreated=0;
FilesEntered=0;
FilesSkipped=0;
filenames = scrapeFilesFromDir(startDir, '*.nex');
unique_filenames = unique(filenames); %FIXME
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12'); %connect to sorted files database
brokenFilesConn = database('BrokenFiles','','','Vendor','MySql','Server','128.95.153.12'); %connect to broken files database
disp('connected to Nex_Paradigm_Sort database and BrokenFiles database');
for i = 1:length(unique_filenames)
    filename = unique_filenames{i}; %FIX LATER
    try
        PID = getparadigmID(filename);
    catch
        BrokenFileAdder(brokenFilesConn, filename, 'No Paradigm ID Available', 'NULL');
        brokenfiles = brokenfiles+1;
        continue;
    end
    if (size(PID,1)>1 || isnan(PID)) %sometimes files have something like a header, but they're screwy and (in reality) broken
        BrokenFileAdder(brokenFilesConn, filename, 'Header Partially Corrupted', 'NULL');
        brokenfiles = brokenfiles+1;
        continue;
    else %continue trying to add the file
        try
            paradigmStringName = paradigmLibrary(PID);
        catch
            BrokenFileAdder(brokenFilesConn, filename, 'PID not in paradigmLibrary', int2str(PID));
            brokenfiles = brokenfiles+1;
            continue;
        end
        %sort out whether a table already exists in the database for that paradigmID. If not, create.
        [tableName, tc] = TableAdder(paradigmStringName, conn);
        TablesCreated = TablesCreated + tc;
        %check that the file isn't in DB already.
        sqlCheckQuery = ['SELECT * FROM ', tableName, ' WHERE fileID=''', filename, ''';'];
        checkTableCursor = exec(conn, sqlCheckQuery);
        fetchCheckTableCursor = fetch(checkTableCursor);
        isFile = fetchCheckTableCursor.Data;
        close(checkTableCursor); close(fetchCheckTableCursor);
        if (~strcmp(isFile,'No Data')) %returned something
            FilesSkipped = FilesSkipped+1;
            continue;
        end
        %Now that there's definitely a table and the file's not in the DB already, try opening the stro
        try
            stro = nex2stro(filename);
        catch
            BrokenFileAdder(brokenFilesConn, filename, 'stro cannot open the file', 'NULL');
            brokenfiles = brokenfiles+1;
            continue;
        end
        if (~isfield(stro,'sum'))  % If there's a error
            BrokenFileAdder(brokenFilesConn, filename, 'no sum field in stro', 'NULL');
            brokenfiles = brokenfiles+1;
            continue;
        else
            [tempRFX, tempRFY, bf] = dbRFSorter(PID, stro);
            if bf
                BrokenFileAdder(brokenFilesConn, filename, 'receptive field error', int2str([tempRFX, tempRFY]));
                brokenfiles = brokenfiles+1;
                continue;
            else
                rfX = tempRFX;
                rfY = tempRFY;
            end
            %all file names are random other than the .nex at the end, so here's my way of trying to pull the date/subject ID out of the name as best I can
            tempdate = filename(isstrprop(filename, 'digit'));
            if ge(length(tempdate), 6)
                tempdate = tempdate(1:6);
                subjectID = filename(1:strfind(filename, tempdate)-1);
                recordingDate = ['20',tempdate(end-1:end), '-', tempdate(1:2), '-', tempdate(3:4)];
                sqlInsertQuery = ['INSERT INTO ', tableName, ' VALUES (NULL, ''', filename,''',''', recordingDate, ''',''',subjectID,''',', rfX,', ', rfY, ',1, NULL, NULL,', int2str(PID), ');'];
            else
                recordingDate = 'NULL';
                sqlInsertQuery = ['INSERT INTO ', tableName, ' VALUES (NULL, ''', filename,''',', recordingDate, ','''',', rfX,', ', rfY, ',1, NULL, NULL,', int2str(PID), ');'];
            end
            %finally add file info to database
            addToTableCursor = exec(conn, sqlInsertQuery);
            FilesEntered = FilesEntered + 1;
            close(addToTableCursor);
        end
    end
end
close(conn); close(brokenFilesConn);
disp('Insert to DBs completed. Files stored in Nex_Paradigm_Sort DB, broken files entered into BrokenFiles DB on Synology server');
disp(['Stats- Tables Created: ', int2str(TablesCreated), ' Files Entered: ', int2str(FilesEntered), ' Files Skipped: ', int2str(FilesSkipped), ' Files Broken: ', int2str(brokenfiles)]);