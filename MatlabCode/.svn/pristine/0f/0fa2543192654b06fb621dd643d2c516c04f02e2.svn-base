%% Section 1 - Walks through entire directory of nex files and creates a list of subdirectories
datapath = 'N:\NexFiles';
brokenfiles = 0;
TablesCreated=0;
FilesEntered=0;
FilesSkipped=0;
%BE IN N:\NEXFILES OR IT WONT RUN
cd(datapath);
% First making a gigantic dirstruct
dirstruct = dir;
subdirs = {};
for i = 1:size(dirstruct,1) % Getting the initial list of directories
    if (dirstruct(i).isdir && ~strncmp(dirstruct(i).name, '.',1))
        cd(dirstruct(i).name);
        subdirs{length(subdirs)+1} = [datapath, filesep,dirstruct(i).name];
        cd('../');
    else
        continue
    end
end
alldone = 0;
startsubdir = 1;
while ~alldone
    nsubdirs = length(subdirs);
    for i = startsubdir:length(subdirs)
        eval(['cd (''',subdirs{i},''')']);
        dirstruct = dir;
        for j = 1:size(dirstruct,1)
            if (dirstruct(j).isdir && ~strncmp(dirstruct(j).name, '.',1));
                datapath = eval('pwd');
                subdirs{length(subdirs)+1} = [datapath, filesep,dirstruct(j).name];
            end
        end
    end
    if length(subdirs) == nsubdirs
        alldone = 1;
    else
        startsubdir = nsubdirs;
    end
end
%% Section 2: OK, now we have a complete list of subdirectories that we can go through, one by one, looking for files to dump into the dbs, dumping them as we go.
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12'); %connect to sorted files database
disp('connected to Nex_Paradigm_Sort database');
totfilelist = {};
for k = 1:length(subdirs)
    cd(subdirs{k}); 
    validFiles = struct2cell(dir('*.nex'));
    validFileNames = validFiles(1,:);
    for l = 1:length(validFileNames)
        filename = validFileNames{l};
        totfilelist{end+1} = filename;
        try
            PID = getparadigmID(filename);
        catch
            BrokenFileAdder(conn, filename, 'No Paradigm ID Available', 'NULL');
            brokenfiles = brokenfiles+1;
            continue;
        end
        if (size(PID,1)>1 || isnan(PID)) %sometimes files have something like a header, but they're screwy and (in reality) broken
            BrokenFileAdder(conn, filename, 'Header Partially Corrupted', 'NULL');
            brokenfiles = brokenfiles+1;
            continue;
        else %continue trying to add the file
            try
                paradigmStringName = paradigmLibrary(PID);
            catch
                BrokenFileAdder(conn, filename, 'PID not in paradigmLibrary', int2str(PID));
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
                BrokenFileAdder(conn, filename, 'stro cannot open the file', 'NULL');
                brokenfiles = brokenfiles+1;
                continue;
            end
            if (~isfield(stro,'sum'))  % If there's a error
                BrokenFileAdder(conn, filename, 'no sum field in stro', 'NULL');
                brokenfiles = brokenfiles+1;
                continue;
            else
                [tempRFX, tempRFY, bf] = dbRFSorter(PID, stro);
                if bf
                    BrokenFileAdder(conn, filename, 'receptive field error', int2str([tempRFX, tempRFY]));
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
end
close(conn); close(brokenFilesConn);
disp('Insert to DBs completed. Files stored in Nex_Paradigm_Sort DB on Synology server');
disp(['Stats- Tables Created: ', int2str(TablesCreated), ' Files Entered: ', int2str(FilesEntered), ' Files Skipped: ', int2str(FilesSkipped), ' Files Broken: ', int2str(brokenfiles)]);