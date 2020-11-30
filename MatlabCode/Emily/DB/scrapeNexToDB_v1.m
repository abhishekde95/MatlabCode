%% Section 1 - Walks through entire directory of nex files and creates a list of subdirectories
datapath = 'N:\NexFiles';
dbname = 'FullParadigmSortedDB';
brokenfiles = {};
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
%% Section 2: OK, now we have a complete list of subdirectories that we can go through, one by one, looking for files to dump into the db, dumping them as we go.
conn = database('FullParadigmSortedDB','admin','adminVector123','Vendor','MySql','Server','localhost'); %connect to database
disp('connected');
outerWaiting = waitbar(0, 'outer run');
for k = 1:length(subdirs)
    try
        cd(subdirs{k});
    catch
        continue;
    end
    dirstruct = dir;
    for l = 1:length(dirstruct)
        waitbar(l/k, outerWaiting, sprintf('%f',[int2str(k), 'out of ', length(subdirs), ' and ', int2str(l), ' out of ', length(dirstruct)]));
        filename = dirstruct(l).name;
        if strcmp(filename, '.') || strcmp(filename, '..') %not real files
            continue;
        elseif(length(filename) < 4) || ~(strcmp(filename(end-3:end),'.nex')) %these files aren't useful in the slightest.
            continue;
        end
        try
            PID = getparadigmID(filename);
        catch
            brokenfiles{end+1} = {filename, 'PID not obtainable'};
            continue;
        end
        if (size(PID,1)>1 || isnan(PID)) %sometimes files have something like a header, but they're screwy and (in reality) broken
            brokenfiles{end+1} = {filename, 'Header partially corrupted'};
            continue;
        else
            %sort out whether a table already exists in the database for that paradigmID.
            try
                paradigmStringName = paradigmLibrary(PID);
            catch
                brokenfiles{end+1} = {filename, 'PID not in paradigmLibrary'};
            end
            tableName = paradigmStringName(1:end-5);
            ifTableCursor = exec(conn, 'SHOW FULL TABLES IN FullParadigmSortedDB;', 10);
            fetchedIfTableCursor = fetch(ifTableCursor);
            checkIfTable = fetchedIfTableCursor.Data;
            close(ifTableCursor); close(fetchedIfTableCursor);
            isTable = strcmpi(tableName, checkIfTable(:,1));
            %If not, create a new one.
            if (sum(isTable)==0)
                %creating a generic table with basic stro information included and nothing else
                %no depth, r values, FP info, c/l/b values, nx/stx, etc.
                %Tables can be modified by hand later.
                uid = 'UID INT NOT NULL AUTO_INCREMENT';
                fID = 'fileID VARCHAR(40) NOT NULL DEFAULT ''''';
                rDate = 'recDate DATE NULL DEFAULT NULL';
                sID = 'subjID VARCHAR(20) NOT NULL DEFAULT ''''';
                rfx = 'rfX INT NULL DEFAULT 0';
                rfy = 'rfY INT NULL DEFAULT 0';
                qual = 'quality BOOL NOT NULL DEFAULT 1';
                notes = 'notes VARCHAR(30) NULL DEFAULT NULL';
                neuron = 'neuron INT UNSIGNED NULL DEFAULT NULL';
                pID = 'paradigmID INT NOT NULL DEFAULT 0';
                sqltablequery = ['CREATE TABLE IF NOT EXISTS ', tableName, ' (', uid,', ',fID,', ',rDate,', ',sID,', ',rfx,', ',rfy,', ',...
                    qual,', ',notes,', ',neuron,', ',pID,', PRIMARY KEY (UID), CONSTRAINT UNIQUE (fileID, UID));'];
                createNewTableCursor = exec(conn, sqltablequery);
                TablesCreated = TablesCreated + 1;
                close(createNewTableCursor);
            else
                %If table exists, check that the file isn't there already.
                sqlCheckQuery = ['SELECT * FROM ', tableName, ' WHERE fileID=''', filename, ''';'];
                checkTableCursor = exec(conn, sqlCheckQuery);
                fetchCheckTableCursor = fetch(checkTableCursor);
                isFile = fetchCheckTableCursor.Data;
                close(checkTableCursor); close(fetchCheckTableCursor);
                if (~strcmp(isFile,'No Data')) %returned something
                    %disp([filename, ' already in ', tableName, ' table']);
                    FilesSkipped = FilesSkipped+1;
                    continue;
                end
            end
            %Now that there's definitely a table and the file's not in the DB already, try opening the stro
            try
                stro = nex2stro(filename);
            catch
                brokenfiles{end+1} = {filename, ' stro does not open the file'};
                continue;
            end
            if (~isfield(stro,'sum'))  % If there's a error
                brokenfiles{end+1} = {filename, ' no sum field in stro'};
                continue
            else %If table exists and file isn't there, add new file to existing table.
                %RF info is named randomly, depending on the paradigm.
                if(PID==212)
                    RFX = stro.sum.exptParams.rfPosX;
                    RFY = stro.sum.exptParams.rfPosY;
                elseif(PID==211)
                    RFX = stro.sum.exptParams.RFPosX;
                    RFY = stro.sum.exptParams.RFPosY;
                elseif (PID==157 || PID==213)
                    RFX = stro.sum.exptParams.stim_x;
                    RFY = stro.sum.exptParams.stim_y;
                elseif (PID==158)
                    RFX = stro.sum.exptParams.targ_x;
                    RFY = stro.sum.exptParams.targ_y;
                else
                    RFX = stro.sum.exptParams.rf_x;
                    RFY = stro.sum.exptParams.rf_y;
                end
                %convert any NaNs to NULLs so SQL can read them
                if isnan(RFX)
                    RFX='NULL';
                else
                    RFX=int2str(RFX);
                end
                if size(RFX,2)>1 || size(RFY,2)>1
                    brokenfiles{end+1} = {filename, 'RFX is stored incorrectly'};
                    continue;
                end
                if isnan(RFY)
                    RFY='NULL';
                else
                    RFY=int2str(RFY);
                end
                %all file names are random other than the MMDDYY000 and the .nex at the end, so here's my way of trying to get around that.
                try
                    getDateIdx = regexp(filename, '\d{6}');
                    stringDate = filename(getDateIdx:getDateIdx+5);
                    recordingDate= int2str(['20',stringDate(end-1:end),'-',stringDate(1:2),'-',stringDate(end-3:end-2)]);
                    subjectID = filename(1:getDateIdx-1);
                catch
                    keyboard;
                    recordingDate = 'NULL';
                    subjectID = filename;
                end
                %dump to db
                sqlInsertQuery = ['INSERT INTO ', tableName, ' VALUES (NULL, ''', filename,''',', recordingDate, ',''',subjectID,''',', RFX,', ', RFY, ',1, NULL, NULL,', int2str(PID), ');'];
                addToTableCursor = exec(conn, sqlInsertQuery);
                if (size(addToTableCursor.Message,1)>0)
                    keyboard;
                end
                disp(addToTableCursor.Message);
                FilesEntered = FilesEntered + 1;
                close(addToTableCursor);
            end
        end
    end
end
close(conn); close(outerWaiting);
disp('disconnected');
disp(['Stats- Tables Created: ', int2str(TablesCreated), ' Files Entered: ', int2str(FilesEntered), ' Files Skipped: ', int2str(FilesSkipped), ' FilesBroken: ', int2str(size(brokenfiles,2))]);
if (size(brokenfiles,1)>0)
    disp('these files could not be inserted into the Database for some reason');
    for q = 1:size(brokenfiles,2)
        disp(brokenfiles{q});
    end
    disp('insert to DB completed');
else
    disp('insert to DB completed with no issues');
end