%used at least by nex2db, if not by more functions. Also useful for
%learning how to add a table to a database from matlab.
function [tableName, TablesCreated] = TableAdder(psn, conn)
tableName = psn(1:end-5);
ifTableCursor = fetch(conn, 'SHOW FULL TABLES IN Nex_Paradigm_Sort;', 10);
isTable = strcmpi(tableName, ifTableCursor(:,1));
if (sum(isTable)==0)
    %creating a generic table with basic stro information included and nothing else
    %Tables can be modified by hand later.
    uid = 'UID INT NOT NULL AUTO_INCREMENT';
    fID = 'fileID VARCHAR(40) NOT NULL DEFAULT ''''';
    rDate = 'recDate DATE NULL DEFAULT NULL';
    sID = 'subjID VARCHAR(20) NOT NULL DEFAULT ''''';
    rfx = 'rfX INT NULL DEFAULT NULL';
    rfy = 'rfY INT NULL DEFAULT NULL';
    qual = 'quality BOOL NOT NULL DEFAULT 1';
    notes = 'notes VARCHAR(30) NULL DEFAULT NULL';
    neuron = 'neuron INT UNSIGNED NULL DEFAULT NULL';
    pID = 'paradigmID INT NOT NULL DEFAULT 0';
    sqltablequery = ['CREATE TABLE IF NOT EXISTS ', tableName, ' (', uid,', ',fID,', ',rDate,', ',sID,', ',rfx,', ',rfy,', ',...
        qual,', ',notes,', ',neuron,', ',pID,', PRIMARY KEY (UID), CONSTRAINT UNIQUE (fileID, UID));'];
    createNewTableCursor = exec(conn, sqltablequery);
    TablesCreated = 1;
    close(createNewTableCursor);
else
    TablesCreated = 0;
end
end