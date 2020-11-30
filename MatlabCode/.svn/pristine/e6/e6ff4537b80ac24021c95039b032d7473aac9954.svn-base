%Called by nex2db. Adds broken files to separate database.
function BrokenFileAdder(conn, filename, reason, info)
    bf_check_query = ['SELECT * FROM BrokenFiles WHERE FileName = ''', filename, ''';'];
    isInDB = fetch(conn, bf_check_query);
    if ~isempty(isInDB)
        if strcmp(info, 'NULL')
            bf_insert_query = ['INSERT INTO BrokenFiles VALUES(''', filename, ''', ''', reason, ''', NULL, NULL);'];
        else
            bf_insert_query = ['INSERT INTO BrokenFiles VALUES(''', filename, ''', ''', reason, ''', ''' , info, ''', NULL);'];
        end
        bf_insert_Cursor = exec(conn, bf_insert_query);
        close(bf_insert_Cursor);
    end
end