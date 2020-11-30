%populates the DB for all LMTF filenames with the name of the monitor used
%during recording. 
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
filenames = fetch(conn, 'SELECT fileID FROM LMTF WHERE recDate < ''2017-02-14'' AND quality = 1;');
for i = 1:length(filenames)  
    file = filenames{i};
    stro = nex2stro(findfile(file));
    monspd = stro.sum.exptParams.mon_spd;
    monspd = reshape(monspd,length(monspd)/3,3);
    tmp = FindMon(monspd);
    update_stmt = sprintf('UPDATE LMTF SET monID = ''%s'' WHERE fileID = ''%s'';', tmp, file);
    exec(conn, update_stmt);
end
d = fetch(conn, 'SELECT COUNT(fileID), subjID, monID, dispLum FROM LMTF WHERE rfX = 50 AND rfY = 0 AND recDate < ''2017-02-14'' AND quality = 1 GROUP BY subjID, monID');
disp(d);