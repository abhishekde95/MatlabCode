function filename_additional = leftoverfiles(filename)
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
allfiles = fetch(conn,'SELECT filename FROM WNthresh');
comments = fetch(conn,'SELECT comments FROM WNthresh');
NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
close(conn);
filename_additional  = allfiles(strcmp(NTmode,'subunit'));
filename_additional(ismember(filename_additional,filename)) = [];
end

