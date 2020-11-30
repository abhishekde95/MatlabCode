%This was what I used to generate the anova comparison for Apollo and Utu
%for the short and long stimulus duration data I collected. We ended up
%using another analysis, but I've kept this around anyway in case it ended
%up being useful.
function anovan_stimdurcomp()
sIDs = {'A', 'U'};
conds = {'RG', 'LUM'};
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
for s = 1:length(sIDs)
   for c = 1:length(conds)
       if strcmp(sIDs{s}, 'A')
           if strcmp(conds{c}, 'RG')
               long = 'std inputs 666 chr t1';
               short = 'std inputs chr 333 take 1';
           else
               long = 'std inputs 666 take 1';
               short = 'std inputs 333 take 1';
           end
       else
           if strcmp(conds{c}, 'RG')
               long = 'std inputs 666 chr';
               short = 'std inputs 333 chr';
           else
               long = 'std inputs 666 lum';
               short = 'std inputs 333 lum';
           end
       end
       long_flist_q = sprintf('SELECT fileID from LMTF WHERE subjID = ''%s'' AND quality = 1 AND notes = ''%s''', sIDs{s}, long);
       short_flist_q = sprintf('SELECT fileID from LMTF WHERE subjID = ''%s'' AND quality = 1 AND notes = ''%s''', sIDs{s}, short);
       long_flist = fetch(conn, long_flist_q);
       short_flist = fetch(conn, short_flist_q);
       data_long = getLMTFrawdata(long_flist);
       data_short = getLMTFrawdata(short_flist);
       long_thresh = log(sqrt(data_long(:,1).^2+data_long(:,2).^2));
       short_thresh = log(sqrt(data_short(:,1).^2+data_short(:,2).^2));
       tot_thresh = [long_thresh; short_thresh];
       dur = [repmat(666,size(data_long,1),1); repmat(333,size(data_short,1),1)];
       TFs = [data_long(:,3); data_short(:,3)];
       anovan(tot_thresh, {TFs, dur});
   end
end
end