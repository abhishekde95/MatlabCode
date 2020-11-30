%part of the same overcomplicated nonsense of "analysis" and several other
%functions. Not really useful for much.
function countMap = sqlCounter()
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
subjs = fetch(conn, 'SELECT DISTINCT subjID FROM LMTF');
countMap = containers.Map();
for s = 1:length(subjs)
   key1 = sprintf('%s_TotLums', subjs{s}); 
   key3 = sprintf('%s_50Lums', subjs{s}); 
   dispLums_q = sprintf('SELECT DISTINCT dispLum FROM LMTF WHERE subjID = ''%s'' AND quality =1', subjs{s});
   dispLums = fetch(conn, dispLums_q);
   key2s = {}; nums = {};
   for d = 1:length(dispLums)
          key2s{end+1} = sprintf('%s_%d', subjs{s}, dispLums{d});
          nums_q = sprintf('SELECT COUNT(fileID) FROM LMTF WHERE subjID = ''%s'' AND dispLum = %d AND quality = 1', subjs{s}, dispLums{d});
          nums{end+1} = fetch(conn, nums_q);
   end
   dl_50_q = sprintf('SELECT DISTINCT dispLum FROM LMTF WHERE subjID = ''%s'' AND quality =1 AND rfX = 50 AND rfY = 0', subjs{s});
   dl_50 = fetch(conn, dl_50_q);
   key4s = {}; num_50s = {};
   for f = 1:length(dl_50)
       key4s{end+1} = sprintf('%s_%d_50', subjs{s}, dl_50{f});
       num_50_q = sprintf('SELECT COUNT(fileID) FROM LMTF WHERE subjID = ''%s'' AND dispLum = %d AND quality = 1 AND rfX = 50 AND rfY = 0', subjs{s}, dispLums{d});
       num_50s{end+1} = fetch(conn, num_50_q);
   end
   keyset = {key1, key2s{:}, key3, key4s{:}};
   valueset = {dispLums, nums{:}, dl_50, num_50s{:}};
   temp_map = containers.Map(keyset, valueset);
   countMap = [countMap; temp_map];
end
close(conn);

k = keys(countMap);
for m = 1:countMap.Count
    temp = k{m};
    disp(temp);
    disp(countMap(temp));
end
end