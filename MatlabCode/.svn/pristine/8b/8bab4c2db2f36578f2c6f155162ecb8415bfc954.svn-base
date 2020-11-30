%5/8/18 - this never ended up getting used for anything, and is probably
%repalced entirely by iterateAndPlotFiles_modularPlusDB.
%General LMTF Analysis: The goal of this is to have a master function which
% calls whatever other LMTF analysis functions exist, as necessary, on
% whichever datasets are deemed necessary. Because this is nebulous, this
% function will probably grow out of control. As of 11/30/16, all the
% program does is collect filenames into a map (via getLMTFfilesMaps()),
% and then converts the file lists, still within the map, to a data matrix
% (via getLMTFrawdata()) EG

function generalLMTFanalysis()
%Determine what kinds of data we already have:
%% DISPLUMs (collapsing across stim durs and ets)
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12'); %connect to sorted files database
dispLums = fetch(conn, 'SELECT DISTINCT dispLum FROM LMTF');
filesSum = {}; filesSum_50only = {}; filesSum_50only_A = {};filesSum_50only_G = {};filesSum_50only_E = {};
for i = 1:length(dispLums)
    if isnan(dispLums{i})
        continue;
    end
    fileReq = sprintf('SELECT fileID FROM LMTF WHERE dispLum = %d', dispLums{i});
    filesSum{i} = fetch(conn, fileReq);
    fileReq_50only = sprintf('SELECT fileID FROM LMTF WHERE dispLum = %d AND rfX = 50 AND rfY = 0', dispLums{i});
    filesSum_50only{i} = fetch(conn, fileReq_50only);
    fileReq_50only_A = sprintf('SELECT fileID FROM LMTF WHERE dispLum = %d AND rfX = 50 AND rfY = 0 AND subjID = ''A''', dispLums{i});
    filesSum_50only_A{i} = fetch(conn, fileReq_50only_A);
    fileReq_50only_G = sprintf('SELECT fileID FROM LMTF WHERE dispLum = %d AND rfX = 50 AND rfY = 0 AND subjID = ''G''', dispLums{i});
    filesSum_50only_G{i} = fetch(conn, fileReq_50only_G);
    fileReq_50only_E = sprintf('SELECT fileID FROM LMTF WHERE dispLum = %d AND rfX = 50 AND rfY = 0 AND subjID = ''E''', dispLums{i});
    filesSum_50only_E{i} = fetch(conn, fileReq_50only_E);
end
close(conn);
g_50_90 = filesSum_50only_G{4};
g_50_52 = filesSum_50only_G{7};
e_50_72 = filesSum_50only_E{1};
e_50_90 = filesSum_50only_E{4};
A_50_72 = filesSum_50only_A{1};
A_50_90 = filesSum_50only_A{4};
A_50_52 = filesSum_50only_A{7};
dispLum_titles = {'g_50_90','g_50_52','e_50_72','e_50_90','A_50_72','A_50_90','A_50_52'};
dispLum_files = {g_50_90, g_50_52, e_50_72, e_50_90, A_50_72, A_50_90, A_50_52};
dispLum_data = {};
for j = 1:length(dispLum_files)
    dispLum_data{j} = getLMTFrawdata(dispLum_files{j});
    binnedDataPlot(dispLum_data{j}, 1,1,1, dispLum_titles{j});
end
gvg_9052 = {dispLum_data{1}, dispLum_data{2}};
eve_9072 = {dispLum_data{3}, dispLum_data{4}};
ava_7290 = {dispLum_data{5}, dispLum_data{6}};
ava_9052 = {dispLum_data{6}, dispLum_data{7}};
avava_729052 = {dispLum_data{5}, dispLum_data{6}, dispLum_data{7}};
gva_9052 = {dispLum_data{1}, dispLum_data{2}, dispLum_data{6}, dispLum_data{7}};
gva_52 = {dispLum_data{2}, dispLum_data{7}};
gva_90 = {dispLum_data{1}, dispLum_data{6}};
eva_7290 = {dispLum_data{3}, dispLum_data{4}, dispLum_data{5}, dispLum_data{6}};
eva_90 = {dispLum_data{4}, dispLum_data{6}};
eva_72 = {dispLum_data{3}, dispLum_data{5}};
avgve_90 = {dispLum_data{1}, dispLum_data{4}, dispLum_data{6}};
bdp_loop = {gvg_9052,eve_9072,ava_7290,ava_9052,avava_729052,gva_9052,gva_52,gva_90,eva_7290,eva_90,eva_72, avgve_90};
bdp_loop_titles = {'g90=b g52=r','e72=b e90=r','a72=b a90=r','a90=b a52=r','a72=b a90=r a52=g','g90=b g52=r a90=g a52=k','g52=b a52=r','g90=b a90=r','e72=b e90=r a72=g a90=k', ...
    'e90=b a90=r','e72=b a72=r', 'g90=b e90=r a90=g'};
for b = 1:length(bdp_loop)
    binnedDataPlot(bdp_loop{b}, 1, 1, 0, bdp_loop_titles{b});
end

%IMPORTANT NOTE: STIMDURS ARE DIRECTLY CORRELATED TO DISPLUMS - DL 52 IS SD
%656, DL 72 IS SD 576, DL 90 IS SD 668
%%  eye tracker: only at dl 90, only at sd 668 only at 50,0
%3 graphs for g, e: et, noet, both
%3 graphs g+e: et, noet, both
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12'); %connect to sorted files database
subjs_et = fetch(conn, 'SELECT DISTINCT subjID FROM LMTF WHERE notes LIKE ''eye%%'' AND rfX = 50 AND rfY = 0 AND quality = 1');
subjs_noet = fetch(conn, 'SELECT DISTINCT subjID FROM LMTF WHERE notes NOT LIKE ''eye%%'' AND rfX = 50 AND rfY = 0 AND quality = 1');
allsubjs = subjs_noet(ismember(subjs_noet, subjs_et));
eye_tracker_map = containers.Map();
for a = 1:length(allsubjs)
    files_query = sprintf('SELECT fileID FROM LMTF WHERE subjID = ''%s'' AND notes LIKE ''eye%%'' AND rfX = 50 AND rfY = 0 AND quality = 1', allsubjs{a});
    files = fetch(conn, files_query);
    data = getLMTFrawdata(files);
    non_files_query = sprintf('SELECT fileID FROM LMTF WHERE subjID = ''%s'' AND notes NOT LIKE ''eye%%'' AND rfX = 50 AND rfY = 0 AND quality = 1', allsubjs{a});
    non_files = fetch(conn, non_files_query);
    non_data = getLMTFrawdata(non_files);
    key_et = sprintf('ET %s', allsubjs{a});
    key_noet = sprintf('No ET %s', allsubjs{a});
    et_keys = {key_et, key_noet};
    et_values = {data, non_data};
    temp_map = containers.Map(et_keys, et_values);
    eye_tracker_map = [eye_tracker_map; temp_map];
end
close(conn);
keyset = keys(eye_tracker_map);
for d = 1:length(keyset)
    binnedDataPlot(eye_tracker_map(keyset{d}), 1,1,1, keyset{d});
end
noet_vals = strfind(keyset, 'No ET');
noetIDX = find(~cellfun(@isempty,noet_vals));
etIDX = find(cellfun(@isempty,noet_vals));
noet_dat = {}; et_dat = {}; tot_dat = {};
for n = 1:length(noetIDX)
    k = keyset{noetIDX(n)};
    noet_dat{end+1} = eye_tracker_map(k);
    tot_dat{end+1} = eye_tracker_map(k);
end
for e = 1:length(etIDX)
    k = keyset{etIDX(e)};
    et_dat{end+1} = eye_tracker_map(k);
    tot_dat{end+1} = eye_tracker_map(k);
end
% naming isn't great but it should still be clear
binnedDataPlot(noet_dat, 1, 0, 1, cell2mat(keyset(noetIDX)));
binnedDataPlot(et_dat, 1, 0, 1, cell2mat(keyset(etIDX)));
binnedDataPlot(tot_dat, 1, 0, 1, cell2mat(keyset));

binnedDataPlot(noet_dat, 1, 1, 0, cell2mat(keyset(noetIDX)));
binnedDataPlot(et_dat, 1, 1, 0, cell2mat(keyset(etIDX)));
binnedDataPlot(tot_dat, 1, 1, 0, cell2mat(keyset));

%% etc
% map = getLMTFfilesMaps();
% full_keyset = keys(map);
% %rewrites file lists with data matrix
% for k = 1:length(full_keyset)
%     current_key = full_keyset{k};
%     data = getLMTFrawdata(map(current_key));
%     map(current_key) = data;
% end
%bdp_titles = {'humans- r = G g = Ab b = Z k = E m = P', 'monkeys- r = A g = N b = S k = F m = U','total files: h = black, m = blue'};
%          if size(subjID) < 2
%              title = bdp_titles
%              dots = 0;
%          else
%              title = subjID;
%              dots = 1;
%          end
%          binnedDataPlot(data{k}, 1, dots, title) %multiplot
%          binnedDataPlot(data{k}, 0, dots, title) %single plot
end
