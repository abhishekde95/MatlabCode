%Three one off functions for examining LMTF data- top displays a list of eye win widths and heights
%for all subjects, second displays grouped list of eyewin widths and
%heights for apollo and utu, and third checks the fpoff - stimoff time
%window to ensure that it's 100-600 ms.

conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
%flist = fetch (conn, 'SELECT fileID FROM LMTF WHERE recDate > ''2016-05-27'' AND recDate < ''2017-02-14'' AND quality = 1 AND subjID IN(''A'', ''N'', ''S'', ''F'', ''U'');');
flist = fetch (conn, 'SELECT fileID, subjID FROM LMTF WHERE recDate < ''2017-02-14'' AND quality = 1;');
disp(size(flist));
eyewin_size = {};
for i = 1:length(flist)
    file = flist{i,1};
    subj = flist{i,2};
    fullPath = findfile(file);
    stro = nex2stro(fullPath);
    width = stro.sum.exptParams.eyewin_w;
    height = stro.sum.exptParams.eyewin_h;
    add = sprintf('UPDATE LMTF SET eyewin_width = %d, eyewin_height = %d WHERE fileID = ''%s'';', width, height, file);
    whadd = exec(conn, add);
    if ~isempty(whadd.Message)
        disp(whadd.Message);
        keyboard;
    end
end
whList = fetch(conn, 'SELECT DISTINCT eyewin_width, eyewin_height, subjID FROM LMTF WHERE recDate < ''2017-02-14'' AND quality = 1;');
disp(whList);
%subjs = fetch(conn, 'SELECT DISTINCT subjID FROM LMTF WHERE recDate < ''2017-02-14'' AND quality = 1;');
%% Section 2 displays grouped eyewin widths and heights for Apollo and Utu, only looking at files for paper
subjs = {'A', 'U'};
for i = 1:length(subjs)
    currs = subjs{i};
    whSubj = sprintf('SELECT DISTINCT eyewin_width, eyewin_height FROM LMTF WHERE recDate < ''2017-02-14'' AND quality = 1 AND subjID = ''%s'';', currs);
    whCount = fetch(conn, whSubj);
    disp(currs);
    disp(whCount);
    fileQ = sprintf('SELECT COUNT(fileID) FROM LMTF WHERE recDate < ''2017-02-14'' AND quality = 1 AND subjID = ''%s'' GROUP BY eyewin_width, eyewin_height', currs);
    whFile = fetch(conn, fileQ);
    disp(whFile);
    whFileMat = cell2mat(whFile);
    numf = sum(whFileMat);
    perc_std = whFileMat(1)/numf;
    perc_diff = sum(whFileMat(2:end))/numf;
    disp(perc_std);
    disp(perc_diff);
end

uniquewhcount = fetch(conn, 'SELECT eyewin_width, eyewin_height, COUNT(fileID) FROM LMTF WHERE recDate < ''2017-02-14'' AND subjID IN(''A'', ''U'') AND quality = 1 GROUP BY eyewin_width, eyewin_height');
disp(uniquewhcount);

%% Section 3: time range for fpoff - stimoff = 100-600ms
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
flist = fetch (conn, 'SELECT fileID, nhp FROM LMTF WHERE recDate < ''2017-02-14'' AND quality = 1;');
eyewin_sizes_h = [];
eyewin_sizes_m = [];
trange = [];
for i = 1:length(flist)
    file = flist{i,1};
    nhp = flist{i,2};
    fullfile = findfile(file);
    stro = nex2stro(fullfile);
    min_time = round(min(stro.trial(:,3) - stro.trial(:,5)), 4);
    max_time = max(stro.trial(:,3) - stro.trial(:,5));
    trange = [trange; min_time max_time];
%     if ~isempty(stro.ras)
%         w = stro.sum.exptParams.eyewin_w;
%         h = stro.sum.exptParams.eyewin_h;
%         if nhp
%             eyewin_sizes_m = [eyewin_sizes_m; w h];
%         else
%             if w == 7 && h == 7 || w == 15 && h == 15
%                 disp(file);
%             end
%             eyewin_sizes_h = [eyewin_sizes_h; w h];
%         end
%     end
end
plot(trange(:,1), trange(:,2));
disp(mean(trange(:,1)));
disp(mean(trange(:,2)));