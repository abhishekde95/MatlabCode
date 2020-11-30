%Called by "analysis". Part of an overly complicated solution to a problem that didn't need any
%of this nonsense. Do not use.
function [actual_raw_data_map, distinct_IDs] = pullDataForAnalysis()
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12'); %connect to sorted files database
%keyset = {'ADL_52', 'ADL_72', 'ADL_90', 'EDL_72', 'EDL_90', 'EET_Y', 'EET_N', 'GET_Y', 'GET_N', 'ALL_H', 'ALL_M'};
%keyset = {'ADL_52', 'ADL_72', 'GET_Y80', 'GET_N80', 'GET_Y5060','GET_N5060','GET_Y50n60', 'GET_N50n60'};

keyset = {'GET_Y', 'GET_N', 'ADL_52', 'ADL_72', 'A2DL_52', 'A2DL_90', 'A3DL_72', 'A3DL_90', 'G1DL_52', 'N1DL_52', 'G2DL_52', 'S2DL_52', 'N3DL_52', 'S3DL_52'};
% DLs = {52, 72, 90}; %doing this here because whatever. Could be better but whatever.
% for d = 1:length(DLs)
%     human_q = sprintf('SELECT DISTINCT subjID FROM LMTF WHERE nhp = 0 AND dispLum = %d', DLs{d});
%     humansList = fetch(conn, human_q);
%     primate_q = sprintf('SELECT DISTINCT subjID FROM LMTF WHERE nhp = 1 AND dispLum = %d', DLs{d});
%     primateList = fetch(conn, primate_q);
%     for h = 1:length(humansList)
%         temp = sprintf('HEACH%d_%s', DLs{d}, humansList{h});
%         keyset{end+1} = temp;
%     end
%     for p = 1:length(primateList)
%         temp = sprintf('MEACH%d_%s', DLs{d}, primateList{p});
%         keyset{end+1} = temp;
%     end
% end
valueset = {}; distinct_IDs = {};
rfx50 = ' AND rfX = 50 AND rfY = 0 AND quality = 1';
for k = 1:length(keyset)
    full_key = keyset{k};
    split_key = strsplit(full_key, '_');
    distinct_IDs{end+1} = split_key{1};
    det_subj = 'subjID = ';
    if strcmp(distinct_IDs{end}, 'HEACH') || strcmp(distinct_IDs{end}, 'MEACH')
        subj_ID = sprintf('''%s''', split_key{2});
    elseif strcmp(distinct_IDs{end}, 'ALL')
        det_subj = '';
        subj_ID = '';
    else
        subj = distinct_IDs{end};
        subj_ID = sprintf('''%s''', subj(1));
    end
    which_data = split_key{2};
    switch which_data
        case 'Y80'
            rfx50 = ' AND rfX = 80 AND rfY = 0 AND quality = 1';
            which_data = ' AND notes LIKE ''eye%%''';
        case 'N80'
            rfx50 = ' AND rfX = 80 AND rfY = 0 AND quality = 1';
            which_data = ' AND notes NOT LIKE ''eye%%''';
        case 'Y5060'
            rfx50 = ' AND rfX = 50 AND rfY = 60 AND quality = 1';
            which_data = ' AND notes LIKE ''eye%%''';
        case 'N5060'
            rfx50 = ' AND rfX = 50 AND rfY = 60 AND quality = 1';
            which_data = ' AND notes NOT LIKE ''eye%%''';
        case 'Y50n60'
            rfx50 = ' AND rfX = 50 AND rfY = -60 AND quality = 1';
            which_data = ' AND notes LIKE ''eye%%''';
        case 'N50n60'
            rfx50 = ' AND rfX = 50 AND rfY = -60 AND quality = 1';
            which_data = ' AND notes NOT LIKE ''eye%%''';
        case '52'
            which_data = ' AND dispLum = 52';
        case '90'
            which_data = ' AND dispLum = 90';
        case '72'
            which_data = ' AND dispLum = 72';
        case 'Y'
            which_data = ' AND notes LIKE ''eye%%''';
        case 'N'
            which_data = ' AND notes NOT LIKE ''eye%%''';
        case 'H'
            which_data = ' nhp = 0';
        case 'M'
            which_data = ' nhp = 1';
        otherwise
            which_data = '';
    end
    files_q = sprintf('SELECT fileID FROM LMTF WHERE %s%s%s%s', det_subj, subj_ID, which_data, rfx50);
    disp(files_q);
    filenames = fetch(conn, files_q);
    temp_data = getLMTFrawdata(filenames);
    valueset{end+1} = temp_data;
end
actual_raw_data_map = containers.Map(keyset, valueset);
distinct_IDs = unique(distinct_IDs);
close(conn);
end