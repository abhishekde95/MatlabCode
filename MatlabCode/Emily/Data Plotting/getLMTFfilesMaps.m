%5/8/18 - not useful for anything. do not use.
%getLMTFfilesMaps: Sets of files get pulled from the database according to 
% different (hardcoded) requirements and are organized into a map by 
% condition. Specifically, each key represents a subject and condition,
% and each value represents the relevant subset of filenames. Returns a
% single map.
% In this way, any combination of datasets can be stored and collected and
% organized in the same fashion, to be passed to the general analysis 
% function. 11/29/16 EG (originally named Data Mountain)
% 11/30/16 added map clearing functionality - removes all keys that have no
% associated values. Also organized so that there would be only one map
% returned, rather than a series of maps.
function general_map = getLMTFfilesMaps()
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12'); %connect to sorted files database
keyset = fetch(conn, 'SELECT DISTINCT subjID FROM LMTF');
keyset_string = sprintf('''%s'',', keyset{:});
valueset_query = sprintf('SELECT fileID FROM LMTF WHERE subjID IN (%s)', keyset_string(1:end-1));
valueset = fetch(conn, valueset_query);
general_map = containers.Map(keyset, valueset);
for m = 1:length(keyset)
    subjID = keyset{m};
    all_files_query = sprintf('CALL all_50_0(''%s'')', subjID);
    af = fetch(conn, all_files_query);
    af_key = sprintf('all files %s', subjID);
    crt_files_query = sprintf('CALL all_50_0_CRT(''%s'')', subjID);
    cf = fetch(conn, crt_files_query);
    cf_key = sprintf('crt files %s', subjID);
    pp_files_query = sprintf('CALL all_50_0_PP(''%s'')', subjID);
    apf = fetch(conn, pp_files_query);
    apf_key = sprintf('all pp files %s', subjID);
    keyset = {af_key, cf_key, apf_key};
    valueset = {af, cf, apf};
    temp_map = containers.Map(keyset, valueset);
    if ~isempty(apf)
        pp_noET_files_query = sprintf('SELECT fileID FROM LMTF WHERE subjID = ''%s'' AND recDate > ''2016-05-27'' AND quality = 1 AND rfX = 50 AND rfY = 0 AND notes NOT LIKE ''eye%%''', subjID);
        pnf = fetch(conn, pp_noET_files_query);
        pnf_key = sprintf('all pp no ET files %s', subjID);
        pp_onlyET_files_query = sprintf('SELECT fileID FROM LMTF WHERE subjID = ''%s'' AND recDate > ''2016-05-27'' AND quality = 1 AND rfX = 50 AND rfY = 0 AND notes LIKE ''eye%%''', subjID);
        pef = fetch(conn, pp_onlyET_files_query);
        pef_key = sprintf('all pp ET files %s', subjID);
%         if sum(ismember(pnf, apf)) == size(pnf,1) || sum(ismember(pef, apf)) == size(pef,1)
%             continue;
%         else
            keyset = {pnf_key, pef_key};
            valueset = {pnf, pef};
            second_map = containers.Map(keyset, valueset);
            temp_map = [temp_map; second_map];
%         end
    end
    general_map = [general_map; map];
end
%cleaning the maps to only contain keys that correspond to values
%i.e. if there is no file list associated with the 
total_keyset = keys(general_map);
for e = 1:general_map.Count
    key = total_keyset{e};
    if isempty(general_map(key))
        remove(general_map, key);
    end
end
close(conn);
end