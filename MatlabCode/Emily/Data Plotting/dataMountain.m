%Not useful for anything. overcomplicated answer to simple problem. 
%Data Mountain: Sets of files get pulled from the database according to 
% different (hardcoded) requirements and are organized into maps by subject 
% and condition. Specifically, each map represents a subject, and each key
% represents a different subset of filenames. 
% For example, E_map contains: 
%   subject name, 
%   list of all filenames corresponding to 50,0 collected on any monitor
%   list of filenames corresponding to 50,0 collected on only CRT monitors
%   list of all filenames corresponding to 50,0 collected on propixx
%   list of filenames corresponding to 50,0 collected on propixx without
%       the eyetracker
%   list of filenames corresponding to 50,0 collected on propixx with
%       the eyetracker
% In this way, any combination of datasets can be stored and collected and
% organized in the same fashion, to be passed further below to the general 
% analysis function. 11/29/16 EG
function dataMountain()
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12'); %connect to sorted files database
all_human_flist = fetch(conn, 'SELECT fileID FROM LMTF WHERE subjID IN(''E'', ''G'', ''Ab'', ''Z'', ''P'') AND quality = 1 AND rfX = 50 AND rfY = 0');
human_map = containers.Map('subj', 'E,G,Ab,Z,P', 'all files', all_human_flist);
all_monkey_flist = fetch(conn, 'SELECT fileID FROM LMTF WHERE subjID IN(''A'', ''N'', ''S'', ''F'', ''U'') AND quality = 1 AND rfX = 50 AND rfY = 0');
monkey_map = containers.Map('subj', 'A,N,S,F,U', 'all files', all_monkey_flist);
total_map = [human_map, monkey_map];
% titles = {'humans- r = G g = Ab b = Z k = E m = P', 'monkeys- r = A g = N b = S k = F m = U','total files: h = black, m = blue'};
E_map = containers.Map('subj', 'E'); A_map = containers.Map('subj', 'A'); G_map = containers.Map('subj', 'G'); U_map = containers.Map('subj', 'U'); P_map = containers.Map('subj', 'P'); 
Z_map = containers.Map('subj', 'Z'); N_map = containers.Map('subj', 'N'); F_map = containers.Map('subj', 'F'); S_map = containers.Map('subj', 'S'); Ab_map = containers.Map('subj', 'Ab');
list = [Ab_map, F_map, S_map, N_map, Z_map, P_map, A_map, E_map, G_map, U_map];
for i = 1:length(list)
    map = list(i);
    subjID = map('subj');
    all_files_query = sprintf('CALL all_50_0(''%s'')', subjID);
    af = fetch(conn, all_files_query);
    crt_files_query = sprintf('CALL all_50_0_CRT(''%s'')', subjID);
    cf = fetch(conn, crt_files_query);
    pp_files_query = sprintf('CALL all_50_0_PP(''%s'')', subjID);
    apf = fetch(conn, pp_files_query);
    first_map = containers.Map('all files', af, 'crt files', cf, 'all pp files', apf);
    map = [map, first_map];
    if ~isempty(map('all pp files'))
        pp_noET_files_query = sprintf('SELECT fileID FROM LMTF WHERE subjID = ''%s'' AND recDate > ''2016-05-27'' AND quality = 1 AND rfX = 50 AND rfY = 0 AND notes NOT LIKE ''eye%''', subjID);
        pnf = fetch(conn, pp_noET_files_query);
        pp_onlyET_files_query = sprintf('SELECT fileID FROM LMTF WHERE subjID = ''%s'' AND recDate > ''2016-05-27'' AND quality = 1 AND rfX = 50 AND rfY = 0 AND notes LIKE ''eye%''', subjID);
        pef = fetch(conn, pp_onlyET_files_query);
        second_map = containers.Map('pp no ET files', pnf, 'pp ET files', pef);
        map = [map; second_map];
    end
    list(i) = map;
end
close(conn);
list = [list total_map];
generalLMTFanalysis(list);
end

function generalLMTFanalysis(list)
for i = 1:length(list)
    
    binnedDataPlot(data{i}, 1, dots, title) %multiplot
    %binnedDataPlot(data{i}, 0, dots, title) %single plot
end
% tf_sum = unique(tfs);
end