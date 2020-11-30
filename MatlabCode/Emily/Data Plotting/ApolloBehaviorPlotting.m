%This was written Pre-Database - so no telling if it's useful. However,
%it's probably better to rewrite this entirely.
%graphing script - insert stro, output plot - not very useful for LMTF
%data, but very helpful when visualizing Apollo's ability levels when
%performing IsoSamp behavior. 
textFileList = {'ApolloIsoTemp.txt'};
startpath = 'C:/NO BACKUP/NexFiles/Emily/Apollo';
for iterator = 1:length(textFileList)
    fileName = textFileList{iterator};
    flist = flatten(fnamesFromTxt2(fullfile(nexfilepath,'nexfilelists', 'Emily', fileName)));
    for f = 1:length(flist)
        filestro = nex2stro(findfile(flist{f},startpath));
        name = strsplit(filestro.sum.fileName, '\');
        %collecting info from file into matrix - correct, Lcc, mcc, tfs
        if filestro.sum.paradigmID == 157 %is lmtf
            compmtx = [filestro.trial(:,16), filestro.trial(:,(11:13))];
            numtrials = num2str(size(compmtx, 1));
            stitch_name = strcat(name{end}, ' lmtf ', numtrials);
        elseif filestro.sum.paradigmID == 107 %is isosamp
            compmtx = [filestro.trial(:,8), filestro.trial(:,(20:21)), filestro.trial(:,17)];
            numtrials = num2str(size(compmtx, 1));
            stitch_name = strcat(name{end}, ' isoSamp ', numtrials);
        else
            disp('not an accepted stro type');
            return; %trying this - test!
        end
        
        %removing duplicate info while also averaging percent of correct responses per l m tf condition
        cut_compmtx = compmtx(:,2:4);
        unique_cut_compmtx = unique(cut_compmtx, 'rows');
        us = size(unique_cut_compmtx, 1);
        combomtx = [];
        for i = 1:us
            idxs = find(ismember(cut_compmtx, unique_cut_compmtx(i,:), 'rows'));
            corr = [];
            for j = 1:size(idxs,1)
                fullrow = compmtx(idxs(j), :);
                corr = [corr, fullrow(1)];
            end
            p_corr = mean(corr);
            combomtx(i, :) = [p_corr, unique_cut_compmtx(i, :)];
        end
        %disp(combomtx);
        
        %plotting info into bubble color/size plot
        ms = combomtx(:,1)*200;
        %fixing markersizes that are size 0
        small_ms = find(ms==0);
        ms(small_ms) = 1;
        figure; colormap(parula(200));
        graph = scatter3(combomtx(:, 2), combomtx(:,3), combomtx(:, 4), ms, ms, 'filled', 'MarkerEdgeColor', 'k');
         colorbar;
        title(stitch_name); xlabel('L'); ylabel('M'); zlabel('TF'); set(gca, 'zscale', 'log');
        drawnow;
    end %loop over one text file
end %loop over all text files