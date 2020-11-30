%this was my attempt to find lmtf files with headers that differed from their
%collected data. Sort of useful but didn't end up getting used much.
function mismatched = matchHeaderToRecords(filenameList)
mismatched = {};    
    for i = 1:length(filenameList)
        file = findfile(filenameList{i},'N:/NexFiles');
        stro = nex2stro(file);
        headerX = abs(stro.sum.exptParams.stim_x);
        headerY = abs(stro.sum.exptParams.stim_y);
        allRecX = stro.trial(:,9);
        recX = abs(allRecX(1));
        allRecY = stro.trial(:,10);
        recY = abs(allRecY(1));
        if headerX == recX && headerY == recY
            continue;
        else
            mismatched{i} = {filenameList{i}, headerX, headerY, recX, recY};
        end
    end
end