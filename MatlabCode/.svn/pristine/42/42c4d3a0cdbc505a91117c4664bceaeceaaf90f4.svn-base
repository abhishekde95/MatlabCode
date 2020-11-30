%comparing diffs of thresh between monk and human for a subset of data
% tfbound - max tf the data can be (only looking at lower slice of tfs)
% min/max bound - min/max LM angle (only looking at subset of color dirs)
% determine via ttest whether two sets of threshs have sig diff
function [val,pval] = mvhdiffs()
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
%all = {'G', 'Ab', 'Z', 'E', 'P', 'A', 'N', 'S', 'F', 'U'};
%all = {'U', 'A', 'S', 'N', 'F', 'Z', 'P', 'G', 'E'};
all = {'U', 'A', 'E', 'G'};
tfboundList = [3, 1, 5, 5, 1];
minboundList = [30, 20, 40, 20, 40];
maxboundList = [60, 70, 50, 70, 50];
for k = 1:length(tfboundList)
tfbound = tfboundList(k);
minbound = minboundList(k);
maxbound = maxboundList(k);
thresh_each = [];
subjErr=[];
for i = 1:length(all)
    if strncmp(all{i},'S',1) || strncmp(all{i},'N',1) || strncmp(all{i},'F',1)
        addition = 'AND monID = ''Dell 4''';
    %elseif strncmp(all{i}, 'U',1) || strncmp(all{i}, 'A', 1) ||strncmp(all{i}, 'E', 1) || strncmp(all{i}, 'G', 1)
    %    addition = 'AND monID = ''ProPixx''';
    else
        addition = '';
    end
    query = sprintf('SELECT fileID FROM LMTF WHERE subjID = ''%s'' AND rfX = 50 AND rfY = 0 AND recDate < ''2017-02-14'' AND quality = 1 %s;', all{i}, addition);
    flist = fetch(conn, query);
    if isempty(flist)
        subjErr = [subjErr, all{i}];
        continue;
    else
        full_data = getLMTFrawdata(flist);
        lowtf_data = full_data([full_data(:,3)<=tfbound],:);
        Lcc = lowtf_data(:,1);
        Mcc = lowtf_data(:,2);
        rad_to_degrees = atan2(Mcc,Lcc)*(180/pi);
        thresh = sqrt(Mcc.^2 + Lcc.^2);
        L=(rad_to_degrees>=minbound & rad_to_degrees<=maxbound);
        avg_thresh = mean(thresh(L));
        if isnan(avg_thresh)
            subjErr = [subjErr, all{i}];
        else
            thresh_each(i) = avg_thresh;
        end
    end
end
%disp(thresh_each);
thresh_each(isnan(thresh_each))=[];
disp(['tfbound ', num2str(tfbound), ', minbound ', num2str(minbound), ', maxbound ', num2str(maxbound),': thresh_each:     ', num2str(thresh_each)]);
disp(['didn''t include the following subjects: ', subjErr]);
[val, pval] = ttest2(thresh_each(1:3), thresh_each(4:end),'tail', 'both');
disp(['pval = ', num2str(pval)]);
if val
    figure;
    plot(thresh_each,'.');
    title([num2str(tfbound),'hz, min ', num2str(minbound), ', max ', num2str(maxbound),', pval ', num2str(pval)]);
end
end
end
%change arbitrary thresholds (3, 30-60) to (5/1, 20-70/40-50) <- four
%combos