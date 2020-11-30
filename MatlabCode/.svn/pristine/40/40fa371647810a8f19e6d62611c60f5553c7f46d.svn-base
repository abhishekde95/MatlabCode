%avg difference of residuals
%comparing thresholds for stimdurs 333 and 666 apollo and utu

function [val,pval] = stimdur_diffs()
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
subjs = {'A', 'U'};
for s = 1:length(subjs)
    disp(subjs{s});
    if s == 1
        all = {'pre 333 training', '333 first run'};
    else
        all = {'std run pre 333', '333 first run'};
    end
    thresh_each = {};
    for i = 1:length(all)
        query = sprintf('SELECT fileID FROM LMTF WHERE subjID = ''%s'' AND rfX = 50 AND rfY = 0 AND quality = 1 AND notes LIKE (''%s'')', subjs{s}, all{i});
        flist = fetch(conn, query);
        disp(size(flist));
        full_data = getLMTFrawdata(flist);
        Lcc = full_data(:,1);
        Mcc = full_data(:,2);
        thresh = sqrt(Mcc.^2 + Lcc.^2);
        %avg_thresh = mean(thresh);
        %thresh_each(i) = avg_thresh;
        thresh_each{i} = thresh;
    end
    disp(thresh_each);
    [val, pval] = ttest2(thresh_each{1}, thresh_each{2},'tail', 'both');
    disp(val);
    disp(pval);
%     if val
%         disp(pval);
%         figure;
%         plot(thresh_each,'.');
%         title([num2str(tfbound),'hz, min ', num2str(minbound), ', max ', num2str(maxbound),', pval ', num2str(pval)]);
%     end
end
end