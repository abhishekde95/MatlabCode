% We ended up not using this analysis in the long run, but potentially
% useful
%compares old 666, new 666, and 333 stimdur for apollo and utu by comparing models
%generated using the various stimdur collected data and a generic subset of
%tfs and thetas. one run r666 - 333, one run o666 - 333, one run o666-r666;
function [preds_1, preds_2, preds_3] = comp_short_long_stimdur()
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
subjs = {'A'}; %A
preds_1 = {};
preds_2 = {};
preds_3 = {};
for s = 1:length(subjs)
    disp(subjs{s});
    %if s == 1
    all = {'second 666 attempt training', 'std inputs 666 take 1', 'std inputs 333 take 1'}; %'propixx'
    %else
    %all = {'std run pre 333', '333 first run', 'postbreak fresh run'};
    %end
    %evaluate model at fixed set of points (tfs and thetas), thetas 0-2pi, tf 1-40
    [tfs,thetas] = meshgrid(linspace(1,40,20),linspace(0,2*pi,20));
    %figure; hold on;
    for p = 1:length(all)
        query = sprintf('SELECT fileID FROM LMTF WHERE subjID = ''%s'' AND rfX = 50 AND rfY = 0 AND quality = 1 AND notes LIKE (''%s'')', subjs{s}, all{p});
        try
            flist = fetch(conn, query);
        catch
            query = 'SELECT fileID FROM LMTF WHERE subjID = ''U'' AND rfX = 50 AND rfY = 0 AND quality = 1 AND (notes LIKE ''postbreak fresh run'' OR notes LIKE ''post resids culling'')';
            flist = fetch(conn, query);
        end
        data = getLMTFrawdata(flist);
        disp(size(data));
        modelParams = LMTFfitting(data, 1);
        if p == 1 %full data
            symbol = '.';
            color = 'black';
        elseif p == 2 %r666 std inputs
            symbol = 'o';
            color = 'green';
        else %333 std inputs
            symbol = 'o';
            color = 'magenta';
        end
        L = data(:,1);
        M = data(:,2);
        TF = data(:,3);
        Loog = logical(data(:,4));
%         
%         plot3(L(~Loog),M(~Loog),TF(~Loog),'ko','MarkerFaceColor',color,'MarkerSize',5, 'Marker', symbol);
%         plot3(-L(~Loog),-M(~Loog),TF(~Loog),'ko','MarkerFaceColor',color,'MarkerSize',5, 'Marker', symbol);
%         plot3(L(Loog),M(Loog),TF(Loog),'ro','MarkerFaceColor','red','MarkerSize',5, 'Marker', symbol);
%         plot3(-L(Loog),-M(Loog),TF(Loog),'ro','MarkerFaceColor','red','MarkerSize',5, 'Marker', symbol);
%         
        % Adjusting the axes so we can see everything
%         lim = max(max(abs([L M])));
%         set(gca,'Xlim',1.1*[-lim lim]);
%         set(gca,'Zscale','log');
%         set(gca,'Ylim',1.1*[-lim lim]);
%         axis square
%         xlabel('L-cone contrast');
%         ylabel('M-cone contrast');
%         zlabel('Temporal Frequency (Hz)');
%         set(gca,'View',[135 20]);
%         axis vis3d
        preds = [];
        preds = LMTF_thresh_from_model(reshape(thetas, [], 1), reshape(tfs, [], 1), modelParams);
        switch all{p}
            case {'second 666 attempt training', 'std run pre 333'}
                preds_1 = {preds, [L, M, TF]};
            case {'std inputs 666 take 1'}
                preds_2 = {preds, [L, M, TF]};
            case {'std inputs 333 take 1'}
                preds_3 = {preds, [L, M, TF]};
            otherwise
                keyboard;
        end
    end
    %plot the difference between his 666 (new and old, separately) thresholds and his 333 thresholds as a function of TF and theta
    %     figure;
    %     plot3(reshape(tfs, [], 1), reshape(thetas, [], 1), [preds_r666-preds_short], '.');
    %     titlestring = sprintf('''%s'' recent 666 - 333', subjs{s});
    %     title(titlestring); xlabel('tfs'); ylabel('thetas'); zlabel('r diffs');
    %     figure;
    %     plot3(reshape(tfs, [], 1), reshape(thetas, [], 1),[preds_o666-preds_short], '.');
    %     titlestring = sprintf('''%s'' old 666 - 333', subjs{s});
    %     title(titlestring); xlabel('tfs'); ylabel('thetas'); zlabel('o diffs');
    %plotting and values from greg earlier
    %> 0 first was further, < 0 first was closer, 0 means identical. is it random or clean? Expect negative, (666-333).
    %     diff1 = log10(preds_r666) - log10(preds_short);
    %     figure;
    %     hist(diff1);
    %     title([subjs{s},' recent 666-333']);
    %     disp([subjs{s}, ' recent 666-333 ', num2str(max(abs(diff1)))]);
    %     diff2 = log10(preds_o666) - log10(preds_r666);
    %     figure;
    %     hist(diff2);
    %     title([subjs{s},' old 666- recent 666']);
    %     disp([subjs{s},' old 666- recent 666 ', num2str(max(abs(diff2)))]);
    %     diff3 = log10(preds_o666) - log10(preds_short);
    %     figure;
    %     hist(diff3);
    %     title([subjs{s},' old 666 - 333']);
    %     disp([subjs{s},' old 666 - 333 ', num2str(max(abs(diff3)))]);
    %     fprintf(['mean "recent 666-333": ', num2str(mean(diff1)), ' ,mean "old 666- recent 666": ', num2str(mean(diff2)), ' ,mean "old 666 - 333": ', num2str(mean(diff3)), '\n']);
end
end