%This is part of lmtfpaperfigures, and was used to generate figure 6.
function [diffs] = boxAndWhisker(X)
subj_names = ['E', 'G', 'U', 'A']; %just hardcoded from my memory from running LMTFpaperfigures last, can make this more flexible in the future
diffs_a = [];
diffs_b = [];
diffs_c = [];
figure;
hold on;
for i = 1:length(X)
    subj = X{i};
    a = subj(:,6)-subj(:,5);
    diffs_a(i,:) = [mean(a), std(a)/sqrt(length(a))];
    b = subj(:,7)-subj(:,6);
    diffs_b(i,:) = [mean(b), std(b)/sqrt(length(b))];
    c = subj(:,8)-subj(:,7);
    diffs_c(i,:) = [mean(c), std(b)/sqrt(length(c))];
end
mean_diffs = [diffs_a(:,1); diffs_b(:,1); diffs_c(:,1)];
ste_diffs =  [diffs_a(:,2); diffs_b(:,2); diffs_c(:,2)];
e = errorbar(mean_diffs, ste_diffs, '.');
for i = 1:length(e.XData)
    switch i
        case {1,5,9}
            label = [subj_names(1), ' '];
        case {2, 6, 10}
            label = [subj_names(2), ' '];
        case {3, 7, 11}
            label = [subj_names(3), ' '];
        case {4, 8, 12}
            label = [subj_names(4), ' '];
        otherwise
            keyboard;
    end
    text(e.XData(i), e.YData(i), label, 'horizontalalignment', 'right', 'color', e.Color);
end
set(gca,'xlim', [0 12.5]);
r = refline(0,0);
set(r, 'color', 'k');
set(r, 'LineStyle', ':');
l1 = line([4.5 4.5],[-.01 .006]);
set(l1, 'color', 'k');
l2 = line([8.5 8.5],[-.01 .006]);
set(l2, 'color', 'k');
set(gca, 'xticklabels', []);
%fix below to label the axes
text(2, -.011, 'trt-rt');
text(6, -.011, 'dtrt - trt');
text(10, -.011, 'ydtrt - dtrt');
end

% %taken from prettypaircomp
% diffs = nan*ones(size(X,2),size(X,2));
% p = nan*ones(size(X,2),size(X,2));
% for i = 1:size(X,2)
%     for j = 1:size(X,2)
%         if strcmp(type,'t-test')
%             diffs(i,j) = nanmean(X(:,i)-X(:,j));
%             [h,p(i,j)] = ttest(X(:,i)-X(:,j));
%         else % type is wilcoxon
%             diffs(i,j) = nanmedian(X(:,i)-X(:,j));
%             p(i,j) = signrank((X(:,i)-X(:,j)));
%         end
%     end
% end
% %L = ['unconstrained''11+2*n''10+3*n(1p1)''10+3*n(theta)''rampy trough''tilted rampy trough''double TRT' 'yoked dTRT'];
% figure;
% %hist(diffs);