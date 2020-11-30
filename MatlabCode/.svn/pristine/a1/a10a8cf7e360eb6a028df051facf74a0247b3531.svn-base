% Comparing which is a better metric: difference in medians vs. median of differences 
% Author - GDLH, 1/20

close all; clearvars;

% Comparing difference in medians versus median of differences
niter = 2000;
n = 20;
data = [];
for i = 1:niter
   % a = normrnd(0,1,n,2);
    a = exprnd(1,n,2);
    data(i,:) = [median(a(:,1))-median(a(:,2)) median(a(:,1)-a(:,2))];
end
bins = linspace(min(data(:)), max(data(:)),30);
figure; axes; hold on;
hist(data,bins); legend({'diff. of med.','med. of diff.'}); set(gca,'Tickdir','out'); axis square; hold off;
