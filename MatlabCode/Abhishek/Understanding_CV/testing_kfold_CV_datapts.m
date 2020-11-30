% To understand the impact of data points of Leave-one-out CV
% Author - Abhishek De, 2/20
close all; clearvars;

x = linspace(0,10,1001);
y = x+randn(size(x));
figure(1),subplot(231); 
plot(x,y,'k.'); hold on; lsline; axis square; set(gca,'Tickdir','out'); xlabel('X'); ylabel('Y'); hold off;

ptstotest = [10 30 50 100 300 500 1000];
pred_error = cell(1,numel(ptstotest)); % cell-array for storing the CV errors
for ii = 1:numel(ptstotest)
    tmperror = [];
    for k = 1:100 % Implementing loops because I 
        ind = randi(numel(x),[ptstotest(ii) 1]);
        x_data = x(ind);
        y_data = y(ind);
        c = cvpartition(x_data,'LeaveOut');
        
        for jj = 1:c.NumTestSets
            slope = x_data(c.training(jj))/y_data(c.training(jj));
            tmperror = [tmperror; sum((y_data(c.test(jj)) - slope*x_data(c.test(jj))).^2)];
        end
        pred_error{ii} = tmperror;
    end
    % plotting the results
    subplot(232); plot(ptstotest(ii)*ones(size(tmperror)),tmperror,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
    subplot(233); histogram(tmperror,'DisplayStyle','stairs','EdgeColor',[ii/numel(ptstotest) 0 0]); hold on;
    subplot(234); plot(ptstotest(ii),median(tmperror),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
    subplot(235); plot(ptstotest(ii),var(tmperror),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
    subplot(236); plot(ptstotest(ii),mean(tmperror),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
end
subplot(232); xlabel('CV points'); ylabel('Error'); axis square; set(gca,'XScale','log','Tickdir','out'); hold off;
subplot(233); xlabel('Error'); ylabel('#'); axis square; set(gca,'YScale','log','Tickdir','out'); hold off;
subplot(234); xlabel('CV points'); ylabel('Median error'); set(gca,'XScale','log','Tickdir','out'); axis square; hold off;
subplot(235); xlabel('CV points'); ylabel('Error variance'); set(gca,'XScale','log','Tickdir','out'); axis square; hold off;
subplot(236); xlabel('CV points'); ylabel('Mean error'); set(gca,'XScale','log','Tickdir','out'); axis square; hold off;

% Results of this simulation:
% Median error decreases with N
% Mean error decreases with N
% Variance of error decreases with N
