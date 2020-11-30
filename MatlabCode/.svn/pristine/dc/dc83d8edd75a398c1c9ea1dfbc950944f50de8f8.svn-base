% Checking non-stationarity: Population Analyses of cells
% Author - Abhishek De, 2/20

close all; clearvars;
plot_counter = 1;
load baselineFRstats.mat % Loading the baseline FRs from all trials for each cell
load RSSE_linearmodel_CV.mat % Robust regression - linear model
load RSSE_quadmodel_CV.mat % Robust regression - quadratic model

% Loading the identities of the cells
load newLUMidx.mat
load newDOidx.mat
load newhardtoclassifyidx.mat
load newSOidx.mat
LUMidx = newLUMidx;
DOidx = newDOidx;
hardtoclassifyidx = [newSOidx; newhardtoclassifyidx]; 

r = []; p = [];
RSSEisoresp_medianofratios = [];
for ii=1:numel(baselineFRstats)
    
    % Computing the Spearman's r
    [tmp_r, tmp_p] = corr((1:numel(baselineFRstats{ii}))',baselineFRstats{ii},'type','Pearson');
    r = [r; tmp_r];
    p = [p; tmp_p];
    
    % computation for calculating median of differences/ratios
    RSSEisoresp_medianofratios = [RSSEisoresp_medianofratios; median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii})];
end

% Plotting the results 
figure(plot_counter); set(gcf,'Name','Assesing non-stationarity');
subplot(321); histogram(p,logspace(-20,0,21),'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'Linewidth',2);
axis square; set(gca,'Tickdir','out','XScale','log','Xlim',[10.^(-20) 1]); xlabel('p val'); ylabel('count'); title('ALL');
subplot(322); histogram(r,linspace(-1,1,21),'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'Linewidth',2); hold on;
histogram(r(p<0.0001),linspace(-1,1,21),'FaceColor',[0 0 0],'EdgeColor',[0 0 0],'Linewidth',2); title('ALL'); 
axis square; set(gca,'Tickdir','out','Xlim',[-1 1]); xlabel('Spearman-r'); ylabel('count'); legend('p>0.0001','p<0.0001');
subplot(323), histogram(abs(r(LUMidx)),linspace(0,1,21),'FaceColor',[0 0 0],'EdgeColor',[1 1 1],'Linewidth',1);
axis square; set(gca,'Tickdir','out','Xlim',[0 1]); xlabel('Spearman-r'); ylabel('count'); title('LUM');
subplot(324), histogram(abs(r(DOidx)),linspace(0,1,21),'FaceColor',[1 0 0],'EdgeColor',[1 1 1],'Linewidth',1);
axis square; set(gca,'Tickdir','out','Xlim',[0 1]); xlabel('Spearman-r'); ylabel('count'); title('DO');
subplot(325), histogram(abs(r(hardtoclassifyidx)),linspace(0,1,21),'FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1],'Linewidth',1);
axis square; set(gca,'Tickdir','out','Xlim',[0 1]); xlabel('Spearman-r'); ylabel('count'); title('HTC');
subplot(326); plot(abs(r(LUMidx)),RSSEisoresp_medianofratios(LUMidx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(abs(r(DOidx)),RSSEisoresp_medianofratios(DOidx),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(abs(r(hardtoclassifyidx)),RSSEisoresp_medianofratios(hardtoclassifyidx),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Xlim',[0 0.6],'Ylim',[0.01 1000],'YScale','log'); xlabel('Spearman-r'); ylabel('NLI');
plot_counter = plot_counter + 1;

