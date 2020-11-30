% A new script for find any correlation between the signal integration within subunits, between the subunits and the isoresponse curve
% Derivative of PopAnalysis_WNthresh_subunit_isoresp.m
% Analysis of the cross validated erros
% Author - Abhishek De, 12/19
close all; clearvars;
plot_counter = 1;

% loading the indices for DO, SO, LUM and hardtoclassifycells
load newDOidx.mat
load newLUMidx.mat
load newSOidx.mat
load newhardtoclassifyidx.mat
DOidx = newDOidx';
LUMidx = newLUMidx';
hardtoclassifyidx = [newhardtoclassifyidx' newSOidx'];
idx = [DOidx LUMidx hardtoclassifyidx];

% Load the isoresponse data
load RSSE_linearmodel_CV.mat % Robust regression
load RSSE_quadmodel_CV.mat
load SSE_linearmodel_CV.mat % Ordinary least squares
load SSE_quadmodel_CV.mat

% Load the integration across the subunits data
load AUROClinsubunits_CV.mat
load AUROCquadsubunits_CV.mat

% Load the integration within the subunit data
load AUROClinS1_CV.mat
load AUROCquadS1_CV.mat
load AUROClinS2_CV.mat
load AUROCquadS2_CV.mat

% Extract the mean and median errors
RSSEisoresp_lin_mean = []; RSSEisoresp_lin_median = []; % Isoresponse data
RSSEisoresp_quad_mean = []; RSSEisoresp_quad_median = [];
Acrosssubunits_lin_mean = []; Acrosssubunits_lin_median = []; % Across subunits data
Acrosssubunits_quad_mean = []; Acrosssubunits_quad_median = [];
Withinsubunits_lin1_mean = []; Withinsubunits_lin1_median = []; % Within subunits data
Withinsubunits_quad1_mean = []; Withinsubunits_quad1_median = [];
Withinsubunits_lin2_mean = []; Withinsubunits_lin2_median = [];
Withinsubunits_quad2_mean = []; Withinsubunits_quad2_median = [];

% For storing median of differences/ ratios
RSSEisoresp_medianofratios = [];
Acrosssubunits_medianofdifferences = [];
Withinsubunits_medianofdifferences = [];

for ii = 1:numel(RSSE_linearmodel)
    
    % Isoresponse data
    RSSEisoresp_lin_mean = [RSSEisoresp_lin_mean; mean(RSSE_linearmodel{ii})];
    RSSEisoresp_lin_median = [RSSEisoresp_lin_median; median(RSSE_linearmodel{ii})];
    RSSEisoresp_quad_mean = [RSSEisoresp_quad_mean; mean(RSSE_quadmodel{ii})];
    RSSEisoresp_quad_median = [RSSEisoresp_quad_median; median(RSSE_quadmodel{ii})];
    
    % Integration across subunits data
    Acrosssubunits_lin_mean = [Acrosssubunits_lin_mean; mean(AUROClinsubunits{ii})];
    Acrosssubunits_lin_median = [Acrosssubunits_lin_median; median(AUROClinsubunits{ii})];
    Acrosssubunits_quad_mean = [Acrosssubunits_quad_mean; mean(AUROCquadsubunits{ii})];
    Acrosssubunits_quad_median = [Acrosssubunits_quad_median; median(AUROCquadsubunits{ii})];
    
    % Integration within subunit data
    Withinsubunits_lin1_mean = [Withinsubunits_lin1_mean; mean(AUROClin1{ii})];
    Withinsubunits_lin1_median = [Withinsubunits_lin1_median; median(AUROClin1{ii})];
    Withinsubunits_quad1_mean = [Withinsubunits_quad1_mean; mean(AUROCquad1{ii})];
    Withinsubunits_quad1_median = [Withinsubunits_quad1_median; median(AUROCquad1{ii})];
    Withinsubunits_lin2_mean = [Withinsubunits_lin2_mean; mean(AUROClin2{ii})];
    Withinsubunits_lin2_median = [Withinsubunits_lin2_median; median(AUROClin2{ii})];
    Withinsubunits_quad2_mean = [Withinsubunits_quad2_mean; mean(AUROCquad2{ii})];
    Withinsubunits_quad2_median = [Withinsubunits_quad2_median; median(AUROCquad2{ii})];
    
    % computation for calculating median of differences/ratios
    RSSEisoresp_medianofratios = [RSSEisoresp_medianofratios; median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii})];
    Acrosssubunits_medianofdifferences = [Acrosssubunits_medianofdifferences; median(AUROCquadsubunits{ii}-AUROClinsubunits{ii})];
    Withinsubunits_medianofdifferences = [Withinsubunits_medianofdifferences; median([median(AUROCquad1{ii}-AUROClin1{ii}) median(AUROCquad2{ii}-AUROClin2{ii})])];
    
end


% Combining the values from each of the subunits 
Withinsubunits_lin_mean = mean([Withinsubunits_lin1_mean Withinsubunits_lin2_mean],2);
Withinsubunits_quad_mean = mean([Withinsubunits_quad1_mean Withinsubunits_quad2_mean],2);
Withinsubunits_lin_median = mean([Withinsubunits_lin1_median Withinsubunits_lin2_median],2);
Withinsubunits_quad_median = mean([Withinsubunits_quad1_median Withinsubunits_quad2_median],2);
figure(plot_counter); set(gcf,'Name','Median values');
subplot(231); plot(RSSEisoresp_medianofratios,100*(Acrosssubunits_medianofdifferences),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 1000],'Ylim',[-2 10]); xlabel('Spatial NLI isoresponse'); ylabel('Across AUROC Quad-Lin');
subplot(232); plot(RSSEisoresp_medianofratios,100*(Withinsubunits_medianofdifferences),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
axis square; set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 1000],'Ylim',[-2 8]); xlabel('Spatial NLI isoresponse'); ylabel('Within AUROC Quad-Lin');
subplot(233); plot(100*(Acrosssubunits_medianofdifferences),100*(Withinsubunits_medianofdifferences),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
axis square; set(gca,'Tickdir','out','Xlim',[-2 10],'Ylim',[-2 8]); xlabel('Across AUROC Quad-Lin'); ylabel('Within AUROC Quad-Lin');
subplot(234); histogram(RSSEisoresp_medianofratios,logspace(-1,3,31),'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); 
set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 1000],'XTick',[0.1 1 10 100 1000]); title('Isoresponse'); xlabel('(Lin/Quad) errors');  axis square;
subplot(235);histogram(100*(Acrosssubunits_medianofdifferences),linspace(-2,8,31),'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-2 10]); xlabel('(Quad-Lin) AUROC'); title('Across subunits'); axis square;
subplot(236);histogram(100*(Withinsubunits_medianofdifferences),linspace(-2,8,31),'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-2 10]); xlabel('(Quad-Lin) AUROC'); title('Within subunits'); axis square;
plot_counter = plot_counter + 1;

% Spearman rank correlation coefficient 
[r1,p1] = corr(RSSEisoresp_medianofratios,100*(Acrosssubunits_medianofdifferences),'type','Spearman');
[r2,p2] = corr(RSSEisoresp_medianofratios,100*(Withinsubunits_medianofdifferences),'type','Spearman');
[r3,p3] = corr(100*(Acrosssubunits_medianofdifferences),100*(Withinsubunits_medianofdifferences),'type','Spearman');

% Pearson's correlation coefficient
[r4,p4] = corr(RSSEisoresp_medianofratios,100*(Acrosssubunits_medianofdifferences),'type','Pearson');
[r5,p5] = corr(RSSEisoresp_medianofratios,100*(Withinsubunits_medianofdifferences),'type','Pearson');
[r6,p6] = corr(100*(Acrosssubunits_medianofdifferences),100*(Withinsubunits_medianofdifferences),'type','Pearson');

%% Breaking the data into simple, DO and hardtoclassify cells
load S1LMS 
load S2LMS
Sconesignal = abs(S1LMS(3,:)) + abs(S2LMS(3,:));
load anglebwvectorsRGB.mat
load anglediffWNchecksubunit.mat
load z_scores.mat

% Comparing the mean and median results across different cell types:
% Looking at the isoresponse data only 
figure(plot_counter); 
subplot(221); plot(RSSEisoresp_lin_median(LUMidx),RSSEisoresp_quad_median(LUMidx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(RSSEisoresp_lin_median(DOidx),RSSEisoresp_quad_median(DOidx),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); plot([0.001 10],[0.001 10],'k');
plot(RSSEisoresp_lin_median(hardtoclassifyidx),RSSEisoresp_quad_median(hardtoclassifyidx),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Xlim',[0.001 10],'Ylim',[0.001 10],'YScale','log','XScale','log'); xlabel('Linear error'); ylabel('Quadratic error'); title('Isoresponse'); hold off;
subplot(222); histogram(RSSEisoresp_medianofratios,logspace(-1,3,31),'EdgeColor',[0 0 0],'Displaystyle','stairs','Linewidth',2); hold on;
histogram(RSSEisoresp_medianofratios(LUMidx),logspace(-1,3,31),'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); hold on;
plot(median(RSSEisoresp_medianofratios(LUMidx)),7,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); title('LUM');
set(gca,'Tickdir','out','Xlim',[0.1 1000],'XScale','log', 'Ylim',[0 35]); axis square; xlabel('Linear error/Quadratic error'); ylabel('#cells');
subplot(223),histogram(RSSEisoresp_medianofratios,logspace(-1,3,31),'EdgeColor',[0 0 0],'Displaystyle','stairs','Linewidth',2); hold on;
histogram(RSSEisoresp_medianofratios(hardtoclassifyidx),logspace(-1,3,31),'FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1]);
plot(median(RSSEisoresp_medianofratios(hardtoclassifyidx)),7,'v','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); title('HTC')
set(gca,'Tickdir','out','Xlim',[0.1 1000],'XScale','log', 'Ylim',[0 35]); axis square; xlabel('Linear error/Quadratic error'); ylabel('#cells');
subplot(224),histogram(RSSEisoresp_medianofratios,logspace(-1,3,31),'EdgeColor',[0 0 0],'Displaystyle','stairs','Linewidth',2); hold on;
histogram(RSSEisoresp_medianofratios(DOidx),logspace(-1,3,31),'FaceColor',[1 0 0],'EdgeColor',[1 1 1]);
plot(median(RSSEisoresp_medianofratios(DOidx)),7,'v','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); title('DO');
set(gca,'Tickdir','out','Xlim',[0.1 1000],'XScale','log', 'Ylim',[0 35]); axis square; xlabel('Linear error/Quadratic error'); ylabel('#cells');
plot_counter = plot_counter + 1;

% Stats: Kruskal Wallis Test
group = [ones(size(LUMidx)) 2*ones(size(DOidx)) 3*ones(size(hardtoclassifyidx))];
p1 = kruskalwallis(RSSEisoresp_medianofratios([LUMidx DOidx hardtoclassifyidx]),group,'off'); % Stats of median of ratios/differences
p8 = kruskalwallis(Withinsubunits_medianofdifferences([LUMidx DOidx hardtoclassifyidx]),group,'off'); % Stats of median of ratios/differences

% Dividing the data into different cell types 
% Looking at the isoresponse, across subunits and within subunits data 
figure(plot_counter); set(gcf,'Name','Median values');
subplot(231); plot(RSSEisoresp_medianofratios(LUMidx),100*(Acrosssubunits_medianofdifferences(LUMidx)),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(RSSEisoresp_medianofratios(DOidx),100*(Acrosssubunits_medianofdifferences(DOidx)),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(RSSEisoresp_medianofratios(hardtoclassifyidx),100*(Acrosssubunits_medianofdifferences(hardtoclassifyidx)),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 1000],'Ylim',[-2 8]); xlabel('Spatial NLI isoresponse'); ylabel('Across AUROC Quad-Lin');
subplot(232); plot(RSSEisoresp_medianofratios(LUMidx),100*(Withinsubunits_medianofdifferences(LUMidx)),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(RSSEisoresp_medianofratios(DOidx),100*(Withinsubunits_medianofdifferences(DOidx)),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(RSSEisoresp_medianofratios(hardtoclassifyidx),100*(Withinsubunits_medianofdifferences(hardtoclassifyidx)),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 1000],'Ylim',[-2 8]); xlabel('Spatial NLI isoresponse'); ylabel('Within AUROC Quad-Lin');
subplot(233); plot(100*(Acrosssubunits_medianofdifferences(LUMidx)),100*(Withinsubunits_medianofdifferences(LUMidx)),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(100*(Acrosssubunits_medianofdifferences(DOidx)),100*(Withinsubunits_medianofdifferences(DOidx)),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(100*(Acrosssubunits_medianofdifferences(hardtoclassifyidx)),100*(Withinsubunits_medianofdifferences(hardtoclassifyidx)),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Xlim',[-2 8],'Ylim',[-2 8]); xlabel('Across AUROC Quad-Lin'); ylabel('Within AUROC Quad-Lin');
subplot(234); histogram(RSSEisoresp_medianofratios(LUMidx),logspace(-1,3,31),'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); hold on;
histogram(RSSEisoresp_medianofratios(DOidx),logspace(-1,3,31),'FaceColor',[1 0 0],'EdgeColor',[1 1 1]);
histogram(RSSEisoresp_medianofratios(hardtoclassifyidx),logspace(-1,3,31),'FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1]);
set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 1000],'XTick',[0.1 1 10 100 1000]); title('Isoresponse'); xlabel('(Lin/Quad) errors');  axis square;
subplot(235);histogram(100*(Acrosssubunits_medianofdifferences(LUMidx)),linspace(-2,8,31),'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); hold on;
histogram(100*(Acrosssubunits_medianofdifferences(DOidx)),linspace(-2,8,31),'FaceColor',[1 0 0],'EdgeColor',[1 1 1]);
histogram(100*(Acrosssubunits_medianofdifferences(hardtoclassifyidx)),linspace(-2,8,31),'FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-2 8]); xlabel('(Quad-Lin) AUROC'); title('Across subunits'); axis square;
subplot(236);histogram(100*(Withinsubunits_medianofdifferences(LUMidx)),linspace(-2,8,31),'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); hold on;
histogram(100*(Withinsubunits_medianofdifferences(DOidx)),linspace(-2,8,31),'FaceColor',[1 0 0],'EdgeColor',[1 1 1]);
histogram(100*(Withinsubunits_medianofdifferences(hardtoclassifyidx)),linspace(-2,8,31),'FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-2 8]); xlabel('(Quad-Lin) AUROC'); title('Within subunits'); axis square;
plot_counter = plot_counter + 1;

% Effect of combined S cone weight
figure(plot_counter); set(gcf,'Name','Effect of S cone signal');
subplot(231); plot(RSSEisoresp_medianofratios,Sconesignal,'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 1000],'Ylim',[0 1.5]); xlabel('Spatial NLI isoresponse'); ylabel('S');
subplot(232); plot(100*(Acrosssubunits_medianofdifferences),Sconesignal,'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
axis square; set(gca,'Tickdir','out','Xlim',[-2 8],'Ylim',[0 1.5]); xlabel('Across AUROC Quad-Lin');  ylabel('S');
subplot(233); plot(100*(Withinsubunits_medianofdifferences),Sconesignal,'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
axis square; set(gca,'Tickdir','out','Xlim',[-2 8],'Ylim',[0 1.5]); xlabel('Within AUROC Quad-Lin');  ylabel('S');
subplot(234); plot(RSSEisoresp_medianofratios(LUMidx),Sconesignal(LUMidx),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(RSSEisoresp_medianofratios(DOidx),Sconesignal(DOidx),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(RSSEisoresp_medianofratios(hardtoclassifyidx),Sconesignal(hardtoclassifyidx),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 1000],'Ylim',[0 1.5]); xlabel('Spatial NLI isoresponse'); ylabel('S');
subplot(235); plot(100*(Acrosssubunits_medianofdifferences(LUMidx)),Sconesignal(LUMidx),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(100*(Acrosssubunits_medianofdifferences(DOidx)),Sconesignal(DOidx),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(100*(Acrosssubunits_medianofdifferences(hardtoclassifyidx)),Sconesignal(hardtoclassifyidx),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Xlim',[-2 8],'Ylim',[0 1.5]); xlabel('Across AUROC Quad-Lin'); ylabel('S');
subplot(236); plot(100*(Withinsubunits_medianofdifferences(LUMidx)),Sconesignal(LUMidx),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(100*(Withinsubunits_medianofdifferences(DOidx)),Sconesignal(DOidx),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(100*(Withinsubunits_medianofdifferences(hardtoclassifyidx)),Sconesignal(hardtoclassifyidx),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Xlim',[-2 8],'Ylim',[0 1.5]); xlabel('Within AUROC Quad-Lin'); ylabel('S');
plot_counter = plot_counter + 1;


% Dependence on z-scores of PC1 from subunit data
figure(plot_counter); set(gcf,'Name','Relationship between WN derived z score of PC1 with other variables');
subplot(221); plot(z_scores(LUMidx),RSSEisoresp_medianofratios(LUMidx),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(z_scores(DOidx),RSSEisoresp_medianofratios(DOidx),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); 
plot(z_scores(hardtoclassifyidx),RSSEisoresp_medianofratios(hardtoclassifyidx),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','YScale','log','Ylim',[0.1 1000]); xlabel('Z-score'); ylabel('NLI');
subplot(222); plot(z_scores(LUMidx),100*(Acrosssubunits_medianofdifferences(LUMidx)),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(z_scores(DOidx),100*(Acrosssubunits_medianofdifferences(DOidx)),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); 
plot(z_scores(hardtoclassifyidx),100*(Acrosssubunits_medianofdifferences(hardtoclassifyidx)),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out'); xlabel('Z-score'); ylabel('Across subunits');
subplot(223); plot(z_scores(LUMidx),100*(Withinsubunits_medianofdifferences(LUMidx)),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(z_scores(DOidx),100*(Withinsubunits_medianofdifferences(DOidx)),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); 
plot(z_scores(hardtoclassifyidx),100*(Withinsubunits_medianofdifferences(hardtoclassifyidx)),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out'); xlabel('Z-score'); ylabel('Within subunits');
subplot(224),histogram(z_scores(LUMidx),linspace(min(z_scores),max(z_scores),21),'EdgeColor',[0 0 0],'Displaystyle','stairs','Linewidth',2); hold on;
plot(median(z_scores(LUMidx)),12,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
histogram(z_scores(DOidx),linspace(min(z_scores),max(z_scores),21),'EdgeColor',[1 0 0],'Displaystyle','stairs','Linewidth',2);
plot(median(z_scores(DOidx)),13,'v','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
histogram(z_scores(hardtoclassifyidx),linspace(min(z_scores),max(z_scores),21),'EdgeColor',[0.5 0.5 0.5],'Displaystyle','stairs','Linewidth',2);
plot(median(z_scores(hardtoclassifyidx)),14,'v','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out'); xlabel('Z-score'); ylabel('Counts');
plot_counter = plot_counter + 1;

% Instead of z-score, I am using the percentile value of the PC1 from the bootstrapped null distribution
load vals.mat
figure(plot_counter); set(gcf,'Name','Relationship between WN derived percentile of PC1 with other variables');
subplot(221); plot(vals(LUMidx),RSSEisoresp_medianofratios(LUMidx),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(vals(DOidx),RSSEisoresp_medianofratios(DOidx),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); 
plot(vals(hardtoclassifyidx),RSSEisoresp_medianofratios(hardtoclassifyidx),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','YScale','log','Ylim',[0.1 1000]); xlabel('Percentile'); ylabel('NLI');
subplot(222); plot(vals(LUMidx),100*(Acrosssubunits_medianofdifferences(LUMidx)),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(vals(DOidx),100*(Acrosssubunits_medianofdifferences(DOidx)),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); 
plot(vals(hardtoclassifyidx),100*(Acrosssubunits_medianofdifferences(hardtoclassifyidx)),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out'); xlabel('Percentile'); ylabel('Across subunits');
subplot(223); plot(vals(LUMidx),100*(Withinsubunits_medianofdifferences(LUMidx)),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(vals(DOidx),100*(Withinsubunits_medianofdifferences(DOidx)),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); 
plot(vals(hardtoclassifyidx),100*(Withinsubunits_medianofdifferences(hardtoclassifyidx)),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out'); xlabel('Percentile'); ylabel('Within subunits');
subplot(224),histogram(vals(LUMidx),linspace(0,100,21),'EdgeColor',[0 0 0],'Displaystyle','stairs','Linewidth',2); hold on;
plot(median(vals(LUMidx)),12,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
histogram(vals(DOidx),linspace(0,100,21),'EdgeColor',[1 0 0],'Displaystyle','stairs','Linewidth',2);
plot(median(vals(DOidx)),13,'v','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
histogram(vals(hardtoclassifyidx),linspace(0,100,21),'EdgeColor',[0.5 0.5 0.5],'Displaystyle','stairs','Linewidth',2);
plot(median(vals(hardtoclassifyidx)),14,'v','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out'); xlabel('Percentile'); ylabel('Counts');
plot_counter = plot_counter + 1;

% Comparing different parameters between the different groups of cells 
% 1) spikes collected during the WN subunit phase
% 2) latency of subunit STA energy
% 3) Number of within-gamut points 
% 4) anglediffWNchecksubunit
% 5) baseline FR
% 6) size of the RF/ subunit STA (in terms of pixels)
load numsubunitspikes.mat
load not_oog_idx_all.mat
load baselineFR.mat
load latencySTAsubunit.mat
load pixelswithinsubunit1.mat
load pixelswithinsubunit2.mat

testpoints_isoresp = [];
for ii = 1:size(not_oog_idx_all,2)
    testpoints_isoresp = [testpoints_isoresp; numel(not_oog_idx_all{1,ii})]; 
end
STAsubunitsize = pixelswithinsubunit1 + pixelswithinsubunit2;

p2 = kruskalwallis(numsubunitspikes([LUMidx DOidx hardtoclassifyidx]),group,'off');
p3 = kruskalwallis(latencySTAsubunit([LUMidx DOidx hardtoclassifyidx]),group,'off');
p4 = kruskalwallis(testpoints_isoresp([LUMidx DOidx hardtoclassifyidx]),group,'off');
p5 = kruskalwallis(anglediffWNchecksubunit([LUMidx DOidx hardtoclassifyidx]),group,'off');
p6 = kruskalwallis(baselineFR([LUMidx DOidx hardtoclassifyidx]),group,'off');
p7 = kruskalwallis(STAsubunitsize([LUMidx DOidx hardtoclassifyidx]),group,'off');

figure(plot_counter); set(gcf,'Name','Differences/Similarities among groups');
subplot(231), boxplot(numsubunitspikes([LUMidx DOidx hardtoclassifyidx]),group); axis square; set(gca,'Tickdir','out','XTick',[1 2 3],'XTicklabel',{'LUM','DO','HTC'}); 
ylabel('num spikes'); title(strcat('p=',num2str(p2,3))); 
subplot(232), boxplot(latencySTAsubunit([LUMidx DOidx hardtoclassifyidx]),group); axis square; set(gca,'Tickdir','out','XTick',[1 2 3],'XTicklabel',{'LUM','DO','HTC'}); 
ylabel('Latency'); title(strcat('p=',num2str(p3,3)));
subplot(233), boxplot(testpoints_isoresp([LUMidx DOidx hardtoclassifyidx]),group); axis square; set(gca,'Tickdir','out','XTick',[1 2 3],'XTicklabel',{'LUM','DO','HTC'}); 
ylabel('Within gamut pts'); title(strcat('p=',num2str(p4,3))); 
subplot(234), boxplot(anglediffWNchecksubunit([LUMidx DOidx hardtoclassifyidx]),group); axis square; set(gca,'Tickdir','out','XTick',[1 2 3],'XTicklabel',{'LUM','DO','HTC'}); 
ylabel('Angle diff check-sub'); title(strcat('p=',num2str(p5,3))); 
subplot(235), boxplot(baselineFR([LUMidx DOidx hardtoclassifyidx]),group); axis square; set(gca,'Tickdir','out','XTick',[1 2 3],'XTicklabel',{'LUM','DO','HTC'}); 
ylabel('Baseline FR'); title(strcat('p=',num2str(p6,3))); 
subplot(236), boxplot(STAsubunitsize([LUMidx DOidx hardtoclassifyidx]),group); axis square; set(gca,'Tickdir','out','XTick',[1 2 3],'XTicklabel',{'LUM','DO','HTC'}); 
ylabel('STA size'); title(strcat('p=',num2str(p7,3))); 
plot_counter = plot_counter + 1;

% Relationship between PC1 zscore and % vals
figure(plot_counter); set(gcf,'Name','Relationship bw z_score and percentile');
subplot(221); plot(z_scores,vals,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
axis square; xlabel('Z score'); ylabel('PC1 percentile'); set(gca,'Tickdir','out','YScale','log','Xlim',[-5 20]); title('ALL');
subplot(222); plot(z_scores(LUMidx),vals(LUMidx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
axis square; xlabel('Z score'); ylabel('PC1 percentile'); set(gca,'Tickdir','out','YScale','log','Xlim',[-5 20]); title('LUM');
subplot(223); plot(z_scores(DOidx),vals(DOidx),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; xlabel('Z score'); ylabel('PC1 percentile'); set(gca,'Tickdir','out','YScale','log','Xlim',[-5 20]); title('DO');
subplot(224); plot(z_scores(hardtoclassifyidx),vals(hardtoclassifyidx),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); 
axis square; xlabel('Z score'); ylabel('PC1 percentile'); set(gca,'Tickdir','out','YScale','log','Xlim',[-5 20]); title('HTC');
plot_counter = plot_counter + 1;

%**********************************************************************%
% A control analyses to check what the results look like when I remove the DO and LUM cells which have high PC1- meaning > 95% eigenvalue
mhardtoclassifyidx = [LUMidx(vals(LUMidx)>=95) DOidx(vals(DOidx)>=95) hardtoclassifyidx];
mLUMidx = LUMidx(vals(LUMidx)<95);
mDOidx = DOidx(vals(DOidx)<95);
mgroup = [ones(size(mLUMidx)) 2*ones(size(mDOidx)) 3*ones(size(mhardtoclassifyidx))];
p9 = kruskalwallis(RSSEisoresp_medianofratios([mLUMidx mDOidx mhardtoclassifyidx]),mgroup,'off'); % Stats of median of ratios/differences
p10 = kruskalwallis(Withinsubunits_medianofdifferences([mLUMidx mDOidx mhardtoclassifyidx]),mgroup,'off'); % Stats of median of ratios/differences
p11 = kruskalwallis(numsubunitspikes([mLUMidx mDOidx mhardtoclassifyidx]),mgroup,'off');
p12 = kruskalwallis(latencySTAsubunit([mLUMidx mDOidx mhardtoclassifyidx]),mgroup,'off');
p13 = kruskalwallis(testpoints_isoresp([mLUMidx mDOidx mhardtoclassifyidx]),mgroup,'off');
p14 = kruskalwallis(anglediffWNchecksubunit([mLUMidx mDOidx mhardtoclassifyidx]),mgroup,'off');
p15 = kruskalwallis(baselineFR([mLUMidx mDOidx mhardtoclassifyidx]),mgroup,'off');
p16 = kruskalwallis(STAsubunitsize([mLUMidx mDOidx mhardtoclassifyidx]),mgroup,'off');

figure(plot_counter); set(gcf,'Name','Results and stats with new classification')
subplot(441); plot(RSSEisoresp_lin_median(mLUMidx),RSSEisoresp_quad_median(mLUMidx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(RSSEisoresp_lin_median(mDOidx),RSSEisoresp_quad_median(mDOidx),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); plot([0.001 10],[0.001 10],'k');
plot(RSSEisoresp_lin_median(mhardtoclassifyidx),RSSEisoresp_quad_median(mhardtoclassifyidx),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Xlim',[0.001 10],'Ylim',[0.001 10],'YScale','log','XScale','log'); xlabel('Linear error'); ylabel('Quadratic error'); title('Isoresponse'); hold off;
subplot(442); histogram(RSSEisoresp_medianofratios,logspace(-1,3,31),'EdgeColor',[0 0 0],'Displaystyle','stairs','Linewidth',2); hold on;
histogram(RSSEisoresp_medianofratios(mLUMidx),logspace(-1,3,31),'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); hold on;
plot(median(RSSEisoresp_medianofratios(mLUMidx)),7,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); title('LUM');
set(gca,'Tickdir','out','Xlim',[0.1 1000],'XScale','log', 'Ylim',[0 35]); axis square; xlabel('Linear error/Quadratic error'); ylabel('#cells');
subplot(443),histogram(RSSEisoresp_medianofratios,logspace(-1,3,31),'EdgeColor',[0 0 0],'Displaystyle','stairs','Linewidth',2); hold on;
histogram(RSSEisoresp_medianofratios(mhardtoclassifyidx),logspace(-1,3,31),'FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1]);
plot(median(RSSEisoresp_medianofratios(mhardtoclassifyidx)),7,'v','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); title('HTC')
set(gca,'Tickdir','out','Xlim',[0.1 1000],'XScale','log', 'Ylim',[0 35]); axis square; xlabel('Linear error/Quadratic error'); ylabel('#cells');
subplot(444),histogram(RSSEisoresp_medianofratios,logspace(-1,3,31),'EdgeColor',[0 0 0],'Displaystyle','stairs','Linewidth',2); hold on;
histogram(RSSEisoresp_medianofratios(mDOidx),logspace(-1,3,31),'FaceColor',[1 0 0],'EdgeColor',[1 1 1]);
plot(median(RSSEisoresp_medianofratios(mDOidx)),7,'v','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); title('DO');
set(gca,'Tickdir','out','Xlim',[0.1 1000],'XScale','log', 'Ylim',[0 35]); axis square; xlabel('Linear error/Quadratic error'); ylabel('#cells');
subplot(445); plot(RSSEisoresp_medianofratios(mLUMidx),100*(Withinsubunits_medianofdifferences(mLUMidx)),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(RSSEisoresp_medianofratios(mDOidx),100*(Withinsubunits_medianofdifferences(mDOidx)),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(RSSEisoresp_medianofratios(mhardtoclassifyidx),100*(Withinsubunits_medianofdifferences(mhardtoclassifyidx)),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 1000],'Ylim',[-2 8]); xlabel('Spatial NLI isoresponse'); ylabel('Within AUROC Quad-Lin');
subplot(446); histogram(100*(Withinsubunits_medianofdifferences(mLUMidx)),linspace(-2,8,31),'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); hold on; 
plot(median(100*(Withinsubunits_medianofdifferences(mLUMidx))),15,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
set(gca,'Tickdir','out','Xlim',[-2 8]); xlabel('(Quad-Lin) AUROC'); title('LUM'); axis square; hold off;
subplot(447); histogram(100*(Withinsubunits_medianofdifferences(mhardtoclassifyidx)),linspace(-2,8,31),'FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1]); hold on;
plot(median(100*(Withinsubunits_medianofdifferences(mhardtoclassifyidx))),15,'v','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); 
set(gca,'Tickdir','out','Xlim',[-2 8]); xlabel('(Quad-Lin) AUROC'); title('HTC'); axis square; hold off;
subplot(448); histogram(100*(Withinsubunits_medianofdifferences(mDOidx)),linspace(-2,8,31),'FaceColor',[1 0 0],'EdgeColor',[1 1 1]); hold on;
plot(median(100*(Withinsubunits_medianofdifferences(mDOidx))),15,'v','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); 
set(gca,'Tickdir','out','Xlim',[-2 8]); xlabel('(Quad-Lin) AUROC'); title('DO'); axis square; hold off;

subplot(449), boxplot(numsubunitspikes([mLUMidx mDOidx mhardtoclassifyidx]),mgroup); axis square; set(gca,'Tickdir','out','XTick',[1 2 3],'XTicklabel',{'LUM','DO','HTC'}); 
ylabel('num spikes'); title(strcat('p=',num2str(p11,3))); 
subplot(4,4,10), boxplot(latencySTAsubunit([mLUMidx mDOidx mhardtoclassifyidx]),mgroup); axis square; set(gca,'Tickdir','out','XTick',[1 2 3],'XTicklabel',{'LUM','DO','HTC'}); 
ylabel('Latency'); title(strcat('p=',num2str(p12,3)));
subplot(4,4,11), boxplot(testpoints_isoresp([mLUMidx mDOidx mhardtoclassifyidx]),mgroup); axis square; set(gca,'Tickdir','out','XTick',[1 2 3],'XTicklabel',{'LUM','DO','HTC'}); 
ylabel('Within gamut pts'); title(strcat('p=',num2str(p13,3))); 
subplot(4,4,12), boxplot(anglediffWNchecksubunit([mLUMidx mDOidx mhardtoclassifyidx]),mgroup); axis square; set(gca,'Tickdir','out','XTick',[1 2 3],'XTicklabel',{'LUM','DO','HTC'}); 
ylabel('Angle diff check-sub'); title(strcat('p=',num2str(p14,3))); 
subplot(4,4,13), boxplot(baselineFR([mLUMidx mDOidx mhardtoclassifyidx]),mgroup); axis square; set(gca,'Tickdir','out','XTick',[1 2 3],'XTicklabel',{'LUM','DO','HTC'}); 
ylabel('Baseline FR'); title(strcat('p=',num2str(p15,3))); 
subplot(4,4,14), boxplot(STAsubunitsize([mLUMidx mDOidx mhardtoclassifyidx]),mgroup); axis square; set(gca,'Tickdir','out','XTick',[1 2 3],'XTicklabel',{'LUM','DO','HTC'}); 
ylabel('STA size'); title(strcat('p=',num2str(p16,3))); 
plot_counter = plot_counter + 1;


%% Some more analyses of the isoresponse and the subunits  
load STAsub1LMS_R.mat
load STAsub2LMS_R.mat
load EigenMatrix1.mat
load EigenMatrix2.mat

NLI = RSSEisoresp_lin_median./RSSEisoresp_quad_median; % Also want to check out the distinctions between linear and nonlinear cells  
Criterion = median(NLI); 

figure(plot_counter); set(gcf,'Name','Rotated subunits of different cells');
for ii = 1:numel(STAsub1LMS_R)
    
    % Plotting the rotated axes of the subunits 
    STAsub1LMS_rotated = STAsub1LMS_R{ii};
    STAsub2LMS_rotated = STAsub2LMS_R{ii};
    
    if any(ismember(newLUMidx,ii)) % Black
         subplot(221); plot3([0; STAsub1LMS_rotated(1,:)],[0;STAsub1LMS_rotated(2,:)],[0; STAsub1LMS_rotated(3,:)],'Color',[0 0 0],'LineWidth',3); hold on;
         subplot(221); plot3([0; STAsub2LMS_rotated(1,:)],[0;STAsub2LMS_rotated(2,:)],[0; STAsub2LMS_rotated(3,:)],'Color',[0 0 0],'LineWidth',3); hold on;
    elseif  any(ismember(newDOidx,ii)) % Red
         subplot(221); plot3([0; STAsub1LMS_rotated(1,:)],[0;STAsub1LMS_rotated(2,:)],[0; STAsub1LMS_rotated(3,:)],'Color',[1 0 0],'LineWidth',3); hold on;
         subplot(221); plot3([0; STAsub2LMS_rotated(1,:)],[0;STAsub2LMS_rotated(2,:)],[0; STAsub2LMS_rotated(3,:)],'Color',[1 0 0],'LineWidth',3); hold on;
    else % Gray 
         subplot(221); plot3([0; STAsub1LMS_rotated(1,:)],[0;STAsub1LMS_rotated(2,:)],[0; STAsub1LMS_rotated(3,:)],'Color',[0.5 0.5 0.5],'LineWidth',3); hold on;
         subplot(221); plot3([0; STAsub2LMS_rotated(1,:)],[0;STAsub2LMS_rotated(2,:)],[0; STAsub2LMS_rotated(3,:)],'Color',[0.5 0.5 0.5],'LineWidth',3); hold on;
    end  
    
end
subplot(221); hold on; line([-1 1],[0 0],[0 0],'Linewidth',3); line([0 0],[-1 1],[0 0],'Linewidth',3); line([0 0],[0 0],[-1 1],'Linewidth',3);
axis equal; axis xy; grid on; xlabel('RL'),ylabel('RM'),zlabel('RS'); title ('All cells');  hold off;

% Plotting the cone weights of all the subunits 
subplot(222); plot(S1LMS(1,newLUMidx),S1LMS(2,newLUMidx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on; 
plot(S2LMS(1,newLUMidx),S2LMS(2,newLUMidx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
plot(S1LMS(1,newDOidx),S1LMS(2,newDOidx),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(S2LMS(1,newDOidx),S2LMS(2,newDOidx),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(S1LMS(1,newhardtoclassifyidx),S1LMS(2,newhardtoclassifyidx),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); 
plot(S2LMS(1,newhardtoclassifyidx),S2LMS(2,newhardtoclassifyidx),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[-1 1],'XTick',-1:0.5:1,'YTick',-1:0.5:1,'Tickdir','out'); 
plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 -1],[0 0],'k'); plot([1 0],[0 -1],'k'); plot([0 -1],[-1 0],'k');  xlabel('L'), ylabel('M');
title('All cells'); 

% Some ROC analyses on the 
Criterions = sort(NLI); Criterions = Criterions(2:end-1);
Y1 = []; Y2 = [];
for ii = 1:numel(Criterions)
    tmp1 = roc(anglebwvectorsRGB(NLI>Criterions(ii)),anglebwvectorsRGB(NLI<=Criterions(ii)));
    tmp2 = roc(anglediffWNchecksubunit(NLI>Criterions(ii)),anglediffWNchecksubunit(NLI<=Criterions(ii)));
    Y1 = [Y1; max([tmp1 1-tmp1])];
    Y2 = [Y2; max([tmp2 1-tmp2])];
end
subplot(223); plot(Criterions,Y1,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); axis square; set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 150]); 
xlabel('Criterion: NLI'); ylabel('Percent Correct ROC'); title('Angle bw RGB subunits');
subplot(224); plot(Criterions,Y2,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); axis square; set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 150]); 
xlabel('Criterion: NLI');  ylabel('Percent Correct ROC'); title('Angle bw check and sub');
plot_counter = plot_counter + 1;

% Calculating inverse M matrix 
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;
M = inv(M');

% Defining NLI : Non-linearity index and looking at the principal directions of the subunits in the 3D RGB space 
nonlin_cells = find(NLI>-1);
figure(plot_counter); set(gcf,'Name','Principal axes of each subunits of non-linear cells');
Anglebw_Princ_LMS1 = []; Anglebw_Princ_LMS2 = [];
for ii = 1:numel(NLI)
    E1 = EigenMatrix1{ii};
    [V1RGB,D1] = eig(E1); [~,P_id1] = sort(diag(1./D1)); % Sorted accordingly
    V1 = M*V1RGB(:,P_id1); % Converting from RGB to LMS
    V1 = V1./repmat(sum(abs(V1),1),[3 1]);
    
    E2 = EigenMatrix2{ii};
    [V2RGB,D2] = eig(E2); [~,P_id2] = sort(diag(1./D2)); % Sorted accordingly
    V2 = M*V2RGB(:,P_id2); % Converting from RGB to LMS
    V2 = V2./repmat(sum(abs(V2),1),[3 1]);
    
    % Determining the angle between the principal vector and the LMS cone weights
    tmp_angle1 = acos(dot(V1(:,3),S1LMS(:,ii))/(norm(V1(:,3))*norm(S1LMS(:,ii))))*180/pi;
    tmp_angle2 = acos(dot(V2(:,3),S2LMS(:,ii))/(norm(V2(:,3))*norm(S2LMS(:,ii))))*180/pi;
    
    Anglebw_Princ_LMS1 = [Anglebw_Princ_LMS1; min([tmp_angle1 180-tmp_angle1])];
    Anglebw_Princ_LMS2 = [Anglebw_Princ_LMS2; min([tmp_angle2 180-tmp_angle2])];
    
    % Plotting the vectors along with their cone weights
    ConeweightsS1 = [V1(:,1)/sum(abs(V1(:,1))) V1(:,2)/sum(abs(V1(:,2))) V1(:,3)/sum(abs(V1(:,3)))];
    ConeweightsS2 = [V2(:,1)/sum(abs(V2(:,1))) V2(:,2)/sum(abs(V2(:,2))) V2(:,3)/sum(abs(V2(:,3)))];
    
    if any(ismember(nonlin_cells,ii))
        subplot(221); plot3([0; V1(1,3)],[0; V1(2,3)],[0; V1(3,3)],'Color',[0 0 0],'LineWidth',2); hold on;
        plot3([0; V2(1,3)],[0; V2(2,3)],[0; V2(3,3)],'Color',[0 0 0],'LineWidth',2); hold on;
        
        subplot(224); plot(ConeweightsS1(1,3),ConeweightsS1(2,3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
        plot(ConeweightsS2(1,3),ConeweightsS2(2,3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
        
        for jj = 1:2
            if D1(P_id1(jj),P_id1(jj))>0
                subplot(222); plot3([0; V1(1,jj)],[0; V1(2,jj)],[0; V1(3,jj)],'Color',[0 1 0],'LineWidth',2); hold on;
                subplot(224); plot(ConeweightsS1(1,jj),ConeweightsS1(2,jj),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
            else
                subplot(223); plot3([0; V1(1,jj)],[0; V1(2,jj)],[0; V1(3,jj)],'Color',[0 1 1],'LineWidth',2); hold on;
                subplot(224); plot(ConeweightsS1(1,jj),ConeweightsS1(2,jj),'o','MarkerFaceColor',[0 1 1],'MarkerEdgeColor',[1 1 1]); hold on;
            end
            
            if D2(P_id2(jj),P_id2(jj))>0
                subplot(222); plot3([0; V2(1,jj)],[0; V2(2,jj)],[0; V2(3,jj)],'Color',[0 1 0],'LineWidth',2); hold on;
                subplot(224); plot(ConeweightsS2(1,jj),ConeweightsS2(2,jj),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
            else
                subplot(223); plot3([0; V2(1,jj)],[0; V2(2,jj)],[0; V2(3,jj)],'Color',[0 1 1],'LineWidth',2); hold on;
                subplot(224); plot(ConeweightsS2(1,jj),ConeweightsS2(2,jj),'o','MarkerFaceColor',[0 1 1],'MarkerEdgeColor',[1 1 1]); hold on;
            end
        end
    end

end

hold on; subplot(221); line([-1 1],[0 0],[0 0],'Linewidth',3); line([0 0],[-1 1],[0 0],'Linewidth',3); line([0 0],[0 0],[-1 1],'Linewidth',3);
axis equal; axis xy; grid on; xlabel('L'),ylabel('M'),zlabel('S'); title('Best response direction'); hold off;
hold on; subplot(222); line([-1 1],[0 0],[0 0],'Linewidth',3); line([0 0],[-1 1],[0 0],'Linewidth',3); line([0 0],[0 0],[-1 1],'Linewidth',3);
axis equal; axis xy; grid on; xlabel('L'),ylabel('M'),zlabel('S'); title('Other response direction'); hold off;
hold on; subplot(223); line([-1 1],[0 0],[0 0],'Linewidth',3); line([0 0],[-1 1],[0 0],'Linewidth',3); line([0 0],[0 0],[-1 1],'Linewidth',3);
axis equal; axis xy; grid on; xlabel('L'),ylabel('M'),zlabel('S'); title('Worst response direction'); hold off;
hold on; subplot(224); axis equal; set(gca,'Xlim',[-1 1],'Ylim',[-1 1],'XTick',-1:0.5:1,'YTick',-1:0.5:1,'Tickdir','out'); 
plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 -1],[0 0],'k'); plot([1 0],[0 -1],'k'); plot([0 -1],[-1 0],'k');  xlabel('L'), ylabel('M');
plot_counter = plot_counter + 1;

% Running a monte carlo simulation to see the distribution of random RGB points in LMS space 
N = 10000;
LMS_pts = []; RGB_pts = [];
for ii = 1:N
    pts_rgb = randn(3,1); 
    pt = M*pts_rgb; % Converting random RGB pts to LMS 
    LMS_pts = [LMS_pts pt./sum(abs(pt))];
    RGB_pts = [RGB_pts pts_rgb./sum(abs(pts_rgb))];
end

figure(plot_counter); set(gcf,'Name','Monte Carlo RGB points in LMS space'); 
subplot(121); plot(RGB_pts(1,:),RGB_pts(2,:),'k.'); hold on;
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[-1 1],'XTick',-1:0.5:1,'YTick',-1:0.5:1,'Tickdir','out'); 
plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 -1],[0 0],'k'); plot([1 0],[0 -1],'k'); plot([0 -1],[-1 0],'k');  xlabel('R'), ylabel('G');
subplot(122); plot(LMS_pts(1,:),LMS_pts(2,:),'k.'); hold on;
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[-1 1],'XTick',-1:0.5:1,'YTick',-1:0.5:1,'Tickdir','out'); 
plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 -1],[0 0],'k'); plot([1 0],[0 -1],'k'); plot([0 -1],[-1 0],'k');  xlabel('L'), ylabel('M');
plot_counter = plot_counter + 1;
