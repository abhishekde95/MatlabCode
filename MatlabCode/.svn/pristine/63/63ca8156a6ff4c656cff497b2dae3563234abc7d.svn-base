%% Script for David Perkel's Neubeh 503 4/27 Homework

clear all
close all
load GaleSampleData.mat

% User-Defined Variables
bins = 50;

% Organize Data and Set Auto-Defined Variables
nonan  = ~isnan(GaleSampleData(1,:));
GSD    = GaleSampleData(:,nonan);
Dprime = nan(1,12);
TPR    = nan(12,bins);
FPR    = nan(12,bins);
dist   = nan(12,bins);
AUC    = nan(1,12);

%% Question 1

% Calculate mean and variance of RS
meanRS = nanmean(GSD);
varRS  = nanvar(GSD);
stdRS  = nanstd(GSD);
n      = sum(~isnan(GSD));
sem    = stdRS./sqrt(n);

% Calculate D'
for n = 1:6
    Dprime(1,n*2-1) = (meanRS(n*3-2)-meanRS(n*3-1))/sqrt(mean([varRS(n*3-2) varRS(n*3-1)]));
    Dprime(1,n*2) = (meanRS(n*3-2)-meanRS(n*3))/sqrt(mean([varRS(n*3-2) varRS(n*3)]));
end

meanRS
Dprime


%% Question 2

% Generate ROC Curves
% 4/28 hw discussion: Using fixed intervals throws away data.  Using each
% datapoint as a new threshold maximizes the data.
for q=1:6
    dist(q*2-1,:) = linspace(min(min(GSD(:,q*3-2:q*3-1))),max(max(GSD(:,q*3-2:q*3-1))),bins);
    for n=1:bins
        TPR(q*2-1,n) = sum(GSD(:,q*3-2) >= dist(q*2-1,n))/sum(~isnan(GSD(:,q*3-2)));
        FPR(q*2-1,n) = sum(GSD(:,q*3-1) >= dist(q*2-1,n))/sum(~isnan(GSD(:,q*3-1)));
    end
    dist(q*2,:) = linspace(min(min(GSD(:,[q*3-2,q*3]))),max(max(GSD(:,[q*3-2,q*3]))),bins);
    for n=1:bins
        TPR(q*2,n) = sum(GSD(:,q*3-2) >= dist(q*2,n))/sum(~isnan(GSD(:,q*3-2)));
        FPR(q*2,n) = sum(GSD(:,q*3) >= dist(q*2,n))/sum(~isnan(GSD(:,q*3)));
    end
    figure(1);hold on; grid on;
    plot(FPR(q*2-1,:),TPR(q*2-1,:),'o-b')
    plot(FPR(q*2,:),TPR(q*2,:),'*-g')    
    title('ROC Curves for BOS vs. revBOS and BOS vs. Conspecific')
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')
    legend('BOS vs. revBOS','BOS vs. conspecific','location','SouthOutside')
    axis([-.1 1.1 -.1 1.1])
end

% 4/27 added ROC, still needs work.
for q=1:6
    allpoints = unique([GSD(:,q*3-2) GSD(:,q*3-1)]);
    thresholds = allpoints(~isnan(allpoints));
    TPR1 = nan(1,length(thresholds));
    FPR1 = TPR1;
    for n=1:length(thresholds)
        thresh = thresholds(n);
        TPR1(n) = sum(GSD(:,q*3-2) >= thresh)/sum(~isnan(GSD(:,q*3-2)));
        FPR1(n) = sum(GSD(:,q*3-1) >= thresh)/sum(~isnan(GSD(:,q*3-1)));
    end
    allpoints = unique([GSD(:,q*3-2) GSD(:,q*3)]);
    thresholds = allpoints(~isnan(allpoints));
    TPR2 = nan(1,length(thresholds));
    FPR2 = TPR2;    
    for n=1:length(thresholds)
        thresh = thresholds(n);
        TPR2(n) = sum(GSD(:,q*3-2) >= thresh)/sum(~isnan(GSD(:,q*3-2)));
        FPR2(n) = sum(GSD(:,q*3) >= thresh)/sum(~isnan(GSD(:,q*3)));
    end
    figure(2);hold on; grid on;
    plot(FPR1,TPR1,'o-b')
    plot(FPR2,TPR2,'*-g')    
    title('ROC Curves for BOS vs. revBOS and BOS vs. Conspecific')
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')
    legend('BOS vs. revBOS','BOS vs. conspecific','location','SouthOutside')
    axis([-.1 1.1 -.1 1.1])
end


% Find area under the curves
FPR = fliplr(FPR);
TPR = fliplr(TPR);
for n = 1:size(FPR,1)
    AUC(n) = trapz(FPR(n,:),TPR(n,:));
end

AUC
    
%% Quesiton 3

% Make Scatter Plot
figure(3); clf; hold on; grid on;
scatter(AUC,Dprime)
title('Dprime vs. AUC')
xlabel('AUC')
ylabel('Dprime')
