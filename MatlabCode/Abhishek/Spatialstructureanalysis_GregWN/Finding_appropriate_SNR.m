% Here, I am checking what would be a good definition of SNR to use 
% Author - Abhishek De, 12/19

close all; clearvars;
%% Approach 1: Plotting R square as a function of SNR (energy of peak frame/energy of noise frame)
% Using the spatial SNR: % See PopAnalysis_createOutputList.m for what is being stored in the 20th column
if ~exist('plot_counter')
    plot_counter = 1;
end
% load Output_ListWN.mat
load Output_ListWN2.mat
load Singleopponent.mat
load modelfits.mat
crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;

Output_List(~Z_cellsofinterest,:) = [];
NLI = cell2mat(Output_List(:,13));
simplecells = NLI<0;
Singleopponent = logical(Singleopponent);
% SNR = Zmax(Z_cellsofinterest)./Z(Z_cellsofinterest,1);
SNR = cell2mat(Output_List(:,20)); 
SNR = SNR(~Singleopponent & simplecells);
% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

RGB_svd = cell2mat(Output_List(~Singleopponent & simplecells,5)');
bkgndlms = cell2mat(Output_List(1,16)');
Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
Mrgbtocc = inv(Mrgbtocc');
conewts_svd = Mrgbtocc * RGB_svd;
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];

load Peraccuracy.mat
load SSE.mat
load Deviation.mat
meanperaccuracy = zeros(size(Peraccuracy));
meanSSE = zeros(size(SSE));
meanR = zeros(size(SSE));
for ii = 1:size(Peraccuracy,1)
    for jj = 1:size(Peraccuracy,2)
        meanperaccuracy(ii,jj) = mean(Peraccuracy{ii,jj});
        meanSSE(ii,jj) = mean(SSE{ii,jj});
        meanR(ii,jj) = mean(cos(Deviation{ii,jj}*pi/180));
    end
end
meanperaccuracy = meanperaccuracy(~Singleopponent & simplecells,:);
meanSSE = meanSSE(~Singleopponent & simplecells,:);
meanR = meanR(~Singleopponent & simplecells,:);

bins = 0:0.1:1.2;
figure(plot_counter); set(gcf,'Name','R=f(SNR): Just the spatial SNR: weighted STA');
subplot(221);  hold on;
plot(log10(SNR(ColorOpponentIds_conewts)),meanR(ColorOpponentIds_conewts,2),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(log10(SNR(Sconedominated_conewts)),meanR(Sconedominated_conewts,2),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]);
plot(log10(SNR(LumIds_conewts)),meanR(LumIds_conewts,2),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[1.8 3.4],'Ylim',[0 1],'XTick',1.8:0.4:3.4,'YTick',0:0.25:1.0); axis square; xlabel('log10(SNR)'); ylabel('R'); title('Gabor'); hold off;
subplot(222); hold on;
plot(log10(SNR(ColorOpponentIds_conewts)),meanR(ColorOpponentIds_conewts,3),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(log10(SNR(Sconedominated_conewts)),meanR(Sconedominated_conewts,3),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]);
 plot(log10(SNR(LumIds_conewts)),meanR(LumIds_conewts,3),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[1.8 3.4],'Ylim',[0 1],'XTick',1.8:0.4:3.4,'YTick',0:0.25:1.0); axis square; xlabel('log10(SNR)'); ylabel('R'); title('DoG'); hold off;

% Kruskal Wallis test for comparing SNR between groups 
data1 = [log10(SNR(LumIds_conewts)); log10(SNR(ColorOpponentIds_conewts)); log10(SNR(Sconedominated_conewts))];
data2 = [meanR(LumIds_conewts,2); meanR(ColorOpponentIds_conewts,2); meanR(Sconedominated_conewts,2)]; % Gabor 
data3 = [meanR(LumIds_conewts,3); meanR(ColorOpponentIds_conewts,3); meanR(Sconedominated_conewts,3)]; % DoG
group = [ones(size(log10(SNR(LumIds_conewts)))); 2*ones(size(log10(SNR(ColorOpponentIds_conewts)))); 3*ones(size(log10(SNR(Sconedominated_conewts))))];
p1 = kruskalwallis(data1,group,'off');
p2 = kruskalwallis(data2,group,'off');
p3 = kruskalwallis(data3,group,'off');
idxs = data1>median(data1);
p4 = kruskalwallis(data3(idxs),group(idxs),'off');
p5 = kruskalwallis(data2(idxs),group(idxs),'off');

% Fitting a saturating exponential to the R_Gabor 
f = @(b,x) b(3)./(1+exp(b(1)*(x-b(2)))); % Objective Function
B1 = lsqcurvefit(f,[-1; 2; 0.5],log10(SNR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts])),meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],2),[-10 -5 eps],[-eps 5 1]);
B2 = lsqcurvefit(f,[-1; 2; 0.5],log10(SNR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts])),meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],3),[-10 -5 eps],[-eps 5 1]);

figure(plot_counter); 
subplot(223); plot(log10(SNR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts])),meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],2),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
hold on; plot(linspace(1.8,3.4,51),f(B1,linspace(1.8,3.4,51)),'k','Linewidth',2);
set(gca,'Tickdir','out','Xlim',[1.8 3.4],'Ylim',[0 1],'XTick',1.8:0.4:3.4,'YTick',0:0.25:1.0); axis square; xlabel('log10(SNR)'); ylabel('R'); title('Gabor'); hold off;
subplot(224); plot(log10(SNR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts])),meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
hold on; plot(linspace(1.8,3.4,51),f(B2,linspace(1.8,3.4,51)),'k','Linewidth',2);
set(gca,'Tickdir','out','Xlim',[1.8 3.4],'Ylim',[0 1],'XTick',1.8:0.4:3.4,'YTick',0:0.25:1.0); axis square; xlabel('log10(SNR)'); ylabel('R'); title('Gabor'); hold off;
plot_counter = plot_counter + 1;

Max_GaborR = B1(3);
Max_DoGR = B2(3);

%% Approach 2: Plotting R square as a function of SNR (energy of peak frame/energy of noise frame)
% Using the complete SNR
if ~exist('plot_counter')
    plot_counter = 1;
end
% load Output_ListWN.mat
load Output_ListWN2.mat
load Singleopponent.mat
load modelfits.mat
crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;

Output_List(~Z_cellsofinterest,:) = [];
NLI = cell2mat(Output_List(:,13));
simplecells = NLI<0;
Singleopponent = logical(Singleopponent);
SNR = Zmax(Z_cellsofinterest);
% SNR = cell2mat(Output_List(:,20)); 
SNR = SNR(~Singleopponent & simplecells);
% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

RGB_svd = cell2mat(Output_List(~Singleopponent & simplecells,5)');
bkgndlms = cell2mat(Output_List(1,16)');
Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
Mrgbtocc = inv(Mrgbtocc');
conewts_svd = Mrgbtocc * RGB_svd;
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];

load Peraccuracy.mat
load SSE.mat
load Deviation.mat
meanperaccuracy = zeros(size(Peraccuracy));
meanSSE = zeros(size(SSE));
meanR = zeros(size(SSE));
for ii = 1:size(Peraccuracy,1)
    for jj = 1:size(Peraccuracy,2)
        meanperaccuracy(ii,jj) = mean(Peraccuracy{ii,jj});
        meanSSE(ii,jj) = mean(SSE{ii,jj});
        meanR(ii,jj) = mean(cos(Deviation{ii,jj}*pi/180));
    end
end
meanperaccuracy = meanperaccuracy(~Singleopponent & simplecells,:);
meanSSE = meanSSE(~Singleopponent & simplecells,:);
meanR = meanR(~Singleopponent & simplecells,:);

bins = 0:0.1:1.2;
figure(plot_counter); set(gcf,'Name','R=f(SNR): Using complete SNR, peak STA frame');
subplot(221);  hold on;
plot(log10(SNR(ColorOpponentIds_conewts)),meanR(ColorOpponentIds_conewts,2),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(log10(SNR(Sconedominated_conewts)),meanR(Sconedominated_conewts,2),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]);
plot(log10(SNR(LumIds_conewts)),meanR(LumIds_conewts,2),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[2.4 3.6],'Ylim',[0 1],'XTick',2.4:0.4:3.6,'YTick',0:0.25:1.0); axis square; xlabel('log10(SNR)'); ylabel('R'); title('Gabor'); hold off;
subplot(222); hold on;
plot(log10(SNR(ColorOpponentIds_conewts)),meanR(ColorOpponentIds_conewts,3),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(log10(SNR(Sconedominated_conewts)),meanR(Sconedominated_conewts,3),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]);
plot(log10(SNR(LumIds_conewts)),meanR(LumIds_conewts,3),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); 
set(gca,'Tickdir','out','Xlim',[2.4 3.6],'Ylim',[0 1],'XTick',2.4:0.4:3.6,'YTick',0:0.25:1.0); axis square; xlabel('log10(SNR)'); ylabel('R'); title('DoG'); hold off;

% Kruskal Wallis test for comparing SNR between groups 
data1 = [log10(SNR(LumIds_conewts)); log10(SNR(ColorOpponentIds_conewts)); log10(SNR(Sconedominated_conewts))];
data2 = [meanR(LumIds_conewts,2); meanR(ColorOpponentIds_conewts,2); meanR(Sconedominated_conewts,2)]; % Gabor 
data3 = [meanR(LumIds_conewts,3); meanR(ColorOpponentIds_conewts,3); meanR(Sconedominated_conewts,3)]; % DoG
group = [ones(size(log10(SNR(LumIds_conewts)))); 2*ones(size(log10(SNR(ColorOpponentIds_conewts)))); 3*ones(size(log10(SNR(Sconedominated_conewts))))];
p1 = kruskalwallis(data1,group,'off');
p2 = kruskalwallis(data2,group,'off');
p3 = kruskalwallis(data3,group,'off');
idxs = data1>median(data1);
p4 = kruskalwallis(data3(idxs),group(idxs),'off');
p5 = kruskalwallis(data2(idxs),group(idxs),'off');

% Fitting a saturating exponential to the R_Gabor 
f = @(b,x) b(3)./(1+exp(b(1)*(x-b(2)))); % Objective Function
B1 = lsqcurvefit(f,[-1; 2; 0.5],log10(SNR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts])),meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],2),[-10 -5 eps],[-eps 5 1]);
B2 = lsqcurvefit(f,[-1; 2; 0.5],log10(SNR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts])),meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],3),[-10 -5 eps],[-eps 5 1]);

figure(plot_counter); 
subplot(223); plot(log10(SNR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts])),meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],2),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
hold on; plot(linspace(2.4,3.6,51),f(B1,linspace(2.4,3.6,51)),'k','Linewidth',2);
set(gca,'Tickdir','out','Xlim',[2.4 3.6],'Ylim',[0 1],'XTick',2.4:0.4:3.6,'YTick',0:0.25:1.0); axis square; xlabel('log10(SNR)'); ylabel('R'); title('Gabor'); hold off;
subplot(224); plot(log10(SNR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts])),meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
hold on; plot(linspace(2.4,3.6,51),f(B2,linspace(2.4,3.6,51)),'k','Linewidth',2);
set(gca,'Tickdir','out','Xlim',[2.4 3.6],'Ylim',[0 1],'XTick',2.4:0.4:3.6,'YTick',0:0.25:1.0); axis square; xlabel('log10(SNR)'); ylabel('R'); title('Gabor'); hold off;
plot_counter = plot_counter + 1;

Max_GaborR = B1(3);
Max_DoGR = B2(3);


%% Approach 3: Plotting R square as a function of SNR (energy of peak frame/energy of noise frame)
% Using SNR from the weighted STA and not just the peak frame: % See PopAnalysis_createOutputList.m for what is being stored in the 22nd column
if ~exist('plot_counter')
    plot_counter = 1;
end
% load Output_ListWN.mat
load Output_ListWN2.mat
load Singleopponent.mat
load modelfits.mat
crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;

Output_List(~Z_cellsofinterest,:) = [];
NLI = cell2mat(Output_List(:,13));
simplecells = NLI<0;
Singleopponent = logical(Singleopponent);
% SNR = Zmax(Z_cellsofinterest);
SNR = cell2mat(Output_List(:,22)); 
SNR = SNR(~Singleopponent & simplecells);
% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

RGB_svd = cell2mat(Output_List(~Singleopponent & simplecells,5)');
bkgndlms = cell2mat(Output_List(1,16)');
Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
Mrgbtocc = inv(Mrgbtocc');
conewts_svd = Mrgbtocc * RGB_svd;
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];

load Peraccuracy.mat
load SSE.mat
load Deviation.mat
meanperaccuracy = zeros(size(Peraccuracy));
meanSSE = zeros(size(SSE));
meanR = zeros(size(SSE));
for ii = 1:size(Peraccuracy,1)
    for jj = 1:size(Peraccuracy,2)
        meanperaccuracy(ii,jj) = mean(Peraccuracy{ii,jj});
        meanSSE(ii,jj) = mean(SSE{ii,jj});
        meanR(ii,jj) = mean(cos(Deviation{ii,jj}*pi/180));
    end
end
meanperaccuracy = meanperaccuracy(~Singleopponent & simplecells,:);
meanSSE = meanSSE(~Singleopponent & simplecells,:);
meanR = meanR(~Singleopponent & simplecells,:);

bins = 0:0.1:1.2;
figure(plot_counter); set(gcf,'Name','R=f(SNR): Weighted SNR: weighted STA');
subplot(221);  hold on;
plot(log10(SNR(ColorOpponentIds_conewts)),meanR(ColorOpponentIds_conewts,2),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(log10(SNR(Sconedominated_conewts)),meanR(Sconedominated_conewts,2),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]);
plot(log10(SNR(LumIds_conewts)),meanR(LumIds_conewts,2),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[2.0 3.6],'Ylim',[0 1],'XTick',2.0:0.4:3.6,'YTick',0:0.25:1.0); axis square; xlabel('log10(SNR)'); ylabel('R'); title('Gabor'); hold off;
subplot(222); hold on;
plot(log10(SNR(ColorOpponentIds_conewts)),meanR(ColorOpponentIds_conewts,3),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(log10(SNR(Sconedominated_conewts)),meanR(Sconedominated_conewts,3),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]);
plot(log10(SNR(LumIds_conewts)),meanR(LumIds_conewts,3),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); 
set(gca,'Tickdir','out','Xlim',[2.0 3.6],'Ylim',[0 1],'XTick',2.0:0.4:3.6,'YTick',0:0.25:1.0); axis square; xlabel('log10(SNR)'); ylabel('R'); title('DoG'); hold off;

% Kruskal Wallis test for comparing SNR between groups 
data1 = [log10(SNR(LumIds_conewts)); log10(SNR(ColorOpponentIds_conewts)); log10(SNR(Sconedominated_conewts))];
data2 = [meanR(LumIds_conewts,2); meanR(ColorOpponentIds_conewts,2); meanR(Sconedominated_conewts,2)]; % Gabor 
data3 = [meanR(LumIds_conewts,3); meanR(ColorOpponentIds_conewts,3); meanR(Sconedominated_conewts,3)]; % DoG
group = [ones(size(log10(SNR(LumIds_conewts)))); 2*ones(size(log10(SNR(ColorOpponentIds_conewts)))); 3*ones(size(log10(SNR(Sconedominated_conewts))))];
p1 = kruskalwallis(data1,group,'off');
p2 = kruskalwallis(data2,group,'off');
p3 = kruskalwallis(data3,group,'off');
idxs = data1>median(data1);
p4 = kruskalwallis(data3(idxs),group(idxs),'off');
p5 = kruskalwallis(data2(idxs),group(idxs),'off');

% Fitting a saturating exponential to the R_Gabor 
f = @(b,x) b(3)./(1+exp(b(1)*(x-b(2)))); % Objective Function
B1 = lsqcurvefit(f,[-1; 2; 0.5],log10(SNR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts])),meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],2),[-10 -5 eps],[-eps 5 1]);
B2 = lsqcurvefit(f,[-1; 2; 0.5],log10(SNR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts])),meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],3),[-10 -5 eps],[-eps 5 1]);

figure(plot_counter); 
subplot(223); plot(log10(SNR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts])),meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],2),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
hold on; plot(linspace(2.0,3.6,51),f(B1,linspace(2.0,3.6,51)),'k','Linewidth',2);
set(gca,'Tickdir','out','Xlim',[2.0 3.6],'Ylim',[0 1],'XTick',2.0:0.4:3.6,'YTick',0:0.25:1.0); axis square; xlabel('log10(SNR)'); ylabel('R'); title('Gabor'); hold off;
subplot(224); plot(log10(SNR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts])),meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
hold on; plot(linspace(2.0,3.6,51),f(B2,linspace(2.0,3.6,51)),'k','Linewidth',2);
set(gca,'Tickdir','out','Xlim',[2.0 3.6],'Ylim',[0 1],'XTick',2.0:0.4:3.6,'YTick',0:0.25:1.0); axis square; xlabel('log10(SNR)'); ylabel('R'); title('Gabor'); hold off;
plot_counter = plot_counter + 1;

Max_GaborR = B1(3);
Max_DoGR = B2(3);

%% Approach 4: Plotting R square as a function of SNR (energy of peak frame/energy of noise frame)
% Using SNR as the ratio of first singular value of the weighted STA and the first singular value of the noise frame 
%See PopAnalysis_createOutputList.m for what is being stored in the 21st column
if ~exist('plot_counter')
    plot_counter = 1;
end
% load Output_ListWN.mat
load Output_ListWN2.mat
load Singleopponent.mat
load modelfits.mat
crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;

Output_List(~Z_cellsofinterest,:) = [];
NLI = cell2mat(Output_List(:,13));
simplecells = NLI<0;
Singleopponent = logical(Singleopponent);
% SNR = Zmax(Z_cellsofinterest);
N = cell2mat(Output_List(:,21)); % Noise 
S = cell2mat(Output_List(:,17)); % Signal
SNR = S(:,1)./N(:,1); 
SNR = SNR(~Singleopponent & simplecells);
% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

RGB_svd = cell2mat(Output_List(~Singleopponent & simplecells,5)');
bkgndlms = cell2mat(Output_List(1,16)');
Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
Mrgbtocc = inv(Mrgbtocc');
conewts_svd = Mrgbtocc * RGB_svd;
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];

load Peraccuracy.mat
load SSE.mat
load Deviation.mat
meanperaccuracy = zeros(size(Peraccuracy));
meanSSE = zeros(size(SSE));
meanR = zeros(size(SSE));
for ii = 1:size(Peraccuracy,1)
    for jj = 1:size(Peraccuracy,2)
        meanperaccuracy(ii,jj) = mean(Peraccuracy{ii,jj});
        meanSSE(ii,jj) = mean(SSE{ii,jj});
        meanR(ii,jj) = mean(cos(Deviation{ii,jj}*pi/180));
    end
end
meanperaccuracy = meanperaccuracy(~Singleopponent & simplecells,:);
meanSSE = meanSSE(~Singleopponent & simplecells,:);
meanR = meanR(~Singleopponent & simplecells,:);

bins = 0:0.1:1.2;
figure(plot_counter); set(gcf,'Name','R=f(SNR): Ratio of singular values, weighted STA');
subplot(221);  hold on;
plot(log10(SNR(ColorOpponentIds_conewts)),meanR(ColorOpponentIds_conewts,2),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(log10(SNR(Sconedominated_conewts)),meanR(Sconedominated_conewts,2),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]);
plot(log10(SNR(LumIds_conewts)),meanR(LumIds_conewts,2),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-0.3 0.7],'Ylim',[0 1],'XTick',-0.3:0.5:0.7,'YTick',0:0.25:1.0); axis square; xlabel('log10(SNR)'); ylabel('R'); title('Gabor'); hold off;
subplot(222); hold on;
plot(log10(SNR(ColorOpponentIds_conewts)),meanR(ColorOpponentIds_conewts,3),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(log10(SNR(Sconedominated_conewts)),meanR(Sconedominated_conewts,3),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]);
plot(log10(SNR(LumIds_conewts)),meanR(LumIds_conewts,3),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); 
set(gca,'Tickdir','out','Xlim',[-0.3 0.7],'Ylim',[0 1],'XTick',-0.3:0.5:0.7,'YTick',0:0.25:1.0); axis square; xlabel('log10(SNR)'); ylabel('R'); title('DoG'); hold off;

% Kruskal Wallis test for comparing SNR between groups 
data1 = [log10(SNR(LumIds_conewts)); log10(SNR(ColorOpponentIds_conewts)); log10(SNR(Sconedominated_conewts))];
data2 = [meanR(LumIds_conewts,2); meanR(ColorOpponentIds_conewts,2); meanR(Sconedominated_conewts,2)]; % Gabor 
data3 = [meanR(LumIds_conewts,3); meanR(ColorOpponentIds_conewts,3); meanR(Sconedominated_conewts,3)]; % DoG
group = [ones(size(log10(SNR(LumIds_conewts)))); 2*ones(size(log10(SNR(ColorOpponentIds_conewts)))); 3*ones(size(log10(SNR(Sconedominated_conewts))))];
p1 = kruskalwallis(data1,group,'off');
p2 = kruskalwallis(data2,group,'off');
p3 = kruskalwallis(data3,group,'off');
idxs = data1>median(data1);
p4 = kruskalwallis(data3(idxs),group(idxs),'off');
p5 = kruskalwallis(data2(idxs),group(idxs),'off');

% Fitting a saturating exponential to the R_Gabor 
f = @(b,x) b(3)./(1+exp(b(1)*(x-b(2)))); % Objective Function
B1 = lsqcurvefit(f,[-1; 2; 0.5],log10(SNR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts])),meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],2),[-10 -5 eps],[-eps 5 1]);
B2 = lsqcurvefit(f,[-1; 2; 0.5],log10(SNR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts])),meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],3),[-10 -5 eps],[-eps 5 1]);

figure(plot_counter); 
subplot(223); plot(log10(SNR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts])),meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],2),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
hold on; plot(linspace(-0.3,0.7,51),f(B1,linspace(-0.3,0.7,51)),'k','Linewidth',2);
set(gca,'Tickdir','out','Xlim',[-0.3 0.7],'Ylim',[0 1],'XTick',-0.3:0.5:0.7,'YTick',0:0.25:1.0); axis square; xlabel('log10(SNR)'); ylabel('R'); title('Gabor'); hold off;
subplot(224); plot(log10(SNR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts])),meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
hold on; plot(linspace(-0.3,0.7,51),f(B2,linspace(-0.3,0.7,51)),'k','Linewidth',2);
set(gca,'Tickdir','out','Xlim',[-0.3 0.7],'Ylim',[0 1],'XTick',-0.3:0.5:0.7,'YTick',0:0.25:1.0); axis square; xlabel('log10(SNR)'); ylabel('R'); title('Gabor'); hold off;
plot_counter = plot_counter + 1;

Max_GaborR = B1(3);
Max_DoGR = B2(3);


%% Approach 5: Plotting R square as a function of SNR (energy of peak frame/energy of noise frame)
% Defining SNR as z instead of z2
if ~exist('plot_counter')
    plot_counter = 1;
end
% load Output_ListWN.mat
load Output_ListWN2.mat
load Singleopponent.mat
load modelfits.mat
crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;

Output_List(~Z_cellsofinterest,:) = [];
NLI = cell2mat(Output_List(:,13));
simplecells = NLI<0;
Singleopponent = logical(Singleopponent);
SNR = sqrt(Zmax(Z_cellsofinterest));
% N = cell2mat(Output_List(:,21)); % Noise 
% S = cell2mat(Output_List(:,17)); % Signal
% SNR = S(:,1)./N(:,1); 
SNR = SNR(~Singleopponent & simplecells);
% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

RGB_svd = cell2mat(Output_List(~Singleopponent & simplecells,5)');
bkgndlms = cell2mat(Output_List(1,16)');
Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
Mrgbtocc = inv(Mrgbtocc');
conewts_svd = Mrgbtocc * RGB_svd;
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];

load Peraccuracy.mat
load SSE.mat
load Deviation.mat
meanperaccuracy = zeros(size(Peraccuracy));
meanSSE = zeros(size(SSE));
meanR = zeros(size(SSE));
for ii = 1:size(Peraccuracy,1)
    for jj = 1:size(Peraccuracy,2)
        meanperaccuracy(ii,jj) = mean(Peraccuracy{ii,jj});
        meanSSE(ii,jj) = mean(SSE{ii,jj});
        meanR(ii,jj) = mean(cos(Deviation{ii,jj}*pi/180));
    end
end
meanperaccuracy = meanperaccuracy(~Singleopponent & simplecells,:);
meanSSE = meanSSE(~Singleopponent & simplecells,:);
meanR = meanR(~Singleopponent & simplecells,:);

bins = 0:0.1:1.2;
figure(plot_counter); set(gcf,'Name','R=f(SNR):peak STA frame but using z');
subplot(221);  hold on;
plot(log10(SNR(ColorOpponentIds_conewts)),meanR(ColorOpponentIds_conewts,2),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(log10(SNR(Sconedominated_conewts)),meanR(Sconedominated_conewts,2),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]);
plot(log10(SNR(LumIds_conewts)),meanR(LumIds_conewts,2),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[1.2 1.8],'Ylim',[0 1],'XTick',1.2:0.3:1.8,'YTick',0:0.25:1.0); axis square; xlabel('log10(SNR)'); ylabel('R'); title('Gabor'); hold off;
subplot(222); hold on;
plot(log10(SNR(ColorOpponentIds_conewts)),meanR(ColorOpponentIds_conewts,3),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(log10(SNR(Sconedominated_conewts)),meanR(Sconedominated_conewts,3),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]);
plot(log10(SNR(LumIds_conewts)),meanR(LumIds_conewts,3),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); 
set(gca,'Tickdir','out','Xlim',[1.2 1.8],'Ylim',[0 1],'XTick',1.2:0.3:1.8,'YTick',0:0.25:1.0); axis square; xlabel('log10(SNR)'); ylabel('R'); title('DoG'); hold off;

% Kruskal Wallis test for comparing SNR between groups 
data1 = [log10(SNR(LumIds_conewts)); log10(SNR(ColorOpponentIds_conewts)); log10(SNR(Sconedominated_conewts))];
data2 = [meanR(LumIds_conewts,2); meanR(ColorOpponentIds_conewts,2); meanR(Sconedominated_conewts,2)]; % Gabor 
data3 = [meanR(LumIds_conewts,3); meanR(ColorOpponentIds_conewts,3); meanR(Sconedominated_conewts,3)]; % DoG
group = [ones(size(log10(SNR(LumIds_conewts)))); 2*ones(size(log10(SNR(ColorOpponentIds_conewts)))); 3*ones(size(log10(SNR(Sconedominated_conewts))))];
p1 = kruskalwallis(data1,group,'off');
p2 = kruskalwallis(data2,group,'off');
p3 = kruskalwallis(data3,group,'off');
idxs = data1>median(data1);
p4 = kruskalwallis(data3(idxs),group(idxs),'off');
p5 = kruskalwallis(data2(idxs),group(idxs),'off');

% Fitting a saturating exponential to the R_Gabor 
f = @(b,x) b(3)./(1+exp(b(1)*(x-b(2)))); % Objective Function
B1 = lsqcurvefit(f,[-1; 2; 0.5],log10(SNR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts])),meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],2),[-10 -5 eps],[-eps 5 1]);
B2 = lsqcurvefit(f,[-1; 2; 0.5],log10(SNR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts])),meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],3),[-10 -5 eps],[-eps 5 1]);

figure(plot_counter); 
subplot(223); plot(log10(SNR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts])),meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],2),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
hold on; plot(linspace(1.2,1.8,51),f(B1,linspace(1.2,1.8,51)),'k','Linewidth',2);
set(gca,'Tickdir','out','Xlim',[1.2 1.8],'Ylim',[0 1],'XTick',1.2:0.3:1.8,'YTick',0:0.25:1.0); axis square; xlabel('log10(SNR)'); ylabel('R'); title('Gabor'); hold off;
subplot(224); plot(log10(SNR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts])),meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
hold on; plot(linspace(1.2,1.8,51),f(B2,linspace(1.2,1.8,51)),'k','Linewidth',2);
set(gca,'Tickdir','out','Xlim',[1.2 1.8],'Ylim',[0 1],'XTick',1.2:0.3:1.8,'YTick',0:0.25:1.0); axis square; xlabel('log10(SNR)'); ylabel('R'); title('Gabor'); hold off;
plot_counter = plot_counter + 1;

Max_GaborR = B1(3);
Max_DoGR = B2(3);

