% Figures for Jude's presentation.
% Author - Abhishek De, 2/20,


%% Figure 1A: Cone weights
close all; clearvars;

if ~exist('plot_counter')
    plot_counter = 1;
end
load Output_List_waveform.mat
load Singleopponent_waveform.mat

crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
Singleopponent_waveform(~Z_cellsofinterest) = [];
NLI = cell2mat(Output_List(:,13));
simplecells = NLI<0;
Singleopponent = Singleopponent_waveform;
Singleopponent = logical(Singleopponent);
% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

RGB_svd = cell2mat(Output_List(simplecells,5)');
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
DO_conewts = [ColorOpponentIds_conewts Sconedominated_conewts];

% First thing I want to plot are the cone weights
figure(plot_counter); subplot(321); hold on;
for ii = 1:numel(LumIds_conewts)
    h(ii) = plot(conewts_svd(1,LumIds_conewts(ii)),conewts_svd(2,LumIds_conewts(ii)),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
    set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
end
for ii = 1:numel(ColorOpponentIds_conewts)
    h(ii) = plot(conewts_svd(1,ColorOpponentIds_conewts(ii)),conewts_svd(2,ColorOpponentIds_conewts(ii)),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
    set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
end
for ii = 1:numel(Sconedominated_conewts)
    h(ii) = plot(conewts_svd(1,Sconedominated_conewts(ii)),conewts_svd(2,Sconedominated_conewts(ii)),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]); hold on;
    set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
end
for ii = 1:numel(Other_conewts)
    h(ii) = plot(conewts_svd(1,Other_conewts(ii)),conewts_svd(2,Other_conewts(ii)),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
    set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
end
plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out');  plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k');  xlabel('L'), ylabel('M'); title('Spatially non-opponent');
set(gcf,'renderer','painters');


%% Figure 1B: Waveform shapes 
if ~exist('plot_counter')
    plot_counter = 1;
end
load Output_List_waveform.mat
load Singleopponent_waveform.mat
load timediff.mat
load peak_waveform.mat
load trough_waveform.mat
load all_waveforms.mat
load burst_index.mat

crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
timediff(~Z_cellsofinterest) = [];
peak_waveform(~Z_cellsofinterest) = [];
trough_waveform(~Z_cellsofinterest) = [];
all_waveforms(~Z_cellsofinterest) = [];
burst_index(~Z_cellsofinterest) = [];
baselineFR = cell2mat(Output_List(:,18));

Singleopponent_waveform(~Z_cellsofinterest) = [];
NLI = cell2mat(Output_List(:,13));
simplecells = NLI<0;
Singleopponent = Singleopponent_waveform;
Singleopponent = logical(Singleopponent);
% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

% Using k-means to divide the data into two groups:
idx = kmeans(timediff,2);
[~,I1] = max([sum(idx==1) sum(idx==2)]); % narrow spiking 
[~,I2] = min([sum(idx==1) sum(idx==2)]); % broad spiking
narrow_wf = find(idx==I1);
broad_wf = find(idx==I2);

narrow_waveform = nan(numel(narrow_wf),1000); count_narrow = 1;
broad_waveform = nan(numel(broad_wf),1000); count_broad = 1;
indices = 1:10:size(narrow_waveform,2);
t = 2.5*(0:1:999)-500;
figure(plot_counter); subplot(322);
for ii = 1:size(all_waveforms,1)
    L = size(all_waveforms{ii},2);
    inds = 1:numel(all_waveforms{ii});
    [~,jj] = min(all_waveforms{ii});
    inds = inds - jj;
    W = all_waveforms{ii};
    W = W - ((min(W)+max(W))/2);
    W = W./max(W);
    if any(ismember(narrow_wf,ii))
        narrow_waveform(count_narrow,200-jj+1:200+L-jj) = W;
        subplot(322); plot(t(indices),narrow_waveform(count_narrow,indices),'Color','c'); hold on;
        count_narrow = count_narrow + 1;
    elseif any(ismember(broad_wf,ii))
        broad_waveform(count_broad,200-jj+1:200+L-jj) = W;
        subplot(322); plot(t(indices),broad_waveform(count_broad,indices),'Color','k'); hold on;
        count_broad = count_broad + 1;
    end
end
subplot(322); set(gca,'Tickdir','out','Xlim',[-200 600],'Ylim',[-1.1 1.1],'YTick',-1:1:1); axis square; 
xlabel('time (us)'); ylabel('Normalized Voltage'); hold off;


%% Figure 1C & 1D
if ~exist('plot_counter')
    plot_counter = 1;
end
load Output_List_waveform.mat
load Singleopponent_waveform.mat
load timediff.mat
load peak_waveform.mat
load trough_waveform.mat
load all_waveforms.mat
load burst_index.mat

crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
timediff(~Z_cellsofinterest) = [];
peak_waveform(~Z_cellsofinterest) = [];
trough_waveform(~Z_cellsofinterest) = [];
all_waveforms(~Z_cellsofinterest) = [];
burst_index(~Z_cellsofinterest) = [];
baselineFR = cell2mat(Output_List(:,18));
Singleopponent_waveform(~Z_cellsofinterest) = [];
NLI = cell2mat(Output_List(:,13));
simplecells = NLI<0;
Singleopponent = Singleopponent_waveform;
Singleopponent = logical(Singleopponent);
% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

RGB_svd = cell2mat(Output_List(simplecells,5)');
bkgndlms = cell2mat(Output_List(1,16)');
Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
Mrgbtocc = inv(Mrgbtocc');
conewts_svd = Mrgbtocc * RGB_svd;
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);
RGB_svdSO = cell2mat(Output_List(Singleopponent & simplecells,5)');
conewts_svdSO = Mrgbtocc * RGB_svdSO;
conewts_svdSO = conewts_svdSO./repmat(sum(abs(conewts_svdSO),1),[3 1]);
conewts_svdSO = conewts_svdSO .* repmat(sign(conewts_svdSO(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

% Next, I am plotting the SVD derived luminance, color-opponent and other cells
ind = find(simplecells);
LUMidx = ind(LumIds_conewts);
COidx = ind(ColorOpponentIds_conewts);
Sconeidx = ind(Sconedominated_conewts);
hardtoclassifyidx = [ind(Other_conewts); find(~simplecells)]; % unclassified and complex cells

% Plotting the histograms of aspect ratios and waveform widths
figure(plot_counter);
subplot(323); histogram(timediff,10:30:700,'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth',2); hold on;
histogram(timediff(LUMidx),10:30:700,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[0 600],'Ylim',[0 70],'XTick',0:200:600,'YTick',0:35:70); ylabel('#Lum cells'); axis square;

subplot(324); histogram(timediff,10:30:700,'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth',2); hold on;
histogram(timediff(COidx),10:30:700,'FaceColor',[1 0 0],'EdgeColor',[1 1 1]);
histogram(timediff(Sconeidx),10:30:700,'DisplayStyle','stairs','EdgeColor',[0 0.5 1.0],'Linewidth',2); 
set(gca,'Tickdir','out','Xlim',[0 600],'Ylim',[0 70],'XTick',0:200:600,'YTick',0:35:70); ylabel('#Color cells'); axis square;

% Figure 1E & 1F

if ~exist('plot_counter')
    plot_counter = 1;
end
load Output_List_waveform.mat
load Singleopponent_waveform.mat
load timediff.mat
load all_waveforms.mat
load burst_index.mat
load Propbursts_returnmap.mat

crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
timediff(~Z_cellsofinterest) = [];
all_waveforms(~Z_cellsofinterest) = [];
burst_index(~Z_cellsofinterest) = [];
Propbursts_returnmap(~Z_cellsofinterest) = [];

figure(plot_counter); subplot(325);
plot(timediff(LUMidx),burst_index(LUMidx),'o','MarkerSize',4,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(timediff(COidx),burst_index(COidx),'o','MarkerSize',4,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(timediff(Sconeidx),burst_index(Sconeidx),'o','MarkerSize',4,'MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Xlim',[100 550],'Ylim',[0 100],'YScale','log','XTick',100:150:550); xlabel('Waveform duration'); ylabel('Burst Index');

subplot(326); plot(timediff(LUMidx),100*Propbursts_returnmap(LUMidx),'o','MarkerSize',4,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on; 
plot(timediff(COidx),100*Propbursts_returnmap(COidx),'o','MarkerSize',4,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(timediff(Sconeidx),100*Propbursts_returnmap(Sconeidx),'o','MarkerSize',4,'MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Xlim',[100 550],'Ylim',[0.01 100],'YScale','log','XTick',100:150:550,'YTick',[0.01 0.1 1 10 100]); xlabel('Waveform duration'); ylabel('Burst Propensity');
set(gcf,'renderer','painters');

% Need to do some ANOVA
group = [ones(size(LUMidx)); 2*ones(size(COidx)); 3*ones(size(Sconeidx))];
data1 = burst_index([LUMidx; COidx; Sconeidx]);
data2 = Propbursts_returnmap([LUMidx; COidx; Sconeidx]);
data3 = timediff([LUMidx; COidx; Sconeidx]);

p1 = kruskalwallis(data1,group,'off');
p2 = kruskalwallis(data2,group,'off');
p3 = kruskalwallis(data3,group,'off');