% Figures for the Waveform paper titled " Chromatic tuning of putative
% inhibitory and excitatory neurons in macaque primary visual cortex."
% Author - Abhishek De, 11/19
 
close all; clearvars;
%% Figure 1: Cone weights
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
DO_conewts = [ColorOpponentIds_conewts Sconedominated_conewts];

% First thing I want to plot are the cone weights
figure(plot_counter); set(gcf,'Name','Cone wts')
subplot(211); hold on;
for ii = 1:numel(LumIds_conewts)
    h(ii) = plot(conewts_svd(1,LumIds_conewts(ii)),conewts_svd(2,LumIds_conewts(ii)),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
    set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
end
for ii = 1:numel(DO_conewts)
    if ii ~=26 | ii ~=30
        h(ii) = plot(conewts_svd(1,DO_conewts(ii)),conewts_svd(2,DO_conewts(ii)),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
        set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
    end
    if ii == 26 | ii ==30
        plot(conewts_svd(1,DO_conewts(ii)),conewts_svd(2,DO_conewts(ii)),'o','MarkerSize',8,'LineWidth',0.5,'MarkerFaceColor','c','MarkerEdgeColor',[1 1 1]); hold on;
        
    end
    
end
for ii = 1:numel(Other_conewts)
    h(ii) = plot(conewts_svd(1,Other_conewts(ii)),conewts_svd(2,Other_conewts(ii)),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
    set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
end
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out'); plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');  xlabel('L'), ylabel('M'); title('Spatially opponent');
RGB_svdSO = cell2mat(Output_List(Singleopponent & simplecells,5)');
conewts_svdSO = Mrgbtocc * RGB_svdSO;
conewts_svdSO = conewts_svdSO./repmat(sum(abs(conewts_svdSO),1),[3 1]);
conewts_svdSO = conewts_svdSO .* repmat(sign(conewts_svdSO(2,:)),[3 1]);

SO_conewts = find(sum(sign(conewts_svdSO(1:2,:)),1)==0);
LNO_conewts = find(sum(sign(conewts_svdSO(1:2,:)),1)==2);
subplot(212); hold on;
for ii = 1:numel(SO_conewts)
    h(ii) = plot(conewts_svdSO(1,SO_conewts(ii)),conewts_svdSO(2,SO_conewts(ii)),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor','c','MarkerEdgeColor',[1 1 1]); hold on;
    set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
end
for ii = 1:numel(LNO_conewts)
    h(ii) = plot(conewts_svdSO(1,LNO_conewts(ii)),conewts_svdSO(2,LNO_conewts(ii)),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor','m','MarkerEdgeColor',[1 1 1]); hold on;
    set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
end
plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out');  plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k');  xlabel('L'), ylabel('M'); title('Spatially non-opponent');
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;

%% Figure 3: Waveform width & peak-to-trough ratio
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

% PLotting the histograms of aspect ratios and waveform widths
figure(plot_counter);
subplot(311); histogram(timediff,10:30:700,'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth',2); hold on;
set(gca,'Tickdir','out','Xlim',[0 600],'Ylim',[0 70],'XTick',0:200:600,'YTick',0:35:70); ylabel('# cells'); axis square;

narrow_waveform = nan(numel(narrow_wf),1000); count_narrow = 1;
broad_waveform = nan(numel(broad_wf),1000); count_broad = 1;
indices = 1:10:size(narrow_waveform,2);
t = 2.5*(0:1:999)-500;
figure(plot_counter); subplot(312);
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
        subplot(312); plot(t(indices),narrow_waveform(count_narrow,indices),'Color','r'); hold on;
        count_narrow = count_narrow + 1;
    elseif any(ismember(broad_wf,ii))
        broad_waveform(count_broad,200-jj+1:200+L-jj) = W;
        subplot(312); plot(t(indices),broad_waveform(count_broad,indices),'Color','k'); hold on;
        count_broad = count_broad + 1;
    end
end
subplot(312); set(gca,'Tickdir','out','Xlim',[-200 600],'Ylim',[-1.1 1.1],'YTick',-1:1:1); axis square; 
xlabel('time (us)'); ylabel('Normalized Voltage'); hold off;

% Further processing for plotting the waveforms: averaged waveforms
NW = nanmean(narrow_waveform,1); NW = NW - ((min(NW)+max(NW))/2); NW = NW/max(NW); % mean narrow waveform
BW = nanmean(broad_waveform,1); BW = BW - ((min(BW)+max(BW))/2); BW = BW/max(BW); % mean broad waveform
errNW = nanstd(narrow_waveform,1)/(sqrt(numel(narrow_waveform))*max(NW));
errBW = nanstd(broad_waveform,1)/(sqrt(numel(broad_waveform))*max(BW));
% subplot(223); errorbar(t(indices),NW(indices),errNW(indices),'Linewidth',2,'Color','k'); hold on; errorbar(t(indices),BW(indices),errBW(indices),'Linewidth',2,'Color','r');  
% set(gca,'Tickdir','out','Xlim',[-200 600],'Ylim',[-1.1 1.1],'YTick',-1:1:1); axis square; 
% xlabel('time (us)'); ylabel('Normalized Voltage'); hold off;
% set(gcf,'renderer','painters'); 



% Further classifying the narrow waveforms 
NW_I = (idx==I1) & burst_index<=3;
NW_E = (idx==I1) & burst_index>3; % chattering cells
BW = idx==I2;

% Analyzing the baseline FR of all the cells
data = [baselineFR(NW_I); baselineFR(NW_E); baselineFR(BW)];
group = [ones(sum(NW_I),1); 2*ones(sum(NW_E),1); 3*ones(sum(BW),1)];
p = kruskalwallis(data,group,'off');
% Plotting the difference between narrow spiking excitatiory 
[p2,h2] = ranksum(baselineFR(NW_I),baselineFR(NW_E));
[p3,h3] = ranksum(baselineFR(NW_E),baselineFR(BW));

figure(plot_counter);
subplot(313); plot(timediff(NW_I),burst_index(NW_I),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on
plot(timediff(NW_E),burst_index(NW_E),'o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 0 0]);
plot(timediff(BW),burst_index(BW),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Xlim',[100 550],'Ylim',[0 100],'YScale','log','XTick',100:150:550); xlabel('Waveform duration'); ylabel('Burst Index');
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;


% Comparing peak to trough ratio of narrow- and broad-spiking neurons
p_pt = ranksum(abs(peak_waveform(idx==I1))./abs(trough_waveform(idx==I1)),abs(peak_waveform(idx==I2))./abs(trough_waveform(idx==I2)));

% Calculating whether there is any correlation between the waveform duration and spontaneous firing rate
baselineFR = cell2mat(Output_List(:,18));
[r1,p1] = corr(timediff,baselineFR,'type','Spearman');

% Trying to fit gaussian mixture model
GMModel1 = fitgmdist(timediff,1);
GMModel2 = fitgmdist(timediff,2);
GMModel3 = fitgmdist(timediff,3);

% Hartigan's dip test
[dip,pdip] = HartigansDipSignifTest(timediff,500);

% Bootstrapping based on Efron and Tibshirani
nboot = 500; h_opt = 25; % h_opt is obtained from the Crossvaldeterminebinsize.m 
[p,occurences] = testbimodality(h_opt,timediff,10,700,nboot);
pval = logspace(-20,-1,51);
y = nboot*binopdf(0,nboot,pval);
p_opt = pval(find(y<occurences(2),1));

%% Figure 3: Waveform width & peak-to-trough ratio
if ~exist('plot_counter')
    plot_counter = 1;
end
load Output_List_waveform.mat
load Singleopponent_waveform.mat
load timediff.mat
load peak_waveform.mat
load trough_waveform.mat

crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
timediff(~Z_cellsofinterest) = [];
peak_waveform(~Z_cellsofinterest) = [];
trough_waveform(~Z_cellsofinterest) = [];

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

RGB_svd = cell2mat(Output_List(~Singleopponent & simplecells,5)');
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
SOidx = find(Singleopponent & simplecells);


% Next, I am plotting the SVD derived luminance, color-opponent and other cells
ind = find(~Singleopponent & simplecells);
LUMidx = ind(LumIds_conewts);
COidx = ind(ColorOpponentIds_conewts);
Sconeidx = ind(Sconedominated_conewts);
DOidx = [COidx; Sconeidx];
LNOidx = SOidx(sum(sign(conewts_svdSO(1:2,:)),1)==2);
hardtoclassifyidx = [ind(Other_conewts); find(~simplecells)]; % unclassified and complex cells
SOidx(sum(sign(conewts_svdSO(1:2,:)),1)==2) = [];
idx = [LUMidx; LNOidx; DOidx; SOidx; hardtoclassifyidx]; % First DO, then simple cells, then SO and then hardtoclassify cells
group = [ones(numel(LUMidx),1); 2*ones(numel(LNOidx),1); 3*ones(numel(DOidx),1); 4*ones(numel(SOidx),1); 5*ones(numel(hardtoclassifyidx),1)];

p1 = kruskalwallis(abs(peak_waveform(idx))./abs(trough_waveform(idx)),group,'off'); % comparing peak/trough
p2 = kruskalwallis(timediff(idx),group,'off'); % comparing trough-peak time

bins = logspace(log10(0.1),log10(10),20);
% PLotting the histograms of aspect ratios and waveform widths
figure(plot_counter);
subplot(511); histogram(timediff(LUMidx),10:30:700,'FaceColor',[0 1 0],'EdgeColor',[1 1 1]); hold on
plot(median(timediff(LUMidx)),19,'v','MarkerFaceColor','g','MarkerSize',6,'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[0 600],'XTick',0:200:600,'Ylim',[0 20],'YTick',[0 20]); ylabel('# cells'); axis square;
subplot(512); hold on;
histogram(timediff(LNOidx),10:30:700,'FaceColor','m','EdgeColor',[1 1 1]);
plot(median(timediff(LNOidx)),1.5,'v','MarkerFaceColor','m','MarkerSize',6,'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[0 600],'Ylim',[0 2],'XTick',0:200:600,'YTick',[0 2]);  ylabel('# cells'); axis square;
subplot(513); hold on;
histogram(timediff(DOidx),10:30:700,'FaceColor',[1 0 0],'EdgeColor',[1 1 1]);
plot(median(timediff(DOidx)),20,'v','MarkerFaceColor','r','MarkerSize',6,'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[0 600],'Ylim',[0 30],'XTick',0:200:600,'YTick',[0 30]);  ylabel('# cells'); axis square;
subplot(514); hold on;
histogram(timediff(SOidx),10:30:700,'FaceColor','c','EdgeColor',[1 1 1]);
plot(median(timediff(SOidx)),9,'v','MarkerFaceColor','c','MarkerSize',6,'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[0 600],'Ylim',[0 10],'XTick',0:200:600,'YTick',[0 10]); xlabel('waveform width (us)'); ylabel('# cells'); axis square;
subplot(515); hold on;
histogram(timediff(hardtoclassifyidx),10:30:700,'FaceColor','k','EdgeColor',[1 1 1]);
plot(median(timediff(hardtoclassifyidx)),14,'v','MarkerFaceColor','k','MarkerSize',6,'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[0 600],'Ylim',[0 15],'XTick',0:200:600,'YTick',[0 15]); xlabel('waveform width (us)'); ylabel('# cells'); axis square;
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;

% Using k-means to divide the data into two groups:
idx = kmeans(timediff,2);
idx = kmeans(timediff,2);
[~,I1] = max([sum(idx==1) sum(idx==2)]); % narrow spiking 
[~,I2] = min([sum(idx==1) sum(idx==2)]); % broad spiking
narrow_wf = find(idx==I1);
broad_wf = find(idx==I2);

% Further classifying the narrow waveforms 
NW_I = (idx==I1) & burst_index<=3;
NW_E = (idx==I1) & burst_index>3; % chattering cells
BW = idx==I2;

% Comparing peak to trough ratio of narrow- and broad-spiking neurons
p_pt = ranksum(abs(peak_waveform(idx==I1))./abs(trough_waveform(idx==I1)),abs(peak_waveform(idx==I2))./abs(trough_waveform(idx==I2)));

% Calculating whether there is any correlation between the waveform duration and spontaneous firing rate
baselineFR = cell2mat(Output_List(:,18));
[r1,p1] = corr(timediff,baselineFR,'type','Spearman');

DO_waveform = nan(numel(DOidx),1000); countDO = 1;
LUM_waveform = nan(numel(LUMidx),1000); countLUM = 1;
LNO_waveform = nan(numel(LNOidx),1000); countLNO = 1;
SO_waveform = nan(numel(SOidx),1000); countSO = 1;
htc_waveform = nan(numel(hardtoclassifyidx),1000); counthtc = 1;
for ii = 1:size(all_waveforms,1)
    L = size(all_waveforms{ii},2);
    inds = 1:numel(all_waveforms{ii});
    [~,jj] = min(all_waveforms{ii});
    inds = inds - jj;
    W = all_waveforms{ii};
    W = W - ((min(W)+max(W))/2);
    if any(ismember(DOidx,ii))
        DO_waveform(countDO,200-jj+1:200+L-jj) = W;
        countDO = countDO + 1;
    elseif any(ismember(LUMidx,ii))
        LUM_waveform(countLUM,200-jj+1:200+L-jj) = W;
        countLUM = countLUM + 1;
    elseif any(ismember(SOidx,ii))
        SO_waveform(countSO,200-jj+1:200+L-jj) = W;
        countSO = countSO + 1;
    elseif any(ismember(LNOidx,ii))
        LNO_waveform(countLNO,200-jj+1:200+L-jj) = W;
        countLNO = countLNO + 1;
    else
        htc_waveform(counthtc,200-jj+1:200+L-jj) = W;
        counthtc = counthtc + 1;
    end
end

% Further processing for plotting the waveforms
DW = nanmean(DO_waveform,1); DW = DW - ((min(DW)+max(DW))/2); DW = DW/max(DW); % mean DO waveform
LW = nanmean(LUM_waveform,1); LW = LW - ((min(LW)+max(LW))/2); LW = LW/max(LW); % mean LUM waveform
LNW = nanmean(LNO_waveform,1); LNW = LNW - ((min(LNW)+max(LNW))/2); LNW = LNW/max(LNW); % mean LNO waveform
SW = nanmean(SO_waveform,1); SW = SW - ((min(SW(100:400))+max(SW(100:400)))/2); SW = SW/max(SW(100:400)); % mean SO waveform
HW = nanmean(htc_waveform,1); HW = HW - ((min(HW)+max(HW))/2); HW = HW/max(HW); % mean htc waveform
errDW = nanstd(DO_waveform,1)/(sqrt(numel(DOidx))*max(DW));
errLW = nanstd(LUM_waveform,1)/(sqrt(numel(LUMidx))*max(LW));
errLNW = nanstd(LNO_waveform,1)/(sqrt(numel(LNOidx))*max(LNW));
errSW = nanstd(SO_waveform,1)/(sqrt(numel(SOidx))*max(SW));
errHW = nanstd(htc_waveform,1)/(sqrt(numel(hardtoclassifyidx))*max(HW));
indices = 1:10:numel(DW);
t = 2.5*(0:1:999)-500;
figure(plot_counter); subplot(211); errorbar(t(indices),HW(indices),errHW(indices),'Linewidth',2,'Color','k'); hold on; errorbar(t(indices),DW(indices),errDW(indices),'Linewidth',2,'Color','r'); errorbar(t(indices),LW(indices),errLW(indices),'Linewidth',2,'Color','g'); 
errorbar(t(indices),SW(indices),errSW(indices),'Linewidth',2,'Color','c'); errorbar(t(indices),LNW(indices),errLNW(indices),'Linewidth',2,'Color','m'); set(gca,'Tickdir','out','Xlim',[-200 600],'Ylim',[-1.1 1.1],'YTick',-1:1:1); axis square; 
xlabel('time (us)'); ylabel('Normalized Voltage'); hold off;

% Bar plot for cells: 
NW_Icount = [sum(ismember(LUMidx,find(NW_I))) sum(ismember(LNOidx,find(NW_I))) sum(ismember(DOidx,find(NW_I))) sum(ismember(SOidx,find(NW_I))) sum(ismember(hardtoclassifyidx,find(NW_I)))];
NW_Ecount = [sum(ismember(LUMidx,find(NW_E))) sum(ismember(LNOidx,find(NW_E))) sum(ismember(DOidx,find(NW_E))) sum(ismember(SOidx,find(NW_E))) sum(ismember(hardtoclassifyidx,find(NW_E)))];
BW_count = [sum(ismember(LUMidx,find(BW))) sum(ismember(LNOidx,find(BW))) sum(ismember(DOidx,find(BW))) sum(ismember(SOidx,find(BW))) sum(ismember(hardtoclassifyidx,find(BW)))];
y = [NW_Icount; NW_Ecount; BW_count];
% subplot(312); bar(y); set(gca,'Tickdir'); axis square;

LUM_count = [sum(ismember(LUMidx,find(NW_I))) sum(ismember(LUMidx,find(NW_E))) sum(ismember(LUMidx,find(BW)))];
LNO_count = [sum(ismember(LNOidx,find(NW_I))) sum(ismember(LNOidx,find(NW_E))) sum(ismember(LNOidx,find(BW)))]; 
DO_count = [sum(ismember(DOidx,find(NW_I))) sum(ismember(DOidx,find(NW_E))) sum(ismember(DOidx,find(BW)))]; 
SO_count = [sum(ismember(SOidx,find(NW_I))) sum(ismember(SOidx,find(NW_E))) sum(ismember(SOidx,find(BW)))];
htc_count = [sum(ismember(hardtoclassifyidx,find(NW_I))) sum(ismember(hardtoclassifyidx,find(NW_E))) sum(ismember(hardtoclassifyidx,find(BW)))];
y = [LUM_count; LNO_count; DO_count; SO_count; htc_count];
subplot(212); bar(y); hold on;
set(gca,'Tickdir','out','Ylim',[0 60],'YTick',0:20:60,'XTick',[1 2 3 4 5],'Xlim',[0 6],'XTicklabel',{'LUM','LNO','DO','SO','htc'}); axis square;
set(gcf,'renderer','painters'); 
plot_counter = plot_counter + 1;


% Need to do some more statistical tests on spatial opponency andn cone-opponency 

%% Figure 4: Control Analyses, cell classification based on cone weights
if ~exist('plot_counter')
    plot_counter = 1;
end

load Output_List_waveform.mat
load Singleopponent_waveform.mat
load timediff.mat
load peak_waveform.mat
load trough_waveform.mat
load all_waveforms.mat

crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
timediff(~Z_cellsofinterest) = [];
peak_waveform(~Z_cellsofinterest) = [];
trough_waveform(~Z_cellsofinterest) = [];
Singleopponent_waveform(~Z_cellsofinterest) = [];
all_waveforms(~Z_cellsofinterest) = [];


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

RGB_svd = cell2mat(Output_List(~Singleopponent & simplecells,5)');
bkgndlms = cell2mat(Output_List(1,16)');
Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
Mrgbtocc = inv(Mrgbtocc');
conewts_svd = Mrgbtocc * RGB_svd;
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
% LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
LumIds_conewts = find(conewts_svd(2,:) + conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)-0.5).^2)<0.3);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts]) = [];
DO_conewts = [ColorOpponentIds_conewts];

% First thing I want to plot are the cone weights
figure(plot_counter); set(gcf,'Name','Control analyses')
subplot(221); hold on;
for ii = 1:numel(LumIds_conewts)
    h(ii) = plot(conewts_svd(1,LumIds_conewts(ii)),conewts_svd(2,LumIds_conewts(ii)),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
    set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
end
for ii = 1:numel(DO_conewts)
    plot(conewts_svd(1,DO_conewts(ii)),conewts_svd(2,DO_conewts(ii)),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor','r','MarkerEdgeColor',[1 1 1]); hold on;
    
end
for ii = 1:numel(Other_conewts)
    h(ii) = plot(conewts_svd(1,Other_conewts(ii)),conewts_svd(2,Other_conewts(ii)),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
    set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
end
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out'); plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');  xlabel('L'), ylabel('M'); title('Spatially opponent');
RGB_svdSO = cell2mat(Output_List(Singleopponent & simplecells,5)');
conewts_svdSO = Mrgbtocc * RGB_svdSO;
conewts_svdSO = conewts_svdSO./repmat(sum(abs(conewts_svdSO),1),[3 1]);
conewts_svdSO = conewts_svdSO .* repmat(sign(conewts_svdSO(2,:)),[3 1]);

% SO_conewts = find(sum(sign(conewts_svdSO(1:2,:)),1)==0);
SO_conewts = find(conewts_svdSO(2,:) - conewts_svdSO(1,:) >thresh & sum(sign(conewts_svdSO(1:2,:)),1)==0 & sqrt((conewts_svdSO(2,:)-0.5).^2 + (conewts_svdSO(1,:)+0.5).^2)<0.3);
LNO_conewts = find(conewts_svdSO(2,:) + conewts_svdSO(1,:) >thresh & sum(sign(conewts_svdSO(1:2,:)),1)==2 & sqrt((conewts_svdSO(2,:)-0.5).^2 + (conewts_svdSO(1,:)-0.5).^2)<0.3);
Other_conewtsSO = 1:size(conewts_svdSO,2); Other_conewtsSO([SO_conewts LNO_conewts]) = [];
% LNO_conewts = find(sum(sign(conewts_svdSO(1:2,:)),1)==2);

subplot(222); hold on;
for ii = 1:numel(SO_conewts)
    h(ii) = plot(conewts_svdSO(1,SO_conewts(ii)),conewts_svdSO(2,SO_conewts(ii)),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor','c','MarkerEdgeColor',[1 1 1]); hold on;
    set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
end
for ii = 1:numel(LNO_conewts)
    h(ii) = plot(conewts_svdSO(1,LNO_conewts(ii)),conewts_svdSO(2,LNO_conewts(ii)),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor','m','MarkerEdgeColor',[1 1 1]); hold on;
    set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
end
for ii = 1:numel(Other_conewtsSO)
    h(ii) = plot(conewts_svdSO(1,Other_conewtsSO(ii)),conewts_svdSO(2,Other_conewtsSO(ii)),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor','k','MarkerEdgeColor',[1 1 1]); hold on;
    set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
end
plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out');  plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k');  xlabel('L'), ylabel('M'); title('Spatially non-opponent');
set(gcf,'renderer','painters');

% Next, I am plotting the SVD derived luminance, color-opponent and other cells
ind = find(~Singleopponent & simplecells);
LUMidx = ind(LumIds_conewts);
DOidx = ind(ColorOpponentIds_conewts);
SOidx = find(Singleopponent & simplecells);
LNOidx = SOidx(LNO_conewts);
hardtoclassifyidx = [ind(Other_conewts); find(~simplecells); SOidx(Other_conewtsSO)]; % unclassified and complex cells
SOidx = SOidx(SO_conewts);
idx = [LUMidx; LNOidx; DOidx; SOidx; hardtoclassifyidx]; % First DO, then simple cells, then SO and then hardtoclassify cells
group = [ones(numel(LUMidx),1); 2*ones(numel(LNOidx),1); 3*ones(numel(DOidx),1); 4*ones(numel(SOidx),1); 5*ones(numel(hardtoclassifyidx),1)];

p1 = kruskalwallis(abs(peak_waveform(idx))./abs(trough_waveform(idx)),group,'off'); % comparing peak/trough
p2 = kruskalwallis(timediff(idx),group,'off'); % comparing trough-peak time


DO_waveform = nan(numel(DOidx),1000); countDO = 1;
LUM_waveform = nan(numel(LUMidx),1000); countLUM = 1;
LNO_waveform = nan(numel(LNOidx),1000); countLNO = 1;
SO_waveform = nan(numel(SOidx),1000); countSO = 1;
htc_waveform = nan(numel(hardtoclassifyidx),1000); counthtc = 1;
for ii = 1:size(all_waveforms,1)
    L = size(all_waveforms{ii},2);
    inds = 1:numel(all_waveforms{ii});
    [~,jj] = min(all_waveforms{ii});
    inds = inds - jj;
    W = all_waveforms{ii};
    W = W - ((min(W)+max(W))/2);
    if any(ismember(DOidx,ii))
        DO_waveform(countDO,200-jj+1:200+L-jj) = W;
        countDO = countDO + 1;
    elseif any(ismember(LUMidx,ii))
        LUM_waveform(countLUM,200-jj+1:200+L-jj) = W;
        countLUM = countLUM + 1;
    elseif any(ismember(SOidx,ii))
        SO_waveform(countSO,200-jj+1:200+L-jj) = W;
        countSO = countSO + 1;
    elseif any(ismember(LNOidx,ii))
        LNO_waveform(countLNO,200-jj+1:200+L-jj) = W;
        countLNO = countLNO + 1;
    else
        htc_waveform(counthtc,200-jj+1:200+L-jj) = W;
        counthtc = counthtc + 1;
    end
end

% Further processing for plotting the waveforms
DW = nanmean(DO_waveform,1); DW = DW - ((min(DW)+max(DW))/2); DW = DW/max(DW); % mean DO waveform
LW = nanmean(LUM_waveform,1); LW = LW - ((min(LW)+max(LW))/2); LW = LW/max(LW); % mean LUM waveform
LNW = nanmean(LNO_waveform,1); LNW = LNW - ((min(LNW)+max(LNW))/2); LNW = LNW/max(LNW); % mean LNO waveform
SW = nanmean(SO_waveform,1); SW = SW - ((min(SW(100:400))+max(SW(100:400)))/2); SW = SW/max(SW(100:400)); % mean SO waveform
HW = nanmean(htc_waveform,1); HW = HW - ((min(HW)+max(HW))/2); HW = HW/max(HW); % mean htc waveform
errDW = nanstd(DO_waveform,1)/(sqrt(numel(DOidx))*max(DW));
errLW = nanstd(LUM_waveform,1)/(sqrt(numel(LUMidx))*max(LW));
errLNW = nanstd(LNO_waveform,1)/(sqrt(numel(LNOidx))*max(LNW));
errSW = nanstd(SO_waveform,1)/(sqrt(numel(SOidx))*max(SW));
errHW = nanstd(htc_waveform,1)/(sqrt(numel(hardtoclassifyidx))*max(HW));
indices = 1:10:numel(DW);
t = 2.5*(0:1:999)-500;
subplot(223); errorbar(t(indices),HW(indices),errHW(indices),'Linewidth',2,'Color','k'); hold on; errorbar(t(indices),DW(indices),errDW(indices),'Linewidth',2,'Color','r'); errorbar(t(indices),LW(indices),errLW(indices),'Linewidth',2,'Color','g'); 
errorbar(t(indices),SW(indices),errSW(indices),'Linewidth',2,'Color','c'); errorbar(t(indices),LNW(indices),errLNW(indices),'Linewidth',2,'Color','m'); set(gca,'Tickdir','out','Xlim',[-200 600],'Ylim',[-1.1 1.1],'YTick',-1:1:1); axis square; 
xlabel('time (us)'); ylabel('Normalized Voltage'); hold off;


LUM_count = [sum(ismember(LUMidx,find(NW_I))) sum(ismember(LUMidx,find(NW_E))) sum(ismember(LUMidx,find(BW)))];
LNO_count = [sum(ismember(LNOidx,find(NW_I))) sum(ismember(LNOidx,find(NW_E))) sum(ismember(LNOidx,find(BW)))]; 
DO_count = [sum(ismember(DOidx,find(NW_I))) sum(ismember(DOidx,find(NW_E))) sum(ismember(DOidx,find(BW)))]; 
SO_count = [sum(ismember(SOidx,find(NW_I))) sum(ismember(SOidx,find(NW_E))) sum(ismember(SOidx,find(BW)))];
htc_count = [sum(ismember(hardtoclassifyidx,find(NW_I))) sum(ismember(hardtoclassifyidx,find(NW_E))) sum(ismember(hardtoclassifyidx,find(BW)))];
y = [LUM_count; LNO_count; DO_count; SO_count; htc_count];
subplot(224); bar(y); set(gca,'Tickdir'); axis square;
set(gcf,'renderer','painters'); 
plot_counter = plot_counter + 1;


%% Additional figures: plotting the STAs of simple, DO, SO and non-opponent cells

if ~exist('plot_counter')
    plot_counter = 1;
end
% load Output_ListWN.mat
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

RGB_svd = cell2mat(Output_List(~Singleopponent & simplecells,5)');
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
SOidx = find(Singleopponent & simplecells);


% Next, I am plotting the SVD derived luminance, color-opponent and other cells
ind = find(~Singleopponent & simplecells);
LUMidx = ind(LumIds_conewts);
COidx = ind(ColorOpponentIds_conewts);
Sconeidx = ind(Sconedominated_conewts);
DOidx = [COidx; Sconeidx];
LNOidx = SOidx(sum(sign(conewts_svdSO(1:2,:)),1)==2);
hardtoclassifyidx = [ind(Other_conewts); find(~simplecells)]; % unclassified and complex cells
SOidx(sum(sign(conewts_svdSO(1:2,:)),1)==2) = [];

% plotting STAs of Simlpe cells
figure(plot_counter); set(gcf,'Name','Simple cells');
for ii = 1:numel(LUMidx)
    tmp_vec_gun = Output_List{LUMidx(ii),2};
    normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
    im = normfactor*tmp_vec_gun + 0.5;
    im = reshape(im,[10 10 3]);
 
    subplot(14,10,ii); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
end
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;

% plotting STAs of DO cells
figure(plot_counter); set(gcf,'Name','DO cells');
for ii = 1:numel(DOidx)
    tmp_vec_gun = Output_List{DOidx(ii),2};
    normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
    im = normfactor*tmp_vec_gun + 0.5;
    im = reshape(im,[10 10 3]);
 
    subplot(14,10,ii); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
end
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;

% plotting STAs of Single-opponent cells
figure(plot_counter); set(gcf,'Name','SO');
for ii = 1:numel(SOidx)
    tmp_vec_gun = Output_List{SOidx(ii),2};
    normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
    im = normfactor*tmp_vec_gun + 0.5;
    im = reshape(im,[10 10 3]);
 
    subplot(14,10,ii); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
end
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;

% plotting STAs of Single-opponent cells
figure(plot_counter); set(gcf,'Name','LNO');
for ii = 1:numel(LNOidx)
    tmp_vec_gun = Output_List{LNOidx(ii),2};
    normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
    im = normfactor*tmp_vec_gun + 0.5;
    im = reshape(im,[10 10 3]);
 
    subplot(14,10,ii); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
end
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;

%% Addtional figure 1: plotting metrics for other the three kinds of neurons 
if ~exist('plot_counter')
    plot_counter = 1;
end
load Output_List_waveform.mat
load Singleopponent_waveform.mat
load timediff.mat
load peak_waveform.mat
load trough_waveform.mat
load all_waveforms.mat
load ISI_peaktime.mat
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
ISI_peaktime(~Z_cellsofinterest) = [];
burst_index(~Z_cellsofinterest) = [];
baselineFR = cell2mat(Output_List(:,18));


Singleopponent_waveform(~Z_cellsofinterest) = [];
NLI = cell2mat(Output_List(:,13));
simplecells = NLI<0;
Singleopponent = Singleopponent_waveform;
Singleopponent = logical(Singleopponent);


% Using k-means to divide the data into two groups:
idx = kmeans(timediff,2);
[~,I1] = max([sum(idx==1) sum(idx==2)]); % narrow spiking 
[~,I2] = min([sum(idx==1) sum(idx==2)]); % broad spiking

NW_I = (idx==I1) & burst_index<=3;
NW_E = (idx==I1) & burst_index>3; % chattering cells
BW = idx==I2;

figure(plot_counter); 
subplot(121); plot(baselineFR(NW_I),burst_index(NW_I),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on
plot(baselineFR(NW_E),burst_index(NW_E),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
plot(baselineFR(BW),burst_index(BW),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','YScale','log'); xlabel('Baseline FR'); ylabel('Burst Index');
subplot(122); plot(ISI_peaktime(NW_I),burst_index(NW_I),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on
plot(ISI_peaktime(NW_E),burst_index(NW_E),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
plot(ISI_peaktime(BW),burst_index(BW),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'XScale','log','YScale','log','Tickdir','out'); axis square; xlabel('ISI peak time'); ylabel('burst_index');
plot_counter = plot_counter + 1;

%% Additional figure 2: Analyzing if there is any bias in terms of cone or spatial opponency 
if ~exist('plot_counter')
    plot_counter = 1;
end
load Output_List_waveform.mat
load Singleopponent_waveform.mat
load timediff.mat
load peak_waveform.mat
load trough_waveform.mat

crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
timediff(~Z_cellsofinterest) = [];
peak_waveform(~Z_cellsofinterest) = [];
trough_waveform(~Z_cellsofinterest) = [];

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

RGB_svd = cell2mat(Output_List(~Singleopponent & simplecells,5)');
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
SOidx = find(Singleopponent & simplecells);


% Next, I am plotting the SVD derived luminance, color-opponent and other cells
ind = find(~Singleopponent & simplecells);
LUMidx = ind(LumIds_conewts);
COidx = ind(ColorOpponentIds_conewts);
Sconeidx = ind(Sconedominated_conewts);
DOidx = [COidx; Sconeidx];
LNOidx = SOidx(sum(sign(conewts_svdSO(1:2,:)),1)==2);
hardtoclassifyidx = [ind(Other_conewts); find(~simplecells)]; % unclassified and complex cells
SOidx(sum(sign(conewts_svdSO(1:2,:)),1)==2) = [];
idx = [LUMidx; LNOidx; DOidx; SOidx; hardtoclassifyidx]; % First DO, then simple cells, then SO and then hardtoclassify cells
group = [ones(numel(LUMidx),1); 2*ones(numel(LNOidx),1); 3*ones(numel(DOidx),1); 4*ones(numel(SOidx),1); 5*ones(numel(hardtoclassifyidx),1)];


% Using k-means to divide the data into two groups:
idx = kmeans(timediff,2);
idx = kmeans(timediff,2);
[~,I1] = max([sum(idx==1) sum(idx==2)]); % narrow spiking 
[~,I2] = min([sum(idx==1) sum(idx==2)]); % broad spiking
narrow_wf = find(idx==I1);
broad_wf = find(idx==I2);

% Further classifying the narrow waveforms 
NW_I = (idx==I1) & burst_index<=3;
NW_E = (idx==I1) & burst_index>3; % chattering cells
BW = idx==I2;

LUM_count = [sum(ismember(LUMidx,find(NW_I))) sum(ismember(LUMidx,find(NW_E))) sum(ismember(LUMidx,find(BW)))];
LNO_count = [sum(ismember(LNOidx,find(NW_I))) sum(ismember(LNOidx,find(NW_E))) sum(ismember(LNOidx,find(BW)))]; 
DO_count = [sum(ismember(DOidx,find(NW_I))) sum(ismember(DOidx,find(NW_E))) sum(ismember(DOidx,find(BW)))]; 
SO_count = [sum(ismember(SOidx,find(NW_I))) sum(ismember(SOidx,find(NW_E))) sum(ismember(SOidx,find(BW)))];
htc_count = [sum(ismember(hardtoclassifyidx,find(NW_I))) sum(ismember(hardtoclassifyidx,find(NW_E))) sum(ismember(hardtoclassifyidx,find(BW)))];
y = [LUM_count; LNO_count; DO_count; SO_count; htc_count];

tot_cells = repmat([numel(LUMidx); numel(LNOidx); numel(DOidx); numel(SOidx); numel(hardtoclassifyidx)],[1 3]);
 

SpatialOp_count = LUM_count + DO_count;
SpatialNOp_count = LNO_count + SO_count;
y = [SpatialOp_count; SpatialNOp_count];
tot_cells = repmat([numel(LUMidx) + numel(DOidx); numel(LNOidx) + numel(SOidx)],[1 3]);
figure(plot_counter); subplot(211); bar(y./tot_cells); hold on;
set(gca,'Tickdir','out','Ylim',[0 1],'YTick',0:0.25:1,'XTick',[1 2],'Xlim',[0 3],'XTicklabel',{'Space O','Space NO'});
set(gca,'Tickdir'); axis square; ylabel('% of cells'); 


ConeOp_count = SO_count + DO_count;
ConeNOp_count = LUM_count + LNO_count;
y = [ConeOp_count; ConeNOp_count];
tot_cells = repmat([numel(SOidx) + numel(DOidx); numel(LUMidx) + numel(LNOidx)],[1 3]);
subplot(212); bar(y./tot_cells); hold on; 
set(gca,'Tickdir','out','Ylim',[0 1],'YTick',0:0.25:1,'XTick',[1 2],'Xlim',[0 3],'XTicklabel',{'CO','CNO'});
axis square; ylabel('% of cells'); 
set(gcf,'renderer','painters'); 
plot_counter = plot_counter + 1;
% Need to do some more statistical tests on spatial opponency andn cone-opponency 

%% Additional figure 3: plotting scatter plots of burst index and time diff for LUM, LNO, DO, SO and unclassified cells

if ~exist('plot_counter')
    plot_counter = 1;
end
load Output_List_waveform.mat
load Singleopponent_waveform.mat
load timediff.mat
load peak_waveform.mat
load trough_waveform.mat
load burst_index.mat

crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
timediff(~Z_cellsofinterest) = [];
peak_waveform(~Z_cellsofinterest) = [];
trough_waveform(~Z_cellsofinterest) = [];
burst_index(~Z_cellsofinterest) = [];

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

RGB_svd = cell2mat(Output_List(~Singleopponent & simplecells,5)');
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
SOidx = find(Singleopponent & simplecells);


% Next, I am plotting the SVD derived luminance, color-opponent and other cells
ind = find(~Singleopponent & simplecells);
LUMidx = ind(LumIds_conewts);
COidx = ind(ColorOpponentIds_conewts);
Sconeidx = ind(Sconedominated_conewts);
DOidx = [COidx; Sconeidx];
LNOidx = SOidx(sum(sign(conewts_svdSO(1:2,:)),1)==2);
hardtoclassifyidx = [ind(Other_conewts); find(~simplecells)]; % unclassified and complex cells
SOidx(sum(sign(conewts_svdSO(1:2,:)),1)==2) = [];
idx = [LUMidx; LNOidx; DOidx; SOidx; hardtoclassifyidx]; % First DO, then simple cells, then SO and then hardtoclassify cells
group = [ones(numel(LUMidx),1); 2*ones(numel(LNOidx),1); 3*ones(numel(DOidx),1); 4*ones(numel(SOidx),1); 5*ones(numel(hardtoclassifyidx),1)];



% Using k-means to divide the data into two groups:
idx = kmeans(timediff,2);
idx = kmeans(timediff,2);
[~,I1] = max([sum(idx==1) sum(idx==2)]); % narrow spiking 
[~,I2] = min([sum(idx==1) sum(idx==2)]); % broad spiking
narrow_wf = find(idx==I1);
broad_wf = find(idx==I2);

% Further classifying the narrow waveforms 
NW_I = (idx==I1) & burst_index<=3;
NW_E = (idx==I1) & burst_index>3; % chattering cells
BW = idx==I2;


bins = logspace(log10(0.1),log10(10),20);
% PLotting the histograms of aspect ratios and waveform widths
figure(plot_counter);
subplot(321); plot(timediff(LUMidx), burst_index(LUMidx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot([290 290],[1 100],'k'); plot([100 290],[3 3],'k'); axis square; set(gca,'Tickdir','out','Ylim',[1 100],'YScale','log','Xlim',[100 550],'XTick',100:150:550); 
xlabel('time diff'); ylabel('burst index'); title('LUM');
subplot(322); plot(timediff(LNOidx), burst_index(LNOidx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot([290 290],[1 100],'k'); plot([100 290],[3 3],'k'); axis square; set(gca,'Tickdir','out','Ylim',[1 100],'YScale','log','Xlim',[100 550],'XTick',100:150:550); 
xlabel('time diff'); ylabel('burst index'); title('Non-opponent');
subplot(323); plot(timediff(DOidx), burst_index(DOidx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot([290 290],[1 100],'k'); plot([100 290],[3 3],'k'); axis square; set(gca,'Tickdir','out','Ylim',[1 100],'YScale','log','Xlim',[100 550],'XTick',100:150:550); 
xlabel('time diff'); ylabel('burst index'); title('DO');
subplot(324); plot(timediff(SOidx), burst_index(SOidx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot([290 290],[1 100],'k'); plot([100 290],[3 3],'k'); axis square; set(gca,'Tickdir','out','Ylim',[1 100],'YScale','log','Xlim',[100 550],'XTick',100:150:550); 
xlabel('time diff'); ylabel('burst index'); title('SO');
subplot(325); plot(timediff(hardtoclassifyidx), burst_index(hardtoclassifyidx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot([290 290],[1 100],'k'); plot([100 290],[3 3],'k'); axis square; set(gca,'Tickdir','out','Ylim',[1 100],'YScale','log','Xlim',[100 550],'XTick',100:150:550); 
xlabel('time diff'); ylabel('burst index'); title('Unclassified');
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;

%% Additional figure: I just want to see whether is any any correlation the previously plotted scatteroplot (burst index vs. timediff) to spatial structure 

if ~exist('plot_counter')
    plot_counter = 1;
end
load Output_List_waveform.mat
load Singleopponent_waveform.mat
load timediff.mat
load peak_waveform.mat
load trough_waveform.mat
load burst_index.mat

crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
timediff(~Z_cellsofinterest) = [];
peak_waveform(~Z_cellsofinterest) = [];
trough_waveform(~Z_cellsofinterest) = [];
burst_index(~Z_cellsofinterest) = [];

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

RGB_svd = cell2mat(Output_List(~Singleopponent & simplecells,5)');
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
SOidx = find(Singleopponent & simplecells);

% Next, I am plotting the SVD derived luminance, color-opponent and other cells
ind = find(~Singleopponent & simplecells);
LUMidx = ind(LumIds_conewts);
COidx = ind(ColorOpponentIds_conewts);
Sconeidx = ind(Sconedominated_conewts);
DOidx = [COidx; Sconeidx];
LNOidx = SOidx(sum(sign(conewts_svdSO(1:2,:)),1)==2);
hardtoclassifyidx = [ind(Other_conewts); find(~simplecells)]; % unclassified and complex cells
SOidx(sum(sign(conewts_svdSO(1:2,:)),1)==2) = [];
idx = [LUMidx; LNOidx; DOidx; SOidx; hardtoclassifyidx]; % First DO, then simple cells, then SO and then hardtoclassify cells
group = [ones(numel(LUMidx),1); 2*ones(numel(LNOidx),1); 3*ones(numel(DOidx),1); 4*ones(numel(SOidx),1); 5*ones(numel(hardtoclassifyidx),1)];

% Using k-means to divide the data into two groups:
idx = kmeans(timediff,2);
idx = kmeans(timediff,2);
[~,I1] = max([sum(idx==1) sum(idx==2)]); % narrow spiking 
[~,I2] = min([sum(idx==1) sum(idx==2)]); % broad spiking
narrow_wf = find(idx==I1);
broad_wf = find(idx==I2);

% Further classifying the narrow waveforms 
NW_I = (idx==I1) & burst_index<=3;
NW_E = (idx==I1) & burst_index>3; % chattering cells
BW = idx==I2;

load Deviation.mat
meanR = zeros(size(Deviation));
for ii = 1:size(Deviation,1)
    for jj = 1:size(Deviation,2)
        meanR(ii,jj) = mean(cos(Deviation{ii,jj}*pi/180));
    end
end

bins = logspace(log10(0.1),log10(10),20);
% PLotting the histograms of aspect ratios and waveform widths
figure(plot_counter);
subplot(321); % Gabor vs DoG comparison
for ii = 1:numel(LUMidx)
    tmp = meanR(LUMidx(ii),2) - meanR(LUMidx(ii),3); 
    plot(timediff(LUMidx(ii)), burst_index(LUMidx(ii)),'o','MarkerSize',10*(0.5+tmp),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]); hold on;
end
plot([290 290],[1 100],'k'); plot([100 290],[3 3],'k'); axis square; set(gca,'Tickdir','out','Ylim',[1 100],'YScale','log','Xlim',[100 550],'XTick',100:150:550); 
xlabel('time diff'); ylabel('burst index'); title('LUM: Gabor vs. DoG');
subplot(323); 
for ii = 1:numel(DOidx)
    tmp = meanR(DOidx(ii),2) - meanR(DOidx(ii),3); 
    plot(timediff(DOidx(ii)), burst_index(DOidx(ii)),'o','MarkerSize',10*(0.5+tmp),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]); hold on;
end

plot([290 290],[1 100],'k'); plot([100 290],[3 3],'k'); axis square; set(gca,'Tickdir','out','Ylim',[1 100],'YScale','log','Xlim',[100 550],'XTick',100:150:550); 
xlabel('time diff'); ylabel('burst index'); title('DO');
subplot(325); 
for ii = 1:numel(hardtoclassifyidx)
    tmp = meanR(hardtoclassifyidx(ii),2) - meanR(hardtoclassifyidx(ii),3); 
    plot(timediff(hardtoclassifyidx(ii)), burst_index(hardtoclassifyidx(ii)),'o','MarkerSize',10*(0.5+tmp),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]); hold on;
end
plot([290 290],[1 100],'k'); plot([100 290],[3 3],'k'); axis square; set(gca,'Tickdir','out','Ylim',[1 100],'YScale','log','Xlim',[100 550],'XTick',100:150:550); 
xlabel('time diff'); ylabel('burst index'); title('Unclassified');

subplot(322); % Gabor vs Crescent comparison
for ii = 1:numel(LUMidx)
    tmp = meanR(LUMidx(ii),2) - meanR(LUMidx(ii),1); 
    plot(timediff(LUMidx(ii)), burst_index(LUMidx(ii)),'o','MarkerSize',10*(0.5+tmp),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]); hold on;
end
plot([290 290],[1 100],'k'); plot([100 290],[3 3],'k'); axis square; set(gca,'Tickdir','out','Ylim',[1 100],'YScale','log','Xlim',[100 550],'XTick',100:150:550); 
xlabel('time diff'); ylabel('burst index'); title('LUM: Gabor vs. Crescent');
subplot(324); 
for ii = 1:numel(DOidx)
    tmp = meanR(DOidx(ii),2) - meanR(DOidx(ii),1); 
    plot(timediff(DOidx(ii)), burst_index(DOidx(ii)),'o','MarkerSize',10*(0.5+tmp),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]); hold on;
end

plot([290 290],[1 100],'k'); plot([100 290],[3 3],'k'); axis square; set(gca,'Tickdir','out','Ylim',[1 100],'YScale','log','Xlim',[100 550],'XTick',100:150:550); 
xlabel('time diff'); ylabel('burst index'); title('DO');
subplot(326); 
for ii = 1:numel(hardtoclassifyidx)
    tmp = meanR(hardtoclassifyidx(ii),2) - meanR(hardtoclassifyidx(ii),1); 
    plot(timediff(hardtoclassifyidx(ii)), burst_index(hardtoclassifyidx(ii)),'o','MarkerSize',10*(0.5+tmp),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]); hold on;
end
plot([290 290],[1 100],'k'); plot([100 290],[3 3],'k'); axis square; set(gca,'Tickdir','out','Ylim',[1 100],'YScale','log','Xlim',[100 550],'XTick',100:150:550); 
xlabel('time diff'); ylabel('burst index'); title('Unclassified');
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;

% Comparing the spatial tuning of NW_I and BW DO cells
[p1,h1] = ranksum(meanR(DOidx(ismember(DOidx,find(NW_E))),2)-meanR(DOidx(ismember(DOidx,find(NW_E))),3),meanR(DOidx(ismember(DOidx,find(BW))),2)-meanR(DOidx(ismember(DOidx,find(BW))),3));

%% For Thesis committee meeting 5

if ~exist('plot_counter')
    plot_counter = 1;
end
load Output_List_waveform.mat
load Singleopponent_waveform.mat
load timediff.mat
load peak_waveform.mat
load trough_waveform.mat

crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
timediff(~Z_cellsofinterest) = [];
peak_waveform(~Z_cellsofinterest) = [];
trough_waveform(~Z_cellsofinterest) = [];

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

RGB_svd = cell2mat(Output_List(~Singleopponent & simplecells,5)');
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
SOidx = find(Singleopponent & simplecells);


% Next, I am plotting the SVD derived luminance, color-opponent and other cells
ind = find(~Singleopponent & simplecells);
LUMidx = ind(LumIds_conewts);
COidx = ind(ColorOpponentIds_conewts);
Sconeidx = ind(Sconedominated_conewts);
DOidx = [COidx; Sconeidx];
LNOidx = SOidx(sum(sign(conewts_svdSO(1:2,:)),1)==2);
hardtoclassifyidx = [ind(Other_conewts); find(~simplecells)]; % unclassified and complex cells
SOidx(sum(sign(conewts_svdSO(1:2,:)),1)==2) = [];
idx = [LUMidx; LNOidx; DOidx; SOidx; hardtoclassifyidx]; % First DO, then simple cells, then SO and then hardtoclassify cells
group = [ones(numel(LUMidx),1); 2*ones(numel(LNOidx),1); 3*ones(numel(DOidx),1); 4*ones(numel(SOidx),1); 5*ones(numel(hardtoclassifyidx),1)];

p1 = kruskalwallis(abs(peak_waveform(idx))./abs(trough_waveform(idx)),group,'off'); % comparing peak/trough
p2 = kruskalwallis(timediff(idx),group,'off'); % comparing trough-peak time

bins = logspace(log10(0.1),log10(10),20);
% PLotting the histograms of aspect ratios and waveform widths
figure(plot_counter);
subplot(211); histogram(timediff([LUMidx;LNOidx]),10:30:700,'FaceColor',[0 1 0],'EdgeColor',[1 1 1]); hold on;
plot([290 290],[0 20],'k');
set(gca,'Tickdir','out','Xlim',[0 600],'XTick',0:200:600,'Ylim',[0 20],'YTick',[0 20]); ylabel('# cells'); axis square;
subplot(212); hold on;
histogram(timediff([DOidx;SOidx]),10:30:700,'FaceColor',[1 0 0],'EdgeColor',[1 1 1]);
plot([290 290],[0 30],'k');
set(gca,'Tickdir','out','Xlim',[0 600],'Ylim',[0 30],'XTick',0:200:600,'YTick',[0 30]);  ylabel('# cells'); axis square;
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;

