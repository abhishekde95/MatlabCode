% Waveform analysis for all the cells along with their baseline FR
% Author - Abhishek De, 10/19

close all; clearvars;
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;
resize_fact = 1;
% Include Gun noise data for statistical tests in order to estimate the RF size
channels = 3;
NPOINTS = 65536;
resize_fact2 = 1;
CHI2CRIT = 0.95; % For flaging a stixel as significant (adding gun and cone noise z-scores, squared)
maxT = 9;
crit = chi2inv(CHI2CRIT,300); % 3 color channels
spikename_options = ['sig001a'; 'sig001b'];
[filename_Lum, spikeIdx_Lum] = fnamesFromTxt2('Lum.txt');
[filename_ColorOpponent, spikeIdx_ColorOpponent] = fnamesFromTxt2('ColorOpponent.txt');

% Loading the Neurothresh files
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
tmp = fetch(conn,'SELECT filename FROM WNthresh'); tmp_filename = tmp.filename;
tmp = fetch(conn,'SELECT NTmode FROM WNthresh'); NTmode = tmp.NTmode;
tmp = fetch(conn,'SELECT spikeidx FROM WNthresh'); spikeidx_NT = tmp.spikeidx;
close(conn);
tmp_filename = tmp_filename(strcmp('subunit',NTmode));
spikeidx_NT = spikeidx_NT(strcmp('subunit',NTmode));
filename_NT = cell(size(tmp_filename));
for kk = 1:size(tmp_filename)
    filename_NT(kk) = {tmp_filename(kk)};
end

% Loading the just the subunit files
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
tmp = fetch(conn,'SELECT filename FROM WNSubunit'); tmp_filename = tmp.filename;
tmp = fetch(conn,'SELECT mode FROM WNSubunit'); subunit_mode = tmp.mode;
tmp = fetch(conn,'SELECT spikeidx FROM WNSubunit'); spikeidx_subunit = tmp.spikeidx;
close(conn);
tmp_filename = tmp_filename(strcmp('STA',subunit_mode));
spikeidx_subunit = spikeidx_subunit(strcmp('STA',subunit_mode));
filename_subunit = cell(size(tmp_filename));
for kk = 1:size(tmp_filename)
    filename_subunit(kk) = {tmp_filename(kk)};
end


Input_List = [filename_Lum; filename_ColorOpponent; filename_NT; filename_subunit];
spikeIdx = [spikeIdx_Lum; spikeIdx_ColorOpponent; spikeidx_NT; spikeidx_subunit];

numcells = numel(Input_List);
files_not_working = [];
count = 1;

load Output_List_waveform.mat
N = size(Output_List,1);
numsubplots = 10;
plot_counter = 1;
figure(plot_counter); set(gcf,'Name','Waveforms');
peak_waveform = []; trough_waveform = []; timediff = [];
peak_hw = []; trough_hw = []; % storing the half-widths
peak_time = []; trough_time = []; % storing the peak and trough times
repolarization_time = []; % repolarization time
color = [0 0 0];
all_waveforms = cell(N,1);
subplot_counter = 1;
removecells = zeros(N,1);
for ii = 1:size(Input_List,1)
    disp(ii);
    flag = 0;
    stro = {};
    for jj = 1:size(Input_List{ii},2)
        try
            tmpstro = nex2stro(findfile(char(Input_List{ii}(jj))));
        catch
            files_not_working = [files_not_working; Input_List{ii}];
            flag = 1;
            break;
        end
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
        if ~any(strcmp(stro.sum.rasterCells,'sig001a_wf'))
            files_not_working = [files_not_working; Input_List{ii}];
            flag = 1;
        end
    end
    if flag
        continue;
    end
    
    samplingrate = stro.sum.waves.storeRates{1};
    stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on');
    fpacqidx = strcmp(stro.sum.trialFields(1,:),'fpacq');
    spikeidx = strcmp(stro.sum.rasterCells(1,:),'sig001a');
    leastcount = 10^6/samplingrate; % in us
    raw_waveform = cell2mat(stro.ras(:,strcmp(stro.sum.rasterCells,'sig001a_wf')));
    idx = mean(raw_waveform,2)<10 & mean(raw_waveform,2)>-10;
    if sum(idx)>20000
        idx(20001:end) = 0; % As of now, this MATLAB along with this machine cannot handle large matrices. Therefore the max number of samples have been restricted to 20,000
    end
    [v,d] = eig(cov(raw_waveform(idx,:)));
    [~,ind] = sort(diag(d));
    vec1 = v(:,ind(end)); vec2 = v(:,ind(end-1));
    Z = linkage([raw_waveform(idx,:)*vec1 raw_waveform(idx,:)*vec2]);
    T = cluster(Z,'maxclust',2);
    cutoff = median([Z(end-2,3) Z(end-1,3)]);
    [~,idx2] = max([sum(T==1) sum(T==2)]);
    waveform = mean(raw_waveform(T==idx2,:),1);
    
%     keyboard;
%     waveform = reversefilter(waveform,1); % reverse filter 
    
    time1 = 0:25:25*(numel(waveform)-1);
    time = 0:2.5:25*(numel(waveform)-1);
    waveform = spline(time1,waveform,time);
    
    % waveform
    [val1,t1] = max(waveform);
    [val2,t2] = min(waveform);
    % Making sure that trough time is before the trough time
    flip = 0;
    if t1 < t2
        waveform = -1*waveform;
        [val1,t1] = max(waveform);
        [val2,t2] = min(waveform);
        flip = 1;
    end
    peak_waveform = [peak_waveform; val1];
    trough_waveform = [trough_waveform; val2];
    timediff = [timediff; abs(time(t1)-time(t2))];
    peak_hw = [peak_hw; range(time(waveform>val1/2))];
    trough_hw = [trough_hw; range(time(waveform<val2/2))];
    peak_time = [peak_time; time(t1)];
    trough_time = [trough_time; time(t2)];
    tmp = time(waveform>val1/2);
    repolarization_time = [repolarization_time; tmp(end)-time(t1)];
    
    A1 = islocalmin(waveform);
    A2 = islocalmax(waveform);
    
    if ~any(A1) | ~any(A2) | sum(A1)>2 | sum(A1)>2
        removecells(count) = 1;
    end
    
    
%     figure(plot_counter); subplot(numsubplots,numsubplots,subplot_counter); plot(time,waveform,'-','color',color,'Linewidth',2); hold on; plot(time(A1),waveform(A1),'r*');
%     plot(time(A2),waveform(A2),'g*');plot([time(t1) time(t1)],[0 val1],'k'); plot([time(t2) time(t2)],[0 val2],'k'); set(gca,'Xlim',[min(time) max(time)],'XTick',0:800:800); drawnow; axis square; hold off;
    
    % Storing the waveforms
    all_waveforms{count} = waveform; % each row is a waveform
    count = count + 1;
    if mod(count,100) == 1 & count > 1
        plot_counter = plot_counter + 1;
        subplot_counter = 0;
    end
    subplot_counter = subplot_counter + 1;

end
plot_counter = plot_counter + 1;

removecells = logical(removecells);
savevariables = 0;
if savevariables 
    save all_waveforms  all_waveforms
    save timediff timediff
    save peak_waveform peak_waveform
    save trough_waveform trough_waveform
    save removecells removecells
    save repolarization_time repolarization_time
end

%% Need to quantify the quality of isolation: comparing singal and noise
figure(plot_counter);
subplot(121); plot(abs(peak_waveform),abs(trough_waveform),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
set(gca,'TickDir','out','Xlim',[0 1],'Ylim',[0 1]); xlabel('Peak'); ylabel('Trough'); title('Waveform vs. noise'); axis square; hold off;
subplot(122); plot(log10(abs(peak_waveform)),log10(abs(trough_waveform)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
set(gca,'TickDir','out','Xlim',[-5 0],'Ylim',[-5 0]); xlabel('log(Peak)'); ylabel('log(Trough)'); title('Waveform vs. noise'); axis square; hold off;
plot_counter = plot_counter + 1;

% Checking whether the spike width is non-unimodal
[dip_stat,pval_dip] = HartigansDipSignifTest(sort(peak_hw+trough_hw),1000);
[dip_stat2,pval_dip2] = HartigansDipSignifTest(sort(timediff),1000);

% Analyses of waveforms
figure(plot_counter);
subplot(331); plot(abs(peak_waveform)./abs(trough_waveform),abs(timediff),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
set(gca,'TickDir','out'); xlabel('Peak/Trough'); ylabel('peak-trough time'); axis square; hold off;
subplot(332); plot(peak_hw,trough_hw,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
set(gca,'TickDir','out'); xlabel('peak half-width'); ylabel('trough half-width'); axis square; hold off;
subplot(333); plot(abs(peak_waveform)./abs(trough_waveform), peak_hw+trough_hw,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
set(gca,'Tickdir','out'); xlabel('peak/trough'); ylabel('peak-hw+trough-hw'); axis square; hold off;
subplot(334); histogram(peak_hw+trough_hw,100:30:700,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); set(gca,'Tickdir','out','Xlim',[100 700]);
xlabel('peak-hw+trough-hw'); ylabel('# cells'); title(strcat('dip p=',num2str(pval_dip,2))); axis square;
subplot(335); plot(peak_time, peak_hw,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
set(gca,'Tickdir','out'); xlabel('peak time'); ylabel('peak half-width'); axis square; hold off;
subplot(336); plot(trough_time, trough_hw,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
set(gca,'Tickdir','out'); xlabel('trough time'); ylabel('trough half-width'); axis square; hold off;
subplot(337), histogram(timediff,10:30:700,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); set(gca,'Tickdir','out','Xlim',[10 700]); xlabel('peak-trough time');
ylabel('# cells'); title(strcat('dip p=',num2str(pval_dip2,2))); axis square;
subplot(338),
for ii = 1:N
    plot(all_waveforms{ii}./max(all_waveforms{ii}),'k');  hold on;
end
set(gca,'Tickdir','out','XTick',[]); axis square; xlabel('time'); hold off;
subplot(339), histogram(log(abs(peak_waveform)./abs(trough_waveform)),10,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); set(gca,'Tickdir','out'); xlabel('peak/trough'); ylabel('# cells'); axis square;
plot_counter = plot_counter + 1;

%% Checking whether there is any pattern that correlates waveform shape with different cells

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
Singleopponent_waveform(~Z_cellsofinterest) = [];
timediff(~Z_cellsofinterest) = [];
peak_waveform(~Z_cellsofinterest) = [];
trough_waveform(~Z_cellsofinterest) = [];
all_waveforms(~Z_cellsofinterest) = [];
peak_hw(~Z_cellsofinterest) = [];
trough_hw(~Z_cellsofinterest) = [];


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
SOidx = find(Singleopponent & simplecells);

% First thing I want to plot are the cone weights
figure(plot_counter); set(gcf,'Name','Cone wts')
subplot(211);plot(conewts_svd(1,LumIds_conewts),conewts_svd(2,LumIds_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(1,ColorOpponentIds_conewts),conewts_svd(2,ColorOpponentIds_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(1,Sconedominated_conewts),conewts_svd(2,Sconedominated_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]);
plot(conewts_svd(1,Other_conewts),conewts_svd(2,Other_conewts),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out'); plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');  xlabel('L'), ylabel('M');
RGB_svdSO = cell2mat(Output_List(Singleopponent & simplecells,5)');
conewts_svdSO = Mrgbtocc * RGB_svdSO;
conewts_svdSO = conewts_svdSO./repmat(sum(abs(conewts_svdSO),1),[3 1]);
conewts_svdSO = conewts_svdSO .* repmat(sign(conewts_svdSO(2,:)),[3 1]);
subplot(212); plot(conewts_svdSO(1,sum(sign(conewts_svdSO(1:2,:)),1)==0),conewts_svdSO(2,sum(sign(conewts_svdSO(1:2,:)),1)==0),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor','c','MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svdSO(1,sum(sign(conewts_svdSO(1:2,:)),1)==2),conewts_svdSO(2,sum(sign(conewts_svdSO(1:2,:)),1)==2),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor','m','MarkerEdgeColor',[1 1 1]); 
plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out');  plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k');  xlabel('L'), ylabel('M');
plot_counter = plot_counter + 1;

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

% Doing a bunch of Kruskal-Wallis Test
p1 = kruskalwallis(abs(peak_waveform(idx))./abs(trough_waveform(idx)),group,'off'); % comparing peak/trough
p2 = kruskalwallis(peak_hw(idx)+trough_hw(idx),group,'off'); % comparing spike width
p3 = kruskalwallis(timediff(idx),group,'off'); % comparing trough-peak time

figure(plot_counter);
subplot(231); boxplot(abs(peak_waveform(idx))./abs(trough_waveform(idx)),group); set(gca,'Tickdir','out','XTicklabel',{'Lum','LNO','DO','SO','oth'}); ylabel('peak/trough');
title(strcat('p=',num2str(p1,2))); axis square;
subplot(232); boxplot(peak_hw(idx)+trough_hw(idx),group); set(gca,'Tickdir','out','XTicklabel',{'Lum','LNO','DO','SO','oth'}); ylabel('peak-hw+trough-hw');
title(strcat('p=',num2str(p2,2))); axis square;
subplot(233); boxplot(timediff(idx),group); set(gca,'Tickdir','out','XTicklabel',{'Lum','LNO','DO','SO','oth'}); ylabel('peak-trough time');
title(strcat('p=',num2str(p3,2))); axis square;
subplot(234),
for ii = 1:numel(all_waveforms)
    inds = 1:numel(all_waveforms{ii});
    [~,jj] = max(all_waveforms{ii});
    inds = inds - jj;
    if any(ismember(DOidx,ii))
        plot(inds,all_waveforms{ii}./max(all_waveforms{ii}),'r');  hold on;
    elseif any(ismember(LUMidx,ii))
        plot(inds,all_waveforms{ii}./max(all_waveforms{ii}),'g');  hold on;
    elseif any(ismember(SOidx,ii))
        plot(inds,all_waveforms{ii}./max(all_waveforms{ii}),'c');  hold on;
    else
        plot(inds, all_waveforms{ii}./max(all_waveforms{ii}),'k');  hold on;
    end
end
set(gca,'Tickdir','out','XTick',[],'Ylim',[-3 1.5]); axis square; hold off;
subplot(235); plot(abs(peak_waveform(DOidx))./abs(trough_waveform(DOidx)),abs(timediff(DOidx)),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(abs(peak_waveform(LUMidx))./abs(trough_waveform(LUMidx)),abs(timediff(LUMidx)),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
plot(abs(peak_waveform(SOidx))./abs(trough_waveform(SOidx)),abs(timediff(SOidx)),'o','MarkerFaceColor','c','MarkerEdgeColor',[1 1 1]);
plot(abs(peak_waveform(hardtoclassifyidx))./abs(trough_waveform(hardtoclassifyidx)),abs(timediff(hardtoclassifyidx)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'TickDir','out'); xlabel('Peak/Trough'); ylabel('peak-trough time'); axis square; hold off;
subplot(236); plot(abs(peak_waveform(DOidx))./abs(trough_waveform(DOidx)),peak_hw(DOidx)+trough_hw(DOidx),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(abs(peak_waveform(LUMidx))./abs(trough_waveform(LUMidx)),peak_hw(LUMidx)+trough_hw(LUMidx),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
plot(abs(peak_waveform(SOidx))./abs(trough_waveform(SOidx)),peak_hw(SOidx)+trough_hw(SOidx),'o','MarkerFaceColor','c','MarkerEdgeColor',[1 1 1]);
plot(abs(peak_waveform(hardtoclassifyidx))./abs(trough_waveform(hardtoclassifyidx)),peak_hw(hardtoclassifyidx)+trough_hw(hardtoclassifyidx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'TickDir','out'); xlabel('Peak/Trough'); ylabel('peak-hw+trough-hw'); axis square; hold off;
plot_counter = plot_counter + 1;

% Looking at the histograms of peak-trough ratio, spike width and trough-peak time for different cells
figure(plot_counter);
subplot(321); histogram(log(abs(peak_waveform)./abs(trough_waveform)),-2:0.2:2,'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth',2); hold on;
histogram(log(abs(peak_waveform(DOidx))./abs(trough_waveform(DOidx))),-2:0.2:2,'FaceColor',[1 0 0],'EdgeColor',[1 1 1]);
histogram(log(abs(peak_waveform(LUMidx))./abs(trough_waveform(LUMidx))),-2:0.2:2,'FaceColor',[0 1 0],'EdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-2 2],'XTick',-2:2:2); xlabel('log(peak/trough)'); ylabel('# cells'); title('DO vs. LUM'); axis square;
subplot(322); histogram(log(abs(peak_waveform)./abs(trough_waveform)),-2:0.2:2,'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth',2); hold on;
histogram(log(abs(peak_waveform(hardtoclassifyidx))./abs(trough_waveform(hardtoclassifyidx))),-2:0.2:2,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
histogram(log(abs(peak_waveform(SOidx))./abs(trough_waveform(SOidx))),-2:0.2:2,'FaceColor','c','EdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-2 2],'XTick',-2:2:2); xlabel('log(peak/trough)'); ylabel('# cells'); title('others'); axis square;
subplot(323); histogram(peak_hw+trough_hw,100:30:700,'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth',2); hold on;
histogram(peak_hw(DOidx)+trough_hw(DOidx),100:30:700,'FaceColor',[1 0 0],'EdgeColor',[1 1 1]);
histogram(peak_hw(LUMidx)+trough_hw(LUMidx),100:30:700,'FaceColor',[0 1 0],'EdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[100 700]); xlabel('peak-hw+trough-hw'); ylabel('# cells'); title('DO vs. LUM'); axis square;
subplot(324); histogram(peak_hw+trough_hw,100:30:700,'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth',2); hold on;
histogram(peak_hw(hardtoclassifyidx)+trough_hw(hardtoclassifyidx),100:30:700,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
histogram(peak_hw(SOidx)+trough_hw(SOidx),100:30:700,'FaceColor','c','EdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[100 700]); xlabel('peak-hw+trough-hw'); ylabel('# cells'); title('others'); axis square;
subplot(325); histogram(timediff,10:30:700,'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth',2); hold on;
histogram(timediff(DOidx),10:30:700,'FaceColor',[1 0 0],'EdgeColor',[1 1 1]);
histogram(timediff(LUMidx),10:30:700,'FaceColor',[0 1 0],'EdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[10 700]); xlabel('peak-trough time'); ylabel('# cells'); title('DO vs. LUM'); axis square;
subplot(326); histogram(timediff,10:30:700,'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth',2); hold on;
histogram(timediff(hardtoclassifyidx),10:30:700,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
histogram(timediff(SOidx),10:30:700,'FaceColor','c','EdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[10 700]); xlabel('peak-trough time'); ylabel('# cells'); title('others'); axis square;
plot_counter = plot_counter+1;

% Trying to compute a average waveform of the cells separated by the the
% DO, LUM and hardtoclassify classes.
figure(plot_counter);
DO_waveform = nan(numel(DOidx),1000); countDO = 1;
LUM_waveform = nan(numel(LUMidx),1000); countLUM = 1;
LNO_waveform = nan(numel(LNOidx),1000); countLNO = 1;
SO_waveform = nan(numel(SOidx),1000); countSO = 1;
htc_waveform = nan(numel(hardtoclassifyidx),1000); counthtc = 1;
for ii = 1:numel(all_waveforms)
    L = size(all_waveforms{ii},2);
    inds = 1:numel(all_waveforms{ii});
    [~,jj] = min(all_waveforms{ii});
    inds = inds - jj;
    W = all_waveforms{ii};
    W = W - ((min(W)+max(W))/2);
    if any(ismember(DOidx,ii))
        DO_waveform(countDO,200-jj+1:200+L-jj) = W;
        subplot(231); plot(inds,W/max(W),'r');  hold on;
        countDO = countDO + 1;
    elseif any(ismember(LUMidx,ii))
        LUM_waveform(countLUM,200-jj+1:200+L-jj) = W;
        subplot(233); plot(inds,W/max(W),'g');  hold on;
        countLUM = countLUM + 1;
    elseif any(ismember(SOidx,ii))
        SO_waveform(countSO,200-jj+1:200+L-jj) = W;
        subplot(232); plot(inds,W/max(W),'c');  hold on;
        countSO = countSO + 1;
    elseif any(ismember(LNOidx,ii))
        LNO_waveform(countLNO,200-jj+1:200+L-jj) = W;
        subplot(234); plot(inds,W/max(W),'m');  hold on;
        countLNO = countLNO + 1;
    else
        htc_waveform(counthtc,200-jj+1:200+L-jj) = W;
        subplot(235); plot(inds,W/max(W),'k');  hold on;
        counthtc = counthtc + 1;
    end
end
subplot(231); set(gca,'TickDir','out','Ylim',[-1 1],'YTick',-1:1:1); ylabel('Normalized Voltage'); title('DO'); axis square;
subplot(232); set(gca,'TickDir','out','Ylim',[-1 1],'YTick',-1:1:1); ylabel('Normalized Voltage'); title('SO'); axis square;
subplot(233); set(gca,'TickDir','out','Ylim',[-1 1],'YTick',-1:1:1); ylabel('Normalized Voltage'); title('LUM'); axis square;
subplot(234); set(gca,'TickDir','out','Ylim',[-1 1],'YTick',-1:1:1); ylabel('Normalized Voltage'); title('LNO'); axis square;
subplot(235); set(gca,'TickDir','out','Ylim',[-1 1],'YTick',-1:1:1); ylabel('Normalized Voltage'); title('others'); axis square;
% Further processing for plotting the waveforms
DW = nanmean(DO_waveform,1); DW = DW - ((min(DW)+max(DW))/2); DW = DW/max(DW); % mean DO waveform
LW = nanmean(LUM_waveform,1); LW = LW - ((min(LW)+max(LW))/2); LW = LW/max(LW); % mean LUM waveform
LNW = nanmean(LNO_waveform,1); LNW = LNW - ((min(LNW)+max(LNW))/2); LNW = LNW/max(LNW); % mean LNO waveform
SW = nanmean(SO_waveform,1); SW = SW - ((min(SW)+max(SW))/2); SW = SW/max(SW); % mean SO waveform
HW = nanmean(htc_waveform,1); HW = HW - ((min(HW)+max(HW))/2); HW = HW/max(HW); % mean htc waveform
errDW = nanstd(DO_waveform,1)/(sqrt(numel(DOidx))*max(DW));
errLW = nanstd(LUM_waveform,1)/(sqrt(numel(LUMidx))*max(LW));
errLNW = nanstd(LNO_waveform,1)/(sqrt(numel(LNOidx))*max(LNW));
errSW = nanstd(SO_waveform,1)/(sqrt(numel(SOidx))*max(SW));
errHW = nanstd(htc_waveform,1)/(sqrt(numel(hardtoclassifyidx))*max(HW));
indices = 1:10:numel(DW);
subplot(236); errorbar(indices,DW(indices),errDW(indices),'Linewidth',2,'Color','r'); hold on; errorbar(indices,LW(indices),errLW(indices),'Linewidth',2,'Color','g'); errorbar(indices,SW(indices),errSW(indices),'Linewidth',2,'Color','c');
errorbar(indices,HW(indices),errHW(indices),'Linewidth',2,'Color','k'); errorbar(indices,LNW(indices),errLNW(indices),'Linewidth',2,'Color','m'); set(gca,'Tickdir','out','Xlim',[100 500],'Ylim',[-1.1 1.1],'YTick',-1:1:1); axis square; ylabel('Normalized Voltage'); hold off;
plot_counter = plot_counter + 1;

% Summary plot
figure(plot_counter); set(gcf,'Name','Summary plot'); 
subplot(221);plot(conewts_svd(1,LumIds_conewts),conewts_svd(2,LumIds_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(1,ColorOpponentIds_conewts),conewts_svd(2,ColorOpponentIds_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(1,Sconedominated_conewts),conewts_svd(2,Sconedominated_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(conewts_svd(1,Other_conewts),conewts_svd(2,Other_conewts),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out'); plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');  xlabel('L'), ylabel('M'); title('Spatially opponent');
subplot(222); plot(conewts_svdSO(1,sum(sign(conewts_svdSO(1:2,:)),1)==0),conewts_svdSO(2,sum(sign(conewts_svdSO(1:2,:)),1)==0),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor','c','MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svdSO(1,sum(sign(conewts_svdSO(1:2,:)),1)==2),conewts_svdSO(2,sum(sign(conewts_svdSO(1:2,:)),1)==2),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor','k','MarkerEdgeColor',[1 1 1]); 
plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out');  plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k');  xlabel('L'), ylabel('M'); title('Spatially non-opponent');
subplot(223); errorbar(indices,DW(indices),errDW(indices),'Linewidth',2,'Color','r'); hold on; errorbar(indices,LW(indices),errLW(indices),'Linewidth',2,'Color','g'); errorbar(indices,SW(indices),errSW(indices),'Linewidth',2,'Color','c');
errorbar(indices,HW(indices),errHW(indices),'Linewidth',2,'Color','k'); set(gca,'Tickdir','out','Xlim',[100 500],'Ylim',[-1.1 1.1],'YTick',-1:1:1); axis square; ylabel('Normalized Voltage'); 
legend('DO','LUM','SO','Others'); hold off;
subplot(224); histogram(timediff(LUMidx),10:30:700,'Normalization','probability','DisplayStyle','stairs','EdgeColor',[0 1 0],'Linewidth',2); hold on;
plot(median(timediff(LUMidx)),0.18,'v','MarkerFaceColor','g','MarkerSize',6,'MarkerEdgeColor',[1 1 1]);
histogram(timediff(DOidx),10:30:700,'Normalization','probability','DisplayStyle','stairs','EdgeColor',[1 0 0],'Linewidth',2);
plot(median(timediff(DOidx)),0.18,'v','MarkerFaceColor','r','MarkerSize',6,'MarkerEdgeColor',[1 1 1]);
histogram(timediff(SOidx),10:30:700,'Normalization','probability','DisplayStyle','stairs','EdgeColor','c','Linewidth',2);
plot(median(timediff(SOidx)),0.18,'v','MarkerFaceColor','c','MarkerSize',6,'MarkerEdgeColor',[1 1 1]);
histogram(timediff(hardtoclassifyidx),10:30:700,'Normalization','probability','DisplayStyle','stairs','EdgeColor','k','Linewidth',2);
plot(median(timediff(hardtoclassifyidx)),0.19,'v','MarkerFaceColor','k','MarkerSize',6,'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Ylim',[0 0.2],'YTick',0:0.1:0.2,'Xlim',[0 700]); xlabel('waveform width (us)'); ylabel('Proportion of cells'); axis square; hold off;
plot_counter = plot_counter + 1;

% PLotting the histograms of aspect ratios and waveform widths
figure(plot_counter);
subplot(521); histogram(timediff,10:30:700,'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth',2); hold on;
histogram(timediff(LUMidx),10:30:700,'FaceColor',[0 1 0],'EdgeColor',[1 1 1]);
plot(median(timediff(LUMidx)),20,'v','MarkerFaceColor','g','MarkerSize',6,'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[0 600],'Ylim',[0 65]); ylabel('# cells'); axis square;
subplot(523); hold on;
histogram(timediff(LNOidx),10:30:700,'FaceColor','m','EdgeColor',[1 1 1]);
plot(median(timediff(LNOidx)),4,'v','MarkerFaceColor','m','MarkerSize',6,'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[0 600],'Ylim',[0 5]);  ylabel('# cells'); axis square;
subplot(525); hold on;
histogram(timediff(DOidx),10:30:700,'FaceColor',[1 0 0],'EdgeColor',[1 1 1]);
plot(median(timediff(DOidx)),20,'v','MarkerFaceColor','r','MarkerSize',6,'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[0 600],'Ylim',[0 30]);  ylabel('# cells'); axis square;
subplot(527); hold on;
histogram(timediff(SOidx),10:30:700,'FaceColor','c','EdgeColor',[1 1 1]);
plot(median(timediff(SOidx)),20,'v','MarkerFaceColor','c','MarkerSize',6,'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[0 600],'Ylim',[0 30]); xlabel('waveform width (us)'); ylabel('# cells'); axis square;
subplot(529); hold on;
histogram(timediff(hardtoclassifyidx),10:30:700,'FaceColor','k','EdgeColor',[1 1 1]);
plot(median(timediff(hardtoclassifyidx)),20,'v','MarkerFaceColor','k','MarkerSize',6,'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[0 600],'Ylim',[0 30]); xlabel('waveform width (us)'); ylabel('# cells'); axis square;
subplot(522); histogram(log(abs(peak_waveform)./abs(trough_waveform)),-2:0.2:2,'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth',2); hold on;
histogram(log(abs(peak_waveform(LUMidx))./abs(trough_waveform(LUMidx))),-2:0.2:2,'FaceColor',[0 1 0],'EdgeColor',[1 1 1]);
plot(median(log(abs(peak_waveform(LUMidx))./abs(trough_waveform(LUMidx)))),30,'v','MarkerFaceColor','g','MarkerSize',6,'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-2 2],'XTick',-2:2:2);  axis square;
subplot(524); hold on;
histogram(log(abs(peak_waveform(LNOidx))./abs(trough_waveform(LNOidx))),-2:0.2:2,'FaceColor','m','EdgeColor',[1 1 1]);
plot(median(log(abs(peak_waveform(LNOidx))./abs(trough_waveform(LNOidx)))),9,'v','MarkerFaceColor','m','MarkerSize',6,'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-2 2],'XTick',-2:2:2,'Ylim',[0 10]); axis square;
subplot(526); hold on;
histogram(log(abs(peak_waveform(DOidx))./abs(trough_waveform(DOidx))),-2:0.2:2,'FaceColor','r','EdgeColor',[1 1 1]);
plot(median(log(abs(peak_waveform(DOidx))./abs(trough_waveform(DOidx)))),39,'v','MarkerFaceColor','r','MarkerSize',6,'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-2 2],'XTick',-2:2:2,'Ylim',[0 40]); axis square;
subplot(528); hold on;
histogram(log(abs(peak_waveform(SOidx))./abs(trough_waveform(SOidx))),-2:0.2:2,'FaceColor','c','EdgeColor',[1 1 1]);
plot(median(log(abs(peak_waveform(SOidx))./abs(trough_waveform(SOidx)))),19,'v','MarkerFaceColor','c','MarkerSize',6,'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-2 2],'XTick',-2:2:2,'Ylim',[0 20]);  axis square;
subplot(5,2,10); hold on;
histogram(log(abs(peak_waveform(hardtoclassifyidx))./abs(trough_waveform(hardtoclassifyidx))),-2:0.2:2,'FaceColor','k','EdgeColor',[1 1 1]);
plot(median(log(abs(peak_waveform(hardtoclassifyidx))./abs(trough_waveform(hardtoclassifyidx)))),19,'v','MarkerFaceColor','k','MarkerSize',6,'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-2 2],'XTick',-2:2:2,'Ylim',[0 20]); xlabel('log(peak/trough)'); axis square;
plot_counter = plot_counter + 1;

%% Some control population analyses looking at the effect of scaling the waveform on
% a)peak-trough ratio
% b)trough time - peak time
% c)peak half-width + trough half-width

var_peaktotrough = [];
var_peakminustroughtime = [];
var_peakhwplustroughhw = [];
for ii = 1:N
    stro = nex2stro(findfile(char(filename(ii))));
    samplingrate = stro.sum.waves.storeRates{1};
    stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on');
    fpacqidx = strcmp(stro.sum.trialFields(1,:),'fpacq');
    spikeidx = strcmp(stro.sum.rasterCells(1,:),'sig001a');
    WVF = mean(cell2mat(stro.ras(:,strcmp(stro.sum.rasterCells,'sig001a_wf')))',2);
    peak_waveform = []; trough_waveform = []; timediff = []; peak_hw = [];
    trough_hw = []; peak_time = []; trough_time = [];
    scaling_factor = logspace(-3,3,11);
    for ii = 1:numel(scaling_factor)
        waveform = WVF*scaling_factor(ii);
        time1 = 0:25:25*(numel(waveform)-1);
        time = 0:2.5:25*(numel(waveform)-1);
        waveform = spline(time1,waveform,time);
        % waveform
        [val1,t1] = max(waveform);
        [val2,t2] = min(waveform);
        % Making sure that peak time is before the trough time
        flip = 0;
        if t1 > t2
            waveform = -1*waveform;
            [val1,t1] = max(waveform);
            [val2,t2] = min(waveform);
            flip = 1;
        end
        peak_waveform = [peak_waveform; val1];
        trough_waveform = [trough_waveform; val2];
        timediff = [timediff; abs(time(t2)-time(t1))];
        peak_hw = [peak_hw; range(time(waveform>val1/2))];
        trough_hw = [trough_hw; range(time(waveform<val2/2))];
        peak_time = [peak_time; time(t1)];
        trough_time = [trough_time; time(t2)];
    end
    var_peaktotrough = [var_peaktotrough; var(peak_waveform./trough_waveform)];
    var_peakminustroughtime = [var_peakminustroughtime; var(timediff)];
    var_peakhwplustroughhw = [var_peakhwplustroughhw; var(peak_hw+trough_hw)];
    
    if ii == 1
        figure(plot_counter);
        subplot(221); plot(scaling_factor,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
        set(gca,'YScale','log'); axis square; ylabel('Scaling factor');
        subplot(222), plot(scaling_factor,abs(peak_waveform./trough_waveform),'o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]); hold on;
        plot(scaling_factor,abs(peak_waveform),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
        plot(scaling_factor,abs(trough_waveform),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
        set(gca,'XScale','log','YScale','log'); axis square; ylabel('peak/trough'); xlabel('scaling factor'); legend('peak/trough','peak','trough'); hold off;
        subplot(223), plot(scaling_factor,timediff,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
        set(gca,'XScale','log'); axis square; ylabel('peak-trough time'); xlabel('scaling factor');
        subplot(224), plot(scaling_factor,peak_hw+trough_hw,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
        set(gca,'XScale','log'); axis square; ylabel('peak-hw+trough-hw'); xlabel('scaling factor');
        plot_counter = plot_counter + 1;
    end
end

figure(plot_counter);
subplot(221); histogram(var_peaktotrough,-0.2:0.02:0.2,'FaceColor',[0 0 0],'EdgeColor',[0 0 0]); axis square; ylabel('#cells');
xlabel('var(peak/trough)'); set(gca,'Tickdir','out','Xlim',[-0.2 0.2],'XTick',-0.2:0.2:0.2);
subplot(222); histogram(var_peakminustroughtime,-0.2:0.02:0.2,'FaceColor',[0 0 0],'EdgeColor',[0 0 0]); axis square; ylabel('#cells');
xlabel('var(peak-trough time)'); set(gca,'Tickdir','out','Xlim',[-0.2 0.2],'XTick',-0.2:0.2:0.2);
subplot(223); histogram(var_peakhwplustroughhw,-0.2:0.02:0.2,'FaceColor',[0 0 0],'EdgeColor',[0 0 0]); axis square; ylabel('#cells');
xlabel('var(peak-hw+trough-hw)'); set(gca,'Tickdir','out','Xlim',[-0.2 0.2],'XTick',-0.2:0.2:0.2);
plot_counter = plot_counter + 1;

