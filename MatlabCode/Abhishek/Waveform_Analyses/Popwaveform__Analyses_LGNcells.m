% Waveform analyses for Greg's LGN data
% Author - Abhishek De, 10/19

close all; clearvars;
conn = database('Nex_Paradigm_Sort','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
tmp_filename = fetch(conn,'SELECT fileID FROM IsoSamp_LGN');
quality = cell2mat(fetch(conn,'SELECT quality FROM IsoSamp_LGN'));
cellclass = fetch(conn,'SELECT cellClass FROM IsoSamp_LGN');
unitidx = fetch(conn,'SELECT spikeCode FROM IsoSamp_LGN');
close(conn);
filename = tmp_filename(quality); % Selecting files in which quality = 1
cellclass = char(cellclass(quality));
unitidx = unitidx(quality);

% Printing the name of the usable LGN filename in a .txt file
printfilenames = 0;
if printfilenames
    fid = fopen('LGNfilenames.txt','w');
    for ii = 1:size(filename,1)
        fprintf(fid,'%s\r\n',string(filename(ii)));
    end
    fclose(fid);
end


plot_counter = 1;
N = numel(filename);
numsubplots = ceil(sqrt(N));
figure(plot_counter); set(gcf,'Name','Waveforms');
peak_waveform = []; trough_waveform = []; timediff = [];
peak_hw = []; trough_hw = []; % storing the half-widths
peak_time = []; trough_time = []; % storing the peak and trough times
peak_noise_waveform = []; trough_noise_waveform = [];
color = [0 0 0];
all_waveforms = cell(N,1);
all_noise_waveforms = cell(N,1);
euclid_waveform = [];
euclid_noise_waveform = [];
baselineFR = [];
varFR = [];
files_waveform_not_included = [];
for ii = 1:numel(filename)
    stro = nex2stro(findfile(char(filename{ii})));
    if ~any(strcmp(stro.sum.rasterCells,unitidx(ii)))
        files_waveform_not_included = [files_waveform_not_included; filename{ii}];
    end
    
    samplingrate = stro.sum.waves.storeRates{1};
    stimonidx = strcmp(stro.sum.trialFields(1,:),'stimon_t');
    fpacqidx = strcmp(stro.sum.trialFields(1,:),'fpacq_t');
    spikeidx = strcmp(stro.sum.rasterCells(1,:),unitidx(ii)); 
    leastcount = 10^6/samplingrate; % in us
    waveform = mean(cell2mat(stro.ras(:,strcmp(stro.sum.rasterCells,strcat(unitidx(ii),'_wf'))))',2);
    std_waveform = std(cell2mat(stro.ras(:,strcmp(stro.sum.rasterCells,strcat(unitidx(ii),'_wf'))))',0,2);
    noise_waveform = mean(cell2mat(stro.ras(:,strcmp(stro.sum.rasterCells,'sig001U_wf')))',2);
    std_noise_waveform = std(cell2mat(stro.ras(:,strcmp(stro.sum.rasterCells,'sig001U_wf')))',0,2);
    time1 = 0:25:25*(numel(waveform)-1);
    time = 0:2.5:25*(numel(waveform)-1);
    waveform = spline(time1,waveform,time);
    std_waveform = spline(time1,std_waveform,time);
    noise_waveform = spline(time1,noise_waveform,time);
    std_noise_waveform = spline(time1,std_noise_waveform,time);
    
    FR = [];
    for jj = 1:size(stro.ras,1)
        FR = [FR; sum(stro.ras{jj,spikeidx}>stro.trial(jj,fpacqidx) & stro.ras{jj,spikeidx}<stro.trial(jj,stimonidx))/(stro.trial(jj,stimonidx) - stro.trial(jj,fpacqidx))];
    end
    baselineFR = [baselineFR; nanmean(FR)];
    varFR = [varFR; nanvar(FR)];
    
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
    
    % Noise_waveform
    if flip
        noise_waveform = -1*noise_waveform;
    end
    
    [val1,~] = max(noise_waveform);
    [val2,~] = min(noise_waveform);
    peak_noise_waveform = [peak_noise_waveform; val1]; 
    trough_noise_waveform = [trough_noise_waveform; val2];
    
    subplot(numsubplots,numsubplots,ii); plot(time,waveform,'-','color',color,'Linewidth',2); hold on; plot(time,noise_waveform,'-','color',[0.5 0.5 0.5],'Linewidth',2);
    plot(time,waveform-std_waveform,'-','color',color,'Linewidth',1); plot(time,waveform+std_waveform,'-','color',color,'Linewidth',1);
    plot(time,noise_waveform-std_noise_waveform,'-','color',[0.5 0.5 0.5],'Linewidth',1); plot(time,noise_waveform+std_noise_waveform,'-','color',[0.5 0.5 0.5],'Linewidth',1);
    line([0 max(time)],[0 0]); set(gca,'Xlim',[min(time) max(time)],'XTick',0:800:800); drawnow; axis square; hold off;
    
    % Storing the waveforms
    all_waveforms{ii} = waveform; % each row is a waveform
    all_noise_waveforms{ii} = noise_waveform;
    euclid_waveform = [euclid_waveform; sqrt(sum(waveform.^2))];
    euclid_noise_waveform = [euclid_noise_waveform; sqrt(sum(noise_waveform.^2))];
end

plot_counter = plot_counter + 1;

% Need to quantify the quality of isolation: comparing singal and noise
figure(plot_counter); 
subplot(221); plot(euclid_waveform,euclid_noise_waveform,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on; % plotting the waveform in PC space
set(gca,'TickDir','out','Xlim',[0 5],'Ylim',[0 5],'XTick',0:1:5,'YTick',0:1:5); line([0 5],[0 5]); xlabel('Waveform'); ylabel('Noise'); title('Euclidean distance'); axis square; hold off;
subplot(222); plot(abs(peak_waveform),abs(trough_waveform),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(abs(peak_noise_waveform),abs(trough_noise_waveform),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); 
set(gca,'TickDir','out','Xlim',[0 1],'Ylim',[0 1]); xlabel('Peak'); ylabel('Trough'); title('Waveform vs. noise'); axis square; hold off; 
subplot(223); plot(log10(abs(peak_waveform)),log10(abs(trough_waveform)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(log10(abs(peak_noise_waveform)),log10(abs(trough_noise_waveform)),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); 
set(gca,'TickDir','out','Xlim',[-5 0],'Ylim',[-5 0]); xlabel('log(Peak)'); ylabel('log(Trough)'); title('Waveform vs. noise'); axis square; hold off; 
plot_counter = plot_counter + 1;

% Checking whether the spike width is non-unimodal
[dip_stat,pval_dip] = HartigansDipSignifTest(peak_hw+trough_hw,1000);
[dip_stat2,pval_dip2] = HartigansDipSignifTest(timediff,1000);

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
subplot(337), histogram(timediff,100:30:700,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); set(gca,'Tickdir','out','Xlim',[100 700]); xlabel('peak-trough time'); 
ylabel('# cells'); title(strcat('dip p=',num2str(pval_dip2,2))); axis square;
subplot(338),
for ii = 1:N
    plot(all_waveforms{ii}./max(all_waveforms{ii}),'k');  hold on;
end
set(gca,'Tickdir','out','XTick',[]); axis square; xlabel('time'); hold off;
subplot(339), histogram(log(abs(peak_waveform)./abs(trough_waveform)),10,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); set(gca,'Tickdir','out'); xlabel('peak/trough'); ylabel('# cells'); axis square;
plot_counter = plot_counter + 1;


% Comparing the spike waveform metrics to that of firing statistics
figure(plot_counter);
subplot(221); plot(baselineFR,varFR,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
xlabel('baseline FR'); ylabel('var FR'); set(gca,'Tickdir','out'); axis square;
subplot(222); plot(peak_hw+trough_hw,baselineFR,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
xlabel('peak-hw+trough-hw'); ylabel('baseline FR'); set(gca,'Tickdir','out'); axis square; 
subplot(223); plot(timediff, baselineFR,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
xlabel('peak-trough time'); ylabel('baseline FR'); set(gca,'Tickdir','out'); axis square; 
subplot(224); plot(abs(peak_waveform)./abs(trough_waveform), baselineFR,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
xlabel('peak/trough'); ylabel('baseline FR'); set(gca,'Tickdir','out'); axis square; 
plot_counter = plot_counter + 1;

%% Comparing waveform stats between parvo and magno 
Parvoidx = find(cellclass(:,1)=='P');
Magnoidx = find(cellclass(:,1)=='M');
Konioidx = find(cellclass(:,1)=='K');
idx = [Parvoidx; Magnoidx; Konioidx]; % First DO, then simple cells, then SO and then hardtoclassify cells
group = [ones(numel(Parvoidx),1);2*ones(numel(Magnoidx),1);3*ones(numel(Konioidx),1)];

% Doing a bunch of Kruskal-Wallis Test
p1 = kruskalwallis(abs(peak_waveform(idx))./abs(trough_waveform(idx)),group,'off'); % comparing peak/trough 
p2 = kruskalwallis(peak_hw(idx)+trough_hw(idx),group,'off'); % comparing spike width
p3 = kruskalwallis(timediff(idx),group,'off'); % comparing trough-peak time
p4 = kruskalwallis(baselineFR(idx),group,'off'); % comparing trough-peak time

figure(plot_counter); 
subplot(331); boxplot(abs(peak_waveform(idx))./abs(trough_waveform(idx)),group); set(gca,'Tickdir','out','XTicklabel',{'P','M','K'}); ylabel('peak/trough'); 
title(strcat('p=',num2str(p1,2))); axis square;
subplot(332); boxplot(peak_hw(idx)+trough_hw(idx),group); set(gca,'Tickdir','out','XTicklabel',{'P','M','K'}); ylabel('peak-hw+trough-hw'); 
title(strcat('p=',num2str(p2,2))); axis square;
subplot(333); boxplot(timediff(idx),group); set(gca,'Tickdir','out','XTicklabel',{'P','M','K'}); ylabel('peak-trough time'); 
title(strcat('p=',num2str(p3,2))); axis square;
subplot(334),
for ii = 1:N
    inds = 1:numel(all_waveforms{ii});
    [~,jj] = max(all_waveforms{ii});
    inds = inds - jj;
    if any(ismember(Parvoidx,ii))
        plot(inds,all_waveforms{ii}./max(all_waveforms{ii}),'r');  hold on;
    elseif any(ismember(Magnoidx,ii))
        plot(inds,all_waveforms{ii}./max(all_waveforms{ii}),'g');  hold on;
    else any(ismember(Konioidx,ii))
        plot(inds, all_waveforms{ii}./max(all_waveforms{ii}),'b');  hold on;
    end
end
set(gca,'Tickdir','out','XTick',[],'Ylim',[-3 1.5]); axis square; hold off;
subplot(335); plot(abs(peak_waveform(Parvoidx))./abs(trough_waveform(Parvoidx)),abs(timediff(Parvoidx)),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(abs(peak_waveform(Magnoidx))./abs(trough_waveform(Magnoidx)),abs(timediff(Magnoidx)),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
plot(abs(peak_waveform(Konioidx))./abs(trough_waveform(Konioidx)),abs(timediff(Konioidx)),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]);
set(gca,'TickDir','out'); xlabel('Peak/Trough'); ylabel('peak-trough time'); axis square; hold off; 
subplot(336); plot(abs(peak_waveform(Parvoidx))./abs(trough_waveform(Parvoidx)),peak_hw(Parvoidx)+trough_hw(Parvoidx),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(abs(peak_waveform(Magnoidx))./abs(trough_waveform(Magnoidx)),peak_hw(Magnoidx)+trough_hw(Magnoidx),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
plot(abs(peak_waveform(Konioidx))./abs(trough_waveform(Konioidx)),peak_hw(Konioidx)+trough_hw(Konioidx),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]);
set(gca,'TickDir','out'); xlabel('Peak/Trough'); ylabel('peak-hw+trough-hw'); axis square; hold off; 
subplot(337); boxplot(baselineFR(idx),group); set(gca,'Tickdir','out','XTicklabel',{'P','M','K'}); ylabel('baselineFR'); 
title(strcat('p=',num2str(p4,2))); axis square;
subplot(338); plot(peak_hw(Parvoidx)+trough_hw(Parvoidx),baselineFR(Parvoidx),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(peak_hw(Magnoidx)+trough_hw(Magnoidx),baselineFR(Magnoidx),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(peak_hw(Konioidx)+trough_hw(Konioidx),baselineFR(Konioidx),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]); hold on;
xlabel('peak-hw+trough-hw'); ylabel('baseline FR'); set(gca,'Tickdir','out'); axis square; hold off;
subplot(339); plot(timediff(Parvoidx), baselineFR(Parvoidx),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on; 
plot(timediff(Magnoidx), baselineFR(Magnoidx),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); 
plot(timediff(Konioidx), baselineFR(Konioidx),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]); 
xlabel('peak-trough time'); ylabel('baseline FR'); set(gca,'Tickdir','out'); axis square; hold off;
plot_counter = plot_counter + 1;

% Looking at the histograms of peak-trough ratio, spike width and trough-peak time for different cells
figure(plot_counter);
subplot(321); histogram(log(abs(peak_waveform)./abs(trough_waveform)),-2:0.2:2,'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth',2); hold on;
histogram(log(abs(peak_waveform(Parvoidx))./abs(trough_waveform(Parvoidx))),-2:0.2:2,'FaceColor',[1 0 0],'EdgeColor',[1 1 1]);
histogram(log(abs(peak_waveform(Magnoidx))./abs(trough_waveform(Magnoidx))),-2:0.2:2,'FaceColor',[0 1 0],'EdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-2 2],'XTick',-2:2:2); xlabel('log(peak/trough)'); ylabel('# cells'); title('P vs. M'); axis square;
subplot(322); histogram(log(abs(peak_waveform)./abs(trough_waveform)),-2:0.2:2,'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth',2); hold on;
histogram(log(abs(peak_waveform(Konioidx))./abs(trough_waveform(Konioidx))),-2:0.2:2,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-2 2],'XTick',-2:2:2); xlabel('log(peak/trough)'); ylabel('# cells'); title('K'); axis square;
subplot(323); histogram(peak_hw+trough_hw,100:30:700,'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth',2); hold on;
histogram(peak_hw(Parvoidx)+trough_hw(Parvoidx),100:30:700,'FaceColor',[1 0 0],'EdgeColor',[1 1 1]);
histogram(peak_hw(Magnoidx)+trough_hw(Magnoidx),100:30:700,'FaceColor',[0 1 0],'EdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[100 700]); xlabel('peak-hw+trough-hw'); ylabel('# cells'); title('P vs. M'); axis square;
subplot(324); histogram(peak_hw+trough_hw,100:30:700,'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth',2); hold on;
histogram(peak_hw(Konioidx)+trough_hw(Konioidx),100:30:700,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[100 700]); xlabel('peak-hw+trough-hw'); ylabel('# cells'); title('K'); axis square;
subplot(325); histogram(timediff,100:30:700,'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth',2); hold on;
histogram(timediff(Parvoidx),100:30:700,'FaceColor',[1 0 0],'EdgeColor',[1 1 1]);
histogram(timediff(Magnoidx),100:30:700,'FaceColor',[0 1 0],'EdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[100 700]); xlabel('peak-trough time'); ylabel('# cells'); title('P vs. M'); axis square;
subplot(326); histogram(timediff,100:30:700,'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth',2); hold on;
histogram(timediff(Konioidx),100:30:700,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[100 700]); xlabel('peak-trough time'); ylabel('# cells'); title('K'); axis square;
plot_counter = plot_counter+1; 

% Trying to compute a average waveform of the cells separated by the the
% DO, LUM and hardtoclassify classes. 
figure(plot_counter);
Parvo_waveform = nan(numel(Parvoidx),1000); countP = 1;
Magno_waveform = nan(numel(Magnoidx),1000); countM = 1;
Konio_waveform = nan(numel(Konioidx),1000); countK = 1;
for ii = 1:N
    L = size(all_waveforms{ii},2);
    inds = 1:numel(all_waveforms{ii});
    [~,jj] = max(all_waveforms{ii});
    inds = inds - jj;
    W = all_waveforms{ii};
    W = W - ((min(W)+max(W))/2);
    if any(ismember(Parvoidx,ii))
        Parvo_waveform(countP,200-jj+1:200+L-jj) = W; 
        subplot(221); plot(inds,W/max(W),'r');  hold on;
        countP = countP + 1;
    elseif any(ismember(Magnoidx,ii))
        Magno_waveform(countM,200-jj+1:200+L-jj) = W;
        subplot(222); plot(inds,W/max(W),'g');  hold on;
        countM = countM + 1;
    else any(ismember(Konioidx,ii))
        Konio_waveform(countK,200-jj+1:200+L-jj) = W;
        subplot(223); plot(inds,W/max(W),'k');  hold on;
        countK = countK + 1;
    end
end
subplot(221); set(gca,'TickDir','out','Ylim',[-1 1],'YTick',-1:1:1); ylabel('Normalized Voltage'); title('Parvo'); axis square;
subplot(222); set(gca,'TickDir','out','Ylim',[-1 1],'YTick',-1:1:1); ylabel('Normalized Voltage'); title('Magno'); axis square;
subplot(223); set(gca,'TickDir','out','Ylim',[-1 1],'YTick',-1:1:1); ylabel('Normalized Voltage'); title('Konio'); axis square;
% Further processing for plotting the waveforms 
DW = nanmean(Parvo_waveform,1); DW = DW - ((min(DW(150:400))+max(DW(150:400)))/2); DW = DW/max(DW(150:400)); % mean Parvo waveform
LW = nanmean(Magno_waveform,1); LW = LW - ((min(LW(150:400))+max(LW(150:400)))/2); LW = LW/max(LW(150:400)); % mean Magno waveform
HW = nanmean(Konio_waveform,1); HW = HW - ((min(HW(150:400))+max(HW(150:400)))/2); HW = HW/max(HW(150:400)); % mean Konio waveform
errDW = nanstd(Parvo_waveform,1)/(sqrt(numel(Parvoidx))*max(DW(150:400))); 
errLW = nanstd(Magno_waveform,1)/(sqrt(numel(Magnoidx))*max(LW(150:400)));
errHW = nanstd(Konio_waveform,1)/(sqrt(numel(Konioidx))*max(HW(150:400)));
subplot(224); errorbar(DW,errDW,'Linewidth',2,'Color','r'); hold on; errorbar(LW,errLW,'Linewidth',2,'Color','g');
errorbar(HW,errHW,'Linewidth',2,'Color','k'); set(gca,'Tickdir','out','Xlim',[150 500],'Ylim',[-1.1 1.1],'YTick',-1:1:1); axis square; ylabel('Normalized Voltage');
legend('Parvo','Magno','Konio'); hold off; 
plot_counter = plot_counter + 1;