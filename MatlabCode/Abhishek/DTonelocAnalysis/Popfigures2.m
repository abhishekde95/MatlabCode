% Developing a second script for making figures for the paper after consulting with Greg and Yasmine.
% Author - Abhishek De, 3/19,
close all; clearvars;

%% Figure 1: Cell counting histology figure
% Using Virusy3-V1-Dlx-ADIHC-20X-001_YS2_forFIJI.tif
if ~exist('plot_counter')
    plot_counter = 1;
end
T_PV = readtable('C:\Users\setup\Google Drive\UW\GRAD_SCHOOL_PAPERS\OPTO_PAPER\cell_counting\prelimResults_PV.csv');
T_ChR2 = readtable('C:\Users\setup\Google Drive\UW\GRAD_SCHOOL_PAPERS\OPTO_PAPER\cell_counting\prelimResults_ChR2.csv');
figure(plot_counter);
plot(T_PV.X,1732-T_PV.Y,'o','MarkerFaceColor',[0 1 0],'MarkerSize',10); hold on;
plot(T_ChR2.X,1732-T_ChR2.Y,'o','MarkerFaceColor',[1 0 0],'MarkerSize',6); axis equal, set(gca,'Xlim',[0 1424],'Ylim',[0 1732],'Tickdir','out');
xlabel('X'), ylabel('Y'); legend('PV','ChR2'); hold off;
plot_counter = plot_counter + 1;

% Calculating 95 % confidence interval
p = 37/41;
err = 1.96*sqrt(p*(1-p)/41);

%% Figure 2
% �	Neurophysiology: One example each of suppression and activation rasters, singles-units along with their waveforms in inset (2 figures)
% �	Scatter plot of all suppression and activation units, classified as single and multi-units (1 figure)
% �	One example of a single-unit from Fixstim interfreq and Gratings displaying a putative �direction selective� inhibitory neuron (2 figures)

if ~exist('plot_counter')
    plot_counter = 1;
end
filename = {'A012319005_wf.nex'; 'A020119004_wf.nex'};
titles = {'activation';'suppresion'};
N = numel(filename);
L = ceil(sqrt(N));
binwidth = .005; % 5 ms bin width
bins = -0.4:binwidth:0.6;
baselinebins = (bins>-0.3 & bins<0);
baselineFR = [];
laserFR = [];
statstestresult = []; % Performing a wilcoxon rank sum test and storing the p value
figure(plot_counter); set(gcf,'Name','Neuronal data');
for jj = 1:N
    stro = nex2stro(findfile(char(filename(jj,:))));
    fponidx  = strcmp(stro.sum.trialFields(1,:),'fpon_t');
    fpacqidx  = strcmp(stro.sum.trialFields(1,:),'fpacq_t');
    fpoffidx  = strcmp(stro.sum.trialFields(1,:),'fpoff_t');
    stimonidx  = strcmp(stro.sum.trialFields(1,:),'stimon_t');
    stimoffidx  = strcmp(stro.sum.trialFields(1,:),'stimoff_t');
    targonidx  = strcmp(stro.sum.trialFields(1,:),'targon_t');
    saccstartidx  = strcmp(stro.sum.trialFields(1,:),'saccstart_t');
    saccendidx  = strcmp(stro.sum.trialFields(1,:),'saccend_t');
    rewidx  = strcmp(stro.sum.trialFields(1,:),'rew_t');
    stimpresentidx  = strcmp(stro.sum.trialFields(1,:),'stimpresent');
    Lcc = strcmp(stro.sum.trialFields(1,:),'lcc');
    Mcc = strcmp(stro.sum.trialFields(1,:),'mcc');
    Scc = strcmp(stro.sum.trialFields(1,:),'scc');
    tf = strcmp(stro.sum.trialFields(1,:),'tf');
    oog = strcmp(stro.sum.trialFields(1,:),'oog');
    stimidx = strcmp(stro.sum.trialFields(1,:),'stim_idx');
    optstim = strcmp(stro.sum.trialFields(1,:),'optstim');
    correct = strcmp(stro.sum.trialFields(1,:),'correct');
    laseron = strcmp(stro.sum.trialFields(1,:),'laseron_t');
    laseroff = strcmp(stro.sum.trialFields(1,:),'laseroff_t');
    analogstrtimeind = strcmp(stro.sum.rasterCells(1,:),'anlgStartTime');
    spikeind = strcmp(stro.sum.rasterCells(1,:),'sig001a');
    heyetraceind = strcmp(stro.sum.rasterCells(1,:),'AD11');
    veyetraceind = strcmp(stro.sum.rasterCells(1,:),'AD12');
    lasertraceind = strcmp(stro.sum.rasterCells(1,:),'AD13');
    correcttrials = stro.trial(:,correct);
    diridxs = stro.trial(:,stimidx);
    stimpresent = logical(stro.trial(:,stimpresentidx));
    LMStriplet = [stro.trial(:,Lcc) stro.trial(:,Mcc) stro.trial(:,Scc)];
    lasertrials = logical(stro.trial(:,optstim));
    samplerate = stro.sum.analog.storeRates{3}; % sample rate at which laser analog pulses are stored in file
    
    % Stimulus absent trials
    idxs = find(~stimpresent);
    count1 = 1;
    count2 = 1;
    PSTHlaser = [];
    tmp_baselineFR = [];
    tmp_laserFR = [];
    tmp_spikecountsbaseline = [];
    tmp_spikecountslaserstimabsent = [];
    
    for ii = 1:numel(idxs)
        ind = idxs(ii);
        analogstartime = stro.ras{ind,analogstrtimeind};
        spiketimes = stro.ras{ind,spikeind};
        t = stro.ras{ind,analogstrtimeind}+[0:1:length(stro.ras{ind,lasertraceind})-1]/samplerate;
        if lasertrials(ind)
            % Pulling out the laseron and off timings from analog traces
            laserontime = t(find(stro.ras{ind,lasertraceind}>0.1,1));
            laserofftime = t(find(stro.ras{ind,lasertraceind}>0.1,1,'last'));
            timedurlaser = laserofftime - laserontime;
            subplot(3,3,3*jj-2); plot(spiketimes(spiketimes>laserontime-0.3 & spiketimes<laserofftime+0.2)-laserontime,count1*ones(size(spiketimes(spiketimes>laserontime-0.3 & spiketimes<laserofftime+0.2))),'k.'); hold on;
            PSTHlaser = [PSTHlaser;  hist(spiketimes-laserontime, bins)/binwidth];
            laserbins = bins>=0 & bins<=timedurlaser;
            tmp_laserFR = [tmp_laserFR; sum(spiketimes>laserontime & spiketimes<laserofftime)/0.3];
            count1 = count1 + 1;
        else
            stimontime = stro.trial(ind,stimonidx);
            tmp_baselineFR = [tmp_baselineFR; sum(spiketimes>stimontime & spiketimes<stimontime+0.3)/0.3 ];
        end
    end
    PSTHlaser = mean(PSTHlaser,1);
    line([0 0],[0 count1]); line([timedurlaser timedurlaser],[0 count1]); xlabel('time'); ylabel('trials'); title(titles{jj});
    set(gca,'Xlim',[-0.15 0.45],'Tickdir','out','XTick',-0.15:0.15:0.45,'Ylim',[0 count1-1]); axis square; hold off;
    
    subplot(3,3,3*jj-1); plot(bins,PSTHlaser,'b','Linewidth',2);
    if jj == 1
        hold on; line([0 0],[0 600]); line([timedurlaser timedurlaser],[0 600]); xlabel('time'); ylabel('FR');
        set(gca,'Xlim',[-0.15 0.45],'Tickdir','out','XTick',-0.15:0.15:0.45,'Ylim',[0 600],'YTick',[0:300:600]); axis square; hold off;
    else
        hold on; line([0 0],[0 100]); line([timedurlaser timedurlaser],[0 100]); xlabel('time'); ylabel('FR');
        set(gca,'Xlim',[-0.15 0.45],'Tickdir','out','XTick',-0.15:0.15:0.45,'Ylim',[0 100],'YTick',[0:50:100]); axis square; hold off;
    end
    
    baselineFR = [baselineFR; mean(tmp_baselineFR)];
    laserFR = [laserFR; mean(tmp_laserFR)];
    samplingrate = stro.sum.waves.storeRates{1};
    leastcount = 10^6/samplingrate; % in us
    waveform = mean(cell2mat(stro.ras(:,3))',2);
    noise_waveform = mean(cell2mat(stro.ras(:,4))',2);
    time = 0:25:25*(numel(waveform)-1);
    [val1,t1] = max(waveform);
    [val2,t2] = min(waveform);
    subplot(3,3,3*jj); plot(time,waveform,'-','color',[0 0 0],'Linewidth',2); hold on; plot(time,noise_waveform,'-','color',[0.5 0.5 0.5],'Linewidth',2);
    plot(time,waveform-std(cell2mat(stro.ras(:,3))',0,2),'-','color',[0 0 0],'Linewidth',1); plot(time,waveform+std(cell2mat(stro.ras(:,3))',0,2),'-','color',[0 0 0],'Linewidth',1);
    plot(time,noise_waveform-std(cell2mat(stro.ras(:,4))',0,2),'-','color',[0.5 0.5 0.5],'Linewidth',1); plot(time,noise_waveform+std(cell2mat(stro.ras(:,4))',0,2),'-','color',[0.5 0.5 0.5],'Linewidth',1);
    line([0 max(time)],[0 0]); set(gca,'Xlim',[min(time) max(time)],'XTick',0:200:800,'Tickdir','out'); axis square; title('Waveform: single unit'); hold off;
    
end

% Population summary of excitation and suppression
if ~exist('plot_counter')
    plot_counter = 1;
end
filename_M = {'M051318008.nex';'M051318009.nex';'M051618004.nex';'M051618005.nex';'M051618006.nex';'M051618007.nex';'M051618008.nex';'M051618010.nex';...
    'M051618011.nex';'M051618012.nex';'M051618013.nex';'M051618014.nex';'M051818004.nex';'M051818005.nex';'M051818006.nex';'M051818007.nex';...
    'M051818010.nex';'M051818011.nex';'M051818012.nex';'M051818013.nex';'M051818014.nex';'M052318001.nex';'M052318002.nex';'M052318003.nex';'M052318004.nex'};
filename_A = {'A011719009.nex';'A011719010.nex';'A011719011.nex';'A011719012.nex';'A012219006.nex';'A012219009.nex';'A012219010.nex';'A012219011.nex';'A012219012.nex';...
    'A012319004.nex';'A012319005.nex';'A012319006.nex';'A013119002.nex';'A013119003.nex';'A013119004.nex';'A013119005.nex';'A013119006.nex';'A020119004.nex';'A020119005.nex';...
    'A020319001.nex';'A020319007.nex';'A020319008.nex';'A020319009.nex';'A020519003.nex';'A020519008.nex';'A020519009.nex';'A021519004.nex';'A021519005.nex';'A021519006.nex';...
    'A021519007.nex';'A021519008.nex'};
filename = [filename_M; filename_A];
Singlevsmultiunit_M = ['M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M'];
Singlevsmultiunit_A = ['M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'S';'S';'S';'M';'S';'S';'M';'S';'M';'M';'M';'S';'M';'M';'M';'M';'M';'S';'S';'S';'S'];
N = numel(filename);
L = ceil(sqrt(N));
binwidth = .005;
bins = -0.4:binwidth:0.6;
baselinebins = (bins>0 & bins<0.3);
baselineFR = [];
laserFR = [];
statstestresult = []; % Performing a wilcoxon rank sum test and storing the p value
ISI = [];
for jj = 1:N
    stro = nex2stro(findfile(char(filename(jj,:))));
    fponidx  = strcmp(stro.sum.trialFields(1,:),'fpon_t');
    fpacqidx  = strcmp(stro.sum.trialFields(1,:),'fpacq_t');
    fpoffidx  = strcmp(stro.sum.trialFields(1,:),'fpoff_t');
    stimonidx  = strcmp(stro.sum.trialFields(1,:),'stimon_t');
    stimoffidx  = strcmp(stro.sum.trialFields(1,:),'stimoff_t');
    targonidx  = strcmp(stro.sum.trialFields(1,:),'targon_t');
    saccstartidx  = strcmp(stro.sum.trialFields(1,:),'saccstart_t');
    saccendidx  = strcmp(stro.sum.trialFields(1,:),'saccend_t');
    rewidx  = strcmp(stro.sum.trialFields(1,:),'rew_t');
    stimpresentidx  = strcmp(stro.sum.trialFields(1,:),'stimpresent');
    Lcc = strcmp(stro.sum.trialFields(1,:),'lcc');
    Mcc = strcmp(stro.sum.trialFields(1,:),'mcc');
    Scc = strcmp(stro.sum.trialFields(1,:),'scc');
    tf = strcmp(stro.sum.trialFields(1,:),'tf');
    oog = strcmp(stro.sum.trialFields(1,:),'oog');
    stimidx = strcmp(stro.sum.trialFields(1,:),'stim_idx');
    optstim = strcmp(stro.sum.trialFields(1,:),'optstim');
    correct = strcmp(stro.sum.trialFields(1,:),'correct');
    laseron = strcmp(stro.sum.trialFields(1,:),'laseron_t');
    laseroff = strcmp(stro.sum.trialFields(1,:),'laseroff_t');
    analogstrtimeind = strcmp(stro.sum.rasterCells(1,:),'anlgStartTime');
    spikeind = strcmp(stro.sum.rasterCells(1,:),'sig001a');
    heyetraceind = strcmp(stro.sum.rasterCells(1,:),'AD11');
    veyetraceind = strcmp(stro.sum.rasterCells(1,:),'AD12');
    lasertraceind = strcmp(stro.sum.rasterCells(1,:),'AD13');
    correcttrials = stro.trial(:,correct);
    diridxs = stro.trial(:,stimidx);
    stimpresent = logical(stro.trial(:,stimpresentidx));
    LMStriplet = [stro.trial(:,Lcc) stro.trial(:,Mcc) stro.trial(:,Scc)];
    lasertrials = logical(stro.trial(:,optstim));
    samplerate = stro.sum.analog.storeRates{3}; % sample rate at which laser analog pulses are stored in file
    
    % Stimulus absent trials
    idxs = find(~stimpresent);
    count1 = 1;
    count2 = 1;
    PSTHlaser = zeros(1,length(bins));
    PSTHbaseline = zeros(1,length(bins));
    PSTHnolaser = zeros(1,length(bins));
    tmp_baselineFR = [];
    tmp_laserFR = [];
    tmp_spikecountsbaseline = [];
    tmp_spikecountslaserstimabsent = [];
    
    for ii = 1:numel(idxs)
        ind = idxs(ii);
        analogstartime = stro.ras{ind,analogstrtimeind};
        spiketimes = stro.ras{ind,spikeind};
        if lasertrials(ind)
            % Pulling out the laseron and off timings from analog traces
            t = stro.ras{ind,analogstrtimeind}+[0:1:length(stro.ras{ind,lasertraceind})-1]/samplerate;
            laserontime = t(find(stro.ras{ind,lasertraceind}>0.1,1));
            laserofftime = t(find(stro.ras{ind,lasertraceind}>0.1,1,'last'));
            timedurlaser = laserofftime - laserontime;
            PSTHlaser = PSTHlaser + (hist(spiketimes-laserontime, bins)/binwidth);
            count1 = count1 + 1;
            laserbins = bins>=0 & bins<=timedurlaser;
            tmp_laserFR = [tmp_laserFR; sum(spiketimes>laserontime & spiketimes<laserofftime)/0.3];
            ISI = [ISI; diff(spiketimes(spiketimes>laserontime & spiketimes<laserofftime))];
            tmp_spikecountsbaseline = [tmp_spikecountsbaseline; sum(spiketimes>laserontime-0.3 & spiketimes<laserontime)];
            tmp_spikecountslaserstimabsent = [tmp_spikecountslaserstimabsent; sum(spiketimes>laserontime & spiketimes<laserofftime)];
        else
            if ~stimpresent(ind)
                stimontime = stro.trial(ind,stimonidx);
                PSTHbaseline = PSTHbaseline + hist(spiketimes-stimontime, bins)/binwidth;
                tmp_baselineFR = [tmp_baselineFR; sum(spiketimes>stimontime & spiketimes<stimontime+0.3)/0.3 ];
            end
        end
    end
    baselineFR = [baselineFR; mean(tmp_baselineFR)];
    laserFR = [laserFR; mean(tmp_laserFR)];
    statstestresult = [statstestresult; ranksum(tmp_spikecountsbaseline,tmp_spikecountslaserstimabsent)];
end

firstletters_filename = char([filename]);
Singlevsmultiunit = [Singlevsmultiunit_M; Singlevsmultiunit_A];
Singleunitidxs = Singlevsmultiunit == 'S';
Maui_idxs = firstletters_filename(:,1)=='M';
Apollo_idxs = firstletters_filename(:,1)=='A';
figure(plot_counter);
subplot(337); plot(baselineFR(~Singleunitidxs),laserFR(~Singleunitidxs),'o','MarkerSize',7,'LineWidth',1.0,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(baselineFR(Singleunitidxs),laserFR(Singleunitidxs),'o','MarkerSize',7,'LineWidth',1.0,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
axis square; xlabel('baselineFR'); ylabel('laserFR'); set(gca,'Xlim',[0 300],'Ylim',[0 300],'TickDir','out','XTick',[0:150:300],'YTick',[0:150:300]); line([0 300],[0 300]); legend('Multi','Single'); hold off;
subplot(338), plot(baselineFR(~Singleunitidxs),laserFR(~Singleunitidxs),'o','MarkerSize',7,'LineWidth',1.0,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(baselineFR(Singleunitidxs),laserFR(Singleunitidxs),'o','MarkerSize',7,'LineWidth',1.0,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
axis square; xlabel('baselineFR'); ylabel('laserFR'); set(gca,'Xlim',[0 30],'Ylim',[0 30],'TickDir','out','XTick',[0:15:30],'YTick',[0:15:30]); line([0 30],[0 30]); legend('Multi','Single'); hold off;
plot_counter = plot_counter + 1;


% Calculating the time of the first spike of the excited cells (Supplementary figure 1 in the new draft)
filename = filename(laserFR>baselineFR);
N = numel(filename);
L = ceil(sqrt(N));
timeoffirstspike = cell(N,1);
meantime = [];
stdtime = [];
groupID = [];
for jj = 1:N
    stro = nex2stro(findfile(char(filename(jj,:))));
    fponidx  = strcmp(stro.sum.trialFields(1,:),'fpon_t');
    fpacqidx  = strcmp(stro.sum.trialFields(1,:),'fpacq_t');
    fpoffidx  = strcmp(stro.sum.trialFields(1,:),'fpoff_t');
    stimonidx  = strcmp(stro.sum.trialFields(1,:),'stimon_t');
    stimoffidx  = strcmp(stro.sum.trialFields(1,:),'stimoff_t');
    targonidx  = strcmp(stro.sum.trialFields(1,:),'targon_t');
    saccstartidx  = strcmp(stro.sum.trialFields(1,:),'saccstart_t');
    saccendidx  = strcmp(stro.sum.trialFields(1,:),'saccend_t');
    rewidx  = strcmp(stro.sum.trialFields(1,:),'rew_t');
    stimpresentidx  = strcmp(stro.sum.trialFields(1,:),'stimpresent');
    Lcc = strcmp(stro.sum.trialFields(1,:),'lcc');
    Mcc = strcmp(stro.sum.trialFields(1,:),'mcc');
    Scc = strcmp(stro.sum.trialFields(1,:),'scc');
    tf = strcmp(stro.sum.trialFields(1,:),'tf');
    oog = strcmp(stro.sum.trialFields(1,:),'oog');
    stimidx = strcmp(stro.sum.trialFields(1,:),'stim_idx');
    optstim = strcmp(stro.sum.trialFields(1,:),'optstim');
    correct = strcmp(stro.sum.trialFields(1,:),'correct');
    laseron = strcmp(stro.sum.trialFields(1,:),'laseron_t');
    laseroff = strcmp(stro.sum.trialFields(1,:),'laseroff_t');
    analogstrtimeind = strcmp(stro.sum.rasterCells(1,:),'anlgStartTime');
    spikeind = strcmp(stro.sum.rasterCells(1,:),'sig001a');
    heyetraceind = strcmp(stro.sum.rasterCells(1,:),'AD11');
    veyetraceind = strcmp(stro.sum.rasterCells(1,:),'AD12');
    lasertraceind = strcmp(stro.sum.rasterCells(1,:),'AD13');
    correcttrials = stro.trial(:,correct);
    diridxs = stro.trial(:,stimidx);
    stimpresent = logical(stro.trial(:,stimpresentidx));
    LMStriplet = [stro.trial(:,Lcc) stro.trial(:,Mcc) stro.trial(:,Scc)];
    lasertrials = logical(stro.trial(:,optstim));
    samplerate = stro.sum.analog.storeRates{3}; % sample rate at which laser analog pulses are stored in file
    
    % Stimulus absent trials
    idxs = find(~stimpresent);
    firstspiketime_tmp = [];
    
    for ii = 1:numel(idxs)
        ind = idxs(ii);
        analogstartime = stro.ras{ind,analogstrtimeind};
        spiketimes = stro.ras{ind,spikeind};
        if lasertrials(ind)
            % Pulling out the laseron and off timings from analog traces
            t = stro.ras{ind,analogstrtimeind}+[0:1:length(stro.ras{ind,lasertraceind})-1]/samplerate;
            laserontime = t(find(stro.ras{ind,lasertraceind}>0.1,1));
            laserofftime = t(find(stro.ras{ind,lasertraceind}>0.1,1,'last'));
            firstspiketime_tmp = [firstspiketime_tmp; spiketimes(find(spiketimes - laserontime>0,1))-laserontime];
        else
        end
    end
    timeoffirstspike{jj} = firstspiketime_tmp;
    meantime = [meantime; mean(firstspiketime_tmp)];
    stdtime = [stdtime; std(firstspiketime_tmp)];
    groupID = [groupID; jj*ones(sum(lasertrials(idxs)),1)];
end

figure(plot_counter); set(gcf,'Name','First spike time');
subplot(121); boxplot(1000*cell2mat(timeoffirstspike)',groupID','PlotStyle','compact'); set(gca,'Xlim',[0 40],'XTick',0:10:40,'Ylim',[0 100],'YTick',0:25:100,'Tickdir','out'); axis square;
xlabel('Site number'); ylabel('First spike time');
subplot(122); histogram(1000*meantime,0:2:50); set(gca,'Xlim',[0 50],'XTick',0:10:50,'Ylim',[0 6],'YTick',0:3:6,'Tickdir','out'); axis square;
ylabel('# of sites'); xlabel('Average first spike time');
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;
%% Figure 3 & Supplementary figure 1
if ~exist('plot_counter')
    plot_counter = 1;
end
filename = {'A012219006.nex';'A012219009.nex';'A012219010.nex';'A012219011.nex';'A012219012.nex';'A013119002.nex';'A013119003.nex';'A013119004.nex';...
    'A013119005.nex';'A013119006.nex';'A020119004.nex';'A020119005.nex';'A020319007.nex';'A020319008.nex';'A020319009.nex';'A020519003.nex';...
    'A020519008.nex';'A020519009.nex'};
Singlevsmultiunit = ['M';'M';'M';'M';'M';'S';'M';'S';'S';'M';'S';'M';'M';'S';'M';'M';'M';'M'];
N = numel(filename);
binwidth = .005;
bins = -0.7:binwidth:1.0;
baselinebins = (bins>0 & bins<0.3);
baselineFR = [];
laserFR = [];
reboundFR = [];
prereboundbaselineFR = [];
reboundFR2 = [];
reboundFR3 = [];
reboundFR4 = [];
statsresult = [];
backtobaseline_time = [];
backtobaseline_time2 = [];
median_spikecount = cell(N,1);
mean_spikecount = cell(N,1);
count  = 1;
for jj = 1:N
    stro = nex2stro(findfile(char(filename(jj,:))));
    fponidx  = strcmp(stro.sum.trialFields(1,:),'fpon_t');
    fpacqidx  = strcmp(stro.sum.trialFields(1,:),'fpacq_t');
    fpoffidx  = strcmp(stro.sum.trialFields(1,:),'fpoff_t');
    stimonidx  = strcmp(stro.sum.trialFields(1,:),'stimon_t');
    stimoffidx  = strcmp(stro.sum.trialFields(1,:),'stimoff_t');
    targonidx  = strcmp(stro.sum.trialFields(1,:),'targon_t');
    saccstartidx  = strcmp(stro.sum.trialFields(1,:),'saccstart_t');
    saccendidx  = strcmp(stro.sum.trialFields(1,:),'saccend_t');
    rewidx  = strcmp(stro.sum.trialFields(1,:),'rew_t');
    stimpresentidx  = strcmp(stro.sum.trialFields(1,:),'stimpresent');
    Lcc = strcmp(stro.sum.trialFields(1,:),'lcc');
    Mcc = strcmp(stro.sum.trialFields(1,:),'mcc');
    Scc = strcmp(stro.sum.trialFields(1,:),'scc');
    tf = strcmp(stro.sum.trialFields(1,:),'tf');
    oog = strcmp(stro.sum.trialFields(1,:),'oog');
    stimidx = strcmp(stro.sum.trialFields(1,:),'stim_idx');
    optstim = strcmp(stro.sum.trialFields(1,:),'optstim');
    correct = strcmp(stro.sum.trialFields(1,:),'correct');
    laseron = strcmp(stro.sum.trialFields(1,:),'laseron_t');
    laseroff = strcmp(stro.sum.trialFields(1,:),'laseroff_t');
    analogstrtimeind = strcmp(stro.sum.rasterCells(1,:),'anlgStartTime');
    spikeind = strcmp(stro.sum.rasterCells(1,:),'sig001a');
    heyetraceind = strcmp(stro.sum.rasterCells(1,:),'AD11');
    veyetraceind = strcmp(stro.sum.rasterCells(1,:),'AD12');
    lasertraceind = strcmp(stro.sum.rasterCells(1,:),'AD13');
    correcttrials = stro.trial(:,correct);
    diridxs = stro.trial(:,stimidx);
    stimpresent = logical(stro.trial(:,stimpresentidx));
    LMStriplet = [stro.trial(:,Lcc) stro.trial(:,Mcc) stro.trial(:,Scc)];
    lasertrials = logical(stro.trial(:,optstim));
    samplerate = stro.sum.analog.storeRates{3}; % sample rate at which laser analog pulses are stored in file
    
    % Stimulus absent trials
    idxs = find(~stimpresent);
    count1 = 1;
    PSTHlaser = [];
    PSTHbaseline =[];
    tmp_baselineFR = [];
    tmp_laserFR = [];
    tmp_reboundFR = [];
    tmp_reboundFR2 = [];
    tmp_reboundFR3 = [];
    tmp_reboundFR4 = [];
    tmp_spikecountsbaseline = [];
    tmp_spikecountslaserstimabsent = [];
    tmp_prereboundbaselineFR = [];
    
    timediff = [];
    for ii = 1:numel(idxs)
        ind = idxs(ii);
        analogstartime = stro.ras{ind,analogstrtimeind};
        if lasertrials(ind)
            t = stro.ras{ind,analogstrtimeind}+[0:1:length(stro.ras{ind,lasertraceind})-1]/samplerate;
            timediff = [timediff; stro.trial(ind,saccstartidx)-t(find(stro.ras{ind,lasertraceind}>0.1,1))];
        end
    end
    [~,~,count2] = unique(timediff);
    
    for ii = 1:numel(idxs)
        ind = idxs(ii);
        analogstartime = stro.ras{ind,analogstrtimeind};
        spiketimes = stro.ras{ind,spikeind};
        if lasertrials(ind)
            % Pulling out the laseron and off timings from analog traces
            t = stro.ras{ind,analogstrtimeind}+[0:1:length(stro.ras{ind,lasertraceind})-1]/samplerate;
            laserontime = t(find(stro.ras{ind,lasertraceind}>0.1,1));
            laserofftime = t(find(stro.ras{ind,lasertraceind}>0.1,1,'last'));
            timedurlaser = laserofftime - laserontime;
            saccadestarttime = stro.trial(ind,saccstartidx)-laserontime;
            if jj == 14 | jj == 11
                figure(plot_counter); subplot(3,2,2*count-1); plot(spiketimes(spiketimes-laserontime<max(bins) & spiketimes-laserontime>min(bins))-laserontime,count2(count1)*ones(size(spiketimes(spiketimes-laserontime<max(bins) & spiketimes-laserontime>min(bins)))),'k.'); hold on;
                plot(saccadestarttime,count2(count1),'v','MarkerSize',4,'color',[1 0 0]);
                plot(stro.trial(ind,fpacqidx)-laserontime,count2(count1),'*','MarkerSize',4,'color',[1 0 0]); axis square;
            end
            PSTHlaser = [PSTHlaser; hist(spiketimes(spiketimes-laserontime<max(bins) & spiketimes-laserontime>min(bins))-laserontime, bins)/binwidth];
            count1 = count1 + 1;
            
            %             laserbins = bins>=0 & bins<=timedurlaser;
            tmp_laserFR = [tmp_laserFR; sum(spiketimes-laserontime>0 & spiketimes-laserontime<=timedurlaser)/timedurlaser];
            tmp_spikecountsbaseline = [tmp_spikecountsbaseline; sum(spiketimes>laserontime-0.3 & spiketimes<laserontime)];
            tmp_spikecountslaserstimabsent = [tmp_spikecountslaserstimabsent; sum(spiketimes>laserontime & spiketimes<laserofftime)];
            
            %             reboundbins = bins>0.30 & bins<=0.35;
            tmp_reboundFR = [tmp_reboundFR; sum(spiketimes-laserontime>=0.30 & spiketimes-laserontime<0.35)/0.050];
            tmp_reboundFR2 = [tmp_reboundFR2; sum(spiketimes-laserontime>=0.35 & spiketimes-laserontime<0.40)/0.050];
            tmp_reboundFR3 = [tmp_reboundFR3; sum(spiketimes-laserontime>=0.40 & spiketimes-laserontime<0.45)/0.050];
            tmp_reboundFR4 = [tmp_reboundFR4; sum(spiketimes-laserontime>=0.45 & spiketimes-laserontime<0.50)/0.050];
            
            %             prereboundbins = bins>=-0.05 & bins< 0;
            tmp_prereboundbaselineFR = [tmp_prereboundbaselineFR; sum(spiketimes-laserontime>=-0.05 & spiketimes-laserontime<0.0)/0.050 ];
        else
            if ~stimpresent(ind)
                stimontime = stro.trial(ind,stimonidx);
                PSTHbaseline = [PSTHbaseline;  hist(spiketimes-stimontime, bins)/binwidth];
                tmp_baselineFR = [tmp_baselineFR; sum(spiketimes>stimontime & spiketimes<=stimontime+0.3)/0.3];
            end
        end
    end
    timewindow = 0.05;
    reference_baselinebin = bins>-1*timewindow & bins<0;
    t = 0.30:binwidth:0.70;
    p = [];
    tmp_median_spikecount = [];
    tmp_mean_spikecount = [];
    
    meanbaseline_spikecount = mean(sum(PSTHlaser(:,reference_baselinebin)*binwidth,2));
    for ii = 1:numel(t)
        t1 = t(ii);
        [p_tmp,h_tmp] = signrank(sum(PSTHlaser(:,reference_baselinebin)*binwidth,2),sum(PSTHlaser(:,bins>t1 & bins<t1+timewindow)*binwidth,2));
        tmp_mean_spikecount = [tmp_mean_spikecount; mean(sum(PSTHlaser(:,bins>t1 & bins<t1+timewindow)*binwidth,2))];
        p = [p; log10(p_tmp)];
        
    end
    t = t + (timewindow/2);
    backtobaseline_time = [backtobaseline_time; t(find(p>log10(0.05),1))];
    
    tmp_time = t(find(mean(PSTHlaser(:,bins>=0.3 & bins<=0.7),1)>0.9*mean(tmp_prereboundbaselineFR),1))-(timewindow/2);
    backtobaseline_time2 = [backtobaseline_time2; tmp_time];
    
    if jj == 14 | jj == 11
        subplot(3,2,2*count-1); set(gca,'Ylim',[0 count1],'Xlim',[-0.3 0.6],'XTick',-0.3:0.3:0.6,'Tickdir','out','YTick',0:15:30); axis square;
        line([0 0],[0 count1]); line([timedurlaser timedurlaser],[0 count1]); line([tmp_time tmp_time],[0 count1],'color',[0 0 0]); hold off;
        subplot(3,2,2*count); plot(bins,mean(PSTHlaser,1),'-','color',[0 0.5 1],'Linewidth',2);
        set(gca,'Ylim',[0 150],'Xlim',[-0.3 0.6],'XTick',-0.3:0.3:0.6,'Tickdir','out','YTick',0:50:150); axis square;
        line([0 0],[0 150]); line([timedurlaser timedurlaser],[0 150]); line([tmp_time tmp_time],[0 150],'color',[0 0 0]); hold off;
        count = count + 1;
    end
    baselineFR = [baselineFR; mean(tmp_baselineFR)];
    laserFR = [laserFR; mean(tmp_laserFR)];
    reboundFR = [reboundFR; mean(tmp_reboundFR)];
    prereboundbaselineFR = [prereboundbaselineFR; mean(tmp_prereboundbaselineFR)]; % 0-50 ms after laser was swicthed was off
    reboundFR2 = [reboundFR2; mean(tmp_reboundFR2)]; % 100-150 ms after the laser was swiched off
    reboundFR3 = [reboundFR3; mean(tmp_reboundFR3)]; % 100-150 ms after the laser was swiched off
    reboundFR4 = [reboundFR4; mean(tmp_reboundFR4)]; % 100-150 ms after the laser was swiched off
    statsresult = [statsresult; signrank(tmp_reboundFR,tmp_prereboundbaselineFR)];
    mean_spikecount{jj} = tmp_mean_spikecount';
    
end
subplot(325); histogram((backtobaseline_time2 - 0.3)*1000,0:25:250,'Normalization','probability','FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
hold on; plot(median((backtobaseline_time2 - 0.3)*1000),0,'kv','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'XTick',0:50:250,'YTick',0:0.15:0.3,'Tickdir','out','Ylim',[0 0.3],'Xlim',[0 250]);
xlabel('time to recover (ms)'); ylabel('Proportion of cells'); axis square

idx = find(logical(t-0.475)==0);
for ii=1:N
    FR = mean_spikecount{ii}/0.050;
    subplot(326); plot(prereboundbaselineFR(ii),mean(FR(1:idx)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
    plot(repmat(prereboundbaselineFR(ii),[idx 1]),FR(1:idx),'k');
end
set(gca,'Tickdir','out','Xlim',[0 34],'Ylim',[0 34],'XTick',0:17:34,'YTick',0:17:34); line([0 34],[0 34]); xlabel('pre-laser FR'); ylabel('post-laser FR'); axis square; hold off;
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;


figure(plot_counter);
subplot(221);plot(prereboundbaselineFR,reboundFR,'o','MarkerSize',6,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
xlabel('Pre laser FR'); ylabel('Post laser FR'); set(gca,'Tickdir','out','Xlim',[0 32],'Ylim',[0 32],'XTick',0:16:32,'YTick',0:16:32);
line([0 32],[0 32]); axis square; title('0-50 ms');
subplot(222);plot(prereboundbaselineFR,reboundFR2,'o','MarkerSize',6,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
xlabel('Pre laser FR'); ylabel('Post laser FR'); set(gca,'Tickdir','out','Xlim',[0 32],'Ylim',[0 32],'XTick',0:16:32,'YTick',0:16:32);
line([0 32],[0 32]); axis square; title('50-100 ms');
subplot(223);plot(prereboundbaselineFR,reboundFR3,'o','MarkerSize',6,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
xlabel('Pre laser FR'); ylabel('Post laser FR'); set(gca,'Tickdir','out','Xlim',[0 32],'Ylim',[0 32],'XTick',0:16:32,'YTick',0:16:32);
line([0 32],[0 32]); axis square; title('100-150 ms');
subplot(224);plot(prereboundbaselineFR,reboundFR4,'o','MarkerSize',6,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
xlabel('Pre laser FR'); ylabel('Post laser FR'); set(gca,'Tickdir','out','Xlim',[0 32],'Ylim',[0 32],'XTick',0:16:32,'YTick',0:16:32);
line([0 32],[0 32]); axis square;  title('150-200 ms');
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;

% Additional analyses for estimating recovery time
figure(plot_counter); set(gcf,'Name','Estimating recovery time');
subplot(221); plot(t,cell2mat(mean_spikecount)'/0.050); set(gca,'Tickdir','out','Ylim',[0 40],'YTick',0:10:40,'Xlim',[0.3 0.5],'XTick',0.3:0.1:0.5); xlabel('time'); ylabel('FR'); axis square
subplot(222); plot(t,repmat(prereboundbaselineFR,[1 numel(t)])'); set(gca,'Tickdir','out','Ylim',[0 40],'YTick',0:10:40,'Xlim',[0.3 0.5],'XTick',0.3:0.1:0.5); xlabel('time'); ylabel('FR'); axis square;
idx = find(logical(t-0.475)==0);

for ii=1:N
    FR = mean_spikecount{ii}/0.050;
    subplot(223); plot(prereboundbaselineFR(ii),FR(1),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
    plot(prereboundbaselineFR(ii),FR(idx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
    plot(repmat(prereboundbaselineFR(ii),[idx-2 1]),FR(2:idx-1),'k');
    subplot(224); plot(prereboundbaselineFR(ii),mean(FR(1:idx)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
    plot(repmat(prereboundbaselineFR(ii),[idx 1]),FR(1:idx),'k');
end
subplot(223); plot(0.5:0.5:32,(0.5:0.5:32)*0.9,'Color',[0.5 0.5 0.5]); plot(0.5:0.5:32,(0.5:0.5:32)*1.1,'Color',[0.5 0.5 0.5]);
set(gca,'Tickdir','out','Xlim',[0 32],'Ylim',[0 32],'XTick',0:16:32,'YTick',0:16:32); line([0 32],[0 32]); xlabel('pre-laser FR'); ylabel('post-laser FR'); axis square; hold off;
subplot(224); plot(0.5:0.5:32,(0.5:0.5:32)*0.9,'Color',[0.5 0.5 0.5]); plot(0.5:0.5:32,(0.5:0.5:32)*1.1,'Color',[0.5 0.5 0.5]);
set(gca,'Tickdir','out','Xlim',[0 32],'Ylim',[0 32],'XTick',0:16:32,'YTick',0:16:32); line([0 32],[0 32]); xlabel('pre-laser FR'); ylabel('post-laser FR'); axis square; hold off;
plot_counter = plot_counter + 1;

% Representing the recovery time (calculated by Greg's method)
backtobaseline_time2 = backtobaseline_time2 - 0.3;
figure(plot_counter); set(gcf,'Name','Recovery time 2');
subplot(131);histogram(backtobaseline_time2*1000,0:20:250,'Normalization','probability','FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
hold on; plot(median(backtobaseline_time2*1000),0,'kv','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'XTick',0:50:250,'YTick',0:0.25:0.5,'Tickdir','out','Ylim',[0 0.5],'Xlim',[0 250]);
xlabel('time to recover (ms): method 2'); ylabel('Proportion of cells'); axis square
subplot(132);histogram((backtobaseline_time-0.3)*1000,0:20:250,'Normalization','probability','FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
hold on; plot(median((backtobaseline_time-0.3)*1000),0,'kv','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'XTick',0:50:250,'YTick',0:0.25:0.5,'Tickdir','out','Ylim',[0 0.5],'Xlim',[0 250]);
xlabel('time to recover (ms):method 1'); ylabel('Proportion of cells'); axis square

[r,p1] = corr(backtobaseline_time-0.3,backtobaseline_time2);
subplot(133); plot(backtobaseline_time-0.3,backtobaseline_time2,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
xlabel('Recovery time 1'); ylabel('Recovery time 2'); set(gca,'Tickdir','out','Xlim',[0 0.3],'Ylim',[0 0.3],'Xtick',0:0.15:0.3,'Ytick',0:0.15:0.3); axis square;
plot_counter = plot_counter + 1;

%% Figure 4
% � Schematic of the SCstimcue task
% �	1 example of SCstimcue session, each from Apollo and Maui along with catch trials (4 figures)
% �	A population figure of probability of saccades in laser vs control trials for �within RF� and �outside RF� from Maui
% (2 figures, Apollo�s population goes in supplementary figure) (2-D histogram or heat-map)

if ~exist('plot_counter')
    plot_counter = 1;
end
filename = ['M042718005.nex';'A012319003.nex'];
figure(plot_counter); set(gcf,'Name','SC stimcue');
for jj = 1:size(filename,1)
    stro = nex2stro(findfile(filename(jj,:)));
    targ_x = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_x'));
    
    targ_y = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_y'));
    uniquetargxy = unique([targ_x targ_y],'rows');
    ntrials = size(stro.trial,1);
    Lcatchtrials = targ_x == 0 & targ_y == 0;
    laseron_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
    laseroff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
    fpon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpon_t'));
    fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
    fix_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpacq_t'));
    targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targon_t'));
    targoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targoff_t'));
    sacinit_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'saccinit_t'));
    sacmade_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'saccmade_t'));
    targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targon_t'));
    targoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targoff_t'));
    
    Llaser = ~isnan(laseron_t);
    samplerate = stro.sum.analog.storeRates{1};
    subplot(2,2,2*jj-1); plot(stro.sum.exptParams.rf_x/10,stro.sum.exptParams.rf_y/10,'o','MarkerSize',12,'MarkerEdgeColor',[0 0 0],'Linewidth',2); hold on;
    subplot(2,2,2*jj); plot(stro.sum.exptParams.rf_x/10,stro.sum.exptParams.rf_y/10,'o','MarkerSize',12,'MarkerEdgeColor',[0 0 0],'Linewidth',2); hold on;
    RFloc = [stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y];
    if norm(RFloc)>100
        plotlimit = 20;
    else
        plotlimit = 10;
    end
    
    for i = 1:size(stro.trial,1)
        if size(stro.sum.rasterCells,2)==5
            x = stro.ras{i,2}*4096/400;
            y = stro.ras{i,3}*4096/400;
            t = stro.ras{i,5}+[0:1:length(x)-1]/samplerate;
        elseif size(stro.sum.rasterCells,2)==6
            x = stro.ras{i,3}*4096/400;
            y = stro.ras{i,4}*4096/400;
            t = stro.ras{i,6}+[0:1:length(x)-1]/samplerate;
        elseif size(stro.sum.rasterCells,2)==7
            x = stro.ras{i,4}*4096/400;
            y = stro.ras{i,5}*4096/400;
            t = stro.ras{i,7}+[0:1:length(x)-1]/samplerate;
        end
        Lt = t>fpoff_t(i) & t < fpoff_t(i) +.3;
        if ~Lcatchtrials(i)
            if ~Llaser(i)
                subplot(2,2,2*jj-1); plot(x(Lt),y(Lt),'color',[0.5 0.5 0.5],'Linewidth',1); hold on;
            else
                subplot(2,2,2*jj-1); plot(x(Lt),y(Lt),'color',[0 0.5 1],'Linewidth',1); hold on;
            end
        else
            if ~Llaser(i)
                subplot(2,2,2*jj); plot(x(Lt),y(Lt),'--','color',[0.5 0.5 0.5],'Linewidth',1); hold on;
            else
                subplot(2,2,2*jj); plot(x(Lt),y(Lt),'--','color',[0 0.5 1],'Linewidth',1); hold on;
            end
        end
    end
    subplot(2,2,2*jj-1); axis square; set(gca,'Xlim',[-1*plotlimit plotlimit],'Ylim',[-1*plotlimit plotlimit],'Tickdir','out','XTick',[-1*plotlimit:5:plotlimit],'YTick',[-1*plotlimit:5:plotlimit]); title(filename(jj,1));
    subplot(2,2,2*jj); axis square; set(gca,'Xlim',[-1*plotlimit plotlimit],'Ylim',[-1*plotlimit plotlimit],'Tickdir','out','XTick',[-1*plotlimit:5:plotlimit],'YTick',[-1*plotlimit:5:plotlimit]); title(strcat(filename(jj,1),':catch trials'));
end
plot_counter = plot_counter + 1;

% Need to calculate timing of events
idx = find(Llaser,1);
x = stro.ras{idx,2}*4096/400; x = x/max(abs(x));
t = stro.ras{idx,5}+[0:1:length(x)-1]/samplerate;
FP = zeros(numel(t),1); FP(t>fpon_t(idx) & t<fpoff_t(idx)) = 1;
Targ = zeros(numel(t),1); Targ(t>targon_t(idx) & t<targoff_t(idx)) = 1;
L = stro.ras{idx,4}*4096/400; L = L/max(abs(L));
figure(plot_counter);
plot(t,FP+10); hold on; plot(t,Targ+7); plot(t,L+4); plot(t,1+x/max(abs(x))); set(gca,'Ylim',[0 12],'Xtick',[],'YTick',[]); axis square;
plot_counter = plot_counter + 1;



%% Figure 5
% �	Schematic of DToneloc task (1 figure)
% �	1 example of DToneloc session (Hits, Misses, CR, FA) from Maui and Apollo monkey along with laser and control staircases (4 figures)
if ~exist('plot_counter')
    plot_counter = 1;
end
load('T_vos1978_Y');
Vlambda = T_vos1978_Y';
HitL = []; HitNL = [];
MissL = []; MissNL = [];
CRL = []; CRNL = [];
FAL = []; FANL = [];
TPRNL = []; TPRL = [];
FPRNL = []; FPRL = [];
alpha = 0.05;
RF = [];
stimabsentL = []; stimabsentNL = [];
stimpresentL = []; stimpresentNL = [];
filename = ['M051418010.nex'; 'A021419009.nex'];
RWmultiplier = [0.7; 0.6];
laserdial = [4.0; 1.8];
count = 1;
figure(plot_counter); set(gcf,'Name','Example Sessions from M & A');
for jj = 1:size(filename,1)
    stro = nex2stro(findfile(char(filename(jj,:))));
    fponidx  = strcmp(stro.sum.trialFields(1,:),'fpon_t');
    fpacqidx  = strcmp(stro.sum.trialFields(1,:),'fpacq_t');
    fpoffidx  = strcmp(stro.sum.trialFields(1,:),'fpoff_t');
    stimonidx  = strcmp(stro.sum.trialFields(1,:),'stimon_t');
    stimoffidx  = strcmp(stro.sum.trialFields(1,:),'stimoff_t');
    targonidx  = strcmp(stro.sum.trialFields(1,:),'targon_t');
    saccstartidx  = strcmp(stro.sum.trialFields(1,:),'saccstart_t');
    saccendidx  = strcmp(stro.sum.trialFields(1,:),'saccend_t');
    rewidx  = strcmp(stro.sum.trialFields(1,:),'rew_t');
    stimpresentidx  = strcmp(stro.sum.trialFields(1,:),'stimpresent');
    Lcc = strcmp(stro.sum.trialFields(1,:),'lcc');
    Mcc = strcmp(stro.sum.trialFields(1,:),'mcc');
    Scc = strcmp(stro.sum.trialFields(1,:),'scc');
    tf = strcmp(stro.sum.trialFields(1,:),'tf');
    oog = strcmp(stro.sum.trialFields(1,:),'oog');
    stimidx = strcmp(stro.sum.trialFields(1,:),'stim_idx');
    optstim = strcmp(stro.sum.trialFields(1,:),'optstim');
    correct = strcmp(stro.sum.trialFields(1,:),'correct');
    laseron = strcmp(stro.sum.trialFields(1,:),'laseron_t');
    laseroff = strcmp(stro.sum.trialFields(1,:),'laseroff_t');
    correcttrials = stro.trial(:,correct);
    diridxs = stro.trial(:,stimidx);
    stimpresent = logical(stro.trial(:,stimpresentidx));
    LMStriplet = [stro.trial(:,Lcc) stro.trial(:,Mcc) stro.trial(:,Scc)];
    lasertrials = logical(stro.trial(:,optstim));
    
    lasertrialidxs = logical(stro.trial(:,optstim));
    colordirchoiceidxs = correcttrials;
    colorstimpresent = logical(stimpresent);
    percentcorrectlasertrials = sum(colordirchoiceidxs & lasertrialidxs)/sum(lasertrialidxs);
    percentcorrectnonlasertrials = sum(colordirchoiceidxs & ~lasertrialidxs)/sum(~lasertrialidxs);
    
    % Laser trials
    Hitlasertrial = sum(colordirchoiceidxs & lasertrialidxs & colorstimpresent); % Hit
    Misslasertrial = sum(~colordirchoiceidxs & lasertrialidxs & colorstimpresent); % Miss
    CRlasertrial = sum(colordirchoiceidxs & lasertrialidxs & ~colorstimpresent); % Correct Reject
    FAlasertrial = sum(~colordirchoiceidxs & lasertrialidxs & ~colorstimpresent); % False Alarm
    CRlasertrialidxs = find(lasertrialidxs & ~colorstimpresent); L1 = numel(CRlasertrialidxs);
    
    % Non-Laser trials
    Hitnonlasertrial = sum(colordirchoiceidxs & ~lasertrialidxs & colorstimpresent); % Hit
    Missnonlasertrial = sum(~colordirchoiceidxs & ~lasertrialidxs & colorstimpresent); % Miss
    CRnonlasertrial = sum(colordirchoiceidxs & ~lasertrialidxs & ~colorstimpresent); % Correct Reject
    FAnonlasertrial = sum(~colordirchoiceidxs & ~lasertrialidxs & ~colorstimpresent); % False Alarm
    CRnonlasertrialidxs = find(~lasertrialidxs & ~colorstimpresent); L2 = numel(CRnonlasertrialidxs);
    
    stimabsentL = [stimabsentL; sum(CRlasertrial) + sum(FAlasertrial)];
    stimabsentNL = [stimabsentNL; sum(CRnonlasertrial) + sum(FAnonlasertrial)];
    stimpresentL = [stimpresentL; sum(Hitlasertrial) + sum(Misslasertrial)];
    stimpresentNL = [stimpresentNL; sum(Hitnonlasertrial) + sum(Missnonlasertrial)];
    
    % Now storing the Hits, Miss, CR and FA for laser and non-laser trials
    HitL = [HitL; Hitlasertrial]; HitNL = [HitNL; Hitnonlasertrial];
    MissL = [MissL; Misslasertrial]; MissNL = [MissNL; Missnonlasertrial];
    CRL = [CRL; CRlasertrial]; CRNL = [CRNL; CRnonlasertrial];
    FAL = [FAL; FAlasertrial]; FANL = [FANL; FAnonlasertrial];
    TPRNL = [TPRNL; Hitnonlasertrial/(Hitnonlasertrial+Missnonlasertrial)]; % True positive ratio, non-laser trial
    FPRNL = [FPRNL; FAnonlasertrial/(FAnonlasertrial+CRnonlasertrial)]; % False positive ratio, non-laser trial
    TPRL = [TPRL; Hitlasertrial/(Hitlasertrial+Misslasertrial)]; % True positive ratio, laser trial
    FPRL = [FPRL; FAlasertrial/(FAlasertrial+CRlasertrial)]; % False positive ratio, laser trial
    
    % Calculating psychophysical detection threshold for non-laser trials
    
    fundamentals = stro.sum.exptParams.fundamentals; % CONE FUNDAMENTALS: L,M,S
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
    mon_spd = stro.sum.exptParams.mon_spd; % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
    M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals
    Lb = M*stro.sum.exptParams.bkgndrgb;
    % Converting all the contrasts to Luminance contrast
    Bkgnd_luminance = stro.sum.exptParams.bkgndrgb'*mon_spd*Vlambda;
    Peak_luminance = (inv(M)*((1+LMStriplet).*(repmat(Lb',[size(LMStriplet,1) 1])))')'*mon_spd*Vlambda;
    Trough_luminance = (inv(M)*((1-LMStriplet).*(repmat(Lb',[size(LMStriplet,1) 1])))')'*mon_spd*Vlambda;
    Contrast_luminance = (Peak_luminance - Trough_luminance)./(Peak_luminance + Trough_luminance);
    LMStripletfile = Contrast_luminance;
    LMScontrastfileNL = Contrast_luminance(colorstimpresent & ~lasertrialidxs,:);
    LMScontrastfileL = Contrast_luminance(colorstimpresent & lasertrialidxs,:);
    % calculating gamut edge luminance contrast
    t = min(((1./stro.sum.exptParams.bkgndrgb)-1));
    OOG_luminance = (stro.sum.exptParams.bkgndrgb*(1+t))' * mon_spd*Vlambda;
    putativeoogcontrast =  (OOG_luminance - Bkgnd_luminance)/Bkgnd_luminance;
    
    
    answers = colordirchoiceidxs(colorstimpresent & ~lasertrialidxs);
    LMScontrastfileNL(LMScontrastfileNL>putativeoogcontrast) = putativeoogcontrast;
    LMScontrastfileL(LMScontrastfileL>putativeoogcontrast) = putativeoogcontrast;
    
    % Calculating psychophysical detection threshold for non-laser trials
    trials = lasertrials;
    % only evaluating the laser trials
    weibullparamsL = [];
    weibullparamsNL = [];
    oogvals = putativeoogcontrast;
    for kk = 1:2
        if kk == 2
            trials = ~trials;
        end
        LMScontrast = [];
        LMScontrast = Contrast_luminance(stimpresent & trials,:);
        answers = correcttrials(stimpresent & trials);
        LMScontrast(LMScontrast>putativeoogcontrast) = putativeoogcontrast;
        
        if kk==1
            contrastL = unique(LMScontrast);
            correctanswersL = zeros(size(contrastL));
            wronganswersL = zeros(size(contrastL));
            trialspercontrastL = zeros(size(contrastL));
            for ss = 1:numel(contrastL)
                trialspercontrastL(ss) = numel(answers(LMScontrast==contrastL(ss)));
                correctanswersL(ss) = sum(answers(LMScontrast==contrastL(ss)));
                wronganswersL(ss) = trialspercontrastL(ss) - correctanswersL(ss);
                percorrectL(ss) = correctanswersL(ss)/trialspercontrastL(ss);
            end
            [aL,bL,gL] = weibullFitforDToneloc(contrastL,[correctanswersL wronganswersL],'mle');
            weibullparamsL = [weibullparamsL; aL bL gL]; % laser trials
            LMScontrastL = LMScontrast;
        else
            contrastNL = unique(LMScontrast);
            correctanswersNL = zeros(size(contrastNL));
            wronganswersNL = zeros(size(contrastNL));
            trialspercontrastNL = zeros(size(contrastNL));
            for ss = 1:numel(contrastNL)
                trialspercontrastNL(ss) = numel(answers(LMScontrast==contrastNL(ss)));
                correctanswersNL(ss) = sum(answers(LMScontrast==contrastNL(ss)));
                wronganswersNL(ss) = trialspercontrastNL(ss) - correctanswersNL(ss);
                percorrectNL(ss) = correctanswersNL(ss)/trialspercontrastNL(ss);
            end
            [aNL,bNL,gNL] = weibullFitforDToneloc(contrastNL,[correctanswersNL wronganswersNL],'mle');
            weibullparamsNL = [weibullparamsNL; aNL bNL gNL]; % laser trials
            LMScontrastNL = LMScontrast;
        end
    end
    
    subplot(2,3,3*jj-2); h = bar([HitNL(jj)/stimpresentNL(jj) HitL(jj)/stimpresentL(jj);MissNL(jj)/stimpresentNL(jj) MissL(jj)/stimpresentL(jj);CRNL(jj)/stimabsentNL(jj) CRL(jj)/stimabsentL(jj);FANL(jj)/stimabsentNL(jj)  FAL(jj)/stimabsentL(jj)]); set(h(2),'FaceColor',[0 0.5 1]); set(h(1),'FaceColor',[0.5 0.5 0.5]);
    set(gca,'XTick',[1 2 3 4],'XTickLabel',{'H','M','CR','FA'},'TickDir','out','Xlim',[0 5],'YTick',0:0.25:1.0); ylabel('Prop of trials'); axis square; title(filename(jj,1));
    subplot(2,3,3*jj-1); plot(LMScontrastfileNL,'color',[0.5 0.5 0.5],'Linewidth',2); hold on; plot(LMScontrastfileL,'color',[0 0.5 1.0],'Linewidth',2); xlabel('trial number'); ylabel('Luminance contrast');
    title(filename(jj,1)); axis square; set(gca,'TickDir','out','Xlim',[0 30],'XTick',[0:15:30]); hold off;
    
    contrastlattice = logspace(log10(0.03),log(1.0),51);
    fitL = gL*(1-exp(-((contrastlattice./aL).^bL)));
    fitNL = gNL*(1-exp(-((contrastlattice./aNL).^bNL)));
    [~,edges] = discretize(contrastNL,logspace(log10(min(contrastNL)-0.0001),log10(max(contrastL)+0.0001),7));
    for ss = 1:numel(edges)-1
        idx = contrastNL>=edges(ss) & contrastNL<edges(ss+1);
        if sum(idx)
            subplot(2,3,3*jj);  plot(edges(ss),sum(correctanswersNL(idx))./(sum(correctanswersNL(idx))+sum(wronganswersNL(idx))),'o','MarkerSize',6,'MarkerFacecolor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
        end
    end
    plot(contrastlattice,fitNL,'-','Linewidth',2,'color',[0.5 0.5 0.5]); hold on;
    for ss = 1:numel(edges)-1
        idx = contrastL>=edges(ss) & contrastL<edges(ss+1);
        if sum(idx)
            subplot(2,3,3*jj);  plot(edges(ss),sum(correctanswersL(idx))./sum((correctanswersL(idx))+sum(wronganswersL(idx))),'o','MarkerSize',6,'MarkerFacecolor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]); hold on;
        end
    end
    plot(contrastlattice,fitL,'-','Linewidth',2,'color',[0 0.5 1.0]); set(gca,'Xlim',[0.03 1],'Tickdir','out','YTick',[0:0.25:1],'Xscale','log','XTick',[0.1 0.3 1],'XTickLabels', {'0.1','0.3','1'}); xlabel('Luminance contrast'); ylabel('Proportion correct'); title('Psychometric function'); axis square; hold off;
    
    [~,p] = equalproptest([HitL(jj) HitNL(jj)],[30 30],0.05);
    disp(p);
end
plot_counter = plot_counter + 1;

% Need to calculate timing of events
idx = 5;%find(lasertrials,1);
x = stro.ras{idx,2}*4096/400;
t = stro.ras{idx,5}+[0:1:length(x)-1]/samplerate;
FP = zeros(numel(t),1); FP(t>stro.trial(idx,fponidx) & t<stro.trial(idx,fpoffidx)) = 1;
Targ = zeros(numel(t),1); Targ(t>stro.trial(idx,targonidx)) = 1;
Stim = zeros(numel(t),1); Stim(t>stro.trial(idx,stimonidx) & t<stro.trial(idx,stimoffidx)) = 1;
L = stro.ras{idx,4}*4096/400; L = L/max(abs(L));
figure(plot_counter);
plot(t,FP+13); hold on; plot(t,Stim+10); plot(t,Targ+7); plot(t,L+4); plot(t,1+x/max(abs(x))); set(gca,'Xlim',[t(1) t(end)],'Ylim',[0 15],'Xtick',[],'YTick',[]); axis square;
plot_counter = plot_counter + 1;

%% Figure 6
% �	DToneloc population: A scatter plot of hits in control vs laser trials at tested locations for Maui and Apollo
% (Total 4 figures: 2 hits control vs laser population plot, 2 mapped RF locations for those session)
if ~exist('plot_counter')
    plot_counter = 1;
end
load('T_vos1978_Y');
Vlambda = T_vos1978_Y';
ratio_controloverlaser = [];
timeperblock = [];
trialsperblock = [];
for mm = 1:1
    if mm == 1
        load filenameoptoM.mat
        load RWmultiplieroptoM.mat
        load laserdialoptoM.mat
        monkeyname = 'Maui';
    else
        load filenameoptoA.mat
        load RWmultiplieroptoA.mat
        load laserdialoptoA.mat
        monkeyname = 'Apollo';
    end
    filename = [filenameopto];
    RWmultiplier = [RWmultiplieropto];
    laserdial = [laserdialopto];
    HitL = []; HitNL = [];
    MissL = []; MissNL = [];
    CRL = []; CRNL = [];
    numdivisionsperfile = 3;
    CRearlylateL = []; % for storing the number of CRs in first half and second half of the block of experiment (laser trials)
    CRearlylateNL = [];% for storing the number of CRs in first half and second half of the block of experiment (no-laser trials)
    FAL = []; FANL = [];
    stimdur = [];
    laserdur = [];
    h_hit = []; p_hit = [];
    h_FA = []; p_FA = [];
    TPRNL = []; TPRL = [];
    FPRNL = []; FPRL = [];
    alpha = 0.05;
    RF = [];
    stimabsentL = []; stimabsentNL = [];
    stimpresentL = []; stimpresentNL = [];
    newfilename = [];
    newRWmultiplier = [];
    newstimsize = [];
    newlaserdial = [];
    weibullparamsL = [];
    weibullparamsNL = [];
    prop_correctatGAMUTL = [];
    prop_correctatGAMUTNL = [];
    DetectThreshL = [];
    oogvals = [];
    count = 1;
    for jj = 1:size(filename,1)
        stro = nex2stro(findfile(char(filename(jj,:))));
        fponidx  = strcmp(stro.sum.trialFields(1,:),'fpon_t');
        fpacqidx  = strcmp(stro.sum.trialFields(1,:),'fpacq_t');
        fpoffidx  = strcmp(stro.sum.trialFields(1,:),'fpoff_t');
        stimonidx  = strcmp(stro.sum.trialFields(1,:),'stimon_t');
        stimoffidx  = strcmp(stro.sum.trialFields(1,:),'stimoff_t');
        targonidx  = strcmp(stro.sum.trialFields(1,:),'targon_t');
        saccstartidx  = strcmp(stro.sum.trialFields(1,:),'saccstart_t');
        saccendidx  = strcmp(stro.sum.trialFields(1,:),'saccend_t');
        rewidx  = strcmp(stro.sum.trialFields(1,:),'rew_t');
        stimpresentidx  = strcmp(stro.sum.trialFields(1,:),'stimpresent');
        Lcc = strcmp(stro.sum.trialFields(1,:),'lcc');
        Mcc = strcmp(stro.sum.trialFields(1,:),'mcc');
        Scc = strcmp(stro.sum.trialFields(1,:),'scc');
        tf = strcmp(stro.sum.trialFields(1,:),'tf');
        oog = strcmp(stro.sum.trialFields(1,:),'oog');
        stimidx = strcmp(stro.sum.trialFields(1,:),'stim_idx');
        optstim = strcmp(stro.sum.trialFields(1,:),'optstim');
        correct = strcmp(stro.sum.trialFields(1,:),'correct');
        laseron = strcmp(stro.sum.trialFields(1,:),'laseron_t');
        laseroff = strcmp(stro.sum.trialFields(1,:),'laseroff_t');
        anlgStartTime = strcmp(stro.sum.rasterCells,'anlgStartTime');
        timeperblock = [timeperblock; stro.ras{end,anlgStartTime}-stro.ras{1,anlgStartTime}];
        trialsperblock = [trialsperblock; size(stro.ras,1)];
        
        correcttrials = stro.trial(:,correct);
        diridxs = stro.trial(:,stimidx);
        stimpresent = logical(stro.trial(:,stimpresentidx));
        LMStriplet = [stro.trial(:,Lcc) stro.trial(:,Mcc) stro.trial(:,Scc)];
        lasertrials = logical(stro.trial(:,optstim));
        if isnan(unique(stro.trial(:,21:22)))
            stimlocations = [stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y];
            RF = [RF; stimlocations];
            newfilename = [newfilename; filename(jj)];
            newRWmultiplier = [newRWmultiplier; RWmultiplier(jj)];
            newlaserdial = [newlaserdial; laserdial(jj)];
            ib = ones(size(lasertrials));
            newstimsize = [newstimsize; stro.sum.exptParams.sigma];
        else
            [stimlocations,~,ib] = unique(stro.trial(:,21:22),'rows');
            RF = [RF; stimlocations];
            L = size(stimlocations,1);
            newfilename = [newfilename; repmat(filename(jj),[L 1])];
            newRWmultiplier = [newRWmultiplier; repmat(RWmultiplier(jj),[L 1])];
            newlaserdial = [newlaserdial; repmat(laserdial(jj),[L 1])];
            newstimsize = [newstimsize; repmat(stro.sum.exptParams.sigma,[L 1])];
        end
        
        for ii = 1:numel(unique(ib))
            ind = logical(ib==ii);
            lasertrialidxs = logical(stro.trial(ind,optstim));
            colordirchoiceidxs = correcttrials(ind);
            colorstimpresent = logical(stimpresent(ind));
            percentcorrectlasertrials = sum(colordirchoiceidxs & lasertrialidxs)/sum(lasertrialidxs);
            percentcorrectnonlasertrials = sum(colordirchoiceidxs & ~lasertrialidxs)/sum(~lasertrialidxs);
            
            % Laser trials
            Hitlasertrial = sum(colordirchoiceidxs & lasertrialidxs & colorstimpresent); % Hit
            Misslasertrial = sum(~colordirchoiceidxs & lasertrialidxs & colorstimpresent); % Miss
            CRlasertrial = sum(colordirchoiceidxs & lasertrialidxs & ~colorstimpresent); % Correct Reject
            FAlasertrial = sum(~colordirchoiceidxs & lasertrialidxs & ~colorstimpresent); % False Alarm
            CRlasertrialidxs = find(lasertrialidxs & ~colorstimpresent); L1 = numel(CRlasertrialidxs);
            
            % Non-Laser trials
            Hitnonlasertrial = sum(colordirchoiceidxs & ~lasertrialidxs & colorstimpresent); % Hit
            Missnonlasertrial = sum(~colordirchoiceidxs & ~lasertrialidxs & colorstimpresent); % Miss
            CRnonlasertrial = sum(colordirchoiceidxs & ~lasertrialidxs & ~colorstimpresent); % Correct Reject
            FAnonlasertrial = sum(~colordirchoiceidxs & ~lasertrialidxs & ~colorstimpresent); % False Alarm
            CRnonlasertrialidxs = find(~lasertrialidxs & ~colorstimpresent); L2 = numel(CRnonlasertrialidxs);
            
            
            stimabsentL = [stimabsentL; sum(CRlasertrial) + sum(FAlasertrial)];
            stimabsentNL = [stimabsentNL; sum(CRnonlasertrial) + sum(FAnonlasertrial)];
            stimpresentL = [stimpresentL; sum(Hitlasertrial) + sum(Misslasertrial)];
            stimpresentNL = [stimpresentNL; sum(Hitnonlasertrial) + sum(Missnonlasertrial)];
            
            % Now storing the Hits, Miss, CR and FA for laser and non-laser trials
            HitL = [HitL; Hitlasertrial]; HitNL = [HitNL; Hitnonlasertrial];
            MissL = [MissL; Misslasertrial]; MissNL = [MissNL; Missnonlasertrial];
            CRL = [CRL; CRlasertrial]; CRNL = [CRNL; CRnonlasertrial];
            FAL = [FAL; FAlasertrial]; FANL = [FANL; FAnonlasertrial];
            TPRNL = [TPRNL; Hitnonlasertrial/(Hitnonlasertrial+Missnonlasertrial)]; % True positive ratio, non-laser trial
            FPRNL = [FPRNL; FAnonlasertrial/(FAnonlasertrial+CRnonlasertrial)]; % False positive ratio, non-laser trial
            TPRL = [TPRL; Hitlasertrial/(Hitlasertrial+Misslasertrial)]; % True positive ratio, laser trial
            FPRL = [FPRL; FAlasertrial/(FAlasertrial+CRlasertrial)]; % False positive ratio, laser trial
            
            fundamentals = stro.sum.exptParams.fundamentals; % CONE FUNDAMENTALS: L,M,S
            fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
            mon_spd = stro.sum.exptParams.mon_spd; % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
            mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
            mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
            M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals
            Lb = M*stro.sum.exptParams.bkgndrgb;
            % Converting all the contrasts to Luminance contrast
            Bkgnd_luminance = stro.sum.exptParams.bkgndrgb'*mon_spd*Vlambda;
            Peak_luminance = (inv(M)*((1+LMStriplet).*(repmat(Lb',[size(LMStriplet,1) 1])))')'*mon_spd*Vlambda;
            Trough_luminance = (inv(M)*((1-LMStriplet).*(repmat(Lb',[size(LMStriplet,1) 1])))')'*mon_spd*Vlambda;
            Contrast_luminance = (Peak_luminance - Trough_luminance)./(Peak_luminance + Trough_luminance);
            %             Contrast_luminance = ((inv(M)*((LMStriplet+1).*(repmat(Lb',[size(LMStriplet,1) 1])))')'*mon_spd*Vlambda-Bkgnd_luminance)/Bkgnd_luminance; % Same as Michelson contrast
            LMScontrastfileNL = Contrast_luminance(colorstimpresent & ~lasertrialidxs,:);
            LMScontrastfileL = Contrast_luminance(colorstimpresent & lasertrialidxs,:);
            % calculating gamut edge luminance contrast
            t = min(((1./stro.sum.exptParams.bkgndrgb)-1));
            OOG_luminance = (stro.sum.exptParams.bkgndrgb*(1+t))' * mon_spd*Vlambda;
            putativeoogcontrast =  (OOG_luminance - Bkgnd_luminance)/Bkgnd_luminance;
            oogvals = [oogvals; putativeoogcontrast];
            
            % Calculating psychophysical detection threshold for non-laser trials
            LMStripletfile = LMStriplet(ind,:);
            trials = lasertrialidxs;
            
            for kk = 1:2
                if kk == 2
                    trials = ~trials;
                end
                LMScontrast = Contrast_luminance(colorstimpresent & trials & ind,:);
                answers = colordirchoiceidxs(colorstimpresent & trials);
                LMScontrast(LMScontrast>putativeoogcontrast) = putativeoogcontrast;
                if kk==1
                    contrastL = unique(LMScontrast);
                    correctanswersL = zeros(size(contrastL));
                    wronganswersL = zeros(size(contrastL));
                    trialspercontrastL = zeros(size(contrastL));
                    for ss = 1:numel(contrastL)
                        trialspercontrastL(ss) = numel(answers(LMScontrast==contrastL(ss)));
                        correctanswersL(ss) = sum(answers(LMScontrast==contrastL(ss)));
                        wronganswersL(ss) = trialspercontrastL(ss) - correctanswersL(ss);
                        percorrectL(ss) = correctanswersL(ss)/trialspercontrastL(ss);
                    end
                    [aL,bL,gL] = weibullFitforDToneloc(contrastL,[correctanswersL wronganswersL],'mle');
                    weibullparamsL = [weibullparamsL; aL bL gL]; % laser trials
                    prop_correctatGAMUTL = [prop_correctatGAMUTL; gL*(1-exp(-((putativeoogcontrast/aL).^bL)))];
                    DetectThreshL = [DetectThreshL; aL];
                else
                    contrastNL = unique(LMScontrast);
                    correctanswersNL = zeros(size(contrastNL));
                    wronganswersNL = zeros(size(contrastNL));
                    trialspercontrastNL = zeros(size(contrastNL));
                    for ss = 1:numel(contrastNL)
                        trialspercontrastNL(ss) = numel(answers(LMScontrast==contrastNL(ss)));
                        correctanswersNL(ss) = sum(answers(LMScontrast==contrastNL(ss)));
                        wronganswersNL(ss) = trialspercontrastNL(ss) - correctanswersNL(ss);
                        percorrectNL(ss) = correctanswersNL(ss)/trialspercontrastNL(ss);
                    end
                    [aNL,bNL,gNL] = weibullFitforDToneloc(contrastNL,[correctanswersNL wronganswersNL],'mle');
                    weibullparamsNL = [weibullparamsNL; aNL bNL gNL]; % laser trials
                    prop_correctatGAMUTNL = [prop_correctatGAMUTNL; gNL*(1-exp(-((putativeoogcontrast/aNL).^bNL)))];
                    
                end
            end
        end
    end
    % Accumulating files over different sessions
    filedates = cell2mat(newfilename);
    datesind = 2:10;
    filedates = str2num(filedates(:,datesind));
    [filedates,~,ib] = unique([filedates RF newRWmultiplier newlaserdial newstimsize],'rows');
    HitLSession = []; HitNLSession = [];
    MissLSession = []; MissNLSession = [];
    CRLSession = []; CRNLSession = [];
    FALSession = []; FANLSession = [];
    stimabsentLSession = [];
    stimabsentNLSession = [];
    stimpresentLSession = [];
    stimpresentNLSession = [];
    p_hits = [];
    p_FA = [];
    p_hitsLFANL = [];
    numfilespersession = [];
    nrows = ceil(sqrt(numel(unique(ib))));
    for ii = 1:numel(unique(ib))
        idx = ib==ii;
        HitLSession = [HitLSession; sum(HitL(idx))];
        HitNLSession = [HitNLSession ;sum(HitNL(idx))];
        MissLSession = [MissLSession; sum(MissL(idx))];
        MissNLSession = [MissNLSession; sum(MissNL(idx))];
        CRLSession = [CRLSession; sum(CRL(idx))];
        CRNLSession = [CRNLSession; sum(CRNL(idx))];
        FALSession = [FALSession; sum(FAL(idx)) ];
        FANLSession = [FANLSession; sum(FANL(idx))];
        stimabsentLSession = [stimabsentLSession; sum(stimabsentL(idx))];
        stimabsentNLSession = [stimabsentNLSession; sum(stimabsentNL(idx))];
        stimpresentLSession = [stimpresentLSession; sum(stimpresentL(idx))];
        stimpresentNLSession = [stimpresentNLSession; sum(stimpresentNL(idx))];
        [~,p1] = equalproptest([sum(HitL(idx)) sum(HitNL(idx))],[sum(stimpresentL(idx)) sum(stimpresentNL(idx))],alpha);
        [~,p2] = equalproptest([sum(FAL(idx)) sum(FANL(idx))],[sum(stimabsentL(idx)) sum(stimabsentNL(idx))],alpha);
        [~,p3] = equalproptest([sum(HitL(idx)) sum(FANL(idx))],[sum(stimpresentL(idx)) sum(stimabsentNL(idx))],alpha);
        p_hits = [p_hits;p1]; p_FA = [p_FA;p2]; p_hitsLFANL = [p_hitsLFANL; p3];
        numfilespersession = [numfilespersession; sum(idx)];
    end
    
    idx = p_hits<1.0;
    HitLSession = HitLSession(idx);
    HitNLSession = HitNLSession(idx);
    MissLSession = MissLSession(idx);
    MissNLSession = MissNLSession(idx);
    CRLSession = CRLSession(idx);
    CRNLSession = CRNLSession(idx);
    FALSession = FALSession(idx);
    FANLSession = FANLSession(idx);
    stimabsentLSession = stimabsentLSession(idx);
    stimabsentNLSession = stimabsentNLSession(idx);
    stimpresentLSession = stimpresentLSession(idx);
    stimpresentNLSession = stimpresentNLSession(idx);
    TPRLSession = HitLSession./(HitLSession+MissLSession);
    TPRNLSession = HitNLSession./(HitNLSession+MissNLSession);
    FPRLSession = FALSession./(FALSession+CRLSession);
    FPRNLSession = FANLSession./(FANLSession+CRNLSession);
    p_hits = p_hits(idx); p_FA = p_FA(idx); p_hitsLFANL = p_hitsLFANL(idx);
    
    figure(plot_counter); set(gcf,'Name',strcat('stats:',monkeyname));
    subplot(2,2,2*mm-1);
    for ii = 1:numel(p_hits)
        if (p_hits(ii)<0.05)
            color = [0 0 0];
        else
            color = [0.5 0.5 0.5];
        end
        px = HitNLSession(ii)./stimpresentNLSession(ii); errorx = sqrt(px*(1-px)/stimpresentNLSession(ii));
        py = HitLSession(ii)./stimpresentLSession(ii); errory = sqrt(py*(1-py)/stimpresentLSession(ii));
        plot(px,py,'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',color,'MarkerEdgeColor',[1 1 1],'Color','k'); hold on;
    end
    line([0 1],[0 1],'Linewidth',1);axis square; xlabel('Hits control'); ylabel('Hits laser'); set(gca,'Xlim',[0 1],'Ylim',[0 1],'Tickdir','out','XTick',[0:0.5:1],'YTick',[0:0.5:1]);
    if mm == 1
        title('Maui');
    else
        title('Apollo');
    end
    [~,p] = corr((HitNLSession./stimpresentNLSession)-(HitLSession./stimpresentLSession),(CRNLSession./stimabsentNLSession)-(CRLSession./stimabsentLSession));
    
    %     bins = logspace(-0.5,log10(30),21);
    %     subplot(2,4,4*mm-2); histogram(prop_correctatGAMUTNL./prop_correctatGAMUTL,bins,'Normalization','probability','FaceColor',[0 0 0],'EdgeColor',[1 1 1]); hold on;
    %     plot(median(prop_correctatGAMUTNL./prop_correctatGAMUTL),0,'kv','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
    %     set(gca,'Ylim',[0 0.15],'Xlim',[0.8 30],'Tickdir','out','XTick',[1 3 10 30],'XScale','log');axis square; xlabel('ratio(prop of correct):control vs laser'); ylabel('probability'); hold off;
    
    ratio_controloverlaser = [ratio_controloverlaser; median(prop_correctatGAMUTNL./prop_correctatGAMUTL)];
    
    %     subplot(2,4,4*mm-1); plot(prop_correctatGAMUTNL,prop_correctatGAMUTL,'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on
    %     set(gca,'Ylim',[0 1],'Xlim',[0 1],'Tickdir','out','XTick',[0:0.5:1],'YTick',[0:0.5:1]); line([0 1],[0 1],'Linewidth',1);axis square; xlabel('prop of correct:control'); ylabel('prop of correct:laser'); hold off;
    % if HitL(end)==0
    %     dprimeL(mm,aa) = norminv(1-(0.5./stimpresentL(end)))-norminv(max([0.5./stimabsentL(end) FAL(end)./stimabsentL(end)]));
    % else
    %     dprimeL(mm,aa) = norminv(min([1-(0.5./stimpresentL(end)) HitL(end)./stimpresentL(end)]))-norminv(max([0.5./stimabsentL(end) FAL(end)./stimabsentL(end)]));
    % end
    dprimeL = norminv(min([1-(0.5./stimpresentLSession) HitLSession./stimpresentLSession],[],2))-norminv(max([0.5./stimabsentLSession FALSession./stimabsentLSession],[],2));
    dprimeNL = norminv(min([1-(0.5./stimpresentNLSession) HitNLSession./stimpresentNLSession],[],2))-norminv(max([0.5./stimabsentNLSession FANLSession./stimabsentNLSession],[],2));
    if any(HitLSession==0)
        ind = HitLSession==0;
        dprimeL(ind) = norminv(0.5./stimabsentLSession(ind))-norminv(max([0.5./stimabsentLSession(ind) FALSession(ind)./stimabsentLSession(ind)],[],2));
    elseif any(HitNLSession==0)
        ind = HitNLSession==0;
        dprimeNL(ind) = norminv(0.5./stimabsentNLSession(ind))-norminv(max([0.5./stimabsentNLSession(ind) FANLSession(ind)./stimabsentNLSession(ind)],[],2));
    elseif any(HitLSession==30)
        ind = HitLSession==30;
        dprimeL(ind) = norminv(1-(0.5./stimpresentLSession(ind)))-norminv(max([0.5./stimabsentLSession(ind) FALSession(ind)./stimabsentLSession(ind)],[],2));
    elseif any(HitNLSession==30)
        ind = HitNLSession==30;
        dprimeNL(ind) = norminv(1-(0.5./stimpresentNLSession(ind)))-norminv(max([0.5./stimabsentNLSession(ind) FANLSession(ind)./stimabsentNLSession(ind)],[],2));
    elseif any(FALSession==0)
        ind = FALSession==0;
        dprimeL(ind) = norminv(min([1-(0.5./stimpresentLSession(ind)) HitLSession(ind)./stimpresentLSession(ind)],[],2))-norminv(0.5./stimabsentLSession(ind));
    elseif any(FANLSession==0)
        ind = FANLSession==0;
        dprimeNL(ind) = norminv(min([1-(0.5./stimpresentNLSession(ind)) HitNLSession(ind)./stimpresentNLSession(ind)],[],2))-norminv(0.5./stimabsentNLSession(ind));
    end
    figure(plot_counter); subplot(2,2,2*mm); plot(dprimeNL,dprimeL,'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on
    set(gca,'Ylim',[-1 4],'Xlim',[-1 4],'Tickdir','out','XTick',[-1:1:4],'YTick',[-1:1:4]); line([-1 4],[-1 4],'Linewidth',1);axis square; xlabel('d prime:control'); ylabel('d prime:laser'); hold off;
    
    % plotting reliability across sessions 
    k = cell2mat(filename);
    [~,ia,ic] = unique(k(:,1:5),'stable','rows');
    mu = []; sigma = [];
    for aa =1:numel(ia)
        figure(plot_counter+1);subplot(2,2,2*mm-1); plot(aa*ones(sum(ic==aa),1),dprimeNL(ic==aa)-dprimeL(ic==aa),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
        mu = [mu; mean(dprimeNL(ic==aa)-dprimeL(ic==aa))];
        sigma = [sigma; std(dprimeNL(ic==aa)-dprimeL(ic==aa))];
    end
    set(gca,'Tickdir','out','Ylim',[-1 4],'YTick',-1:1:4,'Xlim',[1 numel(ia)],'XTick',1:1:numel(ia)); axis square; xlabel('Session'); ylabel('dprimeNL-dprimeL');
    figure(plot_counter+1);subplot(2,2,2*mm); errorbar(1:numel(ia),mu,sigma,'-o','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
    set(gca,'Tickdir','out','Ylim',[-1 4],'YTick',-1:1:4,'Xlim',[1 numel(ia)],'XTick',1:1:numel(ia)); axis square; xlabel('Session'); ylabel('dprimeNL-dprimeL');
end
figure(plot_counter); set(gcf,'renderer','painters');
figure(plot_counter+1); set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;


average_time_per_trials = mean(timeperblock./trialsperblock);

%% Figure 7
% �	DToneloc: effect of laser intensity on hits in control and laser trials (3 figures: 1 control staircase, 1 laser staircase, 1 diff in hits in laser and control trials
% vs laserdial/ laser intensity)(files collected from Apollo on 02/15/19)
if ~exist('plot_counter')
    plot_counter = 1;
end
load('T_vos1978_Y');
Vlambda = T_vos1978_Y';
filename = {'A021519004.nex';'A021519005.nex';'A021519006.nex';'A021519007.nex';'A021519008.nex';'A021519009.nex';'A021519010.nex'};
RWmultiplier = [0.8*ones(numel(filename),1)];
load laserdetails.mat
laserdial = [1.00;0.50;1.50;0.75;2.00;0.25;2.50];
laserdial = spline(laserdetails.dial,laserdetails.laserpower,laserdial);
HitL = []; HitNL = [];
MissL = []; MissNL = [];
CRL = []; CRNL = [];
FAL = []; FANL = [];
TPRNL = []; TPRL = [];
FPRNL = []; FPRL = [];
RF = [];
stimabsentL = []; stimabsentNL = [];
stimpresentL = []; stimpresentNL = [];
HitsaccadeRT = []; MisssaccadeRT = [];
CRsaccadeRT = []; FAsaccadeRT = [];
newfilename = [];
newRWmultiplier = [];
newstimsize = [];
newlaserdial = [];
thresholdsNL = [];
outofgamutptsNL = [];
averaged_staircasecontrastNL = [];
averaged_staircasecontrastL = [];
dprime_diff = [];
dprime_control = [];
staircaseterminationptNL = [];
cgradations = linspace(0.1,1,numel(filename));
[~,~,idx] = unique(laserdial);
cgradations = cgradations(idx);
figure(plot_counter), set(gcf,'Name','Effect of Laser');
for jj = 1:size(filename,1)
    stro = nex2stro(findfile(char(filename(jj,:))));
    fponidx  = strcmp(stro.sum.trialFields(1,:),'fpon_t');
    fpacqidx  = strcmp(stro.sum.trialFields(1,:),'fpacq_t');
    fpoffidx  = strcmp(stro.sum.trialFields(1,:),'fpoff_t');
    stimonidx  = strcmp(stro.sum.trialFields(1,:),'stimon_t');
    stimoffidx  = strcmp(stro.sum.trialFields(1,:),'stimoff_t');
    targonidx  = strcmp(stro.sum.trialFields(1,:),'targon_t');
    saccstartidx  = strcmp(stro.sum.trialFields(1,:),'saccstart_t');
    saccendidx  = strcmp(stro.sum.trialFields(1,:),'saccend_t');
    rewidx  = strcmp(stro.sum.trialFields(1,:),'rew_t');
    stimpresentidx  = strcmp(stro.sum.trialFields(1,:),'stimpresent');
    Lcc = strcmp(stro.sum.trialFields(1,:),'lcc');
    Mcc = strcmp(stro.sum.trialFields(1,:),'mcc');
    Scc = strcmp(stro.sum.trialFields(1,:),'scc');
    tf = strcmp(stro.sum.trialFields(1,:),'tf');
    oog = strcmp(stro.sum.trialFields(1,:),'oog');
    stimidx = strcmp(stro.sum.trialFields(1,:),'stim_idx');
    optstim = strcmp(stro.sum.trialFields(1,:),'optstim');
    correct = strcmp(stro.sum.trialFields(1,:),'correct');
    laseron = strcmp(stro.sum.trialFields(1,:),'laseron_t');
    laseroff = strcmp(stro.sum.trialFields(1,:),'laseroff_t');
    correcttrials = stro.trial(:,correct);
    diridxs = stro.trial(:,stimidx);
    stimpresent = logical(stro.trial(:,stimpresentidx));
    LMStriplet = [stro.trial(:,Lcc) stro.trial(:,Mcc) stro.trial(:,Scc)];
    lasertrials = logical(stro.trial(:,optstim));
    if isnan(unique(stro.trial(:,21:22)))
        stimlocations = [stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y];
        RF = [RF; stimlocations];
        newfilename = [newfilename; filename(jj)];
        newRWmultiplier = [newRWmultiplier; RWmultiplier(jj)];
        newlaserdial = [newlaserdial; laserdial(jj)];
        ib = ones(size(lasertrials));
        newstimsize = [newstimsize; stro.sum.exptParams.sigma];
    else
        [stimlocations,~,ib] = unique(stro.trial(:,21:22),'rows');
        RF = [RF; stimlocations];
        L = size(stimlocations,1);
        newfilename = [newfilename; repmat(filename(jj),[L 1])];
        newRWmultiplier = [newRWmultiplier; repmat(RWmultiplier(jj),[L 1])];
        newlaserdial = [newlaserdial; repmat(laserdial(jj),[L 1])];
        newstimsize = [newstimsize; repmat(stro.sum.exptParams.sigma,[L 1])];
    end
    
    for ii = 1:numel(unique(ib))
        ind = logical(ib==ii);
        lasertrialidxs = logical(stro.trial(ind,optstim));
        colordirchoiceidxs = correcttrials(ind);
        colorstimpresent = logical(stimpresent(ind));
        percentcorrectlasertrials = sum(colordirchoiceidxs & lasertrialidxs)/sum(lasertrialidxs);
        percentcorrectnonlasertrials = sum(colordirchoiceidxs & ~lasertrialidxs)/sum(~lasertrialidxs);
        
        % Laser trials
        Hitlasertrial = sum(colordirchoiceidxs & lasertrialidxs & colorstimpresent); % Hit
        Misslasertrial = sum(~colordirchoiceidxs & lasertrialidxs & colorstimpresent); % Miss
        CRlasertrial = sum(colordirchoiceidxs & lasertrialidxs & ~colorstimpresent); % Correct Reject
        FAlasertrial = sum(~colordirchoiceidxs & lasertrialidxs & ~colorstimpresent); % False Alarm
        CRlasertrialidxs = find(lasertrialidxs & ~colorstimpresent); L1 = numel(CRlasertrialidxs);
        
        % Non-Laser trials
        Hitnonlasertrial = sum(colordirchoiceidxs & ~lasertrialidxs & colorstimpresent); % Hit
        Missnonlasertrial = sum(~colordirchoiceidxs & ~lasertrialidxs & colorstimpresent); % Miss
        CRnonlasertrial = sum(colordirchoiceidxs & ~lasertrialidxs & ~colorstimpresent); % Correct Reject
        FAnonlasertrial = sum(~colordirchoiceidxs & ~lasertrialidxs & ~colorstimpresent); % False Alarm
        CRnonlasertrialidxs = find(~lasertrialidxs & ~colorstimpresent); L2 = numel(CRnonlasertrialidxs);
        
        
        stimabsentL = [stimabsentL; sum(CRlasertrial) + sum(FAlasertrial)];
        stimabsentNL = [stimabsentNL; sum(CRnonlasertrial) + sum(FAnonlasertrial)];
        stimpresentL = [stimpresentL; sum(Hitlasertrial) + sum(Misslasertrial)];
        stimpresentNL = [stimpresentNL; sum(Hitnonlasertrial) + sum(Missnonlasertrial)];
        
        % Now storing the Hits, Miss, CR and FA for laser and non-laser trials
        HitL = [HitL; Hitlasertrial]; HitNL = [HitNL; Hitnonlasertrial];
        MissL = [MissL; Misslasertrial]; MissNL = [MissNL; Missnonlasertrial];
        CRL = [CRL; CRlasertrial]; CRNL = [CRNL; CRnonlasertrial];
        FAL = [FAL; FAlasertrial]; FANL = [FANL; FAnonlasertrial];
        TPRNL = [TPRNL; Hitnonlasertrial/(Hitnonlasertrial+Missnonlasertrial)]; % True positive ratio, non-laser trial
        FPRNL = [FPRNL; FAnonlasertrial/(FAnonlasertrial+CRnonlasertrial)]; % False positive ratio, non-laser trial
        TPRL = [TPRL; Hitlasertrial/(Hitlasertrial+Misslasertrial)]; % True positive ratio, laser trial
        FPRL = [FPRL; FAlasertrial/(FAlasertrial+CRlasertrial)]; % False positive ratio, laser trial
        
        % Calculating psychophysical detection threshold for non-laser trials
        LMStripletfile = LMStriplet(ind,:);
        LMScontrastfileNL = sqrt(sum(LMStripletfile(colorstimpresent & ~lasertrialidxs,:).^2,2));
        LMScontrastfileL = sqrt(sum(LMStripletfile(colorstimpresent & lasertrialidxs,:).^2,2));
        staircaseterminationptNL = [staircaseterminationptNL; LMScontrastfileNL(end)];
        
        fundamentals = stro.sum.exptParams.fundamentals; % CONE FUNDAMENTALS: L,M,S
        fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
        mon_spd = stro.sum.exptParams.mon_spd; % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
        mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
        mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
        M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals
        Lb = M*stro.sum.exptParams.bkgndrgb;
        % Converting all the contrasts to Luminance contrast
        Bkgnd_luminance = stro.sum.exptParams.bkgndrgb'*mon_spd*Vlambda;
        Peak_luminance = (inv(M)*((1+LMStriplet).*(repmat(Lb',[size(LMStriplet,1) 1])))')'*mon_spd*Vlambda;
        Trough_luminance = (inv(M)*((1-LMStriplet).*(repmat(Lb',[size(LMStriplet,1) 1])))')'*mon_spd*Vlambda;
        Contrast_luminance = (Peak_luminance - Trough_luminance)./(Peak_luminance + Trough_luminance);
        LMScontrastfileNL = Contrast_luminance(colorstimpresent & ~lasertrialidxs,:);
        LMScontrastfileL = Contrast_luminance(colorstimpresent & lasertrialidxs,:);
        % calculating gamut edge luminance contrast
        t = min(((1./stro.sum.exptParams.bkgndrgb)-1));
        OOG_luminance = (stro.sum.exptParams.bkgndrgb*(1+t))' * mon_spd*Vlambda;
        putativeoogcontrast =  (OOG_luminance - Bkgnd_luminance)/Bkgnd_luminance;
        
        % A simple criteria for detecting "bad files"
        tmp = stro.trial(ind & stimpresent & ~lasertrials,oog);
        outofgamutptsNL = [outofgamutptsNL; sum(tmp(end-4:end))];
        
        answers = colordirchoiceidxs(colorstimpresent & ~lasertrialidxs);
        LMScontrastfileNL(LMScontrastfileNL>putativeoogcontrast) = putativeoogcontrast;
        LMScontrastfileL(LMScontrastfileL>putativeoogcontrast) = putativeoogcontrast;
        contrast = unique(LMScontrastfileNL);
        correctanswers = zeros(size(contrast));
        wronganswers = zeros(size(contrast));
        trialspercontrast = zeros(size(contrast));
        
        for mm = 1:numel(contrast)
            trialspercontrast(mm) = numel(answers(LMScontrastfileNL==contrast(mm)));
            correctanswers(mm) = sum(answers(LMScontrastfileNL==contrast(mm)));
            wronganswers(mm) = trialspercontrast(mm) - correctanswers(mm);
        end
        [a,~,~] = weibullFitforDToneloc(contrast,[correctanswers wronganswers],'mle');
        thresholdsNL = [thresholdsNL; a];
        
    end
    
    averaged_staircasecontrastL = [averaged_staircasecontrastL; mean(LMScontrastfileL(1:numel(LMScontrastfileL)/2+1))];
    averaged_staircasecontrastNL = [averaged_staircasecontrastNL; mean(LMScontrastfileNL(numel(LMScontrastfileNL)/2+1:end))];
    dprimeL = norminv(min([1-(0.5./stimpresentL(jj)) HitL(jj)./stimpresentL(jj)],[],2)) - norminv(max([0.5./stimabsentL(jj) FAL(jj)./stimabsentL(jj)],[],2));
    dprimeNL = norminv(min([1-(0.5./stimpresentNL(jj)) HitNL(jj)./stimpresentNL(jj)],[],2)) - norminv(max([0.5./stimabsentNL(jj) FANL(jj)./stimabsentNL(jj)],[],2));
    dprime_diff = [dprime_diff; dprimeNL-dprimeL];
    dprime_control = [dprime_control; dprimeNL];
    figure(plot_counter);
    subplot(321); plot(LMScontrastfileL,'-','color',[1-cgradations(jj) 1-cgradations(jj) 1-cgradations(jj)],'Linewidth',2); hold on;
    subplot(322); plot(LMScontrastfileNL,'-','color',[1-cgradations(jj) 1-cgradations(jj) 1-cgradations(jj)],'Linewidth',2); hold on;
    subplot(323); plot(laserdial(jj),HitNL(jj)./stimpresentNL(jj) - HitL(jj)./stimpresentL(jj),'o','MarkerSize',6,'MarkerFaceColor',[1-cgradations(jj) 1-cgradations(jj) 1-cgradations(jj)],'MarkerEdgeColor',[1 1 1]); hold on;
    subplot(324); plot(laserdial(jj),dprime_diff(end),'o','MarkerSize',6,'MarkerFaceColor',[1-cgradations(jj) 1-cgradations(jj) 1-cgradations(jj)],'MarkerEdgeColor',[1 1 1]); hold on;
    
    
end
subplot(321); xlabel('trials'); ylabel('Luminance contrast'); title('laser staircase'); set(gca,'Ylim',[0.15 0.75],'TickDir','out','Xlim',[0 30],'YTick',[0.15:0.15:0.75],'XTick',[0:15:30]); axis square; hold off;
subplot(322); xlabel('trials'); ylabel('Luminance contrast'); title('control staircase'); set(gca,'Ylim',[0.15 0.75],'TickDir','out','Xlim',[0 30],'YTick',[0.15:0.15:0.75],'XTick',[0:15:30]); axis square; hold off;

% plotting the differences in hits between laser and control trials
[model] = sigmoidalfit_AD(laserdial,HitNL./stimpresentNL - HitL./stimpresentL);
[model2] = sigmoidalfit_AD(laserdial,dprime_diff);
x = linspace(0,max(laserdial),51);
fit = model(1)*(x.^model(3))./(x.^model(3)+model(2).^model(3));
fit2 = model2(1)*(x.^model2(3))./(x.^model2(3)+model2(2).^model2(3));
subplot(323);plot(linspace(0,max(laserdial),51),fit,'color',[0 0 0],'Linewidth',2); xlabel('Laser power (mW)'); ylabel('Hits control - Hits laser'); set(gca,'Ylim',[0 1],'Xlim',[0 100],'XTick',[0:20:100],'TickDir','out'); axis square; hold off;
subplot(324);plot(linspace(0,max(laserdial),51),fit2,'color',[0 0 0],'Linewidth',2); xlabel('Laser power (mW)'); ylabel('dprime control-laser'); set(gca,'Ylim',[0 3],'Xlim',[0 100],'XTick',[0:20:100],'TickDir','out'); axis square; hold off;
[~,p] = corr(laserdial,linspace(0,numel(laserdial)-1,numel(laserdial))');

% Verifying whether there is any correlation between laser power and averaged control staicase contrast
[r1,p1] = corr(laserdial,averaged_staircasecontrastNL,'type','Spearman');
[r2,p2] = corr(laserdial,averaged_staircasecontrastL,'type','Spearman');

[r3,p3] = corr(laserdial,averaged_staircasecontrastNL,'type','Pearson');
[r4,p4] = corr(laserdial,averaged_staircasecontrastL,'type','Pearson');
[r5,p5] = corr(laserdial,dprime_control,'type','Pearson');


% �	Effect of hits in control and  laser trials as function of sessions (Fetch analysis)
load('T_vos1978_Y');
Vlambda = T_vos1978_Y';
filenameopto1 = {'M042518002.nex';'M042518004.nex';'M042518005.nex';'M042518008.nex';'M042518006.nex';'M042518007.nex'}; % 4/25
RWmultiplieropto1 =[1.0;1.0;1.0;1.0;1.0;1.0];
laserdialopto1 = [2.5;2.5;2.5;2.5;2.5;2.5];

filenameopto2 = {'M042618003.nex';'M042618004.nex';'M042618005.nex';'M042618006.nex';'M042618007.nex'}; % 4/26, 200 stim, 300 laser
RWmultiplieropto2 =[1.0;1.0;1.0;1.0;1.0];
laserdialopto2 = [2.5;2.5;2.5;2.5;2.5];

filenameopto3 = {'M051418008.nex';'M051418009.nex';'M051418010.nex';'M051418011.nex';'M051418013.nex';'M051418014.nex';'M051418015.nex';'M051418016.nex'}; %5/14, 200 stim, 300 laser
RWmultiplieropto3 =[1.0;1.0;0.7;0.7;1.0;1.0;0.7;0.7];
laserdialopto3 = 4.0*ones(size(RWmultiplieropto3));

filenameopto4 = {'M051518003.nex';'M051518004.nex';'M051518005.nex';'M051518006.nex';'M051518007.nex';'M051518008.nex';'M051518009.nex';'M051518010.nex';'M051518011.nex'}; %5/15, 200 stim, 300 laser
RWmultiplieropto4 =[1.0;0.7;0.5;0.5;0.3;0.6;1.0;0.4;0.7];
laserdialopto4 = 4.0*ones(size(RWmultiplieropto4));

filenameopto5 = {'A020119001.nex';'A020119002.nex';'A020119003.nex';'A020119004.nex';'A020119005.nex';'A020119006.nex';'A020119009.nex';'A020119010.nex';'A020119011.nex';'A020119012.nex';'A020119013.nex'};
RWmultiplieropto5 = [0.8;0.8;0.8;0.8;1.0;0.7;0.8;0.8;0.7;0.8;0.7];
laserdialopto5 = [2.0;2.0;2.0;2.0;2.0;2.0;2.0;2.0;2.0;2.0;0.2];

filenameopto6 = {'A012319007.nex';'A012319008.nex';'A012319009.nex';'A012319010.nex';'A012319011.nex'};
RWmultiplieropto6 = [0.7;0.7;0.7;0.7;0.7];
laserdialopto6 = [1.5;1.5;1.5;1.5;1.5];

filenameopto7 = {'A013019006.nex';'A013019007.nex';'A013019008.nex';'A013019009.nex';'A013019012.nex'};
RWmultiplieropto7 = [0.7;0.7;0.5;0.5;0.7];
laserdialopto7 = [2.2;2.2;2.2;2.2;2.2];

filenameopto8 = {'A020519001.nex';'A020519002.nex';'A020519003.nex';'A020519004.nex';'A020519005.nex';'A020519009.nex'};
RWmultiplieropto8 = [0.8;0.7;0.7;0.8;1.0;0.8];
laserdialopto8 = [1.5;1.5;1.5;1.5;1.5;1.5];

filenameopto9 = {'A013119001.nex';'A013119002.nex';'A013119005.nex';'A013119006.nex';'A013119007.nex';'A013119008.nex'};
RWmultiplieropto9 = [0.7;0.7;0.6;0.6;0.7;1.0];
laserdialopto9 = [1.0;1.0;1.0;1.0;1.0;1.0];

oogvals = [];
nrows = 9;
figure(plot_counter);
RFs = cell(nrows,1);
diffinpropshits_early = cell(nrows,1); % 1-480 trials
diffinpropshits_late = cell(nrows,1); % 480+ trials
diffinpropshits = zeros(nrows,11);
dprime_early = cell(nrows,1); % 1-480 trials
dprime_late = cell(nrows,1); % 480+ trials
dprime = zeros(nrows,11);
timeperblock = [];
trialsperblock = [];
for mm = 1:nrows
    filename = eval(strcat('filenameopto',num2str(mm)));
    HitL = []; HitNL = [];
    MissL = []; MissNL = [];
    CRL = []; CRNL = [];
    FAL = []; FANL = [];
    TPRNL = []; TPRL = [];
    FPRNL = []; FPRL = [];
    RF = [];
    stimabsentL = []; stimabsentNL = [];
    stimpresentL = []; stimpresentNL = [];
    weibullparamsL = [];
    weibullparamsNL = [];
    stimsizeinsigmas = [];
    cgradations = linspace(0.2,1,size(filename,1));
    contrastlattice = linspace(0,1.3,21);
    for aa = 1:size(filename,1)
        stro = nex2stro(findfile(char(filename(aa,:))));
        fponidx  = strcmp(stro.sum.trialFields(1,:),'fpon_t');
        fpacqidx  = strcmp(stro.sum.trialFields(1,:),'fpacq_t');
        fpoffidx  = strcmp(stro.sum.trialFields(1,:),'fpoff_t');
        stimonidx  = strcmp(stro.sum.trialFields(1,:),'stimon_t');
        stimoffidx  = strcmp(stro.sum.trialFields(1,:),'stimoff_t');
        targonidx  = strcmp(stro.sum.trialFields(1,:),'targon_t');
        saccstartidx  = strcmp(stro.sum.trialFields(1,:),'saccstart_t');
        saccendidx  = strcmp(stro.sum.trialFields(1,:),'saccend_t');
        rewidx  = strcmp(stro.sum.trialFields(1,:),'rew_t');
        stimpresentidx  = strcmp(stro.sum.trialFields(1,:),'stimpresent');
        Lcc = strcmp(stro.sum.trialFields(1,:),'lcc');
        Mcc = strcmp(stro.sum.trialFields(1,:),'mcc');
        Scc = strcmp(stro.sum.trialFields(1,:),'scc');
        tf = strcmp(stro.sum.trialFields(1,:),'tf');
        oog = strcmp(stro.sum.trialFields(1,:),'oog');
        stimidx = strcmp(stro.sum.trialFields(1,:),'stim_idx');
        optstim = strcmp(stro.sum.trialFields(1,:),'optstim');
        correct = strcmp(stro.sum.trialFields(1,:),'correct');
        laseron = strcmp(stro.sum.trialFields(1,:),'laseron_t');
        laseroff = strcmp(stro.sum.trialFields(1,:),'laseroff_t');
        correcttrials = stro.trial(:,correct);
        diridxs = stro.trial(:,stimidx);
        stimpresent = logical(stro.trial(:,stimpresentidx));
        LMStriplet = [stro.trial(:,Lcc) stro.trial(:,Mcc) stro.trial(:,Scc)];
        lasertrials = logical(stro.trial(:,optstim));
        anlgStartTime = strcmp(stro.sum.rasterCells,'anlgStartTime');
        timeperblock = [timeperblock; stro.ras{end,anlgStartTime}-stro.ras{1,anlgStartTime}];
        trialsperblock = [trialsperblock; size(stro.ras,1)];
        
        RF = [RF; stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y];
        stimsizeinsigmas = [stimsizeinsigmas; stro.sum.exptParams.sigma];
        
        lasertrialidxs = logical(stro.trial(:,optstim));
        colordirchoiceidxs = correcttrials;
        colorstimpresent = logical(stimpresent);
        percentcorrectlasertrials = sum(colordirchoiceidxs & lasertrialidxs)/sum(lasertrialidxs);
        percentcorrectnonlasertrials = sum(colordirchoiceidxs & ~lasertrialidxs)/sum(~lasertrialidxs);
        
        % Laser trials
        Hitlasertrial = sum(colordirchoiceidxs & lasertrialidxs & colorstimpresent); % Hit
        Misslasertrial = sum(~colordirchoiceidxs & lasertrialidxs & colorstimpresent); % Miss
        CRlasertrial = sum(colordirchoiceidxs & lasertrialidxs & ~colorstimpresent); % Correct Reject
        FAlasertrial = sum(~colordirchoiceidxs & lasertrialidxs & ~colorstimpresent); % False Alarm
        
        % Non-Laser trials
        Hitnonlasertrial = sum(colordirchoiceidxs & ~lasertrialidxs & colorstimpresent); % Hit
        Missnonlasertrial = sum(~colordirchoiceidxs & ~lasertrialidxs & colorstimpresent); % Miss
        CRnonlasertrial = sum(colordirchoiceidxs & ~lasertrialidxs & ~colorstimpresent); % Correct Reject
        FAnonlasertrial = sum(~colordirchoiceidxs & ~lasertrialidxs & ~colorstimpresent); % False Alarm
        
        stimabsentL = [stimabsentL; sum(CRlasertrial) + sum(FAlasertrial)];
        stimabsentNL = [stimabsentNL; sum(CRnonlasertrial) + sum(FAnonlasertrial)];
        stimpresentL = [stimpresentL; sum(Hitlasertrial) + sum(Misslasertrial)];
        stimpresentNL = [stimpresentNL; sum(Hitnonlasertrial) + sum(Missnonlasertrial)];
        
        % Now storing the Hits, Miss, CR and FA for laser and non-laser trials
        HitL = [HitL; Hitlasertrial]; HitNL = [HitNL; Hitnonlasertrial];
        MissL = [MissL; Misslasertrial]; MissNL = [MissNL; Missnonlasertrial];
        CRL = [CRL; CRlasertrial]; CRNL = [CRNL; CRnonlasertrial];
        FAL = [FAL; FAlasertrial]; FANL = [FANL; FAnonlasertrial];
        TPRNL = [TPRNL; Hitnonlasertrial/(Hitnonlasertrial+Missnonlasertrial)]; % True positive ratio, non-laser trial
        FPRNL = [FPRNL; FAnonlasertrial/(FAnonlasertrial+CRnonlasertrial)]; % False positive ratio, non-laser trial
        TPRL = [TPRL; Hitlasertrial/(Hitlasertrial+Misslasertrial)]; % True positive ratio, laser trial
        FPRL = [FPRL; FAlasertrial/(FAlasertrial+CRlasertrial)]; % False positive ratio, laser trial
        
        fundamentals = stro.sum.exptParams.fundamentals; % CONE FUNDAMENTALS: L,M,S
        fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
        mon_spd = stro.sum.exptParams.mon_spd; % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
        mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
        mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
        M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals
        Lb = M*stro.sum.exptParams.bkgndrgb;
        % Converting all the contrasts to Luminance contrast
        Bkgnd_luminance = stro.sum.exptParams.bkgndrgb'*mon_spd*Vlambda;
        Contrast_luminance = ((inv(M)*((LMStriplet+1).*(repmat(Lb',[size(LMStriplet,1) 1])))')'*mon_spd*Vlambda-Bkgnd_luminance)/Bkgnd_luminance;
        LMScontrastfileNL = Contrast_luminance(colorstimpresent & ~lasertrialidxs,:);
        LMScontrastfileL = Contrast_luminance(colorstimpresent & lasertrialidxs,:);
        % calculating gamut edge luminance contrast
        t = min(((1./stro.sum.exptParams.bkgndrgb)-1));
        OOG_luminance = (stro.sum.exptParams.bkgndrgb*(1+t))' * mon_spd*Vlambda;
        putativeoogcontrast =  (OOG_luminance - Bkgnd_luminance)/Bkgnd_luminance;
        
        trials = lasertrials;
        % only evaluating the laser trials
        for kk = 1:2
            if kk == 2
                trials = ~trials;
            end
            LMScontrast = [];
            LMScontrast = Contrast_luminance(stimpresent & trials,:);
            answers = correcttrials(stimpresent & trials);
            LMScontrast(LMScontrast>putativeoogcontrast) = putativeoogcontrast;
            oogvals = [oogvals; putativeoogcontrast];
            
            
            if kk==1
                contrastL = unique(LMScontrast);
                correctanswersL = zeros(size(contrastL));
                wronganswersL = zeros(size(contrastL));
                trialspercontrastL = zeros(size(contrastL));
                for jj = 1:numel(contrastL)
                    trialspercontrastL(jj) = numel(answers(LMScontrast==contrastL(jj)));
                    correctanswersL(jj) = sum(answers(LMScontrast==contrastL(jj)));
                    wronganswersL(jj) = trialspercontrastL(jj) - correctanswersL(jj);
                    percorrectL(jj) = correctanswersL(jj)/trialspercontrastL(jj);
                end
                [aL,bL,gL] = weibullFitforDToneloc(contrastL,[correctanswersL wronganswersL],'mle');
                weibullparamsL = [weibullparamsL; aL bL gL]; % laser trials
                LMScontrastL = LMScontrast;
            else
                contrastNL = unique(LMScontrast);
                correctanswersNL = zeros(size(contrastNL));
                wronganswersNL = zeros(size(contrastNL));
                trialspercontrastNL = zeros(size(contrastNL));
                for jj = 1:numel(contrastNL)
                    trialspercontrastNL(jj) = numel(answers(LMScontrast==contrastNL(jj)));
                    correctanswersNL(jj) = sum(answers(LMScontrast==contrastNL(jj)));
                    wronganswersNL(jj) = trialspercontrastNL(jj) - correctanswersNL(jj);
                    percorrectNL(jj) = correctanswersNL(jj)/trialspercontrastNL(jj);
                end
                [aNL,bNL,gNL] = weibullFitforDToneloc(contrastNL,[correctanswersNL wronganswersNL],'mle');
                weibullparamsNL = [weibullparamsNL; aNL bNL gNL]; % laser trials
                LMScontrastNL = LMScontrast;
            end
        end
        fitL = gL*(1-exp(-((contrastlattice./aL).^bL)));
        fitNL = gNL*(1-exp(-((contrastlattice./aNL).^bNL)));
        diffinpropshits(mm,aa) = HitNL(end)./stimpresentNL(end)-HitL(end)./stimpresentL(end);
        
        if HitNL(end)==0
            dprimeNL(mm,aa) = norminv(1-(0.5./stimpresentNL(end)))-norminv(max([0.5./stimabsentNL(end) FANL(end)./stimabsentNL(end)]));
        else
            dprimeNL(mm,aa) = norminv(min([1-(0.5./stimpresentNL(end)) HitNL(end)./stimpresentNL(end)]))-norminv(max([0.5./stimabsentNL(end) FANL(end)./stimabsentNL(end)]));
        end
        if HitL(end)==0
            dprimeL(mm,aa) = norminv(1-(0.5./stimpresentL(end)))-norminv(max([0.5./stimabsentL(end) FAL(end)./stimabsentL(end)]));
        else
            dprimeL(mm,aa) = norminv(min([1-(0.5./stimpresentL(end)) HitL(end)./stimpresentL(end)]))-norminv(max([0.5./stimabsentL(end) FAL(end)./stimabsentL(end)]));
        end
        
    end
    RFs{mm} = RF;
    diffinpropshits_late{mm} = mean(HitNL(5:end)./stimpresentNL(5:end)-HitL(5:end)./stimpresentL(5:end));
    diffinpropshits_early{mm} = mean(HitNL(1:4)./stimpresentNL(1:4)-HitL(1:4)./stimpresentL(1:4));
    dprimeL_late{mm} = mean(norminv(min([1-(0.5./sum(stimpresentL(5:end))) sum(HitL(5:end))./sum(stimpresentL(5:end))],[],2))-norminv(max([0.5./sum(stimabsentL(5:end)) sum(FAL(5:end))./sum(stimabsentL(5:end))],[],2)));
    dprimeL_early{mm} = mean(norminv(min([1-(0.5./sum(stimpresentL(1:4))) sum(HitL(1:4))./sum(stimpresentL(1:4))],[],2))-norminv(max([0.5./sum(stimabsentL(1:4)) sum(FAL(1:4))./sum(stimabsentL(1:4))],[],2)));
    dprimeNL_late{mm} = mean(norminv(min([1-(0.5./sum(stimpresentNL(5:end))) sum(HitNL(5:end))./sum(stimpresentNL(5:end))],[],2))-norminv(max([0.5./sum(stimabsentNL(5:end)) sum(FANL(5:end))./sum(stimabsentNL(5:end))],[],2)));
    dprimeNL_early{mm} = mean(norminv(min([1-(0.5./sum(stimpresentNL(1:4))) sum(HitNL(1:4))./sum(stimpresentNL(1:4))],[],2))-norminv(max([0.5./sum(stimabsentNL(1:4)) sum(FANL(1:4))./sum(stimabsentNL(1:4))],[],2)));
end
% plot_counter = plot_counter + 5;
average_time_per_trials = mean(timeperblock./trialsperblock);
std_time_per_trials = std(timeperblock./trialsperblock);

files_in_eachsession = [numel(filenameopto1);numel(filenameopto2);numel(filenameopto3);numel(filenameopto4);numel(filenameopto5);numel(filenameopto6);numel(filenameopto7);numel(filenameopto8);numel(filenameopto9)];
avg_trials_withinsession = mean(files_in_eachsession)* mean(trialsperblock);
std_trials_withinsession = std(files_in_eachsession)* mean(trialsperblock);

std_fordiffs = [];
std_fordprimeL = [];
std_fordprimeNL = [];
std_fordprime_diffs = [];
for ii = 1:size(diffinpropshits,2)
    idx = diffinpropshits(:,ii)>0;
    std_fordiffs = [std_fordiffs; std(diffinpropshits(idx,ii))/sqrt(sum(idx))];
    std_fordprimeL = [std_fordprimeL; std(dprimeL(idx,ii))/sqrt(sum(idx))];
    std_fordprimeNL = [std_fordprimeNL; std(dprimeNL(idx,ii))/sqrt(sum(idx))];
    std_fordprime_diffs = [std_fordprime_diffs; std(dprimeNL(idx,ii)-dprimeL(idx,ii))/sqrt(sum(idx))];
end

% Quantifying the effect
% figure(plot_counter); set(gcf,'Name','Fetcsh analysis');
% subplot(241); errorbar(sum(diffinpropshits)./sum(diffinpropshits>0),std_fordiffs,'-o','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); axis square;
% xlabel('trial no.'); set(gca,'Ylim',[0 1],'YTick',[0:0.25:1],'Tickdir','out','Xlim',[1 11],'XTick',[1:2:11]);
% subplot(245); plot(cell2mat(diffinpropshits_early),cell2mat(diffinpropshits_late),'o','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
% plot([0 1],[0 1],'k'); set(gca,'Xlim',[0 1],'Ylim',[0 1],'Tickdir','out','XTick',0:0.25:1,'YTick',0:0.25:1); axis square; xlabel('control-laser 1-480 trials'); ylabel('control-laser 480+ trials');
% subplot(242); errorbar(sum(dprimeL)./sum(dprimeL~=0),std_fordprimeL,'-o','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); axis square;
% xlabel('trial no.'); ylabel('d prime laser'); set(gca,'Tickdir','out','Xlim',[1 11],'XTick',[1:2:11]);
% subplot(246); plot(cell2mat(dprimeL_early),cell2mat(dprimeL_late),'o','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
% plot([-2 2],[-2 2],'k');xlabel('d prime laser 1-480 trials'); ylabel('d prime laser 480+ trials'); set(gca,'Tickdir','out'); axis square;
% set(gca,'Tickdir','out','XTick',-2:2:2,'Xlim',[-2 2],'YTick',-2:2:2,'Ylim',[-2 2]); hold off;
% subplot(243); errorbar(sum(dprimeNL)./sum(dprimeNL~=0),std_fordprimeNL,'-o','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); axis square;
% xlabel('trial no.'); ylabel('d prime no-laser'); set(gca,'Tickdir','out','Xlim',[1 11],'XTick',[1:2:11]);
% subplot(247); plot(cell2mat(dprimeNL_early),cell2mat(dprimeNL_late),'o','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
% plot([-2 2],[-2 2],'k');xlabel('d prime no-laser 1-480 trials'); ylabel('d prime no-laser 480+ trials'); set(gca,'Tickdir','out'); axis square;
% set(gca,'Tickdir','out','XTick',-2:2:2,'Xlim',[-2 2],'YTick',-2:2:2,'Ylim',[-2 2]); hold off;
subplot(325); errorbar(sum(dprimeNL-dprimeL)./sum(dprimeNL~=0),std_fordprime_diffs,'-o','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); axis square;
xlabel('trial no.'); ylabel('d prime diff'); set(gca,'Ylim',[0 3],'YTick',[0:1:3],'Tickdir','out','Xlim',[1 11],'XTick',[1:2:11]);
subplot(326); plot(cell2mat(dprimeNL_early)-cell2mat(dprimeL_early),cell2mat(dprimeNL_late)-cell2mat(dprimeL_late),'o','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot([0 3.0],[0 3.0],'k');xlabel('d prime diff 1-480 trials'); ylabel('d prime diff 480+ trials'); set(gca,'Tickdir','out'); axis square;
set(gca,'Tickdir','out','XTick',0:1.5:3.0,'Xlim',[0 3.0],'YTick',0:1.5:3.0,'Ylim',[0 3.0]); hold off;
plot_counter = plot_counter + 1;

% Some statistical tests on the early vs late trials: diff in prop of hits, d prime diff analyses
p1 = signrank(cell2mat(diffinpropshits_early),cell2mat(diffinpropshits_late));
[~,p2] = ttest(cell2mat(diffinpropshits_early),cell2mat(diffinpropshits_late));
p3 = signrank(cell2mat(dprimeNL_early)-cell2mat(dprimeL_early),cell2mat(dprimeNL_late)-cell2mat(dprimeL_late));
[~,p4] = ttest(cell2mat(dprimeNL_early)-cell2mat(dprimeL_early),cell2mat(dprimeNL_late)-cell2mat(dprimeL_late));

% Correlation coefficient on the averaged data: diff in prop of hits, d prime diff analyses
[r5,p5] = corr((sum(diffinpropshits)./sum(diffinpropshits>0))',(1:1:size(diffinpropshits,2))','type','Spearman');
[r6,p6] = corr((sum(dprimeNL-dprimeL)./sum(dprimeNL~=0))',(1:1:size(diffinpropshits,2))','type','Spearman');


% Doing some regression analyses on indivisual sessions as suggested by Greg
regress_slopes_diff = [];
regress_slopes_dprimediff = [];
for ii = 1:size(diffinpropshits)
    beta = [ones(length(find(diffinpropshits(ii,:)>0)'),1) find(diffinpropshits(ii,:)>0)']\(diffinpropshits(ii,diffinpropshits(ii,:)>0))';
    beta1 = [ones(length(find(diffinpropshits(ii,:)>0)'),1) find(diffinpropshits(ii,:)>0)']\(dprimeNL(ii,diffinpropshits(ii,:)>0)-dprimeL(ii,diffinpropshits(ii,:)>0))';
    regress_slopes_diff = [regress_slopes_diff; beta(2)];
    regress_slopes_dprimediff = [regress_slopes_dprimediff; beta1(2)];
end
p7 = signrank(regress_slopes_diff);
[~,p8] = ttest(regress_slopes_diff);
p9 = signrank(regress_slopes_dprimediff);
[~,p10] = ttest(regress_slopes_dprimediff);

figure(plot_counter); set(gcf,'Name','Regression slopes');
subplot(121); histogram(regress_slopes_diff,-0.1:0.02:0.1,'FaceColor',[0 0 0]); hold on; plot(median(regress_slopes_diff),0,'kv','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
xlabel('Regress slopes diff'); ylabel('# sessions'); axis square; set(gca,'Ylim',[0 3],'YTick',0:1:3);
subplot(122); histogram(regress_slopes_dprimediff,-0.3:0.05:0.3,'FaceColor',[0 0 0]); hold on; plot(median(regress_slopes_dprimediff),0,'kv','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
xlabel('Regress slopes d prime diff'); ylabel('# sessions'); axis square; set(gca,'Ylim',[0 3],'YTick',0:1:3);
plot_counter = plot_counter + 1;



%% Supplementary figure 1
% Population part
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filename = fetch(conn,'SELECT filename FROM SCstimcue');
training = fetch(conn,'SELECT training FROM SCstimcue');
mode = fetch(conn,'SELECT mode FROM SCstimcue');
close(conn);
ind = find(strcmp(training,'no') & strcmp(mode,'SCstimcue'));
L = ceil(sqrt(numel(ind)));
probsacc_withinRF = cell(numel(ind),1);
probsacc_withinRFL = cell(numel(ind),1);
probsacc_withinRFNL = cell(numel(ind),1);
probsacc_outsideRF = cell(numel(ind),1);
probsacc_outsideRFL = cell(numel(ind),1);
probsacc_outsideRFNL = cell(numel(ind),1);
distancefromtarget_withinRFL = cell(numel(ind),1);
distancefromtarget_withinRFNL = cell(numel(ind),1);
distancefromtarget_outsideRFL = cell(numel(ind),1);
distancefromtarget_outsideRFNL = cell(numel(ind),1);
distfromtarget_NLcatch = cell(numel(ind),1);
distfromtarget_Lcatch = cell(numel(ind),1);

for ii = 1:numel(ind)
    targetlocations = [];
    targethitsallNL = []; % non-laser trials
    targethitsallL = []; % laser trials
    stimpresentationNL = [];% non-laser trials
    stimpresentationL = []; % laser trials
    distfromtargetL = []; % distance between final eye position and RF in laser trials
    distfromtargetNL = [];% distance between final eye position and RF in control trials
    RFloc = [];
    fileofinterest = char(filename(ind(ii),:));
    
    stro = nex2stro(findfile(fileofinterest));
    targ_x = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_x'));
    targ_y = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_y'));
    uniquetargxy = unique([targ_x targ_y],'rows');
    ntrials = size(stro.trial,1);
    Lcatchtrials = targ_x == 0 & targ_y == 0;
    laseron_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
    laseroff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
    fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
    sacint_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'saccinit_t'));
    sacmade_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'saccmade_t'));
    targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targon_t'));
    targoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targoff_t'));
    Llaser = ~isnan(laseron_t);
    samplerate = stro.sum.analog.storeRates{1};
    [C,ia,ib] = unique(stro.trial(:,12:13),'rows');
    targetlocations = [targetlocations; C];
    targethitsNL = zeros(size(C,1),1);
    targethitsL = zeros(size(C,1),1);
    SPNL = zeros(size(C,1),1);
    SPL = zeros(size(C,1),1);
    tmp_distfromtargetNLcatch = [];
    tmp_distfromtargetLcatch = [];
    RFloc = [RFloc; stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y];
    if norm(RFloc)>100
        halfwidth = 4.0; % in degrees of visual angle
    else
        halfwidth = 1.5; % in degrees of visual angle
    end
    for i = 1:size(stro.trial,1)
        if size(stro.ras,2)==5
            x = stro.ras{i,2}*4096/400;
            y = stro.ras{i,3}*4096/400;
            t = stro.ras{i,5}+[0:1:length(x)-1]/samplerate;
        elseif size(stro.ras,2)==6
            x = stro.ras{i,3}*4096/400;
            y = stro.ras{i,4}*4096/400;
            t = stro.ras{i,6}+[0:1:length(x)-1]/samplerate;
        elseif size(stro.ras,2)==7
            x = stro.ras{i,4}*4096/400;
            y = stro.ras{i,5}*4096/400;
            t = stro.ras{i,7}+[0:1:length(x)-1]/samplerate;
        end
        if ~Llaser(i) & ~Lcatchtrials(i)
            
            Lt = t>fpoff_t(i) & t < fpoff_t(i)+0.30;
            target = stro.trial(i,[12 13]); % x and y target locations of that trial
            idx = find(sum(target==C,2)==2);
            if sum(x(Lt)>(target(1)/10)-halfwidth & x(Lt)<(target(1)/10)+halfwidth & y(Lt)>(target(2)/10)-halfwidth & y(Lt)<(target(2)/10)+halfwidth)>10
                targethitsNL(idx) = targethitsNL(idx) + 1;
            end
            SPNL(idx) = SPNL(idx) + 1;
            finaleyepos = [x(find(Lt,1,'last')) y(find(Lt,1,'last'))];
            distfromtargetNL(idx) = norm((target/10) - finaleyepos);
        end
        if Llaser(i) & ~Lcatchtrials(i) % Analyzing laser trials
            
            Lt = t>fpoff_t(i) & t < fpoff_t(i)+0.30;
            target = stro.trial(i,[12 13]); % x and y target locations of that trial
            idx = find(sum(target==C,2)==2);
            if sum(x(Lt)>(target(1)/10)-halfwidth & x(Lt)<(target(1)/10)+halfwidth & y(Lt)>(target(2)/10)-halfwidth & y(Lt)<(target(2)/10)+halfwidth)>10
                targethitsL(idx) = targethitsL(idx) + 1;
            end
            analogstartime = stro.ras{i,5};
            spiketimes = stro.ras{i,1};
            laserontime = laseron_t(i);
            laserofftime = laseroff_t(i);
            timedurlaser = laserofftime - laserontime;
            SPL(idx) = SPL(idx) + 1;
            finaleyepos = [x(find(Lt,1,'last')) y(find(Lt,1,'last'))];
            distfromtargetL(idx) = norm((target/10) - finaleyepos);
        end
        % Now calculating results for catch trials
        if ~Llaser(i) & Lcatchtrials(i)
            
            Lt = t>fpoff_t(i) & t < fpoff_t(i)+0.30;
            finaleyepos = [x(find(Lt,1,'last')) y(find(Lt,1,'last'))];
            tmp_distfromtargetNLcatch = [tmp_distfromtargetNLcatch; norm((target/10) - finaleyepos)];
        end
        if Llaser(i) & Lcatchtrials(i)
            
            Lt = t>fpoff_t(i) & t < fpoff_t(i)+0.30;
            finaleyepos = [x(find(Lt,1,'last')) y(find(Lt,1,'last'))];
            tmp_distfromtargetLcatch = [tmp_distfromtargetLcatch; norm((target/10) - finaleyepos)];
        end
    end
    targethitsallNL = [targethitsallNL; targethitsNL];
    targethitsallL = [targethitsallL; targethitsL];
    stimpresentationNL = [stimpresentationNL; SPNL];
    stimpresentationL = [stimpresentationL; SPL];
    withinRFidx = all(C==RFloc,2);
    outsideRFidx = all(C~=RFloc & C~=[0 0],2);
    probsacc_withinRFL{ii} = targethitsallL(withinRFidx)./stimpresentationL(withinRFidx);
    probsacc_withinRFNL{ii} = targethitsallNL(withinRFidx)./stimpresentationNL(withinRFidx);
    probsacc_outsideRFL{ii} = targethitsallL(outsideRFidx)./stimpresentationL(outsideRFidx);
    probsacc_outsideRFNL{ii} = targethitsallNL(outsideRFidx)./stimpresentationNL(outsideRFidx);
    probsacc_withinRF{ii} = (targethitsallNL(withinRFidx)./stimpresentationNL(withinRFidx)) -(targethitsallL(withinRFidx)./stimpresentationL(withinRFidx));
    probsacc_outsideRF{ii} = (targethitsallNL(outsideRFidx)./stimpresentationNL(outsideRFidx)) -(targethitsallL(outsideRFidx)./stimpresentationL(outsideRFidx));
    
    distancefromtarget_withinRFL{ii} = mean(distfromtargetL(withinRFidx));
    distancefromtarget_withinRFNL{ii} = mean(distfromtargetNL(withinRFidx));
    distancefromtarget_outsideRFL{ii} = mean(distfromtargetL(outsideRFidx));
    distancefromtarget_outsideRFNL{ii} = mean(distfromtargetNL(outsideRFidx));
    distfromtarget_NLcatch{ii} = mean(tmp_distfromtargetNLcatch);
    distfromtarget_Lcatch{ii} = mean(tmp_distfromtargetLcatch);
end

k = cell2mat(filename(ind));
monkey_ID = k(:,1);
idx = monkey_ID == 'M';
bins = -1:0.1:1;
% figure(plot_counter); set(gcf,'Name','Population: SCstimcue');
% h1 = histogram(cell2mat(probsacc_withinRF(idx)),bins,'Normalization','probability','Visible','off');val1 = h1.Values;
% h2 = histogram(cell2mat(probsacc_outsideRF(idx)),bins,'Normalization','probability','Visible','off');val2 = h2.Values;
% h3 = histogram(cell2mat(probsacc_withinRF(~idx)),bins,'Normalization','probability','Visible','off');val3 = h3.Values;
% h4 = histogram(cell2mat(probsacc_outsideRF(~idx)),bins,'Normalization','probability','Visible','off');val4 = h4.Values;
% subplot(321); bar([val1;val2]','stacked'); axis square; set(gca,'Tickdir','out','Xlim',[0.5 20.5]);
% subplot(322); bar([val3;val4]','stacked'); axis square; set(gca,'Tickdir','out','Xlim',[0.5 20.5]);
% subplot(323); histogram(cell2mat(probsacc_withinRF(idx)),bins,'Normalization','probability','FaceColor',[0 0 0],'EdgeColor',[1 1 1]); hold on; histogram(cell2mat(probsacc_outsideRF(idx)),bins,'Normalization','probability','FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1]);
% xlabel('P(sacc control - laser)'); ylabel('pdf'); axis square; set(gca,'Tickdir','out','Ylim',[0 0.6],'YTick',0:0.3:0.6); title('M');
% subplot(324); histogram(cell2mat(probsacc_withinRF(~idx)),bins,'Normalization','probability','FaceColor',[0 0 0],'EdgeColor',[1 1 1]); hold on; histogram(cell2mat(probsacc_outsideRF(~idx)),bins,'Normalization','probability','FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1]);
% xlabel('P(sacc control - laser)'); ylabel('pdf'); axis square; set(gca,'Tickdir','out','Ylim',[0 0.6],'YTick',0:0.3:0.6); title('A');
% set(gcf,'renderer','painters');
% subplot(325),plot(cell2mat(distancefromtarget_withinRFL(idx)),cell2mat(distancefromtarget_withinRFNL(idx)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
% plot(cell2mat(distancefromtarget_outsideRFL(idx)),cell2mat(distancefromtarget_outsideRFNL(idx)),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); xlabel('Laser'); ylabel('Control'); axis square;
% set(gca,'Tickdir','out','Xlim',[0 10],'Ylim',[0 10],'XTick',0:5:10,'YTick',0:5:10); title('Maui'); line([0 10],[0 10]); hold off;
% subplot(326),plot(cell2mat(distancefromtarget_withinRFL(~idx)),cell2mat(distancefromtarget_withinRFNL(~idx)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
% plot(cell2mat(distancefromtarget_outsideRFL(~idx)),cell2mat(distancefromtarget_outsideRFNL(~idx)),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); xlabel('Laser'); ylabel('Control'); axis square;
% set(gca,'Tickdir','out','Xlim',[0 20],'Ylim',[0 20],'XTick',0:10:20,'YTick',0:10:20); title('Apollo'); line([0 20],[0 20]); hold off;
% plot_counter = plot_counter + 1;

[p3] = signrank(cell2mat(probsacc_withinRFL(idx)),cell2mat(probsacc_withinRFNL(idx))); % For Maui: within RF
[p4] = signrank(cell2mat(probsacc_outsideRFL(idx)),cell2mat(probsacc_outsideRFNL(idx))); % For Maui: outside RF
[p5] = signrank(cell2mat(probsacc_withinRFL(~idx)),cell2mat(probsacc_withinRFNL(~idx))); % For Apollo: within RF
[p6] = signrank(cell2mat(probsacc_outsideRFL(~idx)),cell2mat(probsacc_outsideRFNL(~idx))); % For Apollo: outside RF

% A new analyses based on Greg's suggestion: difference between the centroids of laser and control target
% figure(plot_counter); set(gcf,'Name','Diff between centroids: eye end points');
% subplot(121); plot(abs(cell2mat(distancefromtarget_withinRFL(idx))-cell2mat(distancefromtarget_withinRFNL(idx))),abs(cell2mat(distancefromtarget_outsideRFL(idx))-cell2mat(distancefromtarget_outsideRFNL(idx))),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
% xlabel('Within RF'); ylabel('Outside RF'); title('Diff in centroids: Maui'); set(gca,'Tickdir','out','Xlim',[0 10],'Ylim',[0 10],'XTick',0:2.5:10,'YTick',0:2.5:10); line([0 10],[0 10]); axis square;
% subplot(122); plot(abs(cell2mat(distancefromtarget_withinRFL(~idx))-cell2mat(distancefromtarget_withinRFNL(~idx))),abs(cell2mat(distancefromtarget_outsideRFL(~idx))-cell2mat(distancefromtarget_outsideRFNL(~idx))),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
% xlabel('Within RF'); ylabel('Outside RF'); title('Diff in centroids: Apollo'); set(gca,'Tickdir','out','Xlim',[0 10],'Ylim',[0 10],'XTick',0:2.5:10,'YTick',0:2.5:10); line([0 10],[0 10]); axis square;
% plot_counter = plot_counter + 1;

p7 = signrank(abs(cell2mat(distancefromtarget_withinRFL(idx))-cell2mat(distancefromtarget_withinRFNL(idx))),abs(cell2mat(distancefromtarget_outsideRFL(idx))-cell2mat(distancefromtarget_outsideRFNL(idx))));
[~,p8] = ttest(abs(cell2mat(distancefromtarget_withinRFL(idx))-cell2mat(distancefromtarget_withinRFNL(idx))),abs(cell2mat(distancefromtarget_outsideRFL(idx))-cell2mat(distancefromtarget_outsideRFNL(idx))));
p9 = signrank(abs(cell2mat(distancefromtarget_withinRFL(~idx))-cell2mat(distancefromtarget_withinRFNL(~idx))),abs(cell2mat(distancefromtarget_outsideRFL(~idx))-cell2mat(distancefromtarget_outsideRFNL(~idx))));
[~,p10] = ttest(abs(cell2mat(distancefromtarget_withinRFL(~idx))-cell2mat(distancefromtarget_withinRFNL(~idx))),abs(cell2mat(distancefromtarget_outsideRFL(~idx))-cell2mat(distancefromtarget_outsideRFNL(~idx))));

% Population analysis for catch trials for Maui and Apollo
pcatch_M = signrank(cell2mat(distfromtarget_Lcatch(idx)),cell2mat(distfromtarget_NLcatch(idx))); % Maui
pcatch_A = signrank(cell2mat(distfromtarget_Lcatch(~idx)),cell2mat(distfromtarget_NLcatch(~idx))); % Apollo

% Trying to work with getSacdata
tmax = 1.0;
nrows = 5;
count = 1;
RFloc = [];
diffcentroidwithin = [];
diffcentroidoutside = [];
Saccadelatency_withinL = [];
Saccadelatency_withinNL = [];
Saccadelatency_outsideL = [];
Saccadelatency_outsideNL = [];
SaccadeAmplitudes_all = [];
SaccadeAmplitudes_catchL = [];
SaccadeAmplitudes_catchNL = [];
SaccadeAmplitudes_catchL = [];
SaccadecatchXY_L = [];
SaccadecatchXY_NL = [];
centroidcatchL = [];
centroidcatchNL = [];
RF = [];
plotfigures = 0;
Catchtrials_exist = [];
Saccadelatency_catchL = [];
Saccadelatency_catchNL = [];
for ii = 1:numel(ind)
    if ii == find(idx==0,1)
        Saccadelatency_withinLidx = numel(Saccadelatency_withinL);
        Saccadelatency_withinNLidx = numel(Saccadelatency_withinNL);
        Saccadelatency_outsideLidx = numel(Saccadelatency_outsideL);
        Saccadelatency_outsideNLidx = numel(Saccadelatency_outsideNL);
        Saccadelatency_catchLidx = numel(Saccadelatency_catchL);
        Saccadelatency_catchNLidx = numel(Saccadelatency_catchNL);
        SaccadeAmplitudes_catchLidx = numel(SaccadeAmplitudes_catchL);
        SaccadeAmplitudes_catchNLidx = numel(SaccadeAmplitudes_catchNL);
    end
    fileofinterest = char(filename(ind(ii),:));
    stro = nex2stro(findfile(fileofinterest));
    fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
    SacData = getSacData(stro,1,Inf,10);
    targ_x = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_x'));
    targ_y = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_y'));
    RFloc = [stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y];
    RF = [RF; RFloc];
    Lcatchtrials = targ_x == 0 & targ_y == 0;
    laseron_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
    laseroff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
    fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
    [C,ia,ib] = unique(stro.trial(:,12:13),'rows');
    SPNL = zeros(size(C,1),1);
    SPL = zeros(size(C,1),1);
    
    Llaser = ~isnan(laseron_t);
    SaccadeAmplitudes = [];
    SaccadeDirection = [];
    SaccadeX = [];
    SaccadeY = [];
    Saccadelatency = [];
    
    if idx(ii)
        plotlimit = 8;
    else
        plotlimit = 15;
    end
    for jj = 1:length(SacData.amplitudes)
        fpofftime = fpoff_t(jj);
        
        if ~isempty(find(SacData.starttimes{jj}-fpofftime>0 & SacData.starttimes{jj}-fpofftime<tmax))
            selected_times_idxs = find(SacData.starttimes{jj}-fpofftime>0 & SacData.starttimes{jj}-fpofftime<tmax);
            %             [~,ix] = max(SacData.amplitudes{jj}(selected_times_idxs)); % contingent on saccade amplitude
            ix = 1; % contingent on first saccade
            SaccadeAmplitudes = [SaccadeAmplitudes; max(SacData.amplitudes{jj}(selected_times_idxs(ix)))];
            SaccadeDirection = [SaccadeDirection; max(SacData.directions{jj}(selected_times_idxs(ix)))];
            SaccadeX = [SaccadeX; SaccadeAmplitudes(end)*cos(SaccadeDirection(end))];
            SaccadeY = [SaccadeY; SaccadeAmplitudes(end)*sin(SaccadeDirection(end))];
            Saccadelatency = [Saccadelatency; SacData.starttimes{jj}(selected_times_idxs(ix))-fpofftime];
        else
            SaccadeAmplitudes = [SaccadeAmplitudes; 0];
            SaccadeDirection = [SaccadeDirection; 0];
            SaccadeX = [SaccadeX; 0];
            SaccadeY = [SaccadeY; 0];
            Saccadelatency = [Saccadelatency; tmax];
        end
    end
    
    SaccadeAmplitudes_all = [SaccadeAmplitudes_all; SaccadeAmplitudes];
    withinRFidx = find(all(C==RFloc,2));
    outsideRFidx = find(all(C~=RFloc & C~=[0 0],2));
    centroidwithinL = [mean(SaccadeX(ib==withinRFidx & Llaser)) mean(SaccadeY(ib==withinRFidx & Llaser))];
    centroidwithinNL = [mean(SaccadeX(ib==withinRFidx & ~Llaser)) mean(SaccadeY(ib==withinRFidx & ~Llaser))];
    diffcentroidwithin = [diffcentroidwithin; sum((centroidwithinL - centroidwithinNL).^2)];
    Saccadelatency_withinL = [Saccadelatency_withinL; Saccadelatency(ib==withinRFidx & Llaser)];
    Saccadelatency_withinNL = [Saccadelatency_withinNL; Saccadelatency(ib==withinRFidx & ~Llaser)];
    
    centroidoutsideL = []; centroidoutsideNL = [];
    outsideL = []; outsideNL = [];
    for jj = 1:numel(outsideRFidx)
        centroidoutsideL = [centroidoutsideL; mean(SaccadeX(ib==outsideRFidx(jj) & Llaser)) mean(SaccadeY(ib==outsideRFidx(jj) & Llaser))];
        centroidoutsideNL = [centroidoutsideNL; mean(SaccadeX(ib==outsideRFidx(jj) & ~Llaser)) mean(SaccadeY(ib==outsideRFidx(jj) & ~Llaser))];
        outsideL = [outsideL; Saccadelatency(ib==outsideRFidx(jj) & Llaser)];
        outsideNL = [outsideNL; Saccadelatency(ib==outsideRFidx(jj) & ~Llaser)];
    end
    diffcentroidoutside = [diffcentroidoutside; mean(sum((centroidoutsideL - centroidoutsideNL).^2),2)];
    Saccadelatency_outsideL = [Saccadelatency_outsideL; outsideL];
    Saccadelatency_outsideNL = [Saccadelatency_outsideNL; outsideNL];
    
    if any(Lcatchtrials & Llaser)
        Catchtrials_exist = [Catchtrials_exist; 1];
        centroidcatchL = [centroidcatchL; norm([mean(SaccadeX(Lcatchtrials & Llaser)) mean(SaccadeY(Lcatchtrials & Llaser))])];
        centroidcatchNL = [centroidcatchNL;  norm([mean(SaccadeX(Lcatchtrials & ~Llaser)) mean(SaccadeY(Lcatchtrials & ~Llaser))])];
        Saccadelatency_catchL = [Saccadelatency_catchL; Saccadelatency(Lcatchtrials & Llaser)];
        Saccadelatency_catchNL = [Saccadelatency_catchNL; Saccadelatency(Lcatchtrials & ~Llaser)];
        SaccadeAmplitudes_catchL = [SaccadeAmplitudes_catchL; SaccadeAmplitudes(Lcatchtrials & Llaser)];
        SaccadeAmplitudes_catchNL = [SaccadeAmplitudes_catchNL; SaccadeAmplitudes(Lcatchtrials & ~Llaser)];
        SaccadecatchXY_L = [SaccadecatchXY_L; SaccadeX(Lcatchtrials & Llaser) SaccadeY(Lcatchtrials & Llaser)];
        SaccadecatchXY_NL = [SaccadecatchXY_NL; SaccadeX(Lcatchtrials & ~Llaser) SaccadeY(Lcatchtrials & ~Llaser)];
    else
        Catchtrials_exist = [Catchtrials_exist; 0];
    end
    
    if plotfigures
        figure(plot_counter),subplot(nrows,2,2*count-1); plot(upsample(SaccadeX(~Llaser & ~Lcatchtrials),2),upsample(SaccadeY(~Llaser & ~Lcatchtrials),2),'Color',[0.5 0.5 0.5]); hold on;
        plot(upsample(SaccadeX(Llaser & ~Lcatchtrials),2),upsample(SaccadeY(Llaser & ~Lcatchtrials),2),'Color',[0 0.5 1.0]); title(fileofinterest); set(gca,'Tickdir','out','Xlim',[-plotlimit plotlimit],'Ylim',[-plotlimit plotlimit],'XTick',-plotlimit:plotlimit:plotlimit,'YTick',-plotlimit:plotlimit:plotlimit); axis square; hold off;
        subplot(nrows,2,2*count);
        if any(~Llaser & Lcatchtrials)
            plot(upsample(SaccadeX(~Llaser & Lcatchtrials),2),upsample(SaccadeY(~Llaser & Lcatchtrials),2),'Color',[0.5 0.5 0.5]); hold on;
        end
        if any(Llaser & Lcatchtrials)
            plot(upsample(SaccadeX(Llaser & Lcatchtrials),2),upsample(SaccadeY(Llaser & Lcatchtrials),2),'Color',[0 0.5 1.0]);
        end
        set(gca,'Tickdir','out','Xlim',[-plotlimit plotlimit],'Ylim',[-plotlimit plotlimit],'XTick',-plotlimit:plotlimit:plotlimit,'YTick',-plotlimit:plotlimit:plotlimit); axis square; hold off;
        count = count + 1;
        if count == nrows+1
            plot_counter = plot_counter + 1;
            count = 1;
        end
    end
end
plot_counter = plot_counter + 1;

% Analyzing trials when target was presented
figure(plot_counter); set(gcf,'Name','getSaccData: Saccade amplitudes and latencies');
subplot(321); plot(diffcentroidwithin(idx)/10,diffcentroidoutside(idx)/10,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
xlabel('Within RF'); ylabel('Outside RF'); title('Diff in centroids: Maui'); set(gca,'Tickdir','out','Xlim',[0 6],'Ylim',[0 6],'XTick',0:2:6,'YTick',0:2:6); line([0 6],[0 6]); axis square;
subplot(322); plot(diffcentroidwithin(~idx)/10,diffcentroidoutside(~idx)/10,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
xlabel('Within RF'); ylabel('Outside RF'); title('Diff in centroids: Apollo'); set(gca,'Tickdir','out','Xlim',[0 6],'Ylim',[0 6],'XTick',0:2:6,'YTick',0:2:6); line([0 6],[0 6]); axis square;
subplot(323); histogram(Saccadelatency_withinNL(1:Saccadelatency_withinNLidx),0:0.025:1.0,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5,'Normalization','probability'); hold on; histogram(Saccadelatency_withinL(1:Saccadelatency_withinLidx),0:0.025:1.0,'FaceColor',[0 0.5 1.0],'Normalization','probability');
xlabel('Latencies (ms)'); ylabel('% saccades'); title('Latency within RF: Maui'); set(gca,'Tickdir','out','Xlim',[0 1],'XTick',0:0.25:1,'Ylim',[0 0.8],'YTick',0.:0.2:0.8); legend('control','laser'); axis square;
subplot(325); histogram(Saccadelatency_outsideNL(1:Saccadelatency_outsideNLidx),0:0.025:1.0,'FaceColor',[0.5 0.5 0.5],'Normalization','probability'); hold on; histogram(Saccadelatency_outsideL(1:Saccadelatency_outsideNLidx),0:0.025:1.0,'FaceColor',[0 0.5 1.0],'Normalization','probability');
xlabel('Latencies (ms)'); ylabel('% saccades');  title('Latency outside RF: Maui'); set(gca,'Tickdir','out','Xlim',[0 1],'XTick',0:0.25:1,'Ylim',[0 0.8],'YTick',0.:0.2:0.8); axis square;
subplot(324); histogram(Saccadelatency_withinNL(Saccadelatency_withinNLidx+1:end),0:0.025:1.0,'FaceColor',[0.5 0.5 0.5],'Normalization','probability'); hold on; histogram(Saccadelatency_withinL(Saccadelatency_withinLidx+1:end),0:0.025:1.0,'FaceColor',[0 0.5 1.0],'Normalization','probability');
xlabel('Latencies (ms)'); ylabel('% saccades'); title('Latency within RF: Apollo'); set(gca,'Tickdir','out','Xlim',[0 1],'XTick',0:0.25:1,'Ylim',[0 0.8],'YTick',0.:0.2:0.8); axis square;
subplot(326); histogram(Saccadelatency_outsideNL(Saccadelatency_outsideNLidx+1:end),0:0.025:1.0,'FaceColor',[0.5 0.5 0.5],'Normalization','probability'); hold on; histogram(Saccadelatency_outsideL(Saccadelatency_outsideNLidx+1:end),0:0.025:1.0,'FaceColor',[0 0.5 1.0],'Normalization','probability');
xlabel('Latencies (ms)'); ylabel('% saccades');  title('Latency outside RF: Apollo'); set(gca,'Tickdir','out','Xlim',[0 1],'XTick',0:0.25:1,'Ylim',[0 0.8],'YTick',0.:0.2:0.8); axis square;
plot_counter = plot_counter + 1;

% Same as before but cutting off at t=0.5; Analyzing trials when target was presented
Saccadelatency_withinNL_mod = Saccadelatency_withinNL; Saccadelatency_withinNL_mod(Saccadelatency_withinNL_mod>0.5)=0.5;
Saccadelatency_withinL_mod = Saccadelatency_withinL; Saccadelatency_withinL_mod(Saccadelatency_withinL_mod>0.5)=0.5;
Saccadelatency_outsideNL_mod = Saccadelatency_outsideNL; Saccadelatency_outsideNL_mod(Saccadelatency_outsideNL_mod>0.5)=0.5;
Saccadelatency_outsideL_mod = Saccadelatency_outsideL; Saccadelatency_outsideL_mod(Saccadelatency_outsideL_mod>0.5)=0.5;
figure(plot_counter); set(gcf,'Name','getSaccData: Saccade amplitudes and latencies');
subplot(321); plot(diffcentroidwithin(idx)/10,diffcentroidoutside(idx)/10,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
xlabel('Within RF'); ylabel('Outside RF'); title('Diff in centroids: Maui'); set(gca,'Tickdir','out','Xlim',[0 5],'Ylim',[0 5],'XTick',0:2.5:5,'YTick',0:2.5:5); line([0 5],[0 5]); axis square;
subplot(322); plot(diffcentroidwithin(~idx)/10,diffcentroidoutside(~idx)/10,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
xlabel('Within RF'); ylabel('Outside RF'); title('Diff in centroids: Apollo'); set(gca,'Tickdir','out','Xlim',[0 5],'Ylim',[0 5],'XTick',0:2.5:5,'YTick',0:2.5:5); line([0 5],[0 5]); axis square;
subplot(323); histogram(Saccadelatency_withinNL_mod(1:Saccadelatency_withinNLidx),0:0.025:0.5,'EdgeColor',[0.5 0.5 0.5],'Normalization','probability','DisplayStyle','stairs'); hold on; histogram(Saccadelatency_withinL_mod(1:Saccadelatency_withinLidx),0:0.025:0.5,'EdgeColor',[0 0.5 1.0],'Normalization','probability','DisplayStyle','stairs');
xlabel('Latencies (ms)'); ylabel('% saccades'); title('Latency within RF: Maui'); set(gca,'Tickdir','out','Xlim',[0 0.5],'XTick',0:0.1:0.5,'Ylim',[0 0.8],'YTick',0.:0.2:0.8); legend('control','laser'); axis square;
subplot(325); histogram(Saccadelatency_outsideNL_mod(1:Saccadelatency_outsideNLidx),0:0.025:0.5,'EdgeColor',[0.5 0.5 0.5],'Normalization','probability','DisplayStyle','stairs'); hold on; histogram(Saccadelatency_outsideL_mod(1:Saccadelatency_outsideNLidx),0:0.025:0.5,'EdgeColor',[0 0.5 1.0],'Normalization','probability','DisplayStyle','stairs');
xlabel('Latencies (ms)'); ylabel('% saccades');  title('Latency outside RF: Maui'); set(gca,'Tickdir','out','Xlim',[0 0.5],'XTick',0:0.1:0.5,'Ylim',[0 0.8],'YTick',0.:0.2:0.8); axis square;
subplot(324); histogram(Saccadelatency_withinNL_mod(Saccadelatency_withinNLidx+1:end),0:0.025:0.5,'EdgeColor',[0.5 0.5 0.5],'Normalization','probability','DisplayStyle','stairs'); hold on; histogram(Saccadelatency_withinL_mod(Saccadelatency_withinLidx+1:end),0:0.025:0.5,'EdgeColor',[0 0.5 1.0],'Normalization','probability','DisplayStyle','stairs');
xlabel('Latencies (ms)'); ylabel('% saccades'); title('Latency within RF: Apollo'); set(gca,'Tickdir','out','Xlim',[0 0.5],'XTick',0:0.1:0.5,'Ylim',[0 0.8],'YTick',0.:0.2:0.8); axis square;
subplot(326); histogram(Saccadelatency_outsideNL_mod(Saccadelatency_outsideNLidx+1:end),0:0.025:0.5,'EdgeColor',[0.5 0.5 0.5],'Normalization','probability','DisplayStyle','stairs'); hold on; histogram(Saccadelatency_outsideL_mod(Saccadelatency_outsideNLidx+1:end),0:0.025:0.5,'EdgeColor',[0 0.5 1.0],'Normalization','probability','DisplayStyle','stairs');
xlabel('Latencies (ms)'); ylabel('% saccades');  title('Latency outside RF: Apollo'); set(gca,'Tickdir','out','Xlim',[0 0.5],'XTick',0:0.1:0.5,'Ylim',[0 0.8],'YTick',0.:0.2:0.8); axis square;
plot_counter = plot_counter + 1;

% Trying the bar graph representation
bins = 0:0.1:0.6;
h1 = histcounts(Saccadelatency_withinNL(1:Saccadelatency_withinNLidx),[bins Inf]); h1 = h1/sum(h1);
h2 = histcounts(Saccadelatency_withinL(1:Saccadelatency_withinLidx),[bins Inf]); h2 = h2/sum(h2);
h3 = histcounts(Saccadelatency_outsideNL(1:Saccadelatency_outsideNLidx),[bins Inf]); h3 = h3/sum(h3);
h4 = histcounts(Saccadelatency_outsideL(1:Saccadelatency_outsideLidx),[bins Inf]); h4 = h4/sum(h4);
h5 = histcounts(Saccadelatency_withinNL(Saccadelatency_withinNLidx+1:end),[bins Inf]); h5 = h5/sum(h5);
h6 = histcounts(Saccadelatency_withinL(Saccadelatency_withinLidx+1:end),[bins Inf]); h6 = h6/sum(h6);
h7 = histcounts(Saccadelatency_outsideNL(Saccadelatency_outsideNLidx+1:end),[bins Inf]); h7 = h7/sum(h7);
h8 = histcounts(Saccadelatency_outsideL(Saccadelatency_outsideLidx+1:end),[bins Inf]); h8 = h8/sum(h8);
figure(plot_counter);
subplot(321); plot(diffcentroidwithin(idx)/10,diffcentroidoutside(idx)/10,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
xlabel('Within RF'); ylabel('Outside RF'); title('Diff in centroids: Maui'); set(gca,'Tickdir','out','Xlim',[0 5],'Ylim',[0 5],'XTick',0:2.5:5,'YTick',0:2.5:5); line([0 5],[0 5]); axis square;
subplot(322); plot(diffcentroidwithin(~idx)/10,diffcentroidoutside(~idx)/10,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
xlabel('Within RF'); ylabel('Outside RF'); title('Diff in centroids: Apollo'); set(gca,'Tickdir','out','Xlim',[0 5],'Ylim',[0 5],'XTick',0:2.5:5,'YTick',0:2.5:5); line([0 5],[0 5]); axis square;
subplot(323); bar(bins,[h1;h2]'); set(gca,'Tickdir','out','Ylim',[0 1],'YTick',0:0.25:1.0,'Xlim',[-0.1 0.7],'XTick',-0.1:0.1:0.7); xlabel('Latency');
title('Latency within RF: Maui'); ylabel('proportion of saccades'); axis square;
subplot(324); bar(bins,[h3;h4]'); set(gca,'Tickdir','out','Ylim',[0 1],'YTick',0:0.25:1.0,'Xlim',[-0.1 0.7],'XTick',-0.1:0.1:0.7); xlabel('Latency');
title('Latency outside RF: Maui'); ylabel('proportion of saccades'); axis square;
subplot(325); bar(bins,[h5;h6]'); set(gca,'Tickdir','out','Ylim',[0 1],'YTick',0:0.25:1.0,'Xlim',[-0.1 0.7],'XTick',-0.1:0.1:0.7); xlabel('Latency');
title('Latency within RF: Apollo'); ylabel('proportion of saccades'); axis square;
subplot(326); bar(bins,[h7;h8]'); set(gca,'Tickdir','out','Ylim',[0 1],'YTick',0:0.25:1.0,'Xlim',[-0.1 0.7],'XTick',-0.1:0.1:0.7); xlabel('Latency');
title('Latency outside RF: Apollo'); ylabel('proportion of saccades'); axis square;
plot_counter = plot_counter + 1;

% Analyzing catch trials: when no target was presented
figure(plot_counter), set(gcf,'Name','getSaccData: Catch trial amplitudes and latencies');
subplot(221); histogram(SaccadeAmplitudes_catchNL(1:SaccadeAmplitudes_catchNLidx),0:1:15,'FaceColor',[0.5 0.5 0.5],'Normalization','probability'); hold on;  histogram(SaccadeAmplitudes_catchL(1:SaccadeAmplitudes_catchLidx),0:1:15,'FaceColor',[0 0.5 1.0],'Normalization','probability');
xlabel('Amplitude'); ylabel('% saccades');  title('Amplitude catch: Maui'); set(gca,'Tickdir','out','Xlim',[0 15],'XTick',0:5:15); axis square;
subplot(222); histogram(SaccadeAmplitudes_catchNL(SaccadeAmplitudes_catchNLidx+1:end),0:1:15,'FaceColor',[0.5 0.5 0.5],'Normalization','probability'); hold on;  histogram(SaccadeAmplitudes_catchL(SaccadeAmplitudes_catchLidx+1:end),0:1:15,'FaceColor',[0 0.5 1.0],'Normalization','probability');
xlabel('Amplitude'); ylabel('% saccades');  title('Amplitude catch: Apollo'); set(gca,'Tickdir','out','Xlim',[0 15],'XTick',0:5:15); axis square;
subplot(223); histogram(Saccadelatency_catchNL(1:Saccadelatency_catchNLidx),0:0.05:1.0,'FaceColor',[0.5 0.5 0.5],'Normalization','probability'); hold on; histogram(Saccadelatency_catchL(1:Saccadelatency_catchNLidx),0:0.05:1.0,'FaceColor',[0 0.5 1.0],'Normalization','probability');
xlabel('Latencies (ms)'); ylabel('% saccades');  title('Latency catch: Maui'); set(gca,'Tickdir','out','Xlim',[0 1],'XTick',0:0.25:1); axis square;
subplot(224); histogram(Saccadelatency_catchNL(Saccadelatency_catchNLidx+1:end),0:0.05:1.0,'FaceColor',[0.5 0.5 0.5],'Normalization','probability'); hold on; histogram(Saccadelatency_catchL(Saccadelatency_catchNLidx+1:end),0:0.05:1.0,'FaceColor',[0 0.5 1.0],'Normalization','probability');
xlabel('Latencies (ms)'); ylabel('% saccades');  title('Latency catch: Apollo'); set(gca,'Tickdir','out','Xlim',[0 1],'XTick',0:0.25:1); axis square;
plot_counter = plot_counter + 1;

% Statistical tests for saccade amplitudes
p1 = signrank(diffcentroidwithin(idx),diffcentroidoutside(idx));
[~,p2] = ttest(diffcentroidwithin(idx),diffcentroidoutside(idx));
p3 = signrank(diffcentroidwithin(~idx),diffcentroidoutside(~idx));
[~,p4] = ttest(diffcentroidwithin(~idx),diffcentroidoutside(~idx));

% Statistical tests for saccade latencies: Maui
p5 = ranksum(Saccadelatency_withinNL(1:Saccadelatency_withinNLidx),Saccadelatency_withinL(1:Saccadelatency_withinLidx));
[~,p6] = ttest2(Saccadelatency_withinNL(1:Saccadelatency_withinNLidx),Saccadelatency_withinL(1:Saccadelatency_withinLidx));
p7 = ranksum(Saccadelatency_outsideNL(1:Saccadelatency_outsideNLidx),Saccadelatency_outsideL(1:Saccadelatency_outsideLidx));
[~,p8] = ttest2(Saccadelatency_outsideNL(1:Saccadelatency_outsideNLidx),Saccadelatency_outsideL(1:Saccadelatency_outsideLidx));
p9 = ranksum(Saccadelatency_withinNL(1:Saccadelatency_withinNLidx),Saccadelatency_outsideNL(1:Saccadelatency_outsideNLidx));
[~,p10] = ttest2(Saccadelatency_withinNL(1:Saccadelatency_withinNLidx),Saccadelatency_outsideNL(1:Saccadelatency_outsideNLidx));

% Statistical tests for saccade latencies: Apollo
p11 = ranksum(Saccadelatency_withinNL(Saccadelatency_withinNLidx+1:end),Saccadelatency_withinL(Saccadelatency_withinLidx+1:end));
[~,p12] = ttest2(Saccadelatency_withinNL(Saccadelatency_withinNLidx+1:end),Saccadelatency_withinL(Saccadelatency_withinLidx+1:end));
p13 = ranksum(Saccadelatency_outsideNL(Saccadelatency_outsideNLidx+1:end),Saccadelatency_outsideL(Saccadelatency_outsideLidx+1:end));
[~,p14] = ttest2(Saccadelatency_outsideNL(Saccadelatency_outsideNLidx+1:end),Saccadelatency_outsideL(Saccadelatency_outsideLidx+1:end));
p15 = ranksum(Saccadelatency_withinNL(Saccadelatency_withinNLidx+1:end),Saccadelatency_outsideNL(Saccadelatency_outsideNLidx+1:end));
[~,p16] = ttest2(Saccadelatency_withinNL(Saccadelatency_withinNLidx+1:end),Saccadelatency_outsideNL(Saccadelatency_outsideNLidx+1:end));

% Population analyses on catch trial centroids
p17 = signrank(centroidcatchL,centroidcatchNL);
[~,p18] = ttest(centroidcatchL,centroidcatchNL);

% Statistical tests for catch amplitudes and latencies: Maui
p21 = ranksum(SaccadeAmplitudes_catchL(1:SaccadeAmplitudes_catchLidx),SaccadeAmplitudes_catchNL(1:SaccadeAmplitudes_catchNLidx));
[~,p22] = ttest2(SaccadeAmplitudes_catchL(1:SaccadeAmplitudes_catchLidx),SaccadeAmplitudes_catchNL(1:SaccadeAmplitudes_catchNLidx));
p23 = ranksum(Saccadelatency_catchL(1:Saccadelatency_catchLidx),Saccadelatency_catchNL(1:Saccadelatency_catchNLidx));
[~,p24] = ttest2(Saccadelatency_catchL(1:Saccadelatency_catchLidx),Saccadelatency_catchNL(1:Saccadelatency_catchNLidx));

% Statistical tests for catch amplitudes and latencies: Apollo
p25 = ranksum(SaccadeAmplitudes_catchL(SaccadeAmplitudes_catchLidx+1:end),SaccadeAmplitudes_catchNL(SaccadeAmplitudes_catchNLidx+1:end));
[~,p26] = ttest2(SaccadeAmplitudes_catchL(SaccadeAmplitudes_catchLidx+1:end),SaccadeAmplitudes_catchNL(SaccadeAmplitudes_catchNLidx+1:end));
p27 = ranksum(Saccadelatency_catchL(Saccadelatency_catchLidx+1:end),Saccadelatency_catchNL(Saccadelatency_catchNLidx+1:end));
[~,p28] = ttest2(Saccadelatency_catchL(Saccadelatency_catchLidx+1:end),Saccadelatency_catchNL(Saccadelatency_catchNLidx+1:end));

% Implementing randomization tests for analyzing catch trials
LatencyData = [Saccadelatency_catchL; Saccadelatency_catchNL];
SaccadeData = [SaccadecatchXY_L; SaccadecatchXY_NL];
Laseridxs = [ones(size(SaccadecatchXY_L,1),1); zeros(size(SaccadecatchXY_L,1),1)];
Monkidxs = [ones(SaccadeAmplitudes_catchLidx,1); zeros(-SaccadeAmplitudes_catchLidx+size(SaccadecatchXY_L,1),1); ones(SaccadeAmplitudes_catchNLidx,1); zeros(-SaccadeAmplitudes_catchNLidx+size(SaccadecatchXY_NL,1),1)];
uniqueconds = unique([Monkidxs],'rows','stable');

% plotting catch trial saccade end points
figure(plot_counter); set(gcf,'Name','Catch trial: saccade end points');
subplot(121); plot(SaccadeData(~Laseridxs & Monkidxs,1),SaccadeData(~Laseridxs & Monkidxs,2),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(SaccadeData(Laseridxs & Monkidxs,1),SaccadeData(Laseridxs & Monkidxs,2),'o','MarkerSize',5,'MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]); xlabel('X'); ylabel('Y'); title('Maui');
set(gca,'Xlim',[-15 15],'Ylim',[-15 15],'Tickdir','out'); axis square; grid on; hold off;
subplot(122); plot(SaccadeData(~Laseridxs & ~Monkidxs,1),SaccadeData(~Laseridxs & ~Monkidxs,2),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(SaccadeData(Laseridxs & ~Monkidxs,1),SaccadeData(Laseridxs & ~Monkidxs,2),'o','MarkerSize',5,'MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]); xlabel('X'); ylabel('Y'); title('Apollo')
set(gca,'Xlim',[-15 15],'Ylim',[-15 15],'Tickdir','out'); axis square; grid on; hold off;
plot_counter = plot_counter + 1;

% Randomization tests on catch trials
tmaxidxs = LatencyData==tmax;
SaccadeData(tmaxidxs,:) = [];
Laseridxs(tmaxidxs) = [];
Monkidxs(tmaxidxs) = [];
niter = 10000;
result = nan*ones(size(uniqueconds,1),2);
data = [SaccadeData Laseridxs];
for condcounter = 1:size(uniqueconds,1)
    L1 =  Monkidxs == uniqueconds(condcounter,1);
    
    tmpdata = data(L1,:);
    % Now we have the data from a single condition ? doing the permutation test.
    
    randomizedteststats = zeros(niter,1);
    for i = 0:niter
        if i == 0
            L = logical(tmpdata(:,3));
        else
            L = L(randperm(length(L)));
        end
        squareddist = sum((mean(tmpdata(L,[1 2]))-mean(tmpdata(~L,[1 2]))).^2);
        if i == 0
            result(condcounter,1) = squareddist;
        else
            randomizedteststats(i) = squareddist;
        end
    end
    result(condcounter,2) = sum(randomizedteststats>=result(condcounter,1))/niter;
    result(condcounter,:)
    figure(plot_counter); subplot(1,2,condcounter); histogram(randomizedteststats,20); hold on; plot(result(condcounter,1),0,'kv','MarkerFaceColor',[0 0 0]); axis square; hold off
    
end
plot_counter = plot_counter + 1;
% Checking out the Saccade Amplitudes from the algorithm
check = 0;
if check
    figure(plot_counter); set(gcf,'Name','Saccade Amplitudes: All');
    histogram(SaccadeAmplitudes_all,0:0.25:10,'FaceColor',[0 0 0]); set(gca,'Tickdir','out','Xlim',[0 10],'XTick',0:2:10); xlabel('Amplitude'); ylabel('# saccades'); axis square;
    plot_counter = plot_counter + 1;
end

%% No figure
% �	RF locations for both the monkeys

if ~exist('plot_counter')
    plot_counter = 1;
end
Monkeyidxs = [];
RF = [];
for mm = 1:2
    if mm == 1
        load filenameoptoM.mat
        load RWmultiplieroptoM.mat
        load laserdialoptoM.mat
        monkeyname = 'Maui';
    else
        load filenameoptoA.mat
        load RWmultiplieroptoA.mat
        load laserdialoptoA.mat
        monkeyname = 'Apollo';
    end
    filename = [filenameopto];
    RWmultiplier = [RWmultiplieropto];
    laserdial = [laserdialopto];
    for jj = 1:size(filename,1)
        stro = nex2stro(findfile(char(filename(jj,:))));
        if isnan(unique(stro.trial(:,21:22)))
            stimlocations = [stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y];
            RF = [RF; stimlocations];
            Monkeyidxs = [Monkeyidxs; filename{jj}(1)];
        else
            [stimlocations,~,ib] = unique(stro.trial(:,21:22),'rows');
            RF = [RF; stimlocations];
            Monkeyidxs = [Monkeyidxs; filename{jj}(1)];
        end
    end
end
RF = RF/10;
figure(plot_counter); plot(RF(Monkeyidxs=='M',1),RF(Monkeyidxs=='M',2),'o','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(RF(Monkeyidxs=='A',1),RF(Monkeyidxs=='A',2),'o','MarkerSize',8,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
plot(0,0,'*','color',[0 0 0]); hold on;
set(gca,'Xlim',[-20 20],'Ylim',[-20 20],'Tickdir','out','XTick',-20:5:20,'YTick',-20:5:20); grid on; legend('Maui','Apollo'); axis square;
xlabel('Degrees'), ylabel('Degrees');
plot_counter = plot_counter + 1;


%% Supplementary figure 2
% �	DToneloc: Retinotopic specificity, figure showing effect in 1 spatial location and not in 2 adjacent locations. (Total 4 figures: 1 figure of RF locations of the randomly
% interleaved stimuli, 3 figures of Hits, Misses, CR, FA for 3 tested locations) (psychometric function)
load('T_vos1978_Y');
Vlambda = T_vos1978_Y';
filename = {'M052318001.nex';'M052318002.nex';'M052318003.nex';'M052318004.nex'}; % 5/23, opto session, tested 3 spatial locations
if ~exist('plot_counter')
    plot_counter = 1;
end
HitsL = []; HitsNL = [];
MissL = []; MissNL = [];
CRL = []; CRNL = [];
FAL = []; FANL = [];
locationstested = [];
timebins = 0:0.02:0.3;
Leftsaccadehists = [];
Rightsaccadehists = [];
stimabsentL = []; stimabsentNL = [];
stimpresentL = []; stimpresentNL = [];
alpha = 0.05;
LMScontrastL = []; LMScontrastNL = [];
answersL = []; answersNL = [];
for jj = 1:numel(filename)
    stro = nex2stro(findfile(char(filename(jj,:)))); % actual opto file
    fponidx  = strcmp(stro.sum.trialFields(1,:),'fpon_t');
    fpacqidx  = strcmp(stro.sum.trialFields(1,:),'fpacq_t');
    fpoffidx  = strcmp(stro.sum.trialFields(1,:),'fpoff_t');
    stimonidx  = strcmp(stro.sum.trialFields(1,:),'stimon_t');
    stimoffidx  = strcmp(stro.sum.trialFields(1,:),'stimoff_t');
    targonidx  = strcmp(stro.sum.trialFields(1,:),'targon_t');
    saccstartidx  = strcmp(stro.sum.trialFields(1,:),'saccstart_t');
    saccendidx  = strcmp(stro.sum.trialFields(1,:),'saccend_t');
    rewidx  = strcmp(stro.sum.trialFields(1,:),'rew_t');
    stimpresentidx  = strcmp(stro.sum.trialFields(1,:),'stimpresent');
    Lcc = strcmp(stro.sum.trialFields(1,:),'lcc');
    Mcc = strcmp(stro.sum.trialFields(1,:),'mcc');
    Scc = strcmp(stro.sum.trialFields(1,:),'scc');
    tf = strcmp(stro.sum.trialFields(1,:),'tf');
    oog = strcmp(stro.sum.trialFields(1,:),'oog');
    stimidx = strcmp(stro.sum.trialFields(1,:),'stim_idx');
    optstim = strcmp(stro.sum.trialFields(1,:),'optstim');
    correct = strcmp(stro.sum.trialFields(1,:),'correct');
    laseron = strcmp(stro.sum.trialFields(1,:),'laseron_t');
    laseroff = strcmp(stro.sum.trialFields(1,:),'laseroff_t');
    stim_x = strcmp(stro.sum.trialFields(1,:),'stim_x');
    stim_y = strcmp(stro.sum.trialFields(1,:),'stim_y');
    correcttrials = stro.trial(:,correct);
    stimpresent = logical(stro.trial(:,stimpresentidx));
    LMStriplet = [stro.trial(:,Lcc) stro.trial(:,Mcc) stro.trial(:,Scc)];
    lasertrials = logical(stro.trial(:,optstim));
    stimlocationsall = [stro.trial(:,stim_x) stro.trial(:,stim_y)];
    [stimlocations,~,ib] = unique(stimlocationsall,'rows');
    locationstested = [locationstested; stimlocations];
    for ii = 1:numel(unique(ib))
        ind = logical(ib==ii);
        LMSdirtripletslaser = LMStriplet(stimpresent & lasertrials & ind,:);
        LMSdirtripletsnonlaser = LMStriplet(stimpresent & ~lasertrials & ind,:);
        LMSdircontrastlaser = sqrt(sum(LMSdirtripletslaser.^2,2));
        LMSdircontrastnonlaser = sqrt(sum(LMSdirtripletsnonlaser.^2,2));
        
        oogidxslaser = logical(stro.trial(stimpresent & lasertrials & ind,oog));
        oogidxsnonlaser = logical(stro.trial(stimpresent & ~lasertrials & ind,oog));
        
        lasertrialnumbers = 1:numel(LMSdircontrastlaser);
        nonlasertrialnumbers = 1:numel(LMSdircontrastnonlaser);
        
        lasertrialidxs = logical(stro.trial(ind,optstim));
        colordirchoiceidxs = correcttrials(ind);
        colorstimpresent = logical(stimpresent(ind));
        percentcorrectlasertrials = sum(colordirchoiceidxs & lasertrialidxs)/sum(lasertrialidxs);
        percentcorrectnonlasertrials = sum(colordirchoiceidxs & ~lasertrialidxs)/sum(~lasertrialidxs);
        
        % Laser trials
        Hitlasertrial = sum(colordirchoiceidxs & lasertrialidxs & colorstimpresent); % Hit
        Misslasertrial = sum(~colordirchoiceidxs & lasertrialidxs & colorstimpresent); % Miss
        CRlasertrial = sum(colordirchoiceidxs & lasertrialidxs & ~colorstimpresent); % Correct Reject
        FAlasertrial = sum(~colordirchoiceidxs & lasertrialidxs & ~colorstimpresent); % False Alarm
        HitsL = [HitsL; Hitlasertrial];
        MissL = [MissL; Misslasertrial];
        CRL = [CRL; CRlasertrial];
        FAL = [FAL; FAlasertrial];
        
        % Non-Laser trials
        Hitnonlasertrial = sum(colordirchoiceidxs & ~lasertrialidxs & colorstimpresent); % Hit
        Missnonlasertrial = sum(~colordirchoiceidxs & ~lasertrialidxs & colorstimpresent); % Miss
        CRnonlasertrial = sum(colordirchoiceidxs & ~lasertrialidxs & ~colorstimpresent); % Correct Reject
        FAnonlasertrial = sum(~colordirchoiceidxs & ~lasertrialidxs & ~colorstimpresent); % False Alarm
        HitsNL = [HitsNL; Hitnonlasertrial];
        MissNL = [MissNL; Missnonlasertrial];
        CRNL = [CRNL; CRnonlasertrial];
        FANL = [FANL; FAnonlasertrial];
        stimabsentL = [stimabsentL; sum(CRlasertrial) + sum(FAlasertrial)];
        stimabsentNL = [stimabsentNL; sum(CRnonlasertrial) + sum(FAnonlasertrial)];
        stimpresentL = [stimpresentL; sum(Hitlasertrial) + sum(Misslasertrial)];
        stimpresentNL = [stimpresentNL; sum(Hitnonlasertrial) + sum(Missnonlasertrial)];
        
        fundamentals = stro.sum.exptParams.fundamentals; % CONE FUNDAMENTALS: L,M,S
        fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
        mon_spd = stro.sum.exptParams.mon_spd; % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
        mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
        mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
        M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals
        Lb = M*stro.sum.exptParams.bkgndrgb;
        % Converting all the contrasts to Luminance contrast
        Bkgnd_luminance = stro.sum.exptParams.bkgndrgb'*mon_spd*Vlambda;
        Peak_luminance = (inv(M)*((1+LMStriplet).*(repmat(Lb',[size(LMStriplet,1) 1])))')'*mon_spd*Vlambda;
        Trough_luminance = (inv(M)*((1-LMStriplet).*(repmat(Lb',[size(LMStriplet,1) 1])))')'*mon_spd*Vlambda;
        Contrast_luminance = (Peak_luminance - Trough_luminance)./(Peak_luminance + Trough_luminance);
        
        
        LMScontrastfileNL = Contrast_luminance(colorstimpresent & ~lasertrialidxs,:);
        LMScontrastfileL = Contrast_luminance(colorstimpresent & lasertrialidxs,:);
        % calculating gamut edge luminance contrast
        t = min(((1./stro.sum.exptParams.bkgndrgb)-1));
        OOG_luminance = (stro.sum.exptParams.bkgndrgb*(1+t))' * mon_spd*Vlambda;
        putativeoogcontrast =  (OOG_luminance - Bkgnd_luminance)/Bkgnd_luminance;
        
        trials = lasertrials;
        for kk = 1:2
            if kk == 2
                trials = ~trials;
            end
            LMScontrast = Contrast_luminance(stimpresent & trials & ind,:);
            answers = correcttrials(stimpresent & trials & ind);
            LMScontrast(LMScontrast>putativeoogcontrast) = putativeoogcontrast;
            if kk == 1
                LMScontrastL = [LMScontrastL; LMScontrast'];
                answersL = [answersL; answers'];
            else
                LMScontrastNL = [LMScontrastNL; LMScontrast'];
                answersNL = [answersNL; answers'];
            end
        end
    end
end
% Averaging out the data multiple files
[locationstested,ia,ib] = unique(locationstested,'rows');
nrows = numel(unique(ib));

% Plotting the RF location wrt to fixation point
cgradations = linspace(0,1,size(locationstested,1));
figure(plot_counter), set(gcf,'Name','Spatial extent: randomly interleaved, multiple locations');
subplot(331); plot(stro.sum.exptParams.fp_x,stro.sum.exptParams.fp_y,'+','Markersize',15,'Linewidth',2);
for ii = 1:numel(unique(ib))
    subplot(331);hold on; plot(locationstested(ii,1), locationstested(ii,2),'o','color',[1-cgradations(ii) 0 cgradations(ii)],'Markersize',15,'Linewidth',2);
    ind = find(ib==ii);
    [~,p1] = equalproptest([sum(HitsL(ind)) sum(HitsNL(ind))],[sum(stimpresentL(ind)) sum(stimpresentNL(ind))],alpha);
    [~,p2] = equalproptest([sum(FAL(ind)) sum(FANL(ind))],[sum(stimabsentL(ind)) sum(stimabsentNL(ind))],alpha);
    subplot(3,3,ii+1); h = bar([sum(HitsL(ind))/(sum(HitsL(ind))+sum(MissL(ind))) sum(HitsNL(ind))/(sum(HitsNL(ind))+sum(MissNL(ind))); sum(MissL(ind))/(sum(HitsL(ind))+sum(MissL(ind))) sum(MissNL(ind))/(sum(HitsNL(ind))+sum(MissNL(ind))); sum(CRL(ind))/(sum(CRL(ind))+sum(FAL(ind))) sum(CRNL(ind))/(sum(CRNL(ind))+sum(FANL(ind))); sum(FAL(ind))/(sum(CRL(ind))+sum(FAL(ind))) sum(FANL(ind))/(sum(CRNL(ind))+sum(FANL(ind)))]); hold on;
    set(h(2),'FaceColor',[0 0.5 1]); set(h(1),'FaceColor',[0.5 0.5 0.5]);
    text(1,max([sum(HitsL(ind));sum(HitsNL(ind))])+20,strcat('p1=',num2str(p1,3))); text(3.5,max([sum(FAL(ind));sum(FANL(ind))])+20,strcat('p2=',num2str(p2,3)));
    ylabel('Prop. of trials'); title(num2str(locationstested(ii,:))); set(gca,'XTick',[1 2 3 4],'Xlim',[0 5],'YTick',[0:0.25:1],'XTickLabel',{'H','M','CR','FA'},'TickDir','Out','XColor',[1-cgradations(ii) 0 cgradations(ii)],'YColor',[1-cgradations(ii) 0 cgradations(ii)]); drawnow;
    
end
subplot(331), hold on; axis square; set(gca,'Xlim',[-100 100],'Ylim',[-100 100],'TickDir','Out'); grid on; xlabel('X'), ylabel('Y'); hold off;
subplot(332); hold on; axis square; hold off;
subplot(333); hold on; axis square; hold off;
subplot(334); hold on; axis square; hold off;

% Need to calculate a psychometric function for all the locations
for ii = 1:numel(ia)
    idx = ia(ii);
    tmplasercontrast = reshape(LMScontrastL(ib==idx,:),[],1); % converting array to a vector
    tmpcontrolcontrast = reshape(LMScontrastNL(ib==idx,:),[],1);
    answerlaser = reshape(answersL(ib==idx,:),[],1);
    answercontrol = reshape(answersNL(ib==idx,:),[],1);
    correctanswersL = []; incorrectanswersL = []; contrastL = [];
    correctanswersNL = []; incorrectanswersNL = []; contrastNL = [];
    [val,ind] = unique(tmplasercontrast);
    for jj = 1:numel(val)
        contrastL = [contrastL; val(jj)];
        correctanswersL = [correctanswersL; sum(answerlaser(tmplasercontrast==val(jj)))];
        incorrectanswersL = [incorrectanswersL; numel(answerlaser(tmplasercontrast==val(jj))) - sum(answerlaser(tmplasercontrast==val(jj)))];
    end
    [aL,bL,gL] = weibullFitforDToneloc(contrastL,[correctanswersL incorrectanswersL],'mle');
    [val,ind] = unique(tmpcontrolcontrast);
    for jj = 1:numel(val)
        contrastNL = [contrastNL; val(jj)];
        correctanswersNL = [correctanswersNL; sum(answercontrol(tmpcontrolcontrast==val(jj)))];
        incorrectanswersNL = [incorrectanswersNL; numel(answercontrol(tmpcontrolcontrast==val(jj))) - sum(answerlaser(tmpcontrolcontrast==val(jj)))];
    end
    [aNL,bNL,gNL] = weibullFitforDToneloc(contrastNL,[correctanswersNL incorrectanswersNL],'mle');
    contrastlattice = logspace(log10(0.03),log10(1.0),51);
    fitL = gL*(1-exp(-((contrastlattice./aL).^bL)));
    fitNL = gNL*(1-exp(-((contrastlattice./aNL).^bNL)));
    figure(plot_counter); subplot(3,3,ii+4); plot(contrastL,correctanswersL./(correctanswersL+incorrectanswersL),'o','MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]); hold on;
    plot(contrastNL,correctanswersNL./(correctanswersNL+incorrectanswersNL),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
    plot(contrastlattice,fitL,'-','Linewidth',2,'color',[0 0.5 1.0]); plot(contrastlattice,fitNL,'-','Linewidth',2,'color',[0.5 0.5 0.5]);
    xlabel('contrast'); ylabel('Prop. correct'); set(gca,'Xlim',[0.03 1.0],'Tickdir','out','YTick',[0:0.25:1],'Xscale','log','XTick',[0.1 0.3 1],'XTickLabels', {'0.1','0.3','1'});
    title(locationstested(ii,:)'); axis square; hold off;
end

plot_counter = plot_counter + 1;


%% Supplementary figure 4
% example of opto-tagging, direction selective cell

if ~exist('plot_counter')
    plot_counter = 1;
end
filename = ['M032518015.nex'];
% Another example of FixStim interfreq is 'A011419005.nex'
for jj = 1:size(filename,1)
    stro = nex2stro(findfile(filename(jj,:)));
    offset = [-.1 .1];  % pre and post time wrt stimon/stimoff
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
    fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
    fpacq_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpacq_t'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
    targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targon_t'));
    targoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targoff_t'));
    opt_stimfreq = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'opt_stimfreq'));
    targ_shown = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_shown'));
    optfreq = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'elec_stimfreq'));
    uniqfreqs = unique(optfreq);
    stimon_t(isnan(stimon_t)) = fpacq_t(isnan(stimon_t))+nanmean(stimon_t-fpacq_t); %
    modulator = stro.sum.exptParams.modulator;
    
    % Hack below
    if (isnan(targ_shown(1)))
        targ_shown(:,1) = 0;
    end
    dur = mode(stimoff_t-stimon_t);
    sync_t = stimon_t;  % Ecode for alignment
    Lspikechans = strncmp(stro.sum.rasterCells,'sig0',4);
    
    figure(plot_counter); set(gcf,'Name','FixStim interfreq');
    annotation('textbox',[.4 .89 .2 .1],'string',['CELL ',filename(jj,:)],'HorizontalAlignment','Center','FitBoxToText','on','EdgeColor','none')
    N = uniqfreqs(2:9)';
    if (any(Lspikechans))
        for whichspike = find(Lspikechans)
            spikes = stro.ras(:,whichspike);
            binwidth = .005; % s
            bins = offset(1):binwidth:dur+offset(2);
            PSTH = zeros(1,length(bins));
            for j = N
                PSTH = zeros(1,length(bins));
                %figure; subplot(2,1,1); hold on;
                subplot(length(N), 2, 2*(find(N==j))-1); hold on;
                L = optfreq == j;
                trlidxs = find(L);
                for counter = 1:sum(L)
                    trlidx = trlidxs(counter);
                    tmpspikes = spikes{trlidx}-sync_t(trlidx);
                    tmpspikes(tmpspikes < offset(1) | tmpspikes > dur+offset(2)) = [];
                    nspikestot = length(tmpspikes);
                    plot([tmpspikes tmpspikes]',[zeros(nspikestot,1) .5*ones(nspikestot,1)]'+counter,'k-');
                    PSTH = PSTH + hist(tmpspikes, bins);
                end
                PSTH = PSTH./(sum(L).*binwidth);
                if (j > 0) % Plotting the time course of optical stimulation
                    if modulator
                        f = 0.5*1000^(1/255)^j;
                        t = linspace(bins(1),bins(end),5e4);
                        y = sin(2*pi*f*t-pi/2);
                        y(t<0) = -1;
                        y(t>dur) = -1;
                        y = y-.5; % for plotting below the spikes
                        plot(t,y,'-','color',[0 0.5 1],'linewidth',2);
                    else
                        secspercycle = 1/unique(j(j > 0));
                        transitions = 0:secspercycle/2:dur;
                        if (ceil(length(transitions)/2) ~= floor(length(transitions)/2))
                            transitions(end+1) = dur; %automatic shutoff
                        end
                        x = [transitions; transitions];
                        x = [offset(1); x(:); max(x(:))+offset(2)];
                        y = [repmat([0 1],1,length(transitions)/2) 0]*2-2;
                        y = [y;y];
                        if (length(x(:)) == length(y(:)))
                            plot(x,y(:)','k-','linewidth',2);
                        end
                    end
                end
                set(gca,'XLim', [-0.1 0.6],'XTick',[-0.1:0.1:0.6],'Ytick',[],'YLim',[-4 sum(L)+5],'Tickdir','out');
                if modulator && j > 0
                    title(['Frequency: ',num2str(round(f)),' Hz']);
                else
                    title(['Frequency: ',num2str(j)]);
                end
                
                % PSTH
                subplot(length(N), 2, 2*(find(N==j))-0); hold on;
                plot(bins,PSTH,'k-','LineWidth',2);
                set(gca,'YLim',[0 600],'YTick',[0:300:600]);
                set(gca,'Xlim',[-0.1 0.6],'Tickdir','out','XTick',[-0.1:0.1:0.6]);
                xlabel('Time (s)','FontSize',12);
                ylabel('Response (sp/s)','FontSize',12);
            end
        end
        set(gcf,'renderer','painters');
        plot_counter = plot_counter + 1;
    end
end

% Plotting the waveform
strowf = nex2stro(findfile('M032518015_wf.nex'));
samplingrate = strowf.sum.waves.storeRates{1};
leastcount = 10^6/samplingrate; % in us
waveform = mean(cell2mat(strowf.ras(:,3))',2);
noise_waveform = mean(cell2mat(strowf.ras(:,4))',2);
time = 0:25:25*(numel(waveform)-1);
[val1,t1] = max(waveform);
[val2,t2] = min(waveform);
figure(plot_counter); plot(time,waveform,'-','color',[0 0 0],'Linewidth',2); hold on; plot(time,noise_waveform,'-','color',[0.5 0.5 0.5],'Linewidth',2);
plot(time,waveform-std(cell2mat(strowf.ras(:,3))',0,2),'-','color',[0 0 0],'Linewidth',1); plot(time,waveform+std(cell2mat(strowf.ras(:,3))',0,2),'-','color',[0 0 0],'Linewidth',1);
plot(time,noise_waveform-std(cell2mat(strowf.ras(:,4))',0,2),'-','color',[0.5 0.5 0.5],'Linewidth',1); plot(time,noise_waveform+std(cell2mat(strowf.ras(:,4))',0,2),'-','color',[0.5 0.5 0.5],'Linewidth',1);
line([0 max(time)],[0 0]); set(gca,'Xlim',[min(time) max(time)],'XTick',0:200:800,'Tickdir','out'); axis square; title('Waveform: single unit'); hold off;
plot_counter = plot_counter + 1;


% Gratings file
filename = ['M032518016.nex'];
GT=nex2stro(findfile(filename)); % Gratings file
orients = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'orient'));
sfs = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'sf'));
tfs = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'tf'));
diams = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'diam'));
protocols = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'protocol'));
stimtype = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'stimtype'));
framerate = GT.sum.exptParams.framerate;
nframes = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'nframes'));
stimon_t = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'stim_on'));
stimoff_t= GT.trial(:,strcmp(GT.sum.trialFields(1,:),'stim_off'));
spikeidxs = find(strncmp(GT.sum.rasterCells,'sig',3));
spikenames = GT.sum.rasterCells(spikeidxs);

figure(plot_counter); set(gcf,'Name',strcat('Gratings:',filename));
for spikeidx = spikeidxs
    Lcc = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'lcont'));
    Mcc = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'mcont'));
    Scc = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'scont'));
    colordirections = unique([Lcc Mcc Scc],'rows');
    
    spikerates = [];
    baselines = [];  baseline_t = 0.5;
    for i = 1:size(GT.trial,1)
        spiketimes = GT.ras{i,spikeidx};
        nspikes = sum(spiketimes > stimon_t(i) & spiketimes < stimoff_t(i));
        spikerates = [spikerates; nspikes./(stimoff_t(i)-stimon_t(i))];
        nspikes = sum(spiketimes > stimon_t(i)-baseline_t & spiketimes < stimon_t(i));
        baselines = [baselines; nspikes./baseline_t];
    end
    
    if any(protocols == 1)
        % Looking at orientation tuning curve(s)
        orienttrials = find(protocols == 1);
        protocolswitches = orienttrials(find(diff(orienttrials) > 1)+1);
        starts = [find(orienttrials == 1,1,'first'); protocolswitches];
        stops = [];
        for i = 1:length(starts)
            if (i == length(starts))
                stops = [stops; max(orienttrials(orienttrials > starts(i)))];
            else
                stops = [stops; max(orienttrials(orienttrials > starts(i) & orienttrials < starts(i+1)))];
            end
        end
        
        for i = 1:length(starts)
            trlidxs = [starts(i):stops(i)];
            x = orients(trlidxs);
            y = spikerates(trlidxs);
            Ltmp = x == min(x);
            y = [y; y(Ltmp)];
            x = [x; x(Ltmp)+2*pi];
            sf = unique(sfs(trlidxs));  % For axis label
            if (length(sf > 1))
                disp('Error: Parameters changed during protocol 1 but "reset" not clicked');
            end
            mu = []; sem = [];
            for j = unique(x)'
                mu(j == unique(x))= mean(y(x==j));
                sem(j == unique(x)) = std(y(x==j))/sqrt(sum(x==j));
            end
        end
        
        
        %         subplot(121); polarplot(unique(x)',mu+sem,'-','Color',[0 0 0]); hold on;
        %         polarplot(unique(x)',mu,'-','Color',[0 0 0],'LineWidth',2); polarplot(unique(x)',mu-sem,'-','Color',[0 0 0]);
        %         polarplot(linspace(0,2*pi,100),repmat(mean(baselines),1,100)); thetaticks(0:45:315); rticks([0 30 60]); title('Orientation tuning');
        
        subplot(121); polarplot(x,y,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
        polarplot(unique(x)',mu,'-','Color',[0 0 0],'LineWidth',2);
        thetaticks(0:45:315); rticks([0 30 60]); title('Orientation tuning');
        
        % Looking at spatial frequency tuning curve(s)
        sftrials = find(protocols == 2);
        protocolswitches = sftrials(find(diff(sftrials) > 1)+1);
        starts = [min(sftrials); protocolswitches];
        stops = [];
        for i = 1:length(starts)
            if (i == length(starts))
                stops = [stops; max(sftrials(sftrials > starts(i)))];
            else
                stops = [stops; max(sftrials(sftrials > starts(i) & sftrials < starts(i+1)))];
            end
        end
        
        for i = 1:length(starts)
            trlidxs = [starts(i):stops(i)];
            x = sfs(trlidxs);
            y = spikerates(trlidxs);
            pp = csape(x,y,'variational');
            xx = linspace(min(sfs),max(sfs),100);
            fit = ppval(pp,xx);
            subplot(122); hold on; plot(x,y,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
            %             plot(xx,fit,'k','Linewidth',2);
            orient = unique(orients(trlidxs)); % for axis labeling
            title(['orientation: ',num2str(orient*180/pi,3),' deg']);
            xlabel('spatial frequency (cyc/deg)');
            ylabel('response (sp/sec)');
            mu = []; sem = [];
            for j = unique(x)'
                mu(j == unique(x))= mean(y(x==j));
                sem(j == unique(x)) = std(y(x==j))/sqrt(sum(x==j));
            end
        end
        %         errorbar(unique(x),mu,sem,'k');
        plot(unique(x)',mu,'k','Linewidth',2);
        %         plot([min(x) max(x)], repmat(mean(baselines),1,2),'k:');
        %         axis tight;
        set(gca,'Tickdir','out','YLim',[0 60],'XLim',[0.3 10],'XTick',[0.3 1 3 10],'XScale','log'); axis square; hold off;
    end
end
plot_counter = plot_counter + 1;


%% No Figure: Waveform analyses
% �	Analysis of single-unit waveforms, correlating suppressed vs activated single units with peak/trough ratio (1 figure: suppressive unit waveform, 1 figure:
% activated unit waveform, 1 figure: histogram of peak/trough ratio segregated by suppressed and activated unit) (Check for the latency of the outlier, suppressed cell, rasters)

if ~exist('plot_counter')
    plot_counter = 1;
end
filename_M = {'M051318008.nex';'M051318009.nex';'M051618004.nex';'M051618005.nex';'M051618006.nex';'M051618007.nex';'M051618008.nex';'M051618010.nex';...
    'M051618011.nex';'M051618012.nex';'M051618013.nex';'M051618014.nex';'M051818004.nex';'M051818005.nex';'M051818006.nex';'M051818007.nex';...
    'M051818010.nex';'M051818011.nex';'M051818012.nex';'M051818013.nex';'M051818014.nex';'M052318001.nex';'M052318002.nex';'M052318003.nex';'M052318004.nex'};
filename_A = {'A011719009.nex';'A011719010.nex';'A011719011.nex';'A011719012.nex';'A012219006.nex';'A012219009.nex';'A012219010.nex';'A012219011.nex';'A012219012.nex';...
    'A012319004.nex';'A012319005.nex';'A012319006.nex';'A013119002.nex';'A013119003.nex';'A013119004.nex';'A013119005.nex';'A013119006.nex';'A020119004.nex';'A020119005.nex';...
    'A020319001.nex';'A020319007.nex';'A020319008.nex';'A020319009.nex';'A020519003.nex';'A020519008.nex';'A020519009.nex';'A021519004.nex';'A021519005.nex';'A021519006.nex';...
    'A021519007.nex';'A021519008.nex';'A021519010.nex'};
filename = [filename_M; filename_A];
Singlevsmultiunit_M = ['M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M'];
Singlevsmultiunit_A = ['M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'S';'S';'S';'M';'S';'S';'M';'S';'M';'M';'M';'S';'M';'M';'M';'M';'M';'S';'S';'S';'S';'M'];
N = numel(filename);
L = ceil(sqrt(N));
binwidth = .005;
bins = -0.4:binwidth:0.6;
baselinebins = (bins>0 & bins<0.3);
baselineFR = [];
laserFR = [];
statstestresult = []; % Performing a wilcoxon rank sum test and storing the p value
ISI = [];
for jj = 1:N
    stro = nex2stro(findfile(char(filename(jj,:))));
    fponidx  = strcmp(stro.sum.trialFields(1,:),'fpon_t');
    fpacqidx  = strcmp(stro.sum.trialFields(1,:),'fpacq_t');
    fpoffidx  = strcmp(stro.sum.trialFields(1,:),'fpoff_t');
    stimonidx  = strcmp(stro.sum.trialFields(1,:),'stimon_t');
    stimoffidx  = strcmp(stro.sum.trialFields(1,:),'stimoff_t');
    targonidx  = strcmp(stro.sum.trialFields(1,:),'targon_t');
    saccstartidx  = strcmp(stro.sum.trialFields(1,:),'saccstart_t');
    saccendidx  = strcmp(stro.sum.trialFields(1,:),'saccend_t');
    rewidx  = strcmp(stro.sum.trialFields(1,:),'rew_t');
    stimpresentidx  = strcmp(stro.sum.trialFields(1,:),'stimpresent');
    Lcc = strcmp(stro.sum.trialFields(1,:),'lcc');
    Mcc = strcmp(stro.sum.trialFields(1,:),'mcc');
    Scc = strcmp(stro.sum.trialFields(1,:),'scc');
    tf = strcmp(stro.sum.trialFields(1,:),'tf');
    oog = strcmp(stro.sum.trialFields(1,:),'oog');
    stimidx = strcmp(stro.sum.trialFields(1,:),'stim_idx');
    optstim = strcmp(stro.sum.trialFields(1,:),'optstim');
    correct = strcmp(stro.sum.trialFields(1,:),'correct');
    laseron = strcmp(stro.sum.trialFields(1,:),'laseron_t');
    laseroff = strcmp(stro.sum.trialFields(1,:),'laseroff_t');
    analogstrtimeind = strcmp(stro.sum.rasterCells(1,:),'anlgStartTime');
    spikeind = strcmp(stro.sum.rasterCells(1,:),'sig001a');
    heyetraceind = strcmp(stro.sum.rasterCells(1,:),'AD11');
    veyetraceind = strcmp(stro.sum.rasterCells(1,:),'AD12');
    lasertraceind = strcmp(stro.sum.rasterCells(1,:),'AD13');
    correcttrials = stro.trial(:,correct);
    diridxs = stro.trial(:,stimidx);
    stimpresent = logical(stro.trial(:,stimpresentidx));
    LMStriplet = [stro.trial(:,Lcc) stro.trial(:,Mcc) stro.trial(:,Scc)];
    lasertrials = logical(stro.trial(:,optstim));
    samplerate = stro.sum.analog.storeRates{3}; % sample rate at which laser analog pulses are stored in file
    
    % Stimulus absent trials
    idxs = find(~stimpresent);
    count1 = 1;
    count2 = 1;
    PSTHlaser = zeros(1,length(bins));
    PSTHbaseline = zeros(1,length(bins));
    PSTHnolaser = zeros(1,length(bins));
    tmp_baselineFR = [];
    tmp_laserFR = [];
    tmp_spikecountsbaseline = [];
    tmp_spikecountslaserstimabsent = [];
    
    for ii = 1:numel(idxs)
        ind = idxs(ii);
        analogstartime = stro.ras{ind,analogstrtimeind};
        spiketimes = stro.ras{ind,spikeind};
        if lasertrials(ind)
            % Pulling out the laseron and off timings from analog traces
            t = stro.ras{ind,analogstrtimeind}+[0:1:length(stro.ras{ind,lasertraceind})-1]/samplerate;
            laserontime = t(find(stro.ras{ind,lasertraceind}>0.1,1));
            laserofftime = t(find(stro.ras{ind,lasertraceind}>0.1,1,'last'));
            timedurlaser = laserofftime - laserontime;
            PSTHlaser = PSTHlaser + hist(spiketimes-laserontime, bins);
            count1 = count1 + 1;
            laserbins = bins>=0 & bins<=timedurlaser;
            tmp_laserFR = [tmp_laserFR; mean(PSTHlaser(laserbins))/timedurlaser];
            ISI = [ISI; diff(spiketimes(spiketimes>laserontime & spiketimes<laserofftime))];
            tmp_spikecountsbaseline = [tmp_spikecountsbaseline; sum(spiketimes>laserontime-0.3 & spiketimes<laserontime)];
            tmp_spikecountslaserstimabsent = [tmp_spikecountslaserstimabsent; sum(spiketimes>laserontime & spiketimes<laserofftime)];
        else
            if ~stimpresent(ind)
                stimontime = stro.trial(ind,stimonidx);
                PSTHbaseline = PSTHbaseline + hist(spiketimes-stimontime, bins);
                tmp_baselineFR = [tmp_baselineFR; mean(PSTHbaseline(baselinebins))/0.3];
            end
        end
    end
    baselineFR = [baselineFR; mean(tmp_baselineFR)];
    laserFR = [laserFR; mean(tmp_laserFR)];
    statstestresult = [statstestresult; ranksum(tmp_spikecountsbaseline,tmp_spikecountslaserstimabsent)];
end

% Including files from FixStim contrast response
filename2 = {'M010718005.nex';'M010718017.nex';'M010818019.nex';'M010818024.nex'};
Singlevsmultiunit_2 = ['M';'S';'M';'S'];
for aa = 1:size(filename2,1)
    stro = nex2stro(findfile(char(filename2(aa,:))));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
    targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targon_t'));
    targoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targoff_t'));
    targ_shown = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_shown'));
    laser_power = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'laser_power'));
    contrasts = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'contrast'));
    fpacq_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpacq_t'));
    analogstrtimeind = strcmp(stro.sum.rasterCells(1,:),'anlgStartTime');
    spikeind = strcmp(stro.sum.rasterCells(1,:),'sig001a');
    heyetraceind = strcmp(stro.sum.rasterCells(1,:),'AD11');
    veyetraceind = strcmp(stro.sum.rasterCells(1,:),'AD12');
    lasertraceind = strcmp(stro.sum.rasterCells(1,:),'AD13');
    uniquecontrasts = unique(contrasts);
    uniquepowers = unique(laser_power);
    sync_t = targon_t;  % Ecode for alignment
    Lspikechans = strncmp(stro.sum.rasterCells,'sig0',4);
    % Plotting the laser Trials
    samplingrate = stro.sum.analog.storeRates{3};
    multiplicationfactor = 0.75;
    spikecounts = cell(1,numel(uniquepowers));
    for jj = 2:numel(uniquepowers)
        L = find(laser_power==uniquepowers(jj));
        laserinterval = [];
        PSTH = zeros(1,length(bins));
        tmp_spikecountsbaseline = [];
        tmp_spikecountslaserstimabsent = [];
        for ii = 1:numel(L)
            ind = L(ii);
            analogstarttime = stro.ras{ind,analogstrtimeind};
            spiketimes = stro.ras{ind,spikeind};
            laserinterval = [laserinterval; stimoff_t(ind)-stimon_t(ind)];
            PSTH = PSTH + hist(spiketimes, bins);
            tmp_spikecountsbaseline = [tmp_spikecountsbaseline; sum(spiketimes>stimon_t(ind)-laserinterval(end) & spiketimes<stimon_t(ind))];
            tmp_spikecountslaserstimabsent = [tmp_spikecountslaserstimabsent; sum(spiketimes>stimon_t(ind) & spiketimes<stimoff_t(ind))];
        end
        baselineFR = [baselineFR; mean(tmp_spikecountsbaseline)/mean(laserinterval)];
        laserFR = [laserFR; mean(tmp_spikecountslaserstimabsent)/mean(laserinterval) ];
        statstestresult = [statstestresult; ranksum(tmp_spikecountsbaseline,tmp_spikecountslaserstimabsent)];
    end
end
firstletters_filename = char([filename;filename2]);
Singlevsmultiunit = [Singlevsmultiunit_M; Singlevsmultiunit_A; Singlevsmultiunit_2];
Singleunitidxs = Singlevsmultiunit == 'S';
Maui_idxs = firstletters_filename(:,1)=='M';
Apollo_idxs = firstletters_filename(:,1)=='A';

% Additional analysis of waveforms
Suppressionvsactivation = zeros(size(baselineFR));
Suppressionvsactivation(baselineFR>laserFR) = 0; % 'Supression'
Suppressionvsactivation(baselineFR<laserFR) = 1; % 'Activation'
idx = find([Singleunitidxs]);
N = ceil(sqrt(numel(idx)));
figure(plot_counter); set(gcf,'Name','Plotting waveforms of isolated units');
peak = []; trough = [];
timediff = [];
newfilename = [filename;filename2];
waveform_all = [];
for ii = 1:numel(idx)
    ind = idx(ii);
    stro = nex2stro(findfile(char(newfilename(ind,:))));
    
    if ~Suppressionvsactivation(ind)
        color = [0.5 0.5 0.5]; % Suppression
    else
        color = [0 0 0]; % Activation
    end
    samplingrate = 40000; %stro.sum.waves.storeRates{1};
    leastcount = 10^6/samplingrate; % in us
    waveform = mean(cell2mat(stro.ras(:,2))',2);
    time1 = 0:25:25*(numel(waveform)-1);
    time = 0:2.5:25*(numel(waveform)-1);
    waveform = spline(time1,waveform,time);
    [val1,t1] = max(waveform); % peak
    [val2,t2] = min(waveform); % trough
    if t1<t2
        tmp = t1;
        t1 = t2;
        t2 = tmp;
        %         waveform = -1*waveform;
    end
    peak = [peak; val1];
    trough = [trough; val2];
    timediff = [timediff; abs(time(t2)-time(t1))];
    figure(plot_counter); subplot(N,N,ii); plot(waveform,'-','color',color,'Linewidth',2); hold on;
    line([t1 t1],[0 val1]); line([t2 t2],[0 val2]);line([0 numel(waveform)],[0 0]); hold off;
    
    waveform_all = [waveform_all; waveform];
end
plot_counter = plot_counter + 1;
% Plotting population stats
peaktotroughratio = peak./abs(trough);
p = ranksum(peaktotroughratio(logical(Suppressionvsactivation(idx))),peaktotroughratio(~Suppressionvsactivation(idx)));
figure(plot_counter); subplot(221); plot(timediff(logical(Suppressionvsactivation(idx))),peaktotroughratio(logical(Suppressionvsactivation(idx))),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(timediff(~(Suppressionvsactivation(idx))),peaktotroughratio(~(Suppressionvsactivation(idx))),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); ylabel('Peak/trough ratio'); xlabel('Spike width (us)');
legend('Activation','Suppression'); axis square; hold off;
subplot(222); plot(peak(logical(Suppressionvsactivation(idx))),abs(trough(logical(Suppressionvsactivation(idx)))),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(peak(~(Suppressionvsactivation(idx))),abs(trough(~(Suppressionvsactivation(idx)))),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); xlabel('Peak amp'); ylabel('trough amp');
legend('Activation','Suppression'); set(gca,'Xlim',[0.02 0.15],'Ylim',[0.02 0.15],'Tickdir','out'); axis square; hold off;
subplot(223); plot(waveform_all(logical(Suppressionvsactivation(idx)),:)','-','color',[0 0 0]); hold on; plot(waveform_all(~logical(Suppressionvsactivation(idx)),:)','-','color',[0.5 0.5 0.5]);hold off;
plot_counter = plot_counter + 1;

%% Supplementary figure 3
% �	DToneloc population: laser control behavior files obtained during behavior, Maui
if ~exist('plot_counter')
    plot_counter = 1;
end
load('T_vos1978_Y');
Vlambda = T_vos1978_Y';
ratio_controloverlaser = [];
for mm = 1:1
    if mm == 1
        load filenamelcM.mat
        load RWmultiplierlcM.mat
        load laserdiallcM.mat
        monkeyname = 'Maui';
    else
        
    end
    filename = [filenamelc];
    RWmultiplier = [RWmultiplierlc];
    laserdial = [laserdiallc];
    HitL = []; HitNL = [];
    MissL = []; MissNL = [];
    CRL = []; CRNL = [];
    numdivisionsperfile = 3;
    CRearlylateL = []; % for storing the number of CRs in first half and second half of the block of experiment (laser trials)
    CRearlylateNL = [];% for storing the number of CRs in first half and second half of the block of experiment (no-laser trials)
    FAL = []; FANL = [];
    stimdur = [];
    laserdur = [];
    h_hit = []; p_hit = [];
    h_FA = []; p_FA = [];
    TPRNL = []; TPRL = [];
    FPRNL = []; FPRL = [];
    alpha = 0.05;
    RF = [];
    stimabsentL = []; stimabsentNL = [];
    stimpresentL = []; stimpresentNL = [];
    newfilename = [];
    newRWmultiplier = [];
    newstimsize = [];
    newlaserdial = [];
    weibullparamsL = [];
    weibullparamsNL = [];
    prop_correctatGAMUTL = [];
    prop_correctatGAMUTNL = [];
    DetectThreshL = [];
    oogvals = [];
    count = 1;
    for jj = 1:size(filename,1)
        stro = nex2stro(findfile(char(filename(jj,:))));
        fponidx  = strcmp(stro.sum.trialFields(1,:),'fpon_t');
        fpacqidx  = strcmp(stro.sum.trialFields(1,:),'fpacq_t');
        fpoffidx  = strcmp(stro.sum.trialFields(1,:),'fpoff_t');
        stimonidx  = strcmp(stro.sum.trialFields(1,:),'stimon_t');
        stimoffidx  = strcmp(stro.sum.trialFields(1,:),'stimoff_t');
        targonidx  = strcmp(stro.sum.trialFields(1,:),'targon_t');
        saccstartidx  = strcmp(stro.sum.trialFields(1,:),'saccstart_t');
        saccendidx  = strcmp(stro.sum.trialFields(1,:),'saccend_t');
        rewidx  = strcmp(stro.sum.trialFields(1,:),'rew_t');
        stimpresentidx  = strcmp(stro.sum.trialFields(1,:),'stimpresent');
        Lcc = strcmp(stro.sum.trialFields(1,:),'lcc');
        Mcc = strcmp(stro.sum.trialFields(1,:),'mcc');
        Scc = strcmp(stro.sum.trialFields(1,:),'scc');
        tf = strcmp(stro.sum.trialFields(1,:),'tf');
        oog = strcmp(stro.sum.trialFields(1,:),'oog');
        stimidx = strcmp(stro.sum.trialFields(1,:),'stim_idx');
        optstim = strcmp(stro.sum.trialFields(1,:),'optstim');
        correct = strcmp(stro.sum.trialFields(1,:),'correct');
        laseron = strcmp(stro.sum.trialFields(1,:),'laseron_t');
        laseroff = strcmp(stro.sum.trialFields(1,:),'laseroff_t');
        correcttrials = stro.trial(:,correct);
        diridxs = stro.trial(:,stimidx);
        stimpresent = logical(stro.trial(:,stimpresentidx));
        LMStriplet = [stro.trial(:,Lcc) stro.trial(:,Mcc) stro.trial(:,Scc)];
        lasertrials = logical(stro.trial(:,optstim));
        if isnan(unique(stro.trial(:,21:22)))
            stimlocations = [stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y];
            RF = [RF; stimlocations];
            newfilename = [newfilename; filename(jj)];
            newRWmultiplier = [newRWmultiplier; RWmultiplier(jj)];
            newlaserdial = [newlaserdial; laserdial(jj)];
            ib = ones(size(lasertrials));
            newstimsize = [newstimsize; stro.sum.exptParams.sigma];
        else
            [stimlocations,~,ib] = unique(stro.trial(:,21:22),'rows');
            RF = [RF; stimlocations];
            L = size(stimlocations,1);
            newfilename = [newfilename; repmat(filename(jj),[L 1])];
            newRWmultiplier = [newRWmultiplier; repmat(RWmultiplier(jj),[L 1])];
            newlaserdial = [newlaserdial; repmat(laserdial(jj),[L 1])];
            newstimsize = [newstimsize; repmat(stro.sum.exptParams.sigma,[L 1])];
        end
        
        for ii = 1:numel(unique(ib))
            ind = logical(ib==ii);
            lasertrialidxs = logical(stro.trial(ind,optstim));
            colordirchoiceidxs = correcttrials(ind);
            colorstimpresent = logical(stimpresent(ind));
            percentcorrectlasertrials = sum(colordirchoiceidxs & lasertrialidxs)/sum(lasertrialidxs);
            percentcorrectnonlasertrials = sum(colordirchoiceidxs & ~lasertrialidxs)/sum(~lasertrialidxs);
            
            % Laser trials
            Hitlasertrial = sum(colordirchoiceidxs & lasertrialidxs & colorstimpresent); % Hit
            Misslasertrial = sum(~colordirchoiceidxs & lasertrialidxs & colorstimpresent); % Miss
            CRlasertrial = sum(colordirchoiceidxs & lasertrialidxs & ~colorstimpresent); % Correct Reject
            FAlasertrial = sum(~colordirchoiceidxs & lasertrialidxs & ~colorstimpresent); % False Alarm
            CRlasertrialidxs = find(lasertrialidxs & ~colorstimpresent); L1 = numel(CRlasertrialidxs);
            
            % Non-Laser trials
            Hitnonlasertrial = sum(colordirchoiceidxs & ~lasertrialidxs & colorstimpresent); % Hit
            Missnonlasertrial = sum(~colordirchoiceidxs & ~lasertrialidxs & colorstimpresent); % Miss
            CRnonlasertrial = sum(colordirchoiceidxs & ~lasertrialidxs & ~colorstimpresent); % Correct Reject
            FAnonlasertrial = sum(~colordirchoiceidxs & ~lasertrialidxs & ~colorstimpresent); % False Alarm
            CRnonlasertrialidxs = find(~lasertrialidxs & ~colorstimpresent); L2 = numel(CRnonlasertrialidxs);
            
            
            stimabsentL = [stimabsentL; sum(CRlasertrial) + sum(FAlasertrial)];
            stimabsentNL = [stimabsentNL; sum(CRnonlasertrial) + sum(FAnonlasertrial)];
            stimpresentL = [stimpresentL; sum(Hitlasertrial) + sum(Misslasertrial)];
            stimpresentNL = [stimpresentNL; sum(Hitnonlasertrial) + sum(Missnonlasertrial)];
            
            % Now storing the Hits, Miss, CR and FA for laser and non-laser trials
            HitL = [HitL; Hitlasertrial]; HitNL = [HitNL; Hitnonlasertrial];
            MissL = [MissL; Misslasertrial]; MissNL = [MissNL; Missnonlasertrial];
            CRL = [CRL; CRlasertrial]; CRNL = [CRNL; CRnonlasertrial];
            FAL = [FAL; FAlasertrial]; FANL = [FANL; FAnonlasertrial];
            TPRNL = [TPRNL; Hitnonlasertrial/(Hitnonlasertrial+Missnonlasertrial)]; % True positive ratio, non-laser trial
            FPRNL = [FPRNL; FAnonlasertrial/(FAnonlasertrial+CRnonlasertrial)]; % False positive ratio, non-laser trial
            TPRL = [TPRL; Hitlasertrial/(Hitlasertrial+Misslasertrial)]; % True positive ratio, laser trial
            FPRL = [FPRL; FAlasertrial/(FAlasertrial+CRlasertrial)]; % False positive ratio, laser trial
            
            fundamentals = stro.sum.exptParams.fundamentals; % CONE FUNDAMENTALS: L,M,S
            fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
            mon_spd = stro.sum.exptParams.mon_spd; % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
            mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
            mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
            M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals
            Lb = M*stro.sum.exptParams.bkgndrgb;
            % Converting all the contrasts to Luminance contrast
            Bkgnd_luminance = stro.sum.exptParams.bkgndrgb'*mon_spd*Vlambda;
            Peak_luminance = (inv(M)*((1+LMStriplet).*(repmat(Lb',[size(LMStriplet,1) 1])))')'*mon_spd*Vlambda;
            Trough_luminance = (inv(M)*((1-LMStriplet).*(repmat(Lb',[size(LMStriplet,1) 1])))')'*mon_spd*Vlambda;
            Contrast_luminance = (Peak_luminance - Trough_luminance)./(Peak_luminance + Trough_luminance);
            %             Contrast_luminance = ((inv(M)*((LMStriplet+1).*(repmat(Lb',[size(LMStriplet,1) 1])))')'*mon_spd*Vlambda-Bkgnd_luminance)/Bkgnd_luminance; % Same as Michelson contrast
            LMScontrastfileNL = Contrast_luminance(colorstimpresent & ~lasertrialidxs,:);
            LMScontrastfileL = Contrast_luminance(colorstimpresent & lasertrialidxs,:);
            % calculating gamut edge luminance contrast
            t = min(((1./stro.sum.exptParams.bkgndrgb)-1));
            OOG_luminance = (stro.sum.exptParams.bkgndrgb*(1+t))' * mon_spd*Vlambda;
            putativeoogcontrast =  (OOG_luminance - Bkgnd_luminance)/Bkgnd_luminance;
            oogvals = [oogvals; putativeoogcontrast];
            
            % Calculating psychophysical detection threshold for non-laser trials
            LMStripletfile = LMStriplet(ind,:);
            trials = lasertrialidxs;
            
            for kk = 1:2
                if kk == 2
                    trials = ~trials;
                end
                LMScontrast = Contrast_luminance(colorstimpresent & trials & ind,:);
                answers = colordirchoiceidxs(colorstimpresent & trials);
                LMScontrast(LMScontrast>putativeoogcontrast) = putativeoogcontrast;
                if kk==1
                    contrastL = unique(LMScontrast);
                    correctanswersL = zeros(size(contrastL));
                    wronganswersL = zeros(size(contrastL));
                    trialspercontrastL = zeros(size(contrastL));
                    for ss = 1:numel(contrastL)
                        trialspercontrastL(ss) = numel(answers(LMScontrast==contrastL(ss)));
                        correctanswersL(ss) = sum(answers(LMScontrast==contrastL(ss)));
                        wronganswersL(ss) = trialspercontrastL(ss) - correctanswersL(ss);
                        percorrectL(ss) = correctanswersL(ss)/trialspercontrastL(ss);
                    end
                    [aL,bL,gL] = weibullFitforDToneloc(contrastL,[correctanswersL wronganswersL],'mle');
                    weibullparamsL = [weibullparamsL; aL bL gL]; % laser trials
                    prop_correctatGAMUTL = [prop_correctatGAMUTL; gL*(1-exp(-((putativeoogcontrast/aL).^bL)))];
                    DetectThreshL = [DetectThreshL; aL];
                else
                    contrastNL = unique(LMScontrast);
                    correctanswersNL = zeros(size(contrastNL));
                    wronganswersNL = zeros(size(contrastNL));
                    trialspercontrastNL = zeros(size(contrastNL));
                    for ss = 1:numel(contrastNL)
                        trialspercontrastNL(ss) = numel(answers(LMScontrast==contrastNL(ss)));
                        correctanswersNL(ss) = sum(answers(LMScontrast==contrastNL(ss)));
                        wronganswersNL(ss) = trialspercontrastNL(ss) - correctanswersNL(ss);
                        percorrectNL(ss) = correctanswersNL(ss)/trialspercontrastNL(ss);
                    end
                    [aNL,bNL,gNL] = weibullFitforDToneloc(contrastNL,[correctanswersNL wronganswersNL],'mle');
                    weibullparamsNL = [weibullparamsNL; aNL bNL gNL]; % laser trials
                    prop_correctatGAMUTNL = [prop_correctatGAMUTNL; gNL*(1-exp(-((putativeoogcontrast/aNL).^bNL)))];
                    
                end
            end
        end
    end
    % Accumulating files over different sessions
    filedates = cell2mat(newfilename);
    datesind = 2:10;
    filedates = str2num(filedates(:,datesind));
    [filedates,~,ib] = unique([filedates RF newRWmultiplier newlaserdial newstimsize],'rows');
    HitLSession = []; HitNLSession = [];
    MissLSession = []; MissNLSession = [];
    CRLSession = []; CRNLSession = [];
    FALSession = []; FANLSession = [];
    stimabsentLSession = [];
    stimabsentNLSession = [];
    stimpresentLSession = [];
    stimpresentNLSession = [];
    p_hits = [];
    p_FA = [];
    p_hitsLFANL = [];
    numfilespersession = [];
    nrows = ceil(sqrt(numel(unique(ib))));
    for ii = 1:numel(unique(ib))
        idx = ii;
        HitLSession = [HitLSession; sum(HitL(idx))];
        HitNLSession = [HitNLSession ;sum(HitNL(idx))];
        MissLSession = [MissLSession; sum(MissL(idx))];
        MissNLSession = [MissNLSession; sum(MissNL(idx))];
        CRLSession = [CRLSession; sum(CRL(idx))];
        CRNLSession = [CRNLSession; sum(CRNL(idx))];
        FALSession = [FALSession; sum(FAL(idx)) ];
        FANLSession = [FANLSession; sum(FANL(idx))];
        stimabsentLSession = [stimabsentLSession; sum(stimabsentL(idx))];
        stimabsentNLSession = [stimabsentNLSession; sum(stimabsentNL(idx))];
        stimpresentLSession = [stimpresentLSession; sum(stimpresentL(idx))];
        stimpresentNLSession = [stimpresentNLSession; sum(stimpresentNL(idx))];
        [~,p1] = equalproptest([sum(HitL(idx)) sum(HitNL(idx))],[sum(stimpresentL(idx)) sum(stimpresentNL(idx))],alpha);
        [~,p2] = equalproptest([sum(FAL(idx)) sum(FANL(idx))],[sum(stimabsentL(idx)) sum(stimabsentNL(idx))],alpha);
        [~,p3] = equalproptest([sum(HitL(idx)) sum(FANL(idx))],[sum(stimpresentL(idx)) sum(stimabsentNL(idx))],alpha);
        p_hits = [p_hits;p1]; p_FA = [p_FA;p2]; p_hitsLFANL = [p_hitsLFANL; p3];
        numfilespersession = [numfilespersession; sum(idx)];
    end
    
    idx = p_hits<1.0;
    HitLSession = HitLSession(idx);
    HitNLSession = HitNLSession(idx);
    MissLSession = MissLSession(idx);
    MissNLSession = MissNLSession(idx);
    CRLSession = CRLSession(idx);
    CRNLSession = CRNLSession(idx);
    FALSession = FALSession(idx);
    FANLSession = FANLSession(idx);
    stimabsentLSession = stimabsentLSession(idx);
    stimabsentNLSession = stimabsentNLSession(idx);
    stimpresentLSession = stimpresentLSession(idx);
    stimpresentNLSession = stimpresentNLSession(idx);
    TPRLSession = HitLSession./(HitLSession+MissLSession);
    TPRNLSession = HitNLSession./(HitNLSession+MissNLSession);
    FPRLSession = FALSession./(FALSession+CRLSession);
    FPRNLSession = FANLSession./(FANLSession+CRNLSession);
    p_hits = p_hits(idx); p_FA = p_FA(idx); p_hitsLFANL = p_hitsLFANL(idx);
    
    figure(plot_counter); set(gcf,'Name',strcat('stats:',monkeyname));
    subplot(1,2,1);
    for ii = 1:numel(p_hits)
        if (p_hits(ii)<0.05)
            color = [0 0 0];
        else
            color = [0.5 0.5 0.5];
        end
        px = HitNLSession(ii)./stimpresentNLSession(ii); errorx = sqrt(px*(1-px)/stimpresentNLSession(ii));
        py = HitLSession(ii)./stimpresentLSession(ii); errory = sqrt(py*(1-py)/stimpresentLSession(ii));
        plot(px,py,'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',color,'MarkerEdgeColor',[1 1 1],'Color','k'); hold on;
    end
    line([0 1],[0 1],'Linewidth',1);axis square; xlabel('Hits control'); ylabel('Hits laser'); set(gca,'Xlim',[0 1],'Ylim',[0 1],'Tickdir','out','XTick',[0:0.5:1],'YTick',[0:0.5:1]);
    if mm == 1
        title('Maui');
    else
        title('Apollo');
    end
    [~,p] = corr((HitNLSession./stimpresentNLSession)-(HitLSession./stimpresentLSession),(CRNLSession./stimabsentNLSession)-(CRLSession./stimabsentLSession));
    ratio_controloverlaser = [ratio_controloverlaser; median(prop_correctatGAMUTNL./prop_correctatGAMUTL)];
    
    dprimeL = norminv(min([1-(0.5./stimpresentLSession) HitLSession./stimpresentLSession],[],2))-norminv(max([0.5./stimabsentLSession FALSession./stimabsentLSession],[],2));
    dprimeNL = norminv(min([1-(0.5./stimpresentNLSession) HitNLSession./stimpresentNLSession],[],2))-norminv(max([0.5./stimabsentNLSession FANLSession./stimabsentNLSession],[],2));
    subplot(1,2,2); plot(dprimeNL,dprimeL,'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on
    set(gca,'Ylim',[1 3],'Xlim',[1 3],'Tickdir','out','XTick',[1:1:3],'YTick',[1:1:3]); line([1 3],[1 3],'Linewidth',1);axis square; xlabel('d prime:control'); ylabel('d prime:laser'); hold off;
    
end
set(gcf,'Name','Laser Control files','renderer','painters');
plot_counter = plot_counter + 1;

p1 = signrank(dprimeNL,dprimeL);
[~,p2] = ttest(dprimeNL,dprimeL);

load laserdetails.mat
laserdial = laserdiallc;
laserdial = spline(laserdetails.dial,laserdetails.laserpower,laserdiallc);

%% No figure : Additonal analyses suggested by Greg
% Comparing control performances between blocks with 0 laser power and some finite laser power
if ~exist('plot_counter')
    plot_counter = 1;
end
load('T_vos1978_Y');
Vlambda = T_vos1978_Y';
filename = {'A021419002.nex';'A021419003.nex';'A021419004.nex';'A021419005.nex';'A021419007.nex';'A021419008.nex';'A021419009.nex'};
RWmultiplier = [1.0;1.0;0.8;0.7;0.8;0.7;0.6];
laserdial = [0;2.4;2.4;2.8;2.4;2.8;1.8];
load laserdetails.mat
laserpower = spline(laserdetails.dial,laserdetails.laserpower,laserdial);

HitL = []; HitNL = [];
MissL = []; MissNL = [];
CRL = []; CRNL = [];
FAL = []; FANL = [];
stimdur = [];
laserdur = [];
h_hit = []; p_hit = [];
h_FA = []; p_FA = [];
TPRNL = []; TPRL = [];
FPRNL = []; FPRL = [];
alpha = 0.05;
RF = [];
stimabsentL = []; stimabsentNL = [];
stimpresentL = []; stimpresentNL = [];
newfilename = [];
newRWmultiplier = [];
newstimsize = [];
newlaserdial = [];
weibullparamsL = [];
weibullparamsNL = [];
prop_correctatGAMUTL = [];
prop_correctatGAMUTNL = [];
DetectThreshL = [];
oogvals = [];
count = 1;
answers_all = [];
LMScontrast_all = [];
for jj = 1:size(filename,1)
    stro = nex2stro(findfile(char(filename(jj,:))));
    fponidx  = strcmp(stro.sum.trialFields(1,:),'fpon_t');
    fpacqidx  = strcmp(stro.sum.trialFields(1,:),'fpacq_t');
    fpoffidx  = strcmp(stro.sum.trialFields(1,:),'fpoff_t');
    stimonidx  = strcmp(stro.sum.trialFields(1,:),'stimon_t');
    stimoffidx  = strcmp(stro.sum.trialFields(1,:),'stimoff_t');
    targonidx  = strcmp(stro.sum.trialFields(1,:),'targon_t');
    saccstartidx  = strcmp(stro.sum.trialFields(1,:),'saccstart_t');
    saccendidx  = strcmp(stro.sum.trialFields(1,:),'saccend_t');
    rewidx  = strcmp(stro.sum.trialFields(1,:),'rew_t');
    stimpresentidx  = strcmp(stro.sum.trialFields(1,:),'stimpresent');
    Lcc = strcmp(stro.sum.trialFields(1,:),'lcc');
    Mcc = strcmp(stro.sum.trialFields(1,:),'mcc');
    Scc = strcmp(stro.sum.trialFields(1,:),'scc');
    tf = strcmp(stro.sum.trialFields(1,:),'tf');
    oog = strcmp(stro.sum.trialFields(1,:),'oog');
    stimidx = strcmp(stro.sum.trialFields(1,:),'stim_idx');
    optstim = strcmp(stro.sum.trialFields(1,:),'optstim');
    correct = strcmp(stro.sum.trialFields(1,:),'correct');
    laseron = strcmp(stro.sum.trialFields(1,:),'laseron_t');
    laseroff = strcmp(stro.sum.trialFields(1,:),'laseroff_t');
    correcttrials = stro.trial(:,correct);
    diridxs = stro.trial(:,stimidx);
    stimpresent = logical(stro.trial(:,stimpresentidx));
    LMStriplet = [stro.trial(:,Lcc) stro.trial(:,Mcc) stro.trial(:,Scc)];
    lasertrials = logical(stro.trial(:,optstim));
    if isnan(unique(stro.trial(:,21:22)))
        stimlocations = [stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y];
        RF = [RF; stimlocations];
        newfilename = [newfilename; filename(jj)];
        newRWmultiplier = [newRWmultiplier; RWmultiplier(jj)];
        newlaserdial = [newlaserdial; laserdial(jj)];
        ib = ones(size(lasertrials));
        newstimsize = [newstimsize; stro.sum.exptParams.sigma];
    else
        [stimlocations,~,ib] = unique(stro.trial(:,21:22),'rows');
        RF = [RF; stimlocations];
        L = size(stimlocations,1);
        newfilename = [newfilename; repmat(filename(jj),[L 1])];
        newRWmultiplier = [newRWmultiplier; repmat(RWmultiplier(jj),[L 1])];
        newlaserdial = [newlaserdial; repmat(laserdial(jj),[L 1])];
        newstimsize = [newstimsize; repmat(stro.sum.exptParams.sigma,[L 1])];
    end
    
    for ii = 1:numel(unique(ib))
        ind = logical(ib==ii);
        lasertrialidxs = logical(stro.trial(ind,optstim));
        colordirchoiceidxs = correcttrials(ind);
        colorstimpresent = logical(stimpresent(ind));
        percentcorrectlasertrials = sum(colordirchoiceidxs & lasertrialidxs)/sum(lasertrialidxs);
        percentcorrectnonlasertrials = sum(colordirchoiceidxs & ~lasertrialidxs)/sum(~lasertrialidxs);
        
        % Laser trials
        Hitlasertrial = sum(colordirchoiceidxs & lasertrialidxs & colorstimpresent); % Hit
        Misslasertrial = sum(~colordirchoiceidxs & lasertrialidxs & colorstimpresent); % Miss
        CRlasertrial = sum(colordirchoiceidxs & lasertrialidxs & ~colorstimpresent); % Correct Reject
        FAlasertrial = sum(~colordirchoiceidxs & lasertrialidxs & ~colorstimpresent); % False Alarm
        CRlasertrialidxs = find(lasertrialidxs & ~colorstimpresent); L1 = numel(CRlasertrialidxs);
        
        % Non-Laser trials
        Hitnonlasertrial = sum(colordirchoiceidxs & ~lasertrialidxs & colorstimpresent); % Hit
        Missnonlasertrial = sum(~colordirchoiceidxs & ~lasertrialidxs & colorstimpresent); % Miss
        CRnonlasertrial = sum(colordirchoiceidxs & ~lasertrialidxs & ~colorstimpresent); % Correct Reject
        FAnonlasertrial = sum(~colordirchoiceidxs & ~lasertrialidxs & ~colorstimpresent); % False Alarm
        CRnonlasertrialidxs = find(~lasertrialidxs & ~colorstimpresent); L2 = numel(CRnonlasertrialidxs);
        
        
        stimabsentL = [stimabsentL; sum(CRlasertrial) + sum(FAlasertrial)];
        stimabsentNL = [stimabsentNL; sum(CRnonlasertrial) + sum(FAnonlasertrial)];
        stimpresentL = [stimpresentL; sum(Hitlasertrial) + sum(Misslasertrial)];
        stimpresentNL = [stimpresentNL; sum(Hitnonlasertrial) + sum(Missnonlasertrial)];
        
        % Now storing the Hits, Miss, CR and FA for laser and non-laser trials
        HitL = [HitL; Hitlasertrial]; HitNL = [HitNL; Hitnonlasertrial];
        MissL = [MissL; Misslasertrial]; MissNL = [MissNL; Missnonlasertrial];
        CRL = [CRL; CRlasertrial]; CRNL = [CRNL; CRnonlasertrial];
        FAL = [FAL; FAlasertrial]; FANL = [FANL; FAnonlasertrial];
        TPRNL = [TPRNL; Hitnonlasertrial/(Hitnonlasertrial+Missnonlasertrial)]; % True positive ratio, non-laser trial
        FPRNL = [FPRNL; FAnonlasertrial/(FAnonlasertrial+CRnonlasertrial)]; % False positive ratio, non-laser trial
        TPRL = [TPRL; Hitlasertrial/(Hitlasertrial+Misslasertrial)]; % True positive ratio, laser trial
        FPRL = [FPRL; FAlasertrial/(FAlasertrial+CRlasertrial)]; % False positive ratio, laser trial
        
        fundamentals = stro.sum.exptParams.fundamentals; % CONE FUNDAMENTALS: L,M,S
        fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
        mon_spd = stro.sum.exptParams.mon_spd; % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
        mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
        mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
        M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals
        Lb = M*stro.sum.exptParams.bkgndrgb;
        % Converting all the contrasts to Luminance contrast
        Bkgnd_luminance = stro.sum.exptParams.bkgndrgb'*mon_spd*Vlambda;
        Peak_luminance = (inv(M)*((1+LMStriplet).*(repmat(Lb',[size(LMStriplet,1) 1])))')'*mon_spd*Vlambda;
        Trough_luminance = (inv(M)*((1-LMStriplet).*(repmat(Lb',[size(LMStriplet,1) 1])))')'*mon_spd*Vlambda;
        Contrast_luminance = (Peak_luminance - Trough_luminance)./(Peak_luminance + Trough_luminance);
        %             Contrast_luminance = ((inv(M)*((LMStriplet+1).*(repmat(Lb',[size(LMStriplet,1) 1])))')'*mon_spd*Vlambda-Bkgnd_luminance)/Bkgnd_luminance; % Same as Michelson contrast
        LMScontrastfileNL = Contrast_luminance(colorstimpresent & ~lasertrialidxs,:);
        LMScontrastfileL = Contrast_luminance(colorstimpresent & lasertrialidxs,:);
        % calculating gamut edge luminance contrast
        t = min(((1./stro.sum.exptParams.bkgndrgb)-1));
        OOG_luminance = (stro.sum.exptParams.bkgndrgb*(1+t))' * mon_spd*Vlambda;
        putativeoogcontrast =  (OOG_luminance - Bkgnd_luminance)/Bkgnd_luminance;
        oogvals = [oogvals; putativeoogcontrast];
        
        % Calculating psychophysical detection threshold for non-laser trials
        LMStripletfile = LMStriplet(ind,:);
        trials = lasertrialidxs;

        if jj==1 
            LMScontrast = Contrast_luminance(colorstimpresent & ind,:);
            answers = colordirchoiceidxs(colorstimpresent);
            LMScontrast(LMScontrast>putativeoogcontrast) = putativeoogcontrast;
            contrastNL = unique(LMScontrast);
            correctanswersNL = zeros(size(contrastNL));
            wronganswersNL = zeros(size(contrastNL));
            trialspercontrastNL = zeros(size(contrastNL));
            percorrectNL = zeros(size(contrastNL));
            for ss = 1:numel(contrastNL)
                trialspercontrastNL(ss) = numel(answers(LMScontrast==contrastNL(ss)));
                correctanswersNL(ss) = sum(answers(LMScontrast==contrastNL(ss)));
                wronganswersNL(ss) = trialspercontrastNL(ss) - correctanswersNL(ss);
                percorrectNL(ss) = correctanswersNL(ss)/trialspercontrastNL(ss);
            end
            [aNL,bNL,gNL] = weibullFitforDToneloc(contrastNL,[correctanswersNL wronganswersNL],'mle'); 
            figure(plot_counter); subplot(223); plot(contrastNL,percorrectNL','o','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
            figure(plot_counter); subplot(224); plot(mean([Contrast_luminance(colorstimpresent & ind & trials,:) Contrast_luminance(colorstimpresent & ind & ~trials,:)],2),'Linewidth',2,'color',[0 0 0]); hold on;
        elseif jj>=2
            LMScontrast = Contrast_luminance(colorstimpresent & ind & ~trials,:);
            LMScontrast(LMScontrast>putativeoogcontrast) = putativeoogcontrast;
            LMScontrast_all = [LMScontrast_all; LMScontrast];
            answers_all = [answers_all; colordirchoiceidxs(colorstimpresent & ~trials)]; 
            figure(plot_counter); subplot(224);  hold on; plot(Contrast_luminance(colorstimpresent & ind & ~trials,:),'Linewidth',2,'color',[0.5 0.5 0.5]);
        end
        if jj == size(filename,1)
            
            contrastNL2 = unique(LMScontrast_all);
            correctanswersNL2 = zeros(size(contrastNL2));
            wronganswersNL2 = zeros(size(contrastNL2));
            trialspercontrastNL2 = zeros(size(contrastNL2));
            percorrectNL2 = zeros(size(contrastNL2));
            for ss = 1:numel(contrastNL2)
                trialspercontrastNL2(ss) = numel(answers_all(LMScontrast_all==contrastNL2(ss)));
                correctanswersNL2(ss) = sum(answers_all(LMScontrast_all==contrastNL2(ss)));
                wronganswersNL2(ss) = trialspercontrastNL2(ss) - correctanswersNL2(ss);
                percorrectNL2(ss) = correctanswersNL2(ss)/trialspercontrastNL2(ss);
            end
            [aNL2,bNL2,gNL2] = weibullFitforDToneloc(contrastNL2,[correctanswersNL2 wronganswersNL2],'mle');
            figure(plot_counter); subplot(223); plot(contrastNL2,percorrectNL2','o','MarkerSize',8,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
        end
    end 
    disp(jj)
end
subplot(224); hold on; xlabel('Trial number'); ylabel('Luminance contrast');  set(gca,'Xlim',[0 30],'XTick',0:15:30,'Ylim',[0 0.8],'YTick',0:0.2:0.8,'Tickdir','out'); axis square;

% Calculating d'
dprimeL = norminv(min([1-(0.5./stimpresentL) HitL./stimpresentL],[],2))-norminv(max([0.5./stimabsentL FAL./stimabsentL],[],2));
dprimeNL = norminv(min([1-(0.5./stimpresentNL) HitNL./stimpresentNL],[],2))-norminv(max([0.5./stimabsentNL FANL./stimabsentNL],[],2));
    

% Bar plots 
bins = 1:0.2:3;
figure(plot_counter);
subplot(221), bar(([mean([HitL(1) HitNL(1)]);HitNL(2:end)]./stimpresentNL)','FaceColor',[0.5 0.5 0.5]); 
set(gca,'Tickdir','out','Xlim',[0 8],'Ylim',[0 1.0],'YTick',0:0.5:1.0); xlabel('Blocks'); ylabel('Prop of Hits'); axis square;
subplot(222); bar([mean([dprimeL(1) dprimeNL(1)]);dprimeNL(2:end)],'FaceColor',[0.5 0.5 0.5]);
set(gca,'Tickdir','out','Xlim',[0 8],'Ylim',[0 3.0],'YTick',0:1.5:3.0); xlabel('Blocks'); ylabel('d prime'); axis square;
contrastlattice = logspace(log10(0.03),log(1.0),51);
fit1 = gNL*(1-exp(-((contrastlattice./aNL).^bNL))); fit2 = gNL2*(1-exp(-((contrastlattice./aNL2).^bNL2)));
subplot(223);plot(contrastlattice,fit1,'Color',[0 0 0],'Linewidth',2); hold on; plot(contrastlattice,fit2,'Color',[0.5 0.5 0.5],'Linewidth',2); xlabel('contrast'); ylabel('Prop. correct'); 
set(gca,'Xlim',[0.03 1.0],'Tickdir','out','YTick',[0:0.25:1],'Xscale','log','XTick',[0.03 0.1 0.3 1],'XTickLabels', {'0.03','0.1','0.3','1'}); axis square; hold off;
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;

% Need to do a Wald's test: Comparing Log-likelihood ratio to the Chi-square cumulative distribution
I = zeros(size([contrastNL; contrastNL2]));
I(1:numel(contrastNL)) = 1;
[aNL_cum,bNL_cum,gNL_cum,~,~,likelihood_cum] = weibullFitforDToneloc([contrastNL; contrastNL2],[[correctanswersNL; correctanswersNL2] [wronganswersNL; wronganswersNL2]],'mle');
[aNL_ind,bNL_ind,gNL_ind,aNL2_ind,~,likelihood_ind] = weibullFitforDToneloc([contrastNL; contrastNL2],[[correctanswersNL; correctanswersNL2] [wronganswersNL; wronganswersNL2]],'mle',[],I);
LLratio = -2*(likelihood_ind-likelihood_cum);
p = 1-chi2cdf(LLratio,1);

fit_cum = gNL_cum*(1-exp(-((contrastlattice./aNL_cum).^bNL_cum)));


