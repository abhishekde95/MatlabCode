% Writing a script for analyzing the effect of laser of on neuronal data (multi-unit)
% Author - Abhishek De, 6/18
close all; clearvars;
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
firstletters_filename = char(filename);
Maui_idxs = firstletters_filename(:,1)=='M';
Apollo_idxs = firstletters_filename(:,1)=='A';
plot_counter = 1;
N = numel(filename); 
L = ceil(sqrt(N));
showstimpresentlasertrials = 0;
showstimabsentnolasertrials = 0;
showstimpresentnolasertrials = 1;
figure(plot_counter); set(gcf,'Name','Stimabsent-laser trials');
if showstimabsentnolasertrials
    figure(plot_counter+1); set(gcf,'Name','Stimabsent-nolaser trials');
end
if showstimpresentlasertrials
    figure(plot_counter+2); set(gcf,'Name','Stimpresent-laser trials');
end
if showstimpresentnolasertrials
    figure(plot_counter+3); set(gcf,'Name','Stimpresent-nolaser trials'); 
end
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
            figure(plot_counter); subplot(L,L,jj); plot(spiketimes(spiketimes>laserontime-0.3 & spiketimes<laserofftime+0.2)-laserontime,count1*ones(size(spiketimes(spiketimes>laserontime-0.3 & spiketimes<laserofftime+0.2))),'k.'); hold on;
            PSTHlaser = PSTHlaser + hist(spiketimes-laserontime, bins)/binwidth;
            count1 = count1 + 1;
            laserbins = bins>=0 & bins<=timedurlaser;
%             tmp_baselineFR = [tmp_baselineFR; mean(PSTHlaser(baselinebins))/0.3];
            tmp_laserFR = [tmp_laserFR; binwidth*mean(PSTHlaser(laserbins))/timedurlaser];
            ISI = [ISI; diff(spiketimes(spiketimes>laserontime & spiketimes<laserofftime))]; 
            tmp_spikecountsbaseline = [tmp_spikecountsbaseline; sum(spiketimes>laserontime-0.3 & spiketimes<laserontime)];
            tmp_spikecountslaserstimabsent = [tmp_spikecountslaserstimabsent; sum(spiketimes>laserontime & spiketimes<laserofftime)];
        else 
            if showstimabsentnolasertrials
                stimontime = stro.trial(ind,stimonidx);
                stimofftime = stro.trial(ind,stimoffidx);
                timedurnolaser = stimofftime - stimontime;
                figure(plot_counter+1); subplot(L,L,jj); plot(spiketimes(spiketimes>stimontime-0.3 & spiketimes<stimofftime+0.2)-stimontime,count2*ones(size(spiketimes(spiketimes>stimontime-0.3 & spiketimes<stimofftime+0.2))),'k.'); hold on;
                PSTHnolaser = PSTHnolaser + hist(spiketimes-stimontime, bins)/binwidth;
                count2 = count2 + 1;
            end
            if ~lasertrials(ind) & ~stimpresent(ind)
%                 keyboard;
                stimontime = stro.trial(ind,stimonidx);
                PSTHbaseline = PSTHbaseline + hist(spiketimes-stimontime, bins)/binwidth;
                tmp_baselineFR = [tmp_baselineFR; binwidth*mean(PSTHbaseline(baselinebins))/0.3];
            end
        end
    end
    baselineFR = [baselineFR; mean(tmp_baselineFR)];
    laserFR = [laserFR; mean(tmp_laserFR)];
    figure(plot_counter); subplot(L,L,jj); plot(bins,PSTHlaser/sum(lasertrials & ~stimpresent),'-','color',[0 0.5 1],'Linewidth',2); 
    line([0 0],[0 count1]); line([timedurlaser timedurlaser],[0 count1]); xlabel('time'); ylabel('trials'); title(char(filename{jj}(1:10))); set(gca,'Xlim',[-0.3 0.5]); hold off;
    if showstimabsentnolasertrials
        figure(plot_counter+1); subplot(L,L,jj); plot(bins,PSTHnolaser/sum(~lasertrials & ~stimpresent),'-','color',[0.5 0.5 0.5],'Linewidth',2);
        line([0 0],[0 count2]); line([timedurnolaser timedurnolaser],[0 count2]); xlabel('time'); ylabel('trials'); title(char(filename{jj}(1:10))); set(gca,'Xlim',[-0.3 0.5]); hold off;
    end
    statstestresult = [statstestresult; ranksum(tmp_spikecountsbaseline,tmp_spikecountslaserstimabsent)];
    
    % stimulus present trials 
    idxs = find(stimpresent);
    count3 = 1;
    count4 = 1;
    PSTHlaser = zeros(1,length(bins));
    PSTHnolaser = zeros(1,length(bins));
    for ii = 1:numel(idxs)
        ind = idxs(ii);
        analogstartime = stro.ras{ind,analogstrtimeind};
        spiketimes = stro.ras{ind,spikeind};
        if lasertrials(ind)
            if showstimpresentlasertrials
                t = stro.ras{ind,analogstrtimeind}+[0:1:length(stro.ras{ind,lasertraceind})-1]/samplerate;
                laserontime = t(find(stro.ras{ind,lasertraceind}>0.1,1));
                laserofftime = t(find(stro.ras{ind,lasertraceind}>0.1,1,'last'));
                timedurlaser = laserofftime - laserontime;
                figure(plot_counter+2); subplot(L,L,jj); plot(spiketimes(spiketimes>laserontime-0.2 & spiketimes<laserofftime+0.2)-laserontime,count3*ones(size(spiketimes(spiketimes>laserontime-0.2 & spiketimes<laserofftime+0.2))),'k.'); hold on;
                PSTHlaser = PSTHlaser + hist(spiketimes-laserontime, bins)/binwidth;
                count3 = count3 + 1;
            end
        else
            if showstimpresentnolasertrials
                stimontime = stro.trial(ind,stimonidx);
                stimofftime = stro.trial(ind,stimoffidx);
                timedurnolaser = stimofftime - stimontime;
                figure(plot_counter+3); subplot(L,L,jj); plot(spiketimes(spiketimes>stimontime-0.2 & spiketimes<stimofftime+0.2)-stimontime,count4*ones(size(spiketimes(spiketimes>stimontime-0.2 & spiketimes<stimofftime+0.2))),'k.'); hold on;
                PSTHnolaser = PSTHnolaser + hist(spiketimes-stimontime, bins)/binwidth;
                count4 = count4 + 1;
            end
        end
    end
    if showstimpresentlasertrials
        figure(plot_counter+2); subplot(L,L,jj); plot(bins,PSTHlaser/sum(lasertrials & stimpresent),'-','color',[0 0.5 1],'Linewidth',2);
        line([0 0],[0 count3]); line([timedurlaser timedurlaser],[0 count3]); xlabel('time'); ylabel('trials'); title(char(filename{jj}(1:10))); set(gca,'Xlim',[-0.2 0.5]); hold off;
    end
    if showstimpresentnolasertrials
        figure(plot_counter+3); subplot(L,L,jj); plot(bins,PSTHnolaser/sum(~lasertrials & stimpresent),'-','color',[0.5 0.5 0.5],'Linewidth',2);
        line([0 0],[0 count4]); line([timedurnolaser timedurnolaser],[0 count4]); xlabel('time'); ylabel('trials'); title(char(filename{jj}(1:10))); set(gca,'Xlim',[-0.2 0.5]); hold off;
    end
end
plot_counter = plot_counter + 4;

Singlevsmultiunit = [Singlevsmultiunit_M; Singlevsmultiunit_A];
Singleunitidxs = Singlevsmultiunit == 'S';
figure(plot_counter); set(gcf,'Name','FR analysis:laser vs baseline'); 
subplot(121); plot(baselineFR(Apollo_idxs),laserFR(Apollo_idxs),'o','MarkerSize',7,'LineWidth',1.0,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(baselineFR(Maui_idxs),laserFR(Maui_idxs),'o','MarkerSize',7,'LineWidth',1.0,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]); hold on; 
axis square; xlabel('baseline FR'); ylabel('laser FR'); set(gca,'Xlim',[0 600],'Ylim',[0 600]); line([0 600],[0 600]); legend('A','M'); hold off;
subplot(122); plot(baselineFR(~Singleunitidxs),laserFR(~Singleunitidxs),'o','MarkerSize',7,'LineWidth',1.0,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(baselineFR(Singleunitidxs),laserFR(Singleunitidxs),'o','MarkerSize',7,'LineWidth',1.0,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on; 
axis square; xlabel('baselineFR'); ylabel('laserFR'); set(gca,'Xlim',[0 600],'Ylim',[0 600],'TickDir','out'); line([0 600],[0 600]); legend('Multi','Single'); hold off;
plot_counter = plot_counter + 1;

%% Marking suppression vs activation indices
Suppressionvsactivation = zeros(size(baselineFR));
Suppressionvsactivation(baselineFR>laserFR) = 0; % 'Supression'
Suppressionvsactivation(baselineFR<laserFR) = 1; % 'Activation'

% Neuronal contrast response functions from stimpresent, control trials
figure(plot_counter); set(gcf,'Name','Spikecounts: Suppression');
figure(plot_counter+1); set(gcf,'Name','Spikecounts: Activation');
L1 = ceil(sqrt(sum(Suppressionvsactivation==0)));
L2 = ceil(sqrt(sum(Suppressionvsactivation==1)));
count1 = 1;
count2 = 1;
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
    lasertraceind = strcmp(stro.sum.rasterCells(1,:),'AD13');
    correcttrials = stro.trial(:,correct);
    diridxs = stro.trial(:,stimidx);
    stimpresent = logical(stro.trial(:,stimpresentidx));
    LMStriplet = [stro.trial(:,Lcc) stro.trial(:,Mcc) stro.trial(:,Scc)];
    lasertrials = logical(stro.trial(:,optstim));
    LMScontrast = sqrt(sum(LMStriplet.^2,2));
    putativeoogcontrast = min(LMScontrast(logical(stro.trial(:,oog))));
    if ~isempty(putativeoogcontrast)
        LMScontrast(LMScontrast>putativeoogcontrast) = putativeoogcontrast;
    end
    % stimulus absent no-laser trials 
    idxs = find(~stimpresent & ~lasertrials);
    spikecounts = [];
    for ii = 1:numel(idxs)
        ind = idxs(ii);
        analogstartime = stro.ras{ind,analogstrtimeind};
        spiketimes = stro.ras{ind,spikeind};
        stimontime = stro.trial(ind,stimonidx);
        stimofftime = stro.trial(ind,stimoffidx);
        stimdur = stimofftime - stimontime;
        spikecounts = [spikecounts; numel(spiketimes(spiketimes>stimontime & spiketimes<stimofftime))/stimdur]; % Storing the avg firing rate
    end
    spikecounts_baseline = mean(spikecounts);
    stderr_baseline = std(spikecounts_baseline)/sqrt(numel(spikecounts_baseline));    
    
    % stimulus present trials 
    idxs = find(stimpresent & ~lasertrials);
    spikecounts = [];
    LMScontrast_stimpresentcontrol = LMScontrast(stimpresent & ~lasertrials);
    for ii = 1:numel(idxs)
        ind = idxs(ii);
        analogstartime = stro.ras{ind,analogstrtimeind};
        spiketimes = stro.ras{ind,spikeind};
        stimontime = stro.trial(ind,stimonidx);
        stimofftime = stro.trial(ind,stimoffidx);
        stimdur = stimofftime - stimontime;
        spikecounts = [spikecounts; numel(spiketimes(spiketimes>stimontime & spiketimes<stimofftime))/stimdur]; % Storing the avg firing rate
    end
    
    [A,~,B] = unique(LMScontrast_stimpresentcontrol);
    avg_spikecounts = zeros(size(A)); std_error = zeros(size(A));
    for i = 1:max(B)
        avg_spikecounts(i) = mean(spikecounts(B==i));
        std_error(i) = std(spikecounts(B==i))/sqrt(sum(B==i));
    end
    A = [0; A];
    avg_spikecounts = [spikecounts_baseline; avg_spikecounts];
    std_error = [stderr_baseline; std_error];
    if ~Suppressionvsactivation(jj)
        color = [0.5 0.5 0.5];
        figure(plot_counter); subplot(L1,L1,count1); errorbar(A, avg_spikecounts, std_error,'o','MarkerFaceColor',color,'MarkerSize',5); title(char(filename{jj}(1:10)));
        count1 = count1 + 1;
    else
        color = [0 0 0];
        figure(plot_counter+1); subplot(L2,L2,count2); errorbar(A, avg_spikecounts, std_error,'o','MarkerFaceColor',color,'MarkerSize',5); title(char(filename{jj}(1:10)));
        count2 = count2 + 1;
    end
end
plot_counter = plot_counter + 2;

%% Isolating the waveforms for single units 

idx = find(Singlevsmultiunit=='S');
N = ceil(sqrt(numel(idx)));
figure(plot_counter); set(gcf,'Name','Plotting waveforms of isolated units');
peak = []; trough = [];
timediff = [];
for ii = 1:numel(idx)
    ind = idx(ii);
    stro = nex2stro(findfile(char(filename(ind,:))));
    if ~Suppressionvsactivation(ind)
        color = [0.5 0.5 0.5]; % Suppression
    else
        color = [0 0 0]; % Activation
    end
    samplingrate = stro.sum.waves.storeRates{1};
    leastcount = 10^6/samplingrate; % in us
    waveform = mean(cell2mat(stro.ras(:,2))',2);
    time = 0:25:25*(numel(waveform)-1);
    [val1,t1] = max(waveform);
    [val2,t2] = min(waveform);
    peak = [peak; val1]; 
    trough = [trough; val2];
    timediff = [timediff; abs(time(t2)-time(t1))];
    figure(plot_counter); subplot(N,N,ii); plot(mean(cell2mat(stro.ras(:,2))',2),'-','color',color,'Linewidth',2); hold on;
    line([t1 t1],[0 val1]); line([t2 t2],[0 val2]);line([0 numel(waveform)],[0 0]); hold off;
end
plot_counter = plot_counter + 1;

% Plotting population stats
peaktotroughratio = peak./abs(trough);
figure(plot_counter); subplot(121); plot(timediff(logical(Suppressionvsactivation(idx))),peaktotroughratio(logical(Suppressionvsactivation(idx))),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(timediff(~(Suppressionvsactivation(idx))),peaktotroughratio(~(Suppressionvsactivation(idx))),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); ylabel('Peak/trough ratio'); xlabel('peak time - trough time (us)'); 
legend('Activation','Suppression'); axis square; hold off;
subplot(122); plot(peak(logical(Suppressionvsactivation(idx))),abs(trough(logical(Suppressionvsactivation(idx)))),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(peak(~(Suppressionvsactivation(idx))),abs(trough(~(Suppressionvsactivation(idx)))),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); xlabel('Peak amp'); ylabel('trough amp'); 
legend('Activation','Suppression'); axis square; hold off;
plot_counter = plot_counter + 1;