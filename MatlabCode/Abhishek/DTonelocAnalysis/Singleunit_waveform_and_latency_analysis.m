% A new script for checking whether there was any relation between the
% waveform widths and latency of spiking 
% Author - Abhishek De, 6/20
close all; clearvars;

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
example_filename_idxs = ismember(filename,{'A012319005.nex'; 'A020119004.nex'});
Singlevsmultiunit_M = ['M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'M'];
Singlevsmultiunit_A = ['M';'M';'M';'M';'M';'M';'M';'M';'M';'M';'S';'S';'S';'M';'S';'S';'M';'S';'M';'M';'M';'S';'M';'M';'M';'M';'M';'S';'S';'S';'S'];
Singlevsmultiunit = [Singlevsmultiunit_M; Singlevsmultiunit_A];
Singleunitidxs = Singlevsmultiunit == 'S';

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
            PSTHlaser = PSTHlaser + hist(spiketimes-laserontime, bins)/binwidth;
            count1 = count1 + 1;
            laserbins = bins>=0 & bins<=timedurlaser;
            tmp_laserFR = [tmp_laserFR; sum(spiketimes>laserontime & spiketimes<laserofftime)/0.3];
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

% Additional analysis of waveforms
Suppressionvsactivation = zeros(size(baselineFR));
Suppressionvsactivation(baselineFR>laserFR) = 0; % 'Supression'
Suppressionvsactivation(baselineFR<laserFR) = 1; % 'Activation'
idx = find(Singleunitidxs);
N = ceil(sqrt(numel(idx)));
figure(plot_counter); set(gcf,'Name','Plotting waveforms of isolated units');
peak = []; trough = [];
timediff = [];
newfilename = filename;
for ii = 1:numel(idx)
    ind = idx(ii);
    stro = nex2stro(findfile(char(newfilename(ind,:))));
    if ~Suppressionvsactivation(ind)
        color = [0.5 0.5 0.5]; % Suppression
    else
        color = [0 0 0]; % Activation
    end
    samplingrate = 40000;
    leastcount = 10^6/samplingrate; % in us
    waveform = mean(cell2mat(stro.ras(:,2))',2);
    time = 0:25:25*(numel(waveform)-1);
    [val1,t1] = max(waveform);
    [val2,t2] = min(waveform);
    peak = [peak; val1]; 
    trough = [trough; val2];
    timediff = [timediff; time(t2)-time(t1)];
    figure(plot_counter); subplot(N,N,ii); plot(mean(cell2mat(stro.ras(:,2))',2),'-','color',color,'Linewidth',2); hold on;
    line([t1 t1],[0 val1]); line([t2 t2],[0 val2]);line([0 numel(waveform)],[0 0]); hold off;
end
plot_counter = plot_counter + 1;
% Plotting population stats

peaktotroughratio = peak./trough;
p = ranksum(peaktotroughratio(logical(Suppressionvsactivation(idx))),peaktotroughratio(~Suppressionvsactivation(idx)));
peaktotroughratio = peak./abs(trough);

figure(plot_counter); subplot(221); plot(abs(timediff(logical(Suppressionvsactivation(idx)))),peaktotroughratio(logical(Suppressionvsactivation(idx))),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(abs(timediff(~(Suppressionvsactivation(idx)))),peaktotroughratio(~(Suppressionvsactivation(idx))),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); ylabel('Peak Amplitude/Trough Amplitude'); xlabel('|Peak time - Trough time (us)|'); 
set(gca,'Tickdir','out','Xlim',[100 400],'Ylim',[0.3 1.0]); legend('Activation','Suppression'); axis square; hold off;
subplot(222); plot(peak(logical(Suppressionvsactivation(idx))),abs(trough(logical(Suppressionvsactivation(idx)))),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(peak(~(Suppressionvsactivation(idx))),abs(trough(~(Suppressionvsactivation(idx)))),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); xlabel('Peak Amplitude'); ylabel('Trough Amplitude'); 
set(gca,'Tickdir','out','Xlim',[0.02 0.14],'Ylim',[0.02 0.14]); legend('Activation','Suppression'); axis square; hold off;


% Latency of activated units 
% Calculating the time of the first spike of the excited cells (Supplementary figure 1 in the new draft)
newfilename = filename(laserFR>baselineFR);
Singleunitidxs_activated = Singleunitidxs(laserFR>baselineFR);
N = numel(newfilename);
L = ceil(sqrt(N));
timeoffirstspike = cell(N,1);
meantime = [];
stdtime = [];
groupID = [];
for jj = 1:N
    stro = nex2stro(findfile(char(newfilename(jj,:))));
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

figure(plot_counter); 
subplot(223); histogram(1000*meantime,0:2:50,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'Linewidth',2); hold on;
histogram(1000*meantime(Singleunitidxs_activated),0:2:50,'FaceColor',[1 0 0],'EdgeColor',[1 1 1],'Linewidth',2);
set(gca,'Xlim',[0 50],'XTick',0:10:50,'Ylim',[0 6],'YTick',0:3:6,'Tickdir','out'); axis square;
ylabel('# of sites'); xlabel('Average first spike time (ms)'); legend('All units','Single units'); title('Latency of activation');
 hold off; 

subplot(224), plot(1000*meantime(Singleunitidxs_activated),abs(timediff(logical(Suppressionvsactivation(idx)))),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Xlim',[0 30],'Ylim',[100 400]); xlabel('Latency to first spike'); ylabel('|Peak time - Trough time (us)|'); legend('Activated units');
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;