% A script for eLife rebuttal
% Author - Abhishek De
% This contains an earlier version of rebuttal figures.
% For more organized version check out Popfigures_Elife_rebuttal.m

%% Neurophysiology: Further characterization of electrophysiological responses
close all; clearvars;
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
visualstimFR_stimp = [];
visualstimbaselineFR_stimp = [];
statstestresult_stimp = [];
ISI = [];
HitsL = []; MissesL = []; CRL = []; FAL = []; stimpresentL = []; stimabsentL = [];
HitsNL = []; MissesNL = []; CRNL = []; FANL = []; stimpresentNL = []; stimabsentNL = [];
dprimeL = []; dprimeNL = [];
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
    
    % Laser trials
    HitsL = [HitsL; sum(stimpresent & lasertrials & correcttrials)];
    MissesL = [MissesL; sum(stimpresent & lasertrials & ~correcttrials)];
    CRL = [CRL; sum(~stimpresent & lasertrials & correcttrials)];
    FAL = [FAL; sum(~stimpresent & lasertrials & ~correcttrials)];
    stimpresentL = [stimpresentL; sum(stimpresent & lasertrials)];
    stimabsentL = [stimabsentL; sum(~stimpresent & lasertrials)];
    
    % No-laser trials
    HitsNL = [HitsNL; sum(stimpresent & ~lasertrials & correcttrials)];
    MissesNL = [MissesNL; sum(stimpresent & ~lasertrials & ~correcttrials)];
    CRNL = [CRNL; sum(~stimpresent & ~lasertrials & correcttrials)];
    FANL = [FANL; sum(~stimpresent & ~lasertrials & ~correcttrials)];
    stimpresentNL = [stimpresentNL; sum(stimpresent & ~lasertrials)];
    stimabsentNL = [stimabsentNL; sum(~stimpresent & ~lasertrials)];
    
    % Calculating the dprimeL and dprimeNL
    dprimeL = [dprimeL; calcdprime(HitsL(end),FAL(end),stimpresentL(end),stimabsentL(end))];
    dprimeNL = [dprimeNL; calcdprime(HitsNL(end),FANL(end),stimpresentNL(end),stimabsentNL(end))];
    
    % Analyses of stimulus absent trials
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
    
    % Analyses of stimulus prsent trials
    idxs = find(stimpresent);
    count1 = 1;
    count2 = 1;
    tmp_visualstimFR = [];
    tmp_visualstimbaselineFR = [];
    
    
    for ii = 1:numel(idxs)
        ind = idxs(ii);
        analogstartime = stro.ras{ind,analogstrtimeind};
        spiketimes = stro.ras{ind,spikeind};
        if ~lasertrials(ind) % looking at the no-laser trials
            if stimpresent(ind)
                stimontime = stro.trial(ind,stimonidx);
                stimofftime = stro.trial(ind,stimoffidx);
                timedur = stimofftime-stimontime;
                tmp_visualstimFR = [tmp_visualstimFR; sum(spiketimes>stimontime & spiketimes<stimofftime)/timedur ];
                
            end
        end
    end
    visualstimFR_stimp = [visualstimFR_stimp; mean(tmp_visualstimFR)];
    statstestresult_stimp = [statstestresult_stimp; ranksum(tmp_visualstimFR,tmp_baselineFR)];
end

firstletters_filename = char([filename]);
Singlevsmultiunit = [Singlevsmultiunit_M; Singlevsmultiunit_A];
Singleunitidxs = Singlevsmultiunit == 'S';
Maui_idxs = firstletters_filename(:,1)=='M';
Apollo_idxs = firstletters_filename(:,1)=='A';

% PLotting in log scale
% Since 0 cannot be represented in log scale, I am changing the zeros to 0.01
baselineFRlog = baselineFR;
baselineFRlog(baselineFRlog==0) = 0.1;
laserFRlog = laserFR;
laserFRlog(laserFRlog==0) = 0.1;
figure(plot_counter);
subplot(221); plot(baselineFRlog(~Singleunitidxs & Maui_idxs),laserFRlog(~Singleunitidxs & Maui_idxs),'s','MarkerSize',8,'LineWidth',1.0,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(baselineFRlog(Singleunitidxs & Maui_idxs),laserFRlog(Singleunitidxs & Maui_idxs),'s','MarkerSize',8,'LineWidth',1.0,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(baselineFRlog(~Singleunitidxs & Apollo_idxs),laserFRlog(~Singleunitidxs & Apollo_idxs),'o','MarkerSize',6,'LineWidth',1.0,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(baselineFRlog(Singleunitidxs & Apollo_idxs),laserFRlog(Singleunitidxs & Apollo_idxs),'o','MarkerSize',6,'LineWidth',1.0,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
axis square; xlabel('baselineFR'); ylabel('laserFR'); set(gca,'Xlim',[0.1 300],'Ylim',[0.1 300],'TickDir','out','XScale','log','YScale','log'); line([0.1 300],[0.1 300],'Color','k'); 
legend('M:multi','A:multi','A:single'); title('Opto stim modulation'); hold off;

% Linking neural activity to behavior
Beh = (dprimeNL-dprimeL)./(dprimeNL+dprimeL);
Neur = abs((laserFR-baselineFR)./(laserFR+baselineFR));
subplot(222); plot(Neur(laserFR>=baselineFR),Beh(laserFR>=baselineFR),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(Neur(laserFR<baselineFR),Beh(laserFR<baselineFR),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); 
legend('Excited','Suppressed'); xlabel('Neuronal laser modulation'); 
ylabel('Behavioral laser modulation'); 
% ylabel('d prime control - laser');
axis square; set(gca,'Tickdir','out','Xlim',[0.5 1]); 
[r,p] = corr(Neur,Beh,'type','Spearman');

% Looking at the stim present trials
visualstimFR_stimplog = visualstimFR_stimp;
visualstimFR_stimplog(visualstimFR_stimplog==0) = 0.1;
visualstimbaselineFR_stimp = baselineFR;
visualstimbaselineFR_stimplog = visualstimbaselineFR_stimp; 
visualstimbaselineFR_stimplog(visualstimbaselineFR_stimplog==0) = 0.1;
subplot(223); plot(visualstimbaselineFR_stimplog(~Singleunitidxs & Maui_idxs & statstestresult_stimp<0.05),visualstimFR_stimplog(~Singleunitidxs & Maui_idxs & statstestresult_stimp<0.05),'s','MarkerSize',8,'LineWidth',1.0,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(visualstimbaselineFR_stimplog(Singleunitidxs & Maui_idxs & statstestresult_stimp<0.05),visualstimFR_stimplog(Singleunitidxs & Maui_idxs & statstestresult_stimp<0.05),'s','MarkerSize',8,'LineWidth',1.0,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(visualstimbaselineFR_stimplog(~Singleunitidxs & Apollo_idxs & statstestresult_stimp<0.05),visualstimFR_stimplog(~Singleunitidxs & Apollo_idxs & statstestresult_stimp<0.05),'o','MarkerSize',6,'LineWidth',1.0,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(visualstimbaselineFR_stimplog(Singleunitidxs & Apollo_idxs & statstestresult_stimp<0.05),visualstimFR_stimplog(Singleunitidxs & Apollo_idxs & statstestresult_stimp<0.05),'o','MarkerSize',6,'LineWidth',1.0,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(visualstimbaselineFR_stimplog(~Singleunitidxs & Maui_idxs & statstestresult_stimp>0.05),visualstimFR_stimplog(~Singleunitidxs & Maui_idxs & statstestresult_stimp>0.05),'s','MarkerSize',8,'LineWidth',1.0,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1]); hold on;
plot(visualstimbaselineFR_stimplog(Singleunitidxs & Maui_idxs & statstestresult_stimp>0.05),visualstimFR_stimplog(Singleunitidxs & Maui_idxs & statstestresult_stimp>0.05),'s','MarkerSize',8,'LineWidth',1.0,'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[1 1 1]); hold on;
plot(visualstimbaselineFR_stimplog(~Singleunitidxs & Apollo_idxs & statstestresult_stimp>0.05),visualstimFR_stimplog(~Singleunitidxs & Apollo_idxs & statstestresult_stimp>0.05),'o','MarkerSize',6,'LineWidth',1.0,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1]); hold on;
plot(visualstimbaselineFR_stimplog(Singleunitidxs & Apollo_idxs & statstestresult_stimp>0.05),visualstimFR_stimplog(Singleunitidxs & Apollo_idxs & statstestresult_stimp>0.05),'o','MarkerSize',6,'LineWidth',1.0,'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[1 1 1]); hold on;
axis square; xlabel('baselineFR'); ylabel('stimpresentFR: no laser'); set(gca,'Xlim',[0.1 100],'Ylim',[0.1 100],'TickDir','out','XScale','log','YScale','log'); line([0.1 100],[0.1 100],'Color','k'); legend('M:multi','A:multi','A:single'); title('Visual stim modulation'); hold off;

subplot(224); plot(visualstimbaselineFR_stimplog(laserFR>=baselineFR & statstestresult_stimp<0.05),visualstimFR_stimplog(laserFR>=baselineFR & statstestresult_stimp<0.05),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(visualstimbaselineFR_stimplog(laserFR<baselineFR & statstestresult_stimp<0.05),visualstimFR_stimplog(laserFR<baselineFR & statstestresult_stimp<0.05),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(visualstimbaselineFR_stimplog(laserFR>=baselineFR & statstestresult_stimp>=0.05),visualstimFR_stimplog(laserFR>=baselineFR & statstestresult_stimp>=0.05),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1]); hold on;
plot(visualstimbaselineFR_stimplog(laserFR<baselineFR & statstestresult_stimp>=0.05),visualstimFR_stimplog(laserFR<baselineFR & statstestresult_stimp>=0.05),'o','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[1 1 1]); hold on;
axis square; xlabel('baselineFR'); ylabel('stimpresentFR: no laser'); set(gca,'Xlim',[0.1 100],'Ylim',[0.1 100],'TickDir','out','XScale','log','YScale','log'); line([0.1 100],[0.1 100],'Color','k'); legend('excited','suppressed'); title('Visual stim modulation'); hold off;
plot_counter = plot_counter + 1;

p1 = ranksum(baselineFR(laserFR>=baselineFR),baselineFR(laserFR<baselineFR)); % Comparing the firing rates of the excited and suppressed sites 

%% Analyses of suppressed and excited sites 
% Calculate the latency of the suppressed sites 
filenameS = filename(laserFR<=baselineFR); % Pulling out files with suppressed sites 
filenameE = filename(laserFR>baselineFR); % Pulling out the files with excited sites 
N = numel(filename);
L = ceil(sqrt(N));
binwidth = .001;
bins = -0.3:binwidth:0.6;
countS = 1;
countE = 1;
ChangeinFR_time = [];
pFR = [];
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
    count = 1;
    PSTH = [];
    for ii = 1:numel(idxs)
        ind = idxs(ii);
        analogstartime = stro.ras{ind,analogstrtimeind};
        spiketimes = stro.ras{ind,spikeind};
        if lasertrials(ind)
            % Pulling out the laseron and off timings from analog traces
            t = stro.ras{ind,analogstrtimeind}+[0:1:length(stro.ras{ind,lasertraceind})-1]/samplerate;
            laserontime = t(find(stro.ras{ind,lasertraceind}>0.1,1));
            laserofftime = t(find(stro.ras{ind,lasertraceind}>0.1,1,'last'));
            timedurlaser = laserofftime-laserontime;
            spiketodisp = spiketimes(spiketimes > laserontime-0.3 & spiketimes<laserofftime+0.3);
            % Plotting the spikes
            if ismember(filename(jj),filenameS)
                figure(plot_counter); subplot(5,4,countS); plot(spiketodisp-laserontime, count*ones(size(spiketodisp)), 'k.'); hold on;
            else
                figure(plot_counter+2); subplot(6,7,countE); plot(spiketodisp-laserontime, count*ones(size(spiketodisp)), 'k.'); hold on;
            end
            count = count + 1;
            PSTH = [PSTH; hist(spiketimes-laserontime, bins)/binwidth];
        else
        end
    end
    PSTH(:,1) = 0; PSTH(:,end) = 0;
    if ismember(filename(jj),filenameS)
        % Suppressed sites 
        % Plotting the spikes
        figure(plot_counter); subplot(5,4,countS);  line([0 0],[0 30]); line([timedurlaser timedurlaser],[0 30]);
        xlabel('time'); ylabel('trials'); title(char(filename{jj}(1:10))); set(gca,'Xlim',[-0.1 0.1],'XTick',[-0.1 0 0.1]); hold off;
        
        % Plotting the PSTH
        figure(plot_counter+1); subplot(5,4,countS); plot(bins,mean(PSTH,1),'Linewidth',2,'Color','k');
        xlabel('time'); ylabel('Ips'); title(char(filename{jj}(1:10))); set(gca,'Xlim',[-0.1 0.1],'XTick',[-0.1 0 0.1]); hold off;
        countS = countS + 1;
    else
        % Excited sites 
        % Plotting the spikes
        figure(plot_counter+2); subplot(6,7,countE);  line([0 0],[0 30]); line([timedurlaser timedurlaser],[0 30]);
        xlabel('time'); ylabel('trials'); title(char(filename{jj}(1:10))); set(gca,'Xlim',[-0.1 0.1],'XTick',[-0.1 0 0.1]); hold off;
        
        % Plotting the PSTH
        figure(plot_counter+3); subplot(6,7,countE); plot(bins,mean(PSTH,1),'Linewidth',2,'Color','k');
        xlabel('time'); ylabel('Ips'); title(char(filename{jj}(1:10))); set(gca,'Xlim',[-0.1 0.1],'XTick',[-0.1 0 0.1]); hold off;
        
        countE = countE + 1;
    end
    
    %%%%%% For calculating the time of excitation and suppression onset
    timewindow = 0.05;
    reference_baselinebin = bins>-1*timewindow & bins<0;
    t = -timewindow:0.001:timewindow;
    p = [];
    
    meanbaseline_spikecount = mean(sum(PSTH(:,reference_baselinebin)*binwidth,2));
    for ii = 1:numel(t)
        t1 = t(ii);
        [p_tmp,h_tmp] = signrank(sum(PSTH(:,reference_baselinebin)*binwidth,2),sum(PSTH(:,bins>t1 & bins<t1+timewindow)*binwidth,2));
        p = [p p_tmp]; 
    end
    t = t+timewindow;
    t2 = t(find(p<0.05,1));
    if isempty(t2)
        t2 = t(end);
    end
    ChangeinFR_time = [ChangeinFR_time; t2];
    pFR = [pFR; p]; % Storing the p-values from the sign rank test 
    
end
plot_counter = plot_counter + 4;

% plotting the results of latency 
figure(plot_counter); set(gcf,'Name','Latency'); 
subplot(121); histogram(ChangeinFR_time(ismember(filename,filenameE)),logspace(-3,-1,31),'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); hold on;
histogram(ChangeinFR_time(ismember(filename,filenameS)),logspace(-3,-1,31),'FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1]);
set(gca,'Tickdir','out','XScale','log'); axis square; xlabel('time (s)'); ylabel('#cells'); title('Latency of FR change: Histogram'); legend('excited','suppressed'); hold off;
subplot(122); plot(t,mean(pFR(ismember(filename,filenameE),:)',2),'Color',[0 0 0],'Linewidth',2); hold on;
plot(t,mean(pFR(ismember(filename,filenameS),:)',2),'Color',[0.5 0.5 0.5],'Linewidth',2); 
plot(t,median(pFR(ismember(filename,filenameE),:)',2),'-.','Color',[0 0 0],'Linewidth',2);
plot(t,median(pFR(ismember(filename,filenameS),:)',2),'-.','Color',[0.5 0.5 0.5],'Linewidth',2);
plot(t,0.05*ones(size(t)),'Color','g','Linewidth',2);
set(gca,'Tickdir','out','XScale','log','YScale','log'); axis square; xlabel('time (s)'); ylabel('p value'); legend('excited: mean','suppressed: mean','excited: median','suppressed: median'); 
title('Latency of FR change: Population'); hold off;
plot_counter = plot_counter + 1;

%% Population Analysis of SCstimcue data
% Figure 4- Supplement 1 
close all;
clearvars;
if ~exist('plot_counter')
    plot_counter = 1;
end
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

[p3] = signrank(cell2mat(probsacc_withinRFL(idx)),cell2mat(probsacc_withinRFNL(idx))); % For Maui: within RF
[p4] = signrank(cell2mat(probsacc_outsideRFL(idx)),cell2mat(probsacc_outsideRFNL(idx))); % For Maui: outside RF
[p5] = signrank(cell2mat(probsacc_withinRFL(~idx)),cell2mat(probsacc_withinRFNL(~idx))); % For Apollo: within RF
[p6] = signrank(cell2mat(probsacc_outsideRFL(~idx)),cell2mat(probsacc_outsideRFNL(~idx))); % For Apollo: outside RF

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
Saccadelatency_catchL = []; Saccadelatency_catchNL = [];
SaccadeAmpwithinL = []; SaccadeAmpwithinNL = [];
SaccadeAmpoutsideL = []; SaccadeAmpoutsideNL = [];
DistancefromtargetwithinL = []; DistancefromtargetwithinNL = [];
DistancefromtargetoutsideL = []; DistancefromtargetoutsideNL = [];
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
    SaccadeAmpwithinL = [SaccadeAmpwithinL; norm(centroidwithinL)];
    SaccadeAmpwithinNL = [SaccadeAmpwithinNL; norm(centroidwithinNL)];
    DistancefromtargetwithinL = [DistancefromtargetwithinL; norm(centroidwithinL - RFloc/10)]; 
    DistancefromtargetwithinNL = [DistancefromtargetwithinNL; norm(centroidwithinNL - RFloc/10)];
    diffcentroidwithin = [diffcentroidwithin; sum((centroidwithinL - centroidwithinNL).^2)];
    Saccadelatency_withinL = [Saccadelatency_withinL; Saccadelatency(ib==withinRFidx & Llaser)];
    Saccadelatency_withinNL = [Saccadelatency_withinNL; Saccadelatency(ib==withinRFidx & ~Llaser)];
    
    centroidoutsideL = []; centroidoutsideNL = [];
    outsideL = []; outsideNL = [];
    tmpSaccadeAmpoutsideL = []; tmpSaccadeAmpoutsideNL = [];
    tmpDistancefromtargetoutsideL = []; tmpDistancefromtargetoutsideNL = [];
    for jj = 1:numel(outsideRFidx)
        centroidoutsideL = [centroidoutsideL; mean(SaccadeX(ib==outsideRFidx(jj) & Llaser)) mean(SaccadeY(ib==outsideRFidx(jj) & Llaser))];
        centroidoutsideNL = [centroidoutsideNL; mean(SaccadeX(ib==outsideRFidx(jj) & ~Llaser)) mean(SaccadeY(ib==outsideRFidx(jj) & ~Llaser))];
        outsideL = [outsideL; Saccadelatency(ib==outsideRFidx(jj) & Llaser)];
        outsideNL = [outsideNL; Saccadelatency(ib==outsideRFidx(jj) & ~Llaser)];
        tmpSaccadeAmpoutsideL = [tmpSaccadeAmpoutsideL; norm(centroidoutsideL(end,:))];
        tmpSaccadeAmpoutsideNL = [tmpSaccadeAmpoutsideNL; norm(centroidoutsideNL(end,:))];
        tmpDistancefromtargetoutsideL = [tmpDistancefromtargetoutsideL; norm(centroidoutsideL(end,:)-(C(outsideRFidx(jj),:)/10))];
        tmpDistancefromtargetoutsideNL = [tmpDistancefromtargetoutsideNL; norm(centroidoutsideNL(end,:)-(C(outsideRFidx(jj),:)/10))];
    end
    diffcentroidoutside = [diffcentroidoutside; mean(sum((centroidoutsideL - centroidoutsideNL).^2),2)];
    Saccadelatency_outsideL = [Saccadelatency_outsideL; outsideL];
    Saccadelatency_outsideNL = [Saccadelatency_outsideNL; outsideNL];
    SaccadeAmpoutsideL = [SaccadeAmpoutsideL; mean(tmpSaccadeAmpoutsideL)];
    SaccadeAmpoutsideNL = [SaccadeAmpoutsideNL; mean(tmpSaccadeAmpoutsideNL)];
    DistancefromtargetoutsideL = [DistancefromtargetoutsideL; mean(tmpDistancefromtargetoutsideL)];
    DistancefromtargetoutsideNL = [DistancefromtargetoutsideNL; mean(tmpDistancefromtargetoutsideNL)];
    
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
subplot(323); histogram(Saccadelatency_withinNL_mod(1:Saccadelatency_withinNLidx),0:0.025:0.5,'EdgeColor',[0.5 0.5 0.5],'FaceColor',[1 1 1],'Normalization','probability'); hold on; histogram(Saccadelatency_withinL_mod(1:Saccadelatency_withinLidx),0:0.025:0.5,'EdgeColor',[0 0.5 1.0],'Normalization','probability','DisplayStyle','stairs');
xlabel('Latencies (ms)'); ylabel('% saccades'); title('Latency within RF: Maui'); set(gca,'Tickdir','out','Xlim',[0 0.5],'XTick',0:0.1:0.5,'Ylim',[0 0.8],'YTick',0.:0.2:0.8); legend('control','laser'); axis square;
subplot(325); histogram(Saccadelatency_outsideNL_mod(1:Saccadelatency_outsideNLidx),0:0.025:0.5,'EdgeColor',[0.5 0.5 0.5],'FaceColor',[1 1 1],'Normalization','probability'); hold on; histogram(Saccadelatency_outsideL_mod(1:Saccadelatency_outsideNLidx),0:0.025:0.5,'EdgeColor',[0 0.5 1.0],'Normalization','probability','DisplayStyle','stairs');
xlabel('Latencies (ms)'); ylabel('% saccades');  title('Latency outside RF: Maui'); set(gca,'Tickdir','out','Xlim',[0 0.5],'XTick',0:0.1:0.5,'Ylim',[0 0.8],'YTick',0.:0.2:0.8); axis square;
subplot(324); histogram(Saccadelatency_withinNL_mod(Saccadelatency_withinNLidx+1:end),0:0.025:0.5,'EdgeColor',[0.5 0.5 0.5],'FaceColor',[1 1 1],'Normalization','probability'); hold on; histogram(Saccadelatency_withinL_mod(Saccadelatency_withinLidx+1:end),0:0.025:0.5,'EdgeColor',[0 0.5 1.0],'Normalization','probability','DisplayStyle','stairs');
xlabel('Latencies (ms)'); ylabel('% saccades'); title('Latency within RF: Apollo'); set(gca,'Tickdir','out','Xlim',[0 0.5],'XTick',0:0.1:0.5,'Ylim',[0 0.8],'YTick',0.:0.2:0.8); axis square;
subplot(326); histogram(Saccadelatency_outsideNL_mod(Saccadelatency_outsideNLidx+1:end),0:0.025:0.5,'EdgeColor',[0.5 0.5 0.5],'FaceColor',[1 1 1],'Normalization','probability'); hold on; histogram(Saccadelatency_outsideL_mod(Saccadelatency_outsideNLidx+1:end),0:0.025:0.5,'EdgeColor',[0 0.5 1.0],'Normalization','probability','DisplayStyle','stairs');
xlabel('Latencies (ms)'); ylabel('% saccades');  title('Latency outside RF: Apollo'); set(gca,'Tickdir','out','Xlim',[0 0.5],'XTick',0:0.1:0.5,'Ylim',[0 0.8],'YTick',0.:0.2:0.8); axis square;
plot_counter = plot_counter + 1;

% Figure for plotting the Actual saccade amplitidues 
figure(plot_counter);
subplot(221); plot(SaccadeAmpoutsideNL(idx),SaccadeAmpoutsideL(idx),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1]); hold on;
plot(SaccadeAmpwithinNL(idx),SaccadeAmpwithinL(idx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
xlabel('Saccade Amp: Control'); ylabel('Saccade Amp: Laser'); title('Maui'); set(gca,'Tickdir','out','Xlim',[0 10],'Ylim',[0 10],'XTick',0:5:10,'YTick',0:5:10); 
plot([0 10],[0 10],'k'); legend('outside RF','inside RF'); axis square;
subplot(222); plot(SaccadeAmpoutsideNL(~idx),SaccadeAmpoutsideL(~idx),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1]); hold on;
plot(SaccadeAmpwithinNL(~idx),SaccadeAmpwithinL(~idx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
xlabel('Saccade Amp: Control'); ylabel('Saccade Amp: Laser'); title('Apollo'); set(gca,'Tickdir','out','Xlim',[0 20],'Ylim',[0 20],'XTick',0:5:20,'YTick',0:5:20); plot([0 20],[0 20],'k'); axis square;
subplot(223); plot(DistancefromtargetoutsideNL(idx),DistancefromtargetoutsideL(idx),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1]); hold on;
plot(DistancefromtargetwithinNL(idx),DistancefromtargetwithinL(idx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
xlabel('Dist from target: Control'); ylabel('Dist from target: Laser'); title('Maui'); set(gca,'Tickdir','out','Xlim',[0 10],'Ylim',[0 10],'XTick',0:5:10,'YTick',0:5:10); 
axis square; set(gca,'XScale','log','YScale','log','Xlim',[0.3 10],'Ylim',[0.3 10],'XTick',[0.3 1 3 10],'YTick',[0.3 1 3 10]); plot([0.3 10],[0.3 10],'k'); legend('outside RF','inside RF');
subplot(224); plot(DistancefromtargetoutsideNL(~idx),DistancefromtargetoutsideL(~idx),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1]); hold on;
plot(DistancefromtargetwithinNL(~idx),DistancefromtargetwithinL(~idx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
xlabel('Dist from target: Control'); ylabel('Dist from target: Laser'); title('Apollo'); set(gca,'Tickdir','out','Xlim',[0 20],'Ylim',[0 20],'XTick',0:5:20,'YTick',0:5:20); 
set(gca,'XScale','log','YScale','log','Xlim',[0.3 20],'Ylim',[0.3 20],'XTick',[0.3 1 3 10 20],'YTick',[0.3 1 3 10 20]); hold on; plot([0.3 20],[0.3 20],'k'); axis square;
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

% Plotting points by session - as asked by the reviewers 
[monkey_ID_session,ia,ib] = unique(k(:,1:7),'rows'); 
new_idx = monkey_ID_session(:,1) == 'M';
DistancefromtargetoutsideNL_session = []; DistancefromtargetoutsideL_session = [];
DistancefromtargetwithinNL_session = []; DistancefromtargetwithinL_session = [];
for ii = 1:size(monkey_ID_session,1)
    % Outside the RF
    DistancefromtargetoutsideNL_session = [DistancefromtargetoutsideNL_session; mean(DistancefromtargetoutsideNL(ib==ii))];
    DistancefromtargetoutsideL_session = [DistancefromtargetoutsideL_session; mean(DistancefromtargetoutsideL(ib==ii))];
    
    % Within the RF
    DistancefromtargetwithinNL_session = [DistancefromtargetwithinNL_session; mean(DistancefromtargetwithinNL(ib==ii))];
    DistancefromtargetwithinL_session = [DistancefromtargetwithinL_session; mean(DistancefromtargetwithinL(ib==ii))];
end

figure(plot_counter); set(gcf,'Name','Same data as before but session wise');
subplot(221); plot(DistancefromtargetoutsideNL_session(new_idx),DistancefromtargetoutsideL_session(new_idx),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1]); hold on;
plot(DistancefromtargetwithinNL_session(new_idx),DistancefromtargetwithinL_session(new_idx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
xlabel('Dist from target: Control'); ylabel('Dist from target: Laser'); title('Maui'); set(gca,'Tickdir','out','Xlim',[0 10],'Ylim',[0 10],'XTick',0:5:10,'YTick',0:5:10); 
axis square; set(gca,'XScale','log','YScale','log','Xlim',[0.3 10],'Ylim',[0.3 10],'XTick',[0.3 1 3 10],'YTick',[0.3 1 3 10]); plot([0.3 10],[0.3 10],'k'); legend('outside RF','inside RF');
subplot(222); plot(DistancefromtargetoutsideNL_session(~new_idx),DistancefromtargetoutsideL_session(~new_idx),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1]); hold on;
plot(DistancefromtargetwithinNL_session(~new_idx),DistancefromtargetwithinL_session(~new_idx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
xlabel('Dist from target: Control'); ylabel('Dist from target: Laser'); title('Apollo'); set(gca,'Tickdir','out','Xlim',[0 20],'Ylim',[0 20],'XTick',0:5:20,'YTick',0:5:20); 
set(gca,'XScale','log','YScale','log','Xlim',[0.3 20],'Ylim',[0.3 20],'XTick',[0.3 1 3 10 20],'YTick',[0.3 1 3 10 20]); hold on; plot([0.3 20],[0.3 20],'k'); axis square;
subplot(223); plot(DistancefromtargetoutsideNL_session(new_idx),DistancefromtargetoutsideL_session(new_idx),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1]); hold on;
plot(DistancefromtargetwithinNL_session(new_idx),DistancefromtargetwithinL_session(new_idx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
xlabel('Dist from target: Control'); ylabel('Dist from target: Laser' ); title('Maui'); set(gca,'Tickdir','out','Xlim',[0 10],'Ylim',[0 10],'XTick',0:5:10,'YTick',0:5:10); 
axis square; set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[0 5 10],'YTick',[0 5 10]); plot([0 10],[0 10],'k'); legend('outside RF','inside RF');
subplot(224); plot(DistancefromtargetoutsideNL_session(~new_idx),DistancefromtargetoutsideL_session(~new_idx),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1]); hold on;
plot(DistancefromtargetwithinNL_session(~new_idx),DistancefromtargetwithinL_session(~new_idx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
xlabel('Dist from target: Control'); ylabel('Dist from target: Laser'); title('Apollo'); set(gca,'Tickdir','out','Xlim',[0 20],'Ylim',[0 20],'XTick',0:5:20,'YTick',0:5:20); 
set(gca,'Xlim',[0 20],'Ylim',[0 20],'XTick',[0 10 20],'YTick',[0 10 20]); hold on; plot([0 20],[0 20],'k'); axis square;
plot_counter = plot_counter + 1;

% Some additional stats
p17 = signrank(DistancefromtargetoutsideNL_session(new_idx),DistancefromtargetoutsideL_session(new_idx)); % Maui - outside 
p18 = signrank(DistancefromtargetwithinNL_session(new_idx),DistancefromtargetwithinL_session(new_idx)); % Maui - within
p19 = signrank(DistancefromtargetoutsideNL_session(~new_idx),DistancefromtargetoutsideL_session(~new_idx)); % Apollo - outside 
p20 = signrank(DistancefromtargetwithinNL_session(~new_idx),DistancefromtargetwithinL_session(~new_idx)); % Apollo - within

%% Continuation of the previous code - plotting the the saccade last point - saccde first point from individual files  
close all;
clearvars;
if ~exist('plot_counter')
    plot_counter = 1;
end
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filename = fetch(conn,'SELECT filename FROM SCstimcue');
training = fetch(conn,'SELECT training FROM SCstimcue');
mode = fetch(conn,'SELECT mode FROM SCstimcue');
close(conn);
ind = find(strcmp(training,'no') & strcmp(mode,'SCstimcue'));
L = ceil(sqrt(numel(ind)));
k = cell2mat(filename(ind));
monkey_ID = k(:,1);
idx = monkey_ID == 'M';

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
plotfigures = 1;
Catchtrials_exist = [];
Saccadelatency_catchL = []; Saccadelatency_catchNL = [];
SaccadeAmpwithinL = []; SaccadeAmpwithinNL = [];
SaccadeAmpoutsideL = []; SaccadeAmpoutsideNL = [];
DistancefromtargetwithinL = []; DistancefromtargetwithinNL = [];
DistancefromtargetoutsideL = []; DistancefromtargetoutsideNL = [];
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
    
    if norm(RFloc/10)<10
        plotlimit = 10;
    else
        plotlimit = 20;
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
    SaccadeAmpwithinL = [SaccadeAmpwithinL; norm(centroidwithinL)];
    SaccadeAmpwithinNL = [SaccadeAmpwithinNL; norm(centroidwithinNL)];
    DistancefromtargetwithinL = [DistancefromtargetwithinL; norm(centroidwithinL - RFloc/10)]; 
    DistancefromtargetwithinNL = [DistancefromtargetwithinNL; norm(centroidwithinNL - RFloc/10)];
    diffcentroidwithin = [diffcentroidwithin; sum((centroidwithinL - centroidwithinNL).^2)];
    Saccadelatency_withinL = [Saccadelatency_withinL; Saccadelatency(ib==withinRFidx & Llaser)];
    Saccadelatency_withinNL = [Saccadelatency_withinNL; Saccadelatency(ib==withinRFidx & ~Llaser)];
    
    centroidoutsideL = []; centroidoutsideNL = [];
    outsideL = []; outsideNL = [];
    tmpSaccadeAmpoutsideL = []; tmpSaccadeAmpoutsideNL = [];
    tmpDistancefromtargetoutsideL = []; tmpDistancefromtargetoutsideNL = [];
    for jj = 1:numel(outsideRFidx)
        centroidoutsideL = [centroidoutsideL; mean(SaccadeX(ib==outsideRFidx(jj) & Llaser)) mean(SaccadeY(ib==outsideRFidx(jj) & Llaser))];
        centroidoutsideNL = [centroidoutsideNL; mean(SaccadeX(ib==outsideRFidx(jj) & ~Llaser)) mean(SaccadeY(ib==outsideRFidx(jj) & ~Llaser))];
        outsideL = [outsideL; Saccadelatency(ib==outsideRFidx(jj) & Llaser)];
        outsideNL = [outsideNL; Saccadelatency(ib==outsideRFidx(jj) & ~Llaser)];
        tmpSaccadeAmpoutsideL = [tmpSaccadeAmpoutsideL; norm(centroidoutsideL(end,:))];
        tmpSaccadeAmpoutsideNL = [tmpSaccadeAmpoutsideNL; norm(centroidoutsideNL(end,:))];
        tmpDistancefromtargetoutsideL = [tmpDistancefromtargetoutsideL; norm(centroidoutsideL(end,:)-(C(outsideRFidx(jj),:)/10))];
        tmpDistancefromtargetoutsideNL = [tmpDistancefromtargetoutsideNL; norm(centroidoutsideNL(end,:)-(C(outsideRFidx(jj),:)/10))];
    end
    diffcentroidoutside = [diffcentroidoutside; mean(sum((centroidoutsideL - centroidoutsideNL).^2),2)];
    Saccadelatency_outsideL = [Saccadelatency_outsideL; outsideL];
    Saccadelatency_outsideNL = [Saccadelatency_outsideNL; outsideNL];
    SaccadeAmpoutsideL = [SaccadeAmpoutsideL; mean(tmpSaccadeAmpoutsideL)];
    SaccadeAmpoutsideNL = [SaccadeAmpoutsideNL; mean(tmpSaccadeAmpoutsideNL)];
    DistancefromtargetoutsideL = [DistancefromtargetoutsideL; mean(tmpDistancefromtargetoutsideL)];
    DistancefromtargetoutsideNL = [DistancefromtargetoutsideNL; mean(tmpDistancefromtargetoutsideNL)];
    
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
    
    if plotfigures & ~idx(ii)
        figure(plot_counter),subplot(5,4,count); plot(upsample(SaccadeX(~Llaser & ~Lcatchtrials),2),upsample(SaccadeY(~Llaser & ~Lcatchtrials),2),'Color',[0.5 0.5 0.5]); hold on;
        plot(upsample(SaccadeX(Llaser & ~Lcatchtrials),2),upsample(SaccadeY(Llaser & ~Lcatchtrials),2),'Color',[0 0.5 1.0]); 
        plot(C(C(:,1)~=0 & C(:,2)~=0,1)/10,C(C(:,1)~=0 & C(:,2)~=0,2)/10,'o','MarkerSize',8,'MarkerEdgeColor',[0 0 0],'Linewidth',1);
        plot(RFloc(1)/10,RFloc(2)/10,'o','MarkerSize',8,'MarkerEdgeColor',[1 0 0],'Linewidth',1);
        title(fileofinterest); set(gca,'Tickdir','out','Xlim',[-plotlimit plotlimit],'Ylim',[-plotlimit plotlimit],'XTick',-plotlimit:plotlimit:plotlimit,'YTick',-plotlimit:plotlimit:plotlimit); axis square; hold off;
        
        % Plotting the catch trials 
        figure(plot_counter+1); subplot(5,4,count);
        if any(~Llaser & Lcatchtrials)
            plot(upsample(SaccadeX(~Llaser & Lcatchtrials),2),upsample(SaccadeY(~Llaser & Lcatchtrials),2),'Color',[0.5 0.5 0.5]); hold on;
        end
        if any(Llaser & Lcatchtrials)
            plot(upsample(SaccadeX(Llaser & Lcatchtrials),2),upsample(SaccadeY(Llaser & Lcatchtrials),2),'Color',[0 0.5 1.0]);
        end
        set(gca,'Tickdir','out','Xlim',[-plotlimit plotlimit],'Ylim',[-plotlimit plotlimit],'XTick',-plotlimit:plotlimit:plotlimit,'YTick',-plotlimit:plotlimit:plotlimit); axis square; hold off;
        count = count + 1;
    end
end
plot_counter = plot_counter + 2;

%% This code is based on Greg's suggestion: 
% Ideally, the time window should start when the saccade starts (by some eye velocity criterion) and end when the saccade ends. 
% getSacData returns the start time and end times (relative to the "anlgStartTime", which is a column in stro.ras) of each saccade. 
% You can use this to pull out the appropriate eye position samples.

close all;
clearvars;
if ~exist('plot_counter')
    plot_counter = 1;
end
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filename = fetch(conn,'SELECT filename FROM SCstimcue');
training = fetch(conn,'SELECT training FROM SCstimcue');
mode = fetch(conn,'SELECT mode FROM SCstimcue');
close(conn);
ind = find(strcmp(training,'no') & strcmp(mode,'SCstimcue'));
L = ceil(sqrt(numel(ind)));
k = cell2mat(filename(ind));
monkey_ID = k(:,1);
idx = monkey_ID == 'M';

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
plotfigures = 1;
Catchtrials_exist = [];
Saccadelatency_catchL = []; Saccadelatency_catchNL = [];
SaccadeAmpwithinL = []; SaccadeAmpwithinNL = [];
SaccadeAmpoutsideL = []; SaccadeAmpoutsideNL = [];
DistancefromtargetwithinL = []; DistancefromtargetwithinNL = [];
DistancefromtargetoutsideL = []; DistancefromtargetoutsideNL = [];
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
    heyetraceind = strcmp(stro.sum.rasterCells(1,:),'AD11');
    veyetraceind = strcmp(stro.sum.rasterCells(1,:),'AD12');
    analogstrtimeind = strcmp(stro.sum.rasterCells(1,:),'anlgStartTime');
    [C,ia,ib] = unique(stro.trial(:,12:13),'rows');
    SPNL = zeros(size(C,1),1);
    SPL = zeros(size(C,1),1);
    samplerate = stro.sum.analog.storeRates{3};
    Llaser = ~isnan(laseron_t);
    
    if norm(RFloc/10)<10
        plotlimit = 10;
    else
        plotlimit = 20;
    end
    
    for jj = 1:length(SacData.amplitudes)

        fpofftime = fpoff_t(jj);
        analogstartime = stro.ras{jj,analogstrtimeind};
        t = stro.ras{jj,analogstrtimeind}+[0:1:length(stro.ras{jj,heyetraceind})-1]/samplerate;
        
        if ~isempty(find(SacData.starttimes{jj}-fpofftime>0 & SacData.starttimes{jj}-fpofftime<tmax))
            selected_times_idxs = find(SacData.starttimes{jj}-fpofftime>0 & SacData.starttimes{jj}-fpofftime<tmax);
            %             [~,ix] = max(SacData.amplitudes{jj}(selected_times_idxs)); % contingent on saccade amplitude
            ix = 1; % contingent on first saccade
            SaccadeStarttime = max(SacData.starttimes{jj}(selected_times_idxs(ix)));
            SaccadeEndTime = max(SacData.endtimes{jj}(selected_times_idxs(ix)));
        else
            SaccadeStarttime = nan;
            SaccadeEndTime = nan;
        end
        tinterest = t>SaccadeStarttime & t<SaccadeEndTime;
        eyex = stro.ras{jj,heyetraceind}(tinterest)*4096/400;
        eyey = stro.ras{jj,veyetraceind}(tinterest)*4096/400;
        
        if plotfigures & ~idx(ii) & ~Lcatchtrials(jj)
            if ~Llaser(jj)
                figure(plot_counter),subplot(5,4,count); plot(eyex,eyey,'Color',[0.5 0.5 0.5]); hold on;
            else
                figure(plot_counter),subplot(5,4,count); plot(eyex,eyey,'Color',[0 0.5 1.0]); hold on;
            end
        end
    end
    
    
    
    if plotfigures & ~idx(ii)
        figure(plot_counter),subplot(5,4,count); hold on;
        plot(C(C(:,1)~=0 & C(:,2)~=0,1)/10,C(C(:,1)~=0 & C(:,2)~=0,2)/10,'o','MarkerSize',8,'MarkerEdgeColor',[0 0 0],'Linewidth',1);
        plot(RFloc(1)/10,RFloc(2)/10,'o','MarkerSize',8,'MarkerEdgeColor',[1 0 0],'Linewidth',1);
        title(fileofinterest); set(gca,'Tickdir','out','Xlim',[-plotlimit plotlimit],'Ylim',[-plotlimit plotlimit],'XTick',-plotlimit:plotlimit:plotlimit,'YTick',-plotlimit:plotlimit:plotlimit); axis square; hold off;
        count = count + 1;
    end
end
plot_counter = plot_counter + 2;



%% Analysis of DToneloc Data
% Plotting the Psychometric data points based on how many times they were presented 
% Figure 5 of the paper

close all; clearvars;
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
    
    subplot(2,3,3*jj-2); h = bar([HitNL(jj)/stimpresentNL(jj) HitL(jj)/stimpresentL(jj); MissNL(jj)/stimpresentNL(jj) MissL(jj)/stimpresentL(jj); CRNL(jj)/stimabsentNL(jj) CRL(jj)/stimabsentL(jj); FANL(jj)/stimabsentNL(jj)  FAL(jj)/stimabsentL(jj)]); set(h(2),'FaceColor',[0 0.5 1]); set(h(1),'FaceColor',[0.5 0.5 0.5]); 
    set(gca,'XTick',[1 2 3 4],'XTickLabel',{'H','M','CR','FA'},'TickDir','out','Xlim',[0 5],'YTick',0:0.25:1.0,'Ylim',[0 1],'YTick',0:0.25:1.0); ylabel('Prop of trials'); axis square; title(filename(jj,1));
    subplot(2,3,3*jj-1); plot(LMScontrastfileNL,'-o','color',[0.5 0.5 0.5],'MarkerSize',3,'Linewidth',2); hold on; plot(LMScontrastfileL,'-o','color',[0 0.5 1.0],'MarkerSize',3,'Linewidth',2); xlabel('trial number'); ylabel('Luminance contrast');
    title(filename(jj,1)); axis square; set(gca,'TickDir','out','Xlim',[0 30],'XTick',[0:15:30],'Ylim',[0.1 1.0],'YTick',[0.1 0.3 1.0],'YScale','log'); hold off;
    
    contrastlattice = logspace(log10(0.03),log(1.0),51);
    fitL = gL*(1-exp(-((contrastlattice./aL).^bL)));
    fitNL = gNL*(1-exp(-((contrastlattice./aNL).^bNL)));
    [~,edges] = discretize(contrastNL,logspace(log10(min(contrastNL)-0.0001),log10(max(contrastL)+0.0001),9));
    for ss = 1:numel(edges)-1
        idx = contrastNL>=edges(ss) & contrastNL<edges(ss+1);
        if sum(idx)
            tot = sum(correctanswersNL(idx)+wronganswersNL(idx));
            subplot(2,3,3*jj);  plot(edges(ss),sum(correctanswersNL(idx))./(sum(correctanswersNL(idx)+wronganswersNL(idx))),'o','MarkerSize',3*sqrt(tot),'MarkerFacecolor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
        end
    end
    plot(contrastlattice,fitNL,'-','Linewidth',2,'color',[0.5 0.5 0.5]); hold on;
    for ss = 1:numel(edges)-1
        idx = contrastL>=edges(ss) & contrastL<edges(ss+1);
        if sum(idx)
            tot = sum(correctanswersL(idx)+wronganswersL(idx));
            subplot(2,3,3*jj);  plot(edges(ss),sum(correctanswersL(idx))./sum((correctanswersL(idx))+sum(wronganswersL(idx))),'o','MarkerSize',3*sqrt(tot),'MarkerFacecolor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]); hold on;
        end
    end
    plot(contrastlattice,fitL,'-','Linewidth',2,'color',[0 0.5 1.0]); set(gca,'Xlim',[0.03 1],'XTick',[0.03 0.1 0.3 1.0],'Tickdir','out','YTick',[0:0.25:1],'Xscale','log','XTick',[0.1 0.3 1],'XTickLabels', {'0.1','0.3','1'}); xlabel('Luminance contrast'); ylabel('Proportion correct'); title('Psychometric function'); axis square; hold off;
    
    [~,p] = equalproptest([HitL(jj) HitNL(jj)],[30 30],0.05);
    disp(p);
end
plot_counter = plot_counter + 1;

%% Retinotopic specificity of the laser stimulation
% Figure 6 in the paper 
close all; clearvars;
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
    subplot(3,3,ii+1); h = bar([sum(HitsNL(ind))/(sum(HitsNL(ind))+sum(MissNL(ind))) sum(HitsL(ind))/(sum(HitsL(ind))+sum(MissL(ind))); sum(MissNL(ind))/(sum(HitsNL(ind))+sum(MissNL(ind))) sum(MissL(ind))/(sum(HitsL(ind))+sum(MissL(ind))) ; sum(CRNL(ind))/(sum(CRNL(ind))+sum(FANL(ind))) sum(CRL(ind))/(sum(CRL(ind))+sum(FAL(ind)));  sum(FANL(ind))/(sum(CRNL(ind))+sum(FANL(ind))) sum(FAL(ind))/(sum(CRL(ind))+sum(FAL(ind)))]); hold on;
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
    
    [~,edges] = discretize(contrastNL,logspace(log10(min(contrastNL)-0.0001),log10(max(contrastL)+0.0001),9));
    figure(plot_counter)
    for ss = 1:numel(edges)-1
        idx = contrastNL>=edges(ss) & contrastNL<edges(ss+1);
        tot = sum(correctanswersNL(idx)+incorrectanswersNL(idx));
        if sum(idx) & tot
            
            subplot(3,3,ii+4); plot(edges(ss),sum(correctanswersNL(idx))./(sum(correctanswersNL(idx)+incorrectanswersNL(idx))),'o','MarkerSize',1.5*sqrt(tot),'MarkerFacecolor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
        end
    end
    for ss = 1:numel(edges)-1
        idx = contrastL>=edges(ss) & contrastL<edges(ss+1);
        tot = sum(correctanswersL(idx)+incorrectanswersL(idx));
        if sum(idx) & tot
            
            subplot(3,3,ii+4); plot(edges(ss),sum(correctanswersL(idx))./(sum(correctanswersL(idx)+incorrectanswersL(idx))),'o','MarkerSize',1.5*sqrt(tot),'MarkerFacecolor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]); hold on;
        end
    end
    plot(contrastlattice,fitL,'-','Linewidth',2,'color',[0 0.5 1.0]); plot(contrastlattice,fitNL,'-','Linewidth',2,'color',[0.5 0.5 0.5]);
    xlabel('contrast'); ylabel('Prop. correct'); set(gca,'Xlim',[0.03 1.0],'XTick',[0.03 0.1 0.3 1.0],'Tickdir','out','YTick',[0:0.25:1],'Xscale','log','XTick',[0.1 0.3 1],'XTickLabels', {'0.1','0.3','1'});
    title(locationstested(ii,:)'); axis square; hold off;
end

plot_counter = plot_counter + 1;
%% Population DToneloc plot
close all; clearvars;
if ~exist('plot_counter')
    plot_counter = 1;
end
load('T_vos1978_Y');
Vlambda = T_vos1978_Y';
ratio_controloverlaser = [];
timeperblock = [];
trialsperblock = [];
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
            
        end
    end
    % Accumulating files over different sessions
    filedates = cell2mat(newfilename);
    datesind = 2:7;
    filedates = str2num(filedates(:,datesind));
%     [filedates,~,ib] = unique([filedates RF newRWmultiplier newlaserdial newstimsize],'rows');
    [SessionID,~,ib] = unique([filedates],'rows'); % This is the line where I change the definition of session
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
    subplot(2,2,mm);
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
    figure(plot_counter); subplot(2,2,mm+2); plot(dprimeNL,dprimeL,'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on
    set(gca,'Ylim',[-1 3],'Xlim',[-1 3],'Tickdir','out','XTick',[-1:1:3],'YTick',[-1:1:3]); line([-1 3],[-1 3],'Linewidth',1);axis square; xlabel('d prime:control'); ylabel('d prime:laser'); hold off;
  
end
figure(plot_counter); set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;
