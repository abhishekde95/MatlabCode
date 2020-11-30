% Population figure for the editor
% Author - Abhishek De, 5/12
% •	Example of excitation along with waveform in inset
% •	Example of suppression along with waveform in inset
% •	1 example SC_stimcue session from Maui
% •	1 DToneloc file with hits, misses, CR and FA, along with psychometric function in the inset
% •	1 population figure from Apollo of hits in laser vs control trials
% •	Maybe a histology figure

%% 
clearvars; close all;
if ~exist('plot_counter')
    plot_counter = 1;
end
filename = {'A012319005.nex';'A020119004.nex'};
% filename = {'A012219011.nex';'A020119004.nex'};
titles = {'activation';'suppression'};
N = numel(filename); 
L = ceil(sqrt(N));
binwidth = .005;
bins = -0.2:binwidth:0.5;
baselinebins = (bins>-0.3 & bins<0);
baselineFR = [];
laserFR = [];
statstestresult = []; % Performing a wilcoxon rank sum test and storing the p value
figure(1);
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
            subplot(2,4,2*jj-1); plot(spiketimes(spiketimes>laserontime-0.35 & spiketimes<laserofftime+0.25)-laserontime,0.5*count1*ones(size(spiketimes(spiketimes>laserontime-0.35 & spiketimes<laserofftime+0.25)))','k.'); hold on;
            PSTHlaser = PSTHlaser + hist(spiketimes-laserontime, bins)/binwidth;
            count1 = count1 + 1;
            laserbins = bins>=0 & bins<=timedurlaser;

        else 
        end
    end
    PSTHlaser = PSTHlaser/(count1-1);
    line([0 0],[0 0.5*count1]); line([timedurlaser timedurlaser],[0 0.5*count1]); 
    set(gca,'Xlim',[-0.15 0.45],'Ylim',[0 0.5*count1],'Tickdir','out','XTick',[-0.15:0.15:0.45]); xlabel('time(s)'); ylabel('Number of trials'); axis square; title(titles{jj}); 
    
    subplot(2,4,2*jj); plot(bins,PSTHlaser,'color',[0 0.5 1],'Linewidth',2); 
    set(gca,'Xlim',[-0.15 0.45],'Tickdir','out','XTick',[-0.15:0.15:0.45],'YTick',[0:100:600]); 
    if max(PSTHlaser(2:end-1))<20
        set(gca,'Ylim',[0 20]);
    end
    xlabel('time(s)'); ylabel('Response (ips)'); axis square; title(titles{jj});     
end

% SC_stimcue example
stro = nex2stro(findfile('M042718005.nex'));
targ_x = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_x'));
targ_y = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_y'));
uniquetargxy = unique([targ_x targ_y],'rows');
ntrials = size(stro.trial,1);
Lcatchtrials = targ_x == 0 & targ_y == 0;
laseron_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
laseroff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
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
subplot(2,4,5); plot(stro.sum.exptParams.rf_x/10,stro.sum.exptParams.rf_y/10,'o','MarkerSize',15,'MarkerEdgeColor',[0 0 0],'Linewidth',3); hold on;
RFloc = [stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y];
if norm(RFloc)>100
    plotlimit = 18;
else
    plotlimit = 8;
end
for i = 1:size(stro.trial,1)
    x = stro.ras{i,2}*4096/400;
    y = stro.ras{i,3}*4096/400;
    t = stro.ras{i,5}+[0:1:length(x)-1]/samplerate;
    Lt = t>fpoff_t(i) & t < fpoff_t(i) +.3;
    if ~Lcatchtrials(i)
        if ~Llaser(i)
            subplot(2,4,5); plot(x(Lt),y(Lt),'color',[0.5 0.5 0.5],'Linewidth',1); hold on;
        else
            subplot(2,4,5); plot(x(Lt),y(Lt),'color',[0 0.5 1],'Linewidth',1); hold on;
        end
    else
    end
end
subplot(2,4,5); axis square; set(gca,'Xlim',[-1*plotlimit plotlimit],'Ylim',[-1*plotlimit plotlimit],'Tickdir','out','XTick',[-8:4:8],'YTick',[-8:4:8]); xlabel('Degrees'); ylabel('Degrees'); title('SCstimcue:Maui');

% DToneloc example file from 
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
RWmultiplier = 0.6;
laserdial = 1.8;
oogvals = [];
stro = nex2stro(findfile('A021419009.nex'));
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

fundamentals = stro.sum.exptParams.fundamentals; % CONE FUNDAMENTALS: L,M,S
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
mon_spd = stro.sum.exptParams.mon_spd; % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals
Lb = M*stro.sum.exptParams.bkgndrgb;
% Converting all the contrasts to Luminance contrast
Bkgnd_luminance = stro.sum.exptParams.bkgndrgb'*mon_spd*Vlambda;
Contrast_luminance = ((inv(M)*((LMStriplet+1).*(repmat(Lb',[120 1])))')'*mon_spd*Vlambda-Bkgnd_luminance)/Bkgnd_luminance;
LMStripletfile = Contrast_luminance;
LMScontrastfileNL = sqrt(sum(LMStripletfile(colorstimpresent & ~lasertrialidxs,:).^2,2));
LMScontrastfileL = sqrt(sum(LMStripletfile(colorstimpresent & lasertrialidxs,:).^2,2));
% calculating gamut edge luminance contrast
t = min(((1./stro.sum.exptParams.bkgndrgb)-1));
OOG_luminance = (stro.sum.exptParams.bkgndrgb*(1+t))' * mon_spd*Vlambda;
putativeoogcontrast =  (OOG_luminance - Bkgnd_luminance)/Bkgnd_luminance;
oogvals = [oogvals; putativeoogcontrast];

% Calculating psychophysical detection threshold for non-laser trials
trials = lasertrials;
% only evaluating the laser trials
weibullparamsL = [];
weibullparamsNL = [];
for kk = 1:2
    if kk == 2
        trials = ~trials;
    end
    LMScontrast = [];
    LMScontrast = sqrt(sum(Contrast_luminance(stimpresent & trials,:).^2,2));
    answers = correcttrials(stimpresent & trials);
    LMScontrast(LMScontrast>putativeoogcontrast) = putativeoogcontrast;
    
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
subplot(2,4,6);h = bar([HitNL HitL; MissNL MissL]); set(h(2),'FaceColor',[0 0.5 1]); set(h(1),'FaceColor',[0.5 0.5 0.5]); ylabel('Number of trials');
set(gca,'XTick',[1 2],'XTickLabel',{'Hits','Misses'},'TickDir','out','Xlim',[0 3],'YTick',[0:10:30]); axis square; title('DToneloc'); hold off;
contrastlattice = logspace(log10(0.03),log10(1.0),51);
fitL = gL*(1-exp(-((contrastlattice./aL).^bL)));
fitNL = gNL*(1-exp(-((contrastlattice./aNL).^bNL)));
[~,edges] = discretize(contrastNL,logspace(log10(min(contrastNL)-0.0001),log10(max(contrastL)+0.0001),7));
for jj = 1:numel(edges)-1
    idx = contrastNL>=edges(jj) & contrastNL<edges(jj+1);
    if sum(idx)
        subplot(2,4,7);  plot(edges(jj),sum(correctanswersNL(idx))./(sum(correctanswersNL(idx))+sum(wronganswersNL(idx))),'o','MarkerSize',8,'MarkerFacecolor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
    end
end
plot(contrastlattice,fitNL,'-','Linewidth',2,'color',[0.5 0.5 0.5]); hold on;
for jj = 1:numel(edges)-1
    idx = contrastL>=edges(jj) & contrastL<edges(jj+1);
    if sum(idx)
        subplot(2,4,7);  plot(edges(jj),sum(correctanswersL(idx))./sum((correctanswersL(idx))+sum(wronganswersL(idx))),'o','MarkerSize',8,'MarkerFacecolor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]); hold on;
    end
end
plot(contrastlattice,fitL,'-','Linewidth',2,'color',[0 0.5 1.0]); set(gca,'Xlim',[0.03 1.0],'Tickdir','out','YTick',[0:0.25:1],'Xscale','log','XTick',[0.03 0.1 0.3 1],'XTickLabels', {'0.03','0.1','0.3','1'}); xlabel('Luminance contrast'); ylabel('Proportion correct'); title('Psychometric function'); axis square; hold off;

%  Population figure for DToneloc from Apollo   
load filenameoptoA.mat
load RWmultiplieroptoA.mat
load laserdialoptoA.mat
monkeyname = 'Apollo';
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
        
        tmpL = []; tmpNL = [];
        for kk = 1:numdivisionsperfile
            tmpL = [tmpL sum(colordirchoiceidxs(CRlasertrialidxs((kk-1)*(L1/numdivisionsperfile)+1:kk*L1/numdivisionsperfile)))];
            tmpNL = [tmpNL sum(colordirchoiceidxs(CRnonlasertrialidxs((kk-1)*(L2/numdivisionsperfile)+1:kk*L2/numdivisionsperfile)))];
        end
        CRearlylateL = [CRearlylateL; tmpL];
        CRearlylateNL = [CRearlylateNL; tmpNL];
        
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
        Contrast_luminance = ((inv(M)*((LMStriplet+1).*(repmat(Lb',[size(LMStriplet,1) 1])))')'*mon_spd*Vlambda-Bkgnd_luminance)/Bkgnd_luminance;
        LMStripletfile = Contrast_luminance;
        LMScontrastfileNL = sqrt(sum(LMStripletfile(colorstimpresent & ~lasertrialidxs,:).^2,2));
        LMScontrastfileL = sqrt(sum(LMStripletfile(colorstimpresent & lasertrialidxs,:).^2,2));
        % calculating gamut edge luminance contrast
        t = min(((1./stro.sum.exptParams.bkgndrgb)-1));
        OOG_luminance = (stro.sum.exptParams.bkgndrgb*(1+t))' * mon_spd*Vlambda;
        putativeoogcontrast =  (OOG_luminance - Bkgnd_luminance)/Bkgnd_luminance;
        LMStripletfile = Contrast_luminance(ind,:);
        LMScontrastfileNL = sqrt(sum(LMStripletfile(colorstimpresent & ~lasertrialidxs,:).^2,2));
        LMScontrastfileL = sqrt(sum(LMStripletfile(colorstimpresent & lasertrialidxs,:).^2,2));
        
        answers = colordirchoiceidxs(colorstimpresent & ~lasertrialidxs);
        
        LMScontrastfileNL(LMScontrastfileNL>putativeoogcontrast) = putativeoogcontrast;
        LMScontrastfileL(LMScontrastfileL>putativeoogcontrast) = putativeoogcontrast;

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
plot_counter = plot_counter + 1;

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

subplot(2,4,8); color = [0 0 0];
for ii = 1:numel(p_hits)
    px = HitNLSession(ii)./stimpresentNLSession(ii); errorx = sqrt(px*(1-px)/stimpresentNLSession(ii));
    py = HitLSession(ii)./stimpresentLSession(ii); errory = sqrt(py*(1-py)/stimpresentLSession(ii));
    plot(px,py,'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',color,'MarkerEdgeColor',[1 1 1],'Color','k'); hold on;
end
line([0 1],[0 1],'Linewidth',1);axis square; xlabel('Proportion of hits (control)'); ylabel('Proportion of hits (laser)'); title('Pop DToneloc: Apollo'); set(gca,'Xlim',[0 1],'Ylim',[0 1],'Tickdir','out','XTick',0:0.25:1,'YTick',0:0.25:1);


