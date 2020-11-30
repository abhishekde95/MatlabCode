% For Yasmine's R01 
% Seeing the effect of laser power on neural excitation and suppression
% Author - Abhishek De, 1/20

close all; clearvars;
if ~exist('plot_counter')
    plot_counter = 1;
end

filename_M = {'M051318008.nex';'M051318009.nex';'M051618004.nex';'M051618005.nex';'M051618006.nex';'M051618007.nex';'M051618008.nex';'M051618010.nex';...
    'M051618011.nex';'M051618012.nex';'M051618013.nex';'M051618014.nex';'M051818004.nex';'M051818005.nex';'M051818006.nex';'M051818007.nex';...
    'M051818010.nex';'M051818011.nex';'M051818012.nex';'M051818013.nex';'M051818014.nex';'M052318001.nex';'M052318002.nex';'M052318003.nex';'M052318004.nex'};
laserdial_M = [1.0; 0.7; 0.7; 0.5; 0.4; 0.3; 0.45; 0.6;...
    0.6; 1.0; 4.0; 4.0; 2.0; 2.0; 2.0; 2.0;...
    2.0; 2.0; 2.0; 2.0; 1.0; 2.0; 2.0; 2.0; 2.0];
filename_A = {'A011719009.nex';'A011719010.nex';'A011719011.nex';'A011719012.nex';'A012219006.nex';'A012219009.nex';'A012219010.nex';'A012219011.nex';'A012219012.nex';...
    'A012319004.nex';'A012319005.nex';'A012319006.nex';'A013119002.nex';'A013119003.nex';'A013119004.nex';'A013119005.nex';'A013119006.nex';'A020119004.nex';'A020119005.nex';...
    'A020319001.nex';'A020319007.nex';'A020319008.nex';'A020319009.nex';'A020519003.nex';'A020519008.nex';'A020519009.nex';'A021519004.nex';'A021519005.nex';'A021519006.nex';...
    'A021519007.nex';'A021519008.nex'};
laserdial_A = [2.4; 2.4; 2.4; 2.4; 2.5; 2.5; 2.7; 2.7; 2.7;...
    2.0; 2.0; 1.3; 1.0; 1.5; 1.5; 1.0; 1.0; 2.0; 2.0;...
    1.3; 2.0; 2.5; 2.0; 1.5; 1.3; 1.5; 1.0; 0.5; 1.5;...
    0.75; 2.0];
filename = [filename_M; filename_A];
laserdial = [laserdial_M; laserdial_A];

load laserdetails
laserdial = spline(laserdetails.dial,laserdetails.laserpower,laserdial);

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
subplot(121); plot(baselineFRlog(~Singleunitidxs & Maui_idxs),laserFRlog(~Singleunitidxs & Maui_idxs),'s','MarkerSize',8,'LineWidth',1.0,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(baselineFRlog(Singleunitidxs & Maui_idxs),laserFRlog(Singleunitidxs & Maui_idxs),'s','MarkerSize',8,'LineWidth',1.0,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(baselineFRlog(~Singleunitidxs & Apollo_idxs),laserFRlog(~Singleunitidxs & Apollo_idxs),'o','MarkerSize',6,'LineWidth',1.0,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(baselineFRlog(Singleunitidxs & Apollo_idxs),laserFRlog(Singleunitidxs & Apollo_idxs),'o','MarkerSize',6,'LineWidth',1.0,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
axis square; xlabel('baselineFR'); ylabel('laserFR'); set(gca,'Xlim',[0.1 300],'Ylim',[0.1 300],'TickDir','out','XScale','log','YScale','log'); line([0.1 300],[0.1 300],'Color','k'); 
legend('M:multi','A:multi','A:single'); title('Opto stim modulation'); hold off;
subplot(122); 
for ii = 1:numel(laserdial)
     plot(baselineFRlog(ii),laserFRlog(ii),'o','MarkerSize',laserdial(ii)/15,'LineWidth',1.0,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
end
axis square; xlabel('baselineFR'); ylabel('laserFR'); set(gca,'Xlim',[0.1 300],'Ylim',[0.1 300],'TickDir','out','XScale','log','YScale','log'); line([0.1 300],[0.1 300],'Color','k'); 
title('Opto stim modulation = f(Laser power)'); hold off;

[r1,p1] = corr(laserdial(laserFR>baselineFR),laserFR(laserFR>baselineFR)-baselineFR(laserFR>baselineFR),'type','Spearman'); % excited sites 
[r2,p2] = corr(laserdial(laserFR<baselineFR),laserFR(laserFR<baselineFR)-baselineFR(laserFR<baselineFR),'type','Spearman'); % suppressed sites 

%% 
close all; clearvars
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
subplot(323);plot(linspace(0,max(laserdial),51),fit,'color',[0 0 0],'Linewidth',2); xlabel('Laser power (mW)'); ylabel('Hits control - Hits laser'); 
title('Impact on behavior'); set(gca,'Ylim',[0 1],'Xlim',[0 100],'XTick',[0:20:100],'TickDir','out'); axis square; hold off;
subplot(324);plot(linspace(0,max(laserdial),51),fit2,'color',[0 0 0],'Linewidth',2); xlabel('Laser power (mW)'); ylabel('dprime control-laser'); 
title('Impact on behavior'); set(gca,'Ylim',[0 3],'Xlim',[0 100],'XTick',[0:20:100],'TickDir','out'); axis square; hold off;
[~,p] = corr(laserdial,linspace(0,numel(laserdial)-1,numel(laserdial))');
plot_counter = plot_counter + 1;

% Checking out the effect of laser power on neuronal modulation
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
            laserontime = stro.trial(ind,10);
            laserofftime =  stro.trial(ind,11);
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

subplot(325); plot(laserdial,(laserFR-baselineFR),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); axis square;
xlabel('Laser power (mW)'); ylabel('Laser FR-Baseline FR'); set(gca,'Tickdir','out'); title('Impact on FR');
subplot(326); plot(laserdial,(laserFR-baselineFR)./(laserFR+baselineFR),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); axis square;
xlabel('Laser power (mW)'); ylabel('Neuronal modulation'); set(gca,'Tickdir','out','YScale','log'); title('Impact on FR');
plot_counter = plot_counter + 1;
