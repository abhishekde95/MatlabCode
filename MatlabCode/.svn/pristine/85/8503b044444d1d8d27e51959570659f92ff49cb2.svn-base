% Some additional info for Elife
% Author - Abhishek De, 10/19

%%  Figure 6: Storing the data from DToneloc experiment
% •	DToneloc population: A scatter plot of hits in control vs laser trials at tested locations for Maui and Apollo
% (Total 4 figures: 2 hits control vs laser population plot, 2 mapped RF locations for those session)
close all; clearvars;
if ~exist('plot_counter')
    plot_counter = 1;
end
load('T_vos1978_Y');
Vlambda = T_vos1978_Y';
ratio_controloverlaser = [];
timeperblock = [];
trialsperblock = [];
T = []; % for storing data in a table
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
    
    % Trying to create a table 
    Hit_laser = HitLSession;
    Miss_laser = MissLSession;
    False_alarm_laser = FALSession;
    Correct_reject_laser = CRLSession;
    Gabor_absent_laser = stimabsentLSession;
    Gabor_present_laser = stimpresentLSession;
    
    Hit_control = HitNLSession;
    Miss_control = MissNLSession;
    False_alarm_control = FANLSession;
    Correct_reject_control = CRNLSession;
    Gabor_absent_control = stimabsentNLSession;
    Gabor_present_control = stimpresentNLSession;

    tmp = char(filename);
    Monkey_ID = tmp(:,1);
    Dates = strcat(tmp(:,2:3),'/',tmp(:,4:5),'/',tmp(:,6:7));
    T = [T; table(Monkey_ID, Dates,Hit_laser,Miss_laser,False_alarm_laser,Correct_reject_laser,Gabor_present_laser,Gabor_absent_laser,Hit_control,Miss_control,False_alarm_control,Correct_reject_control,Gabor_present_control,Gabor_absent_control)];
    
end
figure(plot_counter); set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;

% Storing the files in an excel sheet format
excelfilename = 'DLX_OPTO_data_Gabor2AFC.xlsx';
writetable(T,excelfilename,'Sheet',1);

%% Storing the neurophysiology data
% Figure 2
% •	Neurophysiology: One example each of suppression and activation rasters, singles-units along with their waveforms in inset (2 figures)
% •	Scatter plot of all suppression and activation units, classified as single and multi-units (1 figure)
% •	One example of a single-unit from Fixstim interfreq and Gratings displaying a putative ‘direction selective’ inhibitory neuron (2 figures)
close all; clearvars;
if ~exist('plot_counter')
    plot_counter = 1;
end
binwidth = .005; % 5 ms bin width
bins = -0.4:binwidth:0.6;
baselinebins = (bins>-0.3 & bins<0);

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
spiketimes_cells = [];
baselineFR = [];
laserFR = [];
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
    spiketimes_tmp = [];
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
            tmp_spikecountsbaseline = [tmp_spikecountsbaseline; sum(spiketimes>laserontime-0.3 & spiketimes<laserontime)];
            tmp_spikecountslaserstimabsent = [tmp_spikecountslaserstimabsent; sum(spiketimes>laserontime & spiketimes<laserofftime)];
            spiketimes_tmp = [spiketimes_tmp; {spiketimes - laserontime}];
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
    spiketimes_cells = [spiketimes_cells; {spiketimes_tmp}];
end
tmp = char(filename);
Monkey_ID = tmp(:,1);
T = table(Monkey_ID, spiketimes_cells, baselineFR, laserFR);
save DLX_OPTO_neurophysiology T

%% Storing the SC stimcue data 
% Adapted from Supplementary figure 2: SC stimcue population 
close all; clearvars;
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
files_for_storing = [];
num_targets_for_storing = [];
Saccadelatency_within_laser_for_storing = [];
Saccadelatency_within_control_for_storing = [];
Saccadelatency_outside_laser_for_storing = [];
Saccadelatency_outside_control_for_storing = [];
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
    
    % Storing data in variables so that I can later organize them in a table
    files_for_storing = [files_for_storing; fileofinterest];
    num_targets_for_storing = [num_targets_for_storing; numel([withinRFidx; outsideRFidx])]; 
    Saccadelatency_within_laser_for_storing = [Saccadelatency_within_laser_for_storing; {Saccadelatency(ib==withinRFidx & Llaser)}]; 
    Saccadelatency_within_control_for_storing = [Saccadelatency_within_control_for_storing; {Saccadelatency(ib==withinRFidx & ~Llaser)}]; 
    Saccadelatency_outside_laser_for_storing = [Saccadelatency_outside_laser_for_storing; {outsideL}];
    Saccadelatency_outside_control_for_storing = [Saccadelatency_outside_control_for_storing; {outsideNL}];
    
end
plot_counter = plot_counter + 1;

Monkey_ID = files_for_storing(:,1);
Dates = strcat(files_for_storing(:,2:3),'/',files_for_storing(:,4:5),'/',files_for_storing(:,6:7));
Saccade_latency_outside_laser = Saccadelatency_outside_laser_for_storing;
Saccade_latency_outside_control = Saccadelatency_outside_control_for_storing;
Saccade_latency_within_laser = Saccadelatency_within_laser_for_storing;
Saccade_latency_within_control = Saccadelatency_within_control_for_storing;
Mean_diff_in_saccade_end_pts_btw_laser_and_control_within = diffcentroidwithin/10;
Mean_diff_in_saccade_end_pts_btw_laser_and_control_outside = diffcentroidoutside/10;
Number_of_targets = num_targets_for_storing;

T = table(Monkey_ID,Dates,Number_of_targets,Saccade_latency_outside_laser,Saccade_latency_outside_control,Saccade_latency_within_laser,Saccade_latency_within_laser,Mean_diff_in_saccade_end_pts_btw_laser_and_control_within,Mean_diff_in_saccade_end_pts_btw_laser_and_control_outside); 
save DLX_Saccade_task_data T 

