% Files for cranking out figures for the meeting
% Author - Abhishek De, 3/19
% Here's the list of the figures I intend to print
% Section 1: DToneloc
% 1.1 �	Example single files from Maui and Apollo with their staircases in laser and control files
% 1.2 �	Population data � differences in hits (and CRs) in control and laser trials from Maui and Apollo
% 1.3 �	Effect of laser on Hits in laser and control trials
% 1.4 �	Spatial extent of the effect
% 1.4.1 o	Example with randomly interleaved stim locations
% 1.4.2 o	Example with files tested for different locations where stim loc was fixed for each file
% 1.4.3 o Examples of laser effects on performance as function of time(plot the psychometric function for each of those files and see if there is any transition over time)
% Section 2: SCstimcue
% 2.1 �	Example file of SCstimcue from Maui and Apollo
% 2.2 �	Example file of SCstimcue in each animal during catch trials
% 2.3 �	A population data of probability of making a saccade into the target location in laser and control trial for both Maui and Apollo
% 2.4 �	Maybe a plot about latency differences between control and laser trials
% 2.5 �	A control file where the target locations were the same as DToneloc saccade target locations
% Section 3: Neuronal data
% 3.1 �	One example file of excitation and suppression
% 3.2 �	A population figure of excitation and suppression
% 3.3 �	�	1 example file (each monkey) from FixStim interfreq to show sinusoidal modulation. Included a grating file M032518016.nex. 
close all; clearvars;

%% Section 1.1 Example single files from Maui and Apollo with their staircases in laser and control files
if ~exist('plot_counter')
    plot_counter = 1;
end
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
    LMStripletfile = LMStriplet;
    LMScontrastfileNL = sqrt(sum(LMStripletfile(colorstimpresent & ~lasertrialidxs,:).^2,2));
    LMScontrastfileL = sqrt(sum(LMStripletfile(colorstimpresent & lasertrialidxs,:).^2,2));
    putativeoogcontrast = min(LMScontrastfileL(logical(stro.trial(colorstimpresent & lasertrialidxs,oog))));
    
    answers = colordirchoiceidxs(colorstimpresent & ~lasertrialidxs);
    if ~isempty(putativeoogcontrast)
        LMScontrastfileNL(LMScontrastfileNL>putativeoogcontrast) = putativeoogcontrast;
        LMScontrastfileL(LMScontrastfileL>putativeoogcontrast) = putativeoogcontrast;
    end
    subplot(2,2,2*jj-1);bar([HitL(jj) HitNL(jj); MissL(jj) MissNL(jj); CRL(jj) CRNL(jj); FAL(jj) FANL(jj)]); ylabel('Count');
    set(gca,'XTick',[1 2 3 4],'XTickLabel',{'H','M','CR','FA'},'TickDir','out','Xlim',[0 5]); axis square; title(filename(jj,1));
    subplot(2,2,2*jj); plot(LMScontrastfileNL,'color',[0.5 0.5 0.5],'Linewidth',2); hold on; plot(LMScontrastfileL,'color',[0 0.5 1.0],'Linewidth',2); xlabel('trial number'); ylabel('contrast');
    title(filename(jj,1)); axis square; set(gca,'TickDir','out'); hold off;
end
plot_counter = plot_counter + 1;

%% Section 1.2 Population data � differences in hits (and CRs) in control and laser trials from Maui and Apollo
% Run this script separately for Maui and Apollo
if ~exist('plot_counter')
    plot_counter = 1;
end
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
            LMStripletfile = LMStriplet(ind,:);
            LMScontrastfileNL = sqrt(sum(LMStripletfile(colorstimpresent & ~lasertrialidxs,:).^2,2));
            LMScontrastfileL = sqrt(sum(LMStripletfile(colorstimpresent & lasertrialidxs,:).^2,2));
            putativeoogcontrast = min(LMScontrastfileL(logical(stro.trial(colorstimpresent & lasertrialidxs,oog))));
            
            answers = colordirchoiceidxs(colorstimpresent & ~lasertrialidxs);
            if ~isempty(putativeoogcontrast)
                LMScontrastfileNL(LMScontrastfileNL>putativeoogcontrast) = putativeoogcontrast;
                LMScontrastfileL(LMScontrastfileL>putativeoogcontrast) = putativeoogcontrast;
            end
        end
    end
    plot_counter = plot_counter + 1;
    
    pL = anova1(CRearlylateL,[],'off');
    pNL = anova1(CRearlylateNL,[],'off');
%     figure(plot_counter); set(gcf,'Name',strcat('Effect: early vs late trials:',monkeyname));
%     subplot(221); boxplot(CRearlylateL,'Notch','on'); hold on; text(1,3,strcat('p=',num2str(pL,2)));axis square; ylabel('Correct Rejects'); set(gca,'Tickdir','out'); title('Laser-ANOVA');
%     subplot(222); boxplot(CRearlylateNL,'Notch','on'); hold on; text(1,3,strcat('p=',num2str(pNL,2))); axis square; ylabel('Correct Rejects'); set(gca,'Tickdir','out'); title('Control-ANOVA');
%     xvals = repmat(1:1:numdivisionsperfile,[size(CRearlylateL,1) 1]); xvals = xvals(:);
%     [r1,p1] = corr(xvals,CRearlylateL(:)); [r2,p2] = corr(xvals,CRearlylateNL(:));
%     subplot(223);plot(xvals,CRearlylateL(:),'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on; lsline;
%     text(1.5,25,strcat('r=',num2str(r1,2))); text(1.5,20,strcat('p=',num2str(p1,2)));
%     axis square; xlabel('intervals'); ylabel('CR'); title('Laser-Regression'); set(gca,'Ylim',[0 30],'Tickdir','out'); hold off;
%     subplot(224);plot(xvals,CRearlylateNL(:),'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on; lsline;
%     text(1.5,25,strcat('r=',num2str(r2,2))); text(1.5,20,strcat('p=',num2str(p2,2)));
%     axis square; xlabel('intervals'); ylabel('CR'); title('Control-Regression'); set(gca,'Ylim',[0 30],'Tickdir','out'); hold off;
    
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
    
    figure(plot_counter); set(gcf,'Name',strcat('stats:',monkeyname)); subplot(221);
%     for ii = 1:numel(p_FA)
%         if (p_FA(ii)<0.05)
%             color = [0 0 0];
%         else
%             color = [0.5 0.5 0.5];
%         end
%         px = CRNLSession(ii)/stimabsentNLSession(ii); errorx = sqrt(px*(1-px)/stimabsentNLSession(ii));
%         py = CRLSession(ii)/stimabsentLSession(ii); errory = sqrt(py*(1-py)/stimabsentLSession(ii));
%         errorbar(px,py,errorx,errorx,errory,errory,'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',color,'MarkerEdgeColor',[1 1 1],'Color','k'); hold on;
%     end
%     line([0 1],[0 1],'Linewidth',1);axis square; xlabel('CR control'); ylabel('CR laser'); set(gca,'Xlim',[0 1],'Ylim',[0 1],'Tickdir','out');
    subplot(121);
    for ii = 1:numel(p_hits)
        if (p_hits(ii)<0.05)
            color = [0 0 0];
        else
            color = [0.5 0.5 0.5];
        end
        px = HitNLSession(ii)./stimpresentNLSession(ii); errorx = sqrt(px*(1-px)/stimpresentNLSession(ii));
        py = HitLSession(ii)./stimpresentLSession(ii); errory = sqrt(py*(1-py)/stimpresentLSession(ii));
        errorbar(px,py,errorx,errorx,errory,errory,'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',color,'MarkerEdgeColor',[1 1 1],'Color','k'); hold on;
    end
    line([0 1],[0 1],'Linewidth',1);axis square; xlabel('Hits control'); ylabel('Hits laser'); set(gca,'Xlim',[0 1],'Ylim',[0 1],'Tickdir','out');
    [~,p] = corr((HitNLSession./stimpresentNLSession)-(HitLSession./stimpresentLSession),(CRNLSession./stimabsentNLSession)-(CRLSession./stimabsentLSession));
%     subplot(223),plot((HitNLSession./stimpresentNLSession)-(HitLSession./stimpresentLSession),(CRNLSession./stimabsentNLSession)-(CRLSession./stimabsentLSession),'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1],'Color','k');
%     xlabel('Hits control - Hit laser'); ylabel('CR control - CR laser'); lsline; title(strcat('p=',num2str(p,2))); axis square;
    subplot(122); plot(stro.sum.exptParams.fp_x,stro.sum.exptParams.fp_y,'+','Markersize',10,'Linewidth',2); hold on; plot(RF(:,1),RF(:,2),'o','MarkerSize',8,'Linewidth',2,'MarkerEdgeColor',[0 0 0]);
    set(gca,'Xlim',[-150 150],'Ylim',[-150 150],'Tickdir','out'); grid on; axis square;
    plot_counter = plot_counter + 1;
end

%% Section 1.3 Effect of laser on Hits in laser and control trials (files collected from Apollo on 02/15/19)
if ~exist('plot_counter')
    plot_counter = 1;
end
filename = {'A021519003.nex';'A021519004.nex';'A021519005.nex';'A021519006.nex';'A021519007.nex';'A021519008.nex';'A021519009.nex';'A021519010.nex'};
RWmultiplier = 0.8*ones(numel(filename),1);
laserdial = [2.0;1.0;0.5;1.5;0.75;2.0;0.25;2.5];
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
staircaseterminationptNL = [];
cgradations = linspace(0.1,1,numel(filename));
cgradations = cgradations([6 4 2 5 3 7 1 8]);
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
        putativeoogcontrast = min(LMScontrastfileL(logical(stro.trial(colorstimpresent & lasertrialidxs,oog))));
        staircaseterminationptNL = [staircaseterminationptNL; LMScontrastfileNL(end)];
        
        % A simple criteria for detecting "bad files"
        tmp = stro.trial(ind & stimpresent & ~lasertrials,oog);
        outofgamutptsNL = [outofgamutptsNL; sum(tmp(end-4:end))];
        
        answers = colordirchoiceidxs(colorstimpresent & ~lasertrialidxs);
        if ~isempty(putativeoogcontrast)
            LMScontrastfileNL(LMScontrastfileNL>putativeoogcontrast) = putativeoogcontrast;
            LMScontrastfileL(LMScontrastfileL>putativeoogcontrast) = putativeoogcontrast;
        end
        contrast = unique(LMScontrastfileNL);
        correctanswers = zeros(size(contrast));
        wronganswers = zeros(size(contrast));
        trialspercontrast = zeros(size(contrast));
        
        for mm = 1:numel(contrast)
            trialspercontrast(mm) = numel(answers(LMScontrastfileNL==contrast(mm)));
            correctanswers(mm) = sum(answers(LMScontrastfileNL==contrast(mm)));
            wronganswers(mm) = trialspercontrast(mm) - correctanswers(mm);
        end
        [a,~,~] = weibullFit(contrast,[correctanswers wronganswers],'mle');
        thresholdsNL = [thresholdsNL; a];
        
    end
    CRsaccadeind = (~stimpresent & correcttrials);
    Misssaccadeind = (stimpresent & ~correcttrials);
    Hitsaccadeind = (stimpresent & correcttrials);
    FAsaccadeind = (~stimpresent & ~correcttrials);
    HitsaccadeRT = [HitsaccadeRT; stro.trial(Hitsaccadeind,saccstartidx)-stro.trial(Hitsaccadeind,fpoffidx)];
    MisssaccadeRT = [MisssaccadeRT; stro.trial(Misssaccadeind,saccstartidx)-stro.trial(Misssaccadeind,fpoffidx)];
    CRsaccadeRT = [CRsaccadeRT; stro.trial(CRsaccadeind,saccstartidx)-stro.trial(CRsaccadeind,fpoffidx)];
    FAsaccadeRT = [FAsaccadeRT; stro.trial(FAsaccadeind,saccstartidx)-stro.trial(FAsaccadeind,fpoffidx)];
    
    figure(plot_counter); subplot(221); plot(LMScontrastfileNL,'-','color',[1-cgradations(jj) 1-cgradations(jj) 1-cgradations(jj)],'Linewidth',2); hold on;
    subplot(222); plot(LMScontrastfileL,'-','color',[1-cgradations(jj) 1-cgradations(jj) 1-cgradations(jj)],'Linewidth',2); hold on;
    subplot(223); plot(laserdial(jj),HitNL(jj)./stimpresentNL(jj) - HitL(jj)./stimpresentL(jj),'o','MarkerSize',8,'MarkerFaceColor',[1-cgradations(jj) 1-cgradations(jj) 1-cgradations(jj)],'MarkerEdgeColor',[1 1 1]); hold on;
    subplot(224); plot(jj,laserdial(jj),'o','MarkerSize',8,'MarkerFaceColor',[1-cgradations(jj) 1-cgradations(jj) 1-cgradations(jj)],'MarkerEdgeColor',[1 1 1]); hold on;
    
end
subplot(221); xlabel('trials'); ylabel('contrast'); title('control staircase'); set(gca,'Ylim',[0 1.5],'TickDir','out'); axis square; hold off;
subplot(222); xlabel('trials'); ylabel('contrast'); title('laser staircase'); set(gca,'Ylim',[0 1.5],'TickDir','out'); axis square; hold off;

% plotting the differences in hits between laser and control trials
[model] = sigmoidalfit_AD(laserdial,HitNL./stimpresentNL - HitL./stimpresentL);
x = linspace(0,max(laserdial),51);
fit = model(1)*(x.^model(3))./(x.^model(3)+model(2).^model(3)); 
subplot(223);plot(linspace(0,max(laserdial),51),fit,'color',[0 0 0],'Linewidth',2); xlabel('laser dial'); ylabel('Hits control - Hits laser'); set(gca,'Ylim',[0 1],'Xlim',[0 max(laserdial)],'TickDir','out'); axis square; hold off;
[~,p] = corr(laserdial,linspace(0,7,8)');
subplot(224); xlabel('Chronology'); ylabel('laserdial'); title(strcat('p=',num2str(p,2))); set(gca,'TickDir','out'); axis square; hold off;
plot_counter = plot_counter + 1;

%% Section 1.4.1 Spatial extent of the effect: Example with randomly interleaved stim locations
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
        
        % Need to calculate the saccade RT for the left and right targets
        Leftsaccadeind = (~stimpresent & correcttrials & ind) |  (stimpresent & ~correcttrials & ind);
        Rightsaccadeind = (stimpresent & correcttrials & ind) |  (~stimpresent & ~correcttrials & ind);
        L1 = hist(stro.trial(Leftsaccadeind,saccstartidx)-stro.trial(Leftsaccadeind,fpoffidx),timebins);
        R1 = hist(stro.trial(Rightsaccadeind,saccstartidx)-stro.trial(Rightsaccadeind,fpoffidx),timebins);
        Leftsaccadehists = [Leftsaccadehists; L1];
        Rightsaccadehists = [Rightsaccadehists; R1];
        
        trials = lasertrials;
        for kk = 1:2
            if kk == 2
                trials = ~trials;
            end
            LMScontrast = sqrt(sum(LMStriplet(stimpresent & trials & ind,:).^2,2));
            putativeoogcontrast = min(LMScontrast(logical(stro.trial(stimpresent & trials & ind,oog))));
            answers = correcttrials(stimpresent & trials & ind);
            if ~isempty(putativeoogcontrast)
                LMScontrast(LMScontrast>putativeoogcontrast) = putativeoogcontrast;
            end
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
subplot(221); plot(stro.sum.exptParams.fp_x,stro.sum.exptParams.fp_y,'+','Markersize',15,'Linewidth',2);
for ii = 1:numel(unique(ib))
    subplot(221);hold on; plot(locationstested(ii,1), locationstested(ii,2),'o','color',[1-cgradations(ii) 0 cgradations(ii)],'Markersize',15,'Linewidth',2);
    ind = find(ib==ii);
    [~,p1] = equalproptest([sum(HitsL(ind)) sum(HitsNL(ind))],[sum(stimpresentL(ind)) sum(stimpresentNL(ind))],alpha);
    [~,p2] = equalproptest([sum(FAL(ind)) sum(FANL(ind))],[sum(stimabsentL(ind)) sum(stimabsentNL(ind))],alpha);
    subplot(2,2,ii+1); bar([sum(HitsL(ind)) sum(HitsNL(ind)); sum(MissL(ind)) sum(MissNL(ind)); sum(CRL(ind)) sum(CRNL(ind)); sum(FAL(ind)) sum(FANL(ind))]); hold on;
    text(1,max([sum(HitsL(ind));sum(HitsNL(ind))])+20,strcat('p1=',num2str(p1,3))); text(3.5,max([sum(FAL(ind));sum(FANL(ind))])+20,strcat('p2=',num2str(p2,3)));
    ylabel('Count'); title(num2str(locationstested(ii,:))); set(gca,'XTick',[1 2 3 4],'XTickLabel',{'H','M','CR','FA'},'TickDir','Out','XColor',[1-cgradations(ii) 0 cgradations(ii)],'YColor',[1-cgradations(ii) 0 cgradations(ii)]); drawnow;
    
end
subplot(221), hold on; axis square; set(gca,'Xlim',[-100 100],'Ylim',[-100 100],'TickDir','Out'); grid on; xlabel('X'), ylabel('Y'); hold off;
subplot(222); hold on; axis square; hold off;
subplot(223); hold on; axis square; hold off;
subplot(224); hold on; axis square; hold off;
plot_counter = plot_counter + 1;

%% 1.4.2 Spatial extent of the effect:	Example with files tested for different locations where stim loc was fixed for each file
if ~exist('plot_counter')
    plot_counter = 1;
end
filenameatRF = {'M051818002.nex';'M051818003.nex';'M051818015.nex'};
filenameoutsideRF = {'M051818004.nex';'M051818005.nex';'M051818006.nex';'M051818007.nex';'M051818008.nex';'M051818009.nex';'M051818010.nex';'M051818011.nex';'M051818012.nex';'M051818013.nex';'M051818014.nex'};
filename = [filenameatRF; filenameoutsideRF];
color = ['r';'g';'b';'k';'m';'c';'y'];
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
    uniquediridxs = unique(diridxs);
    lasertrials = logical(stro.trial(:,optstim));
    N = numel(uniquediridxs);
    RF = [RF; stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y];
    for ii = 1:N
        lasertrialidxs = logical(stro.trial(diridxs == ii-1,optstim));
        colordirchoiceidxs = correcttrials(diridxs == ii-1);
        colorstimpresent = logical(stimpresent(diridxs==ii-1));
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
    end
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
    
end

[RFlocations,ia,ib] = unique(RF,'rows');
nrows = ceil(sqrt(numel(unique(ib))));
figure(plot_counter); subplot(nrows,nrows,1);plot(0,0,'+','Markersize',15,'Linewidth',2); hold on;
subplot(nrows,nrows,9);plot(0,0,'+','Markersize',15,'Linewidth',2); hold on;
for ii = 1:numel(unique(ib))
    subplot(nrows,nrows,1);plot(RFlocations(ii,1),RFlocations(ii,2),'o','MarkerSize',7,'LineWidth',0.5,'MarkerFaceColor',color(ii),'MarkerEdgeColor',color(ii));hold on;
    ind = find(ib==ii);
    [~,p1] = equalproptest([sum(HitL(ind)) sum(HitNL(ind))],[sum(stimpresentL(ind)) sum(stimpresentNL(ind))],alpha);
    [~,p2] = equalproptest([sum(FAL(ind)) sum(FANL(ind))],[sum(stimabsentL(ind)) sum(stimabsentNL(ind))],alpha);
    subplot(nrows,nrows,ii+1), bar([sum(HitL(ind)) sum(HitNL(ind)); sum(MissL(ind)) sum(MissNL(ind)); sum(CRL(ind)) sum(CRNL(ind)); sum(FAL(ind)) sum(FANL(ind))]); hold on;
    text(1,max([sum(HitNL(ind)) sum(MissL(ind))])+15,strcat('p1=',num2str(p1,2))); text(3.5,max([sum(FAL(ind)) sum(FANL(ind))])+20,strcat('p2=',num2str(p2,2)));
    ylabel('Count'); set(gca,'XTick',[1 2 3 4],'XTickLabel',{'H','M','CR','FA'},'Ylim',[0 90],'TickDir','out'); title(num2str(RFlocations(ii,:)));axis square; set(gca,'XColor',color(ii),'YColor',color(ii)); hold off;
    
    % Another way
    tmp = sum(HitNL(ind))./sum(stimpresentNL(ind)) - sum(HitL(ind))./sum(stimpresentL(ind));
    subplot(nrows,nrows,9);plot(RFlocations(ii,1),RFlocations(ii,2),'o','MarkerSize',7,'LineWidth',0.5,'MarkerFaceColor',[0.8-tmp 0.8-tmp 0.8-tmp],'MarkerEdgeColor',[1 1 1]);hold on;
    
end
subplot(nrows,nrows,1);grid on; set(gca,'Xlim',[-80 80],'Ylim',[-80 80],'TickDir','out'); axis square; xlabel('X'); ylabel('Y'); title('RF locations');
subplot(nrows,nrows,9);grid on; set(gca,'Xlim',[-80 80],'Ylim',[-80 80],'TickDir','out'); axis square; xlabel('X'); ylabel('Y'); title('Magnitude of the effect');
plot_counter = plot_counter + 1;

%% Section: 1.4.3 - o	Examples of laser effects on performance as function of time 
% (plot the psychometric function for each of those files and see if there is any transition over time).
if ~exist('plot_counter')
    plot_counter = 1;
end

filenameopto1 = {'M042518004.nex';'M042518005.nex';'M042518008.nex';'M042518006.nex';'M042518007.nex'}; % 4/25
RWmultiplieropto1 =[1.0;1.0;1.0;1.0;1.0];
laserdialopto1 = [2.5;2.5;2.5;2.5;2.5];

filenameopto2 = {'M042618003.nex';'M042618004.nex';'M042618005.nex';'M042618006.nex';'M042618007.nex'}; % 4/26, 200 stim, 300 laser
RWmultiplieropto2 =[1.0;1.0;1.0;1.0;1.0];
laserdialopto2 = [2.5;2.5;2.5;2.5;2.5];

filenameopto3 = {'M042718006.nex';'M042718007.nex';'M042718008.nex';'M042718009.nex'}; % 4/27, 200 stim, 300 laser
RWmultiplieropto3 =[1.0;1.0;1.0;1.0];
laserdialopto3 = [2.7;2.7;2.7;2.7];

filenameopto4 = {'M042718012.nex';'M042718013.nex';'M042718014.nex'}; % 4/27, 100 stim, 150 laser
RWmultiplieropto4 =[1.0;1.0;1.0];
laserdialopto4 = [3.2;3.2;3.2];

filenameopto5 = {'A011719009.nex';'A011719010.nex';'A011719011.nex';'A011719012.nex'};
RWmultiplieropto5 = [1.0;1.0;1.0;1.0];
laserdialopto5 = [2.4;2.4;2.4;2.4];
oogvals = [];
nrows = 5;
figure(plot_counter);
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
        
        trials = lasertrials;
        % only evaluating the laser trials
        for kk = 1:2
            if kk == 2
                trials = ~trials;
            end
            LMScontrast = [];
            LMScontrast = sqrt(sum(LMStriplet(stimpresent & trials,:).^2,2));
            putativeoogcontrast = min(LMScontrast(logical(stro.trial(stimpresent & trials,oog))));
            answers = correcttrials(stimpresent & trials);
            if ~isempty(putativeoogcontrast)
                LMScontrast(LMScontrast>putativeoogcontrast) = putativeoogcontrast;
                oogvals = [oogvals; putativeoogcontrast];
            end
            
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
        for jj = 1:numel(contrastL)
            subplot(nrows,5,mm*5-4);  plot(contrastL(jj),correctanswersL(jj)./(correctanswersL(jj)+wronganswersL(jj)),'o','MarkerSize',correctanswersL(jj) + wronganswersL(jj),'MarkerFacecolor',[1-cgradations(aa) 1-cgradations(aa) 1-cgradations(aa)],'MarkerEdgeColor',[1 1 1]); hold on;
        end
        plot(contrastlattice,fitL,'-','Linewidth',2,'color',[1-cgradations(aa) 1-cgradations(aa) 1-cgradations(aa)]);
        subplot(nrows,5,mm*5-3); plot(LMScontrastL,'-','Linewidth',2,'color',[1-cgradations(aa) 1-cgradations(aa) 1-cgradations(aa)]); hold on;
        subplot(nrows,5,mm*5-2); plot(aa,HitL(end)./stimpresentL(end),'o','MarkerFacecolor',[1-cgradations(aa) 1-cgradations(aa) 1-cgradations(aa)],'MarkerEdgeColor',[1 1 1]); hold on;
        subplot(nrows,5,mm*5-1); plot(LMScontrastNL,'-','Linewidth',2,'color',[1-cgradations(aa) 1-cgradations(aa) 1-cgradations(aa)]); hold on;
        fitNL = gNL*(1-exp(-((contrastlattice./aNL).^bNL)));
        for jj = 1:numel(contrastNL)
            subplot(nrows,5,mm*5);  plot(contrastNL(jj),correctanswersNL(jj)./(correctanswersNL(jj)+wronganswersNL(jj)),'o','MarkerSize',correctanswersNL(jj) + wronganswersNL(jj),'MarkerFacecolor',[1-cgradations(aa) 1-cgradations(aa) 1-cgradations(aa)],'MarkerEdgeColor',[1 1 1]); hold on;
        end
        plot(contrastlattice,fitNL,'-','Linewidth',2,'color',[1-cgradations(aa) 1-cgradations(aa) 1-cgradations(aa)]);
    end
    subplot(nrows,5,mm*5-4); hold on; xlabel('contrast'); ylabel('prop of correct'); set(gca,'Tickdir','out','Ylim',[0 1],'Xlim',[0 1.3]); title(strcat('SessionL:',num2str(mm))); hold off;
    subplot(nrows,5,mm*5-3); hold on; xlabel('trial no.'); ylabel('contrast'); set(gca,'Tickdir','out','Ylim',[0 1.5]); 
    if mm == 1
       title('laser staircase');
    end
    hold off;
    subplot(nrows,5,mm*5-2); hold on; xlabel('file no.'); ylabel('phits control-laser'); set(gca,'Tickdir','out','Ylim',[0 1]); hold off;
    subplot(nrows,5,mm*5-1); hold on; xlabel('trial no.'); ylabel('contrast'); set(gca,'Tickdir','out','Ylim',[0 1.5]); 
    if mm == 1
       title('control staircase'); 
    end
    hold off;
    subplot(nrows,5,mm*5); hold on; xlabel('contrast'); ylabel('prop of correct'); set(gca,'Tickdir','out','Ylim',[0 1],'Xlim',[0 1.3]); title(strcat('SessionNL:',num2str(mm))); hold off;
    plot_counter = plot_counter + 1;
end

%% Section 2.1: Example file of SCstimcue from Maui and Apollo
% 2.2: Example file of SCstimcue from Maui and Apollo

if ~exist('plot_counter')
    plot_counter = 1;
end
filename = ['M042718005.nex';'A012319003.nex'];
figure(plot_counter); set(gcf,'Name','SC stimcue');
nrows = size(filename,1);
for jj = 1:size(filename,1)
    stro = nex2stro(findfile(filename(jj,:)));
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
    subplot(nrows,2,2*jj-1); plot(stro.sum.exptParams.rf_x/10,stro.sum.exptParams.rf_y/10,'o','MarkerSize',12,'MarkerEdgeColor',[0 0 0],'Linewidth',2); hold on;
    subplot(nrows,2,2*jj); plot(stro.sum.exptParams.rf_x/10,stro.sum.exptParams.rf_y/10,'o','MarkerSize',12,'MarkerEdgeColor',[0 0 0],'Linewidth',2); hold on;
    RFloc = [stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y];
    if norm(RFloc)>100
        plotlimit = 18;
    else 
        plotlimit = 12;
    end
    
    for i = 1:size(stro.trial,1)
        x = stro.ras{i,2}*4096/400;
        y = stro.ras{i,3}*4096/400;
        t = stro.ras{i,5}+[0:1:length(x)-1]/samplerate;
        Lt = t>fpoff_t(i) & t < fpoff_t(i) +.3;
        if ~Lcatchtrials(i)
            if ~Llaser(i)
                subplot(nrows,2,2*jj-1); plot(x(Lt),y(Lt),'color',[0.5 0.5 0.5],'Linewidth',1); hold on;
            else
                subplot(nrows,2,2*jj-1); plot(x(Lt),y(Lt),'color',[0 0.5 1],'Linewidth',1); hold on;
            end
        else
            if ~Llaser(i)
                subplot(nrows,2,2*jj); plot(x(Lt),y(Lt),'--','color',[0.5 0.5 0.5],'Linewidth',1); hold on;
            else
                subplot(nrows,2,2*jj); plot(x(Lt),y(Lt),'--','color',[0 0.5 1],'Linewidth',1); hold on;
            end
        end
    end
    subplot(nrows,2,2*jj-1); axis square; set(gca,'Xlim',[-1*plotlimit plotlimit],'Ylim',[-1*plotlimit plotlimit],'Tickdir','out'); title(filename(jj,1));
    subplot(nrows,2,2*jj); axis square; set(gca,'Xlim',[-1*plotlimit plotlimit],'Ylim',[-1*plotlimit plotlimit],'Tickdir','out'); title(strcat(filename(jj,1),':catch trials'));
end
plot_counter = plot_counter + 1;

%% Section 2.3: A population data of probability of making a saccade into the target location in a laser and control trials for Maui and Apollo
if ~exist('plot_counter')
    plot_counter = 1;
end
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filename = fetch(conn,'SELECT filename FROM SCstimcue');
training = fetch(conn,'SELECT training FROM SCstimcue');
close(conn);
idx = ismember(filename,'M051718020.nex');
filename(idx) = [];
training(idx) = [];
ind = find(strcmp(training,'no'));
L = ceil(sqrt(numel(ind)));
probsacc_withinRF = cell(numel(ind),1);
probsacc_withinRFL = cell(numel(ind),1);
probsacc_withinRFNL = cell(numel(ind),1);
probsacc_outsideRF = cell(numel(ind),1);
probsacc_outsideRFL = cell(numel(ind),1);
probsacc_outsideRFNL = cell(numel(ind),1);
monkey_ID = [];
ploteachfile = 0;
if ploteachfile
    figure(plot_counter); set(gcf,'Name','Laser trials');
end
for ii = 1:numel(ind)
    targetlocations = [];
    targethitsallNL = []; % non-laser trials
    targethitsallL = []; % laser trials
    stimpresentationNL = [];% non-laser trials
    stimpresentationL = []; % laser trials
    RFloc = [];
    fileofinterest = char(filename(ind(ii),:));
    monkey_ID = [monkey_ID; fileofinterest(1)];
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
    RFloc = [RFloc; stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y];
    if norm(RFloc)>100
        halfwidth = 5.0; % in degrees of visual angle
    else
        halfwidth = 3.0; % in degrees of visual angle
    end
    for i = 1:size(stro.trial,1)
        if ~Llaser(i)
            x = stro.ras{i,2}*4096/400;
            y = stro.ras{i,3}*4096/400;
            t = stro.ras{i,5}+[0:1:length(x)-1]/samplerate;
            Lt = t>fpoff_t(i) & t < fpoff_t(i)+0.30;
            target = stro.trial(i,[12 13]); % x and y target locations of that trial
            idx = find(sum(target==C,2)==2);
            if sum(x(Lt)>(target(1)/10)-halfwidth & x(Lt)<(target(1)/10)+halfwidth & y(Lt)>(target(2)/10)-halfwidth & y(Lt)<(target(2)/10)+halfwidth)>10
                targethitsNL(idx) = targethitsNL(idx) + 1;
            end
            SPNL(idx) = SPNL(idx) + 1;
        end
        if Llaser(i) % Analyzing laser trials
            x = stro.ras{i,2}*4096/400;
            y = stro.ras{i,3}*4096/400;
            t = stro.ras{i,5}+[0:1:length(x)-1]/samplerate;
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
        end
    end
    targethitsallNL = [targethitsallNL; targethitsNL];
    targethitsallL = [targethitsallL; targethitsL];
    stimpresentationNL = [stimpresentationNL; SPNL];
    stimpresentationL = [stimpresentationL; SPL];
    withinRFidx = all(C==RFloc,2);
    probsacc_withinRFL{ii} = targethitsallL(withinRFidx)./stimpresentationL(withinRFidx);
    probsacc_withinRFNL{ii} = targethitsallNL(withinRFidx)./stimpresentationNL(withinRFidx);
    probsacc_outsideRFL{ii} = targethitsallL(~withinRFidx)./stimpresentationL(~withinRFidx);
    probsacc_outsideRFNL{ii} = targethitsallNL(~withinRFidx)./stimpresentationNL(~withinRFidx);
    probsacc_withinRF{ii} = (targethitsallNL(withinRFidx)./stimpresentationNL(withinRFidx)) -(targethitsallL(withinRFidx)./stimpresentationL(withinRFidx));
    probsacc_outsideRF{ii} = (targethitsallNL(~withinRFidx)./stimpresentationNL(~withinRFidx)) -(targethitsallL(~withinRFidx)./stimpresentationL(~withinRFidx));
    timebins = 0:0.01:0.3;
    if ploteachfile
        figure(plot_counter);subplot(L,L,ii);
        plot(RFloc(:,1)/10,RFloc(:,2)/10,'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
        for jj = 1:size(targethitsallL)
            if targethitsallL(jj)
                if sum(targetlocations(jj,:)~=0,2)
                    plot(targetlocations(jj,1)/10,targetlocations(jj,2)/10,'o','MarkerSize',5*targethitsallL(jj)/(stimpresentationNL(jj)+stimpresentationL(jj)),'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); hold on;
                end
            end
        end
        grid on; set(gca,'Xlim',[-18 18],'Ylim',[-18 18],'Tickdir','out'); title(char(filename(ind(ii)))); xlabel('X'); ylabel('Y'); hold off;
    end
end
plot_counter = plot_counter + 1;
idx  = monkey_ID == 'M';
bins = -1:0.1:1;
[~,p1] = kstest2(cell2mat(probsacc_withinRF(idx)),cell2mat(probsacc_outsideRF(idx))); % For Maui
[~,p2] = kstest2(cell2mat(probsacc_withinRF(~idx)),cell2mat(probsacc_outsideRF(~idx))); % For Apollo
figure(plot_counter); set(gcf,'Name','Population: SCstimcue');
subplot(221); histogram(cell2mat(probsacc_withinRF(idx)),bins,'Normalization','probability'); hold on; histogram(cell2mat(probsacc_outsideRF(idx)),bins,'Normalization','probability'); 
xlabel('P(sacc control - laser)'); ylabel('pdf'); axis square; set(gca,'Tickdir','out'); title('M');
subplot(222); histogram(cell2mat(probsacc_withinRF(~idx)),bins,'Normalization','probability'); hold on; histogram(cell2mat(probsacc_outsideRF(~idx)),bins,'Normalization','probability'); 
xlabel('P(sacc control - laser)'); ylabel('pdf'); axis square; set(gca,'Tickdir','out'); title('A');
subplot(223);cdfplot(cell2mat(probsacc_withinRF(idx))); hold on; cdfplot(cell2mat(probsacc_outsideRF(idx))); xlabel('P(sacc control - laser)'); 
ylabel('cdf'); axis square; set(gca,'Tickdir','out'); title(strcat('M: p=',num2str(p1,2))); hold off;
subplot(224);cdfplot(cell2mat(probsacc_withinRF(~idx))); hold on; cdfplot(cell2mat(probsacc_outsideRF(~idx))); xlabel('P(sacc control - laser)'); 
ylabel('cdf'); axis square; set(gca,'Tickdir','out'); title(strcat('A: p=',num2str(p2,2))); hold off;
plot_counter = plot_counter + 1;

[~,p3] = ttest2(cell2mat(probsacc_withinRFL(idx)),cell2mat(probsacc_withinRFNL(idx))); % For Maui: within RF
[~,p4] = ttest2(cell2mat(probsacc_outsideRFL(idx)),cell2mat(probsacc_outsideRFNL(idx))); % For Maui: outside RF
[~,p5] = ttest2(cell2mat(probsacc_withinRFL(~idx)),cell2mat(probsacc_withinRFNL(~idx))); % For Apollo: within RF
[~,p6] = ttest2(cell2mat(probsacc_outsideRFL(~idx)),cell2mat(probsacc_outsideRFNL(~idx))); % For Apollo: outside RF
figure(plot_counter);set(gcf,'Name','Population: Within and outside RF');
bins = 0:0.1:1;
subplot(221); histogram(cell2mat(probsacc_withinRFL(idx)),bins,'Normalization','probability','FaceColor',[0 0.5 1]); hold on; histogram(cell2mat(probsacc_withinRFNL(idx)),bins,'Normalization','probability','FaceColor',[0.5 0.5 0.5]); 
xlabel('P(sacc)'); ylabel('pdf'); axis square; set(gca,'Tickdir','out'); title(strcat('M:within,p=',num2str(p3,2)));
subplot(222);histogram(cell2mat(probsacc_outsideRFL(idx)),bins,'Normalization','probability','FaceColor',[0 0.5 1]); hold on; histogram(cell2mat(probsacc_outsideRFNL(idx)),bins,'Normalization','probability','FaceColor',[0.5 0.5 0.5]); 
xlabel('P(sacc)'); ylabel('pdf'); axis square; set(gca,'Tickdir','out'); title(strcat('M:outside,p=',num2str(p4,2)));
subplot(223); histogram(cell2mat(probsacc_withinRFL(~idx)),bins,'Normalization','probability','FaceColor',[0 0.5 1]); hold on; histogram(cell2mat(probsacc_withinRFNL(~idx)),bins,'Normalization','probability','FaceColor',[0.5 0.5 0.5]); 
xlabel('P(sacc)'); ylabel('pdf'); axis square; set(gca,'Tickdir','out'); title(strcat('A:within,p=',num2str(p5,2)));
subplot(224);histogram(cell2mat(probsacc_outsideRFL(~idx)),bins,'Normalization','probability','FaceColor',[0 0.5 1]); hold on; histogram(cell2mat(probsacc_outsideRFNL(~idx)),bins,'Normalization','probability','FaceColor',[0.5 0.5 0.5]); 
xlabel('P(sacc)'); ylabel('pdf'); axis square; set(gca,'Tickdir','out'); title(strcat('A:outside,p=',num2str(p6,2)));
plot_counter = plot_counter + 1;

[~,tp1] = ttest2(cell2mat(probsacc_withinRFL(idx)),cell2mat(probsacc_withinRFNL(idx)));
[~,tp2] = ttest2(cell2mat(probsacc_outsideRFL(idx)),cell2mat(probsacc_outsideRFNL(idx)));
[~,tp3] = ttest2(cell2mat(probsacc_withinRFL(~idx)),cell2mat(probsacc_withinRFNL(~idx)));
[~,tp4] = ttest2(cell2mat(probsacc_outsideRFL(~idx)),cell2mat(probsacc_outsideRFNL(~idx)));

figure(plot_counter); set(gcf,'Name','A different representation');
subplot(221); plot([cell2mat(probsacc_withinRFL(idx)) cell2mat(probsacc_withinRFNL(idx))]','-o','color',[0 0 0]); set(gca,'Xlim',[0 3]); title('M:within'); ylabel('P(sacc)');
subplot(222); plot([cell2mat(probsacc_outsideRFL(idx)) cell2mat(probsacc_outsideRFNL(idx))]','-o','color',[0 0 0]); set(gca,'Xlim',[0 3]); title('M:outside'); ylabel('P(sacc)');
subplot(223); plot([cell2mat(probsacc_withinRFL(~idx)) cell2mat(probsacc_withinRFNL(~idx))]','-o','color',[0 0 0]); set(gca,'Xlim',[0 3]); title('A:within'); ylabel('P(sacc)');
subplot(224); plot([cell2mat(probsacc_outsideRFL(~idx)) cell2mat(probsacc_outsideRFNL(~idx))]','-o','color',[0 0 0]); set(gca,'Xlim',[0 3]);title('A:outside'); ylabel('P(sacc)');
plot_counter = plot_counter + 1;

%% Section 2.4: Latency differences between laser and control trials (data divided by within and outside RF)
if ~exist('plot_counter')
    plot_counter = 1;
end
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filename = fetch(conn,'SELECT filename FROM SCstimcue');
training = fetch(conn,'SELECT training FROM SCstimcue');
close(conn);
ind = find(strcmp(training,'no'));
L = ceil(sqrt(numel(ind)));
binwidth = .005;
bins = -0.4:binwidth:0.6;
withinRFlatencyL = cell(numel(ind),1);
withinRFlatencyNL = cell(numel(ind),1);
outsideRFlatencyL = cell(numel(ind),1);
outsideRFlatencyNL = cell(numel(ind),1);
for ii = 1:numel(ind)
    targetlocations = [];
    targethitsallNL = []; % non-laser trials
    targethitsallL = []; % laser trials
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
    [C,ia,ib] = unique(stro.trial(:,12:13),'rows');
    RFloc = [RFloc; stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y];
    latency = sacint_t - fpoff_t;
%     latency(isnan(latency)) = 0.5; % Artificially allocating all the NaNs as 0.5 ms 
    withinRFlatencyL{ii} = latency(all(stro.trial(:,[12 13]) == RFloc,2) & Llaser & ~Lcatchtrials);
    withinRFlatencyNL{ii} = latency(all(stro.trial(:,[12 13]) == RFloc,2) & ~Llaser & ~Lcatchtrials);
    outsideRFlatencyL{ii} = latency(all(stro.trial(:,[12 13]) ~= RFloc,2) & Llaser & ~Lcatchtrials);
    outsideRFlatencyNL{ii} = latency(all(stro.trial(:,[12 13]) ~= RFloc,2) & ~Llaser & ~Lcatchtrials);
end

monkeyIdxs = cell2mat(filename(ind));
monkeyIdxs = monkeyIdxs(:,1);
idxs = monkeyIdxs == 'M'; % For Maui
withinRFlaser1 = cell2mat(withinRFlatencyL(idxs)); withinRFlaser1 = withinRFlaser1(~isnan(withinRFlaser1)); 
withinRFcontrol1 = cell2mat(withinRFlatencyNL(idxs)); withinRFcontrol1 = withinRFcontrol1(~isnan(withinRFcontrol1)); 
outsideRFlaser1 = cell2mat(outsideRFlatencyL(idxs)); outsideRFlaser1 = outsideRFlaser1(~isnan(outsideRFlaser1)); 
outsideRFcontrol1 = cell2mat(outsideRFlatencyNL(idxs)); outsideRFcontrol1 = outsideRFcontrol1(~isnan(outsideRFcontrol1)); 
withinRFlaser2 = cell2mat(withinRFlatencyL(~idxs)); withinRFlaser2 = withinRFlaser2(~isnan(withinRFlaser2)); 
withinRFcontrol2 = cell2mat(withinRFlatencyNL(~idxs)); withinRFcontrol2 = withinRFcontrol2(~isnan(withinRFcontrol2)); 
outsideRFlaser2 = cell2mat(outsideRFlatencyL(~idxs)); outsideRFlaser2 = outsideRFlaser2(~isnan(outsideRFlaser2)); 
outsideRFcontrol2 = cell2mat(outsideRFlatencyNL(~idxs)); outsideRFcontrol2 = outsideRFcontrol2(~isnan(outsideRFcontrol2)); 

p1 = ranksum(withinRFlaser1,withinRFcontrol1);
p2 = ranksum(outsideRFlaser1,outsideRFcontrol1);
p3 = ranksum(withinRFlaser2,withinRFcontrol2);
p4 = ranksum(outsideRFlaser2,outsideRFcontrol2);
bins = 0:0.02:0.5;
figure(plot_counter); set(gcf,'Name','Latency analysis');
subplot(221); histogram(withinRFlaser1,bins,'Normalization','probability','FaceColor',[0 0.5 1]); hold on; histogram(withinRFcontrol1,bins,'Normalization','probability','FaceColor',[0.5 0.5 0.5]);
plot(median(withinRFlaser1),0,'v','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]); plot(median(withinRFcontrol1),0,'v','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
xlabel('latency'); ylabel('pdf'); title('M: Within RF'); axis square; set(gca,'Tickdir','out'); hold off;
subplot(222); histogram(outsideRFlaser1,bins,'Normalization','probability','FaceColor',[0 0.5 1]); hold on; histogram(outsideRFcontrol1,bins,'Normalization','probability','FaceColor',[0.5 0.5 0.5]);
plot(median(outsideRFlaser1),0,'v','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]); plot(median(outsideRFcontrol1),0,'v','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
xlabel('latency'); ylabel('pdf'); title('M: Outside RF'); axis square; set(gca,'Tickdir','out'); hold off;
subplot(223); histogram(withinRFlaser2,bins,'Normalization','probability','FaceColor',[0 0.5 1]); hold on; histogram(withinRFcontrol2,bins,'Normalization','probability','FaceColor',[0.5 0.5 0.5]);
plot(median(withinRFlaser2),0,'v','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]); plot(median(withinRFcontrol2),0,'v','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
xlabel('latency'); ylabel('pdf'); title('A: Within RF'); axis square; set(gca,'Tickdir','out'); hold off;
subplot(224); histogram(outsideRFlaser2,bins,'Normalization','probability','FaceColor',[0 0.5 1]); hold on; histogram(outsideRFcontrol2,bins,'Normalization','probability','FaceColor',[0.5 0.5 0.5]);
plot(median(outsideRFlaser2),0,'v','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]); plot(median(outsideRFcontrol2),0,'v','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
xlabel('latency'); ylabel('pdf'); title('A: Outside RF'); axis square; set(gca,'Tickdir','out'); hold off;
plot_counter = plot_counter + 1;

%% Section 2.5: Control file where target locations were same as DToneloc saccade target locations
if ~exist('plot_counter')
    plot_counter = 1;
end
filename = ['M051718020.nex'];
nrows = size(filename,1);
figure(plot_counter); set(gcf,'Name','SC stimcue: DToneloc saccade locations'); 
for jj = 1:size(filename,1)
    stro = nex2stro(findfile(filename(jj,:)));
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
    analogstrtimeind = strcmp(stro.sum.rasterCells(1,:),'anlgStartTime');
    spikeind = strcmp(stro.sum.rasterCells(1,:),'sig001a');
    heyetraceind = strcmp(stro.sum.rasterCells(1,:),'AD11');
    veyetraceind = strcmp(stro.sum.rasterCells(1,:),'AD12');
    lasertraceind = strcmp(stro.sum.rasterCells(1,:),'AD13');
    
    Llaser = ~isnan(laseron_t);
    samplerate = stro.sum.analog.storeRates{1};
    subplot(nrows,2,2*jj-1); hold on;
    subplot(nrows,2,2*jj); hold on;
    for i = 1:size(stro.trial,1)
        x = stro.ras{i,heyetraceind}*4096/400;
        y = stro.ras{i,veyetraceind}*4096/400;
        t = stro.ras{i,analogstrtimeind}+[0:1:length(x)-1]/samplerate;
        Lt = t>fpoff_t(i) & t < fpoff_t(i) +.3;
        if ~Lcatchtrials(i)
            if ~Llaser(i)
                subplot(nrows,2,2*jj-1); plot(x(Lt),y(Lt),'color',[0.5 0.5 0.5],'Linewidth',1); hold on;
            else
                subplot(nrows,2,2*jj-1); plot(x(Lt),y(Lt),'color',[0 0.5 1],'Linewidth',1); hold on;
            end
        else
            if ~Llaser(i)
                subplot(nrows,2,2*jj); plot(x(Lt),y(Lt),'--','color',[0.5 0.5 0.5],'Linewidth',1); hold on;
            else
                subplot(nrows,2,2*jj); plot(x(Lt),y(Lt),'--','color',[0 0.5 1],'Linewidth',1); hold on;
            end
        end
    end
    subplot(nrows,2,2*jj-1); axis square; set(gca,'Xlim',[-5 5],'Ylim',[-5 5],'Tickdir','out'); title(filename(jj,1));
    subplot(nrows,2,2*jj); axis square; set(gca,'Xlim',[-5 5],'Ylim',[-5 5],'Tickdir','out'); title(strcat(filename(jj,1),':catch trials'));
end
plot_counter = plot_counter + 1;

%% Section 3.1 : Neuronal data: One example file of excitation and suppression
if ~exist('plot_counter')
    plot_counter = 1;
end
filename = {'A012319005.nex'; 'A020119004.nex'};
titles = {'activation';'suppresion'};
N = numel(filename); 
L = ceil(sqrt(N));
binwidth = .005;
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
            subplot(2,2,jj); plot(spiketimes(spiketimes>laserontime-0.3 & spiketimes<laserofftime+0.2)-laserontime,count1*ones(size(spiketimes(spiketimes>laserontime-0.3 & spiketimes<laserofftime+0.2))),'k.'); hold on;
            PSTHlaser = PSTHlaser + hist(spiketimes-laserontime, bins);
            count1 = count1 + 1;
            laserbins = bins>=0 & bins<=timedurlaser;
            tmp_baselineFR = [tmp_baselineFR; mean(PSTHlaser(baselinebins))/0.3];
            tmp_laserFR = [tmp_laserFR; mean(PSTHlaser(laserbins))/timedurlaser];
        else 
        end
    end
    baselineFR = [baselineFR; mean(tmp_baselineFR)];
    laserFR = [laserFR; mean(tmp_laserFR)];
    %hold on; subplot(2,2,jj); plot(bins,PSTHlaser,'color',[0 0.5 1],'Linewidth',2); 
    line([0 0],[0 count1]); line([timedurlaser timedurlaser],[0 count1]); xlabel('time'); ylabel('trials'); title(titles{jj}); 
    set(gca,'Xlim',[-0.3 0.5],'Tickdir','out'); axis square; hold off;
end

% Section 3.2 : Population summary of excitation and suppression
% I am merging this with previous section
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
figure(plot_counter); 
subplot(223); plot(baselineFR(Apollo_idxs),laserFR(Apollo_idxs),'o','MarkerSize',7,'LineWidth',1.0,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(baselineFR(Maui_idxs),laserFR(Maui_idxs),'o','MarkerSize',7,'LineWidth',1.0,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on; 
axis square; xlabel('baselineFR'); ylabel('laserFR'); set(gca,'Xlim',[0 90],'Ylim',[0 90],'TickDir','out'); line([0 90],[0 90]); legend('A','M'); hold off;
subplot(224); plot(baselineFR(~Singleunitidxs),laserFR(~Singleunitidxs),'o','MarkerSize',7,'LineWidth',1.0,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(baselineFR(Singleunitidxs),laserFR(Singleunitidxs),'o','MarkerSize',7,'LineWidth',1.0,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on; 
axis square; xlabel('baselineFR'); ylabel('laserFR'); set(gca,'Xlim',[0 90],'Ylim',[0 90],'TickDir','out'); line([0 90],[0 90]); legend('Multi','Single'); hold off;
plot_counter = plot_counter + 1;

% Additional analysis of waveforms
Suppressionvsactivation = zeros(size(baselineFR));
Suppressionvsactivation(baselineFR>laserFR) = 0; % 'Supression'
Suppressionvsactivation(baselineFR<laserFR) = 1; % 'Activation'
idx = find(Singleunitidxs);
N = ceil(sqrt(numel(idx)));
figure(plot_counter); set(gcf,'Name','Plotting waveforms of isolated units');
peak = []; trough = [];
timediff = [];
newfilename = [filename;filename2];
for ii = 1:numel(idx)
    ind = idx(ii);
    stro = nex2stro(findfile(char(newfilename(ind,:))));
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
    timediff = [timediff; time(t2)-time(t1)];
    figure(plot_counter); subplot(N,N,ii); plot(mean(cell2mat(stro.ras(:,2))',2),'-','color',color,'Linewidth',2); hold on;
    line([t1 t1],[0 val1]); line([t2 t2],[0 val2]);line([0 numel(waveform)],[0 0]); hold off;
end
plot_counter = plot_counter + 1;
% Plotting population stats

peaktotroughratio = peak./trough;
p = ranksum(peaktotroughratio(logical(Suppressionvsactivation(idx))),peaktotroughratio(~Suppressionvsactivation(idx)));
peaktotroughratio = peak./abs(trough);

figure(plot_counter); subplot(121); plot(abs(timediff(logical(Suppressionvsactivation(idx)))),peaktotroughratio(logical(Suppressionvsactivation(idx))),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(abs(timediff(~(Suppressionvsactivation(idx)))),peaktotroughratio(~(Suppressionvsactivation(idx))),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); ylabel('Peak Amplitude/Trough Amplitude'); xlabel('|Peak time - Trough time (us)|'); 
set(gca,'Tickdir','out','Xlim',[100 400],'Ylim',[0.3 1.0]); legend('Activation','Suppression'); axis square; hold off;
subplot(122); plot(peak(logical(Suppressionvsactivation(idx))),abs(trough(logical(Suppressionvsactivation(idx)))),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(peak(~(Suppressionvsactivation(idx))),abs(trough(~(Suppressionvsactivation(idx)))),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); xlabel('Peak Amplitude'); ylabel('Trough Amplitude'); 
set(gca,'Tickdir','out','Xlim',[0.02 0.14],'Ylim',[0.02 0.14]); legend('Activation','Suppression'); axis square; hold off;
plot_counter = plot_counter + 1;

%% Section 3.3: 1 example file from FixStim interfreq to show an example files about sinusoidal modulation
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
                set(gca,'XLim', [offset(1) dur+offset(2)],'Ytick',[],'YLim',[-4 sum(L)+5]);
                if modulator && j > 0
                    title(['Frequency: ',num2str(round(f)),' Hz']);
                else
                    title(['Frequency: ',num2str(j)]);
                end
                
                % PSTH
                subplot(length(N), 2, 2*(find(N==j))-0); hold on;
                plot(bins,PSTH,'k-','LineWidth',2);
                set(gca,'YLim',[0 10*ceil(max(PSTH(:)/10))+1]);
                set(gca,'Xlim',[offset(1) dur+offset(2)]);
                xlabel('Time (s)','FontSize',12);
                ylabel('Response (sp/s)','FontSize',12);
            end
        end
        plot_counter = plot_counter + 1;
    end
end

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
        
        
        subplot(121); polarplot(unique(x)',mu+sem,'-','Color',[0 0 0]); hold on;
        polarplot(unique(x)',mu,'-','Color',[0 0 0],'LineWidth',2); polarplot(unique(x)',mu-sem,'-','Color',[0 0 0]);
        polarplot(linspace(0,2*pi,100),repmat(mean(baselines),1,100)); title('Orientation tuning');
        
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
            subplot(122); hold on; plot(x,y,'k.'); plot(xx,fit,'k','Linewidth',2);
            orient = unique(orients(trlidxs)); % for axis labeling
            title(['orientation: ',num2str(orient*180/pi,3),' deg']);
            set(gca,'XScale','log');
            xlabel('spatial frequency (cyc/deg)');
            ylabel('response (sp/sec)');
            set(gca,'YLim',[0 max(spikerates)]);
            mu = []; sem = [];
            for j = unique(x)'
                mu(j == unique(x))= mean(y(x==j));
                sem(j == unique(x)) = std(y(x==j))/sqrt(sum(x==j));
            end
        end
        errorbar(unique(x),mu,sem,'k');
        plot([min(x) max(x)], repmat(mean(baselines),1,2),'k:'); axis tight; 
        set(gca,'Tickdir','out'); axis square; hold off;
    end
end