% Analyzing data files from 1 day (8 files) to observe the effect of laser
% power on differences in hits in laser and control trials
% Author - Abhishek De, 2/19
close all; clearvars;
filename = {'A021519003.nex';'A021519004.nex';'A021519005.nex';'A021519006.nex';'A021519007.nex';'A021519008.nex';'A021519009.nex';'A021519010.nex'};
RWmultiplier = 0.8*ones(numel(filename),1);
laserdial = [2.0;1.0;0.5;1.5;0.75;2.0;0.25;2.5];
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
plot_counter = 1;
RF = [];
stimabsentL = []; stimabsentNL = [];
stimpresentL = []; stimpresentNL = [];
HitsaccadeRT = []; MisssaccadeRT = [];
CRsaccadeRT = []; FAsaccadeRT = [];
smatrix = [];
newfilename = [];
newRWmultiplier = [];
newstimsize = [];
newlaserdial = [];
thresholdsNL = [];
outofgamutptsNL = [];
count = 1;
staircaseterminationptNL = [];
cgradations = linspace(0,1,numel(filename));
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
    
    figure(plot_counter); subplot(221); plot(LMScontrastfileNL,'-','color',[cgradations(jj) 0 1-cgradations(jj)],'Linewidth',2); hold on; 
    subplot(222); plot(LMScontrastfileL,'-','color',[cgradations(jj) 0 1-cgradations(jj)],'Linewidth',2); hold on;
    subplot(223); plot(laserdial(jj),HitNL(jj)./stimpresentNL(jj) - HitL(jj)./stimpresentL(jj),'o','MarkerSize',8,'MarkerFaceColor',[cgradations(jj) 0 1-cgradations(jj)],'MarkerEdgeColor',[1 1 1]); hold on;
    subplot(224); plot(jj,laserdial(jj),'o','MarkerSize',8,'MarkerFaceColor',[cgradations(jj) 0 1-cgradations(jj)],'MarkerEdgeColor',[1 1 1]); hold on;
    
end
subplot(221); xlabel('trials'); ylabel('contrast'); title('control staircase'); set(gca,'Ylim',[0.3 1.3],'TickDir','out'); axis square; hold off;
subplot(222); xlabel('trials'); ylabel('contrast'); title('laser staircase'); set(gca,'Ylim',[0.3 1.3],'TickDir','out'); axis square; hold off;

% plotting the differences in hits between laser and control trials
[model] = sigmoidalfit_AD(laserdial,HitNL./stimpresentNL - HitL./stimpresentL);
fit = 1./(1+exp(model(1)*(linspace(0,max(laserdial),51)-model(2))));
subplot(223);plot(linspace(0,max(laserdial),51),fit,'color',[0 0 0],'Linewidth',2); xlabel('laser dial'); ylabel('Hits control - Hits laser'); set(gca,'Ylim',[0 1],'Xlim',[0 max(laserdial)],'TickDir','out'); axis square; hold off;
[~,p] = corr(laserdial,linspace(0,7,8)');
subplot(224); xlabel('Chronology'); ylabel('laserdial'); title(strcat('p=',num2str(p,2))); set(gca,'TickDir','out'); axis square; hold off;
plot_counter = plot_counter + 1;
