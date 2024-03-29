% DTonelocAnalysis2 - for analysing files with interleaved staircases at
% different spatial locations
% Author - Abhishek De, 5/18
close all; clearvars;
filename = {'M052318001.nex';'M052318002.nex';'M052318003.nex';'M052318004.nex'}; % 5/23, opto session, tested 3 spatial locations
% filename = {'M052518007.nex';'M052518009.nex'}; % 5/25, opto session, tested 3 spatial locations
plot_counter = 1;
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
decisioncriterion = 0;
xvals = linspace(-4,4,101);
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
figure(plot_counter); set(gcf,'Name','Averaged data');
for ii = 1:numel(unique(ib))
    ind = find(ib==ii);
%     keyboard;
    [~,p1] = equalproptest([sum(HitsL(ind)) sum(HitsNL(ind))],[sum(stimpresentL(ind)) sum(stimpresentNL(ind))],alpha);
    [~,p2] = equalproptest([sum(FAL(ind)) sum(FANL(ind))],[sum(stimabsentL(ind)) sum(stimabsentNL(ind))],alpha);
    subplot(nrows,5,1 + 5*(ii-1)); bar([sum(HitsL(ind)) sum(HitsNL(ind)); sum(MissL(ind)) sum(MissNL(ind)); sum(CRL(ind)) sum(CRNL(ind)); sum(FAL(ind)) sum(FANL(ind))]); hold on;
    text(3.5,max([FAL;FANL])+5,strcat('p1=',num2str(p1))); text(3.5,max([FAL;FANL])+10,strcat('p2=',num2str(p2)));
    ylabel('Count'); title(num2str(locationstested(ii,:))); set(gca,'XTick',[1 2 3 4],'XTickLabel',{'H','M','CR','FA'}); drawnow;
    
    Saccadetoleft = sum(Leftsaccadehists(ind,:),1);
    Saccadetoright = sum(Rightsaccadehists(ind,:),1);
    subplot(nrows,5,2+5*(ii-1)); bar(timebins,[Saccadetoleft' Saccadetoright']); ylabel('Count');
    set(gca,'Xlim',[0 0.3]);  drawnow;
    
    cL = LMScontrastL(ind,:);
    answers1 = answersL(ind,:);
    if size(cL,1)< size(cL,2)
        cL = cL';
    end
    contrastL = unique(cL);
    correctanswersL = zeros(size(contrastL));
    wronganswersL = zeros(size(contrastL));
    trialspercontrastL = zeros(size(contrastL));
    for jj = 1:numel(contrastL)
        trialspercontrastL(jj) = numel(answers1(cL==contrastL(jj)));
        correctanswersL(jj) = sum(answers1(cL==contrastL(jj)));
        wronganswersL(jj) = trialspercontrastL(jj) - correctanswersL(jj);
        percorrectL(jj) = correctanswersL(jj)/trialspercontrastL(jj);
        subplot(nrows,5,3+5*(ii-1)); plot(contrastL(jj),percorrectL(jj),'o','MarkerSize',trialspercontrastL(jj)/5,'LineWidth',1.0,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]); hold on;
    end
    [a,b,g] = weibullFit(contrastL,[correctanswersL wronganswersL],'mle');
    pts = linspace(min(contrastL)/1.3,max(contrastL)*1.3,101);
    plot(pts,g +(0.5-g).*exp(-((pts./a).^b)),'b','Linewidth',2); 
    set(gca,'Xscale','log'); xlabel('contrast'); ylabel('Perf'); set(gca,'Xlim',[10^(-1.0) 10^(-0.1)]);
    cNL = LMScontrastNL(ind,:);
    answers2 = answersNL(ind,:);
    if size(cNL,1)< size(cNL,2)
        cNL = cNL';
    end
    contrastNL = unique(cNL);
    correctanswersNL = zeros(size(contrastNL));
    wronganswersNL = zeros(size(contrastNL));
    trialspercontrastNL = zeros(size(contrastNL));
    for jj = 1:numel(contrastNL)
        trialspercontrastNL(jj) = numel(answers2(cNL==contrastNL(jj)));
        correctanswersNL(jj) = sum(answers2(cNL==contrastNL(jj)));
        wronganswersNL(jj) = trialspercontrastNL(jj) - correctanswersNL(jj);
        percorrectNL(jj) = correctanswersNL(jj)/trialspercontrastNL(jj);
        subplot(nrows,5,3+5*(ii-1)); plot(contrastNL(jj),percorrectNL(jj),'o','MarkerSize',trialspercontrastNL(jj)/5,'LineWidth',1.0,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
    end
    [a,b,g] = weibullFit(contrastNL,[correctanswersNL wronganswersNL],'mle');
    pts = linspace(min(contrastNL)/1.3,max(contrastNL)*1.3,101);
    plot(pts,g +(0.5-g).*exp(-((pts./a).^b)),'r','Linewidth',2); hold off;
    set(gca,'Xscale','log'); xlabel('contrast'); ylabel('Perf'); title(num2str(locationstested(ii,:))); set(gca,'Xlim',[10^(-1.0) 10^(-0.3)]);hold off;
    subplot(nrows,5,4+5*(ii-1)); plot(HitsL(ind)-HitsNL(ind),'b','Linewidth',2); hold on; plot(CRL(ind)-CRNL(ind),'g','Linewidth',2); xlabel('time'); ylabel('diff');
     set(gca,'Ylim',[-20 20]);
    
    % Represent the Hits, Miss, CR and FA in terms of gaussian signal and noise distribution
    [mu1L,mu2L,sigma2L] = calcsignalnoisedist(sum(HitsL(ind)),sum(MissL(ind)),sum(CRL(ind)),sum(FAL(ind)),decisioncriterion);
    [mu1NL,mu2NL,sigma2NL] = calcsignalnoisedist(sum(HitsNL(ind)),sum(MissNL(ind)),sum(CRNL(ind)),sum(FANL(ind)),decisioncriterion);
    subplot(nrows,5,5*(ii-1)+5); plot(xvals,0.5*normpdf(xvals,mu1L,1)/max(normpdf(xvals,mu1L,1)),'b','Linewidth',2); hold on; plot(xvals,0.5*normpdf(xvals,mu2L,sigma2L)/max(normpdf(xvals,mu2L,sigma2L)),'b','Linewidth',2);
    plot(xvals,0.5*normpdf(xvals,mu1NL,1)/max(normpdf(xvals,mu1NL,1)),'r','Linewidth',2); plot(xvals,0.5*normpdf(xvals,mu2NL,sigma2NL)/max(normpdf(xvals,mu2NL,sigma2NL)),'r','Linewidth',2); line([decisioncriterion decisioncriterion],[0 1.0]);
    xlabel('decision var'), ylabel('pd');
    
end
plot_counter = plot_counter + 1;

% Plotting the RF location wrt to fixation point
figure(plot_counter), set(gcf,'Name','RF location');
plot(stro.sum.exptParams.fp_x,stro.sum.exptParams.fp_y,'+','Markersize',15,'Linewidth',2);
hold on; plot(locationstested(:,1), locationstested(:,2),'o', 'Markersize',15,'Linewidth',2);
axis square; set(gca,'Xlim',[-100 100],'Ylim',[-100 100]); grid on; xlabel('X'), ylabel('Y'); hold off; plot_counter = plot_counter + 1;