% DTonelocAnalysismultiplelocations - for analysing files with interleaved staircases at
% different spatial locations
% Author - Abhishek De, 5/18
close all; clearvars;
% stro = nex2stro(findfile('M052118001.nex'));
% stro = nex2stro(findfile('M052118002.nex'));
% stro = nex2stro(findfile('M052118003.nex'));
% stro = nex2stro(findfile('M052118004.nex'));
% stro = nex2stro(findfile('M052218002.nex'));
% stro = nex2stro(findfile('M052218003.nex'));
% stro = nex2stro(findfile('M052218004.nex'));
% stro = nex2stro(findfile('M052218005.nex'));
% stro = nex2stro(findfile('M052318001.nex')); % actual opto file
stro = nex2stro(findfile('M052318002.nex')); % actual opto file
% stro = nex2stro(findfile('M052318003.nex')); % actual opto file
% stro = nex2stro(findfile('M052318004.nex')); % actual opto file
% stro = nex2stro(findfile('M052518007.nex')); % actual opto file
% stro = nex2stro(findfile('M052518009.nex')); % actual opto file
plot_counter = 1;
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
%
correcttrials = stro.trial(:,correct);
stimpresent = logical(stro.trial(:,stimpresentidx));
LMStriplet = [stro.trial(:,Lcc) stro.trial(:,Mcc) stro.trial(:,Scc)];
lasertrials = logical(stro.trial(:,optstim));
stimlocationsall = [stro.trial(:,stim_x) stro.trial(:,stim_y)];
[stimlocations,ia,ib] = unique(stimlocationsall,'rows');

nrows = numel(unique(ib));
HitsL = []; HitsNL = [];
MissL = []; MissNL = [];
CRL = []; CRNL = [];
FAL = []; FANL = [];
decisioncriterion = 0; % Am assuming that the ideal observer always has a decision criterion placed at x=0, x>0, look right, x<0, look left
color = ['b'; 'r'];
timebins = 0:0.01:0.3;
xvals = linspace(-3,3,101);
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
    
    figure(plot_counter);subplot(nrows,5,5*(ii-1)+1);plot(lasertrialnumbers,LMSdircontrastlaser,'b','Linewidth',2); hold on;
    plot(lasertrialnumbers(oogidxslaser==1),LMSdircontrastlaser(oogidxslaser==1),'o','MarkerSize',5,'LineWidth',1.0,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]); hold on;
    plot(nonlasertrialnumbers,LMSdircontrastnonlaser,'r','Linewidth',2); hold on;
    plot(nonlasertrialnumbers(oogidxsnonlaser==1),LMSdircontrastnonlaser(oogidxsnonlaser==1),'o','MarkerSize',5,'LineWidth',1.0,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]);
    xlabel('Trial number'); ylabel('contrast'); title(num2str(stimlocations(ii,:))); hold off; drawnow;
    
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
    TPRlasertrial = Hitlasertrial/(Hitlasertrial+FAlasertrial); % True positive ratio
    FPRlasertrial = Misslasertrial/(Misslasertrial+CRlasertrial); % False positive ratio
    HitsL = [HitsL; Hitlasertrial];
    MissL = [MissL; Misslasertrial];
    CRL = [CRL; CRlasertrial];
    FAL = [FAL; FAlasertrial];    
    
    % Non-Laser trials
    Hitnonlasertrial = sum(colordirchoiceidxs & ~lasertrialidxs & colorstimpresent); % Hit
    Missnonlasertrial = sum(~colordirchoiceidxs & ~lasertrialidxs & colorstimpresent); % Miss
    CRnonlasertrial = sum(colordirchoiceidxs & ~lasertrialidxs & ~colorstimpresent); % Correct Reject
    FAnonlasertrial = sum(~colordirchoiceidxs & ~lasertrialidxs & ~colorstimpresent); % False Alarm
    TPRnonlasertrial = Hitnonlasertrial/(Hitnonlasertrial+FAnonlasertrial); % True positive ratio
    FPRnonlasertrial = Missnonlasertrial/(Missnonlasertrial+CRnonlasertrial); % False positive ratio
    HitsNL = [HitsNL; Hitnonlasertrial];
    MissNL = [MissNL; Missnonlasertrial];
    CRNL = [CRNL; CRnonlasertrial];
    FANL = [FANL; FAnonlasertrial]; 
    
    
    subplot(nrows,5,5*(ii-1)+2); bar([Hitlasertrial Hitnonlasertrial; Misslasertrial Missnonlasertrial; CRlasertrial CRnonlasertrial; FAlasertrial FAnonlasertrial]); ylabel('Count');
    title(num2str(stimlocations(ii,:))); set(gca,'XTick',[1 2 3 4],'XTickLabel',{'H','M','CR','FA'}); drawnow;
    
    
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
        contrast = unique(LMScontrast);
        correctanswers = zeros(size(contrast));
        wronganswers = zeros(size(contrast));
        trialspercontrast = zeros(size(contrast));
        
        subplot(nrows,5,(ii-1)*5+3);
        for jj = 1:numel(contrast)
            trialspercontrast(jj) = numel(answers(LMScontrast==contrast(jj)));
            correctanswers(jj) = sum(answers(LMScontrast==contrast(jj)));
            wronganswers(jj) = trialspercontrast(jj) - correctanswers(jj);
            percorrect(jj) = correctanswers(jj)/trialspercontrast(jj);
            plot(contrast(jj),percorrect(jj),'o','MarkerSize',trialspercontrast(jj)/2,'LineWidth',1.0,'MarkerFaceColor',color(kk),'MarkerEdgeColor',color(kk)); hold on;
        end
        [a,b,g] = weibullFit(contrast,[correctanswers wronganswers],'mle');
        pts = linspace(min(contrast)/1.3,max(contrast)*1.3,101);
        plot(pts,g +(0.5-g).*exp(-((pts./a).^b)),color(kk),'Linewidth',2);     
    end
    set(gca,'Xscale','log'); xlabel('contrast'); ylabel('Perf'); title(num2str(stimlocations(ii,:))); set(gca,'Xlim',[10^(-1.0) 10^(-0.3)]); hold off;
    
    % Need to calculate the saccade RT for the left and right targets
    Leftsaccadeind = (~stimpresent & correcttrials & ind) |  (stimpresent & ~correcttrials & ind);
    Rightsaccadeind = (stimpresent & correcttrials & ind) |  (~stimpresent & ~correcttrials & ind);
    subplot(nrows,5,5*(ii-1)+4); histogram(stro.trial(Leftsaccadeind,saccstartidx)-stro.trial(Leftsaccadeind,fpoffidx),timebins);hold on; 
    histogram(stro.trial(Rightsaccadeind,saccstartidx)-stro.trial(Rightsaccadeind,fpoffidx),timebins);
    set(gca,'Xlim',[0 0.3]); xlabel('RT');      
    % Represent the Hits, Miss, CR and FA in terms of gaussian signal and noise distribution
    [mu1L,mu2L,sigma2L] = calcsignalnoisedist(HitsL(end),MissL(end),CRL(end),FAL(end),decisioncriterion);
    [mu1NL,mu2NL,sigma2NL] = calcsignalnoisedist(HitsNL(end),MissNL(end),CRNL(end),FANL(end),decisioncriterion);
    subplot(nrows,5,5*(ii-1)+5); plot(xvals,0.5*normpdf(xvals,mu1L,1)/max(normpdf(xvals,mu1L,1)),'b','Linewidth',2); hold on; plot(xvals,0.5*normpdf(xvals,mu2L,sigma2L)/max(normpdf(xvals,mu2L,sigma2L)),'b','Linewidth',2);
    plot(xvals,0.5*normpdf(xvals,mu1NL,1)/max(normpdf(xvals,mu1NL,1)),'r','Linewidth',2); plot(xvals,0.5*normpdf(xvals,mu2NL,sigma2NL)/max(normpdf(xvals,mu2NL,sigma2NL)),'r','Linewidth',2); line([decisioncriterion decisioncriterion],[0 1.0]);
    xlabel('decision var'), ylabel('pd');  title(num2str(stimlocations(ii,:)));
    
    % Now I want to test if the noise distribution in laser trials could be accounted by just a shift or a shift + change in variance 
    out = comparemodelsnoisedist(CRL(end),FAL(end),CRNL(end),FANL(end),decisioncriterion);
    disp(out);
    
end

plot_counter = plot_counter + 1;

% Plotting the RF location wrt to fixation point
figure(plot_counter), set(gcf,'Name','RF location');
plot(stro.sum.exptParams.fp_x,stro.sum.exptParams.fp_y,'+','Markersize',15,'Linewidth',2);
hold on; plot(stimlocations(:,1), stimlocations(:,2),'o', 'Markersize',15,'Linewidth',2);
set(gca,'Xlim',[-100 100],'Ylim',[-100 100]); grid on; xlabel('X'), ylabel('Y'); axis equal; hold off; plot_counter = plot_counter + 1;