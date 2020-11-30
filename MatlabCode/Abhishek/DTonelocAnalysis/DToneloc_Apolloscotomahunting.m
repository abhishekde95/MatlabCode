% A new script for analysing the leftover visual spaces where he can do the
% detection task
% Author - Abhishek De, 5/2/19
close all; clearvars;

% filename = ['M051418010.nex'; 'A021419009.nex'];
filename = ['A041919001.nex'; 'A041919002.nex'; 'A041919003.nex'; 'A041919004.nex'; 'A041919005.nex'; 'A041919006.nex'; 'A041919007.nex'; 'A041919008.nex'; 'A041919009.nex'; 'A041919010.nex';...
    'A042319001.nex'; 'A042319002.nex'; 'A042319004.nex'; 'A042319005.nex'; 'A042319006.nex'; 'A042319007.nex'; 'A042319008.nex';...
    'A042419001.nex'; 'A042419002.nex'; 'A042419003.nex'; 'A042419004.nex'; 'A042419005.nex'; 'A042419006.nex'; 'A042419007.nex'; 'A042419008.nex';...
    'A042619001.nex'; 'A042619002.nex'; 'A042619003.nex'; 'A042619004.nex'; 'A042619005.nex'; 'A042619006.nex'; 'A042619007.nex'; 'A042619008.nex'; 'A042619009.nex'; 'A042619010.nex'; 'A042619011.nex';...
    'A042919001.nex'; 'A042919002.nex'; 'A042919003.nex'; 'A042919004.nex'; 'A042919005.nex'; 'A042919006.nex'; 'A042919007.nex'; 'A042919008.nex'; 'A042919009.nex'; 'A042919010.nex'; 'A042919011.nex';...
    'A050119001.nex'; 'A050119002.nex'; 'A050119003.nex'; 'A050119004.nex'; 'A050119005.nex'; 'A050119006.nex'; 'A050119007.nex'; 'A050119008.nex'; 'A050119009.nex'; 'A050119010.nex';...
    'A050219002.nex'; 'A050219003.nex'; 'A050219004.nex'; 'A050219005.nex'; 'A050219006.nex'; 'A050219007.nex'; 'A050219008.nex'; 'A050219009.nex'; 'A050219010.nex';'A050219011.nex';'A050219012.nex';'A050219013.nex';'A050219014.nex';...
    'A050319001.nex'; 'A050319002.nex'; 'A050319003.nex'; 'A050319004.nex'; 'A050319005.nex'; 'A050319006.nex'; 'A050319007.nex'; 'A050319008.nex'; 'A050319009.nex'; 'A050319010.nex'];

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
weibullparams = [];
plot_counter = 1;
count = 1;
numplots = ceil(sqrt(size(filename,1)));
contrastlattice = logspace(log10(0.03),log(1.0),51);

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
    Contrast_luminance = ((inv(M)*((LMStriplet+1).*(repmat(Lb',[size(LMStriplet,1) 1])))')'*mon_spd*Vlambda-Bkgnd_luminance)/Bkgnd_luminance;
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
        
    oogvals = putativeoogcontrast;
    LMScontrast = [];
    LMScontrast = Contrast_luminance(stimpresent,:);
    answers = correcttrials(stimpresent);
    LMScontrast(LMScontrast>putativeoogcontrast) = putativeoogcontrast;
    contrastL = unique(LMScontrast);
    correctanswersL = zeros(size(contrastL));
    wronganswersL = zeros(size(contrastL));
    trialspercontrastL = zeros(size(contrastL));
    percorrectL = zeros(size(contrastL));
    for ss = 1:numel(contrastL)
        trialspercontrastL(ss) = numel(answers(LMScontrast==contrastL(ss)));
        correctanswersL(ss) = sum(answers(LMScontrast==contrastL(ss)));
        wronganswersL(ss) = trialspercontrastL(ss) - correctanswersL(ss);
        percorrectL(ss) = correctanswersL(ss)./(correctanswersL(ss)+wronganswersL(ss));
    end
    [aL,bL,gL] = weibullFitforDToneloc(contrastL,[correctanswersL wronganswersL],'mle');
    weibullparams = [weibullparams; aL bL gL]; % laser trials
    LMScontrastL = LMScontrast;
    RF = [RF; stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y];
    fitL = gL*(1-exp(-((contrastlattice./aL).^bL)));
    
    figure(plot_counter); subplot(numplots,numplots,jj); plot(contrastL,percorrectL,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]); hold on;
    plot(contrastlattice,fitL,'color',[0.5 0.5 0.5],'Linewidth',3); set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1]); title(num2str(RF(jj,:))); hold off;
    
    figure(plot_counter+1); subplot(numplots,numplots,jj); h = bar([sum(HitL(end))/(sum(HitL(end))+sum(MissL(end))) sum(HitNL(end))/(sum(HitNL(end))+sum(MissNL(end))); sum(MissL(end))/(sum(HitL(end))+sum(MissL(end))) sum(MissNL(end))/(sum(HitNL(end))+sum(MissNL(end))); sum(CRL(end))/(sum(CRL(end))+sum(FAL(end))) sum(CRNL(end))/(sum(CRNL(end))+sum(FANL(end))); sum(FAL(end))/(sum(CRL(end))+sum(FAL(end))) sum(FANL(end))/(sum(CRNL(end))+sum(FANL(end)))]); hold on; 
    set(h(2),'FaceColor',[0 0.5 1]); set(h(1),'FaceColor',[0.5 0.5 0.5]);
    set(gca,'XTick',[1 2 3 4],'Xlim',[0 5],'YTick',[0:0.5:1],'XTickLabel',{'H','M','CR','FA'},'TickDir','Out'); drawnow;

    
end
plot_counter = plot_counter + 2;
detthreshold = weibullparams(:,1);
detthreshold(detthreshold>putativeoogcontrast) = putativeoogcontrast;
parafoveal = sqrt(sum(RF.^2,2))/10 < 7;
figure(plot_counter); set(gcf,'Name','Detection thresholds');
subplot(121); st = tpaps(RF(parafoveal,:)',detthreshold(parafoveal)'); fnplt(st); hold on; 
 plot3(RF(parafoveal,1),RF(parafoveal,2),detthreshold(parafoveal),'o','MarkerFaceColor',[0 0 0]);
xlabel('X');ylabel('Y');zlabel('Luminance Contrast');
subplot(122); st = tpaps(RF(~parafoveal,:)',detthreshold(~parafoveal)'); fnplt(st); hold on; 
 plot3(RF(~parafoveal,1),RF(~parafoveal,2),detthreshold(~parafoveal),'o','MarkerFaceColor',[0 0 0]);
xlabel('X');ylabel('Y');zlabel('Luminance Contrast');
plot_counter = plot_counter + 1;