% Analysis to investigate it the psychometric function drifts over time, kinda like analysis in Fetch et al; 2018
% Author - Abhishek De, 3/19
close all; clearvars;
load('T_vos1978_Y');
Vlambda = T_vos1978_Y';
filenameopto1 = {'M042518004.nex';'M042518005.nex';'M042518008.nex';'M042518006.nex';'M042518007.nex'}; % 4/25
RWmultiplieropto1 =[1.0;1.0;1.0;1.0;1.0];
laserdialopto1 = [2.5;2.5;2.5;2.5;2.5];

filenameopto2 = {'M042618003.nex';'M042618004.nex';'M042618005.nex';'M042618006.nex';'M042618007.nex'}; % 4/26, 200 stim, 300 laser
RWmultiplieropto2 =[1.0;1.0;1.0;1.0;1.0];
laserdialopto2 = [2.5;2.5;2.5;2.5;2.5];

filenameopto3 = {'M051418008.nex';'M051418009.nex';'M051418010.nex';'M051418011.nex';'M051418013.nex';'M051418014.nex';'M051418015.nex';'M051418016.nex'}; %5/14, 200 stim, 300 laser 
RWmultiplieropto3 =[1.0;1.0;0.7;0.7;1.0;1.0;0.7;0.7];
laserdialopto3 = 4.0*ones(size(RWmultiplieropto3));

filenameopto4 = {'M051518003.nex';'M051518004.nex';'M051518005.nex';'M051518006.nex';'M051518007.nex';'M051518008.nex';'M051518009.nex';'M051518010.nex';'M051518011.nex'}; %5/15, 200 stim, 300 laser 
RWmultiplieropto4 =[1.0;0.7;0.5;0.5;0.3;0.6;1.0;0.4;0.7];
laserdialopto4 = 4.0*ones(size(RWmultiplieropto4));

filenameopto5 = {'A020119001.nex';'A020119002.nex';'A020119003.nex';'A020119004.nex';'A020119005.nex';'A020119006.nex';'A020119009.nex';'A020119010.nex';'A020119011.nex';'A020119012.nex';'A020119013.nex'};
RWmultiplieropto5 = [0.8;0.8;0.8;0.8;1.0;0.7;0.8;0.8;0.7;0.8;0.7];
laserdialopto5 = [2.0;2.0;2.0;2.0;2.0;2.0;2.0;2.0;2.0;2.0;0.2];

filenameopto6 = {'A012319007.nex';'A012319008.nex';'A012319009.nex';'A012319010.nex';'A012319011.nex'};
RWmultiplieropto6 = [0.7;0.7;0.7;0.7;0.7];
laserdialopto6 = [1.5;1.5;1.5;1.5;1.5];

filenameopto7 = {'A013019006.nex';'A013019007.nex';'A013019008.nex';'A013019009.nex';'A013019012.nex'};
RWmultiplieropto7 = [0.7;0.7;0.5;0.5;0.7];
laserdialopto7 = [2.2;2.2;2.2;2.2;2.2];

filenameopto8 = {'A020519002.nex';'A020519003.nex';'A020519004.nex';'A020519005.nex';'A020519009.nex'};
RWmultiplieropto8 = [0.7;0.7;0.8;1.0;0.8];
laserdialopto8 = [1.5;1.5;1.5;1.5;1.5];

oogvals = [];
plot_counter = 1;
nrows = 8;
figure(plot_counter);
RFs = cell(nrows,1);
diffinpropshits_early = cell(nrows,1); % 1-480 trials
diffinpropshits_late = cell(nrows,1); % 480+ trials
diffinpropshits = zeros(nrows,11);
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
        LMScontrastfileNL = Contrast_luminance(colorstimpresent & ~lasertrialidxs,:);
        LMScontrastfileL = Contrast_luminance(colorstimpresent & lasertrialidxs,:);
        % calculating gamut edge luminance contrast
        t = min(((1./stro.sum.exptParams.bkgndrgb)-1));
        OOG_luminance = (stro.sum.exptParams.bkgndrgb*(1+t))' * mon_spd*Vlambda;
        putativeoogcontrast =  (OOG_luminance - Bkgnd_luminance)/Bkgnd_luminance;
        
        trials = lasertrials;
        % only evaluating the laser trials
        for kk = 1:2
            if kk == 2
                trials = ~trials;
            end
            LMScontrast = [];
            LMScontrast = Contrast_luminance(stimpresent & trials,:);
            answers = correcttrials(stimpresent & trials);
            LMScontrast(LMScontrast>putativeoogcontrast) = putativeoogcontrast;
            oogvals = [oogvals; putativeoogcontrast];
            
            
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
            figure(1); subplot(3,3,mm);  plot(contrastL(jj),correctanswersL(jj)./(correctanswersL(jj)+wronganswersL(jj)),'o','MarkerSize',correctanswersL(jj) + wronganswersL(jj),'MarkerFacecolor',[1-cgradations(aa) 1-cgradations(aa) 1-cgradations(aa)],'MarkerEdgeColor',[1 1 1]); hold on;
        end
        plot(contrastlattice,fitL,'-','Linewidth',2,'color',[1-cgradations(aa) 1-cgradations(aa) 1-cgradations(aa)]);
        figure(2); subplot(3,3,mm); plot(LMScontrastL,'-','Linewidth',2,'color',[1-cgradations(aa) 1-cgradations(aa) 1-cgradations(aa)]); hold on;
                figure(4); subplot(3,3,mm); plot(LMScontrastNL,'-','Linewidth',2,'color',[1-cgradations(aa) 1-cgradations(aa) 1-cgradations(aa)]); hold on;
        fitNL = gNL*(1-exp(-((contrastlattice./aNL).^bNL)));
        for jj = 1:numel(contrastNL)
            figure(5); subplot(3,3,mm);  plot(contrastNL(jj),correctanswersNL(jj)./(correctanswersNL(jj)+wronganswersNL(jj)),'o','MarkerSize',correctanswersNL(jj) + wronganswersNL(jj),'MarkerFacecolor',[1-cgradations(aa) 1-cgradations(aa) 1-cgradations(aa)],'MarkerEdgeColor',[1 1 1]); hold on;
        end
        plot(contrastlattice,fitNL,'-','Linewidth',2,'color',[1-cgradations(aa) 1-cgradations(aa) 1-cgradations(aa)]);
        diffinpropshits(mm,aa) = HitNL(end)./stimpresentNL(end)-HitL(end)./stimpresentL(end);
    end
    RFs{mm} = RF;
    diffinpropshits_late{mm} = mean(HitNL(5:end)./stimpresentNL(5:end)-HitL(5:end)./stimpresentL(5:end));
    diffinpropshits_early{mm} = mean(HitNL(1:4)./stimpresentNL(1:4)-HitL(1:4)./stimpresentL(1:4));
    
    
    figure(1); subplot(3,3,mm); hold on; xlabel('contrast'); ylabel('prop of correct'); set(gca,'Tickdir','out','Ylim',[0 1],'Xlim',[0 1.3]); title(strcat('SessionL:',num2str(mm))); hold off;
    figure(2); subplot(3,3,mm); hold on; xlabel('trial no.'); ylabel('contrast'); set(gca,'Tickdir','out','Ylim',[0 0.7]); 
    if mm == 1
       title('laser staircase');
    end
    hold off;
    figure(3); subplot(3,3,mm); figure(3); subplot(3,3,mm); plot(HitNL./stimpresentNL-HitL./stimpresentL,'-o','Linewidth',2,'MarkerFaceColor',[0 0 0]); hold on;
    xlabel('file no.'); ylabel('phits control-laser'); set(gca,'Tickdir','out','Ylim',[0 1]); hold off;
    figure(4); subplot(3,3,mm); hold on; xlabel('trial no.'); ylabel('contrast'); set(gca,'Tickdir','out','Ylim',[0 0.7]); 
    if mm == 1
       title('control staircase'); 
    end
    hold off;
    figure(5); subplot(3,3,mm); hold on; xlabel('contrast'); ylabel('prop of correct'); set(gca,'Tickdir','out','Ylim',[0 1],'Xlim',[0 1.3]); title(strcat('SessionNL:',num2str(mm))); hold off;
    
end
plot_counter = plot_counter + 5;

std_fordiffs = [];
for ii = 1:size(diffinpropshits,2)
    idx = diffinpropshits(:,ii)>0;
    std_fordiffs = [std_fordiffs; std(diffinpropshits(idx,ii))/sqrt(sum(idx))]; 
end


figure(plot_counter); subplot(221);plot(cell2mat(diffinpropshits_early),cell2mat(diffinpropshits_late),'o','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot([0 1],[0 1],'k'); set(gca,'Xlim',[0 1],'Ylim',[0 1],'Tickdir','out','XTick',0:0.25:1,'YTick',0:0.25:1); axis square; xlabel('control-laser 1-480 trials'); ylabel('control-laser 480+ trials');
subplot(222);
for jj = 1:nrows
    text(cell2mat(diffinpropshits_early(jj)),cell2mat(diffinpropshits_late(jj)),num2str(jj)); hold on;
end
plot([0 1],[0 1],'k'); set(gca,'Xlim',[0 1],'Ylim',[0 1],'Tickdir','out','XTick',0:0.25:1,'YTick',0:0.25:1); axis square; xlabel('control-laser 1-480 trials'); ylabel('control-laser 480+ trials');
subplot(223); histogram(cell2mat(diffinpropshits_early)-cell2mat(diffinpropshits_late),-0.3:0.05:0.3,'FaceColor',[0 0 0]); hold on; 
plot(median(cell2mat(diffinpropshits_early)-cell2mat(diffinpropshits_late)),0,'kv','MarkerFaceColor',[0 1 0]); xlabel('early-late'); set(gca,'Tickdir','out'); axis square; hold off;
subplot(224); errorbar(sum(diffinpropshits)./sum(diffinpropshits>0),std_fordiffs,'-o','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); axis square; xlabel('trial no.'); set(gca,'Ylim',[0 1],'Tickdir','out','Xlim',[1 11],'XTick',[1:2:11]);
plot_counter = plot_counter + 1;


