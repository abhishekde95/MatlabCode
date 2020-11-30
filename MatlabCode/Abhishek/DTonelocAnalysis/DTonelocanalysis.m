% Script for analysing DToneloc data
% Author - Abhishek De, 3/18

close all;
clearvars;
% stro = nex2stro(findfile('M032518019.nex')); % actual opto file
% stro = nex2stro(findfile('M032518020.nex'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LMScontrast%%%%%%%%%%%%%%%%%%%%%%%%%
% stro = nex2stro(findfile('M032918001.nex')); % just behavior without opto files
% stro = nex2stro(findfile('M032918002.nex'));
% stro = nex2stro(findfile('M032918003.nex'));
% stro = nex2stro(findfile('M032918004.nex'));
% stro = nex2stro(findfile('M033018002.nex'));
% stro = nex2stro(findfile('M033018003.nex'));

% stro = nex2stro(findfile('M040318001.nex')); % Lum
% stro = nex2stro(findfile('M040318002.nex')); % L-M @20,-40
% stro = nex2stro(findfile('M040318003.nex')); % L-M @20,-50
% stro = nex2stro(findfile('M040318004.nex')); % L-M @20,-50
% stro = nex2stro(findfile('M040318005.nex')); % Lum @20,-60
% stro = nex2stro(findfile('M040318006.nex')); % L-M @20,-60
% stro = nex2stro(findfile('M040318007.nex')); % L+M+S @20,-60
% stro = nex2stro(findfile('M040418001.nex')); % L+M+S, @47,-67, staircase OFF
% stro = nex2stro(findfile('M040418002.nex')); % L+M+S, @47,-67, staircase ON
% stro = nex2stro(findfile('M040418003.nex')); % L+M+S, @47,-67, staircase ON
% stro = nex2stro(findfile('M040418004.nex')); % L+M+S, @47,-67, staircase OFF
% stro = nex2stro(findfile('M040518001.nex')); % luminance, staircase ON @50,-68
% stro = nex2stro(findfile('M040518002.nex')); % luminance, staircase ON
% stro = nex2stro(findfile('M040518003.nex')); % luminance, staircase ON
% stro = nex2stro(findfile('M040518004.nex')); % luminance, staircase ON
% stro = nex2stro(findfile('M040518005.nex')); % luminance, staircase ON
% stro = nex2stro(findfile('M040518006.nex')); % luminance, staircase ON
% stro = nex2stro(findfile('M040518007.nex')); % luminance, staircase ON

% stro = nex2stro(findfile('M040618001.nex')); % luminance, staircase ON
% stro = nex2stro(findfile('M040618002.nex')); % luminance, staircase ON
% stro = nex2stro(findfile('M040618003.nex')); % luminance, staircase ON
% stro = nex2stro(findfile('M040618004.nex')); % luminance, staircase ON
% stro = nex2stro(findfile('M040618005.nex')); % luminance, staircase ON
% stro = nex2stro(findfile('M040618006.nex')); % luminance, staircase ON
% stro = nex2stro(findfile('M040618007.nex')); % luminance, staircase ON

% stro = nex2stro(findfile('M040918001.nex')); % @ 34,-70
% stro = nex2stro(findfile('M040918002.nex')); % @ 50,-65
% stro = nex2stro(findfile('M040918003.nex')); % @ 60,-55
% stro = nex2stro(findfile('M040918004.nex')); % @ 20,-55

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% opto files -4/12
% stro = nex2stro(findfile('M041218001.nex'));
% stro = nex2stro(findfile('M041218003.nex'));
% stro = nex2stro(findfile('M041218004.nex'));
% stro = nex2stro(findfile('M041218006.nex'));
% stro = nex2stro(findfile('M041218007.nex'));
% stro = nex2stro(findfile('M041218008.nex'));
% stro = nex2stro(findfile('M041218009.nex'));

% -4/25
% stro = nex2stro(findfile('M042518002.nex')); % 400(stim dur),300(laser dur)
% stro = nex2stro(findfile('M042518004.nex')); % 200, 300
% stro = nex2stro(findfile('M042518005.nex')); % 200, 300
% stro = nex2stro(findfile('M042518006.nex')); % 200, 200
% stro = nex2stro(findfile('M042518007.nex')); % 200, 150
% stro = nex2stro(findfile('M042518008.nex')); % 200, 300

stro = nex2stro(findfile('M051718014.nex'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
% filename = fetch(conn,'SELECT filename FROM DToneloc');
% training = fetch(conn,'SELECT training FROM DToneloc');
% comments = fetch(conn,'SELECT comments FROM DToneloc');
% close (conn);

plot_counter = 1;
alpha = 0.05;
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
coldirname = {'Lum','L-M','S'};
lasertrials = logical(stro.trial(:,optstim));
N = numel(uniquediridxs);

for ii = 1:N
    
    LMSdirtripletslaser = LMStriplet(stimpresent & lasertrials & diridxs == ii-1,:);
    LMSdirtripletsnonlaser = LMStriplet(stimpresent & ~lasertrials & diridxs == ii-1,:);
    LMSdircontrastlaser = sqrt(sum(LMSdirtripletslaser.^2,2));
    LMSdircontrastnonlaser = sqrt(sum(LMSdirtripletsnonlaser.^2,2));
    
    oogidxslaser = logical(stro.trial(stimpresent & lasertrials & diridxs == ii-1,oog));
    oogidxsnonlaser = logical(stro.trial(stimpresent & ~lasertrials & diridxs == ii-1,oog));
    
    lasertrialnumbers = 1:numel(LMSdircontrastlaser);
    nonlasertrialnumbers = 1:numel(LMSdircontrastnonlaser);
    figure(plot_counter);subplot(2,N,ii);plot(lasertrialnumbers,LMSdircontrastlaser,'b-o','Linewidth',2); hold on;
    plot(lasertrialnumbers(oogidxslaser==1),LMSdircontrastlaser(oogidxslaser==1),'o','MarkerSize',5,'LineWidth',1.0,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]);
    xlabel('Trial number'); ylabel('contrast'); title(coldirname{ii}); hold off; drawnow;
    
    subplot(2,N,ii+N);plot(nonlasertrialnumbers,LMSdircontrastnonlaser,'r-o','Linewidth',2); hold on;
    plot(nonlasertrialnumbers(oogidxsnonlaser==1),LMSdircontrastnonlaser(oogidxsnonlaser==1),'o','MarkerSize',5,'LineWidth',1.0,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]);
    xlabel('Trial number'); ylabel('contrast'); title(coldirname{ii}); hold off; drawnow;
    
end
figure(plot_counter);set(gcf,'Name','Laser ON staircase: stim present + stim absent');
plot_counter = plot_counter + 1;

% Now I will combine stimpresent and stimabsent for each color direction
% and analyze the data
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
    TPRlasertrial = Hitlasertrial/(Hitlasertrial+FAlasertrial); % True positive ratio
    FPRlasertrial = Misslasertrial/(Misslasertrial+CRlasertrial); % False positive ratio
    
    % Non-Laser trials
    Hitnonlasertrial = sum(colordirchoiceidxs & ~lasertrialidxs & colorstimpresent); % Hit
    Missnonlasertrial = sum(~colordirchoiceidxs & ~lasertrialidxs & colorstimpresent); % Miss
    CRnonlasertrial = sum(colordirchoiceidxs & ~lasertrialidxs & ~colorstimpresent); % Correct Reject
    FAnonlasertrial = sum(~colordirchoiceidxs & ~lasertrialidxs & ~colorstimpresent); % False Alarm
    TPRnonlasertrial = Hitnonlasertrial/(Hitnonlasertrial+FAnonlasertrial); % True positive ratio
    FPRnonlasertrial = Missnonlasertrial/(Missnonlasertrial+CRnonlasertrial); % False positive ratio
    
    figure(plot_counter);subplot(1,N,ii); bar([percentcorrectlasertrials percentcorrectnonlasertrials]); ylabel('Per correct');
    title(coldirname{ii}); set(gca,'XTick',[1 2],'XTickLabel',{'L','NL'}); drawnow;
    
    figure(plot_counter+1);subplot(1,N,ii); bar([Hitlasertrial Hitnonlasertrial; Misslasertrial Missnonlasertrial; CRlasertrial CRnonlasertrial; FAlasertrial FAnonlasertrial]); ylabel('Count');
    title(coldirname{ii}); set(gca,'XTick',[1 2 3 4],'XTickLabel',{'H','M','CR','FA'}); legend('L','NL'); drawnow;
    
    figure(plot_counter+2);subplot(1,N,ii); bar([Misslasertrial+CRlasertrial Missnonlasertrial+CRnonlasertrial; Hitlasertrial+FAlasertrial Hitnonlasertrial+FAnonlasertrial]); ylabel('Count');
    title(coldirname{ii}); set(gca,'XTick',[1 2],'XTickLabel',{'Left','Right'}); legend('L','NL'); drawnow;
    
    figure(plot_counter+3); subplot(1,N,ii); plot([0 FPRlasertrial 1],[0 TPRlasertrial 1],'b','Linewidth',2); hold on;
    plot([0 FPRnonlasertrial 1],[0 TPRnonlasertrial 1],'r','Linewidth',2); xlabel('FPR'); ylabel('TPR'); title(coldirname{ii}); hold off;
    
end
figure(plot_counter);set(gcf,'Name','Percentage correct: Laser vs Non-Laser Trials (both stim present and stim absent)');
figure(plot_counter+1);set(gcf,'Name','Hit Miss CR FA: Laser vs Non-Laser');
figure(plot_counter+2);set(gcf,'Name','Leftward and Rightward choices: Laser vs Non-Laser');
figure(plot_counter+3);set(gcf,'Name','ROC Analysis');
plot_counter = plot_counter + 4;

% Looking for the right and wrong choices for stim absent trials,
% separately for laser and non-laser trials, This might be relevant for
% some choice probability related analysis
choiceslasertrialsstimabsent = correcttrials(~stimpresent & lasertrials);
choicesnonlasertrialsstimabsent = correcttrials(~stimpresent & ~lasertrials);
figure(plot_counter);
bar([sum(choiceslasertrialsstimabsent) sum(choicesnonlasertrialsstimabsent); sum(~choiceslasertrialsstimabsent) sum(~choicesnonlasertrialsstimabsent)]); 
set(gca,'XTick',[1 2],'XTickLabel',{'Left','Right'}); legend('L','NL'); title('Stim absent trials'); drawnow;
set(gcf,'Name','Choices: Stimulus absent trials');
plot_counter = plot_counter + 1;

% The next part is for some neural analysis, I specifically want to compare
% the firing rate of a single/multi-unit for the stimulus absent laser and
% non-laser trials
idxs = find(~stimpresent);
FRnolaser = [];
FRlaser = [];
count = 1;
figure(plot_counter);
for ii = 1:numel(idxs)
    ind = idxs(ii);
    analogstartime = stro.ras{ind,5};
    if lasertrials(ind)
%         keyboard;
        laserontime = stro.trial(ind,laseron);
        laserofftime = stro.trial(ind,laseroff);
        spiketimes = stro.ras{ind,1};
        timedurlaser = laserofftime - laserontime;
        FRlaser = [FRlaser; numel(spiketimes(spiketimes>laserontime & spiketimes<laserofftime))/timedurlaser];
        plot(spiketimes(spiketimes>laserontime-0.1 & spiketimes<laserofftime+0.1)-laserontime,count*ones(size(spiketimes(spiketimes>laserontime-0.1 & spiketimes<laserofftime+0.1))),'k.'); hold on;
        count = count + 1;
    end
    fpacqtime = stro.trial(ind,fpacqidx);
    fpofftime = stro.trial(ind,fpoffidx);
    spiketimes = stro.ras{ind,1};
    spiketimes = spiketimes(spiketimes>fpacqtime & spiketimes<fpofftime);
    timedur = fpofftime - fpacqtime;
    FRnolaser = [FRnolaser; numel(spiketimes)/timedur];
end
figure(plot_counter); line([0 0],[0 count]); line([timedurlaser timedurlaser],[0 count]); xlabel('time'); ylabel('trials'); title('Rasters for laser trials'); set(gca,'Xlim',[-0.1 0.4]); hold off;
set(gcf,'Name','Response: stim absent, laser trials'); plot_counter = plot_counter + 1;
figure(plot_counter); set(gcf,'Name','Histogram of baseline FR: stim absent trials'); subplot(121); hist(FRlaser,10); xlabel('FR'); ylabel('Frequency'); title('Laser(FR bt. laseroff & on)');
subplot(122); hist(FRnolaser,10); xlabel('FR'); ylabel('Frequency'); title('Non-laser(FR bt. fpoff & fpacq)'); 
plot_counter = plot_counter + 1;

% Trying to figure out if I can draw a psychometric function from the data points
color = ['b'; 'r'];
trials = lasertrials;
for kk = 1:2
    if kk == 2
        trials = ~trials;
    end
    for ii = 1:N
        LMScontrast = sqrt(sum(LMStriplet(stimpresent & trials & diridxs==ii-1,:).^2,2));
        putativeoogcontrast = min(LMScontrast(logical(stro.trial(stimpresent & trials & diridxs==ii-1,oog))));
        answers = correcttrials(stimpresent & trials & diridxs==ii-1);
        if ~isempty(putativeoogcontrast)
            LMScontrast(LMScontrast>putativeoogcontrast) = putativeoogcontrast;
        end
        contrast = unique(LMScontrast);
        correctanswers = zeros(size(contrast));
        wronganswers = zeros(size(contrast));
        trialspercontrast = zeros(size(contrast));
        
        figure(plot_counter); subplot(2,N,ii+(kk-1)*N);
        for jj = 1:numel(contrast)
            trialspercontrast(jj) = numel(answers(LMScontrast==contrast(jj)));
            correctanswers(jj) = sum(answers(LMScontrast==contrast(jj)));
            wronganswers(jj) = trialspercontrast(jj) - correctanswers(jj);
            percorrect(jj) = correctanswers(jj)/trialspercontrast(jj);
            plot(contrast(jj),percorrect(jj),'o','MarkerSize',trialspercontrast(jj),'LineWidth',1.0,'MarkerFaceColor',color(kk),'MarkerEdgeColor',color(kk)); hold on;
        end
        [a,b,g] = weibullFit(contrast,[correctanswers wronganswers],'mle');
        pts = linspace(min(contrast)/1.3,max(contrast)*1.3,101);
        plot(pts,g +(0.5-g).*exp(-((pts./a).^b)),color(kk),'Linewidth',2);
        set(gca,'Xscale','log'); xlabel('contrast'); ylabel('Perf'); title(coldirname{ii}); set(gca,'Xlim',[10^(-1.5) 10^(-0.1)]); hold off;
    end
end
set(gcf,'Name','Psychometric functions');
plot_counter = plot_counter + 1;

% Plotting the RF location wrt to fixation point
figure(plot_counter), set(gcf,'Name','RF location');
plot(stro.sum.exptParams.fp_x,stro.sum.exptParams.fp_y,'+','Markersize',15,'Linewidth',2);
hold on; plot(stro.sum.exptParams.rf_x, stro.sum.exptParams.rf_y,'o', 'Markersize',15,'Linewidth',2);
set(gca,'Xlim',[-100 100],'Ylim',[-100 100]); grid on; hold off; plot_counter = plot_counter + 1;

% Adding a new section, Want to analyze the saccade trajectories during
% Laser stim absent trials
stimabsentL_idxs = ~colorstimpresent & lasertrialidxs;
stimabsentNL_idxs = ~colorstimpresent & ~lasertrialidxs;
samplerate = stro.sum.analog.storeRates{1};
fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
figure(plot_counter); set(gcf,'Name','Saccade traj: stim absent L & NL');
for ii = 1:size(lasertrialidxs,1)
    if size(stro.ras,2)==5
        x = stro.ras{ii,2}*4096/400;
        y = stro.ras{ii,3}*4096/400;
        t = stro.ras{ii,5}+[0:1:length(x)-1]/samplerate;
    else
        x = stro.ras{ii,1}*4096/400;
        y = stro.ras{ii,2}*4096/400;
        t = stro.ras{ii,4}+[0:1:length(x)-1]/samplerate;
    end
    if stimabsentL_idxs(ii)
        Lt = t>fpoff_t(ii) & t < fpoff_t(ii)+0.30;
        subplot(121); plot(x(Lt),y(Lt),'b'); hold on;
    elseif stimabsentNL_idxs(ii)
        Lt = t>fpoff_t(ii) & t < fpoff_t(ii)+0.30;
        subplot(122); plot(x(Lt),y(Lt),'r'); hold on;
    end
end
subplot(121); xlabel('X degrees'); ylabel('Y degrees'); title('stim absent: Laser');  set(gca,'Xlim',[-5 5],'Ylim',[-5 5]); axis square; grid on; hold off;
subplot(122); xlabel('X degrees'); ylabel('Y degrees'); title('stim absent: Non-Laser');  set(gca,'Xlim',[-5 5],'Ylim',[-5 5]); axis square; grid on; hold off;
plot_counter = plot_counter + 1;

[h1,p1] = equalproptest([Hitlasertrial Hitnonlasertrial],[sum(stimpresent & lasertrials) sum(stimpresent & ~lasertrials)],alpha);
[h2,p2] = equalproptest([FAlasertrial FAnonlasertrial],[sum(~stimpresent & lasertrials) sum(~stimpresent & ~lasertrials)],alpha);

%% Looking at the eye positions for all the staircases
figure(plot_counter);
trials = lasertrials;
color = ['r';'k'];
samplingrate = stro.sum.analog.storeRates{2};
for kk = 1:2
    if kk == 2
        trials = ~trials;
    end
    for ii = 1:numel(uniquediridxs)
        idxs = find((trials & diridxs == ii-1)==1);
        stimsituation = stimpresent(idxs);
        for jj = 1:numel(idxs)
            analogstartime = stro.ras{idxs(jj),5};
            fpofftime = stro.trial(idxs(jj),fpoffidx);
            saccendtime = stro.trial(idxs(jj),saccendidx);
            timepts = analogstartime:(1/samplingrate):analogstartime + (numel(stro.ras{idxs(jj),2})-1)/samplingrate;
            timeptsofinterest = find(timepts>fpofftime & timepts<saccendtime);
            eye_x = stro.ras{idxs(jj),2};
            eye_y = stro.ras{idxs(jj),3};
            if stimsituation(jj)==1
                c = color(1);
            else
                c = color(2);
            end
            subplot(2,3,ii+(kk-1)*numel(uniquediridxs)),plot(eye_x(timeptsofinterest),eye_y(timeptsofinterest),c); hold on;
        end
        xlabel('x'),ylabel('y'); axis equal; title(coldirname{ii}); hold off;
    end 
end
plot_counter = plot_counter + 1;

%% Testing if the staircase is working fine
lasertrialidxs = logical(stro.trial(:,optstim));
startingcontrast = 1;
multiplicativefactor = 0.75;

contrasttrajectorylaser = [startingcontrast];
figure(plot_counter);
for ii = 1:numel(uniquediridxs)
    contrasttrajectorynonlaser = [startingcontrast];
    contrasttrajectorylaser = [startingcontrast];
    trialidxsnonlaser = ~lasertrialidxs & diridxs==ii-1;
    trialidxslaser = lasertrialidxs & diridxs==ii-1;
    trialidxsnonlaserstimpresent = stimpresent(trialidxsnonlaser);
    trialidxslaserstimpresent = stimpresent(trialidxslaser);
    monkeyanswersnonlaser = correcttrials(trialidxsnonlaser);
    monkeyanswerslaser = correcttrials(trialidxslaser);
    correctconsecutivecountnonlaser = 0;
    correctconsecutivecountlaser = 0;
    for jj = 1:numel(monkeyanswersnonlaser)
        if monkeyanswersnonlaser(jj) == 1
            correctconsecutivecountnonlaser = correctconsecutivecountnonlaser + 1;
            if correctconsecutivecountnonlaser == 3
                correctconsecutivecountnonlaser = 0;
                tmp = contrasttrajectorynonlaser(end)*multiplicativefactor;
            else
                tmp = contrasttrajectorynonlaser(end);
            end
        else
            correctconsecutivecountnonlaser = 0;
            tmp = contrasttrajectorynonlaser(end)/multiplicativefactor;
        end
        contrasttrajectorynonlaser = [contrasttrajectorynonlaser; tmp];
    end
    for jj = 1:numel(monkeyanswerslaser)
        if monkeyanswerslaser(jj) == 1
            correctconsecutivecountlaser = correctconsecutivecountlaser + 1;
            if correctconsecutivecountlaser == 3
                correctconsecutivecountlaser = 0;
                tmp = contrasttrajectorylaser(end)*multiplicativefactor;
            else
                tmp = contrasttrajectorylaser(end);
            end
        else
            correctconsecutivecountlaser = 0;
            tmp = contrasttrajectorylaser(end)/multiplicativefactor;
        end
        contrasttrajectorylaser = [contrasttrajectorylaser; tmp];
    end
    contrasttrajectorylaser(end) = [];
    subplot(2,3,ii+3); plot(contrasttrajectorynonlaser(trialidxsnonlaserstimpresent),'r','Linewidth',2); hold on;
    plot(contrasttrajectorynonlaser(trialidxsnonlaserstimpresent),'o','MarkerSize',5,'LineWidth',1.0,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
    xlabel('trial number'); ylabel('contrast'); title(coldirname{ii}); hold off;
    subplot(2,3,ii); plot(contrasttrajectorylaser(trialidxslaserstimpresent),'b','Linewidth',2); hold on;
    plot(contrasttrajectorylaser(trialidxslaserstimpresent),'o','MarkerSize',5,'LineWidth',1.0,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
    xlabel('trial number'); ylabel('contrast'); title(coldirname{ii}); hold off;
end
plot_counter = plot_counter + 1;
