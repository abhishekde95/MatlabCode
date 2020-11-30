% Analyze SCstimcue behavior data to plot Maui's natural scotoma
% Author - Abhishek De, 4/18, SCstimcueAnalysesAD.m
clearvars; close all;
plot_counter = 2;
% filename  = {'M041618005.nex';'M042118002.nex';'M042118003.nex';'M042118005.nex';'M042118006.nex';'M042118007.nex';'M042318001.nex';'M042318002.nex';'M042318003.nex';'M042318004.nex';'M042318005.nex'};
% filename1 = {'M051018002.nex';'M051018003.nex';'M051018004.nex';'M051018005.nex'}; % 5/10, laser control trials
% filename2 = {'M051118001.nex';'M051118002.nex';'M051118003.nex';'M051118004.nex'}; % 5/11, laser control trials
filename3 = {'M051318006.nex';'M051318007.nex'}; %5/13, opto experiment
filename4 = {'M051418004.nex';'M051418005.nex';'M051418007.nex'}; %5/14, opto experiment
filename5 = {'M042518001.nex';'M042518003.nex'}; %4/25, opto experiment
filename6 = {'M051518001.nex';'M051518002.nex'}; %5/15, opto experiment
filename7 = {'M042618002.nex';'M042618008.nex'}; %4/26, opto experiment
filename8 = {'M042718004.nex';'M042718005.nex';'M042718011.nex'}; %4/27, opto experiment
filename9 = {'M051618001.nex';'M051618002.nex'}; %5/16, opto experiment
filename10 = {'M051718013.nex';'M051718017.nex'}; %5/17 opto experiment
filename11 = {'M051718020.nex'}; %5/17 opto experiment, saccade targets same as the DToneloc saccde targets
filename12 = {'M052218001.nex'}; %5/22 beh experiment, saccade targets same as the DToneloc saccde targets
filename13 = {'P053018001.nex';'P053018002.nex';'P053018003.nex';'P053018004.nex';'P053018005.nex'};
filename14 = {'A080718001.nex';'A080718002.nex';'A080718003.nex';'A080718004.nex';'A080718005.nex';'A080718006.nex';'A080718007.nex'}; %8/7 behavior Apollo
filename15 = {'A080818001.nex';'A080818002.nex';'A080818003.nex';'A080818004.nex';'A080818005.nex';'A0808180006'}; %8/8 behavior Apollo
filename16 = {'A080918001.nex';'A080918002.nex';'A080918003.nex';'A080918004.nex'}; %8/9 behavior Apollo 
filename17 = {'A081018001.nex';'A081018002.nex';'A081018003.nex';'A081018004.nex';'A081018005.nex';'A081018006.nex';'A081018007.nex'}; %8/10 behavior Apollo
filename18 = {'A081318001.nex';'A081318002.nex'};
filename19 = {'A081418001.nex';'A081418002.nex'};
filename20 = {'A081618001.nex'};
filename21 = {'P083018001.nex';'P083018002.nex';'P083018003.nex';'P083018004.nex';'P083018005.nex';'P083018006.nex';'P083018007.nex';'P083018008.nex';'P083018009.nex';'P083018010.nex';'P083018011.nex'};
% Files collected after Apollo's surgery
filename22 = {'A111918001.nex';'A111918002.nex';'A111918003.nex';'A111918004.nex';'A111918005.nex';'A111918006.nex'};
filename23 = {'A112018001.nex';'A112018002.nex';'A112018003.nex';'A112018004.nex';'A112018005.nex';'A112018006.nex';'A112018007.nex';'A112018008.nex'};
filename24 = {'A112118001.nex';'A112118002.nex';'A112118003.nex';'A112118004.nex';'A112118005.nex';'A112118006.nex';'A112118007.nex';'A112118008.nex'};
% opto sessions Apollo: Day 1 & 2
filename25 = {'A011419001.nex'; 'A011419003.nex'; 'A011519001.nex'; 'A011519002.nex'; 'A011519003.nex'; 'A011519004.nex'; 'A011519009.nex'; 'A011519010.nex'};
filename26 = {'A011719007.nex'; 'A011719008.nex'};

mode = 0;
if mode == 0
    filename = {'A020119007.nex'};
else
    filename =  filename26;
end
targetlocations = [];
targethitsallNL = []; % non-laser trials
stimpresentationNL = [];% non-laser trials
targethitsallL = []; % laser trials
stimpresentationL = []; % laser trials
halfwidth = 5.0; % in degrees of visual angle
RFloc = [];
eyetrajectories_catchL = {};
eyetrajectories_catchNL = {};
plot_individualfiles = 1;
for ii = 1:numel(filename)
    disp(ii);
    fileofinterest = char(filename(ii,:));
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
    for i = 1:size(stro.trial,1)
        if size(stro.ras,2)==5
            x = stro.ras{i,2}*4096/400;
            y = stro.ras{i,3}*4096/400;
            t = stro.ras{i,5}+[0:1:length(x)-1]/samplerate;
        else
            x = stro.ras{i,1}*4096/400;
            y = stro.ras{i,2}*4096/400;
            t = stro.ras{i,4}+[0:1:length(x)-1]/samplerate;
        end
        Lt = t>fpoff_t(i) & t < fpoff_t(i)+0.30;
        target = stro.trial(i,[12 13]); % x and y target locations of that trial
        idx = find(sum(target==C,2)==2);
        if ~Llaser(i)
            if sum(x(Lt)>(target(1)/10)-halfwidth & x(Lt)<(target(1)/10)+halfwidth & y(Lt)>(target(2)/10)-halfwidth & y(Lt)<(target(2)/10)+halfwidth)>10
                targethitsNL(idx) = targethitsNL(idx) + 1;
            end
            SPNL(idx) = SPNL(idx) + 1;
            if Lcatchtrials(i)
                eyetrajectories_catchNL = [eyetrajectories_catchNL; {x(Lt)} {y(Lt)}];
            else
                if ~isnan(sacmade_t(i))
                    Lsac = t>sacmade_t(i) & t<sacmade_t(i) + 0.040;
                    figure(1), plot(mean(x(Lsac)),mean(y(Lsac)),'b.'); hold on;
                end
            end
        end
        if Llaser(i)
            if sum(x(Lt)>(target(1)/10)-halfwidth & x(Lt)<(target(1)/10)+halfwidth & y(Lt)>(target(2)/10)-halfwidth & y(Lt)<(target(2)/10)+halfwidth)>10
                targethitsL(idx) = targethitsL(idx) + 1;
            end
            SPL(idx) = SPL(idx) + 1;
            if Lcatchtrials(i)
                eyetrajectories_catchL = [eyetrajectories_catchL; {x(Lt)} {y(Lt)}];
            else
                if ~isnan(sacmade_t(i))
                    Lsac = t>sacmade_t(i) & t<sacmade_t(i) + 0.040;
                    figure(1), plot(mean(x(Lsac)),mean(y(Lsac)),'r.'); hold on;
                end
            end
        end
    end
    targethitsallNL = [targethitsallNL; targethitsNL];
    targethitsallL = [targethitsallL; targethitsL];
    stimpresentationNL = [stimpresentationNL; SPNL];
    stimpresentationL = [stimpresentationL; SPL];
    if  plot_individualfiles
        timebins = 0:0.01:0.3;
        figure(plot_counter);
        nrows = ceil(sqrt(numel(unique(ib))));
        for jj = 1:numel(unique(ib))
            idx = logical(ib==jj);
            subplot(nrows,nrows,jj);histogram(sacint_t(idx)-fpoff_t(idx),timebins);
            xlabel('time'); title(num2str(targetlocations(jj,:)));
        end
        plot_counter = plot_counter + 1;
    end
end
plot_counter = plot_counter + 1;
% targetlocations = unique(targetlocations,'rows');
% Proportion of saccades to target locations during non-laser trials
Const = 7;
figure(plot_counter);set(gcf,'Name','Non-Laser trials');
plot(RFloc(:,1)/10,RFloc(:,2)/10,'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
for jj = 1:size(targethitsallNL)
    if targethitsallNL(jj)
        if sum(targetlocations(jj,:)~=0,2)
        plot(targetlocations(jj,1)/10,targetlocations(jj,2)/10,'o','MarkerSize',Const*targethitsallNL(jj)/stimpresentationNL(jj),'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]); hold on;
        end
    end
end
grid on; axis square; set(gca,'Xlim',[-17 17],'Ylim',[-17 17]); title('Non-laser trials'); xlabel('X degrees'); ylabel('Y degrees'); hold off;
plot_counter = plot_counter + 1;
% Proportion of saccades to target locations during laser trials
figure(plot_counter); set(gcf,'Name','Laser trials');
plot(RFloc(:,1)/10,RFloc(:,2)/10,'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
for jj = 1:size(targethitsallL)
    if targethitsallL(jj)
        if sum(targetlocations(jj,:)~=0,2)
        plot(targetlocations(jj,1)/10,targetlocations(jj,2)/10,'o','MarkerSize',Const*targethitsallL(jj)/stimpresentationL(jj),'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); hold on;
        end
    end
end
grid on; axis square; set(gca,'Xlim',[-17 17],'Ylim',[-17 17]); title('Laser trials'); xlabel('X degrees'); ylabel('Y degrees'); hold off;
plot_counter = plot_counter + 1;

% Now plotting the ratio of saccades during laser vs non-laser trials
tmptargethitsallL = targethitsallL;
tmptargethitsallNL = targethitsallNL;
figure(plot_counter);
for jj = 1:size(tmptargethitsallL)
    if tmptargethitsallL(jj)
        if sum(targetlocations(jj,:)~=0,2)
        plot(targetlocations(jj,1)/10,targetlocations(jj,2)/10,'o','MarkerSize',Const*tmptargethitsallL(jj)/(tmptargethitsallNL(jj)+tmptargethitsallL(jj)),'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]); hold on;
        end
    end
end
grid on; axis square; set(gca,'Xlim',[-17 17],'Ylim',[-17 17]); title('L/NL target hits ratio'); xlabel('X degrees'); ylabel('Y degrees'); hold off;
plot_counter = plot_counter + 1;

% Next, I want to analyze the eye trajectories during catch trials
figure(plot_counter); set(gcf,'Name','Catch trials eye trajectories');
for ii = 1:size(eyetrajectories_catchL,1)
    subplot(121); plot(eyetrajectories_catchL{ii,1},eyetrajectories_catchL{ii,2},'b'); hold on
end
title('Laser trials'); grid on; xlabel('X degrees');set(gca,'Xlim',[-17 17],'Ylim',[-17 17]);  ylabel('Y degrees'); axis square; hold off;
for jj = 1:size(eyetrajectories_catchNL,1)
    subplot(122); plot(eyetrajectories_catchNL{jj,1},eyetrajectories_catchNL{jj,2},'r'); hold on
end
title('No-laser trials'); grid on; xlabel('X degrees');set(gca,'Xlim',[-17 17],'Ylim',[-17 17]);  ylabel('Y degrees'); axis square; hold off; 
plot_counter = plot_counter + 1;