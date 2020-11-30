
% Copying this piece of code from Greg's FixStimAnalyses.m
% FixStim_contrastresponse analyses
% Author - Abhishek De, 1/18
close all; clearvars;
% stro = nex2stro(findfile('Junk010818001.nex')); % Junk file for testing 8 different contrast levels with max contrast 1.0
% stro = nex2stro(findfile('M010718005.nex')); % multiunit, excitatory
% stro = nex2stro(findfile('M010718017.nex')); % multiunit, excitatory, good data
stro = nex2stro(findfile('M010818019.nex')); % multiunit,suppressive, good data
% stro = nex2stro(findfile('M010818023.nex')); % single unit, suppressive
% stro = nex2stro(findfile('M010818024.nex')); % single unit, suppressive

% ******************Not so good files**********************
% stro = nex2stro(findfile('M010818016.nex')); % multiunit,hard to tell, excitatory?
% stro = nex2stro(findfile('M010818017.nex')); % multiunit, hard to tell, excitatory?
% *********************************************************
% stro = nex2stro(findfile('M032518012.nex'));

conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filename = fetch(conn,'SELECT filename FROM Fixstim_CR');
training = fetch(conn,'SELECT isolation FROM Fixstim_CR');
comments = fetch(conn,'SELECT comments FROM Fixstim_CR');
close(conn);


plot_counter = 1;
offset = [-.2 .75];
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
lasertraceind = strcmp(stro.sum.rasterCells(1,:),'AD13');
uniquecontrasts = unique(contrasts);
uniquepowers = unique(laser_power);
Lspikechans = strncmp(stro.sum.rasterCells,'sig0',4);

sync_t = targon_t;  % Ecode for alignment
if (any(Lspikechans))
    for whichspike = find(Lspikechans)
        figure(plot_counter);
        spikes = stro.ras(:,whichspike);
        binwidth = .005;
        bins = offset(1):binwidth:offset(2);
        PSTH = zeros(1,length(bins));
        for j = uniquecontrasts'
            peakPSTH = 0;
            for i = uniquepowers'
                PSTH = zeros(1,length(bins));
                subplot(numel(uniquecontrasts),4,4*(find(uniquecontrasts==j)-1)+2*find(uniquepowers==i)-1); hold on;
                L = contrasts == j & laser_power == i;
                trlidxs = find(L);
                targetinterval = [];
                for counter = 1:sum(L)
                    trlidx = trlidxs(counter);
                    tmpspikes = spikes{trlidx}-sync_t(trlidx);
                    targetinterval = [targetinterval; targoff_t(trlidx)-sync_t(trlidx)];
                    tmpspikes(tmpspikes < offset(1) | tmpspikes > offset(2)) = [];
                    nspikestot = length(tmpspikes);
                    plot([tmpspikes tmpspikes]',[zeros(nspikestot,1) .5*ones(nspikestot,1)]'+counter,'k-');
                    PSTH = PSTH + hist(tmpspikes, bins);
                end
                set(gca,'XLim', offset,'Ytick',[],'YLim',[0 sum(L)+1]);
                line([0 0],[0 sum(L)+1],'Color',[1 0 0]); %
                line([mean(targetinterval) mean(targetinterval)],[0 sum(L)+1],'Color',[1 0 0]); %
                title(['Contrast: ',num2str(j),' Power: ',num2str(i)]);
                
                % PSTH
                subplot(numel(uniquecontrasts),4,4*(find(uniquecontrasts==j)-1)+2*find(uniquepowers==i)); hold on;
                plot(bins,PSTH,'k-','LineWidth',2);
                set(gca,'YLim',[0 10]);
                set(gca,'Xlim',offset);
                line([0 0],[0 10],'Color',[1 0 0]);
                line([mean(targetinterval) mean(targetinterval)],[0 10],'Color',[1 0 0]);
                peakPSTH = max([peakPSTH, PSTH]);
            end
        end
        plot_counter = plot_counter + 1;
    end
end

% Plotting the laser Trials
samplingrate = stro.sum.analog.storeRates{3};
multiplicationfactor = 0.2;
spikecounts = cell(1,numel(uniquepowers));
for jj = 1:numel(uniquepowers)
    figure(plot_counter);
    if uniquepowers(jj) > 0
        set(gcf,'Name','Laser Trials');
    else
        set(gcf,'Name','Non-laser Trials');
        
    end
    L = find(laser_power==uniquepowers(jj));
    laserinterval = [];
    targetinterval = [];
    tmp_spikecounts = [];
    for ii = 2:numel(L)
        ind = L(ii);
        reference = targon_t(ind);
        analogstarttime = stro.ras{ind,analogstrtimeind};
        spiketimes = stro.ras{ind,spikeind};
        tmp_spikecounts = [tmp_spikecounts; numel(find((spiketimes>targon_t(ind) & spiketimes<targoff_t(ind))>0))];
        time = linspace(analogstarttime,(analogstarttime + numel(stro.ras{ind,lasertraceind})/samplingrate),numel(stro.ras{ind,lasertraceind}));
        spiketimes = spiketimes - reference;
        time = time - reference;
        laserinterval = [laserinterval; stimoff_t(ind)-stimon_t(ind)];
        targetinterval = [targetinterval; targoff_t(ind)-targon_t(ind)];
        subplot(121); plot(time,multiplicationfactor*ii + stro.ras{ind,lasertraceind},'k'); hold on;
        subplot(122); plot(spiketimes,multiplicationfactor*ii + ones(size(spiketimes)),'k.'); hold on;
    end
    spikecounts{jj} = tmp_spikecounts;
    subplot(121);title('Analog Laser Traces'); xlabel('Time'), ylabel('Trial Number');
    line([0 0],[0 multiplicationfactor*(ii+1)],'Color',[0 1 0]); line([mean(targetinterval) mean(targetinterval)],[0 multiplicationfactor*(ii+1)],'Color',[0 1 0]); hold off;
    if uniquepowers(jj) > 0
        line([nanmean(stimon_t-targon_t) nanmean(stimon_t-targon_t)],[0 multiplicationfactor*(ii+1)],'Color',[0 0 1]); line([nanmean(stimon_t-targon_t)+mean(laserinterval) nanmean(stimon_t-targon_t)+mean(laserinterval)],[0 multiplicationfactor*(ii+1)],'Color',[0 0 1]); hold off;
    end
    subplot(122);title('Rasters'); xlabel('Time'), ylabel('Trial Number');
    line([0 0],[0 multiplicationfactor*(ii+1)],'Color',[0 1 0]); line([mean(targetinterval) mean(targetinterval)],[0 multiplicationfactor*(ii+1)],'Color',[0 1 0]); hold off;
    if uniquepowers(jj) > 0
        line([nanmean(stimon_t-targon_t) nanmean(stimon_t-targon_t)],[0 multiplicationfactor*(ii+1)],'Color',[0 0 1]); line([nanmean(stimon_t-targon_t)+mean(laserinterval) nanmean(stimon_t-targon_t)+mean(laserinterval)],[0 multiplicationfactor*(ii+1)],'Color',[0 0 1]); hold off;
    end
    drawnow; plot_counter = plot_counter + 1;
end

% Plotting a spikecounts for laser vs non-laser trials
figure(plot_counter), set(gcf,'Name','Laser vs Non-laser trials');
histogram(spikecounts{1}); hold on; histogram(spikecounts{2}); legend('Non-Laser','Laser');  
xlabel('spike counts'); ylabel('Num Trials');  hold off;

ROCval = roc(spikecounts{1},spikecounts{2});
if ROCval<0.5
    ROCval = 1- ROCval;
end
disp(ROCval);