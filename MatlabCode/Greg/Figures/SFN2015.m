% Figures for the SFN 2015 poster
% Optogenetic Projection Targeting in the Primate Ocuolomotor System
%
% Contents
% 1) Rasters and PSTHs from SC stimulation experiments.
%
% 2) Magnitude of optogenetic effect as a function of retinal position
%
% 3) Behavioral effect (or lack thereof) from optostim
%%
% Section 1
% Rasters from SC stimulation + PSTHs
% Largely taken from FixStimAnalyses, Section 5
INSET = 0;
stro = nex2stro(findfile('S013013011'));
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

% There was a brief period of time when you couldn't determine if the modulator was used or not
% (the modulator header parameter didn't exist). Account for that below.
[~, toks] = isvalidnexfilename(stro.sum.fileName);
file_date = datenum(str2double([['20' toks{1}{4}] toks{1}(2) toks{1}(3)]));
if file_date >= datenum([2015 5 7]) && file_date <= datenum([2015 5 14])
    modulator = 1;
elseif file_date < datenum([2015 5 7])
    modulator = 0;
else % after May 14, 2015
    modulator = stro.sum.exptParams.modulator;
end

% Hack below 
if (isnan(targ_shown(1)))
    targ_shown(:,1) = 0;
end
dur = mode(stimoff_t-stimon_t);

sync_t = stimon_t;  % Ecode for alignment

spikes = stro.ras(:,strcmp(stro.sum.rasterCells,'sig001a'));
binwidth = .002; % s
PSTHstd = .003; % s
bins = offset(1):binwidth:dur+offset(2);
PSTH = zeros(1,length(bins));
if (INSET)
    uniqfreqs = 250;
    binwidth = .00002; 
    PSTHstd = .00002; % For inset
end
for j = uniqfreqs'
    PSTH = zeros(1,length(bins));
    figure; subplot(2,1,1); hold on;
    L = optfreq == j;
    trlidxs = find(L);
    for counter = 1:sum(L)
        trlidx = trlidxs(counter);
        tmpspikes = spikes{trlidx}-sync_t(trlidx);
        tmpspikes(tmpspikes < offset(1) | tmpspikes > dur+offset(2)) = [];
        nspikestot = length(tmpspikes);
        plot(tmpspikes,zeros(nspikestot,1)+counter/1.5,'k.','MarkerSize',72*.15); % Inches to points conversion
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
            plot(t,y,'k-','linewidth',2);
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
    if (INSET)
        set(gca,'Xlim',[-.01 .04]);
    end
    
    % PSTH
    subplot(2,1,2); hold on;
    %plot(bins,PSTH,'k-','LineWidth',2);
    kernel = exppdf(linspace(0,2,3*PSTHstd/binwidth),1);
    kernel = kernel./sum(kernel);
    plot(bins,conv(PSTH,kernel,'same'),'k-','linewidth',2);
    set(gca,'YLim',[0 10*ceil(max(PSTH(:)/10))+1]);
    set(gca,'Xlim',[offset(1) dur+offset(2)]);
    xlabel('Time (s)','FontSize',12);
    ylabel('Response (sp/s)','FontSize',12);
    if (INSET)
        set(gca,'Xlim',[-.01 .04]);
    end
end


%%
% Section 2
% Magnitude of the optogenetic effect as a function of RF position

SCALEPOINTS = 1; % 1 = scale points to magnitude of optogenetic effect
[fnames, spikeIdx] = fnamesFromTxt2('FixStim.txt');
data = [];
for filecounter = 1:size(fnames,1)
    filename = findfile(char(fnames{filecounter}));
    paradigmID = getparadigmID(filename);
    if paradigmID ~= 105
        error([fnames{counter},' is not a FixStim file']);
    else
        stro = nex2stro(filename);
        rfx = stro.sum.exptParams.rf_x;
        rfy = stro.sum.exptParams.rf_y;
        ras = stro.ras(:,strcmp(stro.sum.rasterCells, 'sig001a'));
        stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
        stimtrials = find(~isnan(stimon_t))';
        tmp = zeros(length(stimtrials),2);
        for i = stimtrials
            spiketimes = ras{i}-stimon_t(i);
            tmp(i == stimtrials,:) = [sum(spiketimes > -.01 & spiketimes < 0) sum(spiketimes > 0 & spiketimes < .01)];
        end
        data = [data; rfx rfy mean(tmp,1) filecounter];
    end
end

L = data(:,1) > 0 & data(:,2) < 0; % Only targets in the appropriate quadrant
data(~L,:) = [];
L = data(:,3) > 0 & data(:,4) > 0;
data(~L,:) = [];

if SCALEPOINTS
    spikecount_ratio = data(:,4)./data(:,3);
    if (isnan(spikecount_ratio) | spikecount_ratio == 0)
        spikecount_ratio = eps;
    elseif (isinf(spikecount_ratio))
        spikecount_ratio = 1;
    end
else
    spikecount_ratio = ones(size(data,1),1);
end


% Plotting a symbol for every experiment
figure; axes; hold on;
unique_rfxys = unique(data(:,[1 2]),'rows');
for i = 1:size(data,1)
    if (SCALEPOINTS)
        h = plot(data(i,1)/10,data(i,2)/10,'k.','MarkerSize',10*spikecount_ratio(i),'MarkerEdgeColor',[spikecount_ratio(i)/max(spikecount_ratio) 0 0]);
    else
        h = plot(data(i,1)/10,data(i,2)/10,'k.','MarkerSize',10*spikecount_ratio(i),'MarkerEdgeColor','black');
    end
end
set(gca,'Xlim',[0 8],'Ylim',[-8 0])
axis equal;
disp(['Max Spike count ratio is ',num2str(max(spikecount_ratio))]);

% Plotting one symbol for each RF location
figure; axes; hold on;
unique_rfxys = unique(data(:,[1 2]),'rows');
for i = 1:size(unique_rfxys,1)
    L = data(:,1) == unique_rfxys(i,1) & data(:,2) == unique_rfxys(i,2);
    scr = max(spikecount_ratio(L));
    h = plot(unique_rfxys(i,1)/10,unique_rfxys(i,2)/10,'k.','MarkerSize',10*scr,'MarkerEdgeColor',[scr/max(spikecount_ratio) 0 0]);    
end
set(gca,'Xlim',[0 8],'Ylim',[-8 0])
axis equal;

% A key
figure; axes; hold on;
set(gca,'Xlim',[0 8],'Ylim',[-8 0])
axis equal;
scr_steps = [1 2 5 max(spikecount_ratio)];
for i = 1:length(scr_steps)
    i
    h = plot(i,-4,'k.','MarkerSize',10*scr_steps(i),'MarkerEdgeColor',[scr_steps(i)/max(spikecount_ratio) 0 0]);
    text(i,-3,num2str(scr_steps(i)));
end

%%
% Section 3
% Fixstim fiber at SC - behavioral effect figure
stro = nex2stro(findfile('S013013008'));
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
fpoff_stim_delay = nanmedian([stimon_t-fpoff_t])
L = isnan(stimon_t);
stimon_t(L) = fpoff_t(L)+fpoff_stim_delay;
H = stro.ras(:,strcmp(stro.sum.rasterCells,'AD11'));
V = stro.ras(:,strcmp(stro.sum.rasterCells,'AD12'));
ep_t = stro.ras(:,strcmp(stro.sum.rasterCells,'anlgStartTime'));
ADfreq = stro.sum.analog.storeRates{1};
offset = [0 .15];
ALIGNEPS = 1;
% Xt and Yt plots
subplotcntr = 0;
figure;  % ({x,y},t)
for i = 0:1  % Target appearing or not
    for j = [0 1]  % H or V
        subplotcntr=subplotcntr+1;
        subplot(2,3,subplotcntr); hold on;
        for k = [0 1] % No stim, optical stim (TTL)
            L = stim_type == k & targ_shown == i;
            trlidxs = find(L);
            for counter = 1:sum(L)
                trlidx = trlidxs(counter);
                h = H{trlidx}*4096/400;
                v = V{trlidx}*4096/400;
                t = [0:length(h)-1]./ADfreq+ep_t{trlidx};
                t = t-stimon_t(trlidx);
                if (ALIGNEPS)
                    idx = find(t.^2 == min(t.^2));
                    h = h-h(idx);
                    v = v-v(idx);
                end
                if (j == 0 & k == 0) % H, no stim
                    plot(t,h,'k-');
                elseif (j == 0 & k == 1) % H, stim
                    plot(t,h,'-','Color',[0 .75 0]);
                elseif (j == 1 & k == 0) % V, no stim
                    plot(t,v,'k-');
                elseif (j == 1 & k == 1) % V, stim
                    plot(t,v,'Color',[0 .75 0]);
                end
            end
        end
        set(gca,'Xlim',offset);
    end
end
subplot(2,3,1)
set(gca,'Ylim',[-5 25]);
% Order of panels:
% 1) No target, horizontal 
% 2) No target, vertical
% 3) Target, horizontal 
% 4) Target, vertical

