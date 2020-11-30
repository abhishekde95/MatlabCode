% Analyzing OFF responses from all the cells during the isoreponse measurement - Are DO cells cone-opponent simple cells?
% Author - Abhishek De, 10/18
close all; clearvars;

load newDOidx.mat
load newLUMidx.mat
load newSOidx.mat
load newhardtoclassifyidx.mat
DOidx = newDOidx';
LUMidx = newLUMidx';
hardtoclassifyidx = [newhardtoclassifyidx' newSOidx'];

% Loading the files 
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filename = fetch(conn,'SELECT filename FROM WNthresh');
NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
close(conn);
filename = filename(strcmp(string(NTmode),"subunit"));
NTmode = NTmode(strcmp(string(NTmode),"subunit"));
spikeidx_NT = spikeidx_NT(strcmp(string(NTmode),"subunit"));

binwidth = 0.010;
N = numel(filename);
t1 = 0.15; t2 = 0.15;
PSTHbins = -t1:binwidth:t2;
PSTH = cell(1,N);
OFFrespindices= find(PSTHbins>=0 & PSTHbins<=0.15); % from STIMOFFCD to STIMOFFCD + 150 ms
afterOFFrespindices= find(PSTHbins>0.15);
eyepossamplingrate = 500; % in Hz
sampligrate = 1000; % in Hz for spikes 
baselineFR = [];
OFF_FR = [];
pOFF = [];
rOFF = [];
plot_counter = 1;
numsubplots = ceil(sqrt(N));
plotrasters = 0;
for ii = 1:N
    ind = ii;
    % Plotting the OFF responses
    fileofinterest = char(filename(ind,:));
    stro = nex2stro(findfile(fileofinterest));
    spikename = 'sig001a';%getSpikenum(stro);
    neurothreshidx = strcmp(stro.sum.trialFields(1,:),'neurothresh');
    maskidx = strcmp(stro.sum.rasterCells(1,:),'subunit_mask');
    basisvecidx = strcmp(stro.sum.rasterCells(1,:),'basis_vec');
    stimoffidx = strcmp(stro.sum.trialFields(1,:),'stim_off');
    latencyidx = strcmp(stro.sum.trialFields(1,:),'latency');
    stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on');
    fpacqidx = strcmp(stro.sum.trialFields(1,:),'fpacq');
    spikestampidx = find(strcmp(stro.sum.rasterCells(1,:),spikename));
    mask_changes = [2];
    all_masks = stro.ras(:,maskidx);
    Fx = @(xi) any(isnan(xi)); % function that finds 'NaN' in a cell array
    inds = find(cellfun(Fx,stro.ras(:,basisvecidx))==0);
    if isempty(inds)
        inds = size(stro.trial,1)-1;
    end
    t_offset = stro.trial(end,latencyidx)/1000;
    neurothreshmode = stro.trial(:,neurothreshidx);
    basisvec_dropidx = inds(end);
    neurothresh_startidx = max(find(neurothreshmode==1,1),basisvec_dropidx+1);
    baselinespikecounts = []; % from stimon - 150ms to stimon
    spikecountsduringstimpresent = []; % from stimon to stimoff
    spikecountsafter = []; % from stimoff to stimoff + 150ms
    idxs = neurothresh_startidx:1:size(stro.trial,1);
    t_offset = stro.trial(end,latencyidx)/1000;
    
    num_dur = [];
    for jj = 1:numel(idxs)
        tmp = cell2mat(stro.ras(idxs(jj),spikestampidx));
        baselinespikecounts = [baselinespikecounts; numel(tmp(tmp<stro.trial(idxs(jj),stimonidx) & tmp>stro.trial(idxs(jj),stimonidx)-0.15))];
        spikecountsduringstimpresent = [spikecountsduringstimpresent; numel(tmp(tmp<stro.trial(idxs(jj),stimoffidx) & tmp>stro.trial(idxs(jj),stimonidx)+t_offset))];
        spikecountsafter = [spikecountsafter; numel(tmp(tmp<stro.trial(idxs(jj),stimoffidx)+0.15 & tmp>stro.trial(idxs(jj),stimoffidx)+t_offset))];
        num_dur = [num_dur; stro.trial(idxs(jj),stimonidx)-stro.trial(idxs(jj),fpacqidx)];    
    end
    [newspikecounts,newidxs] = sort(spikecountsduringstimpresent);
    
    [r,p] = corr(spikecountsduringstimpresent,spikecountsafter,'type','Spearman');
    baselineFR = [baselineFR; mean(baselinespikecounts)/0.15];
    OFF_FR = [OFF_FR; mean(spikecountsafter)/(0.15-t_offset)];
    rOFF = [rOFF; r];
    pOFF = [pOFF; p];
    
    count = 1; spcounts = [];
    % for plotting purposes
    if plotrasters
        analogstridx = find(strcmp(stro.sum.rasterCells(1,:),'anlgStartTime'));
        idxs = idxs(newidxs);
        spcount = [];
        for jj = 1:numel(idxs)
            tmp = cell2mat(stro.ras(idxs(jj),spikestampidx));
            spikes = tmp(tmp<stro.trial(idxs(jj),stimoffidx)+t2 & tmp>stro.trial(idxs(jj),stimonidx)-t1);
            spikes = spikes - stro.trial(idxs(jj),stimonidx);
            figure(plot_counter);subplot(numsubplots,numsubplots,ii); plot(spikes,(count)*ones(1,length(spikes)),'k.'); hold on;
            count = count + 1;
            spcount = [spcount; numel(tmp(tmp<stro.trial(idxs(jj),stimoffidx) & tmp>stro.trial(idxs(jj),stimonidx)))];
        end
        if any(DOidx == ind)
            c = 'r';
        elseif any(LUMidx == ind)
            c = 'g';
        else
            c = 'k';
        end
        set(gca,'XColor',c,'YColor',c,'Xlim',[-0.2 0.45],'Ylim',[0 count]);
        line([0 0],[0 numel(idxs)],'Color',[1 0 0]);
        line([0.3 0.3],[0 numel(idxs)],'Color',[1 0 0]);
        line([-0.15 -0.15],[0 numel(idxs)],'Color',[1 0 0]);
        line([0.45 0.45],[0 numel(idxs)],'Color',[1 0 0]); drawnow; hold off;
    end
    
end
plot_counter = plot_counter + 1;
savevariables = 0;
if savevariables == 1
    save rOFF rOFF
    save pOFF pOFF
    save baselineFR baselineFR
    save OFF_FR OFF_FR
end

%% Doing some more population analyses
group = [ones(numel(DOidx),1); 2*ones(numel(LUMidx),1); 3*ones(numel(hardtoclassifyidx),1)];
idx = [DOidx LUMidx hardtoclassifyidx];
p1 = kruskalwallis(baselineFR(idx),group,'off'); % Comparing the baselineFR 
p2 = kruskalwallis(OFF_FR(idx),group,'off'); % Comparing the OFF FR 
p3 = kruskalwallis(OFF_FR(idx)-baselineFR(idx),group,'off'); % Comparing the difference between the OFF FR and the baseline FR
p4 = signtest(OFF_FR(DOidx)-baselineFR(DOidx)); % Testing the significance of OFF responses for DO cells 
p5 = signtest(OFF_FR(LUMidx)-baselineFR(LUMidx)); % Testing the significance of OFF responses for LUM cells
p6 = signtest(OFF_FR(hardtoclassifyidx)-baselineFR(hardtoclassifyidx)); % Testing the significance of OFF responses for LUM cells

M = max([baselineFR;OFF_FR]);
figure(plot_counter);
subplot(221); plot(baselineFR,OFF_FR,'o','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on; plot([0 M],[0 M],'k'); 
set(gca,'Tickdir','out','Xlim',[0 M],'Ylim',[0 M]); axis square; xlabel('baseline FR'); ylabel('OFF response'); 
subplot(222); plot(baselineFR(DOidx),OFF_FR(DOidx),'o','MarkerSize',8,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on;  
plot(baselineFR(LUMidx),OFF_FR(LUMidx),'o','MarkerSize',8,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
plot(baselineFR(hardtoclassifyidx),OFF_FR(hardtoclassifyidx),'o','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); legend('DO','LUM','htc');
plot([0 M],[0 M],'k'); set(gca,'Tickdir','out','Xlim',[0 M],'Ylim',[0 M]); axis square; xlabel('baseline FR'); ylabel('OFF response'); 
subplot(223); boxplot(baselineFR(idx),group); set(gca,'Tickdir','out','XTicklabel',{'DO','Lum','oth'}); ylabel('baseline FR'); 
title(strcat('p=',num2str(p1,2))); axis square;
subplot(224); boxplot(OFF_FR(idx)-baselineFR(idx),group); set(gca,'Tickdir','out','XTicklabel',{'DO','Lum','oth'}); ylabel('OFF FR - baseline FR'); 
title(strcat('p=',num2str(p3,2))); axis square;
plot_counter = plot_counter + 1;

figure(plot_counter); 
boxplot(rOFF(idx),group); set(gca,'Tickdir','out','XTicklabel',{'DO','Lum','oth'}); ylabel('r OFF'); axis square;
plot_counter = plot_counter + 1;
