% Compute autocorrelogram: based on " A Distinct class of bursting neurons
% with strong gamma synchronization and stimulus selectivity in monkey V1. 
% Author - Abhishek De, 11/19

close all; clearvars;
load fundamentals.mat 
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;
resize_fact = 1;
% Include Gun noise data for statistical tests in order to estimate the RF size
channels = 3;
NPOINTS = 65536;
resize_fact2 = 1;
CHI2CRIT = 0.95; % For flaging a stixel as significant (adding gun and cone noise z-scores, squared)
maxT = 9;
crit = chi2inv(CHI2CRIT,300); % 3 color channels
spikename_options = ['sig001a'; 'sig001b'];
[filename_Lum, spikeIdx_Lum] = fnamesFromTxt2('Lum.txt');
[filename_ColorOpponent, spikeIdx_ColorOpponent] = fnamesFromTxt2('ColorOpponent.txt');

% Loading the Neurothresh files
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
tmp = fetch(conn,'SELECT filename FROM WNthresh'); tmp_filename = tmp;
tmp = fetch(conn,'SELECT NTmode FROM WNthresh'); NTmode = tmp;
tmp = fetch(conn,'SELECT spikeidx FROM WNthresh'); spikeidx_NT = tmp;
close(conn);
tmp_filename = tmp_filename(strcmp('subunit',NTmode));
spikeidx_NT = spikeidx_NT(strcmp('subunit',NTmode));
filename_NT = cell(size(tmp_filename));
for kk = 1:size(tmp_filename)
    filename_NT(kk) = {tmp_filename(kk)};
end

% Loading the just the subunit files 
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
tmp = fetch(conn,'SELECT filename FROM WNSubunit'); tmp_filename = tmp;
tmp = fetch(conn,'SELECT mode FROM WNSubunit'); subunit_mode = tmp;
tmp = fetch(conn,'SELECT spikeidx FROM WNSubunit'); spikeidx_subunit = tmp;
close(conn);
tmp_filename = tmp_filename(strcmp('STA',subunit_mode));
spikeidx_subunit = spikeidx_subunit(strcmp('STA',subunit_mode));
filename_subunit = cell(size(tmp_filename));
for kk = 1:size(tmp_filename)
    filename_subunit(kk) = {tmp_filename(kk)};
end

% Merging all the files in the list
Input_List = [filename_Lum; filename_ColorOpponent; filename_NT; filename_subunit];
spikeIdx = cell2mat([spikeIdx_Lum; spikeIdx_ColorOpponent; spikeidx_NT; spikeidx_subunit]);

plot_counter = 1;
numcells = 9;%numel(Input_List);
files_not_working = [];
count = 1;
Singleopponent_waveform = [];
numsubplots = 3;
maximagesperplot = numsubplots^2;
ISI = [];
counter = 1;
ISI_peaktime = [];
burst_index = [];
for ii = 1:numcells
    flag = 0;
    disp(ii);
    filename = char(Input_List{ii}{1}); % acquiring the filename (1st column) from the List
    
    WN = {};
    for jj = 1:size(Input_List{ii},2)
        try 
            tmpstro = nex2stro(findfile(char(Input_List{ii}(jj))));
        catch 
            files_not_working = [files_not_working; Input_List{ii}];
            flag = 1;
            break;
        end
        if (isempty(WN))
            WN = tmpstro;
        else
            WN = strocat(WN, tmpstro);
        end
        if ~any(strcmp(WN.sum.rasterCells,'sig001a_wf'))
            files_not_working = [files_not_working; Input_List{ii}];
            flag = 1;
        end
    end
    if flag 
        continue;
    end

    framerate = WN.sum.exptParams.framerate;
    nstixperside = WN.sum.exptParams.nstixperside;
    ntrials = length(WN.sum.absTrialNum);
    fponidx = find(strcmp(WN.sum.trialFields(1,:),'fp_on'));
    stimonidx = find(strcmp(WN.sum.trialFields(1,:),'stim_on'));
    stimoffidx = find(strcmp(WN.sum.trialFields(1,:),'all_off'));
    nframesidx = find(strcmp(WN.sum.trialFields(1,:),'num_frames'));
    noisetypeidx = find(strcmp(WN.sum.trialFields(1,:),'noise_type'));
    sigmaidxs = strmatch('sigma',WN.sum.trialFields(1,:));
    hepidx = find(strcmp(WN.sum.rasterCells(1,:),'AD11'));
    vepidx = find(strcmp(WN.sum.rasterCells(1,:),'AD12'));
    maskidx = strcmp(WN.sum.rasterCells(1,:), 'subunit_mask');
    anlgStartTimeidx = find(strcmp(WN.sum.rasterCells(1,:),'anlgStartTime'));
    L = WN.trial(:,noisetypeidx)==1;
    mu1idx = find(strcmp(WN.sum.trialFields(1,:),'mu1'));
    mu2idx = find(strcmp(WN.sum.trialFields(1,:),'mu2'));
    mu3idx = find(strcmp(WN.sum.trialFields(1,:),'mu3'));
    sigma1idx = find(strcmp(WN.sum.trialFields(1,:),'sigma1'));
    sigma2idx = find(strcmp(WN.sum.trialFields(1,:),'sigma2'));
    sigma3idx = find(strcmp(WN.sum.trialFields(1,:),'sigma3'));
    maskidx = strcmp(WN.sum.rasterCells(1,:),'subunit_mask');
    basisvecidx = strcmp(WN.sum.rasterCells(1,:),'basis_vec');
    spikeidx = spikeIdx(ii);
    spikename = spikename_options(spikeidx,:); 
    spikeidx = strcmp(WN.sum.rasterCells(1,:),spikename);
    gammaTable = WN.sum.exptParams.gamma_table;
    gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
    gammaTable1 = interp1(linspace(0,255,256),gammaTable,linspace(0,255,65536), 'spline');
    invgamma = InvertGamma(gammaTable, 0);
    sigmavect = unique(WN.trial(L,[sigma1idx sigma2idx sigma3idx]),'rows')/1000;
    sigmavect(all(any(sigmavect == 0),2),:) = [];
    gausslims = [WN.sum.exptParams.gauss_locut WN.sum.exptParams.gauss_hicut]/1000;
    x = linspace(gausslims(1),gausslims(2),NPOINTS);
    Fx = norminv(x)*sigmavect(1);
    sigmacorrectionfactor = std(Fx)./sigmavect(1);
    muvar = (sigmavect(1)*sigmacorrectionfactor)^2;
    
    
    WN.ras(~L ,:) = []; % modiftying the WN structure
    WN.trial(~L,:) = []; % modiftying the WN structure
    mask_changes = [2 size(WN.trial,1)];
    if any(basisvecidx)
        mask_changes = [2];
        all_masks = WN.ras(:,maskidx);
        Fx = @(xi) any(isnan(xi)); % function that finds 'NaN' in a cell array
        inds = find(cellfun(Fx,WN.ras(:,basisvecidx))==0);
        if isempty(inds)
            inds = size(WN.trial,1)-1;
        end
        last_wntrial =  inds(1)-1;
        for k = 3:last_wntrial
            if isequal(all_masks{k}, all_masks{k-1}) %|| all(all_masks{k} == 0) && any(isnan(all_masks{k-1}))
                continue
            else
                mask_changes = [mask_changes k-1 k]; %#ok<AGROW>
            end
        end
        if mask_changes(end) == last_wntrial
            mask_changes(end) = [];
        else
            mask_changes = [mask_changes  last_wntrial];
        end
        mask_changes = reshape(mask_changes , 2, []);
        mask_changes  = mask_changes(:,1);  
        
        idxs = zeros(size(WN.trial,1),1);
        idxs(mask_changes(2,1)+1:end) = 1;
        idxs = logical(idxs);
        WN.ras(idxs,:) = []; % modiftying the WN structure 
        WN.trial(idxs,:) = []; % modiftying the WN structure
    end
    
    % Next, I want to compute autocorrelogram
    N = [];
    timedur = 2.0;
    binwidth = 0.0005; % 1 milliseconds, sampling at 2000 Hz
    bins = [0:binwidth:timedur-binwidth];
    spiketimes = [];
    for jj = mask_changes(1):mask_changes(end)
        stimontime = WN.trial(jj,stimonidx);
        spikes = WN.ras{jj,spikeidx}-stimontime-0.2;
        spikes = spikes(spikes<=timedur & spikes>=0);
        tmp =  histcounts(spikes,[0:binwidth:timedur]);
        N = [N; tmp];
        spiketimes = [spiketimes; spikes];
        
    end
    PSTH = mean(N,1);
    
    % calculating burst index:based on jude's code; https://snl.salk.edu/~jude/
    halfidx = size(N,2)/2;
    info = fit_isi([N(:,1:halfidx); N(:,halfidx+1:end)], 2, 0);
    burst_index = [burst_index; max(info.all.pfit)./max(info.zall.pfit)];
    
%     [tmp,edges] = histcounts(diff(cell2mat(WN.ras(:,spikeidx))),0:binwidth:0.04); % sampling at 2000 Hz, just the first 40 ms
    edges = 0:binwidth:0.04;
    tmp = info.isibin(1:0.04/binwidth)*info.isicnt(1);
    ISI{counter} = tmp./sum(tmp); % storing the ISI 
    
    
    % Performing cross validation on ISI
    K = 5; 
    c = cvpartition(tmp,'k',K);
    SSE = [];
    xvals = edges(1:end-1);
    polyorder = 5:2:40; % polynomidal order of 5-40
    for pol = 1:numel(polyorder)
        tmp_SSE = [];
        for kk = 1:K
            p = polyfit(xvals(c.training(kk)),tmp(c.training(kk)),polyorder(pol)); % sampling at 2000 Hz
            tmp_SSE = [tmp_SSE; sum((polyval(p,xvals(c.test(kk)))-tmp(c.test(kk))).^2)];
        end
        SSE = [SSE; mean(tmp_SSE)];
    end
    [~,minpolyidx] = min(SSE);
    polyorder_opt = polyorder(minpolyidx);
    p = polyfit(xvals,tmp,polyorder_opt);
    
    % finding peak ISI time
    [~,I] = max(polyval(p,edges));
    ISI_peaktime = [ISI_peaktime; edges(I)];
    
    figure(plot_counter); subplot(numsubplots,numsubplots,count); plot(bins,PSTH,'k'); set(gca,'Tickdir','out','Xlim',[0 timedur]); axis square
    figure(plot_counter+4); subplot(numsubplots,numsubplots,count); bar(edges(1:end-1),tmp,'FaceColor','k'); 
    hold on; plot(edges(1:end-1), polyval(p,edges(1:end-1)),'r','Linewidth',2); set(gca,'Tickdir','out','Xlim',[0 0.04]); axis square; hold off;
    
    count = count + 1;
    counter = counter + 1;
    if count == maximagesperplot+1
        plot_counter = plot_counter + 1;
        count = 1;
    end
end
plot_counter = plot_counter + 1;


%%
savevariables = 0;
if savevariables 
    save ISI_peaktime ISI_peaktime
    save burst_index burst_index
end


% Further comparison
load Output_List_waveform.mat
load timediff.mat


crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
timediff(~Z_cellsofinterest) = [];
baselineFR = cell2mat(Output_List(:,18));
ISI_peaktime(~Z_cellsofinterest) = [];
burst_index(~Z_cellsofinterest) = [];
baselineFR(~Z_cellsofinterest) = [];


idx = kmeans(timediff,2);
[~,I1] = max([sum(idx==1) sum(idx==2)]); % narrow spiking 
[~,I2] = min([sum(idx==1) sum(idx==2)]); % broad spiking


%% Further splitting them into different classes 
NW_I = (idx==I1) & burst_index<=3;
NW_E = (idx==I1) & burst_index>3; % chattering cells
BW = idx==I2;

% Analyzing the baseline FR of all the cells
data = [baselineFR(NW_I); baselineFR(NW_E); baselineFR(BW)];
group = [ones(sum(NW_I),1); 2*ones(sum(NW_E),1); 3*ones(sum(BW),1)];
p = kruskalwallis(data,group,'off');


% Plotting the difference between narrow spiking excitatiory 
[p2,h2] = ranksum(baselineFR(NW_I),baselineFR(NW_E));
[p3,h3] = ranksum(baselineFR(NW_E),baselineFR(BW));

figure(plot_counter); subplot(221); boxplot(data,group); set(gca,'Tickdir','out','YScale','log'); axis square;
subplot(222); plot(timediff(NW_I),burst_index(NW_I),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on
plot(timediff(NW_E),burst_index(NW_E),'o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 0 0]);
plot(timediff(BW),burst_index(BW),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','YScale','log'); xlabel('Waveform duration'); ylabel('Burst Index');
subplot(223); plot(baselineFR(NW_I),burst_index(NW_I),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on
plot(baselineFR(NW_E),burst_index(NW_E),'o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 0 0]);
plot(baselineFR(BW),burst_index(BW),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','YScale','log'); xlabel('Baseline FR'); ylabel('Burst Index');
subplot(224); plot(ISI_peaktime(NW_I),burst_index(NW_I),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on
plot(ISI_peaktime(NW_E),burst_index(NW_E),'o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 0 0]);
plot(ISI_peaktime(BW),burst_index(BW),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'XScale','log','YScale','log','Tickdir','out'); axis square; xlabel('ISI peak time'); ylabel('burst_index');
plot_counter = plot_counter + 1;

