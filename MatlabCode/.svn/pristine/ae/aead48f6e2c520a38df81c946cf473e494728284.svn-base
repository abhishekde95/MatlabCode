% A similar analysis to Popwaveform_ISI_burstindex.m but now looking at the
% auto-correlation function and firing rate during the stimulation period. 
% Author - Abhishek De, 12/19
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
numcells = numel(Input_List);
files_not_working = [];
count = 1;
numsubplots = 2;
maximagesperplot = numsubplots^2;
Firing_rate_stim = [];
counter = 1;
BRI = [];
Propbursts_returnmap = []; % based on return map
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
    prespikeISI = []; postspikeISI = [];
    for jj = mask_changes(1):mask_changes(end)
        stimontime = WN.trial(jj,stimonidx);
        stimofftime = WN.trial(jj,stimoffidx);
        spikes = WN.ras{jj,spikeidx}-stimontime-0.2;
        spikes = spikes(spikes<=timedur & spikes>=0);
        tmp =  histcounts(spikes,[0:binwidth:timedur]);
        N = [N; tmp];
        spiketimes = [spiketimes; spikes];
        
        % Pre- and post-spike time calculation
        spikes2 = WN.ras{jj,spikeidx};
        spikes2 = spikes2(spikes2>=stimontime & spikes2<=stimofftime);
        prespikeISI = [prespikeISI; diff(spikes2(1:end-1))];
        postspikeISI = [postspikeISI; diff(spikes2(2:end))];
    end
    Propbursts_returnmap = [Propbursts_returnmap; sum((prespikeISI>0.1 & postspikeISI<0.004))/numel(prespikeISI>0.1)];
    PSTH = mean(N,1);
    Firing_rate_stim = [Firing_rate_stim; sum(N(:))/(size(N,1)*timedur)];
     
%     figure(plot_counter); subplot(numsubplots,numsubplots,count); plot(prespikeISI,postspikeISI,'k.');
%     set(gca,'Tickdir','out','Xlim',[0.0001 1],'XScale','log','XTick',[0.0001 0.001 0.01 0.1 1.0],'Ylim',[0.0001 1],'YScale','log','YTick',[0.0001 0.001 0.01 0.1 1.0]); axis square

    maxtimelag = 0.1/0.0005;
    tmp_cross_corr = zeros(1,2*maxtimelag+1);
    tmpmat = [];
    
    % Based on analyses from "Anderson EB, Mitchell JF, Reynolds JH.
    % Attention-dependent reductions in burstiness and action-potential height in macaque area V4. Nature Neuroscience."
%     for idx1 = 1:size(N,1)-1
%         for idx2 = idx1+1:size(N,1)
%             [temp,lag] = xcorr(N(idx1,:),N(idx2,:),maxtimelag,'biased');
%             tmp_cross_corr = tmp_cross_corr + temp ;
%             tmpmat = [tmpmat; temp];
%         end
%     end
%     cross_corr_shiftpredictor = tmp_cross_corr/(size(N,1)*(size(N,1)-1)/2);
%     [auto_corr,time_lag] = xcorr(PSTH,PSTH,maxtimelag,'biased');
%     corrected_autocorr = auto_corr-cross_corr_shiftpredictor;
%     corrected_autocorr = 0.5*(corrected_autocorr + fliplr(corrected_autocorr));
%     corrected_autocorr_z = corrected_autocorr./std([tmpmat; fliplr(tmpmat)],0,1); 
%     
%     idx = (0.0005*time_lag>=0.001 & 0.0005*time_lag<=0.004);
%     BRI = [BRI; mean(corrected_autocorr_z(idx))];
    
%     Corrected autocorrelation function
%     figure(plot_counter); subplot(numsubplots,numsubplots,count); plot(0.0005*time_lag(time_lag>0),corrected_autocorr_z(time_lag>0),'k');
%     set(gca,'Tickdir','out','Xlim',[0.001 0.1],'XScale','log','XTick',[0.001 0.01 0.1]); axis square
%     
%     k = cell2mat(WN.ras(:,spikeidx));
%     figure(plot_counter+4); subplot(numsubplots,numsubplots,count); plot(prespikeISI*1000,postspikeISI*1000,'k.'); 
%     set(gca,'Tickdir','out','XScale','log','YScale','log','Xlim',[1 10000],'Ylim',[1 10000]); axis square;
        
    % Actual autocorrelation function
%     figure(plot_counter+4); subplot(numsubplots,numsubplots,count); plot(0.0005*time_lag(time_lag>0),auto_corr(time_lag>0),'k');
%     set(gca,'Tickdir','out','Xlim',[0.001 0.1],'XScale','log','XTick',[0.001 0.01 0.1]); axis square
    
    count = count + 1;
    counter = counter + 1;
    if count == maximagesperplot+1
        plot_counter = plot_counter + 1;
        count = 1;
    end
end
plot_counter = plot_counter + 1;

savevariables = 1; 
if savevariables 
    save Propbursts_returnmap Propbursts_returnmap
end

