% Analyses of waveforms: 1 file at a time
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
tmp = fetch(conn,'SELECT filename FROM WNthresh'); tmp_filename = tmp.filename;
tmp = fetch(conn,'SELECT NTmode FROM WNthresh'); NTmode = tmp.NTmode;
tmp = fetch(conn,'SELECT spikeidx FROM WNthresh'); spikeidx_NT = tmp.spikeidx;
close(conn);
tmp_filename = tmp_filename(strcmp('subunit',NTmode));
spikeidx_NT = spikeidx_NT(strcmp('subunit',NTmode));
filename_NT = cell(size(tmp_filename));
for kk = 1:size(tmp_filename)
    filename_NT(kk) = {tmp_filename(kk)};
end

% Loading the just the subunit files
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
tmp = fetch(conn,'SELECT filename FROM WNSubunit'); tmp_filename = tmp.filename;
tmp = fetch(conn,'SELECT mode FROM WNSubunit'); subunit_mode = tmp.mode;
tmp = fetch(conn,'SELECT spikeidx FROM WNSubunit'); spikeidx_subunit = tmp.spikeidx;
close(conn);
tmp_filename = tmp_filename(strcmp('STA',subunit_mode));
spikeidx_subunit = spikeidx_subunit(strcmp('STA',subunit_mode));
filename_subunit = cell(size(tmp_filename));
for kk = 1:size(tmp_filename)
    filename_subunit(kk) = {tmp_filename(kk)};
end

Input_List = [filename_Lum; filename_ColorOpponent; filename_NT; filename_subunit];
spikeIdx = [spikeIdx_Lum; spikeIdx_ColorOpponent; spikeidx_NT; spikeidx_subunit];

numcells = numel(Input_List);
N = numcells;
plot_counter = 1;
maximagesperplot = 64;
numsubplots = ceil(sqrt(maximagesperplot));
count = 1;
% Loading a file
singlecellanalyses = 0;
for ii = 1:N
    disp(ii);
    flag = 0;
    stro = {};
    for jj = 1:size(Input_List{ii},2)
        try
            tmpstro = nex2stro(findfile(char(Input_List{ii}(jj))));
        catch
            flag = 1;
            break;
        end
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
        if ~any(strcmp(stro.sum.rasterCells,'sig001a_wf'))
            flag = 1;
        end
    end
    
    if flag
        continue;
    end
    
    noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
    basisvecidx = strcmp(stro.sum.rasterCells(1,:),'basis_vec');
    maskidx = strcmp(stro.sum.rasterCells(1,:),'subunit_mask');
    T = stro.trial(:,noisetypeidx)==1;
    stro.ras(~T ,:) = []; % modiftying the WN structure
    stro.trial(~T,:) = []; % modiftying the WN structure
    mask_changes = [2 size(stro.trial,1)];
    if any(basisvecidx)
        mask_changes = [2];
        all_masks = stro.ras(:,maskidx);
        Fx = @(xi) any(isnan(xi)); % function that finds 'NaN' in a cell array
        inds = find(cellfun(Fx,stro.ras(:,basisvecidx))==0);
        if isempty(inds)
            inds = size(stro.trial,1)-1;
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
        
        idxs = zeros(size(stro.trial,1),1);
        idxs(mask_changes(2,1)+1:end) = 1;
        idxs = logical(idxs);
        stro.ras(idxs,:) = []; % modiftying the WN structure
        stro.trial(idxs,:) = []; % modiftying the WN structure
    end
    
    samplingrate = stro.sum.waves.storeRates{1};
    stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on');
    fponidx = strcmp(stro.sum.trialFields(1,:),'fp_on');
    fpacqidx = strcmp(stro.sum.trialFields(1,:),'fpacq');
    spikeidx = strcmp(stro.sum.rasterCells(1,:),'sig001a');
    anlgStartTimeidx = strcmp(stro.sum.rasterCells(1,:),'anlgStartTime');
    
    % Extracting waveforms
    rawsignal_waveform = cell2mat(stro.ras(:,strcmp(stro.sum.rasterCells,'sig001a_wf')));
    idx = mean(rawsignal_waveform,2)<10 & mean(rawsignal_waveform,2)>-10;
    
    if sum(idx)>20000
        idx(20001:end) = 0; % As of now, this MATLAB along with this machine cannot handle large matrices. Therefore the max number of samples have been restricted to 20,000
    end
    
    signal_waveform = mean(rawsignal_waveform(idx,:),1);
    rawnoise_waveform = cell2mat(stro.ras(:,strcmp(stro.sum.rasterCells,'sig001U_wf')));
    noise_waveform = mean(rawnoise_waveform,1);
    
    [v,d] = eig(cov(rawsignal_waveform(idx,:)));
    [~,ind] = sort(diag(d));
    vec1 = v(:,ind(end)); vec2 = v(:,ind(end-1));
    
    if singlecellanalyses
        % Plotting the data in PC space
        figure(plot_counter); subplot(321); plot(rawsignal_waveform(idx,:)*vec1, rawsignal_waveform(idx,:)*vec2,'k.'); axis square;
        
        % plotting the waveform
        subplot(322); plot(signal_waveform,'k','Linewidth',2); axis square; set(gca,'Tickdir','out');
        
        % plotting all the signal waveforms
        subplot(323); plot(rawsignal_waveform(idx,:)','k'); axis square; set(gca,'Tickdir','out');
        
        % visualizing a dendrogram
        Z = linkage([rawsignal_waveform(idx,:)*vec1 rawsignal_waveform(idx,:)*vec2]);
        T = cluster(Z,'maxclust',2);
        cutoff = median([Z(end-2,3) Z(end-1,3)]);
        subplot(324); dendrogram(Z,'ColorThreshold',cutoff); axis square, set(gca,'Tickdir','out','XTick',[]);
        
        % Pruning some of the waveforms after running through the dendrogram
        [~,idx2] = max([sum(T==1) sum(T==2)]);
        [~,idx3] = min([sum(T==1) sum(T==2)]);
        newsignal_waveform1 = mean(rawsignal_waveform(T==idx2,:),1);
        newsignal_waveform2 = mean(rawsignal_waveform(T==idx3,:),1);
        
        subplot(325); plot(rawsignal_waveform(T==idx2,:)*vec1, rawsignal_waveform(T==idx2,:)*vec2,'r.'); hold on;
        plot(rawsignal_waveform(T==idx3,:)*vec1, rawsignal_waveform(T==idx3,:)*vec2,'g.'); axis square; set(gca,'Tickdir','out'); hold off;
        subplot(326); plot(newsignal_waveform1,'r','Linewidth',2); hold on; plot(newsignal_waveform2,'g','Linewidth',2); axis square; set(gca,'Tickdir','out');
        legend(strcat('n=',num2str(sum(T==idx2))),strcat('n=',num2str(sum(T==idx3)))); hold off;
    else
        
        % plotting the waveform
        figure(plot_counter); subplot(numsubplots,numsubplots,count); plot(signal_waveform,'k','Linewidth',2); hold on;
        
        % Pruning some of the waveforms after running through the dendrogram
        try
            Z = linkage([rawsignal_waveform(idx,:)*vec1 rawsignal_waveform(idx,:)*vec2]);
            T = cluster(Z,'maxclust',2);
            cutoff = median([Z(end-2,3) Z(end-1,3)]);
            [~,idx2] = max([sum(T==1) sum(T==2)]);
            [~,idx3] = min([sum(T==1) sum(T==2)]);
            newsignal_waveform1 = mean(rawsignal_waveform(T==idx2,:),1);
            newsignal_waveform2 = mean(rawsignal_waveform(T==idx3,:),1);
        catch
            keyboard;
        end
        
        subplot(numsubplots,numsubplots,count); plot(newsignal_waveform1,'r','Linewidth',2);
        plot(newsignal_waveform2,'g','Linewidth',2); axis square; set(gca,'Tickdir','out'); hold off;
        
        count = count + 1;
        
        if count == maximagesperplot+1
            plot_counter = plot_counter + 1;
            count = 1;
        end
    end
end




