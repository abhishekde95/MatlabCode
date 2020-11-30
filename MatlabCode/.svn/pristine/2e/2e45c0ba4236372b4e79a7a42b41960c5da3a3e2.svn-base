% This script will allow me to see how the RGB weights converge to the STA RGBs during the second phase of the experiment
% Author - Abhishek De, 11/18
close all; clearvars;
plot_counter = 1;

load totspikes.mat
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

numsubunitspikes = totspikes(:,1); % Spikes collected during the WN subunit stimuli
plot_counter = 1;
N = numel(filename);
numsubplots = ceil(sqrt(N));
angledifflastpt = [];
anglediffall = nan(N,max(numsubunitspikes));
for ii = 1:N
    anglediff = [];
    global reversalflagidx stepsizescale stepsize nreversals
    stro = nex2stro(findfile(char(filename(ii,:))));
    spikename = 'sig001a';%getSpikenum(stro);
    maskidx = strcmp(stro.sum.rasterCells(1,:),'subunit_mask');
    spikeidx = strcmp(stro.sum.rasterCells(1,:),spikename);
    nstixperside = stro.sum.exptParams.nstixperside;
    ngammasteps = 2^16; % 65536
    linepredtol = stro.sum.exptParams.linepredtol;
    stepsizescale = stro.sum.exptParams.stepsizescale;
    stepsize = stro.sum.exptParams.stepsize;
    nreversals = stro.sum.exptParams.nreversals;
    reversalflagidx = strcmp(stro.sum.trialFields(1,:),'reversalflag');
    bkgnd_ridx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_r'));
    bkgnd_gidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_g'));
    bkgnd_bidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_b'));
    seedidx = strcmp(stro.sum.trialFields(1,:),'seed');
    correctidx = strcmp(stro.sum.trialFields(1,:),'correct');
    nframesidx = strcmp(stro.sum.trialFields(1,:),'num_frames');
    stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on');
    stimoffidx = strcmp(stro.sum.trialFields(1,:),'stim_off');
    fponidx = strcmp(stro.sum.trialFields(1,:),'fp_on');
    fpacqidx = strcmp(stro.sum.trialFields(1,:),'fpacq');
    weightsidx = strcmp(stro.sum.rasterCells(1,:),'weights');
    muidxs = [find(strcmp(stro.sum.trialFields(1,:),'mu1')), ...
        find(strcmp(stro.sum.trialFields(1,:),'mu2')), ...
        find(strcmp(stro.sum.trialFields(1,:),'mu3'))];
    sigmaidxs = [find(strcmp(stro.sum.trialFields(1,:),'sigma1')), ...
        find(strcmp(stro.sum.trialFields(1,:),'sigma2')), ...
        find(strcmp(stro.sum.trialFields(1,:),'sigma3'))];
    latencyidx = strcmp(stro.sum.trialFields(1,:),'latency');
    basisvecidx = strcmp(stro.sum.rasterCells(1,:),'basis_vec');
    msperframe = 1000/stro.sum.exptParams.framerate;
    basisvecdiridx = strcmp(stro.sum.trialFields(1,:),'weights_idx');
    neurothreshidx = strcmp(stro.sum.trialFields(1,:),'neurothresh'); % when exactly the neurothresh trials started
    targetspikerateidx = strcmp(stro.sum.trialFields(1,:),'targetspikerate');
    ntrials = size(stro.trial,1);
    maxT = 15; % this represents the temporal part in the spatiotemporal receptive field
    xx = linspace(stro.sum.exptParams.gauss_locut/1000, stro.sum.exptParams.gauss_hicut/1000,ngammasteps); % xx represents the probabilities. For more info, have a look at the MATLAB 'norminv' function.
    yy = norminv(xx'); % defining norminv to extract the values for which the cdf values range between gauss_locut and gauss_hicut
        
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
    use_STCOVmex_ST = 1; % 1 - Will use STCOVmex_ST, 0 - Will use STCOVmex
    subunitmasktrials = mask_changes(:,2);
    
    % As of now just calculate the WNsubunit STA
    st_mask = stro.ras{subunitmasktrials(1),maskidx}; % subunit mask
    st_mask(st_mask == 0) = Inf;
    [stIdxs,~,~] = unique(st_mask); % now the Infs map to nsubunits+1
    num_subunits = length(stIdxs)-any(isinf(stIdxs)); % nsubunits, like subunits A and B
    if use_STCOVmex_ST
        STCOV_st('init', {num_subunits 3 maxT});
    else
        STCOVmex('init', {num_subunits 3 maxT});
    end
    cum_rgbs = []; % array for storing all the rgbs frames
    cum_n = [];
    
    for k = subunitmasktrials(1):subunitmasktrials(2)
        nframes = stro.trial(k,nframesidx);
        if (nframes == 0)
            continue;
        end
        seed = stro.trial(k,seedidx);
        mu = stro.trial(k,muidxs)/1000;
        sigma = stro.trial(k,sigmaidxs)/1000;
        
        % org_mask tells u if u have updated the mask or not. If org_mask is non-zero it means at this particular trial
        % u have selected the subunits and need to analyse its computation
        org_mask = stro.ras{k,maskidx};
        if any(org_mask)
            org_mask(org_mask == 0) = Inf;
            [subunitIdxs,~,mask] = unique(org_mask); % now the Infs map to nsubunits+1
            nrandnums_perchannel = length(subunitIdxs)-any(isinf(subunitIdxs)); % nsubunits, like subunits A and B
            mask = [mask; mask+max(mask); mask+2*max(mask)]; %#ok<AGROW>
        else
            nrandnums_perchannel = nstixperside^2; % In this case, it is the 100 pixels which are flickering when no subunits are selected
        end
        
        % assuming Gaussian gun noise only, random number generator
        % routine as a mexfile (getEJrandnums.mexw64)
        invnormcdf = bsxfun(@plus, bsxfun(@times, yy, sigma), mu);
        randnums = getEJrandnums(3*nrandnums_perchannel*nframes, seed); % random numbers for 300 x 9 pixels
        % This is the extracted colors for subunits/pixels using the seed number
        randnums = reshape(randnums, 3*nrandnums_perchannel, nframes);
        for gun = 1:3
            idxs = (1:nrandnums_perchannel)+nrandnums_perchannel*(gun-1);
            randnums(idxs,:) = reshape(invnormcdf(randnums(idxs,:)+1,gun),[length(idxs) nframes]);
        end
        
        rgbs = randnums;
        t_stimon = stro.trial(k, stimonidx);
        spiketimes = (stro.ras{k,spikeidx}-t_stimon)*1000; % observing spiketimes in milliseconds
        frametimes = linspace(0, nframes*msperframe, nframes)+(msperframe/2)';
        spiketimes(spiketimes < maxT*msperframe) = [];
        spiketimes(spiketimes > frametimes(end)) = [];
        n = hist(spiketimes, frametimes);
        if use_STCOVmex_ST
            STCOV_st(rgbs(:),n);
        else
            STCOVmex(rgbs(:),n);
        end
        cum_rgbs = [cum_rgbs rgbs];
        cum_n = [cum_n n];
        
    end
    if use_STCOVmex_ST
        out = STCOV_st('return');
    else
        out = STCOVmex('return');
    end
    STS = out{1};  % A (dimension) x 9(frames) matrix
    clear STCOV_st out;
    STAs = STS/numsubunitspikes(ii);
    [~,whichframe] = max(fliplr(sum(STAs.^2,1)));
    stindices = find(cum_n>0);
    stindices = stindices - whichframe+1;
    
    STA = STAs(:,maxT-whichframe+1);
    for jj = 1:numel(stindices)
        newSTA = mean(cum_rgbs(:,stindices(1:jj)),2);
        anglediff = [anglediff; 180*acos(dot(STA,newSTA)/(norm(STA)*norm(newSTA)))/pi];
    end
    anglediffall(ii,1:numel(anglediff)) = anglediff';
    if any(DOidx == ii)
        c = 'b';
    elseif any(LUMidx == ii)
        c = 'g';
    else
        c = 'k';
    end
%    figure(plot_counter); subplot(numsubplots,numsubplots,ii), plot(anglediff,'Linewidth',2); 
%    set(gca,'XColor',c,'YColor',c); hold off;
end
plot_counter = plot_counter + 2;


% Here I am looking at a few more analyses 
err = nanstd(anglediffall,1)./sqrt(nansum(anglediffall,1));
err_DO = nanstd(anglediffall([DOidx],:),1)./sqrt(nansum(anglediffall([DOidx],:),1));
err_LUM = nanstd(anglediffall([LUMidx],:),1)./sqrt(nansum(anglediffall([LUMidx],:),1));
err_htc = nanstd(anglediffall([hardtoclassifyidx],:),1)./sqrt(nansum(anglediffall([hardtoclassifyidx],:),1));
figure(plot_counter); 
subplot(221); errorbar(nanmean(anglediffall,1),err,'Linewidth', 2); ylabel('angle diff'); xlabel('no. spikes'); axis square; set(gca,'Tickdir','out','Ylim',[0 90],'Xlim',[0 max(numsubunitspikes)]); 
subplot(222); plot(nanmean(anglediffall([DOidx],:),1),'r','Linewidth', 2); hold on; 
plot(nanmean(anglediffall([LUMidx],:),1),'g','Linewidth', 2); plot(nanmean(anglediffall([hardtoclassifyidx],:),1),'k','Linewidth', 2);
ylabel('angle diff'); xlabel('Number of spikes'); axis square; set(gca,'Tickdir','out','Ylim',[0 90],'Xlim',[0 max(numsubunitspikes)]); legend('DO','LUM','htc'); hold off;
subplot(223); errorbar(nanmean(anglediffall([DOidx],:),1),err_DO,'r','Linewidth',1); hold on; 
errorbar(nanmean(anglediffall([LUMidx],:),1),err_LUM,'g','Linewidth',1); errorbar(nanmean(anglediffall([hardtoclassifyidx],:),1),err_htc,'k','Linewidth',1);
ylabel('angle diff'); xlabel('Number of spikes'); axis square; set(gca,'Tickdir','out','Ylim',[0 90],'Xlim',[0 max(numsubunitspikes)]); legend('DO','LUM','htc'); hold off;
plot_counter = plot_counter + 1;