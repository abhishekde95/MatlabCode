% A new population analysis where I am comparing the model residual ratio
% and the "percentile" value of the highest eigenvalue
% Author - Abhishek De, 1/18
close all; clearvars;

% Loading the files 
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filename = fetch(conn,'SELECT filename FROM WNthresh');
NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
close(conn);
filename = filename(strcmp(string(NTmode),"subunit"));
NTmode = NTmode(strcmp(string(NTmode),"subunit"));
spikeidx_NT = spikeidx_NT(strcmp(string(NTmode),"subunit"));

plot_counter = 1;
vals = zeros(numel(filename),1);
z_scores = zeros(numel(filename),1);
numsubunitspikes = zeros(numel(filename),1);
PCS1RGB = zeros(3,numel(filename));
PCS2RGB = zeros(3,numel(filename));
PCS1LMS = zeros(3,numel(filename));
PCS2LMS = zeros(3,numel(filename));
use_STCOVmex_ST = 0; % 1 - Will use STCOVmex_ST, 0 - Will use STCOVmex

for ii = 1:numel(filename)
    disp(ii);
    stro = nex2stro(findfile(char(filename(ii,:))));
    global spikename maskidx spikeidx nstixperside ngammasteps seedidx nframesidx
    global fponidx stimoffidx stimonidx muidxs sigmaidxs basisvecidx fpacqidx latencyidx basisvecdiridx neurothreshidx targetspikerateidx
    global msperframe ntrials maxT xx yy M
    spikename = 'sig001a';%getSpikenum(stro);
    maskidx = strcmp(stro.sum.rasterCells(1,:),'subunit_mask');
    spikeidx = strcmp(stro.sum.rasterCells(1,:),spikename);
    nstixperside = stro.sum.exptParams.nstixperside;
    ngammasteps = 2^16; % 65536
    seedidx = strcmp(stro.sum.trialFields(1,:),'seed');
    nframesidx = strcmp(stro.sum.trialFields(1,:),'num_frames');
    stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on');
    stimoffidx = strcmp(stro.sum.trialFields(1,:),'stim_off');
    fponidx = strcmp(stro.sum.trialFields(1,:),'fp_on');
    fpacqidx = strcmp(stro.sum.trialFields(1,:),'fpacq');
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
    
    % Obtaining the M matrix, code extracted from Greg, fitting a cubic spline
    % using the command 'spline'. 'SplineRaw' only availabe through
    % psychtoolbox which I currently don't have now.
    fundamentals = stro.sum.exptParams.fundamentals; % CONE FUNDAMENTALS: L,M,S
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
    mon_spd = stro.sum.exptParams.mon_spd; % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
    M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals
    M = inv(M');
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
    subunitmasktrials = mask_changes(:,2);
    % As of now just calculate the WNsubunit STA
    st_mask = stro.ras{subunitmasktrials(1),maskidx}; % subunit mask
    st_mask(st_mask == 0) = Inf;
    [stIdxs,~,~] = unique(st_mask); % now the Infs map to nsubunits+1
    num_subunits = length(stIdxs)-any(isinf(stIdxs)); % nsubunits, like subunits A and B
    STCOVmex('init', {num_subunits 3 maxT});
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
        % u have selected the subunits and need to analyse its computation.
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
        % Deleting the spikes taking place before the first 9 frames as I need to look at the 9 preceding frames
        spiketimes(spiketimes < maxT*msperframe) = [];
        % Deleting the spikes that take place after the stimulus was
        % removed as it would imply that these spikes do not result from
        % the stimulus shown on the screen
        spiketimes(spiketimes > frametimes(end)) = [];
        n = hist(spiketimes, frametimes);
        STCOVmex(rgbs(:),n);
        cum_rgbs = [cum_rgbs rgbs];
        cum_n = [cum_n n]; 
    end
    
    out = STCOVmex('return');    
    STS = out{1};  % A (dimension) x 9(frames) matrix
    STCross = out{2};  % A (dimension x frames)x (dimension x frames)  matrix
    nspikes = out{3}; % Number of spikes in the given file
    clear STCOVmex out;
    numsubunitspikes(ii) = nspikes;
    % Coverting the STS and the STCross into STA and STC respectively
    STAs = STS/nspikes;
    [~,whichframe] = max(sum(STAs.^2,1));
    tmp = STS(:,whichframe)*STS(:,whichframe)';
    STCross = reshape(STCross(:,whichframe),[6 6]);
    STCs = (nspikes.*STCross-tmp)/(nspikes*(nspikes-1));
    P = eye(size(STCs)) - STAs(:,whichframe)*inv(STAs(:,whichframe)'*STAs(:,whichframe))*STAs(:,whichframe)'; % Subtracting such that the STA is orthogonal to PC
    STCs = P*STCs*P';
    [tmp,d] = eig(STCs);
    [eig_PC,ind] = sort(diag(d)); % storing all the eigenvalues
    [~,tmp_percentileval,tmpzscore] = permutation_test(stro,mask_changes,plot_counter,eig_PC(end),0,whichframe);
    vals(ii) = tmp_percentileval;
    z_scores(ii) = tmpzscore;
    PC1 = tmp(:,ind(end));
    PCS1RGB(1,ii) = PC1(1); PCS2RGB(1,ii) = PC1(2);
    PCS1RGB(2,ii) = PC1(3); PCS2RGB(2,ii) = PC1(4);
    PCS1RGB(3,ii) = PC1(5); PCS3RGB(3,ii) = PC1(6);
    tmp1 = M*PCS1RGB(:,ii); tmp1 = tmp1/(sum(abs(tmp1))); PCS1LMS(:,ii) = tmp1;
    tmp2 = M*PCS2RGB(:,ii); tmp2 = tmp2/(sum(abs(tmp2))); PCS2LMS(:,ii) = tmp2;
end
savevariables = 0;
if savevariables == 1
    save vals vals
    save z_scores z_scores
    save numsubunitspikes numsubunitspikes
    save PCS1RGB PCS1RGB
    save PCS2RGB PCS2RGB
    save PCS1LMS PCS1LMS
    save PCS2LMS PCS2LMS
end
%%
N = numel(filename_c);
load logresidualratio.mat
load anglebwvectorsRGB.mat
load absS.mat
load vals.mat
load z_scores.mat
% for color cells
prettycorr([vals(1:N),z_scores(1:N),logresidualratio(1:N),anglebwvectorsRGB(1:N),absS(1:N)'],{'eig1prctile','zscores','LogResRat','AngleRGB','absS'});
set(gcf,'Name','Color Cells');
% for luminance cells 
prettycorr([vals(N+1:end),z_scores(N+1:end),logresidualratio(N+1:end),anglebwvectorsRGB(N+1:end),absS(N+1:end)'],{'eig1prctile','zscores','LogResRat','AngleRGB','absS'});
set(gcf,'Name','Luminance Cells');