% Waveform analyses- Fitting GLM to evaluate burstiness-based on Jude Mitchell's suggestion
% Author- Abhishek De, 2/20

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


%% Trying out for 1 cell first

% Essential GLM parameters for Alison's code 
kbasprs.neye = 0; % number of "identity" basis vectors near time of spike;
kbasprs.ncos = 7; % number of raised-cosine vectors to use
kbasprs.b = 10; % how nonlinear to make spacings (larger -> more linear)

%%% basis functions for post-spike kernel
ihbasprs.ncols = 7;  % number of basis vectors for post-spike kernel
ihbasprs.hpeaks = [.1 100];  % peak location for first and last vectors, in ms
ihbasprs.b = 1;  % how nonlinear to make spacings (larger -> more linear)
ihbasprs.absref = 1; % absolute refractory period, in ms

% Other parameters for the GLM
softRect = 0;    % use exponential nonlinearity
plotFlag = 1;    % plot fit
saveFlag = 0;    % save fit to fid, in new folder
maxIter = 1000;  % max number of iterations for fitting, also used for maximum number of function evaluations(MaxFunEvals)
tolFun = 1e-6;  % function tolerance for fitting
L2pen = 0;       % penalty on L2-norm of parameter coefficients

% Some WN paramaters
global maxT
numcells = 2;
plot_counter = 1;
maxT = 15; % this represents the temporal part in the spatiotemporal receptive field
temporal_filters = cell(numcells,1);
postspike_filters = cell(numcells,1);
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
    maskidx = strcmp(WN.sum.rasterCells(1,:), 'subunit_mask');
    anlgStartTimeidx = find(strcmp(WN.sum.rasterCells(1,:),'anlgStartTime'));
    seedidx = strcmp(WN.sum.trialFields(1,:),'seed');
    L = WN.trial(:,noisetypeidx)==1;
    mu1idx = find(strcmp(WN.sum.trialFields(1,:),'mu1'));
    mu2idx = find(strcmp(WN.sum.trialFields(1,:),'mu2'));
    mu3idx = find(strcmp(WN.sum.trialFields(1,:),'mu3'));
    sigma1idx = find(strcmp(WN.sum.trialFields(1,:),'sigma1'));
    sigma2idx = find(strcmp(WN.sum.trialFields(1,:),'sigma2'));
    sigma3idx = find(strcmp(WN.sum.trialFields(1,:),'sigma3'));
    maskidx = strcmp(WN.sum.rasterCells(1,:),'subunit_mask');
    basisvecidx = strcmp(WN.sum.rasterCells(1,:),'basis_vec');
    muidxs = [find(strcmp(WN.sum.trialFields(1,:),'mu1')), ...
        find(strcmp(WN.sum.trialFields(1,:),'mu2')), ...
        find(strcmp(WN.sum.trialFields(1,:),'mu3'))];
    sigmaidxs = [find(strcmp(WN.sum.trialFields(1,:),'sigma1')), ...
        find(strcmp(WN.sum.trialFields(1,:),'sigma2')), ...
        find(strcmp(WN.sum.trialFields(1,:),'sigma3'))];
    msperframe = 1000/WN.sum.exptParams.framerate;
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
    ngammasteps = 2^16; % 65536
    xx = linspace(WN.sum.exptParams.gauss_locut/1000, WN.sum.exptParams.gauss_hicut/1000,ngammasteps); % xx represents the probabilities. For more info, have a look at the MATLAB 'norminv' function.
    yy = norminv(xx');
    
    
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
    
    % Pulling out the stimulus frames and responses- Adapted from WN code
    Stim = []; % for storing stimulus values
    response = []; % for storing the response values
    STCOVmex('init', {nstixperside.^2 3 maxT}); % initializing STCOVmex 
   
    for k = mask_changes(1):mask_changes(2)
        nframes = WN.trial(k,nframesidx);
        if (nframes == 0)
            continue;
        end
        seed = WN.trial(k,seedidx);
        mu = WN.trial(k,muidxs)/1000;
        sigma = WN.trial(k,sigmaidxs)/1000;
        org_mask = WN.ras{k,maskidx};
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
        t_stimon = WN.trial(k, stimonidx);
        spiketimes = (WN.ras{k,spikeidx}-t_stimon)*1000; % observing spiketimes in milliseconds
        frametimes = linspace(0, nframes*msperframe, nframes)+(msperframe/2)';
        spiketimes(spiketimes < maxT*msperframe) = [];
        spiketimes(spiketimes > frametimes(end)) = [];
        n = hist(spiketimes, frametimes);
        STCOVmex(rgbs(:),n);
        
        % Appending the stim and response array
        Stim = [Stim rgbs];
        response = [response n];
    end
    out_gun = STCOVmex('return');
    STS_gun = out_gun{1}; STCross_gun = out_gun{2}; nspikes_gun = out_gun{3}; clear out_gun;
    STAs_gun = STS_gun/nspikes_gun; % This might act as an initial guess 
    
    % Calculating the weighted STA
    peakframe = sum(STAs_gun.^2,1) == max(sum(STAs_gun.^2,1));
    id = find(peakframe==1);
    if id~=1
        peakframe(id-1)= 1;
    end
    if id <=maxT-1
        peakframe(id+1)=1;
    end
    
    STAweights = sqrt(sum(STAs_gun(:,peakframe).^2));
    STAweights = STAweights./sum(STAweights);
    spatiochromatic_STA = STAs_gun(:,peakframe)*STAweights'; % Initial guess for spatiochromatic STA
    temporal_STA = sum(STAs_gun.^2,1); % Initial guess for temporal STA
    
%     spatiochromatic_STA = 2*rand([300 1])-1;
    
    %GLM fitting using Alison Weber's code
    nkt = maxT; % number of ms in stim filter
    kbasprs.kpeaks = [.1 round(nkt/1.2)];  % position of first and last bump (relative to identity bumps)
    dt = msperframe; % in ms
    [s, k, h, dc, prs, kbasis, hbasis] = fit_glm_AD(Stim',response',dt,nkt,kbasprs,ihbasprs,[],softRect,plotFlag,maxIter,tolFun,L2pen,spatiochromatic_STA,temporal_STA);
    
    % Displaying the STA
    sSTA = reshape(spatiochromatic_STA,[10 10 3]);
    sSTA = 0.5 + 0.5*sSTA/(max(abs(sSTA(:)))+eps);
    figure(plot_counter);
    subplot(5,4,4*ii-3); image(sSTA); axis square; set(gca,'Tickdir','out','XTick',[],'YTick',[]); hold off;
    
    % Displaying the GLM spatial filter
    s = reshape(s,[10 10 3]);
    s1 = s;
    s1 = 0.5 + 0.5*s1/(max(abs(s1(:)))+eps);
    figure(plot_counter);
    subplot(5,4,4*ii-2); image(s1); axis square; set(gca,'Tickdir','out','XTick',[],'YTick',[]); hold off;
    
    
    % plotting the results of the pre-spike filter 
    figure(plot_counter);
    subplot(5,4,4*ii-1); plot(flipud(k),'k','Linewidth',2);  
    axis square; set(gca,'Tickdir','out','Xlim',[1 maxT]); hold off;
    
    % plotting the results of the post-spike filter 
    figure(plot_counter);
    subplot(5,4,4*ii); plot(-1*h,'r','Linewidth',2); 
    axis square; set(gca,'Tickdir','out','Xlim',[1 maxT]); hold off;
    
    % Storing the variables 
    temporal_filters{ii} = k;
    postspike_filters{ii} = h;
end
plot_counter = plot_counter + 1;

%% Some analysis on the post-spike filters 

