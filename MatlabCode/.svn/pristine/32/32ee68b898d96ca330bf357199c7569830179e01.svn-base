% Derivative of PopAnalysis_WNthresh_intbwsubunits_crossval.m
% Implementing crossvalidation 
% Author - Abhishek De, 12/19
close all; clearvars;

conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filename = fetch(conn,'SELECT filename FROM WNthresh');
NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
close(conn);
filename = filename(strcmp(string(NTmode),"subunit"));
NTmode = NTmode(strcmp(string(NTmode),"subunit"));
spikeidx_NT = spikeidx_NT(strcmp(string(NTmode),"subunit"));

plot_counter = 1;
nspikes = [];
nrows = 5;
count = 1;
AUROClinsubunits = cell(1,numel(filename)); 
AUROCquadsubunits = cell(1,numel(filename));

for ii = 1:numel(filename)
    disp(ii);
    if count == nrows+1
        count = 1;
        plot_counter = plot_counter + 1;
    end
    global reversalflagidx stepsizescale stepsize nreversals nframesidx maxT
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
    
    if num_subunits>2
        keyboard;
    end
    
    for mode = 1:2 % first round for calculating the STAs and the second round for projecting them onto the 2 subunits within the STAs
        if mode == 1
            if use_STCOVmex_ST
                STCOV_st('init', {num_subunits 3 maxT});
            else
                STCOVmex('init', {num_subunits 3 maxT});
            end
        else
            initargs = define_basis_vec(nrandnums_perchannel,stro,STAs,[],0);
            STPROJmod('init',initargs); % initialising the STPROJmod
        end
             
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
            if mode == 1
                if use_STCOVmex_ST
                    STCOV_st(rgbs(:),n);
                else
                    STCOVmex(rgbs(:),n);
                end  
            else
                STPROJmod(rgbs(:),n);
            end                 
        end
        if mode == 1
            if use_STCOVmex_ST
                out = STCOV_st('return');
            else
                out = STCOVmex('return');
            end
            STS = out{1};  % A (dimension) x 9(frames) matrix
            nspikes = [nspikes; out{3}]; % Number of spikes in the given file
            clear STCOV_st out;
            STAs = STS/nspikes(end);
            STAs = fliplr(STAs);
        else
            out = STPROJmod('return');
            projs = out{1};
            Lspike = logical(out{2});
            clear STPROJmod out;
        end
    end
    
    % Need to implement a cross-validation here 
    C = cvpartition(Lspike,'KFold',10);
    tmp_AUROClin = [];
    tmp_AUROCquad = [];
    for jj = 1:C.NumTestSets
        % Fit a GLM & GQM after the frames have been projected onto the subunits
        mdllin =  fitglm(projs(C.training(jj),:),Lspike(C.training(jj)),'linear','Distribution','binomial','Link','logit');
        mdlquad =  fitglm(projs(C.training(jj),:),Lspike(C.training(jj)),'quadratic','Distribution','binomial','Link','logit');
        
        % Performing some additional analyses on the model fits (GLM and GQM)
        predlin = predict(mdllin,projs(C.test(jj),:)); % perdiction from GLM
        predquad = predict(mdlquad,projs(C.test(jj),:)); % perdiction from GQM
        tmp_AUROClin = [tmp_AUROClin; rocN(predlin(Lspike(C.test(jj))),predlin(~Lspike(C.test(jj))))];
        tmp_AUROCquad = [tmp_AUROCquad; rocN(predquad(Lspike(C.test(jj))),predquad(~Lspike(C.test(jj))))];      
    end
    
    AUROClinsubunits{ii} = tmp_AUROClin;
    AUROCquadsubunits{ii} = tmp_AUROCquad;
    
    % Calculating the actual histogram
    divs = 11;
    X_low = prctile(projs(:,1),5); X_high = prctile(projs(:,1),95);
    Y_low = prctile(projs(:,2),5); Y_high = prctile(projs(:,2),95);
    [Nraw,Xedges,Yedges] = histcounts2(projs(:,1),projs(:,2),linspace(X_low,X_high,divs),linspace(Y_low,Y_high,divs));
    Nspike = histcounts2(projs(Lspike,1),projs(Lspike,2),linspace(X_low,X_high,divs),linspace(Y_low,Y_high,divs));
    
    [X,Y] = meshgrid(linspace(X_low,X_high,divs),linspace(Y_low,Y_high,divs)); 
    resp1 = predict(mdllin,[X(:) Y(:)]); resp1 = reshape(resp1,size(X));
    resp2 = predict(mdlquad,[X(:) Y(:)]); resp2 = reshape(resp2,size(X));
    figure(plot_counter),
    subplot(nrows,3,3*count-2); imagesc(Xedges,Yedges,Nspike'./Nraw'); set(gca,'Tickdir','out','XTick',[X_low X_high],'YTick',[Y_low Y_high]); colormap('gray'); axis square; axis xy; xlabel('S1'); ylabel('S2'); 
    subplot(nrows,3,3*count-1); contour(Xedges,Yedges,resp1); set(gca,'Tickdir','out','XTick',[],'YTick',[]); axis equal; xlabel('S1'); ylabel('S2'); 
    subplot(nrows,3,3*count); contour(Xedges,Yedges,resp2); set(gca,'Tickdir','out','XTick',[],'YTick',[]); axis equal; xlabel('S1'); ylabel('S2');
    
    count = count + 1;
end

savevariables = 0;
if savevariables
    save AUROClinsubunits_CV AUROClinsubunits
    save AUROCquadsubunits_CV AUROCquadsubunits
end

%% Further analyses on the crossvalidated-AUROC 

AUROClin_median = [];
AUROCquad_median = [];
diffofmedians = [];
for ii = 1:numel(AUROClinsubunits)
    AUROClin_median = [AUROClin_median; median(AUROClinsubunits{ii})];
    AUROCquad_median = [AUROCquad_median; median(AUROCquadsubunits{ii})];
    diffofmedians = [diffofmedians; median(AUROCquadsubunits{ii}-AUROClinsubunits{ii})];
end

bins = linspace(-2,8,31);
figure(plot_counter);
subplot(121); plot(100*AUROClin_median,100*AUROCquad_median,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
set(gca,'Tickdir','out','Xlim',[50 100],'Ylim',[50 100]); title('Median ROC'); xlabel('Lin'); ylabel('Quad'); axis square;
subplot(122),histogram(100*(diffofmedians),bins,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
set(gca,'Tickdir','out'); title('Median ROC'); xlabel('Quad-Lin AUROC'); axis square;

plot_counter = plot_counter + 1;
