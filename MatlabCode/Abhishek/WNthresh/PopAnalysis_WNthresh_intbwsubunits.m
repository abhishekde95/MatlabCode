% A new script to quantitatively figure out the interaction between the
% subunits from the Whitenoise subunit data
% Author - Abhishek De, 10/18
close all;
clearvars;
close all; clearvars;
load filename_c.mat % Color cells
load filename_l.mat % Luminance cells
load filenameSTAvsPC1.mat
load filenameSTAsubunits.mat
load S1LMS.mat % cone wts for subunit 1
load S2LMS.mat % cone wts for subunit 2
load linear_modelparams.mat % parameters from the linear fit
load quad_modelparams.mat % parameters from the quadratic fit
load SSE_linearmodel.mat % sum of squared errors from the linear fit
load SSE_quadmodel.mat % sum of squared errors from the quadratic fit
load RFstructure_c.mat % RF color cells
load RFstructure_l.mat % RF luminance cells
load z_scores.mat
filename = [filename_c; filename_l];
fileswithisoresponsecontours = numel([filename_c; filename_l]);
RFstructure =[RFstructure_c; RFstructure_l];
plot_counter = 1;
nspikes = [];
nrows = 5;
count = 1;
mSSElin = []; mSSEquad = [];
Fstatisticnum = []; Fstatisticden = [];
Fprob = []; diffmodelpred = []; LLRprob = []; Devprob = []; sigquadterms = [];
AIClin = []; AICquad = []; BIClin = []; BICquad = [];
for ii = 1:1%numel(filename)
    if count == nrows+1
        count = 1;
        plot_counter = plot_counter + 1;
    end
    if ii < fileswithisoresponsecontours
        index = ii;
    else
        index = ii;
    end
    global reversalflagidx stepsizescale stepsize nreversals nframesidx maxT
    stro = nex2stro(findfile(char(filename(index,:))));
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
    % Fit a GLM & GQM after the frames have been projected onto the subunits
    mdllin =  fitglm(projs,Lspike,'linear','Distribution','binomial','Link','logit');
    mdlquad =  fitglm(projs,Lspike,'quadratic','Distribution','binomial','Link','logit');
    
    % Performing some additional analyses on the model fits (GLM and GQM)
    predlin = predict(mdllin,projs); % perdiction from GLM
    predquad = predict(mdlquad,projs); % perdiction from GQM 
    prederrorlin =  predlin-Lspike; % residuals, GLM 
    prederrorquad =  predquad-Lspike; % residuals, GQM
   
    SSElin = sum(prederrorlin.^2);
    SSEquad = sum(prederrorquad.^2);
   
    mSSElin = [mSSElin; SSElin/numel(prederrorlin)];
    mSSEquad = [mSSEquad; SSEquad/numel(prederrorquad)];
    
    % Here I want to try some model comparisons techniques Ftest ANOVA
    Fstatisticnum = (SSElin - SSEquad)/(mdlquad.NumCoefficients - mdllin.NumCoefficients);
    Fstatisticden = SSEquad/(numel(Lspike) - mdlquad.NumCoefficients);
    Fprob = [Fprob; 1-fcdf(Fstatisticnum/Fstatisticden,mdlquad.NumCoefficients - mdllin.NumCoefficients,numel(Lspike) - mdlquad.NumCoefficients)];
        
    % Trying to compare the prediction differences between the GLM and GQM
    diffmodelpred = [diffmodelpred; sum((predlin-predquad).^2)/numel(predlin)];
    
    % Log-likelihood ratio test (Wilk's theorem)
    LLRprob = [LLRprob; 1-chi2cdf(-2*(mdllin.LogLikelihood - mdlquad.LogLikelihood),mdlquad.NumCoefficients - mdllin.NumCoefficients)];
    
    %Deviance test
    Devprob = [Devprob; 1-chi2cdf(mdllin.Deviance-mdlquad.Deviance,mdlquad.NumCoefficients - mdllin.NumCoefficients)];
    
    % Storing the number of coefficients of quadratic terms which have a p-value < 0.001
    sigquadterms = [sigquadterms; sum(mdlquad.Coefficients.pValue(4:end)<0.001)];
    
    % Calculating AIC
    AIClin = [AIClin; 2*mdllin.NumCoefficients - 2*mdllin.LogLikelihood];
    AICquad = [AICquad; 2*mdlquad.NumCoefficients - 2*mdlquad.LogLikelihood];
    
    % Calculating BIC
    BIClin = [BIClin; mdllin.NumCoefficients*log(numel(predlin)) - 2*mdllin.LogLikelihood];
    BICquad = [BICquad; mdlquad.NumCoefficients*log(numel(predquad)) - 2*mdlquad.LogLikelihood];
    
    % Calculating the actual histogram
    divs = 11;
    X_low = prctile(projs(:,1),1); X_high = prctile(projs(:,1),99);
    Y_low = prctile(projs(:,2),1); Y_high = prctile(projs(:,2),99);
    [Nraw,Xedges,Yedges] = histcounts2(projs(:,1),projs(:,2),linspace(X_low,X_high,divs),linspace(Y_low,Y_high,divs));
    Nspike = histcounts2(projs(Lspike,1),projs(Lspike,2),linspace(X_low,X_high,divs),linspace(Y_low,Y_high,divs));
    
    [X,Y] = meshgrid(linspace(X_low,X_high,divs),linspace(Y_low,Y_high,divs)); 
    resp1 = predict(mdllin,[X(:) Y(:)]); resp1 = reshape(resp1,size(X));
    resp2 = predict(mdlquad,[X(:) Y(:)]); resp2 = reshape(resp2,size(X));
    figure(plot_counter),
    subplot(nrows,3,3*count-2); imagesc(Xedges,Yedges,Nspike'./Nraw'); set(gca,'Tickdir','out','XTick',[X_low X_high],'YTick',[Y_low Y_high]); colormap('gray'); axis square; axis xy; xlabel('S1'); ylabel('S2'); 
    subplot(nrows,3,3*count-1); contour(Xedges,Yedges,resp1); set(gca,'Tickdir','out','XTick',[X_low X_high],'YTick',[Y_low Y_high]); axis equal; xlabel('S1'); ylabel('S2'); 
    subplot(nrows,3,3*count); contour(Xedges,Yedges,resp2); set(gca,'Tickdir','out','XTick',[X_low X_high],'YTick',[Y_low Y_high]); axis equal; xlabel('S1'); ylabel('S2');
    
    count = count + 1;
end

savevariables = 0;
if savevariables == 1
    save Fprob Fprob
    save diffmodelpred diffmodelpred
    save LLRprob LLRprob
    save Devprob Devprob    
    save sigquadterms sigquadterms
    save AIClin AIClin
    save AICquad AICquad
    save BIClin BIClin
    save BICquad BICquad
    save mSSElin mSSElin
    save mSSEquad mSSEquad
end