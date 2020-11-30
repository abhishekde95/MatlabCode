% This is an analysis on the WhiteNoise subunit data suggested by Greg.
% I have implemented the analysis for a single neuron in the
% SingleCellAnalysis_RGB3D.m. This script is an extension of the population
% analyses
% Author - Abhishek De, 3/18
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
% [OCidx, LMidx, LUMidx, SOidx, hardtoclassifyidx] = classifycells(S1LMS,S2LMS);
% idx = [OCidx LMidx LUMidx SOidx hardtoclassifyidx];
% idx = [OCidx LMidx];
filename = [filename_c; filename_l];
fileswithisoresponsecontours = numel([filename_c; filename_l]);
% idx = [idx (1:numel([filenameSTAsubunits; filenameSTAvsPC1]))+numel(idx)];
% filename = [filename_c; filename_l; filenameSTAsubunits; filenameSTAvsPC1];
RFstructure =[RFstructure_c; RFstructure_l];
plot_counter = 1;
nspikes = [];
x = linspace(-0.2,0.2,51);
[X,Y,Z] = meshgrid(x,x,x);
nrows = 5;
count = 1;
mostsensitivedirsubunit1 = [];
mostsensitivedirsubunit2 = [];
lowesteigexplainedvarianceS1 = [];
lowesteigexplainedvarianceS2 = [];
mostsensitivedirGQMfit = [];
midsensitivedirGQMfit = [];
leastsensitivedirGQMfit = [];
mostsensitivedirGLMfit = [];
rratio = log10(SSE_linearmodel./SSE_quadmodel);
[~,ss] = sort(rratio); % sorting from lowest to highest rratio values
conicsection = []; % 0 - ellipsoid, 1 - hyperboloid
subunitnonlinearityindex = zeros(numel(filename),2);
significantlynonlinear = zeros(numel(filename),2);
Fprob1 = []; Fprob2 = [];
diffmodelpred1 = []; diffmodelpred2 = [];
LLRprob1 = []; LLRprob2 = [];    
Devprob1 = []; Devprob2 = [];
sigquadterms1 = []; sigquadterms2 = [];
AIClin1 = []; AICquad1 = []; AIClin2 = []; AICquad2 = []; 
BIClin1 = []; BICquad1 = []; BIClin2 = []; BICquad2 = [];
mSSElin1 = []; mSSEquad1 = []; mSSElin2 = []; mSSEquad2 = [];
for ii = 1:1%numel(filename)
    if count == nrows+1
        count = 1;
        plot_counter = plot_counter + 1;
    end
    if ii < fileswithisoresponsecontours
        %         index = ss(ii);
        index = ii;
    else
        index = ii;
    end
    global reversalflagidx stepsizescale stepsize nreversals
    stro = nex2stro(findfile(char(filename(index,:))));
    spikename = 'sig001a'; %getSpikenum(stro);
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
        % Deleting the spikes taking place before the first 9 frames as I need to look at the 9 preceding frames
        spiketimes(spiketimes < maxT*msperframe) = [];
        % Deleting the spikes that take place after the stimulus was
        % removed as it would imply that these spikes do not result from
        % the stimulus shown on the screen
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
    nspikes = [nspikes; out{3}]; % Number of spikes in the given file
    clear STCOV_st out;
    STAs = STS/nspikes(end);
    basisvec1ST = STAs(1:2:end,:); % temporal phosphor dynamics of basisvec 1
    basisvec2ST = STAs(2:2:end,:); % temporal phosphor dynamics of basisvec 2
    [~,whichframe1] = max(fliplr(sum(basisvec1ST.^2,1)));
    [~,whichframe2] = max(fliplr(sum(basisvec2ST.^2,1)));
    stindices = find(cum_n>0);
    stindices1 = stindices - whichframe1+1;
    stindices2 = stindices - whichframe2+1;
    rawrgb_onbasisvec1 = cum_rgbs(1:2:end,:);
    rawrgb_onbasisvec2 = cum_rgbs(2:2:end,:);
    responseS1 = -1*ones(1,size(rawrgb_onbasisvec1,2)); responseS1(stindices1) = 1;
    responseS2 = -1*ones(1,size(rawrgb_onbasisvec2,2)); responseS2(stindices2) = 1;
    muS1spike = mean(rawrgb_onbasisvec1(:,responseS1==1),2);
    muS2spike = mean(rawrgb_onbasisvec2(:,responseS2==1),2);
    img = zeros(size(mask));
    img(mask==1) = muS1spike(1); img(mask==2) = muS2spike(1);
    img(mask==4) = muS1spike(2); img(mask==5) = muS2spike(2);
    img(mask==7) = muS1spike(3); img(mask==8) = muS2spike(3);
    img = reshape(img,[10 10 3]);
    img = 0.5*img/(max(abs(img(:)))+0.001) + 0.5;
    covS1spike = cov(rawrgb_onbasisvec1(:,responseS1==1)');
    covS1 = cov(rawrgb_onbasisvec1');
    covS2spike = cov(rawrgb_onbasisvec2(:,responseS2==1)');
    covS2 = cov(rawrgb_onbasisvec2');
    p1 = mvnpdf([X(:) Y(:) Z(:)],[0 0 0],covS1); p1 = reshape(p1,size(X));
    p2 = mvnpdf([X(:) Y(:) Z(:)],muS1spike',covS1spike); p2 = reshape(p2,size(X));
    p3 = mvnpdf([X(:) Y(:) Z(:)],[0 0 0],covS2); p3 = reshape(p3,size(X));
    p4 = mvnpdf([X(:) Y(:) Z(:)],muS2spike',covS2spike); p4 = reshape(p4,size(X));
    ratio1 = p2./p1;
    ratio2 = p4./p3;
    val1 = prctile(ratio1(:),[10 50 90]);
    val2 = prctile(ratio2(:),[10 50 90]);
    probsS1 = mvnpdf(rawrgb_onbasisvec1',muS1spike',covS1spike)./mvnpdf(rawrgb_onbasisvec1',[0 0 0],covS1);
    probsS2 = mvnpdf(rawrgb_onbasisvec2',muS2spike',covS2spike)./mvnpdf(rawrgb_onbasisvec2',[0 0 0],covS2);
    binsx = linspace(prctile(probsS1,30),prctile(probsS1,70),11);
    binsy = linspace(prctile(probsS2,30),prctile(probsS2,70),11);
    nraw = hist3([probsS1 probsS2],{binsx,binsy}); nraw(nraw==0) = 1;
    nspike = hist3([probsS1(responseS1==1) probsS2(responseS2==1)],{binsx,binsy});
    
    meanraw = [mean(probsS1) mean(probsS2)];
    covraw = cov([probsS1 probsS2]);
    meanspike = [mean(probsS1(responseS1==1)) mean(probsS2(responseS2==1))];
    covspike = cov([probsS1(responseS1==1) probsS2(responseS2==1)]);
    [binsX binsY] = meshgrid(binsx,binsy);
    resp = mvnpdf([binsX(:) binsY(:)],meanspike,covspike)./mvnpdf([binsX(:) binsY(:)],meanraw,covraw);
    resp = reshape(resp,size(binsX));
    figure(plot_counter), subplot(nrows,5,5*count-4); contour(binsX,binsY,resp); axis square; set(gca,'XTick',[],'YTick',[]);
    
    % Comparing GLM and GQM fits
    responseS1(responseS1==-1) = 0;
    responseS2(responseS2==-1) = 0;
    mdl1lin =  fitglm(rawrgb_onbasisvec1',responseS1','linear','Distribution','binomial','Link','logit');
    mdl1quad =  fitglm(rawrgb_onbasisvec1',responseS1','quadratic','Distribution','binomial','Link','logit');
    mdl2lin =  fitglm(rawrgb_onbasisvec2',responseS2','linear','Distribution','binomial','Link','logit');
    mdl2quad =  fitglm(rawrgb_onbasisvec2',responseS2','quadratic','Distribution','binomial','Link','logit');
    
    % Performing some additional analyses on the model fits (GLM and GQM)
    predlin1 = predict(mdl1lin,rawrgb_onbasisvec1'); % perdiction from GLM subunit 1
    predquad1 = predict(mdl1quad,rawrgb_onbasisvec1'); % perdiction from GQM subunit 1
    predlin2 = predict(mdl2lin,rawrgb_onbasisvec2'); % perdiction from GLM subunit 2
    predquad2 = predict(mdl2quad,rawrgb_onbasisvec2'); % perdiction from GQM subunit 2
    prederrorlin1 =  predlin1-responseS1'; % residuals, GLM subunit 1
    prederrorquad1 =  predquad1-responseS1'; % residuals, GQM subunit 1
    prederrorlin2 =  predlin2-responseS2'; % residuals, GLM subunit 2
    prederrorquad2 =  predquad2-responseS2'; % residuals, GQM subunit 2
    SSElin1 = sum(prederrorlin1.^2);
    SSEquad1 = sum(prederrorquad1.^2);
    SSElin2 = sum(prederrorlin2.^2);
    SSEquad2 = sum(prederrorquad2.^2);
    mSSElin1 = [mSSElin1; SSElin1/numel(prederrorlin1)];
    mSSEquad1 = [mSSEquad1; SSEquad1/numel(prederrorquad1)];
    mSSElin2 = [mSSElin2; SSElin2/numel(prederrorlin2)];
    mSSEquad2 = [mSSEquad2; SSEquad2/numel(prederrorquad2)];
    
    % Here I want to try some model comparisons techniques Ftest ANOVA
    Fstatistic1num = (SSElin1 - SSEquad1)/(mdl1quad.NumCoefficients - mdl1lin.NumCoefficients);
    Fstatistic1den = SSEquad1/(numel(responseS1) - mdl1quad.NumCoefficients);
    Fstatistic2num = (SSElin2 - SSEquad2)/(mdl2quad.NumCoefficients - mdl2lin.NumCoefficients);
    Fstatistic2den = SSEquad2/(numel(responseS1) - mdl2quad.NumCoefficients);
    Fprob1 = [Fprob1; 1-fcdf(Fstatistic1num/Fstatistic1den,mdl1quad.NumCoefficients - mdl1lin.NumCoefficients,numel(responseS1) - mdl1quad.NumCoefficients)];
    Fprob2 = [Fprob2; 1-fcdf(Fstatistic2num/Fstatistic2den,mdl2quad.NumCoefficients - mdl2lin.NumCoefficients,numel(responseS2) - mdl2quad.NumCoefficients)];
    
    % Trying to compare the prediction differences between the GLM and GQM
    diffmodelpred1 = [diffmodelpred1; sum((predlin1-predquad1).^2)/numel(predlin1)];
    diffmodelpred2 = [diffmodelpred2; sum((predlin2-predquad2).^2)/numel(predlin2)];
    
    % Log-likelihood ratio test (Wilk's theorem)
    LLRprob1 = [LLRprob1; 1-chi2cdf(-2*(mdl1lin.LogLikelihood - mdl1quad.LogLikelihood),mdl1quad.NumCoefficients - mdl1lin.NumCoefficients)];
    LLRprob2 = [LLRprob2; 1-chi2cdf(-2*(mdl2lin.LogLikelihood - mdl2quad.LogLikelihood),mdl2quad.NumCoefficients - mdl2lin.NumCoefficients)];
    
    %Deviance test
    Devprob1 = [Devprob1; 1-chi2cdf(mdl1lin.Deviance-mdl1quad.Deviance,mdl1quad.NumCoefficients - mdl1lin.NumCoefficients)];
    Devprob2 = [Devprob2; 1-chi2cdf(mdl2lin.Deviance-mdl2quad.Deviance,mdl2quad.NumCoefficients - mdl2lin.NumCoefficients)];
    
    % Storing the number of coefficients of quadratic terms which have a p-value < 0.001
    sigquadterms1 = [sigquadterms1; sum(mdl1quad.Coefficients.pValue(4:end)<0.001)];
    sigquadterms2 = [sigquadterms2; sum(mdl2quad.Coefficients.pValue(4:end)<0.001)];
    
    % Calculating AIC
    AIClin1 = [AIClin1; 2*mdl1lin.NumCoefficients - 2*mdl1lin.LogLikelihood];
    AICquad1 = [AICquad1; 2*mdl1quad.NumCoefficients - 2*mdl1quad.LogLikelihood];
    AIClin2 = [AIClin2; 2*mdl2lin.NumCoefficients - 2*mdl2lin.LogLikelihood];
    AICquad2 = [AICquad2; 2*mdl2quad.NumCoefficients - 2*mdl2quad.LogLikelihood];
    
    % Calculating BIC
    BIClin1 = [BIClin1; mdl1lin.NumCoefficients*log(numel(predlin1)) - 2*mdl1lin.LogLikelihood];
    BICquad1 = [BICquad1; mdl1quad.NumCoefficients*log(numel(predquad1)) - 2*mdl1quad.LogLikelihood];
    BIClin2 = [BIClin2; mdl2lin.NumCoefficients*log(numel(predlin2)) - 2*mdl2lin.LogLikelihood];
    BICquad2 = [BICquad2; mdl2quad.NumCoefficients*log(numel(predquad2)) - 2*mdl2quad.LogLikelihood];
    
    % Comparing GLM and GQM fits on the projection plane of subunit 1 and
    % subunit 2
    
    if index<=fileswithisoresponsecontours
        if isempty(inds)
            inds = size(stro.trial,1)-1;
        end
        neurothreshmode = stro.trial(:,neurothreshidx);
        basisvec_dropidx = inds(end);
        neurothresh_startidx = max(find(neurothreshmode==1,1),basisvec_dropidx+1); % There is possibility of pressing the Neurothresh button multiple times
        num_targetspikerates = unique(stro.trial(neurothresh_startidx:end,targetspikerateidx));
        t_offset = stro.trial(end,latencyidx)/1000;
        vect = stro.ras{basisvec_dropidx,basisvecidx};
        basisvec_size = nstixperside*nstixperside*3;
        numvect = (numel(vect)/basisvec_size)-1;
        basisvec = cell(1,numvect);
        for jj = 1:numvect
            tmp_vec = vect((jj-1)*basisvec_size+1:basisvec_size*jj) + vect(numvect*basisvec_size+1:basisvec_size*(numvect+1));
            basisvec{jj} = reshape(tmp_vec,[nstixperside nstixperside 3]);
        end
        bkgnd_monitor = reshape(vect(numvect*basisvec_size+1:basisvec_size*(numvect+1)),[nstixperside nstixperside 3]);
        norms = cell(1,numel(num_targetspikerates));
        completed_search_alongdir = cell(1,numel(num_targetspikerates));
        for jj = 1:1 % as of now just focusing on 1 target firing rates
            idxs = find(~isnan(stro.trial(:,correctidx)) & stro.trial(:,targetspikerateidx)==num_targetspikerates(jj));
            idxs(idxs<=neurothresh_startidx) = [];
            different_weights = unique(stro.trial(idxs,basisvecdiridx));
            tmp_norm = [];
            tmp_completed_search_alongdir = [];
            for kk = 1:numel(different_weights)
                idxs1 = find(stro.trial(:,basisvecdiridx) == different_weights(kk));
                idxs1(idxs1<neurothresh_startidx) = [];
                raster_data = stro.ras(idxs1,1);
                tmp_norm = [tmp_norm; stro.ras{idxs1(end),weightsidx}'];
                for mm = 1:size(raster_data,1)
                    tmp = raster_data{mm} ;
                    spikes = tmp(tmp>stro.trial(idxs1(mm),stimonidx)+t_offset & tmp < stro.trial(idxs1(mm),stimoffidx));
                    spikes = spikes - stro.trial(idxs1(mm),stimonidx)-t_offset;
                end
                [flag, gamutViolation] = Check_ifcompletely_probed(stro,basisvec,bkgnd_monitor,idxs1,tmp_norm(end,:));
                % flag = 0, incompletely probed
                % flag = 1, completely probed
                % gamutViolation = 1, out of gamut point
                tmp_completed_search_alongdir = [tmp_completed_search_alongdir; flag gamutViolation];
            end
            norms{jj} = tmp_norm;
            completed_search_alongdir{jj} = tmp_completed_search_alongdir;
        end
        basisvec1 = basisvec{1}-bkgnd_monitor;
        basisvec2 = basisvec{2}-bkgnd_monitor;
        tmp = basisvec1 + basisvec2;
        for jj = 1:1
            tmp = norms{jj};
            completed_dir = completed_search_alongdir{jj};
            probed_dirs = logical(completed_dir(:,1)==1); % only including the directions that have been completely probed
            oog_idx = find(completed_dir(:,1)==1 & completed_dir(:,2)==1); % probed and out of gamut
            not_oog_idx = find(completed_dir(:,1)==1 & completed_dir(:,2)==0);
            fact = 0.5./sqrt(tmp(probed_dirs,1).^2 + tmp(probed_dirs,2).^2); % factor needed to extract unit vector
            [THETA,RHO] = cart2pol(tmp(:,1),tmp(:,2));
            ind = (1:numel(THETA))';
            r = fliplr(linspace(0,1,numel(ind)));
            b = fliplr(r);
        end
        if any(THETA>(3*pi/4))
            allthetas = linspace(-pi,pi,100);
            newtheta = linspace(-pi,pi,101);
        else
            allthetas = linspace(-pi/4,3*pi/4,100);
            newtheta = linspace(-pi/4,3*pi/4,101);
        end
        final_model1 = linear_modelparams(index,:);
        rho1 = 1./(final_model1*[cos(allthetas); sin(allthetas)]);
        LOOGtmp1= rho1<0;
        [x_lin,y_lin] = pol2cart(allthetas(~LOOGtmp1),rho1(~LOOGtmp1));
        final_model3 = quad_modelparams(index,:);
        A = final_model3(1); B = final_model3(2); C = final_model3(3); D = final_model3(4); E = final_model3(5);
        rho3 = [];
        p = [A*cos(newtheta').^2+B*sin(newtheta').^2+C*(cos(newtheta').*sin(newtheta')) D*cos(newtheta')+E*sin(newtheta') -1*ones(numel(newtheta'),1)];
        for kk = 1:size(p,1)
            rts = roots(p(kk,:));
            if all(imag(rts)==0) & ~isempty(rts)
                if all(rts>0)
                    r = 100000; %min(rts);
                elseif all(rts<0)
                    r = 100000;
                else
                    r = max(rts);
                end
            else
                r = 100000;
            end
            rho3 = [rho3; r];
        end
        L = rho3>0 & rho3==real(rho3);
        [x_quad,y_quad] = pol2cart(newtheta(L),rho3(L)');
        [x_orig, y_orig] = pol2cart(THETA,RHO);
    end
    
    figure(plot_counter);
    for jj = 1:3
        fv1 = isosurface(X,Y,Z,ratio1,val1(jj));
        fv2 = isosurface(X,Y,Z,ratio2,val2(jj));
        subplot(nrows,5,5*count-3); plot3(fv1.vertices(:,1),fv1.vertices(:,2),fv1.vertices(:,3),'o','MarkerSize',2,'LineWidth',0.5); hold on;
        subplot(nrows,5,5*count-2); plot3(fv2.vertices(:,1),fv2.vertices(:,2),fv2.vertices(:,3),'o','MarkerSize',2,'LineWidth',0.5); hold on;
        if jj == 2
            [v1,d1] = eig(cov(fv1.vertices)); [~,id1] = sort(diag(d1));
            [v2,d2] = eig(cov(fv2.vertices)); [~,id2] = sort(diag(d2));
            wts1 = v1(:,id1(1)); wts1 = wts1/sum(abs(wts1));
            wts2 = v2(:,id2(1)); wts2 = wts2/sum(abs(wts2));
            mostsensitivedirsubunit1 = [mostsensitivedirsubunit1 wts1];
            mostsensitivedirsubunit2 = [mostsensitivedirsubunit2 wts2];
            lowesteigexplainedvarianceS1 = [lowesteigexplainedvarianceS1 100*d1(id1(1),id1(1))/sum(d1(:))];
            lowesteigexplainedvarianceS2 = [lowesteigexplainedvarianceS2 100*d2(id2(1),id2(1))/sum(d2(:))];
        end
    end
    subplot(nrows,5,5*count-3); plot3(2*[-STAs(1,maxT-whichframe1+1) STAs(1,maxT-whichframe1+1)],2*[-STAs(3,maxT-whichframe1+1) STAs(3,maxT-whichframe1+1)],2*[-STAs(5,maxT-whichframe1+1) STAs(5,maxT-whichframe1+1)],'-k','Linewidth',4);
    if Fprob1(end)<0.01
        subunitnonlinearityindex(ii,1) = 1;
        set(gca,'XColor','g','YColor','g','ZColor','g'); hold on;
        coeff = mdl1quad.Coefficients.Estimate;
        Q1 = [coeff(8) coeff(5) coeff(6); coeff(5) coeff(9) coeff(7); coeff(6) coeff(7) coeff(10)];
        [v1,d1] = eig(Q1); [~,id1] = sort(diag(d1));
        wts1 = v1(:,id1(3)); wts1 = wts1/norm(wts1);
        wtsm1 = v1(:,id1(2)); wtsm1 = wtsm1/norm(wtsm1);
        wtsl1 = v1(:,id1(1)); wtsl1 = wtsl1/norm(wtsl1);
        mostsensitivedirGQMfit = [mostsensitivedirGQMfit; wts1'];
        midsensitivedirGQMfit = [midsensitivedirGQMfit; wtsm1'];
        leastsensitivedirGQMfit = [leastsensitivedirGQMfit; wtsl1'];
        if sum(sign(diag(d1)))<0
            conicsection = [conicsection; 1]; % hyperboloid - narrower than linear tuning
        else
            conicsection = [conicsection; 0]; % ellipsoid - broader than linear tuning
        end
        significantlynonlinear(ii,1) = 1;
    else
        coeff = mdl1lin.Coefficients.Estimate(2:end); wts1 = coeff/norm(coeff);
        mostsensitivedirGLMfit = [mostsensitivedirGLMfit; wts1'];
        significantlynonlinear(ii,1) = 0;
    end
    hold off;
    subplot(nrows,5,5*count-2); plot3(2*[-STAs(2,maxT-whichframe2+1) STAs(2,maxT-whichframe2+1)],2*[-STAs(4,maxT-whichframe2+1) STAs(4,maxT-whichframe2+1)],2*[-STAs(6,maxT-whichframe2+1) STAs(6,maxT-whichframe2+1)],'-k','Linewidth',4);
    if Fprob2(end)<0.01
        subunitnonlinearityindex(ii,2) = 1;
        set(gca,'XColor','g','YColor','g','ZColor','g'); hold on;
        coeff = mdl2quad.Coefficients.Estimate;
        Q2 = [coeff(8) coeff(5) coeff(6); coeff(5) coeff(9) coeff(7); coeff(6) coeff(7) coeff(10)];
        [v2,d2] = eig(Q2); [~,id2] = sort(diag(d2));
        wts2 = v2(:,id2(3)); wts2 = wts2/norm(wts2);
        wtsm2 = v2(:,id2(2)); wtsm2 = wtsm2/norm(wtsm2);
        wtsl2 = v2(:,id2(1)); wtsl2 = wtsl2/norm(wtsl2);
        mostsensitivedirGQMfit = [mostsensitivedirGQMfit; wts2'];
        midsensitivedirGQMfit = [midsensitivedirGQMfit; wtsm2'];
        leastsensitivedirGQMfit = [leastsensitivedirGQMfit; wtsl2'];
        if sum(sign(diag(d2)))<0
            conicsection = [conicsection; 1]; % hyperboloid of one sheet - narrower than linear tuning
        else
            conicsection = [conicsection; 0]; % ellipsoid - broader than linear tuning
        end
        significantlynonlinear(ii,2) = 1;
    else
        coeff = mdl1lin.Coefficients.Estimate(2:end); wts2 = coeff/norm(coeff);
        mostsensitivedirGLMfit = [mostsensitivedirGLMfit; wts2'];
        significantlynonlinear(ii,2) = 0;
    end
    hold off;
    subplot(nrows,5,5*count-1); image(img); set(gca,'XTick',[],'YTick',[]); % plots the STA
    if index<=fileswithisoresponsecontours
        subplot(nrows,5,5*count); plot(x_lin,y_lin,'g','Linewidth',2); hold on; plot(x_quad,y_quad,'k','Linewidth',1);
        plot(x_orig(not_oog_idx), y_orig(not_oog_idx),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'PickableParts','none','MarkerEdgeColor',[0 0 1]);
        axis equal; set(gca,'XLim',[-2 2],'YLim',[-2 2]);
        if rratio(index)>1.0
            set(gca,'XColor','g','YColor','g');
        end
        %         if lowesteigexplainedvarianceS1(end)>5 | lowesteigexplainedvarianceS2(end)>5
        %             set(gca,'XColor','g','YColor','g');
        %         end
        hold off;
    end
    count = count + 1;
end

savevariables = 0;
if savevariables == 1
    save Fprob1 Fprob1
    save Fprob2 Fprob2
    save diffmodelpred1 diffmodelpred1
    save diffmodelpred2 diffmodelpred2
    save LLRprob1 LLRprob1
    save LLRprob2 LLRprob2
    save Devprob1 Devprob1
    save Devprob2 Devprob2
    save sigquadterms1 sigquadterms1
    save sigquadterms2 sigquadterms2
    save AIClin1 AIClin1
    save AICquad1 AICquad1
    save AIClin2 AIClin2
    save AICquad2 AICquad2
    save BIClin1 BIClin1
    save BICquad1 BICquad1
    save BIClin2 BIClin2
    save BICquad2 BICquad2
    save mSSElin1 mSSElin1
    save mSSEquad1 mSSEquad1
    save mSSElin2 mSSElin2
    save mSSEquad2 mSSEquad2
end
%%
load newOCidx.mat
load newLMidx.mat
load newSOidx.mat
load newLUMidx.mat
load newhardtoclassifyidx.mat
OCidx = newOCidx;
LMidx = newLMidx;
SOidx = newSOidx;
LUMidx = newLUMidx;
hardtoclassifyidx = newhardtoclassifyidx;
idx = [OCidx LMidx LUMidx SOidx hardtoclassifyidx];

plot_counter = plot_counter + 1;
[r,p] = corr(rratio,(lowesteigexplainedvarianceS1(1:numel(rratio))+lowesteigexplainedvarianceS2(1:numel(rratio)))');
figure(plot_counter); subplot(221); plot(mostsensitivedirsubunit1(1,:),mostsensitivedirsubunit1(2,:),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); xlabel('R'); ylabel('G');
set(gca,'Xlim',[-1 1],'Ylim',[-1 1]);
hold on ;plot(mostsensitivedirsubunit2(1,:),mostsensitivedirsubunit2(2,:),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); xlabel('R'); ylabel('G');
set(gca,'Xlim',[-1 1],'Ylim',[-1 1]); title('Gun wts'); hold off;
subplot(222); plot(log10(SSE_quadmodel),log10(SSE_linearmodel),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); xlabel('SSE quad'); ylabel('SSE lin');
subplot(223); plot(rratio,lowesteigexplainedvarianceS1(1:numel(rratio))+lowesteigexplainedvarianceS2(1:numel(rratio)),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); hold on;
lsline; text(1,30,strcat('r=',num2str(r))); text(1,35,strcat('p=',num2str(p))); xlabel('rratio'); ylabel('Sum lowest eig var');
subplot(224); plot(mostsensitivedirGQMfit(:,1),mostsensitivedirGQMfit(:,2),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); hold on;
plot(mostsensitivedirGLMfit(:,1),mostsensitivedirGLMfit(:,2),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
plot(midsensitivedirGQMfit(:,1),midsensitivedirGQMfit(:,2),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
plot(leastsensitivedirGQMfit(:,1),leastsensitivedirGQMfit(:,2),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]);
xlabel('R'); ylabel('G'); set(gca,'Xlim',[-1 1],'Ylim',[-1 1]);
plot_counter = plot_counter + 1;

% Obtaining the M matrix, code extracted from Greg, fitting a cubic spline
% using the command 'spline'. 'SplineRaw' only availabe through
% psychtoolbox which I currently don't have now.
fundamentals = stro.sum.exptParams.fundamentals; % CONE FUNDAMENTALS: L,M,S
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
mon_spd = stro.sum.exptParams.mon_spd; % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals
bkgndrgb = unique(bkgnd_monitor(:),'stable');
bkgndlms = M*bkgndrgb; % this is in cone excitations
Mrgbtocc = diag(1./bkgndlms)*M;
M = inv(M');

% Now I will convert the gun weights from the GQM and GLM to cone weights. I do see some pattern.
LMSdirsGQMmostsensitivedir = M*mostsensitivedirGQMfit'; LMSdirsGQMmostsensitivedir = LMSdirsGQMmostsensitivedir./repmat(sqrt(sum(LMSdirsGQMmostsensitivedir.^2,1)),[3 1]);
LMSdirsGQMmidsensitivedir = M*midsensitivedirGQMfit'; LMSdirsGQMmidsensitivedir = LMSdirsGQMmidsensitivedir./repmat(sqrt(sum(LMSdirsGQMmidsensitivedir.^2,1)),[3 1]);
LMSdirsGQMleastsensitivedir = M*leastsensitivedirGQMfit'; LMSdirsGQMleastsensitivedir = LMSdirsGQMleastsensitivedir./repmat(sqrt(sum(LMSdirsGQMleastsensitivedir.^2,1)),[3 1]);
LMSdirsGLMmostsensitivedir = M*mostsensitivedirGLMfit'; LMSdirsGLMmostsensitivedir = LMSdirsGLMmostsensitivedir./repmat(sqrt(sum(LMSdirsGLMmostsensitivedir.^2,1)),[3 1]);
figure(plot_counter); plot(LMSdirsGQMmostsensitivedir(1,:),LMSdirsGQMmostsensitivedir(2,:),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); hold on;
plot(LMSdirsGLMmostsensitivedir(1,:),LMSdirsGLMmostsensitivedir(2,:),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
plot(LMSdirsGQMmidsensitivedir(1,:),LMSdirsGQMmidsensitivedir(2,:),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
plot(LMSdirsGQMleastsensitivedir(1,:),LMSdirsGQMleastsensitivedir(2,:),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]);
xlabel('L'); ylabel('M'); set(gca,'Xlim',[-1 1],'Ylim',[-1 1]);
plot_counter = plot_counter + 1;

% Want to plot the proportion of nonlinear subunits within the nonlinear
% double opponent neurons
OC_ii = find(rratio(1:numel(OCidx))>0.4);
LM_ii = find(rratio(numel(OCidx)+1:numel(OCidx)+numel(LMidx))>0.4);
OC_c = hist(sum(significantlynonlinear(OC_ii,:),2),[0 1 2]);
LM_c = hist(sum(significantlynonlinear(LM_ii,:),2),[0 1 2]);
figure(plot_counter), bar([OC_c; LM_c]','stacked'); legend('OC','LM');
plot_counter = plot_counter + 1;

RFstructure = RFstructure(idx);
CS_ii = find(rratio(1:numel(OCidx)+numel(LMidx))>0.4 & RFstructure(1:numel(OCidx)+numel(LMidx))==1);
Ad_ii = find(rratio(1:numel(OCidx)+numel(LMidx))>0.4 & RFstructure(1:numel(OCidx)+numel(LMidx))==2);
CS_c = hist(sum(significantlynonlinear(CS_ii,:),2),[0 1 2]);
Ad_c = hist(sum(significantlynonlinear(Ad_ii,:),2),[0 1 2]);
figure(plot_counter), bar([CS_c; Ad_c]','stacked'); legend('CS','Adjacent');
plot_counter = plot_counter + 1;

% Here I want to see the how many subunits are non-linear on an average for
% linear and non-linear neurons
iiii = find(rratio>0.5); iiiii = find(rratio<=0.5);
a = hist(sum(significantlynonlinear(iiii,:),2),[0 1 2]); a = a/sum(a);
b = hist(sum(significantlynonlinear(iiiii,:),2),[0 1 2]); b = b/sum(b);
figure(plot_counter), bar([a; b]'); legend('non-linear isoresp','linear isoresp');
plot_counter = plot_counter + 1;

% Plotting the proportion of broadly and narrowly tuned subunits
c = hist(conicsection,[0 1]);
figure(plot_counter), bar(categorical({'broad','narrow'}),c); title('color tuning'); ylabel('cell counts');
plot_counter = plot_counter + 1;

% Want to see how z_scores are correlated with linear and non-linear
% isoresponse contours
z_scores = z_scores(idx);
lin_ii  = find(rratio(1:numel(OCidx)+numel(LMidx))<=0.4);
nonlin_ii  = find(rratio(1:numel(OCidx)+numel(LMidx))>0.4);
figure(plot_counter),plot(ones(size(lin_ii)),log(abs(z_scores(lin_ii))),'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]); hold on;
plot(2*ones(size(nonlin_ii)),log(abs(z_scores(nonlin_ii))),'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]); hold on;
set(gca,'Xlim',[0 3]); hold off; plot_counter = plot_counter + 1;



%% Now plotting in 3D RGB and LMS space

% finding the principal axes of all the points in LMS plane
[v,d] = eig(cov([LMSdirsGQMmostsensitivedir'; LMSdirsGLMmostsensitivedir';LMSdirsGQMmidsensitivedir'; LMSdirsGQMleastsensitivedir']));
[vrgb,drgb] = eig(cov([mostsensitivedirGQMfit; mostsensitivedirGLMfit; midsensitivedirGQMfit; leastsensitivedirGQMfit]));
vlms = M*vrgb;
vlms = vlms./repmat(sqrt(sum(vlms.^2,1)),[3 1]);
figure(plot_counter); subplot(121);plot3(mostsensitivedirGQMfit(:,1),mostsensitivedirGQMfit(:,2),mostsensitivedirGQMfit(:,3),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); hold on;
plot3(mostsensitivedirGLMfit(:,1),mostsensitivedirGLMfit(:,2),mostsensitivedirGLMfit(:,3),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
plot3(midsensitivedirGQMfit(:,1),midsensitivedirGQMfit(:,2),midsensitivedirGQMfit(:,3),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
plot3(leastsensitivedirGQMfit(:,1),leastsensitivedirGQMfit(:,2),leastsensitivedirGQMfit(:,3),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]);
plot3([-vrgb(1,1) vrgb(1,1)],[-vrgb(2,1) vrgb(2,1)],[-vrgb(3,1) vrgb(3,1)],'k','Linewidth',2);
plot3([-vrgb(1,2) vrgb(1,2)],[-vrgb(2,2) vrgb(2,2)],[-vrgb(3,2) vrgb(3,2)],'k','Linewidth',2);
plot3([-vrgb(1,3) vrgb(1,3)],[-vrgb(2,3) vrgb(2,3)],[-vrgb(3,3) vrgb(3,3)],'k','Linewidth',2);
xlabel('R'); ylabel('G'); zlabel('B'); title('RGB'); set(gca,'Xlim',[-1 1],'Ylim',[-1 1],'Zlim',[-1 1]);
subplot(122);plot3(LMSdirsGQMmostsensitivedir(1,:),LMSdirsGQMmostsensitivedir(2,:),LMSdirsGQMmostsensitivedir(3,:),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); hold on;
plot3(LMSdirsGLMmostsensitivedir(1,:),LMSdirsGLMmostsensitivedir(2,:),LMSdirsGLMmostsensitivedir(3,:),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
plot3(LMSdirsGQMmidsensitivedir(1,:),LMSdirsGQMmidsensitivedir(2,:),LMSdirsGQMmidsensitivedir(3,:),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
plot3(LMSdirsGQMleastsensitivedir(1,:),LMSdirsGQMleastsensitivedir(2,:),LMSdirsGQMleastsensitivedir(3,:),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]);
plot3([-v(1,1) v(1,1)],[-v(2,1) v(2,1)],[-v(3,1) v(3,1)],'k','Linewidth',2);
plot3([-v(1,2) v(1,2)],[-v(2,2) v(2,2)],[-v(3,2) v(3,2)],'k','Linewidth',2);
plot3([-v(1,3) v(1,3)],[-v(2,3) v(2,3)],[-v(3,3) v(3,3)],'k','Linewidth',2);
xlabel('L'); ylabel('M');zlabel('S'); title('LMS'); set(gca,'Xlim',[-1 1],'Ylim',[-1 1],'Zlim',[-1 1]);
plot_counter = plot_counter + 1;

SisoRGB = inv(Mrgbtocc)*[0 0 1]';
SisoRGB = SisoRGB/norm(SisoRGB);
RedGreenRGB = inv(Mrgbtocc)*[1 -1 0]';
RedGreenRGB = RedGreenRGB/norm(RedGreenRGB);
% LumRGB = inv(Mrgbtocc)*[1.32 1.63 1.04]';
LumRGB = inv(Mrgbtocc)*[1.0 1.0 0]';
LumRGB = LumRGB/norm(LumRGB);
[v1,d1] = eig(cov([mostsensitivedirGQMfit; mostsensitivedirGLMfit])); [~,id1] = sort(diag(d1));
[v2,d2] = eig(cov([midsensitivedirGQMfit; leastsensitivedirGQMfit])); [~,id2] = sort(diag(d2));
[v3,d3] = eig(cov([LMSdirsGQMmostsensitivedir'; LMSdirsGLMmostsensitivedir'])); [~,id3] = sort(diag(d3));
[v4,d4] = eig(cov([LMSdirsGQMmidsensitivedir'; LMSdirsGQMleastsensitivedir'])); [~,id4] = sort(diag(d4));
v1LMS = M*v1(:,id1(1)); v1LMS = v1LMS/norm(v1LMS);

% Want to see the most-sensitive and mid- & least-sensitive directions
% separately in RGB space. I think there is a pattern
figure(plot_counter), subplot(221), plot3(mostsensitivedirGQMfit(:,1),mostsensitivedirGQMfit(:,2),mostsensitivedirGQMfit(:,3),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); hold on;
plot3(mostsensitivedirGLMfit(:,1),mostsensitivedirGLMfit(:,2),mostsensitivedirGLMfit(:,3),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
plot3([-v1(1,id1(1)) v1(1,id1(1))],[-v1(2,id1(1)) v1(2,id1(1))],[-v1(3,id1(1)) v1(3,id1(1))],'k','Linewidth',2);
plot3([-SisoRGB(1) SisoRGB(1)],[-SisoRGB(2) SisoRGB(2)],[-SisoRGB(3) SisoRGB(3)],'b','Linewidth',2);
plot3([-RedGreenRGB(1) RedGreenRGB(1)],[-RedGreenRGB(2) RedGreenRGB(2)],[-RedGreenRGB(3) RedGreenRGB(3)],'r','Linewidth',2);
plot3([-LumRGB(1) LumRGB(1)],[-LumRGB(2) LumRGB(2)],[-LumRGB(3) LumRGB(3)],'g','Linewidth',2);
xlabel('R'); ylabel('G'); zlabel('B'); title('Most-sensitive RGB'); set(gca,'Xlim',[-1 1],'Ylim',[-1 1],'Zlim',[-1 1]); hold off;
subplot(222), plot3(midsensitivedirGQMfit(:,1),midsensitivedirGQMfit(:,2),midsensitivedirGQMfit(:,3),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]); hold on;
plot3(leastsensitivedirGQMfit(:,1),leastsensitivedirGQMfit(:,2),leastsensitivedirGQMfit(:,3),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]);
plot3([-v2(1,id2(1)) v2(1,id2(1))],[-v2(2,id2(1)) v2(2,id2(1))],[-v2(3,id2(1)) v2(3,id2(1))],'k','Linewidth',2);
xlabel('R'); ylabel('G'); zlabel('B'); title('Mid- & least-sensitive RGB'); set(gca,'Xlim',[-1 1],'Ylim',[-1 1],'Zlim',[-1 1]); hold off;
subplot(223), plot3(LMSdirsGQMmostsensitivedir(1,:),LMSdirsGQMmostsensitivedir(2,:),LMSdirsGQMmostsensitivedir(3,:),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); hold on;
plot3(LMSdirsGLMmostsensitivedir(1,:),LMSdirsGLMmostsensitivedir(2,:),LMSdirsGLMmostsensitivedir(3,:),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
plot3([-v3(1,id3(1)) v3(1,id3(1))],[-v3(2,id3(1)) v3(2,id3(1))],[-v3(3,id3(1)) v3(3,id3(1))],'k','Linewidth',2);
plot3([0 0],[0 0],[-1 1],'b','Linewidth',2);
plot3([-v1LMS(1) v1LMS(1)],[-v1LMS(2) v1LMS(2)],[-v1LMS(3) v1LMS(3)],'g','Linewidth',2);
xlabel('L'); ylabel('M');zlabel('S'); title('Most-sensitive LMS'); set(gca,'Xlim',[-1 1],'Ylim',[-1 1],'Zlim',[-1 1]); hold off;
subplot(224), plot3(LMSdirsGQMmidsensitivedir(1,:),LMSdirsGQMmidsensitivedir(2,:),LMSdirsGQMmidsensitivedir(3,:),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]); hold on;
plot3(LMSdirsGQMleastsensitivedir(1,:),LMSdirsGQMleastsensitivedir(2,:),LMSdirsGQMleastsensitivedir(3,:),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]);
plot3([-v4(1,id4(1)) v4(1,id4(1))],[-v4(2,id4(1)) v4(2,id4(1))],[-v4(3,id4(1)) v4(3,id4(1))],'k','Linewidth',2);
xlabel('L'); ylabel('M');zlabel('S'); title('Mid- & least-sensitive LMS'); set(gca,'Xlim',[-1 1],'Ylim',[-1 1],'Zlim',[-1 1]); hold off;
plot_counter = plot_counter + 1;


figure(plot_counter), plot3(mostsensitivedirGQMfit(:,1),mostsensitivedirGQMfit(:,2),mostsensitivedirGQMfit(:,3),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); hold on;
plot3(mostsensitivedirGLMfit(:,1),mostsensitivedirGLMfit(:,2),mostsensitivedirGLMfit(:,3),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]); hold on;
plot3([-v1(1,id1(1)) v1(1,id1(1))],[-v1(2,id1(1)) v1(2,id1(1))],[-v1(3,id1(1)) v1(3,id1(1))],'k','Linewidth',2);
plot3([-SisoRGB(1) SisoRGB(1)],[-SisoRGB(2) SisoRGB(2)],[-SisoRGB(3) SisoRGB(3)],'b','Linewidth',2);
% plot3([-RedGreenRGB(1) RedGreenRGB(1)],[-RedGreenRGB(2) RedGreenRGB(2)],[-RedGreenRGB(3) RedGreenRGB(3)],'r','Linewidth',2);
% plot3([-LumRGB(1) LumRGB(1)],[-LumRGB(2) LumRGB(2)],[-LumRGB(3) LumRGB(3)],'g','Linewidth',2);
xlabel('R'); ylabel('G'); zlabel('B'); title('Most-sensitive RGB'); set(gca,'Xlim',[-1 1],'Ylim',[-1 1],'Zlim',[-1 1]); hold off;
plot_counter = plot_counter + 1;



%% Some control simulation
[x,y,z] = sphere;
xyzLMS = M*[x(:) y(:) z(:)]';
tmp = xyzLMS./repmat(sqrt(sum(xyzLMS.^2,1)),[3 1]);
figure(plot_counter); subplot(221); plot3(x(:),y(:),z(:),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]); xlabel('R'),ylabel('G'),zlabel('B');
subplot(222); plot3(xyzLMS(1,:),xyzLMS(2,:),xyzLMS(3,:),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]); xlabel('L'),ylabel('M'),zlabel('S');
subplot(223); plot3(tmp(1,:),tmp(2,:),tmp(3,:),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]); xlabel('L'),ylabel('M'),zlabel('S');


