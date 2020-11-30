% A follow up on the PopCellAnalysis_RGB3D.m
% Implementing the cross validation on the the RGB data of each subunits
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
x = linspace(-0.2,0.2,51);
[X,Y,Z] = meshgrid(x,x,x);
nrows = 5;
count = 1;
AUROClin1 = cell(1,numel(filename)); 
AUROCquad1 = cell(1,numel(filename));
AUROClin2 = cell(1,numel(filename));  
AUROCquad2 = cell(1,numel(filename));
EigenMatrix1 = cell(1,numel(filename));
EigenMatrix2 = cell(1,numel(filename));
spikename_options = ['sig001a'; 'sig001b'];

for ii = 1:numel(filename)
    disp(ii);
    if count == nrows+1
        count = 1;
        plot_counter = plot_counter + 1;
    end
    
    global reversalflagidx stepsizescale stepsize nreversals
    stro = nex2stro(findfile(char(filename(ii,:))));
    spikename = 'sig001a';
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
    
    % Need to implement the cross-validation here: 2 cross vals for 2 subunits
    responseS1(responseS1==-1) = 0; 
    responseS1 = logical(responseS1);
    responseS2(responseS2==-1) = 0;
    responseS2 = logical(responseS2);
    C1 = cvpartition(responseS1,'KFold',10);
    C2 = cvpartition(responseS2,'KFold',10);
    tmp_AUROClin1 = []; tmp_AUROCquad1 = [];
    tmp_AUROClin2 = []; tmp_AUROCquad2 = [];

    % Comparing GLM and GQM fits
    for jj = 1:C1.NumTestSets
        mdl1lin =  fitglm(rawrgb_onbasisvec1(:,C1.training(jj))',responseS1(C1.training(jj))','linear','Distribution','binomial','Link','logit');
        mdl1quad =  fitglm(rawrgb_onbasisvec1(:,C1.training(jj))',responseS1(C1.training(jj))','quadratic','Distribution','binomial','Link','logit');
        mdl2lin =  fitglm(rawrgb_onbasisvec2(:,C2.training(jj))',responseS2(C2.training(jj))','linear','Distribution','binomial','Link','logit');
        mdl2quad =  fitglm(rawrgb_onbasisvec2(:,C2.training(jj))',responseS2(C2.training(jj))','quadratic','Distribution','binomial','Link','logit');
        
        % Performing some additional analyses on the model fits (GLM and GQM)
        predlin1 = predict(mdl1lin,rawrgb_onbasisvec1(:,C1.test(jj))'); % perdiction from GLM subunit 1
        predquad1 = predict(mdl1quad,rawrgb_onbasisvec1(:,C1.test(jj))'); % perdiction from GQM subunit 1
        predlin2 = predict(mdl2lin,rawrgb_onbasisvec2(:,C2.test(jj))'); % perdiction from GLM subunit 2
        predquad2 = predict(mdl2quad,rawrgb_onbasisvec2(:,C2.test(jj))'); % perdiction from GQM subunit 2
        
        % Quantifying accuracy using AUROC 
        tmp_AUROClin1 = [tmp_AUROClin1; rocN(predlin1(responseS1(C1.test(jj))),predlin1(~responseS1(C1.test(jj))))];
        tmp_AUROCquad1 = [tmp_AUROCquad1; rocN(predquad1(responseS1(C1.test(jj))),predquad1(~responseS1(C1.test(jj))))];
        tmp_AUROClin2 = [tmp_AUROClin2; rocN(predlin2(responseS2(C2.test(jj))),predlin2(~responseS2(C2.test(jj))))];
        tmp_AUROCquad2 = [tmp_AUROCquad2; rocN(predquad2(responseS2(C2.test(jj))),predquad2(~responseS2(C2.test(jj))))];
    end
    
    AUROClin1{ii} = tmp_AUROClin1;
    AUROCquad1{ii} = tmp_AUROCquad1;
    AUROClin2{ii} = tmp_AUROClin2;
    AUROCquad2{ii} = tmp_AUROCquad2;
    
    % Fitting the full model to extract the eigen directions
    mdl1quad_full =  fitglm(rawrgb_onbasisvec1',responseS1','quadratic','Distribution','binomial','Link','logit');
    mdl2quad_full =  fitglm(rawrgb_onbasisvec2',responseS2','quadratic','Distribution','binomial','Link','logit');
    Q1 = mdl1quad_full.Coefficients.Estimate(5:end);
    Q2 = mdl2quad_full.Coefficients.Estimate(5:end);
    EigenMatrix1{ii} = [Q1(6) Q1(1)/2 Q1(2)/2; Q1(1)/2 Q1(5) Q1(3)/2; Q1(2)/2 Q1(3)/2 Q1(4)];
    EigenMatrix2{ii} = [Q2(6) Q2(1)/2 Q2(2)/2; Q2(1)/2 Q2(5) Q2(3)/2; Q2(2)/2 Q2(3)/2 Q2(4)];
    
    figure(plot_counter);
    for jj = 1:3
        fv1 = isosurface(X,Y,Z,ratio1,val1(jj));
        fv2 = isosurface(X,Y,Z,ratio2,val2(jj));
        subplot(nrows,2,2*count-1); plot3(fv1.vertices(:,1),fv1.vertices(:,2),fv1.vertices(:,3),'o','MarkerSize',2,'LineWidth',0.5); hold on;
        subplot(nrows,2,2*count); plot3(fv2.vertices(:,1),fv2.vertices(:,2),fv2.vertices(:,3),'o','MarkerSize',2,'LineWidth',0.5); hold on;
    end
    subplot(nrows,2,2*count-1); plot3(2*[-STAs(1,maxT-whichframe1+1) STAs(1,maxT-whichframe1+1)],2*[-STAs(3,maxT-whichframe1+1) STAs(3,maxT-whichframe1+1)],2*[-STAs(5,maxT-whichframe1+1) STAs(5,maxT-whichframe1+1)],'-k','Linewidth',4); axis square; hold off;
    subplot(nrows,2,2*count); plot3(2*[-STAs(2,maxT-whichframe1+1) STAs(2,maxT-whichframe1+1)],2*[-STAs(4,maxT-whichframe1+1) STAs(4,maxT-whichframe1+1)],2*[-STAs(6,maxT-whichframe1+1) STAs(6,maxT-whichframe1+1)],'-k','Linewidth',4); axis square; hold off;
    count = count + 1;
    
    keyboard;
end

%%
savevariables = 0; 
if savevariables 
    save AUROClinS1_CV AUROClin1
    save AUROCquadS1_CV AUROCquad1
    save AUROClinS2_CV AUROClin2
    save AUROCquadS2_CV AUROCquad2
    save EigenMatrix1 EigenMatrix1
    save EigenMatrix2 EigenMatrix2
end

AUROClin1_median = [];
AUROCquad1_median = [];
AUROClin2_median = [];
AUROCquad2_median = [];
differenceofmedians1 = [];
differenceofmedians2 = [];

for ii = 1:numel(AUROClin1)
    AUROClin1_median = [AUROClin1_median; median(AUROClin1{ii})];
    AUROCquad1_median = [AUROCquad1_median; median(AUROCquad1{ii})];
    AUROClin2_median = [AUROClin2_median; median(AUROClin2{ii})];
    AUROCquad2_median = [AUROCquad2_median; median(AUROCquad2{ii})];
    differenceofmedians1 = [differenceofmedians1; median(AUROCquad1{ii}-AUROClin1{ii})];
    differenceofmedians2 = [differenceofmedians2; median(AUROCquad2{ii}-AUROClin2{ii})];
end

% Plotting the scatterplots
figure(plot_counter);
subplot(221); plot(AUROClin1_median,AUROCquad1_median,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); set(gca,'Tickdir','out','Xlim',[0.5 1],'Ylim',[0.5 1]); 
axis square; xlabel('Lin S1'); ylabel('Quad S1'); title('Median error');
subplot(222); plot(AUROClin2_median,AUROCquad2_median,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); set(gca,'Tickdir','out','Xlim',[0.5 1],'Ylim',[0.5 1]); 
axis square; xlabel('Lin S2'); ylabel('Quad S2'); title('Median error');


% Plotting the histogram of the differences between the linear and quadratic predictions
bins = linspace(-2,8,31);
subplot(223), histogram(100*(differenceofmedians1),bins,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
set(gca,'Tickdir','out'); title('Median ROC: S1'); xlabel('Quad-Lin AUROC'); axis square;
subplot(224), histogram(100*(differenceofmedians2),bins,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
set(gca,'Tickdir','out'); title('Median ROC: S2'); xlabel('Quad-Lin AUROC'); axis square;
plot_counter = plot_counter + 1;



