% Figures for the cone signal integration paper
% Author - Abhishek De, 08/20
% Right now, writing this script for the 
close all; 
clearvars;
plot_counter = 1;

%% Figure X1: Population plot of the Within signal integration NLI

if ~exist('plot_counter')
    plot_counter = 1;
end

load conewts_svd.mat
load vals.mat
thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];

LUMidx = LumIds_conewts;
DOidx = [ColorOpponentIds_conewts Sconedominated_conewts];
hardtoclassifyidx = [Other_conewts];
hardtoclassifyidx = [hardtoclassifyidx LUMidx(vals(LUMidx)>=95) DOidx(vals(DOidx)>=95)];
LUMidx = LUMidx(vals(LUMidx)<95);
DOidx = DOidx(vals(DOidx)<95);

% Load the integration within the subunit data
load AUROClinS1_CV.mat
load AUROCquadS1_CV.mat
load AUROClinS2_CV.mat
load AUROCquadS2_CV.mat

% Loading the LGN data: calculate the median value 
load ConesignalNLI_LGN.mat

% AN empty list/array for storing the ratio of errors from the GLM/GQM
ConesignalNLI = [];

for ii = 1:numel(AUROCquad1)
    
    % Converting the errors into performances 
    ErrorS1_lin = 1-AUROClin1{ii};
    ErrorS1_quad = 1-AUROCquad1{ii};
    ErrorS2_lin = 1-AUROClin2{ii};
    ErrorS2_quad = 1-AUROCquad2{ii};
    
    tmp = [log10(median(ErrorS1_lin./ErrorS1_quad)) log10(median(ErrorS2_lin./ErrorS2_quad))];
    
    ConesignalNLI = [ConesignalNLI; median(tmp)];
end

indices = [109 31 35];
% Plotting the results for SVD based cone weight classification including the PC1 z-scores 
figure(plot_counter);
subplot(121); histogram(ConesignalNLI(LUMidx),-0.04:0.005:0.1,'Displaystyle','stairs','EdgeColor',[0 0 0]); hold on;
plot(median(ConesignalNLI(LUMidx)),20,'v','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
plot(ConesignalNLI(indices(1)),19,'s','MarkerSize',8,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);

 histogram(ConesignalNLI(DOidx),-0.04:0.005:0.1,'Displaystyle','stairs','EdgeColor',[1 0 0]); hold on;
plot(median(ConesignalNLI(DOidx)),18,'v','MarkerSize',8,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(ConesignalNLI(indices(2)),17,'s','MarkerSize',8,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);

histogram(ConesignalNLI(hardtoclassifyidx),-0.04:0.005:0.1,'Displaystyle', 'stairs','EdgeColor',[0.5 0.5 0.5]); hold on;
plot(median(ConesignalNLI(hardtoclassifyidx)),16,'v','MarkerSize',8,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
plot(ConesignalNLI(indices(3)),15,'s','MarkerSize',8,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-0.03 0.06],'Ylim',[0 20],'YTick',[0 10 20]); xlabel('cone signal NLI'); ylabel('Count'); title('Within subunits'); axis square; hold off;

subplot(122); histogram(ConesignalNLI_LGN,-0.04:0.005:0.1,'Displaystyle','stairs','EdgeColor',[0 0.5 1.0]); hold on;
plot(median(ConesignalNLI_LGN),60,'v','MarkerSize',8,'MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-0.03 0.06],'Ylim',[0 60],'YTick',[0 30 60]); xlabel('cone signal NLI'); ylabel('Count'); title('Within subunits'); axis square; hold off;

set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;

% Comparing spatial NLI across cell types
group = [ones(size(LUMidx)) 2*ones(size(DOidx)) 3*ones(size(hardtoclassifyidx))];
data = ConesignalNLI([LUMidx'; DOidx'; hardtoclassifyidx']); 
p1 = kruskalwallis(data,group,'off');


%% Visualizing the iso-probability of GLM/GQM fits 

if ~exist('plot_counter')
    plot_counter = 1;
end

indices = [109 31 35]; % indices of the example simple, DO and hardtoclassify cell

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
nrows = 3;
count = 1;
AUROClin1 = cell(1,numel(filename)); 
AUROCquad1 = cell(1,numel(filename));
AUROClin2 = cell(1,numel(filename));  
AUROCquad2 = cell(1,numel(filename));
EigenMatrix1 = cell(1,numel(filename));
EigenMatrix2 = cell(1,numel(filename));
spikename_options = ['sig001a'; 'sig001b'];
filename = filename(indices);

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
    val1 = prctile(ratio1(:),[25 50 75]);
    val2 = prctile(ratio2(:),[25 50 75]);
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
        subplot(221); plot3(fv1.vertices(:,1),fv1.vertices(:,2),fv1.vertices(:,3),'o','MarkerSize',2,'LineWidth',0.5, 'MarkerFaceColor',[0 0 1]); hold on;
        subplot(222); plot3(fv2.vertices(:,1),fv2.vertices(:,2),fv2.vertices(:,3),'o','MarkerSize',2,'LineWidth',0.5, 'MarkerFaceColor',[0 0 1]); hold on;
    end
    subplot(221); plot3(3*[-STAs(1,maxT-whichframe1+1) STAs(1,maxT-whichframe1+1)],3*[-STAs(3,maxT-whichframe1+1) STAs(3,maxT-whichframe1+1)],3*[-STAs(5,maxT-whichframe1+1) STAs(5,maxT-whichframe1+1)],'-k','Linewidth',4);
     plot3([0 STAs(1,maxT-whichframe1+1)],[0 STAs(3,maxT-whichframe1+1)],[0 STAs(5,maxT-whichframe1+1)],'ro','Linewidth',4);
    xlabel('R'), ylabel('G'), zlabel('B'); axis square; axis xy; hold off;
    subplot(222); plot3(3*[-STAs(2,maxT-whichframe1+1) STAs(2,maxT-whichframe1+1)],3*[-STAs(4,maxT-whichframe1+1) STAs(4,maxT-whichframe1+1)],3*[-STAs(6,maxT-whichframe1+1) STAs(6,maxT-whichframe1+1)],'-k','Linewidth',4); 
    xlabel('R'), ylabel('G'), zlabel('B'); axis square; axis xy; hold off;
    
    
    subplot(223), plot3(rawrgb_onbasisvec1(1,responseS1==1), rawrgb_onbasisvec1(2,responseS1==1),rawrgb_onbasisvec1(3,responseS1==1),'g.'); hold on;
    plot3(rawrgb_onbasisvec1(1,responseS1==0), rawrgb_onbasisvec1(2,responseS1==0),rawrgb_onbasisvec1(3,responseS1==0),'b.');
    xlabel('R'), ylabel('G'), zlabel('B'); axis square; set(gca,'Xlim',[-0.5 0.5], 'Ylim',[-0.5 0.5],'Zlim',[-0.5 0.5]); hold off;
    
    subplot(224), plot3(rawrgb_onbasisvec2(1,responseS2==1), rawrgb_onbasisvec1(2,responseS2==1),rawrgb_onbasisvec1(3,responseS2==1),'g.'); hold on;
    plot3(rawrgb_onbasisvec2(1,responseS2==0), rawrgb_onbasisvec2(2,responseS2==0),rawrgb_onbasisvec2(3,responseS2==0),'b.');
    xlabel('R'), ylabel('G'), zlabel('B'); axis square; set(gca,'Xlim',[-0.5 0.5], 'Ylim',[-0.5 0.5],'Zlim',[-0.5 0.5]); hold off;

    
    plot_counter = plot_counter + 1;
    count = count + 1;
end


