% Some more analyses based on friendly reviewer's comments
% Author - Abhishek De, 9/20

% Section 1: Smoothing the spatial structure of the STA

% Section 2: Isoresponse data from example DO cell

% Section 3: White noise analysis of spatial integration. Check section 3
% of Calculating_effect_contrast.m

% Section 4: Extracting color-time signal from the example DO cell

% Section 5: Adaptive staircase in an example direction

% Section 6: Regularization of the STA and the spatial weighting function

% Section 7: Relationship between the RGB vectors of the 2 subunits in angular units 

% Section 8 - Finding a lenient criteria for a cone weights (Robustness of cone-weight classification)

close all; clearvars;
plot_counter = 1;

%% Section 1: Smoothing the spatial structure of the STA


if ~exist('plot_counter')
    plot_counter = 1;
end

% Loading all the files 
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filename = fetch(conn,'SELECT filename FROM WNthresh');
NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
close(conn);
filename = filename(strcmp(string(NTmode),"subunit"));
NTmode = NTmode(strcmp(string(NTmode),"subunit"));
spikeIdx = spikeidx_NT(strcmp(string(NTmode),"subunit"));
spikename_options = ['sig001a'; 'sig001b'];


% Include Gun noise data for statistical tests in order to estimate the RF size
channels = 3;
NPOINTS = 65536;
resize_fact2 = 1;
CHI2CRIT = 0.95; % For flaging a stixel as significant (adding gun and cone noise z-scores, squared)
maxT = 15;
crit = chi2inv(CHI2CRIT,300); % 3 color channels

% Loading the relevant file
% filenumber = 31; % red-green
filenumber = 24; % blue-yellow
WN = nex2stro(findfile(char(filename(filenumber))));


framerate = WN.sum.exptParams.framerate;
nstixperside = WN.sum.exptParams.nstixperside;
ntrials = length(WN.sum.absTrialNum);
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

% Calculating the M matrix for files
fundamentals = WN.sum.exptParams.fundamentals;
mon_spd = WN.sum.exptParams.mon_spd;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;



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
spikeidx = spikeIdx(filenumber);
spikename = spikename_options(spikeidx,:);



% Calculating STA and STC for frames which triggered spikes
out_gun = getWhtnsStats(WN,maxT,'STCOVmex',{nstixperside^2,3,maxT},spikename);
STS_gun = out_gun{1}; STCross_gun = out_gun{2}; nspikes_gun = out_gun{3}; clear out_gun;
STAs_gun = STS_gun/nspikes_gun;

% Code for Statistical testing begins here
s_gun = std(STAs_gun(:,1));
STAs_gun_z = STAs_gun./s_gun;

% Spatial map
grandz = zeros([nstixperside nstixperside]);
maxzs = [];
for i = 1:maxT
    tmp_gun = reshape(STAs_gun_z(:,i),[nstixperside nstixperside 3]); % This is the only place in the code where I use data from gun noise
    grandz = grandz+sum(tmp_gun.^2,3);
    maxzs = [maxzs; sum(sum(tmp_gun(:,:,1).^2)) sum(sum(tmp_gun(:,:,2).^2)) sum(sum(tmp_gun(:,:,3).^2))];
end
peakframe1 = max(maxzs(:,1)) == maxzs(:,1);
peakframe2 = max(maxzs(:,2)) == maxzs(:,2);
peakframe3 = max(maxzs(:,3)) == maxzs(:,3);
peakframe = max(sum(maxzs,2)) == sum(maxzs,2);
id = find(peakframe==1);
latency = find(peakframe)*1000/WN.sum.exptParams.framerate;
if id~=1
    peakframe(id-1)= 1;
end
if id <=maxT-1
    peakframe(id+1)=1;
end

STAweights = sqrt(sum(STAs_gun(:,peakframe).^2));
STAweights = STAweights./sum(STAweights);
tmpSTA = STAs_gun(:,peakframe)*STAweights';

[u,s,v] = svd(reshape(tmpSTA,[nstixperside^2 3]));

% Plotting the weighted RGB STA
normfactor = 0.5/(max(abs(tmpSTA(:)))+0.01);
im = normfactor*tmpSTA + 0.5;
im = reshape(im,[10 10 3]);

% Plotting the grayscale spatial weighting function
spatial_weighting = u(:,1);
normfactor = 0.5/(max(abs(spatial_weighting))+0.01);
im1 = normfactor*spatial_weighting + 0.5;
im1 = reshape(im1,[10 10]);

% Resizing the STA and the weighting function
resizeSTA = imresize(reshape(tmpSTA,[nstixperside nstixperside 3]),2);
normfactor = 0.5/(max(abs(resizeSTA(:)))+0.01);
resizeSTA = normfactor*resizeSTA + 0.5;

resizeweightingfunc = imresize(reshape(spatial_weighting,[nstixperside nstixperside]),2);
normfactor = 0.5/(max(abs(resizeweightingfunc(:)))+0.01);
resizeweightingfunc = normfactor*resizeweightingfunc + 0.5;

%%
% Plotting the results
figure(plot_counter);
subplot(421); image(im); axis square; set(gca,'XTick',[],'YTick',[]); title('STA')
subplot(422); imagesc(im1); axis square; set(gca,'XTick',[],'YTick',[]); title('Spatial weighting function'); colormap('gray');
subplot(423); image(resizeSTA); axis square; set(gca,'XTick',[],'YTick',[]); title('Regularized STA')
subplot(424); imagesc(resizeweightingfunc); axis square; set(gca,'XTick',[],'YTick',[]); title('Regularized weighting');
subplot(425); image(imgaussfilt(imresize(im,2),1.0)); axis square; set(gca,'XTick',[],'YTick',[]); title('Gaussian filt STA')
subplot(426); imagesc(imgaussfilt(imresize(im1,2),1.0)); axis square; set(gca,'XTick',[],'YTick',[]); title('Gaussian filt weighting')
subplot(427); image(imresize(imgaussfilt(imresize(im,2),1.0),0.5)); axis square; set(gca,'XTick',[],'YTick',[]); title('Same size STA')
subplot(428); imagesc(imresize(imgaussfilt(imresize(im1,2),1.0),0.5)); axis square; set(gca,'XTick',[],'YTick',[]); title('Same size weighting function')
plot_counter = plot_counter + 1;

%% enahancing the STA
imR = sigmoid(squeeze(im(:,:,1)),5,0.5);
imG = sigmoid(squeeze(im(:,:,2)),1,0.5);
imB = sigmoid(squeeze(im(:,:,3)),5,0.5);
im_new = cat(3,imR,imG,imB);
figure(plot_counter);
image(im_new);axis square; set(gca,'XTick',[],'YTick',[]); title('Same size STA')
plot_counter = plot_counter + 1;

%% Trying few more ways of cleaning the STA

% Just gaussian filtering the STA 
im = imgaussfilt(reshape(tmpSTA,[10 10 3]),1.0);
normfactor = 0.5/(max(abs(im(:)))+0.01);
im = normfactor*im + 0.5;

im_orig = reshape(tmpSTA,[10 10 3]);
normfactor = 0.5/(max(abs(im_orig(:)))+0.01);
im_orig = normfactor*im_orig + 0.5;

rawenergymap = sum(reshape(tmpSTA,[10 10 3]).^2,3);
energymap = imgaussfilt(rawenergymap,1.0);

figure(plot_counter);
subplot(221); image(im_orig); axis square; set(gca,'Tickdir','out','XTick',[],'YTick',[]);
subplot(222); image(im); axis square; set(gca,'Tickdir','out', 'XTick',[],'YTick',[]);
subplot(223); imagesc(rawenergymap); axis square; set(gca,'Tickdir','out', 'XTick',[],'YTick',[]);
subplot(224); imagesc(energymap); axis square; set(gca,'Tickdir','out', 'XTick',[],'YTick',[]);
plot_counter = plot_counter + 1;

tmpSTA3 = reshape(tmpSTA,[10 10 3]);
energymapR = imgaussfilt(abs(squeeze(tmpSTA3(:,:,1))),1);
energymapG = imgaussfilt(abs(squeeze(tmpSTA3(:,:,2))),2);
energymapB = imgaussfilt(abs(squeeze(tmpSTA3(:,:,3))),1);

figure(plot_counter);
subplot(231); imagesc(energymapR); axis square; set(gca,'Tickdir','out','XTick',[],'YTick',[]); 
subplot(232); imagesc(energymapG); axis square; set(gca,'Tickdir','out','XTick',[],'YTick',[]);
subplot(233); imagesc(energymapB); axis square; set(gca,'Tickdir','out','XTick',[],'YTick',[]);
% subplot(224); imagesc(energymapB+energymapG+energymapR); axis square; set(gca,'Tickdir','out','XTick',[],'YTick',[]);
subplot(234); imagesc(squeeze(tmpSTA3(:,:,1))); axis square; set(gca,'Tickdir','out','XTick',[],'YTick',[]);
subplot(235); imagesc(squeeze(tmpSTA3(:,:,2))); axis square; set(gca,'Tickdir','out','XTick',[],'YTick',[]);
subplot(236); imagesc(squeeze(tmpSTA3(:,:,3))); axis square; set(gca,'Tickdir','out','XTick',[],'YTick',[]);
plot_counter = plot_counter + 1;

%% multiplying STA with the energy map 
weighted_im = tmpSTA3.*repmat((energymap/max(energymap(:))),[1 1 3]);
normfactor = 0.5/(max(abs(weighted_im(:)))+0.01);
weighted_im = normfactor*weighted_im + 0.5;

imR = sigmoid(squeeze(weighted_im(:,:,1)),8,0.5);
imG = sigmoid(squeeze(weighted_im(:,:,2)),2,0.5);
imB = sigmoid(squeeze(weighted_im(:,:,3)),8,0.5);
im_new = cat(3,imR,imG,imB);

% Performing the same thing for the spatial weighting function
weighted_sim = reshape(spatial_weighting,[10 10]).*(energymap/max(energymap(:)));
normfactor = 0.5/(max(abs(weighted_sim(:)))+0.01);
weighted_sim = normfactor*weighted_sim + 0.5;
sim = sigmoid(weighted_sim,2,0.5);

% Trying out a different strategy 
eR = abs(squeeze(tmpSTA3(:,:,1).^2)); eR = eR/max(eR(:));
eG = abs(squeeze(tmpSTA3(:,:,2).^2)); eG = eG/max(eG(:));
eB = abs(squeeze(tmpSTA3(:,:,3).^2)); eB = eB/max(eB(:));
STA_orig = tmpSTA3.*cat(3,eR,eG,eB);
normfactor = 0.5/(max(abs(STA_orig(:)))+0.01);
STA_orig = normfactor*STA_orig + 0.5;

imR_orig = sigmoid(squeeze(STA_orig(:,:,1)),10,0.5);
imG_orig = sigmoid(squeeze(STA_orig(:,:,2)),5,0.5);
imB_orig = sigmoid(squeeze(STA_orig(:,:,3)),10,0.5);
im_new_orig = cat(3,imR_orig,imG_orig,imB_orig);

figure(plot_counter);
subplot(221); image(im_new); axis square; set(gca,'XTick',[],'YTick',[]); title('Same size STA')
subplot(222); imagesc(weighted_sim); axis square; set(gca,'XTick',[],'YTick',[]); title('Same size spatial weighting function'); colormap('gray');
subplot(223); image(im_new_orig); axis square; set(gca,'XTick',[],'YTick',[]); title('Same size STA')
plot_counter = plot_counter + 1;



%% Section 2: Isoresponse data from example DO cell 

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

% Classifying cells 
LUMidx = LumIds_conewts;
DOidx = [ColorOpponentIds_conewts Sconedominated_conewts];
hardtoclassifyidx = [Other_conewts];
hardtoclassifyidx = [hardtoclassifyidx LUMidx(vals(LUMidx)>=95) DOidx(vals(DOidx)>=95)];
LUMidx = LUMidx(vals(LUMidx)<95);
DOidx = DOidx(vals(DOidx)<95);

load RSSE_linearmodel_CV.mat % Robust regression
load RSSE_quadmodel_CV.mat

% Calculating the median of differences/ratios
RSSEisoresp_medianofratios = [];
for ii = 1:numel(RSSE_linearmodel)   
    % computation for calculating median of differences/ratios
    RSSEisoresp_medianofratios = [RSSEisoresp_medianofratios; median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii})];   
end

% Loading all the files 
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filename = fetch(conn,'SELECT filename FROM WNthresh');
NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
close(conn);
filename = filename(strcmp(string(NTmode),"subunit"));
NTmode = NTmode(strcmp(string(NTmode),"subunit"));
spikeidx_NT = spikeidx_NT(strcmp(string(NTmode),"subunit"));

% These are files that I am concerned with
load RHO_all.mat
load THETA_all.mat
load not_oog_idx_all.mat
load oog_idx_all.mat
load linear_modelparams.mat
load quad_modelparams.mat
load subunitbasisvec.mat
indices = [24];% example DO 
plot_counter = 3;
for zz = 1:numel(indices)
    stro = nex2stro(findfile(char(filename(indices(zz)))));
    global spikename maskidx spikeidx neurothreshidx nstixperside ngammasteps seedidx nframesidx correctidx
    global fponidx stimoffidx stimonidx muidxs sigmaidxs basisvecidx weightsidx fpacqidx targetspikerateidx basisvecdiridx latencyidx reversalflagidx parentverticesidx
    global msperframe ntrials maxT xx yy M linepredtol stepsizescale stepsize nreversals oogscale
    spikename = 'sig001a';
    maskidx = strcmp(stro.sum.rasterCells(1,:),'subunit_mask');
    spikeidx = strcmp(stro.sum.rasterCells(1,:),spikename);
    basisvecidx = strcmp(stro.sum.rasterCells(1,:),'basis_vec');
    weightsidx = strcmp(stro.sum.rasterCells(1,:),'weights');
    parentverticesidx = strcmp(stro.sum.rasterCells(1,:),'parentvertices');
    nstixperside = stro.sum.exptParams.nstixperside;
    ngammasteps = 2^16; % 65536
    linepredtol = stro.sum.exptParams.linepredtol;
    stepsizescale = stro.sum.exptParams.stepsizescale;
    stepsize = stro.sum.exptParams.stepsize;
    nreversals = stro.sum.exptParams.nreversals;
    oogscale = stro.sum.exptParams.oogscale;
    seedidx = strcmp(stro.sum.trialFields(1,:),'seed');
    nframesidx = strcmp(stro.sum.trialFields(1,:),'num_frames');
    stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on');
    stimoffidx = strcmp(stro.sum.trialFields(1,:),'stim_off');
    fponidx = strcmp(stro.sum.trialFields(1,:),'fp_on');
    fpacqidx = strcmp(stro.sum.trialFields(1,:),'fpacq');
    basisvecdiridx = strcmp(stro.sum.trialFields(1,:),'weights_idx');
    neurothreshidx = strcmp(stro.sum.trialFields(1,:),'neurothresh'); % when exactly the neurothresh trials started
    targetspikerateidx = strcmp(stro.sum.trialFields(1,:),'targetspikerate');
    correctidx = strcmp(stro.sum.trialFields(1,:),'correct');
    muidxs = [find(strcmp(stro.sum.trialFields(1,:),'mu1')), ...
        find(strcmp(stro.sum.trialFields(1,:),'mu2')), ...
        find(strcmp(stro.sum.trialFields(1,:),'mu3'))];
    sigmaidxs = [find(strcmp(stro.sum.trialFields(1,:),'sigma1')), ...
        find(strcmp(stro.sum.trialFields(1,:),'sigma2')), ...
        find(strcmp(stro.sum.trialFields(1,:),'sigma3'))];
    latencyidx = strcmp(stro.sum.trialFields(1,:),'latency');
    reversalflagidx = strcmp(stro.sum.trialFields(1,:),'reversalflag');
    msperframe = 1000/stro.sum.exptParams.framerate;
    ntrials = size(stro.trial,1);
    maxT = 10; % this represents the temporal part in the spatiotemporal receptive field
    xx = linspace(stro.sum.exptParams.gauss_locut/1000, stro.sum.exptParams.gauss_hicut/1000,ngammasteps); % xx represents the probabilities. For more info, have a look at the MATLAB 'norminv' function.
    yy = norminv(xx'); % defining norminv to extract the values for which the cdf values range between gauss_locut and gauss_hicut
    t_offset = stro.trial(end,latencyidx)/1000;
    
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
    
    if isfield(stro.sum.exptParams,'nrepframes')
        if ~isnan(stro.sum.exptParams.nrepframes)
            nvblsperstimupdate = stro.sum.exptParams.nrepframes;
        else
            nvblsperstimupdate = 1;
        end
    else
        nvblsperstimupdate = 1;
    end
    
    mask_changes = reshape(mask_changes , 2, []);
    for ii = 1:2
        idxs = zeros(size(stro.trial,1),1);
        idxs(mask_changes(1,ii):mask_changes(2,ii)) = 1;
        idxs = logical(idxs);
        WN = stro;
        WN.ras(~idxs,:) = []; WN.trial(~idxs,:) = [];
        if ii == 1
            nrandnumsperchannel = nstixperside^2;
        else
            mask = stro.ras{mask_changes(1,ii),maskidx}; % subunit mask
            mask(mask == 0) = Inf;
            [stIdxs,~,~] = unique(mask); % now the Infs map to nsubunits+1
            num_subunits = length(stIdxs)-any(isinf(stIdxs)); % nsubunits, like subunits A and B
            mask(isinf(mask)) = num_subunits + 1;
            nrandnumsperchannel = num_subunits;
        end
        out_gun = getWhtnsStats(WN,maxT,'STCOVmex',{nrandnumsperchannel,3,maxT},spikename);
        STS_gun = out_gun{1}; STCross_gun = out_gun{2}; nspikes_gun = out_gun{3}; clear out_gun;
        STAs_gun = STS_gun/nspikes_gun;
        energy = sum(STAs_gun.^2,1);
        peakframe = energy == max(energy);
        id = find(peakframe==1);
        latency = find(peakframe)*1000/stro.sum.exptParams.framerate;
        if id~=1
            peakframe(id-1)= 1;
        end
        if id <=maxT-1
            peakframe(id+1)=1;
        end
        STAweights = sqrt(sum(STAs_gun(:,peakframe).^2));
        STAweights = STAweights./sum(STAweights);
        tmpSTA = STAs_gun(:,peakframe)*STAweights'; % weighted combination of peak and its adjacent frames
        tmpSTA2 = STAs_gun(:,id); % just the peak frame
        if ii == 1
            [u1,~,v1] = svd(reshape(tmpSTA,[nstixperside^2 3])');
            SpatialRF = reshape(v1(:,1),[nstixperside nstixperside]);
        elseif ii == 2
            G_mask = [mask; mask+max(mask); mask+2*max(mask)];
            tmpSTA = expand_vector(tmpSTA,num_subunits,G_mask,1);
            tmpSTA2 = expand_vector(tmpSTA2,num_subunits,G_mask,1);
        end
        normfactor = 0.5/(max(abs(tmpSTA(:)))+0.01);
        tmpSTA = normfactor*tmpSTA + 0.5;
        tmpSTA = reshape(tmpSTA,[nstixperside nstixperside 3]);
        normfactor = 0.5/(max(abs(tmpSTA2(:)))+0.01);
        tmpSTA2 = normfactor*tmpSTA2 + 0.5;
        tmpSTA2 = reshape(tmpSTA2,[nstixperside nstixperside 3]);
        figure(1),subplot(numel(indices),4,4*(zz-1)+ii); image(tmpSTA); set(gca,'XTick',[],'YTick',[]); axis square;
    end
    
    % Determining when Neurothresh mode was active, plotting the basis vector, working correctly
    neurothreshmode = stro.trial(:,neurothreshidx);
    basisvec_dropidx = inds(end);
    neurothresh_startidx = max(find(neurothreshmode==1,1),basisvec_dropidx+1); % There is possibility of pressing the Neurothresh button multiple times
    num_targetspikerates = unique(stro.trial(neurothresh_startidx:end,targetspikerateidx));
    vect = stro.ras{basisvec_dropidx,basisvecidx};
    basisvec_size = nstixperside*nstixperside*3;
    numvect = (numel(vect)/basisvec_size)-1;
    basisvec = cell(1,numvect);
    % Actual basis vec
    for ii = 1:numvect
        tmp_vec = vect((ii-1)*basisvec_size+1:basisvec_size*ii) + vect(numvect*basisvec_size+1:basisvec_size*(numvect+1));
        basisvec{ii} = reshape(tmp_vec,[nstixperside nstixperside 3]);
    end
    bkgnd_monitor = reshape(vect(numvect*basisvec_size+1:basisvec_size*(numvect+1)),[nstixperside nstixperside 3]);
    tmp1 = unique(basisvec{1}-bkgnd_monitor,'stable');
    tmp2 = unique(basisvec{2}-bkgnd_monitor,'stable');
   
    figure(1); subplot(numel(indices),4,4*zz-1); image(imresize(255*(SpatialRF./(2*max(abs(SpatialRF(:))))+.5),2)); set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255))
    figure(1); subplot(numel(indices),4,4*zz); image(subunitbasisvec{indices(zz)}); set(gca,'XTick',[],'YTick',[]);axis square;
    plot_counter = plot_counter + 1;
    
    ind = indices(zz);
    THETA = THETA_all{1,ind};
    THETA = THETA * pi/180; % converting to radians
    if any(THETA>(135*pi/180))
        allthetas = linspace(-pi,pi,100);
        newtheta = linspace(-pi,pi,101);
    else
        allthetas = linspace(-pi/4,3*pi/4,100);
        newtheta = linspace(-pi/4,3*pi/4,101);
    end
    RHO = RHO_all{1,ind};
    oog_idx = oog_idx_all{1,ind};
    not_oog_idx = not_oog_idx_all{1,ind};
    outofgamut = zeros(size(THETA));
    outofgamut(oog_idx) = 1;
    outofgamut = logical(outofgamut);
    [x_orig, y_orig] = pol2cart(THETA,RHO);

    % Linear model predictions
    rho1 = 1./(linear_modelparams(ind,:)*[cos(allthetas); sin(allthetas)]);
    LOOGtmp1= rho1<0;
    [x_lin,y_lin] = pol2cart(allthetas(~LOOGtmp1),rho1(~LOOGtmp1));
    
    % Quadratic model predictions
    [x_quad,y_quad,rho3] = calc_xyvalues(allthetas, quad_modelparams(ind,:));
    L = rho3>0 & rho3==real(rho3);
    [x_quad2,y_quad2] = pol2cart(newtheta(L),rho3(L)');
    
    % Plotting the isoresponse data
    figure(2); subplot(1,numel(indices),zz);
    for ii = 1:numel(not_oog_idx)
        h(ii) = plot(x_orig(not_oog_idx(ii)), y_orig(not_oog_idx(ii)),'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
        set(h(ii),'ButtonDownFcn',{@dispimage,x_orig(not_oog_idx(ii)),y_orig(not_oog_idx(ii)),basisvec,bkgnd_monitor});%(x_orig(not_oog_idx(ii)),y_orig(not_oog_idx(ii)),basisvec,bkgnd_monitor)});
    end
    plot(x_lin,y_lin,'g','Linewidth',2); hold on; plot(x_quad2,y_quad2,'r','Linewidth',2);
    xlabel('X'), ylabel('Y');
    
    
    if any(outofgamut)
        plot(upsample(x_orig(outofgamut),2),upsample(y_orig(outofgamut),2),'color',[0.5 0.5 0.5]);
    end
    set(gca,'Tickdir','out'); drawnow; axis square; 
    
    set(gca,'Xlim',[-0.8 0.8],'Ylim',[-0.8 0.8],'YTick',-0.8:0.4:0.8,'XTick',[-0.8:0.4:0.8]); hold off;
    
    
    % Now, I am going to plot the
    wts  = [1 0;1 1;0 1;-1 1;-1 0;-1 -1;0 -1;1 -1];
    subplotidxs = [1;2;3;4;6;7;8;9];
    for ii = 1:numel(subplotidxs)
        vec = wts(ii,1)*(basisvec{1}-bkgnd_monitor) + wts(ii,2)*(basisvec{2}-bkgnd_monitor);
        normfactor = 0.5/(max(abs(vec(:)))+0.01);
        vec = normfactor*vec + 0.5;
        figure(3),subplot(numel(indices),size(wts,1),size(wts,1)*(zz-1)+ii); image(vec); set(gca,'XTick',[],'YTick',[]); axis square;
    end
    
    % Plotting the data, linear and non-linear fits in the R-Theta space
    figure(4); subplot(1,numel(indices),zz);
    plot(allthetas(~LOOGtmp1),log10(rho1(~LOOGtmp1)),'g','Linewidth',2); hold on; 
    plot(newtheta(L),log10(rho3(L)),'r','LineWidth',2); hold on; 
    plot(THETA(not_oog_idx), log10(RHO(not_oog_idx)),'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'PickableParts','none','MarkerEdgeColor',[1 1 1]);
    plot(THETA(oog_idx), log10(RHO(oog_idx)),'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'PickableParts','none','MarkerEdgeColor',[1 1 1]); 
    xlabel('Theta'); ylabel('log10(R)');
    
    axis square; hold off;
   
    
    plot_counter = plot_counter + 1;
end

% Obtaining the target firing rates of the example neurons
load TFR.mat
FR = TFR(1,indices);

%% Section 3: White noise analysis of spatial integration. Check section 3
% of Calculating_effect_contrast.m


%% Section 4: Extracting color-time signal from the example DO cell


if ~exist('plot_counter')
    plot_counter = 1;
end

% Loading all the files 
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filename = fetch(conn,'SELECT filename FROM WNthresh');
NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
close(conn);
filename = filename(strcmp(string(NTmode),"subunit"));
NTmode = NTmode(strcmp(string(NTmode),"subunit"));
spikeidx_NT = spikeidx_NT(strcmp(string(NTmode),"subunit"));


% Include Gun noise data for statistical tests in order to estimate the RF size
channels = 3;
NPOINTS = 65536;
resize_fact2 = 1;
CHI2CRIT = 0.95; % For flaging a stixel as significant (adding gun and cone noise z-scores, squared)
maxT = 15;
crit = chi2inv(CHI2CRIT,300); % 3 color channels


WN = nex2stro(findfile(char(filename(24))));
framerate = WN.sum.exptParams.framerate;
nstixperside = WN.sum.exptParams.nstixperside;
ntrials = length(WN.sum.absTrialNum);
stimonidx = find(strcmp(WN.sum.trialFields(1,:),'stim_on'));
stimoffidx = find(strcmp(WN.sum.trialFields(1,:),'all_off'));
nframesidx = find(strcmp(WN.sum.trialFields(1,:),'num_frames'));
noisetypeidx = find(strcmp(WN.sum.trialFields(1,:),'noise_type'));
sigmaidxs = strmatch('sigma',WN.sum.trialFields(1,:));
hepidx = find(strcmp(WN.sum.rasterCells(1,:),'AD11'));
vepidx = find(strcmp(WN.sum.rasterCells(1,:),'AD12'));
maskidx = strcmp(WN.sum.rasterCells(1,:), 'subunit_mask');
anlgStartTimeidx = find(strcmp(WN.sum.rasterCells(1,:),'anlgStartTime'));
mu1idx = find(strcmp(WN.sum.trialFields(1,:),'mu1'));
mu2idx = find(strcmp(WN.sum.trialFields(1,:),'mu2'));
mu3idx = find(strcmp(WN.sum.trialFields(1,:),'mu3'));
sigma1idx = find(strcmp(WN.sum.trialFields(1,:),'sigma1'));
sigma2idx = find(strcmp(WN.sum.trialFields(1,:),'sigma2'));
sigma3idx = find(strcmp(WN.sum.trialFields(1,:),'sigma3'));
maskidx = strcmp(WN.sum.rasterCells(1,:),'subunit_mask');
basisvecidx = strcmp(WN.sum.rasterCells(1,:),'basis_vec');
latencyidx = strcmp(WN.sum.trialFields(1,:),'latency');
neurothreshidx = strcmp(WN.sum.trialFields(1,:),'neurothresh');

L = WN.trial(:,noisetypeidx)==1;
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

% Getting the background rgb/lms
% Calculating the M matrix
fundamentals = WN.sum.exptParams.fundamentals;
mon_spd = WN.sum.exptParams.mon_spd;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

ridx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_r'));
gidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_g'));
bidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_b'));
bkgndRGB = [mode(WN.trial(:,ridx)), mode(WN.trial(:,gidx)), mode(WN.trial(:,bidx))];
bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
bkgndlms = M*bkgndrgb;
Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
Mrgbtocc = inv(Mrgbtocc');

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
   
    
    % Modifying the STRO/WN structure after I have pulled out the
    % basis vector
    idxs = zeros(size(WN.trial,1),1);
    idxs(mask_changes(2,1):mask_changes(2,2)) = 1;
    idxs = logical(idxs);
    WN.ras(~idxs,:) = []; % modiftying the WN structure
    WN.trial(~idxs,:) = []; % modiftying the WN structure
end
spikename = 'sig001a';

% Calculating STA and STC for frames which triggered spikes
out_gun = getWhtnsStats_AD(WN,maxT,'STCOVmex',{2,3,maxT},2,spikename);
STS_gun = out_gun{1}; STCross_gun = out_gun{2}; nspikes_gun = out_gun{3}; clear out_gun;
STAs_gun = STS_gun/nspikes_gun;

% PLotting the color signals for individual subuinits
figure(plot_counter); set(gcf,'Name','Subunit color-time signal');
subplot(121), plot(STAs_gun(1,:),'r','Linewidth',2); hold on; plot(STAs_gun(3,:),'g','Linewidth',2); plot(STAs_gun(5,:),'b','Linewidth',2);
xlabel('time'); ylabel('Amplitude'); axis square; set(gca,'Tickdir','out'); hold off;
subplot(122), plot(STAs_gun(2,:),'r','Linewidth',2); hold on; plot(STAs_gun(4,:),'g','Linewidth',2); plot(STAs_gun(6,:),'b','Linewidth',2);
xlabel('time'); ylabel('Amplitude'); axis square; set(gca,'Tickdir','out'); hold off;
plot_counter = plot_counter + 1;

%% % Section 5: Adaptive staircase in an example direction

if ~exist('plot_counter')
    plot_counter = 1;
end

% Loading all the files 
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filename = fetch(conn,'SELECT filename FROM WNthresh');
NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
close(conn);
filename = filename(strcmp(string(NTmode),"subunit"));
NTmode = NTmode(strcmp(string(NTmode),"subunit"));
spikeidx_NT = spikeidx_NT(strcmp(string(NTmode),"subunit"));

stro = nex2stro(findfile(char(filename(24))));
global spikename maskidx spikeidx neurothreshidx nstixperside ngammasteps seedidx nframesidx correctidx
global fponidx stimoffidx stimonidx muidxs sigmaidxs basisvecidx weightsidx fpacqidx targetspikerateidx basisvecdiridx latencyidx reversalflagidx parentverticesidx
global msperframe ntrials maxT xx yy M linepredtol stepsizescale stepsize nreversals oogscale

spikename = getSpikenum(stro);
maskidx = strcmp(stro.sum.rasterCells(1,:),'subunit_mask');
spikeidx = strcmp(stro.sum.rasterCells(1,:),spikename);
basisvecidx = strcmp(stro.sum.rasterCells(1,:),'basis_vec');
weightsidx = strcmp(stro.sum.rasterCells(1,:),'weights');
parentverticesidx = strcmp(stro.sum.rasterCells(1,:),'parentvertices');
nstixperside = stro.sum.exptParams.nstixperside;
ngammasteps = 2^16; % 65536
linepredtol = stro.sum.exptParams.linepredtol;
stepsizescale = stro.sum.exptParams.stepsizescale;
stepsize = stro.sum.exptParams.stepsize;
nreversals = stro.sum.exptParams.nreversals;
oogscale = stro.sum.exptParams.oogscale;
seedidx = strcmp(stro.sum.trialFields(1,:),'seed');
nframesidx = strcmp(stro.sum.trialFields(1,:),'num_frames');
stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on');
stimoffidx = strcmp(stro.sum.trialFields(1,:),'stim_off');
fponidx = strcmp(stro.sum.trialFields(1,:),'fp_on');
fpacqidx = strcmp(stro.sum.trialFields(1,:),'fpacq');
basisvecdiridx = strcmp(stro.sum.trialFields(1,:),'weights_idx');
neurothreshidx = strcmp(stro.sum.trialFields(1,:),'neurothresh'); % when exactly the neurothresh trials started
targetspikerateidx = strcmp(stro.sum.trialFields(1,:),'targetspikerate');
correctidx = strcmp(stro.sum.trialFields(1,:),'correct');
muidxs = [find(strcmp(stro.sum.trialFields(1,:),'mu1')), ...
    find(strcmp(stro.sum.trialFields(1,:),'mu2')), ...
    find(strcmp(stro.sum.trialFields(1,:),'mu3'))];
sigmaidxs = [find(strcmp(stro.sum.trialFields(1,:),'sigma1')), ...
    find(strcmp(stro.sum.trialFields(1,:),'sigma2')), ...
    find(strcmp(stro.sum.trialFields(1,:),'sigma3'))];
latencyidx = strcmp(stro.sum.trialFields(1,:),'latency');
reversalflagidx = strcmp(stro.sum.trialFields(1,:),'reversalflag');
msperframe = 1000/stro.sum.exptParams.framerate;
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

% Determining when Neurothresh mode was active, plotting the basis vector, working correctly
t_offset = stro.trial(end,latencyidx)/1000;
neurothreshmode = stro.trial(:,neurothreshidx);
basisvec_dropidx = inds(end);
neurothresh_startidx = max(find(neurothreshmode==1,1),basisvec_dropidx+1); % There is possibility of pressing the Neurothresh button multiple times
num_targetspikerates = unique(stro.trial(neurothresh_startidx:end,targetspikerateidx));
vect = stro.ras{basisvec_dropidx,basisvecidx};
basisvec_size = nstixperside*nstixperside*3;
numvect = (numel(vect)/basisvec_size)-1;
basisvec = cell(1,numvect);
figure(plot_counter);
% Actual basis vec
bkgnd_monitor = reshape(vect(numvect*basisvec_size+1:basisvec_size*(numvect+1)),[nstixperside nstixperside 3]);
for ii = 1:numvect
    tmp_vec = vect((ii-1)*basisvec_size+1:basisvec_size*ii) + vect(numvect*basisvec_size+1:basisvec_size*(numvect+1));
    basisvec{ii} = reshape(tmp_vec,[nstixperside nstixperside 3]);  
end

% This section works correctly
% plotting the weights corresponding to the weight directions - useful for
% checking if the direction indexes and the directions are aligned or not
weight_direction_mat = [];
for i=neurothresh_startidx:size(stro.trial,1)
    weight_direction_mat = [weight_direction_mat; stro.trial(i,basisvecdiridx) stro.ras{i,weightsidx}'/norm(stro.ras{i,weightsidx}')];
end
[~,idx] = sort(weight_direction_mat(:,1)); % sorting it according to the increasing order of weight direction indexes
weight_direction_mat1 = weight_direction_mat(idx,:);

% plotting the rasterplots for neurothresh trials
norms = cell(1,numel(num_targetspikerates));
completed_search_alongdir = cell(1,numel(num_targetspikerates));
for jj = 1: numel(num_targetspikerates)
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
        for ii = 1:size(raster_data,1)
            tmp = raster_data{ii} ;
            spikes = tmp(tmp>stro.trial(idxs1(ii),stimonidx)+t_offset & tmp < stro.trial(idxs1(ii),stimoffidx));
            spikes = spikes - stro.trial(idxs1(ii),stimonidx)-t_offset;
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

%%
figure(plot_counter);
cum_flag = []; % to check if I am looking at all the probed directions

t_offset1 = -0.2;
figure(plot_counter);

% Good directions: 14, 24, 25, 26
dir = 14;
for jj = 1: numel(num_targetspikerates)
    tmp_n = [];
    tmp_wts = [];
    tmp_parentvertices = [];
    idxs1 = find(stro.trial(:,basisvecdiridx) == dir);
    
    for jj = 1:numel(idxs1)
        tmp_parentvertices = [tmp_parentvertices; stro.ras{idxs1(jj),parentverticesidx}];
        tmp_wts = [tmp_wts; stro.ras{idxs1(jj),weightsidx}'];
        tmp_n = [tmp_n; norm(stro.ras{idxs1(jj),weightsidx})];
    end
    [flag, gamutViolation] = Check_ifcompletely_probed(stro,basisvec,bkgnd_monitor,idxs1,tmp_wts(end,:));
    cum_flag = [cum_flag; flag];
    if gamutViolation == 1
        c = [0 1 0];
    else
        c = [0 0 1];
    end
    subplot(2,2,1),plot(tmp_n,'-o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); xlabel('Trials'), ylabel('Contrast'); 
    set(gca,'Tickdir','out','YScale','log','Ylim',[0.03 1],'YTick',[0.03 0.1 0.3 1]); axis square;
    img = tmp_wts(end,1)*(basisvec{1}-bkgnd_monitor) + tmp_wts(end,2)*(basisvec{2}-bkgnd_monitor);
    subplot(2,2,2),imagesc(0.5*img/(max(abs(img(:)))+0.01) + 0.5); set(gca,'XTick',[],'YTick',[]); axis image; colormap('gray');
    raster_data = stro.ras(idxs1,1);
    num_dur =[];
    firing_rate = [];
    
    for ii = 1:size(raster_data,1)
        tmp = raster_data{ii} ;
        spikes = tmp(tmp<stro.trial(idxs1(ii),stimoffidx) & tmp>stro.trial(idxs1(ii),stimonidx)+t_offset);
        num_dur = [num_dur; (stro.trial(idxs1(ii), stimoffidx)- stro.trial(idxs1(ii),stimonidx))-t_offset];
        firing_rate = [firing_rate; numel(spikes)/num_dur(end)];
    end
    subplot(2,2,4); plot(tmp_n,firing_rate,'-o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
    xlabel('Contrast'); ylabel('Firing rate '); axis square;
    set(gca,'Tickdir','out','XScale','log')
    hold off;
    subplot(2,2,3),plot(firing_rate,'-o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
    set(gca,'Tickdir','out','Ylim',[0 250],'YTick',[0 50 100 150 200 250]); axis square; xlabel('Trials'), ylabel('FR');
end

if dir == 3
    subplot(2,2,4); hold on; set(gca,'Xlim',[0.03 1],'XTick',[0.03 0.1 0.3 1],'Ylim',[0 250],'YTick',[0 50 100 150 200 250]);
    subplot(2,2,3); hold on; plot([0 20],[20 20],'k'); hold off;
end
plot_counter = plot_counter + 1;

%% Section 6: Regularization of the STA and the spatial weighting function

if ~exist('plot_counter')
    plot_counter = 1;
end

% Loading all the files 
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filename = fetch(conn,'SELECT filename FROM WNthresh');
NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
close(conn);
filename = filename(strcmp(string(NTmode),"subunit"));
NTmode = NTmode(strcmp(string(NTmode),"subunit"));
spikeIdx = spikeidx_NT(strcmp(string(NTmode),"subunit"));
spikename_options = ['sig001a'; 'sig001b'];


% Include Gun noise data for statistical tests in order to estimate the RF size
channels = 3;
NPOINTS = 65536;
resize_fact2 = 1;
CHI2CRIT = 0.95; % For flaging a stixel as significant (adding gun and cone noise z-scores, squared)
maxT = 15;
crit = chi2inv(CHI2CRIT,300); % 3 color channels

% Loading the relevant file
% filenumber = 31; % red-green
filenumber = 24; % blue-yellow
WN = nex2stro(findfile(char(filename(filenumber))));

msperframe = 1000/WN.sum.exptParams.framerate;
framerate = WN.sum.exptParams.framerate;
nstixperside = WN.sum.exptParams.nstixperside;
ngammasteps = 2^16; % 65536
ntrials = length(WN.sum.absTrialNum);
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
muidxs = [find(strcmp(WN.sum.trialFields(1,:),'mu1')), ...
        find(strcmp(WN.sum.trialFields(1,:),'mu2')), ...
        find(strcmp(WN.sum.trialFields(1,:),'mu3'))];
    sigmaidxs = [find(strcmp(WN.sum.trialFields(1,:),'sigma1')), ...
        find(strcmp(WN.sum.trialFields(1,:),'sigma2')), ...
        find(strcmp(WN.sum.trialFields(1,:),'sigma3'))];
maskidx = strcmp(WN.sum.rasterCells(1,:),'subunit_mask');
seedidx = strcmp(WN.sum.trialFields(1,:),'seed');
basisvecidx = strcmp(WN.sum.rasterCells(1,:),'basis_vec');
gammaTable = WN.sum.exptParams.gamma_table;
gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
gammaTable1 = interp1(linspace(0,255,256),gammaTable,linspace(0,255,65536), 'spline');
invgamma = InvertGamma(gammaTable, 0);

gausslims = [WN.sum.exptParams.gauss_locut WN.sum.exptParams.gauss_hicut]/1000;
x = linspace(gausslims(1),gausslims(2),NPOINTS);
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
spikeidx = spikeIdx(filenumber);
spikename = spikename_options(spikeidx,:);

% Storing the stimuli and the corresponding spikes

cum_rgbs = []; % a cell for storing all the rgbs frames
cum_n = [];

for k = 1:size(WN.trial,1)
    nframes = WN.trial(k,nframesidx);
    if (nframes == 0)
        continue;
    end
    seed = WN.trial(k,seedidx);
    mu = WN.trial(k,muidxs)/1000;
    sigma = WN.trial(k,sigmaidxs)/1000;
    
    % org_mask tells u if u have updated the mask or not. If org_mask is non-zero it means at this particular trial
    % u have selected the subunits and need to analyse its computation
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
    % Deleting the spikes taking place before the first 9 frames as I need to look at the 9 preceding frames
    spiketimes(spiketimes < maxT*msperframe) = [];
    % Deleting the spikes that take place after the stimulus was
    % removed as it would imply that these spikes do not result from
    % the stimulus shown on the screen
    spiketimes(spiketimes > frametimes(end)) = [];
    n = hist(spiketimes, frametimes);
    
    cum_rgbs = [cum_rgbs {rgbs}];
    cum_n = [cum_n {n}];
    
end

% Calculating the STA the standard way 
spikeidx = spikeIdx(filenumber);
spikename = spikename_options(spikeidx,:);
out_gun = getWhtnsStats(WN,maxT,'STCOVmex',{nstixperside^2,3,maxT},spikename);
STS_gun = out_gun{1}; nspikes_gun = out_gun{3}; clear out_gun;
STAs_gun = STS_gun/nspikes_gun;

% Calculating the weighted STA 
s_gun = std(STAs_gun(:,1));
STAs_gun_z = STAs_gun./s_gun;

% Spatial map
grandz = zeros([nstixperside nstixperside]);
maxzs = [];
for i = 1:maxT
    tmp_gun = reshape(STAs_gun_z(:,i),[nstixperside nstixperside 3]); % This is the only place in the code where I use data from gun noise
    grandz = grandz+sum(tmp_gun.^2,3);
    maxzs = [maxzs; sum(sum(tmp_gun(:,:,1).^2)) sum(sum(tmp_gun(:,:,2).^2)) sum(sum(tmp_gun(:,:,3).^2))];
end
peakframe1 = max(maxzs(:,1)) == maxzs(:,1);
peakframe2 = max(maxzs(:,2)) == maxzs(:,2);
peakframe3 = max(maxzs(:,3)) == maxzs(:,3);
peakframe = max(sum(maxzs,2)) == sum(maxzs,2);
id = find(peakframe==1);
latency = find(peakframe)*1000/WN.sum.exptParams.framerate;
if id~=1
    peakframe(id-1)= 1;
end
if id <=maxT-1
    peakframe(id+1)=1;
end

STAweights = sqrt(sum(STAs_gun(:,peakframe).^2));
STAweights = STAweights./sum(STAweights);
tmpSTA = STAs_gun(:,peakframe)*STAweights';

% Plotting the weighted RGB STA
normfactor = 0.5/(max(abs(tmpSTA(:)))+0.01);
im = normfactor*tmpSTA + 0.5;
im = reshape(im,[10 10 3]);

% Creating the cosine time filter
nkt = maxT-1; % number of ms in stim filter
neye = 0; % number of "identity" basis vectors near time of spike;
ncos = 5; % number of raised-cosine vectors to use
kpeaks = [.1 round(nkt/1.2)];  % position of first and last bump (relative to identity bumps)
b = 10;
kdt = 1;  % spacing of x axis must be in units of 1

nlin = @(x)log(x+1e-20);
invnl = @(x)exp(x)-1e-20; % inverse nonlinearity

yrnge = nlin(kpeaks+b);  
db = diff(yrnge)/(ncos-1);  % spacing between raised cosine peaks
ctrs = yrnge(1):db:yrnge(2);  % centers for basis vectors
mxt = nkt; %invnl(yrnge(2)+2*db)-b; % maximum time bin
kt0 = [0:kdt:mxt]';
nt = length(kt0);        % number of points in iht
ff = @(x,c,dc)(cos(max(-pi,min(pi,(x-c)*pi/dc/2)))+1)/2; % raised cosine basis vector
kbasis0 = ff(repmat(nlin(kt0+b), 1, ncos), repmat(ctrs, nt, 1), db);

% Concatenate identity-vectors
nkt0 = size(kt0,1);
kbasis = [[eye(neye); zeros(nkt0,neye)] [zeros(neye, ncos); kbasis0]];
% kbasis = flipud(kbasis);  % flip so fine timescales are at the end.
nkt0 = size(kbasis,1);
kbasis = normalizecols(kbasis);

figure(plot_counter);
subplot(121); image(im); axis square; set(gca,'XTick',[], 'YTick',[]); 
subplot(122); plot(kbasis,'Linewidth',2); axis square; set(gca,'Tickdir','out');
plot_counter = plot_counter + 1;

% Passing the stimulus and the spiking sequences to the regularizations function

% Have written more extensive code in the MatlabCode/Abhishek/Check_regularize_STA
V = regularizeSTA(cum_rgbs, cum_n, maxT, tmpSTA, kbasis, flip(STAs_gun,2));

%Don't know why this is not working at all. Am leaving this for the
%time being but will come back to some other day.

%% Section 7: Relationship between the RGB vectors of the 2 subunits in angular units 

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

load S1RGB.mat
load S2RGB.mat

unitvectors = [];

for ii = 1:size(S1RGB,2)
    vec1 = S1RGB(:,ii);
    vec2 = S2RGB(:,ii);
    vec = exp(atan2(norm(cross(vec1,vec2)),dot(vec1,vec2))*i);
    
    % Converting the angle difference into unit vectors
    unitvectors = [unitvectors; real(vec) imag(vec)];
end

% Summary circular stats for all cells
meanvec = mean(unitvectors,1);
stdvec = std(unitvectors,1);
% Calculating the circular mean and circular standard deviation
circmean = atan2(meanvec(2),meanvec(1))*180/pi;
circstd = atan2(stdvec(2),stdvec(1))*180/pi;


% Summary circular stats for Simple cells
meanvec_LUM = mean(unitvectors(LUMidx,:),1);
stdvec_LUM = std(unitvectors(LUMidx,:),1);
% Calculating the circular mean and circular standard deviation
circmean_LUM = atan2(meanvec_LUM(2),meanvec_LUM(1))*180/pi;
circstd_LUM = atan2(stdvec_LUM(2),stdvec_LUM(1))*180/pi;


% Summary circular stats for DO cells
meanvec_DO = mean(unitvectors(DOidx,:),1);
stdvec_DO = std(unitvectors(DOidx,:),1);
% Calculating the circular mean and circular standard deviation
circmean_DO = atan2(meanvec_DO(2),meanvec_DO(1))*180/pi;
circstd_DO = atan2(stdvec_DO(2),stdvec_DO(1))*180/pi;

% Summary circular stats for HTC cells
meanvec_HTC = mean(unitvectors(hardtoclassifyidx,:),1);
stdvec_HTC = std(unitvectors(hardtoclassifyidx,:),1);
% Calculating the circular mean and circular standard deviation
circmean_HTC = atan2(meanvec_HTC(2),meanvec_HTC(1))*180/pi;
circstd_HTC = atan2(stdvec_HTC(2),stdvec_HTC(1))*180/pi;

%% Section 8 - Finding a lenient criteria for a cone weights (Robustness of cone-weight classification)

if ~exist('plot_counter')
    plot_counter = 1;
end

load conewts_svd.mat
load vals.mat

LumIds_conewts = find(conewts_svd(1,:)>0.1 & conewts_svd(3,:)>-0.10 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find((conewts_svd(1,:)<-0.1 & conewts_svd(2,:)>0.1) | (conewts_svd(3,:) <-0.1 & conewts_svd(2,:)>0.1));

LUMidx = LumIds_conewts;
DOidx = [ColorOpponentIds_conewts];
hardtoclassifyidx = 1:size(conewts_svd,2);
LUMidx = LUMidx(vals(LUMidx)<95);
DOidx = DOidx(vals(DOidx)<95);
hardtoclassifyidx([LUMidx DOidx]) = [];

% Checking the correlation with non-linearity indices 
% Load the isoresponse data`
load RSSE_linearmodel_CV.mat % Robust regression
load RSSE_quadmodel_CV.mat
% Load the integration within the subunit data
load AUROClinsubunits_CV.mat
load AUROCquadsubunits_CV.mat

% For storing the Isoresponse NLI
RSSEisoresp_medianofratios = [];

% For storing the white noise NLI
Whitenoise_NLI = [];

indices = [109 24 74];
for ii = 1:numel(AUROClinsubunits) 
  
    % White noise NLI 
    Error_quad = 1-(AUROCquadsubunits{ii});
    Error_lin = 1-(AUROClinsubunits{ii});
    Whitenoise_NLI = [Whitenoise_NLI; log10(median(Error_lin./Error_quad))];
    
    % Isoresponse NLI
    RSSEisoresp_medianofratios = [RSSEisoresp_medianofratios; median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii})];
end


figure(plot_counter);
subplot(324); histogram(RSSEisoresp_medianofratios(LUMidx),logspace(-1,3,31),'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); hold on;
plot(median(RSSEisoresp_medianofratios(LUMidx)),10,'v','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 100],'XTick',[0.1 1 10 100],'Ylim',[0 10],'YTick',[0 5 10]); ylabel('Count'); xlabel('Isoresponse NLI'); axis square; hold off;
subplot(322); histogram(RSSEisoresp_medianofratios(DOidx),logspace(-1,3,31),'FaceColor',[1 0 0],'EdgeColor',[1 1 1]); hold on;
plot(median(RSSEisoresp_medianofratios(DOidx)),20,'v','MarkerSize',8,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 100],'XTick',[0.1 1 10 100],'Ylim',[0 20],'YTick',[0 10 20]); ylabel('Count'); xlabel('Isoresponse NLI'); axis square; hold off;
subplot(326); histogram(RSSEisoresp_medianofratios(hardtoclassifyidx),logspace(-1,3,31),'FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1]); hold on;
plot(median(RSSEisoresp_medianofratios(hardtoclassifyidx)),15,'v','MarkerSize',8,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); 
set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 100],'XTick',[0.1 1 10 100],'Ylim',[0 15],'YTick',[0 5 10 15]); ylabel('Count'); xlabel('Isoresponse NLI'); axis square; hold off;


% Plotting the white noise NLIs on a cell-type basis
figure(plot_counter);
subplot(321); histogram(Whitenoise_NLI(DOidx),-0.02:0.005:0.1,'FaceColor',[1 0 0],'EdgeColor',[1 1 1]); hold on;
plot(median(Whitenoise_NLI(DOidx)),25,'v','MarkerSize',8,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-0.02 0.08],'XTick',-0.02:0.02:0.08,'Ylim',[0 25],'YTick',[0 5 10 15 20 25]); xlabel('Whitenoise NLI'); ylabel('Count'); axis square; hold off;
subplot(323); histogram(Whitenoise_NLI(LUMidx),-0.02:0.005:0.1,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); hold on;
plot(median(Whitenoise_NLI(LUMidx)),10,'v','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-0.02 0.08],'XTick',-0.02:0.02:0.08,'Ylim',[0 10],'YTick',[0 5 10]); xlabel('Whitenoise NLI'); ylabel('Count'); axis square; hold off;
subplot(325); histogram(Whitenoise_NLI(hardtoclassifyidx),-0.02:0.005:0.1,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1]); hold on;
plot(median(Whitenoise_NLI(hardtoclassifyidx)),15,'v','MarkerSize',8,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-0.02 0.08],'XTick',-0.02:0.02:0.08,'Ylim',[0 15],'YTick',[0 5 10 15]); xlabel('Whitenoise NLI'); ylabel('Count'); axis square; hold off;
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;

% Some stats on the whitenoise NLI
[p1,~] = ranksum(Whitenoise_NLI(LUMidx),Whitenoise_NLI(DOidx));
[p2,~] = ranksum(Whitenoise_NLI(DOidx),Whitenoise_NLI(hardtoclassifyidx));
[p3,~] = ranksum(Whitenoise_NLI(LUMidx),Whitenoise_NLI(hardtoclassifyidx));
[median(Whitenoise_NLI(LUMidx)) median(Whitenoise_NLI(DOidx)) median(Whitenoise_NLI(hardtoclassifyidx))]

% Some stats on the isoresponse NLI
[p4,~] = ranksum(log(RSSEisoresp_medianofratios(LUMidx)),log(RSSEisoresp_medianofratios(DOidx)));
[p5,~] = ranksum(log(RSSEisoresp_medianofratios(DOidx)),log(RSSEisoresp_medianofratios(hardtoclassifyidx)));
[p6,~] = ranksum(log(RSSEisoresp_medianofratios(LUMidx)),log(RSSEisoresp_medianofratios(hardtoclassifyidx)));
[p7,~] = ranksum(log(RSSEisoresp_medianofratios([LUMidx DOidx])),log(RSSEisoresp_medianofratios(hardtoclassifyidx)));

[median(log(RSSEisoresp_medianofratios(LUMidx))) median(log(RSSEisoresp_medianofratios(DOidx))) median(log(RSSEisoresp_medianofratios(hardtoclassifyidx)))]

% PLotting the cone weights
figure(plot_counter);
plot(conewts_svd(1,LUMidx), conewts_svd(2,LUMidx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(1,DOidx), conewts_svd(2,DOidx),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(conewts_svd(1,hardtoclassifyidx), conewts_svd(2,hardtoclassifyidx),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
axis equal, set(gca,'Xlim',[-1 1],'Ylim',[0 1]);
plot_counter = plot_counter + 1;
