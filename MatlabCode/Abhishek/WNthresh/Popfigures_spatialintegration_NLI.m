% Same script as Popfigures_spatialintegration.m but with NLI (along with cone weights) used for cell classification as opposed to PC1 
% Author - Abhishek De, 5/20
close all; clearvars;
plot_counter = 1;

%% Figure 1: Impact of linear and non-linear spatial filtering in image processing 

if ~exist('plot_counter')
    plot_counter = 1;
end

I_color = imread('peppers.png');
I_orig = rgb2gray(I_color);
I = im2double(I_orig)-0.5;

% Designing a filter
F = zeros(4,2); S1 = F; S2 = F;
S1(:,1) = 1; S2(:,2) = -1;
F = S1+S2;

% Linear filtered image 
I_linear = conv2(I,F,'same');
I_linear(I_linear<0) = 0;
% I_linear(I_linear>prctile(I_linear(:),90))= 1;
I_linear = I_linear + 0.5-median(I_linear(:));

% Performing the squaring non-linear filtering operation
I_nonlinS1v = conv2(I,S1,'same'); I_S1 = I_nonlinS1v; I_nonlinS1v(I_nonlinS1v<0) = 0; 
I_nonlinS2v = conv2(I,S2,'same'); I_S2 = I_nonlinS2v; I_nonlinS2v(I_nonlinS2v<0) = 0;
I_nonlin = I_nonlinS1v.^2 + I_nonlinS2v.^2; I_nonlin(I_nonlin<0) = 0;
I_nonlin = I_nonlin./max(I_nonlin(:));
I_nonlin(I_nonlin<0) = 0; % Analogous to rectifying non-linearity
% I_nonlin(I_nonlin>prctile(I_nonlin(:),95))= 1;
I_nonlin = I_nonlin + 0.5-median(I_nonlin(:));

% Plotting the figures 
figure(plot_counter); set(gcf,'Name','Analyses of images: grayscale')
subplot(131); imshow(I_orig); axis square; set(gca,'XTick',[],'YTick',[]); title('original pic');
subplot(132); imshow(im2uint8(I_linear)); colormap('gray'); axis square; set(gca,'XTick',[],'YTick',[]); title('Linear');
subplot(133); imshow(im2uint8(I_nonlin)); colormap('gray'); axis square; set(gca,'XTick',[],'YTick',[]); title('Squaring non-linearity');
plot_counter = plot_counter + 1;

%% Figure 2: Staircase procedure from an example cell

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

stro = nex2stro(findfile(char(filename(31))));
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


figure(plot_counter);
cum_flag = []; % to check if I am looking at all the probed directions

t_offset1 = -0.2;
figure(plot_counter);
dir = 3;
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

%% Figure 3: Iso-response data from example LUM, DO and HTC cells 
if ~exist('plot_counter')
    plot_counter = 1;
end

load conewts_svd.mat
load WNsubunitNLI.mat
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
hardtoclassifyidx = [hardtoclassifyidx LUMidx(WNsubunitNLI(LUMidx)>=0) DOidx(WNsubunitNLI(DOidx)>=0)];
LUMidx = LUMidx(WNsubunitNLI(LUMidx)<0);
DOidx = DOidx(WNsubunitNLI(DOidx)<0);

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
indices = [109 31 35];% example LUM, DO and HTC cells
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
   
    figure(1); subplot(numel(indices),4,4*zz-1); image(255*(SpatialRF./(2*max(abs(SpatialRF(:))))+.5)); set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255))
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
    
    figure(3); subplot(1,numel(indices),zz);
    for ii = 1:numel(not_oog_idx)
        h(ii) = plot(x_orig(not_oog_idx(ii)), y_orig(not_oog_idx(ii)),'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]); hold on;
        set(h(ii),'ButtonDownFcn',{@dispimage,x_orig(not_oog_idx(ii)),y_orig(not_oog_idx(ii)),basisvec,bkgnd_monitor});%(x_orig(not_oog_idx(ii)),y_orig(not_oog_idx(ii)),basisvec,bkgnd_monitor)});
    end
    plot(x_lin,y_lin,'k','Linewidth',2); hold on; plot(x_quad2,y_quad2,'r','Linewidth',2); 
    
    if any(outofgamut)
        plot(upsample(x_orig(outofgamut),2),upsample(y_orig(outofgamut),2),'color',[0.5 0.5 0.5]);
    end
    set(gca,'Tickdir','out'); drawnow; axis square; 
    if zz==1 
        set(gca,'Xlim',[-1.6 1.6],'Ylim',[-1.6 1.6],'YTick',-1.6:0.8:1.6,'XTick',[-1.6:0.8:1.6]); hold off;
    elseif zz==2
        set(gca,'Xlim',[-0.6 0.6],'Ylim',[-0.6 0.6],'YTick',-0.6:0.3:0.6,'XTick',[-0.6:0.3:0.6]); hold off;
    elseif zz==3
        set(gca,'Xlim',[-2 2],'Ylim',[-2 2],'YTick',-2:1:2,'XTick',[-2:1:2]); hold off;
    end
    
    % Now, I am going to plot the
    wts  = [1 0;1 1;0 1;-1 1;-1 0;-1 -1;0 -1;1 -1];
    subplotidxs = [1;2;3;4;6;7;8;9];
    for ii = 1:numel(subplotidxs)
        vec = wts(ii,1)*(basisvec{1}-bkgnd_monitor) + wts(ii,2)*(basisvec{2}-bkgnd_monitor);
        normfactor = 0.5/(max(abs(vec(:)))+0.01);
        vec = normfactor*vec + 0.5;
        figure(2),subplot(numel(indices),size(wts,1),size(wts,1)*(zz-1)+ii); image(vec); set(gca,'XTick',[],'YTick',[]); axis square;
    end
    plot_counter = plot_counter + 1;
end

% Obtaining the target firing rates of the example neurons
load TFR.mat
FR = TFR(1,indices);

%% Figure 4: Population analysis of isoresponse curves 
if ~exist('plot_counter')
    plot_counter = 1;
end

load conewts_svd.mat
load WNsubunitNLI.mat
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
hardtoclassifyidx = [hardtoclassifyidx LUMidx(WNsubunitNLI(LUMidx)>=0) DOidx(WNsubunitNLI(DOidx)>=0)];
LUMidx = LUMidx(WNsubunitNLI(LUMidx)<0);
DOidx = DOidx(WNsubunitNLI(DOidx)<0);

% Checking the correlation with non-linearity indices 
% Load the isoresponse data
load RSSE_linearmodel_CV.mat % Robust regression
load RSSE_quadmodel_CV.mat

% For storing median of differences/ratios
RSSEisoresp_medianofratios = [];
RSSEisoresp_lin_median = []; RSSEisoresp_quad_median = []; % Isoresponse data

for ii = 1:numel(RSSE_linearmodel)   
    % computation for calculating median of differences/ratios
    RSSEisoresp_medianofratios = [RSSEisoresp_medianofratios; median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii})];  
    
    
    RSSEisoresp_lin_median = [RSSEisoresp_lin_median; median(RSSE_linearmodel{ii})];
    RSSEisoresp_quad_median = [RSSEisoresp_quad_median; median(RSSE_quadmodel{ii})];
end
RSSEisoresp_medianofratios(RSSEisoresp_medianofratios<0.1) = 0.1;
indices = [109 31 35];
% Plotting the results for SVD based cone weight classification including the PC1 z-scores 
figure(plot_counter);
subplot(411); plot(RSSEisoresp_lin_median(hardtoclassifyidx),RSSEisoresp_quad_median(hardtoclassifyidx),'o','Markersize',4,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(RSSEisoresp_lin_median(LUMidx),RSSEisoresp_quad_median(LUMidx),'o','Markersize',4,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(RSSEisoresp_lin_median(indices(1)),RSSEisoresp_quad_median(indices(1)),'o','Markersize',4,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 1 0]);
plot(RSSEisoresp_lin_median(DOidx),RSSEisoresp_quad_median(DOidx),'o','Markersize',4,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); plot([0.0001 10],[0.0001 10],'k');
plot(RSSEisoresp_lin_median(indices(2)),RSSEisoresp_quad_median(indices(2)),'o','Markersize',4,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 1 0]);
plot(RSSEisoresp_lin_median(indices(3)),RSSEisoresp_quad_median(indices(3)),'o','Markersize',4,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0 1 0]);
axis square; set(gca,'Tickdir','out','Xlim',[0.0001 10],'Ylim',[0.0001 10],'YScale','log','XScale','log','XTick',[0.0001 0.001 0.01 0.1 1 10],'YTick',[0.0001 0.001 0.01 0.1 1 10]); xlabel('Linear error'); ylabel('Quadratic error'); title('Isoresponse'); hold off;
subplot(412); histogram(RSSEisoresp_medianofratios(LUMidx),logspace(-1,3,31),'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); hold on;
plot(median(RSSEisoresp_medianofratios(LUMidx)),15,'v','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
plot(RSSEisoresp_medianofratios(indices(1)),14,'s','MarkerSize',8,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 100],'XTick',[0.1 1 10 100],'Ylim',[0 15],'YTick',[0 5 10 15]); ylabel('Count'); title('Isoresponse'); xlabel('median CV ratio of errors'); axis square; hold off;
subplot(413); histogram(RSSEisoresp_medianofratios(DOidx),logspace(-1,3,31),'FaceColor',[1 0 0],'EdgeColor',[1 1 1]); hold on;
plot(median(RSSEisoresp_medianofratios(DOidx)),10,'v','MarkerSize',8,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(RSSEisoresp_medianofratios(indices(2)),9,'s','MarkerSize',8,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 100],'XTick',[0.1 1 10 100],'Ylim',[0 10],'YTick',[0 5 10]); ylabel('Count'); title('Isoresponse'); xlabel('median CV ratio of errors'); axis square; hold off;
subplot(414); histogram(RSSEisoresp_medianofratios(hardtoclassifyidx),logspace(-1,3,31),'FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1]); hold on;
plot(median(RSSEisoresp_medianofratios(hardtoclassifyidx)),15,'v','MarkerSize',8,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); 
plot(RSSEisoresp_medianofratios(indices(3)),14,'s','MarkerSize',8,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 100],'XTick',[0.1 1 10 100],'Ylim',[0 15],'YTick',[0 5 10 15]); ylabel('Count'); title('Isoresponse'); xlabel('median CV ratio of errors'); axis square; hold off;
plot_counter = plot_counter + 1;

% Comparing spatial NLI across cell types
group = [ones(size(LUMidx)) 2*ones(size(DOidx)) 3*ones(size(hardtoclassifyidx))];
data = RSSEisoresp_medianofratios([LUMidx'; DOidx'; hardtoclassifyidx']); 
p1 = kruskalwallis(data,group);

% Comparing spatial NLI between simple and DO cells
[p2,h] = ranksum(RSSEisoresp_medianofratios(LUMidx),RSSEisoresp_medianofratios(DOidx));

% Some control analyses to check the relationship between spatial structure and non-linearity
% Classifying the spatial structure as 1 (center-surround) or 2 (adjacent) based on selection of subunits
RFstructure = [2 1 2 1 2 2 1 2 1 1 2,...
               2 2 1 2 2 2 2 2 2 1 1,...
               2 2 2 2 2 2 2 2 2 1 2,...
               2 1 2 2 2 2 1 2 2 1 2,...
               2 1 2 1 1 1 1 2 2 2 2,...
               2 2 2 1 1 2 1 2 2 2 1,...
               1 1 1 2 2 2 1 2 2 2 2,...
               2 2 2 1 1 1 1 1 2 2 2,...
               1 1 2 2 1 2 2 2 2 2 2,...
               2 2 2 1 2 1 2 2 1 2 2,...
               2 2 1 2 2 2];
                   
% Checking the link between subunit geometry and spatial NLI
[p3,h] = ranksum(RSSEisoresp_medianofratios(RFstructure==1),RSSEisoresp_medianofratios(RFstructure==2));
           
% Loading the CV errors from Gabor and DoG fits 
load IsoresponseRF_Peraccuracy.mat
load IsoresponseRF_SSE.mat
load IsoresponseRF_Deviation.mat
meanperaccuracy = zeros(size(Peraccuracy));
meanSSE = zeros(size(SSE));
meanR = zeros(size(SSE));
for ii = 1:size(Peraccuracy,1)
    for jj = 1:size(Peraccuracy,2)
        meanperaccuracy(ii,jj) = mean(Peraccuracy{ii,jj});
        meanSSE(ii,jj) = mean(SSE{ii,jj});
        meanR(ii,jj) = mean(cos(Deviation{ii,jj}*pi/180));
    end
end

diffGaborDoG = meanR(:,2)-meanR(:,3); % Difference between Pearson's r (Gabor-DoG)

% [center-surround & Gabor   center-surround & DoG;...
%  Adjacent & Gabor          Adjacent & DoG]
X = [sum(RFstructure'==1 & diffGaborDoG>=0) sum(RFstructure'==1 & diffGaborDoG<0);...
    sum(RFstructure'==2 & diffGaborDoG>=0) sum(RFstructure'==2 & diffGaborDoG<0)];

[h,p4] = fishertest(X);

%% Figure 5: Conceptual model of cone-signal integration
if ~exist('plot_counter')
    plot_counter = 1;
end

x = linspace(-1,1,21);
y = linspace(-1,1,21);
z = linspace(-1,2,21);
[X,Y,Z] = meshgrid(x,y,z);
V1 = 3*X-2*Y+2*Z;
V2 = +1.5*max(0,1-X).^2+1.5*max(0,1-Y).^2+2*max(0,Z).^2;

figure(plot_counter); set(gcf,'Name','Cone signal integration');
subplot(121); p1 = patch(isosurface(X,Y,Z,V1,prctile(V1(:),60))); hold on; set(p1,'EdgeColor','none','FaceAlpha',0.5,'EdgeAlpha',0.5);
axis square; set(gca,'Tickdir','out','XTick',[-1 0 1],'YTick',[-1 0 1],'ZTick',[-1 0.5 2]); 
view(-22.2000, 42.8000); grid on; camlight('headlight'); xlabel('L'); ylabel('M'); zlabel('S');
subplot(122); p2 = patch(isosurface(X,Y,Z,V2,prctile(V2(:),30)));  hold on; set(p2,'EdgeColor','none','FaceAlpha',0.5,'EdgeAlpha',0.5);
axis square; set(gca,'Tickdir','out','XTick',[-1 0 1],'YTick',[-1 0 1],'ZTick',[-1 0.5 2]); 
view(-22.2000, 42.8000); grid on;  xlabel('L'); ylabel('M'); zlabel('S');
plot_counter = plot_counter + 1;

%% Figure 6: Population analysis of cone signal combination within individual subunits

if ~exist('plot_counter')
    plot_counter = 1;
end

load conewts_svd.mat
load WNsubunitNLI.mat
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
hardtoclassifyidx = [hardtoclassifyidx LUMidx(WNsubunitNLI(LUMidx)>=0) DOidx(WNsubunitNLI(DOidx)>=0)];
LUMidx = LUMidx(WNsubunitNLI(LUMidx)<0);
DOidx = DOidx(WNsubunitNLI(DOidx)<0);

% Load the integration within the subunit data
load AUROClinS1_CV.mat
load AUROCquadS1_CV.mat
load AUROClinS2_CV.mat
load AUROCquadS2_CV.mat

% For storing median of differences/ratios
Withinsubunits_medianofdifferences = [];

for ii = 1:numel(RSSE_linearmodel)      
    Withinsubunits_medianofdifferences = [Withinsubunits_medianofdifferences; median([median(AUROCquad1{ii}-AUROClin1{ii}) median(AUROCquad2{ii}-AUROClin2{ii})])];
end

indices = [109 31 35];
% Plotting the results for SVD based cone weight classification including the PC1 z-scores 
figure(plot_counter);
subplot(311); histogram(100*(Withinsubunits_medianofdifferences(LUMidx)),linspace(-2,8,21),'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); hold on;
plot(median(100*(Withinsubunits_medianofdifferences(LUMidx))),20,'v','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
plot(100*(Withinsubunits_medianofdifferences(indices(1))),18,'s','MarkerSize',8,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-2 8],'Ylim',[0 20],'YTick',[0 10 20]); xlabel('median CV GQM-GLM AUROC'); ylabel('Count'); title('Within subunits'); axis square; hold off;
subplot(312); histogram(100*(Withinsubunits_medianofdifferences(DOidx)),linspace(-2,8,21),'FaceColor',[1 0 0],'EdgeColor',[1 1 1]); hold on;
plot(median(100*(Withinsubunits_medianofdifferences(DOidx))),15,'v','MarkerSize',8,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(100*(Withinsubunits_medianofdifferences(indices(2))),12,'s','MarkerSize',8,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-2 8],'Ylim',[0 15],'YTick',[0 5 10 15]); xlabel('median CV GQM-GLM AUROC'); ylabel('Count'); title('Within subunits'); axis square; hold off;
subplot(313); histogram(100*(Withinsubunits_medianofdifferences(hardtoclassifyidx)),linspace(-2,8,21),'FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1]); hold on;
plot(median(100*(Withinsubunits_medianofdifferences(hardtoclassifyidx))),15,'v','MarkerSize',8,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
plot(100*(Withinsubunits_medianofdifferences(indices(3))),14,'s','MarkerSize',8,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-2 8],'Ylim',[0 15],'YTick',[0 5 10 15]); xlabel('median CV GQM-GLM AUROC'); ylabel('Count'); title('Within subunits'); axis square; hold off;
plot_counter = plot_counter + 1;

% Comparing spatial NLI across cell types
group = [ones(size(LUMidx)) 2*ones(size(DOidx)) 3*ones(size(hardtoclassifyidx))];
data = Withinsubunits_medianofdifferences([LUMidx'; DOidx'; hardtoclassifyidx']); 
p1 = kruskalwallis(data,group);

% Loading the LGN data: calculate the median value 
load Withinsubunits_medianofdifferences_LGN.mat


%% Figure 7: Relationship between isoresponse data and cone signal integration

if ~exist('plot_counter')
    plot_counter = 1;
end

load conewts_svd.mat
load WNsubunitNLI.mat
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
hardtoclassifyidx = [hardtoclassifyidx LUMidx(WNsubunitNLI(LUMidx)>=0) DOidx(WNsubunitNLI(DOidx)>=0)];
LUMidx = LUMidx(WNsubunitNLI(LUMidx)<0);
DOidx = DOidx(WNsubunitNLI(DOidx)<0);

% Checking the correlation with non-linearity indices 
% Load the isoresponse data
load RSSE_linearmodel_CV.mat % Robust regression
load RSSE_quadmodel_CV.mat

% Load the integration within the subunit data
load AUROClinS1_CV.mat
load AUROCquadS1_CV.mat
load AUROClinS2_CV.mat
load AUROCquadS2_CV.mat

% For storing median of differences/ratios
RSSEisoresp_medianofratios = [];
Withinsubunits_medianofdifferences = [];

for ii = 1:numel(RSSE_linearmodel)   
    % computation for calculating median of differences/ratios
    RSSEisoresp_medianofratios = [RSSEisoresp_medianofratios; median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii})];
    Withinsubunits_medianofdifferences = [Withinsubunits_medianofdifferences; median([median(AUROCquad1{ii}-AUROClin1{ii}) median(AUROCquad2{ii}-AUROClin2{ii})])];   
end

indices = [109 31 35];
figure(plot_counter); plot(RSSEisoresp_medianofratios(LUMidx),100*(Withinsubunits_medianofdifferences(LUMidx)),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(RSSEisoresp_medianofratios(indices(1)),100*(Withinsubunits_medianofdifferences(indices(1))),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 1 0]); hold on;
plot(RSSEisoresp_medianofratios(DOidx),100*(Withinsubunits_medianofdifferences(DOidx)),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(RSSEisoresp_medianofratios(indices(2)),100*(Withinsubunits_medianofdifferences(indices(2))),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 1 0]); hold on;
plot(RSSEisoresp_medianofratios(hardtoclassifyidx),100*(Withinsubunits_medianofdifferences(hardtoclassifyidx)),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
plot(RSSEisoresp_medianofratios(indices(3)),100*(Withinsubunits_medianofdifferences(indices(3))),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0 1 0]); hold on;
axis square; set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 1000],'XTick',[0.1 1 10 100 1000],'Ylim',[-2 8],'YTick',-2:2:8); xlabel('Spatial NLI isoresponse'); ylabel('Within AUROC Quad-Lin');
plot_counter = plot_counter + 1;

[r,p] = corr(RSSEisoresp_medianofratios,Withinsubunits_medianofdifferences,'type','Spearman');

%% Figure 8: Impact of S-cone input on signal integration: Isoresponse data and cone signal integration

if ~exist('plot_counter')
    plot_counter = 1;
end

load conewts_svd.mat
load WNsubunitNLI.mat
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
hardtoclassifyidx = [hardtoclassifyidx LUMidx(WNsubunitNLI(LUMidx)>=0) DOidx(WNsubunitNLI(DOidx)>=0)];
LUMidx = LUMidx(WNsubunitNLI(LUMidx)<0);
DOidx = DOidx(WNsubunitNLI(DOidx)<0);

% Checking whether there is any correlation between S-cone input and non-linearities 

% Load the isoresponse data
load RSSE_linearmodel_CV.mat % Robust regression
load RSSE_quadmodel_CV.mat

% Load the integration within the subunit data
load AUROClinS1_CV.mat
load AUROCquadS1_CV.mat
load AUROClinS2_CV.mat
load AUROCquadS2_CV.mat

% For storing median of differences/ratios
RSSEisoresp_medianofratios = [];
Withinsubunits_medianofdifferences = [];

for ii = 1:numel(RSSE_linearmodel)   
    % computation for calculating median of differences/ratios
    RSSEisoresp_medianofratios = [RSSEisoresp_medianofratios; median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii})];
    Withinsubunits_medianofdifferences = [Withinsubunits_medianofdifferences; median([median(AUROCquad1{ii}-AUROClin1{ii}) median(AUROCquad2{ii}-AUROClin2{ii})])];   
end

indices = [109 31 35];
figure(plot_counter); subplot(121); plot(RSSEisoresp_medianofratios(LUMidx),abs(conewts_svd(3,LUMidx)),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(RSSEisoresp_medianofratios(indices(1)),abs(conewts_svd(3,indices(1))),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 1 0]); hold on;
plot(RSSEisoresp_medianofratios(DOidx),abs(conewts_svd(3,DOidx)),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(RSSEisoresp_medianofratios(indices(2)),abs(conewts_svd(3,indices(2))),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 1 0]); hold on;
plot(RSSEisoresp_medianofratios(hardtoclassifyidx),abs(conewts_svd(3,hardtoclassifyidx)),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
plot(RSSEisoresp_medianofratios(indices(3)),abs(conewts_svd(3,indices(3))),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0 1 0]); hold on;
axis square; set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 1000],'XTick',[0.1 1 10 100 1000],'Ylim',[0 0.6],'YTick',0:0.2:0.6); xlabel('Spatial NLI isoresponse'); ylabel('S cone input');
subplot(122); plot(100*Withinsubunits_medianofdifferences(LUMidx),abs(conewts_svd(3,LUMidx)),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(100*Withinsubunits_medianofdifferences(indices(1)),abs(conewts_svd(3,indices(1))),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 1 0]); hold on;
plot(100*Withinsubunits_medianofdifferences(DOidx),abs(conewts_svd(3,DOidx)),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(100*Withinsubunits_medianofdifferences(indices(2)),abs(conewts_svd(3,indices(2))),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 1 0]); hold on;
plot(100*Withinsubunits_medianofdifferences(hardtoclassifyidx),abs(conewts_svd(3,hardtoclassifyidx)),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
plot(100*Withinsubunits_medianofdifferences(indices(3)),abs(conewts_svd(3,indices(3))),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0 1 0]); hold on;
axis square; set(gca,'Tickdir','out','Xlim',[-2 8],'XTick',-2:2:8,'Ylim',[0 0.6],'YTick',0:0.2:0.6); xlabel('Within AUROC Quad-Lin'); ylabel('S cone input');
plot_counter = plot_counter + 1;

[r1,p1] = corr(RSSEisoresp_medianofratios,abs(conewts_svd(3,:)'),'type','Spearman');
[r2,p2] = corr(100*Withinsubunits_medianofdifferences,abs(conewts_svd(3,:)'),'type','Spearman');

%% Figure S2: LGN analysis 
% Check the code LGN_analyses/RGB3D_analysesonLGN.m


%% Determining the some additional numbers

% Number of neurons from each Monkey
load Pangu.mat
load Maui.mat
NEURONS_MONKEY1 = numel(Maui);
NEURONS_MONKEY2 = numel(Pangu);

% RF locations 
load RF_loc.mat
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
close(conn);
RF_LOCATIONS = RF_loc(strcmp(string(NTmode),"subunit"),:);
RF_AMP = sqrt(sum(RF_LOCATIONS.^2,2))/10;