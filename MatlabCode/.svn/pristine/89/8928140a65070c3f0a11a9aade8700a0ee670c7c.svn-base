% This is an analysis on the WhiteNoise subunit data suggested by Greg
% Author - Abhishek De, 12/17
close all; clearvars;
load filename_c.mat % Color cells
load filename_l.mat % Luminance cells
load S1LMS.mat % cone wts for subunit 1
load S2LMS.mat % cone wts for subunit 2
[OCidx, LMidx, LUMidx, SOidx, hardtoclassifyidx] = classifycells(S1LMS,S2LMS);
idx = [OCidx LMidx LUMidx SOidx hardtoclassifyidx];
filename = [filename_c; filename_l];
plot_counter = 1;

stro = nex2stro(findfile(char(filename(idx(1),:))));
global spikename maskidx spikeidx nstixperside ngammasteps seedidx nframesidx
global fponidx stimoffidx stimonidx muidxs sigmaidxs basisvecidx fpacqidx latencyidx basisvecdiridx neurothreshidx targetspikerateidx
global msperframe ntrials maxT xx yy M
spikename = getSpikenum(stro);
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
out = STCOV_st('return'); % returns the covariance matrix on frame by frame basis
STS = out{1};  % A (dimension) x 9(frames) matrix
STCross = out{2};  % A (dimension x frames)x (dimension x frames)  matrix
nspikes = out{3}; % Number of spikes in the given file
clear STCOV_st out;
% Coverting the STS and the STCross into STA and STC respectively
STAs = STS/nspikes;
basisvec1ST = STAs(1:2:end,:); % temporal phosphor dynamics of basisvec 1
basisvec2ST = STAs(2:2:end,:); % temporal phosphor dynamics of basisvec 2
tmp = STS(:)*STS(:)';
STCs = (nspikes.*STCross-tmp)/(nspikes*(nspikes-1));
P = eye(size(STCs)) - STAs(:)*inv(STAs(:)'*STAs(:))*STAs(:)'; % Subtracting such that the STA is orthogonal to PC
STCs = P*STCs*P';
[tmp,d] = eig(STCs);
eig_PC = sort(diag(d)); % storing all the eigenvalues
v = real(tmp);
[~, idxs] = sort(diag(d));
v = v(:,idxs);
PC1 = v(:,end);
PC1 = reshape(PC1,[6 maxT]);

[~,whichframe1] = max(fliplr(sum(basisvec1ST.^2,1)));
[~,whichframe2] = max(fliplr(sum(basisvec2ST.^2,1)));
[~,whichframe3] = max(fliplr(sum(PC1(1:2:end,:).^2,1)));
[~,whichframe4] = max(fliplr(sum(PC1(2:2:end,:).^2,1)));

stindices = find(cum_n>0);
stindices1 = stindices - whichframe1+1;
stindices2 = stindices - whichframe2+1;
rawrgb_onbasisvec1 = cum_rgbs(1:2:end,:);
rawrgb_onbasisvec2 = cum_rgbs(2:2:end,:);

min_val = floor(min(cum_rgbs(:))*100)/100;
max_val = ceil(max(cum_rgbs(:))*100)/100;
num_bins = 11;
bin_interval = (max_val-min_val)/num_bins;
nbins1 = min_val:bin_interval:max_val;
nbins = linspace(prctile(cum_rgbs(:),5), prctile(cum_rgbs(:),95),numel(nbins1)+2);
xmin = min(nbins); xmax = max(nbins);
new_nbins = [nbins(1)-mean(diff(nbins)), nbins, mean(diff(nbins)) + nbins(end)];

% This analysis is similar to what Charlie and Greg tried on their 3D
% gaussian cloud to observe how much signal is present

meanrawrgb1 = mean(rawrgb_onbasisvec1,2);
meanrawrgb2 = mean(rawrgb_onbasisvec2,2);
meanspikergb1 = mean(rawrgb_onbasisvec1(:,stindices1),2);
meanspikergb2 = mean(rawrgb_onbasisvec2(:,stindices2),2);
dir1 = meanspikergb1 - meanrawrgb1;
dir2 = meanspikergb2 - meanrawrgb2;
projrawrgb1 = dir1'*rawrgb_onbasisvec1;
projspikergb1 = dir1'*rawrgb_onbasisvec1(:,stindices1);
projrawrgb2 = dir2'*rawrgb_onbasisvec2;
projspikergb2 = dir2'*rawrgb_onbasisvec2(:,stindices2);

% Extracting the basis vectors
neurothreshmode = stro.trial(:,neurothreshidx);
basisvec_dropidx = inds(end);
neurothresh_startidx = max(find(neurothreshmode==1,1),basisvec_dropidx+1); % There is possibility of pressing the Neurothresh button multiple times
num_targetspikerates = unique(stro.trial(neurothresh_startidx:end,targetspikerateidx));
vect = stro.ras{basisvec_dropidx,basisvecidx};
basisvec_size = nstixperside*nstixperside*3;
numvect = (numel(vect)/basisvec_size)-1;
basisvec = cell(1,numvect);
bkgnd_monitor = reshape(vect(numvect*basisvec_size+1:basisvec_size*(numvect+1)),[nstixperside nstixperside 3]);
% plotting the results
figure(plot_counter);
for ii = 1:numvect
    tmp_vec = vect((ii-1)*basisvec_size+1:basisvec_size*ii) + vect(numvect*basisvec_size+1:basisvec_size*(numvect+1));
    basisvec{ii} = reshape(tmp_vec,[nstixperside nstixperside 3]);
    
end
% Some additional calculations
basisvec1mask = org_mask;
basisvec1mask(basisvec1mask == Inf) = 0;
basisvec1mask(basisvec1mask == 2) = 0;
basisvec2mask = org_mask;
basisvec2mask(basisvec2mask == Inf) = 0;
basisvec2mask(basisvec2mask == 1) = 0;
m = basisvec{1}-bkgnd_monitor;
m = squeeze(m(:,:,1));
if ~(m(:)'*basisvec1mask)
    tmp = basisvec{2};
    basisvec{2} = basisvec{1};
    basisvec{1} = tmp;
end

BinEdges = linspace(-0.04,0.04,51);
figure(plot_counter); subplot(521),image(basisvec{1}); set(gca,'Xtick',[],'YTick',[]);
subplot(522),image(basisvec{2}); set(gca,'Xtick',[],'YTick',[]);
subplot(523),plot(basisvec1ST(1,:),'r','Linewidth',2); hold on; plot(basisvec1ST(2,:),'g','Linewidth',2);
plot(basisvec1ST(3,:),'b','Linewidth',2); plot(size(basisvec1ST,2)-whichframe1+1,0,'kv','MarkerFaceColor','g'); hold off;
subplot(524),plot(basisvec2ST(1,:),'r','Linewidth',2); hold on; plot(basisvec2ST(2,:),'g','Linewidth',2);
plot(basisvec2ST(3,:),'b','Linewidth',2); plot(size(basisvec2ST,2)-whichframe2+1,0,'kv','MarkerFaceColor','g'); hold off;
subplot(525),plot(PC1(1,:),'r','Linewidth',2); hold on; plot(PC1(3,:),'g','Linewidth',2); plot(PC1(5,:),'b','Linewidth',2);
plot(size(PC1,2)-whichframe3+1,0,'kv','MarkerFaceColor','r'); hold off;
subplot(526),plot(PC1(2,:),'r','Linewidth',2); hold on; plot(PC1(4,:),'g','Linewidth',2); plot(PC1(6,:),'b','Linewidth',2);
plot(size(PC1,2)-whichframe4+1,0,'kv','MarkerFaceColor','r'); hold off;
subplot(527), histogram(projrawrgb1,'BinEdges',BinEdges,'Normalization','probability'); hold on;
histogram(projspikergb1,'BinEdges',BinEdges,'Normalization','probability','FaceColor','r'); title('Basisvec1'); hold off;
subplot(528), histogram(projrawrgb2,'BinEdges',BinEdges,'Normalization','probability'); hold on;
histogram(projspikergb2,'BinEdges',BinEdges,'Normalization','probability','FaceColor','r'); title('Basisvec2'); hold off;
subplot(529),ecdf(projrawrgb1), hold on, ecdf(projspikergb1), title('Basisvec1'); legend('Raw','Spike','Location','northwest'); hold off;
subplot(5,2,10),ecdf(projrawrgb2), hold on, ecdf(projspikergb2), title('Basisvec2'); legend('Raw','Spike','Location','northwest'); hold off;
plot_counter = plot_counter + 1;

% This is the where the I begin looking for how cone signal combination takes place
% For Basisvec 1
[n_spike1_rg,n_raw1_rg,non_lin1_rg] = calcplanarspikestats(rawrgb_onbasisvec1([1 2],:),stindices1,new_nbins); % RG-plane
[n_spike1_gb,n_raw1_gb,non_lin1_gb] = calcplanarspikestats(rawrgb_onbasisvec1([2 3],:),stindices1,new_nbins); % GB-plane
[n_spike1_rb,n_raw1_rb,non_lin1_rb] = calcplanarspikestats(rawrgb_onbasisvec1([1 3],:),stindices1,new_nbins); % RB-plane

% For Basisvec 2 (only RG plane)
[n_spike2_rg,n_raw2_rg,non_lin2_rg] = calcplanarspikestats(rawrgb_onbasisvec2([1 2],:),stindices2,new_nbins); % RG-plane
[n_spike2_gb,n_raw2_gb,non_lin2_gb] = calcplanarspikestats(rawrgb_onbasisvec2([2 3],:),stindices2,new_nbins); % GB-plane
[n_spike2_rb,n_raw2_rb,non_lin2_rb] = calcplanarspikestats(rawrgb_onbasisvec1([1 3],:),stindices2,new_nbins); % RB-plane

% I am trying out a new analysis suggested by Greg, Am blurring all the points in 3-D space and cutting slices to see if if there is
% any evidence of linear-nonlinear combination of gun/cone signals
numdivisions = 9;
Vert1basisvec1 = linspace(min(rawrgb_onbasisvec1(1,:)),max(rawrgb_onbasisvec1(1,:)),numdivisions);
Vert2basisvec1 = linspace(min(rawrgb_onbasisvec1(2,:)),max(rawrgb_onbasisvec1(2,:)),numdivisions);
Vert3basisvec1 = linspace(min(rawrgb_onbasisvec1(3,:)),max(rawrgb_onbasisvec1(3,:)),numdivisions);
[count1 edges1 mid1 loc1] = histcn(rawrgb_onbasisvec1',Vert1basisvec1,Vert2basisvec1,Vert3basisvec1);
[count2 edges2 mid2 loc2] = histcn(rawrgb_onbasisvec1(:,stindices1)',Vert1basisvec1,Vert2basisvec1,Vert3basisvec1);

Vert1basisvec2 = linspace(min(rawrgb_onbasisvec2(1,:)),max(rawrgb_onbasisvec2(1,:)),numdivisions);
Vert2basisvec2 = linspace(min(rawrgb_onbasisvec2(2,:)),max(rawrgb_onbasisvec2(2,:)),numdivisions);
Vert3basisvec2 = linspace(min(rawrgb_onbasisvec2(3,:)),max(rawrgb_onbasisvec2(3,:)),numdivisions);
[count3 edges3 mid3 loc3] = histcn(rawrgb_onbasisvec2',Vert1basisvec2,Vert2basisvec2,Vert3basisvec2);
[count4 edges4 mid4 loc4] = histcn(rawrgb_onbasisvec2(:,stindices2)',Vert1basisvec2,Vert2basisvec2,Vert3basisvec2);
count1 = imgaussfilt3(count1,2);
count2 = imgaussfilt3(count2,2);
count3 = imgaussfilt3(count3,2);
count4 = imgaussfilt3(count4,2);
ratio1 = count2./count1;
ratio1(isnan(ratio1)) = 0;
ratio2 = count4./count3;
ratio2(isnan(ratio2)) = 0;
% ratio1 = imgaussfilt3(ratio1,2);
% ratio2 = imgaussfilt3(ratio2,2);
plotmode = 0; % 0 - don't plot, 1 - plot
% converting LMS mechs into RGB mechs to see if those directions have anything to do with the planar/non-planar isosurfaces
LMSmechs = [1 -1 0; 1 1 0; 0 0 1];
RGBmechdirs = 5*(inv(M)*LMSmechs'); % multiplying with an arbitrary number in order to scale the lines (cardinal mechanisms);

if plotmode
    for ii = 1:size(count1,3)
        figure(plot_counter), subplot(ceil(sqrt(size(count1,3))),ceil(sqrt(size(count1,3))),ii), imagesc(squeeze(ratio1(ii,:,:)));  set(gca,'XTick',[],'YTick',[]); axis square; axis xy; xlabel('B'),ylabel('G');
        figure(plot_counter+1), subplot(ceil(sqrt(size(count2,3))),ceil(sqrt(size(count2,3))),ii), imagesc(squeeze(ratio1(:,ii,:))); set(gca,'XTick',[],'YTick',[]); axis square; axis xy; xlabel('B'),ylabel('R');
        figure(plot_counter+2), subplot(ceil(sqrt(size(count2,3))),ceil(sqrt(size(count2,3))),ii), imagesc(squeeze(ratio1(:,:,ii))); set(gca,'XTick',[],'YTick',[]); axis square; axis xy; xlabel('G'),ylabel('R');
        figure(plot_counter+3), subplot(ceil(sqrt(size(count1,3))),ceil(sqrt(size(count1,3))),ii), imagesc(squeeze(ratio2(ii,:,:))); set(gca,'XTick',[],'YTick',[]); axis square; axis xy; xlabel('B'),ylabel('G');
        figure(plot_counter+4), subplot(ceil(sqrt(size(count2,3))),ceil(sqrt(size(count2,3))),ii), imagesc(squeeze(ratio2(:,ii,:))); set(gca,'XTick',[],'YTick',[]); axis square; axis xy; xlabel('B'),ylabel('R');
        figure(plot_counter+5), subplot(ceil(sqrt(size(count2,3))),ceil(sqrt(size(count2,3))),ii), imagesc(squeeze(ratio2(:,:,ii))); set(gca,'XTick',[],'YTick',[]); axis square; axis xy; xlabel('G'),ylabel('R');
    end
    figure(plot_counter); set(gcf,'Name','Basisvec1 BG');
    figure(plot_counter+1); set(gcf,'Name','Basisvec1 BR');
    figure(plot_counter+2); set(gcf,'Name','Basisvec1 GR');
    figure(plot_counter+3); set(gcf,'Name','Basisvec2 BG');
    figure(plot_counter+4); set(gcf,'Name','Basisvec2 BR');
    figure(plot_counter+5); set(gcf,'Name','Basisvec2 GR');
    plot_counter = plot_counter + 6;
end

numsurfaces = 3;
cgradient = linspace(0,1,numsurfaces);
fv1 = cell(1,numsurfaces);
fv2 = cell(1,numsurfaces);
isovals1 = linspace(prctile(ratio1(:),25),prctile(ratio1(:),75),numsurfaces);
isovals2 = linspace(prctile(ratio2(:),25),prctile(ratio2(:),75),numsurfaces);
for ii = 1:numsurfaces
    fv1{ii} = isosurface(mid1{1},mid1{2},mid1{3},ratio1,isovals1(ii));
    fv2{ii} = isosurface(mid3{1},mid3{2},mid3{3},ratio2,isovals2(ii));
    %     figure(plot_counter);isosurface(mid1{1},mid1{2},mid1{3},ratio1,isovals1(ii)); drawnow; hold on;
    %     figure(plot_counter+1);isosurface(mid3{1},mid3{2},mid3{3},ratio2,isovals2(ii));drawnow; hold on;
end

responseS1 = -1*ones(1,size(rawrgb_onbasisvec1,2)); responseS1(stindices1) = 1;
responseS2 = -1*ones(1,size(rawrgb_onbasisvec2,2)); responseS2(stindices2) = 1;
covS1spike = cov(rawrgb_onbasisvec1(:,responseS1==1)');
covS1nospike = cov(rawrgb_onbasisvec1(:,responseS1==-1)');
muS1spike = mean(rawrgb_onbasisvec1(:,responseS1==1),2);
muS1nospike = mean(rawrgb_onbasisvec1(:,responseS1==-1),2);
covS2spike = cov(rawrgb_onbasisvec2(:,responseS2==1)');
covS2nospike = cov(rawrgb_onbasisvec2(:,responseS2==-1)');
muS2spike = mean(rawrgb_onbasisvec2(:,responseS2==1),2);
muS2nospike = mean(rawrgb_onbasisvec2(:,responseS2==-1),2);

% I am writing this script to understand if the boundary between two separated gaussian a plane or a curved surface
x = linspace(-0.15,0.15,51);
[X,Y,Z] = meshgrid(x,x,x);
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
figure(plot_counter); set(gcf,'Name','Isosurfaces from Ratio of Gaussians');
for ii = 1:3
    fv1 = isosurface(X,Y,Z,ratio1,val1(ii));
    fv2 = isosurface(X,Y,Z,ratio2,val2(ii));
    fv1lms = M*fv1.vertices'; fv2lms = M*fv2.vertices';
    subplot(121); plot3(fv1.vertices(:,1),fv1.vertices(:,2),fv1.vertices(:,3),'o','MarkerSize',2,'LineWidth',0.5); hold on;
    subplot(122); plot3(fv2.vertices(:,1),fv2.vertices(:,2),fv2.vertices(:,3),'o','MarkerSize',2,'LineWidth',0.5); hold on;
    
    if ii==2
        [v1,d1] = eig(cov(fv1.vertices));
        [v2,d2] = eig(cov(fv2.vertices));
    end
end
subplot(121); plot3(2*[-STAs(1,maxT-whichframe1+1) STAs(1,maxT-whichframe1+1)],2*[-STAs(3,maxT-whichframe1+1) STAs(3,maxT-whichframe1+1)],2*[-STAs(5,maxT-whichframe1+1) STAs(5,maxT-whichframe1+1)],'-k','Linewidth',4);
xlabel('R'), ylabel('G'); zlabel('B'); title('Subunit 1'); axis equal; hold off;
subplot(122); plot3(2*[-STAs(2,maxT-whichframe2+1) STAs(2,maxT-whichframe2+1)],2*[-STAs(4,maxT-whichframe2+1) STAs(4,maxT-whichframe2+1)],2*[-STAs(6,maxT-whichframe2+1) STAs(6,maxT-whichframe2+1)],'-k','Linewidth',4);
xlabel('R'), ylabel('G'); zlabel('B'); title('Subunit 2'); axis equal; hold off;
% subplot(223); plot3([0 0],[0 0],[-3 3],'-k','Linewidth',3); plot3([0 0],[-3 3],[0 0],'-k','Linewidth',3); plot3([-3 3],[0 0],[0 0],'-k','Linewidth',3); xlabel('L'); ylabel('M'); zlabel('S'); axis equal; hold off;
% subplot(224); plot3([0 0],[0 0],[-3 3],'-k','Linewidth',3); plot3([0 0],[-3 3],[0 0],'-k','Linewidth',3); plot3([-3 3],[0 0],[0 0],'-k','Linewidth',3); xlabel('L'); ylabel('M'); zlabel('S'); axis equal; hold off;
plot_counter = plot_counter + 2;

%************************************Fitting and testing GLMs*****************************
% The next step would be to fit a logistic regression model
% Regression
responseS1(responseS1==-1)=0;
responseS2(responseS2==-1)=0;
mdl1lin =  fitglm(rawrgb_onbasisvec1',responseS1','linear','Distribution','binomial','Link','logit');
mdl1quad =  fitglm(rawrgb_onbasisvec1',responseS1','quadratic','Distribution','binomial','Link','logit');
mdl2lin =  fitglm(rawrgb_onbasisvec2',responseS2','linear','Distribution','binomial','Link','logit');
mdl2quad =  fitglm(rawrgb_onbasisvec2',responseS2','quadratic','Distribution','binomial','Link','logit');
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

% Here I want to try some model comparisons techniques Ftest ANOVA
Fstatistic1num = (SSElin1 - SSEquad1)/(mdl1quad.NumCoefficients - mdl1lin.NumCoefficients);
Fstatistic1den = SSEquad1/(numel(responseS1) - mdl1quad.NumCoefficients);
Fstatistic2num = (SSElin2 - SSEquad2)/(mdl2quad.NumCoefficients - mdl2lin.NumCoefficients);
Fstatistic2den = SSEquad2/(numel(responseS1) - mdl2quad.NumCoefficients);
Fprob1 = 1-fcdf(Fstatistic1num/Fstatistic1den,mdl1quad.NumCoefficients - mdl1lin.NumCoefficients,numel(responseS1) - mdl1quad.NumCoefficients);
Fprob2 = 1-fcdf(Fstatistic2num/Fstatistic2den,mdl2quad.NumCoefficients - mdl2lin.NumCoefficients,numel(responseS2) - mdl2quad.NumCoefficients);

% Trying to compare the prediction differences between the GLM and GQM
diffmodelpred1 = sum((predlin1-predquad1).^2)/numel(predlin1);
diffmodelpred2 = sum((predlin2-predquad2).^2)/numel(predlin2);

% Log-likelihood ratio test (Wilk's theorem)
LLRprob1 = 1-chi2cdf(-2*(mdl1lin.LogLikelihood - mdl1quad.LogLikelihood),mdl1quad.NumCoefficients - mdl1lin.NumCoefficients);
LLRprob2 = 1-chi2cdf(-2*(mdl2lin.LogLikelihood - mdl2quad.LogLikelihood),mdl2quad.NumCoefficients - mdl2lin.NumCoefficients);

%Deviance test
Devprob1 = 1-chi2cdf(mdl1lin.Deviance-mdl1quad.Deviance,mdl1quad.NumCoefficients - mdl1lin.NumCoefficients);
Devprob2 = 1-chi2cdf(mdl2lin.Deviance-mdl2quad.Deviance,mdl2quad.NumCoefficients - mdl2lin.NumCoefficients);

% Storing the number of coefficients of quadratic terms which have a p-value < 0.001
sigquadterms1 = sum(mdl1quad.Coefficients.pValue(4:end)<0.001);
sigquadterms2 = sum(mdl2quad.Coefficients.pValue(4:end)<0.001);

% Calculating AIC
AIClin1 = 2*mdl1lin.NumCoefficients - 2*mdl1lin.LogLikelihood;
AICquad1 = 2*mdl1quad.NumCoefficients - 2*mdl1quad.LogLikelihood;
AIClin2 = 2*mdl2lin.NumCoefficients - 2*mdl2lin.LogLikelihood;
AICquad2 = 2*mdl2quad.NumCoefficients - 2*mdl2quad.LogLikelihood;

% Calculating BIC
BIClin1 = mdl1lin.NumCoefficients*log(numel(predlin1)) - 2*mdl1lin.LogLikelihood;
BICquad1 = mdl1quad.NumCoefficients*log(numel(predquad1)) - 2*mdl1quad.LogLikelihood;
BIClin2 = mdl2lin.NumCoefficients*log(numel(predlin2)) - 2*mdl2lin.LogLikelihood;
BICquad2 = mdl2quad.NumCoefficients*log(numel(predquad2)) - 2*mdl2quad.LogLikelihood;

fprintf('GLM and GQM results................................\n');
fprintf('Subunit 1 results................................\n');
fprintf('F-test p=%1.4f\n',Fprob1);
fprintf('Square of prediction differences between GLM and GQM %2.7f\n',diffmodelpred1);
fprintf('LLR test p=%1.4f\n',LLRprob1);
fprintf('Deviance test p=%1.4f\n',Devprob1);
fprintf('Number of significant quadratic terms %d\n',sigquadterms1);
fprintf('AIC lin model %f\n',AIClin1);
fprintf('AIC quad model %f\n',AICquad1);
fprintf('BIC lin model %f\n',BIClin1);
fprintf('BIC quad model %f\n\n',BICquad1);

fprintf('Subunit 2 results................................\n');
fprintf('F-test p=%1.4f\n',Fprob2);
fprintf('Square of prediction differences between GLM and GQM %2.7f\n',diffmodelpred2);
fprintf('LLR test p=%1.4f\n',LLRprob2);
fprintf('Deviance test p=%1.4f\n',Devprob2);
fprintf('Number of significant quadratic terms %d\n',sigquadterms2);
fprintf('AIC lin model %f\n',AIClin2);
fprintf('AIC quad model %f\n',AICquad2);
fprintf('BIC lin model %f\n',BIClin2);
fprintf('BIC quad model %f\n\n',BICquad2);
 
% Classification instead of regression
LDA1 = classify(rawrgb_onbasisvec1',rawrgb_onbasisvec1',responseS1','linear');
QDA1 = classify(rawrgb_onbasisvec1',rawrgb_onbasisvec1',responseS1','quadratic');
LDA2 = classify(rawrgb_onbasisvec2',rawrgb_onbasisvec2',responseS2','linear');
QDA2 = classify(rawrgb_onbasisvec2',rawrgb_onbasisvec2',responseS2','quadratic');

% Calculating percentage correct classification from LDA and QDA
LDA1perf = sum(xor(LDA1,responseS1'))/numel(responseS1);
QDA1perf = sum(xor(QDA1,responseS1'))/numel(responseS1);
LDA2perf = sum(xor(LDA2,responseS2'))/numel(responseS2);
QDA2perf = sum(xor(QDA2,responseS2'))/numel(responseS2);

% Calculating Area under curve (AUC) from ROC
AUClin1 = calcAUC(LDA1,responseS1');
AUCquad1 = calcAUC(QDA1,responseS1');
AUClin2 = calcAUC(LDA2,responseS2');
AUCquad2 = calcAUC(QDA2,responseS2');

fprintf('LDA and QDA results................................\n');
fprintf('Subunit 1 results................................\n');
fprintf('LDA %1.3f\n',LDA1perf);
fprintf('QDA %1.3f\n',QDA1perf);
fprintf('AUC lin %1.3f\n',AUClin1);
fprintf('AUC quad %1.3f\n\n',AUCquad1);

fprintf('Subunit 2 results................................\n');
fprintf('LDA %1.3f\n',LDA2perf);
fprintf('QDA %1.3f\n',QDA2perf);
fprintf('AUC lin %1.3f\n',AUClin2);
fprintf('AUC quad %1.3f\n\n',AUCquad2);

%%
V1lin = reshape(predict(mdl1lin,[X(:) Y(:) Z(:)]),size(X));
V1quad = reshape(predict(mdl1quad,[X(:) Y(:) Z(:)]),size(X));
V2lin = reshape(predict(mdl2lin,[X(:) Y(:) Z(:)]),size(X));
V2quad = reshape(predict(mdl2quad,[X(:) Y(:) Z(:)]),size(X));
figure(plot_counter); set(gcf,'Name','GLMfits');
subplot(221),isosurface(X,Y,Z,V1lin,prctile(V1lin(:),50)); hold on; isosurface(X,Y,Z,V1lin,prctile(V1lin(:),10)); isosurface(X,Y,Z,V1lin,prctile(V1lin(:),90));
plot3(2*[-STAs(1,maxT-whichframe1+1) STAs(1,maxT-whichframe1+1)],2*[-STAs(3,maxT-whichframe1+1) STAs(3,maxT-whichframe1+1)],2*[-STAs(5,maxT-whichframe1+1) STAs(5,maxT-whichframe1+1)],'-k','Linewidth',4);
xlabel('R'); ylabel('G'); zlabel('B'); title('S1 linear'); hold off; axis equal;
subplot(222),isosurface(X,Y,Z,V1quad,prctile(V1quad(:),50)); hold on; isosurface(X,Y,Z,V1quad,prctile(V1quad(:),10)); isosurface(X,Y,Z,V1quad,prctile(V1quad(:),90));
plot3(2*[-STAs(1,maxT-whichframe1+1) STAs(1,maxT-whichframe1+1)],2*[-STAs(3,maxT-whichframe1+1) STAs(3,maxT-whichframe1+1)],2*[-STAs(5,maxT-whichframe1+1) STAs(5,maxT-whichframe1+1)],'-k','Linewidth',4);
xlabel('R'); ylabel('G'); zlabel('B'); title('S1 quad'); hold off; axis equal;
subplot(223),isosurface(X,Y,Z,V2lin,prctile(V2lin(:),50)); hold on; isosurface(X,Y,Z,V2lin,prctile(V2lin(:),10)); isosurface(X,Y,Z,V2lin,prctile(V2lin(:),90));
plot3(2*[-STAs(2,maxT-whichframe2+1) STAs(2,maxT-whichframe2+1)],2*[-STAs(4,maxT-whichframe2+1) STAs(4,maxT-whichframe2+1)],2*[-STAs(6,maxT-whichframe2+1) STAs(6,maxT-whichframe2+1)],'-k','Linewidth',4);
xlabel('R'); ylabel('G'); zlabel('B'); title('S2 linear'); hold off; axis equal;
subplot(224),isosurface(X,Y,Z,V2quad,prctile(V2quad(:),50)); hold on; isosurface(X,Y,Z,V2quad,prctile(V2quad(:),10)); isosurface(X,Y,Z,V2quad,prctile(V2quad(:),90));
plot3(2*[-STAs(2,maxT-whichframe2+1) STAs(2,maxT-whichframe2+1)],2*[-STAs(4,maxT-whichframe2+1) STAs(4,maxT-whichframe2+1)],2*[-STAs(6,maxT-whichframe2+1) STAs(6,maxT-whichframe2+1)],'-k','Linewidth',4);
xlabel('R'); ylabel('G'); zlabel('B'); title('S2 quad'); hold off; axis equal;
plot_counter = plot_counter + 1;



%% I am trying something based on Fred's and Greg's suggestion after my thesis committee meeting
probsS1 = mvnpdf(rawrgb_onbasisvec1',muS1spike',covS1spike)./mvnpdf(rawrgb_onbasisvec1',[0 0 0],covS1);
probsS2 = mvnpdf(rawrgb_onbasisvec2',muS2spike',covS2spike)./mvnpdf(rawrgb_onbasisvec2',[0 0 0],covS2);
binsx = linspace(prctile(probsS1,5),prctile(probsS1,95),15);
binsy = linspace(prctile(probsS2,5),prctile(probsS2,95),15);
nraw = hist3([probsS1 probsS2],{binsx,binsy}); nraw(nraw==0) = 1;
nspike = hist3([probsS1(responseS1==1) probsS2(responseS2==1)],{binsx,binsy});

% Now maybe try projecting all the stimuli onto the STA of the subunits
projontoS1raw = muS1spike'* rawrgb_onbasisvec1;
projontoS1spike = muS1spike'* rawrgb_onbasisvec1(:,responseS1==1);
projontoS2raw = muS2spike'* rawrgb_onbasisvec2;
projontoS2spike = muS2spike'* rawrgb_onbasisvec2(:,responseS2==1);
covprojraw = cov([projontoS1raw; projontoS2raw]'); muraw = [mean(projontoS1raw) mean(projontoS2raw)];
covprojspike = cov([projontoS1spike; projontoS2spike]'); muspike = [mean(projontoS1spike) mean(projontoS2spike)];
[binsX1 binsY1] = meshgrid(linspace(min(projontoS1raw),max(projontoS1raw),51),linspace(min(projontoS2raw),max(projontoS2raw),51));
linresp = mvnpdf([binsX1(:) binsY1(:)],muspike,covprojspike)./mvnpdf([binsX1(:) binsY1(:)],muraw,covprojraw);
linresp = reshape(linresp,size(binsX1));

figure(plot_counter), subplot(221); plot(probsS1,probsS2,'.'); hold on; plot(probsS1(responseS1==1),probsS2(responseS2==1),'r.'); axis square; xlabel('prob S1'); ylabel('prob S2'); hold off;
subplot(222);imagesc(binsx,binsy,nspike./nraw); axis square; axis xy; xlabel('prob S2'); ylabel('prob S1');
subplot(223); plot(projontoS1raw,projontoS2raw,'.'); hold on; plot(projontoS1spike,projontoS2spike,'r.');
axis square; xlabel('proj S1'); ylabel('proj S2'); title('Linearity assumption');hold off;
subplot(224); contour(binsX1,binsY1,linresp); axis square; xlabel('proj S1'); ylabel('proj S2'); title('Linearity assumption');hold off;
plot_counter = plot_counter + 1;

% Using inbuilt MATLAB kernel density estimator
figure(plot_counter);subplot(121),ksdensity([probsS1 probsS2]); hold on; ksdensity([probsS1(responseS1==1) probsS2(responseS2==1)]); hold off;