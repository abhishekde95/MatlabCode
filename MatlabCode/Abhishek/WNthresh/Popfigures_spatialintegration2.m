% A new script for generating figures related to the spatial integration part
% Greg's suggestion was to divide the story into 2 papers,
% after Script for generating figures for the spatial integration paper
% Author - Abhishek De, 6/9/2020
close all; clearvars;
plot_counter = 1;

%% Figure 1: Impact of linear and non-linear spatial filtering in image processing 
% This code has been derived from Abhishek/Physiology_modeling/Edge_processing_linear_nonlinear

if ~exist('plot_counter')
    plot_counter = 1;
end

close all; clearvars;
I_color = imread('peppers.png');

% Designing a filter
F = zeros(4,2); S1 = F; S2 = F;
S1(:,1) = 1; S2(:,2) = -1;
F = S1+S2;

% Forming the M matrix 
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

row  = size(I_color,1); col = size(I_color,2);
I_cone_excitation = M*RGB2XWFormat(im2double(I_color))';
I_cone_excitation = XW2RGBFormat(I_cone_excitation',row,col);
I_Lcone_excitation = I_cone_excitation(:,:,1);
I_Mcone_excitation = I_cone_excitation(:,:,2);
I_Scone_excitation = I_cone_excitation(:,:,3);

I_Lcone = (I_Lcone_excitation-median(I_Lcone_excitation(:)))/(median(I_Lcone_excitation(:)));
I_Mcone = (I_Mcone_excitation-median(I_Mcone_excitation(:)))/(median(I_Mcone_excitation(:)));
I_Scone = (I_Scone_excitation-median(I_Scone_excitation(:)))/(median(I_Scone_excitation(:)));
I_LminusM = I_Lcone - I_Mcone;  
I_MminusL = I_Mcone - I_Lcone; 

% Linear filtering
I_color_filtered = imfilter(I_LminusM,F,'conv');
I_color_filtered(I_color_filtered<0) = 0;
I_color_filtered = I_color_filtered./(1.5*max(I_color_filtered(:)));
I_color_filtered = I_color_filtered + 0.5-median(I_color_filtered(:));

% Nonlinear filtering 
I_LminusM(I_LminusM<0) = 0; I_MminusL(I_MminusL<0) = 0;
% I_color_filteredNL = I_LminusM.^2 + I_MminusL.^2; % Squaring non-linearity
I_color_filteredNL = I_LminusM.^2 + I_MminusL.^2; % The vanilla ReLu non-linearity
I_color_filteredNL(I_color_filteredNL<0) = 0;
I_color_filteredNL = I_color_filteredNL./(1.5*max(I_color_filteredNL(:)));
I_color_filteredNL = I_color_filteredNL + 0.5-median(I_color_filteredNL(:));

% Visualizing the filter
A = 0.5*ones(100,100,3);
RGB_forLminusM = inv(M')*[1;-1;0]; RGB_forLminusM = RGB_forLminusM./(max(abs(RGB_forLminusM)));
RGB_forMminusL = inv(M')*[-1;1;0]; RGB_forMminusL = RGB_forMminusL./(max(abs(RGB_forMminusL)));
A(20:80,20:50,1) = RGB_forMminusL(1); A(20:80,20:50,2) = RGB_forMminusL(2); A(20:80,20:50,3) = RGB_forMminusL(3); % Subunit 1-> L-M
A(20:80,51:80,1) = RGB_forLminusM(1); A(20:80,51:80,2) = RGB_forLminusM(2); A(20:80,51:80,3) = RGB_forLminusM(3); % Subunit 2-> M-L

% Normalizing the L-M  and M-L responses
I_LminusM_norm = I_LminusM - median(I_LminusM(:));
I_LminusM_norm = I_LminusM_norm./(2*max(I_LminusM_norm(:)));
I_LminusM_norm = I_LminusM_norm + 0.5;

I_MminusL_norm = I_MminusL - median(I_MminusL(:));
I_MminusL_norm = I_MminusL./(2*max(I_MminusL(:)));
I_MminusL_norm = I_MminusL_norm + 0.5;

% Plotting the images
plot_counter = 1;
figure(plot_counter); set(gcf,'Name','Analyses of images: RGB');
subplot(321); imagesc(A); axis square; set(gca,'XTick',[],'YTick',[]); title('Filter');
subplot(322); imshow(I_color); axis square; set(gca,'XTick',[],'YTick',[]); title('original pic');
subplot(323); imshow(im2double(I_LminusM_norm)); axis square; set(gca,'XTick',[],'YTick',[]); colormap('gray'); title('Rectified L-M map');
subplot(324); imshow(im2double(I_MminusL_norm)); axis square; set(gca,'XTick',[],'YTick',[]); colormap('gray'); title('Rectified M-L map');
subplot(325); imshow(im2double(I_color_filtered)); axis square; set(gca,'XTick',[],'YTick',[]); colormap('gray'); hold off; title('Linear filter');
subplot(326); imshow(im2double(I_color_filteredNL)); axis square; set(gca,'XTick',[],'YTick',[]); colormap('gray'); hold off; title('Relu non-linearity');
plot_counter = plot_counter + 1;

% Trying out the non-rectified maps
I_LminusMnon = I_Lcone - I_Mcone;  
I_MminusLnon = I_Mcone - I_Lcone;
I_LmminusMnon = I_LminusMnon - median(I_LminusMnon(:));
I_LminusMnon = I_LminusMnon./(2*max(abs(I_LminusMnon(:))));
I_LminusMnon = I_LminusMnon + 0.5;

I_MmminusLnon = I_MminusLnon - median(I_MminusLnon(:));
I_MminusLnon = I_MminusLnon./(2*max(abs(I_MminusLnon(:))));
I_MminusLnon = I_MminusLnon + 0.5;

% PLotting the distribution of the L-, M- and S-cone contrast
figure(plot_counter);
subplot(221); histogram(I_Lcone_excitation(:),0:0.01:0.2,'DisplayStyle','stairs','EdgeColor',[1 0 0]); hold on;
histogram(I_Mcone_excitation(:),0:0.01:0.2,'DisplayStyle','stairs','EdgeColor',[0 1 0]); 
histogram(I_Scone_excitation(:),0:0.01:0.2,'DisplayStyle','stairs','EdgeColor',[0 0 1]); 
plot(median(I_Lcone_excitation(:)),0,'v','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(median(I_Mcone_excitation(:)),0,'v','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
plot(median(I_Scone_excitation(:)),0,'v','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out'); xlabel('cone excitation'); ylabel('Pixel count'); legend('L','M','S');
subplot(222); histogram(I_Lcone(:),-1:0.1:1,'DisplayStyle','stairs','EdgeColor',[1 0 0]); hold on;
histogram(I_Mcone(:),-1:0.1:1,'DisplayStyle','stairs','EdgeColor',[0 1 0]); 
histogram(I_Scone(:),-1:0.1:1,'DisplayStyle','stairs','EdgeColor',[0 0 1]); 
plot(median(I_Lcone(:)),0,'v','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(median(I_Mcone(:)),2500,'v','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
plot(median(I_Scone(:)),5000,'v','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out'); xlabel('cone contrast'); ylabel('Pixel count'); legend('L','M','S');
subplot(223); imshow(im2double(I_LminusMnon)); axis square; colormap('gray'); set(gca,'XTick',[],'YTick',[]); title('L-M map');
subplot(224); imshow(im2double(I_MminusLnon)); axis square; colormap('gray'); set(gca,'XTick',[],'YTick',[]); title('M-L map');

plot_counter = plot_counter + 1;
%% Figure 2-part1: Extracting color-time signal from the example DO cell and some masked random WN stimulus
% Along with its Spatial weighting function and STA

% Rest of the figure was adapted from the SFN and COSYNE poster

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

% Creating some masked WN subunit stimuli 
subunit_mask = reshape(WN.ras{2,maskidx},[nstixperside nstixperside]);
stim = zeros(size(subunit_mask));
subunit1_mask = repmat(subunit_mask==1,[1 1 3]);
subunit2_mask = repmat(subunit_mask==2, [1 1 3]);
bkgnd_mask = repmat(0.5*(subunit_mask==0),[1 1 3]);
figure(plot_counter);
for ii = 1:9
    S1 = cat(3,repmat(rand(1),[10 10]),repmat(rand(1),[10 10]),repmat(rand(1),[10 10]));
    S2 = cat(3,repmat(rand(1),[10 10]),repmat(rand(1),[10 10]),repmat(rand(1),[10 10]));
    subplot(3,3,ii); image(S1.*subunit1_mask +S2.*subunit2_mask + bkgnd_mask); axis square; set(gca,'XTick',[],'YTick',[]);
end


%% A continuation of Figure 2- Code for plotting the checkerboard STA along with the spatial weighting map

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

% Plotting the grayscale spatial weighting function
spatial_weighting = u(:,1);
tmpSTA3 = reshape(tmpSTA,[10 10 3]);
rawenergymap = sum(reshape(tmpSTA,[10 10 3]).^2,3);
energymap = imgaussfilt(rawenergymap,1.0);

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

figure(plot_counter);
subplot(121); image(im_new); axis square; set(gca,'XTick',[],'YTick',[]); title('Same size STA')
subplot(122); imagesc(weighted_sim); axis square; set(gca,'XTick',[],'YTick',[]); title('Same size spatial weighting function'); colormap('gray');
plot_counter = plot_counter + 1;

%% Figure 2-part 2: An exmaple of firing rate map calculated using the WN subunit analysis
% The code is in Calculating_effect_contrast.m
% The code will plot the firing rate surface at a rotated angle slong with
% the isoresponse contours from a GLM fit


%% Figure 2-part3: Determining the spatial interaction between the subunits using WN projection analysis- an analysis that could motivate why isoresponse is needed in the first place 

if ~exist('plot_counter')
    plot_counter = 1;
end

load conewts_svd.mat
load vals.mat
load S1RGB.mat
load S2RGB.mat

% Calculating some metrics for circular stats and discarding SO cells
unitvectors = [];
anglebwvectors = [];
for ii = 1:size(S1RGB,2)
    vec1 = S1RGB(:,ii);
    vec2 = S2RGB(:,ii);
    vec = exp(atan2(norm(cross(vec1,vec2)),dot(vec1,vec2))*i);
    
    % Storing the angle bewteen the vectors
    anglebwvectors = [anglebwvectors; atan2(norm(cross(vec1,vec2)),dot(vec1,vec2))*180/pi];
    
    % Converting the angle difference into unit vectors
    unitvectors = [unitvectors; real(vec) imag(vec)];
end

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

% Considering only the spatially opponent subunits

% Load the integration within the subunit whitenoise analysis data
load AUROClinsubunits_CV.mat
load AUROCquadsubunits_CV.mat

% Load the isoresponse data
load RSSE_linearmodel_CV.mat % Robust regression
load RSSE_quadmodel_CV.mat


% For storing the Isoresponse NLI
Isoresponse_NLI = [];

% For storing the white noise NLI
Whitenoise_NLI = [];

for ii = 1:numel(AUROClinsubunits) 
    
    % White noise NLI 
    Error_quad = 1-(AUROCquadsubunits{ii});
    Error_lin = 1-(AUROClinsubunits{ii});
    Whitenoise_NLI = [Whitenoise_NLI; log10(median(Error_lin./Error_quad))];
    
    % Isoresponse NLI
    Isoresponse_NLI = [Isoresponse_NLI; log10(median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii}))];
end

indices = [109 24 74];

% Plotting the white noise NLIs on a cell-type basis
figure(plot_counter);
subplot(311); histogram(Whitenoise_NLI(DOidx),-0.02:0.005:0.1,'FaceColor',[1 0 0],'EdgeColor',[1 1 1]); hold on;
plot(median(Whitenoise_NLI(DOidx)),15,'v','MarkerSize',8,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(Whitenoise_NLI(indices(2)),18,'s','MarkerSize',8,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-0.02 0.08],'XTick',-0.02:0.02:0.08,'Ylim',[0 20],'YTick',[0 10 20]); xlabel('Whitenoise NLI'); ylabel('Count'); title('Across subunits'); axis square; hold off;
subplot(312); histogram(Whitenoise_NLI(LUMidx),-0.02:0.005:0.1,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); hold on;
plot(median(Whitenoise_NLI(LUMidx)),15,'v','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
plot(Whitenoise_NLI(indices(1)),18,'s','MarkerSize',8,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-0.02 0.08],'XTick',-0.02:0.02:0.08,'Ylim',[0 20],'YTick',[0 10 20]); xlabel('Whitenoise NLI'); ylabel('Count'); title('Across subunits'); axis square; hold off;
subplot(313); histogram(Whitenoise_NLI(hardtoclassifyidx),-0.02:0.005:0.1,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1]); hold on;
plot(median(Whitenoise_NLI(hardtoclassifyidx)),20,'v','MarkerSize',8,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
plot(Whitenoise_NLI(indices(3)),20,'s','MarkerSize',8,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-0.02 0.08],'XTick',-0.02:0.02:0.08,'Ylim',[0 20],'YTick',[0 10 20]); xlabel('Whitenoise NLI'); ylabel('Count'); title('Across subunits'); axis square; hold off;
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;

% Comparing spatial NLI across cell types
group = [ones(size(LUMidx)) 2*ones(size(DOidx)) 3*ones(size(hardtoclassifyidx))];
data = Whitenoise_NLI([LUMidx'; DOidx'; hardtoclassifyidx']); 
p1 = kruskalwallis(data,group','off');

[r2,p2] = corr(Isoresponse_NLI,Whitenoise_NLI,'type','Spearman');

% DO + Simple cells vs. hardtoclassify cells
[p3,h] = ranksum(Whitenoise_NLI([DOidx LUMidx]),Whitenoise_NLI(hardtoclassifyidx));


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

%% Figure 3-part1: Iso-response example data from example LUM, DO and HTC cells 
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
indices = [109 24 74];% example LUM, DO and HTC cells
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
    figure(2); subplot(1,numel(indices),zz);
    for ii = 1:numel(not_oog_idx)
        h(ii) = plot(x_orig(not_oog_idx(ii)), y_orig(not_oog_idx(ii)),'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
        set(h(ii),'ButtonDownFcn',{@dispimage,x_orig(not_oog_idx(ii)),y_orig(not_oog_idx(ii)),basisvec,bkgnd_monitor});%(x_orig(not_oog_idx(ii)),y_orig(not_oog_idx(ii)),basisvec,bkgnd_monitor)});
    end
    plot(x_lin,y_lin,'g','Linewidth',2); hold on; plot(x_quad2,y_quad2,'r','Linewidth',2);
    xlabel('X'), ylabel('Y');
    if zz == 1
        title('Simple');
    elseif zz==2
        title('DO');
    else
        title('Unclassified')
    end
    
    if any(outofgamut)
        plot(upsample(x_orig(outofgamut),2),upsample(y_orig(outofgamut),2),'color',[0.5 0.5 0.5]);
    end
    set(gca,'Tickdir','out'); drawnow; axis square; 
    if zz==1 
        set(gca,'Xlim',[-1.6 1.6],'Ylim',[-1.6 1.6],'YTick',-1.6:0.8:1.6,'XTick',[-1.6:0.8:1.6]); hold off;
    elseif zz==2
        set(gca,'Xlim',[-2.0 2.0],'Ylim',[-2.0 2.0],'YTick',-2.0:1.0:2.0,'XTick',[-2.0:1.0:2.0]); hold off;
    elseif zz==3
        set(gca,'Xlim',[-0.2 0.2],'Ylim',[-0.2 0.2],'YTick',-0.2:0.1:0.2,'XTick',[-0.2:0.1:0.2]); hold off;
    end
    
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
    if zz == 1
        title('Simple');
    elseif zz==2
        title('DO');
    else
        title('Unclassified')
    end
    axis square; hold off;
   
    
    plot_counter = plot_counter + 1;
end

% Obtaining the target firing rates of the example neurons
load TFR.mat
FR = TFR(1,indices);


%% Figure 3-part 2: Population analysis of isoresponse curves 
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
indices = [109 24 74];
% Plotting the results for SVD based cone weight classification including the PC1 z-scores 
figure(plot_counter);
subplot(311); histogram(RSSEisoresp_medianofratios(LUMidx),logspace(-1,3,31),'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); hold on;
plot(median(RSSEisoresp_medianofratios(LUMidx)),15,'v','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
plot(RSSEisoresp_medianofratios(indices(1)),14,'s','MarkerSize',8,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 100],'XTick',[0.1 1 10 100],'Ylim',[0 15],'YTick',[0 5 10 15]); ylabel('Count'); title('Isoresponse'); xlabel('median CV ratio of errors'); axis square; hold off;
subplot(312); histogram(RSSEisoresp_medianofratios(DOidx),logspace(-1,3,31),'FaceColor',[1 0 0],'EdgeColor',[1 1 1]); hold on;
plot(median(RSSEisoresp_medianofratios(DOidx)),10,'v','MarkerSize',8,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(RSSEisoresp_medianofratios(indices(2)),9,'s','MarkerSize',8,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 100],'XTick',[0.1 1 10 100],'Ylim',[0 10],'YTick',[0 5 10]); ylabel('Count'); title('Isoresponse'); xlabel('median CV ratio of errors'); axis square; hold off;
subplot(313); histogram(RSSEisoresp_medianofratios(hardtoclassifyidx),logspace(-1,3,31),'FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1]); hold on;
plot(median(RSSEisoresp_medianofratios(hardtoclassifyidx)),15,'v','MarkerSize',8,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); 
plot(RSSEisoresp_medianofratios(indices(3)),14,'s','MarkerSize',8,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 100],'XTick',[0.1 1 10 100],'Ylim',[0 15],'YTick',[0 5 10 15]); ylabel('Count'); title('Isoresponse'); xlabel('median CV ratio of errors'); axis square; hold off;
plot_counter = plot_counter + 1;

% Comparing spatial NLI across cell types
group = [ones(size(LUMidx)) 2*ones(size(DOidx)) 3*ones(size(hardtoclassifyidx))];
data = RSSEisoresp_medianofratios([LUMidx'; DOidx'; hardtoclassifyidx']); 
p1 = kruskalwallis(data,group,'off');

% Comparing spatial NLI between simple and DO cells
[p2,h] = ranksum(RSSEisoresp_medianofratios(LUMidx),RSSEisoresp_medianofratios(DOidx));

% Comparing spatial NLI between simple + DO cells and other cells
[p5,h] = ranksum(RSSEisoresp_medianofratios([LUMidx DOidx]),RSSEisoresp_medianofratios(hardtoclassifyidx));

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

% Checking how target firing rates of the DO, simple and unclassified cells

load baselineFRstats.mat % Baseline FR rates
load TFR.mat % Target firing rates

figure(plot_counter);
indices = [DOidx LUMidx hardtoclassifyidx];
TFRprctile = [];
for iter = 1:numel(baselineFRstats)
    ii = indices(iter);
    baselineFR = baselineFRstats{ii};
    if ismember(ii,DOidx)
        c = [1 0 0];
    elseif ismember(ii, LUMidx)
        c = [0 0 0];
    else
        c = [0.5 0.5 0.5];
    end
    
    plot([iter iter],[prctile(baselineFR,5) prctile(baselineFR,95)],'color',c); hold on 
    
    TFRprctile = [TFRprctile; invprctile(baselineFR, TFR(1,ii))];
    for jj = 1:2
        if TFR(jj,ii)>0
            plot(iter, TFR(jj,ii), 'o', 'MarkerFaceColor',c,'MarkerEdgeColor',[1 1 1]);
        end
    end
    
end
axis square; set(gca,'Tickdir','out', 'Xlim',[1 numel(baselineFRstats)], 'Ylim', [0 120], 'YTick', 0:30:120);
xlabel('Cell number'); ylabel('Firing rate'); 
set(gcf,'renderer','painters'); hold off

% Doing some stats on the target firing rates 
TFR_DO = TFR(:,DOidx);
TFR_DO = TFR_DO(:); TFR_DO = TFR_DO(TFR_DO>0);

TFR_LUM = TFR(:,LUMidx);
TFR_LUM = TFR_LUM(:); TFR_LUM = TFR_LUM(TFR_LUM>0);

TFR_htc = TFR(:, hardtoclassifyidx);
TFR_htc = TFR_htc(:); TFR_htc = TFR_htc(TFR_htc>0);

group = [ones(size(TFR_DO)); 2*ones(size(TFR_LUM)); 3*ones(size(TFR_htc))];
data = [TFR_DO; TFR_LUM; TFR_htc]; 
p2 = kruskalwallis(data,group,'off');


% Could for used for some presentation purposes 
figure(plot_counter)
histogram(RSSEisoresp_medianofratios(LUMidx),logspace(-1,3,31),'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); hold on;
plot(median(RSSEisoresp_medianofratios(LUMidx)),15,'v','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
plot(RSSEisoresp_medianofratios(indices(1)),14,'s','MarkerSize',8,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
histogram(RSSEisoresp_medianofratios(DOidx),logspace(-1,3,31),'FaceColor',[1 0 0],'EdgeColor',[1 1 1]); hold on;
plot(median(RSSEisoresp_medianofratios(DOidx)),10,'v','MarkerSize',8,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(RSSEisoresp_medianofratios(indices(2)),9,'s','MarkerSize',8,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
histogram(RSSEisoresp_medianofratios(hardtoclassifyidx),logspace(-1,3,31),'FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1]); hold on;
plot(median(RSSEisoresp_medianofratios(hardtoclassifyidx)),15,'v','MarkerSize',8,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); 
plot(RSSEisoresp_medianofratios(indices(3)),14,'s','MarkerSize',8,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 100],'XTick',[0.1 1 10 100],'Ylim',[0 15],'YTick',[0 5 10 15]); ylabel('Count'); title('Isoresponse'); xlabel('median CV ratio of errors'); axis square; hold off;
plot_counter = plot_counter + 1;

%% Figure 3 - final part: plotting the relationship between the whitenoise and isoresponse NLI 

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

% Checking the correlation with non-linearity indices 
% Load the isoresponse data
load RSSE_linearmodel_CV.mat % Robust regression
load RSSE_quadmodel_CV.mat
% Load the integration within the subunit data
load AUROClinsubunits_CV.mat
load AUROCquadsubunits_CV.mat

% For storing the Isoresponse NLI
Isoresponse_NLI = [];

% For storing the white noise NLI
Whitenoise_NLI = [];

indices = [109 24 74];
for ii = 1:numel(AUROClinsubunits) 
  
    % White noise NLI 
    Error_quad = 1-(AUROCquadsubunits{ii});
    Error_lin = 1-(AUROClinsubunits{ii});
    Whitenoise_NLI = [Whitenoise_NLI; log10(median(Error_lin./Error_quad))];
    
    % Isoresponse NLI
    Isoresponse_NLI = [Isoresponse_NLI; log10(median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii}))];
end

figure(plot_counter); hold on;
plot(Isoresponse_NLI(LUMidx), Whitenoise_NLI(LUMidx), 'o', 'MarkerFaceColor', [0 0 0],'MarkerEdgeColor',[1 1 1]);
plot(Isoresponse_NLI(indices(1)), Whitenoise_NLI(indices(1)), 'o', 'MarkerFaceColor', [0 0 0],'MarkerEdgeColor',[0 1 0]);
plot(Isoresponse_NLI(DOidx), Whitenoise_NLI(DOidx), 'o', 'MarkerFaceColor', [1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(Isoresponse_NLI(indices(2)), Whitenoise_NLI(indices(2)), 'o', 'MarkerFaceColor', [1 0 0],'MarkerEdgeColor',[0 1 0]);
plot(Isoresponse_NLI(hardtoclassifyidx), Whitenoise_NLI(hardtoclassifyidx), 'o', 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
plot(Isoresponse_NLI(indices(3)), Whitenoise_NLI(indices(3)), 'o', 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor',[0 1 0]);
set(gca,'Tickdir','out','Xlim',[-1 2],'XTick',[-1 0 1 2],'Ylim',[-0.02 0.08],'YTick',[-0.02:0.02:0.08]); 
xlabel('Isoresponse NLI'); ylabel('WhiteNoise NLI'); axis square; hold off;
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;

% Some basic stats
[r1,p1] = corr(Isoresponse_NLI(LUMidx),Whitenoise_NLI(LUMidx),'type','Spearman');
[r2,p2] = corr(Isoresponse_NLI(DOidx),Whitenoise_NLI(DOidx),'type','Spearman');
[r3,p3] = corr(Isoresponse_NLI(hardtoclassifyidx),Whitenoise_NLI(hardtoclassifyidx),'type','Spearman');
[rc,pc] = corr(Isoresponse_NLI,Whitenoise_NLI,'type','Spearman');

%% Figure 4: Downstream circuitry of simple and DO cells
% No code for this part
% The figure was provided by Greg

%% Figure S1: Stimulus comparison and motivation for the 3 phases of the experiment

% Check the code in Calculating_effect_contrast.m file for stimulus
% comparison. The figure was generated using that code. 

%% Figure S2: An illustration of GQM and GLM analysis 

%% Figure S3: Adaptive staircase procedure from an example cell
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

%
figure(plot_counter);
cum_flag = []; % to check if I am looking at all the probed directions

t_offset1 = -0.2;
figure(plot_counter);
dir = 24; % Good directions: 14, 24, 25, 26
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
    set(gca,'Tickdir','out','YScale','log','Ylim',[0.1 1],'YTick',[0.1 0.3 1]); axis square;
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
    set(gca,'Tickdir','out','Ylim',[0 120],'YTick',[0 30 60 90 120]); axis square; xlabel('Trials'), ylabel('FR');
end


subplot(2,2,4); hold on; set(gca,'Xlim',[0.1 1],'XTick',[0.1 0.3 1],'Ylim',[0 120],'YTick',[0 30 60 90 120]);
subplot(2,2,3); hold on; plot([0 20],[30 30],'k'); hold off;

plot_counter = plot_counter + 1;

%% Figure S4
% Robustness of the isoresponse curves based on the last point
% Plotting the average of the last 5 points as well for each probed direction: checking for the example cells only

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
indices = [109 24 74];% example LUM, DO and HTC cells
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
    
    
    Fx = @(xi) any(isnan(xi)); % function that finds 'NaN' in a cell array
    inds = find(cellfun(Fx,stro.ras(:,basisvecidx))==0);
    if isempty(inds)
        inds = size(stro.trial,1)-1;
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
    
    
    % Plotting the staircase termination and the geometric mean of last 3 staircase
    % termination points
    completed_dir = completed_search_alongdir{1};
    probed_dirs = logical(completed_dir(:,1)==1); % only including the directions that have been completely probed
    oog_idx = completed_dir(:,1)==1 & completed_dir(:,2)==1; % probed and out of gamut
    not_oog_idx = completed_dir(:,1)==1 & completed_dir(:,2)==0;
    figure(1);
    for ii = 1:size(norms{1},1)
        dir = ii;
        for kk = 1: numel(num_targetspikerates)
            tmp_n = [];
            tmp_wts = [];
            tmp_parentvertices = [];
            idxs1 = find(stro.trial(:,basisvecdiridx) == dir);
            if ~isempty(idxs1)
                for jj = 1:numel(idxs1)
                    tmp_parentvertices = [tmp_parentvertices; stro.ras{idxs1(jj),parentverticesidx}];
                    tmp_wts = [tmp_wts; stro.ras{idxs1(jj),weightsidx}'];
                    tmp_n = [tmp_n; norm(stro.ras{idxs1(jj),weightsidx})];
                end
                raster_data = stro.ras(idxs1,1);
                num_dur =[];
                firing_rate = [];
                if (~isempty(tmp_n) & not_oog_idx(ii)) 
                    figure(1), 
                    subplot(3,2,2*zz-1); plot(sign(tmp_wts(end,1))*geomean(abs(tmp_wts(end-4:end,1))),sign(tmp_wts(end,2))*geomean(abs(tmp_wts(end-4:end,2))),'o','MarkerFacecolor',[0 0 0], 'MarkerEdgeColor',[1 1 1]); hold on;
                    subplot(3,2,2*zz); plot(tmp_wts(end,1),tmp_wts(end,2),'o','MarkerFacecolor',[0 0 0], 'MarkerEdgeColor',[1 1 1]); hold on;     
                elseif (~isempty(tmp_n) & oog_idx(ii))
                    figure(1), 
                    subplot(3,2,2*zz-1); plot([0 tmp_wts(end,1)],[0 tmp_wts(end,2)],'k'); hold on;
                    subplot(3,2,2*zz); plot([0 tmp_wts(end,1)],[0 tmp_wts(end,2)],'k'); hold on;     
                end
            end
        end
    end
 
   
    xlabel('X'), ylabel('Y');
    if zz == 1
        title('Simple');
    elseif zz==2
        title('DO');
    else
        title('Unclassified')
    end
    
     
    if zz==1 
        subplot(3,2,2*zz-1); set(gca,'Xlim',[-1.6 1.6],'Ylim',[-1.6 1.6],'YTick',-1.6:0.8:1.6,'XTick',[-1.6:0.8:1.6]); set(gca,'Tickdir','out'); drawnow; axis square; hold off;
        subplot(3,2,2*zz); set(gca,'Xlim',[-1.6 1.6],'Ylim',[-1.6 1.6],'YTick',-1.6:0.8:1.6,'XTick',[-1.6:0.8:1.6]);set(gca,'Tickdir','out'); drawnow; axis square; hold off;
    elseif zz==2
        subplot(3,2,2*zz-1); set(gca,'Xlim',[-2.0 2.0],'Ylim',[-2.0 2.0],'YTick',-2.0:1.0:2.0,'XTick',[-2.0:1.0:2.0]); set(gca,'Tickdir','out'); drawnow; axis square; hold off;
        subplot(3,2,2*zz); set(gca,'Xlim',[-2.0 2.0],'Ylim',[-2.0 2.0],'YTick',-2.0:1.0:2.0,'XTick',[-2.0:1.0:2.0]); set(gca,'Tickdir','out'); drawnow; axis square; hold off;
    elseif zz==3
        subplot(3,2,2*zz-1); set(gca,'Xlim',[-0.2 0.2],'Ylim',[-0.2 0.2],'YTick',-0.2:0.1:0.2,'XTick',[-0.2:0.1:0.2]); set(gca,'Tickdir','out'); drawnow; axis square; hold off;
        subplot(3,2,2*zz); set(gca,'Xlim',[-0.2 0.2],'Ylim',[-0.2 0.2],'YTick',-0.2:0.1:0.2,'XTick',[-0.2:0.1:0.2]); set(gca,'Tickdir','out'); drawnow; axis square; hold off;
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
    
    if zz == 1
        title('Simple');
    elseif zz==2
        title('DO');
    else
        title('Unclassified')
    end
    axis square; hold off;
    
    plot_counter = plot_counter + 1;
end

%% Figure S5- Population analysis of isoresponse curves 
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
indices = [109 24 74];
% Plotting the results for SVD based cone weight classification including the PC1 z-scores 
figure(plot_counter);
plot(RSSEisoresp_lin_median(hardtoclassifyidx),RSSEisoresp_quad_median(hardtoclassifyidx),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(RSSEisoresp_lin_median(LUMidx),RSSEisoresp_quad_median(LUMidx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(RSSEisoresp_lin_median(indices(1)),RSSEisoresp_quad_median(indices(1)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 1 0]);
plot(RSSEisoresp_lin_median(DOidx),RSSEisoresp_quad_median(DOidx),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); plot([0.0001 10],[0.0001 10],'k');
plot(RSSEisoresp_lin_median(indices(2)),RSSEisoresp_quad_median(indices(2)),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 1 0]);
plot(RSSEisoresp_lin_median(indices(3)),RSSEisoresp_quad_median(indices(3)),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0 1 0]);
axis square; set(gca,'Tickdir','out','Xlim',[0.0001 10],'Ylim',[0.0001 10],'YScale','log','XScale','log','XTick',[0.0001 0.001 0.01 0.1 1 10],'YTick',[0.0001 0.001 0.01 0.1 1 10]); xlabel('Linear error'); ylabel('Quadratic error'); title('Isoresponse'); hold off;

%% Figure S6 - Robustness of classification
if ~exist('plot_counter')
    plot_counter = 1;
end

load conewts_svd.mat
load vals.mat

LumIds_conewts = find(conewts_svd(1,:)>0 & conewts_svd(3,:) >0);
ColorOpponentIds_conewts = find(conewts_svd(1,:)<0 | conewts_svd(3,:) <0);

LUMidx = LumIds_conewts;
DOidx = [ColorOpponentIds_conewts];
hardtoclassifyidx = [LUMidx(vals(LUMidx)>=95) DOidx(vals(DOidx)>=95)];
LUMidx = LUMidx(vals(LUMidx)<95);
DOidx = DOidx(vals(DOidx)<95);

% Checking the correlation with non-linearity indices 
% Load the isoresponse data
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

%% Relationship between S-cone signal and the isoresponse NLI

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

% Checking whether there is any correlation between S-cone input and non-linearities 

% Load the isoresponse data
load RSSE_linearmodel_CV.mat % Robust regression
load RSSE_quadmodel_CV.mat

% For storing median of differences/ratios
RSSEisoresp_medianofratios = [];


for ii = 1:numel(RSSE_linearmodel)   
    % computation for calculating median of differences/ratios
    RSSEisoresp_medianofratios = [RSSEisoresp_medianofratios; median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii})];
end

indices = [109 24 74];
figure(plot_counter);  plot(RSSEisoresp_medianofratios(LUMidx),abs(conewts_svd(3,LUMidx)),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(RSSEisoresp_medianofratios(indices(1)),abs(conewts_svd(3,indices(1))),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 1 0]); hold on;
plot(RSSEisoresp_medianofratios(DOidx),abs(conewts_svd(3,DOidx)),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(RSSEisoresp_medianofratios(indices(2)),abs(conewts_svd(3,indices(2))),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 1 0]); hold on;
plot(RSSEisoresp_medianofratios(hardtoclassifyidx),abs(conewts_svd(3,hardtoclassifyidx)),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
plot(RSSEisoresp_medianofratios(indices(3)),abs(conewts_svd(3,indices(3))),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0 1 0]); hold on;
axis square; set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 1000],'XTick',[0.1 1 10 100 1000],'Ylim',[0 0.6],'YTick',0:0.2:0.6); xlabel('Spatial NLI isoresponse'); ylabel('S cone input');
plot_counter = plot_counter + 1;

[r1,p1] = corr(RSSEisoresp_medianofratios,abs(conewts_svd(3,:)'),'type','Spearman');

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


%% Some more analyses that Greg suggested
%  1) To make the isoprobability NLI and isoresponse NLI definition more consistent 
%  2) To split the other cell category into a) significant PC1 b) non-significant PC1

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
% Introducing new categories based on Greg's suggestion
hardtoclassifyidx_woPC1 = [Other_conewts(vals(Other_conewts)<95)];
hardtoclassifyidx_wPC1 = [Other_conewts(vals(Other_conewts)>=95)]; 
hardtoclassifyidx_wPC1_DOsimple = [LUMidx(vals(LUMidx)>=95) DOidx(vals(DOidx)>=95)]; 
LUMidx = LUMidx(vals(LUMidx)<95);
DOidx = DOidx(vals(DOidx)<95);

% Checking the correlation with non-linearity indices 
% Load the isoresponse data
load RSSE_linearmodel_CV.mat % Robust regression
load RSSE_quadmodel_CV.mat

% Load the integration within the subunit data
load AUROClinsubunits_CV.mat
load AUROCquadsubunits_CV.mat

% For storing median of differences/ratios
Acrosssubunits_medianofdifferences = [];
Acrosssubunits_lin_median = []; Acrosssubunits_quad_median = [];

% For storing median of differences/ratios
RSSEisoresp_medianofratios = [];
RSSEisoresp_lin_median = []; RSSEisoresp_quad_median = []; % Isoresponse data

for ii = 1:numel(RSSE_linearmodel)   
    % Isoresponse data - computation for calculating median of differences/ratios
    RSSEisoresp_medianofratios = [RSSEisoresp_medianofratios; median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii})];  
    RSSEisoresp_lin_median = [RSSEisoresp_lin_median; median(RSSE_linearmodel{ii})];
    RSSEisoresp_quad_median = [RSSEisoresp_quad_median; median(RSSE_quadmodel{ii})];
    
    % Storing the WN subunit spatial interaction data 
    Error_quad = 1-(AUROCquadsubunits{ii});
    Error_lin = 1-(AUROClinsubunits{ii});
    
    % A new definition of Acrosssubunits_medianofdifferences
    Acrosssubunits_medianofdifferences = [Acrosssubunits_medianofdifferences; median(Error_lin./Error_quad)];
    Acrosssubunits_lin_median = [Acrosssubunits_lin_median; median(Error_lin)];
    Acrosssubunits_quad_median = [Acrosssubunits_quad_median; median(Error_quad)];
    
end
RSSEisoresp_medianofratios(RSSEisoresp_medianofratios<0.1) = 0.1;
indices = [109 31 74];

% Isoresponse data: Plotting the results for SVD based cone weight classification including the PC1 z-scores 
figure(plot_counter);
subplot(221); plot(RSSEisoresp_lin_median(LUMidx),RSSEisoresp_quad_median(LUMidx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(RSSEisoresp_lin_median(DOidx),RSSEisoresp_quad_median(DOidx),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); plot([0.0001 10],[0.0001 10],'k');
axis square; set(gca,'Tickdir','out','Xlim',[0.0001 10],'Ylim',[0.0001 10],'YScale','log','XScale','log','XTick',[0.0001 0.001 0.01 0.1 1 10],'YTick',[0.0001 0.001 0.01 0.1 1 10]); 
xlabel('median Linear error'); ylabel('median Quadratic error'); legend ('LUM','DO'); title('Isoresponse'); hold off;

subplot(222); plot(RSSEisoresp_lin_median(hardtoclassifyidx_woPC1),RSSEisoresp_quad_median(hardtoclassifyidx_woPC1),'o','MarkerFaceColor',[1 0 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(RSSEisoresp_lin_median(hardtoclassifyidx_wPC1),RSSEisoresp_quad_median(hardtoclassifyidx_wPC1),'o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]); hold on;  
plot(RSSEisoresp_lin_median(hardtoclassifyidx_wPC1_DOsimple),RSSEisoresp_quad_median(hardtoclassifyidx_wPC1_DOsimple),'o','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]); hold on;
axis square; set(gca,'Tickdir','out','Xlim',[0.0001 10],'Ylim',[0.0001 10],'YScale','log','XScale','log','XTick',[0.0001 0.001 0.01 0.1 1 10],'YTick',[0.0001 0.001 0.01 0.1 1 10]); plot([0.0001 10],[0.0001 10],'k'); 
xlabel('median Linear error'); ylabel('median Quadratic error'); legend('HTC woPC1','HTC wPC1','DO simple wPC1'); title('Isoresponse'); hold off;

subplot(223); histogram(log10(RSSEisoresp_medianofratios(LUMidx)),-1:0.1:3,'DisplayStyle','stairs','EdgeColor',[0 0 0],'Linewidth',2); hold on;
histogram(log10(RSSEisoresp_medianofratios(DOidx)),-1:0.1:3,'DisplayStyle','stairs','EdgeColor',[1 0 0],'Linewidth',2); hold on;
plot(log10(median(RSSEisoresp_medianofratios(LUMidx))),14,'v','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
plot(log10(median(RSSEisoresp_medianofratios(DOidx))),15,'v','MarkerSize',8,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-1 3],'XTick',-1.0:0.5:3.0,'Ylim',[0 15],'YTick',[0 5 10 15]); 
ylabel('Count'); title('Isoresponse'); xlabel('Isoresponse NLI'); legend ('LUM','DO'); axis square; hold off;

subplot(224); histogram(log10(RSSEisoresp_medianofratios(hardtoclassifyidx_woPC1)),-1:0.1:3,'DisplayStyle','stairs','EdgeColor',[1 0 0.5],'Linewidth',2); hold on;
histogram(log10(RSSEisoresp_medianofratios(hardtoclassifyidx_wPC1)),-1:0.1:3,'DisplayStyle','stairs','EdgeColor',[1 0.5 0],'Linewidth',2); hold on;
histogram(log10(RSSEisoresp_medianofratios(hardtoclassifyidx_wPC1_DOsimple)),-1:0.1:3,'DisplayStyle','stairs','EdgeColor',[0 0.5 1.0],'Linewidth',2); hold on;
plot(log10(median(RSSEisoresp_medianofratios(hardtoclassifyidx_woPC1))),13,'v','MarkerSize',8,'MarkerFaceColor',[1.0 0 0.5],'MarkerEdgeColor',[1 1 1]); 
plot(log10(median(RSSEisoresp_medianofratios(hardtoclassifyidx_wPC1))),14,'v','MarkerSize',8,'MarkerFaceColor',[1.0 0.5 0],'MarkerEdgeColor',[1 1 1]); 
plot(log10(median(RSSEisoresp_medianofratios(hardtoclassifyidx_wPC1_DOsimple))),15,'v','MarkerSize',8,'MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]); 
set(gca,'Tickdir','out','Xlim',[-1 3],'XTick',-1.0:0.5:3.0,'Ylim',[0 15],'YTick',[0 5 10 15]); 
ylabel('Count'); title('Isoresponse'); xlabel('Isoresponse NLI'); legend('HTC woPC1','HTC wPC1','DO simple wPC1'); axis square; hold off;
plot_counter = plot_counter + 1;

% Comparing spatial NLI across cell types
group = [ones(size(LUMidx)) 2*ones(size(DOidx)) 3*ones(size(hardtoclassifyidx_woPC1)) 4*ones(size(hardtoclassifyidx_wPC1)) 5*ones(size(hardtoclassifyidx_wPC1_DOsimple))];
data = RSSEisoresp_medianofratios([LUMidx'; DOidx'; hardtoclassifyidx_woPC1'; hardtoclassifyidx_wPC1'; hardtoclassifyidx_wPC1_DOsimple']); 
% p1 = kruskalwallis(data,group);


% WN subunit spatial interaction data 
figure(plot_counter);
subplot(221); plot(Acrosssubunits_lin_median(LUMidx),Acrosssubunits_quad_median(LUMidx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(Acrosssubunits_lin_median(DOidx),Acrosssubunits_quad_median(DOidx),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); plot([0.08 1],[0.08 1],'k');
axis square; set(gca,'Tickdir','out','Xlim',[0.08 1],'Ylim',[0.08 1],'YScale','log','XScale','log','XTick',[0.08 0.1 1],'YTick',[0.08 0.1 1]); 
xlabel('median GLM error'); ylabel('median GQM error'); legend ('LUM','DO'); title('Isoprobability'); hold off;

subplot(222); plot(Acrosssubunits_lin_median(hardtoclassifyidx_woPC1),Acrosssubunits_quad_median(hardtoclassifyidx_woPC1),'o','MarkerFaceColor',[1 0 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(Acrosssubunits_lin_median(hardtoclassifyidx_wPC1),Acrosssubunits_quad_median(hardtoclassifyidx_wPC1),'o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]); hold on;  
plot(Acrosssubunits_lin_median(hardtoclassifyidx_wPC1_DOsimple),Acrosssubunits_quad_median(hardtoclassifyidx_wPC1_DOsimple),'o','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]); hold on;
axis square; set(gca,'Tickdir','out','Xlim',[0.08 1],'Ylim',[0.08 1],'YScale','log','XScale','log','XTick',[0.08 0.1 1],'YTick',[0.08 0.1 1]); plot([0.08 1],[0.08 1],'k'); 
xlabel('median GLM error'); ylabel('median GQM error'); legend('HTC woPC1','HTC wPC1','DO simple wPC1'); title('Isoprobability'); hold off;

subplot(223); histogram(log10(Acrosssubunits_medianofdifferences(LUMidx)),-0.02:0.005:0.1,'DisplayStyle','stairs','EdgeColor',[0 0 0],'Linewidth',2); hold on;
histogram(log10(Acrosssubunits_medianofdifferences(DOidx)),-0.02:0.005:0.1,'DisplayStyle','stairs','EdgeColor',[1 0 0],'Linewidth',2); hold on;
plot(log10(median(Acrosssubunits_medianofdifferences(LUMidx))),14,'v','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
plot(log10(median(Acrosssubunits_medianofdifferences(DOidx))),15,'v','MarkerSize',8,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-0.02 0.1],'XTick',-0.02:0.02:0.1,'Ylim',[0 15],'YTick',[0 5 10 15]); 
ylabel('Count'); title('Isoprobability'); xlabel('Isoprobability NLI'); legend ('LUM','DO'); axis square; hold off;

subplot(224); histogram(log10(Acrosssubunits_medianofdifferences(hardtoclassifyidx_woPC1)),-0.02:0.005:0.1,'DisplayStyle','stairs','EdgeColor',[1 0 0.5],'Linewidth',2); hold on;
histogram(log10(Acrosssubunits_medianofdifferences(hardtoclassifyidx_wPC1)),-0.02:0.005:0.1,'DisplayStyle','stairs','EdgeColor',[1 0.5 0],'Linewidth',2); hold on;
histogram(log10(Acrosssubunits_medianofdifferences(hardtoclassifyidx_wPC1_DOsimple)),-0.02:0.005:0.1,'DisplayStyle','stairs','EdgeColor',[0 0.5 1.0],'Linewidth',2); hold on;
plot(log10(median(Acrosssubunits_medianofdifferences(hardtoclassifyidx_woPC1))),13,'v','MarkerSize',8,'MarkerFaceColor',[1.0 0 0.5],'MarkerEdgeColor',[1 1 1]); 
plot(log10(median(Acrosssubunits_medianofdifferences(hardtoclassifyidx_wPC1))),14,'v','MarkerSize',8,'MarkerFaceColor',[1.0 0.5 0],'MarkerEdgeColor',[1 1 1]); 
plot(log10(median(Acrosssubunits_medianofdifferences(hardtoclassifyidx_wPC1_DOsimple))),15,'v','MarkerSize',8,'MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]); 
set(gca,'Tickdir','out','Xlim',[-0.02 0.1],'XTick',-0.02:0.02:0.1,'Ylim',[0 15],'YTick',[0 5 10 15]); 
ylabel('Count'); title('Isoprobability'); xlabel('Isoprobability NLI'); legend('HTC woPC1','HTC wPC1','DO simple wPC1'); axis square; hold off;

plot_counter = plot_counter + 1;


%% A monte-carlo simulation to quantify how DO cells were selected
% Don't remember why I did this - AD, 12/20

iter = 10000;

plot_counter = 1;
thresh = 0.8;
for ii=1:iter
    pts = randn(3,1);
    pts = sign(pts(2))*pts./sum(abs(pts));
    
    flag = pts(2) - pts(1) >thresh & sum(sign(pts(1:2)),1)==0 & sqrt((pts(2,:)-0.5).^2 + (pts(1,:)+0.5).^2)<0.3;
    
    if flag 
        plot(pts(1), pts(2),'k.'); hold on;
    end
end
axis square; set(gca,'Tickdir','out','Xlim',[-1.0 0], 'Ylim',[0 1.0]); hold off;

