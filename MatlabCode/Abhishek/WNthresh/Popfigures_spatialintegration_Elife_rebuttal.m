% Spatial intergration rebuttal analysis for eLife spatial integration
% paper
% Author - Abhishek De, 05/30/2020

close all; clearvars;
plot_counter = 1;

%% Figure 1: Impact of linear and non-linear spatial filtering in image processing 
% This code has been derived from Abhishek/Physiology_modeling/Edge_processing_linear_nonlinear
%*********************************************************************************
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
thresh = 0;
I_LminusM(I_LminusM<thresh) = thresh; I_MminusL(I_MminusL<thresh) = thresh;
% I_color_filteredNL = I_LminusM.^2 + I_MminusL.^2; % Squaring non-linearity
I_color_filteredNL = I_LminusM + I_MminusL; % The vanilla ReLu non-linearity
%I_color_filteredNL(I_color_filteredNL<0) = 0;
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

% Modying the I_color_filtered to enhance the constrast
K = I_color_filtered;
K(K>0.55) = 0.75;
K(K>0.53 & K<0.55) = 0.6;

% Plotting the images
plot_counter = 1;
figure(plot_counter); set(gcf,'Name','Analyses of images: RGB');
subplot(321); imagesc(A); axis square; set(gca,'XTick',[],'YTick',[]); title('Filter');
subplot(322); imshow(I_color); axis square; set(gca,'XTick',[],'YTick',[]); title('original pic');
subplot(323); imshow(im2double(I_LminusM_norm)); axis square; set(gca,'XTick',[],'YTick',[]); colormap('gray'); title('Rectified L-M map');
subplot(324); imshow(im2double(I_MminusL_norm)); axis square; set(gca,'XTick',[],'YTick',[]); colormap('gray'); title('Rectified M-L map');
subplot(325); imshow(im2double(K)); axis square; set(gca,'XTick',[],'YTick',[]); colormap('gray'); hold off; title('Linear filter');
subplot(326); imshow(im2double(I_color_filteredNL)); axis square; set(gca,'XTick',[],'YTick',[]); colormap('gray'); hold off; title('Relu non-linearity');
plot_counter = plot_counter + 1;

%% Trying out the non-rectified maps
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
%*********************************************************************************

if ~exist('plot_counter')
    plot_counter = 1;
end

% Loading all the files
try 
    conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
    filename = fetch(conn,'SELECT filename FROM WNthresh');
    NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
    spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
    close(conn);
    filename = filename(strcmp(string(NTmode),"subunit"));
    NTmode = NTmode(strcmp(string(NTmode),"subunit"));
    spikeIdx = spikeidx_NT(strcmp(string(NTmode),"subunit"));

catch
    csv_filename = '/Users/abhishekde/Desktop/MatlabCode/Abhishek/CSV_PHPmyadmin_files/WNthresh.csv';
    [filename, NTmode, spikeIdx] = get_WNthreshdata_from_csvfile(csv_filename, 'subunit');
    spikeIdx = str2num(cell2mat(spikeIdx));
end


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
try 
    conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
    filename = fetch(conn,'SELECT filename FROM WNthresh');
    NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
    spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
    close(conn);
    filename = filename(strcmp(string(NTmode),"subunit"));
    NTmode = NTmode(strcmp(string(NTmode),"subunit"));
    spikeIdx = spikeidx_NT(strcmp(string(NTmode),"subunit"));

catch
    csv_filename = '/Users/abhishekde/Desktop/MatlabCode/Abhishek/CSV_PHPmyadmin_files/WNthresh.csv';
    [filename, NTmode, spikeIdx] = get_WNthreshdata_from_csvfile(csv_filename, 'subunit');
    spikeIdx = str2num(cell2mat(spikeIdx));
end
    
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
filenumber = 74; % blue-yellow
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

imR = sigmoid(squeeze(weighted_im(:,:,1)),7,0.5);
imG = sigmoid(squeeze(weighted_im(:,:,2)),5,0.5);
imB = sigmoid(squeeze(weighted_im(:,:,3)),3,0.5);
im_new = cat(3,imR,imG,imB);

% Performing the same thing for the spatial weighting function
weighted_sim = reshape(spatial_weighting,[10 10]).*(energymap/max(energymap(:)));
normfactor = 0.5/(max(abs(weighted_sim(:)))+0.01);
weighted_sim = normfactor*weighted_sim + 0.5;
sim = sigmoid(weighted_sim,2,0.5);

figure(plot_counter);
subplot(121); imagesc(im_new); axis square; set(gca,'XTick',[],'YTick',[]); title('Same size STA')
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
load S1RGB_svd.mat
load S2RGB_svd.mat
load angulardifferences_RGB.mat
S1RGB = S1RGB_svd;
S2RGB = S2RGB_svd;
anglebwvectors = angulardifference_RGB;
% SpatiallyOpponent = sum(sign(S1RGB).*sign(S2RGB),1)<3;
SpatiallyOpponent = anglebwvectors'>90;

% Calculating some metrics for circular stats and discarding SO cells
unitvectors = [];
for ii = 1:size(S1RGB,2)
    vec1 = S1RGB(:,ii);
    vec2 = S2RGB(:,ii);
    vec = exp(angulardifference_RGB(ii)*(pi/180)*i);
    
    % Converting the angle difference into unit vectors
    unitvectors = [unitvectors; real(vec) imag(vec)];
end

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

LUMidx = LumIds_conewts;
DOidx = [ColorOpponentIds_conewts Sconedominated_conewts];
hardtoclassifyidx = [Other_conewts];
hardtoclassifyidx = [hardtoclassifyidx LUMidx(vals(LUMidx)>=95) DOidx(vals(DOidx)>=95)];
LUMidx = LUMidx(vals(LUMidx)<95);
DOidx = DOidx(vals(DOidx)<95);

% Considering only the spatially opponent subunits
LUMidx = LUMidx(SpatiallyOpponent(LUMidx));
DOidx = DOidx(SpatiallyOpponent(DOidx));
hardtoclassifyidx = hardtoclassifyidx(SpatiallyOpponent(hardtoclassifyidx));

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

% DO vs. Simple cells
[p4,~] = ranksum(Whitenoise_NLI([DOidx]),Whitenoise_NLI(LUMidx));

% DO vs. hardtoclassify cells, LUM vs. hardtoclassify cells
[p5,~] = ranksum(Whitenoise_NLI([DOidx]),Whitenoise_NLI(hardtoclassifyidx));
[p6,~] = ranksum(Whitenoise_NLI([LUMidx]),Whitenoise_NLI(hardtoclassifyidx));

%*************************************CIRCULAR STATS *******************
% I am still unsure about the math here: Need to double check 
% Summary circular stats for all cells
All = [LUMidx DOidx hardtoclassifyidx];
% Calculating the circular mean and circular standard deviation
circmean = mean(atan2d(unitvectors(All',2),unitvectors(All',1)));
circstd = std(atan2d(unitvectors(All',2),unitvectors(All',1)));
Rmean = mean(cos(angle(unitvectors(All,1)+ i*unitvectors(All,2)))); % Pearson correlation coefficient using angles
Rstd = std(cos(angle(unitvectors(All,1)+ i*unitvectors(All,2)))); 

% Summary circular stats for Simple cells
% Calculating the circular mean and circular standard deviation
circmean_LUM = mean(atan2d(unitvectors(LUMidx',2),unitvectors(LUMidx',1)));
circstd_LUM = std(atan2d(unitvectors(LUMidx',2),unitvectors(LUMidx',1)));
Rmean_LUM = mean(cos(angle(unitvectors(LUMidx,1)+ i*unitvectors(LUMidx,2)))); 
Rstd_LUM = std(cos(angle(unitvectors(LUMidx,1)+ i*unitvectors(LUMidx,2))));

% Summary circular stats for DO cells
% Calculating the circular mean and circular standard deviation
circmean_DO = mean(atan2d(unitvectors(DOidx',2),unitvectors(DOidx',1)));
circstd_DO = std(atan2d(unitvectors(DOidx',2),unitvectors(DOidx',1)));
Rmean_DO = mean(cos(angle(unitvectors(DOidx,1)+ i*unitvectors(DOidx,2)))); 
Rstd_DO = std(cos(angle(unitvectors(DOidx,1)+ i*unitvectors(DOidx,2))));

% Summary circular stats for HTC cells
% Calculating the circular mean and circular standard deviation
circmean_HTC = mean(atan2d(unitvectors(hardtoclassifyidx',2),unitvectors(hardtoclassifyidx',1)));
circstd_HTC = std(atan2d(unitvectors(hardtoclassifyidx',2),unitvectors(hardtoclassifyidx',1)));
Rmean_HTC = mean(cos(angle(unitvectors(hardtoclassifyidx,1)+ i*unitvectors(hardtoclassifyidx,2)))); 
Rstd_HTC = std(cos(angle(unitvectors(hardtoclassifyidx,1)+ i*unitvectors(hardtoclassifyidx,2))));

%% Figure 3-part1: Iso-response example data from example LUM, DO and HTC cells 
% Contains OLD indexes but still valid with the NEW criteria, so not changing this part of the code
%*********************************************************************************

if ~exist('plot_counter')
    plot_counter = 1;
end

load conewts_svd.mat
load vals.mat
thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

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
try 
    % Using the JDBC connection
    conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
    filename = fetch(conn,'SELECT filename FROM WNthresh');
    NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
    spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
    close(conn);
    filename = filename(strcmp(string(NTmode),"subunit"));
    NTmode = NTmode(strcmp(string(NTmode),"subunit"));
    spikeidx_NT = spikeidx_NT(strcmp(string(NTmode),"subunit"));

catch
    csv_filename = '/Users/abhishekde/Desktop/MatlabCode/Abhishek/CSV_PHPmyadmin_files/WNthresh.csv';
    [filename, NTmode, spikeIdx] = get_WNthreshdata_from_csvfile(csv_filename, 'subunit');
    spikeidx_NT = str2num(cell2mat(spikeIdx));
end

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

%% Plotting the cone wts
if ~exist('plot_counter')
    plot_counter = 1;
end

if ~exist('plot_counter')
    plot_counter = 1;
end

indices = [109 24 74];

figure(plot_counter);
subplot(311); bar(conewts_svd(:,indices(2)),'FaceColor','k'); 
set(gca,'Ylim',[-1 1],'Xlim',[0 4],'XTick',[1 2 3],'XTicklabel',{'L','M','S'},'Tickdir','out'); axis square;
subplot(312); bar(conewts_svd(:,indices(1)),'FaceColor','k'); 
set(gca,'Ylim',[-1 1],'Xlim',[0 4],'XTick',[1 2 3],'XTicklabel',{'L','M','S'},'Tickdir','out'); axis square;
subplot(313); bar(conewts_svd(:,indices(3)),'FaceColor','k'); 
set(gca,'Ylim',[-1 1],'Xlim',[0 4],'XTick',[1 2 3],'XTicklabel',{'L','M','S'},'Tickdir','out'); axis square;
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;

%% Figure 3-part 2: Population analysis of isoresponse curves 

if ~exist('plot_counter')
    plot_counter = 1;
end

load conewts_svd.mat
load vals.mat
load vals.mat
load S1RGB_svd.mat
load S2RGB_svd.mat
load angulardifferences_RGB.mat
anglebwvectors = angulardifference_RGB;
S1RGB = S1RGB_svd;
S2RGb = S2RGB_svd;
% SpatiallyOpponent = sum(sign(S1RGB).*sign(S2RGB),1)<3;
SpatiallyOpponent = anglebwvectors'>90;

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

LUMidx = LumIds_conewts;
DOidx = [ColorOpponentIds_conewts Sconedominated_conewts];
hardtoclassifyidx = [Other_conewts];
hardtoclassifyidx = [hardtoclassifyidx LUMidx(vals(LUMidx)>=95) DOidx(vals(DOidx)>=95)];
LUMidx = LUMidx(vals(LUMidx)<95);
DOidx = DOidx(vals(DOidx)<95);

% Considering only the spatially opponent subunits
LUMidx = LUMidx(SpatiallyOpponent(LUMidx));
DOidx = DOidx(SpatiallyOpponent(DOidx));
hardtoclassifyidx = hardtoclassifyidx(SpatiallyOpponent(hardtoclassifyidx));

% Checking the correlation with non-linearity indices 
% Load the isoresponse data
load RSSE_linearmodel_CV.mat % Robust regression
load RSSE_quadmodel_CV.mat

% For storing median of differences/ratios
RSSEisoresp_medianofratios = [];
RSSEisoresp_lin_median = []; RSSEisoresp_quad_median = []; % Isoresponse data
Isoresponse_NLI = [];

for ii = 1:numel(RSSE_linearmodel)   
    % computation for calculating median of differences/ratios
    RSSEisoresp_medianofratios = [RSSEisoresp_medianofratios; median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii})];  
    Isoresponse_NLI = [Isoresponse_NLI; log10(RSSEisoresp_medianofratios(end))];
    
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
data = Isoresponse_NLI([LUMidx'; DOidx'; hardtoclassifyidx']); 
p1 = kruskalwallis(data,group,'off');

% Comparing spatial NLI between simple and DO cells
[p2,h] = ranksum(Isoresponse_NLI(LUMidx),Isoresponse_NLI(DOidx));

% Comparing spatial NLI between simple + DO cells and other cells
[p5,h] = ranksum(Isoresponse_NLI([LUMidx DOidx]),Isoresponse_NLI(hardtoclassifyidx));

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

%% Figure 3 - final part: plotting the relationship between the whitenoise and isoresponse NLI 

if ~exist('plot_counter')
    plot_counter = 1;
end

load conewts_svd.mat
load vals.mat

load S1RGB_svd.mat
load S2RGB_svd.mat
load angulardifferences_RGB.mat
anglebwvectors = angulardifference_RGB;
S1RGB = S1RGB_svd;
S2RGb = S2RGB_svd;
% SpatiallyOpponent = sum(sign(S1RGB).*sign(S2RGB),1)<3;
SpatiallyOpponent = anglebwvectors'>90;

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];


LUMidx = LumIds_conewts;
DOidx = [ColorOpponentIds_conewts Sconedominated_conewts];
hardtoclassifyidx = [Other_conewts];
hardtoclassifyidx = [hardtoclassifyidx LUMidx(vals(LUMidx)>=95) DOidx(vals(DOidx)>=95)];
LUMidx = LUMidx(vals(LUMidx)<95);
DOidx = DOidx(vals(DOidx)<95);

% Considering only the spatially opponent subunits
LUMidx = LUMidx(SpatiallyOpponent(LUMidx));
DOidx = DOidx(SpatiallyOpponent(DOidx));
hardtoclassifyidx = hardtoclassifyidx(SpatiallyOpponent(hardtoclassifyidx));

% Checking the correlation with non-linearity indices 
% Load the isoresponse data
load RSSE_linearmodel_CV.mat % Robust regression
load RSSE_quadmodel_CV.mat
% Load the integration within the subunit data
load AUROClinsubunits_CV.mat
load AUROCquadsubunits_CV.mat


Isoresponse_NLI = []; % For storing the Isoresponse NLI
Whitenoise_NLI = []; % For storing the white noise NLI
JK_error_whitenoise = []; % Storing Jackknife based error estimates
JK_error_isoresponse = []; % Storing Jackknife based error estimates

indices = [109 24 74];
for ii = 1:numel(AUROClinsubunits) 
  
    % White noise NLI 
    Error_quad = 1-(AUROCquadsubunits{ii});
    Error_lin = 1-(AUROClinsubunits{ii});
    Whitenoise_NLI = [Whitenoise_NLI; log10(median(Error_lin./Error_quad))];
    
    % Isoresponse NLI
    Isoresponse_NLI = [Isoresponse_NLI; log10(median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii}))];
    
    % Computing error via Jackknife resampling
    % Jackknife White noise error
    N1 = numel(Error_quad);
    X1 = log10(jackknife(@median, Error_lin./Error_quad));
    error_WN = sqrt(sum((X1-mean(X1)).^2)*(N1-1)/N1);
    JK_error_whitenoise = [JK_error_whitenoise; error_WN];
    
    % Jackknife Isoresponse error
    N2 = numel(RSSE_linearmodel{ii});
    X2 = log10(jackknife(@median, RSSE_linearmodel{ii}./RSSE_quadmodel{ii}));
    error_Iso = sqrt(sum((X2-mean(X2)).^2)*(N2-1)/N2);
    JK_error_isoresponse = [JK_error_isoresponse; error_Iso];
end

% figure(plot_counter); hold on;
% plot(Isoresponse_NLI(LUMidx), Whitenoise_NLI(LUMidx), 'o', 'MarkerFaceColor', [0 0 0],'MarkerEdgeColor',[1 1 1]);
% plot(Isoresponse_NLI(indices(1)), Whitenoise_NLI(indices(1)), 'o', 'MarkerFaceColor', [0 0 0],'MarkerEdgeColor',[0 1 0]);
% plot(Isoresponse_NLI(DOidx), Whitenoise_NLI(DOidx), 'o', 'MarkerFaceColor', [1 0 0],'MarkerEdgeColor',[1 1 1]);
% plot(Isoresponse_NLI(indices(2)), Whitenoise_NLI(indices(2)), 'o', 'MarkerFaceColor', [1 0 0],'MarkerEdgeColor',[0 1 0]);
% plot(Isoresponse_NLI(hardtoclassifyidx), Whitenoise_NLI(hardtoclassifyidx), 'o', 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
% plot(Isoresponse_NLI(indices(3)), Whitenoise_NLI(indices(3)), 'o', 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor',[0 1 0]);
% set(gca,'Tickdir','out','Xlim',[-1 2],'XTick',[-1 0 1 2],'Ylim',[-0.02 0.08],'YTick',[-0.02:0.02:0.08]); 
% xlabel('Isoresponse NLI'); ylabel('WhiteNoise NLI'); axis square; hold off;
% set(gcf,'renderer','painters');
% plot_counter = plot_counter + 1;

% Some basic stats
[r1,p1] = corr(Isoresponse_NLI(LUMidx),Whitenoise_NLI(LUMidx),'type','Spearman');
[r2,p2] = corr(Isoresponse_NLI(DOidx),Whitenoise_NLI(DOidx),'type','Spearman');
[r3,p3] = corr(Isoresponse_NLI(hardtoclassifyidx),Whitenoise_NLI(hardtoclassifyidx),'type','Spearman');
[rc,pc] = corr(Isoresponse_NLI,Whitenoise_NLI,'type','Spearman');


% Plotting the same plot as above but with errors
figure(plot_counter); hold on;
errorbar(Isoresponse_NLI(hardtoclassifyidx), Whitenoise_NLI(hardtoclassifyidx), JK_error_whitenoise(hardtoclassifyidx), JK_error_whitenoise(hardtoclassifyidx), JK_error_isoresponse(hardtoclassifyidx), JK_error_isoresponse(hardtoclassifyidx), 'o', 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1], 'color', [0.5 0.5 0.5], 'MarkerSize', 8);
errorbar(Isoresponse_NLI(LUMidx), Whitenoise_NLI(LUMidx),JK_error_whitenoise(LUMidx), JK_error_whitenoise(LUMidx), JK_error_isoresponse(LUMidx), JK_error_isoresponse(LUMidx), 'o', 'MarkerFaceColor', [0 0 0],'MarkerEdgeColor',[1 1 1], 'color', [0 0 0], 'MarkerSize', 8);
errorbar(Isoresponse_NLI(DOidx), Whitenoise_NLI(DOidx), JK_error_whitenoise(DOidx), JK_error_whitenoise(DOidx), JK_error_isoresponse(DOidx), JK_error_isoresponse(DOidx), 'o', 'MarkerFaceColor', [1 0 0],'MarkerEdgeColor',[1 1 1], 'color', [1 0 0], 'MarkerSize', 8);

errorbar(Isoresponse_NLI(indices(3)), Whitenoise_NLI(indices(3)), JK_error_whitenoise(indices(3)), JK_error_whitenoise(indices(3)), JK_error_isoresponse(indices(3)), JK_error_isoresponse(indices(3)), 'o', 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor',[0 1 0], 'color', [0.5 0.5 0.5], 'MarkerSize', 8);
errorbar(Isoresponse_NLI(indices(1)), Whitenoise_NLI(indices(1)), JK_error_whitenoise(indices(1)), JK_error_whitenoise(indices(1)), JK_error_isoresponse(indices(1)), JK_error_isoresponse(indices(1)), 'o', 'MarkerFaceColor', [0 0 0],'MarkerEdgeColor',[0 1 0], 'color', [0 0 0], 'MarkerSize', 8);
errorbar(Isoresponse_NLI(indices(2)), Whitenoise_NLI(indices(2)), JK_error_whitenoise(indices(2)), JK_error_whitenoise(indices(2)), JK_error_isoresponse(indices(2)), JK_error_isoresponse(indices(2)), 'o', 'MarkerFaceColor', [1 0 0],'MarkerEdgeColor',[0 1 0], 'color', [1 0 0], 'MarkerSize', 8);
set(gca,'Tickdir','out','Xlim',[-4 3],'XTick',-4:3,'Ylim',[-0.04 0.14],'YTick',-0.04:0.02:0.14); 
xlabel('Isoresponse NLI'); ylabel('WhiteNoise NLI'); axis square; hold off;
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;


%% Figure 4: Conceptual model of signal integration within subunits (cone-signal integration)
% The visual effects were enhanced using Illustrator
%*********************************************************************************

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

%% Figure 5: Quantification of signal integration within subunit: Cone signal NLI 
%*********************************************************************************

if ~exist('plot_counter')
    plot_counter = 1;
end

% Laoding data for cone weight calculation  
load conewts_svd.mat
load vals.mat
load S1RGB_svd.mat
load S2RGB_svd.mat
load angulardifferences_RGB.mat
anglebwvectors = angulardifference_RGB;
S1RGB = S1RGB_svd;
S2RGb = S2RGB_svd;
SpatiallyOpponent = anglebwvectors'>90;

% Classifying cells based on cone-weights and PC1 signficance
thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

LUMidx = LumIds_conewts;
DOidx = [ColorOpponentIds_conewts Sconedominated_conewts];
hardtoclassifyidx = [Other_conewts];
hardtoclassifyidx = [hardtoclassifyidx LUMidx(vals(LUMidx)>=95) DOidx(vals(DOidx)>=95)];
LUMidx = LUMidx(vals(LUMidx)<95);
DOidx = DOidx(vals(DOidx)<95);

% Considering only the spatially opponent subunits
LUMidx = LUMidx(SpatiallyOpponent(LUMidx));
DOidx = DOidx(SpatiallyOpponent(DOidx));
hardtoclassifyidx = hardtoclassifyidx(SpatiallyOpponent(hardtoclassifyidx));

% Load the integration within the subunit data
load AUROClinS1_CV.mat
load AUROCquadS1_CV.mat
load AUROClinS2_CV.mat
load AUROCquadS2_CV.mat

% Storing the within subunit NLI
Within_subunit_NLI = [];

for ii = 1:numel(AUROClin1) 
    
    % GLM error for both the subunits
    Error_lin_subunit1 = 1-AUROClin1{ii};
    Error_lin_subunit2 = 1-AUROClin2{ii};
    
    % GQM error for both the subunits
    Error_quad_subunit1 = 1-AUROCquad1{ii};
    Error_quad_subunit2 = 1-AUROCquad2{ii};
    
    % Calculating the median ratio of errors both the subunits in log scale
    median_NLI_subunit1 = log10(median(Error_lin_subunit1./ Error_quad_subunit1));
    median_NLI_subunit2 = log10(median(Error_lin_subunit2./ Error_quad_subunit2));
    
    % Storing the NLIs
    Within_subunit_NLI = [Within_subunit_NLI; median([median_NLI_subunit1 median_NLI_subunit2])];
end


% Example cells from each cell type: Simple, DO and NSNDO (hardtoclassify) cells  
indices = [109 24 74];

% Simple cells
figure(plot_counter);
subplot(311); histogram(Within_subunit_NLI(LUMidx),linspace(-0.02,0.08,21),'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); hold on;
plot(median(Within_subunit_NLI(LUMidx)),20,'v','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
plot(Within_subunit_NLI(indices(1)),19,'s','MarkerSize',8,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-0.02 0.08],'Ylim',[0 20],'YTick',[0 10 20]); xlabel('Within subunit NLI'); ylabel('# cells'); axis square; hold off;

% DO cells
subplot(312); histogram(Within_subunit_NLI(DOidx),linspace(-0.02,0.08,21),'FaceColor',[1 0 0],'EdgeColor',[1 1 1]); hold on;
plot(median(Within_subunit_NLI(DOidx)),15,'v','MarkerSize',8,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(Within_subunit_NLI(indices(2)),12,'s','MarkerSize',8,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-0.02 0.08],'Ylim',[0 15],'YTick',[0 5 10 15]); xlabel('Within subunit NLI'); ylabel('# cells'); axis square; hold off;

% hardtoclassify cells-NSNDO cells
subplot(313); histogram(Within_subunit_NLI(hardtoclassifyidx),linspace(-0.02,0.08,21),'FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1]); hold on;
plot(median(Within_subunit_NLI(hardtoclassifyidx)),20,'v','MarkerSize',8,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
plot(Within_subunit_NLI(indices(3)),18,'s','MarkerSize',8,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-0.02 0.08],'Ylim',[0 20],'YTick',[0 10 20]); xlabel('Within subunit NLI'); ylabel('# cells'); axis square; hold off;
plot_counter = plot_counter + 1;

% Statistical tests
%******************
% Kruskal Wallis test: Comparing 'within subunit NLIs' across cell types
group = [ones(size(LUMidx)) 2*ones(size(DOidx)) 3*ones(size(hardtoclassifyidx))];
data = Within_subunit_NLI([LUMidx'; DOidx'; hardtoclassifyidx']); 
p1 = kruskalwallis(data,group,'off');

% Mann-Whitney U test: DO vs. Simple cells 
[p2,~] = ranksum(Within_subunit_NLI([LUMidx]),Within_subunit_NLI(DOidx));

% Mann-Whitney U test: DO vs. hardtoclassify cells 
[p3,~] = ranksum(Within_subunit_NLI([DOidx]),Within_subunit_NLI(hardtoclassifyidx));

% Mann-Whitney U test: Simple vs. hardtoclassify cells 
[p4,~] = ranksum(Within_subunit_NLI([LUMidx]),Within_subunit_NLI(hardtoclassifyidx));

% Mann-Whitney U test: DO + simple cells vs hardtoclassify cells
[p5,~] = ranksum(Within_subunit_NLI([LUMidx DOidx]),Within_subunit_NLI(hardtoclassifyidx));

% Comparing the cone signal NLIs with simple and DO cells
load ConesignalNLI_LGN.mat

[p6, ~]=  ranksum(Within_subunit_NLI([LUMidx]), ConesignalNLI_LGN);
[p7, ~]=  ranksum(Within_subunit_NLI([DOidx]), ConesignalNLI_LGN);
[p8, ~]=  ranksum(Within_subunit_NLI([DOidx LUMidx]), ConesignalNLI_LGN);
[p9, ~] = ranksum(Within_subunit_NLI([hardtoclassifyidx]), ConesignalNLI_LGN);
[p10, ~]=  signrank(ConesignalNLI_LGN);

%% Figure 6 (OLD): Relationship between cone signal NLI to white noise NLI and Isoresponse NLI- 
%NO LONGER USED AS OF 9/21
%*********************************************************************************

if ~exist('plot_counter')
    plot_counter = 1;
end

% Laoding data for cone weight calculation  
load conewts_svd.mat
load vals.mat
load S1RGB_svd.mat
load S2RGB_svd.mat
load angulardifferences_RGB.mat
anglebwvectors = angulardifference_RGB;
S1RGB = S1RGB_svd;
S2RGb = S2RGB_svd;
SpatiallyOpponent = anglebwvectors'>90;

% Classifying cells based on cone-weights and PC1 signficance
thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

LUMidx = LumIds_conewts;
DOidx = [ColorOpponentIds_conewts Sconedominated_conewts];
hardtoclassifyidx = [Other_conewts];
hardtoclassifyidx = [hardtoclassifyidx LUMidx(vals(LUMidx)>=95) DOidx(vals(DOidx)>=95)];
LUMidx = LUMidx(vals(LUMidx)<95);
DOidx = DOidx(vals(DOidx)<95);

% Considering only the spatially opponent subunits
LUMidx = LUMidx(SpatiallyOpponent(LUMidx));
DOidx = DOidx(SpatiallyOpponent(DOidx));
hardtoclassifyidx = hardtoclassifyidx(SpatiallyOpponent(hardtoclassifyidx));

% Load the integration within the subunit data
load AUROClinS1_CV.mat
load AUROCquadS1_CV.mat
load AUROClinS2_CV.mat
load AUROCquadS2_CV.mat

% Load the integration within the subunit whitenoise analysis data
load AUROClinsubunits_CV.mat
load AUROCquadsubunits_CV.mat

% Load the isoresponse data
load RSSE_linearmodel_CV.mat % Robust regression
load RSSE_quadmodel_CV.mat

% Storing the within subunit NLI
Within_subunit_NLI = [];

% For storing the Isoresponse NLI
Isoresponse_NLI = [];

% For storing the white noise NLI
Whitenoise_NLI = [];

for ii = 1:numel(AUROClin1) 
    
    % GLM error for both the subunits
    Error_lin_subunit1 = 1-AUROClin1{ii};
    Error_lin_subunit2 = 1-AUROClin2{ii};
    
    % GQM error for both the subunits
    Error_quad_subunit1 = 1-AUROCquad1{ii};
    Error_quad_subunit2 = 1-AUROCquad2{ii};
    
    % Calculating the median ratio of errors both the subunits in log scale
    median_NLI_subunit1 = log10(median(Error_lin_subunit1./ Error_quad_subunit1));
    median_NLI_subunit2 = log10(median(Error_lin_subunit2./ Error_quad_subunit2));
    
    % Storing the NLIs
    Within_subunit_NLI = [Within_subunit_NLI; median([median_NLI_subunit1 median_NLI_subunit2])];
    
    
    % White noise NLI 
    Error_quad = 1-(AUROCquadsubunits{ii});
    Error_lin = 1-(AUROClinsubunits{ii});
    Whitenoise_NLI = [Whitenoise_NLI; log10(median(Error_lin./Error_quad))];
    
    % Isoresponse NLI
    Isoresponse_NLI = [Isoresponse_NLI; log10(median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii}))];
end


indices = [109 24 74];
figure(plot_counter); 
% White noise NLI vs. Cone-signal NLI
subplot(121); plot(Whitenoise_NLI(LUMidx),Within_subunit_NLI(LUMidx),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(Whitenoise_NLI(indices(1)),Within_subunit_NLI(indices(1)),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 1 0]); hold on;
plot(Whitenoise_NLI(DOidx),Within_subunit_NLI(DOidx),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(Whitenoise_NLI(indices(2)),Within_subunit_NLI(indices(2)),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 1 0]); hold on;
plot(Whitenoise_NLI(hardtoclassifyidx),Within_subunit_NLI(hardtoclassifyidx),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
plot(Whitenoise_NLI(indices(3)),Within_subunit_NLI(indices(3)),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0 1 0]); hold on;
axis square; set(gca,'Tickdir','out','Xlim',[-0.02 0.08],'XTick',-0.02:0.02:0.08,'Ylim',[-0.02 0.08],'YTick',-0.02:0.02:0.08); 
xlabel('White noise NLI'); ylabel('Cone signal NLI');

% Isoresponse NLI vs. Cone-signal NLI
subplot(122); plot(Isoresponse_NLI(LUMidx),Within_subunit_NLI(LUMidx),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(Isoresponse_NLI(indices(1)),Within_subunit_NLI(indices(1)),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 1 0]); hold on;
plot(Isoresponse_NLI(DOidx),Within_subunit_NLI(DOidx),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(Isoresponse_NLI(indices(2)),Within_subunit_NLI(indices(2)),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 1 0]); hold on;
plot(Isoresponse_NLI(hardtoclassifyidx),Within_subunit_NLI(hardtoclassifyidx),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
plot(Isoresponse_NLI(indices(3)),Within_subunit_NLI(indices(3)),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0 1 0]); hold on;
axis square; set(gca,'Tickdir','out','Xlim',[-1 2],'XTick',-1:1:2,'Ylim',[-0.02 0.08],'YTick',-0.02:0.02:0.08); 
xlabel('Isoresponse NLI'); ylabel('Cone signal NLI');
set(gcf, 'renderer', 'painters');
plot_counter = plot_counter + 1;

[r1,p1] = corr(Whitenoise_NLI([LUMidx, DOidx, hardtoclassifyidx]),Within_subunit_NLI([LUMidx, DOidx, hardtoclassifyidx]),'type','Spearman');
[r2,p2] = corr(Isoresponse_NLI([LUMidx, DOidx, hardtoclassifyidx]),Within_subunit_NLI([LUMidx, DOidx, hardtoclassifyidx]),'type','Spearman');

%% Figure 6: Relationship between S-cone strength and NLIs 

if ~exist('plot_counter')
    plot_counter = 1;
end

% Laoding data for cone weight calculation  
load conewts_svd.mat
load vals.mat
load S1RGB_svd.mat
load S2RGB_svd.mat
load angulardifferences_RGB.mat
anglebwvectors = angulardifference_RGB;
S1RGB = S1RGB_svd;
S2RGb = S2RGB_svd;
SpatiallyOpponent = anglebwvectors'>90;

% Classifying cells based on cone-weights and PC1 signficance
thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

LUMidx = LumIds_conewts;
DOidx = [ColorOpponentIds_conewts Sconedominated_conewts];
hardtoclassifyidx = [Other_conewts];
hardtoclassifyidx = [hardtoclassifyidx LUMidx(vals(LUMidx)>=95) DOidx(vals(DOidx)>=95)];
LUMidx = LUMidx(vals(LUMidx)<95);
DOidx = DOidx(vals(DOidx)<95);

% Considering only the spatially opponent subunits
LUMidx = LUMidx(SpatiallyOpponent(LUMidx));
DOidx = DOidx(SpatiallyOpponent(DOidx));
hardtoclassifyidx = hardtoclassifyidx(SpatiallyOpponent(hardtoclassifyidx)); 

% Load the integration within the subunit data
load AUROClinS1_CV.mat
load AUROCquadS1_CV.mat
load AUROClinS2_CV.mat
load AUROCquadS2_CV.mat

% Load the integration within the subunit whitenoise analysis data
load AUROClinsubunits_CV.mat
load AUROCquadsubunits_CV.mat

% Load the isoresponse data
load RSSE_linearmodel_CV.mat % Robust regression
load RSSE_quadmodel_CV.mat


Within_subunit_NLI = []; % Storing the within subunit NLI
Isoresponse_NLI = []; % For storing the Isoresponse NLI
Whitenoise_NLI = []; % For storing the white noise NLI

JK_error_whitenoise = []; % For storing jacknife whitenoise error
JK_error_isoresponse = []; % For storing jackknife isorespone error
JK_error_withinsubunit = []; % For storing the cone isoresponse NLI error

for ii = 1:numel(AUROClin1) 
    
    % GLM error for both the subunits
    Error_lin_subunit1 = 1-AUROClin1{ii};
    Error_lin_subunit2 = 1-AUROClin2{ii};
    
    % GQM error for both the subunits
    Error_quad_subunit1 = 1-AUROCquad1{ii};
    Error_quad_subunit2 = 1-AUROCquad2{ii};
    
    % Calculating the median ratio of errors both the subunits in log scale
    median_NLI_subunit1 = log10(median(Error_lin_subunit1./ Error_quad_subunit1));
    median_NLI_subunit2 = log10(median(Error_lin_subunit2./ Error_quad_subunit2));
    
    % Storing the NLIs
    Within_subunit_NLI = [Within_subunit_NLI; median([median_NLI_subunit1 median_NLI_subunit2])];
    
    % White noise NLI 
    Error_quad = 1-(AUROCquadsubunits{ii});
    Error_lin = 1-(AUROClinsubunits{ii});
    Whitenoise_NLI = [Whitenoise_NLI; log10(median(Error_lin./Error_quad))];
    
    % Isoresponse NLI
    Isoresponse_NLI = [Isoresponse_NLI; log10(median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii}))];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Computing jackknife errors %%%%
    
    % Jackknife White noise error
    N1 = numel(Error_quad);
    X1 = log10(jackknife(@median, Error_lin./Error_quad));
    error_WN = sqrt(sum((X1-mean(X1)).^2)*(N1-1)/N1);
    JK_error_whitenoise = [JK_error_whitenoise; error_WN];
    
    % Jackknife Isoresponse error
    N2 = numel(RSSE_linearmodel{ii});
    X2 = log10(jackknife(@median, RSSE_linearmodel{ii}./RSSE_quadmodel{ii}));
    error_Iso = sqrt(sum((X2-mean(X2)).^2)*(N2-1)/N2);
    JK_error_isoresponse = [JK_error_isoresponse; error_Iso];
    
    % Jackknife cone isoresponse error 
    N3 = numel(Error_lin_subunit1);
    subunit_1_NLI = Error_lin_subunit1./ Error_quad_subunit1;
    subunit_2_NLI = Error_lin_subunit2./ Error_quad_subunit2;
    X3 = median([log10(jackknife(@median, subunit_1_NLI)) log10(jackknife(@median, subunit_2_NLI))],2);
    error_subunit = sqrt(sum((X3-mean(X3)).^2)*(N3-1)/N3);
    JK_error_withinsubunit = [JK_error_withinsubunit; error_subunit];
end


indices = [109 24 74];
figure(plot_counter); 
markersize= 6;

subplot(311); % Within subunit NLI 
errorbar(Within_subunit_NLI(LUMidx),abs(conewts_svd(3,LUMidx)),JK_error_withinsubunit(LUMidx),'horizontal', 'o','MarkerSize',markersize,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1], 'color', [0 0 0]); hold on;
errorbar(Within_subunit_NLI(DOidx),abs(conewts_svd(3,DOidx)),JK_error_withinsubunit(DOidx),'horizontal','o','MarkerSize',markersize,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1], 'color',[1 0 0]);
errorbar(Within_subunit_NLI(hardtoclassifyidx),abs(conewts_svd(3,hardtoclassifyidx)), JK_error_withinsubunit(hardtoclassifyidx),'o','horizontal','MarkerSize',markersize,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1], 'color',[0.5 0.5 0.5]);

errorbar(Within_subunit_NLI(indices(1)),abs(conewts_svd(3,indices(1))),JK_error_withinsubunit(indices(1)),'horizontal', 'o','MarkerSize',markersize,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 1 0], 'color', [0 0 0]); hold on;
errorbar(Within_subunit_NLI(indices(2)),abs(conewts_svd(3,indices(2))),JK_error_withinsubunit(indices(2)),'horizontal','o','MarkerSize',markersize,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 1 0], 'color',[1 0 0]);
errorbar(Within_subunit_NLI(indices(3)),abs(conewts_svd(3,indices(3))), JK_error_withinsubunit(indices(3)),'o','horizontal','MarkerSize',markersize,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0 1 0], 'color',[0.5 0.5 0.5]);
axis square; set(gca,'Tickdir','out','Xlim',[-0.04 0.08],'XTick',-0.04:0.02:0.08,'Ylim',[0 0.4],'YTick',0:0.1:0.4); 
xlabel('Cone signal NLI'); ylabel('S cone input');

subplot(312); % White noise NLI
errorbar(Whitenoise_NLI(LUMidx),abs(conewts_svd(3,LUMidx)),JK_error_whitenoise(LUMidx),'horizontal', 'o','MarkerSize',markersize,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1], 'color', [0 0 0]); hold on;
errorbar(Whitenoise_NLI(DOidx),abs(conewts_svd(3,DOidx)),JK_error_whitenoise(DOidx),'horizontal','o','MarkerSize',markersize,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1], 'color',[1 0 0]);
errorbar(Whitenoise_NLI(hardtoclassifyidx),abs(conewts_svd(3,hardtoclassifyidx)), JK_error_whitenoise(hardtoclassifyidx),'o','horizontal','MarkerSize',markersize,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1], 'color',[0.5 0.5 0.5]);

errorbar(Whitenoise_NLI(indices(1)),abs(conewts_svd(3,indices(1))),JK_error_whitenoise(indices(1)),'horizontal', 'o','MarkerSize',markersize,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 1 0], 'color', [0 0 0]); hold on;
errorbar(Whitenoise_NLI(indices(2)),abs(conewts_svd(3,indices(2))),JK_error_whitenoise(indices(2)),'horizontal','o','MarkerSize',markersize,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 1 0], 'color',[1 0 0]);
errorbar(Whitenoise_NLI(indices(3)),abs(conewts_svd(3,indices(3))), JK_error_whitenoise(indices(3)),'o','horizontal','MarkerSize',markersize,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0 1 0], 'color',[0.5 0.5 0.5]);
axis square; set(gca,'Tickdir','out','Xlim',[-0.04 0.14],'XTick',-0.04:0.02:0.14,'Ylim',[0 0.4],'YTick',0:0.1:0.4); 
xlabel('White noise NLI'); ylabel('S cone input');

subplot(313); % Isoresponse NLI
errorbar(Isoresponse_NLI(LUMidx),abs(conewts_svd(3,LUMidx)),JK_error_isoresponse(LUMidx),'horizontal', 'o','MarkerSize',markersize,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1], 'color', [0 0 0]); hold on;
errorbar(Isoresponse_NLI(DOidx),abs(conewts_svd(3,DOidx)),JK_error_isoresponse(DOidx),'horizontal','o','MarkerSize',markersize,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1], 'color',[1 0 0]);
errorbar(Isoresponse_NLI(hardtoclassifyidx),abs(conewts_svd(3,hardtoclassifyidx)), JK_error_isoresponse(hardtoclassifyidx),'o','horizontal','MarkerSize',markersize,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1], 'color',[0.5 0.5 0.5]);

errorbar(Isoresponse_NLI(indices(1)),abs(conewts_svd(3,indices(1))),JK_error_isoresponse(indices(1)),'horizontal', 'o','MarkerSize',markersize,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 1 0], 'color', [0 0 0]); hold on;
errorbar(Isoresponse_NLI(indices(2)),abs(conewts_svd(3,indices(2))),JK_error_isoresponse(indices(2)),'horizontal','o','MarkerSize',markersize,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 1 0], 'color',[1 0 0]);
errorbar(Isoresponse_NLI(indices(3)),abs(conewts_svd(3,indices(3))), JK_error_isoresponse(indices(3)),'o','horizontal','MarkerSize',markersize,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0 1 0], 'color',[0.5 0.5 0.5]);
axis square; set(gca,'Tickdir','out','Xlim',[-4 3],'XTick',-4:1:4,'Ylim',[0 0.4],'YTick',0:0.1:0.4); 
xlabel('Isoresponse NLI'); ylabel('S cone input');
set(gcf,'renderer', 'painters');
plot_counter = plot_counter + 1;

X = [LUMidx, DOidx, hardtoclassifyidx];
[r1,p1] = corr(Within_subunit_NLI(X),abs(conewts_svd(3,(X))'),'type','Spearman');
[r2,p2] = corr(Whitenoise_NLI(X),abs(conewts_svd(3,(X))'),'type','Spearman');
[r3,p3] = corr(Isoresponse_NLI(X),abs(conewts_svd(3,(X))'),'type','Spearman');

[r4,p4] = corr(Within_subunit_NLI(LUMidx),abs(conewts_svd(3,(LUMidx))'),'type','Spearman');
[r5,p5] = corr(Within_subunit_NLI(DOidx),abs(conewts_svd(3,(DOidx))'),'type','Spearman');
[r6,p6] = corr(Within_subunit_NLI(hardtoclassifyidx),abs(conewts_svd(3,(hardtoclassifyidx))'),'type','Spearman');

[r7,p7] = corr(Whitenoise_NLI(LUMidx),abs(conewts_svd(3,(LUMidx))'),'type','Spearman');
[r8,p8] = corr(Whitenoise_NLI(DOidx),abs(conewts_svd(3,(DOidx))'),'type','Spearman');
[r9,p9] = corr(Whitenoise_NLI(hardtoclassifyidx),abs(conewts_svd(3,(hardtoclassifyidx))'),'type','Spearman');

[r10,p10] = corr(Isoresponse_NLI(LUMidx),abs(conewts_svd(3,(LUMidx))'),'type','Spearman');
[r11,p11] = corr(Isoresponse_NLI(DOidx),abs(conewts_svd(3,(DOidx))'),'type','Spearman');
[r12,p12] = corr(Isoresponse_NLI(hardtoclassifyidx),abs(conewts_svd(3,(hardtoclassifyidx))'),'type','Spearman');

group = [ones(size(LUMidx)) 2*ones(size(DOidx)) 3*ones(size(hardtoclassifyidx))];
data = abs(conewts_svd([LUMidx'; DOidx'; hardtoclassifyidx'])); 
p13 = kruskalwallis(data,group,'off');
%% Figure 7: Downstream circuitry of simple and DO cells
% No code for this part
% The figure was provided by Greg


%% SUPPLEMENTAL FIGURES


%% Figure 2-Figure Supplement 1: Stimulus comparison and motivation for the 3 phases of the experiment

% Check the code in Calculating_effect_contrast.m file for stimulus
% comparison. The figure was generated using that code. 


%% Figure 2-Figure Supplement 2: Scatterplot of firing rates between the pixel WN and hyperpixel WN noise

if ~exist('plot_counter')
    plot_counter = 1;
end


% Loading all the files 
try 
    % Using the JDBC connection
    conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
    filename = fetch(conn,'SELECT filename FROM WNthresh');
    NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
    spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
    close(conn);
    filename = filename(strcmp(string(NTmode),"subunit"));
    NTmode = NTmode(strcmp(string(NTmode),"subunit"));
    spikeidx_NT = spikeidx_NT(strcmp(string(NTmode),"subunit"));

catch
    csv_filename = '/Users/abhishekde/Desktop/MatlabCode/Abhishek/CSV_PHPmyadmin_files/WNthresh.csv';
    [filename, NTmode, spikeIdx] = get_WNthreshdata_from_csvfile(csv_filename, 'subunit');
    spikeidx_NT = str2num(cell2mat(spikeIdx));
end


% Classifying cells into simple, DO and hardtoclassify cells
load conewts_svd.mat
load vals.mat

load S1RGB_svd.mat
load S2RGB_svd.mat
load angulardifferences_RGB.mat
anglebwvectors = angulardifference_RGB;
S1RGB = S1RGB_svd;
S2RGb = S2RGB_svd;
% SpatiallyOpponent = sum(sign(S1RGB).*sign(S2RGB),1)<3;
SpatiallyOpponent = anglebwvectors'>90;

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

LUMidx = LumIds_conewts;
DOidx = [ColorOpponentIds_conewts Sconedominated_conewts];
hardtoclassifyidx = [Other_conewts];
hardtoclassifyidx = [hardtoclassifyidx LUMidx(vals(LUMidx)>=95) DOidx(vals(DOidx)>=95)];
LUMidx = LUMidx(vals(LUMidx)<95);
DOidx = DOidx(vals(DOidx)<95);

% Considering only the spatially opponent subunits
LUMidx = LUMidx(SpatiallyOpponent(LUMidx));
DOidx = DOidx(SpatiallyOpponent(DOidx));
hardtoclassifyidx = hardtoclassifyidx(SpatiallyOpponent(hardtoclassifyidx));


% Matrix for storing the firing rates (mean and std)
mean_firing_rates = zeros(2, numel(filename));
std_firig_rates = zeros(2, numel(filename));
for ii= 1:numel(filename)
    fileofinterest = char(filename(ii,:));
    disp(fileofinterest);
    stro = nex2stro(findfile(fileofinterest));
    spikename = 'sig001a';%getSpikenum(stro);
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
    
    % Calculate Firing rates in pixel WN and hyperpixel WN
    for index = 1:2
        num_spikes =[];
        num_dur = [];
        for jj = mask_changes(1,index):mask_changes(2,index)
            nframes = stro.trial(jj,nframesidx);
            
            if nframes >0
                t_stimon = stro.trial(jj, stimonidx);
                num_spikes = [num_spikes; numel((stro.ras{jj,spikeidx}-t_stimon)*1000)];

                frametimes = linspace(0, nframes*msperframe, nframes)+(msperframe/2)';
                num_dur = [num_dur; range(frametimes)/1000];
            end

        end
        
        % Storing the mean and standard deviation of firing rates
        fr = num_spikes./num_dur;
        mean_firing_rates(index,ii) = mean(fr);
        %std_firig_rates(index,ii) = sem(fr);
        
        M = mean(fr);
        V = var(fr);
        n = numel(fr);
        std_firig_rates(index,ii) = sqrt(V./n./((log(10)*M).^2));
        
    end
    
end
    

% Comparing the Firing rates from pixel and hyperpixel WN 
FR_max = 1000;
FR_min = 0.3;

figure(plot_counter), hold on;
errorbar(log10(mean_firing_rates(1,LUMidx)), log10(mean_firing_rates(2,LUMidx)),std_firig_rates(2,LUMidx), std_firig_rates(2,LUMidx), std_firig_rates(1,LUMidx), std_firig_rates(1,LUMidx), 'o', 'MarkerFaceColor', [0 0 0],'MarkerEdgeColor',[1 1 1], 'color', [0 0 0]); 
errorbar(log10(mean_firing_rates(1,DOidx)), log10(mean_firing_rates(2,DOidx)),std_firig_rates(2,DOidx), std_firig_rates(2,DOidx),std_firig_rates(1,DOidx), std_firig_rates(1,DOidx), 'o', 'MarkerFaceColor', [1 0 0],'MarkerEdgeColor',[1 1 1], 'color', [1 0 0]); 
errorbar(log10(mean_firing_rates(1,hardtoclassifyidx)), log10(mean_firing_rates(2,hardtoclassifyidx)), std_firig_rates(2,hardtoclassifyidx), std_firig_rates(2,hardtoclassifyidx),  std_firig_rates(1,hardtoclassifyidx), std_firig_rates(1,hardtoclassifyidx),'o', 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1], 'color', [0.5 0.5 0.5]); 
plot(log10([FR_min FR_max]), log10([FR_min FR_max]), 'color', 'k')
set(gca,'Tickdir','out','Xlim',[log10(FR_min) log10(FR_max)],'Ylim',[log10(FR_min) log10(FR_max)], 'XTick', log10([0.3 1 10 100 1000]),'YTick', log10([0.3 1 10 100 1000]), 'XTickLabel', [0.3 1 10 100 1000],'YTickLabel', [0.3 1 10 100 1000]); 
xlabel('Pixel WN firing rates'); ylabel('Hyperpixel WN firing rates'); axis square; hold off;
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;




%% Figure 2-Figure Supplement 3: Cone weights plot
% Greg made this figure


%% Figure 2-Figure Supplement 4: An illustration of GQM and GLM analysis 

%% Figure 2-Figure Supplement 5: Relationship between NLI and PC1 

% Showing that it is possible to make a model neuron that does not have 
% a significant PC1 but appears nonlinear in the white noise NLI analysis.
 
% Simulation 1
% Two linear subunits, half rectified and multiplied.
% PC1 criterion will not be triggered, but white noise NLI will show the
% nonlinearity. This was a kind of nonlinearity we were concerned about.
 
close all; clearvars;


% Computing cone weights
% Forming the M matrix 
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;


n = 1000000;
nbins = 20;
stim = normrnd(0,1,n,6);
kernel1 = [.5 -.5 .2];
kernel2 = -kernel1;
conewts = inv(M)'*kernel1'; 
conewts = conewts/sum(abs(conewts));
 
% Simulation 1
lingen = stim*[[kernel1';0;0;0], [0;0;0;kernel2']];
resp = prod(max(lingen+1,0),2);
resp = (resp-min(resp))./(max(resp)-min(resp));
resp = resp>unifrnd(zeros(n,1),ones(n,1));
 
STA = resp'*stim;
[~,d] = eig(cov(stim(resp,:)));
figure; subplot(2,1,1);
plot(flipud(diag(d)),'ko','MarkerFaceColor','black')
set(gca, 'Ylim', [0 1.2]);
PC_orig = diag(d);

% Permutation test
PC_vals = [];
for iter=1:1000
    resp_mod = circshift(resp,randi(numel(resp)));
    STA_mod = resp_mod'*stim;
    [~,d1] = eig(cov(stim(resp_mod,:)));
    PCs = diag(d1);
    PC_vals = [PC_vals; PCs(end)];
end
PC_prctile1 = (find(sort(PC_vals)>PC_orig(end),1)/1000)*100;
PC_prctile1 = (find(sort(PC_vals)>PC_orig(end),1)/1000)*100;
if isempty(PC_prctile1)
    PC_prctile1 = 100;
end
 


clear allstim spikestim
allstim(:,1)=stim*[STA(1:3),0 0 0]';
allstim(:,2)=stim*[0 0 0 STA(4:6)]';
spikestim(:,1)=stim(resp,:)*[STA(1:3),0 0 0]';
spikestim(:,2)=stim(resp,:)*[0 0 0 STA(4:6)]';
 
bins = [nbins prctile(allstim(:),5) prctile(allstim(:),95)]';


[allstim_hist,~]=hist2(allstim,[bins bins]);
[spikestim_hist,~]=hist2(spikestim,[bins bins]);
subplot(2,1,2);
imagesc(spikestim_hist./allstim_hist*255); colormap(gray(255)); hold on;
contour(spikestim_hist./allstim_hist*255,'LineColor','w'); axis square;
set(gca,'Xtick',[],'Ytick',[]);

% Fit a GLM & GQM after the frames have been projected onto the subunits
mdllin1 =  fitglm(lingen,resp,'linear','Distribution','binomial','Link','logit');
mdlquad1 =  fitglm(lingen,resp,'quadratic','Distribution','binomial','Link','logit');
predlin1 = predict(mdllin1,lingen); % perdiction from GLM
predquad1 = predict(mdlquad1,lingen); % perdiction from GQM
Error_lin1 = 1-rocN(predlin1(resp),predlin1(~resp));
Error_quad1 = 1- rocN(predquad1(resp),predquad1(~resp)); 
WhiteNoise_NLI_1 = log10(Error_lin1/Error_quad1);


%
% Simulation 2
% Conversely: Two subunits with a halfwave rectified response to one color
% channel and a fullwave rectified response to another. Added linearly.
% This cell will have a significant PC1 but the GQM will not fit better
% than the GLM (because the fitting is in 2D).
lingen = [];
lingen(:,1) = stim(:,1).*kernel1(1)+(stim(:,2).*kernel1(2)).^2+stim(:,3).*kernel1(3); % nonlinear subfield
lingen(:,2) = stim(:,[4:6])*kernel2'; % linear subfield
resp = sum(lingen,2);
resp = (resp-min(resp))./(max(resp)-min(resp));
resp = resp>unifrnd(zeros(n,1),ones(n,1));
 
STA = resp'*stim;
[v,d] = eig(cov(stim(resp,:)));
figure; subplot(2,1,1);
plot(flipud(diag(d)),'ko','MarkerFaceColor','black')
PC_orig = diag(d);

% Permutation test
PC_vals = [];
for iter=1:1000
    resp_mod = circshift(resp,randi(numel(resp)));
    STA_mod = resp_mod'*stim;
    [~,d1] = eig(cov(stim(resp_mod,:)));
    PCs = diag(d1);
    PC_vals = [PC_vals; PCs(end)];
end
PC_prctile2 = (find(sort(PC_vals)>PC_orig(end),1)/1000)*100;
if isempty(PC_prctile2)
    PC_prctile2 = 100;
end
 
clear allstim spikestim
allstim(:,1)=stim*[STA(1:3),0 0 0]';
allstim(:,2)=stim*[0 0 0 STA(4:6)]';
spikestim(:,1)=stim(resp,:)*[STA(1:3),0 0 0]';
spikestim(:,2)=stim(resp,:)*[0 0 0 STA(4:6)]';
 
bins = [nbins prctile(allstim(:),5) prctile(allstim(:),95)]';
 
[allstim_hist,~]=hist2(allstim,[bins bins]);
[spikestim_hist,~]=hist2(spikestim,[bins bins]);
subplot(2,1,2);
imagesc(spikestim_hist./allstim_hist*255); colormap(gray(255)); hold on;
contour(spikestim_hist./allstim_hist*255,'LineColor','w'); axis square;
set(gca,'Xtick',[],'Ytick',[]); hold off;

% Fit a GLM & GQM after the frames have been projected onto the subunits
mdllin2 =  fitglm(lingen,resp,'linear','Distribution','binomial','Link','logit');
mdlquad2 =  fitglm(lingen,resp,'quadratic','Distribution','binomial','Link','logit');
predlin2 = predict(mdllin2,lingen); % perdiction from GLM
predquad2 = predict(mdlquad2,lingen); % perdiction from GQM
Error_lin2 = 1-rocN(predlin2(resp),predlin2(~resp));
Error_quad2 = 1- rocN(predquad2(resp),predquad2(~resp)); 
WhiteNoise_NLI_2 = log10(Error_lin2/Error_quad2);



%% Figure 3-Figure Supplement 1: Adaptive staircase procedure from an example cell
if ~exist('plot_counter')
    plot_counter = 1;
end

% Loading all the files 
try 
    % Using the JDBC connection
    conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
    filename = fetch(conn,'SELECT filename FROM WNthresh');
    NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
    spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
    close(conn);
    filename = filename(strcmp(string(NTmode),"subunit"));
    NTmode = NTmode(strcmp(string(NTmode),"subunit"));
    spikeidx_NT = spikeidx_NT(strcmp(string(NTmode),"subunit"));

catch
    csv_filename = '/Users/abhishekde/Desktop/MatlabCode/Abhishek/CSV_PHPmyadmin_files/WNthresh.csv';
    [filename, NTmode, spikeIdx] = get_WNthreshdata_from_csvfile(csv_filename, 'subunit');
    spikeidx_NT = str2num(cell2mat(spikeIdx));
end

stro = nex2stro(findfile(char(filename(24))));
global spikename maskidx spikeidx neurothreshidx nstixperside ngammasteps seedidx nframesidx correctidx
global fponidx stimoffidx stimonidx muidxs sigmaidxs basisvecidx weightsidx fpacqidx targetspikerateidx basisvecdiridx latencyidx reversalflagidx parentverticesidx
global msperframe ntrials maxT xx yy M linepredtol stepsizescale stepsize nreversals oogscale


spikename = 'sig001a';
spikeidx = strcmp(stro.sum.rasterCells(1,:),spikename);
maskidx = strcmp(stro.sum.rasterCells(1,:),'subunit_mask');
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
        raster_data = stro.ras(idxs1,spikeidx);
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
    raster_data = stro.ras(idxs1,spikeidx);
    num_dur =[];
    firing_rate = [];
    
    for ii = 1:size(raster_data,1)
        tmp = raster_data{ii} ;
        %t_offset=0;
        %spikes = tmp(tmp<stro.trial(idxs1(ii),stimoffidx) & tmp>stro.trial(idxs1(ii),stimonidx)+t_offset);
        spikes = sum(tmp<stro.trial(idxs1(ii),stimoffidx) & tmp>(stro.trial(idxs1(ii),stimonidx)+t_offset));
        num_dur = [num_dur; (stro.trial(idxs1(ii), stimoffidx)- stro.trial(idxs1(ii),stimonidx))-t_offset];
        firing_rate = [firing_rate; spikes/num_dur(end)];
    end
    
    subplot(2,2,4); plot(tmp_n,firing_rate,'-o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
    xlabel('Contrast'); ylabel('Firing rate '); axis square;
    set(gca,'Tickdir','out','XScale','log')
    hold off;
    
    subplot(2,2,3),plot(firing_rate,'-o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
    set(gca,'Tickdir','out','Ylim',[0 120],'YTick',[0 30 60 90 120]); axis square; xlabel('Trials'), ylabel('FR');
    
end


subplot(2,2,4); hold on; set(gca,'Xlim',[0.1 1],'XTick',[0.1 0.3 1],'Ylim',[0 120],'YTick',[0 30 60 90 120]);
subplot(2,2,3); hold on; plot([0 20],[num_targetspikerates num_targetspikerates],'k'); hold off;

plot_counter = plot_counter + 1;

%% Figure 3 - Figure Supplement 2: Distribution of baseline and target firing rates of DO, simple and NSNDO cells 

if ~exist('plot_counter')
    plot_counter = 1;
end

% Loading all the files 
try 
    % Using the JDBC connection
    conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
    filename = fetch(conn,'SELECT filename FROM WNthresh');
    NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
    spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
    close(conn);
    filename = filename(strcmp(string(NTmode),"subunit"));
    NTmode = NTmode(strcmp(string(NTmode),"subunit"));
    spikeidx_NT = spikeidx_NT(strcmp(string(NTmode),"subunit"));

catch
    csv_filename = '/Users/abhishekde/Desktop/MatlabCode/Abhishek/CSV_PHPmyadmin_files/WNthresh.csv';
    [filename, NTmode, spikeIdx] = get_WNthreshdata_from_csvfile(csv_filename, 'subunit');
    spikeidx_NT = str2num(cell2mat(spikeIdx));
end


load conewts_svd.mat
load vals.mat
load vals.mat
load S1RGB_svd.mat
load S2RGB_svd.mat
load angulardifferences_RGB.mat
anglebwvectors = angulardifference_RGB;
S1RGB = S1RGB_svd;
S2RGb = S2RGB_svd;
% SpatiallyOpponent = sum(sign(S1RGB).*sign(S2RGB),1)<3;
SpatiallyOpponent = anglebwvectors'>90;

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

LUMidx = LumIds_conewts;
DOidx = [ColorOpponentIds_conewts Sconedominated_conewts];
hardtoclassifyidx = [Other_conewts];
hardtoclassifyidx = [hardtoclassifyidx LUMidx(vals(LUMidx)>=95) DOidx(vals(DOidx)>=95)];
LUMidx = LUMidx(vals(LUMidx)<95);
DOidx = DOidx(vals(DOidx)<95);

% Considering only the spatially opponent subunits
LUMidx = LUMidx(SpatiallyOpponent(LUMidx));
DOidx = DOidx(SpatiallyOpponent(DOidx));
hardtoclassifyidx = hardtoclassifyidx(SpatiallyOpponent(hardtoclassifyidx));

% Checking how target firing rates of the DO, simple and unclassified cells
load baselineFRstats.mat % Baseline FR rates
load TFR.mat % Target firing rates


% Error plot representation
figure(plot_counter);
cell_idx = [DOidx LUMidx hardtoclassifyidx];
[~, indices] = sort(TFR(1,cell_idx));
indices = cell_idx(indices);
TFRprctile = [];
TFRzscore = []; % for storing the z-scores
TFRprctile_phase3 = [];
TFRzscore_phase3 = []; % for storing the z-scores
for iter = 1:numel(indices)
    disp(iter);
    
    ii = indices(iter);
    
    % Pulling out the respective file
    fileofinterest = char(filename(ii,:));
    stro = nex2stro(findfile(fileofinterest));
    spikename = 'sig001a';
    maskidx = strcmp(stro.sum.rasterCells(1,:),'subunit_mask');
    spikeidx = strcmp(stro.sum.rasterCells(1,:),spikename);
    basisvecidx = strcmp(stro.sum.rasterCells(1,:),'basis_vec');
    stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on');
    fpacqidx = strcmp(stro.sum.trialFields(1,:),'fpacq');
    
    mask_changes = [2];
    all_masks = stro.ras(:,maskidx);
    Fx = @(xi) any(isnan(xi)); % function that finds 'NaN' in a cell array
    inds = find(cellfun(Fx,stro.ras(:,basisvecidx))==0);
    if isempty(inds)
        inds = size(stro.trial,1)-1;
    end
    last_wntrial =  inds(1)-1;
    
    % Separately plotting the baseline firing rates for Phases 1 & 2, and Phase 3
    num_spikes =[];
    num_dur = [];
    for jj = 1:size(stro.ras,1)
        tmp = stro.ras{jj,spikeidx} ;
        if ~isnan((stro.trial(jj,stimonidx)- stro.trial(jj,fpacqidx)))
            spikes = tmp(tmp<stro.trial(jj,stimonidx) & tmp>stro.trial(jj,fpacqidx));
            num_spikes = [num_spikes; numel(spikes)];
            num_dur = [num_dur; (stro.trial(jj,stimonidx)- stro.trial(jj,fpacqidx))];
        end
    end
    %baselineFRstats{aa} = num_spikes./num_dur;
    
    baselineFR = num_spikes./num_dur;
    baselineFR_phase12 = baselineFR(1:last_wntrial);
    baselineFR_phase3 = baselineFR(last_wntrial+1:end);
    
    if ismember(ii,DOidx)
        c = [1 0 0];
    elseif ismember(ii, LUMidx)
        c = [0 0 0];
    else
        c = [0.5 0.5 0.5];
    end
    
    % Phase 1 and 2
    subplot(131); plot([iter iter],[prctile(baselineFR_phase12,5) prctile(baselineFR_phase12,95)],'color',c); hold on 
    plot(iter, TFR(1,ii), 'o', 'MarkerFaceColor',c,'MarkerEdgeColor',[1 1 1]);
    
    % Phase 3
    subplot(132); plot([iter iter],[prctile(baselineFR_phase3,5) prctile(baselineFR_phase3,95)],'color',c); hold on 
    plot(iter, TFR(1,ii), 'o', 'MarkerFaceColor',c,'MarkerEdgeColor',[1 1 1]);
    
    
     % All phases combined
    subplot(133); plot([iter iter],[prctile(baselineFR,5) prctile(baselineFR,95)],'color',c); hold on 
    plot(iter, TFR(1,ii), 'o', 'MarkerFaceColor',c,'MarkerEdgeColor',[1 1 1]);
    
    % Prctile for target firing rate as a function of baseline FRs from all
    % phases 
    TFRprctile = [TFRprctile; invprctile(baselineFR, TFR(1,ii))];
    TFRzscore = [TFRzscore; (TFR(1,ii)-mean(baselineFR))/std(baselineFR)];
    
    % Prctile for target firing rate as a function of baseline FRs from all
    % phases 
    TFRprctile_phase3 = [TFRprctile_phase3; invprctile(baselineFR_phase3, TFR(1,ii))];
    TFRzscore_phase3 = [TFRzscore_phase3; (TFR(1,ii)-mean(baselineFR_phase3))/std(baselineFR_phase3)];
  
    
end

subplot(131); axis square; set(gca,'Tickdir','out', 'Ylim',[0 100], 'YTick',0:20:100, 'Xlim', [0 98], 'XTick', [0 98]);
xlabel('Cell number'); ylabel('Firing rate'); title('Phase 1 & 2')

subplot(132); axis square; set(gca,'Tickdir','out', 'Ylim',[0 100], 'YTick',0:20:100, 'Xlim', [0 98], 'XTick', [0 98]);
xlabel('Cell number'); ylabel('Firing rate'); title('Phase 3')

subplot(133); axis square; set(gca,'Tickdir','out', 'Ylim',[0 100], 'YTick',0:20:100, 'Xlim', [0 98], 'XTick', [0 98]);
xlabel('Cell number'); ylabel('Firing rate'); title('Phase 1, 2 & 3')

set(gcf,'renderer','painters'); hold off
plot_counter = plot_counter + 1;

% Doing some stats on the target firing rates 
TFR_DO = TFR(1,DOidx);
TFR_DO = TFR_DO(:); TFR_DO = TFR_DO(TFR_DO>0);

TFR_LUM = TFR(1,LUMidx);
TFR_LUM = TFR_LUM(:); TFR_LUM = TFR_LUM(TFR_LUM>0);

TFR_htc = TFR(1, hardtoclassifyidx);
TFR_htc = TFR_htc(:); TFR_htc = TFR_htc(TFR_htc>0);

group = [ones(size(TFR_DO)); 2*ones(size(TFR_LUM)); 3*ones(size(TFR_htc))];
data = [TFR_DO; TFR_LUM; TFR_htc]; 
p2 = kruskalwallis(data,group,'off');



%% Figure 3?Figure Supplement 3- Population analysis of isoresponse curves: linear error vs. non-linear error 
if ~exist('plot_counter')
    plot_counter = 1;
end

load conewts_svd.mat
load vals.mat
load S1RGB_svd.mat
load S2RGB_svd.mat
S1RGB = S1RGB_svd;
S2RGB = S2RGB_svd;
anglebwvectors = angulardifference_RGB;
% SpatiallyOpponent = sum(sign(S1RGB).*sign(S2RGB),1)<3;
SpatiallyOpponent = anglebwvectors'>90;

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];


LUMidx = LumIds_conewts;
DOidx = [ColorOpponentIds_conewts Sconedominated_conewts];
hardtoclassifyidx = [Other_conewts];
hardtoclassifyidx = [hardtoclassifyidx LUMidx(vals(LUMidx)>=95) DOidx(vals(DOidx)>=95)];
LUMidx = LUMidx(vals(LUMidx)<95);
DOidx = DOidx(vals(DOidx)<95);

% Considering only the spatially opponent subunits
LUMidx = LUMidx(SpatiallyOpponent(LUMidx));
DOidx = DOidx(SpatiallyOpponent(DOidx));
hardtoclassifyidx = hardtoclassifyidx(SpatiallyOpponent(hardtoclassifyidx));

% Checking the correlation with non-linearity indices 
% Load the isoresponse data
load RSSE_linearmodel_CV.mat % Robust regression
load RSSE_quadmodel_CV.mat

% For storing median linear and quadratic errors and their standard errors
RSSEisoresp_lin_median = [];  % Isoresponse data-linear error
RSSEisoresp_quad_median = []; % Isoresponse data-quadratic error

% Jackknife error 
lin_std_error = []; 
quad_std_error = [];

for ii = 1:numel(RSSE_linearmodel)       
    
    RSSEisoresp_lin_median = [RSSEisoresp_lin_median; median(RSSE_linearmodel{ii})];
    RSSEisoresp_quad_median = [RSSEisoresp_quad_median; median(RSSE_quadmodel{ii})];
    
    % Computing estimate of errors
    % Linear error
    N1 = numel(RSSE_linearmodel{ii});
    X1 = jackknife(@median, RSSE_linearmodel{ii});
    lin_std_error = [lin_std_error; sqrt(sum((X1-mean(X1)).^2)*(N1-1)/N1)];
    
    % Quadratic error
    N2 = numel(RSSE_quadmodel{ii});
    X2 = jackknife(@median, RSSE_quadmodel{ii});
    quad_std_error = [quad_std_error; sqrt(sum((X2-mean(X2)).^2)*(N2-1)/N2)];
    
end

% Plotting the results for SVD based cone weight classification including the PC1 z-scores 
indices = [109 24 74]; % example cells
% figure(plot_counter);
% plot(RSSEisoresp_lin_median(hardtoclassifyidx),RSSEisoresp_quad_median(hardtoclassifyidx),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
% plot(RSSEisoresp_lin_median(LUMidx),RSSEisoresp_quad_median(LUMidx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
% plot(RSSEisoresp_lin_median(indices(1)),RSSEisoresp_quad_median(indices(1)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 1 0]);
% plot(RSSEisoresp_lin_median(DOidx),RSSEisoresp_quad_median(DOidx),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); plot([0.0001 10],[0.0001 10],'k');
% plot(RSSEisoresp_lin_median(indices(2)),RSSEisoresp_quad_median(indices(2)),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 1 0]);
% plot(RSSEisoresp_lin_median(indices(3)),RSSEisoresp_quad_median(indices(3)),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0 1 0]);


figure(plot_counter); hold on;
errorbar(RSSEisoresp_lin_median(hardtoclassifyidx), RSSEisoresp_quad_median(hardtoclassifyidx), lin_std_error(hardtoclassifyidx), lin_std_error(hardtoclassifyidx), quad_std_error(hardtoclassifyidx), quad_std_error(hardtoclassifyidx), 'o', 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1], 'color', [0.5 0.5 0.5]);
plot(RSSEisoresp_lin_median(indices(3)), RSSEisoresp_quad_median(indices(3)), 'o', 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor',[0 1 0]);
errorbar(RSSEisoresp_lin_median(LUMidx), RSSEisoresp_quad_median(LUMidx),lin_std_error(LUMidx), lin_std_error(LUMidx), quad_std_error(LUMidx), quad_std_error(LUMidx), 'o', 'MarkerFaceColor', [0 0 0],'MarkerEdgeColor',[1 1 1], 'color', [0 0 0]);
plot(RSSEisoresp_lin_median(indices(1)), RSSEisoresp_quad_median(indices(1)), 'o', 'MarkerFaceColor', [0 0 0],'MarkerEdgeColor',[0 1 0], 'MarkerSize', 8);
errorbar(RSSEisoresp_lin_median(DOidx), RSSEisoresp_quad_median(DOidx), lin_std_error(DOidx), lin_std_error(DOidx), quad_std_error(DOidx), quad_std_error(DOidx), 'o', 'MarkerFaceColor', [1 0 0],'MarkerEdgeColor',[1 1 1], 'color', [1 0 0]);
plot(RSSEisoresp_lin_median(indices(2)), RSSEisoresp_quad_median(indices(2)), 'o', 'MarkerFaceColor', [1 0 0],'MarkerEdgeColor',[0 1 0], 'MarkerSize', 8);
plot([0.0001 10],[0.0001 10],'k');
axis square; 
set(gca,'Tickdir','out','Xlim',[0.0001 10],'Ylim',[0.0001 10],'YScale','log','XScale','log','XTick',[0.0001 0.001 0.01 0.1 1 10],'YTick',[0.0001 0.001 0.01 0.1 1 10]); 
xlabel('Linear error'); ylabel('Quadratic error'); title('Isoresponse'); hold off;
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;



%% Figure 3?Figure Supplement 4 - Robustness of classification
if ~exist('plot_counter')
    plot_counter = 1;
end

load conewts_svd.mat
load S1RGB_svd.mat
load S2RGB_svd.mat
load angulardifferences_RGB.mat
load vals.mat
S1RGB = S1RGB_svd;
S2RGB = S2RGB_svd;
anglebwvectors = angulardifference_RGB;
% SpatiallyOpponent = sum(sign(S1RGB).*sign(S2RGB),1)<3;
SpatiallyOpponent = anglebwvectors'>90;
signS1 = sum(sign(S1RGB),1);
signS2 = sum(sign(S2RGB),1);
hardtoclassifyidx = 1:numel(signS1);

% Projecting onto a Luminance axis
% load fundamentals.mat % Loading the cone fundamentals 
% fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
load('T_vos1978_Y'); % Loading the luminosity function
load T_cones_smj.mat
luminosity = T_vos1978_Y;
cone_fundamentals = T_cones_smj;
u = cone_fundamentals'\luminosity'; % converting luminosity into cone mechanisms
u = u/norm(u); % Normalizing

conewts_mod = conewts_svd./repmat(sqrt(sum(conewts_svd.^2,1)),[3 1]); % Normalizing the cone weights
proj = abs(conewts_mod'*u); % Projecting cone weights onto the luminance cone mechanism 
proj = proj/max(proj);

% vals -> PC1 permuation test results 
LUMidx = find((proj>0.67 & vals<95 & SpatiallyOpponent') ==1)'; % Classifying cells
DOidx = find((proj<0.33 & vals<95 & SpatiallyOpponent') ==1)';
hardtoclassifyidx([LUMidx DOidx]) = [];
hardtoclassifyidx = hardtoclassifyidx(SpatiallyOpponent(hardtoclassifyidx));


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
plot(median(RSSEisoresp_medianofratios(LUMidx)),13,'v','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 100],'XTick',[0.1 1 10 100],'Ylim',[0 15],'YTick',[0 5 10 15]); ylabel('Count'); xlabel('Isoresponse NLI'); axis square; hold off;
subplot(322); histogram(RSSEisoresp_medianofratios(DOidx),logspace(-1,3,31),'FaceColor',[1 0 0],'EdgeColor',[1 1 1]); hold on;
plot(median(RSSEisoresp_medianofratios(DOidx)),8,'v','MarkerSize',8,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 100],'XTick',[0.1 1 10 100],'Ylim',[0 10],'YTick',[0 5 10]); ylabel('Count'); xlabel('Isoresponse NLI'); axis square; hold off;
subplot(326); histogram(RSSEisoresp_medianofratios(hardtoclassifyidx),logspace(-1,3,31),'FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1]); hold on;
plot(median(RSSEisoresp_medianofratios(hardtoclassifyidx)),15,'v','MarkerSize',8,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); 
set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 100],'XTick',[0.1 1 10 100],'Ylim',[0 20],'YTick',[0 5 10 15 20]); ylabel('Count'); xlabel('Isoresponse NLI'); axis square; hold off;


% Plotting the white noise NLIs on a cell-type basis
figure(plot_counter);
subplot(321); histogram(Whitenoise_NLI(DOidx),-0.02:0.005:0.1,'FaceColor',[1 0 0],'EdgeColor',[1 1 1]); hold on;
plot(median(Whitenoise_NLI(DOidx)),8,'v','MarkerSize',8,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-0.02 0.08],'XTick',-0.02:0.02:0.08,'Ylim',[0 10],'YTick',[0 5 10]); xlabel('Whitenoise NLI'); ylabel('Count'); axis square; hold off;
subplot(323); histogram(Whitenoise_NLI(LUMidx),-0.02:0.005:0.1,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); hold on;
plot(median(Whitenoise_NLI(LUMidx)),13,'v','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-0.02 0.08],'XTick',-0.02:0.02:0.08,'Ylim',[0 15],'YTick',[0 5 10 15]); xlabel('Whitenoise NLI'); ylabel('Count'); axis square; hold off;
subplot(325); histogram(Whitenoise_NLI(hardtoclassifyidx),-0.02:0.005:0.1,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1]); hold on;
plot(median(Whitenoise_NLI(hardtoclassifyidx)),15,'v','MarkerSize',8,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-0.02 0.08],'XTick',-0.02:0.02:0.08,'Ylim',[0 20],'YTick',[0 5 10 15 20]); xlabel('Whitenoise NLI'); ylabel('Count'); axis square; hold off;
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;

% Some stats on the whitenoise NLI
[p1,~] = ranksum(Whitenoise_NLI(LUMidx),Whitenoise_NLI(DOidx));
[p2,~] = ranksum(Whitenoise_NLI([LUMidx DOidx]),Whitenoise_NLI(hardtoclassifyidx));

% Some stats on the isoresponse NLI
[p3,~] = ranksum(log(RSSEisoresp_medianofratios(LUMidx)),log(RSSEisoresp_medianofratios(DOidx)));
[p4,~] = ranksum(log(RSSEisoresp_medianofratios([LUMidx DOidx])),log(RSSEisoresp_medianofratios(hardtoclassifyidx)));
[~,p5] = ttest2(log(RSSEisoresp_medianofratios([LUMidx DOidx])),log(RSSEisoresp_medianofratios(hardtoclassifyidx)));

%% Figure 5? Figure Supplement 1: Data from LGN cells RFs were obtained using spike-triggered averaging of checkerboard white noise

% Check the code LGN_analyses/RGB3D_analysesonLGN.m

%% Determining the some additional numbers

% Number of neurons from each Monkey
load Pangu.mat
load Maui.mat
NEURONS_MONKEY1 = numel(Maui);
NEURONS_MONKEY2 = numel(Pangu);

% RF locations 
load RF_loc.mat
RF_AMP = sqrt(sum(RF_loc.^2,2))/10;



