%% Pretty much the same analysis present in Pop_Analysis_ConeVGun2.m
% This script deals mostly with gunnoise based STA
% Script written primarily for population analysis
% Author - Abhishek De, 11/16
% Cells used for this script were taken from ConeVGun2.txt
close all; clearvars;
Input_List = {'N021915003.nex';'N022715002.nex';'N040415002.nex';'N050415001.nex';'N041115001.nex';'N050215001-header-merge.nex';...
              'N050215002.nex';'N050415002.nex';'N060515002.nex';'N060615002.nex';'N060615005.nex';'N060815001.nex';...
              'N060815006.nex';'N061615001.nex';'N061715002.nex';'N061715004.nex';'N062315002.nex';'N063015002.nex';...
              'N063015003.nex';'N070115001.nex';'N070715005.nex';'N071015003.nex';'N071015004.nex';'N071615002.nex';...
              'N071815001.nex';'N072215003.nex';'N072615002.nex';'N072815001.nex';'N073015003.nex';'N073015004.nex';...
              'N080715002.nex';'N080715003.nex';'N081015001.nex';'N081115002.nex';'N081215002.nex';'N081515001.nex';...
              'N081515003.nex';'N081715001.nex';'N081715002.nex';'N082015002.nex';'P050116001.nex';'P052016002.nex';...
              'P052316005.nex';'P052816003.nex';'P052916001.nex';'P072716002.nex';'P072816001.nex';'P072816004.nex';...
              'P080116002.nex';'P080216003.nex';'P080316001.nex';'P080516003.nex';'P080816001.nex';'P080816003.nex';...
              'P081616002.nex';'P032916001.nex';'P033116001.nex';'P042816002.nex';'P051916001.nex';'P052316004.nex';...
              'P052816001.nex';'P072816002.nex';'P080216001.nex';'P080616001.nex';'P081016004.nex';'P081116001.nex';...
              'P081516002.nex';'P081716002.nex';'P081716003.nex';'P082616002.nex';'P081716001.nex';'P082216001.nex';...
              'P082216002.nex';'P081516002.nex';'M113016005.nex';'M113016006.nex';'M120716002.nex';'M120916002.nex';...
              'M121716002.nex';'M121816002.nex';'M121816003.nex';'M121816006.nex';'M121916002.nex';'M121916003.nex';...
              'M120116002.nex';'M120216001.nex';'M120216002.nex';'M120916003.nex';'M121216001.nex';'M121616004.nex';...
              'M121916001.nex';'M122016003.nex';'M122116001.nex';'M122316002.nex';'M122416004.nex';'M122616001.nex';...
              'M122916001.nex';'M123116004.nex';'M010117004.nex';'M010217001.nex';'M010317003.nex';'M010317004.nex';...
              'M010617003.nex';'M010917001.nex';'M122416001.nex';'M122416005.nex';'M122616002.nex';'M122516003.nex';...
              'M122516001.nex';'M122516004.nex';'M010117001.nex';'M010317001.nex';'M010317002.nex';'M010617002.nex'};
Output_List = cell(numel(Input_List),9); % filename, conenoise based STA, gunnoise based STA, RF,
% Cone Weights, Orange-Cyan/Lime_magenta/Luminance, Cone Weights bas   ed on the peakframe, Peakframe eigenveactor from SVD,gabor fit
% Constructing the M matrix
load fundamentals.mat 
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;
resize_fact = 7;

% Include Gun noise data for statistical tests in order to estimate the RF size
channels = 3;

for ii = 1:numel(Input_List)
    filename = Input_List{ii}; % acquiring the filename from the List
    Output_List{ii,1} = filename;
    disp(filename);
    WN=nex2stro(findfile(filename));
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
    basisvecidx = strcmp(WN.sum.rasterCells(1,:),'basis_vec');
    anlgStartTimeidx = find(strcmp(WN.sum.rasterCells(1,:),'anlgStartTime'));
    eyestart_t = [WN.ras{:,anlgStartTimeidx}]';
    eyesampperiod = 1/WN.sum.analog.storeRates{1};
    gammaTable = WN.sum.exptParams.gamma_table;
    gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
    gammaTable1 = interp1(linspace(0,255,256),gammaTable,linspace(0,255,65536), 'spline');
    invgamma = InvertGamma(gammaTable, 0);
    
    % Getting the background rgb/lms
    ridx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_r'));
    gidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_g'));
    bidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_b'));
    bkgndRGB = [mode(WN.trial(:,ridx)), mode(WN.trial(:,gidx)), mode(WN.trial(:,bidx))];
    bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
    bkgndlms = M*bkgndrgb;
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
    
    spikename = getSpikenum(WN,'first'); % Getting the spikes present in the first channel
    spikeidx = find(strcmp(WN.sum.rasterCells(1,:),spikename));
    maxT = 10;
    tmpstro_gun = WN; % creating a temporary stro that will be used for analysing stimulus presented in gun space
    tmpstro_gun.ras(mask_changes(2,1)+1:end,:) = [];
    tmpstro_gun.trial(mask_changes(2,1)+1:end,:) = [];
    out_gun = getWhtnsStats(tmpstro_gun,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, spikename);
    STAs_gun = out_gun{1}; STCs_gun = out_gun{2}; nspikes_gun = out_gun{3};
     
    % Code for Statistical testing begins here 
    CHI2CRIT = .95; % For flaging a stixel as significant (adding gun and cone noise z-scores, squared)
    s_gun = std(STAs_gun(:,1));
    STAs_gun_z = STAs_gun./s_gun;
    
    % Spatial map
    grandz = zeros([nstixperside nstixperside]);
    maxzs = [];
    for i = 1:maxT
        tmp_cone = reshape(STAs_gun_z(:,i),[nstixperside nstixperside 3]);
        tmp_gun = reshape(STAs_gun_z(:,i),[nstixperside nstixperside 3]); % This is the only place in the code where I use data from gun noise
        grandz = grandz+sum(tmp_gun.^2,3);
        maxzs(i) = sum(sum(sum(tmp_gun.^2,3)));
    end
    peakframe = max(maxzs) == maxzs;
    crit = chi2inv(CHI2CRIT,channels*maxT); % 3 color channels 
    L = grandz > crit;
    
    % Now get largest contiguous block
    [i,j] = ind2sub(size(L),find(L));
    ij = [i,j];
    T = clusterdata(ij,'linkage','single','distance','euclidean','criterion','distance','cutoff',sqrt(2));
    clusternmembers = [];
    for k =1:max(T)
        clusternmembers(k) = sum(T == k);
    end
    dominantcluster = find(clusternmembers == max(clusternmembers));
    clustermat = zeros(nstixperside, nstixperside);
    clustermat(sub2ind(size(clustermat),ij(T==dominantcluster(1),1),ij(T==dominantcluster(1),2))) = 1;
    % Then get convex hull
    dominantclusteridxs = ij(T==dominantcluster(1),:);
    K=convhull(dominantclusteridxs(:,1), dominantclusteridxs(:,2));
    tmp = dominantclusteridxs(K,:);
    [x,y] = meshgrid(1:nstixperside,1:nstixperside);
    inRF = reshape(inpolygon(x(:),y(:),dominantclusteridxs(K,2),dominantclusteridxs(K,1)),[nstixperside, nstixperside]);
    Output_List{ii,4} = inRF;
    whichpix = find(inRF);
    tmpSTA = reshape(STAs_gun, [nstixperside^2 3 maxT]);
    Output_List{ii,3} = tmpSTA(:,:,peakframe);
    tmpSTA = permute(tmpSTA, [2 1 3]);
    tmpSTA1 = tmpSTA(:,whichpix,:);
    STAgunmat = reshape(tmpSTA1,[3 length(whichpix)*maxT]);
     
    % Acquiring the cone weights from the STA
    [u,~,v] = svd(STAgunmat);
    if (sum(v(:,1)) < 0)
        u = -u;
    end
    coneweights = inv(M')*u(:,1);
    coneweights_gun = coneweights./sum(abs(coneweights));
    Output_List{ii,5} = coneweights_gun;
    
    % Determining if the neuron is Orange-Cyan/Lime_magenta/Luminance
    l = coneweights_gun(1); m = coneweights_gun(2); s = coneweights_gun(3);
    if (sign(s)*sign(m) == 1) && (sign(s)*sign(l)== -1)
        Output_List{ii,6} = 'Orange-Cyan';
    elseif (sign(s)*sign(m) == -1) && (sign(s)*sign(l)== 1)
        Output_List{ii,6} = 'Lime-Magenta';
    elseif (sign(s)*sign(m) == 1) && (sign(s)*sign(l)== 1)
        Output_List{ii,6} = 'Luminance';
    elseif (sign(s)*sign(m) == -1) && (sign(s)*sign(l)== -1)
        Output_List{ii,6} = 'Blue_Yellow';
    end
    
    % Now calculating the cone weights based on all the pixels present in peakframe
    STAgunmat_pf = squeeze(tmpSTA(:,:,peakframe));
    [u1,~,v1] = svd(STAgunmat_pf);
    coneweights_pf = inv(M')*u1(:,1);
    coneweights_gun_pf = coneweights_pf./sum(abs(coneweights_pf));
    if sign(coneweights_gun_pf(2)) * sign(coneweights_gun(2)) == -1
        Output_List{ii,7} = -1*coneweights_gun_pf;
    else
        Output_List{ii,7} = coneweights_gun_pf;
    end
    Output_List{ii,8} = reshape(v1(:,1),[nstixperside nstixperside]);
    [~,Output_List{ii,9}] = gaborfit_AD(imresize(Output_List{ii,8},resize_fact));
    
end
save('Output_List_Gun.mat','Output_List');

%% Plotting the conenoise STA, gunnoise STA, RF and cone weights 
kk = 1;
R = size(Output_List,1);
num_rows = 10; % Number of cells in a figure
C = 7;
resize_fact2 = 1;
for ii = 1:R
    figure(kk);
    
    % gunnoise based STA
    tmp_vec_gun = Output_List{ii,3};
    normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
    im = normfactor*Output_List{ii,3} + 0.5;
    im = reshape(im,[nstixperside nstixperside 3]);
    subplot(num_rows,C,((ii-(kk-1)*num_rows)-1)*C+1); image(im); set(gca,'XTick',[],'YTick',[]);
    if mod(ii,num_rows) ==1 
        title('Gun');
    end
    
    % RF size
    subplot(num_rows,C,((ii-(kk-1)*num_rows)-1)*C+2); imagesc(Output_List{ii,4});set(gca,'XTick',[],'YTick',[]);
    if mod(ii,num_rows) ==1 
        title('RF');
    end
    
    % Cone Weights
    subplot(num_rows,C,((ii-(kk-1)*num_rows)-1)*C+3); bar([Output_List{ii,5}, Output_List{ii,7}]);
    set(gca,'XTick',[1 2 3],'XTickLabel',{'L','M','S'});
    if mod(ii,num_rows) ==1 
        title('Cone Wts');
    end
    
    % Grayscale Image derived from Peak frame
    subplot(num_rows,C,((ii-(kk-1)*num_rows)-1)*C+4); imagesc(imresize(Output_List{ii,8},resize_fact2)); set(gca,'XTick',[],'YTick',[]);
    if mod(ii,num_rows) ==1 
        title('Peakframe');
    end
    
    % Gabor fit 
    subplot(num_rows,C,((ii-(kk-1)*num_rows)-1)*C+5); imagesc(Output_List{ii,9}); set(gca,'XTick',[],'YTick',[]);
    if mod(ii,num_rows) ==1 
        title('Gabor Fit');
    end
    
    % 2-D Fourier Transform 
    imf = Output_List{ii,8};
    imf = log(abs(fftshift(fft2(imresize(imf,resize_fact2)))));
    subplot(num_rows,C,((ii-(kk-1)*num_rows)-1)*C+6); imagesc(imf); set(gca,'XTick',[],'YTick',[]);
    if mod(ii,num_rows) ==1 
        title('DFT');
    end
    
    % resized STA
    im1 = cat(3,imresize(im(:,:,1),resize_fact),imresize(im(:,:,2),resize_fact),imresize(im(:,:,3),resize_fact));
    subplot(num_rows,C,((ii-(kk-1)*num_rows)-1)*C+7); image(im1); set(gca,'XTick',[],'YTick',[]);
    if mod(ii,num_rows) ==1 
        title('R-STA');
    end
    
    if mod(ii,num_rows) == 0
        kk = kk + 1;
    end
end

%% Sorting the cell classes as Orange-Cyan, Lime-Magenta and Luminance and
% plotting their RF size and cone weights
RF_size = [];
Orange_Cyan_idx = find(strcmp(Output_List(:,6),'Orange-Cyan')==1);
Lime_Magenta_idx = find(strcmp(Output_List(:,6),'Lime-Magenta')==1);
Luminance_idx = find(strcmp(Output_List(:,6),'Luminance')==1);
BY_idx = find(strcmp(Output_List(:,6),'Blue_Yellow')==1);
coneweights = [];
for ii = 1:R
    tmp_RF = Output_List{ii,4};
    RF_size = [RF_size; sum(tmp_RF(:))];
    coneweights = [coneweights; Output_List{ii,5}'];
end
figure(kk+1); subplot(221); plot(ones(numel(Orange_Cyan_idx),1),RF_size(Orange_Cyan_idx),'.b','MarkerSize',25); hold on;
plot(2*ones(numel(Lime_Magenta_idx),1),RF_size(Lime_Magenta_idx),'.r','MarkerSize',25);
plot(3*ones(numel(Luminance_idx),1),RF_size(Luminance_idx),'.g','MarkerSize',25); 
plot(4*ones(numel(BY_idx),1),RF_size(BY_idx),'.m','MarkerSize',25); 
set(gca,'XTick',[1 2 3 4],'XTickLabel',{'Orange-Cyan','Lime-Magenta','Luminance','BY'}); ylabel('RF size'); title('RF sizes'); xlim([0.5 4.5]); hold off;

% Plotting the proportion of neurons in a scatter plot
figure(kk+1); subplot(222);scatter(coneweights(Orange_Cyan_idx,1)-coneweights(Orange_Cyan_idx,2),coneweights(Orange_Cyan_idx,3),'b', 'filled'); hold on;
scatter(coneweights(Lime_Magenta_idx,1)-coneweights(Lime_Magenta_idx,2),coneweights(Lime_Magenta_idx,3),'r','filled');
scatter(coneweights(Luminance_idx,1)-coneweights(Luminance_idx,2),coneweights(Luminance_idx,3),'g','filled');
scatter(coneweights(BY_idx,1)-coneweights(BY_idx,2),coneweights(BY_idx,3),'m','filled');
xlabel('L-M'); ylabel('S'); title('Isoluminant plane');hold off;

% Plotting the cone weights in 3 dimensional L,M,S space
figure(kk+1); subplot(223);scatter3(coneweights(Orange_Cyan_idx,1),coneweights(Orange_Cyan_idx,2),coneweights(Orange_Cyan_idx,3),'b', 'filled'); hold on;
scatter3(coneweights(Lime_Magenta_idx,1),coneweights(Lime_Magenta_idx,2),coneweights(Lime_Magenta_idx,3),'r','filled');
scatter3(coneweights(Luminance_idx,1),coneweights(Luminance_idx,2),coneweights(Luminance_idx,3),'g','filled');
scatter3(coneweights(BY_idx,1),coneweights(BY_idx,2),coneweights(BY_idx,3),'m','filled');
xlabel('L'); ylabel('M'); zlabel('S'); title('Cone constrast space'); hold off; kk = kk+1;