%% Pretty much the same analysis present in Single_cell_Analysis.m
% Script written primarily for population analysis
% Author - Abhishek De, 11/16
% Cells used for this script were taken from ConeVGun2.txt
close all; clearvars;
Input_List = {'K031308006.nex';'K031708002.nex';'K032608005.nex';'K032708003.nex';'K033108004.nex';'K040108001.nex';...
              'K040108002.nex';'K040108003.nex';'K040208001.nex';'K040208003.nex';'K040708003.nex';'K040808001.nex';...
              'K043008001.nex';'K043008002.nex';'K050108003.nex';'K050508004.nex';'K051408003.nex';'K052108002.nex';...
              'K052308001.nex';'K053008003.nex';'K102808004.nex';'K110408006.nex';'K110708002.nex';'K111208001.nex';...
              'K121508002.nex';'K121508002.nex';'K020609002.nex';'K022309001.nex';'S022310008.nex';'S090210009.nex'};
Output_List = cell(numel(Input_List),6); % filename, conenoise based STA, gunnoise based STA, RF, Cone Weights, Orange-Cyan/Lime_magenta/Luminance
% Constructing the M matrix
load fundamentals.mat 
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

% Include Gun noise data for statistical tests in order to estimate the RF size
channels = 6;

for ii = 1:numel(Input_List)
    filename = Input_List{ii}; % acquiring the filename from the List
    Output_List{ii,1} = filename;
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
    
    spikename = getSpikenum(WN,'first'); % Getting the spikes present in the first channel
    spikeidx = find(strcmp(WN.sum.rasterCells(1,:),spikename));
    maxT = 10;
    Lgunnoise = WN.trial(:,noisetypeidx) == 1;
    Lconenoise = WN.trial(:,noisetypeidx) == 2;
    tmpstro_cone = WN; % creating a temporary stro that will be used for analysing stimulus presented in cone space
    tmpstro_gun = WN; % creating a temporary stro that will be used for analysing stimulus presented in gun space
    tmpstro_cone.ras(Lgunnoise,:) = [];
    tmpstro_cone.trial(Lgunnoise,:) = [];
    tmpstro_gun.ras(Lconenoise,:) = [];
    tmpstro_gun.trial(Lconenoise,:) = [];
    out_gun = getWhtnsStats(tmpstro_gun,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, spikename);
    out_cone = getWhtnsStats(tmpstro_cone,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, spikename);
    STAs_gun = out_gun{1}; STCs_gun = out_gun{2}; nspikes_gun = out_gun{3};
    STAs_cone = out_cone{1}; STCs_cone = out_cone{2}; nspikes_cone = out_cone{3};
    
    % Code for Statistical testing begins here 
    CHI2CRIT = .95; % For flaging a stixel as significant (adding gun and cone noise z-scores, squared)
    s_cone = std(STAs_cone(:,1));
    STAs_cone_z = STAs_cone./s_cone;
    s_gun = std(STAs_gun(:,1));
    STAs_gun_z = STAs_gun./s_gun;
    
    % Spatial map
    grandz = zeros([nstixperside nstixperside]);
    maxzs = [];
    for i = 1:maxT
        tmp_cone = reshape(STAs_cone_z(:,i),[nstixperside nstixperside 3]);
        tmp_gun = reshape(STAs_gun_z(:,i),[nstixperside nstixperside 3]); % This is the only place in the code where I use data from gun noise
        grandz = grandz+sum(tmp_gun.^2,3)+sum(tmp_cone.^2,3);
        maxzs(i) = sum(sum(sum(tmp_gun.^2,3)+sum(tmp_cone.^2,3)));
    end
    peakframe = max(maxzs) == maxzs;
    crit = chi2inv(CHI2CRIT,channels*maxT); % 6 = 3 color channels * 2 stimulus sets
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
    clustermat(sub2ind(size(clustermat),ij(T==dominantcluster,1),ij(T==dominantcluster,2))) = 1;
    % Then get convex hull
    dominantclusteridxs = ij(T==dominantcluster,:);
    K=convhull(dominantclusteridxs(:,1), dominantclusteridxs(:,2));
    tmp = dominantclusteridxs(K,:);
    [x,y] = meshgrid(1:nstixperside,1:nstixperside);
    inRF = reshape(inpolygon(x(:),y(:),dominantclusteridxs(K,2),dominantclusteridxs(K,1)),[nstixperside, nstixperside]);
    Output_List{ii,4} = inRF;
    whichpix = find(inRF);
    tmpSTA_gun = reshape(STAs_gun, [nstixperside^2 3 maxT]); Output_List{ii,3} = tmpSTA_gun(:,:,peakframe);
    tmpSTA = reshape(STAs_cone, [nstixperside^2 3 maxT]); Output_List{ii,2} = tmpSTA(:,:,peakframe);
    tmpSTA = permute(tmpSTA, [2 1 3]);
    tmpSTA = tmpSTA(:,whichpix,:);
    STAconemat = reshape(tmpSTA,[3 length(whichpix)*maxT]);
    conesigmas = WN.trial(WN.trial(:,noisetypeidx) == 2, sigmaidxs)/1000;  % cone excitation units
    
    % Acquiring the cone weights from the STA
    [u,s,v] = svd(STAconemat);
    u(:,1) = u(:,1)./mode(abs(conesigmas))';  % Think about this line
    coneweights_cone = u(:,1)./sum(abs(u(:,1)));
    Output_List{ii,5} = coneweights_cone;
    
    % Determining if the neuron is Orange-Cyan/Lime_magenta/Luminance
    l = coneweights_cone(1); m = coneweights_cone(2); s = coneweights_cone(3);
    if (sign(s)*sign(m) == 1) && (sign(s)*sign(l)== -1)
        Output_List{ii,6} = 'Orange-Cyan';
    elseif (sign(s)*sign(m) == -1) && (sign(s)*sign(l)== 1)
        Output_List{ii,6} = 'Lime-Magenta';
    elseif (sign(s)*sign(m) == 1) && (sign(s)*sign(l)== 1)
        Output_List{ii,6} = 'Luminance';
    elseif (sign(s)*sign(m) == -1) && (sign(s)*sign(l)== -1)
        Output_List{ii,6} = 'Blue_Yellow';
    end
end
save('Output_List_ConeVGun2.mat','Output_List');

%% Plotting the conenoise STA, gunnoise STA, RF and cone weights 
kk = 1;
R = size(Output_List,1);
num_rows = 8;
C = 4;
for ii = 1:R
    figure(kk);
    % conenoise based STA 
    tmp_vec_cone = Output_List{ii,2};
    normfactor = 0.5/(max(abs(tmp_vec_cone(:)))+0.01);
    im = normfactor*Output_List{ii,2} + 0.5;
    im = reshape(im,[nstixperside nstixperside 3]);
    subplot(num_rows,C,((ii-(kk-1)*num_rows)-1)*C+1); image(im); set(gca,'XTick',[],'YTick',[]);
    if mod(ii,num_rows) == 1 
        title('Cone');
    end
    
    % gunnoise based STA
    tmp_vec_gun = Output_List{ii,3};
    normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
    im = normfactor*Output_List{ii,3} + 0.5;
    im = reshape(im,[nstixperside nstixperside 3]);
    subplot(num_rows,C,((ii-(kk-1)*num_rows)-1)*C+2); image(im); set(gca,'XTick',[],'YTick',[]);
    if mod(ii,num_rows) == 1 
        title('Gun');
    end
    
    % RF size
    subplot(num_rows,C,((ii-(kk-1)*num_rows)-1)*C+3); imagesc(Output_List{ii,4}); set(gca,'XTick',[],'YTick',[]);
    if mod(ii,num_rows) ==1 
        title('RF');
    end
    
    % Cone Weights
    subplot(num_rows,C,((ii-(kk-1)*num_rows)-1)*C+4); bar(Output_List{ii,5});
    set(gca,'XTick',[1 2 3],'XTickLabel',{'L','M','S'});
    if mod(ii,num_rows) == 1 
        title('Cone Wts');
    end
    if mod(ii,num_rows) == 0
        kk = kk + 1;
    end
end

% Sorting the neuronal classes as Orange-Cyan, Lime-Magenta and Luminance and plotting their RF size and cone weights
RF_size = [];
Orange_Cyan_idx = find(strcmp(Output_List(:,end),'Orange-Cyan')==1);
Lime_Magenta_idx = find(strcmp(Output_List(:,end),'Lime-Magenta')==1);
Luminance_idx = find(strcmp(Output_List(:,end),'Luminance')==1);
BY_idx = find(strcmp(Output_List(:,end),'Blue_Yellow')==1);
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