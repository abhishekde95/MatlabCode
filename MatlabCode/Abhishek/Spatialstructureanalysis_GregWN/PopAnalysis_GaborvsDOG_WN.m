% This is a variant of the code I wrote before and am trying it out for Greg's old WN data
% Comparing the Gabor and DOG fits for the WN data
% Author - Abhishek De, 1/19
close all; clearvars;
load fundamentals.mat 
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;
resize_fact = 1;
% Include Gun noise data for statistical tests in order to estimate the RF size
channels = 3;
NPOINTS = 65536;
resize_fact2 = 1;
spikename_options = ['sig001a'; 'sig001b'];
[filename_Lum, spikeIdx_Lum] = fnamesFromTxt2('Lum.txt');
[filename_LvsM, spikeIdx_LvsM] = fnamesFromTxt2('LvsM.txt');
[filename_ColorOpponent, spikeIdx_ColorOpponent] = fnamesFromTxt2('ColorOpponent.txt');
[filename_BY, spikeIdx_BY] = fnamesFromTxt2('BYcandidates.txt');
Input_List = [filename_Lum; filename_LvsM; filename_ColorOpponent; filename_BY];
spikeIdx = [spikeIdx_Lum; spikeIdx_LvsM; spikeIdx_ColorOpponent; spikeIdx_BY];
cellIds = [repmat({'Lum'},size(filename_Lum)); repmat({'LvsM'},size(filename_LvsM)); repmat({'ColorOpponent'},size(filename_ColorOpponent)); repmat({'BY'},size(filename_BY))];
LumIds = strcmp(cellIds,'Lum');
LvsMIds = strcmp(cellIds,'LvsM');
ColorOpponentIds = strcmp(cellIds,'ColorOpponent');
BYIds = strcmp(cellIds,'BY');

for ii = 1:numel(Input_List)
    disp(ii);
    filename = char(Input_List{ii}{1}); % acquiring the filename (1st column) from the List
    Output_List{ii,1} = filename;
    WN = {};
    for jj = 1:size(Input_List{ii},2)
        tmpstro = nex2stro(findfile(char(Input_List{ii}(jj))));
        if (isempty(WN))
            WN = tmpstro;
        else
            WN = strocat(WN, tmpstro);
        end
    end
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
    ridx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_r'));
    gidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_g'));
    bidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_b'));
    bkgndRGB = [mode(WN.trial(:,ridx)), mode(WN.trial(:,gidx)), mode(WN.trial(:,bidx))];
    bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
    bkgndlms = M*bkgndrgb;
   
    mask_changes = [2 size(WN.trial,1)];
    spikeidx = spikeIdx(ii);
    spikename = spikename_options(spikeidx,:); 
    maxT = 15;
    out_gun = getWhtnsStats(WN,maxT,'STCOVmex',{nstixperside^2,3,maxT},spikename);
    STAs_gun = out_gun{1}; STCs_gun = out_gun{2}; nspikes_gun = out_gun{3};
     
    % Code for Statistical testing begins here 
    CHI2CRIT = .95; % For flaging a stixel as significant (adding gun and cone noise z-scores, squared)
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
    tmpSTA = reshape(STAs_gun, [nstixperside^2 3 maxT]);
    Output_List{ii,2} = squeeze(sum(tmpSTA(:,:,peakframe),3));
    tmpSTA = permute(tmpSTA, [2 1 3]);
    try
        K=convhull(dominantclusteridxs(:,1), dominantclusteridxs(:,2));
        tmp = dominantclusteridxs(K,:);
        [x,y] = meshgrid(1:nstixperside,1:nstixperside);
        inRF = reshape(inpolygon(x(:),y(:),dominantclusteridxs(K,2),dominantclusteridxs(K,1)),[nstixperside, nstixperside]);
        Output_List{ii,3} = inRF;
    catch
        Output_List{ii,3} = zeros(nstixperside,nstixperside);
    end

    STAgunmat_pf = Output_List{ii,2};
    [u1,~,v1] = svd(STAgunmat_pf');
    m = reshape(v1(:,1),[nstixperside nstixperside]);
    Output_List{ii,4} = sign(m).*(abs(m).^(1.0));
    Output_List{ii,5} = u1(:,1); % storing the SVD derived spatial RF
    
    % storing latency
    Output_List{ii,6} = latency;
    
    % Storing the zscore based on spatiotemporal energy
    zscoremeans = STAs_gun./(sigmavect(1)/sqrt(nspikes_gun));
    Output_List{ii,7} = sum(zscoremeans.^2,1);
    
end

save Output_ListWN Output_List