
% IsoSamp population analyses
% Section 0) Constants that are used by > 1 cell
% Section 1) d' as a function of RF position, color, and TF. First try.
% Section 2) d' for individual neurons as a function of eccentricity.
% Section 3) Calculating d' of an ideal observer of the cones with access
% to the cones within a midget cell RF. Comparing to single neuron d'.
% Section 4) Calculating d' of an ideal observer of a population of neurons
% and comparing to the monkey's d' (1.27)
% Section 4.1) Run this after section 4 to get the cone model d' (assuming
% an ideal observer with access to all of the cones affected by the
% stimulus).
% Section 4.2) Comparing the sensitivity of the cone model to the behavioral
% model (ideal observer analysis across all dimensions)
% Section 4.3) Comparing the sensitivity of the cone model to the behavioral
% model (just eccentricity).
% Section 5) Cone weights from white noise files.
% Section 6) Comparing stimuli used in the IsoSamp experiments to
% psychophysical thresholds obtained from the most recent models.
% Section 7) Looking at the noise spectrum of neurons 
% Section 8) Looking at SNR using a K-nearest neighbors algorithm using the
% Victor and Purpura spike distance metric
% Section 9) Comparing cone weights from white noise to IsoSamp d-prime
% Section 10) looking at responses around the time of fixational eye
% movements
% Section 11) Checking the robustness of the efficiency analysis 
% (relative sensitivity of cones and LGN cells) 
% Section 12) Sliding window d' analysis for LGN cells
% Section 13) Single cycle d' analysis
% Section 14) d' and spike rate as a function of window length
% Section 15) d' based on spike counts
% Section 16) Comparing d' to parametric estimate of ROC to non-parametric
% estimate of ROC
% Section 17) Comparing parametric and nonparametric d's
% Section 18) Trying to figure out what values of Q and k work best for the
% nonparametric d' analysis
%%
% Section 0: Constants
CELLTYPE = 'M'; % 'M' for magno or 'P' for parvo
MONKEY = 'Utu';
if CELLTYPE == 'P'
    RFSIZESOURCE = 'Watson'; % "Perry et al."(M/P) "Dacey and Petersen"(P), "Derrington and Lennie"(M/P) "Croner and Kaplan"(M/P) "Watson" (P)
else
    RFSIZESOURCE = 'Derrington and Lennie';
end
MMPERDEG = 0.223; % mm/deg (Perry and Cowey 1985)
DEGPERMM = 1/MMPERDEG; % deg/mm
DPRIMEMETHOD = 1;

% Convention is to express RFs as (equivalent) eccentricity on the *temporal*
% retina (where acuity is not unusually high).
% ecc_to_diam_deg is function that takes eccentricity as as argument (in
% temporal retina equivalent DVA) and returns the RF "diameter"
% which is *two* standard deviations of a Gaussian fit to the center 
% (containing 68% of the probability mass in one dimension and 39% in two).
clear ecc_to_diam_deg
if strcmp(RFSIZESOURCE,'Derrington and Lennie')
    if CELLTYPE == 'M'
        ecc_to_diam_deg = @(rf_r_deg) 10.^(-1.2459+0.0345*rf_r_deg); % temporal retina equivalent
    elseif CELLTYPE == 'P'
        ecc_to_diam_deg = @(rf_r_deg) 10.^(-1.3322+0.0265*rf_r_deg);
    end
elseif strcmp(RFSIZESOURCE, 'Dacey and Petersen')
    if CELLTYPE == 'P'
        ecc_to_diam_deg = @(rf_r_deg) 10.^(-1.6192+0.0457*rf_r_deg); %  Assuming their DF diameters were 1 SD. See RGC_RF_sizes, section 5
    end
elseif strcmp(RFSIZESOURCE, 'Perry et al.')
    if CELLTYPE == 'M'
        ecc_to_diam_deg = @(rf_r_deg) 0.0673+0.0165*rf_r_deg;
    elseif CELLTYPE == 'P'
        ecc_to_diam_deg = @(rf_r_deg) 10.^(-1.7174+0.0434*rf_r_deg);
    end
    % I need to multiply Croner and Kaplan's rc values by 1/sqrt(2) to get
    % them into stds (X eccentricity should be multiplied by 0.8 to average
    % the sizes of temporal RFs (which they report) and nasal RFs (which 
    % are equivalent to temporal RFs if you go out another factor of 1/.61 = 1.64.)
elseif strcmp(RFSIZESOURCE, 'Croner and Kaplan')
    if CELLTYPE == 'M'
        ecc_to_diam_deg = @(rf_r_deg)(0.0121*rf_r_deg+0.0761); % 8/6/18. See RGC_RF_sizes, section 3.
    elseif CELLTYPE == 'P'
        ecc_to_diam_deg = @(rf_r_deg)(0.00350*rf_r_deg+0.0450);% 8/6/18. See RGC_RF_sizes, section 4
    end
    elseif strcmp(RFSIZESOURCE, 'Watson')
    if CELLTYPE == 'P'
        a = 0.9729; % Table 1
        r2 = 1.084; % Table 1
        re = 7.633; % Table 1
        dc_0 = 14804.6; % Cone density of fovea
        rm = 41.03; % See Equation 7    
        ecc_to_diam_deg = @(rf_r_deg)(sqrt(2./(sqrt(3).*... % Equation 9. Distance between adjacent midget RF centers.
            (2*dc_0.*(1+rf_r_deg./rm).^-1.*(a*(1+(rf_r_deg./r2)).^-2+(1-a)*exp(-rf_r_deg./re)))...
            /2))... % Dividing y by 2 to estimate RF size from only ON or OFF mosaics (halved density). 
            *.80); % Dacey and Petersen: monkey midget DFs are .77:.81 x smaller than human midget DFs
        
    end
    % Eq. 10 is giving me 10,000 RFs/deg^2 at 0.1� eccentricity. Should be
    % 3*higher density (3*smaller RF diameters).
end

if ~exist('ecc_to_diam_deg','var')
    error('Unable to assign ecc_to_diam_deg. Check definitions of CELLTYPE and RFSIZESOURCE')
end
    
% how many parvo RFs fit inside a magno RF?
rf_areas = [];
for rf_r_deg = 2:12
     rf_areas= [rf_areas; pi*(0.0683+0.0165*rf_r_deg).^2 pi*(DEGPERMM*8.64*(rf_r_deg*MMPERDEG).^1.04/1000).^2];
end
% order: [magno parvo]
%rf_areas(:,1)./rf_areas(:,2)
% Looks like between 6 and 30

% EJ says midgets should be at 4x the density of parasols. 

% Moving out 1 mm along the nasal retina is like moving 0.61 mm along the
% temporal retina in terms of RF size (nasal retina is like an extended
% fovea). Croner and Kaplan express eccentricity as temporal equivalent. I
% should probably average nasal and temporal RF sizes 
%%
% Section 1
% d' as a function of RF position, color, and TF.
% Luminance sensitivity is represented by a black disk.
% Red-green sensitivity is represented by a red ring.

[filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',{CELLTYPE},'subjID',{MONKEY(1)});
HIGHLOWTHRESHOLD = 5;
data = [];
for i = 1:length(filenames)
    stro = {};
    for j = 1:length(filenames{i})
        stro{j} = nex2stro(char(findfile(filenames{i}(j), fullfile(nexfilepath,'Greg',MONKEY))));
    end
    stro = strocat(stro);
    [uniquestim, dprime] = IsoSampGetDPrime(stro,DPRIMEMETHOD, spikecds(i));
    Lblank = uniquestim(:,1) == 0 & uniquestim(:,2) == 0 & uniquestim(:,3) == 0;
    Lrg = sign(uniquestim(:,1)) ~= sign(uniquestim(:,2));
    Llum = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & ~Lblank;
    LlowTF = uniquestim(:,3) < HIGHLOWTHRESHOLD;
    LhighTF = uniquestim(:,3) >= HIGHLOWTHRESHOLD;
    data = [data; stro.sum.exptParams.rf_x/10, stro.sum.exptParams.rf_y/10, nanmean(dprime(Lrg&LlowTF)) nanmean(dprime(Lrg&LhighTF))  nanmean(dprime(Llum&LlowTF)) nanmean(dprime(Llum&LhighTF))];
end

% Setting up column labels
XIDX = 1;
YIDX = 2;
RGLOWIDX = 3;
RGHIGHIDX = 4;
LUMLOWIDX = 5;
LUMHIGHIDX = 6;
labels = {'RGlow','RGhigh','LUMlow','LUMhigh'};

% Looking at d-prime for high/low frequency RG/LUM stimuli across the
% visual field
figure; subplot(2,1,1); hold on;
for i = 1:size(data,1)
    h = plot(abs(data(i,XIDX)),data(i,YIDX),'ko','MarkerFaceColor','black');
    set(h,'MarkerSize',max(data(i,LUMLOWIDX)*20,1));
    h = plot(abs(data(i,1)),data(i,2),'ro');
    set(h,'MarkerSize',max(data(i,RGLOWIDX)*20,1));
end
title('Low frequencies');
axis equal;
axis square;

subplot(2,1,2); hold on;
for i = 1:size(data,1)
    h = plot(abs(data(i,XIDX)),data(i,YIDX),'ko','MarkerFaceColor','black');
    set(h,'MarkerSize',max(data(i,LUMHIGHIDX)*20,1));
    h = plot(abs(data(i,1)),data(i,2),'ro');
    set(h,'MarkerSize',max(data(i,RGHIGHIDX)*20,1));
end
title('High frequencies');
axis equal;
axis square;
%prettycorr(data(:,[3:end]),labels);
%prettycorr(data,['rf','rfy',labels]);
%title('Individual neuron d-primes');

%%
% Section 2
% Estimates of single neuron d' as a function of eccentricity

MONKEYS={'Apollo','Utu'};
CELLTYPES = {'M','P'};
DPRIMEMETHOD = 1;

data = [];
allfilenames = {};
for monkey = MONKEYS
    for celltype = CELLTYPES
        [filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',celltype,'subjID',{monkey{1}(1)});     
        for i = 1:length(filenames)
            stro = {};
            for j = 1:length(filenames{i})
                stro{j} = nex2stro(char(findfile(filenames{i}(j), fullfile(nexfilepath,'Greg',monkey))));
            end
            stro = strocat(stro);
            rfx = stro.sum.exptParams.rf_x/10;
            rfy = stro.sum.exptParams.rf_y/10;
            %rf_r_deg = sqrt(.8*rfx^2+rfy^2);
            %RF_diam_deg = ecc_to_diam_deg(rf_r_deg);
            [uniquestim, dprime] = IsoSampGetDPrime(stro,DPRIMEMETHOD,spikecds(i));
            
            Lblank = uniquestim(:,1) == 0 & uniquestim(:,2) == 0 & uniquestim(:,3) == 0;
            Llum = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & ~Lblank;
            
            % Getting RF size at rfx, rfy
            [theta,rf_r_deg] = cart2pol(rfx,rfy);
            % This stuff below *must* agree with definitions at the top of
            % the script
            if strcmp(celltype,'M')
                rf_diam_deg = 10.^(-1.2459+0.0345*rf_r_deg); % D&L 1984 temporal retina equivalent
            else
                a = 0.9729; % Table 1
                r2 = 1.084; % Table 1
                re = 7.633; % Table 1
                dc_0 = 14804.6; % Cone density of fovea
                rm = 41.03; % See Equation 7
                rf_diam_deg = sqrt(2./(sqrt(3).*... % Equation 9. Distance between adjacent midget RF centers.
                    (2*dc_0.*(1+rf_r_deg./rm).^-1.*(a*(1+(rf_r_deg./r2)).^-2+(1-a)*exp(-rf_r_deg./re)))...
                    ./2))... % Dividing y by 2 to estimate RF size from only ON or OFF mosaics (halving the density).
                    *.58; % Dividing the distance between adjacent midget RF centers by .58 to get RF radius. Makes midget RFs 1/8 area of parasol.
            end
            
            data = [data; find(strcmp(monkey,MONKEYS)) find(strcmp(celltype,CELLTYPES)) mean(dprime(Llum)) rfx rfy rf_diam_deg];
            allfilenames{length(allfilenames)+1} = filenames{i}(1);
        end
    end
end

colors = [1 1 1; 0 0 0];
hax = [];
figure('Position',[440 100 645 700]);
hax(1) = axes('Units','inches','Position',[.75 .5 2 4],'Tickdir','out');
hax(2) = axes('Units','inches','Position',[5 .5 2 4],'Tickdir','out');
hax(3) = axes('Units','inches','Position',[.75 6 3 3],'Tickdir','out'); 
hax(4) = axes('Units','inches','Position',[5 6 3 3],'Tickdir','out'); 

for monkey_idx = unique(data(:,1))'
    for celltype_idx = unique(data(:,2))'
        L = data(:,1) == monkey_idx & data(:,2) == celltype_idx;
        axes(hax(celltype_idx)); hold on;
        for i = find(L)'
            h = plot(data(i,4),data(i,5),'ko','MarkerFaceColor','none');
            set(h,'MarkerSize',max(.1,data(i,3)./max(data(:,3))*40));
            set(h,'ButtonDownfcn',['disp(''',char(allfilenames{i}),''')']);
            if monkey_idx ==2
                set(h,'Marker','^');
            end
        end
    end
end

for i = 1:2
    axes(hax(i));
    set(gca,'Ylim',[-14 14],'Xlim',[-14 0]);
    plot(0,0,'r+','MarkerSize',10);
end

% Scatterplot of d-prime vs eccentricity
axes(hax(3)); hold on;
set(gca,'Xlim',[2 14]);
xlabel('eccentricity (�)');
ylabel('SNR');
for celltype_idx = unique(data(:,2))'
    for monkey_idx = unique(data(:,1))'
        L = data(:,1) == monkey_idx & data(:,2) == celltype_idx;
        h = plot(sqrt(data(L,4).^2+data(L,5).^2),data(L,3),'ko','MarkerFaceColor',colors(celltype_idx,:));
        if monkey_idx ==2
            set(h,'Marker','^');
        end
        [b,~,~,~,stats] = regress(data(L,3),[ones(sum(L),1) sqrt(data(L,4).^2+data(L,5).^2)]);
        disp(['SNR: ',MONKEYS{monkey_idx},': ',CELLTYPES{celltype_idx},': p = ',num2str(stats(3)),': b = ',num2str(b(2))]);
        %if stats(3) < 0.05
        %    h = plot([2 14],[[1 2]*b [1 14]*b],'k-');
        %    if monkey_idx == 2
        %        set(h,'LineStyle','--');
        %    end
        %end
    end
    L = data(:,2) == celltype_idx;
    [b,~,~,~,stats] = regress(data(L,3),[ones(sum(L),1) sqrt(data(L,4).^2+data(L,5).^2)]);
    disp(['SNR: ',CELLTYPES{celltype_idx},': p = ',num2str(stats(3)),': b = ',num2str(b(2))]);
    plot([2 14],[[1 2]*b [1 14]*b],'k-');
end

% Scatterplot of d-prime/RF size vs eccentricity
axes(hax(4)); hold on;
set(gca,'Xlim',[2 14]);
xlabel('eccentricity (�)');
ylabel('SNR/RFdiam');
for celltype_idx = unique(data(:,2))'
    for monkey_idx = unique(data(:,1))'
        L = data(:,1) == monkey_idx & data(:,2) == celltype_idx;
        h = plot(sqrt(data(L,4).^2+data(L,5).^2),data(L,3)./data(L,6),'ko','MarkerFaceColor',colors(celltype_idx,:));
        if monkey_idx ==2
            set(h,'Marker','^');
        end
        [b,~,~,~,stats] = regress(data(L,3)./data(L,6),[ones(sum(L),1) sqrt(data(L,4).^2+data(L,5).^2)]);
        disp(['SNR/RFdiam: ',MONKEYS{monkey_idx},': ',CELLTYPES{celltype_idx},': p = ',num2str(stats(3)),': b = ',num2str(b(2))]);
        %if stats(3) < 0.05
        %    h = plot([2 14],[[1 2]*b [1 14]*b],'k-');
        %    if monkey_idx == 2
        %        set(h,'LineStyle','--');
        %    end
        %end
    end
    L = data(:,2) == celltype_idx;
    [b,~,~,~,stats] = regress(data(L,3)./data(L,6),[ones(sum(L),1) sqrt(data(L,4).^2+data(L,5).^2)]);
    disp(['SNR/RFdiam: ',CELLTYPES{celltype_idx},': p = ',num2str(stats(3)),': b = ',num2str(b(2))]);
    plot([2 14],[[1 2]*b [1 14]*b],'k-');
end

%%
% Section 2.1
% Asking whether neuronal d-prime depends on eccentricity and/or cone
% weights. Averaging d-primes across TFs?

% Below taken from DTcones_gh.m (lines 834:840)
temporalcoeffs=[150.9676 -1.2220 35.9979 -0.1567 9.9936 -0.0258]; % (Goodchild et al., 1996)
nasalcoeffs=[176.6624 -7.9473 94.4908 -0.3518 18.6761 -0.0236]; % (Packer et al. 1989) See GrantBrainStorming.m section 30
conedensfun = @(coeffs,x)(coeffs(1).*(exp(coeffs(2).*x)))+...
    (coeffs(3).*(exp(coeffs(4).*x)))+...
    (coeffs(5).*(exp(coeffs(6).*x)));

[filenames,spikecds,~,neuronid] = fnamesFromTxt('IsoSamp_LGN','cellClass',{CELLTYPE},'subjID',{MONKEY(1)});

data = [];
for i = 1:length(filenames)
    stro = {};
    for j = 1:length(filenames{i})
        stro{j} = nex2stro(char(findfile(filenames{i}(j), fullfile(nexfilepath,'Greg',MONKEY))));
    end
    stro = strocat(stro);
    rfx = stro.sum.exptParams.rf_x/10;
    rfy = stro.sum.exptParams.rf_y/10;
    %rf_r_deg = sqrt(.8*rfx^2+rfy^2);
    %RF_diam_deg = ecc_to_diam_deg(rf_r_deg);
    [uniquestim, dprime] = IsoSampGetDPrime(stro,DPRIMEMETHOD,spikecds(i));
    
    % Now getting white noise
    [wnfilenames] = fnamesFromTxt('WhiteNoiseLGN_forIS','cellClass',{CELLTYPE},'subjID',{MONKEY(1)},'neuron',{num2str(neuronid(i))});
    stro = {};
    for j = 1:length(wnfilenames)
        stro{j} = nex2stro(char(findfile(wnfilenames{1}{j}, fullfile(nexfilepath,'Greg',MONKEY))));
    end
    stro = strocat(stro);
    spikename = ['sig001',char(double('a'+spikecds(i)-1))];
    nstixperside = stro.sum.exptParams.nstixperside;

    % Reconstructing the M matrix and gamma table
    fundamentals = stro.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = stro.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
 
    maxT = 6;
    out = getWhtnsStats(stro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, spikename);
    STA = out{1}/out{3};
    STA = reshape(STA,[nstixperside.^2  3 maxT]);
 
    noise = sum(STA(:,:,1).^2,2)*maxT;
    energy = sum(sum(STA.^2,2),3);
    whichpix_tmp = energy > max(noise) | energy == max(energy);
    
    if sum(whichpix_tmp) > 1
        [tmpi,tmpj] = ind2sub([nstixperside nstixperside],find(whichpix_tmp));
        ij = [tmpi,tmpj];
        T = clusterdata(ij,'linkage','single','distance','euclidean','criterion','distance','cutoff',sqrt(2));
        
        clusternmembers = [];
        for k =1:max(T)
            clusternmembers(k) = sum(T == k);
        end
        dominantcluster = find(clusternmembers == max(clusternmembers));
        
        clustermat = zeros(nstixperside, nstixperside);
        clustermat(sub2ind(size(clustermat),ij(T==dominantcluster,1),ij(T==dominantcluster,2))) = 1;
        
        if sum(clustermat(:)) > 2 & length(unique(ij(T==dominantcluster,1))) > 1 & length(unique(ij(T==dominantcluster,2))) > 1
            % Then get convex hull
            dominantclusteridxs = ij(T==dominantcluster,:);
            K=convhull(dominantclusteridxs(:,1), dominantclusteridxs(:,2));
            tmp = dominantclusteridxs(K,:);
            [x,y] = meshgrid(1:nstixperside,1:nstixperside);
            inRF = reshape(inpolygon(x(:),y(:),dominantclusteridxs(K,2),dominantclusteridxs(K,1)),[nstixperside, nstixperside]);
        else
            inRF = clustermat;
        end
    else
        inRF = reshape(whichpix_tmp,nstixperside,nstixperside);
    end
    reshape(inRF,nstixperside,nstixperside)
    temporalSTA = squeeze(STA(logical(inRF(:)),:,:));
    if ~ismatrix(temporalSTA)
        temporalSTA = mean(temporalSTA,1);
    end
    [u,~,~] = svd(squeeze(temporalSTA));
    cone_weights = inv(M')*u(:,1);
    prefcolordir = atan2(cone_weights(2),cone_weights(1)); % in LM plane
    
    Lblank = uniquestim(:,1) == 0 & uniquestim(:,2) == 0 & uniquestim(:,3) == 0;
    Llum = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & ~Lblank;
    Lrg = sign(uniquestim(:,1)) ~= sign(uniquestim(:,2));

    % Getting cone density at rfx, rfy
    [theta,r] = cart2pol(rfx,rfy);
    coneDensity = (cos(theta).^2.*(conedensfun(temporalcoeffs,r)+conedensfun(nasalcoeffs,r))/2) + (sin(theta).^2.*conedensfun(temporalcoeffs,r));
    
    data = [data; mean(dprime(Lrg)) rfx rfy coneDensity prefcolordir];
end

figure; axes; hold on;

pref_dir = mod(data(:,5),pi);
plot(pref_dir,data(:,1),'ro');
set(gca,'Xlim',[0 pi]);
plot([pi/4 pi/4],[0 .5])
plot([3*pi/4 3*pi/4],[0 .5])
ylabel('d-prime');
xlabel('pref dir in LM plane');

figure;
subplot(2,2,1); hold on;
for i = 1:size(data,1)
   h = plot(data(i,2),data(i,3),'ko','MarkerFaceColor',[.5 .5 .5]);
   set(h,'MarkerSize',max(.1,data(i,1)./max(data(:,1))*20))
end
plot(0,0,'+','MarkerSize',12);
axis equal;
title([MONKEY,' ',CELLTYPE,' cells (n = ',num2str(length(filenames)),')'])

subplot(2,2,3);
plot(sqrt(data(:,2).^2+data(:,3).^2),data(:,1),'ko','MarkerFaceColor','black');
lsline;
xlabel('eccentricity (�)');
ylabel('SNR');
title('SNR');
[b,bint,~,~,stats] = regress(data(:,1),[ones(size(data,1),1) sqrt(data(:,2).^2+data(:,3).^2)]);
if stats(3) < 0.05
    set(gca,'Xcolor','yellow','Ycolor','yellow');
end

subplot(2,2,2); hold on;
dp_cd = data(:,1).*data(:,4); % d-prime x cone density
for i = 1:size(data,1)
   h = plot(data(i,2),data(i,3),'ko','MarkerFaceColor',[.5 .5 .5]);
   set(h,'MarkerSize',max(.1,dp_cd(i)./max(dp_cd)*20))
end
plot(0,0,'+','MarkerSize',12);
title('SNR * cone density');

axis equal;
subplot(2,2,4);
plot(sqrt(data(:,2).^2+data(:,3).^2),dp_cd,'ko','MarkerFaceColor','black');
lsline;
ylims = get(gca,'Ylim');
set(gca,'Ylim',[0 ceil(ylims(2)/10)*10]);
xlabel('eccentricity (�)');
ylabel('SNR * cone density');
[b,bint,~,~,stats] = regress(dp_cd,[ones(size(data,1),1) sqrt(data(:,2).^2+data(:,3).^2)]);
if stats(3) < 0.05
    set(gca,'Xcolor','yellow','Ycolor','yellow');
end


%%
% Section 3
% Efficiency of input (efficiency = neuron d'/cone ensemble d')
% Estimates of single neuron d' based on estimates of RF size and the cone noise model.
% Only considering cones inside the (average) RF.D
% Interesting if single neuron d' varies with eccentricity but efficiency
% does not.
% 1) Upper panel of per-cell plots show the (approximate) d' of the population of cones 
% inside the RF to the stimuli shown. Lower panel shows "input efficiency":
% d' of neuron/d' of cones in RF. An inefficiency of '1' means no SNR drop
% from cones to neuron (or error in RF size estimate or something else).
% 2) Scatterplot of input efficiency for different stimulus types.
% 3) Boxplot shows distributions of input efficiencies (averaged across
% conditions). Each neuron is counted once per box plot.
% 4&5) Lack of correlation between eccentricity and input efficiency suggests
% that, when you adjust for number of cones in the RF, the relationship
% between neural sensitivity and eccentricity (which is confounded with
% contrast) goes away.

TEMPORONASALSCALEFACTOR = 0.8;

[filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',{CELLTYPE},'subjID',{MONKEY(1)});
fundamentals = load('T_cones_smj10');
params = [];
params.runType = 'isosamp';
params.obsMethod = 'obsMethod_filteredWtFxn';
params.impulseResponse = 'rieke';
params.DTV1_fname = [];
params.DTNT_fname = [];
params.unitTest = false;
params.eqMosaic = false;
params.saveDir = '/Users/greghorwitz/Documents/MATLAB';
params.notes = 'IsoSamp test';
params.parallelOperations = false;
params.eyeType = 'monkey';
params.eyeNumber = 1; % LGN neurons are monocular
params.coneSampRate = 2400;
params.flatPowerSpect = false;
params.enableScones = false;
params.sacamp_deg = 0;
params.sacdur_s = 0;
data = [];
DEBUGGING = [];
for a = 1:length(filenames)
    stro = {};
    for j = 1:length(filenames{a})
        stro{j} = nex2stro(char(findfile(filenames{a}(j), fullfile(nexfilepath,'Greg',MONKEY))));
    end
    stro = strocat(stro);
    sigma_gabor = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sigma'))); % DVA

    [uniquestim, dprime] = IsoSampGetDPrime(stro,DPRIMEMETHOD,spikecds(a));
    
    params.stro = stro;
    cal.gamma = params.stro.sum.exptParams.gamma_table; % Only used for "dtv1". Can ignore.
    spds = params.stro.sum.exptParams.mon_spd;
    spds = SplineSpd([380:4:780]',reshape(spds,length(spds)/3,3),[380:5:780]');
    cal.monSpect = spds(:);
    M = fundamentals.T_cones_smj10*spds;
    cal.Mmtx = M(:);
    cal.frameRate = params.stro.sum.exptParams.framerate;
    cal.bkgndrgb = params.stro.sum.exptParams.bkgndrgb';
    cal.fname = 'test';
    cal.monSpectWavelengths = linspace(380,780,101);
    cal.pixperdeg = params.stro.sum.exptParams.pixperdeg;
    params.monCalFile = cal;
    
    rfx = stro.sum.exptParams.rf_x/10;
    rfy = stro.sum.exptParams.rf_y/10;
    rf_r_deg = sqrt((rfx./TEMPORONASALSCALEFACTOR)^2+rfy^2);
    rf_diam_deg = ecc_to_diam_deg(rf_r_deg); % diameter RF
    RF_STD = rf_diam_deg/2; % 1 standard deviation of Gaussian RF
    params.gab.sd = sqrt(1/((1/RF_STD.^2)+(1/sigma_gabor.^2))); % Product of RF Gaussian and Gabor Gaussian
    DEBUGGING =[DEBUGGING; RF_STD sigma_gabor params.gab.sd]
    [gab, cones, mon, idlob, params] = DTcones_gh(params,0);
    
    conemodeldata = [];
    for i = 1:size(idlob.analyticMean,1) % looping over color direction
        for j = 1:size(idlob.analyticMean(i,:),2) % looping over contrast/TF
            if ~isempty(idlob.analyticMean{i,j})
                tmp_lm_mu = idlob.analyticMean{i,j}([1 2]);
                tmp_lm_var = idlob.analyticVar{i,j}([1 2]);
                
                tf = gab.driftRates{i}(j);
                conemodeldata = [conemodeldata; gab.colorDirs(i,[1 2]) tf tmp_lm_mu tmp_lm_var];
            end
        end
    end
    % columns of data: L, M, TF, mu_L, mu_M, sigma2_L, sigma2_M
    cone_dprimes = [];
    for i = 1:size(uniquestim,1)
        L = sign(conemodeldata(:,1)) == sign(uniquestim(i,1)) &...
            sign(conemodeldata(:,2)) == sign(uniquestim(i,2)) &...
            conemodeldata(:,3) == uniquestim(i,3);
        if sum(L) > 1
            disp('too many condition matches');
            continue
        elseif sum(L) == 0
            cone_dprime = nan;
        else
            cone_dprime = sqrt(conemodeldata(L,[4 5]).^2*(1./conemodeldata(L,[6 7])')); % Not assuming equal L:M
        end
        cone_dprimes = [cone_dprimes; cone_dprime];
    end
    photon_dprimes = IsoSampGetPhotonDPrime (gab.flashTimeProf, mon.frameRate, mon.bkgndlms_Rstar, cat(3,cones.num_L,cones.num_M), 4, uniquestim);

    data{a}.rfxy = [stro.sum.exptParams.rf_x/10, stro.sum.exptParams.rf_y/10];
    data{a}.uniquestim = uniquestim;
    data{a}.neurondprime = dprime;
    data{a}.conedprime = cone_dprimes;
    data{a}.photondprime = photon_dprimes;
end

% Raw data dump
figure; axes; hold on;
clear h;
for a = 1:length(data)
    for whichstim = 1
        if whichstim == 1
            L = sign(data{a}.uniquestim(:,1)) == sign(data{a}.uniquestim(:,2)) & data{a}.uniquestim(:,1) ~= 0;
        else
            L = sign(data{a}.uniquestim(:,1)) ~= sign(data{a}.uniquestim(:,2));
        end
        h(1) = plot(data{a}.uniquestim(L,3),data{a}.neurondprime(L),'k-');
        h(2) = plot(data{a}.uniquestim(L,3),data{a}.conedprime(L),'k:');
        h(3) = plot(data{a}.uniquestim(L,3),data{a}.photondprime(L),'k--');

        if whichstim == 2
            set(h,'Color','red');
        end
    end
end
set(gca,'Xscale','log');

% Plotting dprime as a function of TF


% Preparing for plotting. Computing efficiencies.

tfbinedges = logspace(0,log10(60),15);
tfbincenters = sqrt(tfbinedges(1:end-1).*tfbinedges(2:end)); % geomean
efficiencies = cell(2,length(tfbincenters));
ns = zeros(2,length(tfbincenters));
for i = 1:length(data)
    in_ecc_range = sqrt(data{i}.rfxy*data{i}.rfxy') >= 2 & abs(data{i}.rfxy(1)) <= 14 & abs(data{i}.rfxy(2))<= 8;
    if ~in_ecc_range
        disp('chucking a cell')
        keyboard
    end
    for STIMTYPE = {'LUM','RG'}
        stimtypeidx = find(strcmp(STIMTYPE, {'LUM','RG'}));
        if strcmp(STIMTYPE,'LUM')
            Lcolordir = sign(data{i}.uniquestim(:,1)) == sign(data{i}.uniquestim(:,2)) & data{i}.uniquestim(:,3) ~= 0;
        else strcmp(STIMTYPE,'RG')
            Lcolordir = sign(data{i}.uniquestim(:,1)) ~= sign(data{i}.uniquestim(:,2));
        end
        
        for j = 1:length(tfbincenters)
            Ltf = data{i}.uniquestim(:,3) > tfbinedges(j) & data{i}.uniquestim(:,3) <= tfbinedges(j+1);
            if sum(Ltf&Lcolordir) > 0
                efficiencies{stimtypeidx,j} = [efficiencies{stimtypeidx,j}; mean(data{i}.neurondprime(Ltf&Lcolordir))/mean(data{i}.conedprime(Ltf&Lcolordir))];               
                ns(stimtypeidx,j) = ns(stimtypeidx,j) + 1; % each cell counts as an independent entity
            end
        end
    end
end

m = zeros(size(efficiencies));
sem = zeros(size(efficiencies));
for i = 1:numel(efficiencies)
    m(i) = mean(efficiencies{i});
    sem(i) = sqrt(var(efficiencies{i})./ns(i));
end

colors = {'black','red'};
figure; axes; hold on;
for i = 1:2
    L = ns(i,:) > 2; % minimum 'n'
    h = plot(tfbincenters(L), m(i,L),'k-o','LineWidth',2);
    set(h,'Color',colors{i},'MarkerFaceColor',colors{i});
    h = patch([tfbincenters(L), fliplr(tfbincenters(L))],[m(i,L)+sem(i,L), fliplr(m(i,L)-sem(i,L))],colors{i},'Facealpha',.5);
end
set(gca,'Xscale','log');
ylabel('Input efficiency');
xlabel('Temporal frequency (Hz)');
plot([tfbincenters(1) tfbincenters(end)],[0 0],'k-');
plot([tfbincenters(1) tfbincenters(end)],[1 1],'k-');

% Can I see the expected decrease in efficiency (for color only) with
% eccentricity?
tmpdata = [];
for i = 1:length(data)
    in_ecc_range = sqrt(data{i}.rfxy*data{i}.rfxy') >= 2 & abs(data{i}.rfxy(1)) <= 12 & abs(data{i}.rfxy(2))<= 8;
    if ~in_ecc_range
        continue
    end
    Llum = sign(data{i}.uniquestim(:,1)) == sign(data{i}.uniquestim(:,2)) & data{i}.uniquestim(:,3) ~= 0;
    Lrg = sign(data{i}.uniquestim(:,1)) ~= sign(data{i}.uniquestim(:,2));
    tmpdata = [tmpdata; data{i}.rfxy nanmean(data{i}.neurondprime(Llum)./data{i}.conedprime(Llum))  nanmean(data{i}.neurondprime(Lrg)./data{i}.conedprime(Lrg))];
end
ecc = sqrt(sum(tmpdata(:,1).^2+tmpdata(:,1).^2,2));
prettycorr([ecc,tmpdata(:,[2 3])],{'ecc','lum efficiency','rg efficiency'}); % Nothing yet given how noisy the data are

%%
% Section 4
%
% D-prime of a population of "like-cells" with RFs that completely cover
% the stimulus. Nothing about cones in this analysis. Each row of "data"
% gives "population-scaled" d' values for a single neuron. 
%
% Population scaling is achieved by making a hexagonal grid of RFs on the
% stimulus (�4 std) spaced so that the distance between neighboring RFs 
% is 2x the standard deviation of the a Gaussian approximation to the RF
% center (so neighboring RFs touch at the 1 std boundary).
% The population d' is a scaled version of the single cell d'.
% See StatsStuff.m sections 12&13 for calculations of optimal weight and SNR
% of each LGN neuron.
%
% Here's where "population_scalefactor" comes from:
% The ideal linear decoder uses weights that given by w=inv(S2)*mu
% Where mu is the signal they carry and S2 is the covariance matrix, which
% we are assuming is given y pairwise RF overlap.
% Let Y be the population decision variable.
% Y = w1*X1+w2*X2+...+wn*Xn
% Where X1 is the response of the recorded neuron and X2:Xn are the
% responses of other neurons that are assumed to have the same variance
% as X1 but different means because their RFs aren't perfectly centered
% on the stimulus. wis are optimal linear weights.
% d'(one neuron) = E(X1)/sqrt(Var(X1))
% Assume that Var(Xi) = 1. This gives E(X1) = d'(one neuron).

% E(Y) = w1*E(X1)+w2*E(X2)+...+wn*E(Xn)
% E(Xi) = ai*E(X1). <-- each neuron's mean (signal) is a scaled version of E(X1)
% (This is essentially a linearity assumption)
% E(Y) = w1*a1*E(X1)+w2*a2*E(X1)+...+wn*an*E(X1) <-- ais are scale factors capturing
% differential drive from from stimulus; wis are optimal weights (from the covariance-savvy
% ideal observer)
% E(Y) = w'*a*(single neuron d')
% Var(Y) = Var(sum(wi*Xi))
% Var(Y) = w1^2*Var(X1)+w2^2*Var(X2)+...+wn^2*Var(Xn)+2*w1*w2*Cov(X1,X2)+...
% Var(Y) = w*w'*S2; = w*S2*w'
% d1'(Y) = E(Y)/sqrt(Var(Y)) = d1'(one neuron)*(w'*a)/sqrt(w'*w*S2)
%        = d1'(one neuron)*sqrt(sum(ai.^2))

% RGlow sensitivity is correlated with vertical position (high d' at lowest 
% positions). I think I'm  delivering less contrast to the LVF than the upper
% so I don't think I'm making this happen by underestimating psychophysical
% sensitivity in the LVF. What is going on here? Is this just because 
% contrast is higher farther from the fovea and we have few points in the
% UVF? 
NUMBEROFEYES = 2;
TEMPORONASALSCALEFACTOR = .8; % Changing X eccentricity to account for the fact that 50% of the time the cell had a nasal RF (mean(0.61,1) = 0.8)
% ecc_to_diam_deg should return 1 SD of RF size in equivalent temporal eccentricity.
ONOFFCORRELATION = 0.05;
RFTRUNCATIONINSD = 3; % How far from the mean to go out *in every direction* for estimates of RF overlap and integrated contrast

[filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',{CELLTYPE},'subjID',{MONKEY(1)});
bpdf_vec=@(x,y,mu_x,mu_y,sigma)(exp(-((x-mu_x).^2)./(2*sigma^2)-((y-mu_y).^2)/(2*sigma^2))./(2*pi*sigma^2)); % bivariate normpdf

data = [];
for i = 1:length(filenames)
    stro = {};
    for j = 1:length(filenames{i})
        stro{j} = nex2stro(char(findfile(filenames{i}(j), fullfile(nexfilepath,'Greg',MONKEY))));
    end
    stro = strocat(stro);
    population_scalefactor = IsoSampGetPopulationScaleFactor(stro, ecc_to_diam_deg,TEMPORONASALSCALEFACTOR,RFTRUNCATIONINSD, ONOFFCORRELATION,NUMBEROFEYES);
    [uniquestim, dprime] = IsoSampGetDPrime(stro,DPRIMEMETHOD, spikecds(i));
    data{i}.filenames = filenames{i};
    data{i}.spikecds = spikecds(i);
    data{i}.rfxy = [stro.sum.exptParams.rf_x/10, stro.sum.exptParams.rf_y/10];
    data{i}.uniquestim = uniquestim;
    data{i}.dprime = dprime;
    data{i}.population_scalefactor = population_scalefactor;
end

emptycells = [];
for i = 1:length(data)
    if isempty(data{i})
        keyboard
        data{i}.rf = [nan nan];
        data{i}.uniquestim = nan;
        emptycells = [emptycells i];
    end
end
%data(emptycells) = [];
% Plotting mean�SEM. It would be better to use median � IQR?
tfbinedges = logspace(0,log10(60),10);
tfbincenters = sqrt(tfbinedges(1:end-1).*tfbinedges(2:end)); % geomean
summeddprimes = zeros(2,length(tfbincenters));
summeddprimes2 = zeros(2,length(tfbincenters));
alldprimes = cell(size(summeddprimes));
ns = zeros(2,length(tfbincenters));
for i = 1:length(data)
    Lrg = sign(data{i}.uniquestim(:,1)) ~= sign(data{i}.uniquestim(:,2));
    Llum = sign(data{i}.uniquestim(:,1)) == sign(data{i}.uniquestim(:,2)) & data{i}.uniquestim(:,3) ~= 0;
    in_ecc_range = sqrt(data{i}.rfxy*data{i}.rfxy') >= 2 & abs(data{i}.rfxy(1)) <= 12 & abs(data{i}.rfxy(2))<= 8;
    if ~in_ecc_range
        disp('Chucking a cell because RF is outside of the acceptable range');
        data{i}.rfxy
        continue
    end
    for j = 1:length(tfbincenters)
        Ltf = data{i}.uniquestim(:,3) > tfbinedges(j) & data{i}.uniquestim(:,3) <= tfbinedges(j+1);
        if sum(Ltf&Llum) > 0
            summeddprimes(1,j) = summeddprimes(1,j)+mean(data{i}.dprime(Ltf&Llum).*data{i}.population_scalefactor);
            summeddprimes2(1,j) = summeddprimes2(1,j)+mean(data{i}.dprime(Ltf&Llum).*data{i}.population_scalefactor).^2;
            alldprimes{1,j} = [alldprimes{1,j};mean(data{i}.dprime(Ltf&Llum).*data{i}.population_scalefactor)];
            ns(1,j) = ns(1,j) + 1; % each cell counts as an independent entity
        end
        if sum(Ltf&Lrg) > 0
            summeddprimes(2,j) = summeddprimes(2,j)+mean(data{i}.dprime(Ltf&Lrg).*data{i}.population_scalefactor);
            summeddprimes2(2,j) = summeddprimes2(2,j)+mean(data{i}.dprime(Ltf&Lrg).*data{i}.population_scalefactor).^2;
            alldprimes{2,j} = [alldprimes{2,j};mean(data{i}.dprime(Ltf&Lrg).*data{i}.population_scalefactor)];
            ns(2,j) = ns(2,j) + 1; % each cell counts as an independent entity
        end
        if any(isnan(summeddprimes(:)))
            keyboard
        end
    end
end
m = summeddprimes./ns;
v = (summeddprimes2-(summeddprimes.^2)./ns)./ns;
sem = sqrt(v./ns);

figure; axes; hold on;
L = ns(1,:) > 1;
h = patch([tfbincenters(L), fliplr(tfbincenters(L))],[m(1,L)+sem(1,L), fliplr(m(1,L)-sem(1,L))],'black') ;
set(h,'Facealpha',.5);
plot(tfbincenters(L),m(1,L),'k-','Linewidth',2);

L = ns(2,:) > 1;
h= patch([tfbincenters(L), fliplr(tfbincenters(L))],[m(2,L)+sem(2,L), fliplr(m(2,L)-sem(2,L))],'red'); 
set(h,'Facealpha',.5);
plot(tfbincenters(L),m(2,L),'r-','Linewidth',2);
set(gca,'Xscale','log');
ylabel('d''');
xlabel('Temporal frequency (Hz)');
title([MONKEY,' ',CELLTYPE,' cells (n = ',num2str(ns(1)),')'])
plot([tfbinedges(1) tfbinedges(end)],[1.27 1.27],'k-')

% One tailed t-tests on a bin by bin basis
for mu = [0 1.27]
    t = (m-mu)./sem;
    p = 1-tcdf(t,ns)
end

% RF eccentricity and TF
% Pre-population scalefactor multiplication
tmp = [];
for i = 1:length(data)
    Llum = sign(data{i}.uniquestim(:,1)) == sign(data{i}.uniquestim(:,2)) & data{i}.uniquestim(:,3) ~= 0;
    LTF = data{i}.uniquestim(:,3) <= 40; % Some TFs are only tested near the fovea
    dprimes = data{i}.dprime(Llum&LTF); % Single neuron d'
    TFs = data{i}.uniquestim(Llum&LTF,3);
    b = regress(dprimes,[ones(length(dprimes),1),log10(TFs)]);
    dprime = mean(dprimes);
    ecc = sqrt(sum(data{i}.rfxy.^2));
    tmp = [tmp; ecc dprime b(2)];
end
figprefs; axes;
plot(tmp(:,1),tmp(:,3),'ko');
xlabel('eccentricity'); ylabel('slope of d-prime vs TF');
title([MONKEY,' ',CELLTYPE,' cells (n = ',num2str(ns(1)),')'])
set(gca,'TickDir','out')

% Does relative sensitivity to TF change systematically across the visual field?

%%
% Section 4.1
% Getting the d' of an ideal obsever of all of the cones affected by the
% stimulus for comparison with a modeled LGN population (from above section).
% Run section 4 first (we're using "data" and "NUMBEROFEYES" from there.)
%
% Plotting on the same axes created in Section 4.
%
fundamentals = load('T_cones_smj10');
clear params
params.runType = 'isosamp';
params.obsMethod = 'obsMethod_filteredWtFxn';
params.impulseResponse = 'rieke';
params.DTV1_fname = [];
params.DTNT_fname = [];
params.unitTest = false;
params.eqMosaic = false;
params.saveDir = '/Users/greghorwitz/Documents/MATLAB';
params.notes = 'IsoSamp test';
params.parallelOperations = false;
params.eyeType = 'monkey';
params.eyeNumber = NUMBEROFEYES; % Ideal observer of all cones in both eyes
params.coneSampRate = 2400;
params.flatPowerSpect = false;
params.enableScones = false;
params.sacamp_deg = 0;
params.sacdur_s = 0;

for a = 1:length(filenames)
    if isnan(data{a}.uniquestim)
        continue
    end
    stro = {};
    for j = 1:length(filenames{a})
        stro{j} = nex2stro(char(findfile(filenames{a}(j), fullfile(nexfilepath,'Greg',MONKEY))));
    end
    stro = strocat(stro);
    
    params.stro = stro;
    cal.gamma = params.stro.sum.exptParams.gamma_table; % Only used for "dtv1". Can ignore.
    spds = params.stro.sum.exptParams.mon_spd;
    spds = SplineSpd([380:4:780]',reshape(spds,length(spds)/3,3),[380:5:780]');
    cal.monSpect = spds(:);
    M = fundamentals.T_cones_smj10*spds;
    cal.Mmtx = M(:);
    cal.frameRate = params.stro.sum.exptParams.framerate;
    cal.bkgndrgb = params.stro.sum.exptParams.bkgndrgb';
    cal.fname = 'test';
    cal.monSpectWavelengths = linspace(380,780,101);
    cal.pixperdeg = params.stro.sum.exptParams.pixperdeg;
    params.monCalFile = cal;
    
    [gab, cones, mon, idlob, params] = DTcones_gh(params,0);
    
    tmpdata = [];
    for i = 1:size(idlob.analyticMean,1) % looping over color direction
        for j = 1:size(idlob.analyticMean(i,:),2) % looping over contrast/TF
            if ~isempty(idlob.analyticMean{i,j})
                tmp_lm_mu = idlob.analyticMean{i,j}([1 2]);
                tmp_lm_var = idlob.analyticVar{i,j}([1 2]);                
                tf = gab.driftRates{i}(j);
                tmpdata = [tmpdata; gab.colorDirs(i,[1 2]).*gab.contrasts{i}(j) tf tmp_lm_mu tmp_lm_var];
            end
        end
    end
    
    % Using uniquestim to order the rows of tmpdata and calculating cone
    % dprimes
    data{a}.conedprime = [];  
    for j = 1:size(data{a}.uniquestim,1)
        L = all(abs(data{a}.uniquestim(j,:)-tmpdata(:,[1 2 3]))<1e-10,2);
        if sum(L) ~= 1
            if all(data{a}.uniquestim(j,:) == 0)
                data{a}.conedprime(j) = nan;
            else
                error('sorting problem');
            end
        else
            v = tmpdata(L,[6 7]); % variance
            m = tmpdata(L,[4 5]); % mean
            data{a}.conedprime(j) = sqrt(m.^2*(1./v'));
        end
    end
end

% Ideal observer of the cones as a function of TF
% remember to get rid of points that are too close to the fovea!
tfbinedges = logspace(0,log10(60),10);
tfbincenters = sqrt(tfbinedges(1:end-1).*tfbinedges(2:end)); % geomean
summeddprimes = zeros(2,length(tfbincenters));
summeddprimes2 = zeros(2,length(tfbincenters));
ns = zeros(2,length(tfbincenters));
for i = 1:length(data)
    if isnan(data{i}.uniquestim)
        continue
    end
    in_ecc_range = sqrt(data{i}.rfxy*data{i}.rfxy') >= 2 & abs(data{i}.rfxy(1)) <= 14 & abs(data{i}.rfxy(2))<= 8;
    if ~in_ecc_range
        fprintf('Chucking a cell because RF (%.1d, %.2d) is outside of the acceptable range\n',data{i}.rfxy);
        continue
    end
    Lrg = sign(data{i}.uniquestim(:,1)) ~= sign(data{i}.uniquestim(:,2));
    Llum = sign(data{i}.uniquestim(:,1)) == sign(data{i}.uniquestim(:,2)) & data{i}.uniquestim(:,3) ~= 0;
    for j = 1:length(tfbincenters)
        Ltf = data{i}.uniquestim(:,3) > tfbinedges(j) & data{i}.uniquestim(:,3) <= tfbinedges(j+1);
        if sum(Ltf&Llum) > 0
            summeddprimes(1,j) = summeddprimes(1,j)+mean(data{i}.conedprime(Ltf&Llum));
            summeddprimes2(1,j) = summeddprimes2(1,j)+mean(data{i}.conedprime(Ltf&Llum)).^2;
            ns(1,j) = ns(1,j) + 1; % each cell counts as an independent entity
        end
        if sum(Ltf&Lrg) > 0
            summeddprimes(2,j) = summeddprimes(2,j)+mean(data{i}.conedprime(Ltf&Lrg));
            summeddprimes2(2,j) = summeddprimes2(2,j)+mean(data{i}.conedprime(Ltf&Lrg)).^2;
            ns(2,j) = ns(2,j) + 1; % each cell counts as an independent entity
        end
    end
end
m = summeddprimes./ns;
v = (summeddprimes2-(summeddprimes.^2)./ns)./ns;
sem = sqrt(v./ns);

%figure; axes; hold on;
L = ns(1,:) > 1;
h = patch([tfbincenters(L), fliplr(tfbincenters(L))],[m(1,L)+sem(1,L), fliplr(m(1,L)-sem(1,L))],'black') ;
set(h,'Facealpha',.5);
plot(tfbincenters(L),m(1,L),'k-','Linewidth',2);

L = ns(2,:) > 1;
h= patch([tfbincenters(L), fliplr(tfbincenters(L))],[m(2,L)+sem(2,L), fliplr(m(2,L)-sem(2,L))],'red'); 
set(h,'Facealpha',.5);
plot(tfbincenters(L),m(2,L),'r-','Linewidth',2);
set(gca,'Xscale','log');
ylabel('d''');
xlabel('Temporal frequency (Hz)');
%title(['Cone ideal observer  (n = ',num2str(ns(1)),')'])
plot([tfbinedges(1) tfbinedges(end)],[1.27 1.27],'k-')

% Looking at cone vs LGN sensitivity as a function of RF position
TF_THRESHOLD = 5; % Hz
labels = {'RGlow','RGhigh','LUMlow','LUMhigh'};
figure;
for i = 1:4 % RGlow, RGhigh, LUMlow, LUMhigh
    subplot(2,2,i); hold on;
    for j = 1:length(data)
        h = plot(data{j}.rfxy(1),data{j}.rfxy(2),'ko','MarkerFaceColor','black');
        Lrg = sign(data{j}.uniquestim(:,1)) ~= sign(data{j}.uniquestim(:,2));
        Llum = sign(data{j}.uniquestim(:,1)) == sign(data{j}.uniquestim(:,2)) & data{j}.uniquestim(:,3) ~= 0;
        LlowTF = data{j}.uniquestim(:,3) < TF_THRESHOLD;
        LhighTF = data{j}.uniquestim(:,3) >= TF_THRESHOLD;
        if i == 1
            dprimes = [mean(data{j}.dprime(Lrg&LlowTF)) mean(data{j}.conedprime(Lrg&LlowTF))];
        elseif i == 2
            dprimes = [mean(data{j}.dprime(Lrg&LhighTF)) mean(data{j}.conedprime(Lrg&LhighTF))];
        elseif i == 3
            dprimes = [mean(data{j}.dprime(Llum&LlowTF)) mean(data{j}.conedprime(Llum&LlowTF))];
        else % i == 4
            dprimes = [mean(data{j}.dprime(Llum&LhighTF)) mean(data{j}.conedprime(Llum&LhighTF))];
        end
        if dprimes(1)>dprimes(2) % if LGN is more sensitive
            set(h,'MarkerFaceColor','red')
        end
        set(h,'MarkerSize',abs(diff(dprimes)))
    end
    title(labels{i});
end
set(gcf,'Name','Population inefficiency');

% The cone model is more efficient than the monkey for
% R/G at high eccentricities because, as Kathy Mullen showed, the decrease
% in RG sensitivity with eccentricity is postreceptoral. I'm seeing that
% populations of peripheral P-LGN neurons are relatively inefficient. This
% is unlikely to be due to the behavioral sensitivity model which just
% provides contrasts that were used identically in the LGN recordings and
% the cone model simulations. (Negative numbers in correlation plot mean
% cones have a high d'). Interpretation: post-receptoral loss in RG
% sensitivity occurs before Parvocellular neurons.

%%
% Section 4.2
% Comparing the sensitivity of the cone model to the sensitivity of the
% behavioral model. In DTcones_gh.m the runType = 'isosamp' option requires
% a stro structure so I'm loading one in and modifying it.
% Cone model doesn't return threshold. It must be found from fitting
% a psychometric function. I need a function that takes L, M, TF, x, and y
% and returns the cone model threshold. Like LMTF_thresh_from_model. Takes
% a long time.

% Grid search on x,y,theta, TF
%[xs,ys,thetas,tfs] = ndgrid([2:2:12],[-8:2:8],[pi/4 3*pi/4],[1 5 10]);
[xs,ys,thetas,tfs] = ndgrid(5,0,[pi/4 3*pi/4],[1 5]);
nconditions = numel(xs);

filename = 'A060817004.nex'; % Just to get the basic format of an IsoSamp file
stro = nex2stro(findfile(fullfile(nexfilepath,'Greg','Apollo',filename)));
examplestrotrial = stro.trial(1,:);
stro.trial = [];
stro.ras = {};

% Getting behavioral model parameters
whichmode = 5; % Which model of behavioral performance to use
LMTFdata = load(fullfile(fileparts(which('IsoSampOnline')), 'private', 'data', 'LMTF.mat'));
[~,fname] = fileparts(stro.sum.fileName);
sid = MONKEY(1);
LMTFstruct = getfield(LMTFdata,sid);
behavioral_model_params = getfield(LMTFstruct.legacy,['mode',num2str(whichmode),'params']);

% Setting up structure for DTcones_gh.m
clear params
params.runType = 'isosamp';
params.obsMethod = 'obsMethod_filteredWtFxn';
params.impulseResponse = 'rieke';
params.DTV1_fname = [];
params.DTNT_fname = [];
params.unitTest = false;
params.eqMosaic = false;
params.saveDir = '/Users/greghorwitz/Documents/MATLAB';
params.notes = 'IsoSamp test';
params.parallelOperations = false;
params.eyeType = 'monkey';
params.eyeNumber = 2;
params.coneSampRate = 2400;
params.flatPowerSpect = false;
params.enableScones = false;
params.sacamp_deg = 0;
params.sacdur_s = 0;
fundamentals = load('T_cones_smj10');
spds = stro.sum.exptParams.mon_spd;
spds = SplineSpd([380:4:780]',reshape(spds,length(spds)/3,3),[380:5:780]');
cal.monSpect = spds(:);
M = fundamentals.T_cones_smj10*spds;
cal.Mmtx = M(:);
cal.frameRate = stro.sum.exptParams.framerate;
cal.bkgndrgb = stro.sum.exptParams.bkgndrgb';
cal.fname = 'test';
cal.monSpectWavelengths = linspace(380,780,101);
cal.pixperdeg = stro.sum.exptParams.pixperdeg;
params.monCalFile = cal;

data = []; % ratios of model-predicted thresholds (behavioral/cone)
for condcounter = 1:nconditions
    rfx = xs(condcounter); 
    rfy = ys(condcounter); 
    TF = tfs(condcounter);
    Lccs = logspace(-3,log10(.1),10)*cos(thetas(condcounter));
    Mccs = logspace(-3,log10(.1),10)*sin(thetas(condcounter));
    
    % Putting new values in stro structure to be passed to DTcones
    stro.sum.exptParams.rf_x = rfx*10; % Charlie's cone model wants positions in DVA*10
    stro.sum.exptParams.rf_y = rfy*10; % Charlie's cone model wants positions in DVA*10
    for i = 1:length(Lccs)
        stro.trial(i,:) = examplestrotrial;
        stro.trial(i,strcmp(stro.sum.trialFields(1,:),'stim_l')) = Lccs(i);
        stro.trial(i,strcmp(stro.sum.trialFields(1,:),'stim_m')) = Mccs(i);
        stro.trial(i,strcmp(stro.sum.trialFields(1,:),'tf'))= TF;
    end
    params.stro = stro;

    [gab, cones, mon, idlob, params] = DTcones_gh(params,0);
    v = idlob.analyticVar{1}([1 2]); % variance
    m = []; % mean
    for i = 1:length(idlob.analyticMean); m(i,:) = idlob.analyticMean{i}([1 2]); end
    % Mean of noise distribution is (0,0)
    % Getting d-primes
    cone_dprimes = sqrt(m.^2*(1./v'));
    %cone_dprimes = nan*ones(size(m,1),1); % An alternative calculation that gives the same answer
    %for i = 1:size(m,1)
    %    cone_dprimes(i) = sqrt(m(i,:)*diag(1./v)*m(i,:)'); % no assumption of equal numbers of L/M cones
    %end
    % I want to find a contrast that will make d-prime = 1.27
    if (min(cone_dprimes) > 1.27 | max(cone_dprimes) < 1.27)
        error('Threshold not within bounds. Expand contrast bounds.');
    end
    cone_threshold = interp1(cone_dprimes, sqrt(Lccs.^2+Mccs.^2)', 1.27,'spline');
    
    local_params = LMTF_global_to_local_model(behavioral_model_params, rfx, rfy, whichmode);
    %local_params = stro.sum.exptParams.localparams; % Sanity checking
    pred_behavioral_threshold = LMTF_thresh_from_model(Lccs(1), Mccs(1), TF,local_params);
    [condcounter xs(condcounter) ys(condcounter) thetas(condcounter) tfs(condcounter) pred_behavioral_threshold./cone_threshold]
    data(condcounter) = pred_behavioral_threshold./cone_threshold;
end

% Plotting
if length(data) > 4
    whichtheta = 3*pi/4; whichtf = 1;
    L = thetas == whichtheta & tfs == whichtf;
    figure; axes; hold on;
    plotxs = unique(xs(L));
    plotys = unique(ys(L));
    zlabel('Cone/behavioral (sensitivity)');
    surf(plotxs,plotys,reshape(data(L),length(plotxs),length(plotys))');
    plot3(xs(L),ys(L),data(L),'o');
    set(gca,'Zscale','log');
    title(['theta = ',num2str(whichtheta),', TF = ',num2str(whichtf),', Max ratio: ',num2str(max(data(L))./min(data(L)))]);
    % Large numbers indicate that behavior is worse than cones
    xlabel('X'); ylabel('Y');
end

%%
% Section 4.3
% Comparing behavioral sensitivity to cone density

model_mode = 5; % Which model of behavioral performance to use
LMTFdata = load(fullfile(fileparts(which('IsoSampOnline')), 'private', 'data', 'LMTF.mat'));
[~,fname] = fileparts(stro.sum.fileName);
sid = MONKEY(1);
LMTFstruct = getfield(LMTFdata,sid);
behavioral_model_params = getfield(LMTFstruct.legacy,['mode',num2str(model_mode),'params']);

% Cone model
% Below taken from DTcones_gh.m (lines 834:840)
% More cones in nasal retina
temporalcoeffs=[150.9676 -1.2220 35.9979 -0.1567 9.9936 -0.0258]; % (Goodchild et al., 1996)
nasalcoeffs=[176.6624 -7.9473 94.4908 -0.3518 18.6761 -0.0236]; % (Packer et al. 1989) See GrantBrainStorming.m section 30
conedensfun = @(coeffs,x)(coeffs(1).*(exp(coeffs(2).*x)))+...
    (coeffs(3).*(exp(coeffs(4).*x)))+...
    (coeffs(5).*(exp(coeffs(6).*x)));

xs = 2:14;
ys = -8:8;
beh_mod_preds = zeros(length(ys),length(xs));
cone_mod_preds = zeros(length(ys),length(xs));
for i = 1:length(xs)
    for j = 1:length(ys)
        localparams = LMTF_global_to_local_model(behavioral_model_params, xs(i), ys(j), model_mode);
        beh_mod_preds(j,i) = LMTF_thresh_from_model(pi/4, 5,localparams);
        [theta,r] = cart2pol(xs(i),ys(j));
        coneDensity = (cos(theta).^2.*(conedensfun(temporalcoeffs,r)+conedensfun(nasalcoeffs,r))/2) + (sin(theta).^2.*conedensfun(temporalcoeffs,r));
        cone_mod_preds(j,i) = coneDensity;
    end
end

figure;
subplot(2,2,1);
imagesc(beh_mod_preds);
axis equal;
subplot(2,2,2);
imagesc(1./cone_mod_preds);
axis equal;
subplot(2,2,3);
thresh_times_cone_dens = beh_mod_preds.*cone_mod_preds; % Thresholds are high where cone density is low
imagesc(tmp);
axis equal;
subplot(2,2,4);
plot(thresh_times_cone_dens./min(thresh_times_cone_dens));

figure; axes; hold on;
surf(beh_mod_preds./min(beh_mod_preds))
surf(thresh_times_cone_dens./min(thresh_times_cone_dens));

% multiplying L-M thresholds by cone density does a pretty good job of
% flattening out the threhshold surface. Doesn't work as well for L+M
% (largely because the cone density is very high near the fovea and
% behavioral threshold aren't super low there).
%%
% Section 5
% Cone weights from white noise as a function of eccentricity.
% Only looking at a single pixel over a maximum of three frames.
maxT = 6;
[filenames,spikecds] = fnamesFromTxt('WhiteNoiseLGN_forIS','cellClass',{CELLTYPE},'subjID',{MONKEY(1)});

data = [];
for i = 1:length(filenames)
    stro = {};
    for j = 1:length(filenames{i})
        stro{j} = nex2stro(char(findfile(filenames{i}(j))));
    end
    stro = strocat(stro);
 
    noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
    sigmaidxs = strmatch('sigma',stro.sum.trialFields(1,:));
    spikename = ['sig001',char(double('a'+spikecds(i)-1))];
    nstixperside = stro.sum.exptParams.nstixperside;
    
    % Reconstructing the M matrix and gamma table
    fundamentals = stro.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = stro.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    gammaTable = stro.sum.exptParams.gamma_table;
    gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
    
    % Getting the background rgb/lms
    ridx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_r'));
    gidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_g'));
    bidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_b'));
    bkgndRGB = [mode(stro.trial(:,ridx)), mode(stro.trial(:,gidx)), mode(stro.trial(:,bidx))];
    bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
    bkgndlms = M*bkgndrgb;
    out = getWhtnsStats(stro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, spikename);
    STAs = out{1};
    energy = sum(STAs.^2);
    whichframe = find(energy == max(energy));
    
    energy = sum(reshape(STAs(:,whichframe).^2,nstixperside^2,3),2);
    whichpixel = find(energy == max(energy));
    
    STArgb = STAs(whichpixel+[0, nstixperside^2, 2*nstixperside^2],:);
    energy = sum(STArgb.^2);
    whichframe = find(energy == max(energy));
    whichframes = whichframe+[-1:1];
    whichframes(whichframes < 1 | whichframes > maxT) = [];
    STArgbtrunc = STArgb(:,whichframes);
    [u,~,~] = svd(STArgbtrunc);
    if u(:,1)'*STArgb(:,whichframe) < 0
        u(:,1) = -u(:,1);
    end
    STAlms = inv(M')*u(:,1); % cone excitation differences
    data = [data; u(:,1)' STAlms' stro.sum.exptParams.rf_x/10 stro.sum.exptParams.rf_y/10];
end

normrgbs = data(:,[1:3]) ./repmat(sum(abs(data(:,[1:3])),2),1,3);
normconeweights = data(:,[4:6]) ./repmat(sum(abs(data(:,[4:6])),2),1,3);
figure; axes; hold on;
for i = 1:size(data,1)
    h = plot(data(i,7),data(i,8),'o','MarkerSize',8);
    set(h,'MarkerFaceColor',normrgbs(i,:)/2+.5,'MarkerEdgeColor',normrgbs(i,:)/2+.5)
end
axis equal;
set(gca,'Color',[.5 .5 .5]);
set(gca,'Xlim',[-20 0],'Ylim',[-10 10]);

figure; axes; hold on;
L = normconeweights(:,3) > 0;
plot(normconeweights(~L,1),normconeweights(~L,2),'ko');
plot(normconeweights(L,1),normconeweights(L,2),'ko','MarkerFaceColor','black');
plot([-1 0 1 0 -1],[0 -1 0 1 0],'k-');
axis square;
xlabel('L cone weight');
ylabel('M cone weight');
title('filled symbols indicate positive S-cone weights');


%%
% Section 6
% Comparing the stimuli that were used in the IsoSamp experiment with the
% monkeys' psychophysical detection thresholds estimated from the *most
% recent* model. 
% I'm getting differences in threshold that are < 1.5-fold 
whichmode = 5; % 5 = yoked double-tilted rampy trough 
[filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',{CELLTYPE},'subjID',{MONKEY(1)});
data = [];
for i = 1:length(filenames)
    stro = {};
    for j = 1:length(filenames{i})
        stro{j} = nex2stro(char(findfile(filenames{i}(j), fullfile(nexfilepath,'Greg',MONKEY))));
    end
    stro = strocat(stro);
    rfx = stro.sum.exptParams.rf_x/10;
    rfy = stro.sum.exptParams.rf_y/10;
    Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_l'));
    Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_m'));
    TF = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
    uniquestim = sortrows(unique([Lcc Mcc TF],'rows'),3); % sorting by TF
    uniquestim(all(uniquestim == 0,2),:) = []; % Getting rid of the blank
    
    LMTFdata = load(fullfile(fileparts(which('IsoSampOnline')), 'private', 'data', 'LMTF.mat'));
    [~,fname] = fileparts(stro.sum.fileName);
    sid = fname(1);
    LMTFstruct = getfield(LMTFdata,sid);
    global_params = getfield(LMTFstruct.legacy,['mode',num2str(whichmode),'params']);
    local_params = LMTF_global_to_local_model(global_params, rfx, rfy, whichmode);
    %local_params = stro.sum.exptParams.localparams; % Sanity checking
    pred_r = LMTF_thresh_from_model(uniquestim,local_params);
    tested_r = sqrt(uniquestim(:,1).^2+uniquestim(:,2).^2);
    
    data = [data; log10(tested_r)-log10(pred_r)];
end
figure;
plot(10.^data,'.');
set(gca,'Yscale','log')
title(['mode ',num2str(whichmode)]);
10^mean(data)
10^std(data)
% Positive value mean tested contrasts were too high, which is the case
% for a few of Apollo's files. Fortunately, Apollo's neurons are actually
% slightly *less* sensitive (relative to the cones) than Utu's.
%%
% Section 7
% Noise spectra and other statistics
[filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',{CELLTYPE},'subjID',{MONKEY(1)});
timebins = [0:.005:.660];
DEBUG = 0; % pretend high frequency stimulation trials are blank trials
PLOTZSCORES = 0;
fftdata =[];
for a = 1:length(filenames)
    stro = {};
    for j = 1:length(filenames{a})
        stro{j} = nex2stro(char(findfile(filenames{a}(j), fullfile(nexfilepath,'Greg',MONKEY))));
    end
    stro = strocat(stro);
    Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_l'));
    Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_m'));
    TF = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
    dur = mean(stimoff_t-stimon_t);
    if DEBUG
        Lblank = TF > 30 & TF < 40; % For debugging
    else
        Lblank = Lcc == 0 & Mcc == 0 & TF == 0; 
    end
    tmpfftdata = [];
    for i=find(Lblank)'
        tmp = stro.ras{i,spikecds(a)}-stimon_t(i);
        tmp(tmp<0) = [];
        tmp(tmp>dur) = [];
        n = hist(tmp,timebins);
        tmpfftdata = [tmpfftdata; fftshift(fft(n))];
    end
    fftdata = [fftdata; mean(abs(tmpfftdata).^2*2/length(tmpfftdata))];
    
    [uniquestim, dprime, signal, noise] = IsoSampGetDPrime(stro, DPRIMEMETHOD, spikecds(a));
    zstats = [];
    bins = linspace(-4,4,30);
    counts = zeros(size(bins));
    if PLOTZSCORES
        figure; subplot(1,2,1); hold on;
        for i = 1:size(uniquestim,1)
            plot(uniquestim(i,3),noise{i},'b.');
            if sign(uniquestim(i,1)) == sign(uniquestim(i,2))
                plot(uniquestim(i,3),signal{i},'k.');
            else
                plot(uniquestim(i,3),signal{i},'r.');
            end
            zstats = [zstats; nanmean(signal{i}) nanstd(signal{i}) nanmean(noise{i}) nanstd(noise{i})];
            counts = counts+hist(signal{i}-mean(signal{i}),bins);
        end
        set(gca,'XScale','log');
        xlabel('TF (Hz)'); ylabel('decision variable');
        
        title(stro.sum.fileName(find(stro.sum.fileName == filesep,1,'last')+1:end))
        Llum = sign(uniquestim(:,1))==sign(uniquestim(:,2)) & ~uniquestim(:,1)==0;
        Lrg = sign(uniquestim(:,1))~=sign(uniquestim(:,2));
        plot(uniquestim(Llum,3), zstats(Llum,1),'k-','linewidth',2); % lum signal mean
        plot(uniquestim(Llum,3), zstats(Llum,2),'k:','linewidth',2); % lum signal std
        plot(uniquestim(Llum,3), zstats(Llum,3),'b-','linewidth',2); % noise mean
        plot(uniquestim(Llum,3), zstats(Llum,4),'b:','linewidth',2); % noise std
        plot(uniquestim(Lrg,3), zstats(Lrg,1),'r-','linewidth',2); % rg signal mean
        plot(uniquestim(Lrg,3), zstats(Lrg,2),'r:','linewidth',2); % rg signal std
        subplot(1,2,2); hold on;
        bar(bins,counts);
        xlabel('Mean-subtracted DV'); ylabel('Counts');
        drawnow();
    end
end
deltaT = timebins(2)-timebins(1);
nyquist = 1./(2*deltaT);
figure; subplot(2,1,1);
imagesc(abs(fftdata));
set(gca,'Xtick',[]);
subplot(2,1,2);
freqs = linspace(-nyquist,nyquist,size(fftdata,2));
plot(freqs,10*log10(mean(fftdata)));
set(gca,'Ylim',[-12 -5],'Xlim',[0 100]);
xlabel('Frequency (Hz)'); ylabel('Power (dB)');
title([MONKEY,' ',CELLTYPE]);
%%
% Section 8
% K-nearest neighbors in Jonathan Victor spike-train distance to look for
% signal in response to stimulus
Q = .01; % Cost for moving a spike
signaloffset = 0.05; % seconds past stimon_t to start counting spikes

DEBUG = 0;
[filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',{CELLTYPE},'subjID',{MONKEY(1)});
data = [];
for a = 1:length(filenames)
    stro = {};
    for j = 1:length(filenames{a})
        stro{j} = nex2stro(char(findfile(filenames{a}(j), fullfile(nexfilepath,'Greg',MONKEY))));
    end
    stro = strocat(stro);
    spikeidx = spikecds(a);
    Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_l'));
    Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_m'));
    TF = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
    if DEBUG % Randomizing the Lcc, Mcc, and TF values
        idxs = randperm(length(Lcc));
        Lcc = Lcc(idxs);
        Mcc = Mcc(idxs);
        TF = TF(idxs);
    end
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
    [uniquestim, dprime] = IsoSampGetDPrime(stro,DPRIMEMETHOD,spikeidx);
    
    spikes = {};
    startidxs = zeros(size(stro.ras,1),1);
    endidxs = zeros(size(stro.ras,1),1);
    for i = 1:size(stro.ras,1) % looping over trials
        tmpspikes = stro.ras{i,spikeidx}-stimon_t(i)-signaloffset;
        tmpspikes(tmpspikes<0) = [];
        tmpspikes(tmpspikes>stimoff_t(i)-stimon_t(i)) = [];
        if i == 1
            if isempty(tmpspikes)
                startidxs(i) = 1;
                endidxs(i) = 0;
            else
                startidxs(i) = 1;
                endidxs(i) = length(tmpspikes);
            end
        else
            if isempty(tmpspikes) % Make the end index 1 less than the start index?
                startidxs(i) = startidxs(i-1,:);
                endidxs(i) = startidxs(i,:)-1;
            end
            startidxs(i) = numel(vertcat(spikes{:}))+1;
            endidxs(i) = startidxs(i)+length(tmpspikes)-1;
        end
        spikes{i} = tmpspikes;
    end
    spiketimes = vertcat(spikes{:});
    
    % spkdl is blazingly fast so calculating *all* pairwise distances
    d = spkdl(spiketimes,startidxs',endidxs',Q); % 1/Q ms
    d_mat = reshape(d,length(startidxs),length(startidxs));
    
    Lblank = Lcc == 0 & Mcc == 0;
    nblanktrials = sum(Lblank);
    tmpdata = [];
    for i = 1:size(uniquestim,1)
        Lstim = Lcc == uniquestim(i,1) & Mcc == uniquestim(i,2) & TF == uniquestim(i,3);
        nsignaltrials = sum(Lstim);
        
        ntrials_tot = nblanktrials+nsignaltrials;
        sub_d_mat = [d_mat(Lblank,Lblank), d_mat(Lblank,Lstim);...
                     d_mat(Lstim,Lblank), d_mat(Lstim,Lstim)];
        sub_d_mat = sub_d_mat+diag(nan(ntrials_tot,1));
        
        % 8/1/18 Redoing. Taking average distances to all points.
        whichcategoryassigned = nanmean(sub_d_mat(1:nblanktrials,:)) > nanmean(sub_d_mat(nblanktrials+1:end,:)); % 0 means no stim, 1 means stim
        % If distances from blank trials are greater than distances from
        % stim trials, trials is a stim trial.
        tieidxs = nanmean(sub_d_mat(1:nblanktrials,:)) == nanmean(sub_d_mat(nblanktrials+1:end,:)); 
        if any(tieidxs)
            disp('Got a tie');
            whichcategoryassigned(tieidxs) = unidrnd(2,1,sum(tieidxs))-1;
        end
        %sorted_d_mat = sort(sub_d_mat); % Sorts each column in ascending order
        %k = (nsignaltrials/2)+1-mod(nsignaltrials/2,2); % an odd number of nearest neighbors to avoid ties
        %critical_values = sorted_d_mat(k,:); % K nearest neighbor
        %k_close_points = sub_d_mat<=repmat(critical_values,size(sub_d_mat,1),1);
        %mask = [-1*ones(nblanktrials,ntrials_tot); ones(nsignaltrials,ntrials_tot)];
        %whichcategoryassigned = sum(k_close_points.*mask)>0; % negative (0) means no stim, positive (1) means stim
        %if any(sum(k_close_points.*mask) == 0)
        %    disp('got a tie');
        %    if i ~= 1 & ~DEBUG
        %        %keyboard
        %    end
        %    tieidxs = sum(k_close_points.*mask) == 0;
        %    whichcategoryassigned(tieidxs) = unidrnd(2,1,sum(tieidxs))-1; % flipping a coin
        %end
        
        nmisclassified = sum(whichcategoryassigned(1:nblanktrials) ~= 0) + sum(whichcategoryassigned(nblanktrials+1:ntrials_tot) == 0);
        % Each term in the nmisclassified sum, above, is biased upward because there are more "unlike" trials than "like" trials
        % The expected value of each term in this sum is
        %E_nmis = nblanktrials*hygestat(ntrials_tot-1,nblanktrials-1,k)/k + nsignaltrials*hygestat(ntrials_tot-1,nsignaltrials-1,k)/k;
        %E_cor = nblanktrials*hygestat(ntrials_tot-1,nsignaltrials,k)/k + nsignaltrials*hygestat(ntrials_tot-1,nblanktrials,k)/k;
        %tmpdata = [tmpdata; nmisclassified*E_nmis/(ntrials_tot/2) ncorrectlyclassified*E_cor/(ntrials_tot/2)];% Correcting for the bias(?)
        ncorrectlyclassified = ntrials_tot-nmisclassified;
        tmpdata = [tmpdata; nmisclassified ncorrectlyclassified];
    end
    data{a}.filenames = filenames{a};
    data{a}.spikecds = spikecds(a);
    data{a}.rfxy = [stro.sum.exptParams.rf_x/10, stro.sum.exptParams.rf_y/10];
    data{a}.uniquestim = uniquestim;
    data{a}.prop_correct_classification = tmpdata(:,2)./(tmpdata(:,1)+tmpdata(:,2));
    data{a}.conventional_dprime = dprime;
end

tfbinedges = logspace(0,log10(60),8);
tfbincenters = sqrt(tfbinedges(1:end-1).*tfbinedges(2:end)); % geomean
classification_performance = cell(2,length(tfbincenters));
ns = zeros(2,length(tfbincenters));
for i = 1:length(data)
    in_ecc_range = sqrt(data{i}.rfxy*data{i}.rfxy') >= 2 & abs(data{i}.rfxy(1)) <= 12 & abs(data{i}.rfxy(2))<= 8;
    if ~in_ecc_range
        disp('chucking a cell')
        continue
    end
    for STIMTYPE = {'LUM','RG'}
        stimtypeidx = find(strcmp(STIMTYPE, {'LUM','RG'}));
        if strcmp(STIMTYPE,'LUM')
            Lcolordir = sign(data{i}.uniquestim(:,1)) == sign(data{i}.uniquestim(:,2)) & data{i}.uniquestim(:,3) ~= 0;
        elseif strcmp(STIMTYPE,'RG')
            Lcolordir = sign(data{i}.uniquestim(:,1)) ~= sign(data{i}.uniquestim(:,2));
        else
            error;
        end
        
        for j = 1:length(tfbincenters)
            Ltf = data{i}.uniquestim(:,3) > tfbinedges(j) & data{i}.uniquestim(:,3) <= tfbinedges(j+1);
            if sum(Ltf&Lcolordir) > 0
                classification_performance{stimtypeidx,j} = [classification_performance{stimtypeidx,j};...
                    mean(data{i}.prop_correct_classification(Ltf&Lcolordir))];               
                ns(stimtypeidx,j) = ns(stimtypeidx,j) + 1; % each cell counts as an independent entity
            end
        end
    end
end

m = zeros(size(classification_performance));
sem = zeros(size(classification_performance));

for i = 1:numel(classification_performance)
    m(i) = mean(classification_performance{i});
    sem(i) = sqrt(var(classification_performance{i})./ns(i));
    tmp = [];
    for j = 1:length(classification_performance{i})
        tmp(j) = norminv(classification_performance{i}(j),0,sqrt(2));
    end
end

% Mean and SEMs of "raw" classification percentages
colors = {'black','red'};
figure; axes; hold on;
for i = 1:2
    L = ns(i,:) > 2; % minimum 'n'
    h = plot(tfbincenters(L), m(i,L),'k-o','LineWidth',2);
    set(h,'Color',colors{i},'MarkerFaceColor',colors{i});
    h = patch([tfbincenters(L), fliplr(tfbincenters(L))],[m(i,L)+sem(i,L), fliplr(m(i,L)-sem(i,L))],colors{i},'Facealpha',.5);
end
set(gca,'Xscale','log');
ylabel('Proportion correct classification');
xlabel('Frequency (Hz)');
plot([tfbincenters(1) tfbincenters(end)],[.5 .5],'k-');
axis square;

% Adding a second axis for "non-parametric d-prime".
ylims = get(gca,'Ylim');
a2 = axes('YAxisLocation', 'Right');
axis square;
set(a2, 'color', 'none','XTick', [], 'YLim', ylims)
ticks_in_proportions = linspace(ylims(1),ylims(2),5);
set(a2, 'Ytick', ticks_in_proportions);
set(a2,'Yticklabel',round(norminv(ticks_in_proportions,0,sqrt(2))*100)/100)
ylabel('Equivalent d-prime');

% Comparing misclassification rates between KNN and conventional d-prime estimates
rs = [];
diffs = [];
for a = 1:length(data)
    dprimes = data{a}.conventional_dprime;
    dprime_based_propcorrect = zeros(length(dprimes),1);
    for i = 1:length(dprimes)
        x = linspace(-6,6+dprimes(i),1000);
        dx = x(2)-x(1);
        y = normpdf(x,dprimes(i),1);
        dprime_based_propcorrect(i) = sum(normcdf(x,0,1).*(y./sum(y)));
    end
    L = ~all(data{a}.uniquestim == 0,2);  
    tmp = corrcoef(data{a}.prop_correct_classification(L), dprime_based_propcorrect(L));
    rs(a) = tmp(2,1);
    diffs(a) = median(dprime_based_propcorrect(L)-data{a}.prop_correct_classification(L));
end
figure;
subplot(2,1,1); hold on;
hist(rs); axis square;
xlabel('Correlation btn KNN and dprime correct classification rates');
ylabel('count');
subplot(2,1,2); hold on;
hist(diffs); axis square;
xlabel('dprime - KNN correct classification rates');
ylabel('count');

%%
% Section 9) Comparing cone weights from white noise to IsoSamp d-prime
% Computing cone weights from the single hottest pixel only right now.
HIGHLOWTHRESHOLD = 5;
[filenames_Isosamp,spikecds_Isosamp,~,neuron_id_Isosamp] = fnamesFromTxt('IsoSamp_LGN','cellClass',{CELLTYPE},'subjID',{MONKEY(1)});
[filenames_WN,spikecds_WN,~,neuron_id_WN] = fnamesFromTxt('WhiteNoiseLGN_forIS','cellClass',{CELLTYPE},'subjID',{MONKEY(1)});
if length(neuron_id_Isosamp) ~= length(filenames_Isosamp)
    keyboard
end
if length(neuron_id_WN) ~= length(filenames_WN)
    keyboard
end
valid_ids = intersect(neuron_id_Isosamp, neuron_id_WN);
data = {};
for i = valid_ids'
    Isosamp_idx = find(neuron_id_Isosamp == i);
    WN_idx = find(neuron_id_WN == i);

    % First processing IsoSamp file(s)
    stro_Isosamp = {};
    for j = 1:length(filenames_Isosamp{Isosamp_idx})
        stro_Isosamp{j} = nex2stro(char(findfile(filenames_Isosamp{Isosamp_idx}(j), fullfile(nexfilepath,'Greg',MONKEY))));
    end
    stro_Isosamp = strocat(stro_Isosamp);
    [uniquestim, dprime] = IsoSampGetDPrime(stro_Isosamp,DPRIMEMETHOD,spikecds_Isosamp(Isosamp_idx));

    % First processing white noise file(s)
    stro_WN = {};
    for j = 1:length(filenames_WN{WN_idx})
        stro_WN{j} = nex2stro(char(findfile(filenames_WN{WN_idx}(j))));
    end
    stro_WN = strocat(stro_WN);
    noisetypeidx = find(strcmp(stro_WN.sum.trialFields(1,:),'noise_type'));
    sigmaidxs = strmatch('sigma',stro_WN.sum.trialFields(1,:));
    spikename = ['sig001',char(double('a'+spikecds_WN(WN_idx)-1))];
    nstixperside = stro_WN.sum.exptParams.nstixperside;
    % Reconstructing the M matrix and gamma table
    fundamentals = stro_WN.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = stro_WN.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences

    gammaTable = stro_WN.sum.exptParams.gamma_table;
    gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
    % Getting the background rgb/lms
    ridx = find(strcmp(stro_WN.sum.trialFields(1,:),'bkgnd_r'));
    gidx = find(strcmp(stro_WN.sum.trialFields(1,:),'bkgnd_g'));
    bidx = find(strcmp(stro_WN.sum.trialFields(1,:),'bkgnd_b'));
    bkgndRGB = [mode(stro_WN.trial(:,ridx)), mode(stro_WN.trial(:,gidx)), mode(stro_WN.trial(:,bidx))];
    bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
    bkgndlms = M*bkgndrgb;
    out = getWhtnsStats(stro_WN,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, spikename);
    STAs = out{1};
    energy = sum(STAs.^2);
    whichframe = find(energy == max(energy));
    
    energy = sum(reshape(STAs(:,whichframe).^2,nstixperside^2,3),2);
    whichpixel = find(energy == max(energy));
    
    STArgb = STAs(whichpixel+[0, nstixperside^2, 2*nstixperside^2],:);
    energy = sum(STArgb.^2);
    whichframe = find(energy == max(energy));
    whichframes = whichframe+[-2:2];
    whichframes(whichframes < 1 | whichframes > maxT) = [];
    STArgbtrunc = STArgb(:,whichframes);
    [u,~,~] = svd(STArgbtrunc);
    if u(:,1)'*STArgb(:,whichframe) < 0
        u(:,1) = -u(:,1);
    end
    STAlms = inv(Mrgbtocc')*u(:,1); % weights on cone contrasts
    idx = length(data)+1;
    data{idx}.coneweights = STAlms./sum(abs(STAlms));
    data{idx}.uniquestim = uniquestim;
    data{idx}.dprime = dprime;
end

labels = {'Lum high','Lum low','RG high','RG low'};
figure; axes; hold on;
for i = 1:length(data)
    uniquestim = data{i}.uniquestim;
    dprime = data{i}.dprime;
    Lblank = uniquestim(:,1) == 0 & uniquestim(:,2) == 0 & uniquestim(:,3) == 0;
    Lrg = sign(uniquestim(:,1)) ~= sign(uniquestim(:,2));
    Llum = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & ~Lblank;
    LlowTF = uniquestim(:,3) < HIGHLOWTHRESHOLD;
    LhighTF = uniquestim(:,3) >= HIGHLOWTHRESHOLD;
    for j = 1:4
        subplot(2,2,j); hold on;
        if j == 1
            L = Llum&LhighTF;
        elseif j == 2
            L = Llum&LlowTF;
        elseif j == 3
            L = Lrg&LhighTF;
        else % j == 4
            L = Lrg&LlowTF;
        end
        h = plot(data{i}.coneweights(1),data{i}.coneweights(2),'ko');
        markersize = mean(dprime(L));
        if markersize > 0
            set(h,'MarkerSize',markersize*20);
        else
            set(h,'MarkerSize',5,'Marker','x');
        end
    end
end
for i = 1:4
    subplot(2,2,i);
    plot([0 1 0 -1 0],[-1 0 1 0 -1],'k-');
    axis square;
    title(labels{i});
end

%%
% Section 10) Looking at response modulations in IsoSamp around the time of fixational
% eye movements.

[filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',{CELLTYPE},'subjID',{MONKEY(1)});
data = [];
spiketimes_wrt_sac = {};
PREOFFSET = 0.05; % How long before saccade initiation to collect spikes (s)
POSTOFFSET = 0.2; % How long after saccade initiation to collect spikes (s)
for i = 1:length(filenames)
    stro = {};
    for j = 1:length(filenames{i})
        stro{j} = nex2stro(char(findfile(filenames{i}(j), fullfile(nexfilepath,'Greg',MONKEY))));
    end
    stro = strocat(stro);
    ntrials = size(stro.trial,1);
    Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_l'));
    Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_m'));
    TF = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
    fpacq_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpacq_t'));
    spiketimes = stro.ras(:,strcmp(stro.sum.rasterCells,'sig001a')); % Only looking at first spike
    sacstats = getSacData(stro,1);
    for j = 1:ntrials
        nsacs = size(sacstats.starttimes{j},1);
        L = ones(nsacs,1);
        L(sacstats.starttimes{j} < fpacq_t(j)) = 0;
        L(sacstats.starttimes{j} > stimoff_t(j)) = 0;
        if any(L)
            for k = find(L')
                data = [data; i j Lcc(j) Mcc(j) TF(j) sacstats.amplitudes{j}(k) sacstats.directions{j}(k) (sacstats.starttimes{j}(k) < stimon_t(j) | (Lcc(j) == 0 & Mcc(j) == 0))];
                tmpspiketimes = spiketimes{j};
                Lspikes = tmpspiketimes > sacstats.starttimes{j}(k)-PREOFFSET & tmpspiketimes < sacstats.starttimes{j}(k)+POSTOFFSET;
                spiketimes_wrt_sac{size(data,1)} = tmpspiketimes(Lspikes)-sacstats.starttimes{j}(k);
            end
        end
    end
end
% Column labels
FILEIDX = 1;
SACIDX = 2;
LCCIDX = 3;
MCCIDX = 4;
TFIDX = 5;
AMPIDX = 6;
DIRIDX = 7;
BLANKSCREEN = 8;

deltaT = .005; % To get 100 Hz nyquist
bins = -PREOFFSET:deltaT:POSTOFFSET;
binwidth = bins(2)-bins(1); % s
figure; axes; hold on;
for whichcond = 0:1
    PSTH = zeros(1,length(bins));
    nsacs = 0;
    for i = 1:length(spiketimes_wrt_sac)
        if data(i,BLANKSCREEN) == whichcond
            [n,~] = hist(spiketimes_wrt_sac{i}, bins);
            PSTH = PSTH+n;
            nsacs = nsacs+1;
        end
    end
    meanPSTH = PSTH/binwidth/nsacs;
    plot(bins(2:end-1), meanPSTH(2:end-1));
end
xlabel('time (s)');
ylabel('firing rate (spikes/sec)');
title([MONKEY,' ',CELLTYPE]);

% Looking at power spectrum of msac-triggered PSTH
nyquist = 1./(2*deltaT);
figure; axes; hold on;
freqs = linspace(-nyquist,nyquist,length(meanPSTH(2:end-1)));
tmp = fftshift(abs(fft(meanPSTH(2:end-1))));
plot(freqs,10*log10(tmp),'ko-','LineWidth',1,'MarkerSize',6);
set(gca,'Xlim',[0 nyquist]);
xlabel('Frequency (Hz)');
ylabel('Power (dB)');


%%
% Section 10.1: Frequency of microsaccade occurrence as a function of time.
% Also power spectra.
% Also amplitudes and durations

filenames = fnamesFromTxt('IsoSamp_LGN','subjID',{MONKEY(1)});
PREOFFSET = -0.1; % How long before stimulus on to collect saccades
POSTOFFSET = 0.1; % How long after stimulus off to collect saccades
deltaT = 0.005;
timebins = PREOFFSET:deltaT:.666+POSTOFFSET;
ntimebins = length(timebins);
timebinwidth = timebins(2)-timebins(1);
deltaT = timebins(2)-timebins(1);
nyquist = 1./(2*deltaT);

means_lum = nan*ones(length(filenames),ntimebins);
means_blank = nan*ones(length(filenames),ntimebins);
power_lum = nan*ones(length(filenames),ntimebins);
power_blank = nan*ones(length(filenames),ntimebins);
amp_and_dur = [];
for i = 1:length(filenames)
    stro = {};
    for j = 1:length(filenames{i})
        stro{j} = nex2stro(char(findfile(filenames{i}(j), fullfile(nexfilepath,'Greg',MONKEY))));
    end
    stro = strocat(stro);
    ntrials = size(stro.trial,1);
    Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_l'));
    Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_m'));
    TF = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
    Lblank = Lcc == 0 & Mcc == 0 & TF == 0; 
    Llum = Lcc > 0 & Mcc > 0;
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
    sacstats = getSacData(stro,1);
    for k = 0:1
        if k == 1
            whichtrials = find(Llum)';
        else
            whichtrials = find(Lblank)';
        end
        tmpdata = zeros(length(whichtrials),ntimebins);
        rowcounter = 1;
        for j = whichtrials
            sactimes = sacstats.starttimes{j}-stimon_t(j);
            L = logical(ones(length(sactimes),1));
            L(sactimes<0+PREOFFSET-timebinwidth/2) = 0;
            L(sactimes>.666+POSTOFFSET+timebinwidth/2) = 0;
            [n,x] = hist(sactimes(L), timebins);
            tmpdata(rowcounter,:) = n;
            rowcounter = rowcounter+1;
            amp_and_dur = [amp_and_dur; sacstats.amplitudes{j}(L) sacstats.durations{j}(L)];
        end
        if all(Lcc(whichtrials) == 0) && all(Mcc(whichtrials) == 0)
            means_blank(i,:) = mean(tmpdata);
            power_blank(i,:) = mean(abs(fftshift(fft(tmpdata,[],2))));
        else
            means_lum(i,:) = mean(tmpdata);
            power_lum(i,:) = mean(abs(fftshift(fft(tmpdata,[],2))));
        end
    end    
end

figure; axes; hold on;
for i = 0:1
    if i == 0
        tmp = means_blank./timebinwidth;
    else
        tmp = means_lum./timebinwidth;
    end
    sem = sqrt(var(tmp)./size(means,1));
    mn = mean(tmp);
    patch([timebins, fliplr(timebins)],[mn+sem, fliplr(mn-sem)],[.5 .5 .5],'Facealpha',.5);
    plot(timebins,mn,'k-','LineWidth',2);
    xlabel('Time(s)');
    ylabel('saccades/s');
end

figure; axes; hold on;
for i = 0:1
    if i == 0
        tmp = power_blank;
    else
        tmp = power_lum;
    end
    freqs = linspace(-nyquist,nyquist,size(tmp,2));
    plot(freqs,10*log10(mean(tmp)),'ko-','LineWidth',1,'MarkerSize',6);
    set(gca,'Xlim',[0 nyquist]);
    xlabel('Frequency (Hz)'); ylabel('Power (dB)');
    title(MONKEY);
end

%%
% Section 11: Effect of counting window on the d-primes of LGN
% neurons and cone signals. Option to "population scale"

[filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',{CELLTYPE},'subjID',{MONKEY(1)});
HIGHLOWTHRESHOLD = 5;
TEMPORONASALSCALEFACTOR = 0.8;

POPSCALE = 1;
freqs = logspace(log10(1), log10(60),15);
period = 1./freqs;
offsets = [.333-period'/2, .333+period'/2];
L = offsets(:,1)<0 | offsets(:,2) > .66;
freqs(L) = [];
period(L) = [];
offsets(L,:) = [];

% Setting up for cone model
fundamentals = load('T_cones_smj10');
params = [];
params.runType = 'isosamp';
params.obsMethod = 'obsMethod_filteredWtFxn';
params.impulseResponse = 'rieke';
params.DTV1_fname = [];
params.DTNT_fname = [];
params.unitTest = false;
params.eqMosaic = false;
params.saveDir = '/Users/greghorwitz/Documents/MATLAB';
params.notes = 'IsoSamp test';
params.parallelOperations = false;
params.eyeType = 'monkey';
if POPSCALE
    params.eyeNumber = 2;
else
    params.eyeNumber = 1; % comparing sensitivity of a monocular neuron to cones
end
params.coneSampRate = 2400;
params.flatPowerSpect = false;
params.enableScones = false;
params.sacamp_deg = 0;
params.sacdur_s = 0;

%data = [];
lumdata = [];
for a = 1:length(filenames)
    stro = {};
    for j = 1:length(filenames{a})
        stro{j} = nex2stro(char(findfile(filenames{a}(j), fullfile(nexfilepath,'Greg',MONKEY))));
    end
    stro = strocat(stro);
    uniquestim = IsoSampGetDPrime(stro);    
    Lblank = uniquestim(:,1) == 0 & uniquestim(:,2) == 0 & uniquestim(:,3) == 0;
    Lrg = sign(uniquestim(:,1)) ~= sign(uniquestim(:,2));
    Llum = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & ~Lblank;
    LlowTF = uniquestim(:,3) < HIGHLOWTHRESHOLD;
    LhighTF = uniquestim(:,3) >= HIGHLOWTHRESHOLD;

    for offsetcounter = 1:size(offsets,1)
        timewindoffset = offsets(offsetcounter,:);
        [~, dprime] = IsoSampGetDPrime(stro, DPRIMEMETHOD, spikecds(a), timewindoffset);
        if POPSCALE
            pop_scale_fact = IsoSampGetPopulationScaleFactor(stro,ecc_to_diam_deg);
            dprime=dprime*pop_scale_fact;
        end

        params.temporalIntegrationLims = timewindoffset; % truncating the temporal integration window for an ideal observer of cone currents
        params.stro = stro;
        spds = params.stro.sum.exptParams.mon_spd;
        spds = SplineSpd([380:4:780]',reshape(spds,length(spds)/3,3),[380:5:780]');
        cal.monSpect = spds(:);
        M = fundamentals.T_cones_smj10*spds;
        cal.Mmtx = M(:);
        cal.frameRate = params.stro.sum.exptParams.framerate;
        cal.bkgndrgb = params.stro.sum.exptParams.bkgndrgb'; % intensities
        cal.fname = 'test';
        cal.monSpectWavelengths = linspace(380,780,101);
        cal.pixperdeg = params.stro.sum.exptParams.pixperdeg;
        params.monCalFile = cal;
        rfx = stro.sum.exptParams.rf_x/10;
        rfy = stro.sum.exptParams.rf_y/10;
        rf_r_deg = sqrt((rfx./TEMPORONASALSCALEFACTOR)^2+rfy^2);
        rf_diam_deg = ecc_to_diam_deg(rf_r_deg)/2; % size of RF in SD
        sigma_gabor = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sigma'))); % DVA
        if POPSCALE
            params.gab.sd = sigma_gabor;
        else
            params.gab.sd = rf_diam_deg/2;
        end
        % End of cone model set up
        [gab, cones, mon, idlob, params] = DTcones_gh(params,0);
        
        conemodeldata = nan*ones(size(uniquestim,1),1);
        for i = 1:size(idlob.analyticMean,1) % looping over color direction
            for j = 1:size(idlob.analyticMean(i,:),2) % looping over contrast/TF
                if ~isempty(idlob.analyticMean{i,j})
                    tmp_lm_mu = idlob.analyticMean{i,j}([1 2]);
                    tmp_lm_var = idlob.analyticVar{i,j}([1 2]);
                    tf = gab.driftRates{i}(j);
                    L = uniquestim(:,3) == gab.driftRates{i}(j) & sign(gab.colorDirs(i,1)) == sign(uniquestim(:,1));
                    if (sum(L) == 0)
                        error('cannot find condition');
                    end
                    conemodeldata(L) = [sum(abs(tmp_lm_mu))/sqrt(sum(tmp_lm_var))];
                end
            end
        end
        
        % Need to calculate flashTimeProfile on the basis of timewindoffset
        t = 0:1/mon.frameRate:gab.length;
        nframes = length(t);
        flashTimeProfile = ones(1,nframes);
        ramp = linspace(0,1,nframes/4);
        flashTimeProfile(1:length(ramp)) = ramp;
        flashTimeProfile(end:-1:end-length(ramp)+1) = ramp;
        flashTimeProfile = flashTimeProfile(t>=min(timewindoffset) & t<=max(timewindoffset));
        sigmaInPix = params.gab.sd*cal.pixperdeg;

        photondprime = IsoSampGetPhotonDPrime (flashTimeProfile, mon.frameRate, mon.bkgndlms_Rstar, sigmaInPix, cat(3,cones.num_L,cones.num_M), uniquestim);

       % data(a,offsetcounter,1,:) = [nanmean(dprime(Lrg&LlowTF)) nanmean(dprime(Lrg&LhighTF))  nanmean(dprime(Llum&LlowTF)) nanmean(dprime(Llum&LhighTF))];
       % data(a,offsetcounter,2,:) = [nanmean(conemodeldata(Lrg&LlowTF)) nanmean(conemodeldata(Lrg&LhighTF))  nanmean(conemodeldata(Llum&LlowTF)) nanmean(conemodeldata(Llum&LhighTF))];
        if offsetcounter == 1
            lumdata{a} = [uniquestim(Llum,3) dprime(Llum) conemodeldata(Llum) photondprime(Llum)];
        else
            lumdata{a} = [lumdata{a} dprime(Llum) conemodeldata(Llum) photondprime(Llum)];
        end
    end
end

% Sanity checking cone model dprimes
tmp = lumdata{end};
figure;
timewindowdurs = offsets(:,2)-offsets(:,1);
plot(tmp(:,4:3:end))


% lumdata has n cells (one per neuron)
% within each cell, the first column are temporal frequencies
% Then, each pair of columns is dprime from the neuron followed by d-prime
% from the cone current model.
% The number of pairs of columns are the number of counting windows

tmp = nan*ones(length(lumdata),length(freqs),3); % cell, frequency, [LGN,cone, photon]
for i = 1:length(lumdata) % looping over cells
    singlefilefreqs = lumdata{i}(:,1);
    for j = 1:length(freqs) % only looking at frequencies with at least one full cycle
        L = softEq(singlefilefreqs,freqs(j),1);
        if sum(L) == 1
            tmp(i,j,1) = lumdata{i}(L,3*(j-1)+1+1); % Second "+1" is to skip the first column which is TF
            tmp(i,j,2) = lumdata{i}(L,3*(j-1)+2+1);
            tmp(i,j,3) = lumdata{i}(L,3*(j-1)+3+1);
        elseif sum(L) > 1
            error('too many matches');
        end
    end
end

figprefs; axes;
for i = 1:3 % LGN, cone current, photon absorption
    mn = nanmean(squeeze(tmp(:,:,i)));
    sd = nanstd(squeeze(tmp(:,:,i)));
    Lbins = ~isnan(sd);
    sem = sd./sqrt(sum(~isnan(squeeze(tmp(:,:,i)))));
    h1 = patch([freqs(Lbins), fliplr(freqs(Lbins))],[mn(Lbins)+sem(Lbins), fliplr(mn(Lbins)-sem(Lbins))],[.5 .5 .5],'Facealpha',.5,'LineStyle','none');
    h2 = plot(freqs(Lbins), mn(Lbins),'ko-','linewidth',2,'MarkerFaceColor','black');
    if i == 2
        set(h1,'FaceColor',[1 .2 .2]);
        set(h2,'Color','red','MarkerEdgeColor','red','MarkerFaceColor','red');
    elseif i == 3
        set(h1,'FaceColor',[.2 .2 1]);
        set(h2,'Color','cyan','MarkerEdgeColor','cyan','MarkerFaceColor','cyan');
    end
end
set(gca,'xscale','log','Xlim',[1 60],'Ylim',[-1 6]);
title([MONKEY,' ',CELLTYPE]);
if POPSCALE
    plot([freqs(1) freqs(end)],[1.27 1.27],'k--');
    set(gca,'Ylim',[-5 10]);
end

%%
% Section 12
% Sliding time window analysis of d'
% Makes the most sense for high frequency stimuli. 
% Basically asking if the effects of contrast adaptation are large enough
% to affect d'.
[filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',{CELLTYPE},'subjID',{MONKEY(1)});

countingwindowdur = .4;
signaloffsets = linspace(-countingwindowdur,.6,12)';
signaloffsets = [signaloffsets signaloffsets+countingwindowdur]; % sliding window
%signaloffsets = [zeros(10,1) linspace(.02,.66,10)']; % stretching window
TFBINEDGES = logspace(0,log10(61),8);
TFBINCENTERS = sqrt(TFBINEDGES(1:end-1).*TFBINEDGES(2:end)); % geomean

lumdata = [];
rgdata = [];

for i = 1:length(filenames)
    stro = {};
    for j = 1:length(filenames{i})
        stro{j} = nex2stro(char(findfile(filenames{i}(j), fullfile(nexfilepath,'Greg',MONKEY))));
    end
    stro = strocat(stro);
    pop_scale_fact = IsoSampGetPopulationScaleFactor(stro,ecc_to_diam_deg);
    tmp = [];
    for j = 1:size(signaloffsets,1)
        [uniquestim, dprime] = IsoSampGetDPrime(stro, DPRIMEMETHOD, spikecds(i), signaloffsets(j,:));
        tmp = [tmp,dprime*pop_scale_fact];
    end
    Lblank = uniquestim(:,1) == 0 & uniquestim(:,2) == 0 & uniquestim(:,3) == 0;
    Lrg = sign(uniquestim(:,1)) ~= sign(uniquestim(:,2));
    Llum = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & ~Lblank;
    for j = 1:length(TFBINEDGES)-1
        Ltf = uniquestim(:,3)>=TFBINEDGES(j) & uniquestim(:,3)<TFBINEDGES(j+1);
        lumdata(i,:,j) = mean(tmp(Ltf&Llum,:),1);
        rgdata(i,:,j) = mean(tmp(Ltf&Lrg,:),1);
    end
end

lummeans = squeeze(nanmean(lumdata));
rgmeans = squeeze(nanmean(rgdata));

figure;
subplot(2,2,1)
surface(lummeans);
set(gca,'Ytick',[1:size(signaloffsets,1)],'YtickLabel',round(((signaloffsets(:,2)-signaloffsets(:,1))/2)*100)/100);
set(gca,'Xtick',[1:length(TFBINCENTERS)],'XtickLabel',round(TFBINCENTERS*10)/10);
xlabel('Frequency (Hz)');
ylabel('Middle time (s)');
title([CELLTYPE,' ',MONKEY]);

subplot(2,2,2)
surface(rgmeans);
set(gca,'Ytick',[1:size(signaloffsets,1)],'YtickLabel',round(((signaloffsets(:,2)-signaloffsets(:,1))/2)*100)/100);
set(gca,'Xtick',[1:length(TFBINCENTERS)],'XtickLabel',round(TFBINCENTERS*10)/10);
xlabel('Frequency (Hz)');
ylabel('Middle time (s)');
title(['RG min freq = ',num2str(1./countingwindowdur),' Hz']);

subplot(2,2,3)
plot(lummeans);
set(gca,'Xtick',1:length(signaloffsets),'XtickLabel',round(((signaloffsets(:,2)-signaloffsets(:,1))/2)*10)/10,'Xscale','linear');
xlabel('Middle time (s)');
ylabel('SNR');
hold on;
plot([0 length(signaloffsets)],[1.27 1.27],'k-');

subplot(2,2,4)
plot(rgmeans);
set(gca,'Xtick',1:length(signaloffsets),'XtickLabel',round(((signaloffsets(:,2)-signaloffsets(:,1))/2)*10)/10,'Xscale','linear');
xlabel('Middle time (s)');
ylabel('Normalized SNR');

%%
% Section 13: Ideal observer of a single cycle in the middle of the
% response.
% A reasonable center point appears to be 0.35 (20 ms past peak)
% Didn't work that well: very noisy and dependent on "centertime".

[filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',{CELLTYPE},'subjID',{MONKEY(1)});
TFBINEDGES = logspace(0,log10(61),12);
TFBINCENTERS = sqrt(TFBINEDGES(1:end-1).*TFBINEDGES(2:end)); % geomean

lumdata = [];
rgdata = [];
for i = 1:length(filenames)
    stro = {};
    for j = 1:length(filenames{i})
        stro{j} = nex2stro(char(findfile(filenames{i}(j), fullfile(nexfilepath,'Greg',MONKEY))));
    end
    stro = strocat(stro);
    uniquestim = IsoSampGetDPrime(stro);
    pop_scale_fact = IsoSampGetPopulationScaleFactor(stro,ecc_to_diam_deg);
    Lblank = uniquestim(:,1) == 0 & uniquestim(:,2) == 0 & uniquestim(:,3) == 0;
    Llum = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & ~Lblank;
    Lrg = sign(uniquestim(:,1)) ~= sign(uniquestim(:,2));
    TFs = unique(uniquestim(~Lblank,3));

    tmp_lum = [];  tmp_rg = [];
    for j = 1:size(TFs,1)
        centertime = 0.3;
        signaloffset = centertime+(1/(2*TFs(j)))*[-1 1]; % one cycle-long window centered on "centertime"
        signaloffset(1) = max(signaloffset(1),0);
        signaloffset(2) = min(signaloffset(2),.666);
        [~, dprime] = IsoSampGetDPrime(stro, DPRIMEMETHOD, spikecds(i), signaloffset);
        dprime = dprime*pop_scale_fact;
        LTF = uniquestim(:,3) == TFs(j);
        if any(Llum&LTF)
            tmp_lum = [tmp_lum;dprime(Llum&LTF)];
        else
            tmp_lum = [tmp_lum;nan];
        end
        if any(Lrg&LTF)
            tmp_rg = [tmp_rg;dprime(Lrg&LTF)];
        else
            tmp_rg = [tmp_rg;nan];
        end
    end
    
    for j = 1:length(TFBINEDGES)-1
        Ltf = TFs>=TFBINEDGES(j) & TFs<TFBINEDGES(j+1);
        lumdata(i,j) = nanmean(tmp_lum(Ltf,:),1);
        rgdata(i,j) = nanmean(tmp_rg(Ltf,:),1);
    end
end

% HACK!
if any(any(abs(lumdata)> 1e15))
    disp('Deleting unrealistic dprime estimate!!')
    lumdata(abs(lumdata)> 1e15) = nan;
end

% LUM
s2 = nanvar(lumdata);
sem = sqrt(s2./sum(~isnan(lumdata)));
mn = nanmean(lumdata);
figure; axes; hold on;
h = patch([TFBINCENTERS, fliplr(TFBINCENTERS)],[mn+sem, fliplr(mn-sem)],[.5 .5 .5],'Facealpha',.5);
plot(TFBINCENTERS, mn,'k-');
set(gca,'Xscale','log')
plot([1 60],[1.27 1.27]);
title([MONKEY,' LUM ',CELLTYPE]);

% RG
s2 = nanvar(rgdata);
sem = sqrt(s2./sum(~isnan(rgdata)));
mn = nanmean(rgdata);
Lnonnan = ~isnan(mn);
figure; axes; hold on;
h = patch([TFBINCENTERS(Lnonnan), fliplr(TFBINCENTERS(Lnonnan))],...
    [mn(Lnonnan)+sem(Lnonnan), fliplr(mn(Lnonnan)-sem(Lnonnan))],[.5 .5 .5],'Facealpha',.5);
plot(TFBINCENTERS, mn,'k-');
set(gca,'Xscale','log')
plot([1 60],[1.27 1.27]);
title([MONKEY,' RG ',CELLTYPE]);

%%
% Section 14
% Counting up numbers of spikes and d' in a counting window. testing the
% simple idea that only the first 'n' spikes get through from the LGN to
% the decision maker.
TFBINCENTERS = logspace(log10(1), log10(60), 15);
[filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',{CELLTYPE},'subjID',{MONKEY(1)});
nwindows = 8; % spike counting windows
signaloffsets = [zeros(nwindows,1) linspace(.04,.66,nwindows)']; % stretching window
data_dprime = [];
data_sp_per_s = [];
PLOTALL = 0;
for i = 1:length(filenames)
    stro = {};
    for j = 1:length(filenames{i})
        stro{j} = nex2stro(char(findfile(filenames{i}(j), fullfile(nexfilepath,'Greg',MONKEY))));
    end
    stro = strocat(stro);
    Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_l'));
    Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_m'));
    TF = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
    spikename = ['sig001',char(abs('a')+(spikecds(i)-1))];
    spiketimes = stro.ras(:,strcmp(stro.sum.rasterCells,spikename));

    pop_scale_fact = IsoSampGetPopulationScaleFactor(stro,ecc_to_diam_deg);
    dprimes = [];
    for j = 1:size(signaloffsets,1)
        [uniquestim, dprime] = IsoSampGetDPrime(stro, DPRIMEMETHOD, spikecds(i), signaloffsets(j,:));
        dprimes = [dprimes, dprime*pop_scale_fact];
    end
    Lblank = uniquestim(:,1) == 0 & uniquestim(:,2) == 0 & uniquestim(:,3) == 0;
    Llum = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & ~Lblank;
    Lrg = sign(uniquestim(:,1)) ~= sign(uniquestim(:,2));
    sp_per_s = zeros(size(uniquestim,1),size(signaloffsets,1));
    sp = zeros(size(uniquestim,1),size(signaloffsets,1));
    for j = 1:size(uniquestim,1)
        L = Lcc == uniquestim(j,1) & Mcc == uniquestim(j,2) & TF == uniquestim(j,3);
        spt = spiketimes(L);
        ston_t = stimon_t(L);
        spikevect = [];
        for k = 1:length(ston_t)
           spikevect = [spikevect; spt{k}-ston_t(k)];
        end
        for k = 1:size(signaloffsets,1)
            nspikes = sum(spikevect>signaloffsets(k,1) & spikevect < signaloffsets(k,2));
            sp(j,k) = nspikes./sum(L);
            sp_per_s(j,k) = sp(j,k)/diff(signaloffsets(k,:));
        end
    end
    if PLOTALL
        figure; axes;
        subplot(2,3,1);
        imagesc(sp(Llum,:));
        axis image ij;
        subplot(2,3,2);
        imagesc(sp_per_s(Llum,:));
        axis image ij;
        subplot(2,3,3);
        imagesc(dprimes(Llum,:));
        axis image ij;
        
        subplot(2,3,4);
        contour(1:size(signaloffsets,1),log10(uniquestim(Llum,3)),sp(Llum,:));
        axis ij;
        subplot(2,3,5);
        contour(1:size(signaloffsets,1),log10(uniquestim(Llum,3)),sp_per_s(Llum,:));
        axis ij;
        subplot(2,3,6);
        contour(1:size(signaloffsets,1),log10(uniquestim(Llum,3)),dprimes(Llum,:),1.27);
        axis ij;
        drawnow;
    end
    % taking TFs < 45 Hz (every cell contributes to every bin)
    LTF = uniquestim(:,3) < 45;
    tmp = [dprimes(Llum&LTF,:); nan*ones(size(data_dprime,1)-sum(Llum&LTF),size(data_dprime,2))];
    data_dprime = cat(3,data_dprime,tmp);
    tmp = [sp_per_s(Llum&LTF,:); nan*ones(size(data_dprime,1)-sum(Llum&LTF),size(sp_per_s,2))];
    data_sp_per_s = cat(3,data_sp_per_s,tmp);
end

figure; 
subplot(2,1,1); hold on;
imagesc(nanmean(data_dprime,3));
colormap(gray)
contour(nanmean(data_dprime,3),[1.27 1.27],'w-','LineWidth',2);
title([MONKEY,' ',CELLTYPE,' SNR']);
tmpTF = uniquestim(Llum,3);
set(gca,'Ytick',[1:2:size(uniquestim,1)],'Yticklabel',num2str(tmpTF(1:2:end),2));
set(gca,'Xtick',[1:2:size(signaloffsets,1)],'Xticklabel',num2str(signaloffsets(1:2:end,2),1));

subplot(2,1,2); hold on;
imagesc(nanmean(data_sp_per_s,3));
colormap(gray)
if strcmp(MONKEY,'Apollo')
    contour(nanmean(data_sp_per_s,3),[25 25],'w-','LineWidth',2);
else
    contour(nanmean(data_sp_per_s,3),[18 18],'w-','LineWidth',2);
end
title('Firing rate');
xlabel('Counting window (s)');
set(gca,'Xtick',[1:2:size(signaloffsets,1)],'Xticklabel',num2str(signaloffsets(1:2:end,2),1));
ylabel('TF (Hz)');
set(gca,'Ytick',[1:2:size(uniquestim,1)],'Yticklabel',num2str(tmpTF(1:2:end),2));

%%
% Section 15
% d' based on spike counts
% There's lots of information in the spike counts but slightly less than in
% the modulation.
[filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',{CELLTYPE},'subjID',{MONKEY(1)});
data = [];
for i = 1:length(filenames)
    counts = {};
    stro = {};
    for j = 1:length(filenames{i})
        stro{j} = nex2stro(char(findfile(filenames{i}(j), fullfile(nexfilepath,'Greg',MONKEY))));
    end
    stro = strocat(stro);
    Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_l'));
    Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_m'));
    TF = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
   
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
    spikename = ['sig001',char(abs('a')+(spikecds(i)-1))];
    spiketimes = stro.ras(:,strcmp(stro.sum.rasterCells,spikename));
    uniquestim = sortrows(unique([Lcc, Mcc, TF],'rows'),3);
    
    for j = 1:size(uniquestim,1)
        L = Lcc == uniquestim(j,1) & Mcc == uniquestim(j,2) & TF == uniquestim(j,3);
        st = stimon_t(L);
        sp = spiketimes(L);
        spikecounts = [];
        for k = 1:length(st)
            tmp = sp{k}-st(k);
            spikecounts = [spikecounts; sum(tmp>0 & tmp <.666)];
        end
        counts{j} = spikecounts;
    end
    %Lblank = Lcc == 0 & Mcc == 0 & TF == 0;
    tmpdata = zeros(size(uniquestim,1),1);
    for j = 1:size(uniquestim,1)
        sp1 = counts{1};
        sp2 = counts{j};
        tmpdata(j) = (mean(sp2)-mean(sp1))/sqrt(mean([var(sp1);var(sp2)]));
    end
    data = [data; uniquestim tmpdata i*ones(length(tmpdata),1)];
end

data(:,1) = sign(data(:,1));
data(:,2) = sign(data(:,2));

unique_conds = sortrows(unique(data(:,1:3),'rows'),3);

mn = nan*ones(size(unique_conds,1),1);
sem = nan*ones(size(unique_conds,1),1);
for i = 1:size(unique_conds,1)
    L = sign(data(:,1)) == unique_conds(i,1) & sign(data(:,2)) == unique_conds(i,2) & data(:,3) == unique_conds(i,3);
    mn(i) = nanmean(data(L,4));
    sem(i) = sqrt(nanvar(data(L,4))/sum(~isnan(data(L,4))));
end
Lblank = unique_conds(:,1) == 0 & unique_conds(:,2) == 0 & unique_conds(:,3) == 0;
Llum = sign(unique_conds(:,1)) == sign(unique_conds(:,2)) & ~Lblank;

figure; axes; hold on;
h = patch([unique_conds(Llum,3)', fliplr(unique_conds(Llum,3)')],[mn(Llum)+sem(Llum); flipud(mn(Llum)-sem(Llum))]',[.5 .5 .5],'Facealpha',.5);
plot(unique_conds(Llum,3),mn(Llum),'ko-','MarkerFaceColor','black');
set(gca,'Xscale','log');
title([MONKEY,' ',CELLTYPE]);
ylabel('SNR based on spike counts');
xlabel('Frequency (Hz)');

%%
% Section 16
% Comparing d' to non-parametric ROC between signal and noise distributions
% (Could convert ROC to equavalent d')
% "signal" and "noise" returned by IsoSampGetDprime are monotonically transformed 
% Mahalanobis distances (transform doesn't affect ROC).
[filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',{CELLTYPE},'subjID',{MONKEY(1)});
data = [];
for i = 1:length(filenames)
    stro = {};
    for j = 1:length(filenames{i})
        stro{j} = nex2stro(char(findfile(filenames{i}(j), fullfile(nexfilepath,'Greg',MONKEY))));
    end
    stro = strocat(stro);
    Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_l'));
    Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_m'));
    TF = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
    [uniquestim, dprime, signal, noise] = IsoSampGetDPrime(stro,DPRIMEMETHOD, spikecds(i));
    Lblank = uniquestim(:,1) == 0 & uniquestim(:,2) == 0 & uniquestim(:,3) == 0;
    Lrg = sign(uniquestim(:,1)) ~= sign(uniquestim(:,2));
    Llum = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & ~Lblank;
    tmp = [];
    for j = find(Llum)'
        n = noise{j}(~isnan(noise{j}));
        s = signal{j}(~isnan(signal{j}));
        tmp = [tmp; i uniquestim(j,3) dprime(j) roc(n,s)];
    end
    data =[data; tmp];
    corr(tmp(:,[3 4]),'type','Spearman')
end
figure; axes; hold on;
plot(data(:,3),data(:,4),'o');
xlabel('d prime'); ylabel('non-parametric auROC');
% Plot theoretical curve
mus = linspace(-1,5,100);
ROCareas = zeros(size(mus));
for i = 1:length(mus)
    x = linspace(-6,mus(i)+6,1000);
    ROCareas(i) = normcdf(x,0,1)*(normpdf(x,mus(i),1)./sum(normpdf(x,mus(i),1)))';
end
plot(mus,ROCareas,'k-');
set(gca,'Xlim',[-2 6],'Ylim',[0.2 1]);

%%
% Section 17
% Comparing parametric and nonparametric d's

Q = 6; % Q = 6, k = 9 is close to optimal. See next section.
k = 9;
MONKEYS={'Apollo','Utu'};
CELLTYPES = {'M','P'};
data = cell(2,2);
for monkey = MONKEYS
    for celltype = CELLTYPES
        [filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',celltype,'subjID',{monkey{1}(1)});
        tmp = [];
        for i = 1:length(filenames)
            stro = {};
            for j = 1:length(filenames{i})
                stro{j} = nex2stro(char(findfile(filenames{i}(j), fullfile(nexfilepath,'Greg',monkey))));
            end
            stro = strocat(stro);
            
            [uniquestim, dprime_parametric] = IsoSampGetDPrime(stro,1,spikecds(i));
            [~, dprime_nonparametric] = IsoSampGetKNNDPrime(stro,[Q k],spikecds(i));
            
            listsofar = data{strcmp(monkey,MONKEYS),strcmp(celltype,CELLTYPES)};
            listsofar{length(listsofar)+1} = [uniquestim dprime_parametric dprime_nonparametric];
            data{strcmp(monkey,MONKEYS),strcmp(celltype,CELLTYPES)}=listsofar;
            
        end
    end
end

% processing data into a file x TF x color x para/nonpara 4-D tensor, outdata
TFBINCENTERS = logspace(log10(1), log10(60), 15); % Native bins
binwidth = log10(TFBINCENTERS(2))-log10(TFBINCENTERS(1));
TFBINEDGES = 10.^(linspace(log10(TFBINCENTERS(1))-binwidth/2,log10(TFBINCENTERS(end))+binwidth/2,length(TFBINCENTERS)+1));
outdata = [];
for i = 1:2 % Monkey
    for j = 1:2 % Celltype
        nfiles = length(data{i,j});
        tmpdata = nan*ones(nfiles,length(TFBINCENTERS),2,2); % file, TF, color, para/nonpara
        for k = 1:nfiles
            uniquestim = data{i,j}{k}(:,1:3);
            tmp =  data{i,j}{k}(:,4:end);
            TFidxs = sum(repmat(uniquestim(:,3),1,length(TFBINEDGES))>repmat(TFBINEDGES,size(uniquestim,1),1),2);
            Lblank = uniquestim(:,1) == 0 & uniquestim(:,2) == 0 & uniquestim(:,3) == 0;
            Llum = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & ~Lblank;
            Lrg = sign(uniquestim(:,1)) ~= sign(uniquestim(:,2));
            for l = 1:size(uniquestim,1)
                if Lblank(l)
                    continue
                end
                for m = 1:2 % para/non-para
                    tmpdata(k,TFidxs(l),Llum(l)+2*Lrg(l),m) = tmp(l,m);
                end
            end
        end
        outdata{i,j} = tmpdata;
    end
end

figprefs;
colors = [0 0 0; 1 0 0];
titles = {'Magnocellular','Parvocellular'};
AXWIDTH = 4;
AXMARGIN = 1.5;
AXCOLS = 3.5+[0 AXWIDTH+AXMARGIN];
AXROWS = fliplr(2+([0:3]*(AXWIDTH+AXMARGIN)));
AXROWS(1:2) = AXROWS(1:2) + 1; % Extra margin
for i = 1:2 % Monkey
    for j = 1:2 % Celltype
        tmp = outdata{i,j};
        hax1 = axes('position',[AXCOLS(j) AXROWS(i*2-1) AXWIDTH AXWIDTH]); hold on;
        hax2 = axes('position',[AXCOLS(j) AXROWS(i*2) AXWIDTH AXWIDTH]); hold on;
        for k = 1:2 % lum/color
            for m = 1:2 % parametric/nonparametric
                axes(hax1);
                if m == 2
                    plot([TFBINCENTERS(1) TFBINCENTERS(end)],[0 0],'k:');
                end

                n = sum(~isnan(tmp(:,:,k,m)));
                L = n >= 2; % 2 data points minimum for plotting
                mn(m,:) = nanmean(tmp(:,:,k,m));
                sd(m,:) = nanstd(tmp(:,:,k,m));
                sem(m,:) = sd(m,:)./sqrt(n);
                
                h_patch = patch([TFBINCENTERS(L), fliplr(TFBINCENTERS(L))],[mn(m,L)+sem(m,L), fliplr(mn(m,L)-sem(m,L))],colors(k,:),'Facealpha',.25,'LineStyle','none');
                h_line = plot(TFBINCENTERS(L),nanmean(tmp(:,L,k,m)),'ko-','Color',colors(k,:),'MarkerFaceColor',colors(k,:),'MarkerSize',4);
                
                if m == 2 % nonparametric
                    set(h_line,'LineStyle','--','MarkerFaceColor','none');
                end
            end
            % Trajectory plot
            axes(hax2);
            plot([-2 4],[-2 4],'k:');
            h = plot(mn(1,L),mn(2,L),'k-','Color',colors(k,:),'linewidth',1);
            for l = 1:sum(L)
                plot(mn(1,l)+sd(1,l)*[-1 1],mn(2,l)*[1 1],'-','color',colors(k,:),'linewidth',0.5);
                plot(mn(1,l)*[1 1],mn(2,l)+sd(2,l)*[-1 1],'-','color',colors(k,:),'linewidth',0.5);
               % h = plot(mn(1,l),mn(2,l),'o','MarkerSize',l/2+2,'MarkerEdgeColor','white','MarkerFaceColor',colors(k,:),'LineWidth',.1);
               h = plot(mn(1,l),mn(2,l),'o','MarkerSize',4,'MarkerEdgeColor','white','MarkerFaceColor',colors(k,:),'LineWidth',.1);
            end
        end
        axes(hax1);
        set(gca,'Xscale','log','Ylim',[-1 4],'Xlim',[TFBINEDGES(1) TFBINEDGES(end)],'Xtick',[1 10],'XtickLabel',{'1','10'});
        title([titles{j},': monkey ',num2str(i)]);
        ylabel('Signal-to-noise ratio (d'')');
        xlabel('frequency (Hz)');
        
        axes(hax2);
        xlabel('parametric d''');
        ylabel('non-parametric d''');
        set(gca,'Xlim',[-2 4],'Ylim',[-2 4]);
    end
end

%%
% Section 18 
% Finding the values of k and Q that work best for the KNN sensitivity
% analysis.
% Restricting attention to magnocellular neurons

[filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',{'M'});
stros = {};
nblanks = [];
for i = 1:length(filenames)
    stro = {};
    for j = 1:length(filenames{i})
        stro{j} = nex2stro(char(findfile(filenames{i}(j), fullfile(nexfilepath,'Greg'))));
    end
    stros{i} = strocat(stro);
    Lcc = stros{i}.trial(:,strcmp(stros{i}.sum.trialFields(1,:),'stim_l'));
    Mcc = stros{i}.trial(:,strcmp(stros{i}.sum.trialFields(1,:),'stim_m'));
    nblanks(i) = sum(Lcc == 0 & Mcc == 0);
end

Qs = linspace(4,9,10);
ks = [3:2:min(nblanks) nan]; % nan means use default k
data = [];
for Qidx = 1:length(Qs)
    for kidx = 1:length(ks)
        [Qidx kidx]
        Q = Qs(Qidx);
        k = ks(kidx);
        if isnan(k)
            k = [];
        end
        tmp = nan*ones(length(stros),1);
        for i = 1:length(stros)
            [uniquestim, dprime_nonparametric] = IsoSampGetKNNDPrime(stros{i},[Q k],spikecds(i));
            Lblank = all(uniquestim == 0,2);
            tmp(i) = nanmean(dprime_nonparametric(~Lblank));
        end
        data(Qidx,kidx) = nanmean(dprime_nonparametric(~Lblank));
    end
end
figure;
axes;
if size(data,1) == 1 | size(data,2) == 1
    plot(Qs,data);
else
    imagesc(data);
    xlabel('k'); ylabel('Q');
    set(gca,'Xtick',1:length(ks),'XTickLabel',ks)
    set(gca,'Ytick',1:length(Qs),'YTickLabel',Qs)
end