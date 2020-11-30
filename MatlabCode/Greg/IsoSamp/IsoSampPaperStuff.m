% IsoSampPaperStuff
%
% Section 0) Constants
%
% Section 1) Data justification for formulas for magno and parvo RF sizes
% Largely taken from RGC_RF_sizes.m
%
% Section 2) t-SNE and STAs from whitenoise. Take largely from WNPop.m
% section 13
%
% Section 3) Average single neuron d' (and extrapolated population d') as a
% function of retinal eccentricity. Also as a function of temporal frequency.
% See IsoSampPop.m sections 2 and 4.
%
% Section 4) Example rasters and d-prime curves
%
% Section 5) Input efficiency (one plots/monkey)
%
% Section 6) Population d-prime with cone ideal observer(s)
%     6.1) numbers of neurons in pool
%     6.2) distribution of RF sizes
%     6.3) Proportion of SNR loss as a function of TF
%
% Section 7) Comparison of d-prime from ON and OFF cells.
%
% Section 8) Microsaccade-triggered PSTHs.
%
% Section 9) PSTHs of microsaccade occurrences.
%
% Section 10) Noise spectra + ISI distributions
%
% Section 11) Psychophysical temporal integration
%
% Section 12) Sliding window d' & firing rate analysis.
%
% Section 13) Looking at distribution of standard deviation estimates of
% the "normalized" Mahalanobis distance.
%
% Section 14) Stats showing that d' increases for M-neurons from 1 to 8 Hz.
%
% Section 15) Comparing d' from non-parametric auROC (on Mahalanobis
% distances)
%
% Section 16) Power spectrum of the nominally 1 Hz stimulus
%
% Section 17) Collecting stro structures for a GitHub repository
%
% Section 18) Trying to figure out how large Jiang et al.'s Gaussian
% stimuli were.
%
% Section 19) Trying to convert a 14% contrast threshold prediction error
% rate to d-primes
%
% -----------------------------------------------------------------------
% Material below here is specifically for the color paper. 
% Section 20) Cone weights for M, P, K.
%
% Section 21) Cone weights from white noise predict sensitivity to L+M and
% L-M modulation in IsoSamp.
%
% Section 22) Analysis of trajectories for different cell types in the 
% L+M/L-M sensitivity plane as TF changes.
%
% Section 23) Power spectra from L+M and L-M (making the point that, on
% average, magnocellular neurons have an increased response to L-M but 
% relatively little entrained modulation.
%
% Section 24) Proportion of SNR loss across stages. Identical to
% section 6.3 but for L-M modulations instead of L+M.
%
% Section 25) d' as a function of TF and counting window (small 2-D space
% of start and end times).
% 
%%
% Section 0)
% Constants that should be shared across cells in this script so that I
% don't accidentally switch formulas here and there.
MMPERDEG = 0.223; % mm/deg (Perry and Cowey 1985)
DEGPERMM = 1/MMPERDEG; % deg/mm
DPRIMEMETHOD = 1;
MONKEYS={'Apollo','Utu'};
CELLTYPES = {'M','P'};
TEMPORONASALSCALEFACTOR = .8;
% DEBUGGING
%TEMPORONASALSCALEFACTOR = 1;

ONOFFCORRELATION = 0.05;

RFTRUNCATIONINSD = 2;
HUMAN2MONKPSCALEFACTOR = .80; % From Dacey and Petersen. reasonable range: [.77 .81];

%TFBINEDGES = logspace(0,log10(60),8);
%TFBINCENTERS = sqrt(TFBINEDGES(1:end-1).*TFBINEDGES(2:end)); % geomean

TFBINCENTERS = logspace(log10(1), log10(60), 15); % Native bins
binwidth = log10(TFBINCENTERS(2))-log10(TFBINCENTERS(1));
TFBINEDGES = 10.^(linspace(log10(TFBINCENTERS(1))-binwidth/2,log10(TFBINCENTERS(end))+binwidth/2,length(TFBINCENTERS)+1));
MONKEY2MARKER = '^';
EXAMPLEMAGNOCELL = 'A061517005.nex'; % example magnocell
EXAMPLEPARVOCELL = 'A061617006.nex'; % example parvocell
EXAMPLEPARVOCELL = 'A061517002.nex'; % example parvocell

GETRFFROMSTAMODE = 2;
GETRFFROMSTATHRESH = .95;
MAXT= 6; % n frames back to compute STA

ecc_to_diam_deg_M = @(rf_r_deg) 10.^(-1.2459+0.0345*rf_r_deg); % temporal retina equivalent
a = 0.9729; % Table 1
r2 = 1.084; % Table 1
re = 7.633; % Table 1
dc_0 = 14804.6; % Cone density of fovea
rm = 41.03; % See Equation 7
ecc_to_diam_deg_P = @(x)(sqrt(2./(sqrt(3).*... % Equation 9. Distance between adjacent midget RF centers.
    (2*dc_0.*(1+x./rm).^-1.*(a*(1+(x./r2)).^-2+(1-a)*exp(-x./re)))...
    ./2))... % Dividing y by 2 to estimate RF size from only ON or OFF mosaics (halving the density).
    *HUMAN2MONKPSCALEFACTOR); % Monkey midget RFs are slightly smaller than human midget RFs

bpdf_vec=@(x,y,mu_x,mu_y,sigma)(exp(-((x-mu_x).^2)./(2*sigma^2)-((y-mu_y).^2)/(2*sigma^2))./(2*pi*sigma^2)); % bivariate normpdf

% Stuff for DTcones_gh
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
%%
% Section 1) Magnocellular and parvocellular receptive field sizes

DL_magno_nasal_data = xlsread('RGC_RFsizes_vs_eccentricity','DL1984MagnoNasalData'); % x = eccentricity (deg), y = RF field radius (deg)
DL_magno_temporal_data = xlsread('RGC_RFsizes_vs_eccentricity','DL1984MagnoTemporalData'); % x = eccentricity (deg), y = RF field radius (deg)

% "*0.61" to convert nasal RFs to equivalent temporal eccentricity
% "*2" to convert radius to diameter
% "/sqrt(2)" to go out only 1 SD instead of 1/e

rf_r_deg = logspace(log10(2),log10(18),20);
hobj = []; % handles for legend

figprefs; axes;
hobj(1) = plot(DL_magno_nasal_data(:,1)*.61,DL_magno_nasal_data(:,2)*2/sqrt(2),'ko');
hobj(2) = plot(DL_magno_temporal_data(:,1),DL_magno_temporal_data(:,2)*2/sqrt(2),'k+');
set(gca,'Yscale','log');
y = log10([DL_magno_nasal_data(:,2); DL_magno_temporal_data(:,2)].*2./sqrt(2));
X = [ones(length(DL_magno_nasal_data)+length(DL_magno_temporal_data),1) [DL_magno_nasal_data(:,1)*.61; DL_magno_temporal_data(:,1)]];
b_equiv_temporal = regress(y,X);
% Error checking
if b_equiv_temporal(1) - -1.2459 > 10^-4 || b_equiv_temporal(2)-.0345 > 10^-4 
    error('regression coefficients do not match');
end

hobj(3) = plot(rf_r_deg,10.^([ones(length(rf_r_deg),1) rf_r_deg']*b_equiv_temporal),'k-');
title('Derrington and Lennie 1984 magnocellular');
xlabel('Temporal equivalent eccentricity (deg)');
ylabel('Center diameter (2 SD in deg)');
set(gca,'Yscale','log');

%  midget RGC RF sizes from Watson model and cone RF sizes from Goodchild
%  et al. 1996 model

% First Watson (using data from nasal *visual field* (= temporal retina)).
% Nasal retina has high cone/RGC density.
a = 0.9729; % Table 1
r2 = 1.084; % Table 1
re = 7.633; % Table 1
dc_0 = 14804.6; % Cone density of fovea
rm = 41.03; % See Equation 7
y = 2*dc_0.*(1+rf_r_deg./rm).^-1.*(a*(1+(rf_r_deg./r2)).^-2+(1-a)*exp(-rf_r_deg./re)); % From Equation 8

% trying to back out RF size from RF density.
% y is in RF density = RFs/deg^2
% sqrt(1./y(1))
% (pi*sqrt(3))/6) above to account for the fact that packed circles.
% have slightly lower density than packed squares and therefore RFs
% should be slightly smaller than 1/sqrt(density).
% rfsize = (pi*sqrt(3))/6.*sqrt(1./y);

% Using Equation 9
% Dividing y by 2 to get density of ON (or OFF) cells only
rfsize = sqrt(2./(sqrt(3).*y./2))*HUMAN2MONKPSCALEFACTOR; % Equation 9. Distance between adjacent midget RF centers. 
if ~all(rfsize == ecc_to_diam_deg_P(rf_r_deg))
    error('Disagreement between nominally identical ways of calculating P-cell RF size');
end

hobj(4) = plot(rf_r_deg,rfsize,'k--');

% Now Goodchild et al. 1996
% Also using temporal retina data
temporalcoeffs=[150.9676 -1.2220 35.9979 -0.1567 9.9936 -0.0258]; % (Goodchild et al., 1996 Table 2)
conedensfun = @(coeffs,x)(coeffs(1).*(exp(coeffs(2).*x)))+...
    (coeffs(3).*(exp(coeffs(4).*x)))+...
    (coeffs(5).*(exp(coeffs(6).*x)));
coneDensity = conedensfun(temporalcoeffs,rf_r_deg);
conesPerMM2=coneDensity*1e3;
conesPerDeg = conesPerMM2*MMPERDEG.^2;

% Using same conversion of density to hexagonal space from Watson 2014
% Equation A4. Measuring the cone "RF" as the intercone spacing times 0.58
% to make it comparable to the midget RF size measurements.

rfsize = sqrt(2./(sqrt(3)*conedensfun(temporalcoeffs,rf_r_deg)*1e3*MMPERDEG.^2))*HUMAN2MONKPSCALEFACTOR;
rfsize = sqrt(2./(sqrt(3)*conedensfun(temporalcoeffs,rf_r_deg)*1e3*MMPERDEG.^2));
hobj(5) = plot(rf_r_deg,rfsize,'k:');
labels = {'Magno nasal (Derrington and Lennie)','Magno temporal (Derrington and Lennie)','Magno fit (Derrington and Lennie)',...
    'Midgets (Watson)','Cones (Goodchild et al.)'};
legend(hobj,labels,'Location','northwest')
set(gca,'Xlim',[0 18],'Ylim',[.01 .5],'TickDir','out');

% Agrees reasonably well with Wool et al. 2018 who say that by ~9� the
% average midget has two cones in the RF center.

% Just some quick calculations to figure out how much of the probability
% mass of a bivariate Gaussian is within 2 SDs (or within 2/.58 = 3.44 SDs)
x = linspace(-6,6,1000)';
y = normpdf(x,0,1)*normpdf(x,0,1)';
y = y./sum(y(:));

[xx,yy] = meshgrid(x,x);
mask = (xx.^2+yy.^2) < 1; % going out 1 SD in all directions (=2 SD diameter)
sum(sum(y.*mask))

mask = (xx.^2+yy.^2) < 4; % 4 = 2^2
sum(sum(y.*mask))

mask = (xx.^2+yy.^2) < (2/HUMAN2MONKPSCALEFACTOR)^2; % 3.4^2 = 11.89
sum(sum(y.*mask))

y = normpdf(x,0,1);
y = y./sum(y);
sum(y(abs(x)<1))

%%
% Section 2
% t-SNE plot showing clustering of STAs (just to show the different
% physiological cells types in the data set are the ones people expect to
% see and that the L-ON and M-ON parocellular neurons are a little harder
% to discriminate than other neuronal types.)

conn = database('Nex_Paradigm_Sort','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
whitenoisedata = fetch(conn, 'SELECT fileID, rfX, rfY, neuron, spikeCode FROM WhiteNoiseLGN_forIS WHERE quality = ''1''');
isosampdata = fetch(conn, 'SELECT fileID, rfX, rfY, neuron, cellClass FROM IsoSamp_LGN WHERE quality = ''1''');
close(conn);

uniqueWNneurons = unique([whitenoisedata{:,4}]);
WNfilenames = {};
for neuronid = uniqueWNneurons
    LWN =  [whitenoisedata{:,4}] == neuronid;
    Lisosamp = [isosampdata{:,4}] == neuronid;
    if any(Lisosamp)
        tmp = [];
        for i = find(LWN)    
        	tmp = [tmp, whitenoisedata(i,1)];
        end
        WNfilenames{length(WNfilenames)+1} = tmp;
        
        % Sanity checks
        tmp=[];
        for i = find(Lisosamp)
            tmp=[tmp;isosampdata{i,2:4}];
        end
        for i = find(LWN)
            tmp=[tmp;whitenoisedata{i,2:4}];
        end
        if any(std(tmp))
            error('White noise and IsoSamp entries into database disagree');
        end
    end
end

data=[];
filenamestems = {};
Ms = [];
for i = 1:length(WNfilenames)
    stro=[];
    for j = 1:length(WNfilenames{i})
        stro = strocat(stro,nex2stro(findfile(WNfilenames{i}{j})));
    end
    spikename = whitenoisedata{strcmp(whitenoisedata(:,1),WNfilenames{i}(1)),5};

    nstixperside = stro.sum.exptParams.nstixperside;
    ntrials = length(stro.sum.absTrialNum);
    spikeidx = find(strcmp(stro.sum.rasterCells(1,:),spikename));
    
    % Reconstructing the M matrix and gamma table
    funds = stro.sum.exptParams.fundamentals;
    funds = reshape(funds,[length(funds)/3,3]);
    mon_spd = stro.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    Ms(:,:,i) = funds'*mon_spd;

    out = getWhtnsStats(stro,MAXT,'STCOVmex', {nstixperside^2, 3, MAXT}, spikename);
    STA = out{1}/out{3};
    inRF = getRFfromSTA(STA,GETRFFROMSTAMODE,GETRFFROMSTATHRESH);
    STA = reshape(STA,[nstixperside.^2  3 MAXT]);
    temporalSTA = squeeze(STA(logical(inRF(:)),:,:));
    drawnow
    if ~ismatrix(temporalSTA)
        temporalSTA = mean(temporalSTA,1);
    end
    %temporalSTA = squeeze(STA(whichpix,:,:));
    
    data(:,:,i) = squeeze(temporalSTA);
    filenames{i} =  WNfilenames{i}{1};
end

% --------------
% Setting up clusters by hand and seeing where they fall in the t-SNE plot
% 1) 1 = P, 2 = M, 3 = K
% 2) 0 = OFF, 1 = ON
% 3) 0 = L, 1 = M
% --------------
svdvars = zeros(size(data,3),2);
celltypemat = nan(size(data,3),3);
for j = 1:size(data,3)
    sta = squeeze(data(:,:,j));
    [u,s,v] = svd(sta);
    svdvars(j,:) = [s(1) sum(diag(s))]; 
    L = strcmp(filenames{j},whitenoisedata(:,1));
    neuronidx = whitenoisedata{L,4};
    celltype = isosampdata{[isosampdata{:,4}] == neuronidx,5};
    switch(celltype)
        case {'P'}
            celltypemat(j,1) = 1;
        case {'M'}
            celltypemat(j,1) = 2;
        case {'K'}
            celltypemat(j,1) = 3;
        otherwise
            celltypemat(j,1) = 4;
    end
    rgb = u(:,1);
    temporal = v(:,1);
    if max(abs(rgb)) == -min(rgb)
        rgb = -rgb;
        temporal = -temporal;
    end
    if max(abs(temporal)) == max(temporal)
        celltypemat(j,2) = 1; % ON
    else
        celltypemat(j,2) = 0; % OFF
    end
    M = Ms(:,:,j);
    lm_ratio = sta(1,2)/sta(2,2);
    if abs(lm_ratio-M(2,1)/M(2,2)) < abs(lm_ratio-M(1,1)/M(1,2))
        celltypemat(j,3) = 1; % M
    else
        celltypemat(j,3) = 0; % L
    end
end

hist(svdvars(:,1)./svdvars(:,2))
mean(svdvars(:,1)./svdvars(:,2))
std(svdvars(:,1)./svdvars(:,2))

clusterids = {'Parvo L-cone ON','Parvo M-cone ON','Parvo L-cone OFF','Parvo M-cone OFF','Magno ON','Magno OFF','K'};
colors = [1 0 0; 0 .5 0; 0 1 1; 1 0 1; .75 .75 0; 0 0 0; 0 0 1];

idx = zeros(size(celltypemat,1),1);
nclusters = length(clusterids);
for i = 1:nclusters
    done = 0;
    L = zeros(size(celltypemat,1),1);
    if strfind(clusterids{i},'Parvo')
        L = L | celltypemat(:,1) == 1;
    elseif strfind(clusterids{i},'Magno')
        L = L | celltypemat(:,1) == 2;
    elseif strfind(clusterids{i},'K')
        L = L | celltypemat(:,1) == 3;
        done = 1;
    end
    if ~done
        if strfind(clusterids{i},'ON')
            L = L & celltypemat(:,2) == 1;
        else
            L = L & celltypemat(:,2) == 0;
        end
        if strfind(clusterids{i},'Parvo')
            if strfind(clusterids{i},'L-cone')
                L = L & celltypemat(:,3) == 0;
            elseif strfind(clusterids{i},'M-cone')
                L = L & celltypemat(:,3) == 1;
            end
        end
    end
    idx(logical(L)) = i;
end

% Reducing the dimensionality of each STA via SVD
newdata = zeros(size(data,1)+size(data,2),size(data,3));
for i = 1:size(data,3)
    sta = squeeze(data(:,:,i));
    [u,s,v] = svd(sta);
    if sta(:,2)'*u(:,1) < 0
        u = -u;
        v = -v;
    end
    newdata(:,i) = [u(:,1);v(:,1)];
end
y = tsne(newdata','Algorithm','exact','Distance','euclidean','LearnRate',1,'Perplexity',5,'Exaggeration',1);

figprefs; axes('Position',[1 1 9 9]);
for i = 1:nclusters
    for j = find(idx == i)'
        h = plot(y(j,1),y(j,2),'ko','MarkerFaceColor',colors(i,:));
    end
end
axis equal;

% STAs for each cluster
figprefs;
framestep = stro.sum.exptParams.nrepframes/stro.sum.exptParams.framerate*1000; % ms/frame
t = [0:size(data,2)-1]*framestep;
for i = 1:nclusters
    subplot(ceil(sqrt(nclusters)),ceil(sqrt(nclusters)),i); hold on;
    set(gca,'ColorOrder',[1 0 0; 0 1 0; 0 0 1],'Ycolor',colors(i,:),'Xcolor',colors(i,:));
    for j = find(idx == i)'
        plot(t,squeeze(data(:,:,j))','LineWidth',2);
    end
    set(gca,'Xlim',[t(1) t(end)],'Ylim',[-.045 .045],'Ytick',[],'TickDir','out');
end
xlabel('Time (ms)');
ylabel('STA (A.U.)');

%%
% Section 3
% dprime (single neuron and population) as a function of eccentricity and
% temporal frequency. LUM only.
MARKERSIZESCALEFACTOR = 40; % so points are visible

data = [];
dprime_v_TF = [];
ns = [];
for MONKEY = MONKEYS
    for CELLTYPE = CELLTYPES
        [filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',CELLTYPE,'subjID',{MONKEY{1}(1)});
        for i = 1:length(filenames)
            stro = {};
            for j = 1:length(filenames{i})
                stro{j} = nex2stro(char(findfile(filenames{i}(j), fullfile(nexfilepath,'Greg',MONKEY))));
            end
            stro = strocat(stro);
            examplecell = 0;
            if any(strcmp(EXAMPLEPARVOCELL,filenames{i})) | any(strcmp(EXAMPLEMAGNOCELL,filenames{i})) 
                examplecell = 1;
            end
            if strcmp(CELLTYPE,'M')
                ecc_to_diam_deg = ecc_to_diam_deg_M;
            elseif strcmp(CELLTYPE,'P')
                ecc_to_diam_deg = ecc_to_diam_deg_P;
            else
                error('Unknown celltype');
            end
            population_scalefactor = IsoSampGetPopulationScaleFactor(stro, ecc_to_diam_deg, TEMPORONASALSCALEFACTOR, RFTRUNCATIONINSD, ONOFFCORRELATION, 2);
            rfx = stro.sum.exptParams.rf_x/10;
            rfy = stro.sum.exptParams.rf_y/10;

            [uniquestim, dprime] = IsoSampGetDPrime(stro,DPRIMEMETHOD,spikecds(i));
            Lblank = uniquestim(:,1) == 0 & uniquestim(:,2) == 0 & uniquestim(:,3) == 0;
            Llum = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & ~Lblank;
      
            tmp = nan*ones(2,length(TFBINCENTERS));
            for j = 1:length(TFBINCENTERS)
                Ltf = uniquestim(:,3) > TFBINEDGES(j) & uniquestim(:,3) <= TFBINEDGES(j+1);
                if sum(Ltf&Llum) > 0
                    tmp(1,j) = mean(dprime(Llum&Ltf));
                    tmp(2,j) = 1; % each cell counts as an independent entity
                end
            end
            
            dprime_v_TF = [dprime_v_TF; tmp(1,:)];
            ns = [ns; tmp(2,:)];
            data = [data; find(strcmp(MONKEY,MONKEYS)) find(strcmp(CELLTYPE,CELLTYPES)) mean(dprime(Llum)) population_scalefactor rfx rfy examplecell];     
        end
    end
end

% Plotting the data
xextent = [-14 0];
yextent = [-10 10];
aspect_ratio = (yextent(2)-yextent(1))./(xextent(2)-xextent(1));
colors = [1 1 1; .5 .5 .5];
hax = [];
figprefs;
hax(1) = axes('Position',[3 14 7 7*aspect_ratio]);
hax(2) = axes('Position',[3 2 7 7*aspect_ratio]);
hax(3) = axes('Position',[13 13 7 7]);
hax(4) = axes('Position',[13 2 7 7]);

for monkey_idx = unique(data(:,1))'
    for celltype_idx = unique(data(:,2))'
        L = data(:,1) == monkey_idx & data(:,2) == celltype_idx;
        axes(hax(celltype_idx)); hold on;
        for i = find(L)'
            h = plot(data(i,5),data(i,6),'ko','MarkerFaceColor','none','LineWidth',1);
            set(h,'MarkerSize',max(.1,data(i,3)./max(data(:,3))*MARKERSIZESCALEFACTOR));
            if monkey_idx ==2
                set(h,'Marker',MONKEY2MARKER);
                if MONKEY2MARKER ~= '^'
                    error('Monkey 2 marker has to be a triangle otherwise MarkerSize scaling does not work');
                else
                    set(h,'MarkerSize',get(h,'MarkerSize')*1.15);
                end
            end
            if data(i,7) % Example cell
                set(h,'LineWidth',2);
            end
        end
    end
end

for i = 1:2
    axes(hax(i));
    set(gca,'Ylim',yextent,'Xlim',xextent);
    set(gca,'YTick',[yextent(1):5:yextent(end)]);
    plot(0,0,'r+','MarkerSize',10);
end

% Scale bar. d' = 1
axes(hax(1));
set(gca,'Units','points');
tmp = get(gca,'Position');
ax_width_points = tmp(3);
tmp = get(gca,'Xlim');
ax_width_data = tmp(2)-tmp(1);
plot([0 MARKERSIZESCALEFACTOR]*(ax_width_data/ax_width_points)-13,[-8.5 -8.5],'k-','LineWidth',2);
set(gca,'Units','centimeters');

% Scatterplots of d-prime vs eccentricity
for i = 1:2
    axes(hax(i+2)); hold on;
    set(gca,'Xlim',[2 14]);
    xlabel('Eccentricity (�)');
    for celltype_idx = unique(data(:,2))'
        for monkey_idx = unique(data(:,1))'
            L = data(:,1) == monkey_idx & data(:,2) == celltype_idx;
            if i == 1
                y = data(L,3); % single neuron d-prime
                str = ' single neuron ';
            else
                y = data(L,4).*data(L,3); % Population d-prime
                str = ' neuronal population ';
            end
            h = plot(sqrt(data(L,5).^2+data(L,6).^2),y,'ko','MarkerFaceColor',colors(celltype_idx,:));
            if monkey_idx ==2
                set(h,'Marker',MONKEY2MARKER);
            end
            [b,~,~,~,stats] = regress(y,[ones(sum(L),1) sqrt(data(L,5).^2+data(L,6).^2)]);
            disp([MONKEYS{monkey_idx},': ',CELLTYPES{celltype_idx},': p = ',num2str(stats(3)),': b = ',num2str(b(2)),str]);
        end
        L = data(:,2) == celltype_idx;
        if i == 1
            y = data(L,3); % single neuron d-prime
            str = ' single neuron ';
        else
            y = data(L,4).*data(L,3); % Population d-prime
            str = ' neuronal population ';
        end
        [b,~,~,~,stats] = regress(y,[ones(sum(L),1) sqrt(data(L,5).^2+data(L,6).^2)]);
        disp([CELLTYPES{celltype_idx},': p = ',num2str(stats(3)),': b = ',num2str(b(2)),str]);
        plot([2 14],[[1 2]*b [1 14]*b],'k-');
    end
    if i == 1
        ylabel('Sensitivity (d'')');
        set(gca,'Ytick',[0 1 2],'Xtick',[2 6 10 14])
    else
        ylabel('Population sensitivity (d'')');        
        set(gca,'Ytick',-1:2:8,'YtickLabel',[0:2:8],'Xtick',[2 6 10 14])
    end
end
set(gcf,'Render','painters');

% More statistical tests:
% (1) Difference in slopes of magno and parvo cells' SNR with eccentricity
% (single neuron SNR and neuron population SNR)
y = data(:,3); % single neuron d-prime
I_celltype = data(:,2)-1;
ecc = sqrt(data(:,5).^2+data(:,6).^2);
[b,bint,~,~,stats] = regress(y,[ones(length(y),1) I_celltype ecc ecc.*I_celltype]);
bint(end,:)

% Plot of dprime vs TF for magno and parvo (both monkeys).
% Individual neuron d's averaged. Not using popn scalefactor.
figprefs;
hax = [];
hax(1) = axes('Position',[4 14 7 7],'Tickdir','out'); hold on;
hax(2) = axes('Position',[4 3 7 7],'Tickdir','out'); hold on;
for monkey_idx = 1:2
    for celltype_idx = 1:2
        L = data(:,1) == monkey_idx & data(:,2) == celltype_idx;
        popn_scale_fact = data(L,4);
        %mn = nanmean(dprime_v_TF(L,:).*popn_scale_fact);
        %sd = nanstd(dprime_v_TF(L,:).*popn_scale_fact);

        mn = nanmean(dprime_v_TF(L,:));
        sd = nanstd(dprime_v_TF(L,:));
        sem = sd./sqrt(nansum(ns(L,:)));
        axes(hax(monkey_idx));
        Lbins = ~isnan(mn);
        patch([TFBINCENTERS(Lbins), fliplr(TFBINCENTERS(Lbins))],[mn(Lbins)+sem(Lbins), fliplr(mn(Lbins)-sem(Lbins))],[.5 .5 .5],'Facealpha',.5);
        h = plot(TFBINCENTERS(Lbins), mn(Lbins),'k-o','LineWidth',1);
        if celltype_idx == 2
            set(h,'MarkerEdgeColor','black','MarkerFaceColor',[.5 .5 .5]);
        else
            set(h,'MarkerEdgeColor','black','MarkerFaceColor','white','MarkerSize',6);
        end
    end
    plot([TFBINEDGES(1) TFBINEDGES(end)],[1.27 1.27],'k-')
    set(gca,'Xscale','log','Xlim',[TFBINEDGES(1) TFBINEDGES(end)],'Ylim',[-.25 4],'FontSize',12);
    xlabel('Frequency (Hz)','FontSize',15);
    ylabel('SNR','FontSize',15);
end
% Adding a second Y-axis showing contrast.
for monkey_idx = 1:2
    MONKEY = MONKEYS{monkey_idx};
    [filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','subjID',{MONKEY(1)});
    stro = {};
    for j = 1:length(filenames{end})
        stro{j} = nex2stro(char(findfile(filenames{end}(j), fullfile(nexfilepath,'Greg',MONKEY))));
    end
    stro = strocat(stro);
    modelparams = LMTF_global_to_local_model(stro.sum.exptParams.modelparams, 5, 0, 5); % 5 degrees out on H meridian. mode 5.
    TF = logspace(log10(TFBINCENTERS(1)), log10(TFBINCENTERS(end)),100);
    r = LMTF_thresh_from_model(pi/4*ones(100,1), TF',modelparams);
    % Reconstructing the M matrix and gamma table
    %fundamentals = stro.sum.exptParams.fundamentals;
    %fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    %mon_spd = stro.sum.exptParams.mon_spd;
    %mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    %mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    %M = fundamentals'*mon_spd;
    
    % Getting the background rgb/lms
    % Gamma table is basically a straight line. I dont want to add noise so
    % just using nominal "voltage" values for "intensity".
    % I cold probably just use Lcone contrast and get the same answer since
    % I'm normalizing at the end.
    %bkgndRGB = stro.sum.exptParams.bkgndrgb;
    %bkgndlms = M*bkgndRGB;
  
    %stimlms = repmat([bkgndlms(1) bkgndlms(2)],length(r),1).*repmat((1+r),1,2);
    %bkgndlum = [bkgndlms(1) bkgndlms(2)]*[1.8 1]';
    %contrast = (stimlms*[1.8 1]'-bkgndlum)./bkgndlum;
    %contrast-r % that was a lot of work for nothing. Contrast = r.
    axes(hax(monkey_idx));
    yyaxis right
    h = plot(TF,r,'b:','Linewidth',2,'Color',[.5 .5 .5]);
    set(gca,'Yscale','log','Ylim',[0.04 1],'YtickLabel',{'0.1', '1'},'Xticklabel',{'1','10'},'Ycolor',[.5 .5 .5]);
    ylabel('Contrast');
    title(MONKEY);
end
set(gcf,'Render','painters');

%%
% Section 4
% Rasters of example neurons. Taken from Section 4 of "IsoSampFigures.m"

filename = EXAMPLEMAGNOCELL;
STIMTYPE = 'RG';
if strcmp(STIMTYPE,'RG')
    binwidth = 0.01; % bins for PSTH
else
    binwidth = 0.008; % bins for PSTH
end

spikeNum = 'sig001a';
if iscell(filename)
    stros = {};
    for i = 1:length(filename)
        stros{i} = nex2stro(findfile(filename{i},[nexfilepath,filesep,'Greg',filesep,'Apollo']));
    end
    stro = strocat(stros);
else
    stro = nex2stro(findfile(filename,[nexfilepath,filesep,'Greg',filesep,'Apollo']));
end

Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_l'));
Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_m'));
TF = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'rew_t'));

dur = mean(stimoff_t-stimon_t);
spikeidx = strcmp(stro.sum.rasterCells(1,:),spikeNum);
spikes = stro.ras(:,spikeidx);
[uniquestim, dprime, signal, noise] = IsoSampGetDPrime(stro,1);

Lblank = all(uniquestim == 0,2);
if strcmp(STIMTYPE, 'LUM')
    Lstimtype = (sign(uniquestim(:,1)) == sign(uniquestim(:,2))) & ~Lblank;
elseif strcmp(STIMTYPE, 'RG')
    Lstimtype = (sign(uniquestim(:,1)) ~= sign(uniquestim(:,2))) & ~Lblank;
else
    error('unknown stimtype');
end
Lwhichstimtypes = Lstimtype | Lblank;

% Rasters for all conditions
offset = [-.1 .1];  % pre and post time wrt stimon/stimoff
bins = offset(1):binwidth:dur+offset(2);
%figure('position',[900 80 375 725]); axes('units','centimeters','position',[2 3 9 20]); hold on; % One column
%set(gca,'TickDir','out'); 
figprefs; axes('position',[2 3 7 15]);
counter = 0;
for j = find(Lwhichstimtypes)'
    L = Lcc == uniquestim(j,1) & Mcc == uniquestim(j,2) & TF == uniquestim(j,3);
    psth = zeros(size(bins));
    spikes_and_trialidxs = [];
    for i = find(L)'
        tmpspikes = spikes{i}-stimon_t(i);
        tmpspikes(tmpspikes < offset(1) | tmpspikes > dur+offset(2)) = [];
        spikes_and_trialidxs = [spikes_and_trialidxs; tmpspikes, repmat(i,length(tmpspikes),1)];
    end
    % plotting the PSTH first for ease of viewing
    psth = hist(spikes_and_trialidxs(:,1),bins);
    h = plot(bins,counter+1+psth./sum(L)./binwidth/10,'-','Linewidth',2,'Color',[.5 .5 .5]);
    for i = find(L)'
        tmpspikes = spikes_and_trialidxs(spikes_and_trialidxs(:,2) == i,1);
        nspikestot = length(tmpspikes);
        h = plot([tmpspikes tmpspikes]',[zeros(nspikestot,1) 1*ones(nspikestot,1)]'+counter,'k-','linewidth',1);
        counter = counter + 1;
    end
    if j ~= max(find(Lwhichstimtypes))
        plot([offset(1) dur+offset(2)],counter*[1 1],'k-');
    end
    
    if uniquestim(j,3) == 0
        h = text(-.15,counter-sum(L)/2,'Blank');
        h(2) = text(dur+offset(2)+.15,counter-sum(L)/2,'0.00');
    else
        h = text(-.15,counter-sum(L)/2,num2str(round(uniquestim(j,3)*10)/10));
        % BELOW: Using L-cone contrast == M-cone contrast as a proxy for luminance contrast
        % (see BenardeteKaplanSims.m)
        h(2) = text(dur+offset(2)+.15,counter-sum(L)/2,num2str(round(uniquestim(j,1)*100)/100));
    end
    set(h,'HorizontalAlignment','right','FontSize',12,'FontAngle','Italic');
    %if uniquestim(j,1) == 0 & uniquestim(j,2) == 0
    %    plot([0 dur dur 0 0],[counter counter 0 0 counter],'k-','LineWidth',1);
    %end
    counter = counter + 1;
end
plot([0 0],[0 counter],'k-','LineWidth',.5);
plot([dur dur],[0 counter],'k-','LineWidth',.5);
set(gca,'Xlim',[0+offset(1) dur+offset(2)],'Ylim',[0 counter],'Ytick',[],'Xtick',[0 .2 .4 .6 .8],'Box','off','FontSize',12,'FontAngle','italic');
xlabel('Time (s)');
set(gcf,'Renderer','painters');

% Dprime plot
figprefs; axes;
se = [];
niter = 200; % Boot strapping SEs
for i = find(Lstimtype)'
    tmp = [];
    for j = 1:niter
        n_signal = length(signal{i});
        n_noise = length(noise{i});
        s = signal{i}(unidrnd(n_signal,n_signal,1));
        n = noise{i}(unidrnd(n_noise,n_noise,1));
        pooled_var = ((n_signal-1)*nanvar(s)+(n_noise-1)*nanvar(n))/(n_signal+n_noise-2);
        tmp(j) = (nanmean(s)-nanmean(n))/sqrt(pooled_var);
    end
    se = [se; std(tmp)];
end
patch([uniquestim(Lstimtype,3)', flipud(uniquestim(Lstimtype,3))'],[dprime(Lstimtype)'+se', fliplr(dprime(Lstimtype)'-se')],[.5 .5 .5],'Facealpha',.5);
plot(uniquestim(Lstimtype,3),dprime(Lstimtype),'ko-','LineWidth',1,'MarkerFaceColor','white');
set(gca,'Xscale','log');
set(gca,'FontSize',12,'FontAngle','italic','XtickLabel',{'1','10','100'},'Ylim',[-.5 4]);
xlabel('Temporal frequency (Hz)');
ylabel('SNR');
filename = stro.sum.fileName(find(stro.sum.fileName == filesep,1,'last')+1:end-4);
text(20,-.25,filename,'FontAngle','Italic','FontSize',9);

%%
% Section 5
% Input efficiency (d' of LGN neurons/d' of cone model)

data = [];
for MONKEY = MONKEYS
    for CELLTYPE = CELLTYPES
        [filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',CELLTYPE,'subjID',{MONKEY{1}(1)});
        for a = 1:length(filenames)
            stro = {};
            for j = 1:length(filenames{a})
                stro{j} = nex2stro(char(findfile(filenames{a}(j), fullfile(nexfilepath,'Greg',MONKEY))));
            end
            stro = strocat(stro);
            sigma_gabor = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sigma'))); % DVA
            
            [uniquestim, dprime] = IsoSampGetDPrime(stro,DPRIMEMETHOD,spikecds(a));
            
            params.stro = stro;
            params.eyenumber = 1;
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
            if strcmp(CELLTYPE,'M')
                rf_diam_deg = ecc_to_diam_deg_M(rf_r_deg); % diameter RF
            elseif strcmp(CELLTYPE,'P')
                rf_diam_deg = ecc_to_diam_deg_P(rf_r_deg); % diameter RF
            else
                error('unknown cell type');
            end
            RF_STD = rf_diam_deg/2; % 1 standard deviation of Gaussian RF
            params.gab.sd = sqrt(1/((1/RF_STD.^2)+(1/sigma_gabor.^2))); % Product of RF Gaussian and Gabor Gaussian
            params.gab.sd = RF_STD;
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
            tmp = nan*ones(1,length(TFBINCENTERS));
            for i = 1:length(TFBINCENTERS)
                Lblank = uniquestim(:,1) == 0 & uniquestim(:,2) == 0 & uniquestim(:,3) == 0;
                Llum = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & ~Lblank;
                LTF = uniquestim(:,3) >= TFBINEDGES(i) & uniquestim(:,3) < TFBINEDGES(i+1);
                L = Llum & LTF;
                tmp(i) = mean(dprime(L)./cone_dprimes(L));
            end
            data = [data; find(strcmp(MONKEY,MONKEYS)) find(strcmp(CELLTYPE,CELLTYPES)) tmp];     
        end
    end
end

figprefs;
hax(1) = axes('position',[4 13 9 9]);
hax(2) = axes('position',[4 2 9 9]);
niter = 200;
for monkey_idx = 1:2
    axes(hax(monkey_idx));
    for celltype_idx = 1:2
        L = data(:,1) == monkey_idx & data(:,2) == celltype_idx;
        med = nanmedian(data(L,3:end));
        sem = nan*ones(1,length(med));
        for i = 1:length(sem)
            n = sum(L);
            tmp = data(L,2+i);
            sem(i) = std(nanmedian(tmp(unidrnd(n,n,niter))));
        end
        Lbins = ~isnan(med);
        h_patch = patch([TFBINCENTERS(Lbins), fliplr(TFBINCENTERS(Lbins))],[med(Lbins)+sem(Lbins), fliplr(med(Lbins)-sem(Lbins))],[.5 .5 .5],'Facealpha',.5);
        h_line = plot(TFBINCENTERS(Lbins), med(Lbins),'k-o','LineWidth',2);
        if celltype_idx == 2
            set(h_line,'MarkerSize',3);
        else
            set(h_line,'MarkerEdgeColor','black','MarkerFaceColor','black');
        end
        %if monkey_idx ==2
        %    set(h_line,'Marker',MONKEY2MARKER);
        %end
    end
    set(gca,'Xscale','log','Xlim',[TFBINEDGES(1) TFBINEDGES(end)],'FontSize',15,'FontAngle','italic');
    set(gca,'Ytick',[0 1 2 3],'XTickLabel',{'1','10'},'Ylim',[-.25 3.25]);
    plot([TFBINEDGES(1) TFBINEDGES(end)],[0 0],'k-');
    plot([TFBINEDGES(1) TFBINEDGES(end)],[1 1],'k-');
    ylabel('Efficiency');
    xlabel('Frequency (Hz)');
end
set(gcf,'Renderer','Painters');

%%
% Section 6
% Population d' from LGN neurons and cones (also cone absorptions)
% 8/5/19: allowing analysis to be restricted to single LGN RFs.

POPULATIONSCALING = 1;
STANDARDPARAMS = 1;
COLORDIR = 'RG';

if STANDARDPARAMS
    RFSCALEFACTOR = 1;
    R_WITHIN_OVERRIDE = [];
    if strcmp(COLORDIR,'LUM')
        OFFSETS = []; % Should be set relative to stimon and stimoff?
    elseif strcmp(COLORDIR,'RG')
        OFFSETS = [.12 .786]; % adding an offset to the counting window. This only affects LGN idlObs.
    end
else % New parameters
    % Constants below are for Reviewer #3. Changing parameter values and see what
    % happens.
    RFSCALEFACTOR = 1; % 2, 0.5 Changing RF *diameter* so area/density changes are greater
    R_WITHIN_OVERRIDE = [];
%    RFTRUNCATIONINSD = 2;
    OFFSETS = []; % adding an offset to the counting window. This only affects LGN idlObs.
end

if ~STANDARDPARAMS % Halving or doubling RF diameter
    ecc_to_diam_deg_P_debug = @(x)(RFSCALEFACTOR*sqrt(2./(sqrt(3).*... % Equation 9. Distance between adjacent midget RF centers.
    (2*dc_0.*(1+x./rm).^-1.*(a*(1+(x./r2)).^-2+(1-a)*exp(-x./re)))...
    ./2))... % Dividing y by 2 to estimate RF size from only ON or OFF mosaics (halving the density).
    *HUMAN2MONKPSCALEFACTOR); % Monkey midget RFs are slightly smaller than human midget RFs
    a = 0.9729; % Table 1
    r2 = 1.084; % Table 1
    re = 7.633; % Table 1
    dc_0 = 14804.6; % Cone density of fovea
    rm = 41.03; % See Equation 7

    ecc_to_diam_deg_M_debug = @(rf_r_deg) RFSCALEFACTOR*(10.^(-1.2459+0.0345*rf_r_deg)); % temporal retina equivalent (Derrington and Lennie)
end

data = cell(2,2);
for MONKEY = MONKEYS
    for CELLTYPE = CELLTYPES
        [filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',CELLTYPE,'subjID',{MONKEY{1}(1)});     
        for i = 1:length(filenames)
            tmp = [];
            cal = [];
            stro = {};
            for j = 1:length(filenames{i})
                stro{j} = nex2stro(char(findfile(filenames{i}(j), fullfile(nexfilepath,'Greg',MONKEY))));
            end
            stro = strocat(stro);
            sigmas_n = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sigmas_n')));
            gabor_sigma = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sigma')));

            [tmp.uniquestim, tmp.dprime] = IsoSampGetDPrime(stro,DPRIMEMETHOD,spikecds(i),OFFSETS); % GDLH 12/12/19 Added OFFSET to LGN only!
            
            Lblank = tmp.uniquestim(:,1) == 0 & tmp.uniquestim(:,2) == 0 & tmp.uniquestim(:,3) == 0;
            Llum = sign(tmp.uniquestim(:,1)) == sign(tmp.uniquestim(:,2)) & ~Lblank;
            
            % Now getting population scale factor
            % Getting RF size at rfx, rfy
            if strcmp(CELLTYPE,'M')
                if ~STANDARDPARAMS
                    ecc_to_diam_deg = ecc_to_diam_deg_M_debug;
                else
                    ecc_to_diam_deg = ecc_to_diam_deg_M;
                end
            else
                if ~STANDARDPARAMS
                    ecc_to_diam_deg = ecc_to_diam_deg_P_debug;
                else
                    ecc_to_diam_deg = ecc_to_diam_deg_P;
                end
            end
            if POPULATIONSCALING
                %keyboard
                if ~STANDARDPARAMS
                    [tmp.population_scalefactor, tmp.ncells] = IsoSampGetPopulationScaleFactor(stro, ecc_to_diam_deg, TEMPORONASALSCALEFACTOR,RFTRUNCATIONINSD, ONOFFCORRELATION, 2, R_WITHIN_OVERRIDE)
                else
                    [tmp.population_scalefactor, tmp.ncells] = IsoSampGetPopulationScaleFactor(stro, ecc_to_diam_deg, TEMPORONASALSCALEFACTOR,RFTRUNCATIONINSD, ONOFFCORRELATION, 2)
                end
            else
                tmp.population_scalefactor = 1; tmp.ncells = 1;
            end
            rfx = stro.sum.exptParams.rf_x/10;
            rfy = stro.sum.exptParams.rf_y/10;
            tmp.ecc = sqrt((rfx./TEMPORONASALSCALEFACTOR)^2+rfy^2);

            % Now getting cone ideal observer sensitivity
            % Cut and paste from IsoSampPop section 4.1
            params.stro = stro;
            spds = params.stro.sum.exptParams.mon_spd;
            spds = SplineSpd([380:4:780]',reshape(spds,length(spds)/3,3),[380:5:780]');
            M = fundamentals.T_cones_smj10*spds;
            cal.monSpect = spds(:);
            cal.Mmtx = M(:);
            cal.frameRate = params.stro.sum.exptParams.framerate;
            cal.bkgndrgb = params.stro.sum.exptParams.bkgndrgb';
            cal.fname = 'test';
            cal.monSpectWavelengths = linspace(380,780,101);
            cal.pixperdeg = params.stro.sum.exptParams.pixperdeg;
            params.monCalFile = cal;
            if POPULATIONSCALING
                params.gab.sd = gabor_sigma;
                params.eyeNumber = 2;
            else
                params.gab.sd = ecc_to_diam_deg(tmp.ecc)/2; % overwriting gabor SD with RF SD
                params.gab.sd = params.gab.sd; % DEBUGGING: artificially inflating RF SD by 5x
                params.eyeNumber = 1;
            end
            
            [gab, cones, mon, idlob, params] = DTcones_gh(params,0);
            
            conedata = [];
            for k = 1:size(idlob.analyticMean,1) % looping over color direction
                for j = 1:size(idlob.analyticMean(k,:),2) % looping over contrast/TF
                    if ~isempty(idlob.analyticMean{k,j})
                        tmp_lm_mu = idlob.analyticMean{k,j}([1 2]);
                        tmp_lm_var = idlob.analyticVar{k,j}([1 2]);
                        tf = gab.driftRates{k}(j);
                        conedata = [conedata; gab.colorDirs(k,[1 2]).*gab.contrasts{k}(j) tf tmp_lm_mu tmp_lm_var];
                    end
                end
            end
            
            % Using uniquestim to order the rows of tmpdata and calculating cone
            % dprimes
            tmp.conedprime = [];
            for j = 1:size(tmp.uniquestim,1)
                L = all(abs(tmp.uniquestim(j,:)-conedata(:,[1 2 3]))<1e-10,2);
                if sum(L) ~= 1
                    if all(tmp.uniquestim(j,:) == 0)
                        tmp.conedprime(j) = nan;
                    else
                        error('sorting problem');
                    end
                else
                    v = conedata(L,[6 7]); % variance
                    m = conedata(L,[4 5]); % mean
                    tmp.conedprime(j) = sqrt(m.^2*(1./v'));
                end
            end
            
            % Need to calculate flashTimeProfile on the basis of timewindoffset
            t = 0:1/mon.frameRate:gab.length;
            nframes = length(t);
            flashTimeProfile = ones(1,nframes);
            ramp = linspace(0,1,nframes/4);
            flashTimeProfile(1:length(ramp)) = ramp;
            flashTimeProfile(end:-1:end-length(ramp)+1) = ramp;
            sigmaInPix = params.gab.sd*cal.pixperdeg;
            
            tmp.photondprime = IsoSampGetPhotonDPrime(flashTimeProfile, mon.frameRate, mon.bkgndlms_Rstar, sigmaInPix, cat(3,cones.num_L,cones.num_M), tmp.uniquestim);
            % Adding to the "data" cell array of cell arrays
            listsofar = data{strcmp(MONKEY,MONKEYS),strcmp(CELLTYPE,CELLTYPES)};
            listsofar{length(listsofar)+1} = tmp;
            data{strcmp(MONKEY,MONKEYS),strcmp(CELLTYPE,CELLTYPES)}=listsofar;
        end
    end
end

% Doing the plotting
figure('Position',[440 100 750 700],'DefaultAxesTickDirMode','manual','DefaultAxesTickdir','out','DefaultAxesYcolor','black','DefaultAxesXcolor','black')
set(gcf,'DefaultAxesFontSize',15,'DefaultAxesFontAngle','italic','DefaultAxesUnits','centimeters');
hax = [];
hax(1,1) = axes('position',[2 2 9 9]); hold on;
hax(1,2) = axes('position',[2 14 9 9]); hold on;
hax(2,1) = axes('position',[14 2 9 9]); hold on;
hax(2,2) = axes('position',[14 14 9 9]); hold on;

for monkey_idx = 1:2
    for celltype_idx = 1:2
        monkey_celltype_data = data{monkey_idx,celltype_idx};
        tmp = nan*ones(4,length(TFBINCENTERS),length(monkey_celltype_data)); % rows: LGN dprime, cone dprime, photon dprime, ns
        for i = 1:length(monkey_celltype_data)
            cell_data = monkey_celltype_data{i};
            Lblank = cell_data.uniquestim(:,1) == 0 & cell_data.uniquestim(:,2) == 0 & cell_data.uniquestim(:,3) == 0;
            Llum = sign(cell_data.uniquestim(:,1)) == sign(cell_data.uniquestim(:,2)) & ~Lblank;
            Lrg = sign(cell_data.uniquestim(:,1)) ~= sign(cell_data.uniquestim(:,2)) & ~Lblank;
            
            if strcmp(COLORDIR,'RG')
                LlumORrg = Lrg;
            else
                LlumORrg = Llum;
            end
            for j = 1:length(TFBINCENTERS)
                Ltf = cell_data.uniquestim(:,3) > TFBINEDGES(j) & cell_data.uniquestim(:,3) <= TFBINEDGES(j+1);
                if sum(Ltf&LlumORrg) > 0
                    tmp(1,j,i) = mean(cell_data.dprime(LlumORrg&Ltf))*cell_data.population_scalefactor;
                    tmp(2,j,i) = mean(cell_data.conedprime(LlumORrg&Ltf));
                    tmp(3,j,i) = mean(cell_data.photondprime(LlumORrg&Ltf));
                    tmp(4,j,i) = 1; % each cell counts as an independent entity
                end
            end
        end
        %allnans = all(isnan(squeeze(tmp(4,:,:))),2);
        %tmp = tmp(:,~allnans,:);
        n = nansum(tmp(4,:,:),3);
        
        % Getting rid of TFs with only 1 data point
        tmp = tmp(:,n>1,:);
        n = n(1:size(tmp,2));
        % End of new hacky change
        
        mn_lgn = nanmean(tmp(1,:,:),3);
        sd_lgn = nanstd(tmp(1,:,:),0,3);
        sem_lgn = sd_lgn./sqrt(n);
        mn_cone = nanmean(tmp(2,:,:),3);
        sd_cone = nanstd(tmp(2,:,:),0,3);
        sem_cone = sd_cone./sqrt(n);
        mn_photon = nanmean(tmp(3,:,:),3);
        sd_photon = nanstd(tmp(3,:,:),0,3);
        sem_photon = sd_photon./sqrt(n);
        axes(hax(monkey_idx,celltype_idx));
        TBC = TFBINCENTERS(1:length(n));
        patch([TBC, fliplr(TBC)],[mn_photon+sem_photon, fliplr(mn_photon-sem_photon)],[0 .7 .7],'Facealpha',.5,'LineStyle','none');
        plot(TBC, mn_photon,'c-o','LineWidth',1,'MarkerEdgeColor','black','MarkerFaceColor','cyan');
        patch([TBC, fliplr(TBC)],[mn_cone+sem_cone, fliplr(mn_cone-sem_cone)],[1 .5 .5],'Facealpha',.5,'LineStyle','none');
        plot(TBC, mn_cone,'r-o','LineWidth',1,'MarkerEdgeColor','black','MarkerFaceColor','red');
        patch([TBC, fliplr(TBC)],[mn_lgn+sem_lgn, fliplr(mn_lgn-sem_lgn)],[.5 .5 .5],'Facealpha',.5,'LineStyle','none');
        plot(TBC, mn_lgn,'k-o','LineWidth',1,'MarkerEdgeColor','black','MarkerFaceColor','black');
        plot([TFBINEDGES(1) TFBINEDGES(end)],[1.27 1.27],'k-');
        if POPULATIONSCALING
            set(gca,'Xscale','log','Xlim',[TBC(1) TBC(end)],'Ylim',[-2 40]);
        else
            set(gca,'Xscale','log','Xlim',[TBC(1) TBC(end)],'Ylim',[-2 10]);
        end
        if ~STANDARDPARAMS & strcmp(COLORDIR,'LUM')
            set(gca,'Ylim',[-2 20]); % Zooming in on LGN d's
        end
        xlabel('Frequency (Hz)');
        ylabel('SNR');
        if STANDARDPARAMS
            title([MONKEYS{monkey_idx},' ',CELLTYPES{celltype_idx}])
        else
           title(['RFSCALE = ',num2str(RFSCALEFACTOR),' R within= ',num2str(R_WITHIN_OVERRIDE)]);
        end
        disp([MONKEYS{monkey_idx},' ',CELLTYPES{celltype_idx},' ',num2str(TBC(end)),' Hz, photon d''=',num2str(mn_photon(end))]) 
    end
end
equatesubplotaxeslims;
set(gcf,'Renderer','painters');


%%
% Section 6.1 Getting numbers of neurons for reviewer 1
% Need to run section 6 first.

figprefs; axes;
hold on;
for monkey_idx = 1:2
    for celltype_idx = 1:2
        monkey_celltype_data = data{monkey_idx,celltype_idx};
        for counter = 1:length(monkey_celltype_data)
            h = plot(monkey_celltype_data{counter}.ecc,2*monkey_celltype_data{counter}.ncells,'ko'); % 2* because 2 mosaics
            if monkey_idx == 2
                set(h,'Marker',MONKEY2MARKER);
            end
            if celltype_idx == 2
                set(h,'MarkerEdgeColor','black','MarkerFaceColor',[.5 .5 .5]);
            else
                set(h,'MarkerEdgeColor','black','MarkerFaceColor','white','MarkerSize',6);
            end
        end
    end
end
xlabel('eccentricity (�)')
ylabel('number of LGN neurons')
        
% The window I'm looking at is 1.2 x 1.2 deg = 1.4 deg^2
ecc = linspace(2,18,100);
p = 1011688*(ecc+2.9144).^-2.67; % formula from where?
m = 2620.2*((ecc-1.8322).^2+5.5638).^-0.801;

plot(ecc,p*1.4,'b-'); % 1.4 because I'm counting ncells in a 1.2 x 1.2� window. Remember some P are really K cells.
plot(ecc,m*1.4,'m-');
set(gca,'YScale','log','Ylim',[100 10000],'Xlim',[0 18],'TickDir','out');

%%
% Section 6.2 RF sizes
data =[];
for MONKEY = MONKEYS
    for CELLTYPE = CELLTYPES
        [filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',CELLTYPE,'subjID',{MONKEY{1}(1)});     
        for i = 1:length(filenames)
            stro = {};
            for j = 1:length(filenames{i})
                stro{j} = nex2stro(char(findfile(filenames{i}(j), fullfile(nexfilepath,'Greg',MONKEY))));
            end
            stro = strocat(stro);
            rfx = stro.sum.exptParams.rf_x/10;
            rfy = stro.sum.exptParams.rf_y/10;
            ecc = sqrt((rfx./TEMPORONASALSCALEFACTOR)^2+rfy^2);
            if strcmp(CELLTYPE,'M')
                ecc_to_diam_deg = ecc_to_diam_deg_M;
            else
                ecc_to_diam_deg = ecc_to_diam_deg_P;
            end
            data =[data; ecc_to_diam_deg(ecc)];
        end
    end
end
sigma_gabor = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sigma'))); % DVA

((sigma_gabor*4)/mean(data))^2 % How much larger the Gabor integration window is than the single RF integration window

%%
% Section 6.3
% Discussion figure: proportion information lost as a function of temporal
% frequency (for Reviewer #2). Need to run section 6 first.

OMITPCELLS = 1;
summarydata2 = [];
for j = 1:length(TFBINCENTERS)
    summarydata1 = [];
    for celltype_idx = 1:2
        for monkey_idx = 1:2  
            monkey_celltype_data = data{monkey_idx,celltype_idx};
            ncells = length(monkey_celltype_data);
            tmp = nan*ones(ncells,3); % ncells x [LGN, cone_current, photon absorption]
            for i = 1:ncells
                cell_data = monkey_celltype_data{i};
                Lblank = cell_data.uniquestim(:,1) == 0 & cell_data.uniquestim(:,2) == 0 & cell_data.uniquestim(:,3) == 0;
                Llum = sign(cell_data.uniquestim(:,1)) == sign(cell_data.uniquestim(:,2)) & ~Lblank;
                Ltf = cell_data.uniquestim(:,3) > TFBINEDGES(j) & cell_data.uniquestim(:,3) <= TFBINEDGES(j+1);
                if sum(Ltf&Llum) > 0
                    tmp(i,1) = mean(cell_data.dprime(Llum&Ltf))*cell_data.population_scalefactor;
                    tmp(i,2) = mean(cell_data.conedprime(Llum&Ltf));
                    tmp(i,3) = mean(cell_data.photondprime(Llum&Ltf));
                end
            end
            summarydata1 = [summarydata1; celltype_idx nanmean(tmp)]; % Accumulating averages (across cells, within cell type)
        end
    end
    % setting P-cell d' to nan is mean d' is < 1.27
    if OMITPCELLS | nanmean(summarydata1(:,2)) < 1.27
        summarydata1(summarydata1(:,1) == 2,2) = nan;
        disp(['got here. TF = ', num2str(TFBINEDGES(j))]);
    end
    summarydata2 = [summarydata2; nanmean(summarydata1(:,2:end))];
end

fraction_left = (summarydata2 - 1.27); % greg = 0 means behavioral sensitivity (d' = 1.27)
fraction_left = fraction_left./repmat(fraction_left(:,3),1,3);
figprefs; axes;
patch([TFBINCENTERS,fliplr(TFBINCENTERS)],[fraction_left(:,1)',zeros(1,size(fraction_left,1),1),],[0 0 0]);
patch([TFBINCENTERS,fliplr(TFBINCENTERS)],[fraction_left(:,2)',fliplr(fraction_left(:,1)'),],[.5 .5 .5]);
patch([TFBINCENTERS,fliplr(TFBINCENTERS)],[fraction_left(:,3)',fliplr(fraction_left(:,2)'),],[1 1 1]);
axis tight;
set(gca,'Ylim',[0 1],'Xscale','log','Xlim',[TFBINCENTERS(1) TFBINCENTERS(end)]);
xlabel('Frequency (Hz)');
ylabel('Fractional signal-to-noise ratio');
set(gcf,'Renderer','painters');
text(5,.7,'Phototransduction','FontName','Helvetica','FontAngle','italic');
text(1.8,.25,'Cones \rightarrow LGN','FontName','Helvetica','FontAngle','italic');

%%
% Section 7
% Looking for differences in sensitivity between ON and OFF cells

conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
whitenoisedata = fetch(conn, 'SELECT fileID, rfX, rfY, neuron, spikeCode FROM WhiteNoiseLGN_forIS WHERE quality = ''1''');
isosampdata = fetch(conn, 'SELECT fileID, rfX, rfY, neuron, cellClass FROM IsoSamp_LGN WHERE quality = ''1''');
close(conn);

uniqueWNneurons = unique([whitenoisedata{:,4}]);
WNfilenames = {};
IsoSampfilenames  = {};
for neuronid = uniqueWNneurons
    LWN =  [whitenoisedata{:,4}] == neuronid;
    Lisosamp = [isosampdata{:,4}] == neuronid;
    if any(Lisosamp)
        tmp = [];
        for i = find(LWN)       
        	tmp = [tmp, whitenoisedata(i,1)];
        end
        WNfilenames{length(WNfilenames)+1} = tmp;
        tmp = [];
        for i = find(Lisosamp)       
        	tmp = [tmp, isosampdata(i,1)];
        end
        IsoSampfilenames{length(IsoSampfilenames)+1} = tmp;
        % Sanity checks
        tmp=[];
        for i = find(Lisosamp)
            tmp=[tmp;isosampdata{i,2:4}];
        end
        for i = find(LWN)
            tmp=[tmp;whitenoisedata{i,2:4}];
        end
        if any(std(tmp))
            error('White noise and IsoSamp entries into database disagree');
        end
    end
end
if length(WNfilenames) ~= length(IsoSampfilenames)
   error('different numbers of white noise and isosamp files');
end

data = [];
for i = 1:length(WNfilenames)
    WNstro=[];
    for j = 1:length(WNfilenames{i})
        WNstro = strocat(WNstro,nex2stro(findfile(WNfilenames{i}{j})));
    end
    neuronidx = whitenoisedata{strcmp(whitenoisedata(:,1),WNfilenames{i}(1)),4};
    celltype = isosampdata{[isosampdata{:,4}] == neuronidx,5};
    spikename = whitenoisedata{strcmp(whitenoisedata(:,1),WNfilenames{i}(1)),5};
    spikeidx = find(strcmp(WNstro.sum.rasterCells(1,:),spikename));
    nstixperside = WNstro.sum.exptParams.nstixperside;

    IsoSampstro = [];
    for j = 1:length(IsoSampfilenames{i})
        IsoSampstro = strocat(IsoSampstro,nex2stro(findfile(IsoSampfilenames{i}{j})));
    end

    % Reconstructing the M matrix and gamma table
    fundamentals = WNstro.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = WNstro.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;

    out = getWhtnsStats(WNstro,MAXT,'STCOVmex', {nstixperside^2, 3, maxT}, spikename);
    STA = out{1}/out{3};
    inRF = getRFfromSTA(STA,GETRFFROMSTAMODE,GETRFFROMSTATHRESH);
    STA = reshape(STA,[nstixperside.^2  3 MAXT]);
    temporalSTA = squeeze(STA(logical(inRF(:)),:,:));
    if ~ismatrix(temporalSTA)
        temporalSTA = squeeze(mean(temporalSTA,1));
    end
   
    switch(celltype)
        case {'P'}
            data(i,1) = 1;
        case {'M'}
            data(i,1) = 2;
        case {'K'}
            data(i,1) = 3;
        otherwise
            data(i,1) = 4;
    end
    [u,s,v] = svd(temporalSTA);
    rgb = u(:,1);
    temporal = v(:,1);
    if max(abs(rgb)) == -min(rgb)
        rgb = -rgb;
        temporal = -temporal;
    end
    if max(abs(temporal)) == max(temporal)
        data(i,2) = 1; % ON
    else
        data(i,2) = 0; % OFF
    end
    lm_ratio = temporalSTA(1,2)/temporalSTA(2,2);
    if abs(lm_ratio-M(2,1)/M(2,2)) < abs(lm_ratio-M(1,1)/M(1,2))
        data(i,3) = 1; % M
    else
        data(i,3) = 0; % L
    end
    [uniquestim, dprime] = IsoSampGetDPrime(IsoSampstro,DPRIMEMETHOD,spikeidx);
    
    Lblank = uniquestim(:,1) == 0 & uniquestim(:,2) == 0 & uniquestim(:,3) == 0;
    Llum = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & ~Lblank;
    data(i,4) = mean(dprime(Llum));
end
% Doing statistical comparisons
L_P = data(:,1) == 1; % Parvocellular neurons
L_M = data(:,1) == 2; % Magnocellular neurons
L_OFF = data(:,2) == 0; % OFF neurons
L_ON = data(:,2) == 1; % ON neurons
L_Lcone = data(:,3) == 0; % L-cone dominated
L_Mcone = data(:,3) == 1; % M-cone dominated

[~,p_p_on_off] = ttest2(data(L_P & L_OFF,4), data(L_P & L_ON,4)) % parvo ON vs OFF
[~,p_p_l_m] = ttest2(data(L_P & L_Lcone,4), data(L_P & L_Mcone,4)) % parvo L vs M

[~,p_m_on_off] = ttest2(data(L_M & L_OFF,4), data(L_M & L_ON,4)) % magno ON vs OFF
[~,p_m_p] = ttest2(data(L_M,4), data(L_P,4)) % magno vs parvo


%%
% Section 8
% Microsaccade-triggered PSTHs

PREOFFSET = 0.05; % How long before saccade initiation to collect spikes (s)
POSTOFFSET = 0.2; % How long after saccade initiation to collect spikes (s)
nbins = 30;
bins = linspace(-PREOFFSET,POSTOFFSET,nbins);
binwidth = bins(2)-bins(1); % s

timebins = [-0.1:.05:.766];
timebinwidth = timebins(2)-timebins(1);
%deltaT = timebins(2)-timebins(1);
%nyquist = 1./(2*deltaT);
msacpsth = [];

figprefs;
hax(1,1) = axes('Position',[2 2 4 4]);
hax(1,2) = axes('Position',[2 12 4 4]);
hax(2,1) = axes('Position',[12 2 4 4]);
hax(2,2) = axes('Position',[12 12 4 4]);
for MONKEY = MONKEYS
    for CELLTYPE = CELLTYPES
        fftdata =[];
        msacpsth = [];
        nsacs = 0;
        [filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',CELLTYPE,'subjID',{MONKEY{1}(1)});
        PSTHs = [];
        for a = 1:length(filenames)
            stro = {};
            for j = 1:length(filenames{a})
                stro{j} = nex2stro(char(findfile(filenames{a}(j), fullfile(nexfilepath,'Greg',MONKEY))));
            end
            sumspikes = zeros(1,length(bins));
            stro = strocat(stro);
            Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_l'));
            Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_m'));
            TF = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
            stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
            stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
            dur = mean(stimoff_t-stimon_t);
            Lblank = Lcc == 0 & Mcc == 0 & TF == 0;
            Llum = Lcc > 0 & Mcc > 0;
            Lrg = sign(Lcc) ~= sign(Mcc);
            
            whichtrials = find(Lblank);
            whichtrials = find(Llum);
            
            ntrials = length(whichtrials);
            fpacq_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpacq_t'));
            stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
            spiketimes = stro.ras(:,spikecds(a));
            sacstats = getSacData(stro,1);
            for j = 1:ntrials
                trialidx = whichtrials(j);
                L = ones(length(sacstats.starttimes{trialidx}),1);
                L(sacstats.starttimes{trialidx} < fpacq_t(trialidx)) = 0;
                L(sacstats.starttimes{trialidx} > stimoff_t(trialidx)) = 0;
                if any(L)
                    for k = find(L')
                        tmpspiketimes = spiketimes{trialidx};
                        Lspikes = tmpspiketimes > sacstats.starttimes{trialidx}(k)-PREOFFSET & tmpspiketimes < sacstats.starttimes{trialidx}(k)+POSTOFFSET;
                        [tmp,~] = hist(tmpspiketimes(Lspikes)-sacstats.starttimes{trialidx}(k), bins);
                        sumspikes = sumspikes + tmp;
                        nsacs = nsacs+1;
                    end
                end
            end
            PSTHs = [PSTHs; sumspikes/nsacs/binwidth];
        end
        axes(hax(strcmp(MONKEY,MONKEYS),strcmp(CELLTYPE,CELLTYPES)));
        meanPSTH = mean(PSTHs);
        sdPSTH = sqrt(var(PSTHs));
        semPSTH = sdPSTH./sqrt(length(filenames));
        patch([bins(2:end-1), fliplr(bins(2:end-1))],[meanPSTH(2:end-1)+semPSTH(2:end-1), fliplr(meanPSTH(2:end-1)-semPSTH(2:end-1))],[.5 .5 .5],'Facealpha',.5);
        plot(bins(2:end-1), meanPSTH(2:end-1),'k-','LineWidth',2);
        set(gca,'Ylim',[0 20],'Xlim',[bins(2) bins(end-1)]);
        xlabel('Time (s)');
        ylabel('Firing rate (spikes/s)');
        text(0,15,[MONKEY{1},' ',CELLTYPE{1}],'FontSize',12,'FontAngle','italic');
        plot([0 0],[0 15],'k:');
    end
end
set(gcf,'Render','painters');

%%
% Section 9
% Microsaccade frequency over time in a trial. and comparison between
% monkeys.
% Taken from IsoSampPop, Section 10.1

PREOFFSET = -0.1; % How long before stimulus on to collect saccades
POSTOFFSET = 0.1; % How long after stimulus off to collect saccades
nbins = 20;
bins = linspace(PREOFFSET,.666+POSTOFFSET,nbins);
binwidth = bins(2)-bins(1);
figprefs;
hax = axes('Position',[2 2 7 7]);
hax(2) = axes('Position',[2 12 7 7]);
nsacspersec = [];
for MONKEY = MONKEYS
    [filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','subjID',{MONKEY{1}(1)});
    means = nan*ones(length(filenames),nbins);
    
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
        sacstats = getSacData(stro,1);
        Lblank = Lcc == 0 & Mcc == 0 & TF == 0;
        Llum = Lcc > 0 & Mcc > 0;
        whichtrials = find(Llum | Lblank)';

        tmpdata = zeros(length(whichtrials),nbins);
        rowcounter = 1;
        for j = whichtrials
            sactimes = sacstats.starttimes{j}-stimon_t(j);
            sactimes(sactimes<0+PREOFFSET-binwidth/2) = [];
            sactimes(sactimes>.666+POSTOFFSET+binwidth/2) = [];
            [n,x] = hist(sactimes, bins);
            tmpdata(rowcounter,:) = n;
            rowcounter = rowcounter+1;
        end
        means(i,:) = mean(tmpdata./binwidth);
        nsacs = sum(sum(tmpdata));
        totaltime = length(whichtrials)*(bins(end)-bins(1));
        nsacspersec=[nsacspersec; find(strcmp(MONKEY,MONKEYS)) nsacs/totaltime]; % monkey, sacs per sec
    end
    axes(hax(strcmp(MONKEY, MONKEYS)));
    sem = sqrt(var(means)./size(means,1));
    mn = mean(means);
    patch([bins, fliplr(bins)],[mn+sem, fliplr(mn-sem)],[.5 .5 .5],'Facealpha',.5);
    plot(bins,mn,'k-','LineWidth',2);
    set(gca,'Xlim',[bins(1), bins(end)],'Ylim',[0 4])
    xlabel('Time (s)');
    ylabel('Saccades/s');
end

L = nsacspersec(:,1) == 1;
[mean(nsacspersec(L,2)) mean(nsacspersec(~L,2))]
ranksum(nsacspersec(L,2),nsacspersec(~L,2))
%%
% Section 10
% Noise spectra + ISI distributions

timebins = [0:.005:.660];
deltaT = timebins(2)-timebins(1);
nyquist = 1./(2*deltaT);

figprefs;
hax(1,1) = axes('Position',[2 2 7 7]);
hax(1,2) = axes('Position',[2 12 7 7]);
hax(2,1) = axes('Position',[12 2 7 7]);
hax(2,2) = axes('Position',[12 12 7 7]);
ISIbins = [0:5:400]/1000; % s
ISIdata = zeros(2,2,length(ISIbins)); % Monkey, cell type, bins
for MONKEY = MONKEYS
    for CELLTYPE = CELLTYPES
        fftdata =[];
        [filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',CELLTYPE,'subjID',{MONKEY{1}(1)});
        for a = 1:length(filenames)
            stro = {};
            for j = 1:length(filenames{a})
                stro{j} = nex2stro(char(findfile(filenames{a}(j), fullfile(nexfilepath,'Greg',MONKEY))));
            end
            stro = strocat(stro);
            ntrials = size(stro.trial,1);
            Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_l'));
            Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_m'));
            TF = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
            stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
            stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
            dur = mean(stimoff_t-stimon_t);
            
            Lblank = Lcc == 0 & Mcc == 0 & TF == 0;
            tmpfftdata = [];
            for i=find(Lblank)'
                tmp = stro.ras{i,spikecds(a)}-stimon_t(i);
                tmp(tmp<0) = [];
                tmp(tmp>dur) = [];
                n = hist(tmp,timebins);
                
                tmpfftdata = [tmpfftdata; fftshift(fft(n))];
                isis = diff(timebins(logical(n)));
                isis(isis>ISIbins(end)) = [];
                n = hist(isis,ISIbins);
                ISIdata(strcmp(MONKEY,MONKEYS),strcmp(CELLTYPE,CELLTYPES),:) = ...
                    squeeze(ISIdata(strcmp(MONKEY,MONKEYS),strcmp(CELLTYPE,CELLTYPES),:)) + n';
            end
            fftdata = [fftdata; mean(abs(tmpfftdata).^2*2/length(tmpfftdata))];
        end
        axes(hax(strcmp(MONKEY,MONKEYS),strcmp(CELLTYPE,CELLTYPES)));
        
        freqs = linspace(-nyquist,nyquist,size(fftdata,2));
        plot(freqs,10*log10(mean(fftdata)),'k-','LineWidth',2);
        set(gca,'Ylim',[-11 -5],'Xlim',[0 100]);
        xlabel('Frequency (Hz)'); ylabel('Power (dB)');
        title([MONKEY{1},' ',CELLTYPE{1}]);
    end
end

% Now looking at the ISIs
figprefs;
hax(1,1) = axes('Position',[2 2 7 7]);
hax(1,2) = axes('Position',[2 12 7 7]);
hax(2,1) = axes('Position',[12 2 7 7]);
hax(2,2) = axes('Position',[12 12 7 7]);
for MONKEY = MONKEYS
    for CELLTYPE = CELLTYPES
        axes(hax(strcmp(MONKEY,MONKEYS),strcmp(CELLTYPE,CELLTYPES)));
        totISIs = sum(squeeze(ISIdata(strcmp(MONKEY,MONKEYS),strcmp(CELLTYPE,CELLTYPES),:)));
        h = bar(ISIbins,squeeze(ISIdata(strcmp(MONKEY,MONKEYS),strcmp(CELLTYPE,CELLTYPES),:)));
        set(h,'Edgecolor','black','FaceColor','black');
        set(gca,'Xlim',[ISIbins(1) ISIbins(end)]);
    end
end

cvs = [];
meanrates = [];
for i = 1:2 % MONKEY
    for j = 1:2 % CELLTYPE
        empirical_distn = squeeze(ISIdata(i,j,:))./sum(squeeze(ISIdata(i,j,:)));
        mn = ISIbins*empirical_distn;
        meanrates(i,j) = 1./mn;
        v = ISIbins.^2*empirical_distn-mn^2;
        cvs(i,j) = sqrt(v)/mn;
    end
end
%%
% Section 10.1 
% Periodogram estimate of noise spectra (so I can get better SNR for
% frequencies that I care about, 1-60 Hz).

timebins = [0:.005:.660];
deltaT = timebins(2)-timebins(1);
nyquist = 1./(2*deltaT);

figprefs;
hax(1,1) = axes('Position',[2 2 7 7]);
hax(1,2) = axes('Position',[2 12 7 7]);
hax(2,1) = axes('Position',[12 2 7 7]);
hax(2,2) = axes('Position',[12 12 7 7]);
for MONKEY = MONKEYS
    for CELLTYPE = CELLTYPES
        fftdata =[];
        [filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',CELLTYPE,'subjID',{MONKEY{1}(1)});
        for a = 1:length(filenames)
            stro = {};
            for j = 1:length(filenames{a})
                stro{j} = nex2stro(char(findfile(filenames{a}(j), fullfile(nexfilepath,'Greg',MONKEY))));
            end
            stro = strocat(stro);
            ntrials = size(stro.trial,1);
            Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_l'));
            Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_m'));
            TF = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
            stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
            stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
            dur = mean(stimoff_t-stimon_t);
            
            Lblank = Lcc == 0 & Mcc == 0 & TF == 0;
            tmpfftdata = [];
            for i=find(Lblank)'
                tmp = stro.ras{i,spikecds(a)}-stimon_t(i);
                tmp(tmp<0) = [];
                tmp(tmp>dur) = [];
                n = hist(tmp,timebins);
                tmpfftdata = [tmpfftdata; fftshift(fft(n))];
            end
            fftdata = [fftdata; mean(abs(tmpfftdata).^2*2/length(tmpfftdata))];
        end
        axes(hax(strcmp(MONKEY,MONKEYS),strcmp(CELLTYPE,CELLTYPES)));
        
        freqs = linspace(-nyquist,nyquist,size(fftdata,2));
        plot(freqs,10*log10(mean(fftdata)),'k-','LineWidth',2);
        set(gca,'Ylim',[-11 -5],'Xlim',[0 60]);
        xlabel('Frequency (Hz)'); ylabel('Power (dB)');
        title([MONKEY{1},' ',CELLTYPE{1}]);
    end
end


%%
% Section 11
% Psychophysical temporal integration

stimtype = 'lum'; % 'lum' or 'chr'
stimdurs = [333 666];
tmp_tfs = logspace(log10(1),log10(20),100);

figprefs;
hax = [];
hax(1) = axes('Position',[4 14 7 7]); hold on;
hax(2) = axes('Position',[4 3 7 7]); hold on;
for MONKEY = MONKEYS
    for stimdur = stimdurs
        conn = database('Nex_Paradigm_Sort','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
        comment = ['std inputs ',num2str(stimdur),' ',stimtype];
        query = sprintf('SELECT fileID FROM LMTF WHERE subjID = ''%s'' AND rfX = 50 AND rfY = 0 AND quality = 1 AND notes LIKE (''%s'')',MONKEY{1}(1),comment);
        flist = fetch(conn, query);
        close(conn)
        rawdata = getLMTFrawdata(flist);
        data = [sqrt(rawdata(:,1).^2+rawdata(:,2).^2) rawdata(:,3)];
        % Loading the model
        isosamppath = which('IsoSampOnline');
        isosamppath(find(isosamppath==filesep,1,'last'):end) = [];
        load ([isosamppath,filesep,'private',filesep','data',filesep,'LMTF.mat']);
        rfx = unique(rawdata(:,5))/10;
        rfy = unique(rawdata(:,6))/10;
        model = LMTF_global_to_local_model(getfield(getfield(eval(MONKEY{1}(1)),'legacy'),'mode5params'), rfx, rfy, 5);
        if strcmp(stimtype,'chr')
            pred = @(omega)1./(model(1+6)*abs(((1i*2*pi*10^model(5+6).*omega+1).^-model(3+6))-model(2+6)*((1i*2*pi*10^(model(5+6)+model(6+6)).*omega+1).^-(model(3+6)+model(4+6)))));
        else
            pred = @(omega)1./(model(1)*abs(((1i*2*pi*10^model(5).*omega+1).^-model(3))-model(2)*((1i*2*pi*10^(model(5)+model(6)).*omega+1).^-(model(3)+model(4)))));
        end
        axes(hax(strcmp(MONKEY, MONKEYS)));
        plot(tmp_tfs,pred(tmp_tfs),'k-');
        plot(tmp_tfs,pred(tmp_tfs)*sqrt(2),'r-');
        TFs = unique(data(:,2),'sorted');
        colors = [1 0 0; 0 0 0];
        for i = 1:length(TFs)
            j = find(stimdur == stimdurs);
            L = data(:,2) == TFs(i);
            plot(TFs(i)+.02*(j-1),median(data(L,1)),'s','MarkerFaceColor',colors(j,:),'MarkerEdgeColor',colors(j,:),'MarkerSize',10);
            plot(TFs(i)*[1 1]+.02*(j-1),[prctile(data(L,1),25) prctile(data(L,1),75)],'-','Color',colors(j,:),'Linewidth',1);
            plot(TFs(i)*[.95 1.05]+.02*(j-1),[prctile(data(L,1),25) prctile(data(L,1),25)],'-','Color',colors(j,:),'Linewidth',1);
            plot(TFs(i)*[.95 1.05]+.02*(j-1),[prctile(data(L,1),75) prctile(data(L,1),75)],'-','Color',colors(j,:),'Linewidth',1);
        end
    end
    set(gca,'Yscale','log','Xscale','log');
    
    title(MONKEY);
    set(gca,'Xlim',[TFs(1)*.9 TFs(end)*1.1]);
    if strcmp(stimtype, 'lum')
        set(gca,'Ylim',[.04 .24],'Ytick',[.05 .1 .2],'Xtick',[1 5 20]);
    else
        set(gca,'Ylim',[.02 .2]);
    end
    xlabel('Frequency (Hz)');
    ylabel('Threshold (contrast)'); % Is this OK?
end

%%
% Section 12
% sliding window d' and firing rate analysis.

nwindows = 9; % number of spike counting windows
signaloffsets = [zeros(nwindows,1) linspace(.04,.66,nwindows)'] % stretching window
MAXDPRIME = 10;
MAXFR = 60;

figprefs;
hax = [];
hax(1,1) = axes('position',[3 12 5 7]); hold on;
hax(1,2) = axes('position',[3 3 5 7]); hold on;
hax(2,1) = axes('position',[10 12 5 7]); hold on;
hax(2,2) = axes('position',[10 3 5 7]); hold on;
hax(3,1) = axes('position',[16 12 1 7],'YAxisLocation','right');
hax(3,2) = axes('position',[16 3 1 7],'YAxisLocation','right');
for MONKEY = MONKEYS
    [filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',{'M'},'subjID',{MONKEY{1}(1)});
    data_dprime = nan*ones(length(TFBINCENTERS), nwindows, length(filenames)); % TF, window, cell
    data_sp_per_s = nan*ones(length(TFBINCENTERS), nwindows, length(filenames));
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
        pop_scale_fact = IsoSampGetPopulationScaleFactor(stro,ecc_to_diam_deg_M, TEMPORONASALSCALEFACTOR, RFTRUNCATIONINSD, ONOFFCORRELATION, 2);
        for j = 1:size(signaloffsets,1)
            [uniquestim, dprime] = IsoSampGetDPrime(stro, DPRIMEMETHOD, spikecds(i), signaloffsets(j,:));
            Lblank = uniquestim(:,1) == 0 & uniquestim(:,2) == 0 & uniquestim(:,3) == 0;
            Llum = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & ~Lblank;
            uniquestim = uniquestim(Llum,:);
            dprime = dprime(Llum);
            
            for k = 1:size(uniquestim,1)
                LTF = softEq(uniquestim(k,3),TFBINCENTERS,5);
                if sum(LTF) ~= 1
                    error('Trouble finding TF');
                end
                data_dprime(k,j,i) = dprime(LTF)*pop_scale_fact;
                L = Lcc == uniquestim(k,1) & Mcc == uniquestim(k,2) & TF == uniquestim(k,3);
                spt = spiketimes(L);
                ston_t = stimon_t(L);
                spikevect = [];
                for l = 1:length(ston_t)
                    spikevect = [spikevect; spt{l}-ston_t(l)];
                end
                nspikes = sum(spikevect>signaloffsets(j,1) & spikevect < signaloffsets(j,2));
                data_sp_per_s(k,j,i) = nspikes./sum(L)/diff(signaloffsets(j,:));
            end
        end
    end
    Lallnan = all(all(isnan(data_dprime),2),3);
    axes(hax(strcmp(MONKEY,MONKEYS),1));
    tmp = nanmean(data_dprime(~Lallnan,:,:),3);
    image(tmp*255/MAXDPRIME,'CDataMapping','direct');
    if any(tmp(:) > MAXDPRIME)
        error('increase MAXDPRIME');
    end
    contour(nanmean(data_dprime(~Lallnan,:,:),3),[1.27 1.27],'w-','LineWidth',2); title([MONKEY{1},' SNR']);
    set(gca,'Ytick',[1:2:length(TFBINCENTERS)]','Yticklabel',num2str(TFBINCENTERS(1:2:end)',2));
    set(gca,'Xtick',[1:2:size(signaloffsets,1)],'Xticklabel',num2str(signaloffsets(1:2:end,2),1));
  
    axes(hax(strcmp(MONKEY,MONKEYS),2));
    tmp = nanmean(data_sp_per_s(~Lallnan,:,:),3);
    image(tmp*255/MAXFR,'CDataMapping','scaled');
    if any(tmp(:) > MAXFR)
        error('increase MAXFR');
    end
    if strcmp(MONKEY,'Apollo')
        contour(nanmean(data_sp_per_s(~Lallnan,:,:),3),[17 17],'w-','LineWidth',2);
    else
        contour(nanmean(data_sp_per_s(~Lallnan,:,:),3),[16 16],'w-','LineWidth',2);
    end
    title('Firing rate');
    xlabel('Integration time (s)');
    ylabel('TF (Hz)');
    set(gca,'Ytick',[1:2:length(TFBINCENTERS)]','Yticklabel',num2str(TFBINCENTERS(1:2:end)',2));
    set(gca,'Xtick',[1:2:size(signaloffsets,1)],'Xticklabel',num2str(signaloffsets(1:2:end,2),2));
end
colormap(gray(255));

% making separate colorbars
axes(hax(3,1));
image([1:255]');
set(gca,'Ylim',[0 255])
yticks = 0:2:MAXDPRIME;
set(gca,'XTick',[],'Ytick',yticks*255/MAXDPRIME,'Yticklabel',yticks);
ylabel('Sensitivity (d'')');

axes(hax(3,2));
image([1:255]');
set(gca,'Ylim',[0 255])
yticks = 0:10:MAXFR;
set(gca,'XTick',[],'Ytick',yticks*255/MAXFR,'Yticklabel',yticks);
ylabel('Firing rate (sp/s)');

%%
% Section 13
% Standard deviations of "normalized" Mahalanobis distances

titles = {'Signal mean','Signal SD','Noise mean','Noise SD'};
meanbins = linspace(-3, 8, 25);
sdbins = linspace(0, 3, 25);
for MONKEY = MONKEYS
    for CELLTYPE = CELLTYPES
        data = [];
        [filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',CELLTYPE,'subjID',{MONKEY{1}(1)});
        for i = 1:length(filenames)
            stro = {};
            for j = 1:length(filenames{i})
                stro{j} = nex2stro(char(findfile(filenames{i}(j), fullfile(nexfilepath,'Greg',MONKEY))));
            end
            stro = strocat(stro);

            [uniquestim, dprime, signal, noise] = IsoSampGetDPrime(stro,DPRIMEMETHOD,spikecds(i));
            Lblank = uniquestim(:,1) == 0 & uniquestim(:,2) == 0 & uniquestim(:,3) == 0;
            Llum = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & ~Lblank;
            for j = find(Llum)'
                data = [data; i uniquestim(j,3) nanmean(signal{j}) nanstd(signal{j}) nanmean(noise{j}) nanstd(noise{j})];
            end
        end
        figure;
        for j = 1:4
            subplot(2,2,j);
            if j/2 == floor(j/2)
                bins = sdbins;
            else
                bins = meanbins;
            end
            hist(data(:,j+2),bins);
            title(titles{j});
            set(gca,'Xlim',[bins(1) bins(end)]);
        end
        set(gcf,'Name',[char(MONKEY),' ',char(CELLTYPE)]);
        prettycorr(data(:,2:end),{'TF','Signal mu','Signal SD','Noise mu','Noise SD'});
        set(gcf,'Name',[char(MONKEY),' ',char(CELLTYPE)]);
    end
end

%%
% Section 14
% Comparing d' for magnocellular neurons between 1 and 8 Hz.
% No plot, just stats

for MONKEY = MONKEYS
    for CELLTYPE = CELLTYPES
        data = [];
        [filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',CELLTYPE,'subjID',{MONKEY{1}(1)});
        for i = 1:length(filenames)
            stro = {};
            for j = 1:length(filenames{i})
                stro{j} = nex2stro(char(findfile(filenames{i}(j), fullfile(nexfilepath,'Greg',MONKEY))));
            end
            stro = strocat(stro);
            [uniquestim, dprime] = IsoSampGetDPrime(stro,DPRIMEMETHOD,spikeidx);
            Lblank = uniquestim(:,1) == 0 & uniquestim(:,2) == 0 & uniquestim(:,3) == 0;
            Llum = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & ~Lblank;
            LTF_lo = uniquestim(:,3) == 1;
            LTF_hi = uniquestim(:,3) > 7 & uniquestim(:,3) < 8;
            
            data = [data; find(strcmp(MONKEY,MONKEYS)) dprime(Llum&LTF_lo) dprime(Llum&LTF_hi)];
        end
    end
end

for i = 1:2
    L = data(:,1) == i;
    p = signrank(data(L,2),data(L,3))
    disp([char(MONKEYS{i}),': mean low: ',num2str(mean(data(L,2))),' mean high: ',num2str(mean(data(L,3))),' p = ',num2str(p)]);
end

%%
% Section 15
% Comparing d' to non-parametric auROC values

data = [];
for MONKEY = MONKEYS
    for CELLTYPE = CELLTYPES       
        [filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',CELLTYPE,'subjID',{MONKEY{1}(1)});
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
            for j = find(Lrg)'
                n = noise{j}(~isnan(noise{j}));
                s = signal{j}(~isnan(signal{j}));
                tmp = [tmp; find(strcmp(MONKEY,MONKEYS)) find(strcmp(CELLTYPE, CELLTYPES)) dprime(j) roc(n,s)];
            end
            data =[data; tmp];
        end
    end
end

figprefs;
axes('Position',[3 3 8 8]);
for MONKEY = MONKEYS
    for CELLTYPE = CELLTYPES
        L = data(:,1) == find(strcmp(MONKEY,MONKEYS));
        L = L & data(:,2) == find(strcmp(CELLTYPE,CELLTYPES));
        sum(L)
        h = plot(data(L,3),data(L,4),'ko','MarkerSize',4,'MarkerFaceColor','none');
        
        if find(strcmp(MONKEY,MONKEYS)) == 2
           set(h,'Marker',MONKEY2MARKER); 
        end
        if find(strcmp(CELLTYPE,CELLTYPES)) == 2
            set(h,'MarkerFaceColor','none','MarkerEdgeColor',[.5 .5 .5]);
        end
    end
end


xlabel('Sensitivity (d'')'); ylabel('non-parametric auROC');
% Plot theoretical curve
mus = linspace(-1,5,100);
ROCareas = zeros(size(mus));
for i = 1:length(mus)
    x = linspace(-6,mus(i)+6,1000);
    ROCareas(i) = normcdf(x,0,1)*(normpdf(x,mus(i),1)./sum(normpdf(x,mus(i),1)))';
end
plot(mus,ROCareas,'k-');
set(gca,'Xlim',[-2 6],'Ylim',[0.2 1])

%%
% Section 16 Power spectrum of the 1 Hz stimulus
t = [0:1/240:.666];
nframes = length(t);
ramplength = nframes/4;
f = 1; % Hz
stim = sin(2*pi*t*f).*[linspace(0,1,ramplength) ones(1, 2*ramplength) linspace(1,0,ramplength)];
plot(t,stim);
[Pxx,F] = periodogram(stim,[],[],240);
plot(F,Pxx);

%%
% Section 17
% Collecting stro structures for a GitHub repository

data = {{} {}; {} {}};
for MONKEY = MONKEYS
    for CELLTYPE = CELLTYPES
        [filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',CELLTYPE,'subjID',{MONKEY{1}(1)});
        for i = 1:length(filenames)
            stro = {};
            for j = 1:length(filenames{i})
                stro{j} = nex2stro(char(findfile(filenames{i}(j))));
            end
            stro = strocat(stro);
            listsofar = data{strcmp(MONKEY,MONKEYS),strcmp(CELLTYPE,CELLTYPES)};
            listsofar{length(listsofar)+1} = stro;
            data{strcmp(MONKEY,MONKEYS),strcmp(CELLTYPE,CELLTYPES)}=listsofar;
        end
    end
end

%%
% Section 18 How large were Jiang et al's Gaussians? (same peak and "volume"
% of a 1� radius disk stimulus). Doing this by brute force trial-and-error.

[x,y] = meshgrid(linspace(-3,3,1000));
disk = x.^2+y.^2 <= 1;
disk_energy = sum(disk(:))
sigma = 1/sqrt(2);
gaussian = normpdf(x,0,sigma).*normpdf(y,0,sigma);
gaussian = gaussian./max(gaussian(:)); % constraint 1: identical peak

gaussian_energy = sum(gaussian(:))

%%
% Section 19) Trying to convert a 14% contrast threshold prediction error
% rate to d-primes. Have to assume a slope of a psychometric function

x = logspace(-4,0,100);
alpha = .2;
beta = 3; % This is the key parameter that I don't know 2 is conservative?
% Assuming a greater slope assumes more change in percent correct with
% contrast.
y = 1-0.5*exp(-1*(x/alpha).^beta);
interp1(x,y,alpha)

plot(x,y);
set(gca,'Xscale','log')

lims = alpha*[1-.14 1+.14];
pct_correct_lims = interp1(x,y,lims);
% Does a 14% change in threshold map to a consistent change in percent
% correct? Yes. 

% +14% increase in threshold raises percent correct by 2.4 (81.6 to 84)
% -14% decrease in threshold lowers percent correct by 2.75 (81.6 to 79)

% Now I just need to find d's that correspond to these limits
dprime = 1.27;
1-normcdf(0,dprime/sqrt(2),1) % See StatsStuff Section 11

% Making a function that converts dprimes to percent correct in 2AFC
dprimes = linspace(.75, 2,100)
prct_correct = 1-normcdf(0,dprimes./sqrt(2),1);
dprime_lims = interp1(prct_correct,dprimes,pct_correct_lims);
% [0.89 1.27 1.70] <-- theoreticals � SE

%%
% Section 20) Cone weights for M, P, and K cells
% Option to try non-standard cone fundamentals

USESTANDARDFUNDAMENTALS = true;

[WNfilenames, ~, WNspikenames, WNneuronids] = fnamesFromTxt('WhiteNoiseLGN_forIS','cellClass',{'M','P','K'});
[ISfilenames, ~, ISspikenames, ISneuronids] = fnamesFromTxt('IsoSamp_LGN','cellClass',{'M','P','K'});

if ~USESTANDARDFUNDAMENTALS
    if exist('absorptance') % From IsoSampPaperStuff2, Section 1.1
        funds = absorptance';
        disp('USING "ABSORPTANCE" fundamentals');
    else
        funds = load('T_cones_synthgh4');
        funds = funds.T_cones_synthgh4;
        if size(funds,1) < size(funds,2)
            funds = funds';
        end
        disp('USING "T_CONES_SYNTHGH4" fundamentals');

        %    funds = load('T_cones_ss10');
        %    funds = SplineRaw(S_cones_ss10,funds.T_cones_ss10',[380 5 81]);
    end
else
    disp('USING STANDARD CONE FUNDAMENTALS (SMJ 10�)');
end
 
strippeddownWNneuronids = [];
for neuronid = WNneuronids'
    LWN =  WNneuronids == neuronid;
    Lisosamp = ISneuronids == neuronid;
    if any(Lisosamp)
        strippeddownWNneuronids = [strippeddownWNneuronids;neuronid];
    end
end

% Getting manually curated cell types
conn = database('Nex_Paradigm_Sort','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
cellTypes = [];
for i = 1:length(strippeddownWNneuronids)
    class = fetch(conn, ['SELECT cellClass FROM WhiteNoiseLGN_forIS WHERE quality = ''1'' AND neuron = ',num2str(strippeddownWNneuronids(i))]);
    if istable(class)
        class = table2cell(class);
    end
    cellTypes{i,1} = cell2mat(unique(class));
end
close(conn)

data = [];
Ms = [];
for i = 1:length(strippeddownWNneuronids)
    idx = WNneuronids == strippeddownWNneuronids(i); % index into WNneuronids and WNspikenames
    stro=[];
    fnames = WNfilenames{idx};
    for j = 1:length(fnames)
        stro = strocat(stro,nex2stro(findfile(fnames{j})));
    end
    spikename = WNspikenames(idx,:);

    nstixperside = stro.sum.exptParams.nstixperside;
    ntrials = length(stro.sum.absTrialNum);
    spikeidx = find(strcmp(stro.sum.rasterCells(1,:),spikename));
    
    % Reconstructing the M matrix and gamma table
    if USESTANDARDFUNDAMENTALS % Use the fundamentals used to create the stimulus
        funds = stro.sum.exptParams.fundamentals;
        funds = reshape(funds,[length(funds)/3,3]);
    end
    mon_spd = stro.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = funds'*mon_spd;
    
    ridx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_r'));
    gidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_g'));
    bidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_b'));
    % Values in stro.trial are gun *intensities*, not voltages. 
    % See IsoSamp.d for request and LMTFSlave.m for definition ofgl.bkgndrgb 3/29/20 GDLH
    bkgndrgb = [mode(stro.trial(:,ridx)), mode(stro.trial(:,gidx)), mode(stro.trial(:,bidx))]';

    bkgndlms = M*bkgndrgb; % this is in cone excitations
    Mrgbtocc = diag(1./bkgndlms)*M; % Weight the fundamentals, not the guns (i.e. premult. by diag...)
    Ms(:,:,i) = Mrgbtocc; % Assumes a weighted sum of cone contrasts
    out = getWhtnsStats(stro,MAXT,'STCOVmex', {nstixperside^2, 3, MAXT}, spikename);
    STA = out{1}/out{3};
    inRF = getRFfromSTA(STA,GETRFFROMSTAMODE,GETRFFROMSTATHRESH);
    STA = reshape(STA,[nstixperside.^2  3 MAXT]);
    temporalSTA = squeeze(STA(logical(inRF(:)),:,:));
    if ~ismatrix(temporalSTA)
        temporalSTA = mean(temporalSTA,1);
    end
    %temporalSTA = squeeze(STA(whichpix,:,:));
    
    data(:,:,i) = squeeze(temporalSTA);
end

% Getting cone weights
svdvars = zeros(size(data,3),2);
coneweights = [];
for i = 1:length(strippeddownWNneuronids)
    sta = squeeze(data(:,:,i));
    [u,s,v] = svd(sta);
    svdvars(i,:) = [s(1) sum(diag(s))]; 
    rgb = u(:,1);
    if rgb'*sum(sta,2) < 0
        rgb = -rgb;
    end
    M = Ms(:,:,i);
    coneweights(i,:) = rgb'*inv(M);
end

%figure; axes; plot3(coneweights(:,1),coneweights(:,2),coneweights(:,3),'ko');
normalizedconeweights = coneweights./repmat(sum(abs(coneweights),2), 1, size(coneweights,2));
colors = [0 .5 1; 0 0 0; 1 0 0];
uniquecellTypes = unique(cellTypes);

figprefs;
axes('position',[4.5,8,8,8]); hold on;
for i = 1:length(uniquecellTypes)
    L = strcmp(cellTypes,uniquecellTypes(i));
    for j = find(L')
        h = plot(normalizedconeweights(j,1), normalizedconeweights(j,2),'o','MarkerEdgeColor',colors(i,:),'MarkerFaceColor','none');
        set(h,'MarkerSize',5);
        if normalizedconeweights(j,3) > 0
           set(h,'MarkerFaceColor',colors(i,:),'MarkerEdgeColor','white','MarkerSize',7,'LineWidth',.25);
        end
    end
end
set(gca,'Xlim',[-1 1],'Ylim',[-1 1],'Xtick',[-1 0 1],'Ytick',[-1 0 1]);
plot([0 1 0 -1 0],[-1 0 1 0 -1],'k-');
plot([-1 1],[0 0],'k:');
plot([0 0],[-1 1],'k:');

title('open symbols mean S-')
xlabel('normalized L-cone weight');
ylabel('normalized M-cone weight');
set(gca,'box','off');

% I'm getting a lot of L-M-S parvocells. 
% Also a lot of S-opponency in the magnocells.
% This sounds like a calibration problem. 
% Can I fix it with my custom macaque cone fundamentals? Yes

% Number of magno cells for which S goes with L+M
L = strcmp(cellTypes,{'M'});
a = sum(sign(sum(normalizedconeweights(L,[1 2]),2)) == sign(normalizedconeweights(L,3)))
binocdf(a, sum(L),.5)
binocdf(sum(L)-a, sum(L),.5)

%%
% Section 21
% Looking at L+M and L-M sensitivity for M, P, and K populations.
% Also asking (for parvocellular neurons) whether the cone weights 
% predict d' (L+M and L-M). 
% Also getting eccentricity.
% Need to run Section 20 first.

data = [];
for i = 1:length(strippeddownWNneuronids)
    idx = strippeddownWNneuronids(i);
    filenames = ISfilenames{ISneuronids == idx};
    spikenames = ISspikenames(ISneuronids == idx,:);
     
    stro=[];
    for j = 1:length(filenames)
        stro = strocat(stro,nex2stro(findfile(filenames{j})));
    end
    [uniquestim, dprime] = IsoSampGetDPrime(stro,DPRIMEMETHOD,abs(spikenames(end)-96));
    Lblank = uniquestim(:,1) == 0 & uniquestim(:,2) == 0 & uniquestim(:,3) == 0;
    Llum = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & ~Lblank;
    Lrg = sign(uniquestim(:,1)) ~= sign(uniquestim(:,2)) & ~Lblank;
    
    % For starters, let's average across TFs
    data = [data; mean(dprime(Llum)) mean(dprime(Lrg)) stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y];
end

% Sensitivity to L+M and L-M by cell type
colors = [0 .5 1; 0 0 0; 1 0 0];
uniquecellTypes = unique(cellTypes);
AXWIDTH = 6;
figprefs;
axes('Position',[12 2 AXWIDTH AXWIDTH])
hold on;
for i = 1:length(uniquecellTypes)
    L = strcmp(cellTypes,uniquecellTypes(i));
    plot(data(L,1), data(L,2),'wo','MarkerFaceColor',colors(i,:)','LineWidth',.1)
end
xlabel('Signal-to-noise ratio (d'') L+M');
ylabel('Signal-to-noise ratio (d'') L-M');
plot([0 max(data(:))], [0 max(data(:))],'k-');
bnds = [min(min(data(:,[1,2]))) max(max(data(:,[1,2])))];
set(gca,'Xlim',[-.1 .1]+bnds,'Ylim',[-.1 .1]+bnds)
set(gca,'XTick',[0:2],'YTick',[0:2]);

% Sensitivity by cone weights
axes('Position',[3 2 AXWIDTH AXWIDTH])
hold on;
for whichcelltype = {'M','P','K'}
    L = strcmp(cellTypes,whichcelltype);
    LM = abs(normalizedconeweights(L,1)+normalizedconeweights(L,2));
    plot(LM,data(L,2)-data(L,1),'wo','MarkerFacecolor',colors(strcmp(whichcelltype, uniquecellTypes),:),'LineWidth',.1);
    xlabel('L+M cone weights from white noise');
    ylabel('d''_{L-M} - d''_{L+M} ')
    
    % Regression
    [b,bint, ~, ~, stats] = regress(data(L,2)-data(L,1), [ones(size(LM,1),1) LM]);
    disp([whichcelltype{1},' slope: ',num2str(b(2)),' p:',num2str(stats(3))]);
end
set(gca,'Ylim',[-2 2]);

% Eccentricity (no effect on cone weights or difference in d's)
figure; 
subplot(2,1,1); hold on;
for whichcelltype = {'M','P','K'}
    L = strcmp(cellTypes,whichcelltype);
    LM = abs(normalizedconeweights(L,1)+normalizedconeweights(L,2));
    plot3(data(L,3)/10, data(L,4)/10, LM,'wo','markerfacecolor',colors(strcmp(whichcelltype, uniquecellTypes),:),'MarkerSize',8);
end
xlabel('Horz. (DVA)'); ylabel('Vert. (DVA)'); zlabel('L+M WN');
subplot(2,1,2); hold on;
for whichcelltype = {'M','P','K'}
    L = strcmp(cellTypes,whichcelltype);
    LM = data(L,1)-data(L,2); % IsoSamp d's
    plot3(data(L,3)/10, data(L,4)/10, LM,'wo','markerfacecolor',colors(strcmp(whichcelltype, uniquecellTypes),:),'MarkerSize',8);
end
xlabel('Horz. (DVA)'); ylabel('Vert. (DVA)'); zlabel('d''_{L+M}-d''_{L-M}');

% 2-D plots
figure;
subplot(2,1,1); hold on;
for whichcelltype = {'M','P','K'}
    L = strcmp(cellTypes,whichcelltype);
    LM = abs(normalizedconeweights(L,1)+normalizedconeweights(L,2));
    plot(sqrt((data(L,3)/10).^2+(data(L,4)/10).^2), LM,'wo','markerfacecolor',colors(strcmp(whichcelltype, uniquecellTypes),:),'MarkerSize',8);
end
xlabel('ecc. (DVA)'); ylabel('L+M WN');

subplot(2,1,2); hold on;
for whichcelltype = {'M','P','K'}
    L = strcmp(cellTypes,whichcelltype);
    LM = data(L,1)-data(L,2); % IsoSamp d's
    plot(sqrt((data(L,3)/10).^2+(data(L,4)/10).^2), LM,'wo','markerfacecolor',colors(strcmp(whichcelltype, uniquecellTypes),:),'MarkerSize',8);
end
xlabel('ecc (DVA)'); ylabel('d''_{L+M}-d''_{L-M}');

% Whitenoise: no effect of eccentricity on cone weights
% IsoSamp: effect of eccentricity on d' is expected from fixed stimulus
% size (magno: high L+M d' at large eccentricities, parvo: high L-M d' at
% large eccentricities)
%%    
% Section 22
% Comparing d' for L+M and L-M as a function of TF

% Loading the data
filenames = [];
spikenames = [];
neuronids = [];
celltypemat = [];
celltypes = {'K','M','P'};
for i = 1:length(celltypes)
    [tmpfilenames, ~, tmpspikenames, tmpneuronids] = fnamesFromTxt('IsoSamp_LGN','cellClass',celltypes(i));
    filenames = [filenames; tmpfilenames];
    spikenames = [spikenames; tmpspikenames];
    neuronids = [neuronids; tmpneuronids];
    celltypemat = [celltypemat; repmat(i, size(tmpfilenames,1),1)];
end

data = [];
for i = 1:length(filenames)
    idx = neuronids(i);
    tmpfilenames = filenames{neuronids == idx};
    tmpspikenames = spikenames(neuronids == idx,:);

    stro=[];
    for j = 1:length(tmpfilenames)
        stro = strocat(stro,nex2stro(findfile(tmpfilenames{j})));
    end
    [uniquestim, dprime] = IsoSampGetDPrime(stro,DPRIMEMETHOD,abs(tmpspikenames(end)-96));
    Lblank = uniquestim(:,1) == 0 & uniquestim(:,2) == 0 & uniquestim(:,3) == 0;
    Llum = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & ~Lblank;
    Lrg = sign(uniquestim(:,1)) ~= sign(uniquestim(:,2)) & ~Lblank;

    tmp = nan*ones(2,length(TFBINCENTERS));
    for j = 1:length(TFBINCENTERS)
        LTF = uniquestim(:,3) >= TFBINEDGES(j) & uniquestim(:,3) < TFBINEDGES(j+1);
        tmp(:,j) = [mean(dprime(Llum&LTF));mean(dprime(Lrg&LTF))];
    end
    data(:,:,i) = tmp;
end
% data is [lum, rg] x TFBINCENTERS x neuron

% Plotting
colors = [0 .5 1; 0 0 0; 1 0 0];
figure;
for i = 1:length(TFBINCENTERS)
   subplot(4,4,i); hold on;
   tmp = squeeze(data(:,i,:))';
   for j = 1:max(celltypemat)
       L = celltypemat == j;
       h = plot(tmp(L,1),tmp(L,2),'wo','MarkerFaceColor',colors(j,:),'MarkerSize',6);
   end
   set(gca,'Xlim',[nanmin(data(:)) nanmax(data(:))],'Ylim',[nanmin(data(:)) nanmax(data(:))])
   axis square
   title(num2str(TFBINCENTERS(i)));
   xlabel('d''_{L+M}');
   ylabel('d''_{L-M}');
end

% Trajectory plot
figure; axes; hold on;
plot([0 3.5],[0 3.5],'k:')
plot([0 0],[-1 3.5],'k-')
plot([-1 3.5],[0 0],'k-')

for i = 1:max(celltypemat)
    L = celltypemat == i;
    mn = nanmean(data(:,:,L),3);
    sd = nanstd(data(:,:,L),[],3);
    plot(mn(1,:),mn(2,:),'-','color',colors(i,:),'linewidth',2);
    for j = 1:size(mn,2)
       plot(mn(1,j)+sd(1,j)*[-1 1],mn(2,j)*[1 1],'-','color',colors(i,:),'linewidth',0.5);
       plot(mn(1,j)*[1 1],mn(2,j)+sd(2,j)*[-1 1],'-','color',colors(i,:),'linewidth',0.5);
       plot(mn(1,j),mn(2,j),'o','MarkerSize',j+2,'MarkerEdgeColor','white','MarkerFaceColor',colors(i,:));
    end
end
axis square
xl=get(gca,'Xlim');
yl=get(gca,'Ylim');
set(gca,'Xlim',[min([xl(1), yl(1)]) max([xl(2), yl(2)])]);
set(gca,'Ylim',[min([xl(1), yl(1)]) max([xl(2), yl(2)])]);
set(gca,'Xlim',[-1 3.5],'Ylim',[-1 3.5]);
xlabel('d''_{L+M}');
ylabel('d''_{L-M}');

%%
% Section 23
% Power spectral analysis for L-M responses 
% (characterizing distortion in magnocells)

whichcellclass = 'P';
INDIVIDUALCELLPLOTS = false;

timebins = [0:.005:.660];
deltaT = timebins(2)-timebins(1);
nyquist = 1./(2*deltaT);
[filenames, ~, spikenames, neuronids] = fnamesFromTxt('IsoSamp_LGN','cellClass',{whichcellclass});

fftdata = nan*ones(2,length(TFBINCENTERS),length(filenames),length(timebins)); % colordir x tf x cell x hz

for i = 1:length(filenames)
    idx = neuronids(i);
    tmpfilenames = filenames{neuronids == idx};
    tmpspikename = spikenames(neuronids == idx,:);

    stro=[];
    for j = 1:length(tmpfilenames)
        stro = strocat(stro,nex2stro(findfile(tmpfilenames{j})));
    end
    
    ntrials = size(stro.trial,1);
    Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_l'));
    Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_m'));
    TF = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
    uniquestim = sortrows(unique([Lcc Mcc TF],'rows'),3);
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
    dur = mean(stimoff_t-stimon_t);
    
    tmpfftdata = nan*ones(ntrials,length(timebins));
    for j=1:ntrials
        spikeidx = strcmp(stro.sum.rasterCells,tmpspikename);
        tmp = stro.ras{j,spikeidx}-stimon_t(j);
        tmp(tmp<0) = [];
        tmp(tmp>dur) = [];
        n = hist(tmp,timebins);
        tmpfftdata(j,:) = fftshift(fft(n));
    end
    
    freqs = linspace(-nyquist, nyquist,size(tmpfftdata,2)); % Is this correct?
    Llum = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & uniquestim(:,3) > 0;
    Lrg = sign(uniquestim(:,1)) ~= sign(uniquestim(:,2)) & uniquestim(:,3) > 0;
 
    colors = hot(length(unique(TF)));
    if INDIVIDUALCELLPLOTS
        figure;
    end
    for j = 1:2 % L+M and L-M
        if INDIVIDUALCELLPLOTS
            subplot(2,1,j); hold on;
        end
        if j == 1
            whichconditions = find(Llum);
        else
            whichconditions = find(Lrg);
        end
        for k = 1:length(whichconditions)
            L = Lcc == uniquestim(whichconditions(k),1) & Mcc == uniquestim(whichconditions(k),2) & TF == uniquestim(whichconditions(k),3);
            amp = mean(abs(tmpfftdata(L,:)));
            if sum(softEq(uniquestim(whichconditions(j),3),TFBINCENTERS,3)) ~= 1
                disp('error')
                keyboard
            end
            fftdata(j,softEq(uniquestim(whichconditions(k),3),TFBINCENTERS,3),i,:) = amp; % colordir x tf x cell x hz
            if INDIVIDUALCELLPLOTS
                plot(freqs,amp,'-','color',colors(k,:));
            end
        end
    end
    subplot(2,1,1);
    title(stro.sum.fileName(find(stro.sum.fileName == filesep,1,'last')+1:end))
end

colors = jet(size(fftdata,2));
figure;
for i = 1:size(fftdata,1) % colordir
    subplot(2,1,i); hold on;
    for j = 1:size(fftdata,2) % stim TF
        if sum(~isnan(fftdata(i,j,:,1)),3) > 1 % Don't use conditions with only one observation
            tmp = squeeze(nanmean(fftdata(i,j,:,:),3));
            df = freqs(2)-freqs(1);
            plot(freqs,tmp.^2/(2*df),'color',colors(j,:));
            ylabel('Power spectral density');
            set(gca,'yscale','log');
        end
    end
end
subplot(2,1,1); title([whichcellclass,' cells: L+M']);
subplot(2,1,2); title([whichcellclass,' cells: L-M']);
xlabel('Frequency (Hz)');
%%
% Section 24: Relative SNR loss and different stages of the visual system
% for L-M modulation. Identical to section 6.3 (Figure 9) but for L-M
% modulations. Need to run section 6 first.

AXWIDTH = 7;
AXMARGIN = 2;
figprefs;
for colordiridx = 1:2 % LUM, RG
    if colordiridx == 1
        COLORDIR = 'LUM';
    else
        COLORDIR = 'RG';
    end
    
    summarydata1 = [];
    for j = 1:length(TFBINCENTERS)
        for celltype_idx = 1:2
            for monkey_idx = 1:2
                monkey_celltype_data = data{monkey_idx,celltype_idx};
                ncells = length(monkey_celltype_data);
                tmp = nan*ones(ncells,3); % ncells x [LGN, cone_current, photon absorption]
                for i = 1:ncells
                    cell_data = monkey_celltype_data{i};
                    Lblank = cell_data.uniquestim(:,1) == 0 & cell_data.uniquestim(:,2) == 0 & cell_data.uniquestim(:,3) == 0;
                    Llum = sign(cell_data.uniquestim(:,1)) == sign(cell_data.uniquestim(:,2)) & ~Lblank;
                    Lrg = sign(cell_data.uniquestim(:,1)) ~= sign(cell_data.uniquestim(:,2)) & ~Lblank;
                    if strcmp(COLORDIR, 'RG')
                        L = Lrg;
                    else
                        L = Llum;
                    end
                    Ltf = cell_data.uniquestim(:,3) > TFBINEDGES(j) & cell_data.uniquestim(:,3) <= TFBINEDGES(j+1);
                    if sum(Ltf&L) > 0
                        tmp(i,1) = nanmean(cell_data.dprime(L&Ltf))*cell_data.population_scalefactor;
                        tmp(i,2) = nanmean(cell_data.conedprime(L&Ltf));
                        tmp(i,3) = nanmean(cell_data.photondprime(L&Ltf));
                    end
                end
                if sum(~isnan(tmp(:,1))) <= 1 % Don't take TFs with only one data point
                    tmp = nan*ones(size(tmp));
                end
                summarydata1 = [summarydata1; celltype_idx TFBINCENTERS(j) nanmean(tmp)]; % Accumulating averages (across cells, within cell type, within color direction)
            end
        end
    end
    % Getting the mean cone ideal observer d's across cell types and monkeys
    L_enough_data = logical([]);
    summarydata2 = [];
    for whichtf = unique(summarydata1(:,2)')
        L = summarydata1(:,2) == whichtf;
        %sum(L)  
        summarydata2= [summarydata2; nanmean(summarydata1(L,4)) nanmean(summarydata1(L,5))];
        L_enough_data = [L_enough_data; ~all(isnan(summarydata1(L,3)))]
    end
    
    axes('position',[2+(colordiridx-1)*(AXWIDTH+AXMARGIN), 4, AXWIDTH, AXWIDTH]);
    fraction_left = summarydata2(L_enough_data,:) - 1.27; % 0 means behavioral sensitivity (d' = 1.27)
    fraction_left = fraction_left./repmat(fraction_left(L_enough_data,2),1,2);
    patch([TFBINCENTERS(L_enough_data),fliplr(TFBINCENTERS(L_enough_data))],[fraction_left(:,1)',zeros(1,size(fraction_left,1))],[.75 .75 .75]);
  %  patch([TFBINCENTERS(L_enough_data),fliplr(TFBINCENTERS(L_enough_data))],[fraction_left(:,2)',fliplr(fraction_left(:,1)')],[1 1 1]);
    
    % LGN d's
    for i = 1:2 % 1 = mean of M and P, 2 = dominant cell type only
        if i == 1
            L = ones(size(summarydata1,1),1);
        elseif strcmp(COLORDIR, 'LUM') % can assume i == 2
            L = summarydata1(:,1) == 1; % M cells only
        else
            L = summarydata1(:,1) == 2; % P cells only
        end
        tmp = [];
        for whichtf = unique(summarydata1(:,2)')
            tmp = [tmp; nanmean(summarydata1(L&summarydata1(:,2) == whichtf,3))];
        end
        tmp = (tmp-1.27)./(summarydata2(:,2) - 1.27);
        %tmp(tmp>0) = 0; % DON'T FORGET WE'RE ROUNDING LOW LGN D's UP TO
        %1.27 here! (Still doing it in plotting).
        h = patch([TFBINCENTERS(L_enough_data),fliplr(TFBINCENTERS(L_enough_data))],[tmp(L_enough_data)',zeros(1,sum(L_enough_data))],[.5 .5 .5]);
       % if i == 2
            set(h,'Facealpha',.5);
      %  end
    end

    % set(gca,'Ylim',[0 1],'Xscale','log','Xlim',[TFBINCENTERS(1), TFBINCENTERS(find(L_enough_data,1,'last'))]);
    set(gca,'Ylim',[0 1],'Xscale','log','Xlim',[TFBINCENTERS(1), TFBINCENTERS(end)],'box','on');
    set(gca,'Xticklabel',[1 10]);
    
    xlabel('Frequency (Hz)');
    ylabel('Fractional signal-to-noise ratio')
    set(gcf,'Renderer','painters');
    text(7,.8,'Phototransduction','FontName','Helvetica','FontAngle','italic');
    text(1.8,.35,'Cones \rightarrow LGN','FontName','Helvetica','FontAngle','italic');
   
    if strcmp(COLORDIR,'RG')
        titlestr = ['L-M'];
    else
        titlestr = ['L+M'];
    end
    title(titlestr);
end

% Scaling the X axes so that the plots are comparable
%xlims = get(gca,'Xlim');
%pos = get(gca,'Position');
%set(gca,'Position',[pos(1) pos(2) pos(3)+log(xlims(2)/TFBINCENTERS(end)) pos(4)]);
%%
% Section 25
% Effects of counting window on population sensitivity. The idea here is to
% show that there is no reasonable window within which P-cell population d' 
% is much greater that 1.27.
% Normalize scale across plots. Include contour at 1.27.

starts = linspace(0,.2,6);
stops = .666+starts;
tmp = fullfact([length(starts) length(stops)]);
offsets = [starts(tmp(:,1))' stops(tmp(:,2))'];

data = cell(2,2); % within every cell we have a little data matrix with Lcc, Mcc, TF, offset start, offset stop, and d'
for MONKEY = MONKEYS
    for CELLTYPE = CELLTYPES
        [filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',CELLTYPE,'subjID',{MONKEY{1}(1)});     
        for i = 1:length(filenames)
            stro = {};
            for j = 1:length(filenames{i})
                stro{j} = nex2stro(char(findfile(filenames{i}(j), fullfile(nexfilepath,'Greg',MONKEY))));
            end
            stro = strocat(stro);
            if strcmp(CELLTYPE,'M')
                ecc_to_diam_deg = ecc_to_diam_deg_M;
            else
                ecc_to_diam_deg = ecc_to_diam_deg_P;
            end

            uniquestim = IsoSampGetDPrime(stro,DPRIMEMETHOD,spikecds(i)); % Just getting uniquestim
            population_scalefactor = IsoSampGetPopulationScaleFactor(stro, ecc_to_diam_deg, TEMPORONASALSCALEFACTOR, RFTRUNCATIONINSD, ONOFFCORRELATION, 2);

            tmp = [uniquestim, nan*ones(size(uniquestim,1),size(offsets,1))];
            for j = 1:size(offsets,1)
                [uniquestim, dprime] = IsoSampGetDPrime(stro,DPRIMEMETHOD,spikecds(i),offsets(j,:));
                tmp(:,j+3) = dprime*population_scalefactor;
            end
            
             % Adding to the "data" cell array of cell arrays
            listsofar = data{strcmp(MONKEY,MONKEYS),strcmp(CELLTYPE,CELLTYPES)};
            listsofar{length(listsofar)+1} = tmp;
            data{strcmp(MONKEY,MONKEYS),strcmp(CELLTYPE,CELLTYPES)}=listsofar;
        end
    end
end

MONKEYIDX = 'both'; % 1 = Apollo, 2 = Utu
CELLTYPEIDX = 1; % 1 = magno, 2 = parvo
COLORDIR = 'RG';

if strcmp(MONKEYIDX,'both')
    tmp = cat(2,data{1,CELLTYPEIDX}, data{2,CELLTYPEIDX});
else
    tmp = data{MONKEYIDX,CELLTYPEIDX};
end

ncells = length(tmp);
outdata = nan*ones(length(TFBINCENTERS),length(starts),length(stops),ncells);
for cellcounter = 1:ncells
    uniquestim = tmp{cellcounter}(:,[1:3]);
    Lblank = uniquestim(:,1) == 0 & uniquestim(:,2) == 0 & uniquestim(:,3) == 0;
    if strcmp(COLORDIR,'RG')
        L = sign(uniquestim(:,1)) ~= sign(uniquestim(:,2)) & ~Lblank;
    else
        L = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & ~Lblank;
    end
    TFs = uniquestim(L,3);
    dprimes = tmp{cellcounter}(L,4:end);
    for i = 1:length(TFs)
        TFidx = find(TFs(i) > TFBINEDGES(1:end-1) & TFs(i) < TFBINEDGES(2:end));
        if size(dprimes,1) >= TFidx
            outdata(TFidx,:,:,cellcounter) = reshape(dprimes(TFidx,:),length(starts),length(stops));
        end
    end
end

ns = [];
for i = 1:size(outdata,1)
    ns(i,:,:) = sum(~isnan(outdata(i,:,:,:)),4);
end

minval = min(squeeze(nanmean(outdata(:,:,:,:),4)),[],'all');
maxval = max(squeeze(nanmean(outdata(:,:,:,:),4)),[],'all');
minval = 0;
maxval = 1.27;
AXWIDTH = 3;
AXMARGIN = .4;
figprefs;
for i = 1:size(outdata,1)
    mn = squeeze(nanmean(outdata(i,:,:,:),4));
    if all(all(ns(i,:,:) > 1)) & any(mn(:) < 1.27)
        axes('position',[2+(i-1)*(AXWIDTH+AXMARGIN),7,AXWIDTH,AXWIDTH],'Box','on');
        hold on;
        im = squeeze(reshape(mn,length(starts),length(stops)));
        image((im+minval)/(maxval-minval)*256); colormap(hot(256));
        [c,h] = contour(im,1.27*[1 2],'k-','LineWidth',2);
        if i == 1
            set(gca,'Ytick',1:length(starts),'YtickLabel',num2str(starts',2));
            set(gca,'Xtick',1:length(stops),'XtickLabel',num2str(stops',2));
            ylabel('start time (s)'); xlabel('end time (s)');
        else
            set(gca,'YTicklabel',[],'XTicklabel',[]);
        end
        title([num2str(round(TFBINCENTERS(i)*10)/10),' Hz']);
        set(gca,'XTickLabelRotation',90)
    end
end

% Colorbar
axes('position',[19,7,.3,AXWIDTH]);
set(gca,'Xtick',[],'YAxisLocation','right','Box','on')
set(gca,'Ytick',linspace(1,size(colormap,1),3),'YtickLabel',round(linspace(minval,maxval,3)*10)/10)
image([1:size(colormap,1)]');
ylabel('Signal-to-noise ratio (d'')');

if strcmp(MONKEYIDX,'both')
    monkeyname = 'Both monkeys';
else
    monkeyname = MONKEYS{MONKEYIDX};
end
set(gcf,'Name',[monkeyname,' ',CELLTYPES{CELLTYPEIDX},' ',COLORDIR]);

% starting the spike count window at 0.12 and ending at 0.746 looks reasonable