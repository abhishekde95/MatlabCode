%% Script for analysing the cone weights and the RF size of a neurons to see if there is any relationship between the 
% spatial frequency tuning and non-cardinal color tuning in isoluminant plane 
% List on files are present in NO BACKUP/NexFiles/nexfilelists/Greg/WhiteNoise/ConeVGun2.txt
% Author - Abhishek De, 11/16, Most of the code has been derived from ConeAnalysis.m written by GDLH
close all; clearvars;
% filename = 'K031308006.nex';
% filename = 'K040108003.nex';
filename = 'K102109004.nex';
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

% Reconstructing the M matrix
fundamentals = WN.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = WN.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

% Getting the background rgb/lms
ridx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_r'));
gidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_g'));
bidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_b'));
bkgndRGB = [mode(WN.trial(:,ridx)), mode(WN.trial(:,gidx)), mode(WN.trial(:,bidx))];
bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
bkgndlms = M*bkgndrgb;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stuff for looking at spike-triggered averages, etc.
spikename = getSpikenum(WN);
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
% Plotting the STA and PCs
if (~exist('mask'))
    mask = zeros(nstixperside, nstixperside);   % Pixel mask.  Someday make a tool to make this non-zero.
end
Lmask = logical(repmat(~mask(:),[3 1]));

% Normalizing images
nstixperside = WN.sum.exptParams.nstixperside;
rowidxs = reshape([1:3*nstixperside^2],[nstixperside^2 3]);
maxes = []; mins = [];
imagevectors = [STAs_cone];
for i = 1:3
    maxes = [maxes; max(max(imagevectors(rowidxs(:,i),:)))];
    mins = [mins; min(min(imagevectors(rowidxs(:,i),:)))];
end
potentialnormfactors = [(1-[.5; .5; .5]-eps)./(maxes-[.5; .5; .5]); (-[.5; .5; .5]+eps)./(mins-[.5; .5; .5])];
% 'eps' in above line is a kludge that is required for avoiding out of bounds errors.
potentialnormfactors(potentialnormfactors < 0) = []; % if min > mu or max < mu
normfactor = min(potentialnormfactors);
muvect = reshape(repmat([.5 .5 .5],nstixperside^2,1),nstixperside^2*3,1);

% Plotting the STA derived from cone noise
figure(1);
for i = 1:size(STAs_cone,2)
    STA = normfactor*(STAs_cone(:,i)-muvect)+muvect;
    STA = reshape(STA,[nstixperside nstixperside 3]);
    subplot(4,size(STAs_cone,2),i);
    image(STA);
    set(gca,'XTick',[],'YTick',[]); axis square;
    if (i == 1)
        ylabel('STA');
    end
end
% Cone significance maps from cone noise
figure(1);
rng = [min(STAs_cone(:)) max(STAs_cone(:))];
STAs_cone_z = STAs_cone./sqrt(nspikes_cone);
STAs_gun_z = STAs_gun./sqrt(nspikes_gun);
alpha = [0.01];  % max of 3 critical values
crit = norminv(1-alpha,0,1);
for i = 1:size(STAs_cone,2)
    for j = 1:3
        rowidxs = [(j-1)*nstixperside^2+1, j*nstixperside^2];
        subplot(4,size(STAs_cone,2),i+size(STAs_cone,2)*(j-1)+maxT);
        im = reshape(STAs_cone(rowidxs(1):rowidxs(2),i),[nstixperside nstixperside]);    
        im = (im-rng(1))/(rng(2)-rng(1));  % rescaling
        im = repmat(im, [1 1 3]);
        zim = reshape(STAs_cone_z(rowidxs(1):rowidxs(2),i),[nstixperside nstixperside]);
        for k = 1:length(alpha)
            pmat = logical(abs(zim)>crit(k));
            sigidxs = find(pmat);
            im(sigidxs+nstixperside^2*(k-1)) = .5;
        end
        image(im);
        axis image;
         if i==1 & j==1
            ylabel('L cone');
        elseif i==1 & j==2
            ylabel('M cone');
        elseif i==1 & j==3
            ylabel('S cone');
        end
        set(gca,'XTick',[],'YTick',[]);
    end
end

% Trying to get a RF automatically, pooling across cone noise. Using the first frame of the STA as the noise distribution
CHI2CRIT = .95; % For flaging a stixel as significant (adding gun and cone noise z-scores, squared)
s_cone = std(STAs_cone(:,1));
STAs_cone_z = STAs_cone./s_cone;
s_gun = std(STAs_gun(:,1));
STAs_gun_z = STAs_gun./s_gun;

% Spatial map
grandz = zeros([nstixperside nstixperside]);
maxzs = [];
include_gunnoise_data = 1;
if include_gunnoise_data
    channels = 6;
else
    channels = 3;
end
for i = 1:maxT
    tmp_cone = reshape(STAs_cone_z(:,i),[nstixperside nstixperside 3]);
    if include_gunnoise_data == 1
        tmp_gun = reshape(STAs_gun_z(:,i),[nstixperside nstixperside 3]); % This is the only place in the code where I use data from gun noise
        grandz = grandz+sum(tmp_gun.^2,3)+sum(tmp_cone.^2,3);
        maxzs(i) = sum(sum(sum(tmp_gun.^2,3)+sum(tmp_cone.^2,3)));
        disp('Including gunnoise data\n');
    else
        grandz = grandz+sum(tmp_gun.^2,3);
        maxzs(i) = sum(sum(sum(tmp_gun.^2,3)));
        disp('Including conenoise data\n');
    end
    colormap(gray);
end
peakframe = max(maxzs) == maxzs;
% figure; axes; imagesc(grandz);
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
figure(2); subplot(121);imagesc(inRF); hold on;
for k = 1:size(ij,1) % just to look at the cluster membership
    text(ij(k,2),ij(k,1),num2str(T(k)));
end
hold off;

% Comparing estimates of cone weights from gun noise and cone noise.
% At the moment, not normalizing either one by the background lms since
% that will have the same effect on both estimates and so isn't useful
% for the comparison between them.
% Finding whichpixels
% Right now looking at all frames - probably a bad idea
% Am I taking the correct stixels?
whichpix = find(inRF);
tmpSTA = reshape(STAs_cone, [nstixperside^2 3 maxT]);
tmpSTA = permute(tmpSTA, [2 1 3]);
tmpSTA = tmpSTA(:,whichpix,:);
STAconemat = reshape(tmpSTA,[3 length(whichpix)*maxT]);
conesigmas = WN.trial(WN.trial(:,noisetypeidx) == 2, sigmaidxs)/1000;  % cone excitation units

% Cone weights from cone noise expressed in contrast units
%[u,s,v] = svd(STAconemat(:,whichpix)./repmat(lmscontrasts',[1 length(whichpix)]));%
[u,s,v] = svd(STAconemat);
u(:,1) = u(:,1)./mode(abs(conesigmas))';  % Think about this line
coneweights_cone = u(:,1)./sum(abs(u(:,1)));
figure(2); subplot(122); bar(coneweights_cone); set(gca,'XTick',[1 2 3],'XTickLabel',{'L','M','S'});