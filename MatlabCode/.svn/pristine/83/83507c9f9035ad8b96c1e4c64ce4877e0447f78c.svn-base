% Vision research paper figures
%
% Contents
% --------
% Section 1)
% Population analysis of cone weights as estimated by binary cone noise and
% by Gaussian gun noise
%
% Section 2)
% A few frames of an RGB white noise stimulus and an LMS white noise stimulus
%
% Section 3)
% STAs from RGB and LMS white noise (plotting as raw, RGB, images).
%     3.1) Test script to make sure that I'm rendering the LMS STA corectly
%
% Section 4)
% Stimulus distributions projected onto LM plane
%     4.1) RG plane and LM plane
%
% Section 5)
% Looking at the distribution of binary cone noise stimuli in gun
% space. Does this distribution account for the funny-looking STA 
% of L+M cells (cone noise STA)?
%
% Section 6) 
% Not really a figure, but an analysis. Getting a representaion of the
% white noise stimulus in R*/sec that I can send to Fred and he can compare
% the responses of a linear cone model and one that includes realistic
% adaptation.
%
% Section 7)
% Maximum likelihood estimates from whitenoise data.

%%
% Section 1
% Population comparison between cone weights measured with binary cone
% noise and with Gaussian gun noise. Taken wholesale from WNPop section 1.
USECORRECTWEIGHTTRANSFORMATION = 1;
datapath = nexfilepath;
[filenames,spikecds] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\WhiteNoise\ConeVGun2.txt');
maxT = 8;
CHI2CRIT = .95; % For flaging a stixel as significant (adding gun and cone noise z-scores, squared)
data = [];
for i = 1:size(filenames,1)
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
    M = fundamentals'*mon_spd; % In cone excitations or equivalently cone excitation differences
    gammaTable = stro.sum.exptParams.gamma_table;
    gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
    
    % Getting the background rgb/lms
    ridx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_r'));
    gidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_g'));
    bidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_b'));
    bkgndRGB = [mode(stro.trial(:,ridx)), mode(stro.trial(:,gidx)), mode(stro.trial(:,bidx))];
    bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
    bkgndlms = M*bkgndrgb; % this is in cone excitations
    Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
    % If rgb is represented in deltas from bkgnd then M*rgb is in cone excitation differences. 
    % If rgb is represented in deltas from bkgnd then Mrgbtocc*rgb is in cone contrast. 

    gunsigmas = unique(stro.trial(stro.trial(:,noisetypeidx) == 1, sigmaidxs),'rows')/1000;  % normalized gun intensity units
    conesigmas = unique(stro.trial(stro.trial(:,noisetypeidx) == 2, sigmaidxs),'rows')/1000;  % cone excitation units
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Getting STAs
    for noisetype = 1:2
        tmpstro = stro;
        L = stro.trial(:,noisetypeidx) == noisetype;
        tmpstro.ras(~L,:) = [];
        tmpstro.trial(~L,:) = [];
        out = getWhtnsStats(tmpstro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, spikename);
        if (noisetype == 1)
            STAs_gun = out{1};
            STCs_gun = out{2};
            nspikes_gun = out{3};
        elseif (noisetype == 2)
            % The convention I'm using here is that both gun and cone STAs
            % are represented in the space in which they are
            % (approximately) rotationally symmetric. Neither of these is
            % cone excitation difference space; the cone noise stimulus is
            % distended in the S-cone direction. We well need to
            % incorporate a correction for this in the calculation of the
            % cone weights from the cone noise STA.
            % Note: cone sigmas are in cone excitation differences, not cone contrasts.
            % See WhiteNoise.m lines 309-311.
            STAs_cone = out{1};
            STCs_cone = out{2};
            nspikes_cone = out{3};
        end
    end
    clear tmpstro;

    % Z-scoring both gun and cone noise to help us find the RF.
    % Using the first frame as an estimate of the noise.
    s_gun = std(STAs_gun(:,1));
    STAs_gun_z = STAs_gun./s_gun;
    s_cone = std(STAs_cone(:,1));
    STAs_cone_z = STAs_cone./s_cone;
    
    % Spatial map
    grandz = zeros([nstixperside nstixperside]);
    figure;
    for framecounter = 1:maxT
        tmp_gun = reshape(STAs_gun_z(:,framecounter),[nstixperside nstixperside 3]);
        tmp_cone = reshape(STAs_cone_z(:,framecounter),[nstixperside nstixperside 3]);
        grandz = grandz+sum(tmp_gun.^2,3)+sum(tmp_cone.^2,3);
        subplot(3,3,framecounter); imagesc(sum(tmp_gun.^2,3)+sum(tmp_cone.^2,3)); axis square;
    end
    crit = chi2inv(CHI2CRIT,6*maxT); % 6 = 3 color channels * 2 stimulus sets
    L = grandz > crit;
    % Now get largest contiguous block
    [tmpi,tmpj] = ind2sub(size(L),find(L));
    ij = [tmpi,tmpj];
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

    % Comparing estimates of cone weights
    whichframes = [2, 3, 4, 5];
    whichpix = find(inRF);
    tmpSTA = reshape(STAs_gun, [nstixperside^2 3 maxT]);
    tmpSTA = permute(tmpSTA, [2 1 3]);
    tmpSTA = tmpSTA(:,whichpix,whichframes);
    STAgunmat = reshape(tmpSTA,[3 length(whichpix)*length(whichframes)]);
    tmpSTA = reshape(STAs_cone, [nstixperside^2 3 maxT]);
    tmpSTA = permute(tmpSTA, [2 1 3]);
    tmpSTA = tmpSTA(:,whichpix,whichframes);
    STAconemat = reshape(tmpSTA,[3 length(whichpix)*length(whichframes)]);
    
    % First, the gun noise STA
    [u,~,v] = svd(STAgunmat);
    if (sum(v(:,1)) < 0)
        u = -u;
    end
    
%     % Trying to deal with the fact that different contrasts (standard 
%     % were used in the gun and cone noisdeviations) of the
%     % gun and cone noise in the same linear transformation that we use for
%     % the weighting vector. The output of M*rgb is lms (lights) in cone
%     % excitation differences. Dividing by cone sigmas brings them to cone
%     % contrast units.
%     Mprime = diag(1./conesigmas)*M;
    
    if USECORRECTWEIGHTTRANSFORMATION
        coneweights = inv(Mrgbtocc')*u(:,1); % The right way 
    else
        coneweights = Mrgbtocc*u(:,1); % The WRONG (naive) way. LMS = M*rgb
    end
    coneweights_gun = coneweights./sum(abs(coneweights));
    
    [u,~,~] = svd(STAconemat); % STAconemat is in a RS(ish) space which is not cone exc. difference space nor cone contrast space which is where we're trying to go
    if (coneweights_gun'*u(:,1)< 0)
        u = -u;
    end
    % In the line below we're converting from the RS LMS space to cone
    % contrast space. The light transformation from RS
    % to CC is CC = diag(conesigmas)*RS. Conesigmas are in cone contrast
    % units by default.
    % so the mechanism transformation is CED = diag(1./conesigmas)*RS
    Mweirdconespacetocc = diag(conesigmas');
    u(:,1) = inv(Mweirdconespacetocc')*u(:,1);
    coneweights_cone = u(:,1)./sum(abs(u(:,1)));
    data = [data; coneweights_cone' coneweights_gun'];
end

%% Normalized cone weight plot
figure(1); clf; hold on; box on;
plot(data(:,[1,4])',data(:,[2,5])','k-');
for i = 1:size(data,1)
    h = plot(data(i,1),data(i,2),'bs','MarkerFaceColor','b');
    set(h,'ButtonDownFcn',['disp(''',char(filenames{i}(1)),''');'])
    if (data(i,3) < 0)
        set(h,'MarkerFaceColor','none')
    end
    h = plot(data(i,4),data(i,5),'ro','MarkerFaceColor','r');
    set(h,'ButtonDownFcn',['disp(''',char(filenames{i}(1)),''');'])
    if (data(i,6) < 0)
        set(h,'MarkerFaceColor','none')
    end
end
plot([-1 0 1 0 -1],[0 -1 0 1 0],'k-');
axis square;
xlabel('Normalized L-cone weight');
ylabel('Normalized M-cone weight');
ticks = linspace(-1,1,5);
set(gca,'xtick',ticks','ytick',ticks)


%% Normalized cone weight plot
figure(1); clf; axes; hold on; box on;
plot(data(:,[1,4])',data(:,[2,5])','k-');
for i = 1:size(data,1)
    h = plot(data(i,1),data(i,2),'ks','MarkerFaceColor',[.2 .2 .2]);
    set(h,'ButtonDownFcn',['disp(''',char(filenames{i}(1)),''');'])
    if (data(i,3) < 0)
        set(h,'MarkerFaceColor','none')
    end
    h = plot(data(i,4),data(i,5),'ko','MarkerFaceColor',[.2 .2 .2]);
    set(h,'ButtonDownFcn',['disp(''',char(filenames{i}(1)),''');'])
    if (data(i,6) < 0)
        set(h,'MarkerFaceColor','none')
    end
end
plot([-1 0 1 0 -1],[0 -1 0 1 0],'k-');
axis square;
xlabel('Normalized L-cone weight');
ylabel('Normalized M-cone weight');
ticks = linspace(-1,1,5);
set(gca,'xtick',ticks','ytick',ticks)

% For plotting the example cell
cn_cw = conenoise_coneweights ./ sum(conenoise_coneweights);
gn_cw = gunnoise_coneweights ./ sum(gunnoise_coneweights);
h = plot(cn_cw(1),cn_cw(2),'ro','Markerfacecolor','r');
if (cn_cw(3) < 0)
    set(h,'MarkerFaceColor','none')
end
h = plot(gn_cw(1),gn_cw(2),'bs','Markerfacecolor','b');
if (gn_cw(3) < 0)
    set(h,'MarkerFaceColor','none')
end
plot([cn_cw(1) gn_cw(1)],[cn_cw(2) gn_cw(2)],'k-')

% Format and save figure
set(gcf,'PaperPositionMode','auto')
name = 'Fig5_RV';
if ismac
    print(['/Users/jpatrickweller/Dropbox/VisionResearchPaper/' name],'-depsc')
else ispc
    print(['C:\Users\jpweller\Dropbox\VisionResearchPaper\' name],'-depsc');
end
disp(['Fig ' name ' done.'])

%%

figure(2); axes; hold on;
plot(data(:,1),data(:,4),'ro','MarkerFaceColor','red');
plot(data(:,2),data(:,5),'go','MarkerFaceColor','green');
plot(data(:,3),data(:,6),'bo','MarkerFaceColor','blue');
plot([-1 1],[-1 1],'k-');
set(gca,'Xlim',[-1 1],'Ylim',[-1 1]);
axis square;
xlabel('Normalized cone weights (cone noise)');
ylabel('Normalized cone weights (gun noise)');

cw1 = mkbasis(data(:,[1 2 3])')';
cw2 = mkbasis(data(:,[4 5 6])')';
rs = sum(cw1.*cw2,2);
figure(3); clf; hold on; box on;
xlim([0 1]);
%bins = linspace(-.4, 1,20);
bins = linspace(0,1,14);
bincounts = histc(rs,bins);
bar(bins,bincounts,'histc')
set(gca,'xtick',linspace(0,1,5),'YTick',linspace(0,28,5))
ylim([0 28])

if ~USECORRECTWEIGHTTRANSFORMATION
    fprintf('these are the coneweights that were highly correlated even though we used the incorrect transformation');
    data(rs>.8,:)
end
%%
% Section 2
% A single frame from a Gaussian gun noise stimulus and a binary cone noise stimulus
nframes = 3;
stro = nex2stro(findfile('K102109004.nex'));
nstixperside = stro.sum.exptParams.nstixperside;
spds = reshape(stro.sum.exptParams.mon_spd,length(stro.sum.exptParams.mon_spd)/3,3);
spds = SplineSpd([380:4:780]',spds,[380:5:780]');
funds = reshape(stro.sum.exptParams.fundamentals,length(stro.sum.exptParams.fundamentals)/3,3);
M = funds'*spds;

% Getting rid of the synthimage trials
noisetype = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'noise_type'));
stro.trial(~noisetype,:) = [];
stro.ras(~noisetype,:) = [];
noisetype = noisetype(noisetype>0);
sigmas = stro.trial(:,strncmp(stro.sum.trialFields(1,:),'sigma',5));
RGBsigma = mode(sigmas(noisetype == 1,:))./1000;
LMSsigma = mode(sigmas(noisetype == 2,:))./100;
n = 0;
for i = 1:nframes
    
    RGBimage = reshape(normrnd(repmat([.5 .5 .5],nstixperside^2,1),repmat(RGBsigma,nstixperside^2,1)),[nstixperside, nstixperside, 3]);
    
    % Plot
    n = n+1;
    figure(n); clf; set(gcf,'paperpositionmode','auto')
    %subplot(ceil(sqrt(nframes/2)),ceil(sqrt(nframes/2)),i);
    image(RGBimage);
    axis square;
    set(gca,'Xtick',[],'Ytick',[]);
    
    % Format and save figure
    set(gcf,'PaperPositionMode','auto')
    name = ['Fig1_RV_' num2str(n)];
    if ismac
        print(['/Users/jpatrickweller/Dropbox/VisionResearchPaper/' name],'-depsc')
    else ispc
        print(['C:\Users\jpweller\Dropbox\VisionResearchPaper\' name],'-depsc');
    end
    disp(['Fig ' name ' done.'])
    
    % Now the cone noise panel
    bkgndLMS = M*[.5 .5 .5]';
    randombinaryvect = 2*(unidrnd(2,nstixperside^2,3)-1.5);
    ccimage = repmat(LMSsigma,size(randombinaryvect,1),1).*randombinaryvect;
    LMSimage =  repmat(bkgndLMS',size(randombinaryvect,1),1).*(ccimage+1);
    RGBimage = (inv(M)*LMSimage')'; % cones to guns
    RGBimage = reshape(RGBimage,[nstixperside, nstixperside, 3]);
    
    % Plot
    n = n+1;
    figure(n); clf; 
    image(RGBimage);
    axis square;
    set(gcf,'paperpositionmode','auto')
    set(gca,'Xtick',[],'Ytick',[]);
    
    % Format and save figure
    set(gcf,'PaperPositionMode','auto')
    name = ['Fig1_RV_' num2str(n)];
    if ismac
        print(['/Users/jpatrickweller/Dropbox/VisionResearchPaper/' name],'-depsc')
    else ispc
        print(['C:\Users\jpweller\Dropbox\VisionResearchPaper\' name],'-depsc');
    end
    disp(['Fig ' name ' done.'])
    
end


%%
% Section 3
% Now plotting raw (RGB) STAs.
% Getting STAs

USEMECHANISMTRANSFORMATIONS = 1;
USEPROPORTIONALCCSPACE = 0;
%stro = nex2stro(findfile('K032608005.nex')); % lum DS
%stro = nex2stro(findfile('K040108003.nex')); % 2 color-opponent cells
%stro = nex2stro(findfile('K032708003.nex')); % small lum
%stro = nex2stro(findfile('K110708002.nex')); % lum [3,4] is good.
%stro = nex2stro(findfile('K052108002.nex')); % lum frame 4 is good

stro = nex2stro(findfile('K040708003.nex')); % Excellent lum example frames [2 3 4]

%stro = nex2stro(findfile('S022310008.nex')); % blue
%stro = nex2stro(findfile('K043008002.nex')); % red

%stro = strocat(nex2stro(findfile('K051408005.nex')),nex2stro(findfile('K051408006.nex'))); % green whichframe = 4
%stro = nex2stro(findfile('K052308001.nex')); % lime-violet
%stro = nex2stro(findfile('K040108001.nex')); % Red - in paper

%if ~exist('whichframes','var')
    whichframes = 3;
%   whichframes = [2 3 4]
%else
%    fprintf('using pre-existing whichframes (hopefully from section 1): ');
%    fprintf('%d ',whichframes); fprintf('\n');
%end

noisetype = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'noise_type'));
nstixperside = stro.sum.exptParams.nstixperside;
spds = reshape(stro.sum.exptParams.mon_spd,length(stro.sum.exptParams.mon_spd)/3,3);
spds = SplineSpd([380:4:780]',spds,[380:5:780]');
funds = reshape(stro.sum.exptParams.fundamentals,length(stro.sum.exptParams.fundamentals)/3,3);
M = funds'*spds;
gammaTable = stro.sum.exptParams.gamma_table;
gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
sigmas = stro.trial(:,strncmp(stro.sum.trialFields(1,:),'sigma',5));
LMSsigma = mode(sigmas(noisetype == 2,:))./100;

% Getting the background rgb/lms
ridx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_r'));
gidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_g'));
bidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_b'));
bkgndRGB = [mode(stro.trial(:,ridx)), mode(stro.trial(:,gidx)), mode(stro.trial(:,bidx))];
bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
bkgndlms = M*bkgndrgb; % this is in cone excitations
Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
Mcctorgb = inv(Mrgbtocc);
Mcctopropcc = diag(1./LMSsigma); % For converting real to proportional cone contrast. L=0.09 gets mapped to 1.0 and S = 0.4 gets mapped to 1.0
Mpropcctocc = inv(Mcctopropcc);
Mrgbtopropcc = Mcctopropcc*Mrgbtocc;
Mpropcctorgb = inv(Mrgbtopropcc);

spikename = getSpikenum(stro);
spikeidx = find(strcmp(stro.sum.rasterCells(1,:),spikename));
maxT = 8;
nstixperside = stro.sum.exptParams.nstixperside;
for nt = 1:2
    tmpstro = stro;
    L = noisetype == nt;
    tmpstro.ras(~L,:) = [];
    tmpstro.trial(~L,:) = [];
    out = getWhtnsStats(tmpstro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, spikename);
    if (nt == 1)
        STAs_gun = out{1}/out{3};
        STCs_gun = out{2};
        nspikes_gun = out{3};
    elseif (nt == 2)
        STAs_cone = out{1}/out{3};
        STCs_cone = out{2};
        nspikes_cone = out{3};        
    end
end
tmpstro = [];

% -------------------------
% Now getting cone weights
% -------------------------
CHI2CRIT = .95; % For flaging a stixel as significant (adding gun and cone noise z-scores, squared)
s_gun = std(STAs_gun(:,1));
STAs_gun_z = STAs_gun./s_gun;
s_cone = std(STAs_cone(:,1));
STAs_cone_z = STAs_cone./s_cone;

% Spatial map
grandz = zeros([nstixperside nstixperside]);
maxzs = [];
for i = 1:maxT
    tmp_gun = reshape(STAs_gun_z(:,i),[nstixperside nstixperside 3]);
    tmp_cone = reshape(STAs_cone_z(:,i),[nstixperside nstixperside 3]);
    grandz = grandz+sum(tmp_gun.^2,3)+sum(tmp_cone.^2,3);
    maxzs(i) = sqrt(sum(sum(sum(tmp_gun.^2,3)+sum(tmp_cone.^2,3))))-sqrt(600); % "-600" = expected value of chi square
end
crit = chi2inv(CHI2CRIT,6*maxT); % 6 = 3 color channels * 2 stimulus sets
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

whichpix = find(inRF);
tmpSTA = reshape(STAs_gun, [nstixperside^2 3 maxT]);
tmpSTA = permute(tmpSTA, [2 1 3]);
tmpSTA = tmpSTA(:,whichpix,:);
STAgunmat = reshape(tmpSTA,[3 length(whichpix)*maxT]);
tmpSTA = reshape(STAs_cone, [nstixperside^2 3 maxT]);
tmpSTA = permute(tmpSTA, [2 1 3]);
tmpSTA = tmpSTA(:,whichpix,:);
STAconemat = reshape(tmpSTA,[3 length(whichpix)*maxT]);
% Note: we still haven't
% factored in the different sigmas for the different cones. This is
% critical to do before plotting the STA, otherwise it's as if every cone
% type experienced the same contrast, which is not true.

% First, the gun noise STA
[u,~,v] = svd(STAgunmat);
if (sum(v(:,1)) < 0)
    u = -u;
end
gunweights_gun = u(:,1)./sum(abs(u(:,1)));
if USEMECHANISMTRANSFORMATIONS
    if USEPROPORTIONALCCSPACE
        coneweights = inv(Mrgbtopropcc')*u(:,1); % The correct transformation for mechanisms
    else
        coneweights = inv(Mrgbtocc')*u(:,1);
    end
else
    if USEPROPORTIONALCCSPACE
        coneweights = Mrgbtopropcc*u(:,1); % The incorrect transformation (correct for lights)
    else
        coneweights = Mrgbtocc*u(:,1);        
    end
end
coneweights_gun = coneweights./sum(abs(coneweights));
% If we used the RGB->LMS transformation matrix for lights we would get:
M*u(:,1); % not so far off, really.

% Cone weights from cone noise expressed in proportional cone contrast units
[u,s,v] = svd(STAconemat);
if (coneweights_gun'*u(:,1)< 0)
    u = -u;
end

% Stimuli are RS in *proportional* cone contrast space
% This transforms to regular cone contrast space
if ~USEPROPORTIONALCCSPACE
    u = inv(Mpropcctocc)'*u;
    % Now u is in real cone contrast space (the magnitude of each weight has
    % been divided by the standard deviation of the cone contrast (?)
    % modulation).
end
coneweights_cone = u(:,1)./sum(abs(u(:,1)));
if USEMECHANISMTRANSFORMATIONS
    if USEPROPORTIONALCCSPACE
        gunweights_cone = inv(Mpropcctorgb)'*u(:,1);
    else
        gunweights_cone = inv(Mcctorgb)'*u(:,1);
    end
else
    if USEPROPORTIONALCCSPACE
        gunweights_cone = Mpropcctorgb*u(:,1); % using the naive (incorrect) calculation
    else
        gunweights_cone = Mcctorgb*u(:,1); % using the naive (incorrect) calculation        
    end
end
gunweights_cone = gunweights_cone./sum(abs(gunweights_cone));

% Plotting STAs at peak frame(s)
% First, STA from RGB noise
%subplot(2,2,1);
STA = reshape(STAs_gun(:,whichframes),[nstixperside nstixperside 3 length(whichframes)]);
for i = 1:length(whichframes)
    STA(:,:,:,i) = STA(:,:,:,i)*maxzs(whichframes(i));
end
STA = mean(STA,4); % 4 is the dimension corresponding to frame
% First, RGB noise
% Normalizing images
normfactor = .5/max(abs(STA(:)));
STAs_gun_norm = STA.*normfactor+.5;

figure(1); clf; hold on; axis square; box on;
set(gca,'XTick',[],'YTick',[])
set(gcf,'paperpositionmode','auto','Name','Gun Noise STA');
image(STAs_gun_norm);
set(gca,'Ydir','reverse');
axis tight;
keyboard
% Format and save figure
name = 'Fig1_RV_gunsta';
export_fig(['/Users/jpatrickweller/Dropbox/VisionResearchPaper/' name],'-depsc')
disp(['Fig ' name ' done.'])


% Second, STA from LMS noise
% GDLH confirmed that this looks right 5/16/16
% Formula for making RGBs for cone noise is:
% inv(M)*((M*bkgndrgb) + (sigmalms.*binaryvector))
% Also see Section 3.1, below.

STA = reshape(STAs_cone(:,whichframes),[nstixperside nstixperside 3 length(whichframes)]);
for i = 1:length(whichframes)
    STA(:,:,:,i) = STA(:,:,:,i)*maxzs(whichframes(i)); % temporal weighting if using > 1 frame
end
STA = mean(STA,4); % 4 is the dimension corresponding to frame
STAconemat = reshape(STA,nstixperside^2, 3)';

% *********** 
% Below, need to account for the fact that the contrasts differed between cone types
% Otherwise the STA is plotted as if all three cone types underwent the
% same modulation in cone contrast, which is not true.
% So the LMS STA is in cone contrast, not proportional cone contrast units

% Uncomment out the line below to use cone contrast. Keep it in to use
% proportional cone contrast.
% STAconemat = STAconemat./repmat(LMSsigma',1,nstixperside^2);
% *********** 
STA_cone_RGB = (inv(M)*STAconemat)'; % "STA_cone_RGB" is in RGBs
% This is supposed to be the STA (before the background is added back in)
% of the LMS noise rendered in RGB
normfactor = .5/max(max(abs(STA_cone_RGB)));
STAs_cone_RGB = STA_cone_RGB.*normfactor+.5;

figure(2); clf; hold on; axis square; box on;
set(gca,'XTick',[],'YTick',[])
set(gcf,'paperpositionmode','auto','Name','Cone Noise STA');
image(reshape(STAs_cone_RGB,[nstixperside nstixperside 3]));
set(gca,'Ydir','reverse');
axis tight
%print(['C:\Users\ghorwitz\Documents\MATLAB\' filename 'coneSTA'],'-depsc')

% Format and save figure
name = 'Fig1B_conesta_RV';
export_fig(['/Users/jpatrickweller/Dropbox/VisionResearchPaper/' name],'-depsc')
disp(['Fig ' name ' done.'])


%% Individual cone weights maps
for i = 1:3
    figure(i); clf; hold on; box on;
    im = reshape(STAconemat(i,:)./max(abs(STAconemat(:)))+.5,[nstixperside nstixperside]);
    image(im*255);
    colormap(gray(255));
    set(gca,'Ydir','reverse');
    axis tight square
    axis square;
    set(gca,'XTick',[],'YTick',[]);
end

% Individual gun weight maps
STA = reshape(STAs_gun(:,whichframes),[nstixperside nstixperside 3 length(whichframes)]);
for i = 1:length(whichframes)
    STA(:,:,:,i) = STA(:,:,:,i)*maxzs(whichframes(i)); % temporal weighting if using > 1 frame
end
STA = mean(STA,4); % 4 is the dimension corresponding to frame
STAgunmat = reshape(STA,nstixperside^2, 3)';
figure;
for i = 1:3
    subplot(1,3,i);
    im = reshape(STAgunmat(i,:)./max(abs(STAgunmat(:)))+.5,[nstixperside nstixperside]);
    image(im*255);
    colormap(gray(255));
    set(gca,'Ydir','reverse');
    axis tight
    axis square; box on;
    set(gca,'XTick',[],'YTick',[])
    if (i == 2)
        title('RGB');
    end
end


%subplot(2,2,3)
figure;
bar([coneweights_cone coneweights_gun]);
h = legend({'cone','gun'},'location','NorthEast');
set(gca,'XTick',[1 2 3],'XTickLabel',{'L','M','S'});
title(stro.sum.fileName(find(stro.sum.fileName == filesep,1,'last')+1:end))
print(['C:\Users\ghorwitz\Documents\MATLAB\' filename 'correlation'],'-depsc')

figure;
bar([gunweights_cone gunweights_gun]);
h = legend({'cone','gun'},'location','NorthEast');
set(gca,'XTick',[1 2 3],'XTickLabel',{'R','G','B'});
title(stro.sum.fileName(find(stro.sum.fileName == filesep,1,'last')+1:end))

% % % ---------------------------------------------
% % % Trying to figure out why an off luminance subunit
% % % might have a postive blue component when the STA
% % % of LMS noise is computed in RGB space.
% % % ---------------------------------------------
% deltacones = (fullfact([2,2,2])-1.5)*2.*repmat(LMSsigma,8,1);
% deltargbs = (inv(Mrgbtocc)*deltacones')';
% [deltacones, deltargbs];
% 
% LnegM = sign(deltacones(:,2)) == -1;
% sum(deltargbs(LnegM,:))
% % Interesting. The stimuli in the white cone binary noise set
% % that are -M have a *postive* red component!
% 
% LnegL = sign(deltacones(:,1)) == -1;
% sum(deltargbs(LnegL,:))
% % Both L- and M- stimuli have positive blue components. Something
% % must be wrong with my RGB rendering.
% 
%%
% Section 3.1
% Double checking how I'm computing the RGB STA of the LMS cone noise.
% Need to run the above cell first.

coneweights = [1 1 0];

% Getting the RGBs of each of the 8 stimuli in the cone noise
% deltacones is in cone contrasts. Converting to cone excitation
% differences.
cone_exc_diff = deltacones.*repmat(bkgndlms',size(deltacones,1),1);
% Now coneverting cone excitation differences to RGBs
RGBs = (inv(M)*cone_exc_diff')'
% Pretending that each stimulus is shown once.
dotproducts = deltacones*coneweights';
dotproducts(dotproducts<0) = 0;
STA_RGB = dotproducts'*RGBs;
STA_cone_exc_diff = dotproducts'*cone_exc_diff;

% Transforming the STA (in cone excitation differences) to RGBs so that I
% can render it.
inv(M)*STA_cone_exc_diff'
%%
% Section 4
% Visualizing stimulus distributions projected onto the LM plane
stro = nex2stro(findfile('K102109004.nex')); % lum DS

noisetype = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'noise_type'));
nstixperside = stro.sum.exptParams.nstixperside;
spds = reshape(stro.sum.exptParams.mon_spd,length(stro.sum.exptParams.mon_spd)/3,3);
spds = SplineSpd([380:4:780]',spds,[380:5:780]');
funds = reshape(stro.sum.exptParams.fundamentals,length(stro.sum.exptParams.fundamentals)/3,3);
M = funds'*spds;
gammaTable = stro.sum.exptParams.gamma_table;
gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
gammaTable1 = interp1(linspace(0,255,256),gammaTable,linspace(0,255,65536), 'spline');
invgamma = InvertGamma(gammaTable, 0);

sigmas = stro.trial(:,strncmp(stro.sum.trialFields(1,:),'sigma',5));
LMSsigma = mode(sigmas(noisetype == 2,:))./100;
muidxs = strncmp('mu',stro.sum.trialFields(1,:),2);
sigmaidxs = strncmp('sigma',stro.sum.trialFields(1,:),5);
bkgndidxs = strncmp('bkgnd',stro.sum.trialFields(1,:),5);
bkgndRGB = unique(stro.trial(:,bkgndidxs),'rows');
bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
bkgndlms = M*bkgndrgb;

sigmargb = unique(stro.trial(noisetype == 1, sigmaidxs),'rows')/1000;
sigmascone = unique(stro.trial(noisetype == 2, sigmaidxs),'rows')/100;

% Cone isolating directions in gun space
lconergb = inv(M)*[1 0 0]';
mconergb = inv(M)*[0 1 0]';
sconergb = inv(M)*[0 0 1]';

% Making a gaussian cloud in 3D gun space
nstim = 1000;
rgb = normrnd(zeros(nstim,3),repmat(sigmargb,nstim,1),nstim,3); % centered on zero
Lomit = logical(sum(rgb> repmat(1-bkgndrgb',nstim,1) | rgb < repmat(-bkgndrgb',nstim,1),2));
rgb(Lomit,:) = []; nstim = sum(~Lomit);
%plot(rgb(:,1),rgb(:,2),'ko','MarkerFaceColor','black')
%set(gca,'Xlim',[0 1],'Ylim',[0 1]);
%axis square;

% Looking at stimuli in cone contrast space
lms = (rgb*M')./repmat(bkgndlms',nstim,1); % cone contrasts of Gaussian gun noise stimuli
figure(1); clf; hold on; box on;
plot(lms(:,1),lms(:,2),'ko','MarkerFaceColor','black','MarkerSize',4);
plot(sigmascone(1), sigmascone(2),'cs','MarkerSize',10,'MarkerFaceColor','cyan','MarkerEdgeColor','black')
plot(sigmascone(1), -sigmascone(2),'rs','MarkerSize',10,'MarkerFaceColor','red','MarkerEdgeColor','black')
plot(-sigmascone(1), sigmascone(2),'gs','MarkerSize',10,'MarkerFaceColor','green','MarkerEdgeColor','black')
plot(-sigmascone(1), -sigmascone(2),'ks','MarkerSize',10,'MarkerFaceColor','y','MarkerEdgeColor','black')
axis square; set(gca,'Xlim',[-1 1],'Ylim',[-1 1]);
xlabel('L-cone contrast');
ylabel('M-cone contrast');

% Making a projection matrix to project rgbs orthogonal to the S-cone axis
P = eye(3)-sconergb*inv(sconergb'*sconergb)*sconergb';
% sanity check: P*sconergb
funnyrgbtransform =  mkbasis([P*lconergb, P*mconergb],'RotOrtho');
% sanity check: [lconergb'; mconergb'; sconergb']*funnyrgbtransform

rgbinfunnyspace = rgb*funnyrgbtransform;
figure; axes; hold on;
plot(rgbinfunnyspace(:,1),rgbinfunnyspace(:,2),'ko','MarkerFaceColor','black','MarkerSize',4);
axis square; set(gca,'Xlim',[-.5 .5],'Ylim',[-.5 .5]);

% Getting the RGBs of the cone white noise stimuli
ccs = (fullfact([2,2,2])-1.5)*2.*repmat(sigmascone,8,1)
conediffs = ccs.*repmat(bkgndlms',size(ccs,1),1);
rgbsofconestim = inv(M)*conediffs';
lmsinfunnyspace = rgbsofconestim'*funnyrgbtransform; 
% bottom four rows of lmsinfunnyspace are the same as top four rows since
% they only differ in the S-cone component (which has been projected out).
% This is a good sanity check.
plot(lmsinfunnyspace(1,1),lmsinfunnyspace(1,2),'ks','MarkerSize',10,'MarkerFaceColor','y','MarkerEdgeColor','black')
plot(lmsinfunnyspace(2,1),lmsinfunnyspace(2,2),'rs','MarkerSize',10,'MarkerFaceColor','red','MarkerEdgeColor','black')
plot(lmsinfunnyspace(3,1),lmsinfunnyspace(3,2),'gs','MarkerSize',10,'MarkerFaceColor','green','MarkerEdgeColor','black')
plot(lmsinfunnyspace(4,1),lmsinfunnyspace(4,2),'cs','MarkerSize',10,'MarkerFaceColor','cyan','MarkerEdgeColor','black')

axis square;
xlabel('funny axis 1');
ylabel('funny axis 2');


%%
% Section 4.1
% Visualizing stimulus distributions projected onto the Red-Green plane and the LM plane
stro = nex2stro(findfile('K102109004.nex')); % lum DS

noisetype = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'noise_type'));
nstixperside = stro.sum.exptParams.nstixperside;
spds = reshape(stro.sum.exptParams.mon_spd,length(stro.sum.exptParams.mon_spd)/3,3);
spds = SplineSpd([380:4:780]',spds,[380:5:780]');
funds = reshape(stro.sum.exptParams.fundamentals,length(stro.sum.exptParams.fundamentals)/3,3);
M = funds'*spds;
gammaTable = stro.sum.exptParams.gamma_table;
gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
gammaTable1 = interp1(linspace(0,255,256),gammaTable,linspace(0,255,65536), 'spline');
invgamma = InvertGamma(gammaTable, 0);

sigmas = stro.trial(:,strncmp(stro.sum.trialFields(1,:),'sigma',5));
sigmaidxs = strncmp('sigma',stro.sum.trialFields(1,:),5);
sigmargb = unique(stro.trial(noisetype == 1, sigmaidxs),'rows')/1000;
sigmascone = unique(stro.trial(noisetype == 2, sigmaidxs),'rows')/100;
muidxs = strncmp('mu',stro.sum.trialFields(1,:),2);
bkgndidxs = strncmp('bkgnd',stro.sum.trialFields(1,:),5);
bkgndRGB = unique(stro.trial(:,bkgndidxs),'rows');
bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
bkgndlms = M*bkgndrgb;
Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
Mcctorgb = inv(Mrgbtocc);
Mcctopropcc = diag(1./sigmascone); % For converting real to proportional cone contrast. L=0.09 gets mapped to 1.0 and S = 0.4 gets mapped to 1.0
Mpropcctocc = inv(Mcctopropcc);
Mrgbtopropcc = Mcctopropcc*Mrgbtocc;
Mpropcctorgb = inv(Mrgbtopropcc);

% Cone isolating directions in gun space
lconergb = inv(Mrgbtocc)*[sigmascone(1) 0 0]';
mconergb = inv(Mrgbtocc)*[0 sigmascone(2) 0]';
sconergb = inv(Mrgbtocc)*[0 0 sigmascone(3)]';

% Making a gaussian cloud in 3D gun space
nstim = 1000;
rgb = normrnd(zeros(nstim,3),repmat(sigmargb,nstim,1),nstim,3); % centered on zero
Lomit = logical(sum(rgb> repmat(1-bkgndrgb',nstim,1) | rgb < repmat(-bkgndrgb',nstim,1),2));
rgb(Lomit,:) = []; nstim = sum(~Lomit);

% R/G plane
figure(1); clf; hold on; box on; axis square;
plot(rgb(:,1),rgb(:,2),'o','MarkerFaceColor','k','MarkerEdgeColor','w')
conestimmat = 2*(fullfact([2 2 2])-1.5); % proportional cone contrast
conestimcols = (conestimmat+1)./2;
conenoisergb = zeros(size(conestimmat));
for i = 1:size(conestimmat,1)
    conenoisergb(i,:) =[lconergb mconergb sconergb]*conestimmat(i,:)';
    h = plot(conenoisergb(i,1),conenoisergb(i,2),'s');
    %set(h(i),'MarkerSize',10,'MarkerFaceColor',conestimcols(i,:),'MarkerEdgeColor','black'); % change this line to change the color of the cone noise symbols
    set(h,'MarkerSize',10,'MarkerFaceColor',conestimcols(i,:)); % change this line to change the color of the cone noise symbols
end
ticks = linspace(-.5,.5,5);
set(gca,'XLim',[-.5 .5],'YLim',[-.5 .5],'XTick',ticks,'YTick',ticks);
xlabel('Red phosphor intensity');
ylabel('Green phosphor intensity');

% Format and save figure
set(gcf,'PaperPositionMode','auto')
name = 'Fig3a_RV';
if ismac
    print(['/Users/jpatrickweller/Dropbox/VisionResearchPaper/' name],'-depsc')
else ispc
    print(['C:\Users\jpweller\Dropbox\VisionResearchPaper\' name],'-depsc');
end
disp(['Fig ' name ' done.'])

% LM plane
figure(2); clf; hold on; box on; axis square;
lms = (rgb*M')./repmat(bkgndlms',nstim,1); % cone contrasts of Gaussian gun noise stimuli
conenoiselms = conenoisergb * Mrgbtocc';
plot(lms(:,1),lms(:,2),'ko');
for i = 1:size(conestimmat,1)
    h = plot(conenoiselms(i,1),conenoiselms(i,2),'s');
    set(h,'MarkerSize',10,'MarkerFaceColor',conestimcols(i,:)); % change this line to change the color of the cone noise symbols
end
ticks = linspace(-.8,.8,5);
set(gca,'Xlim',[-.8 .8],'XTick',ticks,'Ylim',[-.8 .8],'YTick',ticks);
xlabel('L-cone contrast');
ylabel('M-cone contrast');

% Format and save figure
set(gcf,'PaperPositionMode','auto')
name = 'Fig3b_RV';
if ismac
    print(['/Users/jpatrickweller/Dropbox/VisionResearchPaper/' name],'-depsc')
else ispc
    print(['C:\Users\jpweller\Dropbox\VisionResearchPaper\' name],'-depsc');
end
disp(['Fig ' name ' done.'])

%%
% Section 5) 
% Examining the locations of the cone noise stimuli in gun space. Can I
% predict the STA of an L+M cell from these.

coneweights = [-0.1341, -0.7464, -0.1194]';% coneweights of cell K032708003 estimated from cone noise
%coneweights = [-0.2859 -0.5587 -0.1553]';  % coneweights of cell K032708003 estimated from gun noise
coneweights = [0.4444, -0.4862, -0.0694]'; %coneweights of cell K043008002 estimated from cone noise

stro = nex2stro(findfile('K032708003.nex')); % lum (just needed for M matrix, sigmas, etc.)

noisetype = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'noise_type'));
nstixperside = stro.sum.exptParams.nstixperside;
spds = reshape(stro.sum.exptParams.mon_spd,length(stro.sum.exptParams.mon_spd)/3,3);
spds = SplineSpd([380:4:780]',spds,[380:5:780]');
funds = reshape(stro.sum.exptParams.fundamentals,length(stro.sum.exptParams.fundamentals)/3,3);
M = funds'*spds;
gammaTable = stro.sum.exptParams.gamma_table;
gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
gammaTable1 = interp1(linspace(0,255,256),gammaTable,linspace(0,255,65536), 'spline');
invgamma = InvertGamma(gammaTable, 0);

sigmas = stro.trial(:,strncmp(stro.sum.trialFields(1,:),'sigma',5));
LMSsigma = mode(sigmas(noisetype == 2,:))./100;
muidxs = strncmp('mu',stro.sum.trialFields(1,:),2);
sigmaidxs = strncmp('sigma',stro.sum.trialFields(1,:),5);
bkgndidxs = strncmp('bkgnd',stro.sum.trialFields(1,:),5);
bkgndRGB = unique(stro.trial(:,bkgndidxs),'rows');
bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
bkgndlms = M*bkgndrgb;
sigmargb = unique(stro.trial(noisetype == 1, sigmaidxs),'rows')/1000;
sigmascone = unique(stro.trial(noisetype == 2, sigmaidxs),'rows')/100;

% Cone isolating directions in gun space
lconergb = inv(M)*[1 0 0]';
mconergb = inv(M)*[0 1 0]';
sconergb = inv(M)*[0 0 1]';
ccs = (fullfact([2,2,2])-1.5)*2.*repmat(sigmascone,8,1);
conediffs = ccs.*repmat(bkgndlms',size(ccs,1),1);
rgbsofconestim = inv(M)*conediffs'; % relative to rgb_background

% Now finding the gun weight vector corresponding to the cone weights
rgbweights = M'*coneweights;
rgbweights = rgbweights./norm(rgbweights);

figure; axes; hold on
plot3([0 rgbweights(1)],[0 rgbweights(2)],[0 rgbweights(3)],'k-','LineWidth',3);
RWScone = [0 0 0]'; % Initialling response-weighted average of the cone stimuli
for i = 1:size(rgbsofconestim,2)
    h = plot3(rgbsofconestim(1,i),rgbsofconestim(2,i),rgbsofconestim(3,i),'o','MarkerEdgeColor','none','MarkerFacecolor',.5+rgbsofconestim(:,i)');
    projection = rgbweights'*rgbsofconestim(:,i);
    set(h,'MarkerSize',6+max(0,projection*200));
    RWScone = RWScone+rgbsofconestim(:,i)*max(0,projection);
end
extremeval = max(abs(rgbsofconestim(:)));
set(gca,'xlim',extremeval*[-1.1 1.1],'ylim',extremeval*[-1.1 1.1],'zlim',extremeval*[-1.1 1.1]);

% Finding the average color of the preferred stimulus
RWScone = RWScone./norm(RWScone);
h = plot3([0 RWScone(1)],[0 RWScone(2)],[0 RWScone(3)],'-','LineWidth',5);
set(h,'Color',RWScone/5+.5);
set(gca,'Color',[.5 .5 .5]);
axis square;


%%
% Section 6
% Getting a representation of the Gaussian white noise stimulus in terms of
% R*/sec. Fred said that he'd take these and feed them through his cone
% model so that we can tell the reviewer how much cone specific adaptation
% matters under the conditions of our experiment.

load rgb2Rstar % This is in Charlie's directory and it appropriate for one of 
% the old Dell monitors. Almost certainly Dell4. 
stro = nex2stro(findfile('K032708003.nex')); % lum (just needed for M matrix, sigmas, etc.)
bkgndRGB = stro.trial(1,strncmp(stro.sum.trialFields(1,:),'bkgnd',5));
allsigmas = stro.trial(:,strncmp(stro.sum.trialFields(1,:),'sigma',5));
noisetype = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'noise_type'));
framerate = stro.sum.exptParams.framerate;
msperframe = 1000/framerate;
gammatable = stro.sum.exptParams.gamma_table;
gammatable = reshape(gammatable, size(gammatable,1)/3,3);
bkgndrgb = [gammatable(bkgndRGB(1)+1,1), gammatable(bkgndRGB(2)+1,2), gammatable(bkgndRGB(3)+1,3)];
bkgndRstar = rgb2Rstar*bkgndrgb';

% Now generating some RGB white noise
nframes = 100;
pbounds = [stro.sum.exptParams.gauss_locut stro.sum.exptParams.gauss_hicut]/1000;
sigmas = mode(allsigmas(noisetype == 1,:))/1000;
sigmabounds = norminv(pbounds,[0 0],[sigmas(1) sigmas(2)]);
rgbs = normrnd(repmat([0 0 0],nframes,1),repmat(sigma,nframes,1), nframes,3);
rgbs(rgbs<sigmabounds(1)) = sigmabounds(1);
rgbs(rgbs>sigmabounds(2)) = sigmabounds(2);
Rstars = rgb2Rstar*(rgbs+repmat(bkgndrgb,nframes,1))';
Rstars = permute(repmat(Rstars',[1 1 round(msperframe)]),[3 1 2]);
Rstars_GWN = reshape(Rstars, size(Rstars,1)*size(Rstars,2),size(Rstars,3));

% Now generating some LMS white noise
sigmas = mode(allsigmas(noisetype == 2,:))/100;
lms = (unidrnd(2,100,3)-1.5)*2;
lmscc = lms.*repmat(sigmas,nframes,1);
% Going from cone contrast to R*/sec
% cc = (A-Ab)/Ab
% A = Ab*(1+cc)
Rstars = repmat(bkgndRstar',nframes,1).*(1+lmscc);
Rstars = permute(repmat(Rstars,[1 1 round(msperframe)]),[3 1 2]);
Rstars_BWN = reshape(Rstars, size(Rstars,1)*size(Rstars,2),size(Rstars,3));

%%
% Section 7 Trying to do maximum likelihood fitting of whitenoise data

%stro = strocat(nex2stro(findfile('K051408005.nex')),nex2stro(findfile('K051408006.nex'))); % green
%whichframe = 4; % The way things are set up right now this can only have one element.

stro = nex2stro(findfile('K040708003.nex')); %  lum example
whichframe = 3;

noisetype = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'noise_type'));
nstixperside = stro.sum.exptParams.nstixperside;
spds = reshape(stro.sum.exptParams.mon_spd,length(stro.sum.exptParams.mon_spd)/3,3);
spds = SplineSpd([380:4:780]',spds,[380:5:780]');
funds = reshape(stro.sum.exptParams.fundamentals,length(stro.sum.exptParams.fundamentals)/3,3);
M = funds'*spds;
gammaTable = stro.sum.exptParams.gamma_table;
gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
sigmas = stro.trial(:,strncmp(stro.sum.trialFields(1,:),'sigma',5));
LMSsigma = mode(sigmas(noisetype == 2,:))./100;

% Getting the background rgb/lms
ridx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_r'));
gidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_g'));
bidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_b'));
bkgndRGB = [mode(stro.trial(:,ridx)), mode(stro.trial(:,gidx)), mode(stro.trial(:,bidx))];
bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
bkgndlms = M*bkgndrgb; % this is in cone excitations
Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences. right multiply by [rgb - bkgndrgb]' =
Mcctopropcc = diag(1./LMSsigma); % For converting real cone contrast to proportional cone contrast
Mpropcctocc = inv(Mcctopropcc);
Mrgbtopropcc = Mcctopropcc*Mrgbtocc; % Make sure this is right

spikename = getSpikenum(stro);
spikeidx = find(strcmp(stro.sum.rasterCells(1,:),spikename));
maxT = 8;
nstixperside = stro.sum.exptParams.nstixperside;
options = optimoptions('fmincon','Display','iter','MaxFunEvals',10^6,'TolX',10^-4,'TolFun',10^-4);

strocellarray = {};
for nt = 1:2
    tmpstro = stro;
    L = noisetype == nt;
    tmpstro.ras(~L,:) = [];
    tmpstro.trial(~L,:) = [];
    out = getWhtnsStats(tmpstro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, spikename);
    if (nt == 1)
        STAs_gun = out{1}/out{3}; % each column is nstixperside x nstixperside x 3
    elseif (nt == 2)
        STAs_cone = out{1}/out{3}; 
    end
    strocellarray{nt} = tmpstro;
end

nframes = [sum(strocellarray{1}.trial(:,12)) sum(strocellarray{2}.trial(:,12))];
estimatednumberofstimuli = max(nframes);
% First, doing ML fitting on gun noise, first to calculate weights in gun
% space, then to calculate weights in proportional cone contrast space.
% At the moment, this code assumes we're estimating only 1 frame-worth of parameters.

basisvectlength = stro.sum.exptParams.nstixperside.^2*3*length(whichframe);
% STPROJmod wants color to change most quickly, then space, then time.
basisvectors = STAs_gun(:,whichframe)./norm(STAs_gun(:,whichframe));

mynorm = @(x) deal([],norm(x(5:end))-1); % out arguments: c, ceq
% Below, development code
% computing Poisson likelihood
%lambda = NRparams(1).*(max(0,out{1}).^NRparams(2)./(NRparams(3).^NRparams(2)+out{1}.^NRparams(2)))+NRparams(4);
% fixing the weights and fitting the other parameters
%ll = -1*sum(out{2}.*log(lambda)-lambda)

% negloglik is an anonymous function that takes a matrix of stimulus projections
% (x(:,1)) and spike counts (x(:,2)), and a vector of Naka-Rushton
% parameters (p = [amplitude, exponent, C50, baseline]).
negloglik = @(p, x) -1*sum(x(:,2).*log(p(1).*(max(x(:,1),0).^p(2)./(p(3).^p(2)+max(x(:,1),0).^p(2)))+p(4)) - ...
    (p(1).*(max(x(:,1),0).^p(2)./(p(3).^p(2)+max(x(:,1),0).^p(2)))+p(4)));
%negloglik(NRparams,[out{1} out{2}])
% mynegloglik is an anonymous function that calls another anonymous
% function (negloglik). The rationale for this make a function that fixes
% the spike counts and stimulus projections and lets us try a bunch of
% different parameters.
out = getWhtnsStats(strocellarray{1},maxT,'STPROJmod', {basisvectors, min(whichframe)-1, estimatednumberofstimuli, [stro.sum.exptParams.nstixperside.^2 3 length(whichframe)] }, spikename);
mynegloglik = @(params) negloglik(params, [out{1} out{2}]); % Keeping the data, [out{1} out{2}], fixed
%mynegloglik(NRparams)

% Doing a quick a dirty analysis to get parameters for the Naka-Rushton
% function. Fixing the projections and firing rate and just varying the
% parameters of the Naka-Rushton function: amplitude, exponent, C50, baseline. 
NRparams = [2, .2, 1, .01]; % Initial guess for gun weights from gun noise. These are appropriate for K051408005/6
NRparams = [5, .7, 5, .1]; % amplitude, exponent, C50, baseline. Initial guess for gun weights from gun noise. These are appropriate for lum cell

[NRparamsout,fv] = fmincon(mynegloglik,NRparams,[],[],[],[],[0 .1 0 0],[20 3 5 1],[],options);

%Sanity check
%lambda = NRparamsout(1).*(max(0,out{1}).^NRparamsout(2)./(NRparamsout(3).^NRparamsout(2)+max(0,out{1}).^NRparamsout(2)))+NRparamsout(4);
%mynegloglik(NRparamsout)

% Fitting the weights and the parameters of the Naka-Rushton simultaneously.
LB = [0 .1 0 0 min(out{1})*ones(1,basisvectlength)];
UB = [20 3 5 1 max(out{1})*ones(1,basisvectlength)];
clear VisResMLerr; % Clear out persistent variables (this function is in VisResMLerr.m)
VisResMLfiterr('init',strocellarray{1},whichframe,spikename,maxT); % Setting up persistent parameters for the fitting
[paramsout_gun_gun,~] = fmincon('VisResMLfiterr',[NRparamsout'; basisvectors],[],[],[],[],LB,UB,mynorm,options); % doing the fitting

% Plotting the ML fit of the gun weights from the gun noise, and the STA
% (for comparison). The convention of the variable names
% ("paramsout_gun_gun") is that the first "gun" refers to the stimulus
% distribution and the second "gun" refers to the space in which the
% weights are calculated (what the weights represent).
figure;
weights = reshape(paramsout_gun_gun(5:end),nstixperside, nstixperside, 3,length(whichframe));
subplot(2,1,1)
image(weights./max(abs(2*weights(:)))+.5);
axis square;
title('ML gun noise, gun weights'); set(gca,'Xtick',[],'Ytick',[])
subplot(2,1,2);
weights = reshape(basisvectors, nstixperside, nstixperside, 3);
image(weights./max(abs(2*weights(:)))+.5);
axis square; set(gca,'Xtick',[],'Ytick',[])
title('STA (initial guess)');

% Now doing ML on proportional cone contrast (gun noise)
LB = [0 .1 0 0 -1*ones(1,basisvectlength)];
UB = [20 3 5 1 ones(1,basisvectlength)];
initialguess = inv(Mrgbtopropcc')*reshape(basisvectors,length(basisvectors)/3,3)'; % This doesn't work very well. Not sure why
initialguess = initialguess';
initialguess = initialguess(:)./norm(initialguess(:));

% Trying this crazy initial guess, which seems to work better for some reason
initialguess = basisvectors./norm(basisvectors);

out = getWhtnsStats(strocellarray{1},maxT,'STPROJmod', {initialguess, min(whichframe)-1, estimatednumberofstimuli, [stro.sum.exptParams.nstixperside.^2 3 length(whichframe)] }, spikename);
mynegloglik = @(params) negloglik(params, [out{1} out{2}]); % Keeping the data, [out{1} out{2}], fixed
[NRparamsout,fv] = fmincon(mynegloglik,NRparams,[],[],[],[],[0 .1 0 0],[20 3 5 1],[],options);
% For luminance cell, the best NR params are [2.5, .7, 2.8, .1]
clear VisResMLerr; % Clear out persistent variables (this function is in VisResMLerr.m)
VisResMLfiterr('init',strocellarray{1},whichframe,spikename,maxT,Mrgbtopropcc);
[paramsout_gun_cone,fv] = fmincon('VisResMLfiterr',[NRparamsout'; initialguess(:)],[],[],[],[],LB,UB,mynorm,options);
weights = reshape(paramsout_gun_cone(5:end),nstixperside, nstixperside, 3,length(whichframe));

figure;
normfact_weights = max(abs(weights(:)));
normfact_initialguess = max(abs(initialguess(:)));
for i = 1:3
    subplot(2,3,i)
    image(255*(weights(:,:,i)./(2*normfact_weights)+.5))
    axis square;
    set(gca,'Xtick',[],'Ytick',[]);
    if i == 2
        title('ML Gun noise, cone weights')
    end
    
    subplot(2,3,i+3)
    image(255*(reshape(initialguess([1:nstixperside^2]+nstixperside^2*(i-1),:),nstixperside,nstixperside)./(2*normfact_initialguess)+.5))
    axis square;
    set(gca,'Xtick',[],'Ytick',[]);
    if i == 2
        title('Initialguess from STA')
    end
end
colormap(gray(255));

% Now doing ML on proportional cone contrast (cone noise)
basisvectors = STAs_cone(:,whichframe)./norm(STAs_cone(:,whichframe));
initialguess = basisvectors;
LB = [0 .1 0 0 -5*ones(1,basisvectlength)];
UB = [20 3 5 1 5*ones(1,basisvectlength)];
NRparamsout = [3 1 5 0.1]; % temporary intial guess. This could be more principled
clear VisResMLerr; % Clear out persistent variables (this function is in VisResMLerr.m)
VisResMLfiterr('init',strocellarray{2},whichframe,spikename,maxT);
[paramsout_cone_cone,fv] = fmincon('VisResMLfiterr',[NRparamsout'; initialguess(:)],[],[],[],[],LB,UB,mynorm,options);
weights = reshape(paramsout_cone_cone(5:end),nstixperside, nstixperside, 3,length(whichframe));

figure;
normfact_weights = max(abs(weights(:)));
normfact_initialguess = max(abs(initialguess(:)));
for i = 1:3
    subplot(2,3,i)
    image(255*(weights(:,:,i)./(2*normfact_weights)+.5))
    axis square;
    set(gca,'Xtick',[],'Ytick',[]);
    if i == 2
        title('ML Cone noise, cone weights')
    end
    
    subplot(2,3,i+3)
    image(255*(reshape(initialguess([1:nstixperside^2]+nstixperside^2*(i-1),:),nstixperside,nstixperside)./(2*normfact_initialguess)+.5))
    axis square;
    set(gca,'Xtick',[],'Ytick',[]);
    if i == 2
        title('Initialguess from STA')
    end
end
colormap(gray(255));

% Now doing ML on gunweights (using cone noise)
basisvectors = STAs_cone(:,whichframe);
initialguess = Mrgbtopropcc'*reshape(basisvectors,length(basisvectors)/3,3)'; % inverse-transpose of light transformation
initialguess = initialguess';
initialguess = initialguess(:)./norm(initialguess);

LB = [0 .1 0 0 -5*ones(1,basisvectlength)];
UB = [20 3 5 1 5*ones(1,basisvectlength)];
clear VisResMLerr; % Clear out persistent variables (this function is in VisResMLerr.m)
VisResMLfiterr('init',strocellarray{2},whichframe,spikename,maxT); % Setting up persistent parameters for the fitting
[paramsout_cone_gun,~] = fmincon('VisResMLfiterr',[paramsout_gun_gun(1:4); initialguess],[],[],[],[],LB,UB,mynorm,options); % doing the fitting
weights = reshape(paramsout_cone_gun(5:end),nstixperside, nstixperside, 3,length(whichframe));

% color figure
figure;
weights = reshape(paramsout_cone_gun(5:end),nstixperside, nstixperside, 3,length(whichframe));
subplot(2,1,1)
image(weights./max(abs(2*weights(:)))+.5);
axis square;
title('ML cone noise, gun weights'); set(gca,'Xtick',[],'Ytick',[])
subplot(2,1,2);
weights = reshape(initialguess, nstixperside, nstixperside, 3);
image(weights./max(abs(2*weights(:)))+.5);
axis square; set(gca,'Xtick',[],'Ytick',[])
title('STA (initial guess)');

% grayscale figure
weights = reshape(paramsout_cone_gun(5:end),nstixperside, nstixperside, 3,length(whichframe));
figure;
for i = 1:3
    subplot(1,3,i)
    image(255*(weights(:,:,i)./max(abs(2*weights(:)))+.5));
    axis square;
    colormap(gray(255))
    set(gca,'Xtick',[],'YTick',[]);
    if i == 2
        title('ML cone noise, gun weights');
    end
end
% Extraction of weights across space. Could to better by eliminating pixels outside of
% the RF
% Cone noise, gun weights
weights = reshape(paramsout_cone_gun(5:end),nstixperside, nstixperside, 3,length(whichframe));
[~,~,v] = svd(reshape(weights,nstixperside^2,3));
conenoise_gunweights = v(:,1)*sign(v(1,1));
% Cone noise, cone weights
weights = reshape(paramsout_cone_cone(5:end),nstixperside, nstixperside, 3,length(whichframe));
[~,~,v] = svd(reshape(weights,nstixperside^2,3));
conenoise_coneweights = v(:,1)*sign(v(1,1));

% Gun noise, gun weights
weights = reshape(paramsout_gun_gun(5:end),nstixperside, nstixperside, 3,length(whichframe));
[~,~,v] = svd(reshape(weights,nstixperside^2,3));
gunnoise_gunweights = v(:,1)*sign(v(1,1));
% Gun noise, cone weights
weights = reshape(paramsout_gun_cone(5:end),nstixperside, nstixperside, 3,length(whichframe));
[~,~,v] = svd(reshape(weights,nstixperside^2,3));
gunnoise_coneweights = v(:,1)*sign(v(1,1));

figure; 
subplot(2,1,1);
bar([conenoise_gunweights, gunnoise_gunweights])
subplot(2,1,2);
bar([conenoise_coneweights, gunnoise_coneweights])
