%%
% Section 1
% Plotting the spectrum of a just barely detectable L-M light (and the
% spectrum of the background)

stro = nex2stro(findfile('K040609001'));
bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b];
guns = reshape(stro.sum.exptParams.mon_spect,81,3);
x = linspace(380,780,81);
bkgndspect = guns*bkgndrgb'
M = reshape(stro.sum.exptParams.m_mtx,3,3);
bkgndlms = M*bkgndrgb';
threshlms = bkgndlms.*(1+[.01 -.01 0]');
(threshlms-bkgndlms)./bkgndlms  % just a sanity check
threshrgb = inv(M)*threshlms;
threshspect = guns*threshrgb;

figure; axes; hold on;
plot(x,bkgndspect);
plot(x,threshspect,'k-');
xlabel('nm');
ylabel('w/sr/m^2');

[thresholds, colorDirs, sfs, QuestTrajectories] = DTquestUnpackGH(stro, 'mode');
thresholdlms = bkgndlms


%%
% Section 2
% Continuing from previous section.  Stealing code from Psychtoolbox's
% "IsomerizationInEyeDemo.m"
whatCalc = 'LivingHumanFovea';
photoreceptors = DefaultPhotoreceptors(whatCalc);
photoreceptors = FillInPhotoreceptors(photoreceptors);
% GH converting watts/sr-m^2-wlinterval to  watts/sr-m^2
radianceWatts = SplineSpd([380 5 81],bkgndspect,[380 1 401]);
load T_xyz1931	
T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,[380 1 401]);
theXYZ = T_xyz*radianceWatts; theLuminance = theXYZ(2);
[nil,pupilAreaMM] = PupilDiameterFromLum(theLuminance,photoreceptors.pupilDiameter.source);
% GH This seems awfully big: 7.5 mm

radianceWatts = SplineSpd([380 5 81],bkgndspect,[380 1 401]);
irradianceWatts = RadianceToRetIrradiance(radianceWatts,[380 1 401], ...
			pupilAreaMM,photoreceptors.eyeLengthMM.value);
irradianceQuanta = EnergyToQuanta([380 1 401],irradianceWatts);
figure; set(gcf,'Position',[100 400 700 300]);
subplot(1,2,1); hold on
set(plot(SToWls([380 1 401]),irradianceQuanta,'r'),'LineWidth',1);
set(title('Light Spectrum'),'FontSize',14);
set(xlabel('Wavelength (nm)'),'FontSize',12);
set(ylabel('Quanta/sec-um^2-wlinterval'),'FontSize',12);
[isoPerConeSecBkgnd,absPerConeSec,photoreceptors] = ...
	RetIrradianceToIsoRecSec(irradianceWatts,[380 1 401],photoreceptors);

radianceWatts = SplineSpd([380 5 81],threshspect,[380 1 401]);
irradianceWatts = RadianceToRetIrradiance(radianceWatts,[380 1 401], ...
			pupilAreaMM,photoreceptors.eyeLengthMM.value);
irradianceQuanta = EnergyToQuanta([380 1 401],irradianceWatts);
subplot(1,2,1); hold on
set(plot(SToWls([380 1 401]),irradianceQuanta,'k'),'LineWidth',1);
[isoPerConeSecThresh,absPerConeSec,photoreceptors] = ...
	RetIrradianceToIsoRecSec(irradianceWatts,[380 1 401],photoreceptors);
%GH above line: what's the difference between the first two output
%arguments?  It depends on the difference between
%photoreceptors.effectiveAbsorbtance and photoreceptors.isomerizationAbsorbtance
subplot(1,2,2);
bar([isoPerConeSecBkgnd, isoPerConeSecThresh])
ylabel('Photoisomerizations/cone/sec');

%%
% Section 3
% For Juan Anguerya's general exam.  Looking at how receptoral and
% post-receptoral noise affect detection contours.
darknoise = 1;

sigmu = 1;
x = [-10:.1:10];
noi = normpdf(x,0,1);
sig = normpdf(x,sigmu,1);
plot(x,noi,'b-')
hold on
plot(x,sig,'g-');
% percent correct in seen/not seen with d' = 1 and optimal criterion
normcdf(sigmu/2,0,1)
% Looks like 69% correct with d'=1 and 84% with d'=2
[X,Y] = meshgrid(x,x);
% l_out = normcdf(X,1,receptornoise);
% m_out = normcdf(Y,1,receptornoise);
% lm_both = l_out .* m_out;
% % P(A or B) = P(A)+P(B)-P(A and B)
% imagesc(l_out+m_out-(l_out .* m_out))
% axis xy

%contour(X,Y,data,3)

%%
% OK, just brute forcing a simulation

x = [-10:.1:10];
l_noise = 1;
m_noise = 1;
lumnoise = 3;
rgnoise = 0;
lumgain = 1;
rggain = 4;
[X,Y] = meshgrid(x,x);
niter = 1000;
data = zeros(size(X));
noise = zeros(1,niter);
s_lum = sqrt(lumgain^2*l_noise^2+m_noise^2+lumnoise^2);
s_rg = sqrt(rggain^2*2*l_noise^2+m_noise^2+rgnoise^2);
% Above equations checked out empirically 9/23/10

tmp = [];
for i = 1:niter
    l = X+normrnd(0,l_noise,size(X));
    m = Y+normrnd(0,m_noise,size(Y));

    lum = lumgain*(l+m)+normrnd(0,lumnoise,size(X));
    rg = rggain*(l-m)+normrnd(0,rgnoise,size(X));
    det = abs(lum) > s_lum | abs(rg) > s_rg;
    %det = (abs(lum)./s_lum).^2+(abs(rg)./s_rg).^2 > 4;  % arbitrary
    data = data+det;
    
    tmp = [tmp; lum(50,50), rg(50,50)];
end

% predicted correlation
denom = sqrt(((l_noise^2+m_noise^2)+(rgnoise^2/rggain^2))*((l_noise^2+m_noise^2)+(lumnoise^2/lumgain^2)));
%(l_noise^2-m_noise^2)/(l_noise^2+m_noise^2)
(l_noise^2-m_noise^2)/denom

plot(tmp(:,1),tmp(:,2),'k.')
corr(tmp)

figure;
imagesc(data>.82*niter)
axis xy;
colormap(gray);
axis square

%%
% Section 4
% Trying to recapitulate the model that proposed by Petri and Fred in their
% manuscript "Nonlinear integration of single photon responses in the inner
% retina sets absolute visual threshold". Evidently, a threshold
% nonlinearity should be sufficient to cause two points of inflection in
% the mean output curve (as input is increased).

x = [0:.01:20];
y = [];
THRESH = 5;
ints = [0:1:4*x(end)];
L = ints < THRESH;
for i = 1:length(x)
    tmppdf = poisspdf(ints,x(i));
    zeroresp = sum(tmppdf(L));
    tmppdf(L) = 0;
    tmppdf(1) = zeroresp
    y(i) = mean(ints*tmppdf');
end
plot(x,y,'k.');
set(gca,'Yscale','log','Xscale','log')
%%
% Section 5
% Making a spinning isosamp data set.

stros = load('catted stros.mat');
stro = stros.stro1;  % Picking one example data set. 3 cpd, Sedna?
trialFields = stro.sum.trialFields(1,:);
stimon_t = stro.trial(:,strcmp(trialFields, 'stimon_t'));
stimoff_t = stro.trial(:,strcmp(trialFields, 'stimoff_t'));
anlgstart_t = stro.ras(:,strcmp(stro.sum.rasterCells, 'anlgStartTime'));
spikes = stro.ras(:,strcmp(stro.sum.rasterCells, 'sig001a'));
sf = stro.trial(:,strcmp(trialFields, 'sf'));
num_stims = min(stro.sum.exptParams.num_stims, size(stro.trial, 1));

lms = [stro.trial(:,strcmp(trialFields, 'stim_l')) ...
    stro.trial(:,strcmp(trialFields, 'stim_m')) ...
    stro.trial(:,strcmp(trialFields, 'stim_s'))]; 

[unqLms,~,unqIdxs] = unique(lms, 'rows');

num_trials = size(spikes, 1);
theseSpikes = cell(size(spikes));

stim_dur = mode(stimoff_t - stimon_t);
offset = [-0.1 stim_dur + 0.1]; % 100 ms window

for trialIdx = 1:num_trials
    tempSpikes = spikes{trialIdx} - stimon_t(trialIdx);
    tempSpikes = tempSpikes(tempSpikes >= offset(1) ...
        & tempSpikes <= offset(2));
    theseSpikes{trialIdx} = tempSpikes;
end

for stimIdx = 1:num_stims
    Lstims = unqIdxs == stimIdx;
    if sum(Lstims) > 0  % number of trials with this stimulus
        spikeTimes = cat(1, theseSpikes{Lstims});
        stim_resp(stimIdx) = numel(spikeTimes)./sum(Lstims);
    end
end


figure(); axes(); hold on;
cmap = jet(64);
for i = 1:num_stims
    plot3(unqLms(i,1), unqLms(i,2), unqLms(i,3), 'o', 'markersize', 5, ...
        'markerfacecolor', cmap(ceil(63*stim_resp(i)/max(stim_resp))+1,:), ...
        'markeredgecolor', 'none');
end
xlabel('L'); ylabel('M'); zlabel('S');
set(gca,'View',[-46 22]);
set(gca,'Xlim',[-4.2772 4.2772]/100,'Ylim',[-4.2772 4.2772]/100,'Zlim',[-42.772 42.772]/100)
axis vis3d;

viewangles = [0:3:520]-46;
viewangles(end) = [];

clear M;
for i = 1:length(viewangles)
    set(gca,'View',[viewangles(i) 22])
    M(i) = getframe(gcf);
end
repeat = 1;     %default = 1
pSearch = 1;    %default = 0
bSearch = 1;    %default = 1
reference = 1;  %default = 0
pixRange = 10;  %default = 10
iFrame = 8;     %default = 8
pFrame = 10;    %default = 10
bFrame = 25;    %default = 25
options = [repeat, pSearch, bSearch, reference, pixRange, iFrame, pFrame, bFrame];
%movie2avi(M,'IsoSampMovie.mpg'); % Doesn't seem to work
if (ispc)
    mpgwrite(M, gray, 'IsoSampMovie.mpg', options);
end

max(stim_resp)./(offset(2)-offset(1)) % maximum firing rate

%%
% making a spinning DTNT movie to go with the IsoSamp movie, above.
% Use Kali, 4 cpd
% most of this was taken from DTNT.m
zacksdata = load('IsoSamp/isosurf_params.mat');
CHEATINGPARAMS = zacksdata.kali.model(:,end);
PLOTSURFACE = 1;
MAKEMOVIE = 0;

[filelist,filelistpath,~] = uigetfile('N:\NexFiles\*.txt','Please select the file list');
if isequal(filelist,0), return; end

% parse every line while ignoring comments
fid = fopen([filelistpath filesep filelist]);
listData = textscan(fid,'%s','CommentStyle','%','Delimiter','\n');
% NOTE: textscan will ignore lines with a comment character (even if it's after the filename)!
fclose(fid);
listData = strtrim(listData{1});

badEntries = cellfun('isempty',...
    regexp(listData,'(^[A-Z]{1,2}(0[1-9]|1[012])(0[1-9]|[12][0-9]|3[01])\d{5}(\.\d)?(\.nex)?$)|(^sf:)'));

blankLines = cellfun('length',listData) == 0;
if sum(badEntries & ~blankLines) ~= 0
    disp('Warning -- the following lines will be ignored:');
    disp(listData(badEntries & ~blankLines));
end
listData(badEntries | blankLines) = [];

sfHeadingIdxs = find(~cellfun('isempty',strfind(listData,'sf:')));
if ~isempty(sfHeadingIdxs)
    listSFs = textscan([listData{sfHeadingIdxs}],'sf:%f');
    listSFs = listSFs{1};
    
    % sort file names by SF
    filenamesBySF = cell(length(listSFs),2);
    for i = 1:length(listSFs)
        if i ~= length(listSFs)
            fns = listData(sfHeadingIdxs(i)+1:sfHeadingIdxs(i+1)-1);
        else
            fns = listData(sfHeadingIdxs(i)+1:end);
        end
        filenamesBySF{i,1} = fns;
        filenamesBySF{i,2} = listSFs(i);
    end
end
% sort list by sf
sfNums = cell2mat(filenamesBySF(:,2));
[sfNums,sortedOrder] = sort(sfNums);
filenamesBySF = filenamesBySF(sortedOrder,:);

% find the indices of the dupes (based on repval from Mathworks FEX)
[B,~,J] = unique(sfNums);
I = 1:length(J);
I2 = find(diff(J) == 0);
B2 = unique([I2(:)' I2(:)'+1]);
[~,~,IR] = unique(J(B2));
if ~isempty(IR)
    POS = I(B2);
    idxsToDelete = [];
    for i = unique(IR)'
        idxsDupes = POS(IR == i);
        temp = filenamesBySF(idxsDupes,1);
        filenamesBySF{idxsDupes(1),1} = vertcat(temp{:});
        idxsToDelete = [idxsToDelete; idxsDupes(2:end)];
    end
    filenamesBySF(idxsToDelete,:) = []; % trim them out
end

LINPREDTOL = 0.3;
INITCONTRASTFACTOR = 4; % multiplier for initial contrast (not starting in the plane of 
% the parent triangle, otherwise contrast would be too low).
[sfidx,okay] = listdlg('PromptString','Select a spatial frequency:',...
    'SelectionMode','single','ListString',cellstr(num2str([filenamesBySF{:,2}]')));
if ~okay, return; end

clear trialspecs;
load ('T_cones_smj10');
load ('T_cones_smj');

% Loading a list of data files, one by one
assumedDataPath = '';
filenames = filenamesBySF{sfidx,1};
for fileidx = 1:size(filenames,1)
    if fileidx == 1
        firstFile = findfile(char(filenames{fileidx}));
        if isempty(firstFile)
            warning('DTNT:fnf','File ''%s'' wasn''t found',char(filenames{fileidx}));
            continue
        else
            [assumedDataPath,~,~] = fileparts(firstFile);
        end
        currFile = firstFile;
    else
        nextFile = [assumedDataPath filesep char(filenames{fileidx})];
        if ~exist(nextFile,'file')
            nextFile = findfile(char(filenames{fileidx}));
            if isempty(nextFile)
                warning('DTNT:fnf','File ''%s'' wasn''t found',char(filenames{fileidx}));
                continue
            else
                [assumedDataPath,~,~] = fileparts(nextFile);
            end
        end
        currFile = nextFile;
    end
    stro = nex2stro(currFile);
    [thresholds, colorDirs, sfs, QuestTrajectories] = DTquestUnpackGH(stro, 'mode');
    figure; set(gcf,'name',char(filenames{fileidx}));
    for i = 1:length(thresholds)
        subplot(2,2,i);
        plot(QuestTrajectories{i},'k.');
        title(num2str(colorDirs(i,:)./norm(colorDirs(i,:))));
    end

    threshs = thresholds';
    lmsmat = mkbasis(colorDirs')';
    
    mon_spd = stro.sum.exptParams.mon_spect;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    M = reshape(stro.sum.exptParams.m_mtx,[3 3]);
    
    % Each row is a color direction for consistency with NeuroThreshOnline.m
    if (fileidx == 1)
        % Making entries into the trialspecs array: round 1, (3 color dirs)
        for i = 1:size(lmsmat,1)
            trialspecs(i).colordir = lmsmat(i,:)./norm(lmsmat(i,:));
            trialspecs(i).predictedthreshold = [];
            trialspecs(i).parentvertices = [];
            trialspecs(i).measuredthreshold = threshs(i);
            trialspecs(i).M = M;
            trialspecs(i).bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b];
        end
        
        % Making entries into the trialspecs array: round 2, (4 color dirs)
        signmat = [(fullfact([2,2])-1.5)*2 ones(4,1)];  % assumes three color dirs
        for i = 1:size(signmat,1)
            idx = length(trialspecs)+1;
            parentvertices = diag(signmat(i,:).*threshs)*lmsmat;
            v = mean(parentvertices);
            trialspecs(idx).colordir = v./norm(v);
            trialspecs(idx).predictedthreshold = norm(v);
            trialspecs(idx).parentvertices = parentvertices;
            trialspecs(idx).measuredthreshold = [];
            trialspecs(idx).M = nan*ones(3);
            trialspecs(i).bkgndrgb = [nan nan nan];
        end
    else  % end first data file (setting up round 1 and 2)
        
        colordirs = reshape([trialspecs.colordir],3,length(trialspecs))';
        for i = 1:size(lmsmat,1)
            L = logical(softEq((lmsmat(i,:)./norm(lmsmat(i,:)))*colordirs', 1,10));
            % Error checking
            if (sum(L) == 0)
                disp(['cannot find color direction in trialspecs: ',num2str(lmsmat(i,:))]);
                continue
            end
            if (sum(L) > 1)
                error('Multiple identical color directions detected.');
                keyboard;
            end
            if (~isempty(trialspecs(L).measuredthreshold))
                error('Already measured a threshold in this direction');
            end
            trialspecs(L).measuredthreshold = threshs(i);
            threshratio = trialspecs(L).measuredthreshold./trialspecs(L).predictedthreshold;
            if (abs(log(threshratio)) > abs(log(1+LINPREDTOL)))
                grandparentvertices = trialspecs(L).parentvertices;
                for j = 1:3
                    idx = length(trialspecs)+1;
                    parentvertices = grandparentvertices;
                    parentvertices(j,:) = trialspecs(L).colordir.*trialspecs(L).measuredthreshold;
                    v = mean(parentvertices);
                    trialspecs(idx).colordir = v./norm(v);
                    trialspecs(idx).predictedthreshold = norm(v);
                    trialspecs(idx).parentvertices = parentvertices;
                    trialspecs(idx).measuredthreshold = [];
                end
            else
                disp(['File #',num2str(fileidx),': Threshold in direction ',num2str(lmsmat(i,:)),' is consistent with linear prediction']);
            end
            if (abs(log(threshratio)) > log(3))
                fprintf('Warning: aberrant threshold point: %s\n',[char(filenames{fileidx}),' (',num2str(trialspecs(L).colordir),')']);
            end
        end
    end
end % Datafile loop

% Plotting
scaled = [];
for i = 1:length(trialspecs)
    if (~isempty(trialspecs(i).measuredthreshold))
        scaled = [scaled; trialspecs(i).colordir.*trialspecs(i).measuredthreshold];
    end
end

figure; axes; hold on;
opengl('software');  % To avoid getframe bug: http://www.mathworks.com/support/bugreports/384622
plot3(scaled(:,1),scaled(:,2),scaled(:,3),'k.')
plot3(-scaled(:,1),-scaled(:,2),-scaled(:,3),'k.')

if (PLOTSURFACE)
    % Plotting the fit
    plotlims = max(abs(scaled));
    [x,y,z] = meshgrid(linspace(-plotlims(1),plotlims(1),150),linspace(-plotlims(2),plotlims(2),150),linspace(-plotlims(3),plotlims(3),50));
    tmp = [x(:) y(:) z(:)];
    %   v = sum(abs(tmp *reshape(initparams(2:end),3,3)).^initparams(1),2);
    fpar = CHEATINGPARAMS;
    v = sum(abs(tmp *reshape(fpar(2:end),3,3)).^fpar(1),2);
    
    fv = isosurface(x,y,z,reshape(v,size(x,1),size(x,2),size(x,3)),1);
    h = patch(fv,'EdgeAlpha',0,'FaceAlpha',.2);
end
xlabel('\DeltaL/L'); ylabel('\DeltaM/M'); zlabel('\DeltaS/S');

set(gca,'View',[-46 22]);
set(gca,'Xlim',[ -4.3460 4.4135],'Ylim',[ -3.3664 3.6412],'Zlim',[-33.8701 36.2056]);
axis vis3d;

viewangles = [0:3:520]-46;
viewangles(end) = [];

if (MAKEMOVIE)
    clear M;
    for i = 1:length(viewangles)
        set(gca,'View',[viewangles(i) 22])
        M(i) = getframe(gcf);
    end
    repeat = 1;     %default = 1
    pSearch = 1;    %default = 0
    bSearch = 1;    %default = 1
    reference = 1;  %default = 0
    pixRange = 10;  %default = 10
    iFrame = 8;     %default = 8
    pFrame = 10;    %default = 10
    bFrame = 25;    %default = 25
    options = [repeat, pSearch, bSearch, reference, pixRange, iFrame, pFrame, bFrame];
    %movie2avi(M,'IsoSampMovie.mpg'); % Doesn't seem to work
    if (ispc)
        mpgwrite(M, gray, 'IsoSampMovie3.mpg', options);
    end
end

%%
% Simulation to look at the impact of using sqrt(a^2+b^2) to discriminate
% 2-D homoskedastic Gaussian distribitions of a and b for which mua = mub+c 
% so a+b is the optimal discriminant function (integrating joint densities 
% over lines of constant a-b).

% a and b represent inputs from the two eyes (which are scattered around 0
% in the noise condition and are scattered around 1 in the signal
% condition)

% dprime of 1.25 is approximately threshold
n = 10000; % Just to confirm analytical distns
munoise = 2*sqrt(2); % Just trying to keep the noise distribution far from 0
% It doesn't make sense for one eye to eve give a negative signal (that
% then got mistaken for a postive signal because of the squaring).
musignal = munoise+1.25;
sigma = 1; % Sigma must = 1 for noncentral chi to work.
noise = normrnd(munoise,sigma,n,2);
signal = normrnd(musignal,sigma,n,2);
PLOTHISTS = 0;

L = [false(size(noise,1),1); true(size(signal,1),1)];

% ---------------------
% The first way: adding
% ---------------------
stat = [sum(noise,2); sum(signal,2)];
bins = linspace(min(stat),max(stat),30);

figure; axes; hold on;
[counts,~] = hist(stat(~L), bins);
[counts,~] = hist(stat(L), bins);
if (PLOTHISTS)
    bar(bins,counts./sum(counts),'r');
    bar(bins,counts./sum(counts),'k');
end
% Gaussian distributions
x = linspace(bins(1),bins(end),200);
normfact = (bins(2)-bins(1))/(x(2)-x(1));
noisepdf = normpdf(x,2*munoise, sigma*sqrt(2));
signalpdf = normpdf(x,2*musignal,sigma*sqrt(2));
plot(x,noisepdf./sum(noisepdf).*normfact,'k-','LineWidth',2);
plot(x,signalpdf./sum(signalpdf).*normfact,'r-','LineWidth',2);

% Getting the probability of error
crit = 2*mean([munoise musignal]);
FAs1 = 1-normcdf(crit,2*munoise, sigma*sqrt(2)) % FAs
misses1 = normcdf(crit,2*musignal, sigma*sqrt(2)) % misses
errorrate1= [.5 .5]*[FAs1; misses1];
title(['a+b: ',num2str(1-errorrate1)]);

% ------------------------
% The second way: distance
% ------------------------
%stat = [sum(noise.^2,2); sum(signal.^2,2)]; % losing the sqrt
%stat = [sum(noise.^2,2)./sigma^2; sum(signal.^2,2)./sigma^2]; % losing the sqrt and dividing by sigma
stat = [sqrt(sum(noise.^2,2)); sqrt(sum(signal.^2,2))]; % Still haven't
%figured out how to work with this one

bins = linspace(min(stat),max(stat),30);
figure; axes; hold on;
[counts,~] = hist(stat(~L), bins);
[counts,~] = hist(stat(L), bins);
if (PLOTHISTS)
    bar(bins,counts./sum(counts),'r');
    bar(bins,counts./sum(counts),'k')
end
% Non-central chi-squared distribution (the analytical solution)
x = linspace(bins(1),bins(end),200);
%x = (linspace(bins(1),bins(end),200)).^2;
normfact = (bins(2)-bins(1))/(x(2)-x(1)); % for getting the height correct
noiselambda = 2*(munoise/sigma).^2;
%noisepdf = ncx2pdf(x,2,noiselambda);
noisepdf = 2.*x.*ncx2pdf((x/sigma^2).^2,2,noiselambda);

signallambda = 2*(musignal/sigma).^2;
%signalpdf = ncx2pdf(x,2,signallambda);
signalpdf = 2.*x.*ncx2pdf((x./sigma^2).^2,2,signallambda);

plot(x,noisepdf./sum(noisepdf).*normfact,'k-','LineWidth',2);
plot(x,signalpdf./sum(signalpdf).*normfact,'r-','LineWidth',2);
errorrate2= [.5 .5]*[FAs2; misses2];
title(['sqrt(a^2+b^2)',num2str(1-errorrate2)]);

% Getting the probability of error
% Using the noncentral chi-squared distribution instead of the noncentral
% chi distribution just because I'm more comfortable with it and the sqrt
% doesn't change the efficiency at all.
ratio = log10(noisepdf./signalpdf);
tmp = ratio.^2;
critguess = x(tmp == min(tmp)).^2;

f = @(x) (log10(ncx2pdf(x,2,noiselambda)./ncx2pdf(x,2,signallambda))).^2; % numerical optimization
crit = fminsearch(@(x) f(x),critguess);

FAs2 = 1-ncx2cdf(crit,2,noiselambda) % FAs
misses2 = ncx2cdf(crit,2,signallambda) % misses


Efficiency = errorrate1/errorrate2
%
% OK, now making a reasonable-looking figure
figure; axes; hold on;
theta = linspace(0,2*pi,200)';
dists = sqrt(chi2inv(linspace(0, 1, 10),2))
for i = 1:length(dists)
    plot(munoise+dists(i)*cos(theta), munoise+dists(i)*sin(theta),'k-','linewidth',2);
end
for i = 1:length(dists)
    plot(musignal+dists(i)*cos(theta), musignal+dists(i)*sin(theta),'r-','linewidth',2);
end
axis equal
set(gca,'Ylim',[0,musignal+3*sigma],'Xlim',[0,musignal+3*sigma]);
%legend({'noise','signal'});
xlabel('Left eye response');
ylabel('Right eye response');

% Ideal observer1
tmp = [1:2*(musignal+3*sigma)];
for i = 1:length(tmp)
    plot([0 tmp(i)],[tmp(i) 0],'k--')
end

% Ideal observer2
tmp = [1:sqrt(2)*(musignal+3*sigma)];
for i = 1:length(tmp)
    plot(tmp(i)*cos(theta),tmp(i)*sin(theta),'k--')
end

%%
% Brainstorming about information loss in midget ganglion cells
% Time and space need to be treated differently (why?)

% Temporal factors:
% (This choice of stimulus may have biased monkey performance against luminance 
% detection by being suboptimal for M cells. Using stimuli with more abrupt
% transitions would presumably increase the importance of the M-pathway).
% Even midget cells go up to higher frequencies than humans can see. How
% where the temporal filters/noise sources are that account for this is
% unknown. Seems that the net SNR for lum and RG are matched behaviorally
% at ~15 Hz (conditioned on the spatial and spectral aspects of the
% stimulus). This probably does not indicate similar SNR in individual
% cells becuase there are more parvo than magnocells. However, how magno
% and parvo cells contribute to the detection of 15 Hz luminance stimuli is
% unclear (refs).

% Spatial factors:
% At the eccentricities we're recording at we think there is ~1 cone/RF
% center. First considering the noiseless GC case.
% If there are equal numbers of GCs and cones, there need not be any loss
% of information. 
% Center/surround organization makes bandpass filtering occur. 
% Some frequencies are reduced to zero (DC)? One can contrive a way to reduce
% some frequencies to zero (e.g. by equally weighting center and surround).
% If there are fewer GCs than cones then there must be a loss
% of information, but this could be a loss of information about stimuli that never
% occur or are behaviorally unimportant (very high frequency stimuli).
% If we add noise after the GCs, then we csn also eliminate information
% that does not get filtered close to but above zero. This surely happens -
% noiseless signal transformation is an idealization. 

%%
% Stuff for Jacob Baudin's general exam
% parameteric formula for an ellipse from Wikipedia:
% X(t)=X_c + a*cos(t)*cos(phi) - b*sin(t)*sin(phi)
% Y(t)=Y_c + a*cos(t)*sin(phi) - b*sin(t)*cos(phi)
% which is where Equation A2 in Lee Mollon Zaidi and Smithson (2009) comes
% from.

theta = 0; % Phase difference between the two components (of the stimulus, which is aligned to the cardinal axes).
phi = 0.25*pi; % Phase delay in the S-cones
t=0:.01:4*pi;

% For recaptiulating Figure 2A
figure; axes; hold on;
plot(sin(t),sin(t-theta+phi),'b.')
plot(sin(t),sin(t-theta+phi-pi),'r.')
axis square

r_cw = sqrt(sin(t).^2+sin(t-theta+phi).^2);
r_ccw = sqrt(sin(t).^2+sin(t-theta+phi-pi).^2);
figure; axes; hold on;
plot(t,r_cw,'b.');
plot(t,r_ccw,'r.');

% r_cw and r_ccw are exactly the same?!

dr_cw = (sin(t).^2+sin(t-theta+phi).^2).^-.5.*(sin(t).*cos(t)+sin(t-theta+phi).*cos(t-theta+phi));
dr_ccw = (sin(t).^2+sin(t-theta+phi-pi).^2).^-.5.*(sin(t).*cos(t)+sin(t-theta+phi-pi).*cos(t-theta+phi-pi));
figure; axes; hold on;
plot(t,dr_cw,'b.');
plot(t,dr_ccw,'r.');
dt = (t(2)-t(1)); % step size
dt_vect = t(1:end-1)+dt/2; 
plot(dt_vect, diff(r_cw)./dt); % Finite difference approximation
% sin(t) is L-M component
% sin(t-theta+phi) is S component



% Average estimated S-cone delay was (10.1+16.1)/2 = 13.1 ms 
% which is 0.25*pi and 0.23*pi radian phase shifts at 10 Hz
% Proof:
% At 10 Hz, 1 cycle (2*pi radians) is 100 ms.
% (.25*pi)/(2*pi) = 0.125
% 0.125*100 = 12.5 ms
% (.23*pi)/(2*pi) = 0.1150
% 0.1150*100 = 11.5 ms

% Question B
% Consider this hypothetical result: the direction of modulation around the hue
% circle (clockwise or counter clockwise) was hardest to perceive when the L-M 
% modulation led the S-cone modulation by 1.7 radians. How would the authors have
% interpreted this result (in terms of delays between (L- and M-cones) and S-cones)

phi = 1.7-pi/2; % phi = 0.1292
(phi/(2*pi))*100 % (phi/(2*pi))*100 = 2.05 ms; the amount S-cone signals are delayed.
phi = 1.7-3*pi/2; % phi = 0.1292
(phi/(2*pi))*100 % (phi/(2*pi))*100 = 2.05 ms; the amount S-cone signals are delayed.


% Do you think the adapting field used in the Lee et al. (2009) study resulted
% in equal cone catches among the L-, M-, and S-cones? Why or why not? 

% Why might the adapting field used in the Lee et al. (2009) study have
% resulted in fewer S-cone R*/s than M- or L-cone R*/s? What affect would you
% expect this difference to have on their results?

%%
% Matlab exercises for the neurostats course 2016 (Fred's teaching)

% I have a hard time tying this stuff to neuroscience. Feel free to adjust.

% First, some set up.
x = 1:9; % All the potential values that X can take
px = log10(1+1./x); % The probabilty that X = x(i)

% Now let's calculate the expected value of X
n = length(x); % The number of elements in x (the 'support' of the distribution)
expected_value = 0; % Initializing a variable that we'll increase iteratively
for i = 1:n  % Start of loop. First time through i=1, second time i=2,... last time i=n
   new_bit = x(i)*px(i); % Weighting x(i) by the probability that X = x(i)...
   expected_value = expected_value + new_bit; % and adding it to the 'expected_value' variable
end % End of the loop.

% Here?s a shorter piece of code that does the same thing:
x*px'
% Confirm that px sums to 1
% Modify both pieces of code to compute E(log10(X)) instead of E(X)

% --------------------------------

% You flip a coin with probability of heads = p.
% if is comes up heads, you flip the coin nine times and 
% count the number of (additional) heads, so you 
% get a random number from 0 to 9 drawn from a 
% binomial distribution. 
% If it comes up tails, you draw a number from a discrete 
% uniform distribution on [0:9].
% What do you expect to get on average?
% Let's solve this two ways.

% The first way is through Monte Carlo simulation
niter = 10000; % I hope this is enough iterations
p = .2; % Gotta pick something
randomdraws = zeros(niter,1); % We'll keep our data in this vector
for i = 1:niter
    coinflip = binornd(1,p); % O = tails, 1 = heads
    if (coinflip == 0)
        randomdraws(i) = binornd(9,p);
    else % (coinflip == 1), there are no other possibilties
        randomdraws(i) = unidrnd(10)-1; % stupid Matlab '1' indexing       
    end
end
mean(randomdraws)

% Now we have to do this for many values of 'p' and keep track the results. 
% (Are the students' Matlab skills up to this task?)

% Here's an alternative (faster and more accurate) route. 
% For each value of p, we'll calulate the PMF and calculate the expected
% value from that.
p = .2; % Gotta pick something
conditionals = zeros(10,2); 
conditionals(:,1) = binopdf([0:9],9,p); % These are lines of code that maybe they can come up with?
conditionals(:,2) = .1; % I'm worried that if we just give them, they won't understand why they work.

pmf = conditionals.*repmat([1-p, p], 10,1);
% Now calculating the expected value
[0:9, 0:9]*pmf(:)

% Now do this for 10000 equally spaced values of p using whichever approach you prefer
all_the_ps = linspace(0,1,1000);
expectedvalues = zeros(length(all_the_ps),1);
for p_idx = 1:length(all_the_ps)
    p = all_the_ps(p_idx);
    conditionals(:,1) = binopdf([0:9],9,p);
    pmf = conditionals.*repmat([1-p, p],10,1);
    expectedvalues(p_idx) = [0:9, 0:9]*pmf(:);
end
plot(all_the_ps,expectedvalues)

% What should you make p if you want this number to be as high as possible?
% If you get a 2, what would be a reasonable guess for p?
% --------------------------------------------------------
% You're doing a flash detection experiment. On a fraction 'p' of the
% trials a flash is shown. On the other fraction of the trials (1-p) no
% flash is shown. 

% When a flash is shown you get a response from a Possion distribution with
% parameter lambda_1. When NO flash is shown you get a response from a 
% Possion distribution with parameter lambda_0.

% Let's say lambda_1 = 3 and lambda_0 = 1;
% First let's look at the pdfs
x = 0:20; % Gotta say something about why we're cutting it off at 20
lambda_1 = 3;
lambda_0 = 1;
pdf1 = poisspdf(x,lambda_1);
pdf0 = poisspdf(x,lambda_0);
figure; axes; hold on;
plot(x,pdf0,'ko-');
plot(x,pdf1,'ro-');

% Let's say you observed a response of 2. Is that more likely to
% be a draw from the flash distribution or the non-flash distribution?

% Let's say the probability of a non-flash is 4x greater than the
% probability of a flash. What is the joint distribution of the response?
% First we make a matrix with the condition distributions in it.
conditionals = [pdf0;pdf1];
% Now we have to multiply the first row by the probability of no flash
% and the second row by the probability of a flash to turn it into a 
% real joint distribution that integrates to 1
% (They should be able to figure out that p_noflash = .8 and p_flash = .2)
pmf = [];
pmf(1,:) = conditionals(1,:) * p_noflash;
pmf(2,:) = conditionals(2,:) * p_flash;

% Let's say you observed a response of 2. 
% What is P(flash|response = 2)? What is P(no_flash|response = 2)?
% Given your answer to the above questions, what should you say 
% when you get a response of 2, "flash" or "no-flash"? What percentage of the
% time will you be correct?

% ---------------------------------------------

% You observe a number of spike counts that are assumed to be draws from a
% Poisson distribution with parameter, lambda.
% They are:
data = [16    17    12     9    13     8    15    10    11    16];

% We want to estimate lambda. We could just calculate xbar.
mean(data)

% Or if you wanted to d it the long way you could type:
sum(data)/length(data)

% This is actually a good estimate (I think the value I actually picked  
% to generate these data was 13). Now let's say that all the values ?10
% got rounded down to zero.

newdata = [16    17    12     0    13     0    15    0    11    16];

% We might still want to estimate lambda, but it's definitely not the mean
% anymore, which is:

mean(newdata)

% It's much lower than before because we've taken a bunch of positive
% numbers and replaced them with zeros. Nevertheless, I generated these
% data by picking a particular value for lambda. It's still soemthing we
% might want to estimate, but mean(x) is going to be biased to small
% values and not a great choice for an estimator. Can we do better?

% We're going to estimate lambda by a brute force application of the method
% of maximum likelihood. First, we're going to do it for the regular data
% (that are draws from a regular, non-truncated Poisson distribution).
% Then, again for the truncated Poisson distribution.

% First just plotting a Poisson distribution on top of the data so you can
% see what we're talking about.
my_lambda_guess = mean(data);
x = 0:30;
y = poisspdf(x,my_lambda_guess);
figure; axes; hold on;
plot(x,y);
plot(data,0,'k*')

% You can think of this Poisson distribition that I plotted as a model
% describing where the data might have come from. As such, we can query the
% model to find out how likely different draws are under the model.
% What's the probability we get a 16 from this distribution?
y(x==16)

% What's the probability we get a 9 from this distribution?
y(x==9)

% What is the probaility that we get the data vector we got from this
% distribution?
% If we assume that all of our data points are independent from each other
% then their joint probability (the probability of getting all of these particular
% draws)  all is just the product of the individual probabilities:

lik = 1;
for i = 1:length(data)
    lik = lik*y(x==data(i));
end
lik

% That's a tiny number. But we shouldn't be too surprised by that. The
% probability of getting *any* particular set of 10 numbers is very small,
% even when eaach the 10 numbers is relatively probable. The reason is that
% there are many, many sets of 10 numbers that you could draw from any
% single Poisson distrubition. If you sum up all of these probabilities,
% across all possible sets of 10 numbers, you'd better get '1'. So each set
% of 10 numbers gets a tiny amount of probability so that across all
% possible sets of 10 numbers the probailities sum to 1.

% That was a bit of a digression. The point is that we have a way of
% asking, for any model we dream up, what is the probability of getting
% data like ours under this model. Now, we can try a bunch of candidate
% models (a bunch of lambdas) and see which describes our data the best.

lambdas = linspace(10,15,100);
lliks = [];
for j = 1:length(lambdas)
    llik = 1;
    y = poisspdf(x,lambdas(j));
    for i = 1:length(data)
        llik = llik*y(x == data(i)); 
    end
    lliks(j) = llik;
end
plot(lambdas,lliks);
xlabel('lambda'); ylabel('likelihood');
% Now, let's find the peak of this curve - the value of lambda that
% maximizes the probaility of getting the data set we got:
lambdas(lliks == max(lliks))

% The answer is mean(data) (within quantization error). So you can look at
% what we just did as a long, labor-intensive way of doing something that 
% you already knew how to do in an easier way. But wait, there's more.

% This technique of maximum-likelihood estimation works in all kinds of
% situations, situations where an easy shortcut to the answer may not exist.
% For example, we can use it to estimate the parameter of a truncated
% Poisson distribution. To do thins, you can use the code I wrote above with
% a tweak; you have to manipulate 'y' (the Poisson PMF) so that all the
% probability mass that's usually on x = 1 to x = 10  gets assigned to
% x = 0. I'll let you figure out how to do this by yourself.

lambdas = linspace(10,15,100);
lliks = [];
for j = 1:length(lambdas)
    llik = 1;
    y = poisspdf(x,lambdas(j));
    % An important line of code goes here
    y(1) = y(1)+sum(y(2:11)); % This line redacted in actual assignment
    for i = 1:length(data)
        llik = llik*y(x == data(i)); 
    end
    lliks(j) = llik;
end
plot(lambdas,lliks);
xlabel('lambda'); ylabel('likelihood');

% What's the maximum likelihood estimate for lambda?
