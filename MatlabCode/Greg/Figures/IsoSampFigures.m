% IsoSamp figures for a Rieke lab meeting (and NYU talk and Tuebingen talk and Rochester talk)
% For a movie of a Gabor see LMTF_figures, section 3.
% For a pan color isoresponse surface movie see SFN2010.m, section 12
% For cone model isoresponse ellipsoid see PBIO2016, section 2
% For contrast-sensitivity functions with nearby data points, see LMTF_figures, section 9
% For STAs, go to WNAnalysis.m
% For cone fundamentals, see CSHL2018.m, Section 1
%
% Section 1) Behavioral L+M (or L-M) temporal contrast-sensitivity function and a
%     theoretical one based on a counter of Poisson photon arrivals.
% Section 2) Behavioral L+M (or L-M) temporal contrast-sensitivity function and a
%     theoretical one based on cone outer segment currents. Need to run section
%     1 before running this one.
% Section 3) Cone model d' as a function of temporal frequency at behavioral
%     detection threshold (e.g. at contrast where monkey's d'=1.27).
% Section 3.1) Gaussians separated by 1.27 sd (d' = 1.27)
% Section 4) Single neuron rasters and d's (and PSTHs)
% Section 5) RF hexagonal grid on the stimulus
% Section 6) Population of neurons (population d' analysis) + cones
% Section 6.1) Like above, but using the Victor and Purpura spike train
% distance metric
% Section 7) RF positions of magno and parvoceullar neurons
% Section 8) Input efficiencies (for back of the carousel)
% Section 9) Single example of putative "S-potential"
% Section 10) Single example of signal and noise distributions (confirming
% that d' is not a crazy metric)
% Section 11) Benardete and Kaplan simulated LGN neuron
% Section 12) temporal impulse response and noise spectrum for cones
% Section 13) Rasters for optogenetic stimulation
% Section 14) Real cone data from Juan
% Section 15) Temporal integration analysis (extra figure)
% Section 16) Rotating ellipsoids so we're looking at them in the LM plane
% Section 16.5) Various spinning neurothresh data
% Section 17) A bunch of Gabors at the monkeys' detection threshold
%%
% Section 1. Example temporal contrast sensitivity function and prediction
% from an ideal observer of photon counts.
% Largely drawn from ForCharlie.m, Section 7.

RFX = 5; % DVA
RFY = 0; % DVA
PLOTDATA = 0;
THETAWEDGE = pi/7; % For plotting raw data
MMPERDEG = 0.223; % mm/deg (Perry and Cowey 1985)
STIMTYPE = 'RG';
SMALLAXES = 0; % xlim 1:30, ylim 1:1000
if strcmp(STIMTYPE,'LUM')
    maxTF = 60;
elseif strcmp(STIMTYPE,'RG')
    maxTF = 30;
end
w = linspace(1,maxTF,120);

bkgndlms_Rstar = [8898 7378 2302]; % per cone per second
% Lifted bkgndlms_Rstar values from DTcones_gh.m (line 569) for Propixx calibration
% done on July 4, 2017. Brainard's "IsomerizatonsInEyeDemo.m" gives
% [9110 6400 750] for the human eye in response to the ProPixx background.
% Pretty close.

% Loading behavioral model data
isosamppath = which('IsoSampOnline');
isosamppath(find(isosamppath==filesep,1,'last'):end) = [];
load ([isosamppath,filesep,'private',filesep','data',filesep,'LMTF.mat']);
% "A" is the structure of LMTF model parameters for Apollo
behavioralmodel = LMTF_global_to_local_model(A.legacy.mode5params, RFX, RFY, 3);
% lum
if strcmp(STIMTYPE,'LUM')
    behavioral_tcsf = @(omega)(behavioralmodel(1)*abs(((1i*2*pi*10^behavioralmodel(5).*omega+1).^-behavioralmodel(3))-behavioralmodel(2)*((1i*2*pi*10^(behavioralmodel(5)+behavioralmodel(6)).*omega+1).^-(behavioralmodel(3)+behavioralmodel(4)))));
    CENTERANGLE = atan2(1,1);
elseif strcmp(STIMTYPE,'RG')
    behavioral_tcsf = @(omega)(behavioralmodel(7)*abs(((1i*2*pi*10^behavioralmodel(11).*omega+1).^-behavioralmodel(9))-behavioralmodel(8)*((1i*2*pi*10^(behavioralmodel(11)+behavioralmodel(12)).*omega+1).^-(behavioralmodel(9)+behavioralmodel(10)))));
    CENTERANGLE = atan2(1,-1);
end
ecc = sqrt((RFX)^2 + (RFY)^2); % in deg
conespermm2 = 150.9*10^3*exp(-1.2*ecc)+35.9*10^3*exp(-.16*ecc)+9.9*10^3*exp(-.03*ecc); % actually cones/mm^2 ~27,000
Sconespermm2 = 2.5*10^3 *exp(-.2*ecc)+1.8*10^3*exp(-.05*ecc);
conespermm2 = conespermm2*2; % Two eyes
Sconespermm2 = Sconespermm2*2;
gaborsddeg = 0.15; % IsoSamp Gabor sigma in DVA
gaborsdmm = gaborsddeg*MMPERDEG;
pixelsizeindeg = 0.01; % This doesn't matter very much
gaborlim = 3; % sds
npix = round(2*gaborlim*gaborsddeg/pixelsizeindeg);
sperframe = 1/240;

pixelsizeinmm = pixelsizeindeg*MMPERDEG;
conesperpixel = pixelsizeinmm.^2*conespermm2;
Sconesperpixel = conesperpixel*(Sconespermm2/conespermm2);
Lconesperpixel = (conesperpixel-Sconesperpixel)/2;
Mconesperpixel = Lconesperpixel;

% making the stimulus/weighting function
envelope = [linspace(0,1,40) ones(1,80) linspace(1,0,40)]; % in frames
nframes = length(envelope);
sinusoid = sin(2*pi*3*linspace(0,nframes*sperframe,nframes)); % sin phase provides zero DC
temporalwtframes = envelope.*sinusoid;
spatialwtpix1D = normpdf(linspace(-gaborlim*gaborsdmm,gaborlim*gaborsdmm,npix),0,gaborsdmm);
spatialwtpix1D = spatialwtpix1D./max(spatialwtpix1D);
spatialwtpix2D = spatialwtpix1D'*spatialwtpix1D;
spacetimewt = repmat(spatialwtpix2D,[1 1 nframes]).*repmat(permute(temporalwtframes,[1 3 2]),[npix npix 1]);

% Getting the noise (no Gabor) distribution
LMScc = [0 0 0]; % cone contrast at peak
L = bkgndlms_Rstar(1)*(1+LMScc(1).*spacetimewt)*Lconesperpixel*sperframe; % L isomerizations each pixel, each frame 
M = bkgndlms_Rstar(2)*(1+LMScc(2).*spacetimewt)*Mconesperpixel*sperframe; % M isomerizations each pixel, each frame 
S = bkgndlms_Rstar(3)*(1+LMScc(3).*spacetimewt)*Sconesperpixel*sperframe; % S isomerizations each pixel, each frame 

Lmn = sum(sum(sum(L.*spacetimewt)));
Lsd = sqrt(sum(sum(sum(L.*spacetimewt.^2)))); % Poisson assumption: var(spacetimewt*X(lambda)) = spacetimewt^2*lambda
Mmn = sum(sum(sum(M.*spacetimewt)));
Msd = sqrt(sum(sum(sum(M.*spacetimewt.^2)))); 
noisedist = [Lmn Mmn; Lsd Msd];
% Finished with noise distribution

ccs = logspace(-3,-2.5,10); % cone contrast (vector norm) at peak.
data = zeros(length(ccs),1);
for i = 1:length(ccs)
    L = bkgndlms_Rstar(1)*(1+cos(pi/4)*ccs(i).*spacetimewt)*Lconesperpixel*sperframe; % L isomerizations each pixel, each frame
    M = bkgndlms_Rstar(2)*(1+sin(pi/4)*ccs(i).*spacetimewt)*Mconesperpixel*sperframe; % M isomerizations each pixel, each frame
    % sin and cos above because ccs is in vector norms, not Lcc = Mcc = cc.
    
    Lmn = sum(sum(sum(L.*spacetimewt)));
    Lsd = sqrt(sum(sum(sum(L.*spacetimewt.^2)))); % Poisson assumption: var(spacetimewt*X(lambda)) = spacetimewt^2*lambda
    Mmn = sum(sum(sum(M.*spacetimewt)));
    Msd = sqrt(sum(sum(sum(M.*spacetimewt.^2))));
    
    signaldist = [Lmn Mmn; Lsd Msd];
    
    v = (signaldist(1,:)-noisedist(1,:))./signaldist(2,:).^2; % linear discriminant vector (is cone contrast vector ?!)
    v = v./norm(v); % <--- This step is important and isn't mentioned in the text
    mu1 = signaldist(1,:)*v';
    mu2 = noisedist(1,:)*v';
    sigmasq = signaldist(2,:).^2*(v.^2)';
    data(i) = (mu1-mu2)./sqrt(sigmasq);
    
    % This is equivalent to above 
    % dprime = (signaldist(1,:)-noisedist(1,:))./signaldist(2,:); % d-prime
    % data(i) = sqrt(dprime*dprime');
end
threshold = interp1(data,ccs,1.27,'spline'); % finding contrast that produces d'=1.27 to match monkey
photon_count_sensitivity = 1/threshold;

sens = behavioral_tcsf(w);
figprefs(1); axes; hold on; set(gcf,'Color',[0 0 0],'InvertHardCopy','off');
set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]);
plot(w,sens,'y-','LineWidth',2);
if PLOTDATA
    Lecc = A.eccs(:,1) == RFX*10 & A.eccs(:,2) == RFY*10;
    rawdata = A.raw{Lecc};
    th = atan2(rawdata(:,2), rawdata(:,1));
    Loog = rawdata(:,4);
    L = th > CENTERANGLE-THETAWEDGE/2 & th < CENTERANGLE+THETAWEDGE/2 & ~Loog;
    contrast = sqrt(rawdata(L,1).^2+rawdata(L,2).^2);
    tf = rawdata(L,3);
    plot(tf,1./contrast,'o','MarkerEdgeColor','yellow','MarkerFaceColor','yellow','MarkerSize',10)
end
set(gca,'Xscale','log','Yscale','log');
plot([w(1) w(end)],[1 1]*photon_count_sensitivity,'w-','LineWidth',2)
set(gca,'Ylim',[.1 1000]);
if SMALLAXES
    set(gca,'Xlim',[1 30],'Ylim',[1 1000]);
end

xlabel('Temporal frequency (Hz)');
ylabel('Sensitivity');
%%
% Section 2
% Just like section 1 above but comparing behavioral sensitivity to cone
% ideal observer's
RFX = 5; % DVA
RFY = 0; % DVA
MMPERDEG = 0.223; % mm/deg (Perry and Cowey 1985)
N_EYES = 2;
SMALLAXES = 0;
% Need to load in a stro structure that we will modify
filename = 'A060817004.nex'; % Just to get the basic format of an IsoSamp file
stro = nex2stro(findfile(fullfile(nexfilepath,'Greg','Apollo',filename)));
examplestrotrial = stro.trial(1,:);
stro.trial = [];
stro.ras = {};
% ----------------------

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
params.eyeNumber = N_EYES; 
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

Lccs = logspace(-3,log10(.3),10);
Mccs = logspace(-3,log10(.3),10);
% Putting new values in stro structure to be passed to DTcones
stro.sum.exptParams.rf_x = RFX*10; % Charlie's cone model wants positions in DVA*10
stro.sum.exptParams.rf_y = RFY*10; % Charlie's cone model wants positions in DVA*10

w = logspace(log10(1),log10(60),20); % temporal frequencies
cone_threshold = nan*ones(length(w),1);
for tf_counter = 1:length(w)
    for i = 1:length(Lccs)
        stro.trial(i,:) = examplestrotrial;
        stro.trial(i,strcmp(stro.sum.trialFields(1,:),'stim_l')) = Lccs(i);
        stro.trial(i,strcmp(stro.sum.trialFields(1,:),'stim_m')) = Mccs(i);
        stro.trial(i,strcmp(stro.sum.trialFields(1,:),'tf'))= w(tf_counter);
    end
    params.stro = stro;
    [gab, cones, mon, idlob, params] = DTcones_gh(params,0);
    
    v = idlob.analyticVar{1}([1 2]); % variance
    m = []; % mean
    for i = 1:length(idlob.analyticMean); m(i,:) = idlob.analyticMean{i}([1 2]); end
    cone_dprimes = sqrt(m.^2*(1./v'));
    % I want to find a contrast that will make d-prime = 1.27
    if (min(cone_dprimes) > 1.27 | max(cone_dprimes) < 1.27)
        error('Threshold not within bounds. Expand contrast bounds.');
    end
    cone_threshold(tf_counter) = interp1(cone_dprimes, sqrt(Lccs.^2+Mccs.^2)', 1.27,'spline');
end
figprefs(1); axes; hold on; set(gcf,'Color',[0 0 0],'InvertHardCopy','off');
set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]);
plot(w,1./cone_threshold,'c-','LineWidth',3);
set(gca,'Xscale','log','Yscale','log','Ylim',[10^-1 10^3]);
if SMALLAXES
    set(gca,'Xlim',[1 30],'Ylim',[1 1000]);
end
%%
% Section 3
% Getting cone model d' for stimuli at the monkey's psychophysical
% threshold.
RFX = 5; % DVA
RFY = 0; % DVA
MMPERDEG = 0.223; % mm/deg (Perry and Cowey 1985)
N_EYES = 2;
w = logspace(log10(1),log10(60),40); % temporal frequencies
theta = pi/4; % Color direction in LM plane

% Loading behavioral model data
isosamppath = which('IsoSampOnline');
isosamppath(find(isosamppath==filesep,1,'last'):end) = [];
load ([isosamppath,filesep,'private',filesep','data',filesep,'LMTF.mat']);
% "A" is the structure of LMTF model parameters for Apollo
behavioralmodel = LMTF_global_to_local_model(A.legacy.mode3params, RFX, RFY, 3);
behavioral_tcsf_lum = @(omega)(behavioralmodel(1)*abs(((1i*2*pi*10^behavioralmodel(5).*omega+1).^-behavioralmodel(3))-behavioralmodel(2)*((1i*2*pi*10^(behavioralmodel(5)+behavioralmodel(6)).*omega+1).^-(behavioralmodel(3)+behavioralmodel(4)))));
behavioral_tcsf_rg = @(omega)(behavioralmodel(7)*abs(((1i*2*pi*10^behavioralmodel(11).*omega+1).^-behavioralmodel(9))-behavioralmodel(8)*((1i*2*pi*10^(behavioralmodel(11)+behavioralmodel(12)).*omega+1).^-(behavioralmodel(9)+behavioralmodel(10)))));
behavioral_tcsf = @(omega, theta) sqrt((behavioral_tcsf_lum(omega).*(cos(theta).*cos(behavioralmodel(13))+(sin(theta).*sin(behavioralmodel(13))))).^2+...
     (behavioral_tcsf_rg(omega)*(cos(theta)./sqrt(2)-(sin(theta)./sqrt(2)))).^2); 
 
% Need to load in a stro structure that we will modify
filename = 'A060817004.nex'; % Just to get the basic format of an IsoSamp file
stro = nex2stro(findfile(fullfile(nexfilepath,'Greg','Apollo',filename)));
examplestrotrial = stro.trial(1,:);
stro.trial = [];
stro.ras = {};
% ----------------------

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
params.eyeNumber = N_EYES; 
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

% Putting new values in stro structure to be passed to DTcones
stro.sum.exptParams.rf_x = RFX*10; % Charlie's cone model wants positions in DVA*10
stro.sum.exptParams.rf_y = RFY*10; % Charlie's cone model wants positions in DVA*10
behavior_threshold = 1./behavioral_tcsf(w, theta);

cone_dprimes = zeros(size(w));
for i = 1:length(w)
    stro.trial(i,:) = examplestrotrial;
    stro.trial(i,strcmp(stro.sum.trialFields(1,:),'stim_l')) = cos(theta)*behavior_threshold(i);
    stro.trial(i,strcmp(stro.sum.trialFields(1,:),'stim_m')) = sin(theta)*behavior_threshold(i);
    stro.trial(i,strcmp(stro.sum.trialFields(1,:),'tf'))= w(i);
end
params.stro = stro;
[gab, cones, mon, idlob, params] = DTcones_gh(params,0);
for i = 1:length(idlob.analyticVar)
    v = idlob.analyticVar{i}([1 2]); % variance
    m = idlob.analyticMean{i}([1 2]);
    cone_dprimes(i) = sqrt(m.^2*(1./v'));
end

figure; axes; hold on; set(gcf,'Color',[0 0 0],'InvertHardCopy','off');
set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]);
if theta == 3*pi/4
    Lw = w<=20;
else
    Lw = true(size(w));
end
plot(w(Lw),cone_dprimes(Lw),'c-','LineWidth',2);
set(gca,'Xscale','log','Yscale','linear','Xlim',[w(find(Lw,1,'first')) w(find(Lw,1,'last'))],'Ylim',[0 max(cone_dprimes)]);
plot([w(find(Lw,1,'first')) w(find(Lw,1,'last'))],[1.27 1.27],'y-','LineWidth',2);
ylabel('D''');
xlabel('Temporal frequency (Hz)');

%%
% Section 3.1
% Normal distributions at various d'
color = [1 1 0];
dprime = 5;
x = linspace(-4,4+dprime,100);
y1 = normpdf(x,0,1);
y2 = normpdf(x,dprime,1);
figure; axes; hold on;
plot(x,y1,'k-','LineWidth',2);
h1 = patch(x',y1',color/1.2,'FaceAlpha',.5,'LineStyle','none');
plot(x',y2','k-','LineWidth',2);
h2 = patch(x,y2,color,'FaceAlpha',.5,'LineStyle','none');
set(gca,'Xlim',[x(1) x(end)]);
set(gca,'Xtick',[]);
set(gcf,'Renderer','Painters');

%%
% Section 4
% Example rasters and d' measurements for a single example neuron
filename = 'A061517005.nex'; % example magnocell
%exampleparvocell = {'A062717002.nex', 'A062717004.nex'}; % example parvocell
%filename = exampleparvocell;

STIMTYPE = 'LUM';

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

% A plot of a single set of rasters?
%LTF = uniquestim(:,3) > 1 & uniquestim(:,3) < 1.6;
%LTF = uniquestim(:,3) < 1.2;
Lblank = all(uniquestim == 0,2);
if strcmp(STIMTYPE, 'LUM')
    Lstimtype = (sign(uniquestim(:,1)) == sign(uniquestim(:,2))) & ~Lblank;
elseif strcmp(STIMTYPE, 'RG')
    Lstimtype = (sign(uniquestim(:,1)) ~= sign(uniquestim(:,2))) & ~Lblank;
else
    error('unknown stimtype');
end
%Lwhichstimtypes = Lstimtype & LTF;
Lwhichstimtypes = Lstimtype;

% Rasters for all conditions except the blank
offset = [-.1 .1];  % pre and post time wrt stimon/stimoff
figure; axes; hold on; set(gcf,'Color',[0 0 0],'InvertHardCopy','off');
set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]); counter = 0;
for j = find(Lwhichstimtypes)'
    L = Lcc == uniquestim(j,1) & Mcc == uniquestim(j,2) & TF == uniquestim(j,3);
    for i = find(L)'
        tmpspikes = spikes{i}-stimon_t(i);
        tmpspikes(tmpspikes < offset(1) | tmpspikes > dur+offset(2)) = [];
        nspikestot = length(tmpspikes);
        plot([tmpspikes tmpspikes]',[zeros(nspikestot,1) 1*ones(nspikestot,1)]'+counter,'y-','linewidth',1);
        counter = counter + 1;
    end
    plot([offset(1) dur+offset(2)],counter*[1 1],'w:');
end
set(gca,'Xlim',[0+offset(1) dur+offset(2)],'Ylim',[0 counter],'Ytick',[],'Xtick',[0 .2 .4 .6 .8],'Box','on');
set(gcf,'Renderer','painters')

% PSTHs
nbins = 30;
figure; axes; hold on; set(gcf,'Color',[0 0 0],'InvertHardCopy','off');
set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]);
tbinedges = linspace(offset(1),dur+offset(2),nbins+1);
tbincenters = tbinedges(1:end-1)+(tbinedges(2)-tbinedges(1))/2;
counter = 0;
for j = find(Lwhichstimtypes)'
    L = Lcc == uniquestim(j,1) & Mcc == uniquestim(j,2) & TF == uniquestim(j,3);
    PSTH = zeros(size(tbinedges)-[0 1]);
    for i = find(L)'
        tmpspikes = spikes{i}-stimon_t(i);
        tmpspikes(tmpspikes < offset(1) | tmpspikes > dur+offset(2)) = [];
        PSTH = PSTH+histcounts(tmpspikes,tbinedges);
    end
    PSTH = PSTH./sum(L);
    plot(tbincenters,counter+PSTH,'y','Linewidth',2);
    counter = counter+1;
end
set(gca,'Xlim',[0+offset(1) dur+offset(2)],'Ylim',[0 counter+1],'Ytick',[],'Xtick',[0 .2 .4 .6 .8],'Box','on');
set(gcf,'Renderer','painters')

% Raster and PSTH for the blank
figure; subplot(2,1,1); hold on; set(gcf,'Color',[0 0 0],'InvertHardCopy','off');
set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]);
counter = 0;
L = Lcc == 0 & Mcc == 0 & TF == 0;
for i = find(L)'
    PSTH = zeros(size(tbinedges)-[0 1]);
    tmpspikes = spikes{i}-stimon_t(i);
    tmpspikes(tmpspikes < offset(1) | tmpspikes > dur+offset(2)) = [];
    nspikestot = length(tmpspikes);
    plot([tmpspikes tmpspikes]',[zeros(nspikestot,1) 1*ones(nspikestot,1)]'+counter,'y-','linewidth',1);
    PSTH = PSTH+histcounts(tmpspikes,tbinedges);
    counter = counter + 1;
end
PSTH = PSTH./sum(L);
plot([offset(1) dur+offset(2)],counter*[1 1],'b:');
set(gca,'Xlim',[0+offset(1) dur+offset(2)],'Ylim',[0 counter],'Ytick',[],'Xtick',[],'Box','on');
set(gcf,'Renderer','painters')
subplot(2,1,2); hold on; set(gcf,'Color',[0 0 0],'InvertHardCopy','off');
set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]);
plot(tbincenters,PSTH,'y','Linewidth',2);
set(gca,'Ylim',[0 1])


closest = @(x,y) find(((x-y).^2) == min((x-y).^2)); % Support function for temporal envelope calculation
bins = linspace(offset(1),offset(2)+dur,1000);
pivotbins = [closest(0,bins) closest(dur/4,bins) closest(3*dur/4,bins) closest(dur,bins)];
temporalenvelope = zeros(size(bins));
temporalenvelope(pivotbins(1):pivotbins(2)) = linspace(0,1,pivotbins(2)-pivotbins(1)+1);
temporalenvelope(pivotbins(2):pivotbins(3)) = ones(1,pivotbins(3)-pivotbins(2)+1);
temporalenvelope(pivotbins(3):pivotbins(4)) = linspace(1,0,pivotbins(4)-pivotbins(3)+1);

% plotting the stimulus
figure; axes; hold on;
counter = 0;
for j = find(Lwhichstimtypes)'
    L = Lcc == uniquestim(j,1) & Mcc == uniquestim(j,2) & TF == uniquestim(j,3);
    counter = counter+sum(L);
    tf = uniquestim(j,3);
    y = cos(bins*2*pi*tf);
    plot(bins,(sum(L)/4)*(y.*temporalenvelope)+counter)
end

% Dprime plot
figure; axes; hold on; set(gcf,'Color',[0 0 0],'InvertHardCopy','off');
set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]);
plot(uniquestim(Lwhichstimtypes,3),dprime(Lwhichstimtypes),'ro-','LineWidth',2,'MarkerFaceColor',[1 .2 .2],'MarkerEdgeColor',[1 .2 .2], 'Color',[1 .2 .2]);
set(gca,'Xscale','log');
xlabel('Temporal frequency (Hz)');
ylabel('d''');
set(gcf,'Renderer','painters')

% histograms for calucating d-prime
lumidxs = find(sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & uniquestim(:,3) > 0)';
for i = lumidxs
    bins = [floor(min([noise{i}]))-1:.5:ceil(max([signal{:}]))+1]
    y1 = hist(signal{i},bins);
    y2 = hist(noise{i},bins);
    figprefs(1);
    bar(bins,[y1;y2]',1.5);
    colormap([1 1 0; .6 .6 .6])
end
%%
% Section 5
% Plotting a hexagonal grid of RF locations over a representation of the
% stimulus. 
%filename = 'A061517005.nex'; CELLTYPE = 'M'; % example magnocell
filename = 'A062717002.nex'; CELLTYPE = 'P'; % example parvocell
RFTRUNCATIONINSD = 3;

if CELLTYPE == 'M'
    ecc_to_diam_deg = @(rf_r_deg) 10.^(-1.2459+0.0345*rf_r_deg); % temporal retina equivalent
elseif CELLTYPE == 'P'
    a = 0.9729; % Table 1
    r2 = 1.084; % Table 1
    re = 7.633; % Table 1
    dc_0 = 14804.6; % Cone density of fovea
    rm = 41.03; % See Equation 7
    ecc_to_diam_deg = @(rf_r_deg)(sqrt(2./(sqrt(3).*... % Equation 9. Distance between adjacent midget RF centers.
        (2*dc_0.*(1+rf_r_deg./rm).^-1.*(a*(1+(rf_r_deg./r2)).^-2+(1-a)*exp(-rf_r_deg./re)))...
        /2))... % Dividing y by 2 to estimate RF size from only ON or OFF mosaics (halved density).
        *.80); % Dacey and Petersen: monkey midget DFs are .77:.81 x smaller than human midget DFs
else
    error('Unknown CELLTYPE');
end
% Support function for calculation of mus and S2
bpdf_vec=@(x,y,mu_x,mu_y,sigma)(exp(-((x-mu_x).^2)./(2*sigma^2)-((y-mu_y).^2)/(2*sigma^2))./(2*pi*sigma^2));
truncated_2D_gaussian = @(x,y,mu_x,mu_y,sigma,trunc) (sqrt((x-mu_x).^2+(y-mu_y).^2)<trunc).*bpdf_vec(x,y,mu_x,mu_y,sigma) + (sqrt((x-mu_x).^2+(y-mu_y).^2)>=trunc).*0;

stro = nex2stro(findfile(filename,[nexfilepath,filesep,'Greg',filesep,'Apollo']));
sigma_gabor = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sigma'))); % DVA
sigmas_n = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sigmas_n')));
rf_ecc = sqrt((stro.sum.exptParams.rf_x/10).^2+(stro.sum.exptParams.rf_y/10).^2);
RF_diam_deg = ecc_to_diam_deg(rf_ecc); % 2 SDs of a Gaussian RF
RFdistance = RF_diam_deg;
RF_STD = RF_diam_deg/2; % 1 standard deviation of Gaussian RF
[x_deg,y_deg] = meshgrid(linspace(-sigma_gabor*2*sigmas_n,sigma_gabor*2*sigmas_n,100));
stim = normpdf(x_deg,0,sigma_gabor).*normpdf(y_deg,0,sigma_gabor);

Rad3Over2 = sqrt(3)/2;
x_centers = [x_deg(1,1):RFdistance:x_deg(1,end)];
closest_to_zero = find(abs(x_centers) == min(abs(x_centers)),1);
x_centers = x_centers-x_centers(closest_to_zero);
x_centers_mat = repmat(x_centers,size(x_centers,2),1);
y_centers = x_centers;
x_centers_mat = x_centers_mat*Rad3Over2;
y_centers_mat = repmat(y_centers',1,size(y_centers,2));
y_centers_mat(:,rem(find(y_centers == 0),2)+1:2:end) = y_centers_mat(:,rem(find(y_centers == 0),2)+1:2:end)+.5*RFdistance;
mus = zeros(numel(x_centers_mat),1);
interRFdistances = zeros(numel(x_centers_mat),numel(x_centers_mat));
for j = 1:numel(x_centers_mat)
    for k = 1:numel(x_centers_mat) % Tabulating distances between RF centers which we'll need for S2
        interRFdistances(j,k) = sqrt((x_centers_mat(j)-x_centers_mat(k)).^2+(y_centers_mat(j)-y_centers_mat(k)).^2);
    end
end

RFtoStimdistances = sqrt(x_centers_mat.^2+y_centers_mat.^2);
uniqueRFtoStimdistances = unique(RFtoStimdistances);

for j = 1:length(uniqueRFtoStimdistances)
    if uniqueRFtoStimdistances(j) > RF_STD*RFTRUNCATIONINSD+sigma_gabor*sigmas_n
        mus(j) = 0;
    else
        overlap_point = @(x,y) truncated_2D_gaussian(x,y,0,0,RF_STD,RF_STD*RFTRUNCATIONINSD).*truncated_2D_gaussian(x,y,0,uniqueRFtoStimdistances(j),sigma_gabor,sigma_gabor*sigmas_n); % Points sigmas_n away from Gabor center are not plotted
        mu=integral2(overlap_point,-RFTRUNCATIONINSD*RF_STD, RFTRUNCATIONINSD*RF_STD, max(-RFTRUNCATIONINSD*RF_STD, uniqueRFtoStimdistances(j)-sigma_gabor*sigmas_n),RFTRUNCATIONINSD*RF_STD,'method','iterated');
        mus(RFtoStimdistances(:) == uniqueRFtoStimdistances(j)) = mu;
    end
    for k = 1:numel(x_centers_mat) % Tabulating distances between RF centers which we'll need for S2
        interRFdistances(j,k) = sqrt((x_centers_mat(j)-x_centers_mat(k)).^2+(y_centers_mat(j)-y_centers_mat(k)).^2);
    end
end

S2 = nan(numel(x_centers_mat));
for dist = unique(interRFdistances)' % This calculation assumes infinitely large RFs?
    overlap_point = @(x,y) bpdf_vec(x,y,0,0,RF_STD).*bpdf_vec(x,y,0,dist,RF_STD);
    S2(interRFdistances==dist)=integral2(overlap_point,RFTRUNCATIONINSD*RF_STD,-RFTRUNCATIONINSD*RF_STD,RFTRUNCATIONINSD*RF_STD,-RFTRUNCATIONINSD*RF_STD);
end

mus = mus./max(mus);
S2 = S2./max(S2(:)); % If two RFs are identical, cov = 1 & cov is proportional to overlap
weights = S2\mus; % these are the ideal weights. No need to normalize to max(weights)
single_mosaic_mean = mus'*weights;
single_mosaic_var=weights'*S2*weights;

population_scalefactor=(mus'*weights)./sqrt(single_mosaic_var); % Single mosaic
population_scalefactor = (2*single_mosaic_mean)/sqrt(2*single_mosaic_var+2*0.05*single_mosaic_var);
population_scalefactor = (2*single_mosaic_mean)/sqrt(2*single_mosaic_var+2*0.05*single_mosaic_var)*sqrt(2); % extra "2" is for # of eyes

figure; axes; hold on;
tmp = RF_STD*[cos(linspace(0,2*pi,200))', sin(linspace(0,2*pi,200))'];
for j = 1:numel(x_centers_mat)
    plot(x_centers_mat(j)+tmp(:,1),y_centers_mat(j)+tmp(:,2),'k-');
end
axis image
set(gca,'Xlim',[x_deg(1), x_deg(end)],'Ylim',[y_deg(1), y_deg(end)]);
set(gcf,'Renderer','painters');

% image of Gabor stimulus
figure; axes; hold on;
imagesc(stim);
colormap(gray);
axis image;
population_scalefactor

%%
% Section 6)
% Distribution "population d's" across an actual population of recorded
% neuron along with cone d's. See sections 4 and 4.1 of IsoSampPop.m

CELLTYPE = 'P';
STIMTYPE = 'RG'; % Which to look at, LUM or RG? Both are computed but only one is plotted.

MMPERDEG = 0.223; % mm/deg (Perry and Cowey 1985)
DEGPERMM = 1/MMPERDEG; % deg/mm
N_EYES = 2; % How many eyes the ideal observer has
ONOFFCORRELATION = .05; % Correlation between ON and OFF cells
RFTRUNCATIONINSD = 3;
H_ECCMULTIPLIER = 0.8; % To account for half of the RFs being in the nasal retina and half in the temporal retina

% Getting IsoSamp files
[filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',{CELLTYPE});

if CELLTYPE == 'M'
    ecc_to_diam_deg = @(rf_r_deg) 10.^(-1.2459+0.0345*rf_r_deg); % D&L 1984 temporal retina equivalent
elseif CELLTYPE == 'P' % Watson  temporal retina equivalent
    HUMAN2MONKPSCALEFACTOR = .80; % From Dacey and Petersen. reasonable range: [.77 .81];
    a = 0.9729; % Table 1
    r2 = 1.084; % Table 1
    re = 7.633; % Table 1
    dc_0 = 14804.6; % Cone density of fovea
    rm = 41.03; % See Equation 7
    ecc_to_diam_deg = @(x)(sqrt(2./(sqrt(3).*... % Equation 9. Distance between adjacent midget RF centers.
        (2*dc_0.*(1+x./rm).^-1.*(a*(1+(x./r2)).^-2+(1-a)*exp(-x./re)))...
        ./2))... % Dividing y by 2 to estimate RF size from only ON or OFF mosaics (halving the density).
        *HUMAN2MONKPSCALEFACTOR); % Monkey midget RFs are slightly smaller than human midget RFs
else
    error('Unknown CELLTYPE');
end

% Support function for calculation of mus and S2
bpdf_vec=@(x,y,mu_x,mu_y,sigma)(exp(-((x-mu_x).^2)./(2*sigma^2)-((y-mu_y).^2)/(2*sigma^2))./(2*pi*sigma^2));

% Setting up for cone model analysis
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
params.eyeNumber = N_EYES; 
params.coneSampRate = 2400;
params.flatPowerSpect = false;
params.enableScones = false;
params.sacamp_deg = 0;
params.sacdur_s = 0;

% Getting neuron population d-primes
data = [];
for i = 1:length(filenames)
    stro = {};
    for j = 1:length(filenames{i})
        stro{j} = nex2stro(char(findfile(filenames{i}(j), fullfile(nexfilepath,'Greg'))));
    end
    stro = strocat(stro);
    rfx = stro.sum.exptParams.rf_x/10;
    rfy = stro.sum.exptParams.rf_y/10;
    rf_r_deg = sqrt((rfx/H_ECCMULTIPLIER)^2+rfy^2);
    RF_diam_deg = ecc_to_diam_deg(rf_r_deg); % SD of Gaussian fit
    RFdistance = RF_diam_deg; % RF centers 2 SDs apart, from EJ. Convention is diameter is 2 SDs
    RF_STD = RF_diam_deg/2; % 1 standard deviation of Gaussian RF

    sigma_gabor = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sigma'))); % DVA
    sigmas_n = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sigmas_n')));
    [x_deg,y_deg] = meshgrid(linspace(-sigma_gabor*2*sigmas_n,sigma_gabor*2*sigmas_n,40));
    Rad3Over2 = sqrt(3)/2;
    x_centers = [x_deg(1,1):RFdistance:x_deg(1,end)];
    closest_to_zero = find(abs(x_centers) == min(abs(x_centers)),1);
    x_centers = x_centers-x_centers(closest_to_zero);
    x_centers_mat = repmat(x_centers,size(x_centers,2),1);
    y_centers = x_centers;
    x_centers_mat = x_centers_mat*Rad3Over2;
    y_centers_mat = repmat(y_centers',1,size(y_centers,2));
    y_centers_mat(:,rem(find(y_centers == 0),2)+1:2:end) = y_centers_mat(:,rem(find(y_centers == 0),2)+1:2:end)+.5*RFdistance;
    mus = zeros(numel(x_centers_mat),1);
    S2 = eye(numel(x_centers_mat));
    interRFdistances = zeros(numel(x_centers_mat),numel(x_centers_mat));
    for j = 1:numel(x_centers_mat)
        overlap_point = @(x,y) bpdf_vec(x,y,0,0,RF_STD).*bpdf_vec(x,y,x_centers_mat(j),y_centers_mat(j),sigma_gabor);
        mus(j)=integral2(overlap_point,RFTRUNCATIONINSD*RF_STD,-RFTRUNCATIONINSD*RF_STD,RFTRUNCATIONINSD*RF_STD,-RFTRUNCATIONINSD*RF_STD); % RF is truncated at "n" SDs
        for k = 1:numel(x_centers_mat) % Tabulating distances between RF centers which we'll need for S2
            interRFdistances(j,k) = sqrt((x_centers_mat(j)-x_centers_mat(k)).^2+(y_centers_mat(j)-y_centers_mat(k)).^2);
        end
    end
    S2 = nan(numel(x_centers_mat));
    for dist = unique(interRFdistances)' % This calculation assumes infinitely large RFs?
        overlap_point = @(x,y) bpdf_vec(x,y,0,0,RF_STD).*bpdf_vec(x,y,0,dist,RF_STD);
        S2(interRFdistances==dist)=integral2(overlap_point,RFTRUNCATIONINSD*RF_STD,-RFTRUNCATIONINSD*RF_STD,RFTRUNCATIONINSD*RF_STD,-RFTRUNCATIONINSD*RF_STD);    
    end
    
    mus = mus./max(mus);
    S2 = S2./max(S2(:)); % If two RFs are identical, cov = 1 & cov is proportional to overlap
    weights = S2\mus; % these are the ideal weights. No need to normalize to max(weights)
    single_mosaic_mean = mus'*weights;
    single_mosaic_var=weights'*S2*weights;
    population_scalefactor = (2*single_mosaic_mean)/sqrt(2*single_mosaic_var+2*abs(ONOFFCORRELATION)*single_mosaic_var)*sqrt(N_EYES)
    
    [uniquestim, dprime] = IsoSampGetDPrime(stro,1,spikecds(i));
    data{i}.rfxy = [stro.sum.exptParams.rf_x/10, stro.sum.exptParams.rf_y/10];
    data{i}.uniquestim = uniquestim;
    data{i}.neurondprime = dprime.*population_scalefactor;
    
    % Now getting the cone dprimes
    params.stro = stro;
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
    for j = 1:size(idlob.analyticMean,1) % looping over color direction
        for k = 1:size(idlob.analyticMean(j,:),2) % looping over contrast/TF
            if ~isempty(idlob.analyticMean{j,k})
                tmp_lm_mu = idlob.analyticMean{j,k}([1 2]);
                tmp_lm_var = idlob.analyticVar{j,k}([1 2]);
                tf = gab.driftRates{j}(k);
                tmpdata = [tmpdata; gab.colorDirs(j,[1 2]).*gab.contrasts{j}(k) tf tmp_lm_mu tmp_lm_var];
            end
        end
    end
    
    % Using uniquestim to order the rows of tmpdata and calculating cone
    % dprimes
    data{i}.conedprimes = [];  
    for j = 1:size(data{i}.uniquestim,1)
        L = all(abs(data{i}.uniquestim(j,:)-tmpdata(:,[1 2 3]))<1e-10,2);
        if sum(L) ~= 1
            if all(data{i}.uniquestim(j,:) == 0)
                data{i}.conedprime(j) = nan;
            else
                error('sorting problem');
            end
        else
            v = tmpdata(L,[6 7]); % variance
            m = tmpdata(L,[4 5]); % mean
            data{i}.conedprime(j) = sqrt(m.^2*(1./v'));
        end
    end
end

% Getting rid of neurons with out-of-range RFs
deletelist = [];
for i = 1:length(data)
    in_ecc_range = sqrt(data{i}.rfxy*data{i}.rfxy') >= 2 & abs(data{i}.rfxy(1)) <= 12 & abs(data{i}.rfxy(2))<= 8;
    if ~in_ecc_range
        disp('Chucking a cell because RF is outside of the acceptable range');
        data{i}.rfxy
        deletelist = [deletelist i];
    end
end
data(deletelist) = [];
deletelist = [];

% Preparing for plotting. Computing conedprimes and neurondprimes.
%if strcmp(STIMTYPE,'LUM')
    tfbinedges = logspace(0,log10(60),10);
%else
%    tfbinedges = logspace(0,log10(30),8);
%end
tfbincenters = sqrt(tfbinedges(1:end-1).*tfbinedges(2:end)); % geomean
neurondprimes = cell(size(tfbincenters));
conedprimes = cell(size(tfbincenters));
ns = zeros(1,length(tfbincenters));
for i = 1:length(data)
    in_ecc_range = sqrt(data{i}.rfxy*data{i}.rfxy') >= 2 & abs(data{i}.rfxy(1)) <= 12 & abs(data{i}.rfxy(2))<= 8;
    if ~in_ecc_range
        continue
    end
    if strcmp(STIMTYPE,'LUM')
        Lcolordir = sign(data{i}.uniquestim(:,1)) == sign(data{i}.uniquestim(:,2)) & data{i}.uniquestim(:,3) ~= 0;
    elseif strcmp(STIMTYPE,'RG')
        Lcolordir = sign(data{i}.uniquestim(:,1)) ~= sign(data{i}.uniquestim(:,2));
    else
        error('Unknown stimtype');
    end
    for j = 1:length(tfbincenters)
        Ltf = data{i}.uniquestim(:,3) > tfbinedges(j) & data{i}.uniquestim(:,3) <= tfbinedges(j+1);
        if sum(Ltf&Lcolordir) > 0
            conedprimes{j} = [conedprimes{j};mean(data{i}.conedprime(Ltf&Lcolordir))];
            neurondprimes{j} = [neurondprimes{j};mean(data{i}.neurondprime(Ltf&Lcolordir))];
            ns(1,j) = ns(1,j) + 1; % each cell counts as an independent entity
        end
    end
end

m = []; sem = []; % neurons; cones
for i = 1:length(neurondprimes)
    m(1,i) = mean(neurondprimes{i});
    sem(1,i) = sqrt(var(neurondprimes{i})./ns(i));
    m(2,i) = mean(conedprimes{i});
    sem(2,i) = sqrt(var(conedprimes{i})./ns(i));
end

% Plotting
figure; axes; hold on; set(gcf,'Color',[0 0 0],'InvertHardCopy','off');
set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]);L = ns(1,:) > 1;
colors = {[1 .2 .2],[0 .75 .75]};
for i = 1:2
    h = patch([tfbincenters(L), fliplr(tfbincenters(L))],[m(i,L)+sem(i,L), fliplr(m(i,L)-sem(i,L))],colors{i});
    set(h,'Facealpha',.5);
    h = plot(tfbincenters(L),m(i,L),'k-','Linewidth',2,'Color',colors{i});
end
set(gca,'Xscale','log');
plot([tfbincenters(1) tfbincenters(end)],[1.27 1.27],'y-','LineWidth',2)
xlabel('Temporal frequency (Hz)');
ylabel('d''');
title([CELLTYPE,' n=(',num2str(size(filenames,1)),')']);
set(gca,'Ylim',[-5 25]);
set(gcf,'Renderer','painters');
if strcmp(STIMTYPE,'RG')
    set(gca,'Xlim',[1 30]);
end
title(['n=',num2str(length(data))],'Color','white')
%%
% Section 7)
% RF positions of magnocellular and parvocellular neurons.

[parvocellnames,~] = fnamesFromTxt('IsoSamp_LGN','cellClass',{'P'});
[magnocellnames,~] = fnamesFromTxt('IsoSamp_LGN','cellClass',{'M'});

exampleparvocell = 'A062717002.nex'; % example parvocell. See section 4.

parvoRFs = [];
magnoRFs = [];
examplecellRF = [];
for i = 1:2
    if i == 1
        filenames = parvocellnames;
    else
        filenames = magnocellnames;
    end
    RFs = [];
    for j = 1:length(filenames)
        if length(filenames{j}) > 1
            filenames{j} = filenames{j}(1);
        end
        stro = nex2stro(char(findfile(filenames{j}, fullfile(nexfilepath,'Greg','Apollo'))));
        RFs = [RFs; stro.sum.exptParams.rf_x/10 stro.sum.exptParams.rf_y/10];
        if strcmp(filenames{j},exampleparvocell)
            examplecellRF = RFs(end,:);
            disp('got here');
        end
    end
    if i == 1
        parvoRFs = RFs;
    else
        magnoRFs = RFs;
    end
end

figure; axes; hold on; set(gcf,'Color',[0 0 0],'InvertHardCopy','off');
set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]);
hp = plot(abs(magnoRFs(:,1)),magnoRFs(:,2),'ko','MarkerFaceColor','white','MarkerSize',6);
hm = plot(abs(parvoRFs(:,1)),parvoRFs(:,2),'rs','MarkerFaceColor','red','MarkerSize',4);
set(gca,'Ylim',max(abs([magnoRFs(:,2); parvoRFs(:,2)]))*[-1 1]);
if ~isempty(examplecellRF)
    plot(abs(examplecellRF(1)),examplecellRF(2), 'y^');
end
axis equal;
axis square;
set(gcf,'Renderer','painters');

%% 
% Section 8
% Input efficiencies. Largely taken from IsoSampPop Section 3
CELLTYPE = 'M';
[filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',{CELLTYPE});

if CELLTYPE == 'M'
    ecc_to_diam_deg = @(rf_r_deg)(0.0080*rf_r_deg+0.05)/sqrt(2);
elseif CELLTYPE == 'P'
    ecc_to_diam_deg = @(rf_r_deg)(0.0032*rf_r_deg+0.02)/sqrt(2);
else
    error('Unknown CELLTYPE');
end

fundamentals = load('T_cones_smj10');
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
for a = 1:length(filenames)
    stro = {};
    for j = 1:length(filenames{a})
        stro{j} = nex2stro(char(findfile(filenames{a}(j), fullfile(nexfilepath,'Greg','Apollo'))));
    end
    stro = strocat(stro);
    [uniquestim, dprime] = IsoSampGetDPrime(stro, 1, spikecds(a));
    
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
    rf_r_deg = sqrt(.8*rfx^2+rfy^2);
    params.gab.sd = ecc_to_diam_deg(rf_r_deg);
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
    % columns of data: L, M, TF, mu_L, mu_M, sigma_L, sigma_M
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
    data{a}.rfxy = [stro.sum.exptParams.rf_x/10, stro.sum.exptParams.rf_y/10];
    data{a}.uniquestim = uniquestim;
    data{a}.neurondprime = dprime;
    data{a}.conedprime = cone_dprimes;
end

% Computing efficiencies.
tfbinedges = logspace(0,log10(60),8);
tfbincenters = sqrt(tfbinedges(1:end-1).*tfbinedges(2:end)); % geomean
efficiencies = cell(2,length(tfbincenters));
ns = zeros(2,length(tfbincenters));
    
for i = 1:length(data)
    in_ecc_range = sqrt(data{i}.rfxy*data{i}.rfxy') >= 2 & abs(data{i}.rfxy(1)) <= 12 & abs(data{i}.rfxy(2))<= 8;
    if ~in_ecc_range
        continue
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

colors = {[1 1 1],[1 .2 .2]};
figure; axes; hold on; set(gcf,'Color',[0 0 0],'InvertHardCopy','off');
set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]);
for i = 1:2
    L = ns(i,:) > 2; % minimum 'n'
    h = plot(tfbincenters(L), m(i,L),'w-o','LineWidth',2);
    set(h,'Color',colors{i},'MarkerFaceColor',colors{i});
    h = patch([tfbincenters(L), fliplr(tfbincenters(L))],[m(i,L)+sem(i,L), fliplr(m(i,L)-sem(i,L))],colors{i},'Facealpha',.5);
end
set(gca,'Xscale','log','Xlim',[1 60]);
ylabel('Input efficiency');
xlabel('Temporal frequency (Hz)');
plot([tfbincenters(1) tfbincenters(end)],[0 0],'w-');
plot([tfbincenters(1) tfbincenters(end)],[1 1],'w-');
set(gcf,'Renderer','painters');
%%
% Section 9
% Putative S-potential + LGN cell
stro_wn1 = nex2stro(findfile(fullfile(nexfilepath,'Greg','Apollo','A061517004')));
stro_isosamp = nex2stro(findfile(fullfile(nexfilepath,'Greg','Apollo','A061517005')));
stro_wn2 = nex2stro(findfile(fullfile(nexfilepath,'Greg','Apollo','A061517006')));

stro = stro_isosamp;
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
dur = mean(stimoff_t-stimon_t);
Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_l'));
Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_m'));
TF = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));

% Cross-correlation plot
spikes1 = stro.ras(:,strcmp(stro.sum.rasterCells, 'sig001a'));
spikes2 = stro.ras(:,strcmp(stro.sum.rasterCells, 'sig001b'));
bins = linspace(0,dur,dur*10000);
data = zeros(1,2*(length(bins)-1)-1);
for i = 1:size(stimon_t,1)
    tmp1 = histcounts(spikes1{i}-stimon_t(i) ,bins);
    tmp2 = histcounts(spikes2{i}-stimon_t(i) ,bins);
    data = data+xcorr(tmp1, tmp2,'coeff');
end
figure; axes; hold on; set(gcf,'Color',[0 0 0],'InvertHardCopy','off','Renderer','painters');
set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]);
binwidth = bins(2)-bins(1);
plot([-1*floor(size(data,2)/2):1:floor(size(data,2)/2)]*binwidth,data./size(stimon_t,1),'wo-','MarkerFaceColor','white','LineWidth',2);
set(gca,'Xtick',[-.005:.001:.005],'XTickLabel',[-5:1:5])
title('negative numbers mean spike 2 leads spike 1');
xlabel('lag (ms)');
ylabel('correlation');
set(gca,'Xlim',[-.005 .005]);

% STAs
stro_wn = strocat(stro_wn1, stro_wn2);
maxT = 5;
figure;
for i = 1:2
    nstixperside = stro_wn.sum.exptParams.nstixperside;
    out = getWhtnsStats(stro_wn,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, stro_wn.sum.rasterCells{i});
    STAs = out{1};
    normfactor = 2*max(abs(STAs(:)));
    for whichframe = 1:maxT
        STA = STAs(:,whichframe)/normfactor+.5;
        STA = reshape(STA,[nstixperside nstixperside 3]);
        subplot(2,maxT,(i-1)*maxT+whichframe);
        image(STA);
        set(gca,'XTick',[],'YTick',[]);
        axis image;
    end
end

% Rasters from IsoSamp
for i = 1:2
    spikes = stro.ras(:,i);
    uniquestim = sortrows(unique([Lcc Mcc TF],'rows'),3);
    Lwhichstimtypes = (sign(uniquestim(:,1)) == sign(uniquestim(:,2))) & ~all(uniquestim == 0,2);
    
    offset = [-.1 .1];  % pre and post time wrt stimon/stimoff
    figure; axes; hold on; set(gcf,'Color',[0 0 0],'InvertHardCopy','off','Renderer','painters');
    set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]);
    counter = 0;
    for j = find(Lwhichstimtypes)'
        L = Lcc == uniquestim(j,1) & Mcc == uniquestim(j,2) & TF == uniquestim(j,3);
        for k = find(L)'
            tmpspikes = spikes{k}-stimon_t(k);
            tmpspikes(tmpspikes < offset(1) | tmpspikes > dur+offset(2)) = [];
            nspikestot = length(tmpspikes);
            plot([tmpspikes tmpspikes]',[zeros(nspikestot,1) 1*ones(nspikestot,1)]'+counter,'y-','linewidth',1);
            counter = counter + 1;
        end
        plot([offset(1) dur+offset(2)],counter*[1 1],'b:');
    end
    set(gca,'Xlim',[0+offset(1) dur+offset(2)],'Ylim',[0 counter],'Ytick',[],'Xtick',[0 .2 .4 .6 .8],'Box','on');
    title(num2str(i));
end

% Overlaid d-primes
[uniquestim, dprime1, signal1, noise1] = IsoSampGetDPrime(stro,1);
% switching colums in ras
tmpstro = stro;
spikes1 = stro.ras(:,1);
spikes2 = stro.ras(:,2);
tmpstro.ras(:,1) = spikes2;
[uniquestim, dprime2, signal2, noise2] = IsoSampGetDPrime(tmpstro,1);

dprimemat = [];
TFs = [];
for j = find(Lwhichstimtypes)'
    dprimemat = [dprimemat, [dprime1(j); dprime2(j)]];
    TFs=[TFs, uniquestim(j,3)];
end

figure; axes; hold on; set(gcf,'Color',[0 0 0],'InvertHardCopy','off','Renderer','painters');
set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]);
h = plot(TFs, dprimemat','LineWidth',2);
set(gca,'Xscale','log');
ylabel('d''');
xlabel('Frequency (Hz)');

% Comparing signal and noise separately for the two cells
s_and_n = [];
for i = 1:length(signal1)
    s_and_n = [s_and_n; nanmean(signal1{i}) nanmean(noise1{i}) nanmean(signal2{i}) nanmean(noise2{i})];
end

figure; axes; hold on;
plot(s_and_n(:,1),'b-');
plot(s_and_n(:,2),'b--','LineWidth',2);
plot(s_and_n(:,3),'r-');
plot(s_and_n(:,4),'r--','LineWidth',2);
legend({'Unit 1 signal','Unit 1 noise','Unit 2 signal','Unit 2 noise'});

%%
% Section 10
% Looking that the distribution of the decision variable (which is the
% Mahalanobis distance from the theoretical mean under a Poisson null
% hypothesis).

stro = nex2stro(findfile(fullfile(nexfilepath,'Greg','Apollo','A061517005')));
figure; set(gcf,'Color',[0 0 0],'InvertHardCopy','off');
axes; hold on;
set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]);

[uniquestim, dprime, signal, noise] = IsoSampGetDPrime(stro, 1, spikecds(a));
zstats = [];
bins = linspace(-4,4,30);
counts = zeros(size(bins));
for i = 1:size(uniquestim,1)
    plot(uniquestim(i,3),noise{i},'.','color',[.25 .25 .25]);
    if sign(uniquestim(i,1)) == sign(uniquestim(i,2))
        plot(uniquestim(i,3),signal{i},'.','color',[.5 .5 0]);
    else
        plot(uniquestim(i,3)*unifrnd(.9,1.1,length(signal{i}),1),signal{i},'.','color',[.5 0 0]);
    end
    zstats = [zstats; nanmean(signal{i}) nanstd(signal{i}) nanmean(noise{i}) nanstd(noise{i})];
    counts = counts+hist(signal{i}-mean(signal{i}),bins);
end
set(gca,'XScale','log');
xlabel('TF (Hz)'); ylabel('decision variable');
set(gca,'Xlim',[.8 100]);

title(stro.sum.fileName(find(stro.sum.fileName == filesep,1,'last')+1:end))
Llum = sign(uniquestim(:,1))==sign(uniquestim(:,2)) & ~uniquestim(:,1)==0;
Lrg = sign(uniquestim(:,1))~=sign(uniquestim(:,2));
plot(uniquestim(Llum,3), zstats(Llum,1),'y:','linewidth',2); % lum signal mean
plot(uniquestim(Llum,3), zstats(Llum,2),'y-','linewidth',2); % lum signal std
plot(uniquestim(Llum,3), zstats(Llum,3),':','linewidth',2,'color',[.5 .5 .5]); % noise mean
plot(uniquestim(Llum,3), zstats(Llum,4),'-','linewidth',2,'color',[.5 .5 .5]); % noise std
plot(uniquestim(Lrg,3), zstats(Lrg,1),'r:','linewidth',2); % rg signal mean
plot(uniquestim(Lrg,3), zstats(Lrg,2),'r-','linewidth',2); % rg signal std
figure; set(gcf,'Color',[0 0 0],'InvertHardCopy','off');
axes; hold on; set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]);
bar(bins,counts,'w');
gaussbins = linspace(bins(1),bins(end),100);
plot(gaussbins,sum(counts)*(length(gaussbins)/length(counts))*normpdf(gaussbins,0,1)./sum(normpdf(gaussbins,0,1)),'m-','LineWidth',2)
xlabel('Mean-subtracted DV'); ylabel('Counts');
%%
% Section 11
% d' from simulated Benardete and Kaplan P-cells (midget ganglion cells and
% average parvocellular neurons

filename = 'A061517002.nex'; % To get baseline firing rate and contrasts. RF=(10.6,-0.5).
whichparameters = 'Figure 6'; % 'Figure 6' or 'Table'. 'Table' is parvo cell, figure 6 is example midget RGC
i = sqrt(-1); % In case i is already defined
ntrials = 5000;

if strcmp(whichparameters,'Figure 6')
    % Taking parameters from Figure 6 legend
    % centers
    center.A.mean = 184.20; % spikes/sec/unit contrast
    center.Hs.mean = 0.69;% s
    center.Ts.mean = 18.61/1000; % s
    center.Nl.mean = 38;
    center.Tl.mean = 1.23/1000;  % s
    center.D.mean = 4.0/1000; % s
    
    % surrounds
    surround.A.mean = 125.33; % spikes/sec/unit contrast
    surround.Hs.mean = .56; % s
    surround.Ts.mean = 33.28/1000; % s
    surround.Nl.mean = 124;
    surround.Tl.mean = .42/1000;  % s
    surround.D.mean = 4.0/1000; % s
elseif strcmp(whichparameters,'Table')
    %
    % ---------- From the table (population) ----------
    % centers (from Table 5)
    center.A.mean = 34.93; % spikes/sec/unit contrast
    center.Hs.mean = .72;
    center.Ts.mean = 46.19/1000; % s
    center.Nl.mean = 27.13;
    center.Tl.mean = 52.48/center.Nl.mean/1000;  % parameterized this way because it's time to peak of STA
    center.D.mean = 3.4/1000; % s (Based on Table 5)
    
    % surrounds (from Table 6)
    surround.A.mean = 37.41; % spikes/sec/unit contrast
    surround.A.mean = 17.41; % actually median

    surround.Hs.mean = .26;
    surround.Ts.mean = 57.31/1000; % s
    surround.Nl.mean = 61.88;
    surround.Tl.mean = 59.97/surround.Nl.mean/1000;  % parameterized this way because it's time to peak of STA
    surround.D.mean = 3.5/1000; % s (Based on Table 3 and Table 6)
else
    error('Unknown parameter set');
end

TF = [0:1:300];
n = length(TF);
omega = TF*2*pi; % rad/s
freqresponse_c = center.A.mean*exp(-i.*omega.*center.D.mean).*...
    (1-(center.Hs.mean./(1+i.*omega.*center.Ts.mean))).*...
    (1./(1+i.*omega.*center.Tl.mean)).^center.Nl.mean;

freqresponse_s = surround.A.mean*exp(-i.*omega.*surround.D.mean).*...
    (1-(surround.Hs.mean./(1+i.*omega.*surround.Ts.mean))).*...
    (1./(1+i.*omega.*surround.Tl.mean)).^surround.Nl.mean;
% Below from http://www.mathworks.com/matlabcentral/newsreader/view_thread/5658
% Getting rid of the highest half of the frequency spectrum
if (rem(n,2))
     freqresponse_c((n+3)/2:n)=[]; centerkernel=ifft([freqresponse_c conj(freqresponse_c((n+1)/2:-1:2))]);
     freqresponse_s((n+3)/2:n)=[]; surroundkernel=ifft([freqresponse_s conj(freqresponse_s((n+1)/2:-1:2))]);
else
    freqresponse_c((n+4)/2:n)=[]; centerkernel=ifft([freqresponse_c conj(freqresponse_c(n/2:-1:2))]);
    freqresponse_s((n+4)/2:n)=[]; surroundkernel=ifft([freqresponse_s conj(freqresponse_s(n/2:-1:2))]);
end
% Note: integral of kernel doesn't depend on delta-T (max TF)
x = [0:1/(TF(end)):1/(TF(end))*(length(centerkernel)-1)]; 
% Due to the bit of Mathworks code above, the highest frequency we're considering is max(TF/2)
% That's why sampling is is steps of 1/(TF(end)) instead of 1/(2*TF(end)).

% Getting rid of "white space"
centerkernel(x>0.5) = [];
surroundkernel(x>0.5) = [];
x(x>0.5) = [];
deltaT = x(2)-x(1);

% Loading an IsoSamp file so I can get the contrasts and baseline firing
% rate
stro = nex2stro(findfile(filename));
Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_l'));
Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_m'));
TF = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
dur = mean(stimoff_t-stimon_t);
spikeidx = strcmp(stro.sum.rasterCells(1,:),getSpikenum(stro));

uniquestim = sortrows(unique([Lcc Mcc TF],'rows'),3); % sorting by TF
uniquestim(all(uniquestim == 0,2),:) = []; % Getting rid of blank
Lblank = Lcc == 0 & Mcc == 0 & TF == 0;

baselinespikerates = [];
for idx = find(Lblank)'
    spiketimes = stro.ras{idx,spikeidx};
    nspikes = sum(spiketimes > stimon_t(idx) & spiketimes < stimon_t(idx)+dur);
    baselinespikerates = [baselinespikerates; nspikes./(stimoff_t(idx)-stimon_t(idx))];
end
baseline = mean(baselinespikerates); % sp/s

nbins = round(.666/deltaT); % time bins for stimulus
temporalenvelope = ones(1,nbins);
temporalenvelope(1:round((1/4)*nbins)) = linspace(0,1,round((1/4)*nbins));
temporalenvelope(end:-1:round((3/4)*nbins)+1) = linspace(0,1,round((1/4)*nbins));
timebins = 0:deltaT:1; % time in s
temporalenvelope = [zeros(1,length(timebins)-length(temporalenvelope)),...
    temporalenvelope];
exampletrial = stro.trial(1,:);
stro.trial = []; stro.ras = [];
for i = 1:size(uniquestim,1)
    tf = uniquestim(i,3);
    centerstim = uniquestim(i,1).*temporalenvelope.*sin(2*pi*tf.*timebins);
    surroundstim = uniquestim(i,2).*temporalenvelope.*sin(2*pi*tf.*timebins);
    centerresp = conv(centerstim,centerkernel,'same');
    surroundresp = conv(surroundstim,surroundkernel,'same');
    
    % getting a firing rate for a Monte Carlo Poisson model
    % (Monte Carlo allows me to use IsoSampgetDPrime.m)
    fr = max(baseline+(centerresp-surroundresp),0);
    start_counting_t = .1;
    stop_counting_t = start_counting_t+.666;
    nspikes = zeros(ntrials,size(timebins,2));
    exampletrial(strcmp(stro.sum.trialFields(1,:),'stim_l')) = uniquestim(i,1);
    exampletrial(strcmp(stro.sum.trialFields(1,:),'stim_m')) = uniquestim(i,2);
    exampletrial(strcmp(stro.sum.trialFields(1,:),'tf')) = uniquestim(i,3);
    exampletrial(strcmp(stro.sum.trialFields(1,:),'stimon_t')) = start_counting_t; % Empirical, this should be more principled!
    exampletrial(strcmp(stro.sum.trialFields(1,:),'stimoff_t')) = stop_counting_t;
    for j = 1:ntrials
        nspikes(j,:) = poissrnd(fr.*deltaT);
        stro.trial = [stro.trial; exampletrial];
        spiketimes = [];
        for k = 1:max(nspikes(j,:))
            spiketimes = [spiketimes timebins(nspikes(j,:)>=k)+deltaT.*unifrnd(-.5,.5,1,sum(nspikes(j,:)>=k))];
        end
        stro.ras{size(stro.ras,1)+1,1} = spiketimes;
        stro.ras{size(stro.ras,1),2} = nan;
    end
end
% adding blank trials
exampletrial(strcmp(stro.sum.trialFields(1,:),'stim_l')) = 0;
exampletrial(strcmp(stro.sum.trialFields(1,:),'stim_m')) = 0;
exampletrial(strcmp(stro.sum.trialFields(1,:),'tf')) = 0;
for j = 1:ntrials
    stro.trial = [stro.trial; exampletrial];
    stro.ras{size(stro.ras,1)+1,1} = unifrnd(start_counting_t,stop_counting_t,1,poissrnd(baseline*.666));
    stro.ras{size(stro.ras,1),2} = nan;
end
[uniquestim, dprime, signal, noise] = IsoSampGetDPrime(stro,1);

figure; axes; hold on; set(gcf,'Color',[0 0 0],'InvertHardCopy','off');
set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]);
for stimtype = 1:2 % non-opponent, opponent
    if (stimtype == 1)
        Lstimtype = sign(uniquestim(:,1)) == sign(uniquestim(:,2));
    else
        Lstimtype = sign(uniquestim(:,1)) ~= sign(uniquestim(:,2));
    end
    h = plot(uniquestim(Lstimtype,3),dprime(Lstimtype),'-','Linewidth',2);
    if (stimtype == 1)
        set(h,'Color','white');
    else
        set(h,'Color','red');
    end
end
set(gca,'Xscale','log');
xlabel('Temporal frequency (Hz)');
ylabel('d''');
set(gca,'Xlim',[1 60]);
title('Benardete and Kaplan model midget cell')
%%
% Section 12
% Temporal impulse response and noise spectrum
% Taken largely from Juan Anguerya's ConeVolution.m and
% makeConePowerSpectrum (inside DTcones_gh.m).
% I'm not sure how Charlie got the units on the noise axes.

Io = 2250; % half-desensitizing value (in R*/sec, from Juan's paper) [old value was 4500]
Ib = 7000;   % from the monitor initalization routine
gain_dark = 0.32;          % from Juan's paper (approximate, and in units of pA/R*)[old value was 0.16]
gainRatio = 1 ./ (1+(Ib./Io));
gainAtBkgnd = gainRatio .* gain_dark % this brigns us to picoamps per R* 

samplingRate = 2000; % Hz
Stim = zeros(floor(samplingRate/5),1);
Stim(1) = 1;
TimeAxis=1:length(Stim);
TimeAxis=TimeAxis./samplingRate;
FilterCoeffs=[0.6745    0.0216    0.0299    0.5311   34.1814];
Filter=ConeEmpiricalDimFlash(FilterCoeffs,TimeAxis);
Filter = Filter./max(Filter); % So that peak is at gainAtBkgnd

% Temporal impulse response function
figure; axes; hold on; set(gcf,'Color',[0 0 0],'InvertHardCopy','off');
set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]);
plot(TimeAxis,gainAtBkgnd*Filter*1000,'w-','LineWidth',2)

%-------------------------
% Now the noise spectrum
% Code taken from Juan's "ConeNoise.m"
%-------------------------

LorentzCoeffs=[0.2   30    2.0    0.05  180    2.5];

% Number of points to be modelled = NoiseLength*samplingRate
% deltaFreq = 1/NoiseLength
FreqAxis_PS=((1:NoiseLength*samplingRate)-1)./(NoiseLength);
Nyquist=samplingRate/2;
FreqAxis_PS=FreqAxis_PS(1:find(FreqAxis_PS<=Nyquist,1,'last'));
ModelNoisePS=lorentzsum_poles(LorentzCoeffs,FreqAxis_PS);

figure; axes; hold on; set(gcf,'Color',[0 0 0],'InvertHardCopy','off');
set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]);
plot(FreqAxis_PS,ModelNoisePS,'w-','LineWidth',2)
set(gca,'YScale','log','Xscale','log')
set(gca,'Xlim',[1 500]);
axis square;

%%
% Section 13
% Rasters for optogenetic stimulation

stro = nex2stro(findfile('F051915014.nex'));

fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targon_t'));
targoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targoff_t'));
opt_stimfreq = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'opt_stimfreq'));
targ_shown = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_shown'));

% Hack below 
if (isnan(targ_shown(1)))
    targ_shown(:,1) = 0;
end



offset = [-.1 .1];  % pre and post time wrt stimon/stimoff
fpacq_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpacq_t'));
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
stimon_t(isnan(stimon_t)) = fpacq_t(isnan(stimon_t))+nanmean(stimon_t-fpacq_t); % Not sur what this hack is for
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
dur = mode(stimoff_t-stimon_t);

sync_t = stimon_t;  % Ecode for alignment
Lspikechans = strncmp(stro.sum.rasterCells,'sig0',4);
optfreq = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'elec_stimfreq'));
uniqfreqs = unique(optfreq);

figure; set(gcf,'Color',[0 0 0],'InvertHardCopy','off','Renderer','painters');
if (any(Lspikechans))
    for whichspike = find(Lspikechans)
        spikes = stro.ras(:,whichspike);
        binwidth = .001; % s
        bins = offset(1):binwidth:dur+offset(2);
        PSTH = zeros(1,length(bins));
        for j = uniqfreqs'
            PSTH = zeros(1,length(bins));
            %figure; subplot(2,1,1); hold on;
            subplot(length(uniqfreqs), 2, 2*(find(uniqfreqs==j))-1); hold on;
            set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]);
            L = optfreq == j;
            trlidxs = find(L);
            for counter = 1:sum(L)
                trlidx = trlidxs(counter);
                tmpspikes = spikes{trlidx}-sync_t(trlidx);
                tmpspikes(tmpspikes < offset(1) | tmpspikes > dur+offset(2)) = [];
                nspikestot = length(tmpspikes);
                plot([tmpspikes tmpspikes]',[zeros(nspikestot,1) .95*ones(nspikestot,1)]'+counter,'y-');
                PSTH = PSTH + hist(tmpspikes, bins);
            end
            PSTH = PSTH./(sum(L).*binwidth);
            if (j > 0) % Plotting the time course of optical stimulation
                if stro.sum.exptParams.modulator
                    f = 0.5*1000^(1/255)^j;
                    t = linspace(bins(1),bins(end),5e4);
                    y = sin(2*pi*f*t-pi/2);
                    y(t<0) = -1;
                    y(t>dur) = -1;
                    y = y-.5; % for plotting below the spikes
                    plot(t,y,'y-','linewidth',2);
                else
                    secspercycle = 1/unique(j(j > 0));                    
                    transitions = 0:secspercycle/2:dur;
                    if (ceil(length(transitions)/2) ~= floor(length(transitions)/2))
                        transitions(end+1) = dur; %automatic shutoff
                    end
                    x = [transitions; transitions];
                    x = [offset(1); x(:); max(x(:))+offset(2)];
                    y = [repmat([0 1],1,length(transitions)/2) 0]*2-2;
                    y = [y;y];
                    if (length(x(:)) == length(y(:)))
                        plot(x,y(:)','y-','linewidth',2);
                    end
                end
            end
            set(gca,'XLim', [offset(1) dur+offset(2)],'Ytick',[],'YLim',[-4 sum(L)+5]);
            if stro.sum.exptParams.modulator && j > 0
                title(['Frequency: ',num2str(round(f)),' Hz'],'Color','white');
            else
                title(['Frequency: ',num2str(j)],'Color','white');
            end
            
            % PSTH
            %subplot(2,1,2); hold on;
            subplot(length(uniqfreqs), 2, 2*(find(uniqfreqs==j))-0); hold on;
            set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]);
            plot(bins,PSTH,'y-','LineWidth',1);
            set(gca,'YLim',[0 10*ceil(max(PSTH(:)/10))+1]);
            set(gca,'Xlim',[offset(1) dur+offset(2)]);
            xlabel('Time (s)','FontSize',12,'Color','white');
            ylabel('Response (sp/s)','FontSize',12,'Color','white');
        end
    end
end

%%
% section 14
% Real cone data from Juan
load('/Users/greghorwitz/Dropbox/Horwitz-Rieke Collaboration/Code/Juan/GaborExample/ExCone_gabor.mat');
subsamples = round(linspace(1,size(ExCone_gabor.Data,2),2000));
for i = 1:3
    figure; axes; hold on; set(gcf,'Color',[0 0 0],'InvertHardCopy','off','Renderer','painters');
    set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]);
    plot(ExCone_gabor.TimeAxis(subsamples),ExCone_gabor.Data(i,subsamples),'-','Color',[.5 .5 .5]);
    set(gca,'Ylim',[-10 10]);
    xlabel('Time (s)');
    ylabel('Current (pA)');
end

figure; axes; hold on; set(gcf,'Color',[0 0 0],'InvertHardCopy','off','Renderer','painters');
set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]);
plot(ExCone_gabor.TimeAxis(subsamples),ExCone_gabor.Mean(subsamples),'-','Color',[.5 .5 .5]);
set(gca,'Ylim',[-10 10]);
xlabel('Time (s)');
ylabel('Current (pA)');
title(['n=',num2str(ExCone_gabor.Mean_n)]);


%%
% Section 15
% Temporal integration control experiments
STIMTYPES = {'RG','LUM'};
conn = database('Nex_Paradigm_Sort','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
colors = [0 1 1; 1 1 0];

for i = 1:length(STIMTYPES)
    figure; set(gcf,'Color',[0 0 0],'InvertHardCopy','off');
    STIMTYPE = STIMTYPES{i};
    if strcmp(STIMTYPE, 'LUM')
        longcomment = 'std inputs 666 take 1';
        shortcomment = 'std inputs 333 take 1';
    else
        longcomment = 'std inputs 666 chr t1';
        shortcomment = 'std inputs chr 333 take 1';
    end
    query = sprintf('SELECT fileID FROM LMTF WHERE subjID = ''A'' AND rfX = 50 AND rfY = 0 AND quality = 1 AND notes LIKE (''%s'')',longcomment);
    flist_long = fetch(conn, query);
    query = sprintf('SELECT fileID FROM LMTF WHERE subjID = ''A'' AND rfX = 50 AND rfY = 0 AND quality = 1 AND notes LIKE (''%s'')',shortcomment);
    flist_short = fetch(conn, query);
    if ~isempty(flist_long)
        data_long = getLMTFrawdata(flist_long);
    else
        data_long= [];
    end
    if ~isempty(flist_short)
        data_short = getLMTFrawdata(flist_short);
    else
        data_short= [];
    end
    
    data = [sqrt(data_long(:,1).^2+data_long(:,2).^2) data_long(:,3) repmat(666,size(data_long,1),1)];
    if ~isempty(data_short)
        data = [data; sqrt(data_short(:,1).^2+data_short(:,2).^2) data_short(:,3) repmat(333,size(data_short,1),1)];
    end
    % Loading the model
    isosamppath = which('IsoSampOnline');
    isosamppath(find(isosamppath==filesep,1,'last'):end) = [];
    load ([isosamppath,filesep,'private',filesep','data',filesep,'LMTF.mat']);
    rfx = unique(data_long(:,5))/10;
    rfy = unique(data_long(:,6))/10;
    
    model = LMTF_global_to_local_model(A.legacy.mode5params, rfx, rfy, 5);
    if strcmp(STIMTYPE,'RG')
        pred = @(omega)1./(model(1+6)*abs(((1i*2*pi*10^model(5+6).*omega+1).^-model(3+6))-model(2+6)*((1i*2*pi*10^(model(5+6)+model(6+6)).*omega+1).^-(model(3+6)+model(4+6)))));
    else
        pred = @(omega)1./(model(1)*abs(((1i*2*pi*10^model(5).*omega+1).^-model(3))-model(2)*((1i*2*pi*10^(model(5)+model(6)).*omega+1).^-(model(3)+model(4)))));
    end
    axes; hold on;
    set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]);
    tmp = logspace(log10(1),log10(20),30);
    plot(tmp,1./pred(tmp),'y-','LineWidth',2);
    plot(tmp,1./(pred(tmp)*sqrt(2)),'c-','LineWidth',2);
    TFs = unique(data(:,2),'sorted');
    durations = unique(data(:,3),'sorted');
    for j = 1:length(TFs)
        for k = 1:length(durations)
            L = data(:,2) == TFs(j) & data(:,3) == durations(k);
            mn = 10.^mean(log10(1./data(L,1)));
            sem = sqrt(var(log10(data(L,1)))/sum(L));
            plot(TFs(j),1./data(L,1),'o','MarkerFaceColor',colors(k,:),'MarkerEdgeColor','black','MarkerSize',4); % Individual data points
            plot(TFs(j),mn,'o','MarkerFaceColor',colors(k,:),'MarkerEdgeColor','black','MarkerSize',14);
            plot([TFs(j) TFs(j)],10.^[log10(mn)+sem log10(mn)-sem],'-','Color',colors(k,:),'LineWidth',4);
        end
    end
    set(gca,'Yscale','log','Xscale','log','Xlim',[.9 30],'Ylim',[2 80]);
    set(gcf,'Renderer','painters');
end
close(conn)

%%
% Section 16
% Rotating ellipsoids (cone model and V1 isoresponse)
% Taken from PBIO2016 section 7
cd ('/Users/greghorwitz/Documents/Manuscripts/Completed/Charlie''s model/For Greg');
load('data_new_LtoM_1.mat');
threshold_pts = gab.colorDirs .* repmat(cones.alpha_analytic,1,size(gab.colorDirs,2));
tmp = [threshold_pts;-threshold_pts];
% getting radii
D = [tmp(:,1) .* tmp(:,1),...
    tmp(:,2) .* tmp(:,2),...
    tmp(:,3) .* tmp(:,3),...
    2*tmp(:,1) .* tmp(:,2),...
    2*tmp(:,1) .* tmp(:,3),...
    2*tmp(:,2) .* tmp(:,3)];
v = (D' * D) \(D' * ones(size(tmp,1),1));
A = [v(1) v(4) v(5);...
    v(4) v(2) v(6);...
    v(5) v(6) v(3)];
[evecs, evals] = eig(A);
evals = diag(evals);
evals = evals([3 2 1]);  % back in LMS order. Fragile
radii = sqrt(1./evals);

[lcc,mcc,scc] = meshgrid(linspace(-radii(1),radii(1),40),...
    linspace(-radii(2),radii(2),40),...
    linspace(-radii(3),radii(3),40));
fn = (lcc./radii(1)).^2+(mcc./radii(2)).^2+(scc./radii(3)).^2;

figure('position', [1165 339 230 454]); axes; hold on;
set(gcf,'Color','none');

p = patch(isosurface(lcc,mcc,scc,fn,1));
set(p,'EdgeColor', 'none', 'FaceAlpha',1,'FaceColor','magenta','Edgealpha',0,'SpecularExponent',1,'SpecularStrength',.1,'DiffuseStrength',1);
lighting gouraud;
axis equal; 
set(gca,'View',[90 0]);
set(gca,'Xtick',[],'Ytick',[],'Ztick',[]);
set(gca,'color','none','visible','off');
axis vis3d;

moviefilename = 'TiltedConeIdlObsThresholds';
if (MAKEMOVIE)
    originalview = get(gca,'View');
    deltaviews = linspace(0,90,60);
    writerObj = VideoWriter(moviefilename);
    open(writerObj);
    for i = 1:length(deltaviews)
        set(gca,'View',originalview+[0 deltaviews(i)]);
        h_light1 = camlight(0,0);
        h_light2 = camlight(20,90);
        h_light3 = camlight(20,-90);

        drawnow;
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
        h_light1.Visible = 'off';
        h_light2.Visible = 'off';
        h_light3.Visible = 'off';

    end
    close(writerObj);
end

%%
% Section 16.5
% Rotating isorepsonse ellipsoid from pan color cell.
% Need to run SFN2010.m section 12 first with these parameters
% "h" is the patch obhect defined in this script.

PLOTMESH = 0;
PLOTFRAME = 0;
PLOTFIRSTTHREE = 0;
MAKEMOVIE = 1;
PLOTPLANE = 0;

SHOWOOG = 0;
WHITEN = 0;
USE10DEG = 1;

filename = 'K072109002.nex'; CELLTYPE = 'pan color'; % elliptical cross section 
filename = 'K082609010.nex'; CELLTYPE = 'pan color'; % Nearly circular cross section
filename = 'K082509007.nex'; CELLTYPE = 'linear';
if strcmp(CELLTYPE,'pan color')
    if strcmp(filename,'K082609010.nex')
        SCALEFACTORS = [.25 .25 .25]; % Good for pancolor cell
    end
    if strcmp(filename,'K072109002.nex')
        SCALEFACTORS = [.5 .5 .5]; % Good for pancolor cell
    end

    PLOTQUAD = 1;
    STARTANGLE = -30;
    ENDANGLE = 360+STARTANGLE;
    ELEVVIEWANGLE = 12;
else
    SCALEFACTORS = [];
    PLOTQUAD = 0;
    PLOTPLANE = 0;
    PLOTFRAME = 1;
    STARTANGLE = 35;
   % ENDANGLE = 494;
   % ELEVVIEWANGLE = 22;
    ELEVVIEWANGLE = 12;
end

%%% NOW RUN SFN2010.m section 12 starting with the nex2stro line
TILT = false;
SPIN = true;
set(gcf,'Color','none');

if strcmp(CELLTYPE,'pan color')
    p = h;
    set(p,'EdgeColor', 'none','Edgealpha',0,'SpecularExponent',1,'SpecularStrength',.1,'DiffuseStrength',1);
    set(gcf,'position', [1165 339 230 454]);
end
axis equal; 
set(gca,'View',[STARTANGLE ELEVVIEWANGLE]);
set(gca,'Xtick',[],'Ytick',[],'Ztick',[]);
set(gca,'color','none','visible','off');
axis vis3d;
set(gca,'Zlim',[-1 1])
axis tight

moviefilename = 'TiltedMonkeyThresholds';
if TILT
     deltaviews = linspace(0,90,60);
     set(gca,'View',[0 ELEVVIEWANGLE]);
end
if SPIN
    deltaviews = linspace(0,360,100);
   % deltaviews = linspace(STARTANGLE,ENDANGLE,60);
end
if (MAKEMOVIE)
    originalview = get(gca,'View');
    writerObj = VideoWriter(moviefilename);
    open(writerObj);
    for i = 1:length(deltaviews)
        if TILT
            set(gca,'View',originalview+[0 deltaviews(i)]);
        end
        if SPIN
            set(gca,'View',originalview+[deltaviews(i) 0]);
        end
        h_light1 = camlight(0,0);
        h_light2 = camlight(20,90);
        h_light3 = camlight(20,-90);

        drawnow;
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
        h_light1.Visible = 'off';
        h_light2.Visible = 'off';
        h_light3.Visible = 'off';

    end
    close(writerObj);
end


%%
% Section 17 
% Gabor movies at the monkeys' detection threshold
RFX = 5; % DVA
RFY = 0; % DVA
STIMTYPE = 'LUM';
filenamestem = ['ThreshGabor',STIMTYPE];
MAKEMOVIE = 1;
if strcmp(STIMTYPE,'LUM')
    maxTF = 15;
    conecoefs = [1 1 0];
elseif strcmp(STIMTYPE,'RG')
    maxTF = 15;
    conecoefs = [1 -1 0];
end

M = [0.0761    0.1524    0.0218;
    0.0275    0.1582    0.0321;
    0.0024    0.0118    0.1220];  % Something standard; I think this is Dell 4
lambda = 1;
bkgndrgb = [.5 .5 .5];
bkgndlms = M*bkgndrgb';
% Loading behavioral model data
isosamppath = which('IsoSampOnline');
isosamppath(find(isosamppath==filesep,1,'last'):end) = [];
load ([isosamppath,filesep,'private',filesep','data',filesep,'LMTF.mat']);
% "A" is the structure of LMTF model parameters for Apollo
behavioralmodel = LMTF_global_to_local_model(A.legacy.mode5params, RFX, RFY, 3);
% lum
if strcmp(STIMTYPE,'LUM')
    behavioral_tcsf = @(omega)(behavioralmodel(1)*abs(((1i*2*pi*10^behavioralmodel(5).*omega+1).^-behavioralmodel(3))-behavioralmodel(2)*((1i*2*pi*10^(behavioralmodel(5)+behavioralmodel(6)).*omega+1).^-(behavioralmodel(3)+behavioralmodel(4)))));
elseif strcmp(STIMTYPE,'RG')
    behavioral_tcsf = @(omega)(behavioralmodel(7)*abs(((1i*2*pi*10^behavioralmodel(11).*omega+1).^-behavioralmodel(9))-behavioralmodel(8)*((1i*2*pi*10^(behavioralmodel(11)+behavioralmodel(12)).*omega+1).^-(behavioralmodel(9)+behavioralmodel(10)))));
end
TFs = linspace(1,maxTF,5);
ccs = 1./behavioral_tcsf(TFs);

figure('units','inches','position',[1 1 4 3]);
axes('units','inches','position',[0 0 4 3],'Visible','off','color',bkgndrgb)

for i = 1:length(ccs) % One loop for each stimulus movie we're going to make

    gaborccs = sqrt(ccs(i).^2/sum(abs(conecoefs))).*conecoefs;
    gaborlms = bkgndlms'.*(gaborccs+1);
    FRAMERATE = 30; % 30 Hz. default for VideoWriter. 
    nframespercycle = FRAMERATE/TFs(i); % (frames/sec)/(cycles/sec) = frames/cycle;
    ncycles = 36/nframespercycle;
    gaborrgb = inv(M)*gaborlms'
    phis = linspace(0,2*pi*ncycles,nframespercycle*ncycles+1);
    phis(end) = [];
    envelope = [linspace(0,1,round(length(phis)/4)) ones(1,round(length(phis)/2))];
    envelope = [envelope, linspace(1,0,length(phis)-length(envelope))];
  
    if (MAKEMOVIE)
        moviefilename = [filenamestem,num2str(round(ncycles))];
        writerObj = VideoWriter(moviefilename);
        open(writerObj);
    end
    for j = 1:length(phis)
        im = DrawGaborEdge(bkgndrgb, gaborrgb-bkgndrgb', [0 0 0], pi, lambda, .15, 1, phis(j), 0, 0, 0, 0, .999, 30);
        image((im-bkgndrgb(1))*envelope(j)+bkgndrgb(1));
        set(gcf,'Color',bkgndrgb);
        set(gca,'XTick',[],'YTick',[],'visible','off');
        axis equal;
        drawnow;
        if (MAKEMOVIE)
            frame = getframe(gca);
            writeVideo(writerObj,frame);
        end
    end
    if (MAKEMOVIE)
        close(writerObj);
    end
end
