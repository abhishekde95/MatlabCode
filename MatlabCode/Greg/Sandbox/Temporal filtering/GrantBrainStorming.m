% Contents
%
% 1: Making time-modulated pulse trains for a microstimulation aim
%
% 2: Playing around with the temporal contrast sensitivity model of
% Watson, A. B. (1986). Temporal sensitivity. Handbook of perception and 
% human performance, 1, 6-1.
%
% 3: Fitting temporal contrast sensitivity functions (6 free
% parameters, difference of two low-pass filters)
%
% 4: Fitting a 2-D surface to LMTF data. A 2-D generalization of the Watson
% 1986 model in which sensitivity as a function of color direction (at a single 
% temporal frequency) is constrained to be an ellipse.
%
% 4.1 Fitting the 2-D Watson surface to LMTF data at individual retinal
% locations.
%
% 4.10 Playing with projections onto two non-orthogonal mechanisms and
% getting the orientation of the threshold ellipse. This will be important
% for an upgrade of tf_fiterr2.m 4/30/17.
%
% 4.11 Switching from a 2-D Watson model in which all of the ellipses are
% identically oriented from the bottom of the funnel to the top to one in
% which the sensitivity of luminance and chromatic mechanisms (not
% necessarily orthogonal) varies across temporal frequencies. This leads to
% ellipses that rotate from the bottom to the top of the funnel. 4/30/17
%
% 4.15 Fitting the 2-D Watson surface across retinal locations assuming
% that only the xi_lum, xi_rg, and theta parameters change with
% eccentricity, and these are allowed to change in an unconstrained way.
%
% 4.16 Using the best-fit parameters of each model as an initial guess for
% each other model. OUTDATED. Works, but is unnecessary and takes a long
% time.
%
% 4.2 Looking at the fitted model parameters as a fuction of retinal
% location using SVD.
%
% 4.21 Fitting a (optionally tilted) rampy trough to the xi_lum and xi_rg
% parameters. GDLH 6/10/17
%
% 4.3 Comparing the ability of one LMTF model (at one retinal position)
% to predict data from another model (at another retinal position).
%
% 4.5 Cross-sections through the parametric and non-parametric surface fits
%
% 5: a rotating movie of the LMTF surface fit (For Rieke lab meeting)
%
% 6: a few drifting Gabor movies (For Rieke lab meeting)
%
% 7: iteratively sampling a sphere (or super-sphere with exponent > 2)
% using the angle between adjacent normal vectors. Relies on convex hull.
% Works well when the surface is convex (exponent > 1) but not when the
% surface has concavities (exponent < 1). LAME.
%
% 8: A misguided attempt to find the boundary of a shape by computing the
% convex hull and slowly whitling away at it until all the points are on
% the boundary. This did not work because sometimes there isn't a union of
% edges from a delaunay triangulation that contains all the points that I
% want. LAME.
%
% 9: Finally realied I can triangulate the points by projecting them onto the unit
% cylinder, convhull, and projecting back. In this iteration, I cycle
% through the edges to find pairs of triangles with very different normal
% vectors. Then I sample in the middle of these triangles.
%
% 10: Like section 9 but instead of finding edges across which the normal
% vectors change by a lot, defining the "surprise" of each point by leaving
% it out of the triangulation and seeing how far off the prediction is.
% Computationally intensive and didn't work very well.
%
% 11: Just sampling the largest faces of the boundary.
%
% 12: Sampling uniformly on the boundary
%
% 13: Sampling the largest face of a threshold surface on the theta, TF domain.
% The rationale here was that we want to sample areas that haven't been
% sampled (e.g. large triangles) and areas where the function value is
% changing a lot (e.g. large triangles). I don't know why this didn't work
% better than it did. Hard ot get it to sample the skinny edge of a
% flattened cylinder.
%
% 14: Hacking around with Gaussian processes to sampling a 1-D temporal
% contrast sensitivity function.
%
% 15: Playing around (manually) with the hyperparameters to see how they
% affect the likelihood of the fit.
%
% 16: Drawing samples from a 2D GP.e
% 
% 17: Implementing a 2D GP fit to data that's periodic in one dimension but
% not the other.
%
% 18: Simulating an LMTF experiment in which we use GP regression to
% figure out where to sample next.
%
% 19: Adding random phase to a sinusoid (Aim 3 backup plan)
%
% 20: Added sinusoidally modulated simulated parvocellular neuronal
% responses with noisy phases to see if that might cancel in V1.
%
% 21: Propixx stuff
%
% 22: Some analysis of ThreePulseMonte code (L-M vs L+M discrimination)
%
% 23: Laminar probe stuff
%
% 24: Playing around with space-time FFTs of gratings with and without eye
% movements.
%
% 25: Looking at observer-to-observer differences in relative sensitivity
% to L+M and L-M (using LMTF data?)
%
% 26: Hacking around with Charlie's cone model (IsoSamp data)
%
% 27: Seeing whether I can get good fits to LMTF data by changing the 
% time constants on the filters instead of the gains. This was Fred's 
% suggestion.
%
% 28: Playing around with threshold ellipses (sensitivity ellipses)?
%
% 29: What does a tilted rampy trough look like in (x,y)? Does it match
% descriptions of cone density? RGC density?
%
% 30: Fitting the Packer et al cone density data
%
% 31: Fitting the Watanabe and Rodieck dendritic field sizes
% -----------------------------------------------------------------
% 32: How would optogenetic silencing of V1 be expected to 
% affect detetion performance under spatial uncertainty?
%
% 33: Just like section 32 but a much simpler model.
%
% 34: Pulling out a bunch of DToneloc files and looking for something to
% provide preliminary data for Aim 2.
%
% -----------------------------------------------------------------
% 35: Superposition of L+M and L-M gratings with different spatial phases
% (in phase = cone isolating = consistent with a material edge(?))
%
% 36: Rasters from a directly excited mDLX5/6 ChR2-expressing neuron
%
% 37: Histogram of optogenetic stimulation latencies
%
% 38: Energy computation between luminance tuned simple cell and
% color-luminance tuned DO cells. Does this predict phase sensitivity?
%% Section 1: Playing around with modulated microstimulation pulse trains

x= zeros(1000,1); % 1 sec at 1 kHz
%y(1:5:end) = 1; % constant frequency 200 Hz train
freq = 100; % in Hz
pulsedensity = (sin(freq*2*pi*[0:1000]/1000)+1)/2;

% Strategy: make a continuoue time non-homogeneous process and then
% quantize time.
carrierfreq = 250;
spiketimes = [1:1000/carrierfreq:1000]/1000; % 200 pulses in 1 s
% warping time
% What is the integral of (sin(ax)+1)/2?
% It's (x-cos(ax)/a)/2
NORMFACT = 2; % scaling the amplitude of the sinewave so we don't get > 1 pulse/time bin
c = freq.*2.*pi;
z = 0.5*(spiketimes-cos(spiketimes.*c)/c);
z = (spiketimes-cos(spiketimes.*c)/c/NORMFACT);
z=z-z(1);

figure; axes; hold on;
plot(spiketimes,[1:length(spiketimes)],'b.-'); % these are evenly spaced
xlabel('time (s)');
ylabel('total # pulses');
plot(z,1:length(z),'m.-');

figure;
axes; hold on;
modulated_pulses = x;
modulated_pulses(round(z*1000)+1) = 1;
steady_pulses = x;
steady_pulses(round(spiketimes*1000)+1) = 1;
plot(repmat(find(modulated_pulses),1,2)',repmat([0 1],length(find(modulated_pulses)),1)','k.-')
set(gca,'Xlim',[100 350]);
set(gca,'Ylim',[0.5 1.5]);
title([num2str(freq),' Hz']);
xlabel('Time (ms)');
[sum(steady_pulses) sum(modulated_pulses)]


%% Section 2: Watson temporal sensitivity model
xi = 1;
zeta = 1;
% filter 1
tau1 = .003;
n1 = 3;
% filter 2
kappa = 1.5;
tau2 = kappa*tau1;
n2 = 4;

% General parameters
nsamps = 1000;
t = linspace(0,1,nsamps); % time in s

% First the impulse response functions. Equation 40 from Watson 1986
% This filter is normalized to fixed area.
h1 = (tau1*factorial(n1-1)).^-1 .*(t./tau1).^(n1-1).*exp(-t./tau1);
h2 = (tau2*factorial(n2-1)).^-1 .*(t./tau2).^(n2-1).*exp(-t./tau2);
h = xi*(h1-zeta*h2);

figure; axes; hold on;
plot(t,h1,'b--');
plot(t,h2,'r--');
plot(t,h,'k-');
xlabel('Time (s)');
ylabel('Amplitude');

% Numerical FFT
Fnum = (1/length(t))*fftshift(fft(h));
figure; axes; hold on;
o = linspace(0,0.5*nsamps/(max(t)-min(t)), length(Fnum(length(Fnum)/2:end)));
plot(o,abs(Fnum(length(Fnum)/2:end)),'b.')
set(gca,'Xscale','log','Yscale','log')

% Analytical FT
omega = logspace(0,2,100);
f1 = (1i*2*pi*tau1.*omega+1).^-n1;
f2 = (1i*2*pi*tau2.*omega+1).^-n2;
fana = xi*(f1-zeta*f2);
plot(omega,abs(fana),'k-');
xlabel('Freq (Hz)');
ylabel('Sensitivity');

%% Section 3: Fitting a 1-D temporal contrast sensitivity function
% Using made-up data. Fitting all 6 of the parameters.

tau1 = .008;
n1 = 9;
tau2 = .005;
n2 = 9;
xi = 4;
zeta = .1;

TFs = [0 .02 .2 1 5 10 25];
lumsens = [10 15 20 30 40 25 10];
rgsens = [50 65 67 70 40 20 5];
figure; axes; hold on;
plot(TFs,lumsens,'ko');
plot(TFs,rgsens,'ro');
options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-8);

% fitting all parameters
for i = 1:2
    if i == 1
        sens = lumsens;
        linestyle = 'k-';
    else
        sens = rgsens;
        linestyle = 'r-';
    end
    
    initparams = [max(sens) .1 n1 n2 tau1 tau2];
    [fpar,fv] = fmincon(@(x) tf_fiterr(x,TFs,sens),initparams,...
        [],[],[],[],[0 0 0 0 0 0],[100 1 20 20 .03 .03],[],options);
    omega = logspace(log10(TFs(2)),log10(TFs(end)),100);
    f1 = (1i*2*pi*fpar(5).*omega+1).^-fpar(3);
    f2 = (1i*2*pi*fpar(6).*omega+1).^-fpar(4);
    f = fpar(1)*(f1-fpar(2)*f2);
    plot(omega,abs(f),linestyle,'linewidth',2);
    set(gca,'XScale','log','YScale','log')
end

%% Section 4: Fitting a 2-D function to LMTF data (threshold as a function of 
% L-cone contrast, M-cone contrast, and temporal frequency.)

flist = flatten(fnamesFromTxt2(fullfile(nexfilepath,'nexfilelists','Greg','LMTF','GregLMTF.txt')));
startpath = '/Volumes/NO BACKUP/NexFiles/Greg/Sedna'
data = [];
for i = 1:length(flist)
    stro = notnex2stro(findfile(flist{i},startpath));
    Lstim_idx = strcmp(stro.sum.trialFields(1,:), 'stim_idx');
    Llcc = strcmp(stro.sum.trialFields(1,:), 'lcc');
    Lmcc = strcmp(stro.sum.trialFields(1,:), 'mcc');
    Ltf = strcmp(stro.sum.trialFields(1,:), 'tf');
    Loog = strcmp(stro.sum.trialFields(1,:), 'oog');
    
    [stim_idxs,init_stim_trial_idxs] = unique(stro.trial(:,Lstim_idx),'last');
    questmodes = stro.trial(init_stim_trial_idxs, Llcc|Lmcc);
    tfs = stro.trial(init_stim_trial_idxs,Ltf);
    
    % Out of gamut checking
    funds = reshape(stro.sum.exptParams.fundamentals,length(stro.sum.exptParams.fundamentals)/3,3);
    if (size(stro.sum.exptParams.mon_spd,1) == 303)
        spds = SplineSpd([380:4:780]',reshape(stro.sum.exptParams.mon_spd,length(stro.sum.exptParams.mon_spd)/3,3),[380:5:780]');    
    else
        spds = SplineSpd([380:2:780]',reshape(stro.sum.exptParams.mon_spd,length(stro.sum.exptParams.mon_spd)/3,3),[380:5:780]');
    end
    M = funds'*spds;
    bkgndrgb = stro.sum.exptParams.bkgndrgb;
    [in_gamut,scalar] = gamutCheck([questmodes zeros(size(questmodes,1),1)]', bkgndrgb, M, 'both');
    questmodes(~in_gamut) = questmodes(~in_gamut).*scalar(~in_gamut);
    data = [data; questmodes tfs ~in_gamut']; % Lcc Mcc TF OOG
end
Loog = logical(data(:,end));

% Stuff we're going to need later
[th,r] = cart2pol(data(:,1),data(:,2));
tf = data(:,3);
LB = [0 0 1 1 .001 .001];
UB = [100 1 20 20 .03 .03];
% Make the filter into a function with a handle
f1 = @(fpar,omega)fpar(1)*abs(((1i*2*pi*fpar(5).*omega+1).^-fpar(3))-fpar(2)*((1i*2*pi*fpar(6).*omega+1).^-fpar(4)));
f2 = @(fpar,omega)fpar(1+6)*abs(((1i*2*pi*fpar(5+6).*omega+1).^-fpar(3+6))-fpar(2+6)*((1i*2*pi*fpar(6+6).*omega+1).^-fpar(4+6)));

% ==========================================
% Now rotating the data a la Patrick to find a rotation that works best

% Trying to switch to fitting threshsolds instead of sensitivities. Sums of
% squared sensitivities look weird when you convert them back to
% thresholds.
x = data(:,1);
y = data(:,2);

initparams = [40 .1 9 3 .005 .002]; % Good for Greg
options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-8);
thetas = linspace(0,pi/2,12);
toterr = [];
fpars = [];
fpar = [initparams, initparams];
for i = 1:length(thetas)
    rotmat = [cos(thetas(i)) -sin(thetas(i)); sin(thetas(i)) cos(thetas(i))];
    % Rotating data clockwise = rotating fit axis counterclockwise.
    xytmp = [x,y]*rotmat;
    [fpar,fv] = fmincon(@(params) tf_fiterr2(params,[xytmp(:,1) xytmp(:,2) tf],Loog),fpar,...
        [],[],[],[],[LB LB],[UB UB],[],options);
    toterr(i) = fv;
    fpars(i,:) = fpar;
end
bestrotidx = find(toterr == min(toterr));
reshape(fpars(bestrotidx,:),[6 2])  % In case the user wants to see the fitted parameters
fpar = fpars(bestrotidx,:);

f1 = @(omega)fpar(1)*abs(((1i*2*pi*fpar(5).*omega+1).^-fpar(3))-fpar(2)*((1i*2*pi*fpar(6).*omega+1).^-fpar(4)));
f2 = @(omega)fpar(1+6)*abs(((1i*2*pi*fpar(5+6).*omega+1).^-fpar(3+6))-fpar(2+6)*((1i*2*pi*fpar(6+6).*omega+1).^-fpar(4+6)));
if (f1(0) < f2(0)) % "If f1 is luminance" then exchange f1 and f2, forcing f1 to be chromatic
    fpar= fpar([[7:12]';[1:6]']);
    f1 = @(omega)fpar(1)*abs(((1i*2*pi*fpar(5).*omega+1).^-fpar(3))-fpar(2)*((1i*2*pi*fpar(6).*omega+1).^-fpar(4)));
    f2 = @(omega)fpar(1+6)*abs(((1i*2*pi*fpar(5+6).*omega+1).^-fpar(3+6))-fpar(2+6)*((1i*2*pi*fpar(6+6).*omega+1).^-fpar(4+6)));
end

% Plotting the raw data 
figure; axes; hold on; 
plot3(x(~Loog),y(~Loog),tf(~Loog),'ko','MarkerFaceColor','black','MarkerSize',5);
plot3(-x(~Loog,1),-y(~Loog),tf(~Loog),'ko','MarkerFaceColor','black','MarkerSize',5);
plot3(x(Loog,1),y(Loog),tf(Loog),'ro','MarkerFaceColor','red','MarkerSize',5);
plot3(-x(Loog,1),-y(Loog),tf(Loog),'ro','MarkerFaceColor','red','MarkerSize',5);
lim = max(max(abs([x y])));
set(gca,'Xlim',1.1*[-lim lim]);
set(gca,'Zscale','log');

set(gca,'Ylim',1.1*[-lim lim]);
axis square
xlabel('L-cone contrast');
ylabel('M-cone contrast');
zlabel('TF (Hz)');
set(gca,'View',[135 20]);
axis vis3d

% Plotting the fit
[xx,yy,zz] = meshgrid(linspace(-max(abs(x)),max(abs(x)),20),...
    linspace(-max(abs(y)),max(abs(y)),20),...
    linspace(min(tf),max(tf),20));
a = abs(f1(zz)).^-1; % chromatic threshold
b = abs(f2(zz)).^-1; % luminance threshold
thtmp = atan2(yy,xx)+thetas(bestrotidx); % clockwise rotation from [L,M] to [a,b]
rtmp = (a.*b)./sqrt((b.*cos(thtmp)).^2+(a.*sin(thtmp)).^2); % radius of ellipse - thank you, Wikipedia

V = sqrt(xx.^2+yy.^2)-rtmp;
FV = isosurface(xx,yy,zz,V,0);
h = patch(FV);
%set(h,'FaceColor','green','EdgeColor','none');
set(h,'FaceColor','green','FaceAlpha',.5','EdgeColor','none','EdgeAlpha',0);
set(gcf,'Renderer','painters');
% Good views 
set(gca,'View',[135 12])
set(gca,'View',[225 12])

% Getting stimuli matched for 1.75 detection threshold for Sedna 5 Hz

tf = 5;
a = abs(f1(tf)).^-1;
b = abs(f2(tf)).^-1;
tmptheta = linspace(0,pi/2,12)';
tmplumrg = [cos(tmptheta) sin(tmptheta)]; % first column is L+M, second is L-M
rotmat = [cos(pi/4) -sin(pi/4);sin(pi/4) cos(pi/4)];
tmplm = tmplumrg*rotmat; % in [L, M]
thtmp = atan2(tmplm(:,2),tmplm(:,1))-thetas(bestrotidx); % clockwise rotation from [L,M] to [a,b]
rtmp = (a.*b)./sqrt((b.*cos(thtmp)).^2+(a.*sin(thtmp)).^2); % radius of ellipse - thank you, Wikipedia
tmplm.*repmat(rtmp,1,2) %

%%
% Section 4.1
% Trying to improve the 2-D model fits both for speed and consistency of
% the parameters (across RF locations).
% Takes a small set of points near the L+M and L-M axes, fits
% 1-D temporal contrast sensitivity functions, and uses these as initial
% guesses for a 2-D fit (inlcuding a fitted theta parameter that is guessed 
% to be pi/4).
% If the FIT flag is set to FALSE this cell skips the fitting and instead
% just plots residuals (model fits and raw data are assumed to be in a
% structure, A, with fields A.models, A.eccs, and A,raw


%textfilepath = fullfile(nexfilepath,'nexfilelists','Greg','LMTF','UtuLMTF.txt');
%textfilepath = fullfile(nexfilepath,'nexfilelists','Greg','LMTF','GregLMTF.txt');
%textfilepath = fullfile(nexfilepath,'nexfilelists','Greg','LMTF','tmpGreg.txt');
%flist = flatten(fnamesFromTxt2(textfilepath));

conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
flist = fetch(conn, 'CALL postPropixxFilenames(''A'')');
flist = fetch(conn, 'SELECT fileID FROM LMTF WHERE subjID=''E'' AND notes IN(''eye tracker'', ''eyetracker'', ''propixx'', ''human button box'')');
close(conn);

startpath = '/Volumes/NO BACKUP/NexFiles';
data = [];
FIT = true; % Set to 1 is you want to do model fitting in this loop. 0 = just plot residuals

for i = 1:length(flist)
    stro = nex2stro(findfile(flist{i},startpath));
    Lstim_idx = strcmp(stro.sum.trialFields(1,:), 'stim_idx');
    Llcc = strcmp(stro.sum.trialFields(1,:), 'lcc');
    Lmcc = strcmp(stro.sum.trialFields(1,:), 'mcc');
    Ltf = strcmp(stro.sum.trialFields(1,:), 'tf');
    Loog = strcmp(stro.sum.trialFields(1,:), 'oog');
    
    [stim_idxs,init_stim_trial_idxs] = unique(stro.trial(:,Lstim_idx),'last');
    questmodes = stro.trial(init_stim_trial_idxs, Llcc|Lmcc);
    tfs = stro.trial(init_stim_trial_idxs,Ltf);
    
    % Out of gamut checking
    funds = reshape(stro.sum.exptParams.fundamentals,length(stro.sum.exptParams.fundamentals)/3,3);
    if (size(stro.sum.exptParams.mon_spd,1) == 303)
        spds = SplineSpd([380:4:780]',reshape(stro.sum.exptParams.mon_spd,length(stro.sum.exptParams.mon_spd)/3,3),[380:5:780]');    
    else
        spds = SplineSpd([380:2:780]',reshape(stro.sum.exptParams.mon_spd,length(stro.sum.exptParams.mon_spd)/3,3),[380:5:780]');
    end
    M = funds'*spds;
    bkgndrgb = stro.sum.exptParams.bkgndrgb;
    [in_gamut,scalar] = gamutCheck([questmodes zeros(size(questmodes,1),1)]', bkgndrgb, M, 'both');
    questmodes(~in_gamut) = questmodes(~in_gamut).*scalar(~in_gamut);
    rfmat = repmat([abs(stro.sum.exptParams.stim_x), stro.sum.exptParams.stim_y],size(questmodes,1),1);
    data = [data;  questmodes tfs ~in_gamut' rfmat]; % Lcc Mcc TF OOG RFx RFy 
end

% Going through the data matrix,  fitting the Watson 2D temporal contrast
% sentivity model for each retinal location, and plotting residuals.
uniqueXYs = unique(data(:,[5 6]),'rows','stable');
if FIT
    models = []; fvs = []; raw = {};
end
for ecc_counter = 1:size(uniqueXYs,1)
    ecc = uniqueXYs(ecc_counter,:);
    Lecc = all(data(:,[5 6]) == repmat(ecc,size(data,1),1),2);
    % Just pulling out the data that's at the right retinal position (ecc).
    data_one_ecc = data(Lecc,:); % From here on out in the loop, using "data_one_ecc" which is the subset of the data at a particular retinal position.
    if (FIT)
        [model, ~, fv] = quickLMTFmodelfit(data_one_ecc(:,[1:4]));
        if (~isnan(model))
           models(:,ecc_counter) = model(:);
           fvs(ecc_counter) = fv;
           raw{ecc_counter} = data_one_ecc;       
        else
            models(:,ecc_counter) = nan;
            fvs(ecc_counter) = nan;
            raw{ecc_counter} = nan;
            continue
        end
    else
        model = A.models(:,ecc_counter); fv = A.fvs(ecc_counter);
    end
    predr = LMTF_thresh_from_model(data_one_ecc(:,[1:3]),model);
    residsPlot(data_one_ecc(:,[1:4]), model,['Error val: ',num2str(fv),' RF: ',num2str(uniqueXYs(ecc_counter,:))]);
    set(gcf,'Position',[5 380 560 420]);
    drawnow;
end
% Cleaning up
L = isnan(models(1,:));
if (sum(L) > 0)
    models(:,L) = [];
    fvs(L) = [];
    raw(L) = [];
    uniqueXYs(L,:) = [];
    olddata = data;
    data = data(ismember(data(:,[5 6]),uniqueXYs,'rows'),:);
end
% Execute the lines below to launch the LTFBrowser
if ~exist('A')
    A = [];
end
A.models = models;
A.eccs = uniqueXYs;
A.raw = raw;
A.fvs = fvs;
LMTFBrowser(A);

%%
% Section 4.10
% In a model with two linear mechanisms that combine with squaring, how
% does the orientation of the isodetection ellipse change with the relative
% sensitivity of the two mechanisms?

mechs = normrnd(0,1,2,2);
%mechs = [1 0; 0 1];
[v,d] = eig(mechs*mechs');
widths = 1./diag(d);
rng = 2*ceil(max(widths));
[xx,yy] = meshgrid(linspace(-rng,rng,100),linspace(-rng,rng,100));

scalefactors = [.5 1 2]; % scaling mechanism 1 only
nscalefactors = length(scalefactors);
figure; axes; hold on;
for i = 1:nscalefactors
    scalemech1mat = [scalefactors(i) 1; scalefactors(i) 1];
    tmpmechs = (mechs.*scalemech1mat); % scaling one of the two mechanisms
    transformedxy = [xx(:) yy(:)]*tmpmechs;
    fv = contour(xx,yy,reshape(sum(transformedxy.^2,2),size(xx)),[1 1]);
    
    [v,d] = eig(tmpmechs*tmpmechs');
    r1 = v(:,1)/sqrt(d(1,1));
    plot([0 r1(1)], [0 r1(2)],'k-');
    r2 = v(:,2)/sqrt(d(2,2));
    plot([0 r2(1)], [0 r2(2)],'k-');
    axis equal;
    drawnow;
   % pause;
end
% Confirmed: changing the relative sensitivity of the two mechanisms
% changes the orientation of the isodetectiod ellipse. Now I just need a
% quick way of figuring out what this orientation is from the two
% mechanisms.

% Below, code to go from "th" (direction of the stimulus
% vector in LM plane) to predicted threshold

mechs = normrnd(0,1,2,2);
[v,d] = eig(mechs*mechs'); % v and d define the detection ellipse
transformedxy = [xx(:) yy(:)]*mechs;
figure; axes; hold on;
fv = contour(xx,yy,reshape(sum(transformedxy.^2,2),size(xx)),[1 1]);

radii = sqrt(1./diag(d));
r1 = v(:,1)*radii(1);
r2 = v(:,2)*radii(2);
th = linspace(0,2*pi,50);

r1theta = atan2(v(2,1),v(1,1)); % angle of the r1 axis (the long axis) of the ellipse relative to L-cone

costerm = bsxfun(@times, radii(2), cos(th-r1theta)).^2;
sinterm = bsxfun(@times, radii(1), sin(th-r1theta)).^2;
f = bsxfun(@rdivide, prod(radii), sqrt(costerm+sinterm));

[x,y] = pol2cart(th,f); % th is with respect to the L-cone axis

plot(x,y,'o')
axis equal

% Finding points on the ellipse by taking the dotproduct of 
% unit vectors onto the mechanisms

unitvectors = [cos(th') sin(th')];
dps = unitvectors*mechs;
sum_squared_dps = sum(dps.^2,2);
[x,y] = pol2cart(th',1./sqrt(sum_squared_dps));

plot(x,y,'m*');

% Now for a given theta (and TF!), how do I get "r"? Is going though unit vectors (above) the
% best way to go? Seems like it is.

%%
% Section 4.11
% Allow ellipses to change orientation with temporal frequency according to
% a "line-element model"(?)

conn = database('Nex_Paradigm_Sort','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
flist = fetch(conn, 'CALL postPropixxFilenames(''E'')');
close(conn);
[data, check] = getLMTFrawdata(flist);
uniqueXYs = unique(data(:,[5 6]),'rows','stable');
for j = 1:size(uniqueXYs,1)
    whichecc = uniqueXYs(j,:);
    Lpos = data(:,5) == whichecc(1) & data(:,6) == whichecc(2);
    data_one_ecc = data(Lpos,:);
    Loog = data_one_ecc(:,4) == 1;

    figure; subplot(2,1,1); hold on;
    plot3(data_one_ecc(~Loog,1),data_one_ecc(~Loog,2),data_one_ecc(~Loog,3),'ko','MarkerSize',5,'MarkerFaceColor','black');
    plot3(-data_one_ecc(~Loog,1),-data_one_ecc(~Loog,2),data_one_ecc(~Loog,3),'ko','MarkerSize',5,'MarkerFaceColor','black');
    plot3(data_one_ecc(Loog,1),data_one_ecc(Loog,2),data_one_ecc(Loog,3),'ro','MarkerSize',5,'MarkerFaceColor','red');
    plot3(-data_one_ecc(Loog,1),-data_one_ecc(Loog,2),data_one_ecc(Loog,3),'ro','MarkerSize',5,'MarkerFaceColor','red');
   
    % Fitting new model and old model
    [initialguess, ~, ~] = quickLMTFmodelfit(data_one_ecc(:,[1:4]));
    options = optimset('Algorithm', 'interior-point', 'MaxFunEvals', 5e4, ...
        'MaxIter', 5e4, 'TolFun', 1e-8, 'Display', 'iter');
    LB = [0   0  1   0 -2  log10(1.1)]; % if delta_n == 1, min(n1) = 2
    UB = [500 1 40   2 -1  log10(2)];
    LBtheta = 0;
    UBtheta = pi/2;

    [model_new,fv_new,~,~,~,~,hessmat] = fmincon(@(params) tf_fiterr2(params,data_one_ecc(:,[1:3]), data_one_ecc(:,4)), initialguess,...
        [],[],[],[],[LB LB LBtheta],[UB UB UBtheta],[],options);

    xlims = get(gca,'Xlim'); ylims = get(gca,'Ylim'); zlims = get(gca,'Zlim');
    [xx,yy,zz] = meshgrid(linspace(.5*xlims(1),.5*xlims(2),60),...
        linspace(.5*ylims(1),.5*ylims(2),60),...
        logspace(log10(1),log10(40),60));
    % Plotting  new surface
    xi_1 = model_new(1);
    zeta_1 = model_new(2);
    n1_1 = model_new(3);
    n2_1 = model_new(3)+model_new(4);
    tau1_1 = 10^model_new(5);
    tau2_1 = 10^(model_new(5)+model_new(6));
    xi_2 = model_new(7);
    zeta_2 = model_new(8);
    n1_2 = model_new(9);
    n2_2 = model_new(9)+model_new(10);
    tau1_2 = 10^model_new(11);
    tau2_2 = 10^(model_new(11)+model_new(12));
    theta = model_new(13);
    f1 = @(omega)xi_1*abs(((1i*2*pi*tau1_1.*omega+1).^-n1_1)-zeta_1*((1i*2*pi*tau2_1.*omega+1).^-n2_1)); % sensitivity
    f2 = @(omega)xi_2*abs(((1i*2*pi*tau1_2.*omega+1).^-n1_2)-zeta_2*((1i*2*pi*tau2_2.*omega+1).^-n2_2));
    tfs = data_one_ecc(:,3);
    a = abs(f1(zz)); % luminance sensitivity
    b = abs(f2(zz)); % chromatic sensitivity
    mechs = [cos(theta) 1/sqrt(2); sin(theta) -1/sqrt(2)];

    % Evaluating the 3-D function (response as a function of L, M, TF)
    % Using isosurface to get funnel
    V = zeros(size(xx));
    for i = 1:numel(xx)
        stimvector = [xx(i) yy(i)];
        V(i) = sqrt(sum((stimvector*mechs*diag([a(i) b(i)])).^2));
    end
    FV = isosurface(xx,yy,zz,V,1);
    h = patch(FV);
    set(h,'FaceColor','blue','FaceAlpha',.5','EdgeColor','none','EdgeAlpha',0);
    set(gcf,'Renderer','painters');
    set(gca,'Zscale','log')
    title(num2str(theta));
    %predr = LMTF_thresh_from_model(data_one_ecc,model_new);
    %r = sqrt(sum(data_one_ecc(:,[1 2]).^2,2));
    %plot(r(~data_one_ecc(:,4)),predr(~data_one_ecc(:,4)),'ko')
    
    % Plotting a few isodetection ellipses
    subplot(2,1,2); hold on;
    for i = 1:size(V,3)
        contour(squeeze(V(:,:,i)),[1 1])
    end 
end

%%
% Section 4.15
% Assuming none of the parameters change across retinal eccentricity except
% xi_LUM, x1_RG and theta. These three parameters are fitted individually for 
% each location; the other ten parameters are assumed to be constant across
% locations. 
% This cell assumes you already have in the workspace:
% "A" (for initial guesses)
% "data" from cell 4.1.

UBtheta = pi/2;
LBtheta = 0;
options = optimset('Algorithm', 'interior-point', 'MaxFunEvals', 5e4, ...
    'MaxIter', 5e4, 'TolFun', 1e-8, 'Display', 'off');
uniqueXYs = A.eccs;
% Need to fit 3 x length(uniqueXYs) parameters
% Need some sort of consensus guess as to the shape of the contrast
% response functions. Evaluating fits on a lattice, scaling to equal
% height, and fitting?
nmodels = size(uniqueXYs,1);
% getting numbers of data points per retinal location
ns = zeros(nmodels,1);
for i = 1:nmodels
   ns(i) = sum(all([data(:,5) == uniqueXYs(i,1), data(:,6) == uniqueXYs(i,2)],2));
end
tfs = logspace(log10(1),log10(60),40);
% Computing weighted (by n) averages of predicted thresholds 
modelfits = zeros(nmodels,length(tfs),2);
for i =1:2 % LUM/RG
    for j = 1:nmodels
        if (i == 1)
            xi = models(1,j);
            zeta = models(2,j);
            n1 = models(3,j);
            n2 = models(3,j)+models(4,j);
            tau1 = 10^models(5,j);
            tau2 = 10^(models(5,j)+models(6,j));
        else
            xi = models(7,j);
            zeta = models(8,j);
            n1 = models(9,j);
            n2 = models(9,j)+models(10,j);
            tau1 = 10^models(11,j);
            tau2 = 10^(models(11,j)+models(12,j));           
        end
        f = @(omega)xi*abs(((1i*2*pi*tau1.*omega+1).^-n1)-zeta*((1i*2*pi*tau2.*omega+1).^-n2));
        modelfits(j,:,i) = f(tfs);
    end
end

%  1D model fits to weighted averages
LB = [0   0  1  0 -3  log10(1.00001)]; % if delta_n == 1, min(n1) = 2
UB = [500 1 40  10 -1  log10(2)];
initparams(1,:) = [20  .9  3 0 -2 log10(1.5)]; % LUM
initparams(2,:) = [100 .1 10 2 -2 log10(1.5)]; % RG

% Plotting slices through individual model fits along the LUM direction (1)
% and along the RG direction (2), weighted averages of these slices + model
% fits to the weighted averages (3).
figure;
wtavg = zeros(2,length(tfs));
params = zeros(2,6); % 2 color directions, 6 parameters for 1D fit
for i = 1:2
    subplot(3,1,i);
    plot(tfs,squeeze(modelfits(:,:,i))');
    wtavg(i,:) = ns'*squeeze(modelfits(:,:,i))./sum(ns);
    params(i,:) = fmincon(@(x) tf_fiterr(x,tfs,wtavg(i,:)),initparams(i,:),...
        [],[],[],[],LB,UB,[],options);
end
subplot(3,1,3); hold on;
plot(tfs,wtavg');
for i = 1:2
    f = @(omega)params(i,1)*abs(((1i*2*pi*10^params(i,5).*omega+1).^-params(i,3))-params(i,2)*((1i*2*pi*10^(params(i,5)+params(i,6)).*omega+1).^-(params(i,3)+params(i,4))));
    plot(tfs,f(tfs),'k:');
end

% Initial guesses for the big minimization problem
hugeparamlist = [params(1,2:end)';params(2,2:end)'];
%  [zeta_LUM, n1_LUM, delta n_LUM, tau1_LUM, kappa_LUM, zeta_RG, n1_RG,...
%  delta n_RG, tau1_RG, kappa_RG];

% Getting initial guesses: just use the best fit parameter
% from the fits to the individual data sets at each retinal location 
for i = 1:nmodels
    hugeparamlist = [hugeparamlist; models(1,i);  models(7,i);  models(13,i)];
end

% Setting up parameter bounds
LB = [0  1  0  -3  log10(1.1) 0  1  0 -3  log10(1.00001)]; 
UB = [1 40  2  -1  log10(2)   1 40  10 -1  log10(2)];
for i = 1:nmodels
   LB = [LB, 0 0 0];
   UB = [UB, 500 500 pi/2];
end
options = optimset(options,'Display','iter','MaxFunEvals',10^6);

% This takes a long time.
[fpar,fv] = fmincon(@(params) tf_fiterr3(params,data,1),hugeparamlist,...
    [],[],[],[],LB,UB,[],options);

% Seeing how the initial guesses changed
ximat = reshape(fpar(11:end),3,nmodels);
initguessmat = reshape(hugeparamlist(11:end),3,nmodels);

figure; 
for i = 1:3
    subplot(3,1,i);
    imagesc([ximat(i,:);initguessmat(i,:)]);
end

% Making a new structure that we can use with LMTFBrowser.
B = A;
for i=1:nmodels
    B.models(:,i) = [ximat(1,i);fpar(1:5); ximat(2,i); fpar(6:10); ximat(3,i)];
end

%%
% Section 4.16
% Fitting each dataset (data at each retinal location) using as an initial 
% guess the best fit from each other data sets. Continue until none of the
% fits improve. No need to repeat fitting individual data sets with models 
% that have already been tried. 
% Workspace needs: uniqueXYs, data, models
% This will take a long time.

options = optimset('Algorithm', 'sqp', 'MaxFunEvals', 5e4,'MaxIter', 5e4, 'TolFun', 1e-8, 'Display', 'off');
LB = [0   0  1   0  -3   log10(1.1)   0    0  1  0 -3   log10(1.1) 0]; % if delta_n == 1, min(n1) = 2
UB = [500 1 40   2  -1   log10(2)    500   1 40  2  -1  log10(2) pi/2];
initialguessidxs = 1:size(uniqueXYs,1);
keepingtrack = 1; % Just need to start with some value in keeping track, it's cleared in two line anyway
waitbar_h = waitbar(0,'Please wait...');
while ~isempty(keepingtrack)
    keepingtrack = [];
    for i = 1:size(uniqueXYs,1) % data comes from retinal location i
        for j = initialguessidxs % initial guess comes from model fit at location j
            fractionalwaythrough = ((i-1)*length(initialguessidxs)+find(j==initialguessidxs))/(size(uniqueXYs,1)*length(initialguessidxs));
            waitbar(fractionalwaythrough,waitbar_h);
            Lecc = all(data(:,[5 6]) == repmat(uniqueXYs(i,:),size(data,1),1),2);
            % Just pulling out the data that's at the right retinal position (ecc).
            [model,fv,~,~,~,~,hessmat] = fmincon(@(params) tf_fiterr2(params,  data(Lecc,[1:3]),data(Lecc,4)), models(:,j),...
                [],[],[],[],LB,UB,[],options);
            if (fv < fvs(i))
                disp('--------------------------------');
                fprintf('Fit to data at (%d, %d) is improved by guessing model params from (%d, %d)\n',uniqueXYs(i,:),uniqueXYs(j,:));
                fprintf('Fitting error decreased from %d to %d\n',fvs(i),fv)
                disp('--------------------------------');
                models(:,i) = model;
                fvs(i) = fv;
                keepingtrack = [keepingtrack; i j];
            end
        end
    end
    % The first column of keepingtrack contains models that were updated in
    % the last round. Use these as initial guesses in nest round.
    if (~isempty(keepingtrack))
        initialguessidxs = unique(keepingtrack(:,1))';
        disp('Improved model fits at retinal locations');
        disp(uniqueXYs(initialguessidxs,:));
    end % otherwise, we're done
end
close (waitbar_h);
disp('Done');
A.models = models;
A.fvs = fvs;
A.eccs = uniqueXYs;

% Checking how many of the parameters hit the boundary
figure; axes; hold on;
plot(models,'.');
plot(LB,'-');
plot(UB,'-');

%%
% 4.2 Looking at the distribution of model parameters across retinal
% space and fitting a model that predicts xi_LUM, xi_RG, and theta
% as a function of retinal position (a quadratic regression in polar
% coordinates).

% Plotting how model parameters change as a function of retinal position
% The regressions should really be weighted regression
models = B.models;
labels = {'xi_{LUM}','zeta_{LUM}','n1_{LUM}','delta_n_{LUM}','tau1_{LUM}','kappa_{LUM}','xi_{RG}','zeta_{RG}','n1_{RG}','delta_n_{RG}','tau1_{RG}','kappa_{RG}','theta'};
Lnonzerovar = diag(cov(models')) > 10^-20;
covmat = cov(models(Lnonzerovar,:)');
prettycorr(models(Lnonzerovar,:)',labels(Lnonzerovar));
[phi,r] = cart2pol(abs(uniqueXYs(:,1))/10,uniqueXYs(:,2)/10);

for i = 1:size(models,1);
    tmpparams = models(i,:);
    if all(tmpparams == tmpparams(1))
        continue;
    end
    figure; subplot(2,2,1); hold on;
    for j = 1:size(uniqueXYs,1)
        %h = plot(uniqueXYs(j,1)/10,uniqueXYs(j,2)/10,'ko','MarkerFaceColor','black','MarkerSize',max(1,10*(tmpparams(j)-min(tmpparams))./(max(tmpparams)-min(tmpparams))));
        plot3(uniqueXYs(j,1)/10,uniqueXYs(j,2)/10,tmpparams(j),'ko','MarkerFaceColor','black','MarkerSize',2);
    end
    xlabel('X (dva)');
    ylabel('Y (dva)');
    axis equal square;
    
    subplot(2,2,2); hold on;
    for j = 1:size(uniqueXYs,1)
        plot3(phi(j),r(j),tmpparams(j),'ko','MarkerFaceColor','black','MarkerSize',2);
    end
    xlabel('angle from horizontal (rad)');
    ylabel('radius (deg)');
    axis equal square;
    [b,bint,~,~,stats] = regress(log10(tmpparams'),[ones(size(uniqueXYs,1),1) [r phi phi.^2]]);
    [phimesh,rmesh] = meshgrid(-pi/2:pi/50:pi/2,0:.1:max(r));
    predprojs = 10.^([ones(size(rmesh(:))), rmesh(:), phimesh(:) phimesh(:).^2]*b);
    surface(phimesh,rmesh,reshape(predprojs,size(rmesh)),'EdgeColor','none','FaceAlpha',.5);
    set(gca,'Zscale','log');
    
    subplot(2,2,3); hold on; axis square;
    parameter_labels = {'radius','phi','phi^2'};
    for k = 1:length(parameter_labels)
       text(0.2,k,[parameter_labels{k},': ',num2str(sign(bint(k+1,1)) == sign(bint(k+1,2)))])
    end
    title('individual significance tests');
    set(gca,'ylim',[0 k+1],'yticklabels',[]);
    subplot(2,2,1);
    [xmesh, ymesh] = pol2cart(phimesh, rmesh);
    surface(xmesh,ymesh,reshape(predprojs,size(rmesh)),'EdgeColor','none','FaceAlpha',.5);
 
    if stats(3) < 0.05
        title([labels{i},' (Significant) '])
    else
        title(labels{i})
    end
end

%%
% Section 4.21
% Trying a new version of the tilted rampy trough model
% b0+b1*r+b2*r*cos(phi)+b3*r*sin(phi)
dirpath = fileparts(which('IsoSampOnline'));
load ([dirpath,filesep,'private\data\LMTF.mat']);
SID = 'U';
s = eval(SID);
eccs = s.eccs;
ximat = log10(s.legacy.mode0models([1 7],:));
[phi,r] = cart2pol(eccs(:,1)/10,eccs(:,2)/10);

figure;
for i = 1:2
    X = [ones(size(eccs,1),1) r r.*cos(2.*phi) r.*sin(2.*phi)];
    b = regress(ximat(i,:)',X)
    [rs,phis] = meshgrid(linspace(2,15,20),linspace(-pi/2,pi/2,20));
    pred = [ones(size(rs(:),1),1) rs(:) rs(:).*cos(2.*phis(:)) rs(:).*sin(2.*phis(:))]*b;
    subplot(2,2,(i-1)*2+1); hold on;
    plot3(r,phi,ximat(i,:),'o');
    surf(rs,phis,reshape(pred,size(rs)));
    set(gca,'View',[90 0]);
    subplot(2,2,(i-1)*2+2);
    [xs,ys] = pol2cart(phis,rs);
    surf(xs,ys,reshape(pred,size(rs)));
    set(gca,'View',[0 90]);
    axis equal;
end
subplot(2,2,1); title([SID,': LUM'],'FontSize',12);
subplot(2,2,3); title([SID,': RG'],'FontSize',12);
%%
% Section 4.3) Taking models from the LMTF.mat structure and using them to 
% predict thresholds measured at other retinal locations.
% NOTE: OUT OF GAMUT POINTS ARE SIMPLY IGNORED IN THIE SCRIPT
% Default is to look at all the residuals. It is is scientifically
% interesting to dissect the differences in the surfaces across retinal
% position. Is it mostly luminance sensitivity that's changing? Color? Low
% TF sensitity? Maybe this is a good rotaton project. OUTDATED.
load /Users/greghorwitz/Desktop/MatlabCode/Zack/IsoSamp/private/data/LMTF.mat

LMTFstruct = A;
neccs = size(LMTFstruct.models,2);
dismat = zeros(neccs,neccs); % dissimilarity matrix
for i = 1:neccs
    for j = 1:neccs
        model = LMTFstruct.models(:,i);
        data_to_fit = LMTFstruct.raw{:,j};
        [~,r_measured] = cart2pol(data_to_fit(:,1),data_to_fit(:,2));
        r_pred = LMTF_thresh_from_model(data_to_fit(:,[1:3]),model);
        % completely ignoring OOG points.
        Loog = data_to_fit(:,4);
        logresid = log10(r_pred(~Loog))-log10(r_measured(~Loog));
        mediansquared= median(logresid.^2); 
        dismat(i,j) = mediansquared;
    end
end


% Averging the pairwise errors across triplets of retinal locations that
% are grouped together by a Delaunay triangulation.
DT = delaunayTriangulation(LMTFstruct.eccs');
figure; axes; hold on; 
selfdistances = diag(dismat);
z = zeros(size(DT.ConnectivityList,1),1)
% Plotting patches
for i = 1:size(DT.ConnectivityList,1)
    nodeidxs = DT.ConnectivityList(i,:);
    points = DT.Points(nodeidxs,:);
    L = ismember(LMTFstruct.eccs',points,'rows');
    subdismat = dismat(find(L),find(L));
    errs = [subdismat(find(~tril(ones(size(subdismat)))));...
        subdismat(find(~triu(ones(size(subdismat)))))];
    meanerr = mean(errs);
    z(i) = (mean(errs) - mean(selfdistances))./std(selfdistances);
    h = patch(points(:,1),points(:,2),z(i));
end
colormap(jet(ceil(max(z))));
cmap = colormap;
colorbar

% Plotting edges
figure; axes; hold on;
edgeorder = [1 2; 1 3; 2 3];
cmap = colormap(jet(20));
z = [];
for i = 1:size(DT.ConnectivityList,1)
    nodeidxs = DT.ConnectivityList(i,:);
    for j = 1:3
        pairidx = nodeidxs(edgeorder(j,:));
        points = DT.Points(pairidx,:);
        meanerr = mean([dismat(pairidx(1),pairidx(2)); dismat(pairidx(2),pairidx(1))]); % The "distance" between the fit of 'i' and data of 'j' & the fit of 'j' and data of 'i'
        z(i) = (meanerr - mean(selfdistances))./std(selfdistances);
        coloridx = max(min(size(cmap,1),round(z(i))),1); % 1<coloridx<20 
        h = plot(points(:,1),points(:,2),'color',cmap(coloridx,:),'linewidth',3);
    end
end

%%
% Section 4.5
% Plotting various cuts through the fitted 2D threshold function

% GP fitting
sqrtn = 40;
[x1,x2] = meshgrid(logspace(log10(1),log10(25),sqrtn),linspace(0, pi,sqrtn));
x_star = [x1(:) x2(:)]; % column order: TF, theta

% Setting up mean, covariance, and likelihood functions + hyperparameters
covfunc = {'covProd',  {{'covMask',{[1 0],'covSEiso'}}, {'covMask',{[0 1],'covPeriodPi'}}}};
meanfunc = @meanConst;
likfunc = @likGauss;

hyp.cov = [0;1;1;-1]; % TF length, TF amplitude, Theta length, Theta amplitude
hyp.mean = -1;
hyp.lik = -3;

% Not using the OOG points to get estimate smoothness, error.
hyp2 = minimize(hyp, @gp, -100,@infExact, meanfunc, covfunc, likfunc, [log10(data(~Loog,3)),th(~Loog)],log10(r(~Loog)));
%hyp2 = hyp;
[th,r] = cart2pol(data(:,1),data(:,2));

% 1-D plots in TF. 
% First the parametric surface
omega = logspace(log10(1),log10(25),100);
if f1(25) < f2(25) % f1 is sensitivity
    ftmp = f1;
    f1 = f2;
    f2 = ftmp;
    clear ftmp;
end
y1 = f1(omega); %lum
y2 = f2(omega);

figure; axes; hold on;
plot(omega,y1,'k-');
plot(omega,y2,'r-');
set(gca,'Xscale','log','Yscale','log','Xlim',[omega(1) omega(end)]);
ylabel('Sensitivity (cone contrast)');
xlabel('Temporal frequency (Hz)');

% Second, superimposing GP contrast sensitivity curves
angle1 = thetas(bestrotidx);
angle2 = mod(thetas(bestrotidx)+pi./2,pi);
y1 = gp(hyp2, @infExact, meanfunc, covfunc, likfunc, [log10(data(~Loog,3)) th(~Loog)],log10(r(~Loog)),[log10(omega'), repmat(angle1,length(omega),1)]);
plot(omega,1./10.^y1,'k--');
y2 = gp(hyp2, @infExact, meanfunc, covfunc, likfunc, [log10(data(~Loog,3)) th(~Loog)],log10(r(~Loog)),[log10(omega'), repmat(angle2,length(omega),1)]);
plot(omega,1./10.^y2,'r--');


% Threshold surface
[m,s2,pm,ps2] = gp(hyp2, @infExact, meanfunc, covfunc, likfunc, [log10(data(~Loog,3)) th(~Loog)],log10(r(~Loog)),[log10(x_star(:,1)), x_star(:,2)]);
% figure; axes; hold on;
% surf(reshape(x1,sqrtn,sqrtn),reshape(x2,sqrtn,sqrtn),reshape(pm,sqrtn,sqrtn))
% plot3(data(~Loog,3),mod(th(~Loog),pi),log10(r(~Loog)),'ko','MarkerSize',10,'MarkerFaceColor','black');
% plot3(data(Loog,3),mod(th(Loog),pi),log10(r(Loog)),'ro','MarkerSize',10,'MarkerFaceColor','red');

% Sensitivity plot 
figure; axes; hold on;
surf(reshape(x1,sqrtn,sqrtn),reshape(x2,sqrtn,sqrtn),reshape(1./10.^pm,sqrtn,sqrtn))
plot3(data(~Loog,3),mod(th(~Loog),pi),1./r(~Loog),'ko','MarkerSize',10,'MarkerFaceColor','black');
plot3(data(Loog,3),mod(th(Loog),pi),1./r(Loog),'ro','MarkerSize',10,'MarkerFaceColor','red');
set(gca,'YTick',[0 pi/4 pi/2 3*pi/4],'YTickLabel',{'L','L+M','M','L-M'});
colormap(hot);
xlabel('TF (Hz)');
zlabel('Sensitivity');
set(gca,'Xscale','log','Zscale','log');
set(gca,'Xlim',[omega(1) omega(end)],'Ylim',[0 pi]);

% Picking a slice to plot in the LM plane
for whichfit = 1:2
    figure; plotcounter = 1;
    nbins = 9;
    bins = logspace(log10(min(data(:,3))), log10(max(data(:,3))),nbins+1);
    for i = 1:length(bins)-1
        TFbounds = [bins(i) bins(i+1)];
        L = data(:,3) >= TFbounds(1) & data(:,3) <= TFbounds(2);
        subplot(ceil(sqrt(nbins)),ceil(sqrt(nbins)),plotcounter); hold on;
        plot(data(L&~Loog,1),data(L&~Loog,2),'k.');
        plot(-data(L&~Loog,1),-data(L&~Loog,2),'k.');
        plot(data(L&Loog,1),data(L&Loog,2),'r.');
        plot(-data(L&Loog,1),-data(L&Loog,2),'r.');
        x_star = [repmat(mean(TFbounds),100,1) linspace(0,2*pi,100)'];
        if (whichfit == 1)  % GP
            [m,s2,pm,ps2] = gp(hyp2, @infExact, meanfunc, covfunc, likfunc, [data(~Loog,3) th(~Loog)],log10(r(~Loog)),x_star);
            [x,y] = pol2cart(linspace(0,2*pi,100)',10.^m);
            set(gcf,'Name','Gaussian Process fit')
        else  % Parametric
            y1 = f1(mean(TFbounds)).^-1; %lum thresholds
            y2 = f2(mean(TFbounds)).^-1;
            thtmp = linspace(0,2*pi,100)'
            rtmp = (y1.*y2)./sqrt((y2.*cos(thtmp)).^2+(y1.*sin(thtmp)).^2); % radius of ellipse - thank you, Wikipedia
            [x,y] = pol2cart(thtmp+thetas(bestrotidx),rtmp);
            set(gcf,'Name','Watson fit')
        end
        plot(x,y,'b-');
        title(['TF: ',num2str(geomean(TFbounds))]);
        axis square;
        if (sum(L) > 0)
            lim = max(max(abs(data(L,[1 2]))));
        else
            lim = max(abs([x;y]));
        end
        if isempty(lim)
           lim = max(abs([x;y]));
        end
        set(gca,'Xlim',lim*[-2 2],'Ylim',lim*[-2 2]);
        plotcounter = plotcounter+1;
    end
end
%%
% Section 5: of the fitted LMTF surface
% Making a movie of a rotating surface

opengl('software');  % To avoid getframe bug: http://www.mathworks.com/support/bugreports/384622
%set(gcf,'Renderer','zbuffer') 
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
mpgwrite(M, gray, 'LMTFMovie2.mpg', options);

%%
% Section 6 
% A few Gabor movies
bkgndrgb = [.5 .5 .5];
gaborrgb = [.3 .3 .1]; % could be better
%gaborrgb = [-.2 .1 -0.01]; % could be better
deltaphi = .1;
edgergb = [0 0 0];
theta = 0;
lambda = 2;
sigma = 1;
gamma = 1;
phi = 0;
xoff = 0;
yoff = 0;
etheta = 0;
edisp = 0;
gausslim = .99;
pixperdeg = 20;
nframes = 150;
phis = phi-[1:nframes]*deltaphi;

im = DrawGaborEdge(bkgndrgb, gaborrgb, edgergb, theta, lambda, sigma, gamma, phi, xoff, yoff, etheta, edisp, gausslim, pixperdeg); 
image(im)
axis image
set(gca,'box','off','Xtick',[],'Ytick',[])
% make a movie
clear M;
for i = 1:nframes
    im = DrawGaborEdge(bkgndrgb, gaborrgb, edgergb, theta, lambda, sigma, gamma, phis(i), xoff, yoff, etheta, edisp, gausslim, pixperdeg);
    image(im);
    axis image;
    set(gca,'box','off','Xtick',[],'Ytick',[]);
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
mpgwrite(M, gray, 'GaborMovie3.mpg', options);


%% -----------------------------------
% Starting here I'm trying to adaptively sample the LMTF surface
%%
% Section 7: Sampling a sphere (or super-sphere) using Delaunay triangulation
% and iteratively 

EXP= .5;
NITER = 40;
NOISE = 0;
% First plotting the sphere
figure; axes; hold on; axis equal;
set(gca,'View',[45 34]);
tmp = linspace(-1.1,1.1,20);
[x,y,z] = meshgrid(tmp,tmp,tmp);
f = isosurface(tmp,tmp,tmp,abs(x.^EXP)+abs(y.^EXP)+abs(z.^EXP),1);
hsphere = patch(f,'FaceAlpha',.2,'EdgeAlpha',0);

% Initial probe stimuli
nstim = 3;
x = normrnd(0,5,nstim,3);
x = x./repmat(sum(abs(x.^EXP),2).^(1/EXP),1,3);
X = [[x(:,1); -x(:,1)],[x(:,2); -x(:,2)],[x(:,3); -x(:,3)]];
for i = 1:NITER
    cla;
    hsphere = patch(f,'FaceAlpha',.2,'EdgeAlpha',0);
    K = convhull(X);
    verts = X(K',:);
    for i = 1:3:size(verts,1)
        h = patch(verts([i:i+2],1),verts([i:i+2],2),verts([i:i+2],3),[0 1 0],'FaceAlpha',.2);
    end
    
    %T = DelaunayTri(X(:,1),X(:,2),X(:,3))
    
    T = TriRep(K,X(:,1),X(:,2),X(:,3));
    P = circumcenters(T);
    %P= mean(cat(3,T.X(T.Triangulation(:,1),:),T.X(T.Triangulation(:,2),:),T.X(T.Triangulation(:,3),:)),3);
    % Above line: P is the midpoint of each triangle. Doesn't work well.
    
    plot3(X(:,1),X(:,2),X(:,3),'m*')
   % plot3(P(:,1),P(:,2),P(:,3),'m*')
    fn = faceNormals(T); % unit normal vectors
    n = neighbors(T); % These are indices into P, which are the coordinates of the face centers
    
    dotprods = zeros(size(n,1),3);
    for i = 1:size(n,1)
        v = fn(n(i,:),:); % unit vectors are on the *rows*
        u = fn(i,:);
        dotprods(i,:) = u*v';
    end
    % Need to put in a check to make sure we don't pick the same point twice!
    newpt = [];
    while (isempty(newpt))
        mindp = min(dotprods(:));
        whichpairs = softEq(dotprods,mindp);
        
        [tmprow, tmpcol] = find(whichpairs);
        %quiver3(P(tmprow,1),P(tmprow,2),P(tmprow,3),fn(tmprow,1),fn(tmprow,2),fn(tmprow,3),0.5, 'color','r');
        if (length(tmprow) ~= 4)
            disp('Should always have four faces tied for smallest dot products');
            if sum(dotprods(:)) == numel(dotprods)
                disp('all centers have already been sampled');
            end
            keyboard
        end
        % finding the pair of faces that "go together"
        % Don't want to group faces on opposite sides of the surface
        candfaces = n(whichpairs);
        tmp = ismember(candfaces, n(candfaces(1),:));
        % first argument is idxs of faces that are members of
        % small-dot-product pairs. Second argument is idxs of faces that are
        % adjacent to the first candidate face.
        if (sum(tmp) ~= 1)
            disp('Should only have one pair of faces with minimal dot product');
        end
        quiver3(P(candfaces(1),1),P(candfaces(1),2),P(candfaces(1),3),fn(candfaces(1),1),fn(candfaces(1),2),fn(candfaces(1),3),0.5, 'color','k');
        quiver3(P(candfaces(tmp),1),P(candfaces(tmp),2),P(candfaces(tmp),3),fn(candfaces(tmp),1),fn(candfaces(tmp),2),fn(candfaces(tmp),3),0.5, 'color','k');
        newpt = [P(candfaces(1),:); P(candfaces(tmp),:)];
        newpt = newpt./repmat(sum(abs(newpt.^EXP),2).^(1/EXP),1,3);
        if (NOISE)
            newpt = newpt.*unifrnd(.99,1.01,size(newpt));
        end
        % Have we already tested this point?
        if (any(all(softEq(repmat(newpt(1,:),size(X,1),1),X),2)))
            dotprods(whichpairs) = 1;
            newpt = []; % So we go to the next lowest dot product
            disp(['Going for the next dot product up: ',num2str(mindp)]);
        end
    end
    X = [X; newpt; -newpt];
    set(gca,'Xlim',[-1 1],'Ylim',[-1 1],'Zlim',[-1 1]);
    drawnow;
    %keyboard;
    % pause
end


figure; axes; hold on;
plot3(X(:,1),X(:,2),X(:,3),'k.');
axis equal


%%
% Section 8)
% Find a triangularization for which every point lies on the free boundary
% Start with the convex hull, Then find a point that lies inside the convex
% hull but close to it. Remove voronoi center that includes that point and a
% two points on the convex hull.
% Sometimes there *is* no edge in the triangulation between pairs points 
% that are neighbors on the true boundary. We need a way to redo the
% triangulation!

npts = 150;
EXP = .4;
x = normrnd(0,1,npts,2);
x = x./repmat(sum(abs(x.^EXP),2).^(1/EXP),1,size(x,2));
%x = x+unifrnd(-.1,.1,size(x))

% Prewhitening
% X is the whitened version, x is the non-whitened version (needed for
% figuring out whether an interior point projects onto a boundary edge) 
Mwht = sqrtm(inv(cov(x)));
X = x*Mwht;

T = DelaunayTri(X(:,1),X(:,2));
modtriidxs = T.Triangulation; % A set of points that I'll modify
K = T.convexHull;

figure; axes; hold on;
for i = 1:size(X,1)
    h(i) = plot(X(i,1),X(i,2),'k.');
    set(h(i),'ButtonDownFcn',['disp(',num2str(i),')']);
end
boundary = [K(1:end-1),K(2:end)];
boundaryplt = [boundary(:,1); boundary(1)];
plot(X(boundaryplt(:,1),1),X(boundaryplt(:,1),2),'g-'); % Yes, this works.
                 

while (size(boundary,1) < size(X,1))
    for Kcounter = 1:length(boundary)
        % Making a big queue of triangles to try and delete
        % As soon as we actually do delete one, we remake the queue
        % (which is inefficient)
        whichtri = vertexAttachments(T,boundary(:,1)); % Which triangles have at least one point on convhull?
        whichverts = [];
        for i = 1:length(whichtri)
            whichverts = [whichverts; T.Triangulation(whichtri{i},:)];
        end
        L1 = ismember(whichverts,boundary); % Indices to verts that lie on conv hull
        v = whichverts(sum(L1,2)==2,:); % Triangles with two vertices on boundary
        v = unique(v,'rows');
    end
    try
        if (~isempty(v))
            Vcounter = 1;
            while(Vcounter <= length(v))
                candidateidxs = v(Vcounter,:);
                edgeidxs = candidateidxs(ismember(candidateidxs,boundary));
                interioridx = candidateidxs(~ismember(candidateidxs,edgeidxs));
                isexterioredge = sum(boundary == repmat(edgeidxs,size(boundary,1),1),2) ==2 | ...
                    sum(boundary == repmat(fliplr(edgeidxs),size(boundary,1),1),2) ==2;
                if any(isexterioredge)
                    notrealedge = 0;
                else
                    notrealedge = 1; % edge between edge points isn't part of the boundary!

                end
                edgepts = x(edgeidxs,:);
                interiorpt = x(interioridx,:);
                %plot(edgepts(:,1),edgepts(:,2),'b+');
                %plot(interiorpt(:,1),interiorpt(:,2),'c+');
                % moving the interior point radially until it hits the line
                % between the egdepts
                a = interiorpt*inv(edgepts);
                interiorscale = sum(a).^-1;
                wt = a(1)./sum(a);
             
                % Sanity check, should lie on line
                % plot(interiorpt(1)*interiorscale,interiorpt(2)*interiorscale,'y*');
                if (wt < 0 | wt > 1)
                    disp('interior point not between edge points')
                elseif (interiorscale < 0)
                    disp('interior point on the wrong side of origin wrt edge points')
                elseif(notrealedge == 1)
                    disp('"edge points" lie on boundary but do not share an edge on boundary')
                else
                    tmpmodtriidxs = modtriidxs;
                    tmpmodtriidxs(all(modtriidxs == repmat(v(Vcounter,:),size(modtriidxs,1),1),2),:) = [];
                    tmpboundary = freeBoundary(TriRep(tmpmodtriidxs,X)); % changed tmpmodtriidxs
                    disp('removing a triangle');
                    % Plotting
                    cla
                    for i = 1:size(X,1)
                        h(i) = plot(X(i,1),X(i,2),'k.');
                        set(h(i),'ButtonDownFcn',['disp(',num2str(i),')']);
                    end
                   % triplot(tmpmodtriidxs,X(:,1),X(:,2),'k-');
                    boundaryplt = [tmpboundary(:,1); tmpboundary(1)];
                    plot(X(boundaryplt(:,1),1),X(boundaryplt(:,1),2),'g-');
                    drawnow;
                %   keyboard
                    
                    Vcounter = length(v)+1; % breaking out of the loop
                    modtriidxs = tmpmodtriidxs; % Debugging
                    boundary = tmpboundary; % Debugging
                    
                end
                Vcounter = Vcounter+1;
            end
        else
            disp('v is empty')
        end
    catch
        disp('hit some kind of error');
        keyboard
    end
end

% 
% 
% % Run this code if you get stuck
% orphans = find(~ismember([1:size(X,1)]',boundary(:,1)));
% Tsub = DelaunayTri(X(orphans,1),X(orphans,2));
% modtriidxs = T.Triangulation; % A set of points that I'll modify
% 
% plot(Tsub.X(:,1),Tsub.X(:,2),'b.');
% K = Tsub.convexHull;
% boundarysub = [K(1:end-1),K(2:end)];
% boundaryplt = [boundarysub(:,1); boundarysub(1)];
% plot(Tsub.X(boundaryplt(:,1),1),Tsub.X(boundaryplt(:,1),2),'g-'); % Yes, this works.
%                  
% 

%%
% Section 9: Sampling in the middle of faces that share an edge and have
% very different normal vectors. Ends up sampling in very straight lines. A
% central problem with this approach seems to be many of the triangles are
% very tall in skinny, spanning the entire vertical extent of the surface.
% Sampling doesn't break up these triangles.

% Set up
fpar =[32.9240 0.4656 14.7960 1.7664  0.0130 0.0139 100.0000  0.4890 1.6281  7.0696 0.0233  0.0077];

f1 = @(fpar,omega)fpar(1)*abs(((1i*2*pi*fpar(5).*omega+1).^-fpar(3))-fpar(2)*((1i*2*pi*fpar(6).*omega+1).^-fpar(4)));
f2 = @(fpar,omega)fpar(1+6)*abs(((1i*2*pi*fpar(5+6).*omega+1).^-fpar(3+6))-fpar(2+6)*((1i*2*pi*fpar(6+6).*omega+1).^-fpar(4+6)));


% Testing
%f1 = @(fpar,omega) 50*ones(size(omega));
%f2 = @(fpar,omega) 5*ones(size(omega));

% Plotting the fit
figure; axes; hold on;
[xx,yy,zz] = meshgrid(linspace(-.3, .3,20),...
    linspace(-.3,.3,20),...
    linspace(1, 25,20));
a = abs(f1(fpar,zz)).^-1;
b = abs(f2(fpar,zz)).^-1;
thtmp = atan2(yy,xx);
rtmp = (a.*b)./sqrt((a.*cos(thtmp)).^2+(b.*sin(thtmp)).^2); % radius of ellipse - thank you, Wikipedia

figure; axes; hold on;
V = sqrt(xx.^2+yy.^2)-rtmp;
FV = isosurface(xx,yy,zz,V,0);
h = patch(FV);
%set(h,'FaceColor','green','EdgeColor','none');
set(h,'FaceColor','green','FaceAlpha',.5','EdgeColor','none','EdgeAlpha',0);

% Sampling for random points.
% Transform x and y to lie on cyl. (x^2+y^2) = 1

% Random starting points
%npts = 4;
%tmprnd = normrnd(0,1,npts,2);
%tmprnd = tmprnd./(sqrt(sum(tmprnd.^2,2))*[1 1]); % unit vectors in LM plane
%omegas = unifrnd(1,25,npts,1);
%rotang = unifrnd(0,2*pi);
%rotmat = [cos(rotang) -sin(rotang); sin(rotang) cos(rotang)];
% Deterministic starting points

tmprnd = [
    1/sqrt(2) 1/sqrt(2);...
    1/sqrt(2) -1/sqrt(2);...
   
    1/sqrt(2) 1/sqrt(2);...
    1/sqrt(2) -1/sqrt(2)];
omegas = [1 1 25 25]';


a = abs(f1(fpar,omegas)).^-1;
b = abs(f2(fpar,omegas)).^-1;
%thtmp = atan2(tmprnd(:,2),tmprnd(:,1));
%rtmp = (a.*b)./sqrt((a.*cos(thtmp)).^2+(b.*sin(thtmp)).^2); % radius of ellipse - thank you, Wikipedia
%[x,y] = pol2cart(thtmp, rtmp);
%X = [x y omegas];  % doing it in polar coordinates
X = [b.*tmprnd(:,1) a.*tmprnd(:,2) omegas];

figure; axes; hold on;
for iter = 1:100
    % -----------------------------------------
    % Transforming all the points to lie on the unit cylinder
    % Need to tranform more? Make sure points are uniform in Y(:,1),
    % Y(:,2)?
    r = sqrt(sum(X(:,[1 2]).^2,2));
    Y=[];
    Y(:,1) = X(:,1)./r;
    Y(:,2) = X(:,2)./r;
    Y(:,3) = X(:,3);
    
 %   Mwht = sqrtm(inv(cov([Y;[-Y(:,1) -Y(:,2) Y(:,3)]])));
 %   Y = Y*Mwht;
    K = convhull([Y;[-Y(:,1) -Y(:,2) Y(:,3)]]);
    verts = [X;[-X(:,1) -X(:,2) X(:,3)]];
    
    cla
    plot3(X(:,1),X(:,2),X(:,3),'ko','MarkerFaceColor','black');
    plot3(-X(:,1),-X(:,2),X(:,3),'ko','MarkerFaceColor','black');
    h = trisurf(K,verts(:,1),verts(:,2),verts(:,3));
    set(h,'FaceColor','green','FaceAlpha',0.2);
    set(gca,'Xlim',[-.3 .3],'Ylim',[-.3 .3],'View',[15 12]);
    
    % Finding normals
    fn = faceNormals(TriRep(K,verts));
    P = mean(cat(3,verts(K(:,1),:),verts(K(:,2),:),verts(K(:,3),:)),3);
    %P = circumcenters(TriRep(K,verts));
 
    %quiver3(P(:,1),P(:,2),P(:,3),fn(:,1),fn(:,2),fn(:,3),0, 'color','r');
    n = neighbors(TriRep(K,verts)); % These are indices into P, which are the coordinates of the face centers
    % Comparing surface normals. This is not space independent. Does it matter?
    dotprods = zeros(size(n,1),3); % dot product of each face normal with normals from the adjacent three faces
    for i = 1:size(n,1)
        v = fn(n(i,:),:); % unit vectors are on the *rows*
        u = fn(i,:);
        dotprods(i,:) = u*v';
        dotprods(i,n(i,:) < i) = nan; % Nan'ing duplicates    
    end
    [sorteddps,dps_idx] = sort(dotprods(:));
    
    [row,col] = ind2sub(size(dotprods),dps_idx);
    newdirs = [];
    prevnormX = [X(:,[1 2])./repmat(sqrt(sum(X(:,[1:2]).^2,2)),1,2) X(:,3)];
    for i = 1:length(row) % loop through here once per edge
       for j = 1:2 % the two sides of the edge
            if size(newdirs,1) == 4
                % disp('got enough new directions')
                break;
            end
            if (j == 1)
                t_idx = row(i); % index into dotprods (or fn)
            else
                t_idx = n(row(i),col(i));
            end
            
            tri = K(t_idx,:);
            %hp = plot3(verts(tri,1),verts(tri,2),verts(tri,3),'ro');
            newnormdir = [P(t_idx,[1 2])./sqrt(sum(P(t_idx,[1 2]).^2)) P(t_idx,3)];
            
%             % Try to sample in the middle of the shared edge
%             sharedidxs = ismember(K(row(i),:),K(n(row(i),col(i)),:));
%             sharedidxs = K(row(i), sharedidxs);
%             newnormdir = mean(verts(sharedidxs,:));
%             newnormdir(:,1:2) = newnormdir(:,1:2)./sqrt(sum(newnormdir(:,[1 2]).^2));
%             if (j == 2)
%                 newnormdir = prevnormX(1,:);
%             end
%             % End of Try to sample in the middle of the shared edge
        
            
            if (softEq(fn(row(i),1),0) & softEq(fn(row(i),2),0))
                % normal is straight up or down
                disp('normal is straight up or down (1)');
            elseif (softEq(fn(n(row(i),col(i)),1),0) & softEq(fn(n(row(i),col(i)),2),0)) % ignoring edges that span the top or bottom
                % normal is straight up or down
                disp('normal is straight up or down (2)');          
            elseif any(abs(prevnormX(:,1:2)*newnormdir(:,1:2)') > .999 & newnormdir(:,3) == prevnormX(:,3))
                % We've already tested this direction
                disp('already tried this direction in a previous round')
            elseif ~isempty(newdirs) &&  any(abs(newdirs(:,1:2)*newnormdir(:,1:2)') > .999 & newnormdir(:,3) == newdirs(:,3))
                % We've already tested this direction
                disp('this direction is already in the queue')
            else
                newdirs = [newdirs; newnormdir];
            end
        end
    end
    % a hack to avoid sampling tons of points arbitrarily near the top and
    % bottom. Just sticking them on the top or bottom.
    newdirs(newdirs(:,3) > .9*max(X(:,3)),3) = max(X(:,3));
    newdirs(newdirs(:,3) < 1.5*min(X(:,3)),3) = min(X(:,3));
    
    % Simulating observer
    a = abs(f1(fpar,newdirs(:,3))).^-1;
    b = abs(f2(fpar,newdirs(:,3))).^-1;
    thtmp = atan2(newdirs(:,2),newdirs(:,1));
    rtmp = (a.*b)./sqrt((a.*cos(thtmp)).^2+(b.*sin(thtmp)).^2); % radius of ellipse - thank you, Wikipedia
    X = [X; [rtmp,rtmp,ones(length(a),1)].*newdirs];
   % pause;
end

%%
% Section 10
% As above, but using iterative Delaunay triangulation to find points that
% are the most surprising, given the other points.
% Critical: The only way to lower the surprise o a boundary point (TF = 1
% or 25) is to sample another boundary point. Sampling interior points will
% not help!

% Set up
fpar =[32.9240 0.4656 14.7960 1.7664  0.0130 0.0139 100.0000  0.4890 1.6281  7.0696 0.0233  0.0077];
fpar(7:end) = fpar(1:6);
f1 = @(fpar,omega)fpar(1)*abs(((1i*2*pi*fpar(5).*omega+1).^-fpar(3))-fpar(2)*((1i*2*pi*fpar(6).*omega+1).^-fpar(4)));
f2 = @(fpar,omega)fpar(1+6)*abs(((1i*2*pi*fpar(5+6).*omega+1).^-fpar(3+6))-fpar(2+6)*((1i*2*pi*fpar(6+6).*omega+1).^-fpar(4+6)));
DEBUG = 0;

% Testing
f1 = @(fpar,omega) 4*ones(size(omega));
f2 = @(fpar,omega) 40*ones(size(omega));

% Deterministic starting points

tmprnd = [
    1/sqrt(2) 1/sqrt(2);...
    1/sqrt(2) -1/sqrt(2);...
    1/sqrt(2) 1/sqrt(2);...
    1/sqrt(2) -1/sqrt(2)];
omegas = [1 1 25 25]';

tmprnd = [
    1 0;...
    0 1;...
    1/sqrt(2) 1/sqrt(2);...
    1/sqrt(2) -1/sqrt(2)];
omegas = [1 1 25 25]';



a = abs(f1(fpar,omegas)).^-1;
b = abs(f2(fpar,omegas)).^-1;
X = [b.*tmprnd(:,1) a.*tmprnd(:,2) omegas]; % Initial points

figure; axes; hold on;

for iter = 1:50
    % -----------------------------------------
    % Transforming all the points to lie on the unit cylinder
    % so that we get the real convex hull
    r = sqrt(sum(X(:,[1 2]).^2,2));
    Y=[];
    Y(:,1) = X(:,1)./r;
    Y(:,2) = X(:,2)./r;
    Y(:,3) = X(:,3);
    fullY = [Y;[-Y(:,1) -Y(:,2) Y(:,3)]];
  %  Mwht = sqrtm(inv(cov(fullY)));
    
    K = convhull(fullY);
    fullX = [X;[-X(:,1) -X(:,2) X(:,3)]];
  
    % Plotting
    plot3(fullX(:,1),fullX(:,2),fullX(:,3),'ko','MarkerFaceColor','black');
    h = trisurf(K,fullX(:,1),fullX(:,2),fullX(:,3));
    set(h,'FaceColor','green','FaceAlpha',.2);
    set(gca,'View',[65,-20]);

   % pause
    delete(h);

    % Need a way to figure out which face of the convex hull a line from
    % the origin passes through
    % Getting predictions based on n-1 points
    surprise = [];
    for i = 1:size(X,1)
        idxs = 1:size(X,1);
        idxs(i)=[];
        tmpY = [Y(idxs,:);[-Y(idxs,1) -Y(idxs,2) Y(idxs,3)]]; % tmpY is missing two points
        tmpK = convhull(tmpY); % pointing to all vertices in full X excet

        % First doing everything in *Y* space to find the triangle that a ray
        % from the Z-axis to outsidept intersects.
        
        % Need a mapping between tmpY and fullY (tmpY is missing two rows!)
        tmpK(tmpK>=i) = tmpK(tmpK>=i)+1;
        tmpK(tmpK>=i+size(X,1)) = tmpK(tmpK>=i+size(X,1))+1; % now tmpK is properly formatted for fullX
        fn = faceNormals(TriRep(tmpK,fullY)); % warnings here are OK
        outsideptY = Y(i,:);
        intersections = []; checks = [];
        % ---------------------------------
        % Finding a point on the convex hull along a ray
        % This still doesn't work. Need a check that intersection inside
        % triangle.
        % ---------------------------------
        for j = 1:size(tmpK,1) % looping over faces
            [I,check]=plane_line_intersect(fn(j,:),fullY(tmpK(j,1),:),[0 0 outsideptY(3)],outsideptY);
            intersections(j,:) = I;
            checks(j,:) = check;
        end
        intersections(checks == 0 | checks == 2,[1 2]) = Inf;
        distances_from_z_axis = sqrt(sum(intersections(:,[1 2]).^2,2));
        whichfaceidx = find(distances_from_z_axis == min(distances_from_z_axis),1); % idx into tmpK
        % Now back to X space  
        v1 = fullX(tmpK(whichfaceidx,1),:)-fullX(tmpK(whichfaceidx,2),:);
        v2 = fullX(tmpK(whichfaceidx,2),:)-fullX(tmpK(whichfaceidx,3),:);
        outsideptX = X(i,:);
        [I,check]=plane_line_intersect(cross(v1,v2),fullX(tmpK(whichfaceidx,1),:),[0 0 outsideptX(3)],outsideptX);
        if (check == 0)
            keyboard
        end
        
        % Debugging ---------------------
        if (DEBUG)
            g(1) = plot3(outsideptX(1),outsideptX(2),outsideptX(3),'r*');
            g(2) = plot3(-outsideptX(1),-outsideptX(2),outsideptX(3),'r*');
            g(3) = plot3(fullX(tmpK(whichfaceidx,1),1),fullX(tmpK(whichfaceidx,1),2),fullX(tmpK(whichfaceidx,1),3),'g*');
            g(4) = plot3(fullX(tmpK(whichfaceidx,2),1),fullX(tmpK(whichfaceidx,2),2),fullX(tmpK(whichfaceidx,2),3),'g*');
            g(5) = plot3(fullX(tmpK(whichfaceidx,3),1),fullX(tmpK(whichfaceidx,3),2),fullX(tmpK(whichfaceidx,3),3),'g*');
            g(6) = plot3(I(1),I(2),I(3),'ks');
            g(7) = plot3(-I(1),-I(2),I(3),'ks');
            h = trisurf(tmpK,fullX(:,1),fullX(:,2),fullX(:,3),'FaceAlpha',.2);
        end
        % -------------------------------
        
        pred = max(sqrt(sum(I([1 2]).^2)),eps); % avoiding Infs
        meas = sqrt(sum(fullX(i,[1 2]).^2,2));
        surprise(i) = abs(log10(meas./pred));
        
        drawnow;
        if (DEBUG)
            disp([pred meas surprise(i)]);
            %  pause
            delete(h);
            delete(g);
        end
    end
   % surprise(X(:,3) == 25) = surprise(X(:,3) == 25)/5; % downweighting boundary surprise
   % surprise(X(:,3) == 1) = surprise(X(:,3) == 1)/5;

   % surprise is associated with each vertex
    cattedsurprise = [surprise'; surprise']; % surprise for each point in fullX
    
    % Getting rid of top and bottom faces
    fn = faceNormals(TriRep(K,fullX));
    nvald = all(fn(:,[1 2]) == 0,2); % top and bottom faces
    tmpK = K(~nvald,:);

    
     % Go to a simplex that has the most surprising vertices
     %grandsurprise = mean(cattedsurprise(tmpK),2); % Mean surprise for a simplex
   %  grandsurprise = max(cattedsurprise(tmpK),[],2);
     %mostsurpriseidx = find(grandsurprise == max(grandsurprise)); % Index into tmpK (which references valid simplices)
  %   mostsurpriseidx = find(grandsurprise>= prctile(grandsurprise,0));
  %   vertsofmostsurprisingfaces = reshape(fullX(tmpK(mostsurpriseidx,:)',:),[3,size(mostsurpriseidx,1),3]);
  %   vertsofmostsurprisingfaces = permute(vertsofmostsurprisingfaces,[2 3 1]); % simplex x (LMTF) x vertices
     
     
   % 1) Picking a point based on the simplices containing the most surprising
   % vertex
    mostsurpriseidx = find(surprise == max(surprise));
    idx = unidrnd(length(mostsurpriseidx)); % pick one at random in case of tie
    mostsurpriseidx = mostsurpriseidx(idx);
    whichsimplicies = vertexAttachments(TriRep(tmpK,fullX),[mostsurpriseidx; mostsurpriseidx+size(fullX,1)/2]);
    whichsimplicies = unique([whichsimplicies{1} whichsimplicies{2}]); % index into tmpK
    whichsimplex = whichsimplicies(1); % biased towards early triangles!
    vertsofmostsurprisingface = fullX(tmpK(whichsimplex,:),:);
   individualsurprises = cattedsurprise(tmpK(whichsimplex,:));
  % parents = squeeze(surroundingverts(:,:,mostsurpriseidx));
  % vertsofmostsurprisingface = parents

    

    %2) 
    %candidates = [];
    %for i = 1:3
   % 	candidiates(i,:) = mean([X(mostsurpriseidx,:); parents(mod([i i+1],3)+1,:)]);
   % end
   % newdirs = candidiates(unidrnd(3),:);
   % alreadydonemat = softEq([X(:,[1 2])./repmat(sqrt(sum(X(:,[1 2]).^2,2)),1,2) X(:,3)], repmat(newdirs,size(X,1),1));
   % if any(all(alreadydonemat,2))
   %     keyboard
   % end
    % Done with picking a point based on the parent triangle of the most 
    % surprising vertex
    
    idx = unidrnd(length(mostsurpriseidx));
    vertsofmostsurprisingface = fullX(tmpK(mostsurpriseidx(idx),:),:);
     %individualsurprises = cattedsurprise(tmpK(mostsurpriseidx(idx),:));
  %  if (iter == 4)
  %      break
  %  end
    
  %  [sortedsurprises,sortedidx] = sort(individualsurprises);
    % if two same-TF verts have much higher surprise than third vert, just
    % go to edge and bisect line segment
   % if (sortedsurprises(3)-sortedsurprises(2) < sortedsurprises(3)/5 & ...
    %        sortedsurprises(2)-sortedsurprises(1) > sortedsurprises(2)/5 & ...
    %        vertsofmostsurprisingface(sortedidx(3),3) == vertsofmostsurprisingface(sortedidx(2),3))
  % if (vertsofmostsurprisingface(sortedidx(3),3) == vertsofmostsurprisingface(sortedidx(2),3))
    %if ((vertsofmostsurprisingface(sortedidx(2),3) == 1 |...
    %        vertsofmostsurprisingface(sortedidx(2),3) == 25) &...
    %       (vertsofmostsurprisingface(sortedidx(3),3) == 1 |...
    %       vertsofmostsurprisingface(sortedidx(3),3) == 25))
    %    individualsurprises(sortedidx) = [0 1 1];
    %    disp('edge');
        
   % else
   %     individualsurprises = [1 1 1]'; % midpoint of face
   %     individualsurprises = unifrnd(0,1,3,1); % Random point
   % end
    %newdirs = sum(repmat(10.^individualsurprises/sum(10.^individualsurprises),1,3).*vertsofmostsurprisingface);
    newdirs = sum(repmat(individualsurprises/sum(individualsurprises),1,3).*vertsofmostsurprisingface);
   %  newdirs = mean(vertsofmostsurprisingface);  % Hard to get to exterior edges this way
    % Random within selected triangle
    tmprnd = unifrnd(0,1,3,1);
    newdirs = sum(repmat(tmprnd/sum(tmprnd),1,3).*vertsofmostsurprisingface);
    
    newdirs([1 2]) = newdirs([1 2])./sqrt(sum(newdirs([1 2]).^2));
    newdirs(newdirs(:,3) > .95*max(X(:,3)),3) = max(X(:,3));
    newdirs(newdirs(:,3) < 2*min(X(:,3)),3) = min(X(:,3));
    
    newdirs
    
    % Simulating observer
    a = abs(f1(fpar,newdirs(:,3))).^-1;
    b = abs(f2(fpar,newdirs(:,3))).^-1;
    thtmp = atan2(newdirs(:,2),newdirs(:,1));
    rtmp = (a.*b)./sqrt((a.*cos(thtmp)).^2+(b.*sin(thtmp)).^2); % radius of ellipse - thank you, Wikipedia
    X = [X; [rtmp,rtmp,ones(length(a),1)].*newdirs];
end
%%
% Section 11:
% As above, but just going for for the center of the largest face on the
% convex hull. This doesn't work well. Unless we calculate areas in some
% transformed space?
% This is a nice demonstration of why we *don't* sample this way.
% 
% Set up
fpar =[32.9240 0.4656 14.7960 1.7664  0.0130 0.0139 100.0000  0.4890 1.6281  7.0696 0.0233  0.0077];

f1 = @(fpar,omega)fpar(1)*abs(((1i*2*pi*fpar(5).*omega+1).^-fpar(3))-fpar(2)*((1i*2*pi*fpar(6).*omega+1).^-fpar(4)));
f2 = @(fpar,omega)fpar(1+6)*abs(((1i*2*pi*fpar(5+6).*omega+1).^-fpar(3+6))-fpar(2+6)*((1i*2*pi*fpar(6+6).*omega+1).^-fpar(4+6)));

% Testing
%f1 = @(fpar,omega) 40*ones(size(omega));
%f2 = @(fpar,omega) 2*ones(size(omega));

% Deterministic starting points

tmprnd = [
    1/sqrt(2) 1/sqrt(2);...
    1/sqrt(2) -1/sqrt(2);...
   
    1/sqrt(2) 1/sqrt(2);...
    1/sqrt(2) -1/sqrt(2)];
omegas = [1 1 25 25]';


a = abs(f1(fpar,omegas)).^-1;
b = abs(f2(fpar,omegas)).^-1;
X = [b.*tmprnd(:,1) a.*tmprnd(:,2) omegas]; % Initial points

figure; axes; hold on;

for iter = 1:50
    % -----------------------------------------
    % Transforming all the points to lie on the unit cylinder
    % so that we get the real convex hull
    r = sqrt(sum(X(:,[1 2]).^2,2));
    Y=[];
    Y(:,1) = X(:,1)./r;
    Y(:,2) = X(:,2)./r;
    Y(:,3) = X(:,3);
    K = convhull([Y;[-Y(:,1) -Y(:,2) Y(:,3)]]);
  
    % Plotting
    fullX = [X;[-X(:,1) -X(:,2) X(:,3)]];
    plot3(fullX(:,1),fullX(:,2),fullX(:,3),'ko','MarkerFaceColor','black');
    h = trisurf(K,fullX(:,1),fullX(:,2),fullX(:,3));
    set(h,'FaceColor','green','FaceAlpha',.2);
    set(gca,'View',[65,-20]);

   % pause
    delete(h);
    areas = [];
    fn = faceNormals(TriRep(K,fullX));
    for i = 1:size(K,1)
        verts = fullX(K(i,:),:);
        a = norm(verts(1,:)-verts(2,:));
        b = norm(verts(2,:)-verts(3,:));
        c = norm(verts(1,:)-verts(3,:));
        s = (a+b+c)/2;
        areas(i) = sqrt(s*(s-a)*(s-b)*(s-c));
    end
    nvald = all(fn(:,[1 2]) == 0,2); % top and bottom faces
    areas(nvald) = 0;
    idx = find(areas == max(areas));
    idx = idx(unidrnd(length(idx))); % pick one at random in case of tie
    bigareavertices = fullX(K(idx,:),:);
    newdirs = mean(bigareavertices);
    
    newdirs([1 2]) = newdirs([1 2])./sqrt(sum(newdirs([1 2]).^2));
   % newdirs(newdirs(:,3) > .95*max(X(:,3)),3) = max(X(:,3));
   % newdirs(newdirs(:,3) < 1.1*min(X(:,3)),3) = min(X(:,3));
    
    newdirs
    
    % Simulating observer
    a = abs(f1(fpar,newdirs(:,3))).^-1;
    b = abs(f2(fpar,newdirs(:,3))).^-1;
    thtmp = atan2(newdirs(:,2),newdirs(:,1));
    rtmp = (a.*b)./sqrt((a.*cos(thtmp)).^2+(b.*sin(thtmp)).^2); % radius of ellipse - thank you, Wikipedia
    X = [X; [rtmp,rtmp,ones(length(a),1)].*newdirs];
end


%%
% Section 12:
% A new idea: sample uniformly on an iteratively updated traingularization
% of the surface. Basically the same as sampling the largest faces.

fpar =[32.9240 0.4656 14.7960 1.7664  0.0130 0.0139 100.0000  0.4890 1.6281  7.0696 0.0233  0.0077];

f1 = @(fpar,omega)fpar(1)*abs(((1i*2*pi*fpar(5).*omega+1).^-fpar(3))-fpar(2)*((1i*2*pi*fpar(6).*omega+1).^-fpar(4)));
f2 = @(fpar,omega)fpar(1+6)*abs(((1i*2*pi*fpar(5+6).*omega+1).^-fpar(3+6))-fpar(2+6)*((1i*2*pi*fpar(6+6).*omega+1).^-fpar(4+6)));


% Testing
%f1 = @(fpar,omega) 200*ones(size(omega));
%f2 = @(fpar,omega) 20*ones(size(omega));

% Deterministic starting points

tmprnd = [
    1/sqrt(2) 1/sqrt(2);...
    1/sqrt(2) -1/sqrt(2);...
   
    1/sqrt(2) 1/sqrt(2);...
    1/sqrt(2) -1/sqrt(2)];
omegas = [1 1 25 25]';


a = abs(f1(fpar,omegas)).^-1;
b = abs(f2(fpar,omegas)).^-1;
X = [b.*tmprnd(:,1) a.*tmprnd(:,2) omegas]; % Initial points

figure; axes; hold on;

for iter = 1:50
    % -----------------------------------------
    % Transforming all the points to lie on the unit cylinder
    % so that we get the real convex hull
    r = sqrt(sum(X(:,[1 2]).^2,2));
    Y=[];
    Y(:,1) = X(:,1)./r;
    Y(:,2) = X(:,2)./r;
    Y(:,3) = X(:,3);
    K = convhull([Y;[-Y(:,1) -Y(:,2) Y(:,3)]]);
  
    % Plotting
    fullX = [X;[-X(:,1) -X(:,2) X(:,3)]];
    plot3(fullX(:,1),fullX(:,2),fullX(:,3),'ko','MarkerFaceColor','black');
    h = trisurf(K,fullX(:,1),fullX(:,2),fullX(:,3));
    set(h,'FaceColor','green','FaceAlpha',.2);
    set(gca,'View',[65,-20]);

    pause
    delete(h);
    areas = [];
    fn = faceNormals(TriRep(K,fullX));
    for i = 1:size(K,1)
        verts = fullX(K(i,:),:);
        a = norm(verts(1,:)-verts(2,:));
        b = norm(verts(2,:)-verts(3,:));
        c = norm(verts(1,:)-verts(3,:));
        s = (a+b+c)/2;
        areas(i) = sqrt(s*(s-a)*(s-b)*(s-c));
    end
    nvald = all(fn(:,[1 2]) == 0,2); % top and bottom faces
    areas(nvald) = 0;
    cumulativearea = [0 cumsum(areas)];
    areathresh = unifrnd(0,sum(areas));
    whichtri = sum(cumulativearea<areathresh);
    bigareavertices = fullX(K(whichtri,:),:);
    weights = unifrnd(0,1,3,1);
    weights = weights./sum(weights);
    newdirs = weights'*bigareavertices;
    
    newdirs([1 2]) = newdirs([1 2])./sqrt(sum(newdirs([1 2]).^2));
       
    % Simulating observer
    a = abs(f1(fpar,newdirs(:,3))).^-1;
    b = abs(f2(fpar,newdirs(:,3))).^-1;
    thtmp = atan2(newdirs(:,2),newdirs(:,1));
    rtmp = (a.*b)./sqrt((a.*cos(thtmp)).^2+(b.*sin(thtmp)).^2); % radius of ellipse - thank you, Wikipedia
    X = [X; [rtmp,rtmp,ones(length(a),1)].*newdirs];
end


figure; axes; hold on;
[theta, meas] = cart2pol(X(:,1),X(:,2));
TF = X(:,3);
plot3(theta,TF,meas,'k.')

%%
% Section 13
% Uniformly on an iteratively updated traingularization
% of the surface but this time in polar coordinates (similar to above, but
% polar)

fpar =[32.9240 0.4656 14.7960 1.7664  0.0130 0.0139 100.0000  0.4890 1.6281  7.0696 0.0233  0.0077];

f1 = @(fpar,omega)fpar(1)*abs(((1i*2*pi*fpar(5).*omega+1).^-fpar(3))-fpar(2)*((1i*2*pi*fpar(6).*omega+1).^-fpar(4)));
f2 = @(fpar,omega)fpar(1+6)*abs(((1i*2*pi*fpar(5+6).*omega+1).^-fpar(3+6))-fpar(2+6)*((1i*2*pi*fpar(6+6).*omega+1).^-fpar(4+6)));


% Testing
f1 = @(fpar,omega) 10*ones(size(omega));
f2 = @(fpar,omega) 1*ones(size(omega));

% Deterministic starting points

initdirs = [
    1/sqrt(2) 1/sqrt(2);...
    1/sqrt(2) -1/sqrt(2);...
   
    1/sqrt(2) 1/sqrt(2);...
    1/sqrt(2) -1/sqrt(2)];
initdirs = [
    0 1;...
    1 0;...
   
    0 1;...
    1 0];

omegas = [1 1 25 25]';

a = abs(f1(fpar,omegas)).^-1;
b = abs(f2(fpar,omegas)).^-1;
thtmp = atan2(initdirs(:,2),initdirs(:,1));
rtmp = (a.*b)./sqrt((a.*cos(thtmp)).^2+(b.*sin(thtmp)).^2); % radius of ellipse - thank you, Wikipedia
X = [initdirs.*repmat(rtmp,1,size(initdirs,2)), omegas]; % Initial points

% Need to put this in the same range as theta
% Need to manipulate meas? Meas is not influencing the areas of triangles
% very much

figure; axes; hold on;
for iter = 1:50
    % -----------------------------------------
    % Working in polar coordinates
    [theta, meas] = cart2pol(X(:,1),X(:,2));
    TF = X(:,3);
    theta = [theta; theta+pi; theta+2*pi];
    theta(theta>2*pi) = -pi+theta(theta>2*pi)-2*pi; % everything should be [-pi:2*pi] now
    meas = [meas; meas; meas];
    TF = [X(:,3); X(:,3); X(:,3)];
    T = delaunay(theta,TF);
    fullpol = [theta TF meas];
    h(1) = trimesh(T,theta,TF,meas);
    h(2) = plot3(theta,TF,meas,'k.');
    drawnow;
    pause;
    delete(h);
    areas = [];
    for i = 1:size(T,1)
        verts = fullpol(T(i,:),:);
        verts(:,2) = (verts(:,2)-min(fullpol(:,2)))./(max(fullpol(:,2))-min(fullpol(:,2)));
        if (~softEq(max(fullpol(:,3)),min(fullpol(:,3))))
            verts(:,3) = log10(verts(:,3));
            verts(:,3) = (verts(:,3)-min(log10(fullpol(:,3))))./(max(log10(fullpol(:,3)))-min(log10(fullpol(:,3))));
            %verts(:,3) = (verts(:,3)-min(fullpol(:,3)))./(max(fullpol(:,3))-min(fullpol(:,3)));
            verts(:,3) = 2*verts(:,3).^2;
        else
            verts(:,3) = 0;
        end
        % verts = fullpol(T(i,:),:);
        % First dimension is *angular* difference
        a = verts(1,:)-verts(2,:);
        angdiff = norm([cos(verts(1))-cos(verts(2)) sin(verts(1))-sin(verts(2))]); % equalivalent to angdiff = acos([cos(verts(1)) sin(verts(1))]*[cos(verts(2)) sin(verts(2))]')
        a(1) = angdiff;

        b = verts(2,:)-verts(3,:);
        angdiff = norm([cos(verts(2))-cos(verts(3)) sin(verts(2))-sin(verts(3))]);
        b(1) = angdiff;
        
        c = verts(1,:)-verts(3,:);
        angdiff = norm([cos(verts(1))-cos(verts(3)) sin(verts(1))-sin(verts(3))]);
        c(1) = angdiff;
        
        A = norm(a); B = norm(b); C = norm(c); % lengths of sides
        s = (A+B+C)/2;
        areas(i) = sqrt(s*(s-A)*(s-B)*(s-C));
        [a;b;c]
    end
    cumulativearea = [0 cumsum(areas)];
    areathresh = unifrnd(0,sum(areas));
    whichtri = sum(cumulativearea<areathresh); % pick a random triangle
    whichtri = find(areas == max(areas),1) % pick the triangle with the greatest area
    bigareavertices = fullpol(T(whichtri,:),:);
    weights = unifrnd(0,1,3,1);
    %weights = [1 1 1]';
    weights = weights./sum(weights);
    newdirs = weights'*bigareavertices;
    newdirs(3) = newdirs(2);
    newdirs(2) = sin(newdirs(1));
    newdirs(1) = cos(newdirs(1));
    newdirs
    newdirs(newdirs(:,3) > .95*max(X(:,3)),3) = max(X(:,3));
    newdirs(newdirs(:,3) < 1.1*min(X(:,3)),3) = min(X(:,3));
    
    % Simulating observer
    a = abs(f1(fpar,newdirs(:,3))).^-1;
    b = abs(f2(fpar,newdirs(:,3))).^-1;
    thtmp = acos(newdirs(:,1));
    rtmp = (a.*b)./sqrt((a.*cos(thtmp)).^2+(b.*sin(thtmp)).^2); % radius of ellipse - thank you, Wikipedia
    X = [X; [rtmp,rtmp,ones(length(a),1)].*newdirs];
end


% Plotting polar
figure; axes; hold on;
[theta, meas] = cart2pol(X(:,1),X(:,2));
plot3(theta,X(:,3),meas,'k.','Markersize',9);
plot3(theta+pi,X(:,3),meas,'k.','Markersize',9);
plot3(theta-pi,X(:,3),meas,'k.','Markersize',9);
set(gca,'XLim',[-pi pi])

% Plotting Cartesian
figure; axes; hold on;
r = sqrt(sum(X(:,[1 2]).^2,2));
Y=[];
Y(:,1) = X(:,1)./r;
Y(:,2) = X(:,2)./r;
Y(:,3) = X(:,3);
K = convhull([Y;[-Y(:,1) -Y(:,2) Y(:,3)]]);

fullX = [X;[-X(:,1) -X(:,2) X(:,3)]];
plot3(fullX(:,1),fullX(:,2),fullX(:,3),'ko','MarkerFaceColor','black');
%h = trisurf(K,fullX(:,1),fullX(:,2),fullX(:,3));
%set(h,'FaceColor','green','FaceAlpha',.2);
set(gca,'View',[65,-20]);

%%
% Section 14: Hacking around with using Gaussian process regression to
% adaptively sample a 1-D temporal frequency threshold (or sensitivity)
% function.


fpar =[32.9240 0.4656 14.7960 1.7664  0.0130 0.0139 100.0000  0.4890 1.6281  7.0696 0.0233  0.0077];

f1 = @(fpar,omega)fpar(1)*abs(((1i*2*pi*fpar(5).*omega+1).^-fpar(3))-fpar(2)*((1i*2*pi*fpar(6).*omega+1).^-fpar(4)));
f2 = @(fpar,omega)fpar(1+6)*abs(((1i*2*pi*fpar(5+6).*omega+1).^-fpar(3+6))-fpar(2+6)*((1i*2*pi*fpar(6+6).*omega+1).^-fpar(4+6)));

omegas = [1 5 15 20 25]';  % gp.m wants everything to be a column vector
sens = f1(fpar,omegas);

meanfunc = @meanConst; hyp.mean = log10(mean(sens)); 
covfunc = @covSEiso; hyp.cov = [.5; -.4]; 
likfunc = @likGauss; hyp.lik = -3;

x_star = linspace(0,25,100)';

figure; axes; hold on;
for i = 1:30
    [m,s2,~,b] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, omegas,log10(sens),x_star);
    f = [m+2*sqrt(b); flipdim(m-2*sqrt(b),1)];
    cla;
    fill([x_star; flipdim(x_star,1)], f, [7 7 7]/8);
    plot(omegas,log10(sens),'k.');
    drawnow;
    pause;
    newstim = x_star(find (b==max(b),1));
    omegas = [omegas; newstim];
    sens = [sens; f1(fpar,newstim)];
end
figure; axes; hold on;
plot(x_star,10.^m,'k.-');
plot(omegas,sens,'ro','MarkerFaceColor','red');


%% Section 15: Playing with hyperparameters
% Let's see how well these hyperparameters work
% We're looking at the negative log likelikihood, which we want to be as
% small as possible.
covfunc = @covSEiso; 
ells = linspace(-2,2,20);
sfs = linspace(-2,2,20);
data = zeros(length(ells),length(sfs));
for i = 1:length(ells)
    for j = 1:length(sfs)
        hyp.cov = [ells(i); sfs(j)]; 
        data(j,i) = gp(hyp, @infExact, meanfunc, covfunc, likfunc, omegas,log10(sens));
        % ells changes with x (on the columns)
    end
end
figure; axes;
imagesc(data);
axis xy;
xlabel('ell'); ylabel('sf');

[i,j] = ind2sub(size(data),find(data(:) == min(data(:))))
data(i,j)
ells(i)
sfs(j)

% Now using the minimize command
% Which doesn't seem to minimize the likelihood over the error sigma
hyp2 = minimize(hyp, @gp, -100,@infExact, meanfunc, covfunc, likfunc, omegas, log10(sens))
[m,s2] = gp(hyp2, @infExact, meanfunc, covfunc, likfunc, omegas,log10(sens),x_star);
figure;
subplot(2,1,1); hold on;
f = [m+2*sqrt(s2); flipdim(m-2*sqrt(s2),1)];
fill([x_star; flipdim(x_star,1)], 10.^f, [7 7 7]/8);
plot(x_star,10.^m,'b-');
plot(omegas,sens,'k.');
ylabel('sensitivity');
set(gca,'Xlim',[x_star(1) x_star(end)]);

subplot(2,1,2); hold on;
f = [m+2*sqrt(s2); flipdim(m-2*sqrt(s2),1)];
fill([x_star; flipdim(x_star,1)], f, [7 7 7]/8);
plot(x_star,m,'b-');
plot(omegas,log10(sens),'k.');
ylabel('log sensitivity');
set(gca,'Xlim',[x_star(1) x_star(end)]);


%%
% Section 16: OK, trying to draw sampling from a 2-D GP.
% Step 1: Can I make a generative GP model that is periodic in one
% dimension and not periodic in the other dimension?

%hypper.cov = [.5;log(2);0]; 
%hypnonper.cov = [.5; -.4]; 

%covProd({@covPeriodic,@covSEiso},hyppr)

%cpr = {@covProd,{@covPeriodic,@covSEiso}}; hyppr = [hypper.cov; hypnonper.cov];
%likfunc = @likGauss; 
%hyppr.lik = 0;
sqrtn = 40;
y = [];
[x1,x2] = meshgrid(linspace(0,1,sqrtn),linspace(0,1,sqrtn));
x = [x1(:) x2(:)];
K = feval(@covSEard,[-2;-1;1],x);
if (rank(K) ~= size(K,1))
    [U,S,V] = svd(K);
    y(:,1) = U*sqrt(S)*normrnd(0,1,sqrtn.^2,1);
else
    y(:,1) = chol(K)'*normrnd(0,1,sqrtn.^2,1);
end

for i = 1:size(y,2)
    figure;
    plot3(x(:,1),x(:,2),y(:,i),'k.');
end

xlabel('x');
ylabel('y');

surf(reshape(x(:,1),sqrtn,sqrtn),reshape(x(:,2),sqrtn,sqrtn),reshape(y,sqrtn,sqrtn))

%%
% Section 17: Fitting real LMTF data with GP regression

% Loading some real data
flist = flatten(findfile(fnamesFromTxt2(fullfile(nexfilepath,'nexfilelists','Greg','LMTF','ZackLMTF.txt'))));
data = [];
for i = 1:length(flist)
    stro = notnex2stro(flist{i});
    Lstim_idx = strcmp(stro.sum.trialFields(1,:), 'stim_idx');
    Llcc = strcmp(stro.sum.trialFields(1,:), 'lcc');
    Lmcc = strcmp(stro.sum.trialFields(1,:), 'mcc');
    Ltf = strcmp(stro.sum.trialFields(1,:), 'tf');
    Loog = strcmp(stro.sum.trialFields(1,:), 'oog');
    
    [stim_idxs,init_stim_trial_idxs] = unique(stro.trial(:,Lstim_idx),'last');
    questmodes = stro.trial(init_stim_trial_idxs, Llcc|Lmcc);
    tfs = stro.trial(init_stim_trial_idxs,Ltf);
    
    % Out of gamut checking (data)
    funds = reshape(stro.sum.exptParams.fundamentals,length(stro.sum.exptParams.fundamentals)/3,3);
    spds = SplineSpd([380:2:780]',reshape(stro.sum.exptParams.mon_spd,length(stro.sum.exptParams.mon_spd)/3,3),[380:5:780]');
    M = funds'*spds;
    bkgndrgb = stro.sum.exptParams.bkgndrgb;
    [in_gamut,scalar] = gamutCheck([questmodes zeros(size(questmodes,1),1)], bkgndrgb, M, 'both');
    questmodes(~in_gamut) = questmodes(~in_gamut).*scalar(~in_gamut);
    data = [data; questmodes tfs ~in_gamut']; % Lcc Mcc TF OOG
end
Loog = logical(data(:,end)); % reusing this variable name
% Below, we're fitting to the OOG points as well as
% the ~OOG points. The rationale is that if we don't use these points 
% the algorithm won't realize that a particular set of stimuli are out of
% and will continue to select stimuli that go OOG. This will have the
% effect of introducing a sharp edge in the fitted data though, which 
% isn't great.
[th,r] = cart2pol(data(:,1),data(:,2)); 
th = mod(th,pi); 

% Setting up a grid of points on which to evaluate the posterior
sqrtn = 40;
[x1,x2] = meshgrid(linspace(1,25,sqrtn),linspace(min(th),max(th),sqrtn));
x_star = [x1(:) x2(:)]; % column order: TF, theta

% Figuring out how large the cone contrast can be, in each theta, 
% before it goes out of gamut.
% If a stimulus is predicted to go out of gamut, don't test it.
funds = reshape(stro.sum.exptParams.fundamentals,size(stro.sum.exptParams.fundamentals,1)/3,3)';
spds= reshape(stro.sum.exptParams.mon_spd,size(stro.sum.exptParams.mon_spd,1)/3,3);
spds = SplineSpd([380:2:780]',spds,[380:5:780]');
M = funds*spds;
bkgndrgb = stro.sum.exptParams.bkgndrgb;
cc = [cos(x2(:,1)), sin(x2(:,1)) zeros(size(x2,1),1)];
[in_gamut,scalar] = gamutCheck(cc, bkgndrgb, M, 'both');
gamutmask = repmat(log10(scalar'),1,size(x2,2));
% end of gamut stuff

% Setting up mean, covariance, and likelihood functions + hyperparameters
covfunc = {'covProd',  {{'covMask',{[1 0],'covSEiso'}}, {'covMask',{[0 1],'covPeriodPi'}}}};
meanfunc = @meanConst;
likfunc = @likGauss;

hyp.cov = [3;0;1;-1];
hyp.mean = .05;
hyp.lik = -1;

% The hyperparameters are pretty reasonable, but here we fine tuned them,
% maximum likelihood style. 
% Not using the OOG points to get estimate smoothness, error.
hyp2 = minimize(hyp, @gp, -100,@infExact, meanfunc, covfunc, likfunc, [data(~Loog,3),th(~Loog)],log10(r(~Loog)));
% Now fitting the data (inefficient, since a lot of the computation here is
% shared with the line above).
[m,s2,pm,ps2] = gp(hyp2, @infExact, meanfunc, covfunc, likfunc, [data(:,3) th],log10(r),x_star);

figure; axes; hold on;
surf(reshape(x1,sqrtn,sqrtn),reshape(x2,sqrtn,sqrtn),reshape(pm,sqrtn,sqrtn))
plot3(data(:,3),th,log10(r),'ko','MarkerSize',10,'MarkerFaceColor','black');

h(1) = surf(reshape(x1,sqrtn,sqrtn),reshape(x2,sqrtn,sqrtn),reshape(pm+2*sqrt(ps2),sqrtn,sqrtn));
h(2) = surf(reshape(x1,sqrtn,sqrtn),reshape(x2,sqrtn,sqrtn),reshape(pm-2*sqrt(ps2),sqrtn,sqrtn));
set(h,'FaceColor',[.7 .7 .7],'Facealpha',.2,'EdgeAlpha',.2)
xlabel('TF'); ylabel('Theta'); zlabel('log_1_0 Threshold');
set(gca,'Zlim',[min(pm(:)+2*sqrt(ps2(:))) max(pm(:)+2*sqrt(ps2(:)))]);

% Adaptive sampling
% Taking 'n' new color directions (not just taking the max)
% Algorithm: Take peak, go down sorted list of s2s until you find one
% that isn't continguous with the first, etc.
ps2(m>gamutmask(:)) = 0; % No even considering stimulus directions that are predicted to leave the gamut. 
figure; axes; hold on;
surf(reshape(x1,sqrtn,sqrtn),reshape(x2,sqrtn,sqrtn),reshape(ps2,sqrtn,sqrtn));
nnewsstim = 4;
mask = zeros(sqrtn,sqrtn);
newstim = []; idx = 1;
[sorteds2,I] = sort(ps2,1,'descend');
while (length(newstim) < nnewsstim & idx < length(ps2))
    [i,j] = ind2sub(size(mask),I(idx));
    neighborhoodcols = [j-1 j j+1];
    neighborhoodrows = [i-1 i i+1];
    neighborhoodcols(neighborhoodcols < 1) = [];
    neighborhoodcols(neighborhoodcols > sqrtn) = [];
    neighborhoodrows(neighborhoodrows < 1) = sqrtn+neighborhoodrows(neighborhoodrows < 1);
    neighborhoodrows(neighborhoodrows > sqrtn) = neighborhoodrows(neighborhoodrows > sqrtn)-sqrtn;
    if (~any(any(mask(neighborhoodrows,neighborhoodcols))))
        newstim = [newstim; x1(i,j) x2(i,j)];
    end
    mask(i,j) = 1;
    idx = idx+1;
end
plot3(newstim(:,1),newstim(:,2),repmat(max(ps2),nnewsstim,1),'ks','MarkerFaceColor','black');

% Plotting the data in LMTF space
figure; axes; hold on;
plot3(data(~Loog,1),data(~Loog,2),data(~Loog,3),'k.');
plot3(-data(~Loog,1),-data(~Loog,2),data(~Loog,3),'k.');
plot3(data(Loog,1),data(Loog,2),data(Loog,3),'r.');
plot3(-data(Loog,1),-data(Loog,2),data(Loog,3),'r.');

% Plotting the final regression fit in LMTF space
[xx,yy,zz] = meshgrid(linspace(-max(abs(data(:,1))),max(abs(data(:,1))),20),...
    linspace(-max(abs(data(:,2))),max(abs(data(:,2))),20),...
    linspace(min(data(:,3)),max(data(:,3)),20));
[th_star,~] = cart2pol(xx,yy);
[th,r] = cart2pol(data(:,1),data(:,2));
m = gp(hyp, @infExact, meanfunc, covfunc, likfunc, [data(:,3) th],log10(r),[zz(:) th_star(:)]);
m = reshape(m,size(xx));
V = sqrt(xx.^2+yy.^2)-10.^m;
FV = isosurface(xx,yy,zz,V,0);
h = patch(FV);
set(h,'FaceColor','blue','FaceAlpha',.5','EdgeColor','none','EdgeAlpha',0);

% Plotting the new stimulus directions
for i = 1:size(newstim,1)
    plot3([-1 1].*max(r).*cos(newstim(i,2)),[-1 1].*max(r).*sin(newstim(i,2)),[1 1].*newstim(i,1),'k-')
end

%%
% Section 18
% Simulating an LMTF experiment using GP regression to adaptively pick
% stimuli to test

% Setting up tthe twisty funnel model (Pretty reasonable parameters - not
% sure where they came from though)
fpar =[32.9240 0.4656 14.7960 1.7664  0.0130 0.0139 100.0000  0.4890 1.6281  7.0696 0.0233  0.0077];
% Apollo at eccentricity = (5,0)
fpar = [19.1948 0.9998 19.9997 1.2534 0.0039 0.0124 55.5097 0.2051 6.3821 19.9731 0.0170 0.0022];

f1 = @(fpar,omega)fpar(1)*abs(((1i*2*pi*fpar(5).*omega+1).^-fpar(3))-fpar(2)*((1i*2*pi*fpar(6).*omega+1).^-fpar(4)));
f2 = @(fpar,omega)fpar(1+6)*abs(((1i*2*pi*fpar(5+6).*omega+1).^-fpar(3+6))-fpar(2+6)*((1i*2*pi*fpar(6+6).*omega+1).^-fpar(4+6)));
sigma = .1; % optional noise on observer's thresholds

% Testing
%f1 = @(fpar,omega) 2*ones(size(omega));
%f2 = @(fpar,omega) 30*ones(size(omega));

% Deterministic starting points
initdirs = [
    1/sqrt(2) 1/sqrt(2);...
    1/sqrt(2) -1/sqrt(2);... 
    1/sqrt(2) 1/sqrt(2);...
    1/sqrt(2) -1/sqrt(2)];
omegas = [1 1 25 25]';
a = abs(f1(fpar,omegas)).^-1;
b = abs(f2(fpar,omegas)).^-1;
thtmp = mod(atan2(initdirs(:,2),initdirs(:,1)),pi);
rtmp = (a.*b)./sqrt((a.*cos(thtmp)).^2+(b.*sin(thtmp)).^2); % radius of ellipse - thank you, Wikipedia
X = [initdirs.*repmat(rtmp,1,size(initdirs,2)), log10(omegas)]; % Initial points
% Note: X = [Lcc, Mcc, log10(TF)]; 

% Setting up hyperparameters for the GP
meanfunc = @meanConst; hyp.mean = -2;
meanfunc = {@meanSum, {@meanConst, @meanLinear}};
likfunc = @likGauss
covfunc = {'covProd',  {{'covMask',{[1 0],{'covMaterniso',1}}}, {'covMask',{[0 1],'covPeriodPi'}}}};
hyp.cov = [2;-.4;0;-1];
hyp.lik = -1;
hyp = struct('cov', [-0.75 .5 1 -1.5]', 'lik', -2.3, 'mean', [-1; .6; 0]);
hyp = struct('cov', [0 .2 .2 -1.8]', 'lik', -2.3, 'mean', [-1; .6; 0]);

% Setting up grid of points on which to evaluate the posterior 
sqrtn = 50;
[x1,x2] = meshgrid(logspace(log10(1),log10(40),sqrtn),linspace(0,pi,sqrtn));
x_star = [log10(x1(:)) x2(:)];

% Figuring out how great the cone contrast has to be in each theta to go out of the gamut
% If a stimulus is predicted to go out of gamut, don't test it.
load('Dell4BitsCal');
cal = cals{end};
load('T_cones_smj10');
funds = T_cones_smj10;
spds = SplineSpd([380:4:780]',cals{end}.P_device,[380:5:780]');
M = funds*spds;
bkgndRGB = round(cals{end}.bgColor*255)+1;
bkgndrgb = diag(cals{end}.gammaTable(bkgndRGB,:));
cc = [cos(x2(:,1)), sin(x2(:,1)) zeros(size(x2,1),1)];
[in_gamut,scalar] = gamutCheck(cc, bkgndrgb, M, 'both');
gamutmask = repmat(log10(scalar'),1,size(x2,2));
% end of gamut stuff

niter = 10;
figure('Position',[440 76 513 722]);
% X is [Lcc, Mcc, log10(TF)]
for iter = 1:niter
    [th,r] = cart2pol(X(:,1),X(:,2));
    th = mod(th,pi);
    [m,s2,pm,ps2] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, [X(:,3) th],log10(r),x_star);
    
    subplot(2,1,1);hold on;
    h(1) = surf(reshape(log10(x1),sqrtn,sqrtn),reshape(x2,sqrtn,sqrtn),reshape(pm,sqrtn,sqrtn));
    h(2) = plot3(X(:,3),th,log10(r),'ko','MarkerSize',10,'MarkerFaceColor','black');
    set(gca,'View',[-37.5 30]);
     
    % Where should we sample next?
    ps2(m>gamutmask(:)) = 0;
    subplot(2,1,2); 
    h(3) = surf(reshape(x1,sqrtn,sqrtn),reshape(x2,sqrtn,sqrtn),reshape(ps2,sqrtn,sqrtn));
    set(gca,'View',[-37.5 30]);
    drawnow;
   % pause
    if (iter ~= niter)
        delete(h);
    end
    
    % Picking the new stimulus direction to test
    TFnTheta = x_star(find(ps2 == max(ps2),1),:);
    newdirs = []; newdirs(1) = cos(TFnTheta(2));
    % Simulating observer
    a = abs(f1(fpar,10^TFnTheta(1))).^-1;
    b = abs(f2(fpar,10^TFnTheta(1))).^-1;
    rtmp = (a.*b)./sqrt((a.*cos(TFnTheta(2))).^2+(b.*sin(TFnTheta(2))).^2); % radius of ellipse - thank you, Wikipedia
    if (~isempty(sigma))
        rtmp = 10.^(log10(rtmp)+normrnd(0,sigma));
    end
    Xupdate = [rtmp*cos(TFnTheta(2)) rtmp*sin(TFnTheta(2)) TFnTheta(1)];
    [flag, scalar] = gamutCheck([Xupdate(1:2) 0],bkgndrgb,M,'both');
    if (~flag)
        Xupdate(1:2) = Xupdate(1:2).*scalar;
    end
    X = [X; Xupdate];
    
    if ~rem(iter,10)
       disp('Updating hyperparamters');
       [th,r] = cart2pol(X(:,1),X(:,2));
       oldhyp = hyp;
       hyp = minimize(hyp, @gp, -100,@infExact, meanfunc, covfunc, likfunc, [X(:,3),th],log10(r));
       disp([oldhyp.cov-hyp.cov]);
       if (hyp.lik < -5) % There's no way the thresholds are this accurate!
           hyp.lik = -5; 
       end
   end
end

figure; axes; hold on;
[th,r] = cart2pol(X(:,1),X(:,2));
th = mod(th,pi);
plot3(th,X(:,3),r,'k.');

figure; axes; hold on;
plot3(X(:,1),X(:,2),10.^X(:,3),'k.');
plot3(-X(:,1),-X(:,2),10.^X(:,3),'k.');

% Plotting the final regression fit
[xx,yy,zz] = meshgrid(linspace(-max(abs(X(:,1))),max(abs(X(:,1))),20),...
    linspace(-max(abs(X(:,2))),max(abs(X(:,2))),20),...
    linspace(min(X(:,3)),max(X(:,3)),20));
[th_star,~] = cart2pol(xx,yy);
[th,r] = cart2pol(X(:,1),X(:,2));
m = gp(hyp, @infExact, meanfunc, covfunc, likfunc, [X(:,3) th],log10(r),[zz(:) th_star(:)]);
m = reshape(m,size(xx));
V = sqrt(xx.^2+yy.^2)-10.^m;
FV = isosurface(xx,yy,10.^zz,V,0); % undoing the log10 on TF for plotting only
h = patch(FV);
%set(h,'FaceColor','green','EdgeColor','none');
set(h,'FaceColor','blue','FaceAlpha',.5','EdgeColor','none','EdgeAlpha',0);

% Plotting the real surface
a = abs(f1(fpar,10.^zz)).^-1; % Simulated observer wants TF, not log10(TF)
b = abs(f2(fpar,10.^zz)).^-1;
thtmp = atan2(yy,xx);
rtmp = (a.*b)./sqrt((a.*cos(thtmp)).^2+(b.*sin(thtmp)).^2); % radius of ellipse - thank you, Wikipedia

V = sqrt(xx.^2+yy.^2)-rtmp;
FV = isosurface(xx,yy,10.^zz,V,0);
h = patch(FV);
%set(h,'FaceColor','green','EdgeColor','none');
set(h,'FaceColor','green','FaceAlpha',.5','EdgeColor','none','EdgeAlpha',0);

%%
% Section 19
% Ading random phase to a sinusoid. Another possibility for Aim 3 if
% microstim doesn't work out for some reason.

Fs = 1000;  % samples/sec
nsec = 1;
t = [1./Fs:1./Fs:nsec]'; % sec
f = 4; % Hz
intendedSD = pi/20;
cutfreq = 10; % Specifying the autocorrelation of the phase random process (Hz)

a = round(cutfreq*nsec); % cutfreq/Fs*Fs*nsec
PSD_int = a/(Fs*nsec).^2/2;
PSD = [(intendedSD.^2/PSD_int)*ones(1,a)'; zeros(1,Fs*nsec-a)'];
PSD = -log(rand(size(PSD))).*PSD; % generate/multiply by negative exponential variates.
A = sqrt(PSD);
A = A.*exp(j*2*pi*rand(size(A)));
phasenoise = real(ifft(A));

[intendedSD std(phasenoise)] % sanity check

figure;
subplot(2,1,1);
plot(t,phasenoise);
ylabel('phase (rad)');
subplot(2,1,2);
y = sin(2*pi*f*t+phasenoise);
plot(t,y,'.-');
ylabel('Amplitude');
xlabel('time (s)');

%%
% Section 20
% Adding simulated parvocellular responses (with random phases) to see
% at what temporal frequency color information is lost.
% Step 1: L-M and M-L no noise
% Step 2: L-M and M-L with Poisson noise
% Step 3: A more realistic parvocellular model
% Never finished this

% Step 1, playing around with addng sinusoids
% number of cells doesn't matter much. With a 10 ms delay 
% there's a subtle drop in signal.
deltaT = .0001;
envelope = [linspace(0,1,round(.667/4/deltaT)) ones(1,round(.667/2/deltaT)) linspace(1,0,round(.667/4/deltaT))]; % in deltaTs
ntimesteps = length(envelope);
tfs = linspace(1,50,10);
ncells = 10;

data = [];
for i = 1:length(tfs)
    tmpdata = zeros(ncells,length(envelope));
    for j = 1:ncells
        delay = unifrnd(0,10)/1000; % in s
        phi = 2*pi*(delay*tfs(i)); % in rad
        sinusoid = sin(2*pi*tfs(i)*linspace(0,ntimesteps*deltaT,ntimesteps)+phi);
        temporalframes = envelope.*sinusoid;
        tmpdata(j,:) = temporalframes;
    end
    data(i,:) = sum(tmpdata);
end
%plot(data')

contrastsens = sqrt(sum(data.^2,2));
figure; plot(tfs,contrastsens,'k.');
xlabel('TF'); ylabel('contrast sensitivity');

%%
% Step 1.5
% Swanson model
Wl = .5;
Wm = 1-Wl;
Rl = 1000;
Rm = 1000;
d = .01; % delay in s
omega = 25;
t = linspace(0,1,1000);
rho = omega*d;
LminusM = Wl*Rl*sin(omega*t)-Wm*Rm*sin(omega*t-rho);
amp = sqrt((Wl*Rl-Wm*Rm*cos(rho)).^2+(Wm*Rm*sin(rho)).^2)

figure; 
subplot(2,1,1)
plot(LminusM)
title(amp)


%%
% Step 2 LGN neurons with Poisson noise brute force simulation
% Parameters I'll need: response amplitude per unit contrast for center and
% surround (as a function of TF). Delay between center and surround.

omega = 30; % Hz
mean = 20; % sp/sec
amp = 10; % sp/sec
CSgainratio = .5;
delay = .01; % s
x = linspace(0,1,1000);
center_resp = amp*sin(2*pi*x*omega);
surround_resp = CSgainratio*amp*sin(2*pi*omega*(x-delay));
% 
% figure; axes; hold on;
% plot(x,center_resp,'b-');
% plot(x,surround_resp,'g-');

% Just taking into account the delay between center and surround
% what does the TF tuning curve look like?
data =[];
tfs = [1:100];
for i = tfs % tf
    center_resp = amp*sin(2*pi*x*i);
    surround_resp = CSgainratio*amp*sin(2*pi*i*(x-delay));
    resp = center_resp-surround_resp;
    data = [data; max(resp)]
end
figure;
plot(tfs,data)


%%
% Step 3
% Parameters for P-cell center responses (LGN)
% From Benardete and Kaplan 1997
rmpath '/Users/greghorwitz/Desktop/MatlabCode/-General-/gpml/util' % Shadows "unwrap" which I need
i = sqrt(-1); % Important!
% centers (expressing all times in seconds)
center.A.mean = 34.93; % spikes/sec/unit contrast
center.A.median = 27.41; % spikes/sec/unit contrast
center.A.sd = 21;
center.Nl.mean = 27.13;
center.Nl.sd = 7.67;
center.Hs.mean = .72;
center.Hs.sd = .13;
center.Ts.mean = 46.19/1000; % s
center.Ts.sd = 24.17/1000; % s
center.D.mean = 3.5/1000; % s (Based on Table 5)
center.D.sd = 2.39/1000; % s
center.Tl.mean = 52.48/center.Nl.mean/1000;  % parameterized this way because it's time to peak of STA
center.Tl.sd = 4.28/center.Nl.sd/1000;

% surrounds
surround.A.mean = 37.41; % spikes/sec/unit contrast
surround.A.median = 17.07; % spikes/sec/unit contrast
surround.A.sd = 54.41;
surround.Nl.mean = 61.88;
surround.Nl.sd = 21.84;
surround.Hs.mean = .26;
surround.Hs.sd = .57;
surround.Ts.mean = 57.31/1000; % s
surround.Ts.sd = 51.67/1000; % s
surround.D.mean = (8.41+1.16)/1000; % s (Based on Table 3 and Table 6)
surround.D.mean = 3.5/1000; % s (Same as center)
surround.D.sd = 1.48/1000; % s
surround.Tl.mean = 59.97/surround.Nl.mean/1000;  % parameterized this way because it's time to peak of STA
surround.Tl.sd = 5.46/surround.Nl.sd/1000;

TF = logspace(log10(0.25),log10(100),101); % Hz
omega = TF*2*pi; % rad/s
% Computing center stuff
freqresponse = center.A.mean*exp(-i.*omega.*center.D.mean).*...
    (1-(center.Hs.mean./(1+i.*omega.*center.Ts.mean))).*...
    (1./(1+i.*omega.*center.Tl.mean)).^center.Nl.mean;
figure;
subplot(2,2,1); hold on;
loglog(TF,abs(freqresponse),'.'); set(gca,'Xlim',[TF(1) TF(end)]);
set(gca,'Xscale','log','Yscale','log')
subplot(2,2,3); hold on;
semilogx(TF,unwrap(angle(freqresponse))/pi,'.'); set(gca,'Xlim',[TF(1) TF(end)],'Xscale','log');

% Computing surround stuff
freqresponse = surround.A.mean*exp(-i.*omega.*surround.D.mean).*...
    (1-(surround.Hs.mean./(1+i.*omega.*surround.Ts.mean))).*...
    (1./(1+i.*omega.*surround.Tl.mean)).^surround.Nl.mean;
subplot(2,2,1);
loglog(TF,abs(freqresponse),'g.'); set(gca,'Xlim',[TF(1) TF(end)]);
subplot(2,2,3);
semilogx(TF,unwrap(angle(freqresponse))/pi,'g.'); set(gca,'Xlim',[TF(1) TF(end)]);

% Now converting back to the time domain

TF = [0:1:20000]; 
n = length(TF);
omega = TF*2*pi; % rad/s
freqresponse_c = center.A.median*exp(-i.*omega.*center.D.mean).*...
    (1-(center.Hs.mean./(1+i.*omega.*center.Ts.mean))).*...
    (1./(1+i.*omega.*center.Tl.mean)).^center.Nl.mean;
freqresponse_s = surround.A.median*exp(-i.*omega.*surround.D.mean).*...
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

% Now I have center and surround filters that will take a time varying
% contrast stimulus and give a response in sp/sec. 
% Still a problem with normalizing!
% We want the integral of the kernel to be something reasonable
% Benardete's deltaT was 1/135
centerkernel = (centerkernel./max(centerkernel))*50; % to put it in sp/sec????
% 50 to make peak response to 0.5 contrast grating = 25 sp/sec 
surroundkernel = surroundkernel*(sum(centerkernel)./sum(surroundkernel))*(1/1.5)
% center gain is 1.5 times surround gain

% From Figure 9: Center gain should be 1.5x surround gain
% p. 183 surround lags center by 6-9 ms (implies a null at > 100 Hz!)
figure;
subplot(2,2,1);
plot(x,centerkernel,'.-'); set(gca,'Xlim',[0 .5]);
subplot(2,2,2);
plot(x,surroundkernel,'.-'); set(gca,'Xlim',[0 .5]);
subplot(2,2,3); hold on;
plot(x, centerkernel+surroundkernel,'r-','LineWidth',2);
plot(x, centerkernel-surroundkernel,'k-','LineWidth',2); set(gca,'Xlim',[0 .5]);
deltaT = x(2)-x(1);

% -----------------------------------
% Brute force time domain simulations
% -----------------------------------

% First, temporal contrast sensitivity assuming one color dimension
% Assuming center response is 1.5x surround response
envelope = [zeros(1,round(.2/deltaT)),...
    linspace(0,1,round(.667/4/deltaT)),...
    ones(1,round(.667/2/deltaT)),...
    linspace(1,0,round(.667/4/deltaT)),...
    zeros(1,round(.2/deltaT))]; % in deltaTs
%envelope = ones(size(envelope));
ntimesteps = length(envelope);
tfs = logspace(log10(1),log10(100),50);
maxcontrast =.01;

RGdata = [];
LUMdata = [];
for i = 1:length(tfs)
    sinusoid = maxcontrast*sin(2*pi*tfs(i)*linspace(0,ntimesteps*deltaT,ntimesteps));
    temporalframes = envelope.*sinusoid;
    c = conv(temporalframes, centerkernel,'valid')*deltaT;
    s = conv(temporalframes, surroundkernel,'valid')*deltaT;
    LUMdata = [LUMdata; c-s];
    RGdata = [RGdata; c+s];
end

figure; 
subplot(2,2,1);
hold on;
imagesc(RGdata);
subplot(2,2,2);
hold on;
plot(tfs,max(abs(RGdata')),'r-','LineWidth',2);
plot(tfs,max(abs(LUMdata')),'k-','LineWidth',2);

%plot(tfs,sqrt(sum(RGdata.^2,2)),'r.-');
%plot(tfs,sqrt(sum(LUMdata.^2,2)),'k.-');
set(gca,'Xlim',[tfs(1) tfs(end)]);
set(gca,'Yscale','log','Xscale','log');


%%
% Putting down some spikes
% contrast = 0.05, tf = 3 or contrast = .5, TF = 30;
% contrast =.2;
% tf = 30;
% baselinerate = 5; % sp/s
% ISOLUM = 1;
% ntrials = 10; % aka number of cells
% 
% sinusoid = contrast*sin(2*pi*tf*linspace(0,ntimesteps*deltaT,ntimesteps));
% temporalframes = envelope.*sinusoid;
% %temporalframes = contrast*ones(1,length(envelope));
% c = conv(temporalframes, centerkernel)*deltaT;
% s = conv(temporalframes, surroundkernel)*deltaT;
% figure; subplot(3,1,1); hold on;
% t = [0:deltaT:deltaT*(length(c)-1)]
% plot(t, c-s); set(gca,'Xlim',[min(t) max(t)]);
% 
% subplot(3,1,2); hold on;
% PSTH = zeros(ntrials,length(c));
% for i = 1:ntrials
%     randnums = unifrnd(0,1,1,length(c));
%     if (ISOLUM == 0)
%         spikes = randnums < (c-s)+(baselinerate*deltaT);
%     else
%         spikes = randnums < (c+s)+(baselinerate*deltaT);
%     end
%     PSTH(i,:) = spikes;
%     plot(t,spikes+(i-1));
%     title(sum(spikes));
% end
% subplot(3,1,3)
% bar(t,sum(PSTH)); set(gca,'Xlim',[min(t) max(t)]);
% 
% % Now simulating effects on a V1 neuron
% Vm = zeros(1,length(c));
% depression = zeros(ntrials,length(c));
% alpha = .995; % time constant of membrane (.99 = 5 msec (100 deltaTs). .995 = 10 ms, .999 = 50 ms, .9995 = 100 ms)
% beta = .1; % < - doesn't do anything
% gamma = 1; %.7
% delta = .999; % .999 time constant for synaptic depression
% for i = 1:length(c)
%    if i>1
%        Vm(i) = sum(exp(-depression(:,i-1)).*beta.*PSTH(:,i))+Vm(i-1)*alpha;
%        %  Vm(i) = (beta*PSTH(i)+Vm(i-1)*alpha)*(1+depression(i-1))^-1;
%        for j = 1:ntrials
%            depression(j,i) = gamma*PSTH(j,i)+depression(j,i-1)*delta;
%        end
%    end
% end
% figure; 
% subplot(2,1,1); plot(t,Vm);
% subplot(2,1,2); plot(t,depression,'k-');

%%
% Trying to implement the Abbott and Chance model
% Ignoring "s" (for now) which has a time constant of 20 seconds.

% trying to increase temporal integration with Tm
% ISOLUM = 1;
% deltaT = 5.0000e-05;  % s
% te = .002; % s
% ti = .010; % s
% td= 0.001; % s  %should be 0.3 (.5 works well)
% tm = .7; % Membrane time constant in ms .03
% V0 = -70; % mV
% Ve = 0;
% Vi = -90;
% Vthresh = -55; % mV
% Vreset = -58;
% nsynexc = 30;
% nsyninh = 10;
% Dexc = zeros(1,nsynexc);
% Dinh = zeros(1,nsyninh);
% d = 0; % 1 mean nothing happens
% gexc = 1.7*ones(1,nsynexc); % 1.5
% ginh = 1.4*ones(1,nsyninh);
% contrast = 0.05;


% Trying to increase temporal integration with syanptic depression
ISOLUM = 0;
deltaT = 5.0000e-05;  % s
te = .002; % s
ti = .010; % s
td= .3; % s  %should be 0.3 (.5 works well)
tm = .03; % Membrane time constant in ms .03
V0 = -70; % mV
Ve = 0;
Vi = -90;
Vthresh = -55; % mV
Vreset = -58;
nsynexc = 8;
nsyninh = 2;
Dexc = zeros(1,nsynexc);
Dinh = zeros(1,nsyninh);
d = 0; % 1 mean nothing happens
gexc = 2.5*ones(1,nsynexc); % 1.5
ginh = .2*ones(1,nsyninh);
contrast = 10;  % .4 for 30 Hz

dur = 5; % s
ntimesteps = dur./deltaT;
t = [0:deltaT:deltaT*(ntimesteps-1)];
V = nan*ones(1,ntimesteps);
spikes = nan*ones(1,ntimesteps);
Ge =  0;
Gi =  0;
V(1) = V0;
baselinerate = 20;
tf = 20;  % Hz

sinusoid = contrast*cos(2*pi*tf*linspace(0,ntimesteps*deltaT,ntimesteps));
envelope = [linspace(0,1,round(.667/4/deltaT)),...
    ones(1,round(.667/2/deltaT)),...
    linspace(1,0,round(.667/4/deltaT))]; % in deltaTs

envelope = [zeros(1,floor((ntimesteps-length(envelope)))/2),...
    envelope,...
    zeros(1,ceil((ntimesteps-length(envelope)))/2)];

temporalframes = envelope.*sinusoid;
%temporalframes = contrast*ones(1,length(envelope));
c = conv(temporalframes, centerkernel*135*deltaT,'same');
s = conv(temporalframes, surroundkernel*135*deltaT,'same');

% Setting up afferent spikes
if (ISOLUM == 0)
    kernel = c-s;
else
    kernel = c+s;
end

exc_t = [];
for i = 1:nsynexc
    tmp = unifrnd(0,1,1,ntimesteps) < (baselinerate+kernel)*deltaT;
    exc_t{i} = find(tmp);
end
exc = zeros(nsynexc,max(cellfun('length',exc_t)));
for i = 1:nsynexc
    exc(i,1:length(exc_t{i})) = exc_t{i};
end

inh_t = [];
for i = 1:nsyninh
    tmp = unifrnd(0,1,1,ntimesteps) < (baselinerate+kernel)*deltaT;
    inh_t{i} = find(tmp);
end
inh = zeros(nsyninh,max(cellfun('length',inh_t)));
for i = 1:nsyninh
    inh(i,1:length(inh_t{i})) = inh_t{i};
end


for i = 1:ntimesteps
    if (i > 1)
        V(i) = V(i-1)+(V0-V(i-1)+Ge*(Ve-V(i-1))+Gi*(Vi-V(i-1)))*(deltaT/tm);
        Ge = Ge-Ge*(deltaT/te);
        Gi = Gi-Gi*(deltaT/ti);
        Dinh = Dinh+(1-Dinh).*(deltaT/td);
        Dexc = Dexc+(1-Dexc).*(deltaT/td);
    end
    if any(exc(:) == i)
        for j = 1:nsynexc
            if any(exc(j,:) == i)
                Ge=Ge+gexc(j)*Dexc(j);
                Dexc(j) = d*Dexc(j);
            end
        end
    end
   
    if any(inh(:) == i)
        for j = 1:nsyninh
            if any(inh(j,:) == i)
                Gi=Gi+ginh(j)*Dinh(j);
                Dinh(j) = d*Dinh(j);
            end
        end
    end
    if (V(i) > Vthresh)
        spikes(i) = 1;
        V(i) = Vreset;
    end
end
figure;
subplot(3,2,1); hold on;
plot(t,V);
plot(t,Vthresh*spikes,'r*')
subplot(3,2,3); hold on;
for j = 1:nsynexc
    plot([exc(j,:); exc(j,:)],[j j+.75],'k-')
end
for j = 1:nsyninh
    plot([inh(j,:); inh(j,:)],[j j+.75]+nsynexc,'k-')
end
set(gca,'Ylim',[1 sum([nsynexc, nsyninh])+.75]);
subplot(3,2,5);
plot(t,temporalframes,'k-')

%%
% TF, Contrast series
TFs = logspace(log10(1),log10(40),10);
%TFs = 3;
contrasts = linspace(0,.2,5);
contrasts = .35;
nreps = 10;
data = [];
dur = 1;

for i = 1:length(TFs)
    for j = 1:length(contrasts)
        for k = 1:nreps
            k
            V = nan*ones(1,ntimesteps);
            spikes = nan*ones(1,ntimesteps);
            contrast = contrasts(j);
            tf = TFs(i);
            Dexc = zeros(1,nsynexc);
            Dinh = zeros(1,nsyninh);

            sinusoid = contrast*sin(2*pi*tf*linspace(0,ntimesteps*deltaT,ntimesteps));
            envelope = [linspace(0,1,round(.667/4/deltaT)),...
                ones(1,round(.667/2/deltaT)),...
                linspace(1,0,round(.667/4/deltaT))]; % in deltaTs
            
            envelope = [zeros(1,floor((ntimesteps-length(envelope)))/2),...
                envelope,...
                zeros(1,ceil((ntimesteps-length(envelope)))/2)];
            
            temporalframes = envelope.*sinusoid;
            %temporalframes = contrast*ones(1,length(envelope));
            c = conv(temporalframes, centerkernel,'same')*deltaT;
            s = conv(temporalframes, surroundkernel,'same')*deltaT;
            
            % Setting up afferent spikes
            if (ISOLUM == 0)
                kernel = c-s;
            else
                kernel = c+s;
            end
            
            exc_t = [];
            for l = 1:nsynexc
                tmp = unifrnd(0,1,1,ntimesteps) < kernel+(baselinerate*deltaT);
                exc_t{l} = find(tmp);
            end
            exc = zeros(nsynexc,max(cellfun('length',exc_t)));
            for l = 1:nsynexc
                exc(l,1:length(exc_t{l})) = exc_t{l};
            end
            
            inh_t = [];
            for l = 1:nsyninh
                tmp = unifrnd(0,1,1,ntimesteps) < kernel+(baselinerate*deltaT);
                inh_t{l} = find(tmp);
            end
            inh = zeros(nsyninh,max(cellfun('length',inh_t)));
            for l = 1:nsyninh
                inh(l,1:length(inh_t{l})) = inh_t{l};
            end
            % Done setting up afferent spikes

            V(1) = V0;
            for step = 1:ntimesteps
                if (step > 1)
                    V(step) = V(step-1)+(V0-V(step-1)+Ge*(Ve-V(step-1))+Gi*(Vi-V(step-1)))*(deltaT/tm);
                    Ge = Ge-Ge*(deltaT/te);
                    Gi = Gi-Gi*(deltaT/ti);
                    Dinh = Dinh+(1-Dinh).*(deltaT/td);
                    Dexc = Dexc+(1-Dexc).*(deltaT/td);
                end
                
                if any(exc(:) == step)
                    for l = 1:nsynexc
                        if any(exc(l,:) == step)
                            Ge=Ge+gexc(l)*Dexc(l);
                            Dexc(l) = d*Dexc(l);
                        end
                    end
                end
                
                if any(inh(:) == step)
                    for l = 1:nsyninh
                        if any(inh(j,:) == step)
                            Gi=Gi+ginh(l)*Dinh(l);
                            Dinh(l) = d*Dinh(l);
                        end
                    end
                end
                
                if (V(step) > Vthresh)
                    spikes(step) = 1;
                    V(step) = Vreset;
                end
            end
            data(i,j,k) = nansum(spikes)./dur
        end
    end
end
mn = mean(data,3);
sd = std(data,0,3);
figure; subplot(2,1,1);
imagesc(mn);
subplot(2,1,2);
errorbar(TFs,mn,sd);
set(gca,'Xscale','log','Xlim',[TFs(1) TFs(end)])
sd.^2./mn
%%
% Section 21: Propixx  stuff
% Need to load propixx_photometer.at which is in Documents/Grants/Horwitz
% R01 competitve renewal

%diff_t = diff(propixx_voltage_trace(:,1));
%propixx_voltage_trace(diff_t <=0,:) = [];

deltaT = mean(diff(propixx_voltage_trace(:,1))); % ms
figure; axes; hold on;

plot(propixx_voltage_trace(:,1),propixx_voltage_trace(:,2),'b-');
t= [propixx_voltage_trace(1,1):deltaT:deltaT*(length(propixx_voltage_trace)-1)]';
v_interp = interp1(propixx_voltage_trace(:,1),propixx_voltage_trace(:,2),t)
plot(t,v_interp,'g-')


figure; axes; hold on;
skipT = 20;
short_t = t(1:end-skipT);
short_v = propixx_voltage_trace([1:end-skipT],2);


n = 10;
plot([0:deltaT:(n*(short_t(end)+deltaT)-deltaT)],repmat(short_v,[n 1]))

%%
% Section 22
% ThreePulseMonte (L-M vs L+M discrimination code)
stro = nex2stro;
rightwardsaccade = stro.trial(:,strcmp('rightward_sacc',stro.sum.trialFields(1,:)));
Lcc = stro.trial(:,strcmp('lcc',stro.sum.trialFields(1,:)));
Mcc = stro.trial(:,strcmp('mcc',stro.sum.trialFields(1,:)));
uniqueLM = unique([Lcc,Mcc],'rows');
n = zeros(size(uniqueLM,1),1);
k = zeros(size(uniqueLM,1),1);

for i = 1:size(uniqueLM,1)
    L = Lcc == uniqueLM(i,1) & Mcc == uniqueLM(i,2);
    n(i) = sum(L); 
    k(i) = sum(L&rightwardsaccade);
end

figure; 
subplot(2,1,1); hold on;
plot(0,0,'y*');
for i = 1:size(uniqueLM,1)
    h = plot(uniqueLM(i,1),uniqueLM(i,2),'k.');
    set(h,'MarkerSize',n(i)/5+1);
end
axis square;

subplot(2,1,2); hold on;
plot(0,0,'y*');
for i = 1:size(uniqueLM,1)
    ratio = k(i)/n(i); % proprotion choices to the right
    [uniqueLM(i,:) ratio]
    h1 = plot(uniqueLM(i,1),uniqueLM(i,2),'ko','MarkerFaceColor','black');
    h2 = plot(uniqueLM(i,1),uniqueLM(i,2),'ro');
    set(h2,'MarkerSize',abs(70*(ratio-.5))+1);
    if (ratio < .5)
        set(h2,'MarkerEdgeColor','green');
    end
    if (n(i) > 5 && abs(ratio-.5) > 2*sqrt(.25/n(i)))
       set(h2,'LineWidth',2);
    end
end
axis square

% Sedna's thresholds for L-M (.049) and L+M (.053) (individual components)
thetas = atan2(uniqueLM(:,1),uniqueLM(:,2))

uniquethetas = unique(thetas);
figure; axes; hold on;
for i = 1:length(uniquethetas)
    L = thetas == uniquethetas(i);
    h = plot(uniquethetas(i),sum(k(L))./sum(n(L)),'ko');
    if (abs(uniquethetas(i)-3*pi/4) < .01 | abs(uniquethetas(i)-pi/4) < .01)
        set(h,'MarkerEdgeColor','red');
    end
end
plot([pi/2 pi/2],[0 1])

%%
% Section 23 laminar probe stuff
stro = nex2stro(findfile('p032015001'));
%stro = nex2stro(findfile('p031915008'));

%if strcmp(stro.sum.fileName(end-13:end-4),'p031915008'); % Reordering columns 
%    stro.ras = stro.ras(:,[9:16, 1:8, 25:32, 17:24 ,33]);
   % stro.sum.rasterCells = stro.sum.rasterCells([9:16, 1:8, 25:32, 17:24 ,33])
%end

stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
fpacq_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_acq'));
fpon_t= stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_on'));

analogstarttime = [stro.ras{:,strcmp( stro.sum.rasterCells,'anlgStartTime')}]';
arate = stro.sum.analog.storeRates{1};
analogL = strncmp(stro.sum.rasterCells,'AD',2);
% getting rid of eye position
analogL(strcmp(stro.sum.rasterCells,'AD11')) = 0;
analogL(strcmp(stro.sum.rasterCells,'AD12')) = 0;
analogL(strcmp(stro.sum.rasterCells,'AD17')) = 0;  % Channel 1 is obviously bad

analogidxs = find(analogL);
ntrials = size(stro.trial,1);

ntsamps = 400;
t_offset = -.1;
data = zeros(length(analogidxs),ntsamps);

figure;
for i = 1:ntrials
    nsamples = size(stro.ras{i,end-2},1);
    t = [0:nsamples-1]/arate+analogstarttime(i);
    t = t-stimon_t(i);
   % t = t-(fpon_t(i)+stimon_t(i));
    Lt = t > t_offset;
    Lt(find(Lt,1)+ntsamps:end) = 0;
    singletrialmat = [stro.ras{i,analogidxs}]';
    data = data+singletrialmat(:,Lt);
    
    flipud(singletrialmat);  % Is this right?
    singletrialmat = singletrialmat(:,Lt); % Getting rid of all times prior to stim on
    imagesc(singletrialmat);
    drawnow;
   % pause;
    cla;
end

figure;
imagesc(-data);
%times = t(Lt);
b1 = 1/arate;
b0 = -t_offset;
xtickpos = [-.1 0 .1 .2 .3 .4];
set(gca,'XTick',(xtickpos+b0)./b1,'XTicklabel',xtickpos);

%kernel = [1 -2 1]'*normpdf(-3:1:3)
%CSD = conv2(data,kernel,'valid');
%imagesc(CSD)


%%
% Section 24
% FFTs of gratings with (and without eye movements). What do these do to
% the spatiotemporal power.

npixels = 200;
nframes = 100;
x = linspace(0,20*pi,npixels);
phi = zeros(nframes,1);
%phi = linspace(0,2*pi,nframes)';
%phi = linspace(0,4*pi,nframes)';
%phi = [zeros(floor(nframes/2),1);pi*ones(ceil(nframes/2),1)];
%phi = [zeros(floor(nframes/2),1);pi/6*ones(ceil(nframes/2),1)];

stim = sin(repmat(x,nframes,1)+repmat(phi,1,npixels));
colormap(gray);
stimfft = fft2(stim);
figure; 
subplot(2,2,1);
imagesc(stim);
ylabel('t');
xlabel('x');
subplot(2,2,3);
imagesc(fftshift(abs(stimfft)));
colormap(gray);
ylabel('omega t');
xlabel('omega x');

subplot(2,2,4);
plot(sum(fftshift(abs(stimfft)),2),linspace(-nframes/2,nframes/2,nframes)');
title(num2str(sum(sum(abs(stimfft)))));

%% How does a translation effect power as a function of TF
% Map out how a translation affects high and low SFs 
% Make a sum of two gratings
% This isn't particularly useful

npixels = 200;
nframes = 200;
x = normrnd(0,1,1,npixels);
stim = zeros(nframes,npixels);
for i = 1:nframes
    stim(i,:) = x;
    x = [normrnd(0,1), x(1:end-1)];
end
sacjump = 0;
stim(ceil(nframes/2):end,:) = stim(ceil(nframes/2):end,[sacjump+1:npixels,1:sacjump]); % throwing in a saccade
figure; 
subplot(3,2,1); 
imagesc(stim); colormap(gray);
stimfft = fft2(stim);
subplot(3,2,3); 
imagesc(fftshift(abs(stimfft)));
subplot(3,2,4);
plot(sum(fftshift(abs(stimfft)),2),linspace(-nframes/2,nframes/2,nframes)');
title(num2str(sum(sum(abs(stimfft)))));
subplot(3,2,5);
plot(linspace(-npixels/2,npixels/2,npixels), sum(fftshift(abs(stimfft))));


%%
% How does the power change with an instantaneous change in spatial phase?
data = [];
phis = linspace(0,2*pi,500);
npixels = 200;
nframes = 100;
x = linspace(0,10*pi,npixels);
for phi = phis
    phi = [zeros(floor(nframes/2),1);phi*ones(ceil(nframes/2),1)];
    stim = sin(repmat(x,nframes,1)+repmat(phi,1,npixels));
    stimfft = fft2(stim);
    data = [data sum(abs(stimfft(:)))];
end
figure; 
plot(phis,data);
% Interesting: the maximum integrated power occurs when the spatial phase
% changes by XXX radians.

%%
% Section 25: How much to L-M to L+M thresholds change across individual
% observers? If there is a lot of variability here this is important
% because that cannot be due to variability at the level of the cones.
filelists = {'ApolloLMTF.txt','FreyaLMTF.txt','NutLMTF.txt','SednaLMTF.txt'};

startpath = '/Volumes/NO BACKUP/NexFiles/Greg/Sedna';
data = [];

for j = 1:length(filelists)
    flist = flatten(fnamesFromTxt2(fullfile(nexfilepath,'nexfilelists','Greg','LMTF',char(filelists{j}))));
    j
    % Iterate over the list of .nex files.
    for i = 1:length(flist)
        stro = notnex2stro(findfile(flist{i},startpath)); % <-- process the information in the nex file and put it into the "stro" structure
        if ~(abs(stro.sum.exptParams.stim_x) == 50 & abs(stro.sum.exptParams.stim_y) <= 10)
            [stro.sum.exptParams.stim_x stro.sum.exptParams.stim_y]
            continue;
        end
        disp(['Got here ',num2str([stro.sum.exptParams.stim_x stro.sum.exptParams.stim_y])])

        % Below, just figuring out what information is in what column
        Lstim_idx = strcmp(stro.sum.trialFields(1,:), 'stim_idx');
        Llcc = strcmp(stro.sum.trialFields(1,:), 'lcc');
        Lmcc = strcmp(stro.sum.trialFields(1,:), 'mcc');
        Ltf = strcmp(stro.sum.trialFields(1,:), 'tf');
        Loog = strcmp(stro.sum.trialFields(1,:), 'oog');
        
        % Getting the threshold points
        [stim_idxs,init_stim_trial_idxs] = unique(stro.trial(:,Lstim_idx),'last');
        questmodes = stro.trial(init_stim_trial_idxs, Llcc|Lmcc);
        tfs = stro.trial(init_stim_trial_idxs,Ltf);
        
        % Out of gamut checking
        funds = reshape(stro.sum.exptParams.fundamentals,length(stro.sum.exptParams.fundamentals)/3,3);
        if (size(stro.sum.exptParams.mon_spd,1) == 303)
            spds = SplineSpd([380:4:780]',reshape(stro.sum.exptParams.mon_spd,length(stro.sum.exptParams.mon_spd)/3,3),[380:5:780]');
        else
            spds = SplineSpd([380:2:780]',reshape(stro.sum.exptParams.mon_spd,length(stro.sum.exptParams.mon_spd)/3,3),[380:5:780]');
        end
        M = funds'*spds;
        bkgndrgb = stro.sum.exptParams.bkgndrgb;
        [in_gamut,scalar] = gamutCheck([questmodes zeros(size(questmodes,1),1)]', bkgndrgb, M, 'both');
        questmodes(~in_gamut) = questmodes(~in_gamut).*scalar(~in_gamut);
        data = [data; repmat(j,length(tfs),1) questmodes tfs ~in_gamut'];
    end
end

% Analysis and plotting
TFthreshs = [1 2];
Loog = logical(data(:,end));
subjectidxs = unique(data(:,1))
figure; axes; hold on;
symbols = {'bs','rd','gv','ko'};
params = [];
for j = subjectidxs'
   if (j == 3 | j == 1)
       continue
   end
    Lsub = data(:,1) == j;
    Ltf = data(:,4) >= TFthreshs(1) & data(:,4) <= TFthreshs(2);
    lm = data(Lsub&Ltf,[2 3]);
    Loog = data(Lsub&Ltf,5);
    x = lm(:,1);
    y = lm(:,2);
    
    initparams = [2*std(x) 2*std(y) acos((x./norm(x))'*y./norm(y))]; %first guess for fminsearch parameters
    options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-5,'Display','none');
    [fpar,fv] = fminsearch(@(params) ellipsefiterr(params,[x y],Loog),initparams, options); %fv not necessary

    % Plotting the fit
    angles = linspace(0,2*pi,100);
    R = (fpar(2)^2-fpar(1)^2)*cos(2*angles-2*fpar(3))+fpar(1)^2+fpar(2)^2;
    Q = sqrt(2)*fpar(1)*fpar(2)*sqrt(R);
    r = Q./R;
    [tmpx,tmpy] = pol2cart(angles,r);

    h = plot([x; -x],[y; -y],symbols{j});
    markercolor = get(h,'Color');
    set(h,'MarkerFaceColor',markercolor);
    h = plot(tmpx,tmpy,'m-','LineWidth',2);
    set(h,'Color',markercolor);
    params = [params; fpar];
end
axis square
axis equal


%%
% Section 26
% Hacking around with Charlie's cone model.
% Simulating responses to stimuli used in an IsoSamp experiment.
% All stimuli equated for psychophysical detection threshold assuming a
% Weibull function. Threshold => p(correct) = 0.8161. See update_Quest.m in
% LMTF/supportfunctions. This implies a d' of 1.27 in a 2AFC
% experiment.

filename = 'A053017002';
%filename = 'A060217005' % Magnocll
%filename = 'A052517002'; % Very close to fovea
fundamentals = load('T_cones_smj10');
params.stro = nex2stro(findfile(filename));
cal.gamma = params.stro.sum.exptParams.gamma_table; % Only used for "dtv1". Can ignore.
spds = params.stro.sum.exptParams.mon_spd;
spds = SplineSpd([380:4:780]',reshape(spds,length(spds)/3,3),[380:5:780]');
cal.monSpect = spds(:);
M = fundamentals.T_cones_smj10*spds;
cal.Mmtx = M(:);
cal.frameRate = params.stro.sum.exptParams.framerate;
cal.bkgndrgb = params.stro.sum.exptParams.bkgndrgb';
cal.fname = filename;
cal.monSpectWavelengths = linspace(380,780,101);
cal.pixperdeg = params.stro.sum.exptParams.pixperdeg;

params.runType = 'isosamp';
params.obsMethod = 'obsMethod_filteredWtFxn';
params.impulseResponse = 'rieke';
params.DTV1_fname = [];
params.DTNT_fname = [];
params.monCalFile = cal;
params.unitTest = false;
params.eqMosaic = false;
params.saveDir = '/Users/greghorwitz/Documents/MATLAB';
params.notes = 'IsoSamp test';
params.parallelOperations = false;
params.eyeType = 'monkey';
params.coneSampRate = 2400;
params.flatPowerSpect = false;
params.enableScones = false;
params.sacamp_deg = 0;
params.sacdur_s = 0;

[gab, cones, mon, idlob, params] = DTcones_gh(params,0);

data = [];
for i = 1:size(idlob.analyticMean,1) % looping over color direction
    for j = 1:size(idlob.analyticMean(i,:),2) % looping over contrast/TF
        if ~isempty(idlob.analyticMean{i,j})
            tmp_lm_mu = idlob.analyticMean{i,j}([1 2]);
            tmp_lm_sigma = sqrt(idlob.analyticVar{i,j}([1 2]));

            tf = gab.driftRates{i}(j);
            data = [data; gab.colorDirs(i,[1 2]) tf tmp_lm_mu tmp_lm_sigma];
        end
    end
end
% columns of data: L, M, TF, mu_L, mu_M, sigma_L, sigma_M
% Need to compute a d'-like measure of SNR
figure; subplot(2,1,1); hold on;
L = data(:,1) > 0;
plot(data(L,3),sqrt(sum(data(L,[4 5]).^2,2)),'ko-');
plot(data(~L,3),sqrt(sum(data(~L,[4 5]).^2,2)),'ro-');
plot(data(L,3),sqrt(sum(data(L,[6 7]).^2,2)),'k:');
plot(data(~L,3),sqrt(sum(data(~L,[6 7]).^2,2)),'r:');
set(gca,'Xscale','log');

subplot(2,1,2); hold on;
plot(data(L,3),sqrt(sum(data(L,[4 5]).^2,2))./sqrt(sum(data(L,[6 7]).^2,2)),'ko-','MarkerFaceColor','black');
plot(data(~L,3),sqrt(sum(data(~L,[4 5]).^2,2))./sqrt(sum(data(~L,[6 7]).^2,2)),'ro-','MarkerFaceColor','red');

set(gca,'Xscale','log');
xlabel('TF (Hz)')
ylabel('SNR (d'')')

subplot(2,1,1); title(filename,'FontSize',16);
%%
% Section 27

% -------------------------------------
% What happens if I *don't* normalize the area under 
% the filter in the Watson 1986 model?
% The ratio of the integrals of h1 and h2 stays fixed.
% Larger taus leads to larger integrals. 
% GDLH 8/29/16
% -------------------------------------

% General parameters
nsamps = 1000;
t = linspace(0,1,nsamps); % time in s
timestep = mean(diff(t));

xi = 1;
zeta = .1;
n1 = 5;
kappa = 1.5;
n2 = 5;

figure;
for tau1 = .001:.001:.01
    tau2 = kappa*tau1;
    h1 = (tau1.^(n1-1)./factorial(n1-1)).*(t./tau1).^(n1-1).*exp(-t./tau1);
    h2 = (tau2.^(n2-1)./factorial(n2-1)).*(t./tau2).^(n2-1).*exp(-t./tau2);
    h = xi*(h1-zeta*h2);
    subplot(2,1,1); hold on;
    plot(t,h,'k-');   
    
    Fnum = (1/length(t))*fftshift(fft(h));
    subplot(2,1,2); hold on;
    o = linspace(0,0.5*nsamps/(max(t)-min(t)), length(Fnum(length(Fnum)/2:end)));
    plot(o,abs(Fnum(length(Fnum)/2:end)),'b-')
    set(gca,'Xscale','log','Yscale','log')
end
xlabel('Time (s)');
ylabel('Amplitude');
% ---------------------------------------
% Numerical test of convolution for a single filter
% Just trying to figure out the analytical form
% of the n-fold convolution of e^-(t/tau)
% ---------------------------------------

f = exp(-t./tau1);
for i = 1:n1
    if i == 1
        numconv = f;
    else 
        numconv = conv(numconv,f)*timestep;
    end
end
numconv = numconv(1:length(t)); % Numerical convolution
figure; axes; hold on;
%plot(t,numconv./sum(numconv));
plot(t,numconv);
h1 = (tau1.^(n1-1)./factorial(n1-1)).*(t./tau1).^(n1-1).*exp(-t./tau1); % analytical convolution
plot(t,h1);
%plot(t,h1./sum(h1));

% Yes, h1 is the analytical formula for the n-fold convolution of 
% exp(-t/tau) (in shape) but I can't get the areas to match perfectly.

% ---------------------------------------
% Now let's take a look at some real data
% And see how TCSFs vary with eccentricity
% ---------------------------------------
load /Users/greghorwitz/Desktop/MatlabCode/Zack/IsoSamp/private/data/LMTF.mat
uniqueXYs = A.eccs;
models = A.legacy.firstroundmodels;

LUM = 0;
modelparamidxoffset = 0;
if ~LUM
    modelparamidxoffset = 6;
end
omega = logspace(0,2,100);
figure; axes; hold on;
for ecc_counter = 1:size(uniqueXYs,1)
    xi = models(1+modelparamidxoffset,ecc_counter);
    zeta = models(2+modelparamidxoffset,ecc_counter);
    n1 = models(3+modelparamidxoffset,ecc_counter);
    n2 = n1+models(4+modelparamidxoffset,ecc_counter);
    tau1 = models(5+modelparamidxoffset,ecc_counter);
    tau2 = tau1+models(6+modelparamidxoffset,ecc_counter);
    
    f1 = (1i*2*pi*10^tau1.*omega+1).^-n1;
    f2 = (1i*2*pi*10^tau2.*omega+1).^-n2;
    f = xi*(f1-zeta*f2);
    h = plot(omega,abs(f),'k-');
    set(h,'LineWidth',2)
    color = [abs(uniqueXYs(ecc_counter,1)/100) 0 abs(uniqueXYs(ecc_counter,2)/100)];
    set(h,'color',color);
    % increasing red = increasing x
    % increasing blue = increasing y
end
set(gca,'XScale','log','YScale','log')

%%
% Section 28
% Potential problem in tf_fiterr2.m. I'm taking a weighted sum of mechanism
% activations weighted by the *threshold* of the mechanism. This doesn't
% make sense. I should be weighting each by the *sensitivity* of the
% mechanism and then inverting at the end to get a threshold.

lm = [cos(linspace(0,2*pi,60))' sin(linspace(0,2*pi,60))'];
theta = 0;
a = 20; % lum sensitivity
b = 25; % rg sensitivity
mechs = [cos(theta) 1/sqrt(2); sin(theta) -1/sqrt(2)];
f = [];
for i = 1:size(lm,1)
    dotprods = lm(i,:)*mechs*diag([a b]);
    f(i) = sqrt(sum(dotprods.^2))^-1; % converting sensitivies to thresholds
end
figure
polar(linspace(0,2*pi,60),f)

%%
% Section 29
% Looking at a tilted rampy trough in (x,y) 
load(fullfile(fileparts(which('IsoSampOnline.m')), 'private', 'data', 'LMTF.mat'))
params = E.legacy.mode5params;
% params (mode 5):  The first eleven elements are the ones that are common to every retinal location:
%         [zeta_LUM, n1_LUM, delta n_LUM, tau1_LUM, kappa_LUM, zeta_RG, n1_RG, delta n_RG, tau1_RG, kappa_RG, theta]
%         After that, come [b0, b1, b2, b3, a0, a1, a3].
regression_coeffs = params(12:15);
[r,phi] = meshgrid(linspace(2,10,40),linspace(-pi/2,pi/2,40));
X = [ones(numel(r),1) r(:) phi(:) phi(:).^2];
sensitivity = reshape(X*regression_coeffs, size(r));
[x,y] = pol2cart(phi, r);
L = double(phi>-pi/2.5 & phi < pi/2.5); % Couldn't do LMTF on the vertical meridian so can't estimate sensitivity there.
L(L==0) = nan;
figure;
surf(x.*L,y.*L,sensitivity.*L);
axis equal
set(gca,'View',[0 90]);
title('Sensitivity');

% A better approximation might have been that isosensitivity contours are ellipses
% Sensitivity decreases with eccentricity
new_sensitivity = 10-sqrt(.5*x.^2+y.^2+.2*x.*y)
figure; subplot(2,1,1);
surf(x,y,new_sensitivity);
axis equal
set(gca,'View',[0 90]);
title('Sensitivity');
subplot(2,1,2); hold on;
surf(phi,r,new_sensitivity,64*ones(size(r))); % yellow
% Can I fit this thing with sensitivity = b0+b1*r+b2*r*cos(2*phi)? <- effect
% of phi decreases at small r

X = [ones(numel(phi),1),r(:),r(:).*cos(2*phi(:)),r(:).*sin(2*phi(:))]; % <---- New tilted rampy trough
[b,bint,~,~,stats] = regress(new_sensitivity(:),X);
pred = X*b;
surf(phi,r,reshape(pred,size(phi)),ones(size(r))); % blue

% Looking at the rampy trough model in (x,y)
figure; subplot(2,1,1);
surf(x,y,sensitivity);
axis equal
set(gca,'View',[0 90]);
title('Sensitivity');
subplot(2,1,2);
surf(phi,r,sensitivity);

% The new model is sensitivity = b0+b1*r+b2*r*cos(2*phi)*b3*r*sin(2*phi)
% What does this look like for various values of the parameters?
b = [-10 -1 .1 -.1]';
surf(x,y,reshape(X*b,size(x)));
axis equal;
set(gca,'View',[0 90]);

basisvector0 = r;
basisvector1 = r.*cos(2*phi);
basisvector2 = r.*sin(2*phi);

figure; 
subplot(2,2,1);
surf(x,y,-basisvector0);
view(0,90); axis equal;
subplot(2,2,2);
surf(x,y,basisvector1);
view(0,90); axis equal;

subplot(2,2,3);
surf(x,y,-basisvector2);
view(0,90); axis equal;

%%
% Section 30
% Fitting a series of exponentials to the Packer et al. (1989) cone density
% data. This will allow me to average nasal and temporal visual fields.
MMPERDEG = 0.223; % mm/deg (Perry and Cowey 1985)

[num,txt,raw]=xlsread('PackerEtAlConeDensities')
nasaldata = num(:,[1 2]);
temporaldata = num(:,[4 5]);
figure; axes; hold on;
plot(nasaldata(:,1)./MMPERDEG,nasaldata(:,2)*1000,'ko');
plot(temporaldata(:,1)./MMPERDEG,temporaldata(:,2)*1000,'bs');
xlabel('ecc (deg)');
ylabel('cones/mm^2');
x = linspace(.1,70,100);

% Plotting the fit from Goodchild et al.
temporalcoeffs=[150.9676 -1.2220 35.9979 -0.1567 9.9936 -0.0258]; % for macaque (Goodchild et al., 1996)

fun = @(coeffs,x)(coeffs(1).*(exp(coeffs(2).*x)))+...
        (coeffs(3).*(exp(coeffs(4).*x)))+...
        (coeffs(5).*(exp(coeffs(6).*x)));
plot(x,1000*fun(temporalcoeffs,x),'b-');
% Now fitting the cone density along the nasal retina
nasalcoeffs = nlinfit(nasaldata(:,1)./MMPERDEG,nasaldata(:,2),fun,temporalcoeffs); % Fitting cones/mm^2*1000 a la Goodchild 
plot(x,1000*fun(nasalcoeffs,x),'k-');

% Now fitting both sets of data to try to get an average curve for stimuli
% along the horizontal meridian
Lnasal = nasaldata(:,1) <= 15; % only fitting out to 15 mm
Ltemporal = temporaldata(:,1) <= 15; % only fitting out to 15 mm
mergedcoeffs = nlinfit([nasaldata(Lnasal,1); temporaldata(Ltemporal,1)]./MMPERDEG,[nasaldata(Lnasal,2); temporaldata(Ltemporal,2)],fun,temporalcoeffs); % Fitting cones/mm^2*1000 a la Goodchild 
plot(x,1000*fun(mergedcoeffs,x),'r-');

% Using nasal and temporal coeffs to predict cone density across the retina
% of each eye (use nasal retina coefficients for veritical meridian).

[rfx,rfy] = meshgrid(linspace(-7,7,100),linspace(-7,7,100));
angles = atan2(rfy,rfx);
rs = sqrt(rfx.^2 + rfy.^2);

Lx = rfx>0;
righteyedensity_h = zeros(size(rfx));
righteyedensity_h(Lx) = cos(angles(Lx)).^2.*fun(temporalcoeffs,rs(Lx));
righteyedensity_h(~Lx) = cos(angles(~Lx)).^2.*fun(nasalcoeffs,rs(~Lx));
righteyedensity_v = sin(angles).^2.*fun(temporalcoeffs,rs);
righteyeconedensity = (righteyedensity_h+righteyedensity_v); % sin^2+cos^2 = 1
surf(rfx,rfy,righteyeconedensity*1000)

% Making a function that will return the (binocular) cone density at a
% given RF location
[phi,r] = cart2pol(rfx,rfy);
conedensity = 1000*((cos(phi).^2.*(fun(temporalcoeffs,r)+fun(nasalcoeffs,r))/2) + (sin(phi).^2.*fun(temporalcoeffs,r)));
surf(rfx,rfy,conedensity)

%%
% Section 31
% Fitting the Watanabe and Rodieck dendritic field sizes
% This will be useful for getting the ratio between midget and parasol cell
% RFs as a function of eccentricity
% Need to load a and b (from "WatanabeAndRodieckDendriticFieldDiameters")
% This is useless because Gauthier et al. 2009 show that RF size is hard to
% relate to DF size.?
% Can I use the ratio of midget:parasol DF size as a proxy for midget:parasol RF size

MMPERDEG = 0.223; % mm/deg (Perry and Cowey 1985)
DEGPERMM = 1/MMPERDEG;

figure;
axes; hold on;
plot(DEGPERMM*a(:,1),a(:,2),'.');
plot(DEGPERMM*b(:,1),b(:,2),'o');
set(gca,'Yscale','log');
ylabel('Dendritic field diameter �m');
xlabel('Nasal equivalent eccentrciity (DVA)');

b_parvo = regress(log10(a(:,2)),[DEGPERMM*a(:,1) ones(size(a,1),1)]);
b_magno = regress(log10(b(:,2)),[DEGPERMM*b(:,1) ones(size(b,1),1)]);

plot([1:30],10.^([1:30]*b_parvo(1)+b_parvo(2)),'r-');
plot([1:30],10.^([1:30]*b_magno(1)+b_magno(2)),'k-');

% Dendritic field area
figure; axes; hold on;
midget_DF_area = pi*(DEGPERMM*a(:,2)./2).^2;
parasol_DF_area = pi*(DEGPERMM*b(:,2)./2).^2;

plot(DEGPERMM*a(:,1),midget_DF_area,'.');
plot(DEGPERMM*b(:,1),parasol_DF_area,'o');
tmp = 10.^([1:30]*b_parvo(1)+b_parvo(2));
plot([1:30],pi*(DEGPERMM*tmp./2).^2,'r-');
tmp = 10.^([1:30]*b_magno(1)+b_magno(2));
plot([1:30],pi*(DEGPERMM*tmp./2).^2,'k-');
set(gca,'Yscale','log')
ylabel('Dendritic field area (DVA^2)');
xlabel('Nasal equivalent eccentrciity (DVA)');

plot([1:30],(10.^([1:30]*b_magno(1)+b_magno(2)))./(10.^([1:30]*b_parvo(1)+b_parvo(2)))) % Ratio of RF diameters
plot([1:30],(10.^([1:30]*b_magno(1)+b_magno(2))).^2./(10.^([1:30]*b_parvo(1)+b_parvo(2))).^2) % Ratio of RF areas

% At large eccentricities, parasol cell DFs are 2.8x midget cell DFs. But
% physiology tells us that their parasol cell RFs are 2x midget cell RFs.
% At small eccentricities (according to Watanabe and Rodieck), parasol cell
% DFs are 5.3x midget cell DFs. This fits with Dacey's assertion that
% midget:parasol ratio maybe as high at 30 in fovea (5.3^2 = 28)

% EJ suggests relying on Croner and Kaplan's RF data but their parasol data
% is almost entirely > 10� and their data are noisy (or maybe parasol cells
% are heterogeneous). They make the point that their data agree with Perry
% et al.'s (1984) dendritic field sizes. So I fit a line to Perry et al.'s data
% and divide by sqrt(2) (that's what I have to do with Croner and Kaplan's
% data to get it into standard deviations).

% To reconstruct Croner and Kaplan Figure 13a

ecc_to_diam_deg = @(rf_r_deg) 0.0683+0.0165*rf_r_deg/sqrt(2); % Perry et al. "/sqrt(2)" because that's what I have to do with Croner and Kaplan's data, and their peripheral data agrees with Perry et al. (more central)
% This converts RF center position (DVA) to the standard deviation of a
% Gaussian RF.
r = linspace(0,30,30)
plot(r,ecc_to_diam_deg(r)/2*sqrt(2)); 
% Dividing by 2 to get from RF diameter to radius and multiplying by sqrt(2)
% to get numbers in "rc" units (comparable to Croner and Kaplan-the
% distance at which sensitivity has dropped by a factor of 1/e).

%%
% Section 32
% Simulating detection using spatial uncertainty, spatial correlations,
% optimal filter. 2D simulation.

npix = 50;
lambda = 10; % spatial correlation parameters
tmp = 1:npix; % change to microns or mms or somthing like that
tmp = tmp-round(length(tmp)/2);
x = repmat(tmp,npix,1);
y = x';
d = squareform(pdist([x(:) y(:)],'euclidean'));
spatialcorrelationmat = exp(-abs(d./lambda));

% testing
overlaps = xcorr2(ones(npix,npix)); % computing the overlap for various shifts in the autocorrelation function.
niter = 100;
runningxcorr = zeros(2*npix-1);
PLOTNOISESAMPLES = 1;
if PLOTNOISESAMPLES
    figure;
end
for i = 1:niter
    a = mvnrnd(zeros(npix.^2,1),spatialcorrelationmat);
    if PLOTNOISESAMPLES & floor(i/10) == i/10
        imagesc(reshape(a,npix,npix));
        drawnow;
    end
    runningxcorr = runningxcorr + xcorr2(reshape(a,npix,npix))./overlaps;
end
figure; 
subplot(2,2,1);
imagesc(runningxcorr./niter);
axis square;
title('Simulation');
subplot(2,2,2); hold on;
plot(runningxcorr(npix,:)./niter,'b-');
axis square;
tmp = 1:2*npix;
tmp = tmp-round(length(tmp)/2);
plot(exp(-abs(tmp/lambda)),'r:');
subplot(2,2,3);
d = sqrt(repmat(tmp,2*npix,1).^2+repmat(tmp',1,2*npix).^2);
imagesc(exp(-(abs(d)/lambda)));
axis square;
title('Theory');

% Finding the optimal filter
% requires defining the signal
sigmastimulus = 3; % 6
tmp = 1:npix;
tmp = tmp-round(length(tmp)/2);
sig1D = normpdf(tmp,0,sigmastimulus);
sig1D = sig1D./max(sig1D);
sig2D = sig1D'*sig1D;
optimalfilter = reshape(pinv(spatialcorrelationmat)*sig2D(:),size(sig2D)); % easier than I thought
optimalfilter = optimalfilter./norm(optimalfilter(:));
% Manually overriding the optimal filter
surround = normpdf(tmp,0,5*sigmastimulus)'*normpdf(tmp,0,5*sigmastimulus);
optimalfilter = sig2D./sum(sum(sig2D))-surround./sum(sum(surround));

figure; 
subplot(2,1,1);
imagesc(optimalfilter);
axis square;
subplot(2,1,2);
plot(optimalfilter(25,:));
axis square;

% Now simulating a psychophysical experiment with a matched filter and
% correlated noise.
OPTORADIUS = 10;
OPTOCENTER = [0 5*sigmastimulus]; % [20 20] is outside, [0 12] is on the edge 
tmp = 1:npix;
tmp = tmp-round(length(tmp)/2);
d = sqrt(repmat(tmp-OPTOCENTER(1),npix,1).^2+repmat(tmp'-OPTOCENTER(2),1,npix).^2);
optostimmask  = d > OPTORADIUS;

stimstrengths = logspace(-2,log10(2),7);
ntrials = 2000;
data =[];
for i = 1:length(stimstrengths)
    noise = mvnrnd(zeros(npix.^2,1),spatialcorrelationmat,ntrials);
    for j = 1:ntrials
        stim = stimstrengths(i)*sig2D+reshape(noise(j,:), size(sig2D));
        stim = stim.*optostimmask;
        dv = stim(:)'*optimalfilter(:);
        data(i,j) = dv;
    end
end
% For simulating a psychometric function
threshold = 0;
mn = sum(data>=threshold,2)./ntrials;
sterr = mn.*(1-mn)/sqrt(ntrials);
%dprime = (mean(data,2)-mean(data(1,:)))./sqrt(mean(var(data,[],2)));
figure; subplot(2,1,1);
imagesc(optimalfilter+optostimmask)
axis square;
title('Matched filter model');
subplot(2,1,2); hold on;
errorbar(stimstrengths,mn,sterr,'o');
%plot(stimstrengths,dprime,'ko-')
set(gca,'Ylim',[0.4 1],'Xlim',[.01 max(stimstrengths)]);
axis square;
set(gca,'Xscale','log')
% -------------------------------------------
% Model spatial uncertainty by
% max stim?
% using several stimulus-sized filters and taking the max
% Could alternatively move the optimal filter around and take the max

stimstrengths = logspace(log10(.2),log10(10),7);
data =[];
for i = 1:length(stimstrengths)
    noise = mvnrnd(zeros(npix.^2,1),spatialcorrelationmat,ntrials);
    for j = 1:ntrials
        stim = stimstrengths(i)*sig2D+reshape(noise(j,:), size(sig2D));
        stim = stim.*optostimmask;
        dv = max(max(conv2(stim,sig2D,'same')));
        data(i,j) = dv;
    end
end
% Distributions of maxes aren't Gaussian so d-prime here is a little weird
% Getting the threshold (no optostim trials)
tmpthresholddata =[];
for j = 1:ntrials
    dv = max(max(conv2(reshape(noise(j,:),size(sig2D)),sig2D,'same')));
    tmpthresholddata = [tmpthresholddata; dv];
end
threshold = median(tmpthresholddata);

figure; subplot(2,1,1);
imagesc(optimalfilter+optostimmask)
axis square;
title('Uncertainty model');
subplot(2,1,2); hold on;

mn = mean(data>threshold,2);
sterr = mn.*(1-mn)/sqrt(ntrials);
errorbar(stimstrengths,mn,sterr,'o');
set(gca,'Ylim',[0.4 1]);
axis square;
set(gca,'Xscale','log')

%%
% Section 33
% As above, but working with a simple, toy model
% a npix x npix grid on the cortical surface
% Bottom line seems to be that to show a marked improvement in the spatial
% uncertainty model with optostim, it's important to have a *lot* of
% irrelevant channels. Also surround needs to be strong.

npix = 5;
signal_template = zeros(npix);
signal_template(ceil(npix/2),ceil(npix/2)) = 1;
noise_template = ones(npix);
center_surround_filter = signal_template-(noise_template./sum(noise_template(:)));
signalamplitudes = [0 logspace(-1.5,1.5,10)];
niter = 5000;
optomask = ones(npix);
optomask(ceil(npix/2)+1:end,:) = 0;
data = [];
for i = 1:length(signalamplitudes)
    for j = 1:niter
        noise = 10+normrnd(0,2.5)*noise_template+normrnd(0,1,size(noise_template));
        sig_plus_noise = signalamplitudes(i)*signal_template+noise;
        
        dv1 = center_surround_filter(:)'*sig_plus_noise(:);
        dv2 = max(sig_plus_noise(:));
        dv3 = sum(sum(center_surround_filter.*sig_plus_noise.*optomask));
        dv4 = max(max(sig_plus_noise.*optomask));
        data = [data; signalamplitudes(i) dv1 dv2 dv3 dv4];
    end
end

% Order of columns
% 1) signal amplitude
% 2) matched filter without optomask
% 3) maximum without optomask
% 4) matched filter with optomask
% 5) maximum with optomask

Lblank = data(:,1) == 0;
uniquestim = unique(data(:,1));
uniquestim(uniquestim == 0) = [];
titles = {'matched filter, no mask','maximum, no mask','matched filter, with mask','maximum with mask'};
metadata = [];
for i = 1:length(uniquestim)
    figure;
    for j = 1:4
        subplot(2,2,j); hold on;
        L = data(:,1) == uniquestim(i);
        hist([data(Lblank,j+1) data(L,j+1)]);
        r = roc(data(Lblank,j+1), data(L,j+1));
        title([titles{j},' roc = ',num2str(r)]);
        metadata(i,j) = r;
    end
end

figure; 
subplot(2,1,1); hold on;
plot(uniquestim,metadata(:,1),'ko-');
plot(uniquestim,metadata(:,3),'bo-');
set(gca,'Xscale','log')
subplot(2,1,2); hold on;
plot(uniquestim,metadata(:,2),'ko-');
plot(uniquestim,metadata(:,4),'bo-');
set(gca,'Xscale','log')


%%
% Section 34. Trying to get preliminary data for Aim 2
conn = database('Abhishek','','','Vendor','MySql','Server','128.95.153.12');
location_query = 'SELECT filename FROM DToneloc WHERE training = ''no''';
filenames = fetch(conn, location_query);
close(conn)
nfp = [nexfilepath,filesep,'Abhishek'];
data = [];

for fileidx = 1:length(filenames)
    disp([fileidx length(filenames)]);
    stro = nex2stro(findfile(filenames{fileidx}, nfp));
    Lstimpresent = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimpresent'));
    Lcorrect = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    Llaser = logical(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'optstim')));    
    lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'lcc'));
    mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'mcc'));
    scc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'scc'));
    oog = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'oog'));
    stim_idx = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_idx'));
    figure; set(gcf,'Name',filenames{fileidx});
    for i = unique(stim_idx)'
        L = stim_idx == i;
        stimuli_laser = [lcc(L&Lstimpresent&Llaser) mcc(L&Lstimpresent&Llaser) scc(L&Lstimpresent&Llaser)];
        stimuli_nolaser = [lcc(L&Lstimpresent&~Llaser) mcc(L&Lstimpresent&~Llaser) scc(L&Lstimpresent&~Llaser)];
        colordir = stimuli_nolaser(1,:)./norm(stimuli_nolaser(1,:));
        % numbers of different types of corrects/incorrects
        nhit_laser = sum(L&Lcorrect&Lstimpresent&Llaser);
        nhit_nolaser = sum(L&Lcorrect&Lstimpresent&~Llaser);
        nFA_laser = sum(L&~Lcorrect&~Lstimpresent&Llaser);
        nFA_nolaser = sum(L&~Lcorrect&~Lstimpresent&~Llaser);
        nmiss_laser = sum(L&~Lcorrect&Lstimpresent&Llaser);
        nmiss_nolaser = sum(L&~Lcorrect&Lstimpresent&~Llaser);
        nCR_laser = sum(L&Lcorrect&~Lstimpresent&Llaser);
        nCR_nolaser = sum(L&Lcorrect&~Lstimpresent&~Llaser);
        subplot(2,3,i+1); hold on;
        title(num2str(colordir))
        plot(sqrt(sum(stimuli_laser.^2,2)),'bo-','MarkerFaceColor','blue');
        plot(sqrt(sum(stimuli_nolaser.^2,2)),'ko-','MarkerFaceColor','black');
        set(gca,'Yscale','log');
        subplot(2,3,i+4)
        h = bar([nhit_laser nhit_nolaser; nmiss_laser nmiss_nolaser; nFA_laser nFA_nolaser; nCR_laser nCR_nolaser]);
        set(h(1),'FaceColor','blue');
        set(h(2),'FaceColor',[.5 .5 .5])
        set(gca,'Xticklabel',{'Hits','Misses','FAs','CRs'},'Xlim',[0 5]);
        ylabel('trial count');
    end
end
   
%%
% Section 35
% Looking a linear superpositions of an L+M grating and an L-M grating with
% different spatial phase offsets.

load Dell4BitsCal.mat
cal = cals{end};
load('T_cones_smj10');
funds = T_cones_smj10;
spds = SplineSpd([380:4:780]',cals{end}.P_device,[380:5:780]');
M = funds*spds;
bkgndRGB = round(cals{end}.bgColor*255)+1;
bkgndrgb = diag(cals{end}.gammaTable(bkgndRGB,:));
invgammaTable = InvertGammaTable(cals{end}.gammaInput, cals{end}.gammaTable, 256);

bkgndlms = M*bkgndrgb;

LMcc = [.3 .3 .3]';
%LMcc = [0 0 0]';
LvMcc = [.04 -.04 0]';
LvMcc = [.01 -.01 0]';

LMS_lm_delta = bkgndlms.*LMcc;
LMS_lvm_delta = bkgndlms.*LvMcc;
rgb_lm_delta = inv(M)*LMS_lm_delta;
rgb_lvm_delta = inv(M)*LMS_lvm_delta;

% L+M in cosine phase
npix = 100;
phaseoffsets = linspace(0,2*pi,9); phaseoffsets(end) = [];
figure;
for phaseoffset = phaseoffsets
    grat_template = repmat(cos(linspace(0,6*pi,npix)),npix,1,3);
    rgb_lm_delta_mat = grat_template.*repmat(permute(rgb_lm_delta,[3 2 1]),npix,npix,1);
    grat_template = repmat(cos(linspace(0,6*pi,npix)+phaseoffset),npix,1,3);
    rgb_lvm_delta_mat = grat_template.*repmat(permute(rgb_lvm_delta,[3 2 1]),npix,npix,1);
    
    im_rgb = repmat(permute(bkgndrgb,[3 2 1]),npix,npix,1)+rgb_lm_delta_mat+rgb_lvm_delta_mat;
    if any(im_rgb(:) > 1 | im_rgb(:) < 0)
        error('rgbs out of range');
    end
    im_RGB = zeros(size(im_rgb));
    for i = 1:3
        im_RGB(:,:,i) = reshape(invgammaTable(round(im_rgb(:,:,i)*255)+1,i),npix,npix,1);
    end
    subplot(ceil(sqrt(length(phaseoffsets))), ceil(sqrt(length(phaseoffsets))), length(get(gcf,'Children'))+1)
    image(im_RGB); axis square;
    title(num2str(phaseoffset));
end
%%
% examining how different amplitudes and phases of L+M and L-M modulate L
% and M cones.
LM_phase = pi/2;
LvM_phase = 0;
LM_amp = 1;
LvM_amp = 1;

LM = LM_amp*sin(linspace(0,6*pi,100)+LM_phase);
LvM = LvM_amp*sin(linspace(0,6*pi,100)+LvM_phase);
L = LM+LvM;
M = LM-LvM;
figure; axes; hold on;
plot(L,'r-');
plot(M,'g-');


%%
% Section 36
% Rasters showing rapid excitation by DLX-ChR2
% Take from DLXGrantFig Section 9

filename = 'A012319005.nex';

stro = nex2stro(findfile(char(filename)));
stimpresentidx  = strcmp(stro.sum.trialFields(1,:),'stimpresent');
Lcc = strcmp(stro.sum.trialFields(1,:),'lcc');
Mcc = strcmp(stro.sum.trialFields(1,:),'mcc');
Scc = strcmp(stro.sum.trialFields(1,:),'scc');
tf = strcmp(stro.sum.trialFields(1,:),'tf');
oog = strcmp(stro.sum.trialFields(1,:),'oog');
optstim = strcmp(stro.sum.trialFields(1,:),'optstim');
correct = strcmp(stro.sum.trialFields(1,:),'correct');
laseron = strcmp(stro.sum.trialFields(1,:),'laseron_t');
laseroff = strcmp(stro.sum.trialFields(1,:),'laseroff_t');
analogstrtimeind = strcmp(stro.sum.rasterCells(1,:),'anlgStartTime');
spikeind = strcmp(stro.sum.rasterCells(1,:),'sig001a');
lasertraceind = strcmp(stro.sum.rasterCells(1,:),'AD13');
stimpresent = logical(stro.trial(:,stimpresentidx));
lasertrials = logical(stro.trial(:,optstim));
samplerate = stro.sum.analog.storeRates{3}; % sample rate at which laser analog pulses are stored in file

% Stimulus absent trials
idxs = find(~stimpresent);
counter = 0;
figprefs;

quickticks = 0:.0005:.002;
slowticks = -.1:.1:.45;

hax = [];
hax(1) = axes('Position',[2 2 4 4]);
hax(2) = axes('Position',[8 2 4 4]);

for ii = 1:numel(idxs) % looping over stimulus absent trial
    ind = idxs(ii);
    if lasertrials(ind)
        analogstartime = stro.ras{ind,analogstrtimeind};
        % Pulling out the laseron and off timings from analog traces
        t = stro.ras{ind,analogstrtimeind}+[0:1:length(stro.ras{ind,lasertraceind})-1]/samplerate;
        laserontime = t(find(stro.ras{ind,lasertraceind}>0.1,1));
        laserofftime = t(find(stro.ras{ind,lasertraceind}>0.1,1,'last'));
        spiketimes = stro.ras{ind,spikeind}-laserontime;
        for i = 1:2
            axes(hax(i));
            if i == 1
                ticks = slowticks;
            else
                ticks = quickticks;
            end
            spiketimes = spiketimes(spiketimes>ticks(1) & spiketimes<ticks(end));
            nspikestot = length(spiketimes);
            plot([spiketimes spiketimes]',[zeros(nspikestot,1) .7*ones(nspikestot,1)]'+counter,'k-','linewidth',1);
        end
        counter = counter + 1;
    end
end
axes(hax(1));
set(gca,'Xlim',[-0.15 0.45],'Tickdir','out','XTick',[-0.15:0.15:0.45]); xlabel('time (s)'); ylabel('Number of trials'); axis square;

axes(hax(2));
set(gca,'Xlim',[quickticks(1) quickticks(end)],'XTick',quickticks,'XtickLabel',quickticks*1000); xlabel('time (ms)'); ylabel('Number of trials'); axis square;
%%
% Section 37: optogenetic stimulation latencies
load '~/Documents/Grants/R01EY018849/2021/Data and Code/latency_E.mat'
binedges = logspace(log10(1),log10(20),15);
bincenters = geomean([binedges(2:end);binedges(1:end-1)])
if min(latency_E) < binedges(1) || max(latency_E) > binedges(end)
    error('Need to extend bin range');
end
[n,x] = histcounts(latency_E,binedges);

figprefs; axes;
h = bar(bincenters,n,2);
set(h,'FaceColor',[.5 .5 .5],'LineWidth',1);
ylabel('# neurons');
xlabel('latency (ms)');
%%
% Section 38
% Simulation of and energy computation across luminance simple cells
% and color-luminance tuned DO cells. Does this simulated color-complex
% cell have spatial phase tuning?

% Only considering L and M cones
x = linspace(0,6*pi,100);
LUMRF1 = repmat(sin(x).*normpdf(x,mean(x),2),2,1); % L; M identical
LUMRF2 = repmat(cos(x).*normpdf(x,mean(x),2),2,1); % L; M

%LUMRF1 = LUMRF1.*repmat([1.1; .2], 1, size(LUMRF1,2)); % Unequal cone weights
%LUMRF2 = LUMRF2.*repmat([1.1; .2], 1, size(LUMRF2,2)); % Unequal cone weights

figure; axes; hold on;
plot(LUMRF1')
plot(LUMRF2')

%COLUMRF1 = LUMRF1; COLUMRF1(2,:) = 0; % L cone only 
%COLUMRF2 = LUMRF2; COLUMRF2(2,:) = 0;

COLUMRF1 = LUMRF1; COLUMRF1(1,:) = 0; % M cone only 
COLUMRF2 = LUMRF2; COLUMRF2(1,:) = 0;

%COLUMRF1 = LUMRF1; COLUMRF1(2,:) = -COLUMRF1(2,:); % L-M 
%COLUMRF2 = LUMRF2; COLUMRF2(2,:) = -COLUMRF2(2,:); % L-M 

%figure; axes; hold on;
%plot(COLUMRF1')
%plot(COLUMRF1')

phases = linspace(0,2*pi,10);
phases(end) = [];

phases = 0:pi/4:7*pi/4;
responses = [];
figure; axes; hold on;
KEEPLUMSTATIONARY = 0; % doesn't matter
for i = 1:length(phases)
    if KEEPLUMSTATIONARY
        LUMstim = cos(x);
        RGstim = cos(x+phases(i));
    else
        LUMstim = cos(x+phases(i));
        RGstim = cos(x);
    end
    
    plot(x,LUMstim,'k-');
    plot(x,RGstim,'r-');
    Lstim = LUMstim/sqrt(2)+RGstim/sqrt(2);
    Mstim = LUMstim/sqrt(2)-RGstim/sqrt(2);
    plot(x,Lstim,'r:');
    plot(x,Mstim,'g:');
    %drawnow;
    %pause;
    %cla
	LUMresp = [Lstim*LUMRF1(1,:)'+Mstim*LUMRF1(2,:)'; Lstim*LUMRF2(1,:)'+Mstim*LUMRF2(2,:)'];
    LUMresp = sqrt(sum(LUMresp.^2));
    COLUMRFresp = [Lstim*COLUMRF1(1,:)'+Mstim*COLUMRF1(2,:)'; Lstim*COLUMRF2(1,:)'+Mstim*COLUMRF2(2,:)'];
    COLUMRFresp = sqrt(sum(COLUMRFresp.^2));
    responses(i,:) = [LUMresp COLUMRFresp];
end

[phases' responses]
% Kind of intuitive. LUM always gives the same response (because it's
% insensitive to the chromatic component) and the responses of the (L-cone
% dominated) color-luminance cells is maximal when the LUM and chromatic 
% are in phase (L+M+L-M = 2L). 

% The only case where response is phase independent is when the luminance
% and chromatic drivers of the complex cell are orthogonal. 

% Under this model, if responses are maximal for phase = 0 they are minimal
% for phase = pi or vice versa. There appears to be no way to construct a
% complex cell this way that responds particularly strongly to e.g. pi/2 or
% to both 0 and pi phase.
figure; axes; hold on;
plot(phases, sqrt(sum(responses.^2,2)));
xlabel('phase (rad)');
ylabel('sim. color-complex cell response');


%% Section 39 
% 3D version of the above simulation
conecolors = [1 0 0; 0 1 0; 0 0 1]; % for plotting
PLOTSTIM = 1;
PLOTMECHS = 0;
EQUATECONTRASTS = 0;

coneweights = [1 1 1;0 0 0];
stimcc = [.5 .5 0; 0 0 .5];
if EQUATECONTRASTS
    respmat = stimcc*coneweights';
    singlestimresponses = sqrt(sum(respmat.^2,2));
    stimcc = stimcc./repmat(singlestimresponses,1,3);
    % Sanity check
    respmat = stimcc*coneweights';
    singlestimresponses = sqrt(sum(respmat.^2,2));
    singlestimresponses
end

x = linspace(0,6*pi,100);
RFtemplate = sin(x).*normpdf(x,mean(x),2); 
RFtemplate(2,:) = cos(x).*normpdf(x,mean(x),2);

RFs = zeros(size(coneweights,2),length(x), 2, size(coneweights,1)); % cones x space x sin/cos x mechanism
for i = 1:size(coneweights,1) % mechanisms
    for j = 1:size(coneweights,2) % cones
        for k = 1:size(RFtemplate,1) % sin/cos
            RFs(j,:,k,i) = RFtemplate(k,:).*coneweights(i,j);
        end
    end
end

% Plotting mechanisms
if PLOTMECHS
    figure;
    for i = 1:size(RFs,4) % mechanisms
        subplot(size(RFs,4),1,i); hold on;
        for j = 1:size(RFs,1) % cones
            for k = 1:size(RFs,3)
                plot(squeeze(RFs(j,:,k,i)),'-','Color',conecolors(j,:));
            end
        end
    end
end

% Simulated experiment

phases = linspace(0,2*pi,80);
if PLOTSTIM
    figure;
end
responses = [];
for i = 1:length(phases)
    stim = zeros(size(stimcc,2),length(x),size(stimcc,1));
    stim(:,:,1) = stimcc(1,:)'*sin(x);
    stim(:,:,2) = stimcc(2,:)'*sin(x+phases(i));
    stim = sum(stim,3);
    if PLOTSTIM
        subplot(length(phases),1,i); hold on; set(gca,'ColorOrder',conecolors);
        plot(stim');
        if i < length(phases)
            set(gca,'Xtick',[]);
        end
    end
    tmp =[];
    for j = 1:size(RFs,3) % sin/cos
        for k = 1:size(RFs,4) % mechanisms
            tmp = [tmp, sum(sum(RFs(:,:,j,k).*stim))];
        end
    end
    responses(i,:) = sqrt(sum(tmp.^2));
end
figure;
subplot(2,1,1); hold on;
plot(phases,responses);
subplot(2,1,2); hold on;
plot(responses.*cos(phases'),responses.*sin(phases'),'k-');
plot(0,0,'k*');
axis equal; 

% under no circumstances, it seems, is the preferred phase anything other
% than 0 or pi. 

% linear neurons are phase sensitive when the two stimulus components both
% act on the neuron. There's one phase where the two components combine
% constructively to drive the neuron. The neuron is inhibited to the
% opposite phase.

% Justfication for isoresponse:
% If I don't equate contrasts in some way, I'm likely not to be able to
% disprove the energy model hypothesis (because if only one component
% drives a response, the  phase of the other component won't).
% Equating contrast actually creates more phase sensitivity?

%%
% Developing intuitions for sums of sine waves of different color directions
% and different phases. Just thinking about stimuli here, not neurons.
% OK to assume orthogonality of lights? Not really.
x = linspace(0,6*pi,100);
phase = pi;
LM = 2*cos(x);
LvM = cos(x+phase);
L = LM+LvM;
M = LM-LvM;

figure; 
subplot(2,2,1); plot(x,LM,'k-');
subplot(2,2,2); plot(x,LvM,'c-');
subplot(2,2,3); plot(x,L,'r-');
subplot(2,2,4); plot(x,M,'g-');

% when two color-orthogonal gratings are in phase, you can just think of
% the vector sum.

% when L+M and L-M gratings are 90� out of phase, L and M components 
% have the same amplitude. Phase varies from 90�, when both stimulus
% components have equal contrast, to 0 (when L+M has greater contrast), to pi
% (when L-M has greater contrast).

% when L+M and L-M gratings are 180� out of phase, think of
% the vector difference.

% Think of the vector sum and difference. Changing the relative phase moves
% the vector that represents the sum (just thinking about amplitudes of L 
% and M components) smoothly from one to the other.

% Consider two gratings. Any superposition can be represented as two
% orthognal vectors in the plane (in color space?) of the two vectors. Try
% a plane. 

% For a linear neuron, can I fix the contrast of one grating and think of
% the other as adding some contrast that is proportional to the phase shift
% between the two components? I can vary the contrast of the second
% grating which should add cos(phi) effective contrast.

%%
% Section 40 
% 3D simulation. (Fixing the contrast of grating 1 and measuring firing rate
% surface (response = f(contrast, phase))for grating 2 leads to
% complicated results. Better to equate response to two lights and fix this
% ratio. In that case, for a linear neuron isoresponse contrast varies with
% the cosine of the phase between the two component gratings. No
% assumptions about angles between lights, mechanism.

coneweights = [.1 .2 0; 2 -2 0];
stimcc = [1 1 1;-.3 .3  1];
respmat = stimcc*coneweights'; % Equating component contrasts for equal responses 
singlestimresponses = sum(respmat.^2,2); % Assuming mechanisms combine with an energy calc.
stimcc = stimcc./sqrt(repmat(singlestimresponses, 1, size(stimcc,2)));

x = linspace(0,6*pi,100);
RFtemplate = sin(x).*normpdf(x,mean(x),2); 
RFtemplate(2,:) = cos(x).*normpdf(x,mean(x),2);

RFs = zeros(size(coneweights,2),length(x), 2, size(coneweights,1)); % cones x space x sin/cos x mechanism
for i = 1:size(coneweights,1) % mechanisms
    for j = 1:size(coneweights,2) % cones
        for k = 1:size(RFtemplate,1) % sin/cos
            RFs(j,:,k,i) = RFtemplate(k,:).*coneweights(i,j);
        end
    end
end

% Simulated experiment

% Varying the phase and contrast (keeping contrast ratio of components
% fixed) to obtain the firing rate surface

phases = linspace(0,2*pi,40);
contrasts = linspace(0,10,300);
%contrasts = [1 1.1];
data =[];

responses = zeros(length(phases), length(contrasts));
for i = 1:length(phases)
    for j = 1:length(contrasts)
        stim = zeros(size(stimcc,2),length(x),size(stimcc,1));
        stim(:,:,1) = contrasts(j)*stimcc(1,:)'*sin(x);
        stim(:,:,2) = contrasts(j)*stimcc(2,:)'*sin(x+phases(i));
        stim = sum(stim,3);
        
        tmp =[];
        for k = 1:size(RFs,3) % sin/cos
            for l = 1:size(RFs,4) % mechanisms
                tmp = [tmp, sum(sum(RFs(:,:,k,l).*stim))];
            end
        end
        responses(i,j) = sqrt(sum(tmp.^2));
        responses(i,j) = sqrt((abs(tmp(1))+abs(tmp(2)))^2+(abs(tmp(3))+abs(tmp(4)))^2);
        data = [data; tmp];
    end
end

% This seems to be working. On a fixed L-cone pedestal, as contrast of M
% increases, response as a function of phase is:
% first constant (second component does nothing yet)
% second cos (M adds a little when in phase with L and subtracts when
% opposite
% third abs(cos) (M is balanced with L and can reduce responses to 0 when
% 180� out of phase)
% fourth cos again (M drives response).
% 
figure;
surface(repmat(phases',1,size(responses,2)),repmat(contrasts,size(responses,1),1),responses)
set(gca,'Xtick',[0:pi/4:2*pi],'XtickLabel',[0:pi/4:2*pi]);
xlabel('phase');
ylabel('contrast');

figure; axes; hold on;
critval = prctile(responses(:),30);
dists = sum(~(responses'>critval));
L = max(responses')<critval;
polar(phases/2,dists,'ko');
polar(phases(L)/2,dists(L),'r*') % Never obtained criterion firing rate



%%
% How does an L+M neuron respond to the superposition of L and M cone
% isolating stimuli with variable phase? 
% The absolute phase of the stimulus wrt the RF makes a big difference.

x = linspace(0,6*pi,100);
RFtemplate = cos(x).*normpdf(x,mean(x),2); 
phases  = linspace(0,2*pi,43); phases(end) = [];
responses = [];
for i  = 1:length(phases) % phase difference between L and M cone-isolating gratings
    for j  = 1:length(phases) % Absolute phase of L cone grating
        stim = [sin(x+phases(i)); sin(x+phases(i)+phases(j))];
        responses(i,j) = sum(stim*RFtemplate'); % L+M
    end
end
responses
% rows are absolute phase difference, colums are relative phase difference
% between L and M.

% Response changes with relative phase (x) like abs(cos(x)) with period pi.
figure; axes; hold on;
%imagesc(responses)
plot(phases, mean(abs(responses)));
plot(phases, abs(cos(phases/2)).*max(mean(abs(responses))))
