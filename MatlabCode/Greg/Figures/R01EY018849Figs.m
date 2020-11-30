% Figures for a R01 competitive renewal
%
% Section 1) LMTF data for a monkey and a human
% Section 2) ProPixx
% Section 3) ChromCat data
% Section 4) Laminar array orientation and TF tuning curves
% Section 5) Isodetection contours in LM plane making the point that cone
% noise cannot be a major contributor
% Section 6) Sedna percent correct as a function of saccade direction.
% stimulus presentation
% Section 7) Cone model percent correct as a function of saccade direction.
% Section 8) LGN data. White noise and Isosamp.

% ------------------------------------
% Figures for a spatial-chromatic R01
% ------------------------------------
% Section 9) Rasters from a luminance complex cells and a color-luminance
% complex cell
% Section 10) Simulation of natural illuminant spectra on Munsell chip
% edge.

%
%%
% Section 1: LMTF data from a human and a monkey

lists = {'NutLMTF.txt','GregLMTF.txt'};
lists = {'ApolloLMTF.txt','NutLMTF.txt'};
data = [];
for listcounter = 1:length(lists)
    flist = flatten(fnamesFromTxt2(fullfile(nexfilepath,'nexfilelists','Greg','LMTF',char(lists{listcounter}))));
    for i = 1:length(flist)
        stro = notnex2stro(findfile(flist{i}));
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
        spdstepsize = 400/(size(stro.sum.exptParams.mon_spd,1)/3-1);
        spds = SplineSpd([380:spdstepsize:780]',reshape(stro.sum.exptParams.mon_spd,length(stro.sum.exptParams.mon_spd)/3,3),[380:5:780]');
        M = funds'*spds;
        bkgndrgb = stro.sum.exptParams.bkgndrgb;
        [in_gamut,scalar] = gamutCheck([questmodes zeros(size(questmodes,1),1)], bkgndrgb, M, 'both');
        questmodes(~in_gamut) = questmodes(~in_gamut).*scalar(~in_gamut);
        data = [data; questmodes tfs ~in_gamut' repmat(listcounter,length(tfs),1)]; % Lcc Mcc TF OOG subject
    end
end
subject = data(:,5);

% Doing the fitting
LB = [0 0 1 1 .001 .001]; % Bounds for constrained search
UB = [100 1 20 20 .03 .03]; % Bounds for constrained search
for i = 1:max(subject)
    L = logical(subject == i);
    x = data(L,1);
    y = data(L,2);
    tf = data(L,3);
    Loog = logical(data(L,4));
    initparams = [40 .1 9 3 .005 .002]; % Good for Greg
    options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-8);
    thetas = linspace(0,pi/2,12);
    toterr = [];
    fpars = [];
    fpar = [initparams, initparams];
    for j = 1:length(thetas)
        rotmat = [cos(thetas(j)) sin(thetas(j)); -sin(thetas(j)) cos(thetas(j))];
        % Rotating data counterclockwise = rotating fit axis clockwise. f1
        % should be luminance (in the first quadrant)
        xytmp = [x,y]*rotmat';
        [fpar,fv] = fmincon(@(params) tf_fiterr2(params,[xytmp(:,1) xytmp(:,2) tf],Loog),fpar,...
            [],[],[],[],[LB LB],[UB UB],[],options);
        toterr(j) = fv;
        fpars(j,:) = fpar;
    end
    bestrotidx = find(toterr == min(toterr));
    reshape(fpars(bestrotidx,:),[6 2])  % In case the user wants to see the fitted parameters
    fpar(i,:) = fpars(bestrotidx,:);
end


% Sensitivity plot
sqrtn = 40;
[x1,x2] = meshgrid(logspace(log10(1),log10(25),sqrtn),linspace(0, pi,sqrtn));
x_star = [x1(:) x2(:)]; % column order: TF, theta
tfticks = [1 2 5 10 20]

figure('Position',[267 70 482 716]);
for i = 1:max(subject)
    
    L = logical(subject == i);
    x = data(L,1);
    y = data(L,2);
    [th,r] = cart2pol(x,y);
    tf = data(L,3);
    Loog = logical(data(L,4));
    f1 = @(omega)fpar(i,1)*abs(((1i*2*pi*fpar(i,5).*omega+1).^-fpar(i,3))-fpar(i,2)*((1i*2*pi*fpar(i,6).*omega+1).^-fpar(i,4)));
    f2 = @(omega)fpar(i,1+6)*abs(((1i*2*pi*fpar(i,5+6).*omega+1).^-fpar(i,3+6))-fpar(i,2+6)*((1i*2*pi*fpar(i,6+6).*omega+1).^-fpar(i,4+6)));
    
    a = abs(f1(x1)).^-1;
    b = abs(f2(x1)).^-1;
    thtmp = x2-thetas(bestrotidx);
    rtmp = (a.*b)./sqrt((a.*cos(thtmp)).^2+(b.*sin(thtmp)).^2); % radius of ellipse - thank you, Wikipedia
    rtmp(rtmp > 1) = nan;
    
    subplot(2,1,i); hold on;
    h = surf(log10(x1),x2,log10(1./rtmp));
    set(h,'Edgecolor','none'); % camlight headlight; lighting phong
    alpha(.7)

    plot3(log10(tf(~Loog)),mod(th(~Loog),pi),log10(1./r(~Loog)),'ko','MarkerSize',10,'MarkerFaceColor','black');
    plot3(log10(tf(Loog)),mod(th(Loog),pi),log10(1./r(Loog)),'ro','MarkerSize',10,'MarkerFaceColor','red');

    set(gca,'YTick',[0 pi/4 pi/2 3*pi/4],'YTickLabel',{'L','L+M','M','L-M'});
    colormap(hot);
    xlabel('TF (Hz)');
    ylabel('Color direction');
    zlabel('Log sensitivity');
    
  %  set(gca,'Xscale','log','Zscale','log');
    set(gca,'Xlim',[min(log10(tf)) max(log10(tf))],'Ylim',[0 pi],'Zlim',log10([1 50]));
    set(gca','View',[140 36]);
    set(gca,'XTick',log10(tfticks),'XtickLabel',tfticks);
    set(gca,'Color',[.5 .5 .5]);
end
set(gcf,'Renderer','OpenGL','Resize','off','PaperPosition',[1 1 5 10]);

%%
% Section 2
load('/Users/greghorwitz/Desktop/MatlabCode/Slave/Monitor Calibration/Monitor data/ProPixx/ProPixx.mat')
cal = cals{end};
load('/Users/greghorwitz/Desktop/MatlabCode/Slave/Monitor Calibration/Monitor data/ProPixx/20140416T123338-TestCalOut.mat')
% do this in case someone executes this analysis code manually using old data (before June 2013)
% where RGBs are 8 bit
if find(RGBs > 1, 1)
    if ~exist('S_from', 'var')
        S_from = [380 400/(size(spdMat, 2)-1) size(spdMat, 2)];
    end
    RGBs = RGBs / 255;
end

% spline all spectra to the sampling lattice of the most recent calibration
S_to = cals{end}.S_device;

maxY = -inf;
minY = inf;

cal = cals{1};
rgbbkgnd = FindModelWeights(cal.P_ambient, cal.P_device);

% get the expected intensities that would result from the input voltages in RGBs
GT_idxs = interp1(cal.gammaInput, (1:size(cal.gammaInput,1))', RGBs(:), 'nearest');
GT_idxs = reshape(GT_idxs, size(RGBs));
GT_idxs = bsxfun(@plus, GT_idxs, 0:size(cal.gammaTable,1):numel(cal.gammaTable)-1);
predicted = bsxfun(@plus, cal.gammaTable(GT_idxs), rgbbkgnd');

% "actual" is a bit of a misnomer since we're finding the coefficients on each of the guns
% by regression.
actual = FindModelWeights(SplineSpd(S_from, spdMat', S_to), ...
    SplineSpd(cal.S_device, cal.P_device, S_to))';

figure(calIdx+prevexistingfigh);
set(gcf, 'DefaultAxesColorOrder', eye(3));
plot(RGBs, actual-predicted, '.', 'markersize', 10);
title(sprintf('%s: Calibration #%d', cals{1}.describe.monitor, calIdx));
ylabel('Actual - Predicted');
xlabel('Requested');

ylims = get(gca, 'YLim');
minY = min([minY ylims(1)]);
maxY = max([maxY ylims(2)]);

%%
% Section 3
% ChromCat data

flist = flatten(fnamesFromTxt2(fullfile(nexfilepath,'nexfilelists','Greg','TPM','SednaChromCat.txt')));
INCLUDEENDPOINTS = 0;

tmpdata = [];
for i = 1:length(flist)
    stro = nex2stro(findfile(flist{i}));
    rightwardsaccade = stro.trial(:,strcmp('rightward_sacc',stro.sum.trialFields(1,:)));
    Lcc = stro.trial(:,strcmp('lcc',stro.sum.trialFields(1,:)));
    Mcc = stro.trial(:,strcmp('mcc',stro.sum.trialFields(1,:)));
    thetas = atan2(Lcc,Mcc);
    uniquethetas = unique(thetas);
    n = zeros(size(uniquethetas,1),1);
    k = zeros(size(uniquethetas,1),1);
    for j = 1:size(uniquethetas,1)
        L = thetas == uniquethetas(j);
        n(j) = sum(L);
        k(j) = sum(L&rightwardsaccade);
    end
    tmpdata = [tmpdata; uniquethetas k n];
end
uniquethetas = unique(thetas);
data = [];
for i = 1:length(uniquethetas)
   L = tmpdata(:,1) == uniquethetas(i);
   data(i,:) = [uniquethetas(i) sum(tmpdata(L,2)) sum(tmpdata(L,3))]
end

if (~INCLUDEENDPOINTS)
   data(1,:) = [];
   data(end,:) = []; 
end

figure; axes; hold on;
phat = data(:,2)./data(:,3);
plot(data(:,1),phat,'ko','MarkerSize',5);
se = sqrt((phat.*(1-phat))./data(:,3))
plot([data(:,1), data(:,1)]', [phat+se, phat-se]','k-')
[alpha,beta,thresh50]=FitWeibYN(data(:,1),data(:,2),data(:,3)-data(:,2),1.6,2)
x = linspace(min(data(:,1)),max(data(:,1)),30);
plot(x,1-exp(-(x/alpha).^beta),'k-')
set(gca,'Ylim',[0 1])

%%
% 3.1 ChromCat figure
listnames = {'fnames_1Hz','fnames_3Hz','fnames_6Hz'};
%fnames_1Hz = {'S102414001','S102414002','S102414003'};
fnames_1Hz = {'S102414001','S102414002','S102414003', 'S102714001','S102714002','S102814001','S102814002','S102814003'};
fnames_3Hz = {'S101714001', 'S101714002'};
fnames_6Hz = {'S102214001','S102214002','S102214003','S102214004'};

data = [];
for listidx = 1:length(listnames)
    filenames = eval(listnames{listidx})
    tmpdata = [];
    for i = 1:length(filenames)
        stro = nex2stro(findfile(filenames{i}));
        rightwardsaccade = stro.trial(:,strcmp('rightward_sacc',stro.sum.trialFields(1,:)));
        Lcc = stro.trial(:,strcmp('lcc',stro.sum.trialFields(1,:)));
        Mcc = stro.trial(:,strcmp('mcc',stro.sum.trialFields(1,:)));
        thetas = atan2(Mcc,Lcc);
        uniquethetas = unique(thetas);
        n = zeros(size(uniquethetas,1),1);
        k = zeros(size(uniquethetas,1),1);
        for j = 1:size(uniquethetas,1)
            L = thetas == uniquethetas(j);
            n(j) = sum(L);
            k(j) = sum(L&rightwardsaccade);
        end
        tmpdata = [tmpdata; uniquethetas k n];
    end
    uniquethetas = unique(tmpdata(:,1));
    for i = 1:length(uniquethetas)
        L = tmpdata(:,1) == uniquethetas(i);
        data = [data; listidx uniquethetas(i) sum(tmpdata(L,2)) sum(tmpdata(L,3))];
    end
end

% Doing the plotting
TFIDX = 1;
THETAIDX = 2;
KIDX = 3;
NIDX = 4;
figure; axes; hold on;
colors = {'red','blue','black'}; % 1 Hz, 3 Hz, then 6 Hz
for listidx = 1:length(listnames)
    L = data(:,TFIDX) == listidx;
    % Getting rid of the end points
    L(data(:,THETAIDX) == -pi/4) = 0;
    L(data(:,THETAIDX) == pi/4) = 0;
    
    phat = (data(L,NIDX)-data(L,KIDX))./data(L,NIDX); % Proportion *left* (achromatic) decisions
    plot(data(L,THETAIDX),phat,'o','MarkerSize',8,'MarkerFaceColor',colors{listidx},'MarkerEdgeColor',colors{listidx});
    se = sqrt((phat.*(1-phat))./data(L,NIDX))
    plot([data(L,THETAIDX), data(L,THETAIDX)]', [phat+se, phat-se]','-','color',colors{listidx})
    [alpha,beta,thresh50]=FitWeibYN(data(L,THETAIDX)+pi/4,data(L,NIDX)-data(L,KIDX),data(L,KIDX),pi/4,2)
    x = linspace(min(data(L,THETAIDX)),max(data(L,THETAIDX)),30);
    plot(x,1-exp(-(((x+pi/4)/alpha)).^beta),'-','color',colors{listidx},'linewidth',2)
end
set(gca,'Ylim',[0 1]);
set(gca,'Xlim',[-pi/4 pi/4]);
xticklabels = [-pi/4:pi/8:pi/4];
set(gca,'XTick',xticklabels,'XTickLabel',{'-pi/4','-pi/8','0','pi/8','pi/4'});
ylabel('Proportion achromatic choices')


%%
% Section 4)
% Laminar probe stuff. LFPs, orientation tuning curves, and TF tuning
% curves

PRINT = 1;
stro = nex2stro(findfile('p032015001'));

if strcmp(stro.sum.fileName(end-13:end-4),'p031915008'); % Reordering columns 
   stro.ras = stro.ras(:,[9:16, 1:8, 25:32, 17:24 ,33]);
   stro.sum.rasterCells = stro.sum.rasterCells([9:16, 1:8, 25:32, 17:24 ,33])
end

stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
fpacq_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_acq'));
fpon_t= stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_on'));

analogstarttime = [stro.ras{:,strcmp( stro.sum.rasterCells,'anlgStartTime')}]';
arate = stro.sum.analog.storeRates{1};
analogL = strncmp(stro.sum.rasterCells,'AD',2);
% getting rid of eye position and channel 1
analogL(strcmp(stro.sum.rasterCells,'AD11')) = 0;
analogL(strcmp(stro.sum.rasterCells,'AD12')) = 0;
analogL(strcmp(stro.sum.rasterCells,'AD17')) = 0;  % Channel 1 is obviously bad

analogidxs = find(analogL);
ntrials = size(stro.trial,1);

ntsamps = 100;
t_offset = -.01;
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
times = t(Lt);
b1 = 1/arate;
b0 = -t_offset;
xtickpos = [-.1 0 .05 .1 .15];
set(gca,'XTick',(xtickpos+b0)./b1,'XTicklabel',xtickpos*1000);
set(gca,'YTick',[1:2:15],'YTickLabel',[1:2:15]*150);
ylabel('Depth (microns)'); xlabel('Time (ms)');
if (PRINT)
    print -dtiff LFP
end

% ---------------------------------------------------------
% Now plotting orientation tuning curves (At the best TF? Across all TFs?)
orients = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'orient'));
sfs = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sf'));
tfs = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
unq_sfs = unique(sfs');
unq_tfs = unique(tfs');
if (length(unq_sfs) > length(unq_tfs))
    fs = sfs;
    unq_fs = unq_sfs;
else
    fs = tfs;
    unq_fs = unq_tfs;
end

unq_orients = unique(orients');
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
stimoff_t= stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
spikeidxs = strncmp(stro.sum.rasterCells,'sig',3);
spikerates = [];
baselines = [];  baseline_t = 0.5;
for spikeidx = find(spikeidxs)
    for i = 1:size(stro.trial,1)
        spiketimes = stro.ras{i,spikeidx};
        nspikes = sum(spiketimes > stimon_t(i) & spiketimes < stimoff_t(i));
        spikerates(i,spikeidx) = nspikes./(stimoff_t(i)-stimon_t(i));
        nspikes = sum(spiketimes > stimon_t(i)-baseline_t & spiketimes < stimon_t(i));
        baselines(i,spikeidx) = nspikes./baseline_t;
    end
end

% Plotting
for spikeidx = find(spikeidxs)
    figure('Position',[219 525 216 420]);
    % Orientation
    data = [];
    for theta = unq_orients
        L = theta == orients;
        data = [data; mean(spikerates(L,spikeidx)) std(spikerates(L,spikeidx)) sum(L)];
    end    
    data = [data(end,:); data([1:end-1],:); ];
    data(end+1,:) = data(1,:);
    or = [unq_orients(end)-2*pi, unq_orients(1:end-1) unq_orients(end)]*180/pi;
    
    subplot(2,1,1); hold on;
    errorbar(or, data(:,1),data(:,2)./sqrt(data(:,3)),'LineWidth',2,'Color','black');
    plot([or(1) or(end)],[mean(baselines(:,spikeidx)) mean(baselines(:,spikeidx))],'k--')
    set(gca,'Xlim',[or(1)-22.5 or(end)+22.5]);
    set(gca,'XTick',[0 90 180 270]);
    ylabel('Spikes/sec');
    xlabel('Orientation (deg)');
    h = title(stro.sum.rasterCells{spikeidx},'FontSize',13);
    yl = get(gca,'Ylim');
    set(gca,'Ylim',[0 yl(2)]);
    
    % Temporal frequency
    data = [];
    for f = unq_fs
        L = fs == f;
        data = [data; mean(spikerates(L,spikeidx)) std(spikerates(L,spikeidx)) sum(L)];
    end    
    subplot(2,1,2); hold on;
    errorbar(unq_fs, data(:,1),data(:,2)./sqrt(data(:,3)),'LineWidth',2,'Color','black');
    plot([unq_fs(1) unq_fs(end)],[mean(baselines(:,spikeidx)) mean(baselines(:,spikeidx))],'k--')
    set(gca,'XScale','log')
    set(gca,'Xlim',[unq_fs(1)*.9 unq_fs(end)*1.1]);
    set(gca,'XTick',[1 2.2 5 11 25]);
    ylabel('Spikes/sec');
    xlabel('Temporal frequency (Hz)');
    set(gca,'Ylim',[0 yl(2)]);
    
    if (PRINT)
        eval(['print -depsc OriTF',num2str(spikeidx)]);
    end
end

%%
% Section 5
% Detection thresholds at low temporal frequency in the LM plane for two
% monkey observers and the cone model ideal observer.

filelists = {'FreyaLMTF.txt','SednaLMTF.txt'};

startpath = '/Volumes/NO BACKUP/NexFiles/Greg/Sedna';
data = [];

for j = 1:length(filelists)
    flist = flatten(fnamesFromTxt2(fullfile(nexfilepath,'nexfilelists','Greg','LMTF',char(filelists{j}))));
    % Iterate over the list of .nex files.
    for i = 1:length(flist)
        stro = notnex2stro(findfile(flist{i},startpath)); % <-- process the information in the nex file and put it into the "stro" structure
        if ~(abs(stro.sum.exptParams.stim_x) == 50 & abs(stro.sum.exptParams.stim_y) <= 10)
            [stro.sum.exptParams.stim_x stro.sum.exptParams.stim_y]
            continue;
        end

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

% Loading the detection thresholds Zack obtained from Charlie's model
% (at 5° eccentricity and 3 Hz)
load '/Users/greghorwitz/Documents/Grants/Horwitz R01 competitive renewal/A1/LM_model_thresholds.mat';
modeldata = [repmat(3,size(color_dirs,1),1) color_dirs(:,[1 2]).*repmat(alphas,1,2) repmat([1 0],size(color_dirs,1),1) ];
data = [data; modeldata];

% Analysis and plotting
TFthreshs = [1 1.2];
Loog = logical(data(:,end));
subjectidxs = unique(data(:,1))
figure; axes; hold on;
symbols = {'bs','rd','kv'};
params = [];
hsym = []; hline = [];
for j = subjectidxs'
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
    
    hsym = plot([x; -x],[y; -y],symbols{j});
    markercolor = get(hsym,'Color');
    set(hsym,'MarkerFaceColor',markercolor);
    if (j == 3)
        delete(hsym);
    end
    hline(j) = plot(tmpx,tmpy,'m-','LineWidth',2);
    set(hline(j),'Color',markercolor);
    params = [params; fpar];
end
axis square;
set(gca,'Xlim',[-.12 .12],'Ylim',[-.12 .12]);
legend(hline,{'Monkey 1','Monkey 2','Cone model'},'location','southeast');

%%
% Section 6
% Percent correct as a function of microsaccade direction.
% Only looking at saccades during the stimulus presentation.
% Get rid of Loogs.

filenames = fnamesFromTxt2('/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/LMTF/ApolloLMTF.txt');
amplitudes = [];
directions = [];
peakv = [];
pathlengths = [];
durations = [];
trialparams = [];
sacduringstim = [];
trialcounter = 1;
for a = 1:size(filenames,1)
    stro = nex2stro(findfile(char(filenames{a})));
    if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
        continue;
    end
    ntrials = size(stro.trial,1);
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
    tf = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
    lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'lcc'));
    mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'mcc'));
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    oog = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'oog'));

    sacstats = getSacData(stro);
    close;
    
    for i = 1:ntrials
        st = sacstats.starttimes{i};
        Lsac = (sacstats.starttimes{i} > stimon_t(i)) & (sacstats.endtimes{i} < stimoff_t(i));
        if any(Lsac)
            amplitudes = [amplitudes; sacstats.amplitudes{i}(Lsac)];
            directions = [directions; sacstats.directions{i}(Lsac)];
            peakv = [peakv; sacstats.peakv{i}(Lsac)];
            durations = [durations; sacstats.durations{i}(Lsac)];
            trialparams = [trialparams; repmat(trialcounter,sum(Lsac),1) repmat(correct(i),sum(Lsac),1) repmat(oog(i),sum(Lsac),1) lcc(sacstats.trialnums{i}(Lsac)) mcc(sacstats.trialnums{i}(Lsac)) tf(sacstats.trialnums{i}(Lsac))];
        end
        trialcounter = trialcounter + 1;
    end
end

[length(unique(trialparams(:,1))) trialcounter]; % Number of trials with microsaccdes during stimulus presentation

% First, % correct for trials with up and down microsaccades
% (trials without microsaccades aren't available).
% Columns of data:
% Correct, oog, mean angle, mean amplitude, tf
trialidxs = unique(trialparams(:,1));
data = nan*ones(length(trialidxs),5);
for i = 1:length(trialidxs)
    L = trialparams(:,1) == trialidxs(i);
    [x,y] = pol2cart(directions(L),amplitudes(L));
    [th,r] = cart2pol(mean(x),mean(y));
    data(i,:) = [mean(trialparams(L,2)), mean(trialparams(L,3)), th, r, mean(trialparams(L,6))]; 
end
Loog = data(:,2) == 1;
olddata = data;
data(Loog,:) = []; % Getting rid of all the out of gamut points
LTF = data(:,5) >= 15; % Needs to be "fast" (i.e. TF > x)
Lup = data(:,3) > 0 & data(:,3) < pi;

% Barplots
ns = [sum(Lup&LTF) sum(~Lup&LTF) sum(Lup&~LTF) sum(~Lup&~LTF)]; % upfast, downfast, upslow, downslow
ps = [sum(Lup&LTF&data(:,1))./ns(1),...
    sum(~Lup&LTF&data(:,1))./ns(2),...
    sum(Lup&~LTF&data(:,1))./ns(3),...
    sum(~Lup&~LTF&data(:,1))./ns(4)];
ses = sqrt((ps.*(1-ps))./ns);
figure; axes; hold on;
h = bar(ps,'k');
plot([1;1]*[1:4],[ps;ps]+[zeros(1,length(ns));ses],'k-','LineWidth',2);
plot([1;1]*[1:4],[ps;ps]-[zeros(1,length(ns));ses],'w-','LineWidth',2);
set(gca,'Ylim',[.5 .95]);

% mean(data(Lup,1))
% mean(data(~Lup,1))
% % Some statistics
% [sum(Lup&data(:,1)), sum(~Lup&data(:,1));...
%  sum(Lup&~data(:,1)), sum(~Lup&~data(:,1))]
% % Looks like there are lots of "saccade down & incorrect" trials
% [orat,~,~,p] = oddsratio(sum(Lup&data(:,1)), sum(~Lup&data(:,1)),...
%                     sum(Lup&~data(:,1)), sum(~Lup&~data(:,1)),0.05)
% 
nbins=8;
bins = linspace(0,2*pi,nbins+1);
bins(end) = [];

% Histogramming data
[~,r_fast] = rose(data(LTF,3),bins);
counts_fast = r_fast([3:4:end]);
[~,r_corfast] = rose(data(data(:,1) == 1&LTF,3),bins);
counts_corfast = r_corfast([3:4:end]);

[~,r_slow] = rose(data(~LTF,3),bins);
counts_slow = r_slow([3:4:end]);
[~,r_corslow] = rose(data(data(:,1) == 1&~LTF,3),bins);
counts_corslow = r_corslow([3:4:end]);

% Standard errors
p = counts_corfast./counts_fast;
n = counts_fast;
sefast = sqrt(p.*(1-p))./sqrt(n);
p = counts_corslow./counts_slow;
n = counts_slow;
seslow = sqrt(p.*(1-p))./sqrt(n);

figure; 
polar([bins bins(1)],[counts_corslow counts_corslow(1)]./[counts_slow counts_slow(1)],'k.-');
hold on;
h = polar([bins bins(1)],[counts_corfast counts_corfast(1)]./[counts_fast counts_fast(1)],'r.-');

for i = 1:length(bins)
    polar([bins(i) bins(i)],(counts_corfast(i)/counts_fast(i))+[-sefast(i) sefast(i)],'r-')
end
for i = 1:length(bins)
    polar([bins(i) bins(i)],(counts_corslow(i)/counts_slow(i))+[-seslow(i) seslow(i)],'k-')
end

% Little inset showing distribution of saccade directions
[t_all,r_all] = rose(data(:,3),15);
t_all = t_all([3:4:end]);
r_all = r_all([3:4:end]);
polar(t_all, r_all./(2*max(r_all)),'k-');
    
% ps = nan*ones(length(bins),1);
% for i = 1:length(bins)
%     [~,p] = equalproptest([counts_corfast(i); counts_corslow(i)],[counts_fast(i); counts_slow(i)],0.05);
%     ps(i) = p;
% end
% Add errorbars?


%%
% Section 7 
% Cone model percent correct as a function of saccade direction (just up vs
% down for now).

%load('/Users/greghorwitz/Documents/Grants/Horwitz R01 competitive renewal/Matlabstuff/LM_eyemovement_thresh_compare.mat'); % 1.1 Hz
load('/Users/greghorwitz/Documents/Grants/Horwitz R01 competitive renewal/Matlabstuff/LM_eyemovement_thresh_compare_20Hz.mat'); % 20 Hz

% emneg.cones.alpha_analytic
% emneg.cones.beta_analytic
% empos.cones.alpha_analytic
% empos.cones.beta_analytic
% 
% for i = 1:size(emneg.gab.contrasts,2)
%     figure; axes; hold on;
%     x = emneg.gab.contrasts{i};
%     y = emneg.idlob.roc_analytic{i};
%     alpha = emneg.cones.alpha_analytic(i);
%     beta = emneg.cones.beta_analytic(i);
%     % Weibull 1-exp-(x/alpha).^beta
%     
%     plot(x,.5*(1-exp(-x/alpha).^beta)+.5)
%     plot(x,y,'k.');
%     set(gca,'Xscale','log');
% end

% Forget about these Weibull fits (which are terrible) just interpolate ROC
% areas.
% Positive numbers mean saccade is *in* the direction of motion (slowing it
% down)
xneg = emneg.gab.contrasts{1};
yneg = emneg.idlob.roc_analytic{1};
xpos = empos.gab.contrasts{1};
ypos = empos.idlob.roc_analytic{1};
figure; axes; hold on;
plot(xneg,yneg,'k.')  % less sensitive
plot(xpos,ypos,'b.')  % more sensitive
set(gca,'Xscale','log');

% Need to chuck out the largest few contrasts because the flattness of
% the psychometric function at saturation is screwing up interp1.
nominalthreshold = (1-exp(-1))/2+.5; % .816
nominalthreshold = .75;
threshold = interp1(ypos(1:20),xpos(1:20),nominalthreshold); % threshold
% Now getting percent correct at that contrast level
pcrt_cor = interp1(xneg(1:20),yneg(1:20),threshold)


%%
% Section 8
% White noise and Isosamp data for an LGN neuron

WN = nex2stro(findfile('A052015003.nex'));
Iso = nex2stro(findfile('A052015005.nex'));

% White noise first
frametime = 1000/WN.sum.exptParams.framerate; % ms
nstixperside = WN.sum.exptParams.nstixperside;
spikename = getSpikenum(WN);
spikeidx = find(strcmp(WN.sum.rasterCells(1,:),spikename));
maxT = 9;
out = getWhtnsStats(WN,maxT,'STAmex', {nstixperside^2, 3, maxT}, spikename);
muvect = reshape(repmat([.5 .5 .5],nstixperside^2,1),nstixperside^2*3,1);
STAs = out{1}./(2*max(abs(out{1}(:)))+eps);
% Plotting WN
figure;
for i = 1:size(STAs,2)
    STA = reshape(STAs(:,i)+muvect,[nstixperside nstixperside 3]);
    subplot(6,size(STAs,2),size(STAs,2)-i+1);
    image(STA);
    set(gca,'XTick',[],'YTick',[]); axis square;
   % title(num2str([i i-1]*frametime,'%6.1f'));
    title(num2str((i-0.5)*frametime,'%6.1f'));
end

% Plotting IsoSamp data
l = Iso.trial(:,strcmp(Iso.sum.trialFields(1,:),'stim_l'));
m = Iso.trial(:,strcmp(Iso.sum.trialFields(1,:),'stim_m'));
tf = Iso.trial(:,strcmp(Iso.sum.trialFields(1,:),'tf'));
stimon_t = Iso.trial(:,strcmp(Iso.sum.trialFields(1,:),'stimon_t'));
stimoff_t = Iso.trial(:,strcmp(Iso.sum.trialFields(1,:),'stimoff_t'));
dur = mean(stimoff_t-stimon_t);
% Counting stimuli and trials
uniquestim = unique([l m tf],'rows');
uniquestim = sortrows(uniquestim,3); % sorting by TF
offset = [dur/3 (2/3)*dur]; 
for i = 1:size(uniquestim,1)
    current_tf = uniquestim(i,3);
    L = l == uniquestim(i,1) & m == uniquestim(i,2) & tf == current_tf;
    tmpspikes = [];
    for j = find(L)'
        tmpspikes = [tmpspikes; Iso.ras{j,spikeidx}-stimon_t(j)];
    end
    tmpspikes(tmpspikes < offset(1) | tmpspikes > offset(2)) = [];
    cos_n_sin = [cos(2*pi*current_tf.*tmpspikes) sin(2*pi*current_tf.*tmpspikes)];
    F1 = sqrt(sum(sum(cos_n_sin).^2))./sum(L)/(offset(2)-offset(1)); % spikes/(trial*sec)
    data(i) = F1;
end

normF1 = data-nanmin(data);
normF1 = normF1./nanmax(normF1);
cjet = jet(255);

figure; axes; hold on;
for j = 1:size(uniquestim,1)
    L = l == uniquestim(j,1) & m == uniquestim(j,2) & tf == uniquestim(j,3);
    if (~isnan(normF1(j)))
        h = plot3(uniquestim(j,1),uniquestim(j,2),uniquestim(j,3),'ko','MarkerSize',10,'MarkerEdgeColor','none');
        set(h,'MarkerFaceColor',[normF1(j) 0 1-normF1(j)]);
    end
end
xlabel('L-cone contrast');
ylabel('M-cone contrast');
zlabel('TF (Hz)');
set(gca,'View',[105 16]);
axis vis3d;
[max(data) min(data)]

%%
% A second IsoSamp LMTF data plot with rasters (A053015005)
%stro = nex2stro(findfile('A053015005.nex'));
stro = nex2stro(findfile('A060115007.nex'));

l = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_l'));
m = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_m'));
tf = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
dur = mean(stimoff_t-stimon_t);
Lchrom_trials = sign(l) ~= sign(m);

% Counting stimuli and trials
uniquestim = unique([l m tf],'rows');
uniquestim = sortrows(uniquestim,3); % sorting by TF

% Plotting rasters
Lspikechans = strncmp(stro.sum.rasterCells,'sig0',4);
offset = [-.5 .5];  % pre and post time wrt stimon/stimoff

figure; 
for i = 1:2
    subplot(2,2,i); hold on; counter = 0;
    spikes = stro.ras(:,1);
    for j = 1:size(uniquestim,1)
        L = l == uniquestim(j,1) & m == uniquestim(j,2) & tf == uniquestim(j,3);
        if i == 1
            L = L&Lchrom_trials;
        else
            L = L&~Lchrom_trials;
        end
        for k = find(L)'
            tmpspikes = spikes{k}-stimon_t(k);
            tmpspikes(tmpspikes < offset(1) | tmpspikes > dur+offset(2)) = [];
            nspikestot = length(tmpspikes);
            plot([tmpspikes tmpspikes]',[zeros(nspikestot,1) 1*ones(nspikestot,1)]'+counter,'k-','linewidth',1);
            counter = counter + 1;
        end
    end
    set(gca,'Xlim',[0 dur],'Ylim',[0 counter],'Ytick',[],'Xtick',[0 .2 .4 .6 .8],'Box','on');
end

% Now the second panel showing the stimulus
bins = linspace(0,dur,100);
temporalenvelope = ones(size(bins));
temporalenvelope(1:round((1/3)*length(bins))) = linspace(0,1,round((1/3)*length(bins)));
temporalenvelope(end:-1:round((2/3)*length(bins))+1) = linspace(0,1,round((1/3)*length(bins)));
contrast = sqrt(l.^2+m.^2);
stimtimecourse = zeros(length(Lchrom_trials), length(bins));

for j = 1:length(Lchrom_trials)
    stimtimecourse(j,:) = contrast(j).*temporalenvelope.*sin(2*pi*tf(j)*bins);
end

%stimtimecourse = sign(stimtimecourse).*stimtimecourse.^2; % non-linearity to augment contrast differences
im = stimtimecourse./(max(abs(stimtimecourse(:)))*2);
im = im*255+127;
for i = 1:2
    if i == 1
        L = Lchrom_trials;
    else
        L = ~Lchrom_trials;        
    end
    imtmp = im(L,:);
    [~,idx] = sort(tf(L));
     subplot(2,2,2+i);

    image(flipud(imtmp(idx,:)));
    xtick = 0:.2:dur;
    set(gca,'XTick',interp1([bins(1) bins(end)],[1 length(bins)],xtick));
    set(gca,'Xticklabel',xtick);
    ytick = round(min(tf(L))):5:round(max(tf(L)));
    set(gca,'YTick',interp1([ytick(1) ytick(end)],[1 sum(L)],ytick));
    set(gca,'Yticklabel',fliplr(ytick));
    colormap(jet(255));
end

% Assuming identical noise for luminance and color, it is sufficient here
% to compare luminance to chromatic *signal* and thereby compare there
% signal to noise ratios.

% Getting a SNR per condition
for i = 1:2 % 1 = chrom, 2 = lum
    for j = 1:size(uniquestim,1)
        L = l == uniquestim(j,1) & m == uniquestim(j,2) & tf == uniquestim(j,3);
        if i == 1
            L = L&Lchrom_trials;
        else
            L = L&~Lchrom_trials;
        end
        for k = find(L)'
            tmpspikes = spikes{k}-stimon_t(k);
            tmpspikes(tmpspikes < offset(1) | tmpspikes > dur+offset(2)) = [];
            nspikestot = length(tmpspikes);
            plot([tmpspikes tmpspikes]',[zeros(nspikestot,1) 1*ones(nspikestot,1)]'+counter,'k-','linewidth',1);
            counter = counter + 1;
        end
    end
end

%%
% Section 9
% Rasters for a luminance-only complex cell and a color-luminance complex cell
% Gratings files
lum_filename = 'S020310001'; % <--luminance only
pan_filename = 'S041510002';
label = {'L+M','L-cone','M-cone','L-M'};
filenames = {lum_filename,pan_filename}
figure;
for fileidx = 1:2
    filename = filenames{fileidx};
    stro = nex2stro(findfile(filename));
    protocols = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'protocol'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoff_t= stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
    spikeidxs = 1;
    Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'lcont'));
    Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'mcont'));
    Scc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'scont'));
    colordirections = [Lcc Mcc Scc];
    uniquecolordirs = [.09 .09 0; 0.1273 0 0; 0 0.1273  0; .09 -.09 0]
   % uniquecolordirs = [.09 .09 0; .09 -.09 0]
    
    Lprotocol = protocols == 4;
    for i = 1:size(uniquecolordirs,1)
        L = Lprotocol &...
            colordirections(:,1) == uniquecolordirs(i,1) & ...
            colordirections(:,2) == uniquecolordirs(i,2) & ...
            colordirections(:,3) == uniquecolordirs(i,3);
        trlidxs = find(L);
        subplot(size(uniquecolordirs,1),2,2*i-1+(fileidx-1)); hold on;
        for j = trlidxs'
            plotidx = find(j==trlidxs);
            sp = stro.ras{j,1}-stimon_t(j);
            plot([sp sp]',[plotidx*ones(length(sp),1) plotidx*ones(length(sp),1)+.5]','k-');
            plot([stimoff_t(j)-stimon_t(j) stimoff_t(j)-stimon_t(j)],[plotidx-.5 plotidx+1],'m-','Linewidth',1);
            plot([0 0],[plotidx-.5 plotidx+1],'m-','Linewidth',1);
        end
        set(gca,'YLim',[.5 plotidx+1],'YTick',[],'XTick',[],'XLim',[-.1 1.1]);
        title(label{i});
    end
    set(gca,'XTick',[0 .5 1]);
    xlabel('Time (s)');
end

%%
% Section 10
% Some kind of figure showing weighted cone-excitation
% differences when Munsell chips are illuminated with natural illuminants.
% Showing a difference between L-M-S and L-M+S weighting functions.
% Implications for double-opponent cells?

load T_cones_smj10
load sur_nickerson.mat
load B_cieday
munsell = sur_nickerson;
wls = [380:5:780];
niter = 100;
orange_cyan_color = [0 .8 .8];
lime_magenta_color = [1 .7 1];
% Coefficients for the daylight spectra
% (from Judd et al. 1964 table 3)
daylight_coefficients = [-1.14 .677;
    -0.784 -0.195;
    -0.293 -0.698;
    0.145 -0.752;
    1.005 -0.378];
nillums = 30;
coeffs = linspace(daylight_coefficients(1,1),daylight_coefficients(end,1),nillums);
daylight_coefficients_interp = [coeffs', interp1(daylight_coefficients(:,1),...
    daylight_coefficients(:,2),...
    coeffs,'spline')'];

illuminants = repmat(B_cieday(:,1),1,nillums)+B_cieday(:,[2 3])*daylight_coefficients_interp';

figure; axes; hold on; set(gca,'TickDir','out'); axis square; set(gca,'Color',[.5 .5 .5]);
for iter = 1:niter
    whichchips(iter,:) = [1 1];
    while whichchips(iter,1) == whichchips(iter,2) % Avoiding same chip on boths sides
        whichchips(iter,:) = unidrnd(size(munsell,2),2,1)';
    end
    reflectance1 = munsell(:,whichchips(iter,1));
    reflectance2 = munsell(:,whichchips(iter,2));
    light1 = illuminants.*reflectance1;
    light2 = illuminants.*reflectance2;
    light0 = repmat(mean(illuminants,2),1,nillums); % a single mean illuminant
    
    lms_1 = T_cones_smj10*light1;
    lms_2 = T_cones_smj10*light2;
    lms_0 = T_cones_smj10*light0;
    
    coneexcitations = cat(3, lms_1',lms_2',lms_0'); % rows = illums, cols = LMS, planes = light1, light2, adaptating light
    coneexcitations = permute(coneexcitations,[2 3 1]);
    blue_orange_signals = []; lime_magenta_signals = [];
    for i = 1:size(illuminants,2)
        l = squeeze(coneexcitations(1,:,i));
        m = squeeze(coneexcitations(2,:,i));
        s = squeeze(coneexcitations(3,:,i));
        l_norm = (l-l(3))./l(3);
        m_norm = (m-m(3))./m(3);
        s_norm = (s-s(3))./s(3);  
        blue_orange_signals = [blue_orange_signals;(s_norm(1)+m_norm(1))/2-l_norm(1) (s_norm(2)+m_norm(2))/2-l_norm(2)];
        lime_magenta_signals = [lime_magenta_signals;(s_norm(1)+l_norm(1))/2-m_norm(1) (s_norm(2)+l_norm(2))/2-m_norm(2)];
    end
    h = plot(blue_orange_signals(:,1),-blue_orange_signals(:,2),'-');
    set(h,'color',orange_cyan_color);
    h = plot(lime_magenta_signals(:,1),-lime_magenta_signals(:,2),'-');
    set(h,'color',lime_magenta_color);

    centered_blue_orange = blue_orange_signals-repmat(mean(blue_orange_signals),size(blue_orange_signals,1),1);
    centered_lime_magenta = lime_magenta_signals-repmat(mean(lime_magenta_signals),size(lime_magenta_signals,1),1);
    blue_orange_var_inline = var(centered_blue_orange*[1/sqrt(2); 1/sqrt(2)]);
    blue_orange_var_orth = var(centered_blue_orange*[1/sqrt(2); -1/sqrt(2)]);
    lime_magenta_var_inline = var(centered_lime_magenta*[1/sqrt(2); 1/sqrt(2)]);
    lime_magenta_var_orth = var(centered_lime_magenta*[1/sqrt(2); -1/sqrt(2)]);
    variances = [variances; blue_orange_var_inline blue_orange_var_orth lime_magenta_var_inline lime_magenta_var_orth];
end
xlabel('Subunit 1 activation','FontSize',15);
ylabel('Subunit 2 activation','FontSize',15);
set(gca,'Xlim',[-.41 .41],'Xtick',[-.4 0 .4],'Xticklabel',['-';'0';'+'])
set(gca,'Ylim',[-.41 .41],'Ytick',[-.4 0 .4],'Yticklabel',['-';'0';'+'])

figure; axes; hold on; set(gca,'TickDir','out'); axis square;
ratio1 = log10(variances(:,2))-log10(variances(:,1));
ratio2 = log10(variances(:,4))-log10(variances(:,3));
bins = linspace(min([ratio1; ratio2]),max([ratio1; ratio2]),20);
[n1,~] = hist(ratio1,bins);
[n2,~] = hist(ratio2,bins);
h = bar(x,[n1;n2]');
set(h(1),'FaceColor',orange_cyan_color);
set(h(2),'FaceColor',lime_magenta_color);
xlabel('log variance ratio','FontSize',15); ylabel('count','FontSize',15);

% And now a single example
chips = [418 428];
reflectance1 = munsell(:,chips(1));
reflectance2 = munsell(:,chips(2));
light1 = illuminants.*reflectance1;
light2 = illuminants.*reflectance2;
light0 = repmat(mean(illuminants,2),1,nillums); % a single mean illuminant

lms_1 = T_cones_smj10*light1;
lms_2 = T_cones_smj10*light2;
lms_0 = T_cones_smj10*light0;

coneexcitations = cat(3, lms_1',lms_2',lms_0'); % rows = illums, cols = LMS, planes = light1, light2, adaptating light
coneexcitations = permute(coneexcitations,[2 3 1]);
blue_orange_signals = []; lime_magenta_signals = [];
for i = 1:size(illuminants,2)
    l = squeeze(coneexcitations(1,:,i));
    m = squeeze(coneexcitations(2,:,i));
    s = squeeze(coneexcitations(3,:,i));
    l_norm = (l-l(3))./l(3);
    m_norm = (m-m(3))./m(3);
    s_norm = (s-s(3))./s(3);
    blue_orange_signals = [blue_orange_signals;(s_norm(1)+m_norm(1))/2-l_norm(1) (s_norm(2)+m_norm(2))/2-l_norm(2)];
    lime_magenta_signals = [lime_magenta_signals;(s_norm(1)+l_norm(1))/2-m_norm(1) (s_norm(2)+l_norm(2))/2-m_norm(2)];
end

% Rendering the munsell surfaces under the illuminants
figure; subplot(2,2,1); hold on;
wls = [380:5:780];
plot(wls,reflectance1,'-','LineWidth',2,'color',[0 0 0]);
plot(wls,reflectance2,'-','LineWidth',2,'color',[.5 .5 .5]);
plot(wls,reflectance1,'-','LineWidth',2,'color',[0 0 0]);
set(gca,'Xlim',[400 700],'Ytick',[0 .25 .5 .75 1]);
ylabel('reflectance'); xlabel('wavelength (nm)');

subplot(2,2,2); hold on;
plot(wls,illuminants(:,1),'-','LineWidth',2,'color',[1 0 0]);
plot(wls,illuminants(:,end),'-','LineWidth',2,'color',[0 0 1]);
set(gca,'Xlim',[400 700],'Ytick',[0 50 100 150 200]);
ylabel('power'); xlabel('wavelength (nm)');

% Suffix convention: illuminant, reflectance
%subplot(2,2,3); hold on;
%plot((coneexcitations(:,1,end)-coneexcitations(:,3,end))./coneexcitations(:,3,end),'Color',[0 0 1]);
%plot((coneexcitations(:,2,end)-coneexcitations(:,3,end))./coneexcitations(:,3,end),'Color',[0 0 1]);
%plot((coneexcitations(:,1,1)-coneexcitations(:,3,1))./coneexcitations(:,3,end),'Color',[1 0 0]);
%plot((coneexcitations(:,2,1)-coneexcitations(:,3,1))./coneexcitations(:,3,end),'Color',[1 0 0]);
%set(gca,'Xlim',[.5 3.5],'Xtick',[1 2 3],'Xticklabel',['L';'M';'S']);
% S-cone change goes with M-cone change in both cases
% As illuminant change from blue to red, S and M drop and L increases.
% But S drops disproportionately. For the less-reflective surface
% especially, S-cone signal drops a lot.


% RGB rendering
load('Dell4BitsCal');
cal = cals{end};
spds = SplineSpd([380:4:780]',cals{end}.P_device,[380:5:780]');
M = T_cones_smj10*spds;
rgb21= inv(M)*coneexcitations(:,1,end);
rgb22= inv(M)*coneexcitations(:,2,end);
rgb11= inv(M)*coneexcitations(:,1,1);
rgb12= inv(M)*coneexcitations(:,2,1);
%rgbbkgnd= inv(M)*coneexcitations(:,3,end);
%rgbbkgnd = rgbbkgnd/2; % Assuming background is a flat reflector

rgbnormfact = max([rgb21;rgb22;rgb11;rgb12]);

rgb21 = max(rgb21./rgbnormfact,0);
rgb22 = max(rgb22./rgbnormfact,0);
rgb11 = max(rgb11./rgbnormfact,0);
rgb12 = max(rgb12./rgbnormfact,0);
%rgbbkgnd = max(rgbbkgnd./rgbnormfact,0);
im = cat(3,[rgb21, rgb22],[rgb11, rgb12]);
subplot(2,2,3);
im = permute(im,[2 3 1]);
image(im)
axis square;
set(gca,'Visible','off');

subplot(2,2,4); set(gca,'TickDir','out','color',[.5 .5 .5]); hold on;
plot(blue_orange_signals(:,1),-blue_orange_signals(:,2),'-','LineWidth',2,'color',orange_cyan_color);
plot(lime_magenta_signals(:,1),-lime_magenta_signals(:,2),'-','LineWidth',2,'color',lime_magenta_color);
axis square;
set(gca,'Ylim',[0 .4],'Xlim',[-.2 .2]);
