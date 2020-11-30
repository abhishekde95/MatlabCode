% Figures for a second IsoSamp paper (about chromatic temporal contrast sensitivity)
% Contents
%
% Section 0) Just setting stuff up. Should be identical to section 0 of
% IsoSampPaperStuff.m except for the identities of the example neurons.
%
% Section 1) Cone weights
%
% Section 1.1) Cone weights with fundamental modified from the Stockman,
% Johnson, and MacLeod 10° fundamentals. 
%
% Section 2) Example rasters and PSTHs for an example P and M cell. Taken
% largely from IsoSampPaperStuff.m section 4.
%
% Section 2.1) d' vs TF curves for single cells. Must run sxn 2 first.
%
% Section 2.2) PSTHs to ask how stationary is the response of the
% course of the trial?
%
% Section 3) Amplitude spectra of FFTs of spike trains
%
% Section 4) Cone weights from white noise predict L-M sensitivity in
% IsoSamp. Also eccentricity and ON vs OFF analysis.
%
% Section 5) How d-prime for both L+M and L-M changes as a function of
% temporal frequency.
%
% Section 5.1) Within neuron analysis of d-prime changes as a function of
% temporal frequency for both L+M and L-M.
%
% Section 6) KNN analysis (nonparametric d')
%
% Section 7) Population parvocellular sensitivity for L-M as a function 
% of spike counting window. Taken from IsoSampPop Section 14.
%
% Section 7.1) Various spike counting windows. Comparing M to P and P to
% d' = 1.27.
%
% Section 8) Population sensitivity analysis comparing M(L+M), M(L-M), 
% P(L+M), and P(L-M) with and without latency offset.
%
% Section 9) Population sensitivity analysis comparing M (L+M) and P(L-M)
% to cone currents.
%
% Section 10) Behavioral contrast sensitivity functions for Apollo and
% Utu.
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

TFBINCENTERS = logspace(log10(1), log10(60), 15); % Native bins
binwidth = log10(TFBINCENTERS(2))-log10(TFBINCENTERS(1));
TFBINEDGES = 10.^(linspace(log10(TFBINCENTERS(1))-binwidth/2,log10(TFBINCENTERS(end))+binwidth/2,length(TFBINCENTERS)+1));
MONKEY2MARKER = '^';
EXAMPLEMAGNOCELL = 'A070717002.nex'; % example magnocell. Different from IsoSampPaperStuff.m
EXAMPLEPARVOCELL = 'A122117002.nex'; % example parvocell. Different from IsoSampPaperStuff.m
EXAMPLEKONIOCELL = 'A062517002.nex'; % strongish responses to L+M

GETRFFROMSTAMODE = 2;
GETRFFROMSTATHRESH = .95;
MAXT = 6; % n frames back to compute STA

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
% Section 1: Cone weights
% Section 20 of IsoSampPaperStuff.m
% (Taken from IsoSampPop.m section 5)
%%
% Section 1.1) Cone weights with fundamental modified from the Stockman,
% Johnson, and MacLeod 10° fundamentals.
% Based on TenDegFundStuff.m, section 7.2

% Removing macular pigment and lens from SMJ 10 deg fundamentals
load('T_cones_smj10');
load('den_mac_ws');
load('den_lens_ws');

%plot(380:5:780,den_lens_ws*1);

macpigtransmittance = 1./(10.^(den_mac_ws*.056)); % Removing mac pig (0.28 = entirely)
lenstransmittance = 1./(10.^(den_lens_ws.*.20));  % Reducing lens (by 20%)?
%macpigtransmittance = 1./(10.^(den_mac_ws*.28)); % Removing mac pig (0.28 = entirely)
%lenstransmittance = 1./(10.^(den_lens_ws.*.50));  % Reducing lens (by 20%)?
%absorptance = T_cones_smj10./repmat(macpigtransmittance',3,1);
absorptance = T_cones_smj10./repmat(lenstransmittance',3,1)./repmat(macpigtransmittance',3,1);
absorptance = absorptance./repmat(max(absorptance,[],2),1,81);
% Dividing because we're assuming that these cone fundamentals already
% include these filters and we're removing them. Low numbers (low
% transmittance) means light is getting absorbed. Dividing by a small
% number increases the corneal sensitivity.

% SMJ assume 10° fundamentals contain .28 WS macular pigment template

figure; 
subplot(2,1,1); hold on;
plot(T_cones_smj10','k-')
plot(absorptance','r-')
%subplot(2,1,2); hold on;
a = T_cones_smj10.*repmat(10.^(.056.*den_mac_ws'+.2.*den_lens_ws'),3,1);
a = a./repmat(max(a,[],2),1,size(a,2));
plot(a','k--','Linewidth',2)

(T_cones_smj10-absorptance)./T_cones_smj10

% How much does a nominal L+M or L-M stimulus change?
load Propixx
cal = cals{end};
spd = SplineSpd([380:4:780]',cal.P_device,[380:5:780]');
Morig = T_cones_smj10*spd
Mnew = absorptance*spd;
rgbbgkgnd= [.5 .5 .5]';

lmsbkgnd = Morig*rgbbgkgnd;
LMstim = lmsbkgnd.*(1+[0.1 0.1 0]')
cc = (LMstim-lmsbkgnd)./lmsbkgnd
RGBstim1 = inv(Morig)*LMstim;

LvMstim = lmsbkgnd.*(1+[0.1 -0.1 0]')
cc = (LvMstim-lmsbkgnd)./lmsbkgnd
RGBstim2 = inv(Morig)*LvMstim;

% Converting back though Mnew
lmsbkgnd = Mnew*rgbbgkgnd;
LMstim1 = Mnew*RGBstim1;
cc1 = (LMstim1-lmsbkgnd)./lmsbkgnd
LMstim2 = Mnew*RGBstim2
cc2 = (LMstim2-lmsbkgnd)./lmsbkgnd

(cc1'-[0.1 0.1 0])./[0.1 0.1 0] % Contrast to the L and M cones reduced by 10% of 10% (1%), S-cone < 1%
(cc2'-[0.1 -0.1 0])./[0.1 -0.1 0] % < M cones are getting 10% more contrast than L, S-cone ~1%


% Now computing cc with new fundamentals






%%
% Section 2: Example rasters and PSTHs
% Rasters of example neurons. Taken from Section 4 of "IsoSampFigures.m"
% Assumes spike 1.
filenames = {EXAMPLEMAGNOCELL, EXAMPLEPARVOCELL, EXAMPLEKONIOCELL};
filenames = {EXAMPLEMAGNOCELL};
STIMTYPES = {'RG','LUM'};
spikeNum = 'sig001a';
AXHEIGHT = 6.5;
AXHMARGIN = 2;
AXVMARGIN = 1;
AXWIDTH = 4;
PSTHSCALEFACTOR = 15;
AXVSTART = 9.3;

axpositions = [3 AXVSTART; 3+AXWIDTH+AXHMARGIN AXVSTART+AXHEIGHT+AXVMARGIN];
axpositions = [3 9 15];

figprefs;
for filenameidx = 1:length(filenames)
    for STIMTYPEidx = 1:length(STIMTYPES)
        STIMTYPE = STIMTYPES{STIMTYPEidx};
        binwidth = 0.008; % bins for PSTH
        filename = char(filenames(filenameidx));
        stro = nex2stro(findfile(filename,[nexfilepath,filesep,'Greg',filesep,'Apollo']));
        
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
        axes('position',[axpositions(filenameidx) AXVSTART+(STIMTYPEidx-1)*(AXHEIGHT+AXVMARGIN) AXWIDTH AXHEIGHT]);
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
            h = plot(bins,counter+1+psth./sum(L)./binwidth/PSTHSCALEFACTOR,'-','Linewidth',2,'Color',[.5 .5 .5]);
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
                h = text(-.12,counter-sum(L)/2,'Blank');
                h(2) = text(dur+offset(2)+.1,counter-sum(L)/2,'0.00');
            else
                h = text(-.12,counter-sum(L)/2,num2str(round(uniquestim(j,3)*10)/10));
                % BELOW: Using L-cone contrast == M-cone contrast as a proxy for luminance contrast
                % (see BenardeteKaplanSims.m)
                h(2) = text(dur+offset(2)+.1,counter-sum(L)/2,num2str(round(abs(uniquestim(j,1))*100)/100));
            end
            set(h,'HorizontalAlignment','right','FontSize',8,'FontAngle','Italic');
            %if uniquestim(j,1) == 0 & uniquestim(j,2) == 0
            %    plot([0 dur dur 0 0],[counter counter 0 0 counter],'k-','LineWidth',1);
            %end
            counter = counter + 1;
        end
        plot([0 0],[0 counter],'k-','LineWidth',.5);
        plot([dur dur],[0 counter],'k-','LineWidth',.5);
        set(gca,'Xlim',[0+offset(1) dur+offset(2)],'Ylim',[0 counter],'Ytick',[],'Xtick',[0 .2 .4 .6 .8],'Box','off','FontSize',8,'FontAngle','italic');
        set(gcf,'Renderer','painters');
      
        if strcmp(stro.sum.fileName(find(stro.sum.fileName == filesep,1,'last')+1:end), EXAMPLEMAGNOCELL)
            titlestr = 'Magno';
        elseif strcmp(stro.sum.fileName(find(stro.sum.fileName == filesep,1,'last')+1:end), EXAMPLEPARVOCELL)
            titlestr = 'Parvo';
        elseif strcmp(stro.sum.fileName(find(stro.sum.fileName == filesep,1,'last')+1:end), EXAMPLEKONIOCELL)
            titlestr = 'Konio';
        else
            error ('unknown neuron');
        end
        if strcmp(STIMTYPE, 'LUM')
            titlestr = [titlestr,' (L+M)'];
        else
            titlestr = [titlestr,' (L-M)'];
        end
        title(titlestr,'FontSize',12);
        
        if STIMTYPEidx == 1
            xlabel('Time (s)','FontSize',10);
        else
            set(gca,'Xticklabel',[]);
        end
        
        if filenameidx == 1
            h(3) = text(-1,AXHEIGHT/2,'Frequency (Hz)','HorizontalAlignment','center','Units','centimeters','Rotation',90,'FontAngle','italic');
        end
        if filenameidx == length(filenames)
            h(4) = text(AXWIDTH+1,AXHEIGHT/2,'Contrast','HorizontalAlignment','center','Units','centimeters','Rotation',-90,'FontAngle','italic');
        end
    end 
end

%%
% Section 2.1
% Single neuron d'. Parvo/Magno x L+M/L-M
% Need to run this after section 2, above

AXYPOS = 2;
axpositions = [3 AXYPOS; 3+AXWIDTH+AXHMARGIN AXYPOS];
axpositions = [3+[0:length(filenames)-1]'*(AXWIDTH+AXHMARGIN), repmat(AXYPOS,length(filenames),1)];
colors = [1 0 1; 0 0 0];

for filenameidx = 1:length(filenames)
    axes('Position',[axpositions(filenameidx,1) axpositions(filenameidx,2) AXWIDTH AXWIDTH])
    plot([1 60], [0 0],'k-','LineWidth',1);
    for STIMTYPEidx = length(STIMTYPES):-1:1
        STIMTYPE = STIMTYPES{STIMTYPEidx};
        filename = char(filenames(filenameidx));
        stro = nex2stro(findfile(filename,[nexfilepath,filesep,'Greg',filesep,'Apollo']));
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
        
        % Dprime plot
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
        patch([uniquestim(Lstimtype,3)', flipud(uniquestim(Lstimtype,3))'],[dprime(Lstimtype)'+se', fliplr(dprime(Lstimtype)'-se')],colors(STIMTYPEidx,:),'Facealpha',.5,'EdgeColor','none');
        plot(uniquestim(Lstimtype,3),dprime(Lstimtype),'ko-','LineWidth',1,'MarkerSize',5,'MarkerEdgeColor','none','MarkerFaceColor',colors(STIMTYPEidx,:),'Color',colors(STIMTYPEidx,:));
    end
    set(gca,'Xlim',[TFBINEDGES(1) TFBINEDGES(end)],'Xscale','log');
    set(gca,'FontSize',8,'FontAngle','italic','XTick',[1 10 100],'XtickLabel',{'1','10','100'},'Ylim',[-2 7],'YTick',[-2:2:6]);
    xlabel('Temporal frequency (Hz)','FontSize',10);
    ylabel('Signal-to-noise ratio (d'')','FontSize',10);
    text(18,-1.5,filename(1:end-4),'FontSize',6);
end

%%
% Section 2.2 
% PSTHs for many cells to look for non-stationary responses over individual
% trials.

maxf = 120; % Maximum frequency in Hz. From section 3, below
timebins = [0:1/(maxf*2):.800];
deltaT = timebins(2)-timebins(1);
MAXNUMBEROFNEURONS = 40; % For memory preallocation. Number of neurons per cell type.

data = nan*ones(3, 2,length(TFBINCENTERS),MAXNUMBEROFNEURONS,length(timebins)); % cell class x colordir x tf x cell x time bins
cellclasses = {'P','M','K'};
colordirections = {'L+M','L-M'};
for cellclassidx = 1:3
    cellclass = cellclasses{cellclassidx};
    [filenames, ~, spikenames, neuronids] = fnamesFromTxt('IsoSamp_LGN','cellClass',{cellclass});

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
        
        tmpdata = nan*ones(ntrials,length(timebins));
        for j=1:ntrials
            spikeidx = strcmp(stro.sum.rasterCells,tmpspikename);
            tmp = stro.ras{j,spikeidx}-stimon_t(j);
            tmp(tmp<timebins(1)) = [];
            tmp(tmp>timebins(end)) = [];
            tmpdata(j,:) = hist(tmp,timebins)./deltaT; % in sp/s
        end
        
        Llum = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & uniquestim(:,3) > 0;
        Lrg = sign(uniquestim(:,1)) ~= sign(uniquestim(:,2)) & uniquestim(:,3) > 0;
               
        for colordiridx = 1:2 % L+M and L-M
            if colordiridx == 1
                whichconditions = find(Llum);
            else
                whichconditions = find(Lrg);
            end
            for k = 1:length(whichconditions)
                L = Lcc == uniquestim(whichconditions(k),1) & Mcc == uniquestim(whichconditions(k),2) & TF == uniquestim(whichconditions(k),3);
                PSTH = mean(tmpdata(L,:));
                if sum(softEq(uniquestim(whichconditions(colordiridx),3),TFBINCENTERS,3)) ~= 1
                    disp('error')
                    keyboard
                end
                if i > MAXNUMBEROFNEURONS
                    error('Increase MAXNUMBEROFNEURONS');
                end
                data(cellclassidx,colordiridx,softEq(uniquestim(whichconditions(k),3),TFBINCENTERS,3),i,:) = PSTH; % colordir x tf x cell x timebin
            end
        end
    end
end

% Plotting
for colordiridx = 1:2
    for cellclassidx = 1:3
        figure; set(gcf,'Name',['Color direction ',colordirections{colordiridx},' Cell class ',cellclasses{cellclassidx}]);
        for TFidx = 1:length(TFBINCENTERS)
            subplot(4,4,TFidx); hold on;
            PSTHs = [];
            for cellcounter = 1:MAXNUMBEROFNEURONS
                tmp = squeeze(data(cellclassidx,colordiridx,TFidx,cellcounter,:))';
                tmp = tmp - mean(tmp); % otherwise hard to define sign of modulation
                if ~all(isnan(tmp))
                    if isempty(PSTHs)
                        if mean(tmp(timebins<.20)) > mean(tmp(timebins>.2)) % trying to get peaks instead of troughs
                            tmp = -tmp;
                        end
                        PSTHs = tmp;
                    else
                        phase_sign = sign(nanmean(PSTHs,1)*tmp');
                        PSTHs = [PSTHs; phase_sign*tmp];
                    end
                end
            end
            if ~isempty(PSTHs)
                title(['TF: ',num2str(TFBINCENTERS(TFidx)),' n = ',num2str(sum(~isnan(PSTHs(:,1))))]);
                plot(timebins, mean(PSTHs,1));
                ylims = get(gca,'Ylim');
                plot([0 0],[ylims(1) ylims(2)],'k-'); % stim on
                plot([.660 .660],[ylims(1) ylims(2)],'k-'); % stim off
                plot(.120+[0 0],[ylims(1) ylims(2)],'k:'); % P-cell counting window start
                plot(.120+[.660 .660],[ylims(1) ylims(2)],'k:'); % P-cell counting window end
                set(gca,'Xlim',[timebins(1) timebins(end)]);
            end
        end
    end
end

% Where is the DC shift in M-cell responses with TF coming from?
% Magno cell #7 shows the effect strongly. 
% Is the issue the tightness of the spiking responses?
colordiridx = 1;  % 1 = L+M, 2 = L-M
cellclassidx = 2; % 2 = Magno
colors = jet(length(TFBINCENTERS));
for cellcounter = 1:15 % n magnocellular neurons
    PSTHs = [];
    for TFidx = 1:length(TFBINCENTERS)
        tmp = squeeze(data(cellclassidx,colordiridx,TFidx,cellcounter,:))';
        if ~all(isnan(tmp))
            if isempty(PSTHs)
                PSTHs = tmp;
            else
                PSTHs = [PSTHs; tmp];
            end
        end
    end
    figure; subplot(1,2,1);
    %imagesc(PSTHs);
    h = plot(PSTHs');
    for i = 1:length(h)
        set(h(i),'Color',colors(i,:))
    end
    
    subplot(1,2,2);
    h = plot(fftshift(abs(fft(PSTHs')'),2)');
    for i = 1:length(h)
        set(h(i),'Color',colors(i,:))
    end
    set(gca,'Yscale','log')
end

% % Trying to figure out where the broadband shift comes from 
% t = [0:100];
% x = rem(t,5)== 0;
% figure; subplot(2,1,1);
% plot(x);
% subplot(2,1,2);
% plot(fftshift(abs(fft(x))));

figure;
scale = 4;
exponent = 1;
t = [0:1000];
x = cos(2*pi*t*.01);
x = x*scale;
x(x>0) = x(x>0).^exponent;
%x(x<0) = 0;
 subplot(2,1,1); hold on;
plot(t,x);
subplot(2,1,2); hold on;
plot(fftshift(abs(fft(x))));
set(gca,'Yscale','log')

%%
% Section 3
% Power spectral densities of spike trains

maxf = 120; % Maximum frequency in Hz
timebins = [0:1/(maxf*2):.660];
deltaT = timebins(2)-timebins(1);
nyquist = 1./(2*deltaT);
MAXNUMBEROFNEURONS = 40; % For memory preallocation. Number of neurons per cell type.

fftdata = nan*ones(2, 2,length(TFBINCENTERS),MAXNUMBEROFNEURONS,length(timebins)); % cell class x colordir x tf x cell x hz
for cellclassidx = 1:2
    if cellclassidx == 1
        cellclass = 'P';
    else
        cellclass = 'M';
    end
    [filenames, ~, spikenames, neuronids] = fnamesFromTxt('IsoSamp_LGN','cellClass',{cellclass});

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
            n = hist(tmp,timebins)./deltaT; % in sp/s
            tmpfftdata(j,:) = fftshift(fft(n)./length(n)); % Dividing by nsamples (https://www.mathworks.com/matlabcentral/answers/162846-amplitude-of-signal-after-fft-operation)
        end
        
        freqs = linspace(-nyquist, nyquist,size(tmpfftdata,2));
        Llum = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & uniquestim(:,3) > 0;
        Lrg = sign(uniquestim(:,1)) ~= sign(uniquestim(:,2)) & uniquestim(:,3) > 0;
        
        colors = hot(length(unique(TF)));
        
        for colordiridx = 1:2 % L+M and L-M
            if colordiridx == 1
                whichconditions = find(Llum);
            else
                whichconditions = find(Lrg);
            end
            for k = 1:length(whichconditions)
                L = Lcc == uniquestim(whichconditions(k),1) & Mcc == uniquestim(whichconditions(k),2) & TF == uniquestim(whichconditions(k),3);
                amp = mean(abs(tmpfftdata(L,:)));
                if sum(softEq(uniquestim(whichconditions(colordiridx),3),TFBINCENTERS,3)) ~= 1
                    disp('error')
                    keyboard
                end
                if i > MAXNUMBEROFNEURONS
                    error('Increase MAXNUMBEROFNEURONS');
                end
                fftdata(cellclassidx,colordiridx,softEq(uniquestim(whichconditions(k),3),TFBINCENTERS,3),i,:) = amp; % colordir x tf x cell x hz
            end
        end
    end
end

AXHEIGHT = 5;
AXWIDTH = 5;
AXHMARGIN = 2;
AXVMARGIN = 2;
MINPOWER = 5;
MAXPOWER = 500;
axpositions = [4 12];

colors = jet(size(fftdata,3));
figprefs;
celltypelabels = {'Parvocellular','Magnocellular'};
colordirlabels = {'L+M','L-M'};
for i = 1:size(fftdata,1) % cell type
    for j = 1:size(fftdata,2) % colordir
        axes('position',[axpositions(1)+(i-1)*(AXWIDTH+AXHMARGIN) axpositions(2)-(j-1)*(AXHEIGHT+AXVMARGIN) AXWIDTH AXHEIGHT])
        hold on;
        for k = 1:size(fftdata,3) % stim TF
            if sum(~isnan(fftdata(i,j,k,:,1)),4) > 1 % Don't use conditions with only one observation
                tmp = squeeze(nanmean(fftdata(i,j,k,:,:),4));
                %df = freqs(2)-freqs(1);
                %plot(freqs,tmp.^2/(2*df),'color',colors(k,:));
                L = freqs >= 0;
                plot(freqs(L),tmp(L).^2,'color',colors(k,:));
                plot(TFBINCENTERS(k),MINPOWER,'^','MarkerEdgeColor','none','MarkerFaceColor',colors(k,:))
            end
        end
        
        set(gca,'yscale','log');
        title([celltypelabels{i},' ',colordirlabels{j}]);
        set(gca,'Xlim',[0 maxf],'Ylim',[MINPOWER MAXPOWER]);
        if i == 1
            ylabel('Power (spikes^2/Hz)');
        end
        if j == 2
            xlabel('Frequency (Hz)');
        end
        set(gca,'Ytick',[1 10 100],'YtickLabel',[1 10 100]);
    end
end
equatesubplotaxeslims;

% Colorbar
axes('position',[axpositions(1)+i*(AXWIDTH+AXHMARGIN)-AXHMARGIN/2 axpositions(2)-(j-1)*(AXHEIGHT+AXVMARGIN) .5 AXHEIGHT])
image([1:size(colors,1)]');
colormap(colors);
set(gca,'YAxisLocation','right');
set(gca,'Ylim',[.5 size(colors,1)+.5],'Xtick',[],'Ytick',[1:2:size(colors,1)],'Yticklabels',round(TFBINCENTERS(1:2:end)*10)/10);
ylabel('Stimulus frequency (Hz)');

% Comparing F2 to F1 magnitude.
% fftdata: [P,M], [L+M, L-M], TF, neuron, SF bin
tmpdata = squeeze(fftdata(2,1,:,:,:));
whichfreqsbin = [];
for i = 1:length(TFBINCENTERS)
    err = (freqs-TFBINCENTERS(i)).^2;
    L1 = err == min(err);
    err = (freqs-TFBINCENTERS(i)*2).^2;
    L2 = err == min(err);
    whichfreqsbin(i,:) = [find(L1,1,'first') find(L2,1,'first')]; % "first" here is a hack. 60 Hz falls directly between FFT bins!
end

% Note, some unique TFBINCENTERS map to the same fft frequency
% so these plots are not independent.
% Looking at F1 vs F2 across neurons with frequency bin
figure; set(gcf,'Name','across cells, within frequency');
data = [];
whichfreqs = TFBINCENTERS > 2; % Skipping the first three bins because of trivial overlap in power (f1 and f2 are in the same FFT bin for 1 Hz)
for i = find(whichfreqs) 
    subplot(4,4,i);
    F1 = tmpdata(i,:,whichfreqsbin(i,1));
    F2 = tmpdata(i,:,whichfreqsbin(i,2));
    plot(F1,F2,'o')
    r = corrcoef([F1(~isnan(F1))' F2(~isnan(F2))']);
    title(['r = ',num2str(r(1,2))]);
    data = [data; r(1,2)];
end
equatesubplotaxeslims(1)
[~,p] = ttest(data);
disp(num2str([mean(data) p]));

% Looking at F1 vs F2 across frequency bin within neuron
ncells = sum(~isnan(tmpdata(1,:,1)));
figure; set(gcf,'Name','across frequency, within cells');
data = [];
for i = 1:ncells
    subplot(4,4,i);
    F1 = diag(squeeze(tmpdata(whichfreqs,i,whichfreqsbin(whichfreqs,1))));
    F2 = diag(squeeze(tmpdata(whichfreqs,i,whichfreqsbin(whichfreqs,2))));
    plot(F1,F2,'o')
    r = corrcoef([F1(~isnan(F1)) F2(~isnan(F2))]);
    title(['r = ',num2str(r(1,2))]);
    data = [data; r(1,2)];
end
equatesubplotaxeslims(1);
[~,p] = ttest(data);
disp(num2str([mean(data) p]));


%%
% Section 4
% Cone weights from white noise predict sensitivity to L-M vs. L+M in
% IsoSamp. Also, eccentricity analysis.

[WNfilenames, ~, WNspikenames, WNneuronids] = fnamesFromTxt('WhiteNoiseLGN_forIS','cellClass',{'M','P','K'});
[ISfilenames, ~, ISspikenames, ISneuronids] = fnamesFromTxt('IsoSamp_LGN','cellClass',{'M','P','K'});

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

% First, getting cone weights
coneweights = zeros(length(strippeddownWNneuronids), 3);
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
    
    funds = stro.sum.exptParams.fundamentals;
    funds = reshape(funds,[length(funds)/3,3]);
    mon_spd = stro.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = funds'*mon_spd;
    
    ridx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_r'));
    gidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_g'));
    bidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_b'));
    bkgndrgb = [mode(stro.trial(:,ridx)), mode(stro.trial(:,gidx)), mode(stro.trial(:,bidx))]';
    bkgndlms = M*bkgndrgb; % this is in cone excitations
    Mrgbtocc = diag(1./bkgndlms)*M;
    M = Mrgbtocc; % Assumes a weighted sum of cone contrasts

    out = getWhtnsStats(stro,MAXT,'STCOVmex', {nstixperside^2, 3, MAXT}, spikename);
    STA = out{1}/out{3};
    inRF = getRFfromSTA(STA,GETRFFROMSTAMODE,GETRFFROMSTATHRESH);
    STA = reshape(STA,[nstixperside.^2  3 MAXT]);
    temporalSTA = squeeze(STA(logical(inRF(:)),:,:));
    if ~ismatrix(temporalSTA)
        temporalSTA = mean(temporalSTA,1);
    end
    sta = squeeze(temporalSTA);
    [u,s,v] = svd(sta);
    rgb = u(:,1);
    if rgb'*sta(:,sum(sta.^2) == max(sum(sta.^2))) < 0
        rgb = -rgb;
    end
    coneweights(i,:) = rgb'*inv(M);
end
normalizedconeweights = coneweights./repmat(sum(abs(coneweights),2),1,3);

% Getting IsoSamp data for L-M and L+M
% averaged across temporal frequency
data = [];
for i = 1:length(strippeddownWNneuronids)
    idx = strippeddownWNneuronids(i);
    filenames = ISfilenames{ISneuronids == idx};
    spikenames = ISspikenames(ISneuronids == idx,:);
     
    if any(strcmp(filenames, EXAMPLEMAGNOCELL)) || any(strcmp(filenames, EXAMPLEPARVOCELL)) || any(strcmp(filenames, EXAMPLEKONIOCELL))
        examplecellflag = 1;
    else 
        examplecellflag = 0;
    end
    
    stro=[];
    for j = 1:length(filenames)
        stro = strocat(stro,nex2stro(findfile(filenames{j})));
    end
    [uniquestim, dprime] = IsoSampGetDPrime(stro,DPRIMEMETHOD,abs(spikenames(end)-96));
    Lblank = uniquestim(:,1) == 0 & uniquestim(:,2) == 0 & uniquestim(:,3) == 0;
    Llum = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & ~Lblank;
    Lrg = sign(uniquestim(:,1)) ~= sign(uniquestim(:,2)) & ~Lblank;
    
    % For starters, let's average across TFs
    data = [data; mean(dprime(Llum)) mean(dprime(Lrg)) stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y examplecellflag];
end

% Plotting
colors = [0 .5 1; 0 0 0; 1 0 0];
uniquecellTypes = unique(cellTypes);
figprefs; axes('Position',[3 3 7 7]); hold on;
LM = abs(normalizedconeweights(:,1)+normalizedconeweights(:,2));
examplecell_L = data(:,5);
plot([0 1],[0 0],'k--','LineWidth',1);
for whichcelltype = {'M','P','K'}
    L = strcmp(cellTypes,whichcelltype);
    plot(LM(L&~examplecell_L),data(L&~examplecell_L,2)-data(L&~examplecell_L,1),'wo','markerfacecolor',colors(strcmp(whichcelltype, uniquecellTypes),:),'MarkerSize',8,'LineWidth',.1);
    plot(LM(L&examplecell_L),data(L&examplecell_L,2)-data(L&examplecell_L,1),'w^','markerfacecolor',colors(strcmp(whichcelltype, uniquecellTypes),:),'MarkerSize',9,'LineWidth',.1);
    % regressions
    %[b,bint, ~, ~, stats] = regress(data(L,2)-data(L,1), [ones(sum(L),1) LM(L)]);
    %disp([whichcelltype{1},' slope d''color-d''lum vs |L+M| cone weights: ',num2str(b(2)),' p:',num2str(stats(3))]);
    [r,p] = corr(data(L,2)-data(L,1), LM(L));
    disp([whichcelltype{1},' Pearson correlation: d''color-d''lum vs |L+M| cone weights: ',num2str(r),' p:',num2str(p)]);
end
xlabel('L+M cone weights from white noise');
ylabel('d''_{L-M} - d''_{L+M} ')
set(gca,'Ylim',[-2 2],'YTick',-2:2);

% d'(L+M) vs d'(L-M)
axes('Position',[13 3 7 7]); hold on;
plot([0 2],[0 2],'k--','LineWidth',1);
for whichcelltype = {'M','P','K'}
    L = strcmp(cellTypes,whichcelltype);
    plot(data(L&~examplecell_L,1), data(L&~examplecell_L,2),'wo','markerfacecolor',colors(strcmp(whichcelltype, uniquecellTypes),:),'MarkerSize',8);
    plot(data(L&examplecell_L,1), data(L&examplecell_L,2),'w^','markerfacecolor',colors(strcmp(whichcelltype, uniquecellTypes),:),'MarkerSize',9);
    [r,p] = corrcoef(data(L,1), data(L,2));
    disp(['d''(L+M) vs. d''(L-M), cell type ',whichcelltype{1},', r=',num2str(r(1,2)),', p=',num2str(p(1,2))]);
end
xlabel('Signal-to-noise ratio (d''_{L+M})');
ylabel('Signal-to-noise ratio (d''_{L-M})')
set(gca,'Xlim',[-.5 2],'Ylim',[-.5 2],'Xtick',[0 1 2],'Ytick',[0 1 2]);

% Regression comparing raw d-prime values (and lum/chrom differences) with
% eccentricity
colordirlabels = {'L+M','L-M'};
for whichcelltype = uniquecellTypes'
    L = strcmp(cellTypes,whichcelltype{1});
    ecc = sqrt(data(L,3).^2+data(L,4).^2)/10;
    % First just straight-up regressions on d' as a function of
    % eccentricity
    
    if strcmp(whichcelltype{1},'M')
        axes('Position',[3 13 7 7]); hold on; plottrue = true;
    elseif strcmp(whichcelltype{1},'P')
        axes('Position',[13 13 7 7]); hold on;plottrue = true;
    else
        plottrue = false;
    end
    
    X = [ones(2*length(ecc),1) [zeros(length(ecc),1); ones(length(ecc),1)] [ecc; ecc] [zeros(length(ecc),1); ecc]];
    [b,bint, ~, ~, stats] = regress([data(L,1); data(L,2)] , X);
    disp([whichcelltype{1},': beta(color_dir x slope(d'' vs ecc)) =',num2str(b(4)),', CI',num2str(bint(4,:))]);
   
    for i = 1:2 % color direction 1 = L+M, 2 = L-M
        [b,~, ~, ~, stats] = regress(data(L,i), [ones(length(ecc),1) ecc]);
        disp([whichcelltype{1},', ',colordirlabels{i},', slope(d'' vs ecc)=',num2str(b(2)),', p=',num2str(stats(3))]);
        
        if plottrue
            if strcmp(colordirlabels(i),'L+M')
                col = [0 0 0];
            elseif strcmp(colordirlabels(i),'L-M')
                col = [1 0 1];
            end
            plot(ecc,data(L,i),'wo','MarkerFaceColor',col,'MarkerSize',8);
            plot([2 14]',[1 2; 1 14]*b,'color',col)
            set(gca,'Xlim',[0 15],'Ylim',[-0.5 2.5]); xlabel('Eccentricity (°)'); ylabel('Signal-to-noise ratio (d'')');
        end
    end

    [b,bint, ~, ~, stats] = regress(data(L,2)-data(L,1), [ones(length(ecc),1) ecc]);
    disp([whichcelltype{1},' slope(d''color-d''lum vs ecc): ',num2str(b(2)),' p:',num2str(stats(3))]);
end

% Comparing  d' values of ON and OFF cell (M and P)
% ON and OFF of the basis of sign of L + M cone weights
colordirlabels = {'L+M','L-M'};
L_ON_cells = sign(sum(normalizedconeweights(:,[1 2]),2)) == 1;
for whichcelltype = {'M','P'}
    L = strcmp(cellTypes,whichcelltype);
    for i=1:2
        [h,p] = ttest2(data(L&L_ON_cells,i), data(L&~L_ON_cells,i));
        disp([whichcelltype{1},': ON vs OFF ',colordirlabels{i},', p = ',num2str(p),' ON mean: ',num2str(mean(data(L&L_ON_cells,i))),' OFF mean: ',num2str(mean(data(L&~L_ON_cells,i)))]);
    end
end
 
% No difference between ON and OFF wrt mean d'
%%
% Section 5
% d-prime as a function of temporal frequency for L+M and L-M, averaged
% across neurons

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
examplecellflags = [];
for i = 1:length(filenames)
    idx = neuronids(i);
    tmpfilenames = filenames{neuronids == idx};
    tmpspikenames = spikenames(neuronids == idx,:);
    if any(strcmp(tmpfilenames, EXAMPLEMAGNOCELL)) || any(strcmp(tmpfilenames, EXAMPLEPARVOCELL)) || any(strcmp(tmpfilenames, EXAMPLEKONIOCELL))
        examplecellflags(i) = 1;
    else
        examplecellflags(i) = 0;
    end
    
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

% Plotting trajectories
colors = [0 .5 1; 0 0 0; 1 0 0]; % could be in section 0. K, M, P

figprefs; axes; hold on;
plot([-1 3],[-1 3],'k--')
plot([0 0],[-1 3],'k-')
plot([-1 3],[0 0],'k-')
for i = 1:max(celltypemat)
    L = celltypemat == i;
    mn = nanmean(data(:,:,L),3);
    sd = nanstd(data(:,:,L),[],3);
    n = sum(~isnan(data(:,:,L)),3);
    Ln = n>1; % not plotting data for conditions with only one cell
    L = Ln(1,:) & Ln(2,:);
    plot(mn(1,L),mn(2,L),'-','color',colors(i,:),'linewidth',2);
    for j = 1:size(mn,2)
        if L(j)
            plot(mn(1,j)+sd(1,j)*[-1 1],mn(2,j)*[1 1],'-','color',colors(i,:),'linewidth',0.1);
            plot(mn(1,j)*[1 1],mn(2,j)+sd(2,j)*[-1 1],'-','color',colors(i,:),'linewidth',0.1);
            plot(mn(1,j),mn(2,j),'o','MarkerSize',j/2+2,'MarkerEdgeColor','white','MarkerFaceColor',colors(i,:));
        end
    end
end
xl=get(gca,'Xlim');
yl=get(gca,'Ylim');
set(gca,'Xlim',[min([xl(1), yl(1)]) max([xl(2), yl(2)])]);
set(gca,'Ylim',[min([xl(1), yl(1)]) max([xl(2), yl(2)])]);
set(gca,'Xlim',[-1 3],'Ylim',[-1 3]);
xlabel('Signal-to-noise ratio (d''_{L+M})');
ylabel('Signal-to-noise ratio (d''_{L-M})');

% dot size calibration legend
axes('position',[10.2 3 1 7]);
for j = find(L)
    plot(.5,j,'o','MarkerSize',j/2+2,'MarkerEdgeColor','white','MarkerFaceColor','black');
    text(1,j,num2str(round(10*TFBINCENTERS(j))/10),'FontName','Helvetica','FontAngle','italic')
end
set(gca,'Ytick',[],'Xtick',[],'Box','on');
text(4,10,'Frequency (Hz)','FontName','Helvetica','FontAngle','italic','Rotation',-90);
set(gca,'Visible','off');

%%
% Section 5.1
% Within-neuron analysis of changes in d' with TF. Need to run section 5
% first to get variables "data" , "celltypemat", and "examplecellflags".
% data is [lum, rg] x TFBINCENTERS x neuron

outdata = [];
celltype_idxs = unique(celltypemat);
for i = celltype_idxs'
    L = celltypemat == i;
    for whichcell = find(L)'
        tmp = [];
        for j = 1:2 % Color direction [lum, rg]
            y = data(j,:,whichcell);
            y = y(~isnan(y));
            X = [ones(length(y),1) log10(TFBINCENTERS(1:length(y)))'];
            [b, bint] = regress(y',X);
            tmp = [tmp, b(2) sign(bint(2,1)) == sign(bint(2,2))];
        end
        outdata = [outdata; tmp];
        if i == 2 & j == 2 & b(2) < 0
            keyboard 
            % What's going on with the magnocell with L-M SNR that decreases with TF?
            % It's A060217005. Big non-time locked responses to high
            % frequency L-M results in a negative d' (many spikes, but
            % smaller F1 component than expected by chance.)
        end
    end
end

figprefs;
axes; hold on;
plot([-1 3],[0 0],'k-');
plot([0 0],[-1 3],'k-');
for i = celltype_idxs'
    L = celltypemat == i;
    for j = find(L)'
        h = plot(outdata(j,1),outdata(j,3),'wo','MarkerFaceColor',colors(i,:),'MarkerSize',7);
        if outdata(j,2) == 0 & outdata(j,4) == 0
            set(h,'MarkerFaceColor',colors(i,:)/3+[1 1 1]*(2/3));
        end
        if examplecellflags(j)
            set(h,'Marker','^','MarkerSize',9);
        end
    end
end
set(gca,'Xlim',[-1 3],'Ylim',[-1 3]);
xlabel('slope of d'' vs. log(TF) L+M');
ylabel('slope of d'' vs. log(TF) L-M');

%%
% Section 6
% KNN analysis. Taken from IsoSampPop section 17

Q = 6; % Q = 6, k = 9 is close to optimal. See next section.
k = 9;
MONKEYS={'Apollo','Utu'};
data = cell(2,2);
for monkey = MONKEYS
    for celltype = {'M'}
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
AXWIDTH = 6;
%celltitles = {'magnocellular','parvocellular'};
axpositions = [2 18 AXWIDTH AXWIDTH; 10 18 AXWIDTH AXWIDTH];
colors = [0 0 0; 1 0 1]; % L+M, L-M
for i = 1:2 % Monkey
    for j = 1:1 % Celltype
        tmp = outdata{i,j};
        axes('position',axpositions(j,:)+[0 -8 0 0]*(i-1));
        plot([TFBINCENTERS(1) TFBINCENTERS(end)],[0 0],'k:');
        hold on;
        for k = 1:2 % lum/color
            for m = 1:2 % parametric/nonparametric
                n = sum(~isnan(tmp(:,:,k,m)));
                L = n >= 2; % 2 data points minimum for plotting
                mn = nanmean(tmp(:,:,k,m));
                sd = nanstd(tmp(:,:,k,m));
                sem = sd./sqrt(n);
                h_patch = patch([TFBINCENTERS(L), fliplr(TFBINCENTERS(L))],[mn(L)+sem(L), fliplr(mn(L)-sem(L))],colors(k,:),'Facealpha',.25,'LineStyle','none');
                h_line = plot(TFBINCENTERS(L),nanmean(tmp(:,L,k,m)),'ko-','Color',colors(k,:),'MarkerFaceColor',colors(k,:),'Markersize',4);
                if m == 2 % nonparametric
                   set(h_line,'LineStyle','--','MarkerFaceColor','none'); 
                end
                %if i == 2
                %    set(h_line,'Marker',MONKEY2MARKER); 
                %end
            end
        end
        set(gca,'Xscale','log','XTick',[1 10],'XTickLabel',[1 10]);
        title(['Monkey: ',num2str(i),' ',celltitles{j}]);
        ylabel('Signal-to-noise ratio (d'')');
        xlabel('Frequency (Hz)');
        set(gca,'Xlim',[TFBINCENTERS(1), TFBINCENTERS(end)],'Ylim',[-1 4],'Ytick',[-1:4]);
    end
end

% Trajectory plot
colors = [0 0 0; 1 0 1];
axpositions = [2 2 AXWIDTH AXWIDTH; 10 2 AXWIDTH AXWIDTH];
for j = 1:1 % Celltype
    axes('Position',axpositions(j,:)); hold on;
    plot([-2 4],[-2 4],'k:');
    for i = 1:2 % Monkey
        tmp = outdata{i,j};
        for k = 1:2 % lum/color
            n = sum(~isnan(tmp(:,:,k,m)));
            L = n >= 2; % 2 data points minimum for plotting
            mn = [nanmean(tmp(:,:,k,1)); nanmean(tmp(:,:,k,2))]';
            sd = [nanstd(tmp(:,:,k,1)); nanstd(tmp(:,:,k,1))]';
            sem = sd./repmat(sqrt(n'),1,2);
            h = plot(mn(L,1),mn(L,2),'k-','Color',colors(k,:),'linewidth',2); 
            for l = 1:sum(L)
                plot(mn(l)+sd(l,1)*[-1 1],mn(l,2)*[1 1],'-','color',colors(k,:),'linewidth',0.1);
                plot(mn(l,1)*[1 1],mn(l,2)+sd(l,2)*[-1 1],'-','color',colors(k,:),'linewidth',0.1);
                h = plot(mn(l,1),mn(l,2),'o','MarkerSize',l/2+2,'MarkerEdgeColor','white','MarkerFaceColor',colors(k,:));
                
                if i == 2
                    set(h,'Marker',MONKEY2MARKER);
                end
            end
        end
    end
    xlabel('parametric d''');
    ylabel('non-parametric d''');
    %title([celltitles{j}]);
    set(gca,'Xlim',[-2 4],'Ylim',[-2 4],'Xtick',[-2:2:4],'Ytick',[-2:2:4]);
end

%%
% Section 7
% Population d' as a function of the begining and ending times of the spike
% counting window. Taken from IsoSampPaperStuff section 25, but only
% looking at parvocellular neurons.

starts = linspace(0,.2,6);
stops = .666+starts;
tmp = fullfact([length(starts) length(stops)]);
offsets = [starts(tmp(:,1))' stops(tmp(:,2))'];

data = cell(2,1); % within every cell we have a little data matrix with Lcc, Mcc, TF, offset start, offset stop, and d'
for MONKEY = MONKEYS
    [filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',{'P'},'subjID',{MONKEY{1}(1)});
    for i = 1:length(filenames)
        stro = {};
        for j = 1:length(filenames{i})
            stro{j} = nex2stro(char(findfile(filenames{i}(j), fullfile(nexfilepath,'Greg',MONKEY))));
        end
        stro = strocat(stro);
        
        uniquestim = IsoSampGetDPrime(stro,DPRIMEMETHOD,spikecds(i)); % Just getting uniquestim
        population_scalefactor = IsoSampGetPopulationScaleFactor(stro, ecc_to_diam_deg_P, TEMPORONASALSCALEFACTOR, RFTRUNCATIONINSD, ONOFFCORRELATION, 2);
        
        tmp = [uniquestim, nan*ones(size(uniquestim,1),size(offsets,1))];
        for j = 1:size(offsets,1)
            [uniquestim, dprime] = IsoSampGetDPrime(stro,DPRIMEMETHOD,spikecds(i),offsets(j,:));
            tmp(:,j+3) = dprime*population_scalefactor;
        end
        
        % Adding to the "data" cell array of cell arrays
        listsofar = data{strcmp(MONKEY,MONKEYS),strcmp(CELLTYPE,CELLTYPES)};
        listsofar{length(listsofar)+1} = tmp;
        data{strcmp(MONKEY,MONKEYS)}=listsofar;
    end
end

MONKEYIDX = 'both'; % 1 = Apollo, 2 = Utu
COLORDIR = 'RG';

if strcmp(MONKEYIDX,'both')
    tmp = cat(2,data{1}, data{2});
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
AMARGIN = .5;
figprefs;
for i = 1:size(outdata,1)
    mn = squeeze(nanmean(outdata(i,:,:,:),4));
    if all(all(ns(i,:,:) > 1)) & any(mn(:) < 1.27)
        axes('position',[2+(i-1)*(AXWIDTH+AMARGIN),7, AXWIDTH,AXWIDTH],'Box','on');
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
axes('position',[19.5,7,.3,AXWIDTH],'Yaxislocation','right');
set(gca,'Xtick',[],'YAxisLocation','right','Box','on')
set(gca,'Ytick',[1 size(colormap,1)],'YtickLabel',[minval maxval])
image([1:size(colormap,1)]');
axis tight;
ylabel('Signal-to-noise ratio (d'')');


%%
% Section 7.1
% Looking effects of spike counting window over a broad range for M cell
% (L+M) and P cells (L-M). Answering two questions: 
% 1) is there a window in which populations of P cells are more sensitive than the monkey?
% 2) is there a window in which M cell population sensitivity to L+M and 
% P cell population sensitivity to L-M are very different?

starts = linspace(0,.3,7);
durations = [.05 .1 .2 .4 .66];
offsets = zeros(length(starts)*length(durations),2);
counter = 1;
for i = 1:length(starts)
    for j = 1:length(durations)
        offsets(counter,:) = [starts(i) starts(i)+durations(j)];
        counter = counter+1;
    end
end

data = cell(2,1); % within every cell we have a little data matrix with Lcc, Mcc, TF, offset start, offset stop, and d'
for CELLTYPE = CELLTYPES
    [filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',CELLTYPE);
    for i = 1:length(filenames)
        stro = {};
        for j = 1:length(filenames{i})
            stro{j} = nex2stro(char(findfile(filenames{i}(j), fullfile(nexfilepath,'Greg'))));
        end
        stro = strocat(stro);
        
        uniquestim = IsoSampGetDPrime(stro,DPRIMEMETHOD,spikecds(i)); % Just getting uniquestim
        if strcmp(CELLTYPE, 'M')
            ecc_to_diam_deg = ecc_to_diam_deg_M;
            disp('M');
        else
            ecc_to_diam_deg = ecc_to_diam_deg_P;
            disp('P');
        end
        
        population_scalefactor = IsoSampGetPopulationScaleFactor(stro, ecc_to_diam_deg, TEMPORONASALSCALEFACTOR, RFTRUNCATIONINSD, ONOFFCORRELATION, 2);
        
        tmp = [uniquestim, nan*ones(size(uniquestim,1),size(offsets,1))];
        for j = 1:size(offsets,1)
            [uniquestim, dprime] = IsoSampGetDPrime(stro,DPRIMEMETHOD,spikecds(i),offsets(j,:));
            tmp(:,j+3) = dprime*population_scalefactor;
        end
        
        % Adding to the "data" cell array of cell arrays
        listsofar = data{strcmp(CELLTYPE,CELLTYPES)};
        listsofar{length(listsofar)+1} = tmp;
        data{strcmp(CELLTYPE,CELLTYPES)}=listsofar;
    end
end

% Plotting one figure per counting window
figprefs;
nplotrows = length(starts);
nplotcols = length(durations);
AXWIDTH = 2.2;
AXMARGIN = .9;
axcols = ([1:nplotcols]-1).*AXWIDTH+([1:nplotcols]*AXMARGIN)+2;
axrows = fliplr(([1:nplotrows]-1).*AXWIDTH+([1:nplotrows]*AXMARGIN)+2);

for i = 1:nplotrows
    for j = 1:nplotcols
        axes('Position',[axcols(j) axrows(i) AXWIDTH AXWIDTH]); hold on;
        for celltypeidx = 1:2 % 1 = M, 2 = P
            dprimes = nan(length(data{celltypeidx}), length(TFBINCENTERS));
            for cellcounter = 1:length(data{celltypeidx})
                tmp_all = data{celltypeidx}{cellcounter};
                uniquestim = tmp_all(:,1:3);
                tmp_popndprimes = tmp_all(:,4:end);
                if celltypeidx == 1
                    Lcond = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & uniquestim(:,1) ~= 0;
                else
                    Lcond = sign(uniquestim(:,1)) ~= sign(uniquestim(:,2));
                end
                Lwindow = offsets(:,1) == starts(i) & offsets(:,2) == starts(i)+durations(j);
                dprimes(cellcounter, 1:sum(Lcond)) = tmp_popndprimes(Lcond,Lwindow);
            end
            if celltypeidx == 1
                color = [0 0 0];
            else
                color = [1 0 1];
            end
            mn = nanmean(dprimes);
            sd = nanstd(dprimes);
            n = sum(~isnan(dprimes));
            L = ~isnan(sd) & n > 1;
            TBC = TFBINCENTERS(L);
            mn = mn(L);
            sem = sd(L)./sqrt(n(L));
            patch([TBC, fliplr(TBC)],[mn+sem, fliplr(mn-sem)],color,'Facealpha',.5,'LineStyle','none');
            plot(TBC, mn,'o-','LineWidth',1,'MarkerEdgeColor','white','MarkerFaceColor',color,'Color',color);
        end
        set(gca,'Xlim',[TFBINCENTERS(1) TFBINCENTERS(end)],'Xtick',[1 10],'XtickLabel',{'1', '10'});
        plot([TFBINCENTERS(1) TFBINCENTERS(end)],[1.27 1.27],'k-'); % Psychophysical d'
        plot([TFBINCENTERS(1) TFBINCENTERS(end)],[0 0],'k-','color',[.5 .5 .5]); % No signal
        
        set(gca,'Xscale','log','FontSize',8);
        if i == 1
            h = title([num2str(durations(j)*1000),' ms '],'FontWeight','normal','FontSize',10);
        end
        if j == 1
            ylabel('Population d''','FontSize',10);
        end
        if j == nplotcols
            h = text(150,5,[num2str(starts(i)*1000),' ms'],'Color',[.5 .5 .5],'rotation',270,...
                'HorizontalAlignment','center','FontAngle','italic','FontSize',10);
        end
        if i == nplotrows
            xlabel('Frequency (Hz)','FontSize',10);
        end
        set(gca,'Ylim',[-2 12]);
    end
end

%%
% Section 8
% Population sensitivity analysis without cones.
% M and P, both color directions, averaged across monkeys

% Optional spike counting offset
OFFSETS = [.12 .786]; % adding an offset to the counting window. This only affects LGN idlObs.
OFFSETS = [];

data = cell(2,2); % [Magno/Parvo] x [L+M/L-M] x [monkey]
for MONKEY = MONKEYS % doesn't hurt to segregate by monkey
    for CELLTYPE = CELLTYPES
        if strcmp(CELLTYPE,'M')
            ecc_to_diam_deg = ecc_to_diam_deg_M;
        else
            ecc_to_diam_deg = ecc_to_diam_deg_P;
        end
        [filenames,spikecds] = fnamesFromTxt('IsoSamp_LGN','cellClass',CELLTYPE,'subjID',{MONKEY{1}(1)});
        for i = 1:length(filenames)
            tmp = [];
            cal = [];
            stro = {};
            for j = 1:length(filenames{i})
                stro{j} = nex2stro(char(findfile(filenames{i}(j), fullfile(nexfilepath,'Greg',MONKEY))));
            end
            stro = strocat(stro);
            [tmp.uniquestim, tmp.dprime] = IsoSampGetDPrime(stro,DPRIMEMETHOD,spikecds(i),OFFSETS); % GDLH 12/12/19 Added OFFSET to LGN only!
            [tmp.population_scalefactor, tmp.ncells] = IsoSampGetPopulationScaleFactor(stro, ecc_to_diam_deg, TEMPORONASALSCALEFACTOR,RFTRUNCATIONINSD, ONOFFCORRELATION, 2);
            listsofar = data{strcmp(MONKEY,MONKEYS),strcmp(CELLTYPE,CELLTYPES)};
            listsofar{length(listsofar)+1} = tmp;
            data{strcmp(MONKEY,MONKEYS),strcmp(CELLTYPE,CELLTYPES)}=listsofar;
        end
    end
end

% Doing the plotting
AXWIDTH = 5;
TFSFORINSET = TFBINCENTERS(TFBINCENTERS < 3);
figprefs;
hax(1) = axes('position',[2 2 AXWIDTH AXWIDTH]); hold on;
hax(2) = axes('position',[2+AXWIDTH+1 2 AXWIDTH AXWIDTH]); hold on;
hforleg = [];
for celltype_idx = 1:2
    celltype_data = cat(2,data{1,celltype_idx},data{2,celltype_idx});
    for COLORDIR = {'RG','LUM'}
        tmp = nan*ones(2,length(TFBINCENTERS),length(celltype_data)); % rows: LGN dprime, ns
        for i = 1:length(celltype_data)
            cell_data = celltype_data{i};
            Lblank = cell_data.uniquestim(:,1) == 0 & cell_data.uniquestim(:,2) == 0 & cell_data.uniquestim(:,3) == 0;
            Llum = sign(cell_data.uniquestim(:,1)) == sign(cell_data.uniquestim(:,2)) & ~Lblank;
            Lrg = sign(cell_data.uniquestim(:,1)) ~= sign(cell_data.uniquestim(:,2)) & ~Lblank;
            
            if strcmp(COLORDIR,'RG')
                LlumORrg = Lrg;
                color = [1 0 1];
            else
                LlumORrg = Llum;
                color = [0 0 0];
            end
            if celltype_idx == 1
                plotstyle = '-^';
            else
                plotstyle = '-o';
            end
            for j = 1:length(TFBINCENTERS)
                Ltf = cell_data.uniquestim(:,3) > TFBINEDGES(j) & cell_data.uniquestim(:,3) <= TFBINEDGES(j+1);
                if sum(Ltf&LlumORrg) > 0
                    tmp(1,j,i) = mean(cell_data.dprime(LlumORrg&Ltf))*cell_data.population_scalefactor;
                    tmp(2,j,i) = 1; % each cell counts as an independent entity
                end
            end
        end
        n = nansum(tmp(2,:,:),3);
        % Getting rid of TFs with 1 or fewer data points
        tmp = tmp(:,n>1,:);
        n = n(1:size(tmp,2));
        
        mn_lgn = nanmean(tmp(1,:,:),3);
        sd_lgn = nanstd(tmp(1,:,:),0,3);
        sem_lgn = sd_lgn./sqrt(n);
        TBC = TFBINCENTERS(1:length(n));
        axes(hax(1));
        patch([TBC, fliplr(TBC)],[mn_lgn+sem_lgn, fliplr(mn_lgn-sem_lgn)],color,'Facealpha',.5,'LineStyle','none');
        h = plot(TBC, mn_lgn,plotstyle,'LineWidth',1,'MarkerEdgeColor','black','MarkerFaceColor',color,'Color',color);
        hforleg = [hforleg; h];
        
        axes(hax(2)); % Inset
        L = ismember(TFBINCENTERS,TFSFORINSET);
        patch([TFSFORINSET, fliplr(TFSFORINSET)],[mn_lgn(L)+sem_lgn(L), fliplr(mn_lgn(L)-sem_lgn(L))],color,'Facealpha',.5,'LineStyle','none');
        plot(TFSFORINSET, mn_lgn(L),plotstyle,'LineWidth',1,'MarkerEdgeColor','black','MarkerFaceColor',color,'Color',color);
  
    end
end
axes(hax(2));
set(gca,'Xscale','log','Xlim',[TFSFORINSET(1)*.9 TFSFORINSET(end)],'Xtick',[1 2 3],'Xticklabel',[1 2 3],'Ylim',[-1 4],'Yticklabel',[-1:4])
xlabel('Frequency (Hz)');
ylabel('Signal-to-noise ratio (d'')');

axes(hax(1));
set(gca,'Xscale','log','Xlim',[TFBINCENTERS(1)*.9 TFBINCENTERS(end)*1.1],'Ylim',[-2 12],'Xtick',[1 10],'XtickLabel',[1 10]);
xlabel('Frequency (Hz)');
if monkey_idx == 1
    ylabel('Signal-to-noise ratio (d'')');
end
title(['offset = ',num2str(OFFSETS)]);
xlims = get(hax(2),'Xlim');
ylims = get(hax(2),'Ylim');
plot(xlims([1,2,2,1,1]),ylims([1,1,2,2,1]),'k:','LineWidth',2); % inset box
hleg = legend(hforleg,{'Magno L-M','Magno L+M','Parvo L-M','Parvo L+M'},'Location','NorthWest');

for i = 1:2
    axes(hax(i));
    plot([TFBINEDGES(1) TFBINEDGES(end)],[1.27 1.27],'k--');
end
set(gcf,'Renderer','painters');

%%
% Section 9
% Comparing population SNR for magnocellular neurons (L+M) and
% parvocellular neurons (L-M) to show that the downstream loss is similar
% but the loss from cones to LGN is different.
% Largely taken from IsoSampPaperStuff.m Section 6.
% 5/11/20 option to turn off the population scaling

POPULATIONSCALING = false;
 
data = cell(2,2); % [Magno L+M and Parvo L-M] x [monkey]
for MONKEY = MONKEYS % doesn't hurt to segregate by monkey
    for CELLTYPE = CELLTYPES
        if strcmp(CELLTYPE,'M')
            ecc_to_diam_deg = ecc_to_diam_deg_M;
            COLORDIR = 'LUM';
            OFFSETS = [];
        else
            ecc_to_diam_deg = ecc_to_diam_deg_P;
            COLORDIR = 'RG';
            OFFSETS = [.12 .786]; % adding an offset to the counting window. This only affects LGN idlObs.
        end
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
            Lrg = sign(tmp.uniquestim(:,1)) ~= sign(tmp.uniquestim(:,2)) & ~Lblank;

            % Now getting population scale factor
            % Getting RF size at rfx, rfy
            if POPULATIONSCALING
                [tmp.population_scalefactor, tmp.ncells] = IsoSampGetPopulationScaleFactor(stro, ecc_to_diam_deg, TEMPORONASALSCALEFACTOR,RFTRUNCATIONINSD, ONOFFCORRELATION, 2);
            else
                tmp.population_scalefactor = 1;
                tmp.ncells = 1;
            end
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
                rf_sigma = ecc_to_diam_deg(sqrt((stro.sum.exptParams.rf_x/10)^2+(stro.sum.exptParams.rf_y/10)^2))/2;
                params.gab.sd = sqrt((gabor_sigma^-2+rf_sigma^-2)^-1);
                params.gab.sd = rf_sigma;
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
            listsofar = data{strcmp(MONKEY,MONKEYS),strcmp(CELLTYPE,CELLTYPES)};
            listsofar{length(listsofar)+1} = tmp;
            data{strcmp(MONKEY,MONKEYS),strcmp(CELLTYPE,CELLTYPES)}=listsofar;
        end
    end
end

% Doing the plotting
%figure('Position',[440 100 750 700],'DefaultAxesTickDirMode','manual','DefaultAxesTickdir','out','DefaultAxesYcolor','black','DefaultAxesXcolor','black')
%set(gcf,'DefaultAxesFontSize',15,'DefaultAxesFontAngle','italic','DefaultAxesUnits','centimeters');

AXWIDTH = 5;
AXWIDTH = 9;
snr_diff_LUM = [];
snr_diff_RG = [];
HATCHDENSITY = 60;
figprefs;
for monkey_idx = 1:2
    hax = axes('position',[6 15-12*(monkey_idx-1) AXWIDTH AXWIDTH]); hold on;
    set(hax,'Xscale','log','Xlim',[TFBINCENTERS(1)*.9 TFBINCENTERS(end)*1.1],'Xtick',[1 10],'XtickLabel',[1 10]);
    if POPULATIONSCALING
        set(hax,'Ylim',[-2 35]);
        plot([TFBINEDGES(1) TFBINEDGES(end)],[1.27 1.27],'k--');
    else
        set(hax,'Ylim',[-1 8]);
    end

    for celltype_idx = 1:2
        if celltype_idx == 1
            COLORDIR = 'LUM';
            color = [0 0 0];
        else
            COLORDIR = 'RG';
            color = [1 0 1];
        end
        %    celltype_data = cat(2,data{1,celltype_idx},data{2,celltype_idx});
        celltype_data = data{monkey_idx,celltype_idx};
        tmp = nan*ones(3,length(TFBINCENTERS),length(celltype_data)); % rows: LGN dprime, cone dprime, ns
        for i = 1:length(celltype_data)
            cell_data = celltype_data{i};
            Lblank = cell_data.uniquestim(:,1) == 0 & cell_data.uniquestim(:,2) == 0 & cell_data.uniquestim(:,3) == 0;
            Llum = sign(cell_data.uniquestim(:,1)) == sign(cell_data.uniquestim(:,2)) & ~Lblank;
            Lrg = sign(cell_data.uniquestim(:,1)) ~= sign(cell_data.uniquestim(:,2)) & ~Lblank;
            
            if strcmp(COLORDIR,'RG')
                LlumORrg = Lrg;
                hatch_angle = 45;
            else
                LlumORrg = Llum;
                hatch_angle = 135;
            end
            for j = 1:length(TFBINCENTERS)
                Ltf = cell_data.uniquestim(:,3) > TFBINEDGES(j) & cell_data.uniquestim(:,3) <= TFBINEDGES(j+1);
                if sum(Ltf&LlumORrg) > 0
                    tmp(1,j,i) = mean(cell_data.dprime(LlumORrg&Ltf))*cell_data.population_scalefactor;
                    tmp(2,j,i) = mean(cell_data.conedprime(LlumORrg&Ltf));
                    tmp(3,j,i) = 1; % each cell counts as an independent entity
                end
            end
        end
        n = nansum(tmp(3,:,:),3);
        % Getting rid of TFs with 1 or fewer data points
        tmp = tmp(:,n>1,:);
        n = n(1:size(tmp,2));
        
        mn_lgn = nanmean(tmp(1,:,:),3);
        sd_lgn = nanstd(tmp(1,:,:),0,3);
        sem_lgn = sd_lgn./sqrt(n);
        mn_cone = nanmean(tmp(2,:,:),3);
        sd_cone = nanstd(tmp(2,:,:),0,3);
        sem_cone = sd_cone./sqrt(n);
        TBC = TFBINCENTERS(1:length(n));
        
        % Patch highlighting the gap between cones and lgn
        h_patch = patch([TBC, fliplr(TBC)],[mn_cone, fliplr(mn_lgn)],color,'FaceAlpha',0,'LineStyle','none');
        hatchfill2(h_patch,'Single','HatchAngle',hatch_angle,'HatchDensity',HATCHDENSITY,'HatchLineWidth',.1,'HatchColor','black');
        % Patch highlighting the gap between lgn and behavior
        if POPULATIONSCALING
            h_patch = patch([TBC, fliplr(TBC)],[mn_lgn, repmat(1.27,1,length(mn_lgn))],color,'FaceAlpha',0,'LineStyle','none');
            hatchfill2(h_patch,'Single','HatchAngle',hatch_angle+45,'HatchDensity',HATCHDENSITY,'HatchLineWidth',.1,'HatchColor','black');
        end
        
        patch([TBC, fliplr(TBC)],[mn_cone+sem_cone, fliplr(mn_cone-sem_cone)],color,'Facealpha',.5,'LineStyle','none');
        plot(TBC, mn_cone,'--^','LineWidth',1,'MarkerEdgeColor','black','MarkerFaceColor',color,'Color',color);
        patch([TBC, fliplr(TBC)],[mn_lgn+sem_lgn, fliplr(mn_lgn-sem_lgn)],color,'Facealpha',.5,'LineStyle','none');
        h = plot(TBC, mn_lgn,'-o','LineWidth',1,'MarkerEdgeColor','black','MarkerFaceColor',color,'Color',color);
        
        % Integrating differences in d' between L+M and L-M
        if celltype_idx == 1
            snr_diff_LUM = [mn_cone-mn_lgn;mn_lgn-1.27];
        else
            snr_diff_RG = [mn_cone-mn_lgn;mn_lgn-1.27];
        end   
    end
    maxTFidx = min(size(snr_diff_LUM,2),size(snr_diff_RG,2));
    diffs = snr_diff_LUM(:,1:maxTFidx)-snr_diff_RG(:,1:maxTFidx); % {cortex, eye} x TF. Each element is LUM-RG 
    %integrated_diffs = sqrt(sum(diffs.^2,2))
    integrated_diffs = sum(abs(diffs),2)
    
    integrated_diffs(1)/integrated_diffs(2)
    
    xlabel('Frequency (Hz)');
    if monkey_idx == 1
        ylabel('Signal-to-noise ratio (d'')');
    end
    title(['Monkey ',num2str(monkey_idx)]); 
       
end
set(gcf,'Renderer','painters');



%%
% Section 10
% Behavioral contrast-sensitivity
w_lum = logspace(log10(1),log10(60),200);
w_rg = w_lum(w_lum<=40);

figprefs; 
axes('position',[4.5,8,8,8]); hold on;
for monkey_idx = 1:length(MONKEYS)
    [filenames, ~, spikenames, neuronids] = fnamesFromTxt('IsoSamp_LGN','subjID',{MONKEYS{monkey_idx}(1)},'cellClass',{'M','P','K'});
    RFs = [];
    for i = 1:length(filenames)
        idx = neuronids(i);
        tmpfilenames = filenames{neuronids == idx};
        tmpspikename = spikenames(neuronids == idx,:);
        stro=[];
        for j = 1:length(tmpfilenames)
            stro = strocat(stro,nex2stro(findfile(tmpfilenames{j})));
        end
        RFs = [RFs;stro.sum.exptParams.rf_x/10 stro.sum.exptParams.rf_y/10];
    end
    % Getting model
    isosamppath = which('IsoSampOnline');
    isosamppath(find(isosamppath==filesep,1,'last'):end) = [];
    load ([isosamppath,filesep,'private',filesep','data',filesep,'LMTF.mat']);
    modelstruct = eval(MONKEYS{monkey_idx}(1));
    
    lumdata = [];
    rgdata = [];
    for i = 1:size(RFs,1)
        behavioralmodel = LMTF_global_to_local_model(modelstruct.legacy.mode5params, RFs(i,1), RFs(i,2), 5);
        lum_tcsf = @(omega)(behavioralmodel(1)*abs(((1i*2*pi*10^behavioralmodel(5).*omega+1).^-behavioralmodel(3))-behavioralmodel(2)*((1i*2*pi*10^(behavioralmodel(5)+behavioralmodel(6)).*omega+1).^-(behavioralmodel(3)+behavioralmodel(4)))));
        rg_tcsf = @(omega)(behavioralmodel(7)*abs(((1i*2*pi*10^behavioralmodel(11).*omega+1).^-behavioralmodel(9))-behavioralmodel(8)*((1i*2*pi*10^(behavioralmodel(11)+behavioralmodel(12)).*omega+1).^-(behavioralmodel(9)+behavioralmodel(10)))));
        lumdata = [lumdata; lum_tcsf(w_lum)];
        rgdata = [rgdata; rg_tcsf(w_rg)];
    end
    colors = [0 0 0; 1 0 1];
    for i = 1:2
        if i == 1
            data = lumdata;
            w = w_lum;
        else
            data = rgdata;
            w = w_rg;
        end
        sd = std(log10(data));
        mn = mean(log10(data));
        patch([w fliplr(w)], 10.^[mn+sd, fliplr(mn-sd)],colors(i,:),'FaceAlpha',.25,'EdgeColor','none')
        h = plot(w, geomean(data),'-','LineWidth',2,'Color',colors(i,:));
        if monkey_idx == 2
            set(h,'LineStyle','--')
        end
    end
end
set(gcf,'Renderer','Painters');
set(gca,'Yscale','log','Xscale','log','Xlim',[min(w_lum) max(w_lum)],'Ylim',[.5 100]);
set(gca,'Xtick',[1 10],'Xticklabel',[1 10]);
set(gca,'ytick',[.1 1 10 100],'Yticklabel',[.1 1 10 100]);
xlabel('Frequency (Hz)');
ylabel('Contrast sensitivity');
