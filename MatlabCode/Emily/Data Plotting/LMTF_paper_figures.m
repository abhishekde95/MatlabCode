%LMTF model fits for paper
%LMTF plane is created separately in illustator.
%% Section 1:
%green fit (traditional funnel)
% See also LMTFpop.m Section 10
SAVEFIGS = 1; % and erase other figs
if SAVEFIGS
    close all
end
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
flist = fetch(conn, 'SELECT fileID FROM LMTF WHERE subjID=''E'' AND monID =''ProPixx'' AND recDate > ''2016-05-27'' AND recDate < ''2017-02-03''AND quality = 1 AND rfX = 50 and rfY = 0');
close(conn);
[~, u_data] = iterateAndPlotFiles_modularPlusDB(flist, 1);
figure(1); set(gcf,'renderer','painters');
if SAVEFIGS
    set(gca,'View',[45,12],'Visible','off','Color','none');
    title([]); export_fig junk1 -png -transparent
    set(gca,'View',[135,12],'Visible','off','Color','none');
    title([]); export_fig junk2 -png -transparent
    
    figure(2); set(gcf,'renderer','painters');
    set(gca,'View',[45,12],'Visible','off','Color','none');
    title([]); export_fig junk3 -png -transparent
    set(gca,'View',[-45,12],'Visible','off','Color','none');
    title([]); export_fig junk4 -png -transparent
end
%% Section 2 Gaussian Process Regression Fit
%not currently used in paper
[u_data, ~, bkgndrgb, MM] = getLMTFrawdata(flist);
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
flist = fetch(conn, 'CALL postPropixxFilenames(''U'')');
close(conn);
fulldata = getLMTFrawdata(flist);
gaussianProcessRegressionFitLMTF(fulldata,0);
gaussianProcessRegressionFitLMTF(fulldata,1);
%% Section 3 mushroom farm DO NOT USE
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
flist = fetch(conn, 'CALL utuEdited');
close(conn);
fulldata = getLMTFrawdata(flist);
lmtfMushroomPlot(fulldata);
%print -dpsc filename -bestfit
%% Section 4 11+2n and meshgrids
%lmtfModelClouds('U', 'firstroundmodels',2); % For analysis
lmtfModelClouds('U', 'mode0models',4); 
% v UNNECESSARY
% For paper figure Make sure to use the same TF as section 5, below (6 Hz)
% ^ UNNECESSARY
%% Section 5 DO NOT USE - DEPRECATED
%11+2n semicircle (commented code was used to create meshgrids)
load(fullfile(fileparts(which('IsoSampOnline.m')), 'private', 'data', 'LMTF.mat'))
TF = 6;
params = U.legacy.mode5params;
b_lum = params(12:15);
b_rg = params([16 17 14 18]);
[x,y] = meshgrid(linspace(0,10,25),linspace(-10,10,60));
[phi,r] = cart2pol(x,y);
X = [ones(numel(r),1), r(:), r(:).*sin(2.*phi(:)), r(:).*cos(2.*phi(:))];
xi_lum = reshape(X*b_lum, size(r)); % This is actually the xi_lum parameter
%[~,~,tmp] = tf_fiterr2([1; params(1:5); 1; params(6:11)],[ones(100,2) linspace(0,15,100)' zeros(100,1)]); % Height of fitted LUM TCSF at "TF" Hz assuming xi_lum = 1
%XiLUMScaleFactor = 1./min(tmp); % sensitivity peak, setting xi_lum to 1. Third output argument from ft_fiterr2 is predicted *threshold*.
[~,~,tmp] = tf_fiterr2([1; params(1:5); 1; params(6:11)],[1 1 TF 0]); % Height of fitted LUM TCSF at "TF" Hz assuming xi_lum = 1
XiLUMScaleFactor = 1/tmp; % sensitivity peak, setting xi_lum to 1. Third output argument from ft_fiterr2 is predicted *threshold*.

xi_rg = reshape(X*b_rg, size(r)); % This is actually the xi_rg parameter
%[~,~,tmp] = tf_fiterr2([1; params(1:5); 1; params(6:11)],[ones(100,1)*[1 -1] linspace(0,15,100)' zeros(100,1)]); % Height of fitted RG TCSF at "TF" Hz assuming xi_rg = 1
%XiRGScaleFactor = 1./min(tmp); % sensitivity peak, setting xi_rg to 1.
[~,~,tmp] = tf_fiterr2([1; params(1:5); 1; params(6:11)],[1 -1 TF 0]); % Height of fitted RG TCSF at "TF" Hz assuming xi_rg = 1
XiRGScaleFactor = 1/tmp; % sensitivity at 5 Hz, setting xi_rg to 1.

log_lum_sens = xi_lum*XiLUMScaleFactor;
log_rg_sens = xi_rg*XiRGScaleFactor;
% Both LUM and RG scaled identically
figure;
h = [];
for i = 1:2
    subplot(1,2,i);
    if (i == 1)
        h(i) = surf(x,y,10.^log_lum_sens);
    else
        h(i) = surf(x,y,10.^log_rg_sens);
    end
    axis tight;
    %axis equal;
    set(h,'EdgeAlpha',0);
    set(gca,'View',[0 90]); 
    lims = [min(10.^[log_lum_sens(:); log_rg_sens(:)]) max(10.^[log_lum_sens(:); log_rg_sens(:)])];
    caxis(lims);
end
set(gcf,'Renderer','painters')
figure; axis; 
h = colorbar;
ticks = 0:5:25;
ticks = ticks(ticks>lims(1) & ticks<lims(2));
set(h,'ticks',ticks/lims(2),'ticklabels',ticks);
%% Section 6 comparison histos parfor
%subj = {'U'};
subj = {'E'};
CVErrs = {size(subj)};
for s = 1:length(subj)
    niter = 600;
    conn = database('Nex_Paradigm_Sort','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
    file_query = sprintf('CALL postPropixxFilenames(''%s'')', subj{s});
    filenames = fetch(conn, file_query);
    close(conn);
    data = getLMTFrawdata(filenames);
    CVerrors = zeros(niter, 8);
    nuniqueXYs = size(unique(data(:,[5 6]),'rows','stable'),1);
    NOOGidxs = Shuffle(find(data(:,4) == 0)); % "Not out of gamut indexes". Don't leave out an OOG point during CV; these rarely change the error
    WARNINGS = {};
    niter = min(niter,length(NOOGidxs));
    fvs = zeros(niter, 8);
    L = false(size(data,1),1);
    %L(NOOGidxs()) = true;
    [input_struct] = LMTF_generate_module_data(data(~L,:), false, {}); %first run just for creating the initial guess input struct
    parfor i=1:niter
        L_inner = L;
        L_inner(NOOGidxs(i)) = true;
        [output_struct, ~] = LMTF_generate_module_data(data(~L_inner,:), false, input_struct); % Here's the workhorse
        %WARNINGS{length(WARNINGS)+1} = warning;
        LXY = all(output_struct.eccs == repmat(data(L_inner,[5 6]),size(output_struct.eccs,1),1),2);
        if sum(LXY) == 0 % If there is not enough data to fit the model at this location, skip this data point
            keyboard;
            CVerrors(i,:) = nan*ones(1,8);
            fvs(i,:) = nan*ones(1,8);
        else
            fvs(i,:) = [sum(output_struct.legacy.firstroundfvs),...
                sum(output_struct.legacy.mode0fvs),...
                sum(output_struct.legacy.mode1fvs),...
                sum(output_struct.legacy.mode1p1fvs),...
                sum(output_struct.legacy.mode2fvs),...
                sum(output_struct.legacy.mode3fvs),...
                sum(output_struct.legacy.mode4fvs),...
                sum(output_struct.legacy.mode5fvs)];
            CVerrors(i,:) = [tf_fiterr2(output_struct.legacy.firstroundmodels(:,LXY),data(L_inner,:)),...
                tf_fiterr2(output_struct.legacy.mode0models(:,LXY),data(L_inner,:)),...
                tf_fiterr2(output_struct.legacy.mode1models(:,LXY),data(L_inner,:)),...
                tf_fiterr2(output_struct.legacy.mode1p1models(:,LXY),data(L_inner,:)),...
                tf_fiterr2(output_struct.legacy.mode2models(:,LXY),data(L_inner,:)),...
                tf_fiterr2(output_struct.legacy.mode3models(:,LXY),data(L_inner,:)),...
                tf_fiterr2(output_struct.legacy.mode4models(:,LXY),data(L_inner,:)),...
                tf_fiterr2(output_struct.legacy.mode5models(:,LXY),data(L_inner,:))];
        end
    end
    CVerrors(all(CVerrors==0,2),:)=[];
    CVErrs{s} = CVerrors;
    outfilename = [subj{s},'_CVout'];
    eval(['save ',outfilename,' fvs CVerrors']);
end
boxAndWhisker(CVErrs);
%% Section 7 LMTF progression
%this is just the gabor patch. Using only Greg's code (below).
isostim = [0.09 -0.09 0; .7 .7 0];

M = [   0.0608    0.1219    0.0175
    0.0220    0.1266    0.0257
    0.0019    0.0095    0.0976];
bkgndrgb = [.5 .5 .5];
bkgndlms = M*bkgndrgb';

edgergb = [0 0 0];
theta = 0;
lambda = 2;
sigma = .5;
gamma = 1;
phi = 0;
xoff = 0;
yoff = 0;
etheta = 0;
edisp = 0;
gausslim = .999;
pixperdeg = 20;

figure;
for i = 1:size(isostim,1)
    gaborrgb = inv(M)*(bkgndlms.*(1+isostim(i,:))') - bkgndrgb;
    im = DrawGaborEdge(bkgndrgb, gaborrgb, edgergb, theta, lambda, sigma, gamma, phi, xoff, yoff, etheta, edisp, gausslim, pixperdeg);
    subplot(ceil(sqrt(size(isostim,1))),ceil(sqrt(size(isostim,1))),i)
    image(im);
    set(gca,'Visible','off')
    axis image
end
%% SECTION 8 RASTERS (FROM GREG)
% Figure showing rasters from one IsoSamp data file
stro = nex2stro(findfile('A061917002.nex'));
Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_l'));
Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_m'));
TF = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'rew_t'));

dur = mean(stimoff_t-stimon_t);
spikeidx = strcmp(stro.sum.rasterCells(1,:),getSpikenum(stro));
spikes = stro.ras(:,spikeidx);
uniquestim = sortrows(unique([Lcc Mcc TF],'rows'),3); % sorting by TF

bins = linspace(0,dur,500);
data = zeros(size(uniquestim,1), length(bins));
temporalenvelope = ones(size(bins));
temporalenvelope(1:round((1/3)*length(bins))) = linspace(0,1,round((1/3)*length(bins)));
temporalenvelope(end:-1:round((2/3)*length(bins))+1) = linspace(0,1,round((1/3)*length(bins)));

for i = 1:size(uniquestim,1)
    contrast = sqrt(uniquestim(i,1).^2+uniquestim(i,2).^2);
    data(i,:) = contrast*temporalenvelope.*sin(2*pi*uniquestim(i,3)*bins);
end
Lstimtype = sign(uniquestim(:,1)) == sign(uniquestim(:,2));

figure;
subplot(1,2,2);
im = 255*(data - min(data(:)))/(max(data(:)) - min(data(:)));
image(flipud(im(Lstimtype,:)));
xtick = 0:.2:dur;
set(gca,'XTick',interp1([bins(1) bins(end)],[1 length(bins)],xtick));
set(gca,'Xticklabel',xtick);
set(gca,'YTick',[]);
colormap(jet(255));
xlabel('time (s)');

% Plotting rasters
offset = [-.1 .1];  % pre and post time wrt stimon/stimoff
subplot(1,2,1); hold on; counter = 0;
for j = find(Lstimtype)'
    L = Lcc == uniquestim(j,1) & Mcc == uniquestim(j,2) & TF == uniquestim(j,3);
    for i = find(L)'
        tmpspikes = spikes{i}-stimon_t(i);
        tmpspikes(tmpspikes < offset(1) | tmpspikes > dur+offset(2)) = [];
        nspikestot = length(tmpspikes);
        plot([tmpspikes tmpspikes]',[zeros(nspikestot,1) 1*ones(nspikestot,1)]'+counter,'k-','linewidth',1);
        counter = counter + 1;
    end
end
set(gca,'Xlim',[0+offset(1) dur+offset(2)],'Ylim',[0 counter],'Ytick',[],'Xtick',[0 .2 .4 .6 .8],'Box','on');
xlabel('time (s)');
%% Section 9 CVErrors
try
    A = load('C:\Users\emily.gelfand\Desktop\MatlabCode\trunk\Emily\A_CVout.mat');
    E = load('C:\Users\emily.gelfand\Desktop\MatlabCode\trunk\Emily\E_CVout.mat');
    U = load('C:\Users\emily.gelfand\Desktop\MatlabCode\trunk\Emily\U_CVout.mat');
    G = load('C:\Users\emily.gelfand\Desktop\MatlabCode\trunk\Emily\G_CVout.mat');
catch
    A = load('/Users/horwitzlab/Desktop/MatlabCode/Emily/A_CVout.mat');
    E = load('/Users/horwitzlab/Desktop/MatlabCode/Emily/E_CVout.mat');
    U = load('/Users/horwitzlab/Desktop/MatlabCode/Emily/U_CVout.mat');
    G = load('/Users/horwitzlab/Desktop/MatlabCode/Emily/G_CVout.mat');
end
subjs = {A, E, U, G};
subj_label = {'A', 'E', 'U', 'G'};
color = {[.37 .24 .6], [.90 .38 .003], [.7 .67 .82], [.99 .72 .39]};
f = figure; axes; hold on;
for s = 1:length(subjs)
    nboot = 200;
    S = subjs{s};
    CVerrors = S.CVerrors;
    fvs = S.fvs;
    CVerrors(all(isnan(CVerrors),2),:) = []; % get rid of nans
    fvs(all(isnan(fvs),2),:) = []; % get rid of nans
    data = [];
    %s_xvals = [s, s+4, s+8];
    for i = 1:3
        m = median(log10(CVerrors(:,end)./CVerrors(:,end-i)));
        btstrp = zeros(nboot,1);
        for j = 1:nboot
            idx = unidrnd(size(CVerrors,1),size(CVerrors,1),1);
            btstrp(j) = median(log10(CVerrors(idx,end)./CVerrors(idx,end-i)));
        end
        subplot(1, 3, i);
        plot(s,m,'s', 'markeredgecolor', color{s}, 'markerfacecolor', color{s});
        text(s,m, subj_label{s}, 'horizontalalignment', 'right', 'color', color{s});
        errorbar(s,m,sqrt(var(btstrp)), 'color', color{s});
    end
xlabelopts = {'yoked v unconstrained','yoked v. lum. only','yoked v. symmetric'};
for i = 1:3
    subplot(1, 3, i);
    plot([0 5],[0 0],'k:');
    set(gca,'Xtick',[], 'ticklength', [0 0], 'xlim', [0 5]);
    set(gca, 'TickDir', 'out');
    set(gca, 'ylim', [-.06 .01]);
    ylabel('log10 error ratio');
    xlabel(xlabelopts{i});
    pbaspect([1 2 1]);
end
end
%% Section 10 sensitivity comparison monkey vs human
senscomp_mvh();
%% Section 11 Rasters take 2
filename = 'A061517005.nex'; % example magnocell
spikeNum = 'sig001a';
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
uniquestim = IsoSampGetDPrime(stro,1);
Lstimtype = (sign(uniquestim(:,1)) == sign(uniquestim(:,2))); % includes blank

% Rasters for all (Lstimtype) conditions
offset = [-.1 .1];  % pre and post time wrt stimon/stimoff
figure; axes; hold on;
set(gca,'TickDir','out');
counter = 0;
for j = find(Lstimtype)'
    L = Lcc == uniquestim(j,1) & Mcc == uniquestim(j,2) & TF == uniquestim(j,3);
    for i = find(L)'
        tmpspikes = spikes{i}-stimon_t(i);
        tmpspikes(tmpspikes < offset(1) | tmpspikes > dur+offset(2)) = [];
        nspikestot = length(tmpspikes);
        plot([tmpspikes tmpspikes]',[zeros(nspikestot,1) 1*ones(nspikestot,1)]'+counter,'k-','linewidth',1);
        counter = counter + 1;
    end
    plot([offset(1) dur+offset(2)],counter*[1 1],'b:');
    h = text(-.12,counter-sum(L)/2,num2str(uniquestim(j,3),2),'FontSize',12,'HorizontalAlignment','right');
end
set(gca,'Xlim',[0+offset(1) dur+offset(2)],'Ylim',[0 counter],'Ytick',[],'Xtick',[],'Box','on');
set(gcf,'Renderer','painters')
%% Section 12 locations tested
%linspace the sizes with linspace(6, 16, 5);
sizes = linspace(6, 16, 5);
figure; hold on;
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
location_query = 'SELECT DISTINCT rfX, rfY, subjID FROM LMTF WHERE recDate > ''2016-05-27'' AND recDate < ''2017-02-03'' AND quality = 1 ORDER BY `LMTF`.`subjID` ASC';
Lsubj_list = fetch(conn, location_query);
close(conn);
plotted = [];
for i = 1:length(Lsubj_list)
    x = Lsubj_list{i,1};
    y = Lsubj_list{i,2};
    if ~isempty(plotted)
        [one,~] = ismember(plotted,[x, y ,1],'rows');
        [two,~] = ismember(plotted,[x, y ,2],'rows');
        [three,~] = ismember(plotted,[x, y ,3],'rows');
        [four,~] = ismember(plotted,[x, y ,4],'rows');
        if sum(one)
            plotted(one,:) = [x, y, 2];
            msize = sizes(2);
        elseif sum(two)
            plotted(two, :) = [x, y, 3];
            msize = sizes(3);
        elseif sum(three)
            plotted(three, :) = [x, y, 4];
            msize = sizes(4);
        elseif sum(four)
            plotted(four, :) = [x, y, 5];
            msize = sizes(5);
        else
            plotted(end+1,:) = [x, y, 1];
            msize = sizes(1);
        end
    else
        plotted(1,:) = [x, y, 1];
        msize = sizes(1);
    end
    plot(x/10, y/10, 'o', 'MarkerFaceColor', 'k', 'markerEdgeColor', [0 0 0], 'MarkerSize', msize);
    plot([0 10], [0 0], '--b');
    pbaspect([1 2 1]);
    set(gca, 'TickDir', 'out');
    set(gca, 'ytick', [-10 -5 0 5 10]);
    set(gca, 'xtick', [0 5 10]);
end
%% section 13 table 1 in paper
% conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
% subjs = {'A', 'U', 'E', 'G'};
% num_dp = {};
% 
% for i = 1:length(subjs)
%    subjID = subjs{i};
%    flist_query = sprintf('SELECT fileID FROM LMTF WHERE subjID = ''%s'' AND quality = 1 AND neuron IS NULL AND monID = ''ProPixx'' AND recDate > ''2016-05-27'' AND recDate < ''2017-06-1''', subjID);
%    flist = fetch(conn, flist_query);
%    numf = length(flist);
%    loc_num_dp = 0;
%    for j = 1:numf
%       currfile = flist{j};
%       locfile = findfile(currfile);
%       stro = nex2stro(locfile);
%       num_stims = stro.sum.exptParams.n_stim;
%       loc_num_dp = loc_num_dp + num_stims;
%    end
%    num_dp{i} = {subjID, loc_num_dp};
%    disp(num_dp{i});
% end

% Just use the LMTF.mat structure!
HIGHLOWTHRESH = 5;

load /Users/greghorwitz/Desktop/MatlabCode/Zack/IsoSamp/private/data/LMTF.mat
subjs = {'A', 'U', 'E', 'G'};
for sub_idx = 1:length(subjs)
    subjID = subjs{sub_idx};
    raws = getfield(eval(subjID),'raw');
    counter = [0 0 0 0];
    for i = 1:length(raws)
        raw = raws{i};
        Lopp = raw(:,1) < 0 | raw(:,2) < 0;
        Lnopp = ~Lopp;
        Lhigh = raw(:,3) > HIGHLOWTHRESH;
        Llow = ~Lhigh;
        %Loog = raw(:,4);
        counter = counter + [sum(Lopp&Llow) sum(Lnopp&Llow) sum(Lopp&Lhigh) sum(Lnopp&Lhigh)]; 
    end
    {subjID sum(counter)}
    counter
end

%% Section 14 rf by date chart
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
subjs = {'A', 'U', 'E', 'G'};

for i = 1:length(subjs)
   subjID = subjs{i};
   query = sprintf('SELECT recDate, rfX, rfY FROM `LMTF` WHERE recDate > ''2016-05-27'' AND recDate < ''2017-02-03'' AND quality = 1 AND subjID = ''%s'' GROUP BY recDate', subjID);
   date_rf_list = fetch(conn, query);
   keyboard;
   dates = datestr(cell2mat(date_rf_list(:,1)));
   date_rf_list = {dates, date_rf_list(:,2), date_rf_list(:,3)};
   xlswrite('daterflist', date_rf_list, subjID);
end
%% Section 15 regression fit to all unconstrained models plus plot
subjs = {'A', 'U', 'E', 'G'};
for i = 1:length(subjs)
   subjID = subjs{i};
   lmtfModelClouds(subjID, 'firstroundmodels', 5);
end