% This code is intended to organize GLMP .nex datafiles into trial, par, 
% and plat structures.

% NOTE: In early datasets, imprecise numbers are passed between matlab and
% rex, resulting in stimuli that lie close to a line, but not exactly on
% one.  I have included a correctiomn for normalized datasets, but no 
% correction for the orignal space.

% TO DO: Provide angle correction for poltheta_orig

% 9/12      Created.    JPW
% 11/16/12  Modified. (Abandoned attempt to correct angles in original
%                datasets space.)  JPW
% 5/15/14   Modified.  Added GLMP structure to help with population
%               analyses.

function [trial par plat] = OrganizeRawGLMPData(rawdata)
global GLMP

fig = 1;
filename = rawdata.sum.fileName(end-13:end-4);

%% Fill in missing fields
if isnan(rawdata.trial(1,strcmp(rawdata.sum.trialFields(1,:),'Scc')))
    rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Scc')) = 0;
end
if isnan(rawdata.trial(1,strcmp(rawdata.sum.trialFields(1,:),'Theta')))
    rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Theta')) = rawdata.sum.exptParams.grating_theta;
end
if all(unique(rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Theta'))) == [-1 1]')
    DM = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Theta'));
    DM(DM==1) = 0;
    DM(DM==-1) = pi;
    rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Theta')) = rawdata.sum.exptParams.grating_theta;
    rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Theta')) = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Theta'))+DM;
    L = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Theta')) > 2*pi;
    rawdata.trial(L,strcmp(rawdata.sum.trialFields(1,:),'Theta')) = rawdata.trial(L,strcmp(rawdata.sum.trialFields(1,:),'Theta'))-2*pi;
end
if isnan(rawdata.trial(1,strcmp(rawdata.sum.trialFields(1,:),'Sigma X')))
    rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Sigma X')) = rawdata.sum.exptParams.sigma_X;
end
if isnan(rawdata.trial(1,strcmp(rawdata.sum.trialFields(1,:),'Sigma Y')))
    rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Sigma Y')) = rawdata.sum.exptParams.sigma_Y;
end
if isnan(rawdata.trial(1,strcmp(rawdata.sum.trialFields(1,:),'Drift Rate')))
    rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Drift Rate')) = rawdata.sum.exptParams.driftrate;
end
if isnan(rawdata.trial(1,strcmp(rawdata.sum.trialFields(1,:),'nstd')))
    rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'nstd')) = rawdata.sum.exptParams.nstd;
end
if isnan(rawdata.trial(1,strcmp(rawdata.sum.trialFields(1,:),'nexp')))
    rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'nexp')) = rawdata.sum.exptParams.nexpanse;
end
if sum(strcmp(rawdata.sum.trialFields(1,:),'Adaptive Paradigm'))==0
    rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'adaptive')) = 0;
end



%% Initiate Trial Structure

trial.datafile = filename;
trial.fpon = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'fpon_t'));
trial.fpacq = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'fpacq_t'));
trial.stimon = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'stimon_t'));
trial.stimoff = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'fpoff_t'));
trial.stimDur = trial.stimoff - trial.fpon;
trial.Lcc_orig = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Lcc'));
trial.Mcc_orig = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Mcc'));
trial.Scc_orig = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Scc'));
trial.stim_orig = [trial.Lcc_orig trial.Mcc_orig];
trial.theta = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Theta'));
trial.sigX = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Sigma X'));
trial.sigY = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Sigma Y'));
trial.driftRate = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Drift Rate'));
trial.nstd = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'nstd'));
trial.nexp = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'nexp'));
trial.tspikes = rawdata.ras(:,1);
[trial.poltheta_orig trial.polrho_orig] = cart2pol(trial.Lcc_orig,trial.Mcc_orig);
if sum(strcmp(rawdata.sum.trialFields(1,:),'Adaptive Paradigm'))==0
    trial.adaptive = zeros(size(trial.fpon));
else
    trial.adaptive = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Adaptive Paradigm'));
end


% Transform stimuli onto normalized grid
% Find max contrast along L+M and L-M directions
LpMrho = max(trial.polrho_orig(trial.poltheta_orig == pi/4,:));
LmMrho = max(trial.polrho_orig(trial.poltheta_orig == -pi/4,:));
% Ideal angles and contrasts
thetas = [pi/4 -pi/4 -3*pi/4 3*pi/4]';
rhos = [LpMrho LmMrho LpMrho LmMrho]';
[x,y] = pol2cart(thetas,rhos);
[normX normY] = pol2cart(thetas,1);
normMat = [[x y]\[normX normY]]';
trial.stim_norm = [normMat*[trial.Lcc_orig trial.Mcc_orig]']';
trial.Lcc_norm = trial.stim_norm(:,1);
trial.Mcc_norm = trial.stim_norm(:,2);

% Correcting thetas (roundoff error)
[trial.poltheta_norm trial.polrho_norm] = cart2pol(trial.stim_norm(:,1),trial.stim_norm(:,2));
nrads = 64;
divs = -pi-pi/(nrads/2):pi/(nrads/2):pi+pi/(nrads/2);
rads = divs(2:2:end);
edges = divs(1:2:end);
[n,binIdx] = histc(trial.poltheta_norm,edges);
[n1,binIdx1] = histc(trial.poltheta_orig,edges);
trial.poltheta_norm = rads(binIdx)';
trial.poltheta_orig = rads(binIdx1)';
if any(trial.poltheta_norm==-pi)
    trial.poltheta_norm(trial.poltheta_norm==-pi)=pi;
    trial.poltheta_orig(trial.poltheta_orig==-pi)=pi;
end        


% % Trying to bin orig pts (theta correction)
% % Get matrix that goes from norm to orig
% idealMat = [[normX normY]\[x y]]';
% 
% % Make ring of cartesian points from polar points in original space
% [edgesX edgesY] = pol2cart(edges,max(trial.polrho_orig));
% [radsX radsY] = pol2cart(rads,max(trial.polrho_orig));
% 
% % Turn ring of normalized cartesian points into ring of original cartesian points
% idealEdgesPts = [idealMat*[edgesX;edgesY]]';
% idealRadsPts = [idealMat*[radsX;radsY]]';
% 
% % Turn ring of orignal cartesian points to ring of polar orginal points
% [idealEdges junk] = cart2pol(idealEdgesPts(:,1),idealEdgesPts(:,2));
% [idealRads junk] = cart2pol(idealRadsPts(:,1),idealRadsPts(:,2));
% 
% % Minor correction for L+M, L-M directions (this sucks)
% idealRads(softEq(idealRads,pi/4)) = pi/4;
% idealRads(softEq(idealRads,3*pi/4)) = 3*pi/4;
% idealRads(softEq(idealRads,-pi/4)) = -pi/4;
% idealRads(softEq(idealRads,-3*pi/4)) = -3*pi/4;
% 
% % Sort transformed points by angle (low to high)
% [idealEdges,Edgesidx] = sort(idealEdges);
% bottom = idealEdges(end) - 2*pi;
% top = idealEdges(1) + 2*pi;
% idealEdges(1) = bottom;
% idealEdges(end) = top;
% [idealRads,Radsidx] = sort(idealRads);
% 
% % histc according to new ideal angles
% [n,binIdx] = histc(trial.poltheta_orig,idealEdges);
% 
% % idealRads = sort(idealRads)
% trial.poltheta_orig = idealRads(binIdx);

% %Plot histogram of angles
% figure(fig); fig=fig+1;
% subplot(1,2,1);hold on; grid on;
% set(gcf,'Name',filename,'NumberTitle','off')
% xlim([-pi-.1 pi+.1]);
% hist(trial.poltheta_norm,edges)
% title('Normalized')
% 
% subplot(1,2,2); hold on; grid on;
% set(gcf,'Name',filename,'NumberTitle','off')
% hist(trial.poltheta_orig,idealEdges)
% xlim([-pi-.1 pi+.1]);
% title('Original')

%% Pull out relevant spikes

% Set up variables
ntrials = size(trial.stim_orig,1);
trial.nspikes = nan(size(rawdata.trial,1),1);
trial.normtspikes = cell(size(rawdata.trial,1),1);
trial.fr = nan(ntrials,1);
trial.baselinensp = nan(ntrials,1);
bins = linspace(0,max(trial.stimDur),50);
PSTH = zeros(1,length(bins));
dt = zeros(ntrials,1);
trial.baselinefr = zeros(ntrials,1);
circDelay = .05;

% Generate spiking profile across all trials
for i = 1:ntrials
    trial.normtspikes{i} = rawdata.ras{i,1} - trial.stimon(i);
    PSTH = PSTH + histc(trial.normtspikes{i}',bins);
end

% Find baseline spikerate
for i = 1:ntrials
    nsp =  sum(rawdata.ras{i,1} > trial.fpacq(i) + circDelay & rawdata.ras{i,1} < trial.stimon(i));
    dt(i) = trial.stimon(i) - trial.fpacq(i);
    trial.baselinensp(i) = nsp;
    trial.baselinefr(i) = nsp./(dt(i));
    %trial.baselinefr(i) = nsp./(dt(i) - circDelay);
end

%Gaussian profile
gauss = @(x,params)(params(2)+params(1)*exp((-(x-params(3)).^2)/params(4).^2));
bguess = mean(PSTH(bins<.25));
aguess = max(PSTH)-bguess;
nsigma = 1;
muguess = bins(find(PSTH == max(PSTH),1,'first'));
sigmaguess = .1; % terrible. Need something better.
fittedparams = fminsearch(@(params)sum((gauss(bins,params)-PSTH).^2),[aguess,bguess,muguess,sigmaguess]);
offset = [fittedparams(3)-nsigma*fittedparams(4) fittedparams(3)+nsigma*fittedparams(4)];

% Now counting up spikes in a window
for i = 1:ntrials
    trial.nspikes(i) = sum(rawdata.ras{i,1} > trial.stimon(i)+offset(1) & rawdata.ras{i,1}<trial.stimon(i)+offset(2));
    dt = offset(2)-offset(1);
    trial.fr(i) = trial.nspikes(i)./dt;
end

trial.GaussianSpikingProfile.a = fittedparams(1);
trial.GaussianSpikingProfile.b = fittedparams(2);
trial.GaussianSpikingProfile.mu = fittedparams(3);
trial.GaussianSpikingProfile.sigma = fittedparams(4);


% Plot spiking profile
figure(fig); fig=fig+1; clf; hold on; grid on;
set(gcf,'Name',filename,'NumberTitle','off')
xlim([0 max(trial.stimDur)])
bar(bins,PSTH);
plot(bins,gauss(bins,fittedparams),'b-.','LineWidth',3)
plot(offset,[fittedparams(1) fittedparams(1)],'k-','linewidth',3);
title('Histogram of Spikes')
xlabel('Time (ms) from Stimulus Onset')
ylabel('# of Spikes')


%% Initiate Par Structure

% By Paradigm Properties (aka by unique stimulus)
[uniqueCond,m,idx] = unique([trial.Lcc_orig trial.Mcc_orig trial.Scc_orig trial.theta trial.sigX trial.sigY trial.driftRate trial.nstd trial.nexp trial.adaptive],'rows');

% Collapse repeated stimuli
par.datafile = filename;
par.Lcc_orig = uniqueCond(:,1);
par.Mcc_orig = uniqueCond(:,2);
par.Scc_orig = uniqueCond(:,3);
par.stim_orig = [par.Lcc_orig par.Mcc_orig];
par.stim_norm = [normMat * par.stim_orig']';
par.Lcc_norm = par.stim_norm(:,1);
par.Mcc_norm = par.stim_norm(:,2);

% Define unique conditions
par.theta = uniqueCond(:,4);
par.sigX = uniqueCond(:,5);
par.sigY = uniqueCond(:,6);
par.driftRate = uniqueCond(:,7);
par.nstd = uniqueCond(:,8);
par.nexp = uniqueCond(:,9);
par.adaptive = uniqueCond(:,10);

% Preallocate Space
par.poltheta_orig = nan(size(uniqueCond,1),1);
par.polrho_orig = nan(size(uniqueCond,1),1);
par.poltheta_norm = nan(size(uniqueCond,1),1);
par.polrho_norm = nan(size(uniqueCond,1),1);
par.frs = cell(size(uniqueCond,1),1);
par.meanfr = nan(size(uniqueCond,1),1);
par.varfrs = nan(size(uniqueCond,1),1);
par.nsamps = nan(size(uniqueCond,1),1);
par.nspikes = cell(size(uniqueCond,1),1);
par.meannspikes = nan(size(uniqueCond,1),1);
par.varnspikes = nan(size(uniqueCond,1),1);
par.baselinefrs = cell(size(uniqueCond,1),1);
par.meanblfr = nan(size(uniqueCond,1),1);
par.normtspikes = cell(size(uniqueCond,1),1);
par.catntspikes = cell(size(uniqueCond,1),1);

for n = 1:size(par.stim_orig,1)
    L = idx == n;
    try
        par.poltheta_orig(n) = unique(trial.poltheta_orig(L));
        par.polrho_orig(n) = unique(trial.polrho_orig(L));
        par.poltheta_norm(n) = unique(trial.poltheta_norm(L));
        par.polrho_norm(n) = unique(trial.polrho_norm(L));
    catch
        keyboard
    end
    par.frs{n} = trial.fr(L)';
    par.meanfr(n) = mean(par.frs{n});
    par.varfrs(n) = var(par.frs{n});
    par.nspikes{n} = trial.nspikes(L)';
    par.nsamps(n) = sum(L);
    par.meannspikes(n) = mean(trial.nspikes(L));
    par.varnspikes(n) = var(trial.nspikes(L));
    par.baselinefrs{n} = trial.baselinefr(L);
    par.meanblfr(n) = mean(par.baselinefrs{n});
    par.normtspikes{n} = trial.normtspikes(L);
    par.catntspikes{n} = cat(1,trial.normtspikes{L});
end


%% Initiate Platform Organization

% By paradigm
[uniqueCond,m,idx] = unique([par.Scc_orig par.theta par.sigX par.sigY par.driftRate par.nstd par.nexp par.adaptive],'rows');
traits = cellstr(['Scc      ';'theta    ';'sigX     ';'sigY     ';...
    'driftRate';'nstd     ';'nexp     ';'adaptive ']);
if size(uniqueCond,1) > 1
    platDiffIdx = find(diff(uniqueCond));
else
    platDiffIdx = [];
end
nplat = size(uniqueCond,1);
plat = cell(nplat,1);

for b = 1:nplat
    L = idx == b;
    plat{b}.datafile = filename;
    traitVals = num2cell(uniqueCond(b,:))';
    if ~isempty(platDiffIdx)
        plat{b}.platDiff = cat(2,traits(platDiffIdx));
    end
    plat{b}.allProp = cat(2,traits,traitVals);
    plat{b}.transMat = normMat;
    plat{b}.par.Lcc_orig = par.Lcc_orig(L,:);
    plat{b}.par.Mcc_orig = par.Mcc_orig(L,:);
    plat{b}.par.Scc_orig = par.Scc_orig(L,:);
    plat{b}.par.stim_orig = par.stim_orig(L,:);
    plat{b}.par.poltheta_orig = par.poltheta_orig(L);
    plat{b}.par.polrho_orig = par.polrho_orig(L);
    plat{b}.par.Lcc_norm = par.Lcc_norm(L);
    plat{b}.par.Mcc_norm = par.Mcc_norm(L);
    plat{b}.par.stim_norm = par.stim_norm(L,:);
    plat{b}.par.poltheta_norm = par.poltheta_norm(L);
    plat{b}.par.polrho_norm = par.polrho_norm(L);
    plat{b}.par.theta = par.theta(L,:);
    plat{b}.par.sigX = par.sigX(L,:);
    plat{b}.par.sigY = par.sigY(L,:);
    plat{b}.par.driftRate = par.driftRate(L,:);
    plat{b}.par.nstd = par.nstd(L,:);
    plat{b}.par.nexp = par.nexp(L,:);
    plat{b}.par.frs = {par.frs{L}}';
    plat{b}.par.meanfr = par.meanfr(L);
    plat{b}.par.nspikes = par.nspikes(L,:);
    plat{b}.par.meannspikes = par.meannspikes(L);
    plat{b}.par.varnspikes = par.varnspikes(L);
    plat{b}.par.varfrs = par.varfrs(L);
    plat{b}.par.nsamps = par.nsamps(L);
    plat{b}.par.baselinefrs = par.baselinefrs(L);
    plat{b}.par.meanblfr = par.meanblfr(L);
    plat{b}.par.normtspikes = par.normtspikes(L);
    plat{b}.par.catntspikes = par.catntspikes(L);
    plat{b}.par.adaptive = par.adaptive(L);
end


% By trial
[uniqueCond,m,idx] = unique([trial.Scc_orig trial.theta trial.sigX trial.sigY trial.driftRate trial.nstd trial.nexp trial.adaptive],'rows');

for b = 1:nplat
    L = idx == b;
    plat{b}.trial.Lcc_orig = trial.Lcc_orig(L);
    plat{b}.trial.Mcc_orig = trial.Mcc_orig(L);
    plat{b}.trial.Scc_orig = trial.Scc_orig(L);
    plat{b}.trial.stim_orig = trial.stim_orig(L,:);
    plat{b}.trial.poltheta_orig = trial.poltheta_orig(L);
    plat{b}.trial.polrho_orig = trial.polrho_orig(L);
    plat{b}.trial.Lcc_norm = trial.Lcc_norm(L);
    plat{b}.trial.Mcc_norm = trial.Mcc_norm(L);
    plat{b}.trial.stim_norm = trial.stim_norm(L,:);
    plat{b}.trial.poltheta_norm = trial.poltheta_norm(L);
    plat{b}.trial.polrho_norm = trial.polrho_norm(L);
    plat{b}.trial.fpon = trial.fpon(L);
    plat{b}.trial.fpacq = trial.fpacq(L);
    plat{b}.trial.stimon = trial.stimon(L);
    plat{b}.trial.stimoff = trial.stimoff(L);
    plat{b}.trial.stimDur = trial.stimDur(L);
    plat{b}.trial.baselinensp = trial.baselinensp(L);
    plat{b}.GaussianSpikingProfile = trial.GaussianSpikingProfile;
    plat{b}.trial.theta = trial.theta(L);
    plat{b}.trial.sigX = trial.sigX(L);
    plat{b}.trial.sigY = trial.sigY(L);
    plat{b}.trial.driftRate = trial.driftRate(L);
    plat{b}.trial.nstd = trial.nstd(L);
    plat{b}.trial.nexp = trial.nexp(L);
    plat{b}.trial.fr = trial.fr(L);
    plat{b}.trial.baselinefr = trial.baselinefr(L);
    plat{b}.trial.normtspikes = trial.normtspikes(L);
    plat{b}.trial.nspikes = trial.nspikes(L);
    plat{b}.trial.tspikes = trial.tspikes(L);
    plat{b}.trial.adaptive = trial.adaptive(L);
end

for b = 1:nplat
    L = idx == b;
    GLMP.datafile = rawdata.sum.fileName(end-13:end-4);
    GLMP.subunit{b}.GLMPIdx = find(L);
    GLMP.subunit{b}.fpon = trial.fpon(L);
    GLMP.subunit{b}.fpacq = trial.fpacq(L);
    GLMP.subunit{b}.stimon = trial.stimon(L);
    GLMP.subunit{b}.stimoff = trial.stimoff(L);
    GLMP.subunit{b}.stimDur = trial.stimDur(L);
    GLMP.subunit{b}.Lcc = trial.Lcc_orig(L);
    GLMP.subunit{b}.Mcc = trial.Mcc_orig(L);
    GLMP.subunit{b}.Scc = trial.Scc_orig(L);
    GLMP.subunit{b}.theta = trial.poltheta_orig(L);
    GLMP.subunit{b}.rho = trial.polrho_orig(L);
    GLMP.subunit{b}.spiketimes = trial.tspikes(L);
    GLMP.subunit{b}.normspiketimes = trial.normtspikes(L);
    GLMP.subunit{b}.nspikes = trial.nspikes(L);
    GLMP.subunit{b}.fr = trial.fr(L);
    GLMP.subunit{b}.blnspikes = trial.baselinensp(L);
    GLMP.subunit{b}.blfr = trial.baselinefr(L);
    
    % Collapsing across like stimuli and computing mean spike rates
    x = GLMP.subunit{b}.Lcc;
    y = GLMP.subunit{b}.Mcc;
    lmcoords = unique([x y],'rows');
    GLMP.subunit{b}.meannspikes = nan(size(lmcoords,1),1);
    GLMP.subunit{b}.meanfr = nan(size(lmcoords,1),1);
    GLMP.subunit{b}.uniqueLcc = lmcoords(:,1);
    GLMP.subunit{b}.uniqueMcc = lmcoords(:,2);
    for i = 1:size(lmcoords,1)
        L = (x == lmcoords(i,1)) & (y == lmcoords(i,2));
        n(i) = sum(L);
        trlidxs = find(L);
        %nspikes = [];
        for j = 1:numel(trlidxs)
            tr = trlidxs(j);
            spiketimes = GLMP.subunit{b}.normspiketimes{tr};
            GLMP.subunit{b}.spiketimes_col{i,j} = spiketimes;
        end
        GLMP.subunit{b}.spiketimes_cat{i} = cat(1,GLMP.subunit{b}.spiketimes_col{i,:});
        GLMP.subunit{b}.meannspikes(i) = mean(GLMP.subunit{b}.nspikes(L));
        GLMP.subunit{b}.meanfr(i) = mean(GLMP.subunit{b}.fr(L));
    end

    [GLMP.subunit{b}.uniquetheta GLMP.subunit{b}.uniquerho] = cart2pol(GLMP.subunit{b}.uniqueLcc,GLMP.subunit{b}.uniqueMcc);

    
end



%% Plot!

for p = 1:numel(plat)
    
    x = plat{p}.par.Lcc_orig;
    y = plat{p}.par.Mcc_orig;
    z = plat{p}.par.meannspikes;
    F = TriScatteredInterp(x,y,z);
    [qx,qy] = meshgrid(min(x):.01:max(x),min(y):.01:max(y));
    qz = F(qx,qy);
    
    % Plot original and transformed stimuli
    figure(fig); clf; fig=fig+1;
    plotTitle = [filename ' Plat # ' num2str(p)];
    set(gcf,'Name',plotTitle,'NumberTitle','off','Position',[138 402 805 349])
    subplot(1,3,1); hold on; grid on;
    plot3(plat{p}.trial.Lcc_orig,plat{p}.trial.Mcc_orig,plat{p}.trial.nspikes,'o')
    plot3(plat{p}.trial.Lcc_orig(1),plat{p}.trial.Mcc_orig(1),plat{p}.trial.nspikes(1),'r*')
    mesh(qx,qy,qz)
    xlabel('Lcc')
    ylabel('Mcc')
    title('Original Dataset')
    axis equal square
    xlim([min(x) max(x)]);
    ylim([min(y) max(y)]);
    
    x = plat{p}.par.Lcc_norm;
    y = plat{p}.par.Mcc_norm;
    z = plat{p}.par.meannspikes;
    F = TriScatteredInterp(x,y,z);
    [qx,qy] = meshgrid(min(x):.1:max(x),min(y):.1:max(y));
    qz = F(qx,qy);
    subplot(1,3,2); hold on; grid on; axis equal
    plot3(plat{p}.trial.stim_norm(:,1),plat{p}.trial.stim_norm(:,2),plat{p}.trial.nspikes,'ko')
    plot3(plat{p}.trial.stim_norm(1,1),plat{p}.trial.stim_norm(1,2),plat{p}.trial.nspikes(1),'r*')
    mesh(qx,qy,qz)
    axis equal square
    xlim([min(x) max(x)]);
    ylim([min(y) max(y)]);
    xlabel('Lcc (Normalized)')
    ylabel('Mcc (Normalized)')
    title('Normalized Dataset')
    
    % Display Stimulus Properties
    subplot(1,3,3); cla
    set(gca,'Visible','off')
    if ~isempty(platDiffIdx)
        propIdx = find(strcmp(plat{p}.platDiff(1),plat{1}.allProp(:,1)));
    end
    propVert = linspace(.8,.2,numel(plat{p}.allProp(:,1)));
    for i = 1:numel(plat{p}.allProp(:,1))
        str = [num2str(cell2mat(plat{p}.allProp(i,1))) ' = ' num2str(cell2mat(plat{p}.allProp(i,2)))];
        if ~isempty(platDiffIdx) && any(i == propIdx)
            text(0,propVert(i),str,'FontName','Courier','FontSize',15,'EdgeColor','red')
        else
            text(0,propVert(i),str,'FontName','Courier','FontSize',15)
        end
    end

    
    
end
