%% This code is intended to calculate population likelihood fits across the data

% Created   4/11/12     JPW
% Making a variant for cell FR variability...   7/26/12     JPW
% Changing code to look at likelihood fits
%    across the population (and saving as new script)  8/12     JPW


clear all
close all

library = 'Users/jpatrickweller/Documents/MATLAB/GLMPDatafiles/';
n=1;

% Bring in raw data
%datafile{n} = [char(library) 'S110211007.nex']; n = n+1;% 1 Symmetric Trench (Most Data)
%datafile{n} = [char(library) 'S110411012_2.nex']; n = n+1;% 2
%datafile{n} = [char(library) 'S110911003.nex']; n = n+1;% 3 Symmetric Trench
%datafile{n} = [char(library) 'S110911007.nex']; n = n+1;% 4 Symmetric Trench
%datafile{n} = [char(library) 'S111011003.nex']; n = n+1;% 5 Asymmetric Trench (Second to Most Data)
%datafile{n} = [char(library) 'S112311006.nex']; n = n+1;% 6 Asymmetric Trench
%datafile{n} = [char(library) 'S113011002.nex']; n = n+1;% 7 Complex Cell
%datafile{n} = [char(library) 'S120511008.nex']; n = n+1;% 8 Symmetric Trench

% Radial Grid
% (Mostly) Luminance
datafile{n} = [char(library) 'S121211003.nex']; n = n+1;% 9 Symmetric Trench
datafile{n} = [char(library) 'S121211006.nex']; n = n+1;% 10 Symmetric L-M Trench
datafile{n} = [char(library) 'S121511004.nex']; n = n+1;% 11 Symmetric L-M Trench (Small Dataset)

% (Mostly) Chromatic
datafile{n} = [char(library) 'S041912004.nex']; n = n+1;% 12 Chromatic Cell? (Bug: Theta Unknown)

% Ellipsoidal
datafile{n} = [char(library) 'S120811003.nex']; n = n+1;% 13 Ellipsoidal Cell? Extremely noisy.
datafile{n} = [char(library) 'S120911003.nex']; n = n+1;% 14 Ellipsoidal Cell? Extremely noisy.
%datafile{n} = [char(library) 'S042012002.nex']; n = n+1;% 15 Ellipsoidal Cell! (Bug: Theta Unknown)
datafile{n} = [char(library) 'A060112003.nex']; n = n+1;% 16 Ellispoidal Cell! (Apollo's First Cell!)
datafile{n} = [char(library) 'A062612003.nex']; n = n+1;% 17 Ellispoidal Cell!
%datafile{n} = [char(library) 'A062612004.nex']; n = n+1;% 18 Ellipsoidal (Continuation of 062612003. Not getting all waveforms?)
datafile{n} = [char(library) 'A070312003.nex']; n = n+1;% 19 Ellipsoidal Cell (Tiny Dataset)
datafile{n} = [char(library) 'A070312005.nex']; n = n+1;% 20 Ellipsoidal Cell
datafile{n} = [char(library) 'A071012004.nex']; n = n+1;% 21 Ellipsoidal Cell (10 repeats)
datafile{n} = [char(library) 'A071212005.nex']; n = n+1;% 22 Ellipsoidal Cell
datafile{n} = [char(library) 'A071612002.nex']; n = n+1;% 23 Ellipsoidal Cell
datafile{n} = [char(library) 'A071612005.nex']; n = n+1;% 24 Ellipsoidal Cell (Same cell as A071612002)

% Plateaus: S-Cone
datafile{n} = [char(library) 'S011312002.nex']; n = n+1;% 25 Ellipse!
datafile{n} = [char(library) 'S020112002.nex']; n = n+1;% 26 Symmetric L-M Trench (check iso)
datafile{n} = [char(library) 'S020212007.nex']; n = n+1;% 27 Asymmetric L-M Trench

% Plateaus: 2 Directions of Motion
datafile{n} = [char(library) 'S021012007.nex']; n = n+1;% 28 Symmetric L-M Trench
datafile{n} = [char(library) 'S032312002.nex']; n = n+1;% 29 Luminance Cell
datafile{n} = [char(library) 'S032912004.nex']; n = n+1;% 30 Symmetric L-M Trench
datafile{n} = [char(library) 'S040412002.nex']; n = n+1;% 31 Chromatic Cell!!
datafile{n} = [char(library) 'S041012002.nex']; n = n+1;% 32 Ellipsoidal Cell!!
%datafile{n} = [char(library) 'S040612002.nex']; n = n+1;% 33 Great iso, but not much structure...


% Bugs...
%rawdata = nex2stro([char(library) 'S032812004.nex']);% Dropped 2 headers - must trim orignal file
%rawdata = nex2stro([char(library) 'S032812006.nex']);% Grating File!

% Setting up some variables
fig = 1;
slope = nan(n-1,1);
sig2gauss = 8;

ndatafiles = size(datafile,2);

% Begin grand loop...
for loop = 1:ndatafiles
    rawdata = nex2stro(datafile{loop});
    disp(['Now processing ' datafile{loop}])
    
    %% Organize Raw Data
    
    % Fill in missing fields
    if isnan(rawdata.trial(1,strcmp(rawdata.sum.trialFields(1,:),'Scc')))
        rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Scc')) = 0;
    end
    if isnan(rawdata.trial(1,strcmp(rawdata.sum.trialFields(1,:),'Theta')))
        rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Theta')) = rawdata.sum.exptParams.theta;
    end
    if all(unique(rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Theta'))) == [-1 1]')
        DM = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Theta'));
        DM(DM==1) = 0;
        DM(DM==-1) = pi;
        rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Theta')) = rawdata.sum.exptParams.theta;
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
    
    
    %% Organize Data by Trial
    
    trial.fpon = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'fpon_t'));
    trial.fpacq = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'fpacq_t'));
    trial.stimon = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'stimon_t'));
    trial.stimoff = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'fpoff_t'));
    trial.duration = trial.stimoff - trial.fpacq;
    trial.Lcc = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Lcc'));
    trial.Mcc = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Mcc'));
    trial.Scc = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Scc'));
    [trial.opoltheta trial.opolrho] = cart2pol(trial.Lcc,trial.Mcc);
    trial.stim = [trial.Lcc trial.Mcc];
    trial.theta = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Theta'));
    trial.sigX = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Sigma X'));
    trial.sigY = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Sigma Y'));
    trial.driftRate_cps = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Drift Rate'));
    trial.driftRate_dvaps = trial.driftRate_cps .* trial.sigY .* sig2gauss;
    trial.nstd = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'nstd'));
    trial.nexp = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'nexp'));
    trial.tspikes = rawdata.ras(:,1);
    
    % Transform stimuli onto normalized grid
    LpMrho = max(trial.opolrho(trial.opoltheta == pi/4,:));
    LmMrho = max(trial.opolrho(trial.opoltheta == -pi/4,:));
    thetas = [pi/4 -pi/4 -3*pi/4 3*pi/4]';
    rhos = [LpMrho LmMrho LpMrho LmMrho]';
    [x,y] = pol2cart(thetas,rhos);
    normMat = [0 1; 1 0; 0 -1; -1 0];
    transMat = [[x y]\normMat]';
    trial.tstim = rndofferr([transMat * [trial.Lcc trial.Mcc]']',3);
    trial.LmM = trial.tstim(:,1);
    trial.LpM = trial.tstim(:,2);
    
    % Correcting thetas (roundoff error)
    [trial.tpoltheta trial.tpolrho] = cart2pol(trial.tstim(:,1),trial.tstim(:,2));
    divs = -pi-pi/64:pi/64:pi+pi/64;
    rads = divs(2:2:end);
    edges = divs(1:2:end);
    [n,binIdx] = histc(trial.tpoltheta,edges);
    trial.tpoltheta = rads(binIdx)';
    if any(trial.tpoltheta==-pi)
        trial.tpoltheta(trial.tpoltheta==-pi)=pi;
    end
    
    % Plot original and transformed stimuli
    figure(fig); clf; fig=fig+1;
    subplot(1,2,1); hold on; grid on;
    plot(trial.Lcc,trial.Mcc,'o')
    plot(trial.Lcc(1),trial.Mcc(1),'r*')
    xlabel('Lcc')
    ylabel('Mcc')
    title('Original Dataset')
    axis equal square
    subplot(1,2,2); hold on; grid on; axis equal
    plot(trial.tstim(:,1),trial.tstim(:,2),'ko')
    plot(trial.tstim(1,1),trial.tstim(1,2),'r*')
    axis equal square
    xlabel('L-M')
    ylabel('L+M')
    title('Transformed Dataset')
    
    
    %% Pull out relevant spikes
    
    % Set up variables
    ntrials = size(trial.stim,1);
    trial.nspikes = nan(size(rawdata.trial,1),1);
    trial.normtspikes = cell(size(rawdata.trial,1),1);
    trial.fr = nan(ntrials,1);
    bins = linspace(0,max(trial.duration),50);
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
        trial.baselinefr(i) = nsp./(dt(i) - circDelay);
    end
    
    % Plot spiking profile
    
    gauss = @(x,params)(params(2)+params(1)*exp((-(x-params(3)).^2)/params(4).^2));
    bguess = mean(PSTH(bins<.25));
    aguess = max(PSTH)-bguess;
    muguess = bins(find(PSTH == max(PSTH),1,'first'));
    sigmaguess = .1; % terrible. Need something better.
    fittedparams = fminsearch(@(params)sum((gauss(bins,params)-PSTH).^2),[aguess,bguess,muguess,sigmaguess]);
    offset = [fittedparams(3)-fittedparams(4) fittedparams(3)+fittedparams(4)];
    
    % Plot figure
    %figure(fig); axes; hold on; fig=fig+1;
    %bar(bins,PSTH);
    %plot(bins,gauss(bins,[aguess,bguess,muguess,sigmaguess]),'m.')
    %plot(bins,gauss(bins,fittedparams),'b-.','LineWidth',3)
    %plot(offset,[fittedparams(1) fittedparams(1)],'k-','linewidth',3);
    
    % Now counting up spikes in a window
    for i = 1:ntrials
        trial.nspikes(i) = sum(rawdata.ras{i,1} > trial.stimon(i)+offset(1) & rawdata.ras{i,1}<trial.stimon(i)+offset(2));
        dt = offset(2)-offset(1);
        trial.fr(i) = trial.nspikes(i)./dt;
    end
    
    
    %% Organizing Data
    
    % By Paradigm Properties (aka by unique stimulus)
    [uniqueCond,m,idx] = unique([trial.Lcc trial.Mcc trial.Scc trial.theta trial.sigX trial.sigY trial.driftRate_cps trial.nstd trial.nexp],'rows');
    
    % Collapse repeated stimuli
    par.Lcc = uniqueCond(:,1);
    par.Mcc = uniqueCond(:,2);
    par.Scc = uniqueCond(:,3);
    par.stim = [par.Lcc par.Mcc];
    [par.opoltheta par.opolrho] = cart2pol(par.Lcc,par.Mcc);
    par.tstim = rndofferr([transMat * par.stim(:,1:2)']',3);
    par.LmM = par.tstim(:,1);
    par.LpM = par.tstim(:,2);
    [par.opoltheta par.opolrho] = cart2pol(par.Lcc,par.Mcc);
    [par.tpoltheta par.tpolrho] = cart2pol(par.tstim(:,1),par.tstim(:,2));
    par.theta = uniqueCond(:,4);
    par.sigX = uniqueCond(:,5);
    par.sigY = uniqueCond(:,6);
    par.driftRate_cps = uniqueCond(:,7);
    par.driftRate_dvaps = par.driftRate_cps .* par.sigY .* sig2gauss;
    par.nstd = uniqueCond(:,8);
    par.nexp = uniqueCond(:,9);
    
    % Preallocate Space
    par.frs = cell(size(uniqueCond,1),1);
    par.meanfr = nan(size(uniqueCond,1),1);
    par.varfrs = nan(size(uniqueCond,1),1);
    par.meannspikes = nan(size(uniqueCond,1),1);
    par.varnspikes = nan(size(uniqueCond,1),1);
    par.baselinefrs = cell(size(uniqueCond,1),1);
    par.meanblfr = nan(size(uniqueCond,1),1);
    par.normtspikes = cell(size(uniqueCond,1),1);
    par.catntspikes = cell(size(uniqueCond,1),1);
    
    for n = 1:size(par.stim,1)
        L = idx == n;
        par.frs{n} = trial.fr(L)';
        par.meanfr(n) = mean(par.frs{n});
        par.varfrs(n) = var(par.frs{n});
        par.nspikes{n} = trial.nspikes(L)';
        par.meannspikes(n) = mean(par.nspikes{n});
        par.varnspikes(n) = var(par.nspikes{n});
        par.baselinefrs{n} = trial.baselinefr(L);
        par.meanblfr(n) = mean(par.baselinefrs{n});
        par.normtspikes{n} = trial.normtspikes(L);
        par.catntspikes{n} = cat(1,trial.normtspikes{L});
    end
    
    
    %% Platform Organization
    
    % By paradigm
    [uniqueCond,m,idx] = unique([par.Scc par.theta par.sigX par.sigY par.driftRate_cps par.nstd par.nexp],'rows');
    nplat = size(uniqueCond,1);
    plat = struct([]);
    
    for b = 1:nplat
        L = idx == b;
        plat(b).par.Lcc = par.Lcc(L,:);
        plat(b).par.Mcc = par.Mcc(L,:);
        plat(b).par.Scc = par.Scc(L,:);
        plat(b).par.stim = par.stim(L,:);
        plat(b).par.opoltheta = par.opoltheta(L);
        plat(b).par.opolrho = par.opolrho(L);
        plat(b).par.tstim = par.tstim(L,:);
        plat(b).par.tpoltheta = par.tpoltheta(L);
        plat(b).par.tpolrho = par.tpolrho(L);
        plat(b).par.LmM = par.LmM(L);
        plat(b).par.LpM = par.LpM(L);
        plat(b).par.theta = par.theta(L,:);
        plat(b).par.sigX = par.sigX(L,:);
        plat(b).par.sigY = par.sigY(L,:);
        plat(b).par.driftRate_cps = par.driftRate_cps(L,:);
        plat(b).par.driftRate_dvaps = par.driftRate_dvaps(L,:);
        plat(b).par.nstd = par.nstd(L,:);
        plat(b).par.nexp = par.nexp(L,:);
        plat(b).par.frs = {par.frs{L}}';
        plat(b).par.meanfr = par.meanfr(L);
        plat(b).par.varfrs = par.varfrs(L);
        plat(b).par.baselinefrs = par.baselinefrs(L);
        plat(b).par.meanblfr = par.meanblfr(L);
        plat(b).par.normtspikes = par.normtspikes(L);
        plat(b).par.catntspikes = par.catntspikes(L);
    end
    
    
    % By trial
    [uniqueCond,m,idx] = unique([trial.Scc trial.theta trial.sigX trial.sigY trial.driftRate_cps trial.nstd trial.nexp],'rows');
    
    for b = 1:nplat
        L = idx == b;
        plat(b).trial.Lcc = trial.Lcc(L);
        plat(b).trial.Mcc = trial.Mcc(L);
        plat(b).trial.Scc = trial.Scc(L);
        plat(b).trial.stim = trial.stim(L,:);
        plat(b).trial.opoltheta = trial.opoltheta(L);
        plat(b).trial.opolrho = trial.opolrho(L);
        plat(b).trial.tstim = trial.tstim(L,:);
        plat(b).trial.tpoltheta = trial.tpoltheta(L);
        plat(b).trial.tpolrho = trial.tpolrho(L);
        plat(b).trial.LmM = trial.LmM(L);
        plat(b).trial.LpM = trial.LpM(L);
        plat(b).trial.theta = trial.theta(L);
        plat(b).trial.sigX = trial.sigX(L);
        plat(b).trial.sigY = trial.sigY(L);
        plat(b).trial.driftRate_cps = trial.driftRate_cps(L);
        plat(b).trial.driftRate_dvaps = trial.driftRate_dvaps(L);
        plat(b).trial.nstd = trial.nstd(L);
        plat(b).trial.nexp = trial.nexp(L);
        plat(b).trial.fr = trial.fr(L);
        plat(b).trial.baselinefr = trial.baselinefr(L);
        plat(b).trial.normtspikes = trial.normtspikes(L);
        plat(b).trial.nspikes = trial.nspikes(L);
        plat(b).trial.tspikes = trial.tspikes(L);
    end
    
    
    %% 2-D Fit Surfaces To All Data
    
    disp('Fitting 2D surface to all of the data...')
    
    % Generating an initial guess
    for p = 1:nplat
        
        % Organize rotations
        allangs = unique(plat(p).trial.tpoltheta);
        angs = unique(plat(p).trial.tpoltheta(plat(p).trial.tpoltheta <= -pi/2+2*eps));
        nrots = numel(angs);
        
        % Generating variables
%         if loop == 1 & p == 1 
%             gofAllData = nan(nrots,nplat,ndatafiles);
%             likelihoodratio0 = nan(nplat,ndatafiles);
%         end
        
        % Rotate through angles
        for rot = 1:nrots
            
            vlb = [0    0 0.001 0.001 1 1 0];
            vub = [200 200    1     1 6 6 100];
            options = optimset('MaxFunEvals',5000,'MaxIter',5000,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
            sigmaguess = 0.1;
            maxpred = []; f = [];
            titles = {'L+M','L-M'};
            for i = 1:2
                if (i == 1)
                    L = plat(p).trial.tpoltheta == angs(rot) | softEq(plat(p).trial.tpoltheta,angs(rot)+pi);
                else
                    L = softEq(plat(p).trial.tpoltheta,angs(rot)+pi/2) | softEq(plat(p).trial.tpoltheta,angs(rot)+3*pi/2);
                end
                projs = plat(p).trial.tpolrho(L);
                topfr = max(plat(p).trial.nspikes(L));
                if (topfr == 0)
                    f(i,:) = [0 0 1 1 2 2 0];
                elseif isempty(topfr)
                    f(i,:) = [0 0 1 1 2 2 0];
                else
                    params0 = [topfr, topfr, sigmaguess, sigmaguess, 2, 2, 0];  % need to constrain B across color directions
                    f(i,:) = fmincon('MyFitNakaRushtonFun',params0,[],[],[],[],vlb,vub,[],options,projs,plat(p).trial.nspikes(L),'asymmetric');
                end
                
                % May want to use all data to generate 1-d nakarushtons, which
                % in turn generate initial guesses...
                pred = MyComputeNakaRushton(f(i,:),projs,'asymmetric');
                
            end
            a = mean(f(1,[1 2]));
            b = mean(f(2,[1 2]));
            params0 = max(plat(p).trial.nspikes);
            params0(2) = a./sum([a b]); %this doesn't work so well
            params0(2) = .5;
            params0(3) = .25;
            params0(4) = 2;
            params0(5) = (f(1,2)./f(1,4))./(f(1,1)./f(1,3)); if (isnan(params0(5))) params0(5) = 1; end;
            params0(6) = (f(2,2)./f(2,4))./(f(2,1)./f(2,3)); if (isnan(params0(6))) params0(6) = 1; end;
            params0(7) = mean(f(:,end));
            
            vlb = [10   0 0.0001 1  0  0  0];
            vub = [1000 1     1  6 10 10 20];
            options = optimset('MaxFunEvals',5000,'MaxIter',5000,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
            rotMat = [cos(angs(rot)) -sin(angs(rot)); sin(angs(rot)) cos(angs(rot))];
            [f0,fval] = fmincon('MyFitNakaRushtonFun',params0,[],[],[],[],vlb,vub,[],options,plat(p).trial.tstim*rotMat,plat(p).trial.nspikes,'surface1');
            gofAllData{rot,p,loop} = -fval;
            [x,y] = meshgrid(linspace(min(plat(p).trial.LpM),max(plat(p).trial.LpM),50), linspace(min(plat(p).trial.LmM),max(plat(p).trial.LmM),50));
            surface = MyComputeNakaRushton(f0,[x(:) y(:)],'surface1');
            
            %         % Plot Surface
            %         figure(fig); fig=fig+1; hold on; grid on;
            %         surf(x,y,reshape(surface,size(x)))
            %         axis([-1 1 -1 1])
            %         tempRotPts = plat(p).trial.tstim*rotMat;
            %         plot3(tempRotPts(:,1),tempRotPts(:,2),plat(p).trial.nspikes,'ko')
            %         xlabel('L-M')
            %         ylabel('L+M')
            %         zlabel('Number of Spikes')
            %         title(['Fit to all the data (n = ',num2str(size(plat(p).trial.LpM,1)),')']);
            
        end
        
        likelihoodratio0(p,loop) = min(cat(1,gofAllData{:,p,loop}))/max(cat(1,gofAllData{:,p,loop}));
        
        %     % Plot goodness of fit
        %     figure(fig); fig=fig+1; hold on; grid on;
        %     plot(angs+pi/2,-gofAllData,'o-')
        %     title('Fit Using All Data')
        %     xlabel('Rotation (Radians)')
        %     ylabel('Log Likelihood')
        
    end
    
    
    %% Building 2-D Surfaces with 2 Vectors
    
    disp('Using 2 orthagonal vectors to create 2D surfaces...')
    
    % Generating an initial guess
    for p = 1:nplat
        
        % Organize Rotations
        allangs = unique(plat(p).trial.tpoltheta);
        angs = unique(plat(p).trial.tpoltheta(plat(p).trial.tpoltheta <= -pi/2+2*eps));
        nrots = numel(angs);
        
        % Generating variables
%         if loop == 1 & p == 1
%             gof2Ax = nan{nrots,nplat,ndatafiles);
%             likelihoodratio1 = nan(nplat,ndatafiles);
%         end
        
        % Rotate through angles
        for rot = 1:nrots
            
            vlb = [0    0 0.001 0.001 1 1 0];
            vub = [200 200    1     1 6 6 100];
            options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
            sigmaguess = 0.1;
            maxpred = []; f = [];
            
            for i = 1:2
                if (i == 1)
                    L = plat(p).trial.tpoltheta == angs(rot) | softEq(plat(p).trial.tpoltheta,angs(rot)+pi);
                else
                    L = softEq(plat(p).trial.tpoltheta,angs(rot)+pi/2) | softEq(plat(p).trial.tpoltheta,angs(rot)+3*pi/2);
                end
                projs = plat(p).trial.tpolrho(L);
                topfr = max(plat(p).trial.nspikes(L));
                if (topfr == 0)
                    f(i,:) = [0 0 1 1 2 2 0];
                else
                    params1 = [topfr, topfr, sigmaguess, sigmaguess, 2, 2, 0];  % need to constrain B across color directions
                    f(i,:) = fmincon('MyFitNakaRushtonFun',params1,[],[],[],[],vlb,vub,[],options,projs,plat(p).trial.nspikes(L),'asymmetric');
                end
                
                % May want to use all data to generate 1-d nakarushtons, which
                % in turn generate initial guesses...
                pred = MyComputeNakaRushton(f(i,:),projs,'asymmetric');
                
            end
            
            % Redefine the index as all 4 directions
            L = plat(p).trial.tpoltheta == angs(rot) | softEq(plat(p).trial.tpoltheta,angs(rot)+pi)...
                | softEq(plat(p).trial.tpoltheta,angs(rot)+pi/2) | softEq(plat(p).trial.tpoltheta,angs(rot)+3*pi/2);
            
            a = mean(f(1,[1 2]));
            b = mean(f(2,[1 2]));
            params1 = max(plat(p).trial.nspikes);
            %params1(2) = a./sum([a b]); %this doesn't work so well
            params1(2) = .5;
            params1(3) = .25;
            params1(4) = 2;
            params1(5) = (f(1,2)./f(1,4))./(f(1,1)./f(1,3)); if (isnan(params1(5))) params1(5) = 1; end;
            params1(6) = (f(2,2)./f(2,4))./(f(2,1)./f(2,3)); if (isnan(params1(6))) params1(6) = 1; end;
            params1(7) = mean(f(:,end));
            
            vlb = [10   0 0.0001 1  0  0  0];
            vub = [1000 1     1  6 10 10 20];
            options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
            rotMat = [cos(angs(rot)) -sin(angs(rot)); sin(angs(rot)) cos(angs(rot))];
            rotPts = plat(p).trial.tstim*rotMat;
            [f1] = fmincon('MyFitNakaRushtonFun',params1,[],[],[],[],vlb,vub,[],options,rotPts(L,:),plat(p).trial.nspikes(L),'surface1');
            [x,y] = meshgrid(linspace(min(plat(p).trial.LpM),max(plat(p).trial.LpM),50), linspace(min(plat(p).trial.LmM),max(plat(p).trial.LmM),50));
            surface = MyComputeNakaRushton(f1,[x(:) y(:)],'surface1');
            
            % Evaluate goodness of fit
            gof2Ax{rot,p,loop} = -MyFitNakaRushtonFun(f1,rotPts,plat(p).trial.nspikes,'surface1');
            
            %         % Plot Surface
            %         figure(fig); fig=fig+1; hold on; grid on;
            %         surf(x,y,reshape(surface,size(x)))
            %         axis([-1 1 -1 1])
            %         plot3(rotPts(L,1),rotPts(L,2),plat(p).trial.nspikes(L),'k*')
            %         plot3(rotPts(~L,1),rotPts(~L,2),plat(p).trial.nspikes(~L),'go');
            %         xlabel('L-M')
            %         ylabel('L+M')
            %         zlabel('Number of Spikes')
            %         title(['Fit to Only 2 Vectors: ',num2str(angs(rot)+pi/2),' and ',num2str(angs(rot)+pi)])
            
            
        end
        
        likelihoodratio1(p,loop) = min(cat(1,gof2Ax{:,p,loop}))/max(cat(1,gof2Ax{:,p,loop}));
        
        % Plot Goodness of Fit
        figure(fig); fig=fig+1; hold on; grid on;
        plot(angs+pi/2,cat(1,gof2Ax{:,p,loop}),'o-')
        title('Fit From Two Vectors')
        xlabel('Rotation (Radians)')
        ylabel('Log Likelihood')
        
        
    end
    
end

figure(fig); fig = fig+1; hold on; grid on;
plot(1:ndatafiles, likelihoodratio0(1,:), 'b')
plot(1:ndatafiles, likelihoodratio0(2,:), 'b--')
plot(1:ndatafiles, likelihoodratio1(1,:), 'g')
plot(1:ndatafiles, likelihoodratio1(2,:), 'g--')
xlabel('Datafile')
ylabel('min(GOF)/max(GOF)')



