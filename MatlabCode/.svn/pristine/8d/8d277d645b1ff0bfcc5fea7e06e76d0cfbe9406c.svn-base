function GridLMPlaneOnline_MJ()
global udpCom tcounter trial datastruct

disp('In GridLMPlaneOnline_MJ')

udpCom.sock = [];
udpCom.port = 6665;
udpCom.rexip = '192.168.1.120';
udpCom.plexip = '192.168.1.121';
udpCom.mac = '192.168.1.122';

C = GridLMPlaneCodes();% defined in a separate m file
f
p = InitPStruct(0, C.HDRCOMPLETECD);
s = InitPlex();
initGUI();
initTrialStruct();
initDataStruct();


% Setting up variables for the online GUI
% Some initializations
TS = TrialSpecCodes;
currentspiketally = [];
stimon_t = [];
fix_t = [];


% The main loop
socketOpen = 0;
while ~socketOpen
    [udpCom.sock, socketOpen] = pnetStart(udpCom.port);
end
plxServer = InitPlex();

bounceOut = false;
while ~bounceOut
    if CheckForESCKey() || dealWithMsgs(udpCom.sock)
        bounceOut = true;
    end
    
    [n, eventList] = PL_GetTS(s);
    if n
        p = ProcessEventList(p, eventList);
    end
    
    % do your magic hither
    
    if (any(p.events == C.FPACQCD))
        fix_t = p.times(find(p.events == C.FPACQCD,1));
        abortflag = 0;
    end
    if (any(p.events == C.STIMONCD))
        stimon_t = p.times(find(p.events == C.STIMONCD,1));
        nspikes = sum(currentspiketally > 0); % Spike times are relative to fix_t
        %gl.trialspecs{gl.paramsidx,TS.BASELINE} =...
        %    [nspikes/(stimon_t-fix_t)*1000];
        %[gl.trialspecs{gl.paramsidx,TS.BASELINE}; nspikes/(stimon_t-fix_t)*1000];
        p.lastprocessed_t = p.times(find(p.events == C.STIMONCD,1));
    end
    if (any(p.spikes{1}) && ~isempty(fix_t) && ~abortflag)
        currentspiketally = [currentspiketally; p.spikes{1}-fix_t];
        PlotRaster(currentspiketally);
        p.spikes{1} = [];
    end
    if (any(p.events == C.ABORTCD))
        abortflag = 1;
        currentspiketally = [];
        p.lastprocessed_t = p.times(find(p.events == C.ABORTCD,1,'last'));
    end
    if (any(p.events == C.FPOFFCD))
        
        stimoff_t = p.times(find(p.events == C.FPOFFCD,1));
        trial.fpacq(tcounter) = fix_t;
        trial.stimon(tcounter) = stimon_t;
        trial.stimoff(tcounter) = stimoff_t;
        trial.duration(tcounter) = trial.stimoff(tcounter)-stimon_t;
        trial.tspikes{tcounter} = currentspiketally + fix_t;
        trial.normtspikes{tcounter} = currentspiketally;
        trial.nspikes(tcounter) = sum(currentspiketally > (stimon_t-fix_t) & currentspiketally < (stimoff_t-fix_t));
        trial.fr(tcounter) = trial.nspikes(tcounter)/trial.duration(tcounter)*1000;
        
        % Spike Historgram
        PlotSpikeHist
        
        % Clear out fields
        currentspiketally = [];
        stimon_t = [];
        fix_t = [];        
        p.lastprocessed_t = stimoff_t;
        
        % Organize, Save, and Plot Data
        organizeParStruct();
        plat = organizePlatStruct();
        plotMJGUI();
        plotJPWGUI(plat);
        
        % Update tcounter before next trial begins
        tcounter = tcounter + 1;
        
        disp('Sending plexdone message...')
        pnet(udpCom.sock, 'write', 'plexdone>> >>');
        pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
        disp('plexdone message sent')
        
    end
    p = CleanUpEvents(p);
end

PL_Close(plxServer);
end

 
function PlotRaster(spiketally)
global trial

    if any(~isnan(trial.duration))
        xmax = max(trial.duration);
    else
        xmax = 3000;
    end
    
    figure(1);
    subplot(2,9,[1 2]); cla; hold on;
    axis([-500 xmax -.5 .5])
    axis square
    plot([spiketally, spiketally]',[-.5; .5]*ones(1,length(spiketally)),'k-');
    title('Incoming Spikes')
    
end


function PlotSpikeHist
global trial tcounter

    % Set up variables
    bins = linspace(-500,max(trial.duration),75);
    PSTH = zeros(1,length(bins));
    circDelay = .05;

    % Generate spiking profile across all trials
    for i = 1:tcounter
        if ~isempty(trial.normtspikes{i})
            PSTH = PSTH + histc(trial.normtspikes{i}',bins);
        end
    end

    % Find baseline spikerate
    nsp =  sum(trial.tspikes{tcounter} > trial.fpacq(tcounter) + circDelay & trial.tspikes{tcounter} < trial.stimon(tcounter));
    dt = trial.stimon(tcounter) - trial.fpacq(tcounter);
    trial.baselinefr(tcounter) = nsp./(dt - circDelay);

    % Fit Gaussian
    gauss = @(x,params)(params(2)+params(1)*exp((-(x-params(3)).^2)/params(4).^2));
    bguess = mean(PSTH(bins<.25));
    aguess = max(PSTH)-bguess;
    muguess = bins(find(PSTH == max(PSTH),1,'first'));
    sigmaguess = .1; % terrible. Need something better.
    fittedparams = fminsearch(@(params)sum((gauss(bins,params)-PSTH).^2),[aguess,bguess,muguess,sigmaguess]);
    offset = [fittedparams(3)-fittedparams(4) fittedparams(3)+fittedparams(4)];

    % Plot spiking profile
    figure(1);
    subplot(2,9,[10 11]); cla; hold on;
    xlim([-500 max(trial.duration)])
    bar(bins,PSTH);
    plot(bins,gauss(bins,fittedparams),'b-.','LineWidth',3);
    plot(offset,[fittedparams(1) fittedparams(1)],'k-','linewidth',3);
    xlabel('Normalized time (ms)')
    ylabel('Spike Count')
    title('Spike Histogram')
    axis square

    % Now counting up spikes in a window
%     for i = 1:tcounter
%         trial.nspikes(i) = sum(trial.normtspikes{i} > offset(1) & trial.normtspikes{i} < offset(2));
%         dt = offset(2) - offset(1);
%         trial.fr(i) = trial.nspikes(i)./dt;
%     end

end


function allDone = dealWithMsgs(socket)
global udpCom datastruct
allDone = false;
msgSize = pnet(socket, 'readpacket', 250, 'noblock');
if ~msgSize, return; end

message = pnet(socket, 'read', msgSize, 'char');
if ~isempty(message)
    [message,allDone] = parseMsg(message,socket);
    if strncmp(message, 'return', 6)
        if (~isempty(datastruct))
            % New file name
            %MJdir = 'N:\NexFiles\Patrick\MJdata\';
            MJdir = 'C:\Documents and Settings\jpatrickweller\My Documents\MATLAB\';
            MJfilenameBase = ['MJout',datestr(date,'yyyymmdd')];
            files = dir([MJdir 'MJout*']);
            if ~isempty(files)
                filenames = cat(1,files.name);
                newsuffix = max(str2num(filenames(:,end-6:end-4)))+1;
                MJfilename = [MJdir MJfilenameBase sprintf('%03d',newsuffix)];
            else
                MJfilename = [MJdir MJfilenameBase '001'];
            end
            disp('saving:');
            MJfilename
            save(MJfilename,'datastruct')
        end
    end
    evalMsg(message);
else
    allDone = true;
end
end


% here we implement the logic/function calls associated with the requested
% variables within the message
function [message,doneFlag] = parseMsg(message, socket)
global udpCom

doneFlag = false;

matches = strtrim(regexp(message, ',', 'split'));
if strncmp(matches{1}, 'sendToRex', 9)
    varToSend = matches{2}; % the second match contains the requested variable(s)
    if strfind(varToSend, 'stimColor')
        keyboard;
    elseif strfind(varToSend, 'some.other.variable')
        % do something else
        %     else
        %         eval(message)
    elseif strncmp(message, 'return', 6) % not sure what catastrophe this prevents
        stk = dbstack();
        if ~strcmp(stk(end).name, mfilename)
            doneFlag = true;
        end
    end
end

end

function evalMsg(message)
global udpCom gl %#ok<NUSED> eval needs access to gl
%disp('In evalMsg')
try
    eval(message);
catch exception
    fprintf('Trouble with message: "%s"\n', message);
    disp(getReport(exception));
end
end

function [stim] = getStimParams(theta,sigX,sigY,DR,nstd,nexp)
global trial tcounter

disp('In getStimParams')

%Set up variables
stimDim = 2;
np = 25;

% New MJ parameters (for 2D surface)
maxsupport = 1;
minsupport = -1;
stimrange = linspace(minsupport, maxsupport, np);
[xx, yy] = meshgrid(stimrange, stimrange);
support = [xx(:) yy(:)]; % 625 points on the input grid

% number of initial stimuli
numInitData = stimDim*10;

% initial stimuli selected from LH design
x_lh = lhsdesign(numInitData, stimDim);
xInit = bsxfun(@plus, bsxfun(@times, (maxsupport - minsupport), x_lh), minsupport);

% Old MJ parameters (for 2D surface)
%support_x = linspace(-1, 1, np);
%support_y = linspace(-1, 1, np);
%[xx, yy] = meshgrid(support_x, support_y);
%support = [xx(:) yy(:)];
%numInitData = stimDim*10;
%whichMethod = 1;
stim_range_min = min(support);
stim_range_max = max(support);
    
if tcounter == 1
    x_lh = lhsdesign(numInitData, stimDim);
    trial.xInit = bsxfun(@plus, bsxfun(@times, (stim_range_max - stim_range_min), x_lh), stim_range_min);
end

if tcounter <= numInitData
    nextX = trial.xInit(tcounter,:);
else   
    [nextX, totData] = computeNextStim_ALalgorithm_NewModel([trial.Lcc_norm trial.Mcc_norm], trial.nspikes, stimDim, support, numInitData, maxsupport, minsupport);
    %[nextX, totData] = computeNextStim_ALalgorithm([trial.Lcc_norm trial.Mcc_norm], trial.nspikes, stimDim, support, numInitData, whichMethod);
end

% Some defaults
Lcc = nextX(1);
Mcc = nextX(2);
Scc = 0;

trial.Lcc_norm(tcounter) = Lcc;
trial.Mcc_norm(tcounter) = Mcc;
trial.nstim(tcounter,:) = [Lcc Mcc];
trial.Scc(tcounter) = Scc;
trial.theta(tcounter) = theta;
trial.sigX(tcounter) = sigX;
trial.sigY(tcounter) = sigY;
trial.driftRate(tcounter) = DR;
trial.nstd(tcounter) = nstd;
trial.nexp(tcounter) = nexp;

transMat = [.08 .08; -.08 .08];
LM = transMat * [Lcc Mcc]';
trial.Lcc_orig(tcounter) = LM(1);
trial.Mcc_orig(tcounter) = LM(2);
trial.LmM(tcounter) = trial.nstim(tcounter,1);
trial.LpM(tcounter) = trial.nstim(tcounter,2);
trial.ostim(tcounter,:) = LM;

stim = [trial.Lcc_orig(tcounter) trial.Mcc_orig(tcounter) trial.Scc(tcounter) trial.theta(tcounter)...
    trial.sigX(tcounter) trial.sigY(tcounter) trial.driftRate(tcounter) trial.nstd(tcounter) trial.nexp(tcounter)];

% Display stimulus parameters, etc.
fprintf('Stimuli already tested = %i\n',tcounter)
fprintf('Current Stimulus Parameters:\n Lcc = %g\n Mcc = %g\n Scc = %g\n Theta = %g\n Sigma X = %g\n Sigma Y = %g\n Drift Rate = %g\n nStd = %g\n nExpanse = %g\n\n', stim);

end


% function [nextX, totData] = computeNextStim_ALalgorithm(x, r, stimDim, support, numInitData, whichMethod)
% global datastruct
% tic;
% 
% thTrial = length(r);
% r(r==0) = 0.1;
% 
% if thTrial == numInitData
%     datastruct.x = x;
%     datastruct.r = r;
%     datastruct.support = support;
%     datastruct.ndim = stimDim;
%     datastruct.norm_mat_support = form_normMat(support, support);  % squared distance
%     
%     g = @(t) log(exp(t)+1);
%     ginv = @(t) log(exp(t)-1);
%     datastruct.g = g;
%     datastruct.ginv = ginv;
%     datastruct.finit = datastruct.ginv(datastruct.r+0.1);
%     
%     datastruct.thrsh_detH = 0.01;
%     datastruct.detH = 1;
%     
%     datastruct.norm_mat = form_normMat(datastruct.x, datastruct.x);  % squared distance
%     datastruct.norm_mat_Kstar = form_normMat(datastruct.support, datastruct.x);
%     
%     load fvar_logexp1_lin.mat;
%     load fmean_logexp1_lin.mat;
%     
%     datastruct.fvar_logexp1 = fvar_logexp1_lin;
%     datastruct.fmean_logexp1 = fmean_logexp1_lin;
%     
%     datastruct.whichMethod = whichMethod;
%     
%     load datastruct2Dtable;
%     
%     datastruct.fxx = datastruct2Dtable.fxx;
%     datastruct.sigxx = datastruct2Dtable.sigxx;
%     datastruct.val = datastruct2Dtable.val;
%     
% else
%     datastruct.x = x;
%     datastruct.r = r;
%     
%     try
%         sqrDist_new = form_normMat(datastruct.x(end,:), datastruct.x);
%         datastruct.norm_mat = [[datastruct.norm_mat; sqrDist_new(1:end-1)] sqrDist_new'];
%         
%         normMat_Kstar_new = form_normMat(datastruct.support, datastruct.x(end,:));
%         datastruct.norm_mat_Kstar = [datastruct.norm_mat_Kstar normMat_Kstar_new];
%         datastruct.finit = [datastruct.finit; datastruct.ginv(r(end)+0.1)];
%     catch
%         datastruct.norm_mat = form_normMat(datastruct.x, datastruct.x);  % squared distance
%         datastruct.norm_mat_Kstar = form_normMat(datastruct.support, datastruct.x);
%         datastruct.finit = datastruct.ginv(datastruct.r+0.1);
%     end
%         
% end
% 
% % update data structure for the rest
% datastruct.nstim = length(datastruct.r);
% maxNumTrial = 100; % maximum number of trials for updating hyperparameters
% 
% if ((rem(thTrial, 10)==0) &&(datastruct.detH>datastruct.thrsh_detH))&&(thTrial<=maxNumTrial)
% 
%     % optimize hyperparameters with analytic form
%     ovrscl_1 = mean(datastruct.r)/2; % overall scale
%     lngthscl_1 = max(max(support))-min(min(support))/2; % variance
%     prs0 = [mean(datastruct.r); ovrscl_1; lngthscl_1];
%     datastruct.K = abs(prs0(2))*exp(-.5/abs(prs0(3)).*datastruct.norm_mat);
%     datastruct.muf = abs(prs0(1));
%     
%     if thTrial == numInitData % update hyperparameters using the initial data
%         [prs, fmapFinal, aFinal, WFinal, sqrtLFinal, neglogev0, detH] = computeFmapAndUpdateTheta(prs0, datastruct);
%         datastruct.detH=detH;
%         datastruct.prs = prs;
%         datastruct.finit = fmapFinal;
%         datastruct.muf = prs(1);
% 
%     else % check if we need to update hyperparameters
%         neglogev0  = updateFmapGivenK(datastruct);
%         if datastruct.neglogev<neglogev0
%             prs = datastruct.prs;
%             datastruct.K = abs(prs(2))*exp(-.5/abs(prs(3)).*datastruct.norm_mat);
%             datastruct.muf = abs(prs(1));
%             [neglogev, fmapFinal, aFinal, WFinal, sqrtLFinal]  = updateFmapGivenK(datastruct);
%             datastruct.neglogev = neglogev;
%             datastruct.finit = fmapFinal;
% 
%         else
%             [prs, fmapFinal, aFinal, WFinal, sqrtLFinal, neglogev0, detH] = computeFmapAndUpdateTheta(prs0, datastruct);
%             datastruct.detH=detH;
%             datastruct.prs = prs;
%             datastruct.finit = fmapFinal;
%             datastruct.muf = prs(1);
%         end
%     end
%  
% else
% 
%     prs = datastruct.prs;
%     datastruct.K = abs(prs(2))*exp(-.5/abs(prs(3)).*datastruct.norm_mat);
%     datastruct.muf = abs(prs(1));
%     [neglogev, fmapFinal, aFinal, WFinal, sqrtLFinal]  = updateFmapGivenK(datastruct); 
% 
%     datastruct.neglogev = neglogev;
%     datastruct.finit = fmapFinal;
% 
% end
% 
% datastruct.Kstar = prs(2)*exp(-.5/prs(3).*datastruct.norm_mat_Kstar);
% [predictiveMean, predictiveVar, nextX] = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal, datastruct.fvar_logexp1, datastruct.whichMethod);
% 
% toc;
% 
% qz = datastruct.fmean_logexp1(predictiveMean, sqrt(predictiveVar));
% 
% datastruct.lambFinal = qz;
% datastruct.aFinal = aFinal;
% datastruct.WFinal = WFinal;
% datastruct.sqrtLFinal = sqrtLFinal;
% 
% totData = datastruct;
% 
% % tic;
% % save totData totData; 
% % toc;
% 
% end


function out = TrialSpecCodes
out.TRIALCOUNT = 1;
out.SPIKES = 2;
out.BASELINE = 3;
out.DIAM = 4;
out.SF = 5;
out.TF = 6;
out.CC = 7;
out.ORIENT = 8;
out.PHI = 9;
out.NCYCLES = 10;
out.STIMTYPE = 11;
out.PROTOCOL = 12;
end


function p = CleanUpEvents(p)

if (p.processnowcode ~= 0)
    pncodeidx = find(p.events == p.processnowcode,1);
    p.lastprocessed_t = p.times(pncodeidx);
end
%disp('In CleanUpEvents')
%p.times
%p.lastprocessed_t
L = p.times <= p.lastprocessed_t;

p.times(L) = [];
p.events(L) = [];
for i = 1:length(p.spikes)
    spiketimevect = p.spikes{i};
    p.spikes{i} = [spiketimevect(spiketimevect > p.lastprocessed_t)];
end
p.processnowflag = 0;
end


function initGUI()

close all

end


function plotMJGUI()
global tcounter trial datastruct

if isempty(datastruct)
    
    figure(1);
    subplot(2,9,[4 5]); cla; hold on; grid on;
    plot(trial.xInit(1:tcounter, 1), trial.xInit(1:tcounter,2), 'mo', 'LineWidth',2,...
        'MarkerEdgeColor','y',...
        'MarkerFaceColor',[.2 0.2 .6],...
        'MarkerSize',8);
    axis square
    set(gca, 'xlim', [-1 1], 'ylim', [-1 1]);
    xlabel('Normalized L-M')
    ylabel('Normalized L+M')
    title('Stimulus Selection')
    
else
        
    np = sqrt(size(datastruct.support,1));
    supportside = linspace(-1, 1, np);
    [xx, yy] = meshgrid(supportside, supportside);
    
    % Chosen Stim
    figure(1);
    subplot(2,9,[4 5]); cla; hold on; grid on;
    contour(xx, yy, reshape(datastruct.lambFinal, length(xx), []), 5); axis image; axis xy; title('estimate');
    plot(datastruct.x(1:end, 1), datastruct.x(1:end,2), 'mo', 'LineWidth',2,...
        'MarkerEdgeColor','y',...
        'MarkerFaceColor',[.2 0.2 .6],...
        'MarkerSize',8);
    set(gca, 'xlim', [-1 1], 'ylim', [-1 1]);
    xlabel('Norubpmalized L-M')
    ylabel('Normalized L+M')
    title('Stimulus Selection')
    
    % FR Contour Surface
    subplot(2,9,[13 14]); cla; hold on; grid on;
    contour(xx, yy, reshape(datastruct.lambFinal, length(xx), []), 5); axis image; axis xy; title('estimated firing map');
    set(gca, 'xlim', [-1 1], 'ylim', [-1 1]);
    xlabel('Normalized L-M')
    ylabel('Normalized L+M')
    title('Surface Contour')

end

end


function plotJPWGUI(plat)
global tcounter

nplat = numel(plat);

    for p = 1:nplat
    
        % 3-D Stem Plot
        figure(1)
        subplot(2,9,[7 8 9 16 17 18]); cla; hold on;
        maxresp = max(plat(p).par.meanfr);
        if maxresp == 0
             maxresp = 5;
        end
        grid on; axis([-1 1 -1 1 0 maxresp])
        axis([-1 1 -1 1])
        xlabel('L-M')
        ylabel('L+M')
        zlabel('Firing Rate (sp/s)')
        title('3D Stem Plot')
        stem3(plat(p).par.LmM,plat(p).par.LpM,plat(p).par.meanfr,'fill','MarkerFaceColor','c');
        if tcounter == 1
            set(gca,'view',[-75 38]);
        else
            a = get(gca,'View');
            set(gca,'View',a);
        end
    
    end

end


function initTrialStruct()
global trial tcounter

    trial.Lcc_norm = nan(2,1);
    trial.Mcc_norm = nan(2,1);
    trial.nstim = nan(2,2);
    %trial.npolTheta = nan(2,1);
    %trial.npolRho = nan(2,1);
    trial.Lcc_orig = nan(2,1);
    trial.Mcc_orig = nan(2,1);
    trial.ostim = nan(2,2);
    %trial.opolTheta = nan(2,1);
    %trial.opolRho = nan(2,1);
    trial.Scc = nan(2,1);
    trial.LmM = nan(2,1);
    trial.LpM = nan(2,1);
    trial.theta = nan(2,1);
    trial.sigX = nan(2,1);
    trial.sigY = nan(2,1);
    trial.driftRate = nan(2,1);
    trial.nstd = nan(2,1);
    trial.nexp = nan(2,1);
    trial.fpacq = nan(2,1);
    trial.stimon = nan(2,1);
    trial.stimoff = nan(2,1);
    trial.duration = nan(2,1);
    trial.tspikes = cell(2,1);
    trial.nspikes = nan(2,1);
    trial.fr = nan(2,1);
    trial.baselinefr = nan(2,1);
    trial.normtspikes = cell(2,1);

    tcounter = 1;
end


function organizeParStruct()
global trial par

    % By Paradigm Properties (aka by unique stimulus)
    [uniqueCond,m,idx] = unique([trial.Lcc_orig trial.Mcc_orig trial.Scc trial.theta trial.sigX trial.sigY trial.driftRate trial.nstd trial.nexp],'rows');
    
    if any(isnan(uniqueCond(:,1)))
        nanIdx = find(isnan(uniqueCond(:,1)));
        uniqueCond(nanIdx,:) = [];
        m(nanIdx) = [];
        idx(nanIdx) = [];
    end
    
    % Collapse repeated stimuli
    par.Lcc_orig = uniqueCond(:,1);
    par.Mcc_orig = uniqueCond(:,2);
    par.Scc = uniqueCond(:,3);
    par.ostim = [par.Lcc_orig par.Mcc_orig];
    %par.opolTheta = trial.opolTheta(m);
    %par.opolRho = trial.opolRho(m);
    %[par.opolTheta par.opolRho] = cart2pol(par.Lcc,par.Mcc);
    %par.nstim = rndofferr([transMat * par.stim(:,1:2)']',3);
    par.Lcc_norm = trial.Lcc_norm(m);
    par.Mcc_norm = trial.Mcc_norm(m);
    par.nstim = [par.Lcc_norm par.Mcc_norm];
    %par.npolTheta = trial.npolTheta(m);
    %par.npolRho = trial.npolRho(m);
    par.LmM = par.nstim(:,1);
    par.LpM = par.nstim(:,2);
    %[par.opolTheta par.opolRho] = cart2pol(par.Lcc,par.Mcc);

    % Correcting thetas (roundoff error)
    %[par.npolTheta par.npolrho] = cart2pol(par.nstim(:,1),par.nstim(:,2));
%     divs = -pi-pi/64:pi/64:pi+pi/64;
%     rads = divs(2:2:end);
%     edges = divs(1:2:end);
%     [n,binIdx] = histc(par.npolTheta,edges);
%     par.npolTheta = rads(binIdx)';
%     if any(trial.npolTheta==-pi)
%         trial.npolTheta(trial.npolTheta==-pi)=pi;
%     end

    par.theta = uniqueCond(:,4);
    par.sigX = uniqueCond(:,5);
    par.sigY = uniqueCond(:,6);
    par.driftRate = uniqueCond(:,7);
    par.nstd = uniqueCond(:,8);
    par.nexp = uniqueCond(:,9);

    % Preallocate Space
    par.frs = cell(size(uniqueCond,1),1);
    par.meanfr = nan(size(uniqueCond,1),1);
    par.varfrs = nan(size(uniqueCond,1),1);
    %par.baselinefrs = cell(size(uniqueCond,1),1);
    %par.meanblfr = nan(size(uniqueCond,1),1);
    par.normtspikes = cell(size(uniqueCond,1),1);
    par.catntspikes = cell(size(uniqueCond,1),1);

    for n = 1:size(par.nstim,1)
        L = idx == n;
        par.frs{n} = trial.fr(L)';
        par.meanfr(n) = mean(par.frs{n});
        par.varfrs(n) = var(par.frs{n});
        %par.baselinefrs{n} = trial.baselinefr(L);
        %par.meanblfr(n) = mean(par.baselinefrs{n});
        par.normtspikes{n} = trial.normtspikes(L);
        par.catntspikes{n} = cat(1,trial.normtspikes{L});
    end
end


function plat = organizePlatStruct()
global trial par 
    
    % By paradigm
    [uniqueCond,m,idx] = unique([par.Scc par.theta par.sigX par.sigY par.driftRate par.nstd par.nexp],'rows');
    nplat = size(uniqueCond,1);
    plat(1:nplat) = struct();

    for p = 1:nplat
        L = idx == p;
        plat(p).par.Lcc_orig = par.Lcc_orig(L,:);
        plat(p).par.Mcc_orig = par.Mcc_orig(L,:);
        plat(p).par.Scc = par.Scc(L,:);
        plat(p).par.ostim = par.ostim(L,:);
        %plat(p).par.opolTheta = par.opolTheta(L);
        %plat(p).par.opolRho = par.opolRho(L);
        plat(p).par.Lcc_norm = par.Lcc_norm(L);
        plat(p).par.Mcc_norm = par.Mcc_norm(L);
        plat(p).par.nstim = par.nstim(L,:);
        %plat(p).par.npolTheta = par.npolTheta(L);
        %plat(p).par.npolRho = par.npolRho(L);
        plat(p).par.LmM = par.LmM(L);
        plat(p).par.LpM = par.LpM(L);
        plat(p).par.theta = par.theta(L,:);
        plat(p).par.sigX = par.sigX(L,:);
        plat(p).par.sigY = par.sigY(L,:);
        plat(p).par.driftRate = par.driftRate(L,:);
        plat(p).par.nstd = par.nstd(L,:);
        plat(p).par.nexp = par.nexp(L,:);
        plat(p).par.frs = {par.frs{L}}';
        plat(p).par.meanfr = par.meanfr(L);
        plat(p).par.varfrs = par.varfrs(L);
        %plat(p).par.baselinefrs = par.baselinefrs(L);
        %plat(p).par.meanblfr = par.meanblfr(L);
        plat(p).par.normtspikes = par.normtspikes(L);
        plat(p).par.catntspikes = par.catntspikes(L);
    end


    % By trial
    [uniqueCond,m,idx] = unique([trial.Scc trial.theta trial.sigX trial.sigY trial.driftRate trial.nstd trial.nexp],'rows');

    for p = 1:nplat
        L = idx == p;
        plat(p).trial.Lcc_orig = trial.Lcc_orig(L);
        plat(p).trial.Mcc_orig = trial.Mcc_orig(L);
        plat(p).trial.Scc = trial.Scc(L);
        plat(p).trial.ostim = trial.ostim(L,:);
        %plat(p).trial.opolTheta = trial.opolTheta(L);
        %plat(p).trial.opolRho = trial.opolRho(L);
        plat(p).trial.Lcc_norm = trial.Lcc_norm(L);
        plat(p).trial.Mcc_norm = trial.Mcc_norm(L);
        plat(p).trial.nstim = trial.nstim(L,:);
        %plat(p).trial.npolTheta = trial.npolTheta(L);
        %plat(p).trial.npolRho = trial.npolRho(L);
        plat(p).trial.LmM = trial.LmM(L);
        plat(p).trial.LpM = trial.LpM(L);
        plat(p).trial.theta = trial.theta(L);
        plat(p).trial.sigX = trial.sigX(L);
        plat(p).trial.sigY = trial.sigY(L);
        plat(p).trial.driftRate = trial.driftRate(L);
        plat(p).trial.nstd = trial.nstd(L);
        plat(p).trial.nexp = trial.nexp(L);
        plat(p).trial.fr = trial.fr(L);
        %plat(p).trial.baselinefr = trial.baselinefr(L);
        %plat(p).trial.normtspikes = trial.normtspikes(L);
        plat(p).trial.nspikes = trial.nspikes(L);
        plat(p).trial.tspikes = trial.tspikes(L);
    end
end


function initDataStruct()
global datastruct

datastruct = [];

end



