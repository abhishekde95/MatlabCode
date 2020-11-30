function GridLMPlaneOnline_Hybrid()
global udpCom tcounter trial datastruct

disp('In GridLMPlaneOnline_Hybrid')

udpCom.sock = [];
udpCom.port = 6665;
udpCom.rexip = '192.168.1.120';
udpCom.plexip = '192.168.1.121';
udpCom.mac = '192.168.1.122';

C = GridLMPlaneCodes();% defined in a separate m file

p = InitPStruct(0, C.HDRCOMPLETECD);
s = InitPlex();
resetVars();
initGUI();
initTrialStruct();
initParStruct();
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
        if tcounter>1
            plat = organizePlatStruct();
        end
        if mod(tcounter,2) == 1
            plotMJGUI();
        else
            plotJPWGUI();
        end
        
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
global trial tcounter

    if any(~isnan(trial.duration))
        xmax = max(trial.duration);
    else
        xmax = 3000;
    end
   
    % Plot spiking profile
    if mod(tcounter,2) == 1
        figure(1);
        subplot(2,9,[1 2]); cla; hold on;
    else
        figure(2);
        subplot(2,5,[1 2]); cla; hold on;
    end
    posSP = spiketally(spiketally > 0);
    %axis([-500 xmax -.5 .5])
    axis([0 xmax -.5 .5])
    axis square
    %plot([spiketally, spiketally]',[-.5; .5]*ones(1,length(spiketally)),'k-');
    plot([posSP, posSP]',[-.5; .5]*ones(1,length(posSP)),'k-');
    title('Incoming Spikes')

end


function PlotSpikeHist
global trial tcounter

    % Set up variables
    bins = linspace(-500,max(trial.duration),75);
    PSTH = zeros(1,length(bins));
    circDelay = .05;

    % Generate spiking profile across all trials - different for GLMP and ALMP
    if mod(tcounter,2) == 1
        L = 1:2:tcounter;
    else
        L = 2:2:tcounter;
    end
    for i = L
        if ~isempty(trial.normtspikes{i})
            PSTH = PSTH + histc(trial.normtspikes{i}',bins);
        end
    end

    % Find baseline spikerate
    nsp =  sum(trial.tspikes{tcounter} > trial.fpacq(tcounter) + circDelay & trial.tspikes{tcounter} < trial.stimoff(tcounter));
    dt = trial.stimon(tcounter) - trial.fpacq(tcounter);
    trial.baselinefr(tcounter) = nsp./(dt - circDelay);

    % Fit Gaussian
%     gauss = @(x,params)(params(2)+params(1)*exp((-(x-params(3)).^2)/params(4).^2));
%     bguess = mean(PSTH(bins<.25));
%     aguess = max(PSTH)-bguess;
%     muguess = bins(find(PSTH == max(PSTH),1,'first'));
%     sigmaguess = .1; % terrible. Need something better.
%     fittedparams = fminsearch(@(params)sum((gauss(bins,params)-PSTH).^2),[aguess,bguess,muguess,sigmaguess]);
%     offset = [fittedparams(3)-fittedparams(4) fittedparams(3)+fittedparams(4)];

    
    % Plot spiking profile
    if mod(tcounter,2)==1
        figure(1);
        subplot(2,9,[10 11]); cla; hold on;
        xlim([-500 max(trial.duration)])
        bar(bins,PSTH);
        %plot(bins,gauss(bins,fittedparams),'b-.','LineWidth',3);
        %plot(offset,[fittedparams(1) fittedparams(1)],'k-','linewidth',3);
        xlabel('Normalized time (ms)')
        ylabel('Spike Count')
        title('Spike Histogram')
        axis square
    else
        figure(2);
        subplot(2,5,[6 7]); cla; hold on;
        xlim([-500 max(trial.duration)])
        bar(bins,PSTH);
        %plot(bins,gauss(bins,fittedparams),'b-.','LineWidth',3);
        %plot(offset,[fittedparams(1) fittedparams(1)],'k-','linewidth',3);
        xlabel('Normalized time (ms)')
        ylabel('Spike Count')
        title('Spike Histogram')
    end
        

    % Now counting up spikes in a window
%     for i = 1:tcounter
%         trial.nspikes(i) = sum(trial.normtspikes{i} > offset(1) & trial.normtspikes{i} < offset(2));
%         dt = offset(2) - offset(1);
%         trial.fr(i) = trial.nspikes(i)./dt;
%     end

end


function allDone = dealWithMsgs(socket)
global udpCom datastruct tcounter

    allDone = false;
    msgSize = pnet(socket, 'readpacket', 250, 'noblock');
    if ~msgSize, return; end

    message = pnet(socket, 'read', msgSize, 'char');
    if ~isempty(message)
        [message,allDone] = parseMsg(message,socket);
        if strncmp(message, 'return', 6)
            if (~isempty(datastruct)) && mod(tcounter,2)==1
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

    try
        eval(message);
    catch exception
        fprintf('Trouble with message: "%s"\n', message);
        disp(getReport(exception));
    end
    
end


function [stim] = getStimParams(theta,sigX,sigY,DR,nstd,nexp)
global tcounter

    tcounter = tcounter+1;

    if mod(tcounter,2) == 1
        [params] = ALMP_StimSelection(theta,sigX,sigY,DR,nstd,nexp);
        adaptive = 1;
    else
        [params] = GLMP_StimSelection(theta,sigX,sigY,DR,nstd,nexp);
        adaptive = 0;
    end
    
    stim = [params adaptive];
    
end


function [stim] = ALMP_StimSelection(theta,sigX,sigY,DR,nstd,nexp)
global trial tcounter

    disp('Selecting Next Stimulus Using ALMP Paradigm...')

    %Set up variables
    stimDim = 2;
    np = 25;

    % New MJ parameters (for 2D surface)
     maxsupport = 1;
     minsupport = -1;
     stimrange = linspace(minsupport, maxsupport, np);
     [xx, yy] = meshgrid(stimrange, stimrange);
     support = [xx(:) yy(:)]; % 625 points on the input grid
    
    % Trying new grid (radial) JPW
%    thetas = linspace(pi/32,2*pi,64)';
%    rhos = linspace(.05,1,15)';
%    polarIdx = fullfact([numel(thetas),numel(rhos)]);
%    [x,y] = pol2cart(thetas(polarIdx(:,1)),rhos(polarIdx(:,2)));
%    x = rndofferr(x,3);
%    y = rndofferr(y,3);
%    support = [x y]; %640 points on the input grid

    % number of initial stimuli
    numInitData = stimDim*10;

    % initial stimuli selected from LH design
    %x_lh = lhsdesign(numInitData, stimDim);
    %xInit = bsxfun(@plus, bsxfun(@times, (maxsupport - minsupport), x_lh), minsupport);

    stim_range_min = min(support);
    stim_range_max = max(support);
    
    if tcounter == 1 || tcounter == 2
        x_lh = lhsdesign(numInitData, stimDim);
        trial.xInit(1:2:numInitData*2,:) = bsxfun(@plus, bsxfun(@times, (stim_range_max - stim_range_min), x_lh), stim_range_min);
    end

    ALMP_StimIdx = strcmp(trial.paradigm,'ALMP');
    if tcounter <= numInitData*2
        nextX = trial.xInit(tcounter,:);
    else   
        [nextX, totData] = computeNextStim_ALalgorithm_NewModel([trial.LmM_norm(ALMP_StimIdx) trial.LpM_norm(ALMP_StimIdx)], trial.nspikes(ALMP_StimIdx), stimDim, support, numInitData, maxsupport, minsupport);
        %[nextX, totData] = computeNextStim_ALalgorithm([trial.Lcc_norm trial.Mcc_norm], trial.nspikes, stimDim, support, numInitData, whichMethod);
    end

    % Some defaults
    %Lcc = nextX(1);
    %Mcc = nextX(2);
    LmM = nextX(1);
    LpM = nextX(2);
    Scc = 0;

    %trial.Lcc_norm(tcounter) = Lcc;
    %trial.Mcc_norm(tcounter) = Mcc;
    trial.LmM_norm(tcounter) = LmM;
    trial.LpM_norm(tcounter) = LpM;
    %trial.stim_norm(tcounter,:) = [Lcc Mcc];
    trial.stim_norm(tcounter,:) = [LmM LpM];
    trial.Scc_orig(tcounter) = Scc;
    trial.theta(tcounter) = theta;
    trial.sigX(tcounter) = sigX;
    trial.sigY(tcounter) = sigY;
    trial.driftRate(tcounter) = DR;
    trial.nstd(tcounter) = nstd;
    trial.nexp(tcounter) = nexp;

    transMat = [.08 .08; -.08 .08];
    %LM = transMat * [Lcc Mcc]';
    LM = transMat * [LmM LpM]';
    trial.Lcc_orig(tcounter) = LM(1);
    trial.Mcc_orig(tcounter) = LM(2);
    %trial.LmM_norm(tcounter) = trial.stim_norm(tcounter,1);
    %trial.LpM_norm(tcounter) = trial.stim_norm(tcounter,2);
    trial.stim_orig(tcounter,:) = LM;
    [polTheta_orig,polRho_orig] = cart2pol(trial.Lcc_orig(tcounter),trial.Mcc_orig(tcounter));
    trial.polTheta_orig(tcounter) = rndofferr(polTheta_orig,3);
    trial.polRho_orig(tcounter) = rndofferr(polRho_orig,3);
    [polTheta_norm,polRho_norm] = cart2pol(trial.LmM_norm(tcounter),trial.LpM_norm(tcounter));
    trial.polTheta_norm(tcounter) = rndofferr(polTheta_norm,3);
    trial.polRho_norm(tcounter) = rndofferr(polRho_norm,3);
    trial.paradigm(tcounter) = cellstr('ALMP');

    stim = [trial.Lcc_orig(tcounter) trial.Mcc_orig(tcounter) trial.Scc_orig(tcounter) trial.theta(tcounter)...
        trial.sigX(tcounter) trial.sigY(tcounter) trial.driftRate(tcounter) trial.nstd(tcounter) trial.nexp(tcounter)];

    % Display stimulus parameters, etc.
    fprintf('Stimuli already tested = %i\n',tcounter-1)
    fprintf('Current Stimulus Parameters:\n Lcc = %.15f\n Mcc = %.15f\n Scc = %g\n Theta = %g\n Sigma X = %g\n Sigma Y = %g\n Drift Rate = %g\n nStd = %g\n nExpanse = %g\n\n', stim);

end

function [stim] = GLMP_StimSelection(theta,sigX,sigY,DR,nstd,nexp)
global currentQueue rnd testedLM thetaspace rhospace trial tcounter

    disp('Selecting Next Stimulus Using GLMP Paradigm...')
    %tcounter = tcounter + 1;
    
    % General variables
    nTestFr = 3;
    sisterAngs = [-pi -pi/2 0 pi/2]';
    sisterLM = nan(4,2);
    
    % Incoming variables
    in.theta = theta;
    in.sigX = sigX;
    in.sigY = sigY;
    in.DR = DR;
    in.nstd = nstd;
    in.nexp = nexp;
   
    % Plateaus
    var.Scc = [0]'; %Scc amounts
    var.theta = [0];  %Additive
    var.sigX = [1]'; %Scalar
    var.sigY = [1]'; %Scalar
    var.driftRate = [1]'; %Scalar
    var.nstd = [1]'; %Scalar
    var.nexp = [1]'; %Scalar

    % Take next in currentQueue, or make a new queue if currentQueue is empty
    if isempty(currentQueue)
        if isempty(rnd)
            rnd = 1;
        else
            rnd = rnd + 1;
        end

        % Construct Polar Grid
        if rnd == 1
            thetaspace = pi/4;
            rhospace = .5;
        elseif rnd>1 && mod(rnd,2)~=0
            thetaspace = thetaspace * .5;
        elseif rnd>1 && mod(rnd,2)==0
            rhospace = rhospace * .5;
        end
    
        theta = 0:thetaspace:pi/2;
        theta(theta==pi/2) = [];
        rho = 0:rhospace:1;
        rho(rho==0) = [];
        
        PolRhoIdx = fullfact([numel(theta) numel(rho)]);
        
        polTheta = theta(PolRhoIdx(:,1))';
        polRho = rho(PolRhoIdx(:,2))';
        polLM = [polTheta polRho];
        untestedLM = polLM(~ismember(polLM,testedLM,'rows'),:);
        
        % Initiate Plat Organization
        PlatIdx = fullfact([numel(var.Scc) numel(var.theta) numel(var.sigX)...
            numel(var.sigY) numel(var.driftRate) numel(var.nstd) numel(var.nexp)]);
        out(:,1) = var.Scc(PlatIdx(:,1));
        out(:,2) = var.theta(PlatIdx(:,2)) + in.theta; %Additive
        out(:,3) = var.sigX(PlatIdx(:,3)) .* in.sigX; %Scalar
        out(:,4) = var.sigY(PlatIdx(:,4)) .* in.sigY; %Scalar
        out(:,5) = var.driftRate(PlatIdx(:,5)) .* in.DR; %Scalar
        out(:,6) = var.nstd(PlatIdx(:,6)) .* in.nstd; %Scalar
        out(:,7) = var.nexp(PlatIdx(:,7)) .* in.nexp; %Scalar
        
        
        for n = 1:size(untestedLM,1)
            
            choice = randi(size(untestedLM,1));
            choiceTheta = untestedLM(choice,1);
            choiceRho = untestedLM(choice,2);
            
            sisterLM(:,1) = deal(choiceTheta) + sisterAngs;
            sisterLM(:,2) = deal(choiceRho);

            StimFamIdx = fullfact([size(sisterLM,1) size(out,1)]);
            stimFam = [sisterLM(StimFamIdx(:,1),:) out(StimFamIdx(:,2),:)];
            
            %Randomize Order of StimFam for each presentation
            tempParIdx = [];
            for q = 1:nTestFr
                order = randperm(size(stimFam,1))';
                tempParIdx = cat(1,tempParIdx,order);
                currentQueue = cat(1,currentQueue,stimFam(order,:));
            end

            % Delete last choice so it is not picked again
            untestedLM(choice,:) = [];
            
        end
        testedLM = cat(1, testedLM, unique(currentQueue(:,1:2),'rows'));
    end
    
    % Transform from polar to cartesian coordinates
    polTheta = currentQueue(1,1);
    polRho = currentQueue(1,2);
    [LmM,LpM] = pol2cart(polTheta,polRho);
    LmM = rndofferr(LmM,3);
    LpM = rndofferr(LpM,3);
    currentQueue(1,1) = LmM;
    currentQueue(1,2) = LpM;
        
    % Rotate and Scale Normalized Grid for Monitor
    transMat = [.08 .08; -.08 .08]; %This matrix transforms from a normalized cone opponent space to an elongated cone opponent space with vertices [(.45,.55),(.55,.45),(-.45,-.55),(-.55,-.45)] calculated by GH and JPW
    currentQueue(1,1:2) = transMat * [LmM LpM]';
    
    % Fill in trial fields
    trial.LmM_norm(tcounter) = LmM;
    trial.LpM_norm(tcounter) = LpM;
    %trial.stim_norm(tcounter,:) = [trial.Lcc_norm(tcounter) trial.Mcc_norm(tcounter)];
    trial.polTheta_norm(tcounter) = polTheta;
    trial.polRho_norm(tcounter) = polRho;
    trial.Lcc_orig(tcounter) = currentQueue(1,1);
    trial.Mcc_orig(tcounter) = currentQueue(1,2);
    trial.stim_orig(tcounter,:) = [trial.Lcc_orig(tcounter) trial.Mcc_orig(tcounter)];
    [polTheta_orig,polRho_orig] = cart2pol(currentQueue(1,1),currentQueue(1,2));
    trial.polTheta_orig(tcounter) = rndofferr(polTheta_orig,3);
    trial.polRho_orig(tcounter) = rndofferr(polRho_orig,3);
    trial.Scc_orig(tcounter) = currentQueue(1,3);
    %trial.LmM_norm(tcounter) = trial.stim_norm(tcounter,1);
    %trial.LpM_norm(tcounter) = trial.stim_norm(tcounter,2);
    trial.theta(tcounter) = currentQueue(1,4);
    trial.sigX(tcounter) = currentQueue(1,5);
    trial.sigY(tcounter) = currentQueue(1,6);
    trial.driftRate(tcounter) = currentQueue(1,7);
    trial.nstd(tcounter) = currentQueue(1,8);
    trial.nexp(tcounter) = currentQueue(1,9);
    trial.paradigm(tcounter) = cellstr('GLMP');


    fprintf('Round = %i \n Stimuli already tested = %i\n Stimuli in current queue = %i \n',rnd,tcounter-1,size(currentQueue,1)-1)
    fprintf('Current Stimulus Parameters:\n Lcc = %g\n Mcc = %g\n Scc = %g\n Theta = %g\n Sigma X = %g\n Sigma Y = %g\n Drift Rate = %g\n nStd = %g\n nExpanse = %g\n\n', currentQueue(1,:))
   
    % Return first in queue, then delete it
    stim = currentQueue(1,:);
    currentQueue(1,:) = [];
    
end


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
global tcounter trial datastruct par

    if isempty(datastruct)
    
        figure(1);
        subplot(2,9,[4 5]); cla; hold on; grid on;
        plot(trial.xInit(1:2:tcounter, 1), trial.xInit(1:2:tcounter,2), 'mo', 'LineWidth',2,...
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
        figure(1)
        subplot(2,9,[13 14]); cla; hold on; grid on;
        contour(xx, yy, reshape(datastruct.lambFinal, length(xx), []), 5); axis image; axis xy; title('estimated firing map');
        set(gca, 'xlim', [-1 1], 'ylim', [-1 1]);
        xlabel('Normalized L-M')
        ylabel('Normalized L+M')
        title('Surface Contour')
        
    end
    
    % 3-D Stem Plot
    figure(1)
    subplot(2,9,[7 8 9 16 17 18]); cla; hold on;
    maxresp = max(par.ALMP.meanfr);
    if maxresp == 0
        maxresp = 5;
    end
    grid on; axis([-1 1 -1 1 0 maxresp])
    axis([-1 1 -1 1])
    xlabel('L-M')
    ylabel('L+M')
    zlabel('Firing Rate (sp/s)')
    title('3D Stem Plot')
    stem3(par.ALMP.LmM_norm,par.ALMP.LpM_norm,par.ALMP.meanfr,'fill','MarkerFaceColor','c');
    if tcounter == 1
        set(gca,'view',[-75 38]);
    else
        a = get(gca,'View');
        set(gca,'View',a);
    end

end


function plotJPWGUI()
global tcounter par

    figure(2)
    subplot(2,5,[4 5 9 10]); cla; hold on;
    maxresp = max(par.GLMP.meanfr);
    if maxresp == 0
        maxresp = 5;
    end
    grid on; axis([-1 1 -1 1 0 maxresp])
    axis([-1 1 -1 1])
    stem3(par.GLMP.LmM_norm,par.GLMP.LpM_norm,par.GLMP.meanfr,'fill','MarkerFaceColor','c');
    xlabel('Normalized L-M')
    ylabel('Normalized L+M')
    zlabel('Firing Rate (sp/s)')
    title('3D Stem Plot')
    if tcounter == 2
        set(gca,'view',[-75 38]);
    else
        a = get(gca,'View');
        set(gca,'View',a);
    end

end

function resetVars()
global tcounter rnd testedLM thetaspace rhospace currentQueue

tcounter = 0;
rnd = [];
testedLM = [];
thetaspace = [];
rhospace = [];
currentQueue = [];

end

function initTrialStruct()
global trial

    %trial.Lcc_norm = nan(2,1);
    %trial.Mcc_norm = nan(2,1);
    %trial.stim_norm = nan(2,2);
    trial.polTheta_norm = nan(2,1);
    trial.polRho_norm = nan(2,1);
    trial.Lcc_orig = nan(2,1);
    trial.Mcc_orig = nan(2,1);
    trial.stim_orig = nan(2,2);
    trial.polTheta_orig = nan(2,1);
    trial.polRho_orig = nan(2,1);
    trial.Scc_orig = nan(2,1);
    trial.LmM_norm = nan(2,1);
    trial.LpM_norm = nan(2,1);
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
    trial.paradigm = cell(2,1);

end


function organizeParStruct()
global trial par tcounter

    str = {'ALMP' 'GLMP'};
    if mod(tcounter,2) == 1
        i = 1;
    else
        i = 2;
    end

    L = strcmp(trial.paradigm,str{i});

    % By Paradigm Properties (aka by unique stimulus)
    [uniqueCond,m,idx] = unique([trial.Lcc_orig(L) trial.Mcc_orig(L) trial.Scc_orig(L) trial.theta(L) trial.sigX(L) trial.sigY(L) trial.driftRate(L) trial.nstd(L) trial.nexp(L)],'rows');

    if any(isnan(uniqueCond(:,1)))
        nanIdx = find(isnan(uniqueCond(:,1)));
        uniqueCond(nanIdx,:) = [];
        m(nanIdx) = [];
        idx(nanIdx) = [];
    end

    % Collapse repeated stimuli
    par.(str{i}).Lcc_orig = uniqueCond(:,1);
    par.(str{i}).Mcc_orig = uniqueCond(:,2);
    par.(str{i}).Scc_orig = uniqueCond(:,3);
    par.(str{i}).stim_orig = [par.(str{i}).Lcc_orig par.(str{i}).Mcc_orig];
    temp = trial.polTheta_orig(L);
    par.(str{i}).polTheta_orig = temp(m);
    temp = trial.polRho_orig(L);
    par.(str{i}).polRho_orig = temp(m);
    %temp = trial.Lcc_norm(L);
    %par.(str{i}).Lcc_norm = temp(m);
    %temp = trial.Mcc_norm(L);
    %par.(str{i}).Mcc_norm = temp(m);
    temp = trial.polTheta_norm(L);
    par.(str{i}).polTheta_norm = temp(m);
    temp = trial.polRho_norm(L);
    par.(str{i}).polRho_norm = temp(m);
    temp = trial.LmM_norm(L);
    par.(str{i}).LmM_norm = temp(m);
    %par.(str{i}).LmM_norm = par.(str{i}).LmM_norm;
    temp = trial.LpM_norm(L);
    par.(str{i}).LpM_norm = temp(m);
    %par.(str{i}).LpM_norm = par.(str{i}).LpM_norm;
    %par.(str{i}).stim_norm = [par.(str{i}).LmM_norm par.(str{i}).LpM_norm];
    par.(str{i}).theta = uniqueCond(:,4);
    par.(str{i}).sigX = uniqueCond(:,5);
    par.(str{i}).sigY = uniqueCond(:,6);
    par.(str{i}).driftRate = uniqueCond(:,7);
    par.(str{i}).nstd = uniqueCond(:,8);
    par.(str{i}).nexp = uniqueCond(:,9);

    % Preallocate Space
    par.(str{i}).frs = cell(size(uniqueCond,1),1);
    par.(str{i}).meanfr = nan(size(uniqueCond,1),1);
    par.(str{i}).varfrs = nan(size(uniqueCond,1),1);
    par.(str{i}).baselinefrs = cell(size(uniqueCond,1),1);
    par.(str{i}).meanblfr = nan(size(uniqueCond,1),1);
    par.(str{i}).normtspikes = cell(size(uniqueCond,1),1);
    par.(str{i}).catntspikes = cell(size(uniqueCond,1),1);

    for n = 1:size(par.(str{i}).stim_orig,1)
        LL = idx == n;
        temp = trial.fr(L);
        par.(str{i}).frs{n} = temp(LL)';
        par.(str{i}).meanfr(n) = mean(par.(str{i}).frs{n});
        par.(str{i}).varfrs(n) = var(par.(str{i}).frs{n});
        temp = trial.baselinefr(L);
        par.(str{i}).baselinefrs{n} = temp(LL)';
        par.(str{i}).meanblfr(n) = mean(par.(str{i}).baselinefrs{n});
        temp = trial.normtspikes(L);
        par.(str{i}).normtspikes{n} = temp(LL);
        temp = trial.normtspikes(L);
        par.(str{i}).catntspikes{n} = cat(1,temp{LL});
    end
    
end


function plat = organizePlatStruct()
global trial par tcounter

    str = {'ALMP' 'GLMP'};

    for i = 1:2
        
        % By paradigm
        [uniqueCond,m,idx] = unique([par.(str{i}).Scc_orig par.(str{i}).theta par.(str{i}).sigX...
            par.(str{i}).sigY par.(str{i}).driftRate par.(str{i}).nstd par.(str{i}).nexp],'rows');
        nplat = size(uniqueCond,1);
        
        for p = 1:nplat
            L = idx == p;
            plat(p).par.(str{i}).Lcc_orig = par.(str{i}).Lcc_orig(L,:);
            plat(p).par.(str{i}).Mcc_orig = par.(str{i}).Mcc_orig(L,:);
            plat(p).par.(str{i}).Scc_orig = par.(str{i}).Scc_orig(L,:);
            plat(p).par.(str{i}).stim_orig = par.(str{i}).stim_orig(L,:);
            plat(p).par.(str{i}).polTheta_orig = par.(str{i}).polTheta_orig(L);
            plat(p).par.(str{i}).polRho_orig = par.(str{i}).polRho_orig(L);
            %plat(p).par.(str{i}).Lcc_norm = par.(str{i}).Lcc_norm(L);
            %plat(p).par.(str{i}).Mcc_norm = par.(str{i}).Mcc_norm(L);
            %plat(p).par.(str{i}).stim_norm = par.(str{i}).stim_norm(L,:);
            plat(p).par.(str{i}).polTheta_norm = par.(str{i}).polTheta_norm(L);
            plat(p).par.(str{i}).polRho_norm = par.(str{i}).polRho_norm(L);
            plat(p).par.(str{i}).LmM_norm = par.(str{i}).LmM_norm(L);
            plat(p).par.(str{i}).LpM_norm = par.(str{i}).LpM_norm(L);
            plat(p).par.(str{i}).theta = par.(str{i}).theta(L,:);
            plat(p).par.(str{i}).sigX = par.(str{i}).sigX(L,:);
            plat(p).par.(str{i}).sigY = par.(str{i}).sigY(L,:);
            plat(p).par.(str{i}).driftRate = par.(str{i}).driftRate(L,:);
            plat(p).par.(str{i}).nstd = par.(str{i}).nstd(L,:);
            plat(p).par.(str{i}).nexp = par.(str{i}).nexp(L,:);
            plat(p).par.(str{i}).frs = {par.(str{i}).frs{L}}';
            plat(p).par.(str{i}).meanfr = par.(str{i}).meanfr(L);
            plat(p).par.(str{i}).varfrs = par.(str{i}).varfrs(L);
            plat(p).par.(str{i}).baselinefrs = par.(str{i}).baselinefrs(L);
            plat(p).par.(str{i}).meanblfr = par.(str{i}).meanblfr(L);
            plat(p).par.(str{i}).normtspikes = par.(str{i}).normtspikes(L);
            plat(p).par.(str{i}).catntspikes = par.(str{i}).catntspikes(L);
        end
        
        
        % By trial
        LL = strcmp(trial.paradigm,str{i});
        [uniqueCond,m,idx] = unique([trial.Scc_orig(LL) trial.theta(LL) trial.sigX(LL)...
            trial.sigY(LL) trial.driftRate(LL) trial.nstd(LL) trial.nexp(LL)],'rows');
        
        for p = 1:nplat
            L = idx == p;
            temp = trial.Lcc_orig(LL);
            plat(p).trial.(str{i}).Lcc_orig = temp(L);
            temp = trial.Mcc_orig(LL);
            plat(p).trial.(str{i}).Mcc_orig = temp(L);
            temp = trial.Scc_orig(LL);
            plat(p).trial.(str{i}).Scc_orig = temp(L);
            plat(p).trial.(str{i}).stim_orig = [plat(p).trial.(str{i}).Lcc_orig plat(p).trial.(str{i}).Mcc_orig];
            temp = trial.polTheta_orig(LL);
            plat(p).trial.(str{i}).polTheta_orig = temp(L);
            temp = trial.polRho_orig(LL);
            plat(p).trial.(str{i}).polRho_orig = temp(L);
            %temp = trial.Lcc_norm(LL);
            %plat(p).trial.(str{i}).Lcc_norm = temp(L);
            %temp = trial.Mcc_norm(LL);
            %plat(p).trial.(str{i}).Mcc_norm = temp(L);
            %plat(p).trial.(str{i}).nstim = [plat(p).trial.(str{i}).Lcc_norm plat(p).trial.(str{i}).Mcc_norm];
            temp = trial.polTheta_norm(LL);
            plat(p).trial.(str{i}).polTheta_norm = temp(L);
            temp = trial.polRho_norm(LL);
            plat(p).trial.(str{i}).polRho_norm = temp(L);
            temp = trial.LmM_norm(LL);
            plat(p).trial.(str{i}).LmM_norm = temp(L);
            temp = trial.LpM_norm(LL);
            plat(p).trial.(str{i}).LpM_norm = temp(L);
            temp = trial.theta(LL);
            plat(p).trial.(str{i}).theta = temp(L);
            temp = trial.sigX(LL);
            plat(p).trial.(str{i}).sigX = temp(L);
            temp = trial.sigY(LL);
            plat(p).trial.(str{i}).sigY = temp(L);
            temp = trial.driftRate(LL);
            plat(p).trial.(str{i}).driftRate = temp(L);
            temp = trial.nstd(LL);
            plat(p).trial.(str{i}).nstd = temp(L);
            temp = trial.nexp(LL);
            plat(p).trial.(str{i}).nexp = temp(L);
            temp = trial.fr(LL);
            plat(p).trial.(str{i}).fr = temp(L);
            temp = trial.baselinefr(LL);
            plat(p).trial.(str{i}).baselinefr = temp(L);
            temp = trial.normtspikes(LL);
            plat(p).trial.(str{i}).normtspikes = temp(L);
            temp = trial.nspikes(LL);
            plat(p).trial.(str{i}).nspikes = temp(L);
            temp = trial.tspikes(LL);
            plat(p).trial.(str{i}).tspikes = temp(L);
        end
    end

end


function initDataStruct()
global datastruct

    datastruct = [];

end

function initParStruct()
global par

par = [];

end




