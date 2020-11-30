function GridLMPlaneOnline()

% 2/13/13  Added parameters to specify 'flash' and 'drift'.    JPW
% 3/22/13  Removed 'flash' and 'drift' parameters.             JPW

global udpCom currentQueue rnd testedLM gridspace trial tcounter

disp('In GridLMPlaneOnline')

currentQueue = [];
rnd = [];
testedLM = [];
gridspace = [];
udpCom.sock = [];
udpCom.port = 6665;
udpCom.rexip = '192.168.1.120';
udpCom.plexip = '192.168.1.121';
udpCom.mac = '192.168.1.122';

C = GridLMPlaneCodes();% defined in a separate m file

p = InitPStruct(0, C.HDRCOMPLETECD);
s = InitPlex();
initGUI();
initTrialStruct();

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
        p.lastprocessed_t = p.times(find(p.events == C.STIMONCD,1));
    end
    if (any(p.spikes{1}) && ~isempty(fix_t) && ~abortflag)
        currentspiketally = [currentspiketally; p.spikes{1}-fix_t]; %Spike times are relative to fix_t
        PlotRaster(currentspiketally);
        p.spikes{1} = []; 
    end
    if (any(p.events == C.ABORTCD))
        abortflag = 1;
        currentspiketally = [];
        %ClearRaster;
        p.lastprocessed_t = p.times(find(p.events == C.ABORTCD,1,'last'));
    end
    if (any(p.events == C.FPOFFCD))
        stimoff_t = p.times(find(p.events == C.FPOFFCD,1));

        trial.fpacq(tcounter) = fix_t;
        trial.stimon(tcounter) = stimon_t;
        trial.stimoff(tcounter) = p.times(find(p.events == C.FPOFFCD,1));    
        trial.duration(tcounter) = stimoff_t-stimon_t;
        trial.tspikes{tcounter} = currentspiketally + stimon_t;
        trial.normtspikes{tcounter} = currentspiketally;
        
        % Spike Historgram
        PlotSpikeHist
        
        % Clear out fields
        currentspiketally = [];
        stimon_t = [];
        fix_t = [];        
        p.lastprocessed_t = stimoff_t;
        
        % Organize and Plot data
        organizeParStruct();
        plat = organizePlatStruct();
        plotGUI(plat);
        
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
    subplot(2,5,[1 2]); cla; hold on;
    axis([-500 xmax -.5 .5])
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
    subplot(2,5,[6 7]); cla; hold on;
    xlim([-500 max(trial.duration)])
    bar(bins,PSTH);
    plot(bins,gauss(bins,fittedparams),'b-.','LineWidth',3);
    plot(offset,[fittedparams(1) fittedparams(1)],'k-','linewidth',3);
    xlabel('Normalized time (ms)')
    ylabel('Spike Count')
    title('Spike Histogram')

    % Now counting up spikes in a window
    for i = 1:tcounter
        trial.nspikes(i) = sum(trial.normtspikes{i} > offset(1) & trial.normtspikes{i} < offset(2));
        dt = offset(2) - offset(1);
        trial.fr(i) = trial.nspikes(i)./dt;
    end

end


function allDone = dealWithMsgs(socket)
    global udpCom
    allDone = false;
    msgSize = pnet(socket, 'readpacket', 250, 'noblock');
    if ~msgSize, return; end

    message = pnet(socket, 'read', msgSize, 'char');
    if ~isempty(message)
        [message,allDone] = parseMsg(message,socket);
        evalMsg(message);
    else
        allDone = true;
    end
end


% here we implement the logic/function calls associated with the requested
% variables within the message
function [message,doneFlag] = parseMsg(message, socket)
    global udpCom
    %disp('In parseMsg')
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
        end
    elseif strncmp(message, 'return', 6) % not sure what catastrophe this prevents
        stk = dbstack();
        if ~strcmp(stk(end).name, mfilename)
            doneFlag = true;
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
    global currentQueue rnd testedLM thetaspace rhospace trial tcounter

    disp('Selecting Next Stimulus')
    tcounter = tcounter + 1;
    
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
        PlatIdx = fullfact([numel(var.Scc) numel(var.theta)...
            numel(var.sigX) numel(var.sigY) numel(var.driftRate)...
            numel(var.nstd) numel(var.nexp) numel(var.presType)]);
        out(:,1) = var.Scc(PlatIdx(:,1));
        out(:,2) = var.theta(PlatIdx(:,2)) + in.theta; %Additive
        out(:,3) = var.sigX(PlatIdx(:,3)) .* in.sigX; %Scalar
        out(:,4) = var.sigY(PlatIdx(:,4)) .* in.sigY; %Scalar
        out(:,5) = var.driftRate(PlatIdx(:,5)) .* in.DR; %Scalar
        out(:,6) = var.nstd(PlatIdx(:,6)) .* in.nstd; %Scalar
        out(:,7) = var.nexp(PlatIdx(:,7)) .* in.nexp; %Scalar
        out(:,8) = var.presType(PlatIdx(:,8));
        
        
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
        testedLM = cat(1,testedLM,unique(currentQueue(:,1:2),'rows'));
    end
    
    % Transform from polar to cartesian coordinates
    polTheta = currentQueue(1,1);
    polRho = currentQueue(1,2);
    [Lcc,Mcc] = pol2cart(polTheta,polRho);
    Lcc = rndofferr(Lcc,3);
    Mcc = rndofferr(Mcc,3);
    currentQueue(1,1) = Lcc;
    currentQueue(1,2) = Mcc;
        
    % Rotate and Scale Normalized Grid for Monitor
    %transMat = [.08 .08; -.08 .08]; %This matrix transforms from a normalized cone opponent space to an elongated cone opponent space with vertices [(.45,.55),(.55,.45),(-.45,-.55),(-.55,-.45)] calculated by GH and JPW
    transMat = [.08 0; 0 .08];
    currentQueue(1,1:2) = transMat * [Lcc Mcc]';
    
    % Fill in trial fields
    trial.Lcc_norm(tcounter) = Lcc;
    trial.Mcc_norm(tcounter) = Mcc;
    trial.nstim(tcounter,:) = [trial.Lcc_norm(tcounter) trial.Mcc_norm(tcounter)];
    trial.npolTheta(tcounter) = polTheta;
    trial.npolRho(tcounter) = polRho;
    trial.Lcc_orig(tcounter) = currentQueue(1,1);
    trial.Mcc_orig(tcounter) = currentQueue(1,2);
    trial.ostim(tcounter,:) = [trial.Lcc_orig(tcounter) trial.Mcc_orig(tcounter)];
    [opolTheta,opolRho] = cart2pol(currentQueue(1,1),currentQueue(1,2));
    trial.opolTheta(tcounter) = rndofferr(opolTheta,3);
    trial.opolRho(tcounter) = rndofferr(opolRho,3);
    trial.Scc(tcounter) = currentQueue(1,3);
    trial.LmM(tcounter) = trial.nstim(tcounter,1);
    trial.LpM(tcounter) = trial.nstim(tcounter,2);
    trial.theta(tcounter) = currentQueue(1,4);
    trial.sigX(tcounter) = currentQueue(1,5);
    trial.sigY(tcounter) = currentQueue(1,6);
    trial.driftRate(tcounter) = currentQueue(1,7);
    trial.nstd(tcounter) = currentQueue(1,8);
    trial.nexp(tcounter) = currentQueue(1,9);

    fprintf('Round = %i \n Stimuli already tested = %i\n Stimuli in current queue = %i \n',rnd,tcounter,size(currentQueue,1))
    fprintf('Current Stimulus Parameters:\n Lcc = %.15f\n Mcc = %.15f\n Scc = %g\n Theta = %g\n Sigma X = %g\n Sigma Y = %g\n Drift Rate = %g\n nStd = %g\n nExpanse = %g\n\n', currentQueue(1,:))
   
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


function plotGUI(plat)
global tcounter

nplat = numel(plat);

for p = 1:nplat
    
    % Heat Map
%     figure(10+p);
%     subplot(2,1,1); cla; hold on;
%     minresp = min(plat(p).par.meanfr);
%     maxresp = max(plat(p).par.meanfr);
%     cmap = colormap(jet(64));
%     %polstim = [plat(p).par.npolTheta plat(p).par.npolRho];
%     %[sortedplat,sortorder] = sort(plat(p).par.polstim(1:numel(plat(p).par.meanfr),2),1,'descend');
%     %[sortedplat,sortorder] = sort(polstim,1,'descend');
%     [sortedplat,sortorder] = sort(plat(p).par.npolRho,1,'descend');
%     
%     for i = 1:size(plat(p).par.meanfr,1)
%         n = sortorder(i);
%         if ~isempty(plat(p).par.meanfr(n))
%             cmapidx = ceil((plat(p).par.meanfr(n)-minresp+1)/(maxresp-minresp+1)*63)+1;
%             h = polar(plat(p).par.npolTheta(n),plat(p).par.npolRho(n),'ko'); hold on;
%             set(h,'MarkerFaceColor',cmap(cmapidx,:))
%         end
%     end
%     xlabel('Normalized L-M')
%     ylabel('Normalized L+M')
%     title('Heat Map')
%     grid on; axis tight
    
    
    % 3-D Stem Plot
    figure(1)
    subplot(2,5,[4 5 9 10]); cla; hold on;
    maxresp = max(plat(p).par.meanfr);
    if maxresp == 0
        maxresp = 5;
    end
    grid on; axis([-1 1 -1 1 0 maxresp])
    axis([-1 1 -1 1])
    stem3(plat(p).par.LmM,plat(p).par.LpM,plat(p).par.meanfr,'fill','MarkerFaceColor','c');
    xlabel('Normalized L-M')
    ylabel('Normalized L+M')
    zlabel('Firing Rate (sp/s)')
    title('3D Stem Plot')
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
    trial.npolTheta = nan(2,1);
    trial.npolRho = nan(2,1);
    trial.Lcc_orig = nan(2,1);
    trial.Mcc_orig = nan(2,1);
    trial.opolTheta = nan(2,1);
    trial.opolRho = nan(2,1);
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

    tcounter = 0;
end


function organizeParStruct()
global trial par

    % By Paradigm Properties (aka by unique stimulus)
    [uniqueCond,m,idx] = unique([trial.Lcc_orig trial.Mcc_orig trial.Scc...
        trial.theta trial.sigX trial.sigY trial.driftRate trial.nstd trial.nexp],'rows');
    
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
    par.opolTheta = trial.opolTheta(m);
    par.opolRho = trial.opolRho(m);
    par.Lcc_norm = trial.Lcc_norm(m);
    par.Mcc_norm = trial.Mcc_norm(m);
    par.nstim = [par.Lcc_norm par.Mcc_norm];
    par.npolTheta = trial.npolTheta(m);
    par.npolRho = trial.npolRho(m);
    par.LmM = par.nstim(:,1);
    par.LpM = par.nstim(:,2);

    % Correcting thetas (roundoff error)
    divs = -pi-pi/64:pi/64:pi+pi/64;
    rads = divs(2:2:end);
    edges = divs(1:2:end);
    [n,binIdx] = histc(par.npolTheta,edges);
    par.npolTheta = rads(binIdx)';
    if any(trial.npolTheta==-pi)
        trial.npolTheta(trial.npolTheta==-pi)=pi;
    end

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
    par.baselinefrs = cell(size(uniqueCond,1),1);
    par.meanblfr = nan(size(uniqueCond,1),1);
    par.normtspikes = cell(size(uniqueCond,1),1);
    par.catntspikes = cell(size(uniqueCond,1),1);

    for n = 1:size(par.nstim,1)
        L = idx == n;
        par.frs{n} = trial.fr(L)';
        par.meanfr(n) = mean(par.frs{n});
        par.varfrs(n) = var(par.frs{n});
        par.baselinefrs{n} = trial.baselinefr(L);
        par.meanblfr(n) = mean(par.baselinefrs{n});
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
        plat(p).par.opolTheta = par.opolTheta(L);
        plat(p).par.opolRho = par.opolRho(L);
        plat(p).par.Lcc_norm = par.Lcc_norm(L);
        plat(p).par.Mcc_norm = par.Mcc_norm(L);
        plat(p).par.nstim = par.nstim(L,:);
        plat(p).par.npolTheta = par.npolTheta(L);
        plat(p).par.npolRho = par.npolRho(L);
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
        plat(p).par.baselinefrs = par.baselinefrs(L);
        plat(p).par.meanblfr = par.meanblfr(L);
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
        plat(p).trial.opolTheta = trial.opolTheta(L);
        plat(p).trial.opolRho = trial.opolRho(L);
        plat(p).trial.Lcc_norm = trial.Lcc_norm(L);
        plat(p).trial.Mcc_norm = trial.Mcc_norm(L);
        plat(p).trial.nstim = trial.nstim(L,:);
        plat(p).trial.npolTheta = trial.npolTheta(L);
        plat(p).trial.npolRho = trial.npolRho(L);
        plat(p).trial.LmM = trial.LmM(L);
        plat(p).trial.LpM = trial.LpM(L);
        plat(p).trial.theta = trial.theta(L);
        plat(p).trial.sigX = trial.sigX(L);
        plat(p).trial.sigY = trial.sigY(L);
        plat(p).trial.driftRate = trial.driftRate(L);
        plat(p).trial.nstd = trial.nstd(L);
        plat(p).trial.nexp = trial.nexp(L);
        plat(p).trial.fr = trial.fr(L);
        plat(p).trial.baselinefr = trial.baselinefr(L);
        plat(p).trial.normtspikes = trial.normtspikes(L);
        plat(p).trial.nspikes = trial.nspikes(L);
        plat(p).trial.tspikes = trial.tspikes(L);
    end
end






