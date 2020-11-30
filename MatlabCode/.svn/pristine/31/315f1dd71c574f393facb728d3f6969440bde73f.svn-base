function GLMSDetectionOnline()

% 2/10/13  New paradigm created.    JPW
% 3/2/13 Started adding GUI stuff  GDLH
% 7/14  Detection paradigm created JPW

global udpCom gl GLMP

disp('Starting GLMS Detection...')

%User-Defined Variables
%gl.datafile = 'N060613003.nex'; % +L-M Responsive (#29) (did I gather any detection data?)
%gl.datafile = 'N061513002.nex'; % -L+M Responsive (#34) *
%gl.datafile = 'N062713004.nex'; % RG Chromatic (#36) *
%gl.datafile = 'N112713002.nex'; % +L-M Responsive *
%gl.datafile = 'N102413002.nex'; % +L-M Chromatic (#43) *
%gl.datafile = 'N110514001.nex'; % +/- Luminance *
%gl.datafile = 'N011215002.nex'; % +L+M Horseshoe cell *
%gl.datafile = 'N112514001.nex'; % -L+M Horseshoe Cell *
gl.datafile = 'N040215002.nex';

gl.whichsub = 1;
gl.nPres = 1; % Repeats entire queue in new random order once all stimuli are tested nPres times.

% Non-User-Defined Variables
gl.library = 'C:\Documents and Settings\JPatrickWeller\My Documents\Dropbox\Patrick\GLMS Data\';
gl.CurrentQueue = [];
gl.gridX = [];
gl.gridY = [];
gl.Lans = 0;
gl.Rans = 0;
gl.thetaspace = pi/4;
gl.rhospace = .5;
abortflag = 0;
gl.alreadyshown = [];
udpCom.sock = [];
udpCom.port = 6665;
udpCom.rexip = '192.168.1.120';

C = GLMSDetectionCodes();% defined in a separate m file
p = InitPStruct(0, C.HDRCOMPLETECD);
s = InitPlex();
initGUI();
ConstructQueue();
PlotPsychPred();
SetUpPsych();
InitStatsStruct();

codestruct = {...
        C.ABORTCD, @abortfn;...
        C.EOTCD, @eotfn;...
        C.FPACQCD, @fpacqfn};
        
% Handshaking
pnet(udpCom.sock, 'write', 'plexready>> >>');
pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
    
% The main loop
socketOpen = 0;
while ~socketOpen
    [udpCom.sock, socketOpen] = pnetStart(udpCom.port);
end
plxServer = InitPlex();

bounceOut = false;
while ~bounceOut
    if CheckForESCKey() || dealWithMsgs()
        bounceOut = true;
    end
    
    [n, eventList] = PL_GetTS(s);
    if n
        p = ProcessEventList(p, eventList);
    end
    
    if any(p.events == p.headercode)
        GetHeader(p);
        a = get(figure(1),'UserData');
        p = RemoveOldEvents(p, p.headercode);
    end
    
    % Dealing with each code, one by one
    for i = 1:size(codestruct,1)
        if (any(p.events == codestruct{i,1}))
            feval(codestruct{i,2})
        end
    end
    
    p = CleanUpEvents(p);
    drawnow;
end % big while loop


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Event Triggered Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function fpacqfn %FPACQ
        abortflag = 0;     
        p.lastprocessed_t = p.times(find(p.events == C.FPACQCD,1));
    end

    function abortfn % ABORT
        abortflag = 1;
        
        % Return stimulus to queue and randomize order
        newOrder = randperm(size(gl.CurrentQueue,1))';
        gl.CurrentQueue = gl.CurrentQueue(newOrder,:);
        p.lastprocessed_t = p.times(find(p.events == C.ABORTCD,1,'last'));
        
    end

    function eotfn % EOT
        if ~abortflag
                        
            % Pull out variables
            l = GetValsIfPossible(C.LCCCD, p, 'double');
            m = GetValsIfPossible(C.MCCCD, p, 'double');
            rfcorr = GetValsIfPossible(C.RFCORRCD, p, 'int');
            anscorr = GetValsIfPossible(C.CORRECTRESPCD, p, 'int');
                        
            % Pull out user variables
            a = get(figure(1),'UserData');
            conpanel = get(a.axes.conpanel,'UserData');
            
            % Update plots
            PlotPsych(l,m,rfcorr,anscorr)
            
            %Update L/R counters
            if sign(GLMP.rf_x) == 1
                str1 = 'Rans';
                str2 = 'Lans';
            else
                str1 = 'Lans';
                str2 = 'Rans';
            end
            if (rfcorr && ~anscorr)
                gl.(str2) = gl.(str2) + 1;
            elseif (~rfcorr && ~anscorr)
                gl.(str1) = gl.(str1) + 1;
            end      
            axes(conpanel.axes.LRcounter); cla; hold on;
            bar([gl.Lans gl.Rans'])
            
            % Update Queue and Queue Counter
            gl.CurrentQueue = gl.CurrentQueue(2:end,:);
            set(conpanel.uicontrols.stimcounter,'string',size(gl.CurrentQueue,1))
            
            % Save figure data
            set(a.axes.conpanel,'UserData',conpanel)
            set(figure(1),'UserData',a);
        end
        
        % Handshaking
        p.lastprocessed_t =  p.times(find(p.events == C.EOTCD,1));
        pnet(udpCom.sock, 'write', 'plexready>> >>');
        pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
    end

PL_Close(plxServer);

end % end of main function

function GetHeader(p)
global gl

    C = GridLMSubunitCodes;
    
    gammaTable = GetValsIfPossible([C.GAMMATABLECD C.GAMMATABLECD], p, 'double');
    fundamentals =  GetValsIfPossible([C.FUNDAMENTALSCD C.FUNDAMENTALSCD], p, 'double');
    monSpd =  GetValsIfPossible([C.MONSPDCD C.MONSPDCD], p, 'double');
    a = get(figure(1),'UserData');
    a.stats.msperframe = 1000/GetValsIfPossible(C.FRAMERATECD, p, 'double');
    a.stats.gammaTable = reshape(gammaTable, length(gammaTable)/3,3);
    a.stats.invgammaTable = InvertGamma(a.stats.gammaTable,1);
    a.stats.monSpd = reshape(monSpd, length(monSpd)/3,3);
    a.stats.fundamentals = reshape(fundamentals, length(fundamentals)/3,3);
    a.stats.bkgndrgb = GetValsIfPossible([C.BKGNDRGBCD C.BKGNDRGBCD], p, 'double');
    P_device = SplineSpd(linspace(380,780,size(a.stats.monSpd,1))',a.stats.monSpd,(380:5:780)');
    a.stats.M = a.stats.fundamentals'*P_device;
    a.stats.invM = inv(a.stats.M);
    a.stats.lmsbinaryrgbmat = nan*ones(4,3);    
    set(figure(1),'UserData',a);
    
    gl.gotHeader = 1;

end

function InitStatsStruct()
global gl

    a = get(figure(1),'UserData');
    a.stats.npixperstix = []; % Number of pixels per side of stixel  
    a.stats.nstixperside = []; % Number of stixels per side of stimulus
    a.stats.msperframe = [];  % ms per frame
    a.stats.gammaTable = []; % gamma tables
    a.stats.invgammaTable = []; % inverse gamma tables
    a.stats.monSpd = []; % monitor phosphor spectra
    a.stats.fundamentals = []; % cone fundamentals
    a.stats.M = [];     % guns to cones matrix
    a.stats.invM = [];  % cones to guns matrix
    set(figure(1),'UserData',a);
    
    gl.gotHeader = 0;
    
end

function SetUpPsych()
global gl

    a = get(figure(1),'UserData');
    glmppanel = get(a.axes.glmppanel,'UserData');
    
    allstim = unique(gl.CurrentQueue(:,1:2),'rows');
    
    % Set up RF's
    for n = 1:2
        str = ['RF' num2str(n)];
        glmppanel.(str).Lcc = allstim(:,1);
        glmppanel.(str).Mcc = allstim(:,2);
        glmppanel.(str).resps = cell(size(glmppanel.(str).Lcc));
        glmppanel.(str).meanresp = nan(size(glmppanel.(str).Lcc));
    end
    
    % Set up AllData structure
    glmppanel.AllData.Lcc = allstim(:,1);
    glmppanel.AllData.Mcc = allstim(:,2);
    glmppanel.AllData.resps = cell(size(glmppanel.(str).Lcc));
    glmppanel.AllData.meanresp = nan(size(glmppanel.(str).Lcc));
    
    set(a.axes.glmppanel,'UserData',glmppanel)
    set(gcf,'UserData',a)

end

function allDone = dealWithMsgs()
global udpCom GLMP gl
%GLMP included so that the RF is queriable 

    allDone = false;
    msgSize = pnet(udpCom.sock, 'readpacket', 250, 'noblock');
    if ~msgSize, return; end

    message = pnet(udpCom.sock, 'read', msgSize, 'char');
    if ~isempty(message)
        if strncmp(message, 'return', 6) % this detects if the script was called from the command line
            stk = dbstack();
            if ~strcmp(stk(end).name, mfilename)
                allDone = true;
            end
        end
    
        try
            eval(message);
        catch exception
            fprintf('Trouble with message: "%s"\n', message);
            disp(getReport(exception));
        end
    else
        allDone = true;
    end
    
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
        p.spikes{i} = spiketimevect(spiketimevect > p.lastprocessed_t);
    end
    p.processnowflag = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stimulus Selection Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [stim] = getStimParams()
global gl 

    if isempty(gl.CurrentQueue)
        
        % Make a new queue
        ConstructQueue();
        
    end
        
    % Take the top stimulus in the queue
    stim = gl.CurrentQueue(1,:);
    
    % Update Lcc/Mcc Val
    a = get(figure(1),'UserData');
    conpanel = get(a.axes.conpanel,'UserData');
    val = [num2str(rndofferr(stim(1),3)) '  ' num2str(rndofferr(stim(2),3))];
    set(conpanel.uicontrols.LMVal,'string',val);
    
    % Save Variables
    set(a.axes.conpanel,'UserData',conpanel);
    set(figure(1),'UserData',a);
    
end


function ConstructQueue()
global gl

    nx = 2^(floor(gl.nRnds/2));
    if mod(gl.nRnds,2) == 0
        thetas = shiftdim(0:(gl.thetaspace/sqrt(nx)):(2*pi-(gl.thetaspace/sqrt(nx))));
    else
        thetas = shiftdim(0:(gl.thetaspace/nx):(2*pi-(gl.thetaspace/nx)));
    end
    rhos = shiftdim(gl.rhospace/nx:gl.rhospace/nx:1);
    rfCorrect = shiftdim([0 1]);
    
    % Totally random order
    FFIdx = fullfact([numel(thetas) numel(rhos) numel(rfCorrect)]);
    stimParamList = repmat([thetas(FFIdx(:,1)) rhos(FFIdx(:,2)) rfCorrect(FFIdx(:,3))],gl.nPres,1);
    tempqueue = stimParamList(randperm(size(stimParamList,1)),:);
    
    % Transform from polar to cartesian coordinates
    [tempLcc,tempMcc] = pol2cart(tempqueue(:,1),tempqueue(:,2));
    
    % Scale Cone Contrast Units for Monitor
    scale = gl.lumcc*gl.colcc./sqrt((gl.colcc.*cos(tempqueue(:,1)-pi/4)).^2 ...
        +(gl.lumcc.*sin(tempqueue(:,1)-pi/4)).^2);
    Lcc = tempLcc .* scale;
    Mcc = tempMcc .* scale;
    
    % Add the other stimulus parameters
    Scc = zeros(size(Lcc));
    DVA = ones(size(Lcc)) * gl.DVAPerStixel;
    stimDur = ones(size(Lcc)) * gl.StimDuration;
    nstix = ones(size(Lcc)) * gl.nStixelGrid;
    RFCorrect = tempqueue(:,3);
    
    % Final queue order
    gl.CurrentQueue = [Lcc Mcc Scc DVA stimDur nstix RFCorrect];
    
    % Update dataset counter
    a = get(figure(1),'UserData');
    conpanel = get(a.axes.conpanel,'UserData');
    val = str2double(get(conpanel.uicontrols.ndatasets,'string'));
    set(conpanel.uicontrols.ndatasets,'string',val+1);
    set(a.axes.conpanel,'UserData',conpanel);
    set(figure(1),'UserData',a)
    
end


function gridx = getX()
global gl

    gridx = gl.gridX;
    
end

function gridy = getY()
global gl

    gridy = gl.gridY;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting Up GUI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function initGUI()
global gl GLMP DN

    % Unpack Nex File
    rawdata = nex2stro([char(gl.library) gl.datafile]);
    [GLMP,DN] = OrganizeRawGLMSData(rawdata);
    
    % Set Up Stimulus Parameters
    nstim = numel(GLMP.subunit{gl.whichsub}.uniqueLcc);
    gl.nRnds = 1;
    while nstim > 16
        nstim = nstim/2;
        gl.nRnds = gl.nRnds + 1;
    end
    lumL = GLMP.subunit{gl.whichsub}.uniquetheta == -pi/4;
    colL = GLMP.subunit{gl.whichsub}.uniquetheta == pi/4;
    gl.colcc = max(GLMP.subunit{gl.whichsub}.uniquerho(lumL));
    gl.lumcc = max(GLMP.subunit{gl.whichsub}.uniquerho(colL));
    gl.maxcc = max(GLMP.subunit{gl.whichsub}.uniquerho);
    gl.gridX = GLMP.subunit{gl.whichsub}.gridX(1,:);
    gl.gridY = GLMP.subunit{gl.whichsub}.gridY(1,:);
    gl.nStixelGrid = DN.NStixGrid(1);
    gl.StimDuration = GLMP.stimDur(1);
    gl.DVAPerStixel = DN.DVAPerStix(1);
    
    % Set up GUI Pannel
    close all;
    figure(1); clf;
    set(figure(1),'DefaultAxesUnits','normalized','position',[100 100 1000 700],'UserData',[])
    a = get(figure(1),'UserData');
    
    % Set up subpannels
    a.axes.conpanel = uipanel('Title','Control Panel','Position',[.01 .02 .25 .96],'BorderType','beveledin','BackgroundColor',[.9 .9 .9]);
    a.axes.glmppanel = uipanel('Title','GLMP','Position',[.275 .02 .715 .96],'BorderType','beveledin','BackgroundColor',[.9 .9 .9]);
    
    % Set Up Control Panel
    conpanel = get(a.axes.conpanel,'UserData');
    
    % L/M value of stimulus
    conpanel.uicontrols.LMVal = uicontrol('Parent',a.axes.conpanel,'Units','normalized',...
        'style','edit','string','~','FontSize',12,...
        'Position',[.2 .87 .6 .05],'BackgroundColor',[.8 .8 .8]);
    conpanel.labels.LMVal = uicontrol('Parent',a.axes.conpanel,'Units','normalized',...
        'style','text','string','Lcc/Mcc Value:','FontSize',14,...
        'HorizontalAlignment','center',...
        'Position',[.1 .92 .8 .05],'BackgroundColor',[.9 .9 .9]);
    
    % # of Stimuli in queue
    conpanel.uicontrols.stimcounter = uicontrol('Parent',a.axes.conpanel,'Units','normalized',...
        'style','edit','string','~','FontSize',12,...
        'Position',[.2 .72 .6 .05],'BackgroundColor',[.8 .8 .8]);
    conpanel.labels.stimcounter = uicontrol('Parent',a.axes.conpanel,'Units','normalized',...
        'style','text','string','# of Stimuli in Queue:','FontSize',14,...
        'HorizontalAlignment','center',...
        'Position',[.05 .77 .9 .05],'BackgroundColor',[.9 .9 .9]);
    
    % Dataset counter
    conpanel.uicontrols.ndatasets = uicontrol('Parent',a.axes.conpanel,'Units','normalized',...
        'style','edit','string',0,'FontSize',12,...
        'Position',[.2 .57 .6 .05],'BackgroundColor',[.8 .8 .8]);
    conpanel.labels.ndatasets = uicontrol('Parent',a.axes.conpanel,'Units','normalized',...
        'style','text','string','Dataset #:','FontSize',14,...
        'HorizontalAlignment','center',...
        'Position',[.1 .62 .8 .05],'BackgroundColor',[.9 .9 .9]);
    
    % L/R Counter
    labels{1} = 'Left';
    labels{2} = 'Right';
    labels = labels';
    conpanel.axes.LRcounter = axes('parent',a.axes.conpanel,'units','normalized',...
        'position',[.125 .05 .775 .38],'Xtick',[1 2],'xlim',[0 3],...
        'XTickLabel',labels,'ygrid','on');
    conpanel.labels.LRcounter = uicontrol('Parent',a.axes.conpanel,'Units','normalized',...
        'style','text','string','Incorrect Choices:','FontSize',14,...
        'HorizontalAlignment','center',...
        'Position',[.1 .455 .8 .05],'BackgroundColor',[.9 .9 .9]);

    % GLMP Panel
    % Neurometric
    glmppanel = get(a.axes.glmppanel,'UserData');
    glmppanel.neuro.axes = axes('Parent',a.axes.glmppanel,'Units','normalized',...
        'Position',[.1 .55 .35 .35],'CameraPosition',[-.9 -.9 5],...
        'XLim',[-gl.maxcc gl.maxcc],'YLim',[-gl.maxcc gl.maxcc],'ZLim',[0 1.1],...
        'XGrid','on','YGrid','on','ZGrid','on');
    title('Neurometric Data','FontSize',14,'FontWeight','bold'); hold on;
    xlabel('Lcc'); ylabel('Mcc'); zlabel('Predicted Percentage Correct')
    
    % Psychometric
    glmppanel.AllData.axes = axes('Parent',a.axes.glmppanel,'Units','normalized',...
        'Position',[.55 .55 .35 .35],'CameraPosition',[-.9 -.9 5],...
        'XLim',[-gl.maxcc gl.maxcc],'YLim',[-gl.maxcc gl.maxcc],'ZLim',[0 1.1],...
        'XGrid','on','YGrid','on','ZGrid','on');
    title('All Psychometric Data','FontSize',14,'FontWeight','bold'); hold on;
    xlabel('Lcc'); ylabel('Mcc'); zlabel('Percentage Correct')
    
    %RF1
    glmppanel.RF1.axes = axes('Parent',a.axes.glmppanel,'Units','normalized',...
        'Position',[.1 .1 .35 .35],'CameraPosition',[-.9 -.9 5],...
        'XLim',[-gl.maxcc gl.maxcc],'YLim',[-gl.maxcc gl.maxcc],'ZLim',[0 1.1],...
        'XGrid','on','YGrid','on','ZGrid','on');
    title('RF1 Data','FontSize',14,'FontWeight','bold'); hold on;
    xlabel('Lcc'); ylabel('Mcc'); zlabel('Percentage Correct')
    
    %RF2
    glmppanel.RF2.axes = axes('Parent',a.axes.glmppanel,'Units','normalized',...
        'Position',[.55 .1 .35 .35],'CameraPosition',[-.9 -.9 5],...
        'XLim',[-gl.maxcc gl.maxcc],'YLim',[-gl.maxcc gl.maxcc],'ZLim',[0 1.1],...
        'XGrid','on','YGrid','on','ZGrid','on');
    title('RF2 Data','FontSize',14,'FontWeight','bold'); hold on;
    xlabel('Lcc'); ylabel('Mcc'); zlabel('Percentage Correct')
    
    % Save User Data
    set(a.axes.conpanel,'UserData',conpanel)
    set(a.axes.glmppanel,'UserData',glmppanel)
    set(figure(1),'UserData',a)

end

function PlotPsychPred()
global gl GLMP

    % Load figure variables
    a = get(figure(1),'UserData');
    glmppanel = get(a.axes.glmppanel,'UserData');
    
    % Preallocate space
    AUC = nan(numel(GLMP.subunit{gl.whichsub}.uniqueLcc),1);
    
    % Run through each stimulus
    for s = 1:numel(GLMP.subunit{gl.whichsub}.uniqueLcc)
        
        idx = GLMP.subunit{gl.whichsub}.uniqueIdx{s};
        spHist = GLMP.subunit{gl.whichsub}.fr(idx)';
        blHist = GLMP.subunit{gl.whichsub}.blfr';
        thresh = unique([spHist blHist]);
        TPR = nan(1,length(thresh));
        FPR = TPR;
        
        for n=1:length(thresh)
            TPR(n) = sum(spHist >= thresh(n))/length(spHist);
            FPR(n) = sum(blHist >= thresh(n))/length(blHist);
        end
        
        AUC(s) = trapz(fliplr(FPR),fliplr(TPR));
        
    end
    
    % Linear Interpolation
    x = GLMP.subunit{gl.whichsub}.uniqueLcc;
    y = GLMP.subunit{gl.whichsub}.uniqueMcc;
    z =  AUC;
    F = TriScatteredInterp(x,y,z);
    [qx,qy] = meshgrid(min(x):.01:max(x),min(y):.01:max(y));
    
    % Plot figure
    axes(glmppanel.neuro.axes); cla; hold on;
    glmppanel.neuro.pts = plot3(x,y,z,'k*');
    glmppanel.neuro.surf = surfc(qx,qy,F(qx,qy));
    
    % Save User Data
    set(a.axes.glmppanel,'UserData',glmppanel)
    set(figure(1),'UserData',a)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Real-Time Updating of GUI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PlotPsych(Lcc,Mcc,rf,correct)
    
    % Load figure data
    a = get(figure(1),'UserData');
    glmppanel = get(a.axes.glmppanel,'UserData');
    
    if rf == 0
        rf = 2; % RF correct = RF1;
    end
    
    % Index/plot responses
    structs{1} = ['RF' num2str(rf)];
    structs{2} = 'AllData';
    
    for n = 1:2
        str = structs{n};
        stimL = ismember([glmppanel.(str).Lcc glmppanel.(str).Mcc],[Lcc Mcc],'rows');
        idx = find(stimL);
        glmppanel.(str).resps{idx} = [glmppanel.(str).resps{idx} correct];
        glmppanel.(str).meanresp(idx) = nanmean(glmppanel.(str).resps{idx});
        
        % Redraw surface for this RF
        x = glmppanel.(str).Lcc;
        y = glmppanel.(str).Mcc;
        z = glmppanel.(str).meanresp;
        F = TriScatteredInterp(x,y,z);
        [qx,qy] = meshgrid(min(x):.01:max(x),min(y):.01:max(y));
        
        % Plot figure
        axes(glmppanel.(str).axes); cla; hold on;
        glmppanel.(str).pts = plot3(x,y,z,'k*');
        glmppanel.(str).surf = surfc(qx,qy,F(qx,qy));
    end

    % Save figure data
    set(a.axes.glmppanel,'UserData',glmppanel)
    set(figure(1),'UserData',a)
        
end

