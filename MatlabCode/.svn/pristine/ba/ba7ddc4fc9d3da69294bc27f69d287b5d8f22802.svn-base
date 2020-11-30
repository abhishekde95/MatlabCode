function GridLMSubunitOnline_old()

% 2/10/13  New paradigm created.    JPW
% 3/2/13 Started adding GUI stuff  GDLH

global udpCom gl

disp('Starting GridLMSubunitOnline...')

%User-Defined Variables
gl.seed = 1;
gl.nStixelGrid = 15;
gl.DVAPerStixel = .075;
gl.StimDuration = .2; %In second
gl.nPres =5;
gl.nframesback = 10;
gl.lumcc = .7;
gl.colcc = .09;
gl.templumcc = gl.lumcc;
gl.tempcolcc = gl.colcc;
gl.maxcc = max(gl.lumcc,gl.colcc);
gl.reset = 0;


% Non-User-Defined Variables
gl.CurrentQueue = [];
gl.epoch = 1; % 1 = Dense Noise RF mapping, 2 = GridLMPlane Testing of Subunits
gl.rnd = 1; % for "phases" within a given epoch
gl.gridX = [];
gl.gridY = [];
gl.thetaspace = pi/4;
gl.rhospace = .5;
abortflag = 0;
gl.startnextepoch = 0;
gl.pixelmask.subunit1 = zeros(gl.nStixelGrid, gl.nStixelGrid);
gl.pixelmask.subunit2 = zeros(gl.nStixelGrid, gl.nStixelGrid);
currentspiketally = [];
stimon_t = [];
gl.stimdur = [];
gl.alreadyshown = [];
udpCom.sock = [];
udpCom.port = 6665;
udpCom.rexip = '192.168.1.120';

C = GridLMSubunitCodes();% defined in a separate m file
p = InitPStruct(0, C.HDRCOMPLETECD);
s = InitPlex();
initGUI();
InitStatsStruct();

codestruct = {...
        C.FPACQCD, @fpacqfn;...
        C.STIMONCD, @stimonfn;...
        C.ABORTCD, @abortfn;...
        C.REWCD, @stimofffn;...
        C.EOTCD, @eotfn};

    
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
    
    % Accumulating spikes
    if (any(p.spikes{1}))
        currentspiketally = [currentspiketally; p.spikes{1}];
        p.spikes{1} = [];
        if ~isempty(stimon_t) && gl.epoch == 1
            figure(1);
            a = get(figure(1),'UserData');
            conpanel = get(a.axes.conpanel,'UserData');
            spikecount = get(conpanel.uicontrols.spikecounter,'UserData') + 1;
            set(conpanel.uicontrols.spikecounter,'String',spikecount,'UserData',spikecount);
        end
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
        currentspiketally = [];
        p.lastprocessed_t = p.times(find(p.events == C.FPACQCD,1));
    end

    function stimonfn % STIMON
        stimon_t = p.times(find(p.events == C.STIMONCD,1));
        p.lastprocessed_t = p.times(find(p.events == C.STIMONCD,1));
    end

    function abortfn % ABORT
        abortflag = 1;
        currentspiketally = [];
        stimon_t = [];
        p.lastprocessed_t = p.times(find(p.events == C.ABORTCD,1,'last'));
    end

    function stimofffn % REW  % This may not be neccesary
        p.lastprocessed_t = p.times(find(p.events == C.REWCD,1));
    end

    function eotfn % EOT
        if ~abortflag
            l = GetValsIfPossible(C.LCCCD, p, 'float');
            m = GetValsIfPossible(C.MCCCD, p, 'float');
            x = GetValsIfPossible([C.GRIDXCD C.GRIDXCD], p, 'int');
            y = GetValsIfPossible([C.GRIDYCD C.GRIDYCD], p, 'int');
            nstix = GetValsIfPossible(C.NSTIXGRIDCD, p, 'int');
            epoch = GetValsIfPossible(C.EPOCHCD, p, 'int');
            if isempty(epoch)
                disp('Epoch not being sent over correctly...')
                keyboard
            end
            
            %Dense Noise
            if epoch == 1 
                a = get(figure(1),'UserData');
                conpanel = get(a.axes.conpanel,'UserData');
                if get(conpanel.uicontrols.resetbutton,'Value') == 0 %no reset
                    a.stats.nstixperside = nstix;
                    gl.nframes = GetValsIfPossible(C.NFRAMESCD, p, 'int');
                    bkgndrgb = a.stats.bkgndrgb;
                    stimcount = get(conpanel.uicontrols.stimcounter,'UserData') + gl.nframes;
                    set(conpanel.uicontrols.stimcounter,'String',stimcount,'UserData',stimcount);
                    set(figure(1),'UserData',a);
                    plotnow = UpdateSTX(gl.seed, gl.nframes, bkgndrgb, currentspiketally-stimon_t);
                    if (plotnow)
                        PlotSTA();
                        drawnow;
                    end
                else %Reset button has been pressed
                    fieldnames = {'pLpM','mLmM','pLmM','mLpM'};
                    labels = {'+L+M' '-L-M' '+L-M' '-L+M'};
                    gl.templumcc = str2double(get(conpanel.axes.lumdisplay,'string'));
                    gl.tempcolcc = str2double(get(conpanel.axes.coldisplay,'string'));
                    dnpanel = get(a.axes.dnpanel,'UserData');
                    
                    %Reset DN figs
                    for f = 1:numel(fieldnames)
                        for nfr = 1:gl.nframesback
                            a.stats.conenoise.spike.(fieldnames{f}).STS = zeros([3*gl.nStixelGrid^2 gl.nframesback]);
                            a.stats.conenoise.spike.(fieldnames{f}).nspikes = 0;
                            hax = dnpanel.axes.(fieldnames{f})(nfr);
                            im = reshape(a.stats.conenoise.spike.(fieldnames{f}).STS(:,n),gl.nStixelGrid,gl.nStixelGrid,3);
                            
                            %Redraw subunits
                            im(:,:,1) = im(:,:,1) + gl.pixelmask.subunit1 * .25;
                            im(:,:,2) = im(:,:,2) + gl.pixelmask.subunit2 * .25;
                            
                            % Actual image
                            h = image(im,'parent',hax);
                            set(h,'ButtonDownFcn',@PixelMaskCallback)
                            b.RGBmask = im;
                            set(hax,'UserData',b,'Xtick',[],'YTick',[])
                            
                            % Image Labels
                            if nfr == 1
                                dnpanel.labels.(fieldnames{f}) = text(-.2,.5,(labels{f}),'Rotation',90,'Units','normalized',...
                                    'HorizontalAlignment','center','Parent',hax);
                            end
                            if f == 1
                                t_label = [num2str(-(nfr-1) * a.stats.msperframe,'%.1f') ' ms'];
                                str = ['frame_' num2str(nfr)];
                                dnpanel.labels.(str) = text(.5,1.1,t_label,'Units','normalized',...
                                    'HorizontalAlignment','center','Parent',hax);
                            end
                        end
                    end
                    set(conpanel.uicontrols.spikecounter,'String',0,'UserData',0)
                    set(conpanel.uicontrols.resetbutton,'Value',0,'BackgroundColor',[1 1 1])
                    set(a.axes.conpanel,'UserData',conpanel)
                    set(figure(1),'UserData',a)
                end
                
            %GLMP
            elseif epoch == 2 
                if (~abortflag)
                    gl.presdur = GetValsIfPossible(C.STIMDURCD, p, 'float'); % in seconds
                    a = get(figure(1),'UserData');
                    conpanel = get(a.axes.conpanel,'UserData');
                    stimcount = get(conpanel.uicontrols.stimcounter,'UserData') + 1;
                    set(conpanel.uicontrols.stimcounter,'String',stimcount,'UserData',stimcount);
                    set(conpanel.uicontrols.spikecounter,'String',[num2str(numel(currentspiketally-stimon_t)/gl.presdur) ' spikes/s']);
                    set(figure(1),'UserData',a);
                    PlotGLMPMaps(currentspiketally-stimon_t,l,m,x,y)
                    gl.alreadyshown = cat(1,gl.alreadyshown,[l m]);
                end
            end
            currentspiketally = [];
            stimon_t = [];
            
        end
        
        % Handshaking
        p.lastprocessed_t =  p.times(find(p.events == C.EOTCD,1));
        pnet(udpCom.sock, 'write', 'plexdone>> >>');
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
    a.stats.murgb = [nan nan nan];
    a.stats.conenoise.spike.nspikes = 0;
    a.stats.conenoise.mu = [0 0 0];
    a.stats.conenoise.sigma = [7 7 7];
    a.stats.lmsbinaryrgbmat = nan*ones(4,3);
    
    fieldnames = {'pLpM','mLmM','pLmM','mLpM'};
    for f = 1:numel(fieldnames)
        a.stats.conenoise.spike.(fieldnames{f}).STS = zeros([3*gl.nStixelGrid^2 gl.nframesback]);
        a.stats.conenoise.spike.(fieldnames{f}).nspikes = 0;
    end
    set(figure(1),'UserData',a);
end

function InitStatsStruct()

    a = get(figure(1),'UserData');
    a.stats.npixperstix = 0; % Number of pixels per side of stixel  
    a.stats.nstixperside = 0; % Number of stixels per side of stimulus
    a.stats.msperframe = 0;  % ms per frame
    a.stats.gammaTable = []; % gamma tables
    a.stats.invgammaTable = []; % inverse gamma tables
    a.stats.monSpd = []; % monitor phosphor spectra
    a.stats.fundamentals = []; % cone fundamentals
    a.stats.M = [];     % guns to cones matrix
    a.stats.invM = [];  % cones to guns matrix
    a.stats.gotHeader = 0;
    set(figure(1),'UserData',a);
    
end

function allDone = dealWithMsgs()
global udpCom

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

    % Determine Epoch
    if gl.startnextepoch == 1
        gl.epoch = 2;
        gl.rnd = 1;
        gl.startnextepoch = 0;
    end

    %%%%% Epoch 1 %%%%%
    if gl.epoch == 1
        
        if gl.rnd == 1 %first time though
           
            fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
            fprintf('Epoch 1: Performing Dense Noise\n');
            fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n');
      
            % Setting up Dense Noise grid
            gl.masterGrid1D_C = -(gl.nStixelGrid-1)/2:(gl.nStixelGrid-1)/2;
            [gl.masterGrid_C_c,gl.masterGrid_C_r] = meshgrid(gl.masterGrid1D_C,fliplr(gl.masterGrid1D_C));

            gl.rnd = 2;
        
        else
            
            [~,gl.seed] = getEJrandnums3(gl.nStixelGrid.^2*gl.nframes,gl.seed);

        end
        
        % Stimulus variables
        gl.gridX = {gl.masterGrid_C_c(:)'};
        gl.gridY = {gl.masterGrid_C_r(:)'};
        stim = [gl.templumcc gl.tempcolcc 0 gl.DVAPerStixel gl.StimDuration gl.nStixelGrid gl.seed gl.epoch];
        
        
    %%%%%% Epoch 2 %%%%% 
    elseif gl.epoch == 2
        
        if isempty(gl.CurrentQueue)
            GLMPPhase();
        end
        
        %%%%%% Send Over 1st Entry in Current Queue %%%%%%
        if any(gl.CurrentQueue(:))
            stim = gl.CurrentQueue(1,:);
            gl.CurrentQueue(1,:) = [];
            
            a = get(figure(1),'UserData');
            conpanel = get(a.axes.conpanel,'UserData');
            set(conpanel.uicontrols.stimqueue,'string',size(gl.CurrentQueue,1))
            
        else
            disp('Nothing in Queue!')
            keyboard
        end
        
    end

end

function GLMPPhase()
global gl

        if gl.rnd == 1
            fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
            fprintf('Epoch 2: Performing GLMP on Subunits\n');
            fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n');
            
            UpdateGUIepoch2();
        end
        
        fprintf('Selecting Round %u Stimuli...\n',gl.rnd)
        [Lcc,Mcc] = ChooseLMStimuli();
        PairWithSubunits(Lcc,Mcc);
        fprintf('Round %u Stimuli Loaded into Queue. \n\n',gl.rnd);
        gl.rnd = gl.rnd+1;
        
end


function [Lcc,Mcc] = ChooseLMStimuli()
global gl

        % Construct Polar Grid
        if gl.rnd == 1
            thetas = shiftdim(0:gl.thetaspace:2*pi-gl.thetaspace);
            rhos = shiftdim(gl.rhospace:gl.rhospace:1);
            
        elseif mod(gl.rnd,2) == 1
            thetas = shiftdim(gl.thetaspace/2:gl.thetaspace:2*pi);
            rhos = shiftdim(gl.rhospace:gl.rhospace:1);
            gl.thetaspace = gl.thetaspace/2;
            
        elseif mod(gl.rnd,2) == 0
            thetas = shiftdim(gl.thetaspace:gl.thetaspace:2*pi);
            rhos = shiftdim(gl.rhospace/2:gl.rhospace:1);
            gl.rhospace = gl.rhospace/2;
            
        else
            disp('Problems with rounds in GLMP...')
            keyboard
            
        end
        
        % Totally random order
        PolRhoIdx = fullfact([numel(thetas) numel(rhos)]);
        rhothetalist = repmat([thetas(PolRhoIdx(:,1)) rhos(PolRhoIdx(:,2))],gl.nPres,1);
        tempqueue = rhothetalist(randperm(size(rhothetalist,1)),:);
        
        % Transform from polar to cartesian coordinates
        [tempLcc,tempMcc] = pol2cart(tempqueue(:,1),tempqueue(:,2));
        scale = gl.lumcc*gl.colcc./sqrt((gl.colcc.*cos(tempqueue(:,1)-pi/4)).^2 ...
            +(gl.lumcc.*sin(tempqueue(:,1)-pi/4)).^2);
        
    % Scale Cone Contrast Units for Monitor
    Lcc = tempLcc .* scale;
    Mcc = tempMcc .* scale;
    
end

function PairWithSubunits(Lcc, Mcc)
global gl

    subnames = {'subunit1' 'subunit2' 'subunit3'};
    whichsubs = find([sum(sum(gl.pixelmask.subunit1))...
        sum(sum(gl.pixelmask.subunit2))]);
    %if numel(whichsubs) == 1
        tempx = cell(numel(whichsubs),1);
        tempy = cell(numel(whichsubs),1);
    %else
    %    tempx = cell(numel(whichsubs)+1,1);
    %    tempy = cell(numel(whichsubs)+1,1);
    %end
    
    for s = 1:numel(whichsubs)
        [tempr,tempc] = ind2sub(size(gl.pixelmask.subunit1),find(gl.pixelmask.(subnames{whichsubs(s)})));
        tempx{s} = tempc' - max(gl.masterGrid1D_C) - 1;
        tempy{s} = (tempr' - max(gl.masterGrid1D_C) - 1) * -1;
        
    end
    %tempx{end} = cat(2,tempx{:});
    %tempy{end} = cat(2,tempy{:});
    
    LMxyIdx = fullfact([numel(Lcc) numel(tempx)]);
    LMxyIdx_rand = LMxyIdx(randperm(size(LMxyIdx,1)),:);
    
    % Set up queue
    param_LM = [Lcc(LMxyIdx_rand(:,1)) Mcc(LMxyIdx_rand(:,1))];
    param_S = zeros(size(LMxyIdx,1),1);
    param_nstixelgrid = repmat(gl.nStixelGrid,size(LMxyIdx,1),1);
    param_dvaperstixel = repmat(gl.DVAPerStixel,size(LMxyIdx,1),1);
    param_stimdur = repmat(gl.StimDuration,size(LMxyIdx,1),1);
    param_seed = repmat(gl.seed,size(LMxyIdx,1),1);
    param_epoch = repmat(gl.epoch,size(LMxyIdx,1),1);
    gl.CurrentQueue = cat(2,param_LM,param_S,param_dvaperstixel, param_stimdur,...
        param_nstixelgrid,param_seed,param_epoch);
    gl.gridX = tempx(LMxyIdx_rand(:,2));
    gl.gridY = tempy(LMxyIdx_rand(:,2));

end

function gridx = getX()
global gl

    % Send over first row in coordinate queue, then delete
    if iscell(gl.gridX(1)) 
        gridx = gl.gridX{1};
        gl.gridX = gl.gridX(2:end);
    else
        gridx = gl.gridX(1,:);
        gl.gridX(1,:) = [];
    end
    
end

function gridy = getY()
global gl

    % Send over first row in coordinate queue, then delete
    if iscell(gl.gridY(1)) %In Epoch 3
        gridy = gl.gridY{1};
        gl.gridY = gl.gridY(2:end);
    else
        gridy = gl.gridY(1,:);
        gl.gridY(1,:) = [];
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting Up / Updating GUI's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function initGUI()
    global gl

    % GUI variables
    fieldnames = {'pLpM','mLmM','pLmM','mLpM'};
    labels = {'+L+M' '-L-M' '+L-M' '-L+M'};
    initBGcol = [.6 .6 .6];
    panelsize = 75;
    panelXpos = linspace(45,825,10);
    panelYpos = linspace(265,15,4);
    
    % Set up GUI Pannel
    close all;
    figure(1);
    set(figure(1),'DefaultAxesUnits','normalized','position',[100 100 1200 800],'UserData',[])
    a = get(figure(1),'UserData');
    
    % Set up subpannels
    a.axes.conpanel = uipanel('Title','Control Panel','Position',[.015 .515 .185 .480],'BorderType','beveledin','BackgroundColor',[.9 .9 .9]);
    a.axes.dnpanel = uipanel('Title','Dense Noise','Position',[.21 .515 .775 .480],'BorderType','beveledin','BackgroundColor',[.9 .9 .9]);
    a.axes.glmppanel = uipanel('Title','GLMP','Position',[.015 .015 .970 .490],'BorderType','beveledin','BackgroundColor',initBGcol);
    
    
    % Set Up Control Panel
    % Subunit Selector
    conpanel = get(a.axes.conpanel,'UserData');
    conpanel.axes.selector = uibuttongroup('Parent',a.axes.conpanel,'Units','normalized',...
        'visible','off','Position',[.05 .75 .9 .2],'Title','Stimulus Selector',...
        'BorderType','etchedout','Visible','on','BackgroundColor',[.9 .9 .9],...
        'Title','Subunit Selection','TitlePosition','centertop');
    conpanel.uicontrols.sub1 = uicontrol('Parent',conpanel.axes.selector,'Units','normalized',...
        'Style','radiobutton','String','Subunit 1','FontSize',13,...
        'pos',[.25 .6 .5 .3],'BackgroundColor',[.9 .9 .9]);
    conpanel.uicontrols.sub2 = uicontrol('parent',conpanel.axes.selector,'Units','normalized',...
        'Style','radiobutton','String','Subunit 2','FontSize',13,...
        'pos',[.25 .2 .5 .3],'BackgroundColor',[.9 .9 .9]);
    set(conpanel.axes.selector,'SelectedObject',conpanel.uicontrols.sub1,'SelectionChangeFcn',@switchstim);
    
    % Stimulus Sliders
    conpanel.axes.sliders = uipanel('Parent',a.axes.conpanel,'Position',[.05 .325 .9 .4],'Title','White Noise Stimulus Parameters',...
        'BorderType','etchedout','BackgroundColor',[.9 .9 .9]);
    conpanel.uicontrols.lumslider = uicontrol('Parent',conpanel.axes.sliders,'Units','normalized',...
        'Style','slider','min',0,'max',gl.lumcc,'value',gl.lumcc,'SliderStep',[1/18 1/18],...
        'pos',[.05 .25 .9 .1],'BackgroundColor',[.9 .9 .9],'Callback',@adjustlumcc);
    conpanel.labels.lumslider = uicontrol('Parent',conpanel.axes.sliders,'Units','normalized',...
        'style','edit','string','Luminance Stixel Value','FontSize',8,...
        'HorizontalAlignment','center','Position',[.05 .375 .7 .15],'BackgroundColor',[.9 .9 .9]);
    conpanel.axes.lumdisplay = uicontrol('Parent',conpanel.axes.sliders,'Units','normalized',...
        'style','edit','string',gl.templumcc,'FontSize',10,...
        'HorizontalAlignment','center','Position',[.75 .375 .2 .15],'BackgroundColor',[.9 .9 .9]);
    conpanel.uicontrols.colslider = uicontrol('Parent',conpanel.axes.sliders,'Units','normalized',...
        'Style','slider','min',0,'max',gl.colcc,'value',gl.colcc,'SliderStep',[1/9 1/9],...
        'pos',[.05 .6 .9 .1],'BackgroundColor',[.9 .9 .9],...
        'Callback',@adjustcolcc);
    conpanel.labels.colslider = uicontrol('Parent',conpanel.axes.sliders,'Units','normalized',...
        'style','edit','string','Chromatic Stixel Value','FontSize',8,...
        'HorizontalAlignment','center','Position',[.05 .725 .7 .15],'BackgroundColor',[.9 .9 .9]);
    conpanel.axes.coldisplay = uicontrol('Parent',conpanel.axes.sliders,'Units','normalized',...
        'style','edit','string',gl.tempcolcc,'FontSize',10,...
        'HorizontalAlignment','center','Position',[.75 .725 .2 .15],'BackgroundColor',[.9 .9 .9]);
    conpanel.uicontrols.resetbutton = uicontrol('Parent',conpanel.axes.sliders,'Units','normalized',...
        'position',[.05 .05 .075 .101],'style','togglebutton','Callback',@switchstim);
    conpanel.labels.resetbutton = uicontrol('Parent',conpanel.axes.sliders,'Units','normalized',...
        'style','edit','string','Reset Stimulus and Statistics','FontSize',8,...
        'HorizontalAlignment','center','Position',[.15 .05 .8 .1],'BackgroundColor',[.9 .9 .9]);
    
    % Epoch Switch
    conpanel.uicontrols.startnextepoch = uicontrol('Parent',a.axes.conpanel,'Units','pixels',...
        'style','togglebutton','Min',0,'Max',1,'Value',0,...
        'HorizontalAlignment','left',...
        'Position',[18 18 20 20],...
        'Callback',@startnextepochCallback);
    conpanel.labels.startnextepoch = uicontrol('Parent',a.axes.conpanel,'Units','normalized',...
        'style','edit','string','Start GridLMPlane','FontSize',10,...
        'HorizontalAlignment','center',...
        'Position',[.2 .05 .7 .05],'BackgroundColor',[.9 .9 .9]);
    
    % Counters
    b.spikecounter = 0;
    b.nstimleft = [];
    b.nstimpresented = 0;
    conpanel.uicontrols.spikecounter = uicontrol('Parent',a.axes.conpanel,'Units','normalized',...
        'style','edit','string','0','FontSize',12,...
        'Position',[.2 .16 .6 .08],'BackgroundColor',[.8 .8 .8],...
        'UserData',b.spikecounter);
    conpanel.labels.spikecounter = uicontrol('Parent',a.axes.conpanel,'Units','normalized',...
        'style','text','string','Number of Spikes Collected','FontSize',10,...
        'HorizontalAlignment','center',...
        'Position',[.1 .25 .8 .05],'BackgroundColor',[.9 .9 .9]);
    conpanel.uicontrols.stimqueue = uicontrol('Parent',a.axes.conpanel,'Units','normalized',...
        'style','edit','string','~','FontSize',12,...
        'Position',[.2 .36 .6 .08],'BackgroundColor',[.8 .8 .8],...
        'Visible','off','UserData',b.nstimleft);
    conpanel.labels.stimqueue = uicontrol('Parent',a.axes.conpanel,'Units','normalized',...
        'style','text','string','Number of Stimuli in Queue','FontSize',10,...
        'HorizontalAlignment','center',...
        'Visible','off','Position',[.1 .45 .8 .05],'BackgroundColor',[.9 .9 .9]);
    conpanel.uicontrols.stimcounter = uicontrol('Parent',a.axes.conpanel,'Units','normalized',...
        'style','edit','string','0','FontSize',12,...
        'Position',[.2 .56 .6 .08],'BackgroundColor',[.8 .8 .8],...
        'Visible','off','UserData',b.nstimpresented);
    conpanel.labels.stimcounter = uicontrol('Parent',a.axes.conpanel,'Units','normalized',...
        'style','text','string','Number of Frames Presented','FontSize',10,...
        'HorizontalAlignment','center',...
        'Visible','off','Position',[.1 .65 .8 .05],'BackgroundColor',[.9 .9 .9]);
    
    %Save User Data
    set(a.axes.conpanel,'UserData',conpanel)
    
    
    % Set Up Dense Noise Panel
    dnpanel = get(a.axes.dnpanel,'UserData');
        
    % Set up Dense Noise figures
    for ypos = 1:numel(panelYpos)
        field = fieldnames{ypos};
        for xpos = 1:numel(panelXpos)
            dnpanel.axes.(field)(xpos) = axes('Parent',a.axes.dnpanel,'Units','pixels',...
                'Xtick',[],'Ytick',[],'box','on',...
                'Position',[panelXpos(xpos) panelYpos(ypos) panelsize panelsize]);
            if xpos == 1
                figure(1);
                dnpanel.labels.(field) = text(-.2,.5,(labels{ypos}),'Rotation',90,'Units','normalized',...
                    'HorizontalAlignment','center');
            end
        end
    end
    
    % Save User Data
    set(a.axes.dnpanel,'UserData',dnpanel)
    
    % GLMP Panel
    figure(1);
    glmppanel = get(a.axes.glmppanel,'UserData');
    glmppanel.axes.sub1 = axes('Parent',a.axes.glmppanel,'Units','normalized',...
        'OuterPosition',[.05 .05 .3 .9],...
        'CameraPosition',[-0.399233 -0.576688 199.22],...
        'XLim',[-gl.maxcc gl.maxcc],'YLim',[-gl.maxcc gl.maxcc],'ZLim',[0 25],...
        'XGrid','on','YGrid','on','ZGrid','on',...
        'box','on');
    set(glmppanel.axes.sub1,'Color',[.3 .3 .3])
    xlabel('Lcc'); ylabel('Mcc'); zlabel('Firing Rate')
    title('Subunit 1','FontSize',12,'FontWeight','bold')
    glmppanel.axes.sub2 = axes('Parent',a.axes.glmppanel,'Units','normalized',...
        'OuterPosition',[.35 .05 .3 .9],...
        'CameraPosition',[-0.399233 -0.576688 199.22],...
        'XLim',[-gl.maxcc gl.maxcc],'YLim',[-gl.maxcc gl.maxcc],'ZLim',[0 25],...
        'XGrid','on','YGrid','on','ZGrid','on',...
        'box','on');    
    set(glmppanel.axes.sub2,'Color',[.3 .3 .3])
    xlabel('Lcc'); ylabel('Mcc'); zlabel('Firing Rate')
    title('Subunit 2','FontSize',12,'FontWeight','bold')
    glmppanel.axes.sub3 = axes('Parent',a.axes.glmppanel,'Units','normalized',...
        'OuterPosition',[.65 .05 .3 .9],...
        'CameraPosition',[-0.399233 -0.576688 199.22],...
        'XLim',[-gl.maxcc gl.maxcc],'YLim',[-gl.maxcc gl.maxcc],'ZLim',[0 25],...
        'XGrid','on','YGrid','on','ZGrid','on',...
        'box','on');
    set(glmppanel.axes.sub3,'Color',[.3 .3 .3])
    xlabel('Lcc'); ylabel('Mcc'); zlabel('Firing Rate')
    title('Subunits 1 and 2','FontSize',12,'FontWeight','bold')

    % Save User Data
    set(a.axes.glmppanel,'UserData',glmppanel)
    set(figure(1),'UserData',a)
    
    function adjustlumcc(~,~)
        newcc = get(conpanel.uicontrols.lumslider,'value');
        set(conpanel.axes.lumdisplay,'string',newcc)
        set(a.axes.conpanel,'UserData',conpanel)
    end

    function adjustcolcc(~,~)
        newcc = get(conpanel.uicontrols.colslider,'value');
        set(conpanel.axes.coldisplay,'string',newcc)
        set(a.axes.conpanel,'UserData',conpanel)
    end

    function switchstim(~,~)
        disp('Resetting statistics to prepare for new stimulus...')
        set(conpanel.uicontrols.resetbutton,'BackgroundColor',[1 0 0])
    end

    function startnextepochCallback(~,~)
        a = get(figure(1),'UserData');
        conpanel = get(a.axes.conpanel,'UserData');
        if any(any(gl.pixelmask.subunit1)) || any(any(gl.pixelmask.subunit2))
            
            if (get(conpanel.uicontrols.startnextepoch,'Value')) == 1;
                set(conpanel.uicontrols.startnextepoch,'BackgroundColor',[0 1 0])
                gl.startnextepoch = 1;
                set(conpanel.uicontrols.startnextepoch,'Enable','off')
            else
                set(conpanel.uicontrols.startnextepoch,'BackgroundColor',[0.92549 0.913725 0.847059])
                gl.startnextepcoh = 0;
            end
        else
            fprintf('\n***** Must Select Subunits Before Proceeding to GLMP! *****\n\n')
        end
    end

end

function UpdateGUIepoch2()

    a = get(figure(1),'UserData');
    conpanel = get(a.axes.conpanel,'UserData');
    glmppanel = get(a.axes.glmppanel,'UserData');
    
    set(conpanel.uicontrols.startnextepoch,'Enable','off','BackgroundColor',[.3 .3 .3])
    set(conpanel.labels.startnextepoch,'string','GLMP in Progress...','BackgroundColor',[.8 .8 .8]);
    set(a.axes.dnpanel,'BackgroundColor',[.8 .8 .8]);
    set(a.axes.glmppanel,'BackgroundColor',[.9 .9 .9]);
    set(conpanel.axes.selector,'BackgroundColor',[.8 .8 .8])
    set(conpanel.uicontrols.sub1,'Enable','off','BackgroundColor',[.8 .8 .8])
    set(conpanel.uicontrols.sub2,'Enable','off','BackgroundColor',[.8 .8 .8])
    set(conpanel.axes.selector,'SelectedObject',[]);
    set(conpanel.labels.stimcounter,'string','Number of Stimuli Presented')
    set(conpanel.uicontrols.stimcounter,'string','0','UserData',0)
    set(conpanel.labels.spikecounter,'string','Response to Last Stimulus')
    set(conpanel.uicontrols.spikecounter,'string','0','UserData',0)

    %Set up User Data fields
    p.LM = [];
    p.nspikes = cell([0 0]);
    p.meannspikes = [];
    p.frs = cell([0 0]);
    p.meanfrs = [];
    set(glmppanel.axes.sub1,'Color',[1 1 1],'UserData',p)
    set(glmppanel.axes.sub2,'Color',[1 1 1],'UserData',p)
    set(glmppanel.axes.sub3,'Color',[1 1 1],'UserData',p)
        
    
end

function selcbk(~,~)
    %Can execute a function for selecting new radio button here
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotnow = UpdateSTX(seed, nframes, bkgndrgb, spiketimes)
global gl 

    a = get(figure(1),'UserData');
    stats = a.stats;
    
    nframesback = gl.nframesback;
    stats.murgb = bkgndrgb;
    stats.conenoise.mu = 0;
    if gl.templumcc == 0 && gl.tempcolcc > 0
        colordirlms = [1 -1 0; -1 1 0];
    elseif gl.templumcc > 0 && gl.tempcolcc == 0
        colordirlms = [1 1 0; -1 -1 0];
    else
        colordirlms = [1 1 0; -1 -1 0; 1 -1 0; -1 1 0];
    end
    randnums = getEJrandnums(stats.nstixperside^2*nframes, seed);
    randnums = mod(randnums, size(colordirlms,1))+1;
    randnums_orig = colordirlms(randnums,:);
    
    %Seperate +/- lum/chrom conditions...
    pLpML = find(randnums_orig(:,1) == 1 & randnums_orig(:,2) == 1);
    mLmML = find(randnums_orig(:,1) == -1 & randnums_orig(:,2) == -1);
    pLmML = find(randnums_orig(:,1) == 1 & randnums_orig(:,2) == -1);
    mLpML = find(randnums_orig(:,1) == -1 & randnums_orig(:,2) == 1);

    idxnames = {'pLpML','mLmML','pLmML','mLpML'};
    fieldnames = {'pLpM','mLmM','pLmM','mLpM'};
    for f = 1:numel(fieldnames)
        
        randnums = zeros(size(randnums_orig));
        randnums(eval(idxnames{f}),:) = randnums_orig(eval(idxnames{f}),:);
        
        randnums = reshape(randnums, [stats.nstixperside^2, nframes, 3]);
        randnums = permute(randnums,[1 3 2]);
        randnums = reshape(randnums, [stats.nstixperside^2*3, nframes]);
        % Each column should be Ls followed by Ms followed by Ss for each frame.
        
        if (any(spiketimes) && nframes ~= 0)
            try
                frametimes = linspace(0, nframes*stats.msperframe, nframes)+(stats.msperframe/2)';
                spiketimes(spiketimes < nframesback*stats.msperframe) = [];
                spiketimes(spiketimes > frametimes(end)) = [];
            catch
                disp('Spiketimes error...')
                keyboard
            end
            [n,x] = hist(spiketimes, frametimes);
            %STCOVmex('init',{stats.nstixperside^2,3,nframesback});
            %STCOVmex(randnums(:), n);
            %out = STCOVmex('return');
            STAmex('init',{stats.nstixperside^2,3,nframesback});
            STAmex(randnums(:), n);
            out = STAmex('return');
            stats.conenoise.spike.(fieldnames{f}).STS = stats.conenoise.spike.(fieldnames{f}).STS+out{1};
            stats.conenoise.spike.(fieldnames{f}).nspikes = stats.conenoise.spike.(fieldnames{f}).nspikes+out{3};
        end
        
        a.stats = stats;
        set(figure(1),'UserData',a);
        
    end
    plotnow = 1;
    
end


function PlotSTA()
global gl
    
    a = get(figure(1),'UserData');
    stats = a.stats;
    dnpanel = get(a.axes.dnpanel,'UserData');
 
    maxSTAval = nan(1,4);
    minSTAval = nan(1,4);
    fieldnames = {'pLpM' 'mLmM' 'pLmM' 'mLpM'};
    labels = {'+L+M' '-L-M' '+L-M' '-L+M'};
    for fn = 1:numel(fieldnames)
        
        STS = stats.conenoise.spike.(fieldnames{fn}).STS;
        n = stats.conenoise.spike.(fieldnames{fn}).nspikes;
        
        STAs.(fieldnames{fn}) = STS./n;
        
        maxSTAval(fn) = max(max(abs(STAs.(fieldnames{fn}))));
        minSTAval(fn) = min(min(abs(STAs.(fieldnames{fn})(1:gl.nStixelGrid^2,:))));
        
    end
    
    maxval = max(maxSTAval);
    minval = min(minSTAval);
    for fn = 1:numel(fieldnames)
        
        for n = 1:gl.nframesback
            
            hax = dnpanel.axes.(fieldnames{fn})(n);
            im = reshape(STAs.(fieldnames{fn})(:,n),gl.nStixelGrid,gl.nStixelGrid,3);
            im(:,:,3) = im(:,:,1);
            im = ((abs(im)-minval)./(maxval-minval)) * .75;
            
            %Redraw subunits
            im(:,:,1) = im(:,:,1) + gl.pixelmask.subunit1 * .25;
            im(:,:,2) = im(:,:,2) + gl.pixelmask.subunit2 * .25;
            
            % Actual image
            h = image(im,'parent',hax);
            set(h,'ButtonDownFcn',@PixelMaskCallback)
            b.RGBmask = im;
            set(hax,'UserData',b,'Xtick',[],'YTick',[])
            
            % Image Labels
            if n == 1
                dnpanel.labels.(fieldnames{fn}) = text(-.2,.5,(labels{fn}),'Rotation',90,'Units','normalized',...
                    'HorizontalAlignment','center','Parent',hax);
            end
            if fn == 1
                t_label = [num2str(-(n-1) * stats.msperframe,'%.1f') ' ms'];
                str = ['frame_' num2str(n)];
                dnpanel.labels.(str) = text(.5,1.1,t_label,'Units','normalized',...
                    'HorizontalAlignment','center','Parent',hax);
            end
            
        end
        
    end
    
    % Nested callback for clicking on RF image
    function PixelMaskCallback(~,~)
        
        try
            a = get(figure(1),'UserData');
            conpanel = get(a.axes.conpanel,'UserData');
            dnpanel = get(a.axes.dnpanel,'UserData');
            
            whichpt = get(gca,'CurrentPoint');
            whichpt = round(whichpt(1,[1 2]));
            whichpt = min([whichpt; [size(gl.pixelmask.subunit1,1) size(gl.pixelmask.subunit1,2)]]);
            
            whichsubunit = get(conpanel.axes.selector,'SelectedObject');
            
            if whichsubunit == conpanel.uicontrols.sub1
                subID = 'subunit1';
                colorscheme = [.25 0 0];
            elseif whichsubunit == conpanel.uicontrols.sub2
                subID = 'subunit2';
                colorscheme = [0 .25 0];
            elseif isempty(whichsubunit)
                fprintf('Must Select a Subunit before selecting stixels!')
                
            else
                disp('Something wrong with subunit indexing...')
                keyboard
            end
            
            gl.pixelmask.(subID)(whichpt(2),whichpt(1)) = ~gl.pixelmask.(subID)(whichpt(2), whichpt(1));
            
            fieldnames = {'pLpM','pLmM','mLpM','mLmM'};
            labels = {'+L+M' '-L-M' '+L-M' '-L+M'};
            
            for f = 1:numel(fieldnames)
                for j = 1:gl.nframesback
                    hax = dnpanel.axes.(fieldnames{f})(j);
                    b = get(hax,'UserData');
                    im = b.RGBmask;
                    if (gl.pixelmask.(subID)(whichpt(2), whichpt(1))) % Selecting a stixel
                        im(whichpt(2),whichpt(1),:) = shiftdim(im(whichpt(2),whichpt(1),:)) + colorscheme';
                        
                    else % Deselecting a stixel
                        im(whichpt(2),whichpt(1),:) = shiftdim(im(whichpt(2),whichpt(1),:)) - colorscheme';
                        
                    end
                    h = image(im,'parent',hax);
                    set(h,'ButtonDownFcn',@PixelMaskCallback);
                    b.RGBmask = im;
                    set(hax,'UserData',b,'Xtick',[],'YTick',[]);
                    
                    % Image Labels
                    if j == 1
                        dnpanel.labels.(fieldnames{f}) = text(-.2,.5,(labels{f}),'Rotation',90,'Units','normalized',...
                            'HorizontalAlignment','center','Parent',hax);
                    end
                    if f == 1
                        t_label = [num2str(-(j-1) * a.stats.msperframe,'%.1f') ' ms'];
                        str = ['frame_' num2str(j)];
                        dnpanel.labels.(str) = text(.5,1.1,t_label,'Units','normalized',...
                            'HorizontalAlignment','center','Parent',hax);
                    end
                end
            end
        catch
            disp('Pixelmask is not working...')
            keyboard
        end
        drawnow;
    end
    
end

function PlotGLMPMaps(spiketimes, l, m, x, y)
global gl
    
    a = get(figure(1),'UserData');
    glmppanel = get(a.axes.glmppanel,'UserData');

    nspikes = numel(spiketimes);

    % Find which subunit was flashed
    c = x + max(gl.masterGrid1D_C) + 1;
    r = abs(y - max(gl.masterGrid1D_C) - 1);
    rc = sort([r c],1);
    [sub1r,sub1c] = find(gl.pixelmask.subunit1);
    sub1rc = sort([sub1r sub1c],1);
    [sub2r,sub2c] = find(gl.pixelmask.subunit2);
    sub2rc = sort([sub2r sub2c],1);
    sub3rc = sort(cat(1,sub1rc,sub2rc),1);
    
    %keyboard
    if isequal(rc,sub1rc)
        hax = glmppanel.axes.sub1;
        plotcol = 'r';
    elseif isequal(rc,sub2rc)
        hax = glmppanel.axes.sub2;
        plotcol = 'g';
    elseif isequal(rc,sub3rc)
        hax = glmppanel.axes.sub3;
        plotcol = 'b';
    else
        disp('Subuninit ID not working...')
        keyboard
    end
    
    p = get(hax,'UserData');
    
    try
        stimIdx = find(ismember(p.LM,[l m],'rows'));
        if isempty(stimIdx)
            stimIdx = size(p.LM,1) + 1;
            p.LM(stimIdx,1:2) = [l m];
            p.nspikes{stimIdx,1} = nspikes;
            p.meannspikes(stimIdx,1) = nspikes;
            p.frs{stimIdx,1} = p.nspikes{stimIdx}/gl.presdur;
            p.meanfrs(stimIdx,1) = p.meannspikes(stimIdx)/gl.presdur;
        else
            p.nspikes{stimIdx} = cat(2,p.nspikes{stimIdx},nspikes);
            p.meannspikes(stimIdx) = mean(p.nspikes{stimIdx,:});
            p.frs{stimIdx} = p.nspikes{stimIdx}./gl.presdur;
            p.meanfrs(stimIdx) = p.meannspikes(stimIdx)/gl.presdur;
        end
        maxfr = max(p.meanfrs);
        
    catch
        disp('Something wrong in PlotGLMP...')
        keyboard
    end
    
    if maxfr == 0
        set(hax,'UserData',p)
    else
        set(hax,'UserData',p,'ZLim',[0 maxfr*1.1])
    end
    axes(hax);
    tempview = get(hax,'View');
    cla; hold on;
    stem3(p.LM(:,1),p.LM(:,2),p.meanfrs,[plotcol '*'])
    set(hax,'View',tempview)
    drawnow
        
end

