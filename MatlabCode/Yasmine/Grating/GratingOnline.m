function GratingOnline()

%% CHARACTERIZING STIMULUS TUNING WITH GRATINGS
% YASMINE EL-SHAMAYLEH (5/19/14)

% Written by GDLH 3/3/08
% Edited by YES for V1-V2 pathway 5/17/14 

% Online analysis stuff for the grating paradigm.
% Unlike WhiteNoiseOnline.m, this function keeps going 
% around the main loop as quickly as possible without
% waiting for any particular code to be dropped.  This
% will be important if we ever want to do anything in
% "real time" (ie within a trial).

% NEW LIST OF PROTOCOLS (YES)
% 0) Null protocol               (all done)
% 11) Orientation (coarse)       (move to 12)
% 12) Spatial frequency          (move to 13)
% 13) Temporal frequency         (move to 14)
% 14) Aperture size              (move to 15)
% 15) Orientation (fine)         (move to 16)
% 16) Classical color tuning     (move to 17)
% 17) Contrast response function with/out optical stim on half trials

% STIMULUS TYPES (I'M ONLY USING 0,2 HERE)
% 0 = grating
% 1 = gabor
% 2 = grating + optical stim

% YES: REPLACED COLOR RADIO BUTTON WITH OSTIM RADIO BUTTON (FOR OPTICAL STIM)

%% BASIC CODE 

% Global Variables
global gl;
global udpCom;                          % Needed by "DealWithMessage()" and associated subfunctions
udpCom.sock = [];
udpCom.port = 6665;
udpCom.rexip = '192.168.1.120';

% Some initializations
SetUpFigure;
C = GratingCodes;                       % Script that defines a bunch of constants
TS = TrialSpecCodes;                    % Another script that defines a bunch of constants
p = InitPStruct(0, C.FIXXCD);           % Structure of the events from each trial
s = InitPlex();                         % Link to Plexon datastream
InitGLStruct();                         % Structure that holds the statistics
[sock, Success] = pnetStart(6665);      % UDP communication with REX
currentspiketally = [];
stimon_t = [];
fix_t = [];

% The main loop
while (~gl.timetoleave)
    gl.timetoleave = CheckForESCKey();
    
    % Check for a message from REX
    msgSize = pnet(udpCom.sock, 'readpacket', 200, 'noblock');
    if(msgSize)
        DealWithMessage(msgSize);
    end
    [n, eventList] = PL_GetTS(s);
    if (n > 0)
        p = ProcessEventList(p, eventList);
    end
    if (~gl.gotheader)
        p = GetHeader(p); % "SetUpTrials" called from here if possible
    end
    if (any(p.events == C.FPACQCD))
        fix_t = p.times(find(p.events == C.FPACQCD,1));
        abortflag = 0;
    end
    if (any(p.events == C.STIMONCD))
        stimon_t = p.times(find(p.events == C.STIMONCD,1));
        nspikes = sum(currentspiketally > 0); % Spike times are relative to fix_t
        gl.trialspecs{gl.paramsidx,TS.BASELINE} =...
            [gl.trialspecs{gl.paramsidx,TS.BASELINE}; nspikes/(stimon_t-fix_t)*1000];
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
        ClearRaster;
        p.lastprocessed_t = p.times(find(p.events == C.ABORTCD,1,'last'));
    end
    if (any(p.events == C.STIMOFFCD))
        stimoff_t = p.times(find(p.events == C.STIMOFFCD,1));
        stimdur = stimoff_t-stimon_t;
        nspikes = sum(currentspiketally > stimon_t-fix_t);
        gl.trialspecs{gl.paramsidx,TS.SPIKES} = [gl.trialspecs{gl.paramsidx,TS.SPIKES}; nspikes/stimdur*1000];
        currentspiketally = [];
        stimon_t = [];
        fix_t = [];
        if (gl.protocol == 16)                      %
            PlotColorTuning;
        elseif (gl.protocol == 17)
            PlotContrastResponse;
        else
            PlotTuningCurve;
        end
        ClearRaster;
        UpdateTrialCountdown;
        p.lastprocessed_t = stimoff_t;
    end
    p = CleanUpEvents(p);
end % big while loop
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UDP COMMUNICATION FUNCTIONS BELOW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DealWithMessage

function DealWithMessage(msgSize)
    global gl;
    global udpCom;
    
    message = pnet(udpCom.sock, 'read', msgSize, 'char');
    if (strncmp(message,'return',6))
        a = dbstack; % Check whether called from another function or from command line 
        if (~strcmp(a(end).name, mfilename))
            gl.timetoleave = 1;
        end
    end
    try
        if (strcmp(message,'sendToRex(udpCom, gl.gratingparams, ''double'', message);'))
            PickGratingParams;
        end
        eval(message);
    catch
        fprintf('Unknown message: %s\n', message);
        error = lasterror;
        disp(error.message);
        disp(error.identifier);
        for i = 1:length(error.stack)
            disp(error.stack(i));
        end 
    end
end

%% InitGLStruct
 
function InitGLStruct()
    global gl;
%     if (~isfield(gl,'preforient'))
%         gl.preforient = nan;
%     end
%     if (~isfield(gl,'prefsf'))
%         gl.prefsf = nan;
%     end
%     if (~isfield(gl,'preftf'))          
%         gl.preftf = nan;
%     end
%     if (~isfield(gl,'prefdiam'))
%         gl.prefdiam = nan;
%     end
%     if (~isfield(gl,'prefccdir'))
%         gl.prefccdir = [0 0 0];         
%     end
    
    gl.preforient = NaN;
    gl.prefsf = NaN;
    gl.preftf = NaN;
    gl.prefdiam = NaN;
    gl.prefccdir = [0 0 0];
    
    gl.splineargmax = [];
    gl.disp.gammaTable = [];
    gl.disp.monSpd = [];
    gl.disp.fundamentals = [];
    gl.disp.bkgndrgb = [];    
    gl.ntrlspercond = [];
    gl.protocol = [];
    gl.gratingparamcols = [];   % which columns of gl.trialspecs contain the grating parameters
    gl.gratingparams = [];      % For sending grtaing parameters to REX.  Updated every trial.
    gl.trialspecs = {};         % A structure containing all the grating parameters used in a given protocol. Updated each time the protocol is changed.
    gl.paramsidx = [];          % index into gl.trialspecs indicating the current trial(from which gl.gratingparams is extracted)
    gl.gotheader = 0;
    gl.timetoleave = 0;         % Flag when to leave this .m file
    gl.stimon_t = [];
    gl.prefisolumcolors = [0 0 0];
end

%% GetHeader

% Getting the header data from the ecodes
% Destructively modifies the p structure to update p.lastprocessed_t.
% Calls SetUpTrials if the complete header has been obtained. 

function p = GetHeader(p)
    global gl;
    
    C = GratingCodes;
    
    [out, p] = GetValsIfPossible([C.GAMMATABLECD C.GAMMATABLECD], p);
    if ~isempty(out)
        disp('Got gammatable');
        gl.disp.gammaTable = reshape(out, length(out)/3,3);
        gl.disp.invgammaTable = InvertGamma(gl.disp.gammaTable,1);
    end

    [out, p] = GetValsIfPossible([C.MONSPDCD C.MONSPDCD], p);
    if ~isempty(out)
        disp('Got monSPD');
        gl.disp.monSpd = reshape(out, length(out)/3,3);
    end

    [out, p] = GetValsIfPossible([C.FUNDAMENTALSCD C.FUNDAMENTALSCD], p);
    if ~isempty(out)
        disp('Got fundamentals');
        gl.disp.fundamentals = reshape(out, length(out)/3,3);
        if (size(gl.disp.monSpd, 1) == 101)
            P_device = SplineSpd([380:4:780]',gl.disp.monSpd,[380:5:780]');
        else
            P_device = SplineSpd([380:2:780]',gl.disp.monSpd,[380:5:780]');
        end
        gl.disp.M = gl.disp.fundamentals'*P_device;
        gl.disp.invM = inv(gl.disp.M);
    end

    [out, p] = GetValsIfPossible(C.FRAMERATECD, p, 'double');
    if (~isempty(out))
        disp('Got framerate');
        gl.disp.msperframe = 1000/out;
    end
    
    [out, p] = GetValsIfPossible([C.BKGNDRGBCD C.BKGNDRGBCD], p, 'double');
    if ~isempty(out)
        disp('Got background rgb');
        gl.disp.bkgndrgb = out;
    end

    headerstruct = {...                                         %
	C.INITPROTOCOLCD, 'gl.protocol',1,'int';...                 
    C.NTRLSPERCONDCD, 'gl.ntrlspercond',1,'int';...                           
    C.INITDIAMCD,'gl.DEFAULTDIAM',100,'int';...
    C.INITSFCD,'gl.DEFAULTSF',100,'int';...
    C.INITTFCD,'gl.DEFAULTTF',100,'int';...                                
    C.INITORIENTCD,'gl.DEFAULTORIENT',180/pi,'int';...
    C.INITPHASECD,'gl.DEFAULTPHASE',100,'int';...
    C.INITLCONTCD,'gl.DEFAULTCONECONTRASTS(1)',100,'int';...
    C.INITMCONTCD,'gl.DEFAULTCONECONTRASTS(2)',100,'int';...
    C.INITSCONTCD,'gl.DEFAULTCONECONTRASTS(3)',100,'int';...
    C.INITNCYCLESCD,'gl.DEFAULTNCYCLES',1,'int';...
    C.NCYCLESPLATCD,'gl.NCYCLESPLATCD',1,'float';...
    C.NCYCLESRAMPCD,'gl.NCYCLESRAMPCD',1,'float';...
    C.NCONTRASTSCD,'gl.NCONTRASTS',1,'int';...
    C.NGABORREPSCD,'gl.NGABORREPS',1,'int';...
    C.MAXCONTSCALEFACTCD,'gl.MAXCONTSCALEFACTOR',1,'float';...
    C.MINCONTSCALEFACTCD,'gl.MINCONTSCALEFACTOR',1,'float';...
    C.STIMDURCD,'gl.STIMDUR',100,'int';...
    C.OPTOFLAGCD,'gl.OPTOFLAG',1,'int'};
    
    for i = 1:size(headerstruct,1)
        [out, p] = GetValsIfPossible(headerstruct{i,1}, p, headerstruct{i,4});
        if ~isempty(out)
            eval([headerstruct{i,2},' = out/',num2str(headerstruct{i,3}),';']);
        end
    end

    % Checking to see if we have the entire header
    if ~isempty(gl.disp.bkgndrgb) && ~isempty(gl.disp.monSpd) && ~isempty(gl.MINCONTSCALEFACTOR) && ~isempty(gl.disp.fundamentals)
        gl.gotheader = 1; % For now assuming if we got this, we got everything
        SetUpTrials(gl.protocol);
    end
end

%% SetUpTrials

% Setting up the gl structure according to the protocol specified by REX. 
% gl.ntrlspercond = 0 means: "do an ANOVA to determine whether we should do another
% block of trials or (if p < 0.05) move to the next protocol".  If gl.ntrlspercond == 0 
% then set gl.trialspecs{:,1} = 1 so that we do one trial in each condition
% before doing the ANOVA again.

% YES: NOT USING ANOVA CODE

function SetUpTrials(protocol)
    global gl;
    
    a = get(gcf,'UserData');
    TS = TrialSpecCodes;
    
    % Defining "constants"
    gl.gratingparamcols = find(~ismember(fieldnames(TS),{'TRIALCOUNT','SPIKES','BASELINE'}));  % Columns of gl.trialspecs containing the grating parameters
    gl.trialspecs = {};
    gl.protocol = protocol;
    
    if gl.protocol==1       %DEFAULT PROTOCOL IN GRATINGS.D IS 1, SO SWITCHING TO 11
        protocol = 11;
        gl.protocol =11;
    else
    end
    
    switch (protocol)
        
        case 11, %ORI COARSE
            orients = [0:pi/4:7*pi/4];                                                 
            ntrialtypes = length(orients);
        
        case 12, %SF
            %sfs = [0.3 0.5004 0.8348 1.3925 2.3228 3.8746 6.4633];              
            sfs = [0.5004 0.8348 1.3925 2.3228 3.8746 6.4633];
            ntrialtypes = length(sfs);
        
        case 13, %TF
            %tfs = [0.3 0.5004 0.8348 1.3925 2.3228 3.8746 6.4633 10.7814 17.9845 30];    %YES: TFS FROM MOVSHON LAB (GOOD FOR V1 AND MT)
            tfs = [0.5004 0.8348 1.3925 2.3228 3.8746 6.4633 10.7814 17.9845 30];    %YES: TFS FROM MOVSHON LAB (GOOD FOR V1 AND MT)
            ntrialtypes = length(tfs);                                                  
        
        case 14, %SIZE
            %apsizes = linspace(0,2,9);                                                  %YES: CODE GETS STUCK AT SIZE 0!!!                                        
            apsizes =linspace(0.1,2.5,9);                                                       %YES: NEW RANGE BASED ON MY ANTIDROMIC DATA
            ntrialtypes = length(apsizes);
        
        case 15, %ORI FINE                      %YES: REVERTED BACK TO COARSE ON 1/2/13
            %orients = [0:pi/8:15*pi/8];
            orients = [0:pi/4:7*pi/4];  
            ntrialtypes = length(orients);
            
        case 16, %COLOR
            colordirections = [
                0.0900    0.0900         0
                0.0900   -0.0900         0
                0         0         0.6364
                0.1273    0              0
                0         0.1273         0
                0.0636   -0.0636   -0.4500
                0.0636   -0.0636    0.4500
                0.0636    0.0636    0.4500
                0.0636    0.0636   -0.4500];    
            colordirections = colordirections*.9;                                       % Ugly hack to keep things in gamut on rig 1
            ntrialtypes = size(colordirections,1);                                      
        
        case 17, %CONTRAST RESPONSE FUNCTION                                                               
            %colordirections = [0.0000 0.0078 0.0156 0.0312 0.0624 0.1249 0.2499 0.4999 1]'*[1 1 1]; %YES: CONTRASTS FROM MOVSHON LAB    
            colordirections = [0.0000 0.0078 0.0156 0.0312 0.0624 0.1249 0.2499 0.4999 0.74]'*[1 1 1]; %YES: 74% max to keep things in gamut on rig 1
            
            if gl.OPTOFLAG 
                ntrialtypes = 2*size(colordirections,1);    %YES: OPTICAL STIM ON HALF OF TRIALS
            else
                ntrialtypes = size(colordirections,1);      %YES: NO OPTICAL STIM
            end

        otherwise
            ntrialtypes = 0;
    end

    
    % Setting up default values
    for i = 1:ntrialtypes
        gl.trialspecs{i,TS.TRIALCOUNT} = max([gl.ntrlspercond 1]);              % No fewer than 1 trial/cond
        gl.trialspecs{i,TS.SPIKES} = [];                                        % spike counts
        if (~isnan(gl.prefdiam))
            gl.trialspecs{i,TS.DIAM} = gl.prefdiam;                             % diam (deg)
        else
            gl.trialspecs{i,TS.DIAM} = gl.DEFAULTDIAM;                  
        end
        if (~isnan(gl.prefsf))
            gl.trialspecs{i,TS.SF} = gl.prefsf;                                 % SF (cycles/deg)
        else
            gl.trialspecs{i,TS.SF} = gl.DEFAULTSF;                      
        end
        if (~isnan(gl.preftf))
            gl.trialspecs{i,TS.TF} = gl.preftf;                                 % TF (cycles/sec);
        else
            gl.trialspecs{i,TS.TF} = gl.DEFAULTTF;                      
        end                         
        if ~isempty(gl.prefccdir) & any(gl.prefccdir)
            gl.trialspecs{i,TS.CC} = gl.prefccdir;                              % cone contrasts            
        else
            gl.trialspecs{i,TS.CC} = gl.DEFAULTCONECONTRASTS;           
        end
        if (~isnan(gl.preforient))
            gl.trialspecs{i,TS.ORIENT} = gl.preforient;                         % orientation (rad)
        else
            gl.trialspecs{i,TS.ORIENT} = gl.DEFAULTORIENT;              
        end
        gl.trialspecs{i,TS.PHI} = gl.DEFAULTPHASE;                              % phase (rad)               
        gl.trialspecs{i,TS.STIMTYPE} = 0;                                       % 0 = grating, 1 = gabor,2 = grating+optical stim
        gl.trialspecs{i,TS.PROTOCOL} = protocol;                                % which protocol
    end
    
    
    if (protocol == 11)                                                          %YES: ORI TUNING COARSE
        disp('Setting up trial parameters for orientation tuning'); 
        for i = 1:ntrialtypes
            gl.trialspecs{i,TS.ORIENT} = orients(i);                            
        end
        
    elseif (protocol == 12)                                                      %YES: SF TUNING
        disp('Setting up trial parameters for SF tuning');
        for i = 1:ntrialtypes
            gl.trialspecs{i,TS.SF} = sfs(i);                                   
        end
        
    elseif (protocol == 13)                                                      %YES: TF TUNING
        disp('Setting up trial parameters for TF tuning');
        for i = 1:ntrialtypes
            gl.trialspecs{i,TS.TF} = tfs(i);                                    
        end
        
    elseif (protocol == 14)                                                      %YES: SIZE TUNING
        disp('Setting up trial parameters for aperture size tuning');
        for i = 1:ntrialtypes
            gl.trialspecs{i,TS.DIAM} = apsizes(i);                              
        end
   
    elseif (protocol == 15)                                                      %YES: ORI TUNING FINE
        disp('Setting up trial parameters for orientation tuning'); 
        for i = 1:ntrialtypes
            gl.trialspecs{i,TS.ORIENT} = orients(i);                            
        end
        
    elseif (protocol == 16)                                                      %YES: COLOR TUNING
        disp('Setting up trial parameters for color tuning');
        for i = 1:ntrialtypes
            gl.trialspecs{i,TS.CC} = colordirections(i,:);                     
        end
        
    elseif (protocol == 17)                                                      %YES: CONTRAST RESPONSE FUNCTION
        disp('Setting up trial parameters for contrast-response function');
        i=1;
        if gl.OPTOFLAG                                                           %YES: OPTOSTIM ON HALF OF TRIALS
            ks = [0 2];
            ntrialtypes = ntrialtypes/2;                                         %YUCK! Change this
        else
            ks = 0;
        end
        for j = 1:ntrialtypes
            for k =ks                                                       
                gl.trialspecs{i,TS.CC} = colordirections(j,:);
                gl.trialspecs{i,TS.STIMTYPE} = k;                               
                i = i+1;
            end
        end
        %set(a.uicontrols.reset,'Enable','off');                                 
        
    elseif (protocol == 0)
        disp('Early termination request received');
    end
    
    % Setting up ncycles to show
    for i = 1:size(gl.trialspecs,1)
        gl.trialspecs{i,TS.NCYCLES} = gl.trialspecs{i,TS.TF}.* gl.STIMDUR ;       %YES: COMPUTING NCYCLES FROM STIMDUR  
    end
        
end

%% protocol [CHANGED PROTOCOL DETAILS]

% Called when we've completed a protocol, this function figures out which
% protocol to run next.

function protocol = PickNewProtocol
    global gl;
    uicontrols = getfield(get(gcf,'UserData'),'uicontrols');

    %YES: PROTOCOL#11 IS ORI COARSE
    if (gl.protocol == 11)
        preforient = gl.splineargmax;
        if (isnan(gl.preforient))
            gl.preforient = preforient;
        else
        end
        disp(['new preferred orientation is ',num2str(preforient*180/pi)]);
        gl.preforient = preforient;
        if (isnan(gl.prefsf))
            protocol = 12;
        else
            protocol = 16;
        end

        %YES: PROTOCOL#12 IS SF
    elseif (gl.protocol == 12)
        prefsf = gl.splineargmax;
        if (isnan(gl.prefsf))
            gl.prefsf = prefsf;
        else
        end
        disp(['new preferred sf is ',num2str(prefsf)]);
        gl.prefsf = prefsf;
        protocol = 13;

        %YES: PROTOCOL#13 IS TF
    elseif (gl.protocol == 13)
        preftf = gl.splineargmax;
        if (isnan(gl.preftf))
            gl.preftf = preftf;
        else
        end
        disp(['new preferred tf is ',num2str(preftf)]);
        gl.preftf = preftf;
        protocol = 14;

        %YES: PROTOCOL#14 IS SIZE
    elseif (gl.protocol == 14)
        gl.prefdiam = gl.splineargmax;
        disp(['new preferred aperture size is ',num2str(gl.prefdiam)]);
        protocol = 15;

        %YES: PROTOCOL#15 IS ORI FINE
    elseif (gl.protocol == 15)
        gl.preforient = gl.splineargmax;
        disp(['new preferred orientation is ',num2str(gl.preforient*180/pi)]);
        protocol = 16;

        %YES: PROTOCOL#16 IS COLOR
    elseif (gl.protocol == 16)
        protocol = 17;

        %YES: PROTOCOL#17 IS CONTRAST WITH/OUT OPTICAL STIM
    elseif (gl.protocol == 17)
        protocol = 0;

        %YES: PROTOCOL#0 IS ALL DONE
    elseif (gl.protocol == 0)
        protocol = 0;
    end

    disp(['Switching to protocol ',num2str(protocol)]);
    gl.splineargmax = [];
    end

%% PickGratingParams

% Loading the appropriate parameters into gl.gratingparams
% Order of grating parameters:
% [diam, sf, tf, lcont, mcont, scont, orient, phase, ncycles, stimtype, protocol].
% See bottom of code.
%
% This function gets called whenever REX sends a message
% requesting gl.gratingparams.

function PickGratingParams()
    global gl;
    global udpCom;
    TS = TrialSpecCodes;
    
    if (~gl.gotheader)
        keyboard
        error('We have not yet received the header!');                                                            
    end

    if (all([gl.trialspecs{:,TS.TRIALCOUNT}] == 0) && gl.ntrlspercond == 0)
        % If we're out of trials and REX has specified that we should do an ANOVA
        % to figure out whether it's time to advance to the next protocol.
        p = anova1([gl.trialspecs{:,TS.SPIKES}],[],'off');
        if (isnan(p) || p > 0.05) % Do another block 
            [gl.trialspecs{:,TS.TRIALCOUNT}] = deal(1);
        end
    end
    if (isempty(gl.trialspecs) || all([gl.trialspecs{:,TS.TRIALCOUNT}] == 0))
        % Out of trials either because the ANOVA was significant or because
        % we're using a fixed number of trials (as opposed to ANOVA).
        protocol = PickNewProtocol;
        if (protocol == 0)   % we're done, send preferred parameters
            gl.gratingparams = [gl.prefdiam; gl.prefsf; gl.preftf; 0; 0; 0; gl.preforient; 0; 0; 0];                                               
            gl.gratingparams(isnan(gl.gratingparams)) = -1;                                             %YES: BETTER TO SEND -1 (WAS 0)? DOES THIS WORK?
            return;
        else
            SetUpTrials(protocol);
        end
    end

    % Picking a randomly chosen trial
    ntrialsleft = [gl.trialspecs{:,TS.TRIALCOUNT}]';       
    choices = find(ntrialsleft == max(ntrialsleft));                                                                    
    gl.paramsidx = choices(unidrnd(length(choices)));
    gl.gratingparams = [gl.trialspecs{gl.paramsidx, gl.gratingparamcols}];
    gl.trialspecs{gl.paramsidx,TS.TRIALCOUNT} = gl.trialspecs{gl.paramsidx,TS.TRIALCOUNT}-1;
    
    %alert Rex that the next trial's params have been chosen
    pnet(udpCom.sock, 'write', 'paramSetupComplete>> >>');
    pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING AND ANALYSIS FUNCTIONS BELOW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SetUpFigure [YES: ADDED SOME THINGS]

function SetUpFigure()
   
    figure(1);
    a = get(gcf,'UserData');                    % In case there's already things there (eg from whitenoise)
    set(gcf,'position',[300 200 600 400])
    set(gcf,'DefaultAxesUnits','pixels')
    clf;

    h1 = axes('position',[50 50 300 200]);
    ylabel('sp/sec');

    h2 = axes('position',[50 300 300 50]);
    set(h2,'Ylim',[-1 1],'Xlim',[0 1000],'NextPlot','add');
    xlabel('Time (ms)');

    for i = 1:3
        h3(i) = axes('position',[50+(i-1)*120 100 80 80],'Visible','off',...
            'XTick',[],'YTick',[],'XColor',[1 1 1],'YColor',[1 1 1]);
    end

    % UICONTROLS
    uicontrol('style','text','Position',[380 40 86 20],'String','Trial countdown:');
    a.uicontrols.ntrialstextbox = uicontrol('style','text','Position',[460 40 30 20],'String',num2str(0));
    uicontrol('style','text','Position',[380 65 86 20],'String','Protocol number:');
    a.uicontrols.protocoltextbox = uicontrol('style','text','Position',[460 65 30 20],'String',num2str(0));
    a.uicontrols.changecol = uicontrol('style','popupmenu','string','Achrom|LvsM|S|L-M+S|L-M-S','Callback',@ChangeCol,'Position',[370 330 80 20]);
    a.uicontrols.sftextbox = uicontrol('style','text','Position',[460 350 80 12],'String','SF');
    a.uicontrols.changesf = uicontrol('style','slider','Min',0.5,'Max',6,'Callback',@ChangeSF,'SliderStep',[.01 .1],'Value',2,'Position',[460 330 80 20]);
    a.uicontrols.tftextbox = uicontrol('style','text','Position',[460 300 80 12],'String','TF');
    a.uicontrols.changetf = uicontrol('style','slider','Min',0.5,'Max',18,'Callback',@ChangeTF,'SliderStep',[.01 .1],'Value',2,'Position',[460 280 80 20]);
    a.uicontrols.diamtextbox = uicontrol('style','text','Position',[460 250 80 12],'String','Diam');
    a.uicontrols.changediam = uicontrol('style','slider','Min',0.1,'Max',3,'Callback',@ChangeDiam,'SliderStep',[.01 .1],'Value',1,'Position',[460 230 80 20]);
    a.uicontrols.reset = uicontrol('style','pushbutton','String','Reset','Position',[370 350 80 20],'Callback',@ResetCallback);
    
%     uicontrol('style','text','Position',[500 65 95 20],'String','Optical Stim?');
%     a.uicontrols.optoradio = uicontrol('style','radio','Position',[535 51 15 10],'Value',1);
    
    a.axeshandles.tuningcurve = h1;
    a.axeshandles.raster = h2;
    a.axeshandles.colorspace = h3;
    set(gcf,'UserData',a);

    % ChangeCol
        function ChangeCol(ev,h)
            global gl
            TS = TrialSpecCodes;

            uicontrols = getfield(get(gcf,'UserData'),'uicontrols');
            whichcolidx = get(uicontrols.changecol,'Value');
            whichcolnames = get(uicontrols.changecol,'String');
            whichcolname = whichcolnames(whichcolidx,:);
            switch deblank(whichcolname)
                case 'Achrom'
                    cc = [.6 .6 .6];
                case 'LvsM'
                    cc = [.09 -.09 0];
                case 'S'
                    cc = [0 0 .6364];
                case 'L-M-S'
                    cc = [0.0636 -0.0636 -0.45];
                case 'L-M+S'
                    cc = [0.0636 -0.0636 0.45];
                otherwise
                    disp('Unknown color choice');
                    keyboard
            end
            cc = cc*.9;
            for i = 1:size(gl.trialspecs,1)
                gl.trialspecs{i,TS.CC} = cc;
            end
            gl.DEFAULTCONECONTRASTS = cc;
            gl.prefccdir = cc;
        end

    %ChangeSF
        function ChangeSF(ev,h)
            global gl
            TS = TrialSpecCodes;

            uicontrols = getfield(get(gcf,'UserData'),'uicontrols');
            sf = get(uicontrols.changesf,'Value');
            for i = 1:size(gl.trialspecs,1)
                gl.trialspecs{i,TS.SF} = sf;
            end
            gl.DEFAULTSF = sf;
            gl.prefsf = sf;
        end

    %ChangeTF                                                     
        function ChangeTF(ev,h)
            global gl
            TS = TrialSpecCodes;

            uicontrols = getfield(get(gcf,'UserData'),'uicontrols');
            tf = get(uicontrols.changetf,'Value');
            for i = 1:size(gl.trialspecs,1)
                gl.trialspecs{i,TS.TF} = tf;
                gl.trialspecs{i,TS.NCYCLES} = tf.* gl.STIMDUR; 
            end
            gl.DEFAULTTF = tf;
            gl.preftf = tf;
        end

    %ChangeDiam
        function ChangeDiam(ev,h)
            global gl
            TS = TrialSpecCodes;

            uicontrols = getfield(get(gcf,'UserData'),'uicontrols');
            diam = get(uicontrols.changediam,'Value');
            for i = 1:size(gl.trialspecs,1)
                gl.trialspecs{i,TS.DIAM} = diam;
            end
            gl.DEFAULTDIAM = diam;
            gl.prefdiam = diam;
        end

    %ResetCallback
        function ResetCallback(ev,h)
            global gl

            gl.preforient = nan;
            gl.prefsf = nan;
            gl.preftf = nan;
            gl.prefdiam = nan;
            gl.prefccdir = [0 0 0];
            SetUpTrials(1);
        end
end

%% PlotTuningCurve

% Plots a tuning curve with a spline fit
% and updates gl.splineargmax (for efficiency).
% Replotting every single point each time we come through here (not adding to the plot).

function PlotTuningCurve()
    global gl;
    TS = TrialSpecCodes;
    
    MS = 15; %MarkerSize
    LW = 1;  %LineWidth
    spikerates = [];
    baselines = [];
    orientations  = [];
    sfs  = [];
    tfs = [];                                                                      
    diams = [];                                                              
    fit = [];                                                                      
    x = [];
    
   for i = 1:size(gl.trialspecs,1)
        spikerate = gl.trialspecs{i,TS.SPIKES};
        spikerates = [spikerates; spikerate];
        orientations = [orientations; repmat(gl.trialspecs{i,TS.ORIENT}, length(spikerate),1)];
        sfs = [sfs; repmat(gl.trialspecs{i,TS.SF}, length(spikerate),1)];
        tfs = [tfs; repmat(gl.trialspecs{i,TS.TF}, length(spikerate),1)];          
        diams = [diams; repmat(gl.trialspecs{i,TS.DIAM}, length(spikerate),1)];
        baselines = [baselines; gl.trialspecs{i,TS.BASELINE}];
    end
    meanbaseline = mean(baselines);
    
    

    
    
    
    allorientations = unique([gl.trialspecs{:,TS.ORIENT}]);
    allsfs = unique([gl.trialspecs{:,TS.SF}]);
    alltfs = unique([gl.trialspecs{:,TS.TF}]);                                      
    alldiams = unique([gl.trialspecs{:,TS.DIAM}]);
    
    % Setting up axes
    axesh = getfield(get(gcf,'UserData'),'axeshandles');
    set(axesh.tuningcurve,'Visible','on');
    set(axesh.colorspace,'Visible','off');
    for i = 1:length(axesh.colorspace)
        delete(get(axesh.colorspace(i),'Children'));
    end
    axes(axesh.tuningcurve);
    h = get(axesh.tuningcurve,'Children');
    delete(h);
    hold on;
    box off;
    
    if (length(allorientations) > 2)
        set(gca,'XTick',[0:90/4:360],'XTickLabel',{'0','','45','','90','','135','','180','','225','','270','','315','','360'},'XLim',[0 360],'XScale','linear'); 
        xlabel('orientation (deg)');
        parametervals = orientations;
        L = parametervals == min(parametervals);
        parametervals = [parametervals; parametervals(L)+2*pi];  % computing in rads, plotting in deg
        spikerates = [spikerates; spikerates(L)];
        splineopt = 'periodic';
        plotcoef = 180/pi;
           
    elseif (length(allsfs) > 2)
        set(gca,'XTick',[.1 1 10],'xticklabel',[.1 1 10],'XLim',[0.1 10],'XScale','log');
        xlabel('sf (cycles/deg)');
        parametervals = sfs;
        splineopt = 'variational';
        plotcoef = 1;
        
    elseif (length(alltfs) > 2)
        set(gca,'XTick',[.1 1 10],'xticklabel',[.1 1 10],'XLim',[0.1 40],'XScale','log');
        xlabel('tf (cycles/sec)');
        parametervals = tfs;
        splineopt = 'variational';
        plotcoef = 1;
        
    elseif (length(alldiams) > 2)
        set(gca,'XTick',[.01 .1 1 3],'xticklabel',[.01 .1 1 3],'XLim',[0.05 3],'XScale','log');
        xlabel('diameter (deg)');
        parametervals = diams;
        splineopt = 'not-a-knot';
        plotcoef = 1;
        
    else
        return
    end

    %YES: SPLINE FITTING HAPPENS HERE
    if (length(unique(parametervals)) > 2)  % Compute the spline
        x = linspace(min(parametervals), max(parametervals),100);
        pp = csape(parametervals,spikerates,splineopt);
        fit = ppval(pp,x);
        gl.splineargmax = min(x(fit == max(fit))); % min to avoid ties
        plot(plotcoef*x,fit,'k-','LineWidth',LW);
        plot(plotcoef*gl.splineargmax,max(fit),'r*','MarkerSize',MS-5);
        plot(plotcoef*[min(parametervals) max(parametervals)], [meanbaseline meanbaseline], 'k:');
    end
    
    plot(plotcoef*parametervals, spikerates,'k.','MarkerSize',MS);
    
%     maxRho=[];
%     if spikerates ~=0
%         maxRho = max([max(spikerates) max(fit)])
%         set(gca,'ylim',[0 ceil((maxRho*10)/100)*10]);
%     else
%     end
    
    hold off;   
end

%% PlotColorTuning

function PlotColorTuning
    global gl
    TS = TrialSpecCodes;

    % Setting up axes
    axesh = getfield(get(gcf,'UserData'),'axeshandles');
    MS = 15; %MarkerSize
    spikerates = [];
    baselines = [];
    ccs = [];
    for i = 1:size(gl.trialspecs,1)
        if (gl.trialspecs{i,TS.PROTOCOL} == 16)                                  
            spikerate = gl.trialspecs{i,TS.SPIKES};
            spikerates = [spikerates; spikerate];
            ccs = [ccs; repmat(gl.trialspecs{i,TS.CC}, length(spikerate),1)];
            baselines = [baselines; gl.trialspecs{i,TS.BASELINE}];
        end
    end
    meanbaseline = mean(baselines);
 
    set(axesh.tuningcurve,'Visible','off');
    set(axesh.colorspace,'Visible','on');
    delete(get(axesh.tuningcurve,'Children'));
    for i = 1:length(axesh.colorspace)
        delete(get(axesh.colorspace(i),'Children'));
    end

    % Plotting
    basis = [1 1 0; 1 -1 0; 0 0 1];
    DKLs = (basis*ccs')';
    coordinates = sign(DKLs);
    peakfr = max([spikerates; .1]); % .1 in case all spikerates = 0
    labels = {'L+M','L-M';'S','L-M';'S','L+M'};
    permmat = [3 1 2; 1 3 2; 2 3 1];
    for i = 1:3
        axes(axesh.colorspace(i));
        set(gca,'XLim',[-peakfr peakfr],'YLim',[-peakfr peakfr]);
        hold on;
        text(0,1.5*peakfr,labels{i,1},'HorizontalAlignment','center');
        text(1.1*peakfr,.2*peakfr,labels{i,2});
        L = coordinates(:,permmat(i,1)) == 0;
        theta = atan2(coordinates(:,permmat(i,2)),coordinates(:,permmat(i,3)));
        [x,y] = pol2cart(theta(L), spikerates(L));
        plot(x, y,'k.','MarkerSize',MS)
        plot(-x,-y,'k.','MarkerSize',MS)
    end
end

%% computePrefColors [CHANGED PROTOCOL#]

% Compute the prefered cardinal and intermediate colors
% Also the color direction (in or out of isoluminant plane) that 
% evoked the greatest response.

function [prefisolumcolors, prefccdir] = computePrefColors()
    global gl
    TS = TrialSpecCodes;

    spikerates = [];
    baselines = [];
    ccs = [];
    for i = 1:size(gl.trialspecs,1)
        if (gl.trialspecs{i,TS.PROTOCOL} == 16)                                     
            spikerate = gl.trialspecs{i,TS.SPIKES};
            spikerates = [spikerates; spikerate];
            ccs = [ccs; repmat(gl.trialspecs{i,TS.CC}, length(spikerate),1)];
            baselines = [baselines; gl.trialspecs{i,TS.BASELINE}];
        end
    end
 %   meanbaseline = mean(baselines);
    prefccdir = [];
    
    uniqueccs = unique(ccs,'rows');
    for a = 1:size(uniqueccs,1)
        tList = ismember(ccs, uniqueccs(a,:), 'rows');
        allRates(a,1) = mean(spikerates(tList)); % -baselines(tList)); % Commented out for now
    end
    prefccdir = uniqueccs(find(allRates==max(allRates),1),:);
    
    % Now getting prefisolumdir
    cardinals = [1 -1 0; 0 0 1];
    intermediates = [1 -1 -1; 1 -1 1];
    trialColors = sign(ccs);
    for a = 1:size(cardinals, 1);
        [tList, dummy] = ismember(trialColors, cardinals(a,:), 'rows');
        cardinalRates(a,1) = mean(spikerates(tList)-baselines(tList));
        cardinalCCS(a,:) = unique(ccs(tList,:), 'rows');
    end
    for a = 1:size(intermediates, 1);
        [tList, dummy] = ismember(trialColors, intermediates(a,:), 'rows');
        intermediateRates(a,1) = mean(spikerates(tList)-baselines(tList));
        intermediateCCS(a,:) = unique(ccs(tList,:),'rows');
    end
    
    %determine the perfered cardinal and intermediate colors. Send them to
    %rex as a vector where the first triplet is the "true" prefered color
    %(on the basis of maximum firing rate)
    [val, cardInd] = max(cardinalRates);
    [val, intInd] = max(intermediateRates);
    if cardinalRates(cardInd) > intermediateRates(intInd)
        prefisolumcolors(1,1:3) = cardinalCCS(cardInd,:);
        prefisolumcolors(1,4:6) = intermediateCCS(intInd,:);
    else
        prefisolumcolors(1,1:3) = intermediateCCS(intInd,:);
        prefisolumcolors(1,4:6) = cardinalCCS(cardInd,:);
    end
    disp(['prefcolor:', num2str(prefccdir)])
end

%% PlotContrastResponse

function PlotContrastResponse
    global gl
    TS = TrialSpecCodes;

    % Setting up axes
    axesh = getfield(get(gcf,'UserData'),'axeshandles');
    MS = 15; %MarkerSize
    LW = 1;  %LineWidth
    spikerates = [];
    baselines = [];
    ccs = [];
    ks = [];
    for i = 1:size(gl.trialspecs,1)
        if (gl.trialspecs{i,TS.PROTOCOL} == 17)
            spikerate = gl.trialspecs{i,TS.SPIKES};
            spikerates = [spikerates; spikerate];
            ccs = [ccs; repmat(gl.trialspecs{i,TS.CC}, length(spikerate),1)];
            ks = [ks; repmat(gl.trialspecs{i,TS.STIMTYPE}, length(spikerate),1)]; %YES: ADDED
            baselines = [baselines; gl.trialspecs{i,TS.BASELINE}];
        end
    end
    meanbaseline = mean(baselines);

    allccs = unique(ccs,'rows');   
    allccs = allccs(:,1);
    allks = unique(ks);           

    set(axesh.tuningcurve,'Visible','on')
    delete(get(axesh.tuningcurve,'Children'));
    set(axesh.colorspace,'Visible','off');
    for i = 1:length(axesh.colorspace)
        delete(get(axesh.colorspace(i),'Children'));
    end
    axes(axesh.tuningcurve); hold on;
    set(gca,'XTick',[.001 .01 0.1 1],'xticklabel',[0.001 0.01 0.1 1]*100,'XLim',[0.005 1.05],'XScale','log');
    xlabel('contrast (%)');
   
    % PARSE BY OSTIM CONDITION
    ind_ostim_off = find(ks==0);
    ind_ostim_on = find(ks==2);
    
    splineopt = 'not-a-knot';
    fit =[];
    x =[];
    fit1 =[];
    fit2 = [];
    plotcoef = 1;
    parametervals = ccs(:,1);
    plot(plotcoef*[.01 max(parametervals)], [meanbaseline meanbaseline], 'k:'); hold on;
    
    if length(allks)==1
        plot(ccs(ind_ostim_off,1), spikerates(ind_ostim_off),'k.','MarkerSize',MS);
        if (length(allccs) > 2)  % Compute the spline
            x = linspace(min(parametervals), max(parametervals),100);
            pp = csape(parametervals,spikerates,splineopt);
            fit = ppval(pp,x);
            gl.splineargmax = min(x(fit == max(fit))); % min to avoid ties
            plot(plotcoef*x,fit,'k-','LineWidth',LW); hold on;
            fit1 = fit;
        end
        
    elseif length(allks)==2
        plot(ccs(ind_ostim_off,1), spikerates(ind_ostim_off),'k.','MarkerSize',MS);
        plot(ccs(ind_ostim_on,1), spikerates(ind_ostim_on),'r.','MarkerSize',MS);
        if (length(ccs(ind_ostim_off,1)) > 2)  % Compute the spline
            x = linspace(min(parametervals(ind_ostim_off)), max(parametervals(ind_ostim_off)),100);
            pp = csape(parametervals(ind_ostim_off),spikerates(ind_ostim_off),splineopt);
            fit = ppval(pp,x);
            gl.splineargmax = min(x(fit == max(fit)));
            plot(plotcoef*x,fit,'k-','LineWidth',LW); hold on;
            fit1 = fit;
        end
        
        if (length(ccs(ind_ostim_on,1)) > 2)  % Compute the spline
            x = linspace(min(parametervals(ind_ostim_on)), max(parametervals(ind_ostim_on)),100);
            pp = csape(parametervals(ind_ostim_on),spikerates(ind_ostim_on),splineopt);
            fit = ppval(pp,x);
            gl.splineargmax = min(x(fit == max(fit))); 
            plot(plotcoef*x,fit,'r-','LineWidth',LW); hold on;
            fit2 = fit;
        end
    end
    
%     if spikerates ~=0
%         maxRho = max([max(spikerates) max([fit1 fit2])]);
%         set(gca,'ylim',[0 ceil((maxRho*10)/100)*10]);
%     else
%     end
    hold off;

end

%% PlotRaster [UNCHANGED]

function PlotRaster(spiketally)
    axesh = getfield(get(gcf,'UserData'),'axeshandles');
    axes(axesh.raster);
    plot([spiketally, spiketally]',[-.5; .5]*ones(1,length(spiketally)),'k-');
    set(gca,'ytick',[],'yticklabel',[]);
end

%% ClearRaster [UNCHANGED]

function ClearRaster
    axesh = getfield(get(gcf,'UserData'),'axeshandles');
    h = get(axesh.raster,'Children');  % Should be able to do this with cla but can't
    delete(h);
end

%% UpdateTrialCountdown

function UpdateTrialCountdown
    global gl;
    uicontrols = getfield(get(gcf,'UserData'),'uicontrols');
    TS = TrialSpecCodes;
    ntrials = sum([gl.trialspecs{:,TS.TRIALCOUNT}]);
    set(uicontrols.ntrialstextbox,'String',num2str(ntrials));
    set(uicontrols.protocoltextbox,'String',num2str(gl.trialspecs{gl.paramsidx,TS.PROTOCOL}));
end

% Here are the columns of the gl.trialspecs array:
% 1: Number of trials in this condition more to do
% 2: Driven spike counts
% 3: Baseline spike counts
% 4: Diameter
% 5: SF
% 6: TF
% 7: Cone contrasts
% 8: Orientation
% 9: Spatial phase
% 10: Number of cycles
% 11: Protocol number

%% UpdateTrialCountdown

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