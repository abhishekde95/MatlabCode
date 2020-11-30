function p = WhiteNoiseOnlinethresh_AD()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Online spike triggered averaging and covariance
% for use with the WhiteNoise.d and WhiteNoise.m.
% GDLH 12/9/07
%
% Upgrade to work with multiple spike channels and
% binary cone noise.
% GDLH 11/3/08
% 
% added the closed loop online features to the existing code - Abhishek
% 04/16
% Searching along equispaced angles - the only adaptive algo is the
% staircase - Abhishek 12/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Global Variables
global subunitcmap
subunitcmap = [0.89412 0.10196 0.1098;0.21569 0.49412 0.72157;0.30196 0.68627 0.2902;0.59608 0.30588 0.63922;1 0.49804 0;1 1 0.2;0.65098 0.33725 0.15686;0.96863 0.50588 0.74902;0.6 0.6 0.6];
global gl % 1/16 - Abhishek created a new global structure
global udpCom;  % Needed by "DealWithM1essage()" and associated subfunctions
udpCom.sock = [];
udpCom.port = 6665;
udpCom.rexip = '192.168.1.120';

C = WhiteNoisethreshCodes;        % Script that defines a bunch of constants
DEFAULTNFRAMES = 10;  % Number of frames to go back before each spike
% Some initializations
p = InitPStruct(C.EOTCD, C.FIXXCD); % Structure of the events from each trial
s = InitPlex();     % Link to Plexon datastream
InitStatsStruct();  % Structure that holds the statistics
InitGLStruct(); % Global data structure, Abhishek - 1/16
SetUpFig(DEFAULTNFRAMES);     % Setting up the figure.  I'm going to try to keep data in the 'UserData' field.
[udpCom.sock, success] = pnetStart(udpCom.port);  % UDP communication with REX
if ~success, return; end
stopnow = 0;    % Flag to break out of the endless while loop (hit ESC to set to 1)
codestruct = {...
    C.FPONCD, @fponfn;...
    C.STIMONCD, @stimonfn;...
    C.FPACQCD, @fpacqfn;...
    C.ABORTCD, @abortfn;...
    C.CORRECTCD, @correctfn;...
    C.EOTCD, @eotfn;...
    };
fpon_t = [];
fpacq_t = [];  % code dropped at the prestimdel (STATE in REX)
stimon_t = [];
stimoff_t = [];
correct_t = [];

% The main loop
while (~stopnow)
    stopnow = CheckForESCKey();
    
    % Check for a message from REX
    msgSize = pnet(udpCom.sock, 'readpacket', 200, 'noblock');
    
    if(msgSize)
        stopnow = DealWithMessage(msgSize);
    end
    
    [n, eventList] = PL_GetTS(s);
    if (n > 0)
        p = ProcessEventList(p, eventList);
        UpdateQueueText(p);
    end
    if (~p.processnowflag)  % Set by ProcessEventList()
        continue;
    end
    
    % In case this is the first trial
    if any(p.events == p.headercode)
        p = RemoveOldEvents(p, p.headercode);
        GetHeader(p);
        init_subunit_axes();
        ClearStats; % Just fills up the structures in 'Userdata' with NaNs
        UpdatePixelMask(nan,nan);
    end
    
    % If we haven't received the header yet go back to the
    % beginning and reset the 'p' structure.
    stats = getfield(get(gcf,'UserData'),'stats');
    if (~stats.gotHeader)
        p = InitPStruct(C.EOTCD, C.FIXXCD);
        continue;
    end
    % By this point, we know we've got the header and p.events contains at
    % least 1 EOTCD.
    
    % gotfulltrial = 0 if there is an error, 1 - if successful trial
    a = get(gcf,'UserData');
    b = get(a.axeshandles.synthImage,'UserData');
    
    [seed, nframes, mu, sigma, bkgndrgb, noisetype, nt, gotfulltrial] = GetTrialParams(p);
    
    if isempty(b.weightstosend) && nt~=1 % only do this in WN trials
        if (gotfulltrial == 0)
            p = CleanUpEvents(p); % Remove a trial with EOT but without data (e.g. fixation breaks)
            continue;
        end
    end
   
    % By this time we know, the trial was a good trial 
    % Dealing with each code one by one - Abhishek 1/16
    % Do the calculations,. proceed to this place only if gotfulltrial = 1;
    
    for i=1:size(codestruct,1)
        if(any(p.events == codestruct{i,1}))
            feval(codestruct{i,2})
        end
    end
    p = CleanUpEvents(p);
    
end % while (~stopnow) , big while loop
    function fponfn
        fpon_t = p.times(find(p.events == C.FPONCD,1));
    end

    function stimonfn
        stimon_t = p.times(find(p.events == C.STIMONCD,1));
    end

    function fpacqfn
        fpacq_t = p.times(find(p.events == C.FPACQCD,1));
    end

    function abortfn
        disp('ABORT code issued');
        if ~isempty(gl.trialspecidx)
            % Enter the if statment if u have received an ABORT code and
            % you happen to in Neurothresh mode. 'norms' and 'stepsize' get
            % updated in 'calculateweightsrequest' whereas 'trialnum' and
            % 'spikerates' get updated in correctfn. So in case of ABORT
            % code, undo the changes you have to the 'norms' and
            % 'stepsize'.
                gl.trialspecs(gl.trialspecidx).norms(end) = [];
                gl.trialspecs(gl.trialspecidx).stepsize(end) = [];
                disp([gl.trialspecidx gl.trialspecs(gl.trialspecidx).weights]);
        end
    end

    function correctfn
        % Enter this function if u have acquired a CORRECT code. 
        % Note - CORRECT code only gets dropped in valid Neurothresh trials
        correct_t = p.times(find(p.events == C.CORRECTCD,1)); 
        stimoff_t = p.times(find(p.events == C.STIMOFFCD,1));
        a = get(gcf,'UserData');
        if get(a.uicontrols.NT,'Value') && ~isempty(stimoff_t)
            t_offset = str2num(get(a.uicontrols.latency,'String'));
            dur = stimoff_t-stimon_t - t_offset;
            spiketimes = p.spikes{get(a.uicontrols.whichspike,'value')}; % for the time being, am just concentrating on the first channel
            spiketally = [spiketimes(spiketimes <stimoff_t & spiketimes>(stimon_t+t_offset))];
            nspikes = numel(spiketally);
            spikerate = nspikes/(dur/1000);
            set(a.uicontrols.spikerate,'String',[num2str(spikerate)]);
            set(gcf,'UserData',a);
            update_NT_raster_plot(spiketally-stimon_t, t_offset);
            
            % for modifying the fields inside structure 'gl'. Below is the list of
            % the fields that are modified inside this function
            % spikerates
            % trial_num
            % reversal_num
            % done
            % measuredthreshold
            disp(['This was a correct trial']);
            disp([gl.trialspecidx gl.trialspecs(gl.trialspecidx).weights]);
            gl.trialspecs(gl.trialspecidx).trial_num = gl.trialspecs(gl.trialspecidx).trial_num + 1; % Initial value of trial_num field is 0
            gl.trialspecs(gl.trialspecidx).spikerates = [gl.trialspecs(gl.trialspecidx).spikerates; spikerate];
            if gl.trialspecs(gl.trialspecidx).trial_num >= 2
                thresh = gl.trialspecs(gl.trialspecidx).targetspikerate;
                rev = sign(gl.trialspecs(gl.trialspecidx).spikerates(end)-thresh)* sign(gl.trialspecs(gl.trialspecidx).spikerates(end-1)-thresh);
                if rev == -1
                    gl.trialspecs(gl.trialspecidx).reversal_num = gl.trialspecs(gl.trialspecidx).reversal_num + 1;
                    if gl.trialspecs(gl.trialspecidx).reversal_num == gl.initparams.nreversals
                        gl.trialspecs(gl.trialspecidx).done = 1;
                        gl.trialspecs(gl.trialspecidx).pending = 0;
                        gl.trialspecs(gl.trialspecidx).measuredthreshold = gl.trialspecs(gl.trialspecidx).norms(end);
                    end
                end
            end
        end  
    end

    function eotfn
        % EOT gets dropped in both WhiteNoise and Neurothresh trials
        if nt==0
            % Enter this condition only if its in WhiteNoise mode 
            if (stats.gotHeader)
                EnableResetButton(0);
                plotnow = UpdateSTX(p, seed, nframes, mu, sigma, bkgndrgb, noisetype);
                if (plotnow)
                    PlotSTA(noisetype);
                    [n_subunits,~] = check_num_subunits();
                    if (~ n_subunits)
                        UpdatePixelMask(nan,nan);
                    end
                    drawnow;
                end
                EnableResetButton(1);
                a = get(gcf,'UserData');
                b = get(a.axeshandles.synthImage,'UserData');
                dur = stimon_t - fpacq_t;
                spiketimes = p.spikes{get(a.uicontrols.whichspike,'value')}; % for the time being, am just concentrating on the first channel
                spiketally = [spiketimes(spiketimes <stimon_t & spiketimes >fpacq_t)];
                nspikes = numel(spiketally);
                b.spiketimes = [b.spiketimes; {spiketimes}];
                b.baseline = [b.baseline ;nspikes/(dur/1000)];
                set(a.axeshandles.synthImage,'UserData',b);
                set(gcf,'UserData',a);
                update_baseline_hist(b.baseline);
            end
        end
    end
PL_Close(s);
end

function SetUpTrialSpecs(weights, targetspikerates, parentvertices, parentoog, predthresh, remainder)
% Adds a new trialspecs right before you start probing in a particular stimulus direction
if (nargin<6)
     remainder = -10;
end
if (nargin<5)
     predthresh = [];
end
if (nargin <4)
    parentoog = [];
end
if (nargin <3)
    parentvertices = [];
end
global gl
idx = length(gl.trialspecs) + 1;
gl.trialspecs(idx).trial_num = 0;
gl.trialspecs(idx).spikerates = [];
gl.trialspecs(idx).norms = [];
gl.trialspecs(idx).reversal_num = 0;
gl.trialspecs(idx).stepsize = [];
gl.trialspecs(idx).done = 0;  % either 0 or 1
gl.trialspecs(idx).active = 1; % either 0 or 1
gl.trialspecs(idx).measuredthreshold = [];
gl.trialspecs(idx).predictedthreshold = predthresh;
gl.trialspecs(idx).targetspikerate = targetspikerates;
gl.trialspecs(idx).parentvertices = parentvertices;
gl.trialspecs(idx).parentoutofgamut = parentoog;
gl.trialspecs(idx).outofgamut = 0;
gl.trialspecs(idx).weights = weights;  % contains the weights to be used by the basis vectors
gl.trialspecs(idx).pending = 1; % Introducing new fields to limit the number of searches in a round
gl.trialspecs(idx).remainder = remainder; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting the header from the ecodes, Looping 
function GetHeader(p)
C = WhiteNoisethreshCodes;
headerstruct = {...
    [C.GAMMATABLECD C.GAMMATABLECD],'double', @gammatablefn;...
    [C.MONSPDCD C.MONSPDCD],'double', @monspdfn;...
    [C.FUNDAMENTALSCD C.FUNDAMENTALSCD],'double',@fundamentalsfn;...
    C.NPIXPERSTIXCD,'int',@npixperstixfn;...
    C.NSTIXPERSIDECD,'int',@nstixpersidefn;...
    C.FRAMERATECD,'double',@nframeratefn;...
    C.GAUSSLOCUTCD,'int',@gausslocutfn;...
    C.GAUSSHICUTCD,'int',@gausshicutfn;...
    C.LINPREDTOLCD,'float',@linpredtolfn;...
    C.STEPSIZECD,'float',@stepsizefn;...
    C.SCALECD,'float',@scalefn;...
    C.NREVERSALSCD,'int',@nreversalsfn;...
    C.OOGSCALECD,'float',@oogscalefn;...
    };

a=get(gcf,'UserData');
global gl;
for i=1:size(headerstruct,1)
    out = GetVal(headerstruct{i,1}, p.events, headerstruct{i,2});
    if(~isempty(out))
        feval(headerstruct{i,3},out);
    end
end
P_device = SplineSpd(linspace(380,780,size(a.stats.monSpd,1))',a.stats.monSpd, (380:5:780)');
a.stats.M = a.stats.fundamentals'*P_device;
a.stats.invM = inv(a.stats.M);
a.stats.gotHeader = 1;
set(gcf,'UserData',a);

% Local nested functions
    function gammatablefn(out)
        disp('Got gammatable');
        gammaTable = codes2num(out);
        a.stats.gammaTable = reshape(gammaTable, length(gammaTable)/3,3);
        a.stats.invgammaTable = InvertGamma(a.stats.gammaTable,1);
    end

    function monspdfn(out)
        monSpd = codes2num(out);
        a.stats.monSpd = reshape(monSpd, length(monSpd)/3,3);
    end

    function fundamentalsfn(out)
        fundamentals = codes2num(out);
        a.stats.fundamentals = reshape(fundamentals, length(fundamentals)/3,3);
    end

    function npixperstixfn(out)
        a.stats.npixperstix = out;
    end

    function nstixpersidefn(out)
        a.stats.nstixperside = out;
    end

    function nframeratefn(out)
        a.stats.msperframe = 1000/out;
    end

    function gausslocutfn(out)
        a.stats.gausslocut = out/1000;
    end

    function gausshicutfn(out)
        a.stats.gausshicut = out/1000;
    end

    function linpredtolfn(out)
        gl.initparams.tolerance = out;
    end

    function scalefn(out)
        gl.initparams.stepsizescale = out;
    end

    function stepsizefn(out)
        gl.initparams.stepsize = out;
    end

    function nreversalsfn(out)
        gl.initparams.nreversals = out;
    end

    function oogscalefn(out)
        gl.initparams.oogscale = out;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting the trial parameters from the ecodes
function [seed, nframes, mu, sigma, bkgndrgb, noisetype, nt, gotfulltrial] = GetTrialParams(p)
C = WhiteNoisethreshCodes;
try
    seed = 0;
    nframes = 0;
    mu = 0;
    sigma = 0;
    bkgndrgb = 0;
    noisetype = 0;
    nt = 0;
    gotfulltrial = -1;
    
    trialstruct = {...
        C.SEEDCD,'long',@seedfn;...
        C.NFRAMESCD,'int',@nframesfn;...
        C.NEUROTHRESHCD,'int',@ntfn;...
        [C.MU1CD C.MU2CD C.MU3CD],'int',@mufn;...
        [C.SIGMA1CD C.SIGMA2CD C.SIGMA3CD],'int',@sigmafn;...
        C.NOISETYPECD,'int',@noisetypefn;...
        [C.BKGNDRCD C.BKGNDGCD C.BKGNDBCD],'int',@bkgndrfn;...
        };
    
    %            keyboard
    for i=1:size(trialstruct,1)
        out = GetVal(trialstruct{i,1}, p.events, trialstruct{i,2});
        if(~isempty(out))
            feval(trialstruct{i,3},out);
        end
    end
    gotfulltrial = 1;
    firststimon = find(p.events == C.STIMONCD, 1);
    firststimoff = find(p.events == C.STIMOFFCD, 1);
    if (firststimon > find(p.events == p.processnowcode, 1))
        firststimon = [];
    end
    if(isempty(firststimon))  % at least one STIMONCD
        error('no STIMONCD found');  %GDLH This is a problem!!
    end
    if all(bkgndrgb==0)
        disp('All background rgb is 0');
        gotfulltrial = 0; % This should fix the problem
% keyboard;
    end
    
catch
    seed = nan;
    nframes = nan;
    mu = nan;
    sigma = nan;
    bkgndrgb = nan;
    noisetype = nan;
    nt = -1;
    gotfulltrial = 0;
    
end
    function seedfn(out)
        seed = out;
    end

    function nframesfn(out)
        nframes = out;
    end

    function ntfn(out)
        nt = out;
    end

    function mufn(out)
        mu = out/1000;
    end

    function sigmafn(out)
        sigma = out/1000;
    end

    function noisetypefn(out)
        noisetype = out;
    end

    function bkgndrfn(out)
        bkgndrgb = out;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializing the statistics that
% will be filled when we get the header.
function InitStatsStruct()
a = get(gcf,'UserData');
a.stats.npixperstix = 0; % Number of pixels per side of stixel
a.stats.nstixperside = 0; % Number of stixels per side of stimulus
a.stats.msperframe = 0;  % ms per frame
a.stats.gammaTable = []; % gamma tables
a.stats.invgammaTable = []; % inverse gamma tables
a.stats.monSpd = []; % monitor phosphor spectra
a.stats.fundamentals = []; % cone fundamentals
a.stats.M = [];     % guns to cones matrix
a.stats.invM = [];  % cones to guns matrix
a.stats.gausslocut = 0; % low critical value for Gaussian cut off
a.stats.gausshicut = 0; % high critical value for Gaussian cut off
a.stats.gotHeader = 0;
set(gcf,'UserData',a);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializing the GL structure: to store the data for NeuroThresh trials
function InitGLStruct()
global  gl
gl.initparams.stepsize = [];
gl.initparams.stepsizescale = [];
gl.initparams.nreversals = [];
gl.initparams.oogscale = [];
gl.initparams.tolerance = []; 
gl.numtargetspikerate = 2;
gl.targetspikerates = [];
gl.trialspecidx = []; % index into gl.trialspecs indicating the current trial
gl.trialspecs = []; % A structure containing all the information about the trial
gl.timetoleave = 0;
gl.wtsstopsendingflag = 0;
gl.current_rem = 0;
gl.max_rem = 0;
gl.dirs_to_be_used = cell(1,gl.numtargetspikerate);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clearing out and preallocating space for
% statistics that don't require the header to
% be dropped (e.g. STS, STCross, and nspikes)
function ClearStats()
a = get(gcf,'UserData');
nframesback = str2double(get(a.uicontrols.nframes,'String'));
a.stats.murgb = [nan nan nan];
for i = 1:2
    a.stats.gunnoise.spike{i}.nspikes = 0;
    a.stats.conenoise.spike{i}.nspikes = 0;
end
a.stats.gunnoise.mu = [nan nan nan];
a.stats.gunnoise.sigma = [nan nan nan];
a.stats.conenoise.mu = [nan nan nan];
a.stats.conenoise.sigma = [nan nan nan];
a.stats.lmsbinaryrgbmat = nan*ones(8,3);
set(gcf,'UserData',a);
b = get(a.axeshandles.synthImage,'UserData');
b.weightstosend = {};
set(a.axeshandles.synthImage,'UserData', b);
if (a.stats.nstixperside)
    for i = 1:2  % preallocating space
        a.stats.gunnoise.spike{i}.STS = zeros([3*a.stats.nstixperside^2 nframesback]);
        a.stats.gunnoise.spike{i}.STCross = zeros([(3*a.stats.nstixperside^2)^2 nframesback]);
        a.stats.conenoise.spike{i}.STS = zeros([3*a.stats.nstixperside^2 nframesback]);
        a.stats.conenoise.spike{i}.STCross = zeros([(3*a.stats.nstixperside^2)^2 nframesback]);
    end
else 
    for i = 1:2
        a.stats.gunnoise.spike{i}.STS = zeros(size(a.stats.gunnoise.spike{i}.STS));
        a.stats.gunnoise.spike{i}.STCross = zeros(size(a.stats.gunnoise.spike{i}.STCross));
        a.stats.conenoise.spike{i}.STS = zeros(size(a.stats.conenoise.spike{i}.STS));
        a.stats.conenoise.spike{i}.STCross = zeros(size(a.stats.conenoise.spike{i}.STCross));
    end
end
set(gcf,'UserData',a);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Disable reset buttons during plotting
% and enable them afterward.
% arg = 0 disable, arg = 1 enable
function EnableResetButton(arg)
controls = getfield(get(gcf,'UserData'),'uicontrols');
if (arg)
    set(controls.reset,'Enable','on');
    set(controls.nframes,'Enable','on');
else
    set(controls.reset,'Enable','inactive');
    set(controls.nframes,'Enable','inactive');
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Disable the pixelmask ButtonDownFcn
% during UpdatePixelMask()
function EnablePixelMask(arg)
axeshandles = getfield(get(gcf,'UserData'),'axeshandles');
h = get(axeshandles.pixelmask,'Children');
if (~isempty(h))
    if verLessThan('matlab', '8.4.0.150421') % < R2014b - Graphics changes
        if (arg)
            set(h(1),'HitTest','on')
        else
            set(h(1),'HitTest','off')
        end
    else % >= R2014b; HitTest -> PickableParts
        if (arg)
            set(h(1),'PickableParts','visible')
        else
            set(h(1),'PickableParts','none')
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the figure for plotting
%   Subfunctions include:
%   SetUpAxes
%   InitSynthImage
function SetUpFig(DEFAULTNFRAMES)
global subunitcmap
gcf = figure(1);
set(gcf,'DefaultAxesUnits','pixels')
set(gcf,'position',[200 320 1200 400]);
set(gcf,'ButtonDownFcn','drawnow');  % in case user drags figure
clf;

a = get(gcf,'UserData');
a.uicontrols.reset = uicontrol('style','pushbutton','Callback',@ResetCallback, 'string','RESET','Position',[20 20 60 20]);
a.uicontrols.nframes = uicontrol('style','Edit','string',num2str(DEFAULTNFRAMES),'Callback',@ResetCallback,'Position',[855 10 40 20]);
a.uicontrols.eventqueuelength = uicontrol('style','text','string','Nan','Position',[15 50 70 30]);
a.uicontrols.STAslider = uicontrol('style','slider','Callback',@UpdateSynthImage,'Min',-1,'Max',1,'Value',0,'Position',[780 170 100 10]);
a.uicontrols.PCslider = uicontrol('style','slider','Callback',@UpdateSynthImage,'Min',-1,'Max',1,'Value',0,'Position',[780  95 100 10]);
a.uicontrols.STAcoefedit = uicontrol('style','edit','string','0','Position',[810 181 40 15],'Callback',@EditCallback);
a.uicontrols.PCcoefedit = uicontrol('style','edit','string','0','Position',[810 79 40 15],'Callback',@EditCallback);
a.uicontrols.alphaslider = uicontrol('style','slider','Callback',@UpdatePixelMask,'Min',0,'Max',.05,'SliderStep',[.01 .1],'Value',0,'Position',[20 185 60 10]);
a.uicontrols.alphatext = uicontrol('style','text','String','0','Position',[20 197 60 12]);
a.uicontrols.whichsigtest = uicontrol('style','popupmenu','Callback',@UpdatePixelMask,'String',{'Mean','Var','STA_r','STA_g','STA_b'},'Position',[25 220 50 15]);
a.uicontrols.imagebattery = uicontrol('style','pushbutton','string','ImageBattery','Callback',@SynthImageBatteryCallback,'Position',[780 200 100 20]);
a.uicontrols.imbatteryncontrasts = uicontrol('style','edit','string','4','Position',[850 220 20 20]);
a.uicontrols.imbatterytype = uicontrol('style','popupmenu','string',{'Standard','L vs M','FixCols'},'Position',[780 220 70 20]);
a.uicontrols.imbatteryclear = uicontrol('style','pushbutton','string','Clr Q','Callback',@SynthImageBatteryClear,'Position',[850 181 30 20]);
a.uicontrols.fitgabor = uicontrol('style','pushbutton','string','Gabr','Callback',@FitGaborCallback,'Position',[780 181 30 20]);
a.uicontrols.contrastnorm = uicontrol('style','radio','Position',[10,170,10,10],'Value',1);
a.uicontrols.smallpc = uicontrol('style','radio','Position',[85,125,10,10]);
a.uicontrols.whichspike = uicontrol('style','popupmenu','String',{'1','2'},'Position',[90 20 40 20]);
a.uicontrols.projortho = uicontrol('style','popupmenu','String',{'None','PC','STA'},'Position',[140 20 60 20]);
a.uicontrols.whichnoisetype = uicontrol('style','popupmenu','String',{'gun','cone'},'Position',[210 20 60 20]);
a.uicontrols.whichcovmex = uicontrol('style','popupmenu','String',{'STCOVmex','STCOV_st'},'Position',[280 20 80 20]);

SetUpAxes();
set(gcf,'UserData',a);
drawnow;

%%%%%%%%%%%%%%%%%%%%
% Setting up the plotting axes
    function SetUpAxes()
        AXESWIDTH = 50;
        %a = get(gcf,'UserData');
        nframes = str2double(get(a.uicontrols.nframes,'String'));
        h1 = []; h2 = [];
        figpos = get(gcf,'Position');
        for i = 1:nframes
            h = axes('position',[(1+1/nframes)*AXESWIDTH*(i-1)+100 170 AXESWIDTH AXESWIDTH]);
            set(gca,'XTick',[],'YTick',[],'Box','on');
            axis image;
            h1 = [h1; h];  % STA
            
            h = axes('position',[(1+1/nframes)*AXESWIDTH*(i-1)+100 85 AXESWIDTH AXESWIDTH]);
            set(h,'XTick',[],'YTick',[],'Box','on');
            axis image;
            h2 = [h2; h];  % PC
        end
        a.axeshandles.STA = h1;
        a.axeshandles.PC = h2;
        a.axeshandles.text = axes('position',[400 50 1 1]); % For nspikes
        set(gca,'Visible','off');
        a.axeshandles.synthImage = axes('position',[805 115 AXESWIDTH AXESWIDTH]); % For the SynthImage
        InitSynthImage(a.axeshandles.synthImage);
        a.axeshandles.pixelmask = axes('position',[25 130 AXESWIDTH AXESWIDTH]); % For the pixelmask
        Initpixelmask(a.axeshandles.pixelmask);
        % subunit axes
        a.axeshandles.subunit = axes('Position',[678 130-.25*AXESWIDTH 1.5*AXESWIDTH 1.5*AXESWIDTH]);
        InitSubunit(a.axeshandles.subunit)
        a.axeshandles.chromatic_STA_image = axes('position',[25 270 2*AXESWIDTH 2*AXESWIDTH]);
        Initchromatic_STA_image(a.axeshandles.chromatic_STA_image);
        a.axeshandles.temporal_chromatic_STA = axes('position',[25+2*AXESWIDTH+20 270 9*AXESWIDTH 2*AXESWIDTH]);
        Inittemporal_chromatic_STA(a.axeshandles.temporal_chromatic_STA);
        a.axeshandles.baseline_hist = axes('position',[25+11*AXESWIDTH+50 270 3*AXESWIDTH 2*AXESWIDTH]);
        Initbaselinehist(a.axeshandles.baseline_hist);
        a.axeshandles.NT_Raster_plot = axes('position',[25+15*AXESWIDTH+30 270 3*AXESWIDTH 2*AXESWIDTH]);
        InitNT_Raster_plot(a.axeshandles.NT_Raster_plot);
        a.axeshandles.text_latency = axes('position',[400 50 1 1]); text(-10,5,['Latency']);
        a.axeshandles.text_TFR = axes('position',[400 50 1 1]); text(100,5,['TFR']);
        a.axeshandles.text_SR = axes('position',[400 50 1 1]); text(185,5,['Spikerate']);
        a.axeshandles.spatialstructure = axes('position',[905 175 1.4*AXESWIDTH 1.4*AXESWIDTH]); set(gca,'XTick',[],'YTick',[]);
        Init_spatial_structure(a.axeshandles.spatialstructure);
        a.axeshandles.conewts = axes('position',[905 85 1.4*AXESWIDTH 1.4*AXESWIDTH]); set(gca,'XTick',[],'YTick',[]);
        Init_cone_wts(a.axeshandles.conewts);
        a.axeshandles.dir_TFR1 = axes('position',[990 230 3*AXESWIDTH 3*AXESWIDTH]); set(gca,'XTick',[],'YTick',[]);
        a.axeshandles.dir_TFR2 = axes('position',[990 30 3*AXESWIDTH 3*AXESWIDTH]); set(gca,'XTick',[],'YTick',[]);

        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Initialize the tmp image
        function InitSubunit(axesh)
            global gl
            maxsubunits = size(subunitcmap, 1);
            subunit_pos = get(axesh, 'Position');
            axes_w = subunit_pos(3);
            slider_w = axes_w; slider_h = 15;
            slsubunit_pos = [subunit_pos(1) subunit_pos(2)+axes_w+5 slider_w slider_h];
            updatemask_h = 22;
            updatemask_pos = [subunit_pos(1) subunit_pos(2)-updatemask_h-5 axes_w updatemask_h];
            text_w = 22; text_h = text_w;
            clearmask_pos = [updatemask_pos(1) updatemask_pos(2)-updatemask_h-30 updatemask_pos(3:4)];
            nudgestim_pos = [updatemask_pos(1) updatemask_pos(2)-updatemask_h-5 updatemask_pos(3:4)];
            NT_pos = [updatemask_pos(1) updatemask_pos(2)-updatemask_h-49 updatemask_pos(3:4)];
            spikerate_disp = [updatemask_pos(1)-100 updatemask_pos(2)-updatemask_h-49 updatemask_pos(3:4)];
            latency_disp = [updatemask_pos(1)-300 updatemask_pos(2)-updatemask_h-49 updatemask_pos(3:4)];
            currenttargetFR_disp = [updatemask_pos(1)-200 updatemask_pos(2)-updatemask_h-49 updatemask_pos(3:4)];
            num_targetspikerate_pos = [updatemask_pos(1)+90 updatemask_pos(2)-updatemask_h-30 updatemask_pos(3)/2 updatemask_pos(4)];
            des_spikerate_disp = [updatemask_pos(1)+ 140 updatemask_pos(2)-updatemask_h-30 updatemask_pos(3)/2 updatemask_pos(4)];
            des_spikerate_disp2 = [updatemask_pos(1)+ 190 updatemask_pos(2)-updatemask_h-30 updatemask_pos(3)/2 updatemask_pos(4)];
            simul_dirs_to_probe = [updatemask_pos(1)+ 240 updatemask_pos(2)-updatemask_h-30 updatemask_pos(3)/2 updatemask_pos(4)];
            baseisvectypepos = [updatemask_pos(1)+90 updatemask_pos(2)-updatemask_h-60 updatemask_pos(3) updatemask_pos(4)];
            
            set(axesh,'UserData',struct('mask',[]));
            a.uicontrols.txsubunit = uicontrol('Style','edit','ForegroundColor',subunitcmap(1,:),'FontSize',12,'FontWeight','bold','String','1','Enable','inactive','Position',[slider_w/2-text_w/2+slsubunit_pos(1) slsubunit_pos(2)+slider_h+5 text_w text_h]);
            a.uicontrols.slsubunit = uicontrol('Style','slider','Callback',{@update_subunit_text,a.uicontrols.txsubunit},'Min',1,'Max',9,'SliderStep',[1/(maxsubunits-1) 1/(maxsubunits-1)],'Value',1,'Position',slsubunit_pos);
            a.uicontrols.updatemask = uicontrol('style','togglebutton','string','Update Mask','Position',updatemask_pos,'Min',0,'Max',1,'Value',0);
            a.uicontrols.clearmask = uicontrol('style','pushbutton','string','Clear Mask','Position',clearmask_pos,'Callback',{@clear_subunit_mask,a.uicontrols.updatemask,axesh});
            a.uicontrols.nudgestim = uicontrol('style','pushbutton','string','Nudge stim','Position',nudgestim_pos,'Callback',@nudge_stim);
            a.uicontrols.NT = uicontrol('style','togglebutton','string','NeuroThresh','Position',NT_pos,'Callback',@NT_gen_stim); % Added by Abhishek on 11/15
            a.uicontrols.spikerate = uicontrol('style','text','string','0','Position',spikerate_disp);
            a.uicontrols.des_spikerate = uicontrol('style','Edit','string','0','Position',des_spikerate_disp);
            a.uicontrols.des_spikerate2 = uicontrol('style','Edit','string','0','Position',des_spikerate_disp2);
            a.uicontrols.simul_dirs_to_probe = uicontrol('style','Edit','string','2','Position',simul_dirs_to_probe);
            a.uicontrols.latency = uicontrol('style','Edit','string','0','Position',latency_disp);
            a.uicontrols.currenttargetFR_disp = uicontrol('style','Edit','string','0','Position',currenttargetFR_disp);
            a.uicontrols.num_targetspikerate = uicontrol('style','Edit','string',num2str(gl.numtargetspikerate),'Position', num_targetspikerate_pos);
            a.uicontrols.basisvectype = uicontrol('style','popupmenu','string',{'STAvsPC','subunits'},'Position',baseisvectypepos);
            
            set(axesh, 'XTick', [], 'YTick', [], 'Box', 'on');
            axis(axesh, 'image');
            grid(axesh, 'off');
            
        end 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the synthetic image
function InitSynthImage(axesh)
set(axesh,'XTick',[],'YTick',[],'Box','on');
global gl
axis image;
b.STA = [];
b.PC = [];
b.gabor = [];
b.localimage = [];
b.localweights = [];
b.weightstosend = {};
b.muimage = [];
b.baseline = [];
b.spiketimes = [];
b.basisvectosend = [];
b.basissendingflag = 0;
gl.wtsstopsendingflag = 0;
b.wts_NT = [];
set(axesh,'UserData',b);
end

%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the pixel mask (sig. tests and mask for PC
% calculations).
function Initpixelmask(axesh)
axes(axesh);
set(axesh,'XTick',[],'YTick',[],'Box','on');
axis image;
c = struct('pixelmask',[]);
set(axesh,'UserData',c);
end

%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the image whose pixel if selected broadcasts the temporal chromatic variation of that pixel 
function Initchromatic_STA_image(axesh)
axes(axesh);
set(axesh,'XTick',[],'YTick',[],'Box','on');
axis image;
c.img = [];
c.select_units = [];
set(axesh,'UserData',c);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the cone weights for the subnits
function Init_cone_wts(axesh)
axes(axesh);
set(axesh,'XTick',[],'YTick',[],'Box','on');
axis image
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the spatial strucure STA of the checkerboard noise 
function Init_spatial_structure(axesh)
axes(axesh);
set(axesh,'XTick',[],'YTick',[],'Box','on');
axis image
end

%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the axes plot which will show the temporal chromatic STA variation 
function Inittemporal_chromatic_STA(axesh)
axes(axesh);
set(axesh,'XTick',[],'YTick',[],'Box','on');
% axis (axesh,'square');
c = struct('temporal_chromatic_STA',[]);
set(axesh,'UserData',c);
end

function Initbaselinehist(axesh)
axes(axesh);
set(axesh,'XTick',[],'YTick',[],'Box','on');
c = struct('baseline',[]);
set(axesh,'UserData',c);
end

function InitNT_Raster_plot(axesh)
axes(axesh);
set(axesh,'XTick',[],'YTick',[],'Box','on');
c = struct('Raster',[]);
set(axesh,'UserData',c);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Force the STAslider and PCslider to the edit values
function EditCallback(~, ~)
controls = getfield(get(gcf,'UserData'),'uicontrols');
STAcoef = str2double(get(controls.STAcoefedit,'string'));
if (STAcoef <= -1)
    STAcoef = -1+eps;
elseif (STAcoef >= 1)
    STAcoef = 1-eps;
end

PCcoef = str2double(get(controls.PCcoefedit,'string'));
if (PCcoef <= -1)
    PCcoef = -1+eps;
elseif (PCcoef >= 1)   % Abhishek made a change
    PCcoef = 1-eps;
end
set(controls.STAslider,'Value',STAcoef);
set(controls.PCslider,'Value',PCcoef);
UpdateSynthImage;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callback invoked when we click the reset button or change the number of
% frames to look at.
function ResetCallback(~, ~)
a = get(gcf, 'UserData');
SetUpFig(str2double(get(a.uicontrols.nframes,'String')));
a = get(gcf, 'UserData'); % this isn't a typo, the previous function modifies UserData
UpdatePixelMask(nan,nan);
SynthImageBatteryClear(nan,nan)
ClearStats(); % Just the things that aren't in the header
set(a.uicontrols.updatemask, 'Value', 0); % don't use (send to REX) the current mask
init_subunit_axes(); % clear the subunit mask
paint_subunit_selection(a.axeshandles.subunit); % delete the mask's representation on the plot
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callback invoked when we click on one of the STA or PC images
function ImageCallback(~, ~, axishandle, v)
EnableResetButton(0);
handles = getfield(get(gcf,'UserData'),'axeshandles');
b = get(handles.synthImage,'UserData');
if (any(ismember(handles.STA, axishandle)))
    hs = handles.STA;
elseif (any(ismember(handles.PC, axishandle)))
    hs = handles.PC;
end
data = [];
idx = WhichFrameSelected(hs);
if (isempty(idx))  % Nothing selected yet
    set(axishandle,'Xcolor',[1 1 0],'Ycolor',[1 1 0]);
    data = v;
else
    if (idx == find(axishandle == hs))  % Deselecting
        set(axishandle,'Xcolor',[0 0 0],'Ycolor',[0 0 0]);
    else   % Selecting one and deselecting another
        set(hs(idx),'Xcolor',[0 0 0],'Ycolor',[0 0 0]);
        set(axishandle,'Xcolor',[1 1 0],'Ycolor',[1 1 0]);
        data = v;
    end
end
if (any(ismember(handles.STA, axishandle)))
    b.STA = data;
    update_subunit_axes(data);
    update_chromatic_STA_image_axes(data);
elseif (any(ismember(handles.PC, axishandle)))
    b.PC = data;
    update_subunit_axes(data);
    update_chromatic_STA_image_axes(data);    % Abhishek made change 08/2015
end
set(handles.synthImage,'UserData',b);
UpdateSynthImage(nan,nan);
[n_subunits,~] = check_num_subunits() ; % 0 corresponds to no subunits
if (~n_subunits)
    UpdatePixelMask(nan,nan);
end
EnableResetButton(1);
end

% Return the index of the selected frame
function idx = WhichFrameSelected(handles)
idx = [];
for i = 1:length(handles)
    if all(get(handles(i),'Xcolor') == [1 1 0])
        idx = [idx i];
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callback invoked when you click the NeuroThresh Button
% Prepare an image to be sent to REX
function NT_gen_stim(~,~)
disp('Neurothresh Mode ON');
a = get(gcf,'UserData');
global gl
mu = a.stats.murgb;
nstixperside = a.stats.nstixperside;
muimage = cat(3, repmat(mu(1),[nstixperside nstixperside]),...
    repmat(mu(2),[nstixperside nstixperside]),...
    repmat(mu(3),[nstixperside nstixperside]));
uicontrols = getfield(get(gcf,'UserData'),'uicontrols');
STAvalue = get(uicontrols.STAslider,'value');
PCvalue = get(uicontrols.PCslider,'value');
b = get(a.axeshandles.synthImage,'UserData');
b.muimage = muimage;
gl.numtargetspikerate = str2num(get(a.uicontrols.num_targetspikerate,'string'));
b.weightstosend = {}; % clearing off stored weights in the cell 
STAimage = b.STA;
if (isempty(STAimage))
    STAimage = muimage;
end

basisvectypestrings = get(a.uicontrols.basisvectype,'String');
basisvectypeval = get(a.uicontrols.basisvectype,'Value');
[num_subunits, Subunit_mask] = check_num_subunits();
if strcmp(basisvectypestrings(basisvectypeval),'subunits') & (num_subunits==2)
    val = [-sqrt(1/2) sqrt(1/2); 0 1; sqrt(1/2) sqrt(1/2); 1 0; sqrt(1/2) -sqrt(1/2); 0 -1; -sqrt(1/2) -sqrt(1/2); -1 0]; % new starting directions, A bit of 2nd and 4th quadrant included- Abhishek 11/16
    tmp1 = create_subunit_basisvec(nstixperside,STAimage(:),Subunit_mask,1,b.muimage(:))-b.muimage;
    tmp2 = create_subunit_basisvec(nstixperside,STAimage(:),Subunit_mask,2,b.muimage(:))-b.muimage;

elseif strcmp(basisvectypestrings(basisvectypeval),'subunits') & (num_subunits>2) % for more than 2 subunits, Abhishek - 2/17
    val = [-sqrt(1/2) sqrt(1/2); 0 1; sqrt(1/2) sqrt(1/2); 1 0; sqrt(1/2) -sqrt(1/2)]; 
    tmp1 = create_subunit_basisvec(nstixperside,STAimage(:),Subunit_mask,1:num_subunits-1,b.muimage(:))-b.muimage;
    tmp2 = create_subunit_basisvec(nstixperside,STAimage(:),Subunit_mask,num_subunits,b.muimage(:))-b.muimage;
else 
    PCimage = b.PC;
    if (isempty(PCimage))
        PCimage = muimage;
    end
    val = [1 0; 0 1; -1 0; 0 -1]; % new starting directions STAvsPC search directions, 360 degree search implemented 12/16
    tmp1 = STAimage-b.muimage;
    tmp2 = PCimage-b.muimage;
end

for j = 1:gl.numtargetspikerate
    
    if str2num(get(a.uicontrols.des_spikerate,'String'))==0 
        gl.targetspikerates(j) = mean(b.baseline) + 3*j*std(b.baseline);
        % Minimum targetspikerate is 20 spikes/sec.
        if gl.targetspikerates(j) < 20 % a bad hack to prevent the targetfiringrate from being 0
            gl.targetspikerates(j) = 20*j;
        end
    else
        % Not sure if this would be the right thing to do in case of multiple targetspikerates
        if j == 1
            gl.targetspikerates(j) = str2num(get(a.uicontrols.des_spikerate,'String'));
        elseif j == 2
            gl.targetspikerates(j) = str2num(get(a.uicontrols.des_spikerate2,'String'));
        end
    end
    
    for i = 1:size(val,1)
        if (norm(val(i,:)) > 0)
            weights = val(i,:)/norm(val(i,:)); % normalizing
        else
            weights = [0 0];
        end
        SetUpTrialSpecs(weights, gl.targetspikerates(j))
        b.weightstosend{i} = weights;
    end
end
b.basisvectosend{1} = tmp1/norm(tmp1(:));
b.basisvectosend{2} = tmp2/norm(tmp2(:));
b.basissendingflag = 1;

% Introduce the calculations for determining the latency and offset
nframesback = str2double(get(a.uicontrols.nframes,'String'));
[n_subunits,~] = check_num_subunits(); % 0 corresponds to no subunits
if (n_subunits)
    nrandnums_perchannel = n_subunits;
else
    nrandnums_perchannel = a.stats.nstixperside^2;
end
whichspike = get(a.uicontrols.whichspike,'Value');
if (get(a.uicontrols.whichnoisetype,'Value') == 1) % gun noise
    STS = a.stats.gunnoise.spike{whichspike}.STS;
     n = a.stats.gunnoise.spike{whichspike}.nspikes;
else % cone noise
    STS = a.stats.conenoise.spike{whichspike}.STS;
     n = a.stats.gunnoise.spike{whichspike}.nspikes;
end

% Latency caluclation - Assumes subunit WhiteNoise has been used before
% clicking the NT button
tmp_basis = STS;
frametimes = [1:nframesback]*a.stats.msperframe;
STAenergy =  sum(tmp_basis.^2);
t_offset = frametimes(STAenergy==max(STAenergy));
if isempty(t_offset) 
    % in case the conditions for calculating t_offset fails
    % You wouldn't need this condition for the actual neuron
    % Effective when using for just testing without a neuron
    t_offset = 40;
end
disp(val);
set(a.uicontrols.latency,'String',num2str(t_offset));
set(a.axeshandles.synthImage,'Userdata',b);
set(gcf,'UserData',a);   
 
end

function out_im = create_subunit_basisvec(nstixperside,img,Subunit_mask,i,bkgnd)
% The aim of this function is create a basis vector out of an input STA
% image and the selected subunit whose index is denoted by the letter 'i'. The
% output image is stored in 'out_im' which has the just 1 subunit
% highlighted while the other subunit and the background painted in background 
% color

img_r = img(1:nstixperside.^2);
img_g = img(nstixperside.^2+1:2*nstixperside.^2);
img_b = img(2*nstixperside.^2+1:3*nstixperside.^2);
vec_r = bkgnd(1:nstixperside.^2);
vec_g = bkgnd(nstixperside.^2+1:2*nstixperside.^2);
vec_b = bkgnd(2*nstixperside.^2+1:3*nstixperside.^2);

if numel(i)==1
    idx = find(Subunit_mask(:) == i);
    vec_r(idx) = img_r(idx);
    vec_g(idx) = img_g(idx);
    vec_b(idx) = img_b(idx);
else
    % Adding some more functionality - Abhishek, 2/17
    for jj = i
        idx = find(Subunit_mask(:) == i(jj));
        vec_r(idx) = img_r(idx);
        vec_g(idx) = img_g(idx);
        vec_b(idx) = img_b(idx);
    end
end
out_im = reshape([vec_r; vec_g; vec_b], [nstixperside nstixperside 3]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callback invoked when you click  the synthetic image.
% Prepare an image to be sent to REX.
function SynthImageCallback(~, ~, axeshandle)
b = get(axeshandle,'UserData');

if (isempty(b.weightstosend))
    b.weightstosend = {b.localweights};
else
    b.weightstosend = [b.weightstosend {b.localweights}];
end
set(axeshandle,'UserData',b);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callback invoked when you click the ImageBatteryButton
% Prepare a series of images to send to REX.  The images
% are maintained in the UserData field of the synthImage axes
% as a cells array of NxNx3 matrices that are sent to REX
% one after the other, as "imagerequest()" messages come in.
%
% If the type of image battery is selected to be "Standard"
% we show a 1-D or 2-D (depending on whether an STA and/or
% PC1 image is selected) linear image set.

function SynthImageBatteryCallback(~, ~)
a = get(gcf,'UserData');
handles = a.axeshandles;
ncontrastlevels = str2double(get(a.uicontrols.imbatteryncontrasts,'string'));
b = get(handles.synthImage,'UserData');

mu = a.stats.murgb;  % mean in [0:1] intensity units
nstixperside = a.stats.nstixperside;
muimage = cat(3, repmat(mu(1),[nstixperside nstixperside]),...
    repmat(mu(2),[nstixperside nstixperside]),...
    repmat(mu(3),[nstixperside nstixperside]));
b.muimage = muimage;
imtypestrings = get(a.uicontrols.imbatterytype,'String');
imtypeval = get(a.uicontrols.imbatterytype,'Value');
STAimage = b.STA;
PCimage = b.PC;
STAlim = get(a.uicontrols.STAslider,'value');
PClim = get(a.uicontrols.PCslider,'value');
STAcontrasts = linspace(-STAlim,STAlim,ncontrastlevels);
PCcontrasts = linspace(-PClim,PClim,ncontrastlevels);

if (isempty(STAimage) || STAlim == 0)
    STAimage = muimage;
    STAcontrasts = 0;
end
if (isempty(PCimage) || PClim == 0)
    PCimage = muimage;
    PCcontrasts = 0;
end

tmp1 = STAimage-b.muimage;
tmp2 = PCimage-b.muimage;
b.basisvectosend{1} = tmp1/norm(tmp1(:));
b.basisvectosend{2} = tmp2/norm(tmp2(:));
b.basissendingflag = 1;

% Calling the appropriate nested subfunction
if (strcmp(imtypestrings(imtypeval),'Standard'))
    weights = StandardImageBattery;
end

b.weightstosend = [b.weightstosend weights(randperm(length(weights)))];
set(a.axeshandles.synthImage,'UserData',b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nested function for doing a standard image battery
    function weights = StandardImageBattery
        weights = cell(0);
        for STAvalue = STAcontrasts
            for PCvalue = PCcontrasts
                w = [STAvalue PCvalue];
                w = w/norm(w);
                weights = [weights {w}];
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compresses an image with subunits such that the image could be represented by fewer numbers
function img = compress_vector(images,num_subunits,subunit_mask)
subunit_mask(subunit_mask == 0) = num_subunits + 1;
img = cell(0);
a = get(gcf,'UserData');
new_mask = subunit_mask(:);
[~,ind,~] = unique(new_mask);
for ii=1:size(images,2)
    tmp_img = images{1,ii};
    for gun = 1:3
        vec = tmp_img(:,:,gun);
        tmp_img(:,:,gun) = zeros(a.stats.nstixperside,a.stats.nstixperside);
        vec_tmp = vec(:);
        tmp_img(1:num_subunits+1,1,gun) = vec_tmp(ind) ;
    end
    img = [img {tmp_img}];
end
end

% Callback function that gets invoked when you click the "Clr Q"
% "Clear Queue" button.  It deletes all the images remaining in the
% queue.  Should be useful if you accidentally request the wrong
% image playback.
function SynthImageBatteryClear(~,~)
global gl
a = get(gcf,'UserData');
b = get(a.axeshandles.synthImage,'UserData');
b.localimage = [];
b.weightstosend = [];
b.localweights = [];
b.muimage = [];
b.baseline = [];
b.spiketimes = [];
b.basisvectosend = [];
b.wts_NT = [];
NT_val = get(a.uicontrols.NT,'Value');
set(a.uicontrols.latency,'String','0');
set(a.uicontrols.des_spikerate,'String','0');
set(a.uicontrols.des_spikerate2,'String','0');
set(a.uicontrols.currenttargetFR_disp,'String',0);
set(a.uicontrols.spikerate,'String','0');
set(a.uicontrols.NT,'Value',0);
set(a.axeshandles.synthImage,'UserData',b);
set(gcf,'UserData',a);
 
gl.numtargetspikerate = 2; % By default the algorithm is designed to measure iso-response curve for 2 different target firing rates
gl.targetspikerates = [];
gl.trialspecidx = []; % index into gl.trialspecs indicating the current trial
gl.trialspecs = []; % A structure containing all the information about the trial
gl.timetoleave = 0;
gl.wtsstopsendingflag = 0;
end

% Destructively modifies the 'pixelmask' field in the axes
% UserData to mark the user's mouse clicks.  Doesn't do any
% plotting though - that's taken care of in UpdatePixelMask.
function PixelMaskCallback(~,~)
a = get(gcf,'UserData');
%stats = a.stats;
[num_units, ~] = check_num_subunits();
if  ~(get(a.uicontrols.NT,'Value') | num_units)
    axeshandles = a.axeshandles;
    b = get(axeshandles.pixelmask,'UserData');
    if (~isempty(b.pixelmask))
        whichpt = get(gca,'CurrentPoint');
        whichpt = round(whichpt(1,[1 2]));
        whichpt = min([whichpt; size(b.pixelmask)]);
        b.pixelmask(whichpt(2), whichpt(1)) = ~b.pixelmask(whichpt(2), whichpt(1));
        set(axeshandles.pixelmask,'UserData',b);
        UpdatePixelMask(nan,nan);
    end
end
end

% Updates the Basline histogram after each successfull WhiteNoise Trial
function update_baseline_hist(baseline_firing_rates)
a = get(gcf, 'UserData');
ax = a.axeshandles.baseline_hist;
axes(ax);
if ~isempty(baseline_firing_rates)
    hist(baseline_firing_rates);hold on;
    ylim = get(ax,'Ylim');
    avg = mean(baseline_firing_rates); 
    plot(avg,ylim(2),'kv','MarkerFacecolor','r'); hold on;
    plot(avg + std(baseline_firing_rates),ylim(2),'kv','MarkerFacecolor','g');
    grid(ax, 'off'); hold off;
end
set(ax,'UserData',baseline_firing_rates);
end

% Updates the Raster plot after each successful Neurothresh Trial
function update_NT_raster_plot(spiketimes, offset)
a = get(gcf, 'UserData');
ax = a.axeshandles.NT_Raster_plot;
axes(ax);
plot([spiketimes, spiketimes]',[-.5; .5]*ones(1,length(spiketimes)),'k-');
line([0 0],[-1 1],'Color',[1 0 0])
line([offset offset],[-1 1],'Color',[1 0 0])
set(gca,'Xlim',[0 310])
set(ax, 'UserData', spiketimes);
grid(ax, 'off');
end

% Updates the contents stored in the subunit axes
function update_subunit_axes(STA)
a = get(gcf,'UserData');
ax = a.axeshandles.subunit;
UD = get(ax,'UserData');
if ~isempty(STA)
    image(STA, 'Parent', ax, 'ButtonDownFcn', @update_subunit_mask);
end
set(ax, 'UserData', UD);
set(ax, 'XTick', [], 'YTick', [], 'Box', 'on');
axis(ax, 'image');
grid(ax, 'off');
paint_subunit_selection(ax);
end

function update_cone_wts(data,ax1)
% I am introducing an axis which gives a bar graph of cone wts of subunits
[num_units, mask] = check_num_subunits();
if num_units>0 & ~isempty(data)
    a = get(gcf,'UserData');
    mask = mask(:);
    dataR = squeeze(data(:,:,1)); dataR = dataR(:);
    dataG = squeeze(data(:,:,2)); dataG = dataG(:);
    dataB = squeeze(data(:,:,3)); dataB = dataB(:);
    LMS = [];
    for ii = 1:num_units
        tmp = [unique(dataR(mask==ii))-unique(dataR(mask==0)); unique(dataG(mask==ii))-unique(dataR(mask==0)); unique(dataB(mask==ii))-unique(dataR(mask==0))];
        tmp = a.stats.invM'*tmp;
        tmp = tmp./norm(tmp);
        LMS = [LMS tmp];
    end
    bar(LMS,'Parent',ax1);
    set(ax1, 'XTick',[1 2 3], 'XTickLabel',{'L','M','S'});
    grid(ax1, 'off');
end
end

function update_spatial_structure(data,ax1)
% I am introducing an axis which gives a bar graph of cone wts of subunits
[num_units, ~] = check_num_subunits();
if num_units==0 & ~isempty(data)
    % Enter only if in the phase 1 of the experiment: Checkerboard noise
    % stage
    a = get(gcf,'UserData');
    dataR = squeeze(data(:,:,1)); dataR = dataR(:);
    dataG = squeeze(data(:,:,2)); dataG = dataG(:);
    dataB = squeeze(data(:,:,3)); dataB = dataB(:);
    tmpdata = [dataR dataG dataB];
    [v,~,~] = svd(tmpdata);
    img = reshape(v(:,1),a.stats.nstixperside,a.stats.nstixperside);
    imagesc(sign(img).*abs(img),'Parent',ax1); colormap('gray');
    set(ax1, 'XTick',[], 'YTick',[]);
    grid(ax1, 'off'); axis square;
end
end

% Updates the chromatic_STA_image axes 
function update_chromatic_STA_image_axes(STA)
a = get(gcf, 'UserData');
ax = a.axeshandles.chromatic_STA_image;
UD = get(ax, 'UserData');
UD.img = STA;
if ~isempty(STA)
    image(STA, 'Parent', ax,'ButtonDownFcn', @update_chromatic_STA_image_mask);
end
set(ax, 'UserData', UD);
set(ax, 'XTick', [], 'YTick', [], 'Box', 'on');
axis(ax, 'image');
grid(ax, 'off');
update_cone_wts(STA,a.axeshandles.conewts);
update_spatial_structure(STA,a.axeshandles.spatialstructure);
end

% Updates the chromatic_STA_image axes
function update_chromatic_STA_image_mask(himage, ~) % buttondown callback for subunit image
global subunitcmap
if verLessThan('matlab', '8.4.0.150421') % R2014b - Graphics changes
    HT = {'HitTest' 'off'}; % old method
else
    HT = {'PickableParts' 'none'}; % new method
end

axs = get(himage, 'Parent');
st = get(gcf, 'SelectionType');
UD = get(axs, 'UserData');
whichpt = get(axs, 'CurrentPoint');
whichpt = round(whichpt(1,1:2));
if strcmp(st, 'normal') % left click
    UD.select_units(UD.select_units == 1) = nan;
    delete(findobj(get(axs, 'children'), 'type', 'rectangle'));
    pos = [whichpt(1)-.5 whichpt(2)-.5 1 1]';
    rectangle('Position', pos, 'FaceColor', subunitcmap(1,:), HT{:});
    UD.select_units(whichpt(1), whichpt(2)) = 1;
    plot_waveform(whichpt(1),whichpt(2)); % for plotting the waveforms
elseif strcmp(st, 'alt') % right click
    UD.select_units(whichpt(1), whichpt(2)) = nan;
    delete(findobj(get(axs, 'children'), 'type', 'rectangle'));
end
set(axs, 'UserData', UD);
end

% plots the temporal chromatic signal of the selected pixel
function plot_waveform(x,y)
a = get(gcf,'UserData');
nframesback = str2double(get(a.uicontrols.nframes,'String'));
b = zeros(3,nframesback);
for jj = 1:nframesback
    tmp_ax_handle = get(a.axeshandles.STA(jj),'UserData');
    b(:,jj) = squeeze(tmp_ax_handle(y,x,:));
end 
if ~isempty(b)
    axes(a.axeshandles.temporal_chromatic_STA)
    plot(b(1,:)-a.stats.murgb(1),'--ro','Linewidth',2); hold on;
    plot(b(2,:)-a.stats.murgb(2),'--go','Linewidth',2); hold on;
    plot(b(3,:)-a.stats.murgb(3),'--bo','Linewidth',2); hold off;
    set(a.axeshandles.temporal_chromatic_STA,'xtickLabel',[0:-1:-nframesback+1]);
end
set(a.axeshandles.temporal_chromatic_STA, 'UserData', b);

end

function paint_subunit_selection(ax)
UD = get(ax, 'UserData');
paint_rectangles(UD.mask_display, ax);
set(ax, 'UserData', UD);
end

function paint_rectangles(mask, ax)
global subunitcmap
if verLessThan('matlab', '8.4.0.150421') % R2014b - Graphics changes
    HT = {'HitTest' 'off'}; % old method
else
    HT = {'PickableParts' 'none'}; % new method
end
axes(ax);
delete(findobj(get(ax, 'children'), 'type', 'rectangle'));
if any(~isnan(mask(:)))
    maxsubunits = size(subunitcmap, 1);
    for n = 1:maxsubunits
        [I,J] = find(mask == n);
        for pos = [J-.5 I-.5 ones(length(I), 2)]'
            rectangle('Position', pos, 'FaceColor', subunitcmap(n,:), HT{:});
        end
    end
end
end

% Initialises the contents of the subunit axes
function init_subunit_axes()
a = get(gcf, 'UserData');
UD = get(a.axeshandles.subunit, 'UserData');
UD.mask_rex = nan(a.stats.nstixperside);
UD.mask_display = nan(a.stats.nstixperside);
UD.nudge = [];
set(a.axeshandles.subunit, 'UserData', UD);
UD1 = get(a.axeshandles.chromatic_STA_image, 'UserData');
UD1.select_units = nan(a.stats.nstixperside);
set(a.axeshandles.chromatic_STA_image, 'UserData', UD1);
end

function update_subunit_text(h, ~, htext)
global subunitcmap
n = round(get(h, 'Value'));
set(h, 'Value', n);
set(htext, 'String', n, 'ForegroundColor', subunitcmap(n,:));
end

% Callback function when click 'Clear Mask'
function clear_subunit_mask(~, ~, hupdatemask, haxsub)
set(hupdatemask, 'Value', 0);
init_subunit_axes();
% Abhishek made a change - 09/2015
UD = get(gcf,'UserData');
nframesback = str2double(get(UD.uicontrols.nframes,'String'));
for i = 1:2
    UD.stats.gunnoise.spike{i}.STS = zeros([3*UD.stats.nstixperside^2 nframesback]);
    UD.stats.gunnoise.spike{i}.STCross = zeros([(3*UD.stats.nstixperside^2)^2 nframesback]);
    UD.stats.gunnoise.spike{i}.nspikes = 0;
    UD.stats.conenoise.spike{i}.STS = zeros([3*UD.stats.nstixperside^2 nframesback]);
    UD.stats.conenoise.spike{i}.STCross = zeros([(3*UD.stats.nstixperside^2)^2 nframesback]);
    UD.stats.conenoise.spike{i}.nspikes = 0;
end
set(gcf,'UserData',UD);
paint_subunit_selection(haxsub);
end

% Callback function when you click 'nudge stim'
function nudge_stim(~,~)
a = get(gcf, 'UserData');
UD = get(a.axeshandles.subunit, 'UserData');
x = -(a.stats.nstixperside-1)/2:(a.stats.nstixperside-1)/2;
normmask = (UD.mask_display >= 1)./sum(UD.mask_display(:) >= 1); % normalized to sum to 1
UD.nudge = [sum(normmask)*x' fliplr(x)*sum(normmask,2)]; % +x is right, +y is up. In stixels.
set(a.axeshandles.subunit, 'UserData',UD);
end

% buttondown callback for subunit image
function update_subunit_mask(himage, ~) 
axsubunit = get(himage, 'Parent');
st = get(gcf, 'SelectionType');
UD = get(axsubunit, 'UserData');
whichpt = get(axsubunit, 'CurrentPoint');
whichpt = round(whichpt(1,1:2));
whichpt = min(max(1, whichpt), size(UD.mask_display, 2));

if strcmp(st, 'normal') % left click
    UD.mask_display(whichpt(2), whichpt(1)) = current_subunit;
elseif strcmp(st, 'alt') % right click
    UD.mask_display(whichpt(2), whichpt(1)) = NaN;
end
set(axsubunit, 'UserData', UD);
paint_subunit_selection(axsubunit);

    function n = current_subunit
        a = get(gcf, 'UserData');
        n = get(a.uicontrols.slsubunit, 'Value');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
% end of callbacks
%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function that returns the number of subunits and the entire mask; is
% meaningful when there are subunits
function [num_units, Subunit_mask] = check_num_subunits()
% Abhishek added this function 08/2015,
a = get(gcf,'UserData');
SB = get(a.axeshandles.subunit, 'Userdata');
Subunit_mask = SB.mask_rex;
Subunit_mask(isnan(Subunit_mask)) = 0;
any_bkg_stixs = any(Subunit_mask(:) == 0);
num_units = numel(unique(Subunit_mask(:))) - any_bkg_stixs ; % 0 corresponds to no subunits
end

% Update STS, STCross, and nspikes
function plotnow = UpdateSTX(p, seed, nframes, mu, sigma, bkgndrgb, noisetype)
C = WhiteNoisethreshCodes;
a = get(gcf,'UserData');
stats = a.stats;
controls = a.uicontrols;
nframesback = str2double(get(controls.nframes,'String'));
gammaTable = stats.gammaTable;
invgamma = stats.invgammaTable;
nstixperside = stats.nstixperside;
controls = a.uicontrols;
[n_subunits,~] = check_num_subunits() ; % 0 corresponds to no subunits
if (n_subunits)
    nrandnums_perchannel = n_subunits;
else
    nrandnums_perchannel = nstixperside^2;
end

%*******************************************************************%
%     M = stats.M;
%     invM = stats.invM;
NGAMMASTEPS = size(invgamma,1);
stats.murgb = [gammaTable(bkgndrgb(1)+1,1);...
    gammaTable(bkgndrgb(2)+1,2);...
    gammaTable(bkgndrgb(3)+1,3)];
% GDLH 11/12/11 Doing this on every trial now.  Previously
% only if stats.murgb was Nan.  For some reason murgb was
% getting set to all zeros and that was causing an error in
% the plotting.
if (noisetype == 1)  % gun noise
    x = linspace(stats.gausslocut, stats.gausshicut, NGAMMASTEPS)';
    stats.gunnoise.mu = mu;          % Normalized intensity units
    stats.gunnoise.sigma = sigma;     % Normalized intensity units
    invnormcdf = zeros(NGAMMASTEPS, 3);
    for gun = 1:3
        invnormcdf(:,gun) = norminv(x)*sigma(gun)+mu(gun);
    end
    randnums = getEJrandnums(3*nrandnums_perchannel*nframes, seed);
    randnums = reshape(randnums, [nrandnums_perchannel*3, nframes]);
    for i = 1:3
        idxs = (1:nrandnums_perchannel)+nrandnums_perchannel*(i-1);
        randnums(idxs,:) = reshape(invnormcdf(randnums(idxs,:)+1,i),[length(idxs),nframes]);
    end
else % noisetype == 2 cone noise
    stats.conenoise.mu = mu;
    stats.conenoise.sigma = sigma;
    colordirlms = sign(fullfact([2,2,2])-1.5);
    randnums = getEJrandnums(nrandnums_perchannel*nframes, seed);
    randnums = mod(randnums, 8)+1;
    randnums = colordirlms(randnums,:);
    randnums = reshape(randnums, [nrandnums_perchannel, nframes, 3]);
    randnums = permute(randnums,[1 3 2]);
    randnums = reshape(randnums, [nrandnums_perchannel*3, nframes]);
    % Each column should be Ls followed by Ms followed by Ss for each frame.
end

t_stimon = p.times(find(p.events == C.STIMONCD, 1));
nspikesthistrial = [0 0];
for i = 1:2
    spiketimes = p.spikes{i}-t_stimon;
    if (any(spiketimes))
        frametimes = linspace(0, nframes*stats.msperframe, nframes)+(stats.msperframe/2)';
        spiketimes(spiketimes < nframesback*stats.msperframe) = [];
        spiketimes(spiketimes > frametimes(end)) = [];
        n = hist(spiketimes, frametimes);
        STCOVmex('init',{nrandnums_perchannel,3,nframesback});
        STCOVmex(randnums(:), n);
        out = STCOVmex('return');
        clear STCOVmex;        
        if (noisetype == 1)
            stats.gunnoise.spike{i}.STS = stats.gunnoise.spike{i}.STS+out{1};
            stats.gunnoise.spike{i}.STCross = stats.gunnoise.spike{i}.STCross+out{2};
            stats.gunnoise.spike{i}.nspikes = stats.gunnoise.spike{i}.nspikes+out{3};
        else
            stats.conenoise.spike{i}.STS = stats.conenoise.spike{i}.STS+out{1};
            stats.conenoise.spike{i}.STCross = stats.conenoise.spike{i}.STCross+out{2};
            stats.conenoise.spike{i}.nspikes = stats.conenoise.spike{i}.nspikes+out{3};
        end
    end
    nspikesthistrial(i) = length(spiketimes);
end
a.stats = stats;
set(gcf,'UserData',a);

plotnow = 0;
if (noisetype == get(a.uicontrols.whichnoisetype,'Value'))
    whichspike = get(a.uicontrols.whichspike,'Value');
    L1 = (noisetype == 1 & stats.gunnoise.spike{whichspike}.nspikes > 20);
    L2 = (noisetype == 2 & stats.conenoise.spike{whichspike}.nspikes > 20);
    if ((L1 | L2) && nspikesthistrial(whichspike) > 1)
        plotnow = 1;
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UDP communication functions below...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stopnow = DealWithMessage(msgSize)
global udpCom;
stopnow = 0;

message = pnet(udpCom.sock, 'read', msgSize, 'char');
if (strncmp(message,'return',6))
    a = dbstack;  % Check whether called from another function or from command line
    if (~strcmp(a(end).name, mfilename))
        stopnow = 1;
    end
end
try
    eval(message);
catch ME
    fprintf('Unknown message: %s\n', message);
    disp(getReport(ME));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Send an image to REX.  We have no idea when
% we will get here relative to the request since part of the
% main loop take a very long time and we poll for messages
% infrequently.


    % This function is called from 'calculateweights request' when all the
    % done flags of the current directions have been set to 1. This function finds new directions which the program 
    % will perform the adatptive staircase procedure.
    function out = PickNewWeights_old % Searches in all four quadrants, do it for STA vs PC
        global gl;
        out = 1;
        Lact = logical([gl.trialspecs.active]);
        oog_scaling_factor = gl.initparams.oogscale;
        if (all(isempty([gl.trialspecs(logical(Lact)).parentvertices])))
            % Getting ready for the 2nd round
            targetspikerates = unique([gl.trialspecs(logical(Lact)).targetspikerate]);
            labels = [1 2 3 4 1]; % ad-hoc implementation, need to automate this
            for ii = 1:length(targetspikerates)
                for jj = 1:numel(labels)-1
                    parentvertices = [labels(jj) labels(jj+1)];
                    tmp_oog = [gl.trialspecs(labels(jj)).outofgamut gl.trialspecs(labels(jj+1)).outofgamut];
                    tmp_vec = [gl.trialspecs(labels(jj)).weights*gl.trialspecs(labels(jj)).measuredthreshold;...
                        gl.trialspecs(labels(jj+1)).weights*gl.trialspecs(labels(jj+1)).measuredthreshold];
                    parentoog = any(tmp_oog>0);
                    if isequal(tmp_oog,[0 0])
                        weights = [1 1]*tmp_vec/2;
                    elseif isequal(tmp_oog,[1 0])
                        weights = [oog_scaling_factor 1]*tmp_vec/2;
                    elseif isequal(tmp_oog,[0 1])
                        weights = [1 oog_scaling_factor]*tmp_vec/2;
                    else
                        weights = [oog_scaling_factor oog_scaling_factor]*tmp_vec/2;
                    end
                    predthresh = norm(weights);
                    weights = weights/predthresh; % normalizing
                    SetUpTrialSpecs(weights, targetspikerates(ii),parentvertices, parentoog, predthresh);
                end
                labels = labels + 4;
            end
            out = 2;
            [gl.trialspecs(logical(Lact)).active] = deal(0);
        
        else 
            % Getting ready for the 3rd round
            labels = find(Lact==1);
            for ii = 1:numel(labels)
                grandparentvertices = gl.trialspecs(labels(ii)).parentvertices;
                threshratio = gl.trialspecs(labels(ii)).measuredthreshold/gl.trialspecs(labels(ii)).predictedthreshold;
                if (abs(log(threshratio)) > abs(log(1+gl.initparams.tolerance))) % split if not within the bounds of linear prediction
                    for jj = 1:2
                        parentvertices = [labels(ii) grandparentvertices(jj)];
                        tmp_oog = [gl.trialspecs(parentvertices(1)).outofgamut gl.trialspecs(parentvertices(2)).outofgamut];
                        if ~all(tmp_oog==1)
                            parentoog = any(tmp_oog>0);
                            tmp_vec = [gl.trialspecs(parentvertices(1)).weights*gl.trialspecs(parentvertices(1)).measuredthreshold;...
                                gl.trialspecs(parentvertices(2)).weights*gl.trialspecs(parentvertices(2)).measuredthreshold];
                            if isequal(tmp_oog,[0 0])
                                weights = [1 1]*tmp_vec/2;
                            elseif isequal(tmp_oog,[1 0])
                                weights = [oog_scaling_factor 1]*tmp_vec/2;
                            elseif isequal(tmp_oog,[0 1])
                                weights = [1 oog_scaling_factor]*tmp_vec/2;
                            end
                            predthresh = norm(weights);
                            weights = weights/predthresh; % normalizing
                            SetUpTrialSpecs(weights,gl.trialspecs(labels(ii)).targetspikerate,parentvertices,parentoog,predthresh);
                        end
                    end
                end
            end
            out = 3;
            [gl.trialspecs(logical(Lact)).active] = deal(0);
        end
    end
    function out = PickNewWeights_previous_func% Abhishek 07/16, Search in first quadrant, do it for subunits
        global gl;
        out = 1;
        Lact = logical([gl.trialspecs.active]);
        oog_scaling_factor = gl.initparams.oogscale;
        if (all(isempty([gl.trialspecs(logical(Lact)).parentvertices])))
            % Getting ready for the 2nd round
            targetspikerates = unique([gl.trialspecs(logical(Lact)).targetspikerate]);
            a = get(gcf,'UserData');
            basisvectypestrings = get(a.uicontrols.basisvectype,'String');
            basisvectypeval = get(a.uicontrols.basisvectype,'Value');
            if strcmp(basisvectypestrings(basisvectypeval),'subunits')
                labels = [1 2 3 4 5]; % ad-hoc implementation, need to automate thi
            else
                labels = [1 2 3];
            end
            
            for ii = 1:length(targetspikerates)
                for jj = 1:numel(labels)-1
                    parentvertices = [labels(jj) labels(jj+1)];
                    tmp_oog = [gl.trialspecs(labels(jj)).outofgamut gl.trialspecs(labels(jj+1)).outofgamut];
                    tmp_vec = [gl.trialspecs(labels(jj)).weights;...
                        gl.trialspecs(labels(jj+1)).weights];
                    parentoog = any(tmp_oog>0);
                    weights = [1 1]*tmp_vec/2;
                    predthresh = norm(weights);
                    weights = weights/predthresh; % normalizing
                    SetUpTrialSpecs(weights, targetspikerates(ii),parentvertices, parentoog, predthresh);
                end
                labels = labels + numel(labels);
            end
            out = 2;
            [gl.trialspecs(logical(Lact)).active] = deal(0);
        
        else 
            % Getting ready for the 3rd round
            labels = find(Lact==1);
            for ii = 1:numel(labels)
                grandparentvertices = gl.trialspecs(labels(ii)).parentvertices;
                for jj = 1:2
                    parentvertices = [labels(ii) grandparentvertices(jj)];
                    tmp_oog = [gl.trialspecs(parentvertices(1)).outofgamut gl.trialspecs(parentvertices(2)).outofgamut];
                    parentoog = any(tmp_oog>0);
                    tmp_vec = [gl.trialspecs(parentvertices(1)).weights;...
                        gl.trialspecs(parentvertices(2)).weights];
                    weights = [1 1]*tmp_vec/2;
                    predthresh = norm(weights);
                    weights = weights/predthresh; % normalizing
                    SetUpTrialSpecs(weights,gl.trialspecs(labels(ii)).targetspikerate,parentvertices,parentoog,predthresh);
                end 
            end
            out = 3;
            [gl.trialspecs(logical(Lact)).active] = deal(0);
        end
    end

   % This is the current PickNewWeights I am using, stick to this 
    function out = PickNewWeights % Abhishek 12/16, Optimising search along different directions in basis space
        global gl;
        out = 1;
        Lact = logical([gl.trialspecs.active]);
        oog_scaling_factor = gl.initparams.oogscale;
        if (all(isempty([gl.trialspecs(logical(Lact)).parentvertices])))
            % Getting ready for the 2nd round
            [gl.trialspecs(logical(Lact)).pending] = deal(0);
            targetspikerates = unique([gl.trialspecs(logical(Lact)).targetspikerate]);
            a = get(gcf,'UserData');
            basisvectypestrings = get(a.uicontrols.basisvectype,'String');
            basisvectypeval = get(a.uicontrols.basisvectype,'Value');
            if strcmp(basisvectypestrings(basisvectypeval),'subunits')
                oogpts = [gl.trialspecs.outofgamut];
                labels = [1 2 3 4 5 6 7 8]; % ad-hoc implementation, need to automate thi
                if all(oogpts(6:8)==0) % i.e. probe in the third quadrant if only all the intial searches are within gamut 
                    subtractnumber = 0;
                else 
                    subtractnumber = 3; 
                end
            else
                labels = [1 2 3 4 1]; % for STA vs PC1 search
                subtractnumber = 0;
            end
            simul_dirs_to_probe = str2num(get(a.uicontrols.simul_dirs_to_probe,'string'));
            for ii = 1:length(targetspikerates)
                dirs_to_be_used = [];
%                 keyboard;
                rem = mod(1:numel(labels)-1-subtractnumber,floor((numel(labels)-1-subtractnumber)/simul_dirs_to_probe));
                rem = (circshift(rem',ii-1))';
                for jj = 1:numel(labels)-1 - subtractnumber
                    parentvertices = [labels(jj) labels(jj+1)];
                    tmp_oog = [gl.trialspecs(labels(jj)).outofgamut gl.trialspecs(labels(jj+1)).outofgamut];
                    tmp_vec = [gl.trialspecs(labels(jj)).weights;...
                        gl.trialspecs(labels(jj+1)).weights];
                    parentoog = any(tmp_oog>0);
                    weights = [0.5 0.5]*tmp_vec/2;
                    predthresh = norm(weights);
                    weights = weights/predthresh; % normalizing
                    SetUpTrialSpecs(weights, targetspikerates(ii),parentvertices, parentoog, predthresh, rem(jj));
                    if rem(jj) == 0
                        gl.trialspecs(end).pending = 0;
                    else
                        gl.trialspecs(end).active = 0;
                    end
                    dirs_to_be_used = [dirs_to_be_used; length(gl.trialspecs)];
                end
                gl.dirs_to_be_used{ii} = dirs_to_be_used;
                labels = labels + (numel([gl.trialspecs.outofgamut])/length(targetspikerates));
                if strcmp(basisvectypestrings(basisvectypeval),'subunits')
                    if all(oogpts(8*ii+6:8*ii)==0)
                        subtractnumber = 0;
                    else
                        subtractnumber = 3;
                    end
                end
            end
            gl.max_rem = max(rem);
            gl.current_rem = 0;
            out = 2;
            [gl.trialspecs(logical(Lact)).active] = deal(0);
            plot_weight_vector();
        
        else 
            % Few checks before getting into the further rounds
            
            pending_idxs = logical([gl.trialspecs.pending] & ~[gl.trialspecs.active] & ~[gl.trialspecs.done] );
            just_finished_idxs = logical(~[gl.trialspecs.pending] & [gl.trialspecs.active] & [gl.trialspecs.done]);
            if ~isempty(find(pending_idxs))
%                 keyboard
                [gl.trialspecs(just_finished_idxs).active] = deal(0);
                % Check if the current remainder used is more than the max
                % remainder
                if gl.current_rem<gl.max_rem
                    gl.current_rem = gl.current_rem + 1;
                    dirs_to_probe_idxs = find([gl.trialspecs.remainder] == gl.current_rem);
                    [gl.trialspecs(dirs_to_probe_idxs).pending] = deal(0);
                    [gl.trialspecs(dirs_to_probe_idxs).active] = deal(1);
                    plot_weight_vector();
                else
                    PickNewWeights();
                end   
            else
                % Getting ready for the 3rd round
                [gl.trialspecs(logical(Lact)).active] = deal(0);
                targetspikerates = unique([gl.trialspecs.targetspikerate]);
                a = get(gcf,'UserData');
                simul_dirs_to_probe = str2num(get(a.uicontrols.simul_dirs_to_probe,'string'));
                for kk = 1: length(targetspikerates)
                    labels = gl.dirs_to_be_used{kk};
                    rem = mod(1:2*numel(labels),floor(numel(labels)*(2/simul_dirs_to_probe))); % complicated formula
                    rem = (circshift(rem',kk-1))';
                    dirs_to_be_used = [];
                    for ii = 1:numel(labels)
                        grandparentvertices = gl.trialspecs(labels(ii)).parentvertices;
                        for jj = 1:2
                            parentvertices = [labels(ii) grandparentvertices(jj)];
                            tmp_oog = [gl.trialspecs(parentvertices(1)).outofgamut gl.trialspecs(parentvertices(2)).outofgamut];
                            parentoog = any(tmp_oog>0);
                            tmp_vec = [gl.trialspecs(parentvertices(1)).weights;...
                                gl.trialspecs(parentvertices(2)).weights];
                            weights = [0.5 0.5]*tmp_vec/2;
                            predthresh = norm(weights);
                            weights = weights/predthresh; % normalizing
                            SetUpTrialSpecs(weights,gl.trialspecs(labels(ii)).targetspikerate,parentvertices,parentoog,predthresh,rem(2*(ii-1)+jj));
                            if rem(2*(ii-1)+jj) == 0
                                gl.trialspecs(end).pending = 0;
                            else
                                gl.trialspecs(end).active = 0;
                            end
                            dirs_to_be_used = [dirs_to_be_used; length(gl.trialspecs)];
                        end
                    end
                    gl.dirs_to_be_used{kk} = [];
                    gl.dirs_to_be_used{kk} = dirs_to_be_used;
                    gl.max_rem = max(rem);
                    gl.current_rem = 0;
                end
                plot_weight_vector();
            end
            out = 3;
        end
    end

    function basissendingflagrequest()
        % Sends the basissendingflag to REX
        a = get(gcf,'UserData');
        b = get(a.axeshandles.synthImage,'UserData');
        sendToRex(udpCom, b.basissendingflag, 'int', 'basissendingflagrequest();');
    end

    function basisvecrequest()  
        % Sends the basis vectors to REX when the basissending flag is set
        % to 1
        a = get(gcf,'UserData');
        b = get(a.axeshandles.synthImage,'UserData');
        if b.basissendingflag
            nstixperside = a.stats.nstixperside;
            im = zeros(nstixperside,nstixperside,3,3);
            im(:,:,:,1) = b.basisvectosend{1};
            im(:,:,:,2) = b.basisvectosend{2};
            im(:,:,:,3) = b.muimage;
            sendToRex(udpCom, im, 'double', 'basisvecrequest();');
            b.basissendingflag = 0;
            set(a.axeshandles.synthImage,'UserData',b);
            set(gcf,'UserData',a);
        end
    end
 
    function latencyrequest()
        % sends the latency to the REX so that it gets dropped into the file
        a = get(gcf,'UserData');
        t_offset = str2num(get(a.uicontrols.latency,'String'));
        if ~isempty(t_offset)
            sendToRex(udpCom, t_offset , 'double', 'latencyrequest();');
        end
    end

    function targetspikeraterequest()
        % sends the targetspikerate to the REX so that it gets dropped into the file
        global gl;
        if ~isempty(gl.trialspecidx)
            sendToRex(udpCom, gl.trialspecs(gl.trialspecidx).targetspikerate , 'double', 'targetspikeraterequest();');
        end
    end
    
    function calculateweightsrequest()
        % Is called the inittrial state in REX, initiates the calculation
        % of the weights (direction) that is eventually passed to the REX
        % through weightsrequest function
        a = get(gcf,'UserData');
        b = get(a.axeshandles.synthImage,'UserData');
        if ~isempty(b.weightstosend)
            nstixperside = a.stats.nstixperside;
            if get(a.uicontrols.NT,'Value');
                % Enter the Neurothresh mode
                % fields 'stepsize' and 'norms are being changed in this section
                global gl
                idxs = find([gl.trialspecs.active] & ~[gl.trialspecs.done]);
                if (isempty(idxs))
                    
                    out = PickNewWeights();% subunits & STA vs PC
                    idxs = find([gl.trialspecs.active] & ~[gl.trialspecs.done]);
                end
                
                if (isempty(idxs))
                    gl.wtsstopsendingflag = 1;
                    gl.trialspecidx = [];
                    set(a.axeshandles.synthImage,'UserData',b); % debugging
                    set(gcf,'UserData',a);
                else
                    % creating the image using the weights
                    gl.trialspecidx = idxs(unidrnd(length(idxs)));
                    %                     [gl.trialspecidx gl.trialspecs(gl.trialspecidx).targetspikerate]
                    weights = gl.trialspecs(gl.trialspecidx).weights;
                    if (abs(weights(1))>0 | abs(weights(2))>0)
                        im = weights(1)*(b.basisvectosend{1}) + weights(2)*(b.basisvectosend{2});
                    else
                        im = zeros(nstixperside,nstixperside,3);
                        gl.trialspecs(gl.trialspecidx).weights = [0 0];
                    end
                    
                    trial = gl.trialspecs(gl.trialspecidx).trial_num; % how many trials for gl.trialspecidx have been completed
                    thresh = gl.trialspecs(gl.trialspecidx).targetspikerate;
                    
                    % Abhishek - 07/16, Added currenttargetFR_disp
                    set(a.uicontrols.currenttargetFR_disp,'String',num2str(thresh));
                        
                    s1=0; s2=0;
                    if (trial < 2)
                        stepsize = gl.initparams.stepsize;
                        gl.trialspecs(gl.trialspecidx).stepsize = [gl.trialspecs(gl.trialspecidx).stepsize; gl.initparams.stepsize];
                    else
                        s2 = gl.trialspecs(gl.trialspecidx).spikerates(end)- thresh;
                        s1 = gl.trialspecs(gl.trialspecidx).spikerates(end-1)- thresh;
                        stim_rev_sgn = sign(s2)-sign(s1);
                        if(abs(stim_rev_sgn))
                            stepsize = gl.trialspecs(gl.trialspecidx).stepsize(end)*gl.initparams.stepsizescale;
                        else
                            stepsize = gl.trialspecs(gl.trialspecidx).stepsize(end);
                        end
                        gl.trialspecs(gl.trialspecidx).stepsize = [gl.trialspecs(gl.trialspecidx).stepsize; stepsize];
                    end
                    
                    if (trial >= 1) % Atlest 1 trial has been completed
                        try
                            if (gl.trialspecs(gl.trialspecidx).spikerates(end) > thresh)
                                norm = gl.trialspecs(gl.trialspecidx).norms(end) * (1-stepsize); % Decreasing contrast
                            elseif (gl.trialspecs(gl.trialspecidx).spikerates(end) <= thresh)
                                norm = gl.trialspecs(gl.trialspecidx).norms(end) * (1+stepsize); % Increasing contrast
                            end
                        catch
                            disp('Error in calculateweightsrequest');
                            keyboard
                        end
                    else
                        if isempty(gl.trialspecs(gl.trialspecidx).predictedthreshold)
                            norm = 0.5; % Value of the norm to be used for the first round
                        else
                            norm = gl.trialspecs(gl.trialspecidx).predictedthreshold;
                        end
                    end
                    gl.trialspecs(gl.trialspecidx).norms = [gl.trialspecs(gl.trialspecidx).norms; norm];
                    im = im*gl.trialspecs(gl.trialspecidx).norms(end) + b.muimage;
                    if (any(im(:)>1) | any(im(:)<0))
                        disp(['OOG detected! idx:', num2str(gl.trialspecidx)]);
                        gl.trialspecs(gl.trialspecidx).norms(end) = [];
                        gl.trialspecs(gl.trialspecidx).stepsize(end) = [];
                        gl.trialspecs(gl.trialspecidx).done = 1;
                        if trial > 0
                            gl.trialspecs(gl.trialspecidx).measuredthreshold = gl.trialspecs(gl.trialspecidx).norms(end);
                        end
                        gl.trialspecs(gl.trialspecidx).outofgamut = 1;
                        calculateweightsrequest();
                    end
                    set(a.axeshandles.synthImage,'UserData',b);
                    set(gcf,'UserData',a);
                end
            end 
        end
    end

    function parentverticesrequest()
        % Send the reversalflag to REX so that it gets dropped into the file
        global gl
        if ~isempty(gl.trialspecidx)
            sendToRex(udpCom, gl.trialspecs(gl.trialspecidx).parentvertices, 'integer', 'parentverticesrequest();');
        end
    end

    function reversalflagrequest()
        % Send the reversalflag to REX so that it gets dropped into the file
        global gl
        if ~isempty(gl.trialspecidx)
            sendToRex(udpCom, gl.trialspecs(gl.trialspecidx).reversal_num, 'integer', 'reversalflagrequest();');
        end
    end

    function weightsidxrequest()
        % Sends the weightsidxrequest to REX so that it gets dropped into the file
        a = get(gcf,'UserData');
        global gl
        if get(a.uicontrols.NT,'Value')
            if ~isempty(gl.trialspecidx)
                sendToRex(udpCom, gl.trialspecidx, 'integer', 'weightsidxrequest();');
            end
        else
            sendToRex(udpCom, gl.trialspecidx, 'integer', 'weightsidxrequest();');
        end   
    end

    function weightsrequest()
        % Send the weights to the REX that gets dropped into the file.
        global gl
        a = get(gcf,'UserData');
        b = get(a.axeshandles.synthImage,'UserData');
        if get(a.uicontrols.NT,'Value')
            % Neurothresh Mode
            if ~isempty(gl.trialspecidx)
                b.wts_NT = gl.trialspecs(gl.trialspecidx).weights * gl.trialspecs(gl.trialspecidx).norms(end);
                sendToRex(udpCom,b.wts_NT , 'double', 'weightsrequest();');  
            else
                if (gl.wtsstopsendingflag)
                    sendToRex(udpCom,[10000 10000], 'double', 'weightsrequest();');
                    disp('Last trial ');
                end
            end
        elseif ~get(a.uicontrols.NT,'Value') & ~isempty(b.weightstosend)
            % SynthImageBattery mode
            disp('Not NT mode');
            b.wts_NT = b.weightstosend{1}; % norm value is 1;
            b.weightstosend(1) = [];
            sendToRex(udpCom,b.wts_NT, 'double', 'weightsrequest();');
        end
        set(a.axeshandles.synthImage,'UserData',b);
        set(gcf,'UserData',a);
    end


    function mask = maskrequest()
        % Send the stixel mask to REX on request.
        UD = get(gcf, 'UserData');
        if get(UD.uicontrols.updatemask, 'Value')
            UDsubunit = get(UD.axeshandles.subunit, 'UserData');
            mask = UDsubunit.mask_display;
            set(UD.uicontrols.updatemask, 'Value',0)
            tmp = all(isnan(mask(:)));
            if (~tmp)
                UDsubunit.mask_rex = mask;
                set(UD.axeshandles.subunit, 'UserData', UDsubunit);
                % Abhishek made a change to this function
                nframesback = str2double(get(UD.uicontrols.nframes,'String'));
                [nrandnums_perchannel,~] =  check_num_subunits(); % 0 corresponds to no subunits
                disp('Mask Update on');
                for i = 1:2
                    UD.stats.gunnoise.spike{i}.STS = zeros([3*nrandnums_perchannel nframesback]);
                    UD.stats.gunnoise.spike{i}.nspikes = 0;
                    UD.stats.conenoise.spike{i}.STS = zeros([3*nrandnums_perchannel nframesback]);
                    UD.stats.conenoise.spike{i}.nspikes = 0;
                    UD.stats.gunnoise.spike{i}.STCross = zeros([(3*nrandnums_perchannel)^2 nframesback]);
                    UD.stats.conenoise.spie{i}.STCross = zeros([(3*nrandnums_perchannel)^2 nframesback]);
                end
                %InitSynthImage(UD.axeshandles.synthImage)
                b = get(UD.axeshandles.synthImage,'UserData');
                b.STA = []; b.PC = []; b.localimage = [];
                set(UD.axeshandles.synthImage,'UserData',b);
                set(gcf,'Userdata',UD);
            else
                disp('No stixel selected');
            end
            
            %***********************************************************  
        else
            mask = getfield(get(UD.axeshandles.subunit, 'UserData'), 'mask_rex');
        end
        mask(isnan(mask)) = 0;
        sendToRex(udpCom, mask(:), 'integer', 'maskrequest();');
    end


    function nudgerequest()
        % Send a displacement of the stimulus [(x,y) in stixels] to REX on request.
        a = get(gcf, 'UserData');
        axeshandles = a.axeshandles;
        b = get(axeshandles.subunit,'UserData');
        if (~isempty(b.nudge))
            sendToRex(udpCom, b.nudge, 'double', 'nudgerequest();'); % double to allow nudging in fractions of stixels
            b.nudge = [];
            set(a.axeshandles.subunit, 'UserData', b);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting functions below...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateQueueText(p)
% Updates the Events and the spikes that gets displayed after the end of each trial
uicontrols = getfield(get(gcf,'UserData'),'uicontrols');
nevents = size(p.events,1);
nspikes = [size(p.spikes{1},1), size(p.spikes{2},1)];  % '{1}' and '{2}' are spike numbers
if (nspikes(2) == 0)
    str = ['E: ',num2str(nevents),'   S:',num2str(nspikes(1))];
else
    str = ['E:',num2str(nevents),' S:',num2str(nspikes(1)),', ',num2str(nspikes(2))];
end
set(uicontrols.eventqueuelength,'string',str);
drawnow;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting the active and the completed directions - Abhishek 12/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_weight_vector()
global gl
a = get(gcf,'UserData');
active_idxs = logical([gl.trialspecs.active] & ~[gl.trialspecs.done]);
finished_idxs = logical(~[gl.trialspecs.active] & [gl.trialspecs.done]);
finished_idxs_oog = logical(~[gl.trialspecs.active] & [gl.trialspecs.done] & [gl.trialspecs.outofgamut]);
targetspikerates_all = [gl.trialspecs.targetspikerate];
targetspikerates = unique(targetspikerates_all);
ax = [a.axeshandles.dir_TFR1; a.axeshandles.dir_TFR2];
for ii = 1:numel(targetspikerates)
    axes(ax(ii)); grid off;
    % Use compass statement to plot the weight vectors
    % For finished points
    finished = finished_idxs & logical(targetspikerates_all == targetspikerates(ii));    
    finished_weights = reshape([gl.trialspecs(finished).weights],2,[]);
    measuredthreshold = [gl.trialspecs(finished).measuredthreshold];
    isoresp_pts =  repmat(measuredthreshold,[2 1]).* finished_weights;
    plot(isoresp_pts(1,:),isoresp_pts(2,:),'bo','MarkerSize',2,'MarkerFaceColor','b'); hold on;
    
    % For oog points
    finishedoog = finished_idxs_oog & logical(targetspikerates_all == targetspikerates(ii));    
    finished_weightsoog = reshape([gl.trialspecs(finishedoog).weights],2,[]);
    measuredthresholdoog = [gl.trialspecs(finishedoog).measuredthreshold];
    isoresp_ptsoog =  repmat(measuredthresholdoog,[2 1]).* finished_weightsoog;
    plot(isoresp_ptsoog(1,:),isoresp_ptsoog(2,:),'go','MarkerSize',2,'MarkerFaceColor','g'); hold on;
    
    % For active points
    active = active_idxs & logical(targetspikerates_all == targetspikerates(ii));
    active_weights = reshape([gl.trialspecs(active).weights],2,[]);
    plot(active_weights(1,:),active_weights(2,:),'ro','MarkerSize',3,'MarkerFaceColor','r');
    hold off;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create and display a linear combination of STA and PC.
% Invoked when we adjust the STA or PC slider, or when
% we select or deselect at STA or PC frame.  Since these
% are all Callbacks, they can't co-occur accidentally.
function UpdateSynthImage(~,~)
% Updates the SynthImage 
a = get(gcf,'UserData');
mu = a.stats.murgb;
nstixperside = a.stats.nstixperside;
muimage = cat(3, repmat(mu(1),[nstixperside nstixperside]),...
    repmat(mu(2),[nstixperside nstixperside]),...
    repmat(mu(3),[nstixperside nstixperside]));

uicontrols = getfield(get(gcf,'UserData'),'uicontrols');
STAvalue = get(uicontrols.STAslider,'value');
PCvalue = get(uicontrols.PCslider,'value');
axeshandles = getfield(get(gcf,'UserData'),'axeshandles');
b = get(axeshandles.synthImage,'UserData');
STAimage = b.STA;
PCimage = b.PC;

if (isempty(STAimage))
    STAimage = muimage;
end
if (isempty(PCimage))
    PCimage = muimage;
end

if (abs(STAvalue) > 0.01 || abs(PCvalue) > 0.01)
    w = [STAvalue PCvalue];
    w = w/(eps + abs(STAvalue)+abs(PCvalue));
    im = w(1)*(STAimage-muimage) + w(2)*(PCimage-muimage) + muimage;
else
    w = [0 0];
    im = muimage;
end

b.localweights = w;
axes(axeshandles.synthImage);
b.localimage = im;
try
    h = image(im);
    set(h,'ButtonDownFcn',{@SynthImageCallback, axeshandles.synthImage});
catch
    disp('replay image error');
    image(muimage);
end
set(axeshandles.synthImage, 'XTick',[],'YTick',[]);
set(a.uicontrols.STAcoefedit,'string',num2str(STAvalue,2));
set(a.uicontrols.PCcoefedit,'string',num2str(PCvalue,2));
set(axeshandles.synthImage,'Userdata',b);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the frames of the STA and PC
function PlotSTA(noisetype)
a = get(gcf,'UserData');
stats = getfield(get(gcf,'UserData'),'stats');
controls = getfield(get(gcf,'UserData'),'uicontrols');
axeshandles = getfield(get(gcf,'UserData'),'axeshandles');
pixelmask = getfield(get(axeshandles.pixelmask,'UserData'),'pixelmask');
nframesback = str2double(get(controls.nframes,'String'));

[n_subunits,~] = check_num_subunits(); % 0 corresponds to no subunits
if (n_subunits)
    nrandnums_perchannel = n_subunits;
else
    nrandnums_perchannel = stats.nstixperside^2;
end

mu = stats.murgb;
% if (all(mu == 0))
%     disp('No background rgb set');
%     keyboard
% end
muvect = reshape(repmat(mu',nrandnums_perchannel,1),nrandnums_perchannel*3,1);
mumat = repmat(muvect,1,nframesback);
whichspike = get(controls.whichspike,'Value');

if (get(controls.whichnoisetype,'Value') == 1) % gun noise
    STS = stats.gunnoise.spike{whichspike}.STS;
    STCross = stats.gunnoise.spike{whichspike}.STCross;
    n = stats.gunnoise.spike{whichspike}.nspikes;
    noisetypestr = 'gun';
else % cone noise
    STS = stats.conenoise.spike{whichspike}.STS;
    STCross = stats.conenoise.spike{whichspike}.STCross;
    n = stats.conenoise.spike{whichspike}.nspikes;
    noisetypestr = 'cone';
end
STAs = STS./n;

po_strings = get(controls.projortho,'String');
po_val = get(controls.projortho,'Value');
PROJORTHO = po_strings(po_val);
opts.issym = 'true';
opts.disp = 0;
opts.isreal = 'true';

axes(axeshandles.text(1));
cla;
text(-200,10,['nspikes(',noisetypestr,' ',num2str(whichspike),') = ',num2str(n)]);

% Iterating to get PCs
PCs = zeros(size(STAs));
if n_subunits
    an_mask = [];
    an_mask = zeros(nrandnums_perchannel, 1);   % Pixel mask. Someday make a tool to make this non-zero.
    Lmask = logical(repmat(~an_mask(:),[3 1]));
    SB = get(a.axeshandles.subunit, 'Userdata');
    Subunit_mask = SB.mask_rex; % declaring Subunit_mask here
    Subunit_mask(isnan(Subunit_mask)) = 0;
else
    Lmask = logical(repmat(pixelmask(:),3,1));
end

if (sum(Lmask) == 0)
    disp('No pixels selected!  Not plotting anything.');
    return;
end

% Don't look at the smallest PC until we have enough spikes.
if (n < 3*sum(Lmask))
    set(controls.smallpc,'Value',0);
end
SMALLPC = get(controls.smallpc,'Value');


for ik = 1:nframesback
    STA = STAs(:,ik);
    tmp = STS(:,ik)*STS(:,ik)';
    STCframe = (n.*STCross(:,ik)-tmp(:))/(n*(n-1));
    STCframe = reshape(STCframe,[sqrt(size(STCframe,1)) sqrt(size(STCframe,1))]);
    
    subSTCframe = STCframe(Lmask, Lmask);
    if (strcmp(PROJORTHO,'PC'))
        % Projecting the STC orthogonal to the (masked) STA
        subSTA = STA(Lmask);
        P = eye(size(subSTCframe))-subSTA*inv(subSTA'*subSTA)*subSTA';
        subSTCframe = P*subSTCframe*P';
    end
    if (isempty(get(axeshandles.PC(ik),'UserData')))  % Initial eigenvector guess
        set(axeshandles.PC(ik),'UserData',normrnd(0,1,3*nrandnums_perchannel,1));
    end
    initialguess = get(axeshandles.PC(ik),'UserData');  % Use previous answer as initial guess
    opts.v0 = initialguess(Lmask);
    
    if (SMALLPC)
        if (strcmp(PROJORTHO,'PC')) % Take penultimate small PC if projecting out STA
            try
                [tmp,d] = eigs(subSTCframe,2,'SR', opts);
                d = diag(d);
                tmp = tmp(:,d == max(d));
            catch
                tmp = zeros(length(tmp),1);
            end
        else
            [tmp,~] = eigs(subSTCframe,1,'SM', opts);
        end
    else
        [tmp,~] = eigs(subSTCframe,1,'LM', opts);
    end
    v = double(Lmask);
    v(Lmask) = tmp;  % Inserting the 0's for the masked pixels
    
    if (initialguess'*v < 0)  % Trying to keep the sign from flipping
        v = -v;
    end
    PCs(:,ik) = v;
    set(axeshandles.PC(ik),'UserData',v);  % Before any normalization
    
end

if (strcmp(PROJORTHO,'STA'))  % Projecting STA orthogonal to PC1
    for iii = 1:nframesback
        proj = STAs(:,iii)'*PCs(:,iii);
        STAs(:,iii) = STAs(:,iii)-(proj.*PCs(:,iii));
    end
end

if (strcmp(noisetypestr,'cone'))
    STAs = ConvertlmsTorgb(STAs, stats.M, stats.murgb, stats.conenoise.sigma);
    PCs = ConvertlmsTorgb(PCs, stats.M, stats.murgb, stats.conenoise.sigma);
end

noiseSTAsd = mean(std(reshape(STAs(:,1),[nrandnums_perchannel, 3])));

if(n_subunits)
    % Normalize STAs and PCs separately
    normfactor = FindImNormFactor(STAs, mu,  nrandnums_perchannel);
    STAs = normfactor*STAs+mumat;
    normfactor = FindImNormFactor(PCs, mu,  nrandnums_perchannel);
    PCs = normfactor*PCs+mumat;
else
    if get(controls.contrastnorm,'Value')    % Yoke the noise contrasts
        PCs = NormToEdge(PCs, noiseSTAsd);
        normfactor = FindImNormFactor([STAs PCs], mu, nrandnums_perchannel);
        STAs = normfactor*STAs+mumat;
        PCs = normfactor*PCs+mumat;
    else
        normfactor = FindImNormFactor(STAs, mu,  nrandnums_perchannel);
        STAs = normfactor*STAs+mumat;
        normfactor = FindImNormFactor(PCs, mu,  nrandnums_perchannel);
        PCs = normfactor*PCs+mumat;
    end
end


% Debugging - looking at orthogonalization btn STAs and PCs
% (STAs-mumat)'*(PCs-mumat)
for jj = 1:nframesback
    axes(axeshandles.STA(jj));
    xcolor = get(gca,'Xcolor'); ycolor = get(gca,'Ycolor');
    cla;
    if (n_subunits)
        STAframe = expand_vector(STAs(:,jj),nrandnums_perchannel,Subunit_mask(:));
        STAframe = reshape(STAframe, [stats.nstixperside,stats.nstixperside,3]);
    else
        STAframe = reshape(STAs(:,jj), [stats.nstixperside,stats.nstixperside,3]);
    end
    h = image(STAframe);
    set(h,'ButtonDownFcn',{@ImageCallback, gca, STAframe});
    set(gca,'XTick',[],'YTick',[],'Xcolor',xcolor,'Ycolor',ycolor);
    set(gca,'UserData',STAframe);
    title(num2str(-(jj-1)));
    
    if all(xcolor == [1 1 0])
        update_subunit_axes(STAframe);
    end
    axes(axeshandles.PC(jj));
    xcolor = get(gca,'Xcolor'); ycolor = get(gca,'Ycolor');
    cla;
    if (n_subunits)
        PCframe = expand_vector(PCs(:,jj),nrandnums_perchannel,Subunit_mask(:));
        PCframe = reshape(PCframe,[stats.nstixperside, stats.nstixperside, 3]);
    else
        PCframe = reshape(PCs(:,jj),[stats.nstixperside, stats.nstixperside, 3]);
    end
    h = image(PCframe);
    set(h,'ButtonDownFcn',{@ImageCallback, gca, PCframe});
    set(gca, 'XTick',[],'YTick',[],'Xcolor',xcolor,'Ycolor',ycolor);

    axis image;
end

    function out = NormToEdge(in, noisestd)
        % Normalizing the PCs individually by making sure that the standard
        % deviation of the pixels around the edge is the same as the standard
        % deviation of the STA in the first frame which is presumably just
        % noise.  Practically, what this means is the the edge pixels should
        % not be impinging on the RF.
        % Taking the mean in noiseSTAsd because we're assuming equal
        % standard deviations for each gun for now.
        stats = getfield(get(gcf,'UserData'),'stats');
        template = reshape(1:stats.nstixperside^2,stats.nstixperside,stats.nstixperside);
        edgepixels = [template(:,1); template(1,2:end-1)'; template(end,2:end-1)'; template(:,end)];
        edgepixelidxs = [edgepixels; edgepixels+stats.nstixperside^2; edgepixels+2*(stats.nstixperside^2)];
        PCelements = in(edgepixelidxs,:);
        PCelements(all(PCelements == 0,2),:) = [];  % getting rid of the masked pixels
        PCsds = std(PCelements);    % One std calculated per PC
        out = in.*repmat(noisestd./PCsds, size(in,1),1);
    end

    function reformed_vec = expand_vector(orig_vec,subunits,mask)
        % This is a function that expands the (3*number of subunits) dimensional vector into the 300 dimensional vector using mask
        a = get(gcf,'UserData');
        n_elements = a.stats.nstixperside^2;
        reformed_vec = zeros(3*n_elements,1);
        for ii = 1:subunits
            vec_r  = mask;
            vec_g  = mask;
            vec_b  = mask;
            vec_r(~(vec_r == ii)) = 0; vec_r(vec_r == ii) = orig_vec(ii,1); % R
            vec_g(~(vec_g == ii)) = 0; vec_g(vec_g == ii) = orig_vec(subunits+ii,1); %G
            vec_b(~(vec_b == ii)) = 0; vec_b(vec_b == ii) = orig_vec(2*subunits+ii,1); % B
            vec = [vec_r;vec_g;vec_b];
            reformed_vec = reformed_vec + vec;
            clear vec vec_r vec_g vec_b;
        end
        
        bkg_r = zeros(n_elements,1); bkg_g = zeros(n_elements,1); bkg_b = zeros(n_elements,1);
        bkg_r(mask(:)==0) = a.stats.murgb(1);
        bkg_g(mask(:)==0) = a.stats.murgb(2);
        bkg_b(mask(:)==0) = a.stats.murgb(3);
        reformed_vec = reformed_vec + [bkg_r; bkg_g; bkg_b];
        clear bkg_r bkg_g bkg_b
    end

% Converting from a matrix of l,m,s differences to rgbs intensities for plotting
    function out = ConvertlmsTorgb(in, M, bkgndrgb, sigma)
        
        invM = inv(M);
        nstix = size(in,1)/3;
        bkgndlms = M*bkgndrgb;
        for i = 1:3
            idxs = (1:nstix)+nstix*(i-1);
            in(idxs,:) = in(idxs,:)*sigma(i)+bkgndlms(i);
        end
        giganticM = zeros(nstix,nstix);
        for i = 1:3
            idx1 = (1:nstix)+nstix*(i-1);
            for j = 1:3
                idx2 = (1:nstix)+nstix*(j-1);
                giganticM(idx1,idx2) = diag(repmat(invM(i,j), nstix, 1));
            end
        end
        out = giganticM*in;
        for i = 1:3
            idxs = (1:nstix)+nstix*(i-1);
            out(idxs,:) = out(idxs,:)-bkgndrgb(i);
        end
    end
end

function normfactor = FindImNormFactor(imagevectors, mu, nrandnums_perchannel)
% Finds a factor which normalizes a set of images so that none of the
% pixels are < 0 or > 1.  The "imagevectors" argument should be a
% matrix in which each column is an image, the first N/3 rows are
% the red, the second N/3 rows are the green, and the third N/3
% rows are the blue. The "mu" argument should be a 3x1 vector of
% mean gun intensities (around which we should scale).
stats = getfield(get(gcf,'UserData'),'stats');
rowidxs = reshape(1:3*nrandnums_perchannel,[nrandnums_perchannel 3]);

maxes = []; mins = [];
for i = 1:3
    maxes = [maxes; max(max(imagevectors(rowidxs(:,i),:)))];
    mins = [mins; min(min(imagevectors(rowidxs(:,i),:)))];
end
potentialnormfactors = abs([(1-mu-eps)./maxes; (-mu+eps)./mins;]);
% 'eps' in above line is a kludge that is required for avoiding
% out of bounds errors.
potentialnormfactors(potentialnormfactors < 0) = []; % if min > mu or max < mu
normfactor = min(potentialnormfactors);
end
% Potential problem here - both this function and PixelMaskCallback
% destructively modify the UserData field of the pixelmask axes.
% There could be a collision, and I'm not entirely sure what would happen.
% Trying to fix this - GDLH 2/19/08
function UpdatePixelMask(~,~)
stats = getfield(get(gcf,'UserData'),'stats');
axeshandles = getfield(get(gcf,'UserData'),'axeshandles');
controls = getfield(get(gcf,'UserData'),'uicontrols');
alpha = get(controls.alphaslider,'Value');
set(controls.alphatext,'String',num2str(alpha));
sigteststrings = get(controls.whichsigtest,'String');
whichsigtest = sigteststrings(get(controls.whichsigtest,'Value'));
whichspike = get(controls.whichspike,'Value');
noisestrings = get(controls.whichnoisetype,'String');
whichnoise = noisestrings(get(controls.whichnoisetype,'Value'));
hidx = WhichFrameSelected(axeshandles.STA);
b = get(axeshandles.pixelmask,'UserData');
if (isempty(b.pixelmask) && (stats.nstixperside > 0))  % Initializing image
    %EnablePixelMask(0); % not sure this is important
    b.pixelmask = ones(stats.nstixperside, stats.nstixperside);
    set(axeshandles.pixelmask,'UserData',b);
    %EnablePixelMask(1); % not sure this is important
    return;
end
if (~isempty(hidx))
    if (strcmp(whichnoise, 'gun'))
        if (isfield(stats,'gunnoise'))
            n = stats.gunnoise.spike{whichspike}.nspikes;
        else
            n = 0;
        end
        % Calculating a variance normalization factor to
        % compensate for the fact that our Gaussian is truncated.
        % Should really have one correction factor per dimension...
        NPOINTS = 65536;
        x = linspace(stats.gausslocut,stats.gausshicut,NPOINTS);
        Fx = norminv(x)*stats.gunnoise.sigma(1);
        sigmacorrectionfactors = std(Fx)./stats.gunnoise.sigma;
        %%%%%%%%%
        STA = stats.gunnoise.spike{whichspike}.STS(:,hidx)/n;
        STCross = stats.gunnoise.spike{whichspike}.STCross(:,hidx);
        STCross = reshape(STCross,[3*stats.nstixperside^2, 3*stats.nstixperside^2]);
        STS2 = diag(STCross);
        correctedsigmavect = reshape(repmat((stats.gunnoise.sigma.*sigmacorrectionfactors)',stats.nstixperside^2,1),stats.nstixperside^2*3,1);
    else
        if (isfield(stats,'conenoise'))
            n = stats.conenoise.spike{whichspike}.nspikes;
        else
            n = 0;
        end
        STA = stats.conenoise.spike{whichspike}.STS(:,hidx)/n;
        STCross = stats.conenoise.spike{whichspike}.STCross(:,hidx);
        STCross = reshape(STCross,[3*stats.nstixperside^2, 3*stats.nstixperside^2]);
        STS2 = diag(STCross);
        correctedsigmavect = reshape(repmat([1 1 1]',stats.nstixperside^2,1),stats.nstixperside^2*3,1);
        
    end
    
    pmat = [];
    if (strcmp(whichsigtest,'Mean'))
        Z = sqrt(n)*STA./correctedsigmavect;
        Z2 = reshape(Z.^2,[stats.nstixperside stats.nstixperside 3]);
        Z2 = sum(Z2,3);
        pmat = chi2cdf(Z2,3);
        %[mean(Z) std(Z)]
    elseif (strncmp(whichsigtest,'STA_',4'))
        % Significance test for each phosphor channel - Abhishek, 12/17
        Z = sqrt(n)*STA./correctedsigmavect;
        Z = reshape(Z,[stats.nstixperside stats.nstixperside 3]);
        whichgun = whichsigtest{end}(end);
%         keyboard
        switch(whichgun)
            case 'r'
                Z = squeeze(Z(:,:,1));
            case 'g'
                Z = squeeze(Z(:,:,2));
            case 'b'
                Z = squeeze(Z(:,:,3));
        end
        pmat = normcdf(Z,0,1);

    elseif (strcmp(whichsigtest,'Var'))
        % Doing the variance test
        s2 = STS2./correctedsigmavect.^2;
        s2mat = reshape(s2,[stats.nstixperside stats.nstixperside 3]);
        s2mat = sum(s2mat,3);
        pmat = chi2cdf(s2mat, 3*n);
        %[mean(s2)/n std(s2)]
        
        % Below is historic stuff that uses the non-central chisquare
        % distribution.  It also works, but is a lot less efficient.
        %s2 = STS2./correctedsigmavect.^2;
        %lambda = n*(muvect./correctedsigmavect).^2;
        %pmat = ncx2cdf(s2,n,lambda)  % Amazingly, this seems to work
        %s2mat = reshape(s2,[stats.nstixperside stats.nstixperside 3]);
        %s2mat = sum(s2mat,3);
        %lambdamat = reshape(lambda,[stats.nstixperside stats.nstixperside 3]);
        %lambdamat = sum(lambdamat,3);
        %pmat = ncx2cdf(s2mat,3*n,lambdamat);
    end
    im = repmat(b.pixelmask-0.5,[1 1 3])/4;
    if (~isempty(pmat))
        im(:,:,1) = im(:,:,1) - ((pmat < alpha)-0.5)/5;
        im(:,:,1) = im(:,:,1) + ((pmat > 1-alpha)-0.5)/5;
    end
    im = im+0.5;
    
    % This is weird, but making an image is erasing the contents of the
    % axes UserField!  It's ugly too.  Get to the bottom of this.
    EnablePixelMask(0); % not sure this is important
    b = get(axeshandles.pixelmask,'UserData');
    h = image(im,'parent',axeshandles.pixelmask);
    set(h,'ButtonDownFcn',@PixelMaskCallback);
    set(axeshandles.pixelmask,'XTick',[],'YTick',[],'Box','on');
    set(axeshandles.pixelmask,'UserData',b);
    EnablePixelMask(1); % not sure this is important
    set(controls.alphatext,'String',['p = ',num2str(alpha,2)]);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a list of all the subfields of the structure
% that is maintained in the UserData field of the plotting
% figure:
%
% a
%   stats
%       npixperstix
%       nstixperside
%       msperframe
%       gammaTable
%       invgammaTable
%       monSpd
%       gausslocut
%       gausshicut
%       gotHeader
%       gunnoise.mu
%       gunnoise.sigma
%       gunnoise.spike{i}.STS
%       gunnoise.spike{i}.STCross
%       gunnoise.spike{i}.nspikes
%   uicontrols
%       reset
%       nframes
%       eventqueuelength
%       STAslider
%       PCslider
%   axeshandles
%       STA
%       PC (UserData of axes contain eignevectors)
%       text
%       synthImage
%       pixelmask
%       subunit
%       
% gl
%       initparams.stepsize 
%       initparams.stepsizescale 
%       initparams.nreversals 
%       initparams.oogscale 
%       initparams.tolerance  
%       numtargetspikerate 
%       targetspikerates 
%       trialspecidx 
%       trialspecs 
%       timetoleave
%       wtsstopsendingflag
%       reversalflag
