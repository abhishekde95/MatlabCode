function p = WhiteNoiseOnline_AD()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Online spike triggered averaging and covariance
% for use with the WhiteNoise.d and WhiteNoise.m.
% GDLH 12/9/07
%
% Upgrade to work with multiple spike channels and
% binary cone noise.
% GDLH 11/3/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global subunitcmap
subunitcmap = [0.89412 0.10196 0.1098;0.21569 0.49412 0.72157;0.30196 0.68627 0.2902;0.59608 0.30588 0.63922;1 0.49804 0;1 1 0.2;0.65098 0.33725 0.15686;0.96863 0.50588 0.74902;0.6 0.6 0.6];

C = WhiteNoiseCodes;        % Script that defines a bunch of constants
DEFAULTNFRAMES = 10;  % Number of frames to go back before each spike
% Some initializations
p = InitPStruct(C.EOTCD, C.FIXXCD); % Structure of the events from each trial
s = InitPlex();     % Link to Plexon datastream
InitStatsStruct();  % Structure that holds the statistics
SetUpFig(DEFAULTNFRAMES);     % Setting up the figure.  I'm going to try to keep data in the 'UserData' field.
% Global Variables
global udpCom;  % Needed by "DealWithMessage()" and associated subfunctions
    udpCom.sock = [];
    udpCom.port = 6665;
    udpCom.rexip = '192.168.1.120';

[udpCom.sock, success] = pnetStart(udpCom.port);  % UDP communication with REX
if ~success, return; end
stopnow = 0;    % Flag to break out of the endless while loop (hit ESC to set to 1)

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
        a = get(gcf,'UserData');
        init_subunit_mask();
        %  Abhishek made a change - removed this if condition
%         if (~isfield(a.stats,'gunnoise') && ~isfield(a.stats,'conenoise'))  % ClearStats only the first time we get the header
        ClearStats; % Just fills up the structures in 'Userdata' with NaNs
%         end
        UpdatePixelMask(nan,nan);
    end
    % If we haven't received the header yet go back to the
    % begining and reset the 'p' structure.
    stats = getfield(get(gcf,'UserData'),'stats');
    if (~stats.gotHeader)
        p = InitPStruct(C.EOTCD, C.FIXXCD);
        continue;
    end
    % By this point, we know we've got the header and p.events contains at
    % least 1 EOTCD.
    [seed, nframes, mu, sigma, bkgndrgb, noisetype, gotfulltrial] = GetTrialParams(p); 
    if (~gotfulltrial)
        p = CleanUpEvents(p); % Remove a trial with EOT but without data (e.g. fixation breaks)
        continue;
    end
    % Do the calculations
    if (stats.gotHeader)
        EnableResetButton(0);
        plotnow = UpdateSTX(p, seed, nframes, mu, sigma, bkgndrgb, noisetype);
        if (plotnow)
            PlotSTA(noisetype);
            % Abhishek made change, introduced some new variables ****************%
            [n_subunits,~] = check_num_subunits();
            if (~ n_subunits)
                UpdatePixelMask(nan,nan);
            end
            drawnow;
        end
        EnableResetButton(1);
        p = CleanUpEvents(p);
    end
    
end % while (~stopnow)
PL_Close(s);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting the header from the ecodes
function GetHeader(p)
    C = WhiteNoiseCodes;

    gammaTable = codes2num(GetVal([C.GAMMATABLECD C.GAMMATABLECD], p.events));
    monSpd = codes2num(GetVal([C.MONSPDCD C.MONSPDCD], p.events));
    fundamentals = codes2num(GetVal([C.FUNDAMENTALSCD C.FUNDAMENTALSCD], p.events));

    a = get(gcf,'UserData');
    a.stats.npixperstix = GetVal(C.NPIXPERSTIXCD, p.events);
    a.stats.nstixperside = GetVal(C.NSTIXPERSIDECD, p.events);
    a.stats.msperframe = 1000/GetVal(C.FRAMERATECD, p.events, 'double');
    a.stats.gausslocut = GetVal(C.GAUSSLOCUTCD, p.events)/1000;
    a.stats.gausshicut = GetVal(C.GAUSSHICUTCD, p.events)/1000;
    a.stats.gammaTable = reshape(gammaTable, length(gammaTable)/3,3);
    a.stats.invgammaTable = InvertGamma(a.stats.gammaTable,1);
    a.stats.monSpd = reshape(monSpd, length(monSpd)/3,3);
    a.stats.fundamentals = reshape(fundamentals, length(fundamentals)/3,3);
    P_device = SplineSpd(linspace(380,780,size(a.stats.monSpd,1))', ...
        a.stats.monSpd, (380:5:780)');
    a.stats.M = a.stats.fundamentals'*P_device;
    a.stats.invM = inv(a.stats.M);
    a.stats.gotHeader = 1;
    set(gcf,'UserData',a);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting the trial parameters from the ecodes
function [seed, nframes, mu, sigma, bkgndrgb, noisetype, gotfulltrial] = GetTrialParams(p)
    C = WhiteNoiseCodes;
    try
        seed = GetVal(C.SEEDCD, p.events, 'long');
        nframes = GetVal(C.NFRAMESCD, p.events);
        mu = GetVal([C.MU1CD C.MU2CD C.MU3CD], p.events)/1000;
        sigma = GetVal([C.SIGMA1CD C.SIGMA2CD C.SIGMA3CD], p.events)/1000;
        noisetype = GetVal(C.NOISETYPECD, p.events);
        bkgndrgb = GetVal([C.BKGNDRCD C.BKGNDGCD C.BKGNDBCD], p.events);
        firstsynthimage = find(p.events == C.SYNTHIMAGECD, 1);
        if (~isempty(firstsynthimage) && (firstsynthimage < find(p.events == p.processnowcode, 1)));
            error('This is a synthetic image trial');
        end
        firststimon = find(p.events == C.STIMONCD, 1);
        if (firststimon > find(p.events == p.processnowcode, 1));
            firststimon = [];
        end
        if(isempty(firststimon));  % at least one STIMONCD
            error('no STIMONCD found');  %GDLH This is a problem!!
        end
        gotfulltrial = 1;
    catch
        seed = nan;
        nframes = nan;
        mu = nan;
        sigma = nan;
        bkgndrgb = nan;
        noisetype = nan;
        gotfulltrial = 0;
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
    b.imagestosend = {};
    set(a.axeshandles.synthImage,'UserData', b);
    if (a.stats.nstixperside)
        for i = 1:2  % preallocating space
            a.stats.gunnoise.spike{i}.STS = zeros([3*a.stats.nstixperside^2 nframesback]);
            a.stats.gunnoise.spike{i}.STCross = zeros([(3*a.stats.nstixperside^2)^2 nframesback]);
            a.stats.conenoise.spike{i}.STS = zeros([3*a.stats.nstixperside^2 nframesback]);
            a.stats.conenoise.spike{i}.STCross = zeros([(3*a.stats.nstixperside^2)^2 nframesback]);
        end
    else % don't know how much to preallocate yet
        for i = 1:2
            a.stats.gunnoise.spike{i}.STS = zeros(size(a.stats.gunnoise.spike{i}.STS));
            a.stats.gunnoise.spike{i}.STCross = zeros(size(a.stats.gunnoise.spike{i}.STCross));
            a.stats.conenoise.spike{i}.STS = zeros(size(a.stats.conenoise.spike{i}.STS));
            a.stats.conenoise.spike{i}.STCross = zeros(size(a.stats.conenoise.spike{i}.STCross));
            % a hack to avoid problems when "reset" is clicked during a header
            % drop: if nstixperside hasn't been dropped yet, we still want
            % to clear the contents of STS and STCross.
        end
    end
    % STCOVmex('init',{a.stats.nstixperside^2,3,nframesback});  %
    % Unnecessary?  GDLH 5/21/09
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
%       SetUpAxes
%       InitSynthImage
function SetUpFig(DEFAULTNFRAMES)
    global subunitcmap
    gcf = figure(1);
    set(gcf,'DefaultAxesUnits','pixels')
    set(gcf,'position',[100 200 900 250]);
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
    a.uicontrols.whichsigtest = uicontrol('style','popupmenu','Callback',@UpdatePixelMask,'String',{'Mean','Var'},'Position',[25 220 50 15]);
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
        a.axeshandles.text = axes('position',[figpos(3)-150 40 1 1]); % For nspikes
        set(gca,'Visible','off');
        a.axeshandles.synthImage = axes('position',[805 115 AXESWIDTH AXESWIDTH]); % For the SynthImage
        InitSynthImage(a.axeshandles.synthImage);
        a.axeshandles.pixelmask = axes('position',[25 130 AXESWIDTH AXESWIDTH]); % For the pixelmask
        Initpixelmask(a.axeshandles.pixelmask);
        % subunit axes
        a.axeshandles.subunit = axes('Position',[678 130-.25*AXESWIDTH 1.5*AXESWIDTH 1.5*AXESWIDTH]);
        InitSubunit(a.axeshandles.subunit) 
%         init_subunit_mask();
        
        function InitSubunit(axesh)
            maxsubunits = size(subunitcmap, 1);
            subunit_pos = get(axesh, 'Position');
            axes_w = subunit_pos(3);
            slider_w = axes_w; slider_h = 15;
            slsubunit_pos = [subunit_pos(1) subunit_pos(2)+axes_w+5 slider_w slider_h];
            updatemask_h = 22;
            updatemask_pos = [subunit_pos(1) subunit_pos(2)-updatemask_h-5 axes_w updatemask_h];
            text_w = 22; text_h = text_w;
            clearmask_pos = [updatemask_pos(1) updatemask_pos(2)-updatemask_h-50 updatemask_pos(3:4)];
            nudgestim_pos = [updatemask_pos(1) updatemask_pos(2)-updatemask_h-5 updatemask_pos(3:4)];
            
            set(axesh,'UserData',struct('mask',[]));
            a.uicontrols.txsubunit = uicontrol('Style','edit','ForegroundColor',subunitcmap(1,:),'FontSize',12,'FontWeight','bold','String','1','Enable','inactive','Position',[slider_w/2-text_w/2+slsubunit_pos(1) slsubunit_pos(2)+slider_h+5 text_w text_h]);
            a.uicontrols.slsubunit = uicontrol('Style','slider','Callback',{@update_subunit_text,a.uicontrols.txsubunit},'Min',1,'Max',9,'SliderStep',[1/(maxsubunits-1) 1/(maxsubunits-1)],'Value',1,'Position',slsubunit_pos);
            a.uicontrols.updatemask = uicontrol('style','togglebutton','string','Update Mask','Position',updatemask_pos,'Min',0,'Max',1,'Value',0);
            a.uicontrols.clearmask = uicontrol('style','pushbutton','string','Clear Mask','Position',clearmask_pos,'Callback',{@clear_subunit_mask,a.uicontrols.updatemask,axesh});
            a.uicontrols.nudgestim = uicontrol('style','pushbutton','string','Nudge stim','Position',nudgestim_pos,'Callback',@nudge_stim);

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
axis image;
b.STA = [];
b.PC = [];
b.gabor = [];
b.localimage = [];
b.imagestosend = {};
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
init_subunit_mask(); % clear the subunit mask
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
    elseif (any(ismember(handles.PC, axishandle)))
        b.PC = data;
        update_subunit_axes(data);    % Abhishek made change 08/2015
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
% Callback invoked when you click  the synthetic image.
% Prepare an image to be sent to REX.
function SynthImageCallback(~, ~, axeshandle)
    b = get(axeshandle,'UserData');
    if (isempty(b.imagestosend))
        b.imagestosend = {b.localimage};
    else
        b.imagestosend = [b.imagestosend {b.localimage}];
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
%
% If the type of image battery is selected to be "L vs M"
% we show a series of images with the same spatial weighting
% function as the STA+PC but with parametrically varied L and M
% cone contrasts (and zero S-cone contrast).  Stimuli are shown
% in 8 directions in the LM plane (0,45,90...)
%
% If the type of image battery is selected to be "FixCols"
% we show a 1-D or 2-D (depending on whether an STA and/or
% PC1 image is selected) image set.  This type of image bettery
% differs from 'standard' in that we replace the STA with the
% the first SVD component, and we force the color direction of the PC
% to be "luminance" which I'm defining here as pure radiance modulation.

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

    % Calling the appropriate nested subfunction
    if (strcmp(imtypestrings(imtypeval),'Standard'))
       images = StandardImageBattery;
       figure(2);
       set(gcf,'position',[1150 500 100*numel(PCcontrasts) 100*numel(STAcontrasts)]);
       for ii = 1:numel(PCcontrasts)* numel(STAcontrasts)

           subplot(numel(STAcontrasts),numel(PCcontrasts),ii),image(images{1,ii});

           set(gca,'XTick',[],'YTick',[],'Box','on');
       end
%        close(gcf);
       figure(1);
    elseif (strcmp(imtypestrings(imtypeval),'L vs M'))
       images = LvsMImageBattery;
    elseif (strcmp(imtypestrings(imtypeval),'FixCols'))
       images = FixColsImageBattery;
    end
    
    % Abhishek made a change - 09/2015
    [n_subunits,subunit_mask] = check_num_subunits() ; % 0 corresponds to no subunits
    if (n_subunits)
        images = compress_vector(images, n_subunits, subunit_mask);
    end
    
    b.imagestosend = [b.imagestosend images(randperm(length(images)))];
    set(a.axeshandles.synthImage,'UserData',b);
 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nested function for doing a standard image battery
    function images = StandardImageBattery
        images = cell(0);
        for STAvalue = STAcontrasts
            for PCvalue = PCcontrasts
                im = STAvalue*(STAimage-muimage)+PCvalue.*(PCimage-muimage);
                im = im/(eps + max(STAcontrasts) + max(PCcontrasts));
                im = im + muimage;
                images = [images {im}];
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nested function for doing an L vs M image battery
    function images = LvsMImageBattery

        im = STAlim*(STAimage-muimage)+PClim.*(PCimage-muimage);
        [u,~,~] = svd(reshape(im,[nstixperside^2 3]));
        spatialweightingfunction = u(:,1);
        spatialweightingfunction = spatialweightingfunction./max(abs(spatialweightingfunction));
        thetas = linspace(0,2*pi-pi/4,8)';  % Having this '8' hardcoded is ugly and inflexible

        contrasts = linspace(0,STAlim/10,ncontrastlevels);
        conecontrasts = [];
        for i = 1:ncontrastlevels
            if (all(contrasts(i) == 0))
                conecontrasts = [conecontrasts; 0 0 0];
            else
                conecontrasts = [conecontrasts; contrasts(i)*[cos(thetas) sin(thetas) zeros(length(thetas),1)]];
            end
        end
        bkgndlms = a.stats.M*mu;
        bkgndlmsmat = repmat(bkgndlms',size(conecontrasts,1),1);
        coneexcitations = (1+conecontrasts).*bkgndlmsmat;
        rgbs = (a.stats.invM*coneexcitations')';
        rgbs = rgbs-repmat(mu',size(rgbs,1),1);

        images = cell(0);
        for i = 1:size(rgbs,1)
            im = spatialweightingfunction*rgbs(i,:);
            im = reshape(im,[nstixperside nstixperside 3])+muimage;
            images = [images {im}];
        end

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nested function for doing a fixed color image battery
    % GDLH: This is outdated - can get rid of this function
    function images = FixColsImageBattery

        [u,s,v] = svd(reshape(STAimage-muimage,[nstixperside^2 3]));
        STAim = u(:,1)*s(1)*v(:,1)';
        STAim = reshape(STAim,[nstixperside nstixperside 3])+muimage;


        [u,~,~] = svd(reshape(PCimage-muimage,[nstixperside^2 3]));
        spatialweightingfunction = u(:,1)./max(abs(u(:,1)));
        PCim = spatialweightingfunction*([.2 .2 .2].*mu');
        PCim = reshape(PCim,[nstixperside nstixperside 3])+muimage;

        images = cell(0);
        for STAvalue = STAcontrasts
            for PCvalue = PCcontrasts
                im = STAvalue*(STAim-muimage)+PCvalue.*(PCim-muimage)+muimage;
                images = [images {im}];
            end
        end
    end
end

function img = compress_vector(images,num_subunits,subunit_mask)
   subunit_mask(subunit_mask == 0) = num_subunits + 1;
   img = cell(0);
%    new_mask = [subunit_mask(:); (num_subunits+1)+subunit_mask(:); 2*(num_subunits+1)+subunit_mask(:)];
%    [~,a,~] = unique(new_mask);
%    for ii = 1:size(images,2)
%        tmp_img = images{1,ii};
%        vec = tmp_img(:);
%        vec_tmp = vec(a);
%        vec_tmp = reshape(vec_tmp,[num_subunits+1 1 3]);
%        img = [img {vec_tmp}];
%    end
% keyboard
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
    a = get(gcf,'UserData');
    b = get(a.axeshandles.synthImage,'UserData');
    b.imagestosend = [];
    set(a.axeshandles.synthImage,'UserData',b);
end

% Destructively modifies the 'pixelmask' field in the axes
% UserData to mark the user's mouse clicks.  Doesn't do any
% plotting though - that's taken care of in UpdatePixelMask.
function PixelMaskCallback(~,~)
    a = get(gcf,'UserData');
    %stats = a.stats;
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

% Fit a Gabor function to the synthimage
% Fitting might take some time and it might
% happen while some other function is rewriting
% the 'UserData' field.  Careful.  Disable this
% button during plotting.
function FitGaborCallback(h,~)
    a = get(gcf,'UserData');
    noisetype = get(a.uicontrols.whichnoisetype,'Value');
    b = get(a.axeshandles.synthImage,'UserData');
    contrast = get(a.uicontrols.STAslider,'Value');
    im = b.localimage;
    npix = size(im,1);
    im = reshape(im,[npix*npix,3]);
    im = im-repmat(a.stats.murgb',[npix*npix,1]);
    if (noisetype == 2)
        im = (a.stats.M*im')./repmat(a.stats.conenoise.sigma,1,npix*npix);
        im(isinf(im)) = 0;  % In case we have a divide by zero.
        im = im';
    end
    [u,~,v] = svd(im');
    if (max(abs(v(:,1))) ~= max(v(:,1)))
        v = -v;
        u = -u;
    end
    disp('Color: ');
    u(:,1)
    if (noisetype == 2) % Converting color direction back to rgb.  Is this correct?
        u(:,1) = a.stats.invM*u(:,1);
        u(:,1) = u(:,1)./norm(u(:,1));
    end
    im = reshape(v(:,1),[npix npix]);
    out = gaborfit(im);
    b = get(a.axeshandles.synthImage,'UserData'); % In case it's changed

    if (out.exitflag)
        b.gabor.rgb = u(:,1)';
        b.gabor.theta = out.theta;
        b.gabor.lambda = out.lambda;
        b.gabor.phi = out.phi;
        b.gabor.sigma = out.sigma;
        b.gabor.gamma = out.gamma;
        b.gabor.xoffset = out.xoffset;
        b.gabor.yoffset = out.yoffset;
        b.gabor.exitflag = out.exitflag;
        PlotGabor(a.axeshandles.synthImage, b.gabor, npix, a.stats.murgb, contrast);
    else
        disp('Gabor fit did not converge');
    end
    set(a.axeshandles.synthImage,'UserData', b);

    % Imbedded function for plotting the gabor
    function PlotGabor(axesh, struct, nsubpix, bkgndrgb, contrast)
        axes(axesh);
        pos = get(axesh,'Position');
        npix = pos(4);
        pixperstix = npix/nsubpix;
        margin = 1-(1/pixperstix);
        interval = linspace(1-margin,nsubpix+margin,npix)-ceil(median(1:nsubpix));
        [X, Y] = meshgrid(interval,interval);
        X = X-struct.xoffset;
        Y = Y+struct.yoffset;
        xprime = X.*cos(-struct.theta)+Y.*sin(-struct.theta);
        yprime = -X.*sin(-struct.theta)+Y.*cos(-struct.theta);
        fittedgabor = exp(-(xprime.^2+struct.gamma.^2.*yprime.^2)./(2.*struct.sigma.^2)).*cos(2.*pi.*yprime./struct.lambda-struct.phi);

        colorgabor = struct.rgb'*reshape(fittedgabor,1,npix^2)/2.5*contrast;
        colorgabor = colorgabor+repmat(bkgndrgb,1,npix^2);
        colorgabor = reshape(colorgabor,[3, npix, npix]);
        colorgabor = permute(colorgabor,[2 3 1]);
        h = image(colorgabor);
        set(h,'ButtonDownFcn',{@UpdateSynthImage});
        set(axesh, 'XTick',[],'YTick',[]);
    end

end

function update_subunit_axes(STA)
a = get(gcf, 'UserData');
ax = a.axeshandles.subunit;
UD = get(ax, 'UserData');
if ~isempty(STA)
    image(STA, 'Parent', ax, 'ButtonDownFcn', @update_subunit_mask);
end
set(ax, 'UserData', UD);
set(ax, 'XTick', [], 'YTick', [], 'Box', 'on');
axis(ax, 'image');
grid(ax, 'off');
paint_subunit_selection(ax);
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

function init_subunit_mask()
a = get(gcf, 'UserData');
UD = get(a.axeshandles.subunit, 'UserData');
UD.mask_rex = nan(a.stats.nstixperside);
UD.mask_display = nan(a.stats.nstixperside);
UD.nudge = [];
set(a.axeshandles.subunit, 'UserData', UD);
end

function update_subunit_text(h, ~, htext)
global subunitcmap
n = round(get(h, 'Value'));
set(h, 'Value', n);
set(htext, 'String', n, 'ForegroundColor', subunitcmap(n,:));
end

function clear_subunit_mask(~, ~, hupdatemask, haxsub)
set(hupdatemask, 'Value', 0);
init_subunit_mask();
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

function nudge_stim(~,~)
a = get(gcf, 'UserData');
UD = get(a.axeshandles.subunit, 'UserData');
x = -(a.stats.nstixperside-1)/2:(a.stats.nstixperside-1)/2;
normmask = (UD.mask_display >= 1)./sum(UD.mask_display(:) >= 1); % normalized to sum to 1
UD.nudge = [sum(normmask)*x' fliplr(x)*sum(normmask,2)]; % +x is right, +y is up. In stixels.
set(a.axeshandles.subunit, 'UserData',UD);
end

function update_subunit_mask(himage, ~) % buttondown callback for subunit image
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
    C = WhiteNoiseCodes;
    a = get(gcf,'UserData');
    stats = a.stats;
    controls = a.uicontrols;
    nframesback = str2double(get(controls.nframes,'String'));
    gammaTable = stats.gammaTable;
    invgamma = stats.invgammaTable;
    nstixperside = stats.nstixperside;
    controls = a.uicontrols;
    % Abhishek made change, introduced some new variables ****************%
    
    [n_subunits,~] = check_num_subunits() ; % 0 corresponds to no subunits
    if (n_subunits)
        nrandnums_perchannel = n_subunits;
        if get(controls.updatemask,'Value')
           disp('Update button Pressed');
            for i = 1:2
                stats.gunnoise.spike{i}.STS = zeros([3*nrandnums_perchannel nframesback]);
                stats.gunnoise.spike{i}.STCross = zeros([(3*nrandnums_perchannel)^2 nframesback]);
                stats.gunnoise.spike{i}.nspikes = 0;
                stats.conenoise.spike{i}.STS = zeros([3*nrandnums_perchannel nframesback]);
                stats.conenoise.spike{i}.STCross = zeros([(3*nrandnums_perchannel)^2 nframesback]);
                stats.conenoise.spike{i}.nspikes = 0;
            end
        end
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
        %rgbmat = [];
        %for i = 1:size(colordirlms,1)
        %    lms = (M*stats.murgb)+(colordirlms(i,:)'.*sigma);
        %    rgbmat = [rgbmat; (invM*lms)'];  % these are rgb intensities
        %end  % GDLH not sure if this is necessary
        %stats.lmsbinaryrgbmat = rgbmat+repmat(stats.murgb',size(rgbmat,1),1); %GDLH not sure if this is necessary
        randnums = getEJrandnums(nrandnums_perchannel*nframes, seed);
        randnums = mod(randnums, 8)+1;
        randnums = colordirlms(randnums,:);
        randnums = reshape(randnums, [nrandnums_perchannel, nframes, 3]);
        randnums = permute(randnums,[1 3 2]);
        randnums = reshape(randnums, [nrandnums_perchannel*3, nframes]);
        % Each column should be Ls followed by Ms followed by Ss for each frame.
    end

    t_stimon = p.times(find(p.events == C.STIMONCD, 1));
    nspikesthistrial = [0 0 0];
    for i = 1:2
        spiketimes = p.spikes{i}-t_stimon;
        if (any(spiketimes))
            frametimes = linspace(0, nframes*stats.msperframe, nframes)+(stats.msperframe/2)';
            spiketimes(spiketimes < nframesback*stats.msperframe) = [];
            spiketimes(spiketimes > frametimes(end)) = [];
            n = hist(spiketimes, frametimes);
            % Doing the accumlation here rather than in
            % STCOVmex  GDLH 11/3/08
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
        if ((L1 | L2) & nspikesthistrial(whichspike) > 1)
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
    function imagerequest() %#ok<DEFNU>
        a = get(gcf,'UserData');
        axeshandles = a.axeshandles;
        b = get(axeshandles.synthImage,'UserData');
        if (~isempty(b.imagestosend))
%             [n_subunits,~] = check_num_subunits();
            im = b.imagestosend{1};
            b.imagestosend(1) = [];
            set(axeshandles.synthImage,'UserData',b);
            nstixperside = a.stats.nstixperside;
            invgammaTable = a.stats.invgammaTable;
            for gun = 1:3
                plane = invgammaTable(round(im(:,:,gun)*65535)+1,gun); % goes from 0:1 intensity to 0:1 voltage 
                im(:,:,gun) = reshape(plane,[nstixperside nstixperside]);
            end
            im = round(im*65535)+1; % goes from 0:1 voltage to 1:65536 voltage
            sendToRex(udpCom, im, 'integer', 'imagerequest();');
        end
    end

% Send the parameters for the best-fitting gabor to REX on request.
    function gaborrequest() %#ok<DEFNU>
        a = get(gcf,'UserData');
        axeshandles = a.axeshandles;
        b = get(axeshandles.synthImage,'UserData');
        
        if (isfield(b.gabor,'theta'))
            if (~isempty(b.gabor.theta))
                paramvector = [b.gabor.rgb...
                    b.gabor.theta b.gabor.lambda b.gabor.phi b.gabor.sigma...
                    b.gabor.gamma b.gabor.xoffset b.gabor.yoffset];
                b.gabor.theta = [];
                set(axeshandles.synthImage,'UserData',b);
            else
                paramvector = [];
            end
            sendToRex(udpCom, paramvector, 'double', 'gaborrequest();');
        end
    end

% Send the stixel mask to REX on request.
    function mask = maskrequest() %#ok<DEFNU>
        UD = get(gcf, 'UserData');
        if get(UD.uicontrols.updatemask, 'Value')
            UDsubunit = get(UD.axeshandles.subunit, 'UserData');
            mask = UDsubunit.mask_display;
            set(UD.uicontrols.updatemask, 'Value', 0);
            UDsubunit.mask_rex = mask;
            set(UD.axeshandles.subunit, 'UserData', UDsubunit);
            
            % Abhishek made a change to this function
            nframesback = str2double(get(UD.uicontrols.nframes,'String'));
            [nrandnums_perchannel,~] =  check_num_subunits(); % 0 corresponds to no subunits
            disp('Update on');
            for i = 1:2
                UD.stats.gunnoise.spike{i}.STS = zeros([3*nrandnums_perchannel nframesback]);
                UD.stats.gunnoise.spike{i}.STCross = zeros([(3*nrandnums_perchannel)^2 nframesback]);
                UD.stats.gunnoise.spike{i}.nspikes = 0;
                UD.stats.conenoise.spike{i}.STS = zeros([3*nrandnums_perchannel nframesback]);
                UD.stats.conenoise.spike{i}.STCross = zeros([(3*nrandnums_perchannel)^2 nframesback]);
                UD.stats.conenoise.spike{i}.nspikes = 0;
            end
            
            InitSynthImage(UD.axeshandles.synthImage)
            set(gcf,'Userdata',UD);
            %***********************************************************
            
        else
            mask = getfield(get(UD.axeshandles.subunit, 'UserData'), 'mask_rex');
        end
        mask(isnan(mask)) = 0;
        sendToRex(udpCom, mask(:), 'integer', 'maskrequest();');
    end

% Send a displacement of the stimulus [(x,y) in stixels] to REX on request.
    function nudgerequest() %#ok<DEFNU>
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create and display a linear combination of STA and PC.
% Invoked when we adjust the STA or PC slider, or when
% we select or deselect at STA or PC frame.  Since these
% are all Callbacks, they can't co-occur accidentally.
function UpdateSynthImage(~,~)
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
        im = STAvalue*(STAimage-muimage) + PCvalue.*(PCimage-muimage);
        im = im/(abs(STAvalue)+abs(PCvalue));
        im = im + muimage;       
    else
        im = muimage;
    end
    
    
%     disp(max(STAimage(:)));
%     disp(min(STAimage(:)));
%     disp(max(PCimage(:)));
%     disp(min(PCimage(:)));
    
    b.localimage = im;
    axes(axeshandles.synthImage);
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
    
    % Abhishek made change, introduced some new variables ****************%
    [n_subunits,~] = check_num_subunits(); % 0 corresponds to no subunits
    if (n_subunits)
        nrandnums_perchannel = n_subunits;
    else
        nrandnums_perchannel = stats.nstixperside^2;
    end
    %*******************************************************************%
    
    mu = stats.murgb;
    if (all(mu == 0))
        disp('No background rgb set');
        keyboard
    end
%     sigma = stats.gunnoise.sigma;

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
    text(0,0,['nspikes(',noisetypestr,' ',num2str(whichspike),') = ',num2str(n)]);

    % Iterating to get PCs
    PCs = zeros(size(STAs));
    mask = logical(repmat(pixelmask(:),3,1));
    % Abhishek made change - introduced 3 new variables from compute_basis_vec
    an_mask = [];
    an_mask = zeros(nrandnums_perchannel, 1);   % Pixel mask. Someday make a tool to make this non-zero.
    Lmask = logical(repmat(~an_mask(:),[3 1]));
    SB = get(a.axeshandles.subunit, 'Userdata');
    Subunit_mask = SB.mask_rex; % declaring Subunit_mask here 
    Subunit_mask(isnan(Subunit_mask)) = 0;
    %*********************************************************************
    
    % Just bail out if no pixels are selected - to avoid a crash
    if (sum(mask) == 0)
        disp('No pixels selected!  Not plotting anything.');
        return;
    end

    % Don't look at the smallest PC until we have enough spikes.
    if (n < 3*sum(mask))
        set(controls.smallpc,'Value',0);
    end
    SMALLPC = get(controls.smallpc,'Value');

    for ik = 1:nframesback
        tmp = STS(:,ik)*STS(:,ik)';
        STCframe = (n.*STCross(:,ik)-tmp(:))/(n*(n-1));
        STCframe = reshape(STCframe,[sqrt(size(STCframe,1)) sqrt(size(STCframe,1))]);
        STA = STAs(:,ik);

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
    
    % Abhishek made change here - 08/2015
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
        % Abhishek made this change - 08/2015
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

    

    % Converting from a matrix of l,m,s differences to rgbs intensities
    % for plotting
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
