function Grating
% Grating.m
%
%    Slave program for displaying drifing gratings.
% GDLH 3/1/08
%
%    Doing all the heavy computations in ShowGrating instead of DrawGrating
% GDLH 4/18/16

% Communication variables
    global udpCom;
    udpCom.port = 6665;
    udpCom.sock = nan;
    udpCom.rexip = '192.168.1.120';

    % Keypress constants
    global KEY
    KbName('UnifyKeyNames');
    KEY.ESC = KbName('escape');
    
    % Variables that are destructively modified by subfunctions
    global gl;
    gl.mondistcm = 0;
    gl.screenWidthcm = 0;
    gl.screenWidthpix = 0;
    gl.screenHeightpix = 0;
    gl.screenCenterXpix = 0;
    gl.screenCenterYpix = 0;
    gl.pixperdeg = 0;
    gl.bkgndRGB = [0 0 0];
    gl.bkgndrgb = [0 0 0];
    gl.windowPtr = 0;
    gl.fliprequest = 0;
    gl.framerate = 0;
    gl.framecounter = 0;
    gl.framecountermax = 0;
    gl.timetoleave = 0;
    gl.vpixx = IsVPixx();
    gl.ccmode = 1;
    gl.tex = []; % precomputing textures GDLH 4/18/16
    
    gl.ep.x = 0;
    gl.ep.y = 0;
    
    gl.fp.on = 0;
    gl.fp.x = 0;
    gl.fp.y = 0;
    gl.fp.size = 0;
    gl.fp.rgb = [0 0 0];
    
    gl.grating.on = 0;
    gl.grating.drawrect = [0 0 0 0];
    gl.grating.diam = 1;
    gl.grating.sf = 1;
    gl.grating.tf = 1;
    gl.grating.conecontrast = [0 0 0];
    gl.grating.rgb = [0 0 0];
    gl.grating.phase = 0;
    gl.grating.orient = 0;
    
    gl.screenflash.on = 0;
    gl.screenflash.rgb = [0 0 0];
    
    gl.cal.gammaTable = [];
    gl.cal.invgammaTable = [];
    gl.cal.monSpd = [];
    gl.cal.fundamentals = [];
    gl.cal.M = zeros(3);
    gl.cal.invM = zeros(3);

    [udpCom.sock, Success] = pnetStart(udpCom.port);
    if ~(Success)
        return
    end
    pnet(udpCom.sock, 'setreadtimeout', 0);
    pnet(udpCom.sock, 'setwritetimeout', 0);
    disp('In Grating');

    messageIsAvailable = 0;
    while ~messageIsAvailable  % endless loop.
        messageIsAvailable = pnet(udpCom.sock, 'readpacket', 1000, 'no block');
        if (messageIsAvailable)
            DealWithMessage(messageIsAvailable);
            messageIsAvailable = 0;   
        end     
        if (gl.windowPtr > 0)
            if (gl.grating.on)
                DrawGrating();
            end
            if (gl.fp.on)
                DrawFP();
            end
            if (gl.screenflash.on)
                DrawScreenFlash();
            end
            if (gl.fliprequest)
                DoFlip();
            else
                gl.timetoleave = check_for_quit;
            end
        end
        if (gl.timetoleave)
            close_textures();
            return;
        end
        [keyisdown,secs,keycode] = KbCheck;
        if (keyisdown && keycode(KEY.ESC))
           ShowCursor;
           close_textures();
           pnet(udpCom.sock, 'close');
           Screen('CloseAll');
           return;
        end
    end
end
%%
% check_for quit
function yn = check_for_quit()
    yn = false;
    [keyisdown,nil,keycode] = KbCheck(-1); %#ok<ASGLU>
    if keyisdown
        yn = keycode(KbName('ESCAPE'));
    end
end

%%
% close_textures()
% Now that we're calculating the grating textures in ShowGrating
% we have to remember to clear them out at the end of each trial and
% when we exit this program.
function close_textures()
    global gl
    if isfield(gl,'tex')
        Screen('Close', gl.tex);
    end
end
%%
% DoFlip() taken from WhiteNoise
% Should these be a single function somewhere?
function DoFlip()
    global gl;
    global udpCom;

    if (gl.framecounter == gl.framecountermax) & gl.screenflash.on
        Screen('FillRect', gl.windowPtr, gl.bkgndRGB/255); % The wrong color
    end%  Need to do this before the flip
    Screen('Flip',gl.windowPtr);
    if (gl.grating.on || gl.screenflash.on)
        if(gl.framecounter == 1)
            pnet(udpCom.sock, 'write', 'MACSTIMON>> >>');
            pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
        elseif (gl.framecounter == gl.framecountermax)
            pnet(udpCom.sock, 'write', 'MACSTIMOFF>> >>');
            pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
            close_textures();
            gl.grating.on = 0;
            gl.screenflash.on = 0;
        end
        gl.framecounter = gl.framecounter + 1;
    end
    gl.fliprequest = 0;
end

%%
% InitDisplay below was taken verbatim from WhiteNoiseGH.
% Should these be a single function somewhere?

function InitDisplay(mondist, screenwidth, calfilename, varargin)
    global gl;
    
    load(calfilename);
    cal = cals{end};
    gl.mondistcm = mondist;
    gl.screenWidthcm = screenwidth;
    gl.bkgndRGB = round(255*cal.bgColor)';
    gl.cal.gammaTable = cal.gammaTable;
    gl.cal.invgammaTable = InvertGammaTable(cal.gammaInput, gl.cal.gammaTable, 2^16);
    gl.cal.monSpd = cal.P_device;
    gl.stim.on = 0;
    gl.fp.on = 0;

    if nargin > 3 && ~isempty(varargin{1})
        fundfilename = varargin{1};
        s = load(fundfilename);
        fns = fieldnames(s);
        P_device = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls([380 5 81]));
        gl.cal.fundamentals = eval(['s.',fns{1}])';
        gl.cal.M = gl.cal.fundamentals'*P_device;
        gl.cal.invM = inv(gl.cal.M);
    end
    
    if nargin > 4 && ~isempty(varargin{2})
        gl.ccmode = varargin{2};
    end
    
    gl.bkgndrgb = [gl.cal.gammaTable(gl.bkgndRGB(1)+1,1),...
                   gl.cal.gammaTable(gl.bkgndRGB(2)+1,2),...
                   gl.cal.gammaTable(gl.bkgndRGB(3)+1,3)]';

    % Gamma correction is done in software so just set hardware gamma
    % lookup to the unity line.
    clut = repmat(linspace(0,1,256),3,1)';
    Screen('LoadNormalizedGammaTable', gl.windowPtr, clut);
           
    % start up the imaging pipeline
    if ~isempty(Screen('Windows'))
        close_textures();
        gl.windowPtr = max(Screen('Windows'));
        if gl.ccmode
            Screen('FillRect', gl.windowPtr, cal.bgColor);
        else
            Screen('FillRect', gl.windowPtr, gl.bkgndRGB);
        end 
    else
        PsychImaging('PrepareConfiguration');
        if gl.ccmode > 0 && gl.vpixx
            PsychImaging('AddTask', 'General', 'EnableDataPixxC48Output', gl.ccmode); % in mode '1': every 2nd column of pixels is ignored
            gl.windowPtr = PsychImaging('OpenWindow', 0, cal.bgColor);
        elseif gl.ccmode > 0
            PsychImaging('AddTask', 'General', 'EnableBits++Color++Output', gl.ccmode);
            gl.windowPtr = PsychImaging('OpenWindow', 0, cal.bgColor);
        else
            gl.windowPtr = PsychImaging('OpenWindow', 0, gl.bkgndRGB);
        end
    end

    gl.framerate = 1/Screen('GetFlipInterval', gl.windowPtr);
    [screenwidthpix, screenheightpix]  = Screen('WindowSize', gl.windowPtr);
    gl.screenWidthpix = screenwidthpix; % using Bits++ in Colour mode each pixel has a 1x2 aspect ratio
    gl.screenHeightpix = screenheightpix;  % but Psychophysicstoolbox doesn't (need to) know about this
    gl.screenCenterXpix = screenwidthpix/2;
    gl.screenCenterYpix = screenheightpix/2;
    
    pixpercm = gl.screenWidthpix/gl.screenWidthcm;
    theta = atan2(gl.screenWidthcm/2, gl.mondistcm)*180/pi;
    cmperdeg = gl.screenWidthcm/(2*theta);
    gl.pixperdeg = pixpercm*cmperdeg;
    
    HideCursor;
end

%%
% "Show" functions are called from REX.  They set up the appropriate fields 
% in the "gl" structure (and set the "...on" toggle field to 1). 
function ShowFP(x, y, size, fpr, fpg, fpb)
    global gl;

    gl.fp.x = x/10;
    gl.fp.y = y/10;
    gl.fp.size = size/10;
    gl.fp.rgb = [fpr, fpg, fpb]/255;
    
    fpsizeinpix = round(gl.pixperdeg*gl.fp.size);

    gl.fp.drawrect = [gl.screenCenterXpix+(gl.fp.x*gl.pixperdeg)-floor(fpsizeinpix/2)+1 ...
                gl.screenCenterYpix-(gl.fp.y*gl.pixperdeg)-floor(fpsizeinpix/2)+1 ...
                gl.screenCenterXpix+(gl.fp.x*gl.pixperdeg)+ceil(fpsizeinpix/2)...
                gl.screenCenterYpix-(gl.fp.y*gl.pixperdeg)+ceil(fpsizeinpix/2)];
    
    gl.fp.on = 1;
    
end

%%
% "Draw" functions get called on every screen refresh so long as the
% corresponding graphical object is to be displayed.  These are the
% functions that actually draw the object to the screen (and they make a
% fliprequest so that the drawn object will actually appear).
function DrawFP()
global gl;

Screen('Fillrect', gl.windowPtr, gl.fp.rgb, gl.fp.drawrect);

gl.fliprequest = 1;
end

%%
function HideFP()
global gl;

gl.fp.on = 0;
gl.fliprequest = 1;

end

%%
% Computes the textures that will be displayed in DrawGrating();
% Previously much of this computation was done in DrawGrating 
% (i.e. during the interframe interval) which precluded the high
% framerates afforded by the ProPixx.
% Updated 4/18/16 by GDLH.
function ShowGrating(rfx, rfy, diam, sf, tf,lcont, mcont, scont,  orient, phase, ncycles)
global gl;

gl.grating.x = rfx/10;
gl.grating.y = rfy/10;
gl.grating.diam = diam;
gl.grating.sf = sf;
gl.grating.tf = tf;
gl.grating.conecontrast = [lcont mcont scont];
gl.grating.orient = orient;
gl.grating.phase = phase;
if (tf == 0)  % If TF is zero, interpret ncycles as number of milliseconds.
    %(formerly frames, but I don't think this has ever been used). GDLH 4/18/16
    gl.framecountermax = floor(ncycles/1000*gl.framerate);
else
    gl.framecountermax = floor(ncycles/tf*gl.framerate);
end

if (gl.grating.diam == 0)  % Why would we ever need this?
    return;
end

x = (gl.grating.x+gl.ep.x)*gl.pixperdeg;
y = (gl.grating.y+gl.ep.y)*gl.pixperdeg;

% Calculating the drawing rectangle
if gl.ccmode
    stimsizeinpix = round(gl.pixperdeg*gl.grating.diam/2);  % /2 counting in double-width pixels
    gl.grating.drawrect = round([gl.screenCenterXpix+x-stimsizeinpix ...
        gl.screenCenterYpix-y-stimsizeinpix ...
        gl.screenCenterXpix+x+stimsizeinpix ...
        gl.screenCenterYpix-y+stimsizeinpix]);
    if(rem(gl.grating.drawrect(1), 2)) %if the rectangle starts on an odd pixel
        gl.grating.drawrect(1) = gl.grating.drawrect(1) - 1;
        gl.grating.drawrect(3) = gl.grating.drawrect(3) - 1;
    end
else
    stimsizeinDWpix = round(gl.pixperdeg*gl.grating.diam);  % counting in single-width pixels
    gl.grating.drawrect = round([gl.screenCenterXpix+x-stimsizeinpix ...
        gl.screenCenterYpix-y-stimsizeinpix ...
        gl.screenCenterXpix+x+stimsizeinpix ...
        gl.screenCenterYpix-y+stimsizeinpix]);
end

gl.ep.x = 0; % resetting
gl.ep.y = 0;

% Calculating grating rgbs
bkgndlms = gl.cal.M*gl.bkgndrgb;
gratinglms = bkgndlms.*[1+gl.grating.conecontrast'];
gl.grating.rgb = gl.cal.invM*gratinglms;
if (any(gl.grating.rgb > 1) | any(gl.grating.rgb < 0))
    gl.grating.rgb = gl.bkgndrgb;
end

% Creating an aperture template
if gl.ccmode
    [x,y] = meshgrid(linspace(-1,1,stimsizeinpix), linspace(-1,1,2*stimsizeinpix));
else
    [x,y] = meshgrid(linspace(-1,1,stimsizeinpix), linspace(-1,1,stimsizeinpix));
end
% Only half the number of columns as rows since we're using colour mode
aperture = (x.^2+y.^2) <= 1+2/stimsizeinpix;  % correction for binning

% copied from DrawStim
pixpercycle = gl.pixperdeg/gl.grating.sf;

xinc = (2*pi/pixpercycle)*cos((pi/2)-gl.grating.orient);
yinc = (2*pi/pixpercycle)*sin((pi/2)-gl.grating.orient);

if gl.ccmode % "*2" so that we can count in single-wide pixels 
    [xramp, yramp] = meshgrid(xinc*([0:2:2*stimsizeinpix-1]), yinc*([0:2*stimsizeinpix-1]));
else
    [xramp, yramp] = meshgrid(xinc*([0:1:stimsizeinpix-1]), yinc*([0:stimsizeinpix-1]));
end

phaseinc = gl.grating.tf*2*pi/gl.framerate;
for nth_frame = 1:gl.framecountermax
    phase = mod(gl.grating.phase+(nth_frame*phaseinc), 2*pi); % Debug
    a = cos(xramp+yramp+phase);
    a = a.*aperture;
    if gl.ccmode
        im = zeros(stimsizeinpix*2, stimsizeinpix, 3);
    else
        im = zeros(stimsizeinpix, stimsizeinpix, 3);
    end
    for plane = 1:3
        tmp = a.*(gl.grating.rgb(plane)-gl.bkgndrgb(plane))+gl.bkgndrgb(plane);
        tmp = round(tmp*size(gl.cal.invgammaTable,1)-1)+1;
        tmp = gl.cal.invgammaTable(tmp, plane);
        %tmp = round(tmp*(size(gl.cal.invgammaTable,1)-1));
        if gl.ccmode
            im(:,:,plane) = reshape(tmp, 2*stimsizeinpix, stimsizeinpix);
        else
            im(:,:,plane) = reshape(tmp, stimsizeinpix, stimsizeinpix);
        end
    end
    gl.tex(nth_frame)=Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(im,[],gl.ccmode),[],[],2);
end

gl.framecounter = 1;
gl.grating.on = 1;
end
%%
function DrawGrating()
    global gl;

    Screen('DrawTexture',gl.windowPtr,gl.tex(gl.framecounter), [], gl.grating.drawrect,[],0);
    
    gl.fliprequest = 1;
end
%%
% Changing the color of the background for nframes frames
function ShowFlash(lcont, mcont, scont, nframes)
global gl;
   
gl.screenflash.on = 1;
bkgndlms = gl.cal.M*gl.bkgndrgb;
lms = bkgndlms.*(1+[lcont mcont scont]');
gl.screenflash.rgb = (gl.cal.invM*lms);

if (any(gl.screenflash.rgb > 1) | any(gl.screenflash.rgb < 0))
    gl.screenflash.rgb = gl.bkgndrgb;
end
gl.framecounter = 1;
gl.framecountermax = nframes;
end

%%
function DrawScreenFlash()
global gl

    Screen('FillRect',gl.windowPtr, gl.screenflash.rgb, [0 0 gl.screenWidthpix gl.screenHeightpix]);
    gl.fliprequest = 1;
end

%%  
function AllOff()
    global gl;
    close_textures(); 
    gl.grating.on = 0;
    gl.fp.on = 0;
    gl.screenflash.on = 0;
    gl.fliprequest = 1;
end

%%
function eyepos(x, y)
    global gl;

    gl.ep.x = x/40;
    gl.ep.y = y/40;
end

%%
function DealWithMessage(msgSize)
    global udpCom;
    global gl;  % needs to be here for functions evaled by this one
    message = pnet(udpCom.sock, 'read', msgSize, 'char');
    if (strncmp(message,'return',6))
        a = dbstack;  % Check whether called from another function or from command line 
        if (~strcmp(a(end).name, mfilename))
            gl.timetoleave = 1;
        end
    end
    try 
        eval(message)
    catch
        fprintf('Ignoring uninterpretable message: "%s"\n',message);
        sca;
        error = lasterror;
        disp(error.message);
        disp(error.identifier);
        disp(error.stack);
    end
end