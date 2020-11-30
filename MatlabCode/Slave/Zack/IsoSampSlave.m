%#ok<*DEFNU>
function IsoSampSlave()
global udpCom gl
udpCom.port = 6665;
udpCom.sock = nan;
udpCom.rexip = '192.168.1.120';

KbName('UnifyKeyNames');

[udpCom.sock, okay] = pnetStart(udpCom.port);
if ~okay, return; end

pnet(udpCom.sock, 'setreadtimeout', 0);
pnet(udpCom.sock, 'setwritetimeout', 0);
fprintf('In %s\n', mfilename);

gl.vpixx = IsVPixx();
gl.ccmode = 1; % in mode '1': every 2nd column of pixels is ignored

gl.fliprequest = false;
gl.framecounter = 0;
gl.framecountermax = 0;

gl.fp.on = false;
gl.fp.x = 0;
gl.fp.y = 0;
gl.fp.size = 0;
gl.fp.rgb = [0 0 0];
gl.fp.drawrect = [0 0 0 0];

gl.targ.on = false;
gl.targ.rgb = []; % relative to background; set in InitDisplay
gl.targ.drawrects = zeros(2,4);

gl.stim.on = false;
gl.stim.template = [];

timetoleave = false;
while ~timetoleave
    timetoleave = deal_with_messages();
    
    if gl.stim.on
        DrawStim();
    end
    if gl.fp.on
        DrawFP();
    end
    if gl.targ.on
        DrawTargs();
    end
    if gl.fliprequest
        DoFlip();
    end
    
    if check_for_quit()
        sca();
        pnet(udpCom.sock, 'close');
        return
    end
end

function yn = check_for_quit()
yn = false;
[keyisdown,nil,keycode] = KbCheck(-1); %#ok<ASGLU>
if keyisdown
    yn = keycode(KbName('ESCAPE'));
end

function DoFlip()
global gl

Screen('Flip', gl.windowPtr);
if gl.stim.on
    if ~gl.framecounter
        message_REX('MACSTIMON');
    elseif gl.framecounter == gl.framecountermax - 1
        message_REX('MACSTIMOFF');
        gl.stim.on = false;
    end
    gl.framecounter = gl.framecounter + 1;
end
gl.fliprequest = false;

function InitDisplay(mondist, screenwidth, varargin)
global gl

calfilename = varargin{1};
load(calfilename);
cal = cals{end}; %#ok<USENS> it comes from loading the above file
gl.mondistcm = mondist;
gl.screenWidthcm = screenwidth;
gl.bkgndRGB = round(255 * cal.bgColor)';
gl.cal.gammaTable = cal.gammaTable;
gl.cal.monSpd = cal.P_device;
gl.cal.invgammaTable = InvertGammaTable(cal.gammaInput, cal.gammaTable, 2^16);

if nargin > 3 && ~isempty(varargin{2})
    fundfilename = varargin{2};
    s = load(fundfilename);
    fns = fieldnames(s);
    gl.cal.fundamentals = s.(fns{1})';
    wavelength_spacing = s.(fns{2});
    P_device = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls(wavelength_spacing));
    gl.cal.M = gl.cal.fundamentals'*P_device;
    gl.cal.invM = inv(gl.cal.M);
end

% start up the imaging pipeline
if ~isempty(Screen('Windows'))
    gl.windowPtr = max(Screen('Windows'));
    Screen('FillRect', gl.windowPtr, cal.bgColor);
else
    PsychImaging('PrepareConfiguration');
    if gl.vpixx
        PsychImaging('AddTask', 'General', 'EnableDataPixxC48Output', gl.ccmode);
    else
        PsychImaging('AddTask', 'General', 'EnableBits++Color++Output', gl.ccmode);
    end
    gl.windowPtr = PsychImaging('OpenWindow', 0, cal.bgColor);
    
    [gl.screenWidthpix, gl.screenHeightpix] = Screen('WindowSize', gl.windowPtr);
    gl.screenCenterXpix = gl.screenWidthpix/2;
    gl.screenCenterYpix = gl.screenHeightpix/2;
    
    pixpercm = gl.screenWidthpix/gl.screenWidthcm;
    theta = atan2(gl.screenWidthcm/2, gl.mondistcm)*180/pi;
    cmperdeg = gl.screenWidthcm/2/theta;
    gl.pixperdeg = pixpercm*cmperdeg;
    
    gl.framerate = 1/Screen('GetFlipInterval', gl.windowPtr);
end

gl.bkgndrgb = [gl.cal.gammaTable(gl.bkgndRGB(1) + 1,1)
    gl.cal.gammaTable(gl.bkgndRGB(2) + 1,2)
    gl.cal.gammaTable(gl.bkgndRGB(3) + 1,3)];

gl.targ.rgb = (gl.bkgndRGB-30)/255;
gl.stim.on = false;
gl.fp.on = false;
gl.targ.on = false;

HideCursor();
message_REX('DISPLAYINIT');

function PrepareGabor(flash_time, theta, sf, phi, sigma, ggamma, tf, nsigmas)
global gl

lambda = 1 / sf;
nframes = ceil(gl.framerate * flash_time / 1000);
ramplength = ceil(nframes / 4);
ramp = linspace(0, 1, ramplength);
plateau = ones(1, nframes - 2 * ramplength);
temporalprofile = [ramp plateau fliplr(ramp)];
stimsizeindeg = sigma * nsigmas; % This is half stim size (we we go out +/- nsigmas)
stimsizeinpix = round(stimsizeindeg * gl.pixperdeg);  % full stim size in doublewide pixels
[x,y] = meshgrid(stimsizeindeg * linspace(-1, 1, stimsizeinpix), ...
    stimsizeindeg * linspace(-1, 1, 2 * stimsizeinpix));
% x and y are in dva
X = x * cos(-theta) + y * sin(-theta);
Y =-x * sin(-theta) + y * cos(-theta);

deltaphase = tf * 2 * pi / gl.framerate;
phases = phi + (0:nframes-1)*deltaphase;
phases = reshape(phases, [1 1 nframes]);
temporalprofile = reshape(temporalprofile, [1 1 nframes]);
expterm = bsxfun(@times, exp(-(X.^2 + ggamma^2 * Y.^2) / 2 / sigma^2), temporalprofile);
costerm = cos(bsxfun(@plus, 2 * pi * Y / lambda, phases));
gl.stim.template = expterm .* costerm;
message_REX('GABORREADY');

function ShowFP(x, y, size, fpr, fpg, fpb)
global gl

gl.fp.x = x / 10;
gl.fp.y = y / 10;
gl.fp.size = size / 10;
gl.fp.rgb = [fpr fpg fpb]/255;

xx = gl.fp.x * gl.pixperdeg;
yy = gl.fp.y * gl.pixperdeg;

fpsizeinpix = round(gl.pixperdeg * gl.fp.size);

gl.fp.drawrect = [gl.screenCenterXpix + xx - floor(fpsizeinpix / 2) + 1, ...
    gl.screenCenterYpix - yy - floor(fpsizeinpix / 2) + 1, ...
    gl.screenCenterXpix + xx + ceil(fpsizeinpix / 2), ...
    gl.screenCenterYpix - yy + ceil(fpsizeinpix / 2)];
gl.fp.on = true;

function DrawFP()
global gl
Screen('FillRect', gl.windowPtr, gl.fp.rgb, gl.fp.drawrect);
gl.fliprequest = true;

function HideFP()
global gl
gl.fp.on = false;
gl.fliprequest = true;

function ShowStim(stimx, stimy, l, m, s)
global gl

stimx = stimx / 10;
stimy = stimy / 10;
stimconecontrast = [l m s];

% Creating the drawing window
x = stimx * gl.pixperdeg;
y = stimy * gl.pixperdeg;
stimsizeinpix = size(gl.stim.template, 2); % in doublewide pix

gl.drawrect = round([gl.screenCenterXpix + x - stimsizeinpix ...
    gl.screenCenterYpix - y - stimsizeinpix ...
    gl.screenCenterXpix + x + stimsizeinpix ...
    gl.screenCenterYpix - y + stimsizeinpix]);

if rem(gl.drawrect(1), 2) %if the rectangle starts on an odd pixel
    gl.drawrect(1) = gl.drawrect(1) - 1;
    gl.drawrect(3) = gl.drawrect(3) - 1;
end

assert(gl.drawrect(3) - gl.drawrect(1) == 2 * stimsizeinpix ...
    && gl.drawrect(4) - gl.drawrect(2) == 2 * stimsizeinpix, ...
    'Incorrect draw window size!');

% Getting RGBs
bkgndlms = gl.cal.M * gl.bkgndrgb;
gl.stim.rgb = gl.cal.invM * (bkgndlms .* (1 + stimconecontrast'));

% Checking for out of gamut errors
% Right now just squeezing it back into the gamut - need to send a
% message to REX that this was done.
lims = min([gl.bkgndrgb'; 1 - gl.bkgndrgb']);
if any(abs(gl.stim.rgb - gl.bkgndrgb) > lims')
    scalefactor = max(abs(gl.stim.rgb - gl.bkgndrgb) ./ lims');  % +(1/NGAMMASTEPS);
    fprintf('Squeezing back into gamut by a factor of %f!\n', scalefactor);
    gl.stim.rgb = (gl.stim.rgb - gl.bkgndrgb) ./ scalefactor + gl.bkgndrgb;
end

gl.framecounter = 0;
gl.framecountermax = size(gl.stim.template, 3);
gl.stim.on = true;

function DrawStim()
global gl

NGAMMASTEPS = size(gl.cal.invgammaTable, 1);
frame = gl.stim.template(:,:,gl.framecounter + 1);
img = zeros(size(gl.stim.template, 1), size(gl.stim.template, 2), 3);

for plane = 1:3
    tmp = frame * (gl.stim.rgb(plane) - gl.bkgndrgb(plane)) + gl.bkgndrgb(plane);
    tmp = round(tmp * (NGAMMASTEPS - 1)) + 1;
    tmp = gl.cal.invgammaTable(tmp,plane);
    img(:,:,plane) = reshape(tmp, size(img(:,:,1)));
end

tex = Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(img, [], gl.ccmode), [], [], 2);
Screen('DrawTexture', gl.windowPtr, tex, [], gl.drawrect, [], 0);
Screen('Close', tex);
gl.fliprequest = true;

function HideStim()
global gl
gl.stim.on = false;
gl.fliprequest = true;

function ShowTargs(correct_x, correct_y, targsize)
global gl

xx = correct_x/10*gl.pixperdeg;
yy = correct_y/10*gl.pixperdeg;
fpsizeinpix = round(gl.pixperdeg*targsize/10);

gl.targ.drawrects(1,:) = [
    gl.screenCenterXpix+xx-floor(fpsizeinpix/2)+1, ...
    gl.screenCenterYpix-yy-floor(fpsizeinpix/2)+1, ...
    gl.screenCenterXpix+xx+ceil(fpsizeinpix/2), ...
    gl.screenCenterYpix-yy+ceil(fpsizeinpix/2)];
gl.targ.drawrects(2,:) = gl.targ.drawrects(1,:) + [-2*xx 0 -2*xx 0];
gl.targ.on = true;

function DrawTargs()
global gl
Screen('FillRect', gl.windowPtr, gl.targ.rgb, gl.targ.drawrects');
gl.fliprequest = true;

function AllOff()
global gl
gl.stim.on = false;
gl.targ.on = false;
gl.fp.on = false;
gl.fliprequest = true;

function message_REX(message)
global udpCom
pnet(udpCom.sock, 'write', [message '>> >>']);
pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);

function leave = deal_with_messages()
global udpCom gl %#ok<NUSED>

leave = false;
bytes_avail = pnet(udpCom.sock, 'readpacket', 1000, 'noblock');
if bytes_avail
    message = pnet(udpCom.sock, 'read', bytes_avail, 'char');
    if strncmp(message, 'return', 6)
        stk = dbstack('-completenames'); % Check whether called from another function or from command line
        if ~strcmp(stk(end).name, mfilename)
            leave = true;
        end
    end
    try
        eval(message);
    catch ME
        fprintf('Trouble with message: "%s"\n', message);
        disp(getReport(ME));
    end
end
