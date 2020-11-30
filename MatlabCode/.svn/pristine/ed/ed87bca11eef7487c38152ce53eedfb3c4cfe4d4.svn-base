%#ok<*DEFNU>
function ThreePulseMonteSlave()
% Communication variables
global udpCom gl
udpCom.port = 6665;
udpCom.sock = nan;
udpCom.rexip = '192.168.1.120';

KbName('UnifyKeyNames');

gl.fliprequest = false;
gl.windowPtr = 0;
gl.framecounter = 0;
gl.framecountermax = 0;

gl.vpixx = IsVPixx();
gl.ccmode = 1;

gl.fp.on = false;
gl.fp.x = 0;
gl.fp.y = 0;
gl.fp.size = 0;
gl.fp.rgb = [0 0 0];
gl.fp.drawrect = [0 0 0 0];

gl.vis.on = false;
gl.vis.x = 0;
gl.vis.y = 0;
gl.vis.size = 0;
gl.vis.rgb = [0 0 0];
gl.vis.drawrect = [0 0 0 0];

gl.stim.on = false;

gl.targ.on = false;
gl.targ.x = 0;
gl.targ.y = 0;
gl.targ.size = 0;
gl.targ.rgb = [0 0 0];
gl.targ.drawrects = zeros(2,4);
gl.targ.wrong_contrast = 1;

[udpCom.sock, Success] = pnetStart(udpCom.port);
if ~Success, return, end

pnet(udpCom.sock, 'setreadtimeout', 0);
pnet(udpCom.sock, 'setwritetimeout', 0);
fprintf('In %s\n', mfilename);

timetoleave = false;
while ~timetoleave
    timetoleave = deal_with_messages();
    if gl.windowPtr > 0
        if gl.fp.on
            DrawFP();
        end
        if gl.vis.on
            DrawVis();
        end
        if gl.targ.on
            DrawTargs();
        end
        if gl.stim.on
            DrawStim();
        end
        if gl.fliprequest
            DoFlip();
        end
    end
    if check_for_quit
        close_textures();
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

function close_textures()
global gl
if isfield(gl, 'tex')
    Screen('Close', gl.tex);
end

function DoFlip()
global gl

Screen('Flip', gl.windowPtr);
if gl.stim.on
    gl.framecounter = gl.framecounter + 1;
    if gl.framecounter == 1
        message_REX('MACSTIMON');
    elseif gl.framecounter >= gl.framecountermax
        message_REX('MACSTIMOFF');
        gl.stim.on = false;
    end
end
gl.fliprequest = false;

function InitDisplay(mondistcm, screenwidthcm, varargin) % varargin{1} and {2} are defined to be the calibration and fundamentals
global gl

calfilename = varargin{1};
load(calfilename);
cal = cals{end}; %#ok<USENS>
gl.mondistcm = mondistcm;
gl.screenWidthcm = screenwidthcm;
gl.bkgndRGB = round(255*cal.bgColor)';
gl.cal.gammaTable = cal.gammaTable;
gl.cal.monSpd = cal.P_device;

gl.cal.invgammaTable = InvertGammaTable(cal.gammaInput, gl.cal.gammaTable, 2^16);

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
    close_textures();
    gl.windowPtr = max(Screen('Windows'));
    Screen('FillRect', gl.windowPtr, gl.bkgndRGB/255);
else
    PsychImaging('PrepareConfiguration');
    if gl.vpixx
        PsychImaging('AddTask', 'General', 'EnableDataPixxC48Output', gl.ccmode); % in mode '1': every 2nd column of pixels is ignored
    else
        PsychImaging('AddTask', 'General', 'EnableBits++Color++Output', gl.ccmode);
    end
    gl.windowPtr = PsychImaging('OpenWindow', 0, gl.bkgndRGB/255);

    [gl.screenWidthpix, gl.screenHeightpix] = Screen('WindowSize', gl.windowPtr);
    gl.screenCenterXpix = gl.screenWidthpix/2;
    gl.screenCenterYpix = gl.screenHeightpix/2;

    pixpercm = gl.screenWidthpix/screenwidthcm;
    theta = atan2(screenwidthcm/2, mondistcm)*180/pi;
    cmperdeg = screenwidthcm/2/theta;
    gl.pixperdeg = pixpercm*cmperdeg;

    gl.framerate = 1/Screen('GetFlipInterval', gl.windowPtr);
end

gl.targ.on = false;
gl.fp.on = false;
gl.vis.on = false;
gl.stim.on = false;

gl.bkgndrgb = [gl.cal.gammaTable(gl.bkgndRGB(1)+1,1);...
    gl.cal.gammaTable(gl.bkgndRGB(2)+1,2);...
    gl.cal.gammaTable(gl.bkgndRGB(3)+1,3)];

HideCursor();
message_REX('DISPLAYINIT');

% "Show" functions are called from REX.  They set up the appropriate fields
% in the "gl" structure (and set the "...on" toggle field to 1).
function ShowFP(x, y, size, fprgb)
global gl

if nargin < 4 || isempty(fprgb)
    gl.fp.rgb = [0 0 0];
else
    gl.fp.rgb = fprgb/255;
end

gl.fp.x = x/10;
gl.fp.y = y/10;
gl.fp.size = size/10;

xx = gl.fp.x*gl.pixperdeg;
yy = gl.fp.y*gl.pixperdeg;

fpsizeinpix = round(gl.pixperdeg*gl.fp.size);

gl.fp.drawrect = [gl.screenCenterXpix + xx - floor(fpsizeinpix/2) + 1, ...
    gl.screenCenterYpix - yy - floor(fpsizeinpix/2) + 1, ...
    gl.screenCenterXpix + xx + ceil(fpsizeinpix/2)    , ...
    gl.screenCenterYpix - yy + ceil(fpsizeinpix/2)];

gl.fp.on = true;

% "Draw" functions get called on every screen refresh so long as the
% corresponding graphical object is to be displayed.  These are the
% functions that actually draw the object to the screen (and they make a
% fliprequest so that the drawn object will actually appear).
function DrawFP()
global gl

Screen('FillRect', gl.windowPtr, gl.fp.rgb, gl.fp.drawrect);
gl.fliprequest = true;

function HideFP()
global gl

gl.fp.on = false;
gl.fliprequest = true;

function ShowVis(x, y, size, visrgb)
global gl

gl.vis.rgb = visrgb/255;
gl.vis.x = x/10;
gl.vis.y = y/10;
gl.vis.size = size/10;

xx = gl.vis.x*gl.pixperdeg;
yy = gl.vis.y*gl.pixperdeg;

vissizeinpix = round(gl.pixperdeg*gl.vis.size);

gl.vis.drawrect = [gl.screenCenterXpix + xx - floor(vissizeinpix/2) + 1, ...
    gl.screenCenterYpix - yy - floor(vissizeinpix/2) + 1, ...
    gl.screenCenterXpix + xx + ceil(vissizeinpix/2)    , ...
    gl.screenCenterYpix - yy + ceil(vissizeinpix/2)];

gl.vis.on = true;

function DrawVis()
global gl

Screen('FillRect', gl.windowPtr, gl.vis.rgb, gl.vis.drawrect);
gl.fliprequest = true;

function HideVis()
global gl

gl.stim.on = false;
gl.vis.on = false;
gl.fliprequest = true;

function PrepareGabor(flash_time, tf, theta, sf, phi, ggamma, sigma, nsigmas)
global gl

lambda = 1 / sf;
nframes = ceil(gl.framerate * flash_time / 1000);
ramplength = ceil(nframes / 4);
ramp = linspace(0, 1, ramplength);
plateau = ones(1, nframes - 2 * ramplength);
temporalprofile = [ramp plateau fliplr(ramp)];
nframes = length(temporalprofile);
stimsizeindeg = sigma * nsigmas; % This is half stim size (we go out +/- nsigmas)
stimsizeinpix = round(stimsizeindeg * gl.pixperdeg); % full stim size in doublewide pixels
[x,y] = meshgrid(stimsizeindeg * linspace(-1, 1, stimsizeinpix), ...
    stimsizeindeg * linspace(-1, 1, 2 * stimsizeinpix));
% x and y are in dva
X = x * cos(-theta) + y * sin(-theta);
Y =-x * sin(-theta) + y * cos(-theta);

temporalprofile = reshape(temporalprofile, [1 1 nframes]);
expterm = bsxfun(@times, exp(-(X.^2 + ggamma^2 * Y.^2) / 2 / sigma^2), temporalprofile);
deltaphase = 2*pi*tf/gl.framerate;
phase = reshape(phi + (0:nframes-1)*deltaphase, [1 1 nframes]);
costerm = cos(bsxfun(@plus, 2*pi*Y/lambda, phase));
gl.stim.template = expterm .* costerm;
message_REX('GABORREADY');

function ShowStim(stimx, stimy, stimconecontrast)
global gl

stimx = stimx/10;
stimy = stimy/10;

% Creating the drawing window
x = stimx*gl.pixperdeg;
y = stimy*gl.pixperdeg;
[frame_size(1),frame_size(2),gl.framecountermax] = size(gl.stim.template);

gl.drawrect = round([gl.screenCenterXpix+x-frame_size(2), ...
    gl.screenCenterYpix-y-frame_size(2), ...
    gl.screenCenterXpix+x+frame_size(2), ...
    gl.screenCenterYpix-y+frame_size(2)]);

if rem(gl.drawrect(1), 2) %if the rectangle starts on an odd pixel
    gl.drawrect(1) = gl.drawrect(1) - 1;
    gl.drawrect(3) = gl.drawrect(3) - 1;
end

assert(gl.drawrect(3)-gl.drawrect(1) == 2*frame_size(2) ...
    && gl.drawrect(4)-gl.drawrect(2) == 2*frame_size(2), ...
    'Incorrect draw window size!');

% Getting RGBs
bkgndlms = gl.cal.M*gl.bkgndrgb;
gl.stim.rgb = gl.cal.invM*(bkgndlms.*(1+stimconecontrast'));

% Checking for out of gamut errors
lims = min([gl.bkgndrgb'; 1-gl.bkgndrgb']);
if any(abs(gl.stim.rgb-gl.bkgndrgb) > lims')
    message_REX('OOG');
    scalefactor = max(abs(gl.stim.rgb-gl.bkgndrgb)./lims'); % +(1/NGAMMASTEPS);
    gl.stim.rgb = (gl.stim.rgb-gl.bkgndrgb)./scalefactor+gl.bkgndrgb;
end

img = zeros([frame_size 3]);
NGAMMASTEPS = size(gl.cal.invgammaTable, 1);
gl.tex = zeros(1,gl.framecountermax);
for nth_frame = 1:gl.framecountermax
    frame = gl.stim.template(:,:,nth_frame);
    for plane = 1:3
        int_to_volt = round((frame*(gl.stim.rgb(plane)-gl.bkgndrgb(plane))+gl.bkgndrgb(plane)) * (NGAMMASTEPS-1))+1;
        img(:,:,plane) = reshape(gl.cal.invgammaTable(int_to_volt,plane), frame_size);
    end
    gl.tex(nth_frame) = Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(img, [], gl.ccmode), [], [], 2);
end

gl.framecounter = 0;
gl.stim.on = true;

function DrawStim()
global gl

Screen('DrawTexture', gl.windowPtr, gl.tex(gl.framecounter+1), [], gl.drawrect, [], 0);
% textures are closed in AllOff
gl.fliprequest = true;

function ShowTargs(correct_x, correct_y, targsize, wrong_targ_contrast, color_targets)
global gl

xx = correct_x/10*gl.pixperdeg;
yy = correct_y/10*gl.pixperdeg;
fpsizeinpix = round(gl.pixperdeg*targsize/10);

gl.targ.drawrects(1,:) = [
    gl.screenCenterXpix+xx-floor(fpsizeinpix/2)+1, ...
    gl.screenCenterYpix-yy-floor(fpsizeinpix/2)+1, ...
    gl.screenCenterXpix+xx+ceil(fpsizeinpix/2), ...
    gl.screenCenterYpix-yy+ceil(fpsizeinpix/2)];
gl.targ.drawrects(2,:) = gl.targ.drawrects(1,:) + [-2*xx 2*yy -2*xx 2*yy];
gl.targ.wrong_contrast = wrong_targ_contrast;
gl.targ.color_targets = color_targets;
gl.targ.on = true;

function DrawTargs()
global gl

if gl.targ.color_targets
    if gl.targ.drawrects(1) < gl.screenCenterXpix
        targ_colors = [1 1 0; 0 0 1]';
    else
        targ_colors = [0 0 1; 1 1 0]';
    end
    Screen('FillRect', gl.windowPtr, targ_colors, gl.targ.drawrects');
else
    Screen('FillRect', gl.windowPtr, gl.targ.rgb, gl.targ.drawrects(1,:)); % correct target
    Screen('FillRect', gl.windowPtr, gl.bkgndRGB/255*(1-gl.targ.wrong_contrast), gl.targ.drawrects(2,:)); % wrong target
end
gl.fliprequest = true;

function AllOff()
global gl

close_textures();
gl.stim.on = false;
gl.targ.on = false;
gl.fp.on = false;
gl.vis.on = false;
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
