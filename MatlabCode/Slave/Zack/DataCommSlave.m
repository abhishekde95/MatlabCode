%#ok<*DEFNU>
function DataCommSlave()
% Communication variables
global udpCom gl
udpCom.port = 6665;
udpCom.sock = nan;
udpCom.rexip = '192.168.1.120';

KbName('UnifyKeyNames');

[udpCom.sock, success] = pnetStart(udpCom.port);
if ~success, return, end

pnet(udpCom.sock, 'setreadtimeout', 0);
pnet(udpCom.sock, 'setwritetimeout', 0);
fprintf('In %s\n', mfilename);

gl.vpixx = IsVPixx();
gl.ccmode = 1;

gl.fliprequest = 0;
gl.timetoleave = 0;
gl.windowPtr = 0;

gl.fp.on = 0;
gl.fp.x = 0;
gl.fp.y = 0;
gl.fp.size = 0;
gl.fp.rgb = [0 0 0];
gl.fp.drawrect = [0 0 0 0];

gl.targ.on = 0;
gl.targ.x = 0;
gl.targ.y = 0;
gl.targ.size = 0;
gl.targ.rgb = [0 0 0];
gl.targ.drawrect = [0 0 0 0];

timetoleave = false;
while ~timetoleave
    timetoleave = deal_with_messages();
    
    if gl.targ.on
        DrawTarg();
    end
    if gl.fp.on
        DrawFP();
    end
    if gl.fliprequest
        DoFlip();
    end
    
    if check_for_quit()
        sca();
        pnet('closeall');
        return
    end
end

function yn = check_for_quit()
yn = false;
[keyisdown,nil,keycode] = KbCheck(-1); %#ok<ASGLU>
if keyisdown
    yn = keycode(KbName('ESCAPE'));
end

% DoFlip() taken from WhiteNoise
function DoFlip() % tidy this up to work with targ
global gl

Screen('Flip', gl.windowPtr);
gl.fliprequest = 0;

% varargin{1} and {2} are expected to be the calibration and fundamentals filename.
function InitDisplay(mondistcm, screenwidthcm, varargin)
global gl

assert(nargin > 2 && ~isempty(varargin{1}), ...
    'You need to provide a calibration file name!');

calfilename = varargin{1};
cal = load(calfilename);
cal = cal.cals{end};
gl.bkgndRGB = round(255*cal.bgColor)';
gl.cal.gammaTable = cal.gammaTable;
gl.cal.monSpd = cal.P_device;
gl.cal.invgammaTable = InvertGammaTable(cal.gammaInput, gl.cal.gammaTable, 2^16);

if nargin > 3 && ~isempty(varargin{2})
    fundfilename = varargin{2};
    s = load(fundfilename);
    fns = fieldnames(s);
    gl.cal.fundamentals = s.(fns{1})';
    P_device = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls([380 5 81]));
    gl.cal.M = gl.cal.fundamentals'*P_device;
    gl.cal.invM = inv(gl.cal.M);
end

% start up the imaging pipeline
if ~isempty(Screen('Windows'))
    gl.windowPtr = max(Screen('Windows'));
    Screen('FillRect', gl.windowPtr, gl.bkgndRGB/255);
else
    PsychImaging('PrepareConfiguration');
    if gl.vpixx
        PsychImaging('AddTask', 'General', 'EnableDataPixxC48Output',   gl.ccmode);
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

gl.bkgndrgb = [gl.cal.gammaTable(gl.bkgndRGB(1)+1,1);...
    gl.cal.gammaTable(gl.bkgndRGB(2)+1,2);...
    gl.cal.gammaTable(gl.bkgndRGB(3)+1,3)];

gl.fp.on = false;
gl.targ.on = false;
HideCursor();
message_REX('DISPLAYINIT');

function InitHeaderVars()
global gl
gl.hd.i = int8(-112);
gl.hd.d = exp(-1);
gl.hd.f = single(-pi);
gl.hd.l = int32(-89457349);
gl.hd.c = 'J';
gl.hd.s = int16(-23784);
gl.hd.us = uint16(33338);
gl.hd.ul = uint32(123456789);
gl.hd.ui = uint32(987654321);
gl.hd.ll = -int64(2)^56+29354209347;
gl.hd.ull = uint64(2)^61-239278469283746;

old_state = rng(11111111);
gl.hd.ai = randi([intmin('int8') intmax('int8')], 20, 1, 'int8');
gl.hd.ad = randn(768, 1);
gl.hd.af = randn(722, 1, 'single');
gl.hd.al = randi([intmin('int32') intmax('int32')], 72, 1, 'int32');
gl.hd.ac = '  this   is quite a string! ! ';
gl.hd.as = randi([intmin('int16') intmax('int16')], 52, 1, 'int16');
gl.hd.aus = randi(intmax('uint16'), 123, 1, 'uint16');
gl.hd.aul = randi(intmax('uint32'), 191, 1, 'uint32');
gl.hd.aui = randi(intmax('uint32'), 222, 1, 'uint32');
gl.hd.all = typecast(randi([intmin('int32') intmax('int32')], 2*435, 1, 'int32'), 'int64');
gl.hd.aull = typecast(randi(intmax('uint32'), 2*511, 1, 'uint32'), 'uint64');
rng(old_state);
message_REX('DISPLAYINIT'); % lazy... reuse this message

function InitTrialVars()
global gl
gl.tp.i = int8(-112);
gl.tp.d = exp(-1);
gl.tp.f = single(-pi);
gl.tp.l = int32(-89457349);
gl.tp.c = 'J';
gl.tp.s = int16(-23784);
gl.tp.us = uint16(33338);
gl.tp.ul = uint32(123456789);
gl.tp.ui = uint32(987654321);
gl.tp.ll = -int64(2)^56+29354209347;
gl.tp.ull = uint64(2)^61-239278469283746;

old_state = rng(22222222);
gl.tp.ai = randi([intmin('int8') intmax('int8')], 20, 1, 'int8');
gl.tp.ad = rand(800,1);
gl.tp.af = rand(799, 1, 'single');
gl.tp.al = randi([intmin('int32') intmax('int32')], 72, 1, 'int32');
gl.tp.ac = '  this   is quite a string! ! ';
gl.tp.as = randi([intmin('int16') intmax('int16')], 52, 1, 'int16');
gl.tp.aus = randi(intmax('uint16'), 123, 1, 'uint16');
gl.tp.aul = randi(intmax('uint32'), 191, 1, 'uint32');
gl.tp.aui = randi(intmax('uint32'), 222, 1, 'uint32');
gl.tp.all = typecast(randi([intmin('int32') intmax('int32')], 2*435, 1, 'int32'), 'int64');
gl.tp.aull = typecast(randi(intmax('uint32'), 2*511, 1, 'uint32'), 'uint64');
rng(old_state);
message_REX('DISPLAYINIT'); % lazy... reuse this message

% "Show" functions are called from REX.  They set up the appropriate fields
% in the "gl" structure (and set the "...on" toggle field to 1).
function ShowFP(x, y, size, fpr, fpg, fpb)
global gl

gl.fp.x = x/10;
gl.fp.y = y/10;
gl.fp.size = size/10;
gl.fp.rgb = [fpr fpg fpb]/255;

xx = gl.fp.x*gl.pixperdeg;
yy = gl.fp.y*gl.pixperdeg;

fpsizeinpix = round(gl.pixperdeg*gl.fp.size);

gl.fp.drawrect = [gl.screenCenterXpix + xx - floor(fpsizeinpix/2) + 1, ...
    gl.screenCenterYpix - yy - floor(fpsizeinpix/2) + 1, ...
    gl.screenCenterXpix + xx + ceil(fpsizeinpix/2)    , ...
    gl.screenCenterYpix - yy + ceil(fpsizeinpix/2)];

gl.fp.on = 1;

% "Draw" functions get called on every screen refresh so long as the
% corresponding graphical object is to be displayed.  These are the
% functions that actually draw the object to the screen (and they make a
% fliprequest so that the drawn object will actually appear).
function DrawFP()
global gl

Screen('FillRect', gl.windowPtr, gl.fp.rgb, gl.fp.drawrect);
gl.fliprequest = 1;

function HideFP()
global gl

gl.fp.on = 0;
gl.fliprequest = 1;

function ShowTarg(x, y, size, fpr, fpg, fpb)
global gl

gl.targ.x = x/10;
gl.targ.y = y/10;
gl.targ.size = size/10;
gl.targ.rgb = [fpr fpg fpb]/255;

xx = gl.targ.x*gl.pixperdeg;
yy = gl.targ.y*gl.pixperdeg;

targsizeinpix = round(gl.pixperdeg*gl.targ.size);

gl.targ.drawrect = ...
    [gl.screenCenterXpix + xx - floor(targsizeinpix/2) + 1 ...
    gl.screenCenterYpix - yy - floor(targsizeinpix/2) + 1 ...
    gl.screenCenterXpix + xx + ceil(targsizeinpix/2)     ...
    gl.screenCenterYpix - yy + ceil(targsizeinpix/2)];

gl.targ.on = 1;

function DrawTarg()
global gl

Screen('FillRect', gl.windowPtr, gl.targ.rgb, gl.targ.drawrect);
gl.fliprequest = 1;

function HideTarg()
global gl

gl.targ.on = 0;
gl.fliprequest = 1;

function AllOff()
global gl

gl.targ.on = 0;
gl.fp.on = 0;
gl.fliprequest = 1;

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
