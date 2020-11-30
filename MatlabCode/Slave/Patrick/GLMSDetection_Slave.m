%#ok<*DEFNU>
function GLMSDetection_Slave
    % Communication variables
    global udpCom KEY gl
    udpCom.port = 6665;
    udpCom.sock = nan;
    udpCom.rexip = '192.168.1.120';

    % Keypress constants
    KbName('UnifyKeyNames');
    KEY.ESC = KbName('escape');

    % Variables that are destructively modified by subfunctions
    % Because we are guaranteed to be running this code after having
    % already
    % setup the global gl structure (via the white noise paradigm)
    % there's no sense in reassigning everything.
    gl.usecolourmode = 1;  % Are we using the Bits++

    gl.fliprequest = 0;
    gl.timetoleave = 0;

    gl.fp.on = 0;
    gl.fp.x = 0;
    gl.fp.y = 0;
    gl.fp.size = 0;
    gl.fp.rgb = [0 0 0];
    gl.fp.drawrect = [0 0 0 0];

    gl.vpixx = IsVPixx();
    gl.ccmode = 1;

    gl.targ.on = 0;

    gl.stim.framecounter = 0;
    gl.stim.on = 0;
    gl.stim.drawrect = [0 0 0 0];
    gl.stim.textureidx = [];

    gl.gridX = [];
    gl.gridY = [];

    [udpCom.sock, Success] = pnetStart(udpCom.port);
    if ~Success, return, end

    pnet(udpCom.sock, 'setreadtimeout', 0);
    pnet(udpCom.sock, 'setwritetimeout', 0);

    % Handshaking (to indicate slave program is initiated)
    pnet(udpCom.sock, 'write', 'macready>> >>');
    pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);

    messageIsAvailable = 0;
    while ~messageIsAvailable  % endless loop.
        messageIsAvailable = pnet(udpCom.sock, 'readpacket', 1000, 'no block');

        if messageIsAvailable
            DealWithMessage(messageIsAvailable);
            messageIsAvailable = 0;
        end

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

        if gl.timetoleave, return, end

        [keyisdown,junk,keycode] = KbCheck();
        if keyisdown && keycode(KEY.ESC)
            ShowCursor();
            pnet(udpCom.sock, 'close');
            sca;
            gl.timetoleave = 1;
        end
    end
end


% DoFlip() taken from WhiteNoise
function DoFlip() % tidy this up to work with targ
    global gl

    Screen('Flip', gl.windowPtr, 0, 0);
%     if gl.stim.on == 1 && gl.stim.framecounter == gl.stim.nframes
%           Handshaking
%           pnet(udpCom.sock, 'write', 'macready>> >>');
%           pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
%           disp('send1')
%     end
    if gl.stim.on == 1
        gl.stim.framecounter = gl.stim.framecounter + 1;
    end
    gl.fliprequest = 0;

end

function InitDisplay(mondist, screenwidth, varargin) % varargin{1} and {2} are defined to be the calibration and fundamentals
    global gl udpCom

    [junk,calfilename] = fileparts(varargin{1});
    calfilename = [calfilename '.mat'];
    
    if (gl.vpixx && ~strcmp(calfilename, 'ProPixx.mat')) || (~gl.vpixx && strcmp(calfilename, 'ProPixx.mat'))
        warning(['There is a mismatch between the requested calibration structure and whether' ...
            ' the ProPixx is physically hooked up to this Slave machine']);
        gl.timetoleave = 1;
    end

    load(calfilename);
    cal = cals{end}; %#ok<USENS>
    gl.mondistcm = mondist;
    gl.screenWidthcm = screenwidth;
    gl.bkgndRGB = round(255*cal.bgColor)';
    gl.cal.gammaTable = cal.gammaTable;
    gl.cal.monSpd = cal.P_device;
    gl.cal.invgammaTable = InvertGammaTable(cal.gammaInput, gl.cal.gammaTable, 2^16);

    if nargin > 3 && ~isempty(varargin{2})
        [junk,fundfilename] = fileparts(varargin{2});
        fundfilename = [fundfilename '.mat'];
        s = load(fundfilename);
        fns = fieldnames(s);
        gl.cal.fundamentals = s.(fns{1})';
        wavelength_spacing = s.(fns{2});
        P_device = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls(wavelength_spacing));
        gl.cal.M = gl.cal.fundamentals'*P_device;
        gl.cal.invM = inv(gl.cal.M);
    end

    gl.bkgndrgb = [gl.cal.gammaTable(gl.bkgndRGB(1)+1,1),...
                   gl.cal.gammaTable(gl.bkgndRGB(2)+1,2),...
                   gl.cal.gammaTable(gl.bkgndRGB(3)+1,3)]';

    % start up the imaging pipeline
    if ~isempty(Screen('Windows'))
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
    end

    gl.framerate = Screen('NominalFrameRate', gl.windowPtr, 1);
    [screenwidthpix, screenheightpix] = Screen('WindowSize', gl.windowPtr);
    gl.screenWidthpix = screenwidthpix;
    gl.screenHeightpix = screenheightpix;
    gl.screenCenterXpix = screenwidthpix/2;
    gl.screenCenterYpix = screenheightpix/2;

    pixpercm = gl.screenWidthpix/gl.screenWidthcm;
    theta = atan2(gl.screenWidthcm/2, gl.mondistcm)*180/pi;
    cmperdeg = gl.screenWidthcm/(2*theta);
    gl.pixperdeg = pixpercm*cmperdeg;

    gl.fp.on = 0;
    gl.targ.on = 0;  % probably want to change this
    gl.targ.rgb = (gl.bkgndRGB-30)/255;
    gl.fp.on = 0;
    gl.stim.on = 0;


    HideCursor();

    % Handshaking (to indicate display function is initiated)
    pnet(udpCom.sock, 'write', 'macready>> >>');
    pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);

end

% "Show" functions are called from REX.  They set up the appropriate fields
% in the "gl" structure (and set the "...on" toggle field to 1).
function ShowFP(x, y, size, fpr, fpg, fpb)
    global gl

    %disp('Show FP')

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

end

% "Draw" functions get called on every screen refresh so long as the
% corresponding graphical object is to be displayed.  These are the
% functions that actually draw the object to the screen (and they make a
% fliprequest so that the drawn object will actually appear).
function DrawFP()
    global gl

    %disp('Draw FP')

    Screen('FillRect', gl.windowPtr, gl.fp.rgb, gl.fp.drawrect);
    gl.fliprequest = 1;

end

function HideFP()
    global gl

    %disp('Hide FP')

    gl.fp.on = 0;
    gl.fliprequest = 1;

end

function DrawTargs()
    global gl

    Screen('FillRect', gl.windowPtr, gl.targ.rgb, gl.targ.drawrects(1,:));
    Screen('FillRect', gl.windowPtr, gl.targ.rgb, gl.targ.drawrects(2,:));
    gl.fliprequest = 1;

end

function sendX(gridX)
    global gl udpCom
    gl.gridX = gridX;

    % Handshaking (gridx received)
    pnet(udpCom.sock, 'write', 'macready>> >>');
    pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);

end

function sendY(gridY)
    global gl udpCom
    gl.gridY = -gridY;

    % Handshaking (gridy received)
    pnet(udpCom.sock, 'write', 'macready>> >>');
    pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);

end


function PrepareStim(nStixelGrid, dvaPerStixel, stimDur, Lcc, Mcc, Scc, rfX, rfY, RFcorrect)
    global gl udpCom

    gl.fp.x = 0;
    gl.fp.y = 0;
    gl.stim.nStixelGrid = nStixelGrid;
    gl.stim.dvaPerStixel = dvaPerStixel;
    gl.stim.stimDur = stimDur;

    % Create Stimulus Space
    windowsize = round(gl.pixperdeg * dvaPerStixel * nStixelGrid / 2);
    colidx = -windowsize:windowsize; % This ensures window size will always be symmetric about the middle pixel, and therefore have an odd # of pixels, but
    rowidx = -windowsize:windowsize;

    %Create box
    box = zeros(numel(rowidx),numel(colidx));

    % Define spatial profile
    gl.stim.pixperstix = round(gl.pixperdeg*dvaPerStixel);
    pixgridX = gl.gridX * gl.pixperdeg * dvaPerStixel;
    pixgridY = gl.gridY * gl.pixperdeg * dvaPerStixel;
    gridpixX = [];
    gridpixY = [];
    for n = 1:numel(pixgridX)
        tempY = round(pixgridY(n)-gl.stim.pixperstix/2) : round(pixgridY(n)+gl.stim.pixperstix/2);
        tempX = round(pixgridX(n)-gl.stim.pixperstix/2) : round(pixgridX(n)+gl.stim.pixperstix/2);
        [temppixX, temppixY] = meshgrid(tempX,tempY);
        gridpixY = cat(1,gridpixY,temppixY(:));
        gridpixX = cat(1,gridpixX,temppixX(:));
    end
    stimgridX = gridpixX + colidx(end) + 1;
    stimgridY = gridpixY + rowidx(end) + 1;
    boxIdx = sub2ind(size(box),stimgridY,stimgridX);
    box(boxIdx) = 1;
    stimspatialprofile =  box(:,1:2:end);

    % Defining a drawing rectangle where the stimulus will be shown
    % [left, top, right, bottom]
    rfcenterpixel_x = round(rfX * gl.pixperdeg); % in pixels
    rfcenterpixel_y = round(rfY * gl.pixperdeg);
    boxleftedge = colidx(1); % in pixels
    boxrightedge = colidx(end);
    boxupperedge = rowidx(1);
    boxloweredge = rowidx(end);
    gl.stim.drawrect = ([...
        (gl.screenCenterXpix + rfcenterpixel_x + boxleftedge - 1) ...
        (gl.screenCenterYpix - rfcenterpixel_y + boxupperedge - 1) ...
        (gl.screenCenterXpix + rfcenterpixel_x + boxrightedge) ...
        (gl.screenCenterYpix - rfcenterpixel_y + boxloweredge) ]);

    if ~rem(gl.stim.drawrect(1), 2)
        gl.stim.drawrect(3) = gl.stim.drawrect(3)+1;
    else
        gl.stim.drawrect(1) = gl.stim.drawrect(1)-1;
    end

    % Creating an NxMx3 image to be displayed.
    im = zeros(size(stimspatialprofile,1), size(stimspatialprofile,2), 3); % preallocating space
    bkgndlms = gl.cal.M * gl.bkgndrgb;
    stimulusrgb = gl.cal.invM*(bkgndlms'.* (1 + [Lcc Mcc Scc]))'; % computing rgbs at the peak

    % Temporal parameters
    gl.stim.nframes = round(gl.framerate * stimDur);

    % Find rgb values for each pixel
    for gun = 1:3

        % rescaling stimspatialprofile to go between bkgndrgb and stimulusrgb (instead of between 0 and 1)
        tmp = (stimulusrgb(gun)-gl.bkgndrgb(gun))*(stimspatialprofile)+gl.bkgndrgb(gun);

        % Rounding the intensity values at each pixel so that we can use
        % these values as indices into the inverse gamma table
        tmp = round(tmp*size(gl.cal.invgammaTable,1))+1;
        if max(max(tmp)) > size(gl.cal.invgammaTable,1)
            disp('Out of gammut!')
            sca; keyboard;
        end

        % converting intensities to voltages
        tmp = gl.cal.invgammaTable(tmp(:), gun);
        im(:,:,gun) = reshape(tmp, size(stimspatialprofile,1), size(stimspatialprofile,2));

        if RFcorrect == 0

            im(:,:,gun) = fliplr(im(:,:,gun));

        end

    end

    % Make movie
    gl.stim.textureidx = Screen('MakeTexture', gl.windowPtr, TranslateToColourModeMex(im, [], gl.ccmode), [], [], 2);

    gl.stim.framecounter = 0;

    pnet(udpCom.sock, 'write', 'macready>> >>');
    pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);

end


%%
% Taking random numbers, putting them through the combined gaussgamma
% function and creating an image.  The order of elements is:
% red(1,1) red(2,1)... red(n,1)... red(n,n) green(1,1)...
function DrawStim()
global gl udpCom


    i = gl.stim.framecounter;
    Screen('DrawTexture',gl.windowPtr,gl.stim.textureidx,[],gl.stim.drawrect,[],0);

    if i < gl.stim.nframes
        gl.fliprequest = 1;
    else
        gl.stim.on = 0;
        Screen('Close',gl.stim.textureidx)
        gl.stim.textureidx = [];

        % Send Macready cue to plexon
        pnet(udpCom.sock, 'write', 'macready>> >>');
        pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
        disp('send1')
    end

end


function displayStim()
    global gl
    gl.stim.on = 1;
end


function ShowTargs(x, y, targsize)
global gl

xx = x/10*gl.pixperdeg;
yy = y/10*gl.pixperdeg;
fpsizeinpix = round(gl.pixperdeg*targsize/10);

gl.targ.drawrects(1,:) = [
    gl.screenCenterXpix+xx-floor(fpsizeinpix/2)+1, ...
    gl.screenCenterYpix-yy-floor(fpsizeinpix/2)+1, ...
    gl.screenCenterXpix+xx+ceil(fpsizeinpix/2), ...
    gl.screenCenterYpix-yy+ceil(fpsizeinpix/2)];
gl.targ.drawrects(2,:) = gl.targ.drawrects(1,:) + [-2*xx -0 -2*xx 0];
gl.targ.on = true;
end


function AllOff()
global gl

    gl.stim.on = 0;
    gl.targ.on = 0;
    gl.fp.on = 0;
    gl.fliprequest = 1;

    if ~isempty(gl.stim.textureidx)
        textureidx = unique(gl.stim.textureidx);
        for n = 1:numel(textureidx)
            Screen('Close',textureidx)
        end
    end

    gl.stim.textureidx = [];

end

function DealWithMessage(msgSize)
    global udpCom gl
    message = pnet(udpCom.sock, 'read', msgSize, 'char');

    if strncmp(message, 'return', 6)
        a = dbstack;  % Check whether called from another function or from command line
        if ~strcmp(a(end).name, mfilename)
            gl.timetoleave = 1;
        end
    end

    try
        eval(message);
        disp(message);
    catch ME
        fprintf('Ignoring uninterpretable message: "%s"\n', message);
        disp(getReport(ME));
        sca; keyboard;
    end
end
