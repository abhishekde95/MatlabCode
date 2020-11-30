%#ok<*DEFNU>
function GridLMSubunit_Slave
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
global gl udpCom

Screen('Flip', gl.windowPtr, 0, 0);
if gl.stim.on == 1 && gl.stim.framecounter == 0
    % Handshaking
    pnet(udpCom.sock, 'write', 'macdone>> >>');
    pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
end
if gl.stim.on == 1
    gl.stim.framecounter = gl.stim.framecounter + 1;
end
gl.fliprequest = 0;
end

function InitDisplay(mondist, screenwidth, varargin) % varargin{1} and {2} are defined to be the calibration and fundamentals
global gl
[nil,calfilename] = fileparts(varargin{1});
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
    [nil,fundfilename] = fileparts(varargin{2});
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

HideCursor();
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

function sendX(gridX)
global gl udpCom
gl.gridX = gridX;

% Handshaking
pnet(udpCom.sock, 'write', 'macdone>> >>');
pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);

end

function sendY(gridY)
global gl udpCom
gl.gridY = -gridY;

% Handshaking
pnet(udpCom.sock, 'write', 'macdone>> >>');
pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);

end


function PrepareStim(nStixelGrid, dvaPerStixel, stimDur, cc1, cc2, cc3, seed, epoch)
global gl udpCom

gl.stim.nStixelGrid = nStixelGrid;
gl.stim.dvaPerStixel = dvaPerStixel;
gl.stim.stimDur = stimDur;
gl.stim.epoch = epoch;
gl.stim.seed = seed;
gl.stim.pixperstix = round(gl.pixperdeg*gl.stim.dvaPerStixel);

if gl.stim.epoch == 1 || gl.stim.epoch == 2

    % Create Stimulus Space
    windowsize = round(gl.pixperdeg * dvaPerStixel * nStixelGrid / 2);
    colidx = -windowsize:windowsize; % This ensures window size will always be symmetric about the middle pixel, and therefore have an odd # of pixels, but
    rowidx = -windowsize:windowsize;

    %Create box
    box = zeros(numel(rowidx),numel(colidx));

    % Define spatial profile
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
    rfcenterpixel_x = round(gl.bar.xy(1) * gl.pixperdeg); % in pixels
    rfcenterpixel_y = round(gl.bar.xy(2) * gl.pixperdeg);
    boxleftedge = colidx(1); % in pixels
    boxrightedge = colidx(end);
    boxupperedge = rowidx(1);
    boxloweredge = rowidx(end);
    gl.stim.drawrect = ([...
        (gl.screenCenterXpix + rfcenterpixel_x + boxleftedge-1) ...
        (gl.screenCenterYpix - rfcenterpixel_y + boxupperedge-1) ...
        (gl.screenCenterXpix + rfcenterpixel_x + boxrightedge) ...
        (gl.screenCenterYpix - rfcenterpixel_y + boxloweredge) ]);

    if ~rem(gl.stim.drawrect(1), 2)
        gl.stim.drawrect(3) = gl.stim.drawrect(3)+1;
    else
        gl.stim.drawrect(1) = gl.stim.drawrect(1)-1;
    end

elseif gl.stim.epoch == 3

    gl.stim.pixperstix = round(gl.pixperdeg*dvaPerStixel);
    pixgridX = gl.gridX * gl.pixperdeg * dvaPerStixel;
    pixgridY = gl.gridY * gl.pixperdeg * dvaPerStixel;
    gridpixX = [];
    gridpixY = [];
    for n = 1:numel(pixgridX) %converting from grid #'s to pixel #'s
        tempY = round(pixgridY(n)-gl.stim.pixperstix/2) : round(pixgridY(n)+gl.stim.pixperstix/2);
        tempX = round(pixgridX(n)-gl.stim.pixperstix/2) : round(pixgridX(n)+gl.stim.pixperstix/2);
        [temppixX, temppixY] = meshgrid(tempX,tempY);
        gridpixY = cat(1,gridpixY,temppixY(:));
        gridpixX = cat(1,gridpixX,temppixX(:));
    end
    coeff = princomp([gridpixX gridpixY]);
    PC_1 = coeff(:,1);
    PC_end = [PC_1(2); -PC_1(1)];
    theta = mod(atan2(PC_1(2),-PC_1(1)),pi); % Preferred direction. 0 means horizonal
    proj = [gridpixX gridpixY] * PC_end;
    halfperiod = range(proj); % in pixels
    gl.stim.sf = 1/(halfperiod * 2 / gl.pixperdeg); % in cycles/deg
    gl.stim.theta = theta;
    gl.stim.sigma = .4;
    gl.stim.nframesplateau = 26;
    gl.stim.nframesramp = 12;
    gl.stim.phi = 0;
    gl.stim.gamma = 1;
    gl.stim.tf = 3;
    gl.stim.nsigmas = 3;
    lambda = 1/gl.stim.sf;

    ramp = linspace(0,1,gl.stim.nframesramp);
    plateau = ones(1, gl.stim.nframesplateau);
    temporalprofile = [ramp plateau fliplr(ramp)];
    nframes = numel(temporalprofile);
    stimsizeindeg = gl.stim.sigma*gl.stim.nsigmas; % This is half stim size (we we go out +/- nsigmas)
    stimsizeinpix = round(stimsizeindeg*gl.pixperdeg);  % full stim size in doublewide pixels
    [x,y] = meshgrid(stimsizeindeg*linspace(-1,1,stimsizeinpix), stimsizeindeg*linspace(-1,1,2*stimsizeinpix));
    % x and y are in dva

    xprime = x*cos(-gl.stim.theta) + y*sin(-gl.stim.theta);
    yprime = -x*sin(-gl.stim.theta) + y*cos(-gl.stim.theta);
    deltaphase = gl.stim.tf * 2 * pi / gl.framerate;
    gabor = zeros(2*stimsizeinpix,stimsizeinpix,nframes);
    for a = 1:length(temporalprofile)
        phase = gl.stim.phi+(a-1)*deltaphase;
        gabor(:,:,a) = temporalprofile(a)*exp(-(xprime.^2 + gl.stim.gamma.^2 .* yprime.^2)./ (2*gl.stim.sigma^2)).*cos(2*pi*yprime./lambda+phase);
    end
    gl.stim.template = gabor;
    gl.stim.nframes = nframes;

    % Creating the drawing window
    x = (gl.bar.xy(1) * gl.pixperdeg) + (mean(gl.gridX) * gl.stim.pixperstix); % in pixels
    y = (gl.bar.xy(2) * gl.pixperdeg) - (mean(gl.gridY) * gl.stim.pixperstix);
    stimsizeinpix = size(gl.stim.template,2); % in doublewide pix
    gl.drawrect = round([gl.screenCenterXpix+x-stimsizeinpix ...
        gl.screenCenterYpix-y-stimsizeinpix ...
        gl.screenCenterXpix+x+stimsizeinpix...
        gl.screenCenterYpix-y+stimsizeinpix]);

    if(rem(gl.drawrect(1), 2)) %if the rectangle starts on an odd pixel
        gl.drawrect(1) = gl.drawrect(1) - 1;
        gl.drawrect(3) = gl.drawrect(3) - 1;
    end

    if (gl.drawrect(3)-gl.drawrect(1) ~= 2*stimsizeinpix || gl.drawrect(4)-gl.drawrect(2) ~= 2*stimsizeinpix)
        sca;
        keyboard;
        % Draw window is the wrong size
    end

end

% Create the actual movie and store
if epoch == 1 %Dense Noise

    lumcc = cc1;
    colcc = cc2;
    Scc = cc3;

    LccLum = sqrt((lumcc^2)/2);
    MccLum = LccLum;
    LccCol = sqrt((colcc^2)/2);
    MccCol = LccCol;
    cc = [LccLum MccLum Scc; -LccLum -MccLum Scc; LccCol -MccCol Scc; -LccCol MccCol Scc];

    gl.stim.mu = [gl.cal.gammaTable(gl.bkgndRGB(1)+1,1),...
        gl.cal.gammaTable(gl.bkgndRGB(2)+1,2),...
        gl.cal.gammaTable(gl.bkgndRGB(3)+1,3)];
    mulms = gl.cal.M*gl.stim.mu';
    stimlms = cc.*repmat(mulms',size(cc,1),1) + repmat(mulms',size(cc,1),1);
    rgbmat = gl.cal.invM*stimlms';

    if any(rgbmat(:) > 1) || any(rgbmat(:) < 0)
        disp('The requested colors are out of gamut');
        sca; keyboard
    end

    % Converting to voltages using bsxfun to get linear indices
    idxs = bsxfun(@plus, round(65535*rgbmat')+1, 0:size(gl.cal.invgammaTable,1):numel(gl.cal.invgammaTable)-1);
    gl.lmsbinaryrgbmat = gl.cal.invgammaTable(idxs);

    % Preparing a template filled with stixel indices for putting the random
    % numbers into the image matrix.
    stixidxs = reshape(1:nStixelGrid^2,nStixelGrid,nStixelGrid);
    pixextents = ones(2*gl.stim.pixperstix,gl.stim.pixperstix); % takes into account the aspect ratio in colour mode
    gl.stim.textureidx = kron(stixidxs,pixextents);

elseif epoch == 2 %GLMP

    % Creating an NxMx3 image to be displayed.
    im = zeros(size(stimspatialprofile,1), size(stimspatialprofile,2), 3); % preallocating space
    bkgndlms = gl.cal.M * gl.bkgndrgb;
    L = cc1;
    M = cc2;
    S = cc3;
    stimulusrgb = gl.cal.invM*(bkgndlms'.* (1 + [L M S]))'; % computing rgbs at the peak

    % Temporal parameters
    gl.stim.nframes = round(gl.framerate * stimDur);

    % Find rgb values for each pixel
    for gun = 1:3

        % rescaling stimspatialprofile to go between bkgndrgb and stimulusrgb (instead of between 0 and 1)
        tmp = (stimulusrgb(gun)-gl.bkgndrgb(gun))*(stimspatialprofile)+gl.bkgndrgb(gun);

        % Rounding the intensity values at each pixel so that we can use
        % these values as indices into the inverse gamma table
        try
        tmp = round(tmp*size(gl.cal.invgammaTable,1))+1;
        catch
            keyboard
        end
        if max(max(tmp)) > size(gl.cal.invgammaTable,1)
            disp('Out of gammut!')
            sca; keyboard;
        end

        % converting intensities to voltages
        tmp = gl.cal.invgammaTable(tmp(:), gun);
        im(:,:,gun) = reshape(tmp, size(stimspatialprofile,1), size(stimspatialprofile,2));

    end

    % Make movie
    gl.stim.textureidx = Screen('MakeTexture', gl.windowPtr, TranslateToColourModeMex(im, [], gl.ccmode), [], [], 2);

elseif epoch == 3
    % This is for testing Gabors (like is NeuroThresh), and comparing
    % to GLMS responses

    L = cc1;
    M = cc2;
    S = cc3;
    stimconecontrast = [L M S];

    % Getting RGBs
    bkgndlms = gl.cal.M*gl.bkgndrgb;
    gl.stim.rgb = gl.cal.invM*(bkgndlms.*(1+stimconecontrast'));

end

gl.stim.framecounter = 0;
pnet(udpCom.sock, 'write', 'macdone>> >>');
pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);

end

%%
% Taking random numbers, putting them through the combined gaussgamma
% function and creating an image.  The order of elements is:
% red(1,1) red(2,1)... red(n,1)... red(n,n) green(1,1)...
function DrawStim()
global gl udpCom
if gl.stim.epoch == 1 %Dense Nosie
    imgsize = [2*gl.stim.nStixelGrid*round(gl.stim.pixperstix) gl.stim.nStixelGrid*round(gl.stim.pixperstix) 3];
    [randnums,gl.stim.seed] = getEJrandnums3(gl.stim.nStixelGrid^2,gl.stim.seed);
    randnums = mod(randnums, size(gl.lmsbinaryrgbmat,1))+1;
    rgbs = gl.lmsbinaryrgbmat(randnums,:);
    img = reshape(rgbs(gl.stim.textureidx,:), imgsize);
    %gl.stim.drawrect = [100 100 220 220]; % take out
    tex = Screen('MakeTexture', gl.windowPtr, img, [], [], 2);
    Screen('DrawTexture', gl.windowPtr, tex, [], gl.stim.drawrect, [], 0);
    Screen('Close', tex);
    gl.fliprequest = 1;

elseif gl.stim.epoch == 2 %GLMP

    i = gl.stim.framecounter;
    Screen('DrawTexture',gl.windowPtr,gl.stim.textureidx,[],gl.stim.drawrect,[],0);

    if i < gl.stim.nframes
        gl.fliprequest = 1;
    else
        gl.stim.on = 0;
        Screen('Close',gl.stim.textureidx)
        gl.stim.textureidx = [];
        pnet(udpCom.sock, 'write', 'macdone>> >>');
        pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
    end

elseif gl.stim.epoch == 3 % Large Gabor
    NGAMMASTEPS = size(gl.cal.invgammaTable,1);
    im = zeros(size(gl.stim.template,1), size(gl.stim.template,2), 3);
    frame = gl.stim.template(:,:,gl.stim.framecounter+1);
    for plane = 1:3
        tmp = frame.*(gl.stim.rgb(plane)-gl.bkgndrgb(plane)) + gl.bkgndrgb(plane);
        tmp = round(tmp*(NGAMMASTEPS-1))+1;
        try
            tmp = gl.cal.invgammaTable(tmp, plane);
        catch
            sca
            keyboard
        end
        im(:,:,plane) = reshape(tmp, size(im(:,:,1)));
    end
    tex=Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(im, [], gl.ccmode), [], [], 2);
    Screen('DrawTexture',gl.windowPtr,tex, [], gl.drawrect, [], 0);
    Screen('Close',tex);

    if gl.stim.framecounter+1 < gl.stim.nframes
        gl.fliprequest = 1;
    else
        gl.stim.on = 0;
        %Screen('Close',gl.stim.textureidx)
        %gl.stim.textureidx = [];
        pnet(udpCom.sock, 'write', 'macdone>> >>');
        pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
        disp('send done message')
    end

else
    disp('Problem with Epochs...')
    sca; keyboard;

end

end


function displayStim()
global gl
gl.stim.on = 1;
end



function AllOff()
global gl

gl.stim.on = 0;
gl.fp.on = 0;
gl.fliprequest = 1;

% if ~isempty(gl.stim.textureidx)
%     textureidx = unique(gl.stim.textureidx);
%     for n = 1:numel(textureidx)
%         Screen('Close',textureidx)
%     end
% end

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
