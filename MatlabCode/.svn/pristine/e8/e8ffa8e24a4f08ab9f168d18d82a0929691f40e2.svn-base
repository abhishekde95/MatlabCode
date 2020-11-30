%% Code created by GDLH as a stimulus paradigm 9/20/11
% Edited to make a bar stimulus by JPW Sept_2011

clear all
close all

%% User-Defined Variables

stim.cc = [0.2 0.2 0];
stim.gaborSigma = 10; %Std in 1/10ths degree VA
stim.nstd = 4;
stim.gaborGamma = 1; %Ratio between x and y std
stim.gaborLambda = stim.gaborSigma/10*stim.nstd; %Spatial period (DVA/cycle)
stim.gaborTheta = pi/4; %Orientation in radians
stim.flashDeltaPhase = 0;
stim.ncycles = 12;
stim.driftRate = 3; %cycles/s


%% Auto-Defined Variables

gaborPhase = 0;
USECOLOURMODE = 0
MONDIST = 100; % cm
STIMSIZEINDEG = 4;
gl.screenWidthcm = 40; % cm
stimpos = [-3 0]; % in degrees of visual angle relative to the center of the screen

load ('T_cones_smj10.mat');
fundamentals = T_cones_smj10;
%load('/MatlabCode/Monitor Calibration/Monitor data/Dell 4/Dell4BitsCal');
load('rig1mon')
% cal = cals{end};
% P_device = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls([380 5 81]));
%M = fundamentals*P_device;
M = fundamentals*monspd;

gl.bkgndRGB = round(255*cal.bgColor)';
gl.bkgndRGB = round(255*cal.bkgndrgb)';
gl.cal.gammaTable = cal.gammaTable;
gl.cal.invgammaTable = InvertGamma(gl.cal.gammaTable, USECOLOURMODE);
gl.cal.monSpd = cal.P_device;
gl.mondistcm = MONDIST;

bkgndrgb = [gl.cal.gammaTable(gl.bkgndRGB(1)+1,1); ...
    gl.cal.gammaTable(gl.bkgndRGB(2)+1,2); ...
    gl.cal.gammaTable(gl.bkgndRGB(3)+1,3)];
bkgndlms = M*bkgndrgb;

try
    % ---- This ia all boilerplate stuff that sets up for displaying ---
    % Note: "SCREEN" is a multi-purpose function that is at the heart of
    % the Psychophysics Toolbox. Type "help screen" at the command prompt for
    % more information.
    gl.windowPtr = Screen('OpenWindow',0, 255*cal.bgColor);
    [screenwidthpix, screenheightpix]  = Screen('WindowSize', gl.windowPtr);
    gl.screenCenterXpix = screenwidthpix/2;
    gl.screenCenterYpix = screenheightpix/2;
    gl.framerate = Screen('FrameRate', gl.windowPtr, 1);
    pixpercm = screenwidthpix/gl.screenWidthcm;
    cmperdeg = gl.screenWidthcm/(2*atan2(gl.screenWidthcm/2, gl.mondistcm)*180/pi);
    gl.pixperdeg = pixpercm*cmperdeg;  % Number os pixels per degree of visual angle

    stim.sizeinpix = round(gl.pixperdeg*(stim.gaborSigma/10 * stim.nstd));
    stim.flashDeltaPhase = stim.driftRate * 2 * pi * (1/gl.framerate);
    stim.nframes = stim.ncycles * (1/stim.driftRate) * gl.framerate;

    if isinteger(stim.sizeinpix)
        halfSize = stim.sizeinpix/2;
        row = -halfSize:halfSize-1; %subtract one so that you don't overlap the texture window;
        col = -halfSize:halfSize-1;
    else
        halfsize = (stim.sizeinpix-1)/2;
        row = -halfsize:halfsize;
        col = -halfsize:halfsize;
    end

    [X, Y] = meshgrid(row, col);
    xprime = X .* cos(-stim.gaborTheta) + Y .* sin(-stim.gaborTheta);
    yprime = -X .* sin(-stim.gaborTheta) + Y .* cos(-stim.gaborTheta);

    gaussian = exp(-(xprime.^2 + stim.gaborGamma.^2 .* yprime.^2)./ (2.* (gl.pixperdeg*stim.gaborSigma/10).^2));

    % Drawing loop
    for i=1:stim.nframes

        gabor = gaussian .* cos(2 .* pi .* yprime ./ (gl.pixperdeg *stim.gaborLambda) + gaborPhase);
        stimspatialprofile = gabor;

        % Loading a linear normalized gamma table
        % We don't want the graphics card to deal with gamma compensation, we
        % have to do it ourselves if we want to take advantage of the Bits++
        clut = repmat(linspace(0,1,256),3,1)';
        Screen('LoadNormalizedGammaTable', 0, clut);

        % Defining a drawing rectangle where the stimulus will be shown
        % [left, bottom, right, top]
        drawrect = round([gl.screenCenterXpix+(stimpos(1)*gl.pixperdeg)-stim.sizeinpix/2 ...
            gl.screenCenterYpix-(stimpos(2)*gl.pixperdeg)-stim.sizeinpix/2 ...
            gl.screenCenterXpix+(stimpos(1)*gl.pixperdeg)+stim.sizeinpix/2 ...
            gl.screenCenterYpix-(stimpos(2)*gl.pixperdeg)+stim.sizeinpix/2]);
        % Shifting the drawing rectangle if it happens to start on an
        % odd-numbered pixel. This is only important when using the Bits++
        if(rem(drawrect(1), 2))
            drawrect(1) = drawrect(1) - 1;
            drawrect(3) = drawrect(3) - 1;
        end

        % Creating an NxMx3 image to be displayed.
        im = zeros(stim.sizeinpix, stim.sizeinpix, 3); % preallocating space
        stimulusrgb = inv(M)*(bkgndlms'.*(1+stim.cc))'; % computing rgbs at the peak
        for gun = 1:3
            % rescaling stimspatialprofile to go between bkgndrgb and stimulusrgb (instead of between 0 and 1)
            tmp = (stimulusrgb(gun)-bkgndrgb(gun))*stimspatialprofile+bkgndrgb(gun);
            % Rounding the intensity values at each pixel so that we can use
            % these values as indices into the inverse gamma table
            tmp = round(tmp*size(gl.cal.invgammaTable,1))+1;
            % Checking for out of gamut errors
            if (any(tmp(:) > size(gl.cal.invgammaTable,1)) | any(tmp(:) < 1))
                [gun min(tmp(:)) max(tmp(:))]
                error('Out of gamut error');
            end
            % converting intensities to voltages
            tmp = gl.cal.invgammaTable(tmp, gun);
            tmp = round(tmp*(size(gl.cal.invgammaTable,1)-1));
            im(:,:,gun) = reshape(tmp, stim.sizeinpix, stim.sizeinpix);
        end
        stimtex=Screen('MakeTexture', gl.windowPtr, im);

        gaborPhase = gaborPhase + stim.flashDeltaPhase;
        Screen('DrawTexture',gl.windowPtr,stimtex,[],drawrect,[],0);
        Screen('Flip',gl.windowPtr);
    end
    Screen('Close',stimtex);
    Screen('CloseAll');
catch
    sca;
    rethrow(lasterror);
end
