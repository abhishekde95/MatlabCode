function cal = CalibrateAmbDrvr(cal,USERPROMPT,whichMeterType,blankOtherScreen)
% cal =  CalibrateAmbDrvr(cal,USERPROMPT,whichMeterType,blankOtherScreen)
%
% This script does the work for monitor ambient calibration.

% 4/4/94		dhb		Wrote it.
% 8/5/94		dhb, ccc	More flexible interface.
% 9/4/94		dhb		Small changes.
% 10/20/94	dhb		Add bgColor variable.
% 12/9/94   ccc   Nine-bit modification
% 1/23/95		dhb		Pulled out working code to be called from elsewhere.
%						dhb		Make user prompting optional.
% 1/24/95		dhb		Get filename right.
% 12/17/96  dhb, jmk  Remove big bug.  Ambient wasn't getting set.
% 4/12/97   dhb   Update for new toolbox.
% 8/21/97		dhb		Don't save files here.
%									Always measure.
% 4/7/99    dhb   NINEBIT -> NBITS
%           dhb   Handle noMeterAvail, RADIUS switches.
% 9/22/99   dhb, mdr  Make boxRect depend on boxSize, defined up one level.
% 12/2/99   dhb   Put background on after white box for aiming.
% 8/14/00   dhb   Call to CMETER('Frequency') only for OS9.
% 8/20/00   dhb   Remove bits arg to SetColor.
% 8/21/00   dhb   Remove RADIUS arg to MeasMonSpd.
% 9/11/00   dhb   Remove syncMode code, any direct refs to CMETER.
% 9/14/00   dhb   Use OpenWindow to open.
%           dhb   Made it a function.
% 7/9/02    dhb   Get rid of OpenWindow, CloseWindow.
% 9/23/02   dhb, jmh  Force background to zero when measurements come on.
% 2/26/03   dhb   Tidy comments.
% 4/1/03    dhb   Fix ambient averaging.
% 8/19/12   dhb   Add codelet suggested by David Jones to clean up at end.  See comment in CalibrateMonSpd.
% 8/19/12   mk    Rewrite setup and clut code to be able to better cope with all
%                 the broken operating systems / drivers / gpus and to also
%                 support DataPixx/ViewPixx devices.

global g_usebitspp g_usecolorpp

% If the global flag for using Bits++ is empty, then it hasn't been
% initialized and default it to 0.
if isempty(g_usebitspp)
    g_usebitspp = 0;
end

% Ditto for the Color++/C48 flag
if isempty(g_usecolorpp)
    g_usecolorpp = 0;
end

if nargin == 0
    return
elseif nargin < 2 || isempty(USERPROMPT)
    USERPROMPT = 1;
elseif nargin < 3 || isempty(whichMeterType)
    CMCheckInit();
elseif nargin < 4 || isempty(blankOtherScreen)
    blankOtherScreen = 0;
end

% User prompt
if USERPROMPT
    if cal.describe.whichScreen == 0
        fprintf('Hit any key to proceed past this message and display a box.\n');
        fprintf('Focus radiometer on the displayed box.\n');
        fprintf('Once meter is set up, hit any key - you will get %g seconds\n',...
            cal.describe.leaveRoomTime);
        fprintf('to leave room.\n');
        KbStrokeWait(-1);
    else
        fprintf('Focus radiometer on the displayed box.\n');
        fprintf('Once meter is set up, hit any key - you will get %g seconds\n',...
            cal.describe.leaveRoomTime);
        fprintf('to leave room.\n');
    end
end

% Blank other screen, if requested:
if blankOtherScreen
    % We simply open an onscreen window with black background color:
    Screen('OpenWindow', cal.describe.whichBlankScreen, 0);
end

% Setup screen to be measured
% ---------------------------

% If g_usecolorspp == 0, then
% prepare imaging pipeline for Bits+ Bits++ CLUT mode, or DataPixx/ViewPixx
% L48 CLUT mode (which is pretty much the same). If such a special output
% device is used, the Screen('LoadNormalizedGammatable', win, clut, 2);
% command uploads 'clut's into the device at next Screen('Flip'), taking
% care of possible graphics driver bugs and other quirks:

% If g_usecolorpp > 0, then we leave the CLUT alone, and let the Color++/C48
% shaders directly paint a central box with high-precision color.
PsychImaging('PrepareConfiguration');

if g_usebitspp == 1
    if g_usecolorpp
        PsychImaging('AddTask', 'General', 'EnableBits++Color++Output', 1);
    else
        % Setup for Bits++ CLUT mode. This will automatically load proper
        % identity gamma tables into the graphics hardware and into the Bits+:
        PsychImaging('AddTask', 'General', 'EnableBits++Bits++Output');
    end
end

if g_usebitspp == 2
    if g_usecolorpp
        PsychImaging('AddTask', 'General', 'EnableDataPixxC48Output', 1);
    else
        % Setup for DataPixx/ViewPixx etc. L48 CLUT mode. This will
        % automatically load proper identity gamma tables into the graphics
        % hardware and into the device:
        PsychImaging('AddTask', 'General', 'EnableDataPixxL48Output');
    end
end

% Open the window:
[window, screenRect] = PsychImaging('OpenWindow', cal.describe.whichScreen, cal.bgColor(:)');
if cal.describe.whichScreen == 0
    HideCursor();
end

if ~g_usecolorpp
    % Zero the color lookup table and set the the surround color in the first position
    theClut = zeros(256,3);
    theClut(1,:) = cal.bgColor(:)';
    % The second position of the clut is the color of aiming square
    theClut(2,:) = [1 1 1];
    Screen('LoadNormalizedGammaTable', window, theClut, 2*(g_usebitspp > 0));
    Screen('FillRect', window, 0); % 0th offset -> LUT(1,:)
end

% Draw a box in the center of the screen:
if ~isfield(cal.describe, 'boxRect')
    boxRect = [0 0 cal.describe.boxSize cal.describe.boxSize];
    boxRect = CenterRect(boxRect,screenRect);
else
    boxRect = cal.describe.boxRect;
end

if g_usecolorpp
    Screen('FillRect', window, [1 1 1], boxRect);
else
    Screen('FillRect', window, 1, boxRect); % 1st offset -> LUT(2,:)
end
Screen('Flip', window);

% Put up a white square for focusing and wait for user
if USERPROMPT == 1
    KbStrokeWait(-1);
    fprintf('Pausing for %d seconds ...', cal.describe.leaveRoomTime);
    WaitSecs(cal.describe.leaveRoomTime);
    fprintf(' done\n');
end

% Start timing
tic;

ambient = zeros(cal.describe.S(3), 1);
for a = 1:cal.describe.nAverage
    if ~g_usecolorpp
        Screen('FillRect', window, 1, boxRect);
        Screen('Flip', window, 0, 1);
    end
    
    % Measure ambient
    if g_usecolorpp
        ambient = ambient + MeasMonSpd(window, [0 0 0]', cal.describe.S, 0, whichMeterType, boxRect);
    else
        ambient = ambient + MeasMonSpd(window, [0 0 0]', cal.describe.S, 0, whichMeterType, theClut);
    end
end
ambient = ambient / cal.describe.nAverage;

% Close the screen, restore cluts:
if ~g_usecolorpp && g_usebitspp
    % Load identity clut on Bits++ / DataPixx et al.:
    BitsPlusPlus('LoadIdentityClut', window);
    Screen('Flip', window, 0, 1);
elseif ~g_usebitspp
    LoadIdentityClut(window);
end

% Show hidden cursor:
if cal.describe.whichScreen == 0
    ShowCursor();
end

% Close all windows:
Screen('CloseAll');

% Report time:
fprintf('CalibrateAmbDrvr measurements took %g minutes\n', toc/60);

% Update structure
Smon = cal.describe.S;
Tmon = WlsToT(Smon);
cal.P_ambient = ambient;
cal.T_ambient = Tmon;
cal.S_ambient = Smon;
