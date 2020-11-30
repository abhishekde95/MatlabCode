function number_of_frames = Jacob480Demo(frames, rotfreq, offsetStepSize, ...
    maxSteps, bkgnd)
% Based on DatapixxShowGabor.m
%
% Draws a Benham's top stimulus at 480 Hz.
% Propixxx controler should be in 120 Hz, 1920 x 1080 pixels mode.
% Set this up by running "vputils", selecting the PP controler (option 3),
% and changing the "edid" to options 3 and 6 (both of which run at 120 Hz).
% Then, go to System Preferences -> Displays and make sure that the display
% is set to 120 Hz and 1920x1080. The 480 frame rate is achieved by
% painting four images per 1920x1080 screen and using
% Datapixx('SetPropixxDlipSequenceProgram',2) to scan across them once
% every 120 Hz frame.
%
% Blue rotation offset is specified in radians. Negative values means blue is
% lagging behind red and green. Positve values means blue precedes red and green.
% For some reason, the images plotted in the figure are inaccurate.
% Stimulus rotates clockwise.
dark = 0;
keyboardSamplingRate = 4;

Datapixx('Open');
Datapixx('SetPropixxDlpSequenceProgram', 2); % 2 for 480, 5 for 1440 Hz, 0 for normal
Datapixx('RegWrRd');
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
screenid = max(Screen('Screens'));
win = PsychImaging('OpenWindow', screenid, 255*bkgnd);
[screenwidth, screenheight]=Screen('WindowSize', win);

if nargin < 1 || isempty(frames)
    frames = 4800;
end
if nargin < 2 || isempty(rotfreq)
    rotfreq = 10; % Hz
end
rotinc = 360/(480/rotfreq); % rotation increment in degrees
if nargin < 3 || isempty(rotfreq)
    bluerotoffset = 0; % rad
end

% Get refresh rate
ifi = Screen('GetFlipInterval', win);

% Making sure ProPixx is if the correct configuration
if (abs(ifi-1/120) > 10^-5 || screenheight ~= 1080 || screenwidth ~= 1920)
    Screen('CloseAll');
    error('Screen set up incorrectly. Need to be in 120 Hz mode with screen size 1920 x 1080.')
end

% Making the texture
tic
% Plotting the arcs
% nrings = [4 3 3];
% ring_theta = [0 pi/3; pi/3 2*pi/3; 2*pi/3 pi];
% start_ring_r = .1; % innermost ring
% ringthickness = .04;
% ringmargin = .04;
% m = zeros(256,256,3);
% [tmpx, tmpy] = meshgrid(1:size(m,1),1:size(m,2));
% tmpr = sqrt((tmpx(:)-(size(m,1)+1)/2).^2+(tmpy(:)-(size(m,2)+1)/2).^2);
% r = tmpr./(size(m,1)/2);
% theta = atan2(tmpy(:)-(size(m,2)+1)/2, tmpx(:)-(size(m,1)+1)/2);
% addedtheta = [0 0 bluerotoffset];
% for gun = 1:3
%     ring_r = start_ring_r;
%     % Generate a hemidisk
%     tmpm = bkgnd*ones(size(m,1),size(m,2));
%     tmpm(mod(theta+addedtheta(gun),2*pi) > pi & r <= 1) = dark; % lightness of solid hemidisk
%     for i = 1:length(nrings) % Looping over the three wedges
%         for j = 1:nrings(i) % looping over the individual arcs in each wedge
%             ring_r = ring_r + ringthickness+ringmargin;
%             L = r > ring_r-ringthickness/2 &...
%                 r < ring_r+ringthickness/2 &...
%                 is_angle_between(theta+addedtheta(gun), ring_theta(i,1), ring_theta(i,2));
%             tmpm(L) = dark;
%         end
%     end
%     m(:,:,gun) = tmpm;
% end
tic
allowableBlueOffsets = (-maxSteps:maxSteps);
precomputedFrames = arrayfun( ...
    @(x) GenerateFrame(x, bkgnd, dark), allowableBlueOffsets * offsetStepSize, ...
    'UniformOutput', false);
offsetValues = arrayfun(@(x) x, allowableBlueOffsets, ...
    'UniformOutput', false);
disp(offsetValues);
frameLookup = containers.Map(offsetValues, precomputedFrames);
disp(['Precomputing frames took: ' num2str(toc) ' seconds.']);
currOffset = 0;

% m = GenerateFrame(bluerotoffset, bkgnd, dark);
figure;
forPlotting = frameLookup(currOffset);
for i = 1:3
    subplot(3,1,i);
    imagesc(forPlotting(:,:,i));
    axis square;
    axis ij;
    colormap(gray);
end

%Transform the gabors into texture
gabortex=Screen('MakeTexture', win, frameLookup(currOffset), [], [], 2);
texrect = Screen('Rect', gabortex);
dstRect = CenterRectOnPoint(texrect, screenwidth/4, screenheight/4);
rotAngle = 0;

% Initial flip
vbl = Screen('Flip', win);
number_of_frames = 1;

% Stimulus
keyboardCheckInterval = 1 / keyboardSamplingRate;
tic;
while number_of_frames < frames
    if toc > keyboardCheckInterval
        [isKey, ~, keyCode, ~] = KbCheck();
        if isKey
            keyName = KbName(keyCode);
            if strcmp(keyName, 'LeftArrow') || strcmp(keyName, 'RightArrow')
                if strcmp(keyName, 'LeftArrow')
                    currOffset = max(currOffset - 1, -maxSteps);
                else
                    currOffset = min(currOffset + 1, maxSteps);
                end
                disp(currOffset)
                tic;
                
                gabortex=Screen( ...
                    'MakeTexture', win, frameLookup(currOffset), [], [], 2);
                texrect = Screen('Rect', gabortex);
                dstRect = CenterRectOnPoint(texrect, screenwidth/4, screenheight/4);
            end
        end
    end
    % FRAME 1
    if number_of_frames > frames
        vbl = Screen('Flip', win, vbl + 0.5 * ifi);
        break;
    end
    Screen('DrawTexture', win, gabortex, [], dstRect+[0 0 0 0], rotAngle, [], [], [], [], []);
    number_of_frames = number_of_frames + 1;
    
    % Animate
    rotAngle = mod(rotAngle + rotinc,360);
    
    % FRAME 2
    if number_of_frames > frames
        vbl = Screen('Flip', win, vbl + 0.5 * ifi);
        break;
    end
    Screen('DrawTexture', win, gabortex, [], dstRect+[screenwidth/2 0 screenwidth/2 0], rotAngle, [], [], [], [], []);
    number_of_frames = number_of_frames + 1;
    
    % Animate
    rotAngle = mod(rotAngle + rotinc,360);
    
    % FRAME 3
    if number_of_frames > frames
        vbl = Screen('Flip', win, vbl + 0.5 * ifi);
        break;
    end
    Screen('DrawTexture', win, gabortex, [], dstRect+[0 screenheight/2 0 screenheight/2], rotAngle, [], [], [], [], []);
    number_of_frames = number_of_frames + 1;
    
    % Animate
    rotAngle = mod(rotAngle + rotinc,360);
    
    % FRAME 4
    if number_of_frames > frames
        vbl = Screen('Flip', win, vbl + 0.5 * ifi);
        break;
    end
    
    Screen('DrawTexture', win, gabortex, [], dstRect+[screenwidth/2 screenheight/2 screenwidth/2 screenheight/2], rotAngle, [], [], [], [], []);
    number_of_frames = number_of_frames + 1;
    
    % Animate
    rotAngle = mod(rotAngle + rotinc,360);
    
    % Update the four Qudrants
    vbl = Screen('Flip', win, vbl + 0.5 * ifi);
    
end

% Blank screen
Screen('Flip', win);

% Exit on escape
KbName('UnifyKeyNames');
escapeKey = KbName('ESCAPE');
while 1
    % Check the state of the keyboard.
    [ keyIsDown, ~, keyCode ] = KbCheck;
    if keyIsDown
        if keyCode(escapeKey)
            Screen('CloseAll')
            break;
        end
    end
end

%Restore PROPixx State
Datapixx('SetPropixxDlpSequenceProgram', 0);
Datapixx('RegWrRd');
Datapixx('close');
clear all;

    function L = is_angle_between(in, angle1, angle2)
        in = mod(in, 2*pi);
        L = in > angle1 & in < angle2;
        L = L | in > angle1+2*pi & in < angle2*2*pi;
        
        %        rAngle = mod(mod(angle2-angle1, 2*pi)+2*pi,2*pi);
        %        if (rAngle > pi)
        %            tmp = angle1;
        %            angle1 = angle2;
        %            angle2 = tmp;
        %        end
        %        if (angle1 <= angle2)
        %            L = in >=angle1 & in <= angle2;
        %        else
        %            L = in >=angle1 | in <= angle2;
        %        end
    end

    function m = GenerateFrame(bluerotoffset, bkgnd, dark)
        
        % Plotting the arcs
        nrings = [4 3 3];
        ring_theta = [0 pi/3; pi/3 2*pi/3; 2*pi/3 pi];
        start_ring_r = .1; % innermost ring
        ringthickness = .04;
        ringmargin = .04;
        m = zeros(256,256,3);
        [tmpx, tmpy] = meshgrid(1:size(m,1),1:size(m,2));
        tmpr = sqrt((tmpx(:)-(size(m,1)+1)/2).^2+(tmpy(:)-(size(m,2)+1)/2).^2);
        r = tmpr./(size(m,1)/2);
        theta = atan2(tmpy(:)-(size(m,2)+1)/2, tmpx(:)-(size(m,1)+1)/2);
        addedtheta = [0 0 bluerotoffset];
        for gun = 1:3
            ring_r = start_ring_r;
            % Generate a hemidisk
            tmpm = bkgnd*ones(size(m,1),size(m,2));
            tmpm(mod(theta+addedtheta(gun),2*pi) > pi & r <= 1) = dark; % lightness of solid hemidisk
            for ii = 1:length(nrings) % Looping over the three wedges
                for j = 1:nrings(ii) % looping over the individual arcs in each wedge
                    ring_r = ring_r + ringthickness+ringmargin;
                    L = r > ring_r-ringthickness/2 &...
                        r < ring_r+ringthickness/2 &...
                        is_angle_between(theta+addedtheta(gun), ring_theta(ii,1), ring_theta(ii,2));
                    tmpm(L) = dark;
                end
            end
            m(:,:,gun) = tmpm;
        end
    end
end