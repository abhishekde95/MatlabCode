function SLMTF_KM() %#ok<*ASGLU,*NASGU>
% Making changes to the paradigm
% INtroducing some changes to the SLMTF so as to compare the Kathy Mullen
% study
% Modified by Abhishek on 10/17 to run SLMTF
global gl
SUBJECT_IDX = 'Abhishek'; % 1: Greg, 2: Abhishek, 3: Emily
STIM_DURATION = 1000; % ms
STIM_POSITION = 0; % tenths of degrees
FPSIZE = 2; % tenths of degrees
NCONTRASTS = 50;
TIME_POSITION = [1 2];
DIST_FROM_SCREEN = 100; % 70 -ProPixx, 100- Dell8, ViewPixx
SCREEN_WIDTH = 36;  % 52 - ViewPixx, 36 - Dell 8, 51 - ProPixx

% for the audio player, couldn't get the beep work on this comp, so used this function
Fs = 10000;
t = linspace(0,0.5,6000);
y = sin(3000*t);
InitializePsychSound;
pahandle = PsychPortAudio('Open', [], [], 0, Fs,1);
PsychPortAudio('FillBuffer', pahandle, y);
cleanup = onCleanup(@() cleanup(pahandle)); %#ok<NODEF>

gl.subject_idx = SUBJECT_IDX;
gl.stim_position = STIM_POSITION;
gl.disp.ccmode = 1;
gl.dist_from_screen = DIST_FROM_SCREEN;
gl.stim.needs_update = true;
gl.is_response_pending = false;
gl.stim.scale_lattice = logspace(-2.0, 0, NCONTRASTS); % defining the contrasts here
gl.fp.xyoffset = [10 0];
gl.fp.size = FPSIZE;
gl.trials = 0; % To keep a track of trial number
gl.MAX_TRIALS = 200; % Denotes the maximum number of trials in a block of an experiment

KbName('UnifyKeyNames');
keys.show_stim = KbName('space');
keys.quit = KbName('ESCAPE');
keys.choose_left = KbName('LeftArrow');
keys.choose_right = KbName('RightArrow');
setup_keyboard_capture(keys);

Screen('Preference', 'SkipSyncTests', 2);
InitDisplay(DIST_FROM_SCREEN, SCREEN_WIDTH, '/Users/horwitzlab/Desktop/SlaveCode/Monitor Calibration/Monitor data/Dell 2/Dell2BitsCal.mat', 'T_cones_smj10.mat');
Screen('TextSize', gl.pwin, 20);

fpsizepix = round(FPSIZE*gl.disp.pixperdeg/10);
fpdrawrect = [gl.disp.center_pix(1)-floor(fpsizepix/2)+1 ...
    gl.disp.center_pix(2)-floor(fpsizepix/2)+1 ...
    gl.disp.center_pix(1)+ceil(fpsizepix/2) ...
    gl.disp.center_pix(2)+ceil(fpsizepix/2)];
fpdrawrect = fpdrawrect +repmat(gl.fp.xyoffset,[1 2])*gl.disp.pixperdeg;

gl.stim.tf = 0;
gl.stim.sp = 0;
gl.stim.sf = 0;
gl.stim.shown = cell(gl.MAX_TRIALS,1);
gl.stim.next_2_be_shown_idx = [];
gl.stim.sf_shown = zeros(gl.MAX_TRIALS,1);
gl.stim.phase_shown = zeros(gl.MAX_TRIALS,1);
gl.trial_counter = (gl.MAX_TRIALS/(4*numel(gl.stim.sf)))*ones(1,4*numel(gl.stim.sf));

gl.correct_ans_stepsize =  ones(1, 4);
setup_color_directions(gl.stim.sf);


% Note - GABOR_TEMPLATE is a 2D cell
[GABOR_TEMPLATE] = Prepare_multiple_Gabor_templates(STIM_DURATION,gl.stim.tf,gl.stim.sp,gl.stim.sf);
% for the initial trial

leave = false; nframes = 0; framecounter = 1;
textures = []; drawrects = [];

while true
    
    if ~gl.stim.on
        [key_pressed, nil, keycode] = KbCheck(-1);
        if key_pressed
            KbReleaseWait(-1);
            leave = process_keys(keys, keycode);
        end
    end
    
    % Begins after you have estimates the starting points for all the color directions
    if gl.stim.on
        if gl.stim.needs_update
            which_interval = randi(max(TIME_POSITION));  % selecting the interval in which u want the chromatic stimuli to be present
            blank_interval = find(TIME_POSITION~=which_interval);
            
            Screen('Close', textures(:));
            curr_stim_pos = STIM_POSITION;
            gl.which_interval = which_interval;
            ind = 1;
            
            sf_idx = randi(numel(gl.stim.sf));
            sp_idx = randi(size(GABOR_TEMPLATE,2));
            textures = zeros(size(GABOR_TEMPLATE{sf_idx,sp_idx},3),2);
            drawrects = zeros(4,2);
            
            gl.stim.sf_shown(gl.trials+1) = gl.stim.sf(sf_idx);
            gl.stim.phase_shown(gl.trials+1) = gl.stim.sp(sp_idx);
            [drawrects(:,blank_interval), textures(:,blank_interval)] = ShowStim(GABOR_TEMPLATE{sf_idx,sp_idx}, curr_stim_pos, 0, ...
                gl.next_scale_value(1) * gl.stim.cc_blank);  % returning the 1st stimuli
            [drawrects(:,which_interval), textures(:,which_interval), nframes] = ShowStim(GABOR_TEMPLATE{sf_idx,sp_idx}, curr_stim_pos, 0, ...
                gl.next_scale_value(1) * gl.stim.cc_dt_chrom(:,1));  % returning the 2nd stimuli
            
            framecounter = 1;
            gl.stim.needs_update = false;
            PsychPortAudio('Start', pahandle, 1, 0, 0);
        end
        DrawStimuli(textures(:,ind), drawrects(:,ind), framecounter - (floor(framecounter/nframes)*(nframes-1)));
    end
    
    if (STIM_POSITION ~= gl.fp.xyoffset(1)+1)
        %if the stimulus is at the center of the screen and the offset is 0, then don't display the FP
        Screen('FillRect', gl.pwin, [0 0 0], fpdrawrect);
    end
    Screen('Flip', gl.pwin);
    
    if gl.stim.on
        if framecounter == 2*nframes && ind == 2
            Screen('Close', textures(:));
            gl.stim.on = false;
            gl.is_response_pending = true;
            framecounter = 1;
        elseif framecounter == nframes && ind == 1
            ind = 2; % This triggers the display of the second stimulus (second interval)
        elseif mod(framecounter,nframes) == round(0.25*nframes) && ind == 2
            PsychPortAudio('Start', pahandle, 1, 0, 0);
        end
        framecounter = framecounter+1;
    end
    
    
    if leave || gl.trials == gl.MAX_TRIALS
        Screen('Close', textures(:));
        break
    end
end

function cleanup(pahandle) %#ok<DEFNU>
RestrictKeysForKbCheck([]);
% Stop playback:
PsychPortAudio('Stop', pahandle);
% Close the audio device:
PsychPortAudio('Close', pahandle);
ListenChar(0);
Screen('Preference', 'SkipSyncTests', 0);
sca();

function setup_keyboard_capture(keys)
ListenChar(2);
accepted_keys = cell2mat(struct2cell(keys));
RestrictKeysForKbCheck(accepted_keys);

function setup_color_directions(sf)
global gl
% Taking into account the number of different spatial frequencies present
N = numel(sf);
L_chrom = repmat([0.0636; 0.0636; -0.0636; -0.0636],[N 1]);
M_chrom = repmat([-0.0636; -0.0636; 0.0636; 0.0636],[N 1]);
S_chrom = repmat([-0.45; 0.45; 0.45; -0.45],[N 1]);
L_achrom = 0; M_achrom = 0; S_achrom = 0;
L = [L_chrom; L_achrom];
M = [M_chrom; M_achrom];
S = [S_chrom; S_achrom];

blank_idx = [zeros(numel(L_chrom),1); 1];
blank_idx = logical(blank_idx);
gl.color_ID = cell(1,N*4);
stim_dir1 = 'orange_sf';
stim_dir2 = 'magenta_sf';
stim_dir3 = 'cyan_sf';
stim_dir4 = 'lime_sf';

for ii = 1:N
    gl.color_ID{4*ii-3} = strcat(stim_dir1,num2str(ii));
    gl.color_ID{4*ii-2} = strcat(stim_dir2,num2str(ii));
    gl.color_ID{4*ii-1} = strcat(stim_dir3,num2str(ii));
    gl.color_ID{4*ii} = strcat(stim_dir4,num2str(ii));
end

gl.orig_color_ID = gl.color_ID;
color_directions = [L M S];
gl.stim.min_scale = gl.stim.scale_lattice(1);
gl.stim.max_scale = gl.stim.scale_lattice(end);

% pull out 1L+1M and out of gamut directions
cc_at_detection = color_directions';
gl.stim.cc_blank = cc_at_detection(:,blank_idx);
gl.stim.cc_dt_chrom = cc_at_detection(:,1:N*4);

% make counters and the vectors holding historical scale factor information
if gl.fp.xyoffset(1) == 0
    factor = 0.5; % [low sf; high sf];
else
    factor = 0.7;
end
gl.correct_count = zeros(1, size(gl.stim.cc_dt_chrom, 2)); % to hold a temporary count of the correct choices made for a stimuli of a particular contrast and a color direction

for ii = 1:N
    start_contrast_idx = ceil(factor(ii)*length(gl.stim.scale_lattice));
    gl.prev_scale_idx(4*ii-3:4*ii) = start_contrast_idx;
    gl.next_scale_value(4*ii-3:4*ii) = gl.stim.scale_lattice(start_contrast_idx);
end

% introducing 2 new variables to store global correct and incorrect counts
gl.cum_correct_count = zeros(numel(gl.stim.scale_lattice), size(gl.stim.cc_dt_chrom, 2));
gl.cum_incorrect_count = zeros(numel(gl.stim.scale_lattice), size(gl.stim.cc_dt_chrom, 2));

shuffle_chromatic_stimuli();

function leave = process_keys(keys, keycode)
global gl
leave = false;
gl.stim.needs_update = true;

if keycode(keys.show_stim)
    % Enter this 'if' condition when u press a space bar i.e. the user instructs the computer to start the next trial.
    gl.stim.on = true;
elseif keycode(keys.quit)
    % Enter this 'if' condition when u press an 'ESC' key
    gl.stim.needs_update = false;
    leave = true;
elseif gl.is_response_pending && (keycode(keys.choose_left) || keycode(keys.choose_right))
    % Enter this 'if' condition is the response from the subject is pending and the subject has pressed one of the arrow keys
    
    evaluate_trial(keycode(keys.choose_left));
    % shuffle the chromatic stimuli here
    gl.is_response_pending = false;
    gl.trials = gl.trials + 1;
    gl.stim.shown{gl.trials} = gl.color_ID{1};
    gl.trial_counter(gl.stim.next_2_be_shown_idx(end)) = gl.trial_counter(gl.stim.next_2_be_shown_idx(end)) - 1;
    if gl.trials < gl.MAX_TRIALS
        shuffle_chromatic_stimuli();
    end
    
end

function evaluate_trial(chose_first_int)
global gl
if chose_first_int && gl.which_interval ==1 || ~chose_first_int && gl.which_interval ==2
    % Enter this 'if' condition if u have answered correctly
    gl.cum_correct_count(gl.prev_scale_idx(1),1) = gl.cum_correct_count(gl.prev_scale_idx(1),1) + 1;
    gl.correct_count(1) = gl.correct_count(1)+1;
    if gl.correct_count(1) == gl.correct_ans_stepsize(1)
        % Decrease the contrast if u have answered correctly 3 times for the stimulus shown at a given contrast.
        gl.prev_scale_idx(1) = max(1, gl.prev_scale_idx(1)-1);
        gl.next_scale_value(1) = gl.stim.scale_lattice(gl.prev_scale_idx(1));
        gl.correct_count(1) = 0;
    end
else
    % Jump here if u have made a mistake
     gl.cum_incorrect_count(gl.prev_scale_idx(1),1) = gl.cum_incorrect_count(gl.prev_scale_idx(1),1) + 1;
     gl.prev_scale_idx(1) = min(length(gl.stim.scale_lattice), gl.prev_scale_idx(1)+1);
     gl.next_scale_value(1) = gl.stim.scale_lattice(gl.prev_scale_idx(1));
     gl.correct_ans_stepsize(1) = 3;
end

function shuffle_chromatic_stimuli()
% Shuffle the indexes of the chromatic stimuli
global gl
shuffle_idxs = randperm(size(gl.stim.cc_dt_chrom, 2));
col_idx = find_col_idx(shuffle_idxs(1));
while gl.trial_counter(col_idx) == 0
    shuffle_idxs = circshift(shuffle_idxs,[0 1]);
    col_idx = find_col_idx(shuffle_idxs(1));
end

gl.stim.next_2_be_shown_idx = [gl.stim.next_2_be_shown_idx; col_idx];
gl.color_ID = gl.color_ID(shuffle_idxs);
gl.prev_scale_idx = gl.prev_scale_idx(shuffle_idxs);
gl.next_scale_value = gl.next_scale_value(shuffle_idxs);
gl.correct_count = gl.correct_count(shuffle_idxs);
gl.cum_correct_count = gl.cum_correct_count(:,shuffle_idxs);
gl.cum_incorrect_count = gl.cum_incorrect_count(:,shuffle_idxs);
gl.stim.cc_dt_chrom = gl.stim.cc_dt_chrom(:,shuffle_idxs);
gl.correct_ans_stepsize = gl.correct_ans_stepsize(shuffle_idxs);

function ii = find_col_idx(idx)
global gl
out = 0;
ii = 0;
while (out == 0)
    ii = ii + 1;
    out = strcmp(gl.orig_color_ID(ii),gl.color_ID(idx));
end

function InitDisplay(mondistcm, screenwidthcm, varargin)
global gl

calfilename = varargin{1};
cal = load(calfilename); cal = cal.cals{end};
gl.disp.bkgndRGB = round(255*cal.bgColor)';
gl.calib.gammatable = cal.gammaTable;
gl.calib.invgammatable = InvertGammaTable(cal.gammaInput, cal.gammaTable, 2^16);
gl.disp.screen_width = screenwidthcm;

if nargin > 3 && ~isempty(varargin{2})
    fundfilename = varargin{2};
    s = load(fundfilename);
    fns = fieldnames(s);
    gl.calib.fundamentals = s.(fns{1})';
    wavelength_spacing = s.(fns{2});
    P_device = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls(wavelength_spacing));
    gl.calib.M = gl.calib.fundamentals' * P_device;
    gl.calib.invM = inv(gl.calib.M);
end

% start up the imaging pipeline
if ~isempty(Screen('Windows'))
    gl.pwin = max(Screen('Windows'));
    Screen('FillRect', gl.pwin, cal.bgColor);
else
    PsychImaging('PrepareConfiguration');
    if IsVPixx()
        PsychImaging('AddTask', 'General', 'EnableDataPixxC48Output', gl.disp.ccmode);
    else
        PsychImaging('AddTask', 'General', 'EnableBits++Color++Output', gl.disp.ccmode);
    end
    gl.pwin = PsychImaging('OpenWindow', 0, cal.bgColor);
end

gl.disp.refreshrate = Screen('FrameRate', gl.pwin, 1);
[gl.disp.dim_pix(1), gl.disp.dim_pix(2)] = Screen('WindowSize', gl.pwin);
gl.disp.center_pix(1) = gl.disp.dim_pix(1)/2;
gl.disp.center_pix(2) = gl.disp.dim_pix(2)/2;

gl.disp.pixpercm = gl.disp.dim_pix(1)/screenwidthcm;
theta = atan2(screenwidthcm/2, mondistcm)*180/pi;
cmperdeg = screenwidthcm/2/theta;
gl.disp.pixperdeg = gl.disp.pixpercm*cmperdeg;

gl.disp.bkgndrgb = [gl.calib.gammatable(gl.disp.bkgndRGB(1)+1,1)
    gl.calib.gammatable(gl.disp.bkgndRGB(2)+1,2)
    gl.calib.gammatable(gl.disp.bkgndRGB(3)+1,3)];

gl.stim.on = false;
HideCursor();

function template = Prepare_multiple_Gabor_templates(STIM_DURATION,tf,sp,sf)
N = numel(sf);
M = numel(sp);
template = cell(1,N);
for ii = 1:N
    for jj = 1:M
    template{ii,jj} = PrepareGabor(STIM_DURATION, tf,sp(jj),sf(ii));
    end
end

function template = PrepareGabor(stim_duration, tf, sp, sf)
% tf -> temporal frequency
% sp -> spatial phase
global gl

sigma = 0.5;
nsigmas = 2;
theta = 0*pi/2;
phi = sp;
ggamma = 1/1.8;

gl.gabor.sigma = sigma;
gl.gabor.nsigmas = nsigmas;
gl.gabor.theta = pi/2;
gl.gabor.ggamma = ggamma;

nframestotal = floor(stim_duration*gl.disp.refreshrate/1000);
nframesramp = floor(nframestotal/4);
nframesplateau = nframestotal - 2*nframesramp;

ramp = linspace(0, 1, nframesramp);
plateau = ones(1, nframesplateau);
temporalprofile = [ramp plateau fliplr(ramp)];
nframes = length(temporalprofile);
stimsizeindeg_x = sigma*nsigmas;
stimsizeinpix_x = round(stimsizeindeg_x*gl.disp.pixperdeg);
stimsizeindeg_y = sigma*nsigmas/ggamma;
stimsizeinpix_y = round(stimsizeindeg_y*gl.disp.pixperdeg);
[x,y] = meshgrid(stimsizeindeg_x*linspace(-1, 1, stimsizeinpix_x), ...
    stimsizeindeg_y*linspace(-1, 1, 2 * stimsizeinpix_y));

X = x*cos(-theta) + y*sin(-theta);
Y =-x*sin(-theta) + y*cos(-theta);

deltaphase = 2*pi*tf/gl.disp.refreshrate;
phases = phi + (0:nframes-1)*deltaphase;
phases = reshape(phases, [1 1 nframes]);
temporalprofile = reshape(temporalprofile, [1 1 nframes]);
expterm = bsxfun(@times, exp(-(X.^2+ ggamma.^2*Y.^2)/2/sigma^2), temporalprofile);
costerm = cos(bsxfun(@plus, 2*pi*Y*sf, phases));
template = expterm .* costerm;

function [drawrect, textures, nframes] = ShowStim(template, stimx, stimy, stimconecontrast)
global gl

stimx = stimx/10;
stimy = stimy/10;

x = stimx*gl.disp.pixperdeg;
y = stimy*gl.disp.pixperdeg;
stimsizeinpix_y  = size(template,1)/2;
stimsizeinpix_x = size(template,2);



drawrect = round([gl.disp.center_pix(1) + x - stimsizeinpix_x ...
    gl.disp.center_pix(2) - y - stimsizeinpix_y ...
    gl.disp.center_pix(1) + x + stimsizeinpix_x ...
    gl.disp.center_pix(2) - y + stimsizeinpix_y]);

if rem(drawrect(1), 2)
    drawrect([1 3]) = drawrect([1 3]) - 1;
end

%assert(drawrect(3) - drawrect(1) == 2*stimsizeinpix ...
 %   && drawrect(4) - drawrect(2) == 2*stimsizeinpix, ...
  %  'Incorrect draw window size!');

bkgndlms = gl.calib.M * gl.disp.bkgndrgb;
gl.stim.rgb = gl.calib.invM * (bkgndlms .* (1 + stimconecontrast));
%sca;
%keyboard;

nframes = size(template, 3);
textures = zeros(nframes, 1);

NGAMMASTEPS = size(gl.calib.invgammatable, 1);
for ii = 1:nframes
    frame = template(:,:,ii);
    img = zeros(size(template, 1), size(template, 2), 3);

    for plane = 1:3
        tmp = frame * (gl.stim.rgb(plane) - gl.disp.bkgndrgb(plane)) + gl.disp.bkgndrgb(plane);
        tmp = round(tmp * (NGAMMASTEPS - 1)) + 1;
        tmp = gl.calib.invgammatable(tmp,plane);
        img(:,:,plane) = reshape(tmp, size(img(:,:,1)));
    end
    textures(ii) = Screen('MakeTexture', gl.pwin, TranslateToColourMode(img, [], gl.disp.ccmode), [], [], 2);
end

function DrawStimuli(textures, drawrects, frame_num)
global gl
Screen('DrawTextures', gl.pwin, textures(frame_num,:), [], drawrects, [], 0);

