function LMTF_2stim() %#ok<*ASGLU,*NASGU>
    % Making changes to the paradigm
    global gl
    cleanup = onCleanup(@() cleanup()); %#ok<NODEF>

    SUBJECT_IDX = 2; % 1: Greg, 2: Zack, 3: Leah
    STIM_DURATION = 1000; % ms
    STIM_POSITION = 30; % tenths of degrees
    FPSIZE = 2; % tenths of degrees
    NCOLORDIRS = 2;
    NCONTRASTS = 35;
    DETECTION_ONLY = false;

    gl.disp.ccmode = 1;
    gl.stim.needs_update = true;
    gl.is_response_pending = false;

    model_params = load('DTNTmodelparamsGZL.mat');
    model_params = model_params.fpars(:,SUBJECT_IDX);

    KbName('UnifyKeyNames');
    keys.show_stim = KbName('space');
    keys.quit = KbName('ESCAPE');
    keys.choose_left = KbName('LeftArrow');
    keys.choose_right = KbName('RightArrow');
    setup_keyboard_capture(keys);

    Screen('Preference', 'SkipSyncTests', 2);
    InitDisplay(100, 52, 'ProPixx.mat', 'T_cones_smj10.mat');
    Screen('TextSize', gl.pwin, 20);

    fpsizepix = round(FPSIZE*gl.disp.pixperdeg/10);
    fpdrawrect = [gl.disp.center_pix(1)-floor(fpsizepix/2)+1 ...
        gl.disp.center_pix(2)-floor(fpsizepix/2)+1 ...
        gl.disp.center_pix(1)+ceil(fpsizepix/2) ...
        gl.disp.center_pix(2)+ceil(fpsizepix/2)];

    setup_color_directions(NCOLORDIRS, NCONTRASTS, model_params);

    gl.stim.tf = 0;
    GABOR_TEMPLATE = PrepareGabor(STIM_DURATION, gl.stim.tf);
    stim_positions = STIM_POSITION*[-1 1];

    leave = false; nframes = 0; framecounter = 1;
    textures = []; drawrects = [];
    while true
        if ~gl.stim.on
            [key_pressed, nil, keycode] = KbCheck(-1);
            if key_pressed
                KbReleaseWait(-1);
                leave = process_keys(keys, keycode);
            end
            DrawFormattedText(gl.pwin, sprintf('%0.3f %0.3f @ %0.3fx DT, TF: %d Hz', ...
                gl.stim.cc_dt_chrom(1:2,1), gl.next_scale_value(1), gl.stim.tf(1)), ...
                'center', ceil(.8*gl.disp.dim_pix(2)), 2);
        end

        if gl.stim.on
            if gl.stim.needs_update
                Screen('Close', textures(:));
                curr_stim_pos = stim_positions(randperm(2));
                gl.chromatic_pos = curr_stim_pos(1);

                [drawrects(:,1), textures(:,1), nframes] = ShowStim(GABOR_TEMPLATE, curr_stim_pos(1), 0, ...
                    gl.next_scale_value(1) * gl.stim.cc_dt_chrom(:,1));
                if ~DETECTION_ONLY
                    [drawrects(:,2), textures(:,2)] = ShowStim(GABOR_TEMPLATE, curr_stim_pos(2), 0, ...
                        gl.next_scale_value(1) * gl.stim.cc_dt_achrom);
                end
                framecounter = 1;
                gl.stim.needs_update = false;
            end
            DrawStimuli(textures, drawrects, framecounter);
        end

        Screen('FillRect', gl.pwin, [0 0 0], fpdrawrect);
        Screen('Flip', gl.pwin);

        if gl.stim.on
            framecounter = framecounter+1;
            if framecounter >= nframes
                Screen('Close', textures(:));
                gl.stim.on = false;
                gl.is_response_pending = true;
                framecounter = 1;
            end
        end

        if leave
            Screen('Close', textures(:));
            break
        end
    end

function cleanup() %#ok<DEFNU>
    RestrictKeysForKbCheck([]);
    ListenChar(0);
    Screen('Preference', 'SkipSyncTests', 0);
    sca();

function setup_keyboard_capture(keys)
    ListenChar(2);
    accepted_keys = cell2mat(struct2cell(keys));
    RestrictKeysForKbCheck(accepted_keys);

function setup_color_directions(NCOLORDIRS, NCONTRASTS, model_params)
    global gl
    color_angles = linspace(0, pi, NCOLORDIRS+1)' + 135; % -L-M and L-M direction
    color_angles(1) = [];

    if rem(NCOLORDIRS, 4) % add 1L+1M angle manually
        color_angles(end+1) = pi/4;
    end
    achromatic_idx = abs(color_angles - pi/4) < 2*eps;

    color_directions = [cos(color_angles) sin(color_angles) ones(length(color_angles), 1)]; % Adding the S component here
    [th,phi,nil] = cart2sph(color_directions(:,1), color_directions(:,2), color_directions(:,3));
    wholesum = sum(abs([cos(phi).*cos(th) cos(phi).*sin(th) sin(phi)]*reshape(model_params(2:end), [3 3])).^model_params(1), 2);
    scales_to_dtsurf = (1./wholesum).^(1/model_params(1))/100;

    cc_at_detection = bsxfun(@times, color_directions, scales_to_dtsurf)';

    [nil,max_scales_to_dtsurf] = cellfun(@(x) gamutCheck(x, gl.disp.bkgndrgb, gl.calib.M, 'both'), ...
        num2cell(cc_at_detection, 1), 'unif', 0);
    max_scales_to_dtsurf = [max_scales_to_dtsurf{:}];

    assert(max_scales_to_dtsurf(achromatic_idx) > 1, ...
        'The predicted detection threshold at 1L+1M is out of gamut?!');

    min_scale_to_dtsurf = min(max_scales_to_dtsurf(max_scales_to_dtsurf > 1));
    gl.stim.scale_lattice = linspace(.75, min_scale_to_dtsurf, NCONTRASTS);
    gl.stim.min_scale = gl.stim.scale_lattice(1);
    gl.stim.max_scale = gl.stim.scale_lattice(end);

    % pull out 1L+1M and out of gamut directions
    gl.stim.cc_dt_achrom = cc_at_detection(:,achromatic_idx);
    gl.stim.cc_dt_chrom = cc_at_detection(:,max_scales_to_dtsurf' > 1 & ~achromatic_idx);

    % make counters and the vectors holding historical scale factor information
    gl.correct_count = zeros(1, size(gl.stim.cc_dt_chrom, 2));
    gl.incorrect_count = zeros(1, size(gl.stim.cc_dt_chrom, 2));
    middle = ceil(length(gl.stim.scale_lattice)/2);
    gl.prev_scale_idx = middle + zeros(1, size(gl.stim.cc_dt_chrom, 2));
    gl.next_scale_value = gl.stim.scale_lattice(middle) + zeros(1, size(gl.stim.cc_dt_chrom, 2));

    shuffle_chromatic_stimuli();

function leave = process_keys(keys, keycode)
    global gl
    leave = false;
    gl.stim.needs_update = true;

    if keycode(keys.show_stim)
        gl.stim.on = true;
    elseif keycode(keys.quit)
        gl.stim.needs_update = false;
        leave = true;
    elseif gl.is_response_pending && (keycode(keys.choose_left) || keycode(keys.choose_right))
        evaluate_trial(keycode(keys.choose_left));
        shuffle_chromatic_stimuli();
        gl.is_response_pending = false;
    end

function evaluate_trial(chose_left)
    global gl
    if chose_left && gl.chromatic_pos < 0 || ~chose_left && gl.chromatic_pos > 0
        gl.correct_count(1) = gl.correct_count(1)+1;
        if gl.correct_count(1) == 3 % decrease contrast
            gl.prev_scale_idx(1) = max(1, gl.prev_scale_idx(1)-2);
            gl.next_scale_value(1) = gl.stim.scale_lattice(gl.prev_scale_idx(1));
            gl.correct_count(1) = 0;
            gl.incorrect_count(1) = 0;
        end
    else
        gl.incorrect_count(1) = gl.incorrect_count(1)+1;
        if gl.incorrect_count(1) > 3
            fprintf('    LMTF: %0.5f %0.5f finished at %0.5fx DT\n', ...
                gl.stim.cc_dt_chrom(1:2,1), gl.next_scale_value(1));
            gl.prev_scale_idx(1) = [];
            gl.next_scale_value(1) = [];
            gl.correct_count(1) = [];
            gl.incorrect_count(1) = [];
            gl.stim.cc_dt_chrom(:,1) = [];
            assert(~isempty(gl.stim.cc_dt_chrom), 'All done.');
        else
            gl.prev_scale_idx(1) = min(length(gl.stim.scale_lattice), gl.prev_scale_idx(1)+1);
            gl.next_scale_value(1) = gl.stim.scale_lattice(gl.prev_scale_idx(1));
        end
    end

function shuffle_chromatic_stimuli()
    global gl

    shuffle_idxs = randperm(size(gl.stim.cc_dt_chrom, 2));
    gl.prev_scale_idx = gl.prev_scale_idx(shuffle_idxs);
    gl.next_scale_value = gl.next_scale_value(shuffle_idxs);
    gl.correct_count = gl.correct_count(shuffle_idxs);
    gl.incorrect_count = gl.incorrect_count(shuffle_idxs);
    gl.stim.cc_dt_chrom = gl.stim.cc_dt_chrom(:,shuffle_idxs);

function InitDisplay(mondistcm, screenwidthcm, varargin)
    global gl

    calfilename = varargin{1};
    cal = load(calfilename); cal = cal.cals{end};
    gl.disp.bkgndRGB = round(255*cal.bgColor)';
    gl.calib.gammatable = cal.gammaTable;
    gl.calib.invgammatable = InvertGammaTable(cal.gammaInput, cal.gammaTable, 2^16);

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

function template = PrepareGabor(stim_duration, tf)
    global gl

    sf = 5;
    sigma = .15;
    nsigmas = 2;
    theta = pi/2;
    phi = 0;
    ggamma = 1;

    nframestotal = floor(stim_duration*gl.disp.refreshrate/1000);
    nframesramp = floor(nframestotal/4);
    nframesplateau = nframestotal - 2*nframesramp;

    ramp = linspace(0, 1, nframesramp);
    plateau = ones(1, nframesplateau);
    temporalprofile = [ramp plateau fliplr(ramp)];
    nframes = length(temporalprofile);
    stimsizeindeg = sigma*nsigmas;
    stimsizeinpix = round(stimsizeindeg*gl.disp.pixperdeg);
    [x,y] = meshgrid(stimsizeindeg*linspace(-1, 1, stimsizeinpix), ...
        stimsizeindeg*linspace(-1, 1, 2 * stimsizeinpix));

    X = x*cos(-theta) + y*sin(-theta);
    Y =-x*sin(-theta) + y*cos(-theta);

    deltaphase = 2*pi*tf/gl.disp.refreshrate;
    phases = phi + (0:nframes-1)*deltaphase;
    phases = reshape(phases, [1 1 nframes]);
    temporalprofile = reshape(temporalprofile, [1 1 nframes]);
    expterm = bsxfun(@times, exp(-(X.^2 + ggamma^2*Y.^2)/2/sigma^2), temporalprofile);
    costerm = cos(bsxfun(@plus, 2*pi*Y*sf, phases));
    template = expterm .* costerm;

function [drawrect, textures, nframes] = ShowStim(template, stimx, stimy, stimconecontrast)
    global gl

    stimx = stimx/10;
    stimy = stimy/10;

    x = stimx*gl.disp.pixperdeg;
    y = stimy*gl.disp.pixperdeg;
    stimsizeinpix = size(template, 2);

    drawrect = round([gl.disp.center_pix(1) + x - stimsizeinpix ...
        gl.disp.center_pix(2) - y - stimsizeinpix ...
        gl.disp.center_pix(1) + x + stimsizeinpix ...
        gl.disp.center_pix(2) - y + stimsizeinpix]);

    if rem(drawrect(1), 2)
        drawrect([1 3]) = drawrect([1 3]) - 1;
    end

    assert(drawrect(3) - drawrect(1) == 2*stimsizeinpix ...
        && drawrect(4) - drawrect(2) == 2*stimsizeinpix, ...
        'Incorrect draw window size!');

    bkgndlms = gl.calib.M * gl.disp.bkgndrgb;
    gl.stim.rgb = gl.calib.invM * (bkgndlms .* (1 + stimconecontrast));

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
