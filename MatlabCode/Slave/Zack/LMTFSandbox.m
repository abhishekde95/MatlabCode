function out = LMTFSandbox()
    global gl
    cleanup = onCleanup(@() cleanup()); %#ok<NODEF>
    STIM_DURATION = 1000; % ms
    NCOLORDIRS = 24;
    NCONTRASTS = 30+1;

    gl.disp.ccmode = 1;
    gl.stim.needs_update = true;
    gl.stim.framecounter = 0;
    gl.stim.framecountermax = 0;
    gl.saved_data = [];
    out = [];

    KbName('UnifyKeyNames');
    keys.show_stim = KbName('space');
    keys.quit = KbName('ESCAPE');
    keys.cycle_tf = KbName('t');
    keys.dec_contrast = KbName('q');
    keys.inc_contrast = KbName('w');
    keys.ccw_colordir = KbName('a');
    keys.cw_colordir = KbName('s');
    keys.save_cc = KbName('ENTER');
    setup_keyboard_queue(keys);

    Screen('Preference', 'SkipSyncTests', 2);
    InitDisplay(100, 37, 'Dell4BitsCal.mat', 'T_cones_smj10.mat');
    Screen('TextSize', gl.pwin, 20);

    gl.stim.color_dirs = [  sin(linspace(0, 2*pi, NCOLORDIRS+1))
                            cos(linspace(0, 2*pi, NCOLORDIRS+1))
                            zeros(1, NCOLORDIRS+1)
                        ];
    gl.stim.color_dirs(:,end) = [];

    scale_factors = zeros(length(gl.stim.color_dirs), 1);
    for ii = 1:length(scale_factors)
        [nil,scale_factors(ii)] = gamutCheck(gl.stim.color_dirs(:,ii), gl.disp.bkgndrgb, gl.calib.M, 'both');
    end
    gl.stim.scales = cell2mat(arrayfun(@(x) {linspace(x, 0, NCONTRASTS)'}, scale_factors)');
    gl.stim.cc = gl.stim.scales(1,1) * gl.stim.color_dirs(:,1);

    gl.stim.tf = [3 15 25]; % TODO: values based on refresh rate?

    leave = false;
    textures = []; drawrect = [];
    while true
        if ~gl.stim.on
            [key_pressed, nil, nil, release_time] = KbQueueCheck();
            if key_pressed
                leave = process_keys(keys, release_time);
            end

            DrawFormattedText(gl.pwin, sprintf('CDIR: [%0.3f %0.3f], CONTRAST: %0.3f, TF: %d Hz', ...
                gl.stim.color_dirs(1:2,1), norm(gl.stim.cc), gl.stim.tf(1)), ...
                'center', ceil(.8*gl.disp.dim_pix(2)), 1);
        end

        if gl.stim.on
            if gl.stim.needs_update
                Screen('Close', textures);
                template = PrepareGabor(STIM_DURATION, gl.stim.tf(1));
                [drawrect, textures] = ShowStim(template, 0, 0, gl.stim.cc);
                gl.stim.needs_update = false;
            end
            DrawStim(textures, drawrect);
        end

        Screen('Flip', gl.pwin);

        if gl.stim.on
            gl.stim.framecounter = gl.stim.framecounter + 1;
            if gl.stim.framecounter >= gl.stim.framecountermax
                Screen('Close', textures);
                textures = [];
                gl.stim.on = false;
            end
        end

        if leave
            Screen('Close', textures);
            out = gl.saved_data;
            break
        end
    end

function cleanup() %#ok<DEFNU>
    KbQueueStop();
    KbQueueFlush();
    KbQueueRelease();
    ListenChar(0);
    Screen('Preference', 'SkipSyncTests', 0);
    sca();

function setup_keyboard_queue(keys)
    ListenChar(2);
    keylist = zeros(1, 256);
    keylist(cell2mat(struct2cell(keys))) = 1;
    KbQueueCreate([], keylist);
    KbQueueStart();

function leave = process_keys(keys, release_time)
    global gl
    leave = false;
    gl.stim.needs_update = true;

    if release_time(keys.show_stim)
        gl.stim.on = true;
    elseif release_time(keys.quit)
        gl.stim.needs_update = false;
        leave = true;
    elseif release_time(keys.save_cc)
        gl.stim.needs_update = false;
        gl.saved_data = [gl.saved_data struct('cc', gl.stim.cc, 'tf', gl.stim.tf(1))];
    elseif release_time(keys.cycle_tf)
        gl.stim.tf = gl.stim.tf([2:end 1]);
    elseif release_time(keys.dec_contrast)
        gl.stim.scales = gl.stim.scales([2:end 1],:);
    elseif release_time(keys.inc_contrast)
        gl.stim.scales = gl.stim.scales([end 1:end-1],:);
    elseif release_time(keys.cw_colordir)
        gl.stim.color_dirs = gl.stim.color_dirs(:,[2:end 1]);
        gl.stim.scales = gl.stim.scales(:,[2:end 1]);
    elseif release_time(keys.ccw_colordir)
        gl.stim.color_dirs = gl.stim.color_dirs(:,[end 1:end-1]);
        gl.stim.scales = gl.stim.scales(:,[end 1:end-1]);
    end

    if gl.stim.needs_update
        gl.stim.cc = gl.stim.scales(1,1) * gl.stim.color_dirs(:,1);
    end

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

    sf = 1;
    sigma = .4;
    nsigmas = 3;
    theta = 0;
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

function [drawrect, textures] = ShowStim(template, stimx, stimy, stimconecontrast)
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

    gl.stim.framecountermax = size(template, 3);
    textures = zeros(1, gl.stim.framecountermax);

    NGAMMASTEPS = size(gl.calib.invgammatable, 1);
    for ii = 1:gl.stim.framecountermax
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
    gl.stim.framecounter = 0;

function DrawStim(textures, drawrect)
    global gl
    Screen('DrawTexture', gl.pwin, textures(gl.stim.framecounter+1), [], drawrect, [], 0);
