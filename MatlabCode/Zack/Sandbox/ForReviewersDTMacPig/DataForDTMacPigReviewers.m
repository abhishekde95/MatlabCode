function DataForDTMacPigReviewers(expts)
if nargin < 1 || isempty(expts), expts = [1 2 3]; end
if any(expts < 1 | expts > 3), return; end

if any(expts == 1 | expts == 2)
    s = load('T_cones_smj10.mat');
    cone_sensitivities = s.T_cones_smj10';
    S_cones = s.S_cones_smj10;
    wavelengths = linspace(S_cones(1), ...
        S_cones(1)+S_cones(2)*(S_cones(3)-1), ...
        size(cone_sensitivities, 1))'; %#ok<NASGU>
end

if any(expts == 3)
    scotcal = load('Dell4blackbkgnd.mat');
    scotcal = scotcal.cals{end};
    nd_transmittance = load('wrattennd1.mat');
    nd_transmittance = nd_transmittance.wratten_nd1(:,2);
    s = load('T_rods');
    scotopic_lum_eff = s.T_rods';
    S_rods = s.S_rods;
    wavelengths = linspace(S_rods(1), S_rods(1)+S_rods(2)*(S_rods(3)-1), ...
        size(scotopic_lum_eff, 1))'; %#ok<NASGU>
    scot_P_device = SplineSpd(scotcal.S_device, scotcal.P_device, S_rods);
    new_mon_spd = scot_P_device.*nd_transmittance(:,ones(1,size(scot_P_device,2))).^6;
    scotlum = 1700*scotopic_lum_eff'*new_mon_spd*5/2;
end

for expt = expts(:)'
    switch expt
        case 1
            textfiles_path = nexfilepath('nexfilelists', 'Greg', 'DTEM');
            text_filenames = dir(fullfile(textfiles_path, '*MacPig.txt'));
        case 2
            textfiles_path = nexfilepath('nexfilelists', 'Greg', 'DTNT');
            text_filenames = dir(fullfile(textfiles_path, '*15Hz.txt'));
        case 3
            textfiles_path = nexfilepath('nexfilelists', 'Greg', 'DTEM');
            text_filenames = dir(fullfile(textfiles_path, '*DTscot.txt'));
    end
    text_files = strcat(textfiles_path, filesep, {text_filenames.name});
    % exclude data that are not part of this study
    text_files(~cellfun('isempty', regexp(text_files, 'Pangu|Sanity'))) = [];
    
    trials(1:length(text_files)) = struct('subject', '', 'data', []);
    monitor_spectra = [];
    background_cc = [];
    subject_counter = 0;
    for text_file = text_files
        subject_counter = subject_counter+1;
        nominal_filenames = flatten(fnamesFromTxt2(text_file{1}));
        valid_names = isvalidnexfilename(nominal_filenames);
        nex_filepaths = findfile(nominal_filenames(valid_names));
        [~,sample_filename] = fileparts(nex_filepaths{1});
        
        if any(strcmp(sample_filename(1), {'G' 'L' 'Z'}))
            subjectID = ['Human ' sample_filename(1)];
        else
            subjectID = ['Monkey ' sample_filename(1)];
        end
        
        trials(subject_counter).subject = subjectID;
        for nex_file = nex_filepaths
            try
                stro = nex2stro(nex_file{1});
            catch
                stro = oldnex2stro(nex_file{1});
            end
            
            ntrials = size(stro.trial, 1);
            bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b]';
            correct = logical(stro.trial(:,strcmp(stro.sum.trialFields(1,:), 'correct')));
            cdir_idx = stro.trial(:,strcmp(stro.sum.trialFields(1,:), 'color_dir'));
            
            if expt ~= 3
                modes = stro.trial(:,strcmp(stro.sum.trialFields(1,:), 'quest_thresh'))/100;
                flash_x = stro.trial(:,strcmp(stro.sum.trialFields(1,:), 'flash_x'))/10;
                spd = stro.sum.exptParams.mon_spect;
                M = reshape(stro.sum.exptParams.m_mtx,3,3);
                bkgndcc = M*bkgndrgb;
                color_dirs = reshape(stro.sum.exptParams.RF_colors, [], 3)';
                color_dirs(softEq([0 0 0], color_dirs, [], 'rows'), :) = [];
            else
                modes = stro.trial(:,strcmp(stro.sum.trialFields(1,:), 'quest_mode'));
                flash_x = stro.trial(:,strcmp(stro.sum.trialFields(1,:), 'stim_x'))/10;
                reds = logical(stro.trial(:,strcmp(stro.sum.trialFields(1,:), 'stim_r')));
                greens = logical(stro.trial(:,strcmp(stro.sum.trialFields(1,:), 'stim_g')));
                blues = logical(stro.trial(:,strcmp(stro.sum.trialFields(1,:), 'stim_b')));
                
                scotlum_trial = modes;
                scotlum_trial(greens) = scotlum_trial(greens)*scotlum(2);
                scotlum_trial(blues) = scotlum_trial(blues)*scotlum(3);
                cdir_idx(greens) = 1;
                cdir_idx(blues) = 2;
                
                % remove red trials
                scotlum_trial(reds) = [];
                correct(reds) = [];
                flash_x(reds) = [];
                cdir_idx(reds) = [];
                
                spd = scot_P_device(:);
                bkgndcc = [0;0;0];
            end
            
            spd_idx = all(softEq(monitor_spectra, spd(:,ones(size(monitor_spectra, 2), 1))), 1);
            bkgnd_idx = all(softEq(background_cc, bkgndcc(:,ones(size(background_cc, 2), 1))), 1);
            
            if any(spd_idx)
                spd_idx = find(spd_idx, 1);
            else
                monitor_spectra = [monitor_spectra spd]; %#ok<AGROW>
                spd_idx = size(monitor_spectra, 2);
            end
            
            if any(bkgnd_idx)
                bkgnd_idx = find(bkgnd_idx, 1);
            else
                background_cc = [background_cc bkgndcc]; %#ok<AGROW>
                bkgnd_idx = size(background_cc, 2);
            end
            
            if expt ~= 3
                ccs = bsxfun(@times, color_dirs(cdir_idx,:), modes./sqrt(sum(color_dirs(cdir_idx,:).^2, 2)));
                trials(subject_counter).data = [trials(subject_counter).data; ...
                    ccs flash_x correct spd_idx(ones(ntrials, 1)) bkgnd_idx(ones(ntrials, 1))];
            else
                trials(subject_counter).data = [trials(subject_counter).data; ...
                    scotlum_trial flash_x correct cdir_idx spd_idx(ones(size(correct)))];
            end
        end
    end
    
    if expt ~= 3
        save(fullfile(fileparts(mfilename('fullpath')), sprintf('Experiment%d.mat', expt)), ...
            'trials', 'cone_sensitivities', 'monitor_spectra', 'background_cc', 'wavelengths');
    else
        save(fullfile(fileparts(mfilename('fullpath')), sprintf('Experiment%d.mat', expt)), ...
            'trials', 'scotopic_lum_eff', 'monitor_spectra', 'nd_transmittance', 'wavelengths');
    end
    clear trials
end
