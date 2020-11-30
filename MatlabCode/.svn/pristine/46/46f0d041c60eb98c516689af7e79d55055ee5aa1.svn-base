% This function loads in stimuli/thresholds from IsoSamp/LMTF NEX files and
% appends those to a running matrix of all previous stimuli/thresholds.

function [data,isostim,rfxy] = extract_nex_data(data, fnames, start_at, isostim)
global gl

nfiles = length(fnames);
if start_at > nfiles
    data = gl.data;
    isostim = gl.isostim;
    rfxy = gl.rfxy;
    return
end

[new_data,new_isostim,rfxy] = get_LMTF_stimuli(fnames, start_at, nfiles, gl.NSTIMPERROUND);
if (size(unique(rfxy,'rows'),1) > 1)
    warning('Different files have different stimulus locations! Defering to user-specified location.');
    rfxy = [0 0];
else
    rfxy = rfxy(1,:);
end

data = [data; canonical_form(new_data)];
isostim = [isostim; canonical_form(new_isostim)];
if ~isempty(isostim) % there could be duplicates if more than one IsoSamp file was used
    isostim = unique(isostim, 'rows', 'stable');
end

function [new_data,isostim,rfxy] = get_LMTF_stimuli(fnames, start_at, nfiles, NSTIMPERROUND)
global gl
isostim = [];
rfxy = [];
% we expect NSTIMPERROUND stimuli per file (although there may be fewer in some early files)
new_data = zeros(NSTIMPERROUND*(nfiles-start_at+1), 3);
last_filled_row = 0;
for roundnum = start_at:nfiles
    stro = nex2stro(fnames{roundnum});
    [~,cur_file] = fileparts(fnames{roundnum});
    
    % handle an IsoSamp replay experiment
    if stro.sum.paradigmID == 107 && strcmp('LMTF', char(stro.sum.exptParams.model_str)')
        Llcc = strcmp(stro.sum.trialFields(1,:), 'stim_l');
        Lmcc = strcmp(stro.sum.trialFields(1,:), 'stim_m');
        Ltf = strcmp(stro.sum.trialFields(1,:), 'tf');
        all_stim = [stro.trial(:,Llcc|Lmcc) stro.trial(:,Ltf)];
        isostim = [isostim; unique(all_stim,'rows','stable')]; %#ok<AGROW>
        isostim(isostim(:,1) == 0 & isostim(:,2) == 0,:) = []; % remove the blank
        inc_rf = [stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y];
        
        % cache the RF location and enforce that all subsequent RFs match
        if isempty(gl.isosamp_rf) || all(inc_rf == gl.isosamp_rf)
            gl.isosamp_rf = inc_rf;
            fprintf('IsoSamp file %s has RF = (%d,%d) \n', cur_file, gl.isosamp_rf);
        else
            error('IsoSamp file %s has a different RF location than the previous IsoSamp file(s)', cur_file);
        end
    elseif stro.sum.paradigmID == 157
        Lstim_idx = strcmp(stro.sum.trialFields(1,:), 'stim_idx');
        Llcc = strcmp(stro.sum.trialFields(1,:), 'lcc');
        Lmcc = strcmp(stro.sum.trialFields(1,:), 'mcc');
        Ltf = strcmp(stro.sum.trialFields(1,:), 'tf');
        
        [~,final_stim_trial_idxs] = unique(stro.trial(:,Lstim_idx), 'last');
        quest_modes = stro.trial(final_stim_trial_idxs, Llcc|Lmcc);
        tfs = stro.trial(final_stim_trial_idxs,Ltf);
        nstim = size(quest_modes, 1);
        new_data((1:nstim)+last_filled_row,:) = [quest_modes tfs];
        last_filled_row = last_filled_row + nstim;
        rfxy = [rfxy; stro.sum.exptParams.stim_x stro.sum.exptParams.stim_y]; % in tenths of deg
        % LMTF stimulus location must match IsoSamp RFs
        if ~isempty(gl.isosamp_rf)
            stim_loc = [stro.sum.exptParams.stim_x stro.sum.exptParams.stim_y];
            if stim_loc(2) ~= gl.isosamp_rf(2) || abs(stim_loc(1)) ~= abs(gl.isosamp_rf(1))
                error('LMTF file %s has a different stimulus location than the IsoSamp file(s)', cur_file);
            end
        end
    else
        warning('%s isn''t an LMTF or IsoSamp-LMTF file! Skipping...', cur_file);
    end
    
end
new_data(last_filled_row+1:end,:) = [];


function out = canonical_form(data)
if isempty(data), out = []; return; end
% pin the final QUEST modes to the edge of the current gamut
[in_gamut,scalars] = gamut_extent(data(:,1:2));
data(~in_gamut,1:2) = data(~in_gamut,1:2).*scalars([1 1],~in_gamut)';
% put the data in the expected form for later
[th,r] = cart2pol(data(:,1), data(:,2));
th = mod(th, pi);
out = [log10(data(:,3)) th log10(r) ~in_gamut'];
