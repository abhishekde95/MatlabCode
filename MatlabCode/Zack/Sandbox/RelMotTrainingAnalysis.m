function out = RelMotTrainingAnalysis(stro, plots)
if nargin < 1 || isempty(stro)
    stro = nex2stro;
end

if nargin < 2 || isempty(plots)
    plots = 0;
end

targxidx = strcmp(stro.sum.trialFields(1,:),'targ_x');
targyidx = strcmp(stro.sum.trialFields(1,:),'targ_y');
correct = logical(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct')));
hepidx = strcmp(stro.sum.rasterCells(1,:),'AD11');
vepidx = strcmp(stro.sum.rasterCells(1,:),'AD12');
t_eyestart = cell2mat(stro.ras(:,strcmp(stro.sum.rasterCells(1,:),'anlgStartTime')));
t_fpacq = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpacq_t'));
t_sac = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targoff_t'));
rf = [stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y];
targ_w = 2 * stro.sum.exptParams.targwin_x / 10;
targ_h = 2 * stro.sum.exptParams.targwin_y / 10;
fp_w = 2 * stro.sum.exptParams.eyewin_x / 10;
fp_h = 2 * stro.sum.exptParams.eyewin_y / 10;
ntrials = size(correct,1);

targets = stro.trial(:, targxidx | targyidx);
target_ids = get_target_ids(targets, rf);
ntargets = length(unique(target_ids));

if plots
    MV_TO_DEG = 4096 / 400;

    H_eye = stro.ras(:,hepidx)';
    V_eye = stro.ras(:,vepidx)';
    
    for i = 1:length(H_eye)
        t = (0:length(H_eye{i})-1) ./ stro.sum.analog.storeRates{1}+t_eyestart(i);
        L = t >= t_fpacq(i) & t <= t_sac(i);
        H_eye{i} = H_eye{i}(L);
        V_eye{i} = V_eye{i}(L);
    end
    
    max_trace_length = max(cellfun(@length, H_eye));
    H_eye = cell2mat(cellfun(@(x) cat(1, x, nan(max_trace_length-length(x), 1)), ...
        H_eye, 'unif', 0));
    V_eye = cell2mat(cellfun(@(x) cat(1, x, nan(max_trace_length-length(x), 1)), ...
        V_eye, 'unif', 0));

    H_eye = H_eye * MV_TO_DEG;
    V_eye = V_eye * MV_TO_DEG;

    draw_targets = unique(targets, 'rows') / 10;
    % figure; axis(5*[-1 1 -1 1]); hold on;
    % plot(H_eye(:,~correct), V_eye(:,~correct), 'r:');
    % plot(H_eye(:,correct), V_eye(:,correct), 'g-', 'linewidth', 1.5);

    figure; axis(5*[-1 1 -1 1]); hold on; set(gca, 'xtick', [], 'ytick', [], 'xcolor', [1 1 1], 'ycolor', [1 1 1]);
    patchline(H_eye(:,correct), V_eye(:,correct), 'linestyle', '-', 'edgecolor', 'g', 'linewidth', 2, 'edgealpha', 0.15);

    for target = draw_targets'
        rectangle('pos', [target(1)-targ_w/2 target(2)-targ_h/2 targ_w targ_h]);
    end
    rectangle('pos', [stro.sum.exptParams.fp_x/10-fp_w/2 stro.sum.exptParams.fp_y/10-fp_h/2 fp_w fp_h]);
end

%#ok<*NOPRT>
correct_trials_per_target_stim = hist(target_ids(correct & stro.trial(:,10)), ntargets);
ntrials_per_target_stim = hist(target_ids(logical(stro.trial(:,10))), ntargets);
correct_trials_per_target_nostim = hist(target_ids(correct & ~stro.trial(:,10)), ntargets);
ntrials_per_target_nostim = hist(target_ids(~logical(stro.trial(:,10))), ntargets);
no_stim_frac = correct_trials_per_target_nostim ./ ntrials_per_target_nostim
stim_frac = correct_trials_per_target_stim ./ ntrials_per_target_stim

correct_trials_per_target = hist(target_ids(correct), ntargets);
ntrials_per_target = hist(target_ids, ntargets);
[h,p] = equalproptest(correct_trials_per_target, ntrials * ones(1,4), 0.05);

out{1} = [h p];
out{2} = 1-binocdf(sum(correct_trials_per_target), ntrials, .25);
out{3} = correct_trials_per_target./ntrials_per_target;
out{4} = 1-binocdf(correct_trials_per_target, ntrials_per_target, .25);
end

% This function maps the target locations to a numeric identifier 1 - 4. The
% Cartesian quadrant of `rf` maps to id 1 and the rest are assigned
% anti-clockwise from that.
function ids = get_target_ids(targets, rf)
rf_quad = get_cartesian_quadrant(rf);
targ_quads = get_cartesian_quadrant(targets);
% `rf_quad - 1` is the number of clockwise rotations necessary to map the rf's
% quadrant as #1
ids = mod(targ_quads - 1 - (rf_quad - 1), 4) + 1;

    function quads = get_cartesian_quadrant(points)
        signs = sign(points); signs(signs > 0) = 0;
        quads = xor(signs(:,1), signs(:,2)) - 2*signs(:,2) + 1;
    end
end
