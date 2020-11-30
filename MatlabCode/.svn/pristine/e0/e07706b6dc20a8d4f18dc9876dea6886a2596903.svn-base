% This script was inspired by Greg and Charlie's QUEST unpacking code. In brief,
% the goal is to return the contrast presented on every trial of one particular
% LMTF experiment given via the first argument as a stro. The script will also
% produce a figure, which you can optionally turn off by passing 'false' as the
% second argument.

function presented_contrast = LMTF_QUEST_unpack(stro, show_plots)
if stro.sum.paradigmID ~= 157, return; end

if nargin < 2 || isempty(show_plots)
    show_plots = true;
end

spectra = reshape(stro.sum.exptParams.mon_spd, [], 3);
fundamentals = reshape(stro.sum.exptParams.fundamentals, [], 3);
P_device = SplineSpd(linspace(380, 780, size(spectra, 1))', spectra, ...
    linspace(380, 780, size(fundamentals, 1))');
M = fundamentals'*P_device;

Lstim_idx = strcmp(stro.sum.trialFields(1,:), 'stim_idx');
Llcc = strcmp(stro.sum.trialFields(1,:), 'lcc');
Lmcc = strcmp(stro.sum.trialFields(1,:), 'mcc');
Ltf = strcmp(stro.sum.trialFields(1,:), 'tf');
Loog = strcmp(stro.sum.trialFields(1,:), 'oog');

[stim_idxs,init_stim_trial_idxs] = unique(stro.trial(:,Lstim_idx), 'first');
nstim = length(stim_idxs);
% on which trials is stimulus j presented?
stim_trial_positions = bsxfun(@eq, stro.trial(:,Lstim_idx), stim_idxs');
assert(all(sum(stim_trial_positions) == size(stro.trial, 1)/nstim), ...
    'The stimuli weren''t equally distributed!');

stimuli = stro.trial(:,Llcc|Lmcc);
requested_contrast = sqrt(sum(stimuli.^2, 2)); % QUEST modes, not necessarily the presented contrast
tfs = stro.trial(init_stim_trial_idxs,Ltf);

unit_LM = bsxfun(@rdivide, stimuli(init_stim_trial_idxs,:), ...
    requested_contrast(init_stim_trial_idxs));

[~,max_cc_scales] = gamutCheck([unit_LM zeros(size(unit_LM,1), 3-size(unit_LM,2))]', ...
    stro.sum.exptParams.bkgndrgb, M, 'both');

% the presented contrast is pinned to the edge of the gamut; do that here too
presented_contrast = requested_contrast;
for ii = 1:nstim
    Loog_stim = stro.trial(:,Lstim_idx) == stim_idxs(ii) & stro.trial(:,Loog);
    presented_contrast(Loog_stim) = max_cc_scales(ii);
end

if show_plots
    figure();
    [L,M] = arrayfun(@(x) rat(x,5e-2), abs(unit_LM(:,1)./unit_LM(:,2)));
    L = L.*sign(unit_LM(:,1));
    M = M.*sign(unit_LM(:,2));
    ncols = ceil(sqrt(nstim));
    nrows = ceil(nstim/ncols);
    for ii = 1:nstim
        subplot(nrows,ncols,ii);
        semilogy(requested_contrast(stim_trial_positions(:,ii)), '-k'); hold on
        semilogy(presented_contrast(stim_trial_positions(:,ii)), '*k');
        title(sprintf('(#%d) %dL : %dM @ %5.2f Hz', stim_idxs(ii), L(ii), M(ii), tfs(ii)));
        if mod(ii-1, ncols) == 0 % 1st plot of every subplot row
            ylabel('cone contrast');
        end
        if ii > (nrows-1)*ncols % every subplot along the bottom row
            xlabel('stimulus presentation');
        else
            set(gca, 'xticklabel', []);
        end
    end
end
