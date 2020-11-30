% REX calls this function every trial to get the next stimulus. Aside from
% init_LMTF, this is the only function REX calls directly in all of LMTF.

function out = next_stimulus()
global gl

stim_idx = gl.conditions(gl.cond_idx,1);
hemifield = gl.conditions(gl.cond_idx,2);

% get contrast estimate from kth QUEST function
k = QUEST_idx_from_condition(gl.cond_idx);
contrast_scale = QUEST_estimate(k);

out = [contrast_scale*gl.stimuli(stim_idx,1:2) gl.stimuli(stim_idx,3) stim_idx hemifield];
