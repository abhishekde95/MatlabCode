% This function maps a condition index to the corresponding kth QUEST function.
% There are n QUEST functions, but cond_idx is in 1:2*n (each stimulus is
% presented to both hemifields, hence the *2). This mod trick works because the
% order of the stimuli in the conditions array is the same in positions 1:n and
% n+1:2*n (and the QUEST functions were initialized in the same 1:n order).

function k = QUEST_idx_from_condition(cond_idx)
global gl
k = mod(cond_idx-1, length(gl.Q.funcs))+1;
