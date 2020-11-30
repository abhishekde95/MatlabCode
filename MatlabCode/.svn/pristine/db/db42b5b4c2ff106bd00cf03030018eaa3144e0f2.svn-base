function update_QUEST(correct)
global gl
k = QUEST_idx_from_condition(gl.cond_idx);

Qbeta = 4;
% QUEST's estimate may have been OOG. The slave presents a stimulus at the
% gamut's edge. We replicate that truncation here.
contrast_scale = min(gl.max_cc_scales(k), QUEST_estimate(k));
likelihood = 1-exp(-(contrast_scale./gl.Q.ranges{k}).^Qbeta)/2;
if ~correct
    likelihood = 1-likelihood;
end
gl.Q.funcs{k} = gl.Q.funcs{k}+log(likelihood);
