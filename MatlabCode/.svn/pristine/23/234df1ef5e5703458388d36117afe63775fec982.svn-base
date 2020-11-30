% This annonyingly terse function simply returns Normal PDFs given the previously
% computed lattice (the x values), mu, and sigma.

function Q = init_QUEST_with_priors(ranges, mu, sd)
Q = cellfun(@(x,m,s) {log(normpdf(x,m,s))}, ranges, num2cell(mu(:)), num2cell(sd(:)));
