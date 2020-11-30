function hyp = default_hyperparameters()
% these defaults are the mean values from G, S, and Z.
% i.e., log10(mean(10.^[. . .]))
%hyp = struct('cov', [-0.75 .5 1 -1.5]', 'lik', -2.3, 'mean', -1);
hyp = struct('cov', [0 .2 .2 -1.8]', 'lik', -2.3, 'mean', [-1; .6; 0]); % GDLH 4/15/16 Trying a new set of default hyperparameters
verify_hyperparameters(hyp); % make sure there are no mistakes with the line above
% if there are issues, this function will halt the program
