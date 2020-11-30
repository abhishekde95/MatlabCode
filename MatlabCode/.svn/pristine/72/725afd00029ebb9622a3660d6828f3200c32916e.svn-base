% This function is insurance that any changes to the structure of the
% hyperparameters within the paradigm are reflected properly in the stro data we
% load in.

function verify_hyperparameters(hyp)
D = 2; % 2 input dimensions

assert(isstruct(hyp), 'This function expects a gpml-like hyperparameter structure');
hyp_vals = struct2cell(orderfields(hyp));
nhyp_vals = cellfun(@length, hyp_vals);
assert(length(nhyp_vals) == 3, 'This function expects hyperparameters for cov, lik, and mean only');

% make sure we have the correct number of hyperparameters expected by the
% specified covariance, likelihood, and mean functions (defined as MATLAB
% functions in private/{cov,lik,mean}func.m).
cov = covfunc();
lik = likfunc();
meanf = meanfunc();

ncovhyp = str2func(['@(D)' feval(cov{:})]);
assert(ncovhyp(D) == nhyp_vals(1), ...
    'The number of cov hyperparameters in the STRO doesn''t match the number required by covfunc.m');
nlikhyp = str2func(['@(D)' feval(lik{:})]);
assert(nlikhyp(D) == nhyp_vals(2), ...
    'The number of lik hyperparameters in the STRO doesn''t match the number required by likfunc.m');
nmeanhyp = str2func(['@(D)' feval(meanf{:})]);
assert(nmeanhyp(D) == nhyp_vals(3), ...
    'The number of mean hyperparameters in the STRO doesn''t match the number required by meanfunc.m');
