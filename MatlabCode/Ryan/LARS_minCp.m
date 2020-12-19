function [betaParse, minCpidx] = LARS_minCp(betaAll, Cp)
%%
% This function finds beta whose step number corresponds to the minimum 
% value of Cp and combines them all into an output matrix. It also returns 
% a vector with the step number chosen
%
% Inputs:
%   'betaAll' [nfeatures, nsteps, noffset] - Contains results of LARS at
%       every step.
%   'Cp' [nsteps] - Vector containing the confidence metric for each step.
%
% Outputs:
%   'betaParse' [nfeatures, noffset] - Contains beta coeffiecients
%       corresponding to the minimum Cp value
%   'minCpidx' [noffset] - Vector containing the step number corresponding
%       to the minimum Cp value.
%
%%
% Find index of Cp minimum for each offset
[~, minCpidx] = min(Cp,[], 1);

% Initialize output array
betaParse = zeros(size(betaAll,1), length(minCpidx));


for ii = 1:length(minCpidx)
    betaParse(:,ii) = betaAll(:,minCpidx(ii),ii); 
end