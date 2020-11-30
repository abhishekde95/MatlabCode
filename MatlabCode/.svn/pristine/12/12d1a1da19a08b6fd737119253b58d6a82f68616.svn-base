function [stims,torex] = DTNT(target_nstims, subject_params, varargin)
% load in parameters from the subject's struct
models = subject_params.models;
domains = subject_params.domain;
sfs = subject_params.sfs;

% use parameters from REX (check IsoSamp.d's GENSTIMS_REQ)
sf = varargin{1};
tf = varargin{3};
restrict = varargin{4};

[~,sfidx] = min(abs(sf-sfs));
model = models(:,sfidx);
domain = domains(:,:,sfidx);

lms = DTNT_subsample(target_nstims, model, domain, false);
lms_restricted = restrict_samples(lms, restrict);

prop_retained = size(lms_restricted, 1) / size(lms, 1);
if ~softEq(prop_retained, 1, 5)
    lms = DTNT_subsample(ceil(target_nstims/prop_retained), model, domain, false);
    lms_restricted = restrict_samples(lms, restrict);
end

in_gamut = gamut_extent(lms_restricted / 100); % need / 100 until we fit models using cc not %cc
lms_restricted = lms_restricted(in_gamut,:) / 100; % THANKS CHARLIE

stims = [lms_restricted tf+zeros(size(lms_restricted,1),1)];
torex = reshape(model,1,[]);
