% [stims,torex] = LMTF(target_nstims, subject_params, varargin)
%
% INPUT
%      target_nstims: the number of stimuli requested (this is an upper
%          bound; the function probably returns fewer because some selected
%          stimuli are going out of gamut). 
%       subject_params: structure containing fields "model" which is 17x1
%          and "domain" which is 6x1, which are, in order: min(Lcc),
%          min(Mcc), min(TF), max(Lcc), max(Mcc), max(TF).
%
% OUTPUT
%       stims: nx4 matrix of stimuli to be presented in the format [Lcc,
%           Mcc, Scc, TF]
%       torex: The contents of subject_params.model
%
% varargin as of April 2016:
% Index | Contents
%     1 | gGabor.sf (spatial frequency, unused)
%     2 | gGabor.theta (orientation, unused)
%     3 | gGabor.tf (temporal frequency, unused)
%     4 | gPublic.r[0], gPublic.r[1], gPublic.th[0], gPublic.th[1],
%     gPublic.z[0], gPublic.z[1] (limits on r, theta, and z - "restrict")
%     5 | gPublic.rf_x, gPublic.rf_y (receptive field)
%     6 | gPublic.special_case
%     7 | gPublic.blanks (show a blank stimulus or not)
%
% Updating to the 18-parameter "tilted rampy trough" model:
%     zeta_lum = fpar(1)
%     n_lum = fpar(2)
%     delta_n_lum = fpar(3)
%     logtau_lum = fpar(4)
%     logkappa_lum = fpar(5)
%     zeta_rg = fpar(6)
%     n_rg = fpar(7)
%     delta_n_rg = fpar(8)
%     logtau_rg  = fpar(9)
%     logkappa_rg = fpar(10)
%     theta = fpar(11)
%     xi_lum = fpar(12)+fpar(13)*r+fpar(14)*phi+fpar(15)*phi^2;
%     xi_rg = fpar(16)+fpar(17)*r+fpar(18)*phi^2;
%
% GDLH 12/1/16
%
% GDLH 1/5/17: Updating tilted rampy trough model so that log10 gain is
% linear along the horizontal meridian.
%     log_xi_lum = fpar(12)+fpar(13)*r+fpar(14)*phi+fpar(15)*phi^2;
%     log_xi_rg = fpar(16)+fpar(17)*r+fpar(18)*phi^2;
%
% GDLH 4/16/17: Switching to the yoked tilted rampy trough model.


function [stims,torex] = LMTF(target_nstims, subject_params, varargin)
global gl
% load in parameters from the subject's struct

model = subject_params.legacy.mode5params;
domain = subject_params.domain; % Do I really need this? 

% use parameters from REX (check IsoSamp.d's GENSTIMS_REQ and genStimSet())
restrict = varargin{4};
rf = varargin{5};
special_case = varargin{6};
blank_stim = varargin{7};

[phi,r] = cart2pol(abs(rf(1))/10,rf(2)/10); % important to abs(x). This is the convention of LMTF and we don't want to square angles > pi/2
% GDLH 4/16/17 Yoked, tilted rampy trough
%log_xi_lum = model(12)+model(13)*r+model(14)*phi+model(15)*phi^2;
%log_xi_rg = model(16)+model(17)*r+model(14)*phi+model(18)*phi^2;
% GDLH 7/5/17 Updated yoked, tilted rampy trough
log_xi_lum = model(12)+model(13)*r+model(14)*r*sin(2*phi)+model(15)*r*cos(2*phi);
log_xi_rg = model(16)+model(17)*r+model(14)*r*sin(2*phi)+model(18)*r*cos(2*phi);

localmodel = [10.^log_xi_lum; model(1:5); 10.^log_xi_rg; model(6:11)];

domain = restrict_LMTF_domain(domain, localmodel, restrict); % Do I need this?
domain = reshape(domain ,3 ,2);
% Note: special cases 1-4 sample L+M and L-M, regardless of the theta
% parameter of the model fit.
switch special_case
    case 1 % Lum only
        L = cos(pi/4);
        M = sin(pi/4);
        TF = logspace(log10(domain(3,1)), log10(domain(3,2)), target_nstims)';
        r = LMTF_thresh_from_model(L, M, TF, localmodel);
        lmtf = [L*r, M*r, TF];
    case 2 % RG only
        L = cos(3*pi/4);
        M = sin(3*pi/4);
        TF = logspace(log10(domain(3,1)), log10(domain(3,2)), target_nstims)';
        r = LMTF_thresh_from_model(L, M, TF, localmodel);
        lmtf = [L*r, M*r, TF];
    case 3 % Lum & RG
        L1 = cos(pi/4);
        M1 = sin(pi/4);
        L2 = cos(3*pi/4);
        M2 = sin(3*pi/4);
        TF = logspace(log10(domain(3,1)), log10(domain(3,2)), ceil(target_nstims/2))';
        r1 = LMTF_thresh_from_model(L1, M1, TF, localmodel);
        r2 = LMTF_thresh_from_model(L2, M2, TF, localmodel);
        lmtf = [[L1*r1;L2*r2], [M1*r1;M2*r2], [TF;TF]];
    case 4 % Lum & RG both polarities
        L1 = cos(pi/4);
        M1 = sin(pi/4);
        L2 = cos(3*pi/4);
        M2 = sin(3*pi/4);
        TF = logspace(log10(domain(3,1)), log10(domain(3,2)), ceil(target_nstims/4))';
        r1 = LMTF_thresh_from_model(L1, M1, TF, localmodel);
        r2 = LMTF_thresh_from_model(L2, M2, TF, localmodel);
        lmtf = [[L1*r1;-L1*r1;L2*r2;-L2*r2], [M1*r1;-M1*r1;M2*r2;-M2*r2], [TF;TF;TF;TF]];
    otherwise
        lmtf = LMTF_subsample(target_nstims, localmodel, domain, false); 
end

% just scale r since we took care of theta and TF above
lmtf_restricted = restrict_samples(lmtf, [restrict(1:2) zeros(1,4)]);

in_gamut = gamut_extent(lmtf_restricted(:,1:2));
lmtf_restricted = lmtf_restricted(in_gamut,:);
% We're eliminating stimuli that are outside of the gamut which is why the
% number of stimuli requested is different from the number returned. Should
% be fixed.

% stims = [L M S=0 TF]
stims = [lmtf_restricted(:,1:2) zeros(size(lmtf_restricted,1),1) lmtf_restricted(:,3)];

if blank_stim, stims = [zeros(size(stims(1,:))); stims]; end

torex = localmodel(:)';