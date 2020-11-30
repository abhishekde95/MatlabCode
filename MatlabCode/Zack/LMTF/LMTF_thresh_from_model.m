function r = LMTF_thresh_from_model(varargin)
%function r = LMTF_thresh_from_model(input)
%
% Pass in L, M, TF, model_params --OR-- LM_theta (whichever theta), TF,
% model_params --OR-- [L,M,TF],model_params.
% In the first case, L and M must have the same
% number of components and with the same dimensions. You can pass in an unequal
% number of TFs and color directions, but the non-singleton dimensions of both
% must be the same (i.e., you may need to transpose one of the inputs or flatten
% them. Here is the mode that describes detection threshold in each of 
% two orthogonal direction in the LM plane:

%       f1 = (1i*2*pi*tau1.*tf+1).^-n1;
%       f2 = (1i*2*pi*tau2.*tf+1).^-n2;
%       f = 1./abs(xi*(f1-zeta*f2));
%
% This function works on on the "local model", which consists of
% 13 parameters (6 for the first color direction, 6 for the second color
% direction, and theta, whichis the angle of the first color direction 
% in the LM plane). It does not accept the full 17 parameter model that
% includes a parametric description of how the parameters xi_lum, xi_rg, 
% and theta vary over retinal location. 
%
% GDLH 5/5/16

fpar = varargin{end};
if nargin == 2 % [L,M,TF],model
    LM_angle = atan2(varargin{1}(:,2), varargin{1}(:,1));
    TF = varargin{1}(:,3);
elseif nargin == 3 % [theta,TF],model
    LM_angle = varargin{1};
    TF = varargin{2};
elseif nargin == 4 % L, M, TF, model
    L = varargin{1};
    M = varargin{2};
    LM_angle = atan2(M, L);
    TF = varargin{3};
end
if length(LM_angle) == 1
   LM_angle = repmat(LM_angle,length(TF),1); 
end

% Order of model parameters
xi_1 = fpar(1);
zeta_1 = fpar(2);
n1_1 = fpar(3);
n2_1 = fpar(3)+fpar(4); % convention: n2 = n1+delta_n
tau1_1 = 10^fpar(5);
tau2_1 = 10^(fpar(5)+fpar(6)); % convention: tau2 = kappa*tau1
xi_2 = fpar(7);
zeta_2 = fpar(8);
n1_2 = fpar(9);
n2_2 = fpar(9)+fpar(10);
tau1_2 = 10^fpar(11);
tau2_2 = 10^(fpar(11)+fpar(12));
theta = fpar(13);

%f1 = @(omega)fpar(1)*abs(((1i*2*pi*fpar(5).*omega+1).^-fpar(3))-fpar(2)*((1i*2*pi*fpar(5)*fpar(6).*omega+1).^-fpar(4)));
%f2 = @(omega)fpar(1+6)*abs(((1i*2*pi*fpar(5+6).*omega+1).^-fpar(3+6))-fpar(2+6)*((1i*2*pi*fpar(5+6)*fpar(6+6).*omega+1).^-fpar(4+6)));
% Reparametrizing. 4/23/16 GDLH

f1 = @(omega)xi_1*abs(((1i*2*pi*tau1_1.*omega+1).^-n1_1)-zeta_1*((1i*2*pi*tau2_1.*omega+1).^-n2_1));
f2 = @(omega)xi_2*abs(((1i*2*pi*tau1_2.*omega+1).^-n1_2)-zeta_2*((1i*2*pi*tau2_2.*omega+1).^-n2_2));
a = abs(f1(TF)); % lum sensitivities
b = abs(f2(TF)); % rg sensitivies

% --- isodetection ellipse orientation varies with TF 
% See also tf_fiterr2.m

dotprods = zeros(numel(LM_angle),2);
dotprods(:,1) = bsxfun(@times,a(:),(cos(LM_angle(:)).*cos(theta))+(sin(LM_angle(:)).*sin(theta)));
dotprods(:,2) = bsxfun(@times,b(:),(cos(LM_angle(:))./sqrt(2))-(sin(LM_angle(:))./sqrt(2)));
r = sqrt(sum(dotprods.^2,2)).^-1;
r = reshape(r,size(TF));
