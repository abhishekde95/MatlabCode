function [model, hessmat, fv, logresids] = quickLMTFmodelfit(varargin)

% [model, hessian, fv, logresids] = quickLMTFmodelfit(L, M, TF, Loog) -OR-
% [model, hessian, fv, logresids] = quickLMTFmodelfit([L, M, TF, Loog])
%
% INPUT
%   (L, M, TF, Loog) either passed as four arguments or an n x 4 matrix
%
% OUTPUT
%   model: 13-element vector containing the parameters from the 2D Watson
%   hessmat: The Hessian matrix of the fitted parameters
%   fv: The minimized fitting error (summed squared log10 residuals)
%   logresids: log10 residuals
%
% Here's what the various model parameters actually mean:
% (parameters 1-6 apply to the LUM temporal filter)
%    1) xi_1: overall gain factor
%    2) zeta_1: gain factor for the subtractive component (filter #2)
%    3) n_1: number of cascaded low pass stages for the excitatory component
%    4) delta_n1: number of cascaded low pass stages for the subtractive
%    component (n1_1+delta_n1)
%    5) tau_1: time constant for the excitatory component
%    6) kappa_1: kappa_1*tau_1 is the time constant for the subtractive component
% (parameters 7-12 apply to the RG temporal filter)
%    7) xi_2: overall gain factor
%    8) zeta_2: gain factor for the subtractive component (filter #2)
%    9) n_2: number of cascaded low pass stages for the excitatory component
%    10) delta_n2: number of cascaded low pass stages for the subtractive component
%    11) tau1_2: time constant for the excitatory component
%    12) kappa_2: kappa_2*tau_2 is the time constant for the subtractive component
%    13) theta: angle (counter-clockwise from horizonal) of the LUM
%    mechanism
%
% Algorithm:
%     Fits two 1-D CRFs to data in a wedge spanning L+M and L-M and uses these
% as initial guesses to the 2-D fit, as implemented by fmincon calling tf_fiterr2.m
%
% GDLH 4/21/16
%
% NOTE: Change call to LMTF_thresh_from_model once convention of model
% parameter passing is changed (last parameter should be theta - don't
% treat this parameter as an entirely different animal).
PLOT1D = false; % <-- For diagnosing problems that might emerge from poor initial 1D fits
% Setting up constants
MINNPTS = 6; % minimum number of points to use for 1D fits
% Watson eq. 40 contains "n-1"s so "1" is a good lower bound for n1 and n2.
LB = [0   0  1   0 -2  log10(1.1)]; % if delta_n == 1, min(n1) = 2
UB = [500 1 40   2 -1  log10(2)];

UBtheta = pi/2;
LBtheta = 0;
options = optimset('Algorithm', 'interior-point', 'MaxFunEvals', 5e4, ...
    'MaxIter', 5e4, 'TolFun', 1e-8, 'Display', 'off');
initparams(1,:) = [20  .9  3 0 -2.5 log10(1.5)]; % LUM
initparams(2,:) = [100 .1 10 0 -2.5 log10(1)]; % RG

% Processing input
if (size(varargin{1},1) < MINNPTS)
   warning('Too few data points to fit model. Aborting.');
   model = nan; hessmat = nan; fv = nan; logresids = nan;
   return
end
if nargin == 1
    if size(varargin{1},2) ~= 4
        error('quickLMTFmodelfit expects (L,M,TF,LOOG) either as a matrix or 4 separate arguments');
    else
        Lcc = varargin{1}(:,1);
        Mcc = varargin{1}(:,2);
        tf = varargin{1}(:,3);
        Loog = varargin{1}(:,4);
    end
elseif nargin == 4
    Lcc = varargin{1};
    Mcc = varargin{2};
    tf = varargin{3};
    Loog = varargin{4};
else
    error('quickLMTFmodelfit expects (L,M,TF,LOOG) either as a matrix or 4 separate arguments');
end

[th,r] = cart2pol(Lcc, Mcc);% L,M
th = mod(th,pi);
angles = [pi/4 3*pi/4]; % For starting guesses
fpar = [];
L = logical([]);
for i = 1:2%
    err = abs(mod(th-angles(i),pi));
    L(:,i) = abs(err) < .17; % ~20° wedge on L+M
    if (sum(L(:,i)) < MINNPTS) % use the points inside the wedge if there are enough of them, otherwise take the 10 points closest to pi/4 (or 3*pi/4)
        disp('Too few points. Widening wedge');
        [~,erridx] = sort(err,1,'ascend');
        L(:,i) = err<=err(erridx(MINNPTS));
        % th(L) % sanity check
    end
    
    fpar(:,i) = fmincon(@(params) tf_fiterr(params, tf(L(:,i)), 1./r(L(:,i))), initparams(i,:), ...
        [], [], [], [], LB, UB, [], options);
end

if (PLOT1D)
    figure;
    for i = 1:2
        subplot(2,1,i); hold on;
        plot(tf(L(:,i)), 1./r(L(:,i)),'ko');
        f1 = @(omega)fpar(1,i)*abs(((1i*2*pi*10^fpar(5,i).*omega+1).^-fpar(3,i))-fpar(2,i)*((1i*2*pi*10^(fpar(5,i)+fpar(6,i)).*omega+1).^-(fpar(3,i)+fpar(4,i))));
        tfs = logspace(log10(1),log10(25),100);
        plot(tfs,f1(tfs),'k-');
    end
end

% Now doing the 2D fit, Keeping theta fixed at pi/4  <--- This step does
% essentially nothing.
%[model,~,~,~,~,~,~] = fmincon(@(params) tf_fiterr2(params,[Lcc,Mcc,tf],Loog), [fpar(:)' pi/4],...
%    [],[],[],[],[LB LB pi/4],[UB UB pi/4],[],options);

model = [fpar(:)' pi/4];
% Now doing the 2D fit, allowing theta to rotate
[model,fv,~,~,~,~,hessmat] = fmincon(@(params) tf_fiterr2(params,[Lcc,Mcc,tf],Loog), model,...
    [],[],[],[],[LB LB LBtheta],[UB UB UBtheta],[],options);

if (nargout == 4)
    pred_r = LMTF_thresh_from_model(Lcc,Mcc,tf,model);
    logresids = log10(r) - log10(pred_r);
end
