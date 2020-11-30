function [err, thresholds, preds] = tf_fiterr2(fittedparams, data, Loog)
% tf_fiterr2
% Function to fit a 2-D surface to LMTF data.
% Fitting log sensitivities (equivalent to fitting log thresholds)
%
% INPUT arguments: (params, data, LOOG) or (params, data)
%     params: xi, zeta, n1, delta n, tau1, kappa
%     data: an nx3 matrix, the columns of which are [L, M, TF]
%         or an nx4 matrix [L, M, TF, Loog]
%     Loog: a binary vector which is '1' for stimulus directions that
%     went out of gamut.
% In data(:,[1 2]), out of gamut points should indicate the edge of the
% monitor gamut. These points contribute to the errors of the fit if the
% predicted threshold is lower than the actual threshold.
%
% OUTPUT is the sum absolute log radial error.
%
% The model for contrast thresholds in a single color direction is given
% by:
%       f1 = (1i*2*pi*tau1.*tf+1).^-n1;
%       f2 = (1i*2*pi*tau2.*tf+1).^-n2;
%       f = 1./abs(xi*(f1-zeta*f2));
% where tau1 and tau2 are the time constants of two hypothesized component
% filters. n1 and n2 are the "order" of the filters (the number of
% 1st order, exponential filters that are cascaded to create one of the
% component filters).
% The above equations are used twice, independently, to fit temporal
% contrast thresholds in two orthogonal directions (given by data(:,1) and
% data(:,2), rotated by ). These two 1-D threshold functions (f and g) are combined to 
% predict thresholds in any stimulus direction by:
%     (f.*g)./sqrt((g.*cos(th)).^2+(f.*sin(th)).^2)
% where f and g are functions of temporal frequency and "th" is the angle
% in the LM plane.
% GDLH 4/10/14

if (nargin < 3)
    if (size(data,2) >= 4) % allowing Loog to be the fourth column of data
        Loog = data(:,4);
    else
        Loog = zeros(size(data,1),1);
        disp('Loog not passed in. Assuming all points were in gamut.');
    end
end
npars = length(fittedparams);
if npars == 6  % radially symmetric funnel - this was more of a test of the fitting procedure than a useful model
    xi_1 = fittedparams(1);
    zeta_1 = fittedparams(2);
    n1_1 = fittedparams(3);
    n2_1 = fittedparams(3)+fittedparams(4);
    tau1_1 = 10^fittedparams(5);
    tau2_1 = 10^(fittedparams(5)+fittedparams(6)); % convention: tau2 = kappa*tau1
elseif npars == 12 || npars == 13 % Two independent filters on orthogonal axes
    xi_1 = fittedparams(1);
    zeta_1 = fittedparams(2);
    n1_1 = fittedparams(3);
    n2_1 = fittedparams(3)+fittedparams(4);
    tau1_1 = 10^fittedparams(5);
    tau2_1 = 10^(fittedparams(5)+fittedparams(6));
    xi_2 = fittedparams(7);
    zeta_2 = fittedparams(8);
    n1_2 = fittedparams(9);
    n2_2 = fittedparams(9)+fittedparams(10);
    tau1_2 = 10^fittedparams(11);
    tau2_2 = 10^(fittedparams(11)+fittedparams(12)); % convention: tau2 = kappa*tau1
else
    keyboard
    error('Fittedparams is the wrong size');
end
if npars == 13 % fitting the rotation too
    rotangle = fittedparams(13);
else
    rotangle = pi/4;
end

tfs = data(:,3);

f1 = @(omega)xi_1*abs(((1i*2*pi*tau1_1.*omega+1).^-n1_1)-zeta_1*((1i*2*pi*tau2_1.*omega+1).^-n2_1));
f2 = @(omega)xi_2*abs(((1i*2*pi*tau1_2.*omega+1).^-n1_2)-zeta_2*((1i*2*pi*tau2_2.*omega+1).^-n2_2));
%a = 1./abs(f1(tfs)); % lum thresholds
%b = 1./abs(f2(tfs)); % rg thresholds
a = abs(f1(tfs)); % lum sensitivities
b = abs(f2(tfs)); % rg sensitivies

[th,r] = cart2pol(data(:,1),data(:,2)); % order of arguments: Lcc, Mcc
% Right now th and r are in LM coordinates.
% a and b are no longer the length of the two axes of the ellipse;
% they are gains on the two mechanisms.

%tic;
%f = zeros(size(r));
%mechs = [cos(rotangle) 1./sqrt(2) ; sin(rotangle) -1./sqrt(2)];
%for i = 1:length(f)
%    dotprods = [cos(th(i)) sin(th(i))]*mechs*diag([a(i) b(i)]);
%    f(i) = sqrt(sum(dotprods.^2))^-1; % converting sensitivies to thresholds
%end
%t0 = toc;
%
% % Alternative way of computing "f" without the loop. Speeds up execution.

% tic
dotprods = zeros(size(r,1),2);
dotprods(:,1) = bsxfun(@times,a,(cos(th).*cos(rotangle))+(sin(th).*sin(rotangle)));
dotprods(:,2) = bsxfun(@times,b,(cos(th)./sqrt(2))-(sin(th)./sqrt(2)));
f_alt = sqrt(sum(dotprods.^2,2)).^-1;
% t1 = toc;
% 
% if any(abs(f_alt-f)./min(f,f_alt)> 10^-14)
%     disp('algorithms for computing r do not agree!')
%     keyboard
% else
%     disp(['Speed up using new algorithm: ',num2str(t0-t1),' s'])
% end
f = f_alt;

% Uncomment line below to compare data and predictions online
% [r f]

logres = log10(r) - log10(f);
% Remove residuals if (1) the search went out of gamut (the mode of the QUEST function 
% was outside of the monitor gamut) and (2) the predicted threshold (f) is GREATER than the
% actual threshold (r).
logres(Loog & (log10(f) > log10(r))) = 0;
err = sum(abs(logres));
if (nargout > 1)
    thresholds = r;
    preds = f; 
end
%err = sum(logres.^2);

end
