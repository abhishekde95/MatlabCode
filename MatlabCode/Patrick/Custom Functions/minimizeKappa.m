function [f] = minimizeKappa(kappa,mu,x)

idx = mu == 0;
% if sum(idx) > 0
%     keyboard
% end
sig = sqrt(mu + (kappa .* mu.^2));
p = ((sig.^2 - mu) ./ sig.^2);
r = (mu.^2 ./ (sig.^2-mu));
tempf = nan(size(mu));
tempf(idx) = 0;
r(isinf(r)) = 0;
p(isnan(p)) = 0;
try
    tempf(~idx) = -gammaln(r(~idx)+x(~idx)) + gammaln(r(~idx))...
        - (log(p(~idx)).*x(~idx)) - (log(1-p(~idx)).*r(~idx)); % from wikipedia
catch
    keyboard
end
f = sum(tempf);


end