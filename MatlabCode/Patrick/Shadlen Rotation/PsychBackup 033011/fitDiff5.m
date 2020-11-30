function err = fitDiff5(theta, data)

cohs   = data(:,1);
t1_obs = data(:,2);
t1_se  = data(:,3);
t2_obs = data(:,4);
t2_se  = data(:,5);


[t1_pred,t2_pred,p_pred]=calcDiff5(cohs, theta);

n1_obs = data(:,6);
n_total = data(:,7);

%n_pred = n_total.*p_pred;
%n_se = sqrt(n_total.*p_pred.*(1-p_pred));

N = size(data,1);

err = 0;

t1_ind = ~isnan(t1_pred) & ~isnan(t1_se);
t2_ind = ~isnan(t2_pred) & ~isnan(t2_se);

% Cost function of RT
err = err + sum((t1_obs(t1_ind) - t1_pred(t1_ind)).^2./(2.*t1_se(t1_ind).^2)) + N./2.*log(2.*pi) + sum(log(t1_se(t1_ind)));
err = err + sum((t2_obs(t2_ind) - t2_pred(t2_ind)).^2./(2.*t2_se(t2_ind).^2)) + N./2.*log(2.*pi) + sum(log(t2_se(t2_ind)));

% Cost function of choice probability
err = err - sum(gammaln(n_total+1) - gammaln(n1_obs+1) - gammaln(n_total-n1_obs+1)...
    + n1_obs.*log(p_pred*eps) + (n_total-n1_obs).*log(1-p_pred+eps));


%sum(gammaln(n+1) - gammaln(n_obs+1) - gammaln(n-n_obs+1) + n_obs.*log(ppred+eps) + (n-n_obs).*log(1-ppred+eps));

%sum(log(prod(1:n_total)./(prod(1:n_obs).*prod(1:(n_total-n_obs))))...
%    + n_obs.*log(p_pred)+(n_obs-n_pred).*log(1.0000001-p_pred));
