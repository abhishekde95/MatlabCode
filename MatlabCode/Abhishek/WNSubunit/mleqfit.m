function lik = mleqfit(input,xdata,raw_data,spike_data)
% Maximum likelihood estimate of the 1-D non-linearity

% options.MaxIter = 10000000;
% options.MaxFunEvals = 1e7;
eps =  1e-5;
a = input(1);
b = input(2);
c = input(3);

% Assuming that the relation between the firing rate and the projection values has an
% exponential relation
prob = a*xdata.^2 + b*xdata + c;
prob(prob>=1-eps) = 1-eps;  % Ugly hack to avoid log10(0) = -Inf errors.
prob(prob<eps) = eps;

% Assuming the error values for a given x value is Bernoulli distributed
logP = log10([prob(:), 1-prob(:)]);
lik = spike_data'.*logP(:,1) + (raw_data-spike_data)'.*logP(:,2);
lik = (-1)*sum(lik);

end





