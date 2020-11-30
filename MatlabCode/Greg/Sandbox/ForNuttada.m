% trying to figure out the SD of (a/(a+b))

mu1 = 20;
sigma1 = 1;
mu2 = 20;
sigma2 = 1;
n = 5;

a = normrnd(mu1, sigma1, n,1);
b = normrnd(mu2, sigma2, n,1);
sd_estimates = [];
%%
% approach 1: brute force
niter = 10000;
data = zeros(niter, 1);
for i = 1:niter
    aprime = normrnd(mu1, sigma1, n,1);
    bprime = normrnd(mu2, sigma2, n,1);
    data(i) = mean(aprime)/(mean(aprime)+mean(bprime));
end
figure; axes; hold on;
hist(data)
title(['Brute force simulated SD = ',num2str(std(data))]);
sd_estimates(1) = std(data);

%%
% approach 2: nonparametric bootstrap

niter = 2000;
data = zeros(niter, 1);
for i = 1:niter
    aprime = a(unidrnd(n, n, 1));
    bprime = b(unidrnd(n, n, 1));
    data(i) = mean(aprime)/(mean(aprime)+mean(bprime));
end
figure; axes; hold on;
hist(data)
title(['Non-parametric bootstrap SD = ',num2str(std(data))]);
sd_estimates(2) = std(data);

%%
% approach 3: parametric bootstrap

niter = 2000;
data = zeros(niter, 1);
for i = 1:niter
    aprime = normrnd(mean(a), std(a), n,1);
    bprime = normrnd(mean(b), std(b), n,1);
    data(i) = mean(aprime)/(mean(aprime)+mean(bprime));
end
figure; axes; hold on;
hist(data)
title(['Parametric bootstrap SD = ',num2str(std(data))]);
sd_estimates(3) = std(data);

%%
% approach 3: linear approximation

% d_da = partical derivative of a/(a+b) wrt to 'a'
d_da = mean(b)./(mean(a)+mean(b))^2;
d_db = -mean(a)./(mean(a)+mean(b))^2;
sr = sqrt((d_da)^2*(var(a)/n)+(d_db)^2*(var(b)/n))
disp(['SD of a/(a+b) by first order Taylor series expansion: ',num2str(sr)]);
sd_estimates(4) = sr;

% Another formula
sqrt((std(a)/mean(a))^2+(std(b)/mean(b))^2)


%%
% Summary
figure;
bar(sd_estimates);

