% Section 1
% Code to compute the count distribution for a renewal process with
% Gamma-distributed intervals.
%
% Section 2
% A test of the formula to compute bounds on r23 given r12 and r23.
%
% Section 3
% Playing around with Cauchy random variables and confidence intervals for
% the ratio of two independent normals.
%
% Section 4
% Playing around with Fisher information. Trying to estimate a Poisson
% parameter.
%
% Section 5
% Simulating the Benjamini and Hochberg procedure to control the false
% discovery rate.
%
% Section 6
% Making a prediction ellipse for a bivariat Gaussian distribution.
% 
% Section 7
% Mean and covariance of a Poisson spike train
% consisting of 'n' spikes projected onto two (not necessarily orthogonal)
% basis vectors.
%
% Section 7.1
% Mean and covariance of a non-homogeneous (sinusoidal) Poisson spike trains
% onto sine and cosine basis vectors. Simulation and close-form.
%
% Section 8
% Playing around the James-Stein estimator
%
% Section 9
% trying to reconstruct Paul Glimcher's model of 2AFC performance.
%
% Section 10
% How are threshold, d' (at a constant contrast), and viewing duration
% related for an ideal observer?
%
% Section 11
% Working on ideal observer (matched filter) for correlated inputs.
%
% Section 12
% As above, but calculating "population d'" for an ideal observer of
% correlated inputs where correlation is proportional to overlap of
% circular RFs on a hexagonal grid.
%
% Section 13 Hexagonal grid of (truncated) Gaussian RFs and quantifiction
% of overlap
%
% Section 14 Experimenting with the Purpura and Victor spike distance
% stuff.
%
% Section 15 Trying to fit an ROC curve by maximum likelihood.
%
% Section 16 Comparing difference in medians to median of differences

%%
% Section 1
%
% Code to compute the count distribution for a renewal process with
% Gamma-distributed intervals.

% User-defined parmeters:
t = 50;  % Counting window length
a = 5; % First parameter of Gamma distn. (mean of exponentials).
b = 5; % Second parameter of Gamma (number of exponentials convolved).
maxspikes = 5*t/(a*b); 
% Largest number of spikes that we're going to consider as possible in a
% window of length 't'.  No just a big number - no principled reason for
% this choice.

%----------------------------
% Computing 1-CDF of count distribution using fact that:
% P(N(t)>=k) = P(Wk<=t)
% where N(t) is count in a window of length t
% and Wk is waiting time to kth renewal (in this case gamma with parameter
% k*b.
x = zeros(maxspikes,1);
for k = 1:maxspikes
    x(k) = gamcdf(t,k*b,a);
end
x = [1;x];  % P(N(t) >= 0) = 1
% x(k) is now P(N(t) > k)
countcdf = 1-x;
countpdf = diff(countcdf);
bar(countpdf);  % here it is

% As a sanity check we comparing this
% distribution against the Poisson (set b = 1 for Poisson process)
if (b == 1)
    pdf = poisspdf([0:maxspikes],t/a);
    hold on;
    plot(pdf,'y*');
end

%%
% Section 2

% Assume three random variables.  If we fix the correlation between
% variables 1 and 2 and we fix the correlation between variables 1 and 3,
% what can we say about the allowable correlations between variables 2 and
% 3?

r12 = -.9;
r13 = 0;
% Here are the bounds on r23 (the correlation between rvs 2 and 3)
bnds = [cos(acos(r12)-acos(r13)) cos(acos(r12)+acos(r13))];

% Now let's test it to make sure i haven't gotten it wrong.
% Building a series of correlation matrices with the values of r12 and r13
% given above.  Then, trying many different values for r23 and seeing if we
% get a legitimate correlation matrix (postive semi-definite).
cormat = eye(3);
cormat(1,2) = r12; cormat(2,1) = r12;
cormat(1,3) = r13; cormat(3,1) = r13;

testrs = linspace(-1,1,100);
success = nan*ones(1,length(testrs));
for i = 1:length(testrs)
    cormat(2,3) = testrs(i); cormat(3,2) = testrs(i);
    success(i) = det(cormat)>=0;
end

figure; axes; hold on;
plot(testrs,success);
plot(bnds,[1 1],'m*');

%%
% Section 3
% Playing around with Cauchy random variables and confidence intervals for
% the ratio of two independent normals.

muX = 10;
muY = 40;
sigmaX = 10;
sigmaY = 4;

% First, a Monte Carlo simulation
niter = 10000;
X = normrnd(muX, sigmaX, niter, 1);
Y = normrnd(muY, sigmaY, niter, 1);
Z = X./Y;
figure;
[n,x] = hist(Z,100);
bar(x,n./max(n));
hold on;

% Now using the formulas from Kamerund 1978
z = linspace(min(Z),max(Z),1000);
w = (sigmaY/sigmaX)*z;  % w is a scaled version of z
s = 1./sqrt(w.^2+1);
k = (muX/sigmaX.*w+muY/sigmaY).*s.^2;
M = -.5*(muY/sigmaY.*w-muX/sigmaX).^2.*s.^2;
Q = k.*s.*sqrt(2*pi).*(1-2.*normcdf(-k./s))+(2*s.^2.*exp(-k.^2./(2*s.^2)));

fg = 1/(2*pi).*Q.*exp(M);
fz = (sigmaY/sigmaX)*fg;
fz = fz./sum(fz);
plot(z,fz./max(fz),'m-');

% CDF
Fz = cumsum(fz);
CI = interp1(Fz,z,[.025 .975])


%%
% Section 4
% Playing around with the Fisher information associated with a Poisson
% parameter.

lambda = 2.1;
x = [0:lambda*10];
y1 = poisspdf(x,lambda);
%deltalambda = 1.0000001;
%y2 = poisspdf(x,lambda+log(deltalambda));
%dlogy = ((log(y1)-log(y2))./log(deltalambda)).^2;

deltalambda = .0000001;
y2 = poisspdf(x,lambda+deltalambda);
dlogy = ((log(y1)-log(y2))./deltalambda).^2;

figure; axes; hold on;
plot(x,y1-y2);
%plot(x,y1,'k-');
%plot(x,y2,'r-');

plot(x,dlogy,'k.');

% Expected value not using MC sampling
EJ = y1*dlogy';
1./EJ
lambda

%%
% Section 5 
% The Benjamini and Hochberg procedure
m = 100; % Total number of hypotheses tested
mo = 95; % true null hypotheses
n = 5; % number of samples in a t-test
mu = 1; % Under Ha
sigma = 1;
niter = 1000;
alpha = 0.05;

assert(m >= mo);

data = zeros(niter,m);
for i = 1:niter
   for j = 1:m
       x = normrnd(0,sigma,n,1);
       if j > mo
           x = x+mu; % alternative hypothesis is true
       end
       [h,p] = ttest(x);
       data(i,j) = p;
   end
end

q_star = alpha;
% calculating q values
proportion_false_discoveries = [];
proportion_misses = [];
for i = 1:niter
   % How many true null hypotheses did we reject?
   v = sum(data(i,[1:mo]) <= alpha); % ideally v = 0
   s = sum(data(i,[mo+1:end]) <= alpha); % ideally s = m-mo+1
   % How many alternative hypotheses hypotheses did we fail to reject?
   u = sum(data(i,[1:mo]) > alpha);
   t = sum(data(i,[mo+1:end]) > alpha);
   
   [ps,idx] = sort(data(i,:));
   % Note: data(i,idx) == ps
   L = ps <= ([1:m]./m)*q_star;
   hypotheses_rejected = idx(L);
   if (sum(hypotheses_rejected) == 0)
       FDR = 0;
   else 
       FDR = sum(hypotheses_rejected <= mo)/length(hypotheses_rejected);
   end
   hypotheses_accepted = idx(~L);
   MISSRATE = sum(hypotheses_accepted > mo)/length(hypotheses_accepted);
   
   proportion_false_discoveries(i,:) = [FDR v./(v+s)];
   proportion_misses(i,:) = [MISSRATE t/(u+t)];
end
% In what fraction of "significant" effects was the null hypothesis true?
figure;
subplot(2,2,1);
hist(proportion_false_discoveries(:,1))
title(num2str(mean(proportion_false_discoveries(:,1))));
subplot(2,2,2);
hist(proportion_false_discoveries(:,2))
title(num2str(mean(proportion_false_discoveries(:,2))));
% In what fraction of "non-sigificant effects" was the alternative
% hypothesis true?
subplot(2,2,3);
hist(proportion_misses(:,1))
title(num2str(mean(proportion_misses(:,1))));
subplot(2,2,4);
hist(proportion_misses(:,2))
title(num2str(mean(proportion_misses(:,2))));

% Output is a 2x2 figure
% Upper left: Proportion of "significant" effects (based on q-value) for
% which the null hypothesis was true. (proportion type 1 errors)
% Upper right: Proportion of "significant" effects (based p-value) for
% which the null hypothesis was true. (proportion type 1 errors)
% Lower left: Proportion of "non-significant" effects (based on p-value) for
% which the null hypothesis was false. (proportion type II errors)
% Lower right:  Proportion of "non-significant" effects (based on q-value) for
% which the null hypothesis was false. (proportion type II errors)

%%
% Section 6
% Prediction ellipse, T^2 test, 
% first two moments of Poisson spike train projected onto two basis vectors
% (useful for IsoSamp analyses).

% Generating fake data
mu = [3 5]';
S2 = [2 1; 1 3];
n = 100;
r = mvnrnd(mu,S2,n);
figure; axes; hold on;
plot(r(:,1),r(:,2),'bo');
plot(mu(1),mu(2),'r*')
plot(mean(r(:,1)),mean(r(:,2)),'b*');
S2hat = cov(r);
tmp = [cos(linspace(0,2*pi,100))' sin(linspace(0,2*pi,100))'];

% Theoretical prediction ellipse
[v,d] = eig(S2);
crit = sqrt(chi2inv(.95,2));
xy = tmp*crit*sqrt(d)*v;
plot(xy(:,1)+mu(1),xy(:,2)+mu(2),'r:')

% Empirical prediction ellipse
[v,d] = eig(S2hat);
%crit = sqrt(chi2inv(.95,2));
crit = sqrt((2*(n-1)/(n-2))*finv(.95,2,n-2)); % <--- Hotelling's T^2! (since we're estimating S2)
xy = tmp*crit*sqrt(d)*v;
plot(xy(:,1)+mean(r(:,1)),xy(:,2)+mean(r(:,2)),'b:');

% How many points are outside of the prediction ellipse?
% First theoretical ellipse
nreject = [0 0]; % Theoretical, Empirical
for i = 1:n
    T2 = (r(i,:)-mu')*inv(S2)*(r(i,:)-mu')';
    crit = chi2inv(.95,2);
    if T2 > crit
       plot(r(i,1),r(i,2),'y*');
       nreject(1) = nreject(1)+1;
    end
    
    T2 = (r(i,:)-mean(r))*inv(S2hat)*(r(i,:)-mean(r))';
    crit = (2*(n-1)/(n-2))*finv(.95,2,n-2);
    if T2 > crit
       plot(r(i,1),r(i,2),'bo');
       nreject(2) = nreject(2)+1;
    end
end
title(num2str(nreject));
%%
% Section 7
% Distribition of Poisson spike trains projected onto basis vectors.
% Needed for IsoSamp analysis.
% 8/2/19 I made a mistake in the original IsoSamp submission. 
% Improved formulas are below (and confimatory simulations)
nspikes = 3; % sp
omega = 1; % Hz
t = .666; % s
niter = 50000;

% simulation
spiketimes = unifrnd(0,t,nspikes,niter);
projs = [sum(cos(2.*pi.*omega.*spiketimes))' sum(sin(2.*pi.*omega.*spiketimes))'];

% Theoretical means
mux = nspikes*sin(2.*pi.*omega.*t)/(2.*pi.*omega.*t);
muy = nspikes*(1-cos(2.*pi.*omega.*t))/(2.*pi.*omega.*t);

disp('Means')
[mean(projs);[mux muy]]

disp('Empirical cov')
cov(projs)

disp('Theoretical cov')
% Variance of the cos term
% Z1...Zn ~ Unif(0,t)
% Yi = cos(2*pi*omega*Zi)
% X = sum(Yi)

% X = sum(cos(2*pi*omega*Z)) % Where the sum is over nspikes (assumed known)
% E(X^2) = E(sum(cos(2*pi*omega*Zi)^2))% Here the sum is across n identical terms
%        + E(sum(cos(2*pi*omega*Zi)*cos(2*pi*omega*Zj))) 
% where Zi and Zj are different spikes and the sum is over n(n-1) identical terms

% Term 1: 
% E(sum(cos(2*pi*omega*Zi)^2))
% nspikes*E(cos(2*pi*omega*Z)^2) % where Z ~ Unif(0,t) (all Zi's are IID)
% E(g(X)) = integral(g(z)f(z)dz) where g(z) = cos(2*pi*omega*z)^2 and f(z) = (1/t)
% nspikes*[(sin(4*pi*omega*t)/(8*pi*omega*t)+(1/2)] <--- final answer

% Term 2:
% E(sum(cos(2*pi*omega*Zi)*cos(2*pi*omega*Zj))
% n(n-1)*integrate(cos(2*pi*omega*z1)*cos(2*pi*omega*z2)/t^2) dz1 dz2)
% Above: 1/t^2 is the probability density of f(z1,z2) where z1, z2 iid unif(0,t)
% n(n-1)/t^2*integrate(cos(2*pi*omega*z1) dz1)*integrate(cos(2*pi*omega*z2) dz2)
% Above: breaking the double integral into two pieces (integrating along z1
% and then along z2)
% n(n-1)/t^2*((sin(2*pi*omega*t)/(2*omega*pi))^2


% Variance of the cos term
term1 = (nspikes/8)*(sin(4*pi*omega*t)/(pi*omega*t)+4);
% term 1 is the first term of E(X^2). I can't ignore the cross terms
% (term 2)
term2 = (nspikes*(nspikes-1)/t^2)*sin(2*pi*omega*t)^2/(2*pi*omega)^2;
a = term1+term2-mux^2;

% Covariance between sin and cos
term1 = nspikes*(sin(2*pi*omega*t)^2/(4*pi*omega*t));
term2 = (nspikes*(nspikes-1)/t^2)*sin(pi*omega*t)^3*cos(pi*omega*t)/(pi*omega)^2;
b = term1+term2-mux*muy;

% Variance of the sine term
term1 = (nspikes/t)*((t/2)-(sin(4*pi*omega*t)/(8*pi*omega)));
term2 = (nspikes*(nspikes-1)/t^2)*sin(pi*omega*t)^4/(pi*omega)^2;
c = term1+term2-muy^2;

[a b; b c]

%% Trying a sinusoidal modulation of spike rate (by warping time)
% Section 7.1

% Making a transformation that maps a raised cosine to a uniform function of the
% same integral.
dur = .66;
nspikes = 20;
ntrials = 1000;
freq = 1.75;
PLOTRASTERS = 1;

data = [];
phis = linspace(0,2*pi,5);
phis(end) = [];
for i = 1:length(phis)
    phi = phis(i);
    % Lame (discrete) way for now
    x = linspace(0,dur,10000); % for the interpolation
    y = 0.5*cos(2.*pi.*freq.*x+phi)+0.5;
    y_cum = (sin(2.*pi.*freq.*x+phi)+(2.*pi.*freq.*x)-sin(phi))./(4*pi*freq); % Not quite right
    y_cum = y_cum*dur/y_cum(end); % hack?
    %plot(x,y_cum)
    warpedspiketimes = interp1(y_cum,x,unifrnd(0,dur,nspikes,ntrials),'spline'); % inverting the integral of 1+cos(x)
    if PLOTRASTERS
        figure; axes; hold on;
        plot(warpedspiketimes,repmat(1:ntrials,size(warpedspiketimes,1),1),'k.');
    end
    % Empirical
    projs = [sum(cos(2.*pi.*freq.*warpedspiketimes))' sum(sin(2.*pi.*freq.*warpedspiketimes))'];
    
    empiricalS2 = cov(projs);
    data = [data; empiricalS2(:)' mean(projs)];
end

% Here's what the theory says:
mu_noise = nspikes*[sin(2*pi*freq*dur)+0 -cos(2*pi*freq*dur)+1]/(2*pi*freq*dur);
term1 = (nspikes/8)*(sin(4*pi*freq*dur)/(pi*freq*dur)+4);
term2 = (nspikes*(nspikes-1)/dur^2)*sin(2*pi*freq*dur)^2/(2*pi*freq)^2;
crossprods(1,1) = term1+term2-mu_noise(1)^2;
term1 = (nspikes*sin(2*pi*freq*dur)^2/(4*pi*freq*dur));
term2 = (nspikes*(nspikes-1)/dur^2)*sin(pi*freq*dur)^3*cos(pi*freq*dur)/(pi*freq)^2;
crossprods(2,1) = term1+term2-prod(mu_noise);
crossprods(1,2) = crossprods(2,1);
term1 = (nspikes/dur)*((dur/2)-(sin(4*pi*freq*dur)/(8*pi*freq)));
term2 = (nspikes*(nspikes-1)/dur^2)*sin(pi*freq*dur)^4/(pi*freq)^2;
crossprods(2,2) = term1+term2-mu_noise(2)^2;
S2 = crossprods;
% 
% % Looking at the covariance
% figure; axes; hold on;
% plot(phis, data(:,1),'k-'); % var x
% plot(phis, data(:,3),'b-'); % cov(x,y)
% plot(phis, data(:,4),'r-'); % var y
% plot([phis(1) phis(end)],S2(1,1)*[1 1],'k-');
% plot([phis(1) phis(end)],S2(1,2)*[1 1],'b--');
% plot([phis(1) phis(end)],S2(2,2)*[1 1],'r--');

% simulated variance is consistently lower than expected under the null
% hypothesis of no modulation.

% Ploting mean + 1 SD contours in x,y plane
figure; axes; hold on;
tmp = linspace(0,2*pi,100);  
plot(data(:,5),data(:,6),'bo');

for i = 1:length(phis)
   S2tmp = reshape(data(i,[1:4]),2,2);
   mntmp = data(i,[5 6]);
   error_ellipse(S2tmp,mntmp,'conf',.5,'style','blue');
end

% Theory assuming unmodulated firing rate
plot(mu_noise(1),mu_noise(2),'k+');
error_ellipse(S2,mu_noise,'conf',.5,'style','black');

axis equal;
%%
% Q: Can I represent time-shifted sinusoids from linear combinations of
% sin and cos basis vectors that modulate more slowly than the stimulus
% duration?
% A: Yes.

tf = 1;
dur = .66; % duration of stimulus
basisvects = [cos(linspace(0,2*pi*tf*dur,nbins+1))' sin(linspace(0,2*pi*tf*dur,nbins+1))']; 
basisvects(end,:) = [];

phi = linspace(0,2*pi,20); % phase shift
figure; axes; hold on;
for i = 1:length(phi)
   stim = cos(linspace(0,2*pi*tf*dur,nbins+1)+phi(i))';
   stim(end) = [];
   b = regress(stim, basisvects)
   pred = basisvects*b;
   resid = pred-stim;
   plot(resid);
end

%% 
% How big is the correction for non-orthogonal basis vectors?
% Not negligible and not monotonically decreasing with TF. A good idea to
% keep it in.

dur = .666;
out = [];
nspikes = 24;
for tf = logspace(log10(1),log10(50),500)
    mu_cont = nspikes*[sin(2*pi*tf*dur)+0 -cos(2*pi*tf*dur)+1]/(2*pi*tf*dur); % Integral of cos(ax) is sin(ax)/a. Integral of sin(ax) is -cos(ax)/a
    cp_cont(1,1) = (1/dur)*((dur/2)+sin(4*pi*tf*dur)/(8*pi*tf)) - (mu_cont(1)/nspikes)^2;
    cp_cont(1,2) = (1/dur)*(-1*cos(4*pi*tf*dur)+1)/(8*pi*tf)-mu_cont(1)*mu_cont(2)/nspikes^2;
    cp_cont(2,1) = cp_cont(1,2);
    cp_cont(2,2) = (1/dur)*((dur/2)-sin(4*pi*tf*dur)/(8*pi*tf)) - (mu_cont(2)/nspikes)^2;
    S2 = cp_cont*nspikes;
    out = cat(3,out,S2);
end

out1 = [];
for i = 1:size(out,3)
    %det(squeeze(out(:,:,i)));
    %squeeze(out(1,2,i))
    [v,d] = eig(squeeze(out(:,:,i)));
    out1 = [out1; diag(d)'];
end

figure; axes; hold on;
plot(out1(:,2)./out1(:,1),'k.-');
ylabel('Ratio of eigenvalues');
xlabel('TF (Hz)');

%%
% Analytically looking at the covariance between the sin and cos
% projections to see how it depends on the counting window
t = .66;
omega = linspace(1,6,100);
mux = sin(2.*pi.*omega.*t)/(2.*pi.*omega.*t); % Ignoring 'n' which cancels in caluclation of 'b' (the covariance btwn sin and cos)
muy = (1+cos(2.*pi.*omega.*t))/(2.*pi.*omega.*t); % Ignoring 'n' which cancels in caluclation of 'b' (the covariance btwn sin and cos)
b = (t/2+(1-cos(4.*pi.*omega.*t)))./(8.*pi.*omega.*t)-(mux*muy);

plot(omega,b)

%%
% Section 8
% Playing around with the James-Stein estimator
sigma2 = 10000;
mu = [10 10 10 10]';
d = length(mu);
niter = 200;
data =[];
for i = 1:niter
    tmp = normrnd(mu,sqrt(sigma2)*ones(d,1));
    est_ls = tmp;
    shrinkage_factor = 1-(d-2)*sigma2/(tmp'*tmp);
    est_js = shrinkage_factor*tmp;
    data = [data; sum((est_ls-mu).^2) sum((est_js-mu).^2)];
end

figure; axes; hold on;
plot(data(:,1),data(:,2),'.');
[minval maxval] = bounds(data);
plot([minval maxval],[minval maxval],'k:');
xlabel('SSE (LS estimator)');
ylabel('SSE (JS estimator)');

%%
% Section 9
% Trying to reconstruct Paul Glimcher's argument for why the area under the
% ROC is not the percent correct in a 2AFC task.

Rmu = 1.27;
Lmu = 0;
sigma = 1;
x = linspace(-4,5,5000);
yr = normpdf(x,Rmu,sigma);
yl = normpdf(x,Lmu,sigma);

figure; axes; hold on;
plot(x,yr);
plot(x,yl);

% numerical 
normcdf(x,Lmu,sigma)*(normpdf(x,Rmu,sigma)./sum(normpdf(x,Rmu,sigma)))'
% Brute force simulation
niter=1000;
count = 0;
for i = 1:niter
    if normrnd(Rmu,sigma) > normrnd(Lmu,sigma)
        count = count+1;
    end
end
count/niter
    
% Trying Paul's 2-D version
[xx,yy] = meshgrid(x,x);
gaussian1 = normpdf(x,Rmu,sigma)'*normpdf(x,Lmu,sigma);
gaussian2 = normpdf(x,Lmu,sigma)'*normpdf(x,Rmu,sigma);
figure; axes; hold on;
imagesc(gaussian1);
plot([0 length(x)],[0 length(x)]);
% According to Paul's model you get a draw either from the 1st bivariate
% distribution and compare it to a criterion that is the main diagonal
% Because of symmetry I can just consider one of the two distribitions 
% to get the error rate
threshold = xx<yy;
sum(sum(gaussian1./sum(gaussian1(:)).*threshold)) % quantization error?
% Another way of getting at this number is to look at a 1-D gaussian
% centered on 1/sqrt(2) with a threshold at 0
1-normcdf(0,1/sqrt(2),sigma)
% Looks like d' = 76% correct either way

% We can change Rmu to 1.27 (instead of 1) and get the right answer (81.5%
% correct) using Paul's method.
1-normcdf(0,Rmu/sqrt(2),sigma)

%%
% Section 10
% How is viewing duration related to d' and threshold of an ideal observer?

noise_mu = 0;
sigma= 1; % STD in a single time step
durations = [1 5 10]; % duration
signal_mu = logspace(-2,1,40); % contrast
dprimes = []; % dprimes
prctcor = []; % 2AFC percent correct
for d = 1:length(durations)
    for m = 1:length(signal_mu)
        dprimes(d,m) = sqrt(durations(d))*(signal_mu(m)-noise_mu)./sigma;
        tmp = linspace(noise_mu-(4*sigma*sqrt(durations(d))),signal_mu(m)+(4*sigma*sqrt(durations(d))),10000);
        binwidth = tmp(2)-tmp(1);
        prctcor(d,m) = binwidth*sum(normpdf(tmp,noise_mu,sigma./sqrt(durations(d))).*(1-normcdf(tmp,signal_mu(m),sigma./sqrt(durations(d)))))
    end
end
% Getting thresholds from psychometric functions
thresh = [];
for d = 1:length(durations)
    thresh(d) = interp1(prctcor(d,:),signal_mu,1-0.5*exp(-1));
end
% rows are durations, columns are contrasts

% psychometric functions as a function of viewing duration
figure;
plot(signal_mu,prctcor');
set(gca,'Xscale','log');
xlabel('contrast');
ylabel('% correct');
%looks like they translate to the left in log contrasts

% Plotting dprime as a function of viewing duration
figure;
subplot(2,2,1);
plot(signal_mu,dprimes');
set(gca,'Xscale','linear');
set(gca,'Yscale','linear');
xlabel('contrast');
ylabel('d prime');
set(gca,'Yscale','log','Xscale','log');

% There's a linear relationship between contrast and dprime. The slope
% depends on the viewing duration. 
% For a particular contrast, as the duration grows, dprime grows
% but at a diminishing rate. This diminishing rate shows up as a shallow
% slope (0.5) in the duration-d-prime relationship.

% The relationship between viewing duration and d' is linear in log-log.
subplot(2,2,2);
plot(durations,dprimes,'.-');
set(gca,'Yscale','log','Xscale','log');
xlabel('duration');
ylabel('d prime');

subplot(2,2,3);
plot(durations,thresh,'.-');
set(gca,'Yscale','log','Xscale','log');
xlabel('duration');
ylabel('threshold');

subplot(2,2,4); hold
plot(thresh,dprimes,'.-');
set(gca,'Yscale','log','Xscale','log');
xlabel('threshold');
ylabel('dprime');
% Each line is a different stimulus duration  

% So log threshold (or equivalently, log sensitivity) is linearly related
% to log d-prime.

(log10(thresh(2))-log10(thresh(1)))/(log10(durations(2))-log10(durations(1)))
(log10(thresh(3))-log10(thresh(1)))/(log10(durations(3))-log10(durations(1)))
% Slope is -.5 on log-log plot. Good. Doubling the stimulus duration
% decreases the threshold by a factor of sqrt(2).

%%
% Section 11: Observer of correlated inputs
% Weighting each neuron in proportion to it's mean, but not
% using the covariance structure across the neurons to make an
% even more ideal ideal-observer.

niter = 10000;
mu = [1;.2; .5]; % multiple signal sources
s2 = 1; % This can stay fixed at '1'
c = 0.1; 
S2 = toeplitz([s2 c c],[s2 c c]);
[v,d] = eig(S2);

% Generating fake data
signals = v*sqrt(d)*normrnd(zeros(size(mu,1),niter),ones(size(mu,1),niter));
signals = (signals+repmat(mu,1,niter))';
%figure; axes; 
%plot(signals(:,1),signals(:,2),'.')
%title(num2str(cov(signals)));

% incorrect calculation of d-prime
%dprime = mean()/sqrt(s2)

% Adjusting the weights of the ideal observer
%w = inv(S2)*mu; % the adjusted weights
%t = signals*w;
%dprime = mean(t./std(t))

% The ideal observer of correlated inputs can perform very well even if
% only one of the inputs is carrying signal because it can subtract off an
% independent estimate of the noise (provided by a correlated channel with
% no signal. This is not what I'm looking for. I want to just add up signal
% and noise of a bunch of correlated neurons

t = signals*mu; % Weighting each neuron by it's SNR (it's mu)
mean(t./std(t))
%mean(t./sqrt(s2*ones(size(mu'))*mu))

% Weighting each neuron by it's SNR (it's mu)
tmp = (mu*mu').*S2;
sum(mu.^2)./sqrt(sum(tmp(:)))

% Great, now I just need to calculate a covariance matrix for the array of
% LGN neurons
% Formula for the overlap between pairs of circles
% from http://jwilson.coe.uga.edu/EMAT6680Su12/Carreras/EMAT6690/Essay2/essay2.html
% d = distance between circle centers
% r = circle radius (assumed same for both circles)
r = .2; d = .4;
A = 2*r.^2.*acos(d./(2*r))-(d./2).*sqrt(4*r^2-d.^2)
%%
% Section 12
% Now trying to apply this to a hexagonal grid of circular RFs
% The covariance is proportional to degree of overlap: if Y1=X1+X2 and
% Y2=X2+X3) then Var(Y1)=Var(X1)+Var(X2) and Cov(Y1,Y2) = Var(X2).
INDIVIDUALPLOTS = true;
INFINITELYLARGESTIMULUS = false; % Every cell gets mu=1. Shows that number of cells really matters. Not useful.
sigma_gabor = 0.15; % DVA
sigmas_n = 3;
data = [];
RFcenterdistances = linspace(.12,.24,10);
%RFcenterdistances = [.12 .16 .18 .2]
RFdiams = .24;
for RFcenterdistances_idx = 1:length(RFcenterdistances)
    for RFdiams_idx = 1:length(RFdiams)
        RF_centerdist_deg = RFcenterdistances(RFcenterdistances_idx); % Distance between RF centers in DVA
        RF_diam_deg = RFdiams(RFdiams_idx);
        [x_deg,y_deg] = meshgrid(linspace(-sigma_gabor*3*sigmas_n,sigma_gabor*3*sigmas_n,40));
        Rad3Over2 = sqrt(3)/2;
        x_centers = [x_deg(1,1):RF_centerdist_deg:x_deg(1,end)];
        % Below is for equating the number of neurons (which does not
        % equate the extent of the infinite stimulus that is covered).
%         if INFINITELYLARGESTIMULUS
%             if length(x_centers) < 12
%                 error('too few x_centers')
%             end
%             if length(x_centers) > 12
%                 x_centers = x_centers(1:12);
%             end
%         end
        size(x_centers)
        closest_to_zero = find(abs(x_centers) == min(abs(x_centers)),1);
        x_centers = x_centers-x_centers(closest_to_zero);
        x_centers_mat = repmat(x_centers,size(x_centers,2),1);
        y_centers = x_centers;
        x_centers_mat = x_centers_mat*Rad3Over2;
        y_centers_mat = repmat(y_centers',1,size(y_centers,2));
        y_centers_mat(:,rem(find(y_centers == 0),2)+1:2:end) = y_centers_mat(:,rem(find(y_centers == 0),2)+1:2:end)+.5*RF_centerdist_deg;
        mus = zeros(numel(y_centers_mat),1);
        S2 = eye(numel(y_centers_mat)); % All neurons are assumed to have unit variance
        tmp = RF_diam_deg/2*[cos(linspace(0,2*pi,100))' sin(linspace(0,2*pi,100))'];
        if INDIVIDUALPLOTS
            figure; subplot(2,2,1); hold on;
        end
        for j = 1:numel(x_centers_mat)
            RFlocations = RF_diam_deg*[-.5:.01:.5];
            [x,y] = meshgrid(RFlocations,RFlocations);
            RFmask = sqrt(x.^2+y.^2) < RF_diam_deg/2;
            mus(j) = sum(sum(RFmask.*normpdf(x_centers_mat(j)+x,0,sigma_gabor).*normpdf(y_centers_mat(j)+y,0,sigma_gabor)));
            if INFINITELYLARGESTIMULUS
                mus(j) = 1;
            end
            for k = 1:numel(x_centers_mat) % Looking for RF overlap, brute-force
                r = RF_diam_deg/2;
                d = sqrt((x_centers_mat(j)-x_centers_mat(k))^2+(y_centers_mat(j)-y_centers_mat(k))^2);
                if 2*r-d > 1e-12 % Without this condition we get negative or imaginary covariances
                    S2(j,k) = 2*r.^2.*acos(d./(2*r))-(d./2).*sqrt(4*r^2-d.^2);
                end
                if S2(j,k) < 0
                    keyboard
                end
                S2(k,j) = S2(j,k);
            end
            if INDIVIDUALPLOTS
                plot(x_centers_mat(j)+tmp(:,1),y_centers_mat(j)+tmp(:,2),'b-');
            end
        end
        mus = mus./max(mus);
        S2 = S2./max(S2(:)); % If two RFs are identical, cov = 1 & cov is proportional to overlap
        weights = S2\mus; % these are the ideal weights. No need to normalize to max(weights)
        if INDIVIDUALPLOTS
            axis square; axis equal;
            subplot(2,2,2);
            imagesc(S2);
            axis square;
            subplot(2,2,3);
            imagesc(reshape(weights,size(x_centers_mat)))
            axis square;
        end
        % Version 1, assuming no correlations between neurons (way wrong)
        % population_scalefactor = sum(mus.^2)./sqrt(sum(mus.^2))
        
        % version 2, ideal observer assumes no correlation between neurons (just
        % adds them up, weighted by their mu) but there are correlations among
        % them.
        tmp = (mus*mus').*S2;
        sum(mus.^2)./sqrt(sum(tmp(:)));
        
        % version 3, ideal observer knows about the correlation between neurons and
        % adjusts the weights accordingly.
        tmp = (weights*weights').*S2;
        %(mus'*weights)./sqrt(sum(tmp(:)))
        data(RFcenterdistances_idx,RFdiams_idx) = (mus'*weights)./sqrt(sum(tmp(:)));
        data(RFcenterdistances_idx,RFdiams_idx) = sum(mus.^2)./sqrt(sum(sum((mus*mus').*S2))); % Non-ideal observer

        %data = [data; (mus'*weights) sqrt(sum(tmp(:)))];
    end
end
figure;
if length(RFdiams) > 1
    plot(RFdiams,data*sqrt(2),'.-');
    ylabel('Population scale factor');
    xlabel('RF diameter (',num2str(RFcenterdistances),'° spacing)');
end
if length(RFcenterdistances) > 1
    plot(RFcenterdistances,data*sqrt(2),'.-');
    ylabel('Population scale factor');
    xlabel(['RF spacing (',num2str(RFdiams),'° diameter)']);
end
% OK, I need to get to the bottom of this. What's going on? Does the weird
% dependance of population scale factor on RF spacing have to do with the 
% It seems to have to do with changes in the number of neurons.
% What if I keep the signal strength flat and the 

%%
% Section 13
% Computing the covariance of neurons from the overlap of Gaussian RFs

sigma_gabor = 0.15; % DVA
sigmas_n = 4; % truncation for Gabor
rc = .06; % Croner and Kaplan "rc" measure
RFsigma = rc/sqrt(2); % 1 SD NOTE: rc in Croner and Kaplan is not standard deviation, it's std*sqrt(2)
RFdistance = RFsigma*2; % 2 SD
sigma_gabor = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sigma'))); % DVA
sigmas_n = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sigmas_n')));
[x_deg,y_deg] = meshgrid(linspace(-sigma_gabor*2*sigmas_n,sigma_gabor*2*sigmas_n,40));
Rad3Over2 = sqrt(3)/2;
x_centers = [x_deg(1,1):RFdistance:x_deg(1,end)];
closest_to_zero = find(abs(x_centers) == min(abs(x_centers)),1);
x_centers = x_centers-x_centers(closest_to_zero);
x_centers_mat = repmat(x_centers,size(x_centers,2),1);
y_centers = x_centers;
x_centers_mat = x_centers_mat*Rad3Over2;
y_centers_mat = repmat(y_centers',1,size(y_centers,2));
y_centers_mat(:,rem(find(y_centers == 0),2)+1:2:end) = y_centers_mat(:,rem(find(y_centers == 0),2)+1:2:end)+.5*RFdistance;
mus1 = zeros(numel(x_centers_mat),1); % Product of Gaussians. All yield basically the same answers.
mus2 = zeros(numel(x_centers_mat),1); % Integrated contrast inside truncated Gaussian envelope.
figure; axes; hold on;
tmp = RFsigma.*[cos(linspace(0,2*pi,100))' sin(linspace(0,2*pi,100))']; % For plotting

bpdf_vec=@(x,y,mu_x,mu_y,sigma)(exp(-((x-mu_x).^2)./(2*sigma^2)-((y-mu_y).^2)/(2*sigma^2))./(2*pi*sigma^2)); % 2D Gaussian PDF
truncated_2D_gaussian = @(x,y,mu_x,mu_y,sigma,trunc) (sqrt((x-mu_x).^2+(y-mu_y).^2)<trunc).*bpdf_vec(x,y,mu_x,mu_y,sigma) + (sqrt((x-mu_x).^2+(y-mu_y).^2)>=trunc).*0;

interRFdistances = zeros(numel(x_centers_mat),numel(x_centers_mat));
for j = 1:numel(x_centers_mat)
    plot(x_centers_mat(j)+tmp(:,1), y_centers_mat(j)+tmp(:,2),'k-'); % RFs should meet at 1 sigma boundary
    % Formulas used below are from https://www.cs.nyu.edu/~roweis/notes/gaussid.pdf
    % See IsoSampAnalysis section 6.5
    b = [x_centers_mat(j); y_centers_mat(j)];
    %mus1(j) = exp(-.5*(b'*b)/(RF_diam_deg^2+sigma_gabor^2)); %https://math.stackexchange.com/questions/1260015/normalizing-factor-for-product-of-gaussian-densities-interpretation-with-bayes
    %mus1 is similar to mus2

    % Taking weighted, summed stimulus contrast. 
    % 0,0 is middle of RF. Lattice for counting is in terms of RF.
    % Could do this analyticall but numerically makes it easier to truncate.
    %overlap_point = @(x,y) bpdf_vec(x,y,0,0,RFsigma).*bpdf_vec(x,y,x_centers_mat(j),y_centers_mat(j),sigma_gabor);
    %mus2(j)=integral2(overlap_point,-4*RFsigma,4*RFsigma,-4*RFsigma,4*RFsigma);
    for k = 1:numel(x_centers_mat)
        interRFdistances(j,k) = sqrt((x_centers_mat(j)-x_centers_mat(k)).^2+(y_centers_mat(j)-y_centers_mat(k)).^2);
    end
end

% filling in S2
RFTRUNCATIONINSD = 1.5; % How many SDs is one RF. Points that are beyond this aren't inlcude in either RF
S2 = nan(numel(x_centers_mat));
for dist = unique(interRFdistances)'
    overlap_point = @(x,y) truncated_2D_gaussian(x,y,0,0,RFsigma,RFsigma*RFTRUNCATIONINSD).*truncated_2D_gaussian(x,y,0,dist,RFsigma,RFsigma*RFTRUNCATIONINSD);
    S2(interRFdistances==dist)=integral2(overlap_point,-RFTRUNCATIONINSD*RFsigma,RFTRUNCATIONINSD*RFsigma, dist-RFTRUNCATIONINSD*RFsigma,RFTRUNCATIONINSD*RFsigma,'method','iterated');

    %y_int_limit1 = @(x) -1*RFTRUNCATIONINSD/2*cos(x*pi/RFTRUNCATIONINSD);
    %y_int_limit2 = @(x) RFTRUNCATIONINSD/2*cos(x*pi/RFTRUNCATIONINSD);
   % if (dist >= 2*RFsigma*RFTRUNCATIONINSD)
   %     S2(interRFdistances==dist)=0;
   % else
    %    S2(interRFdistances==dist)=integral2(overlap_point,-RFTRUNCATIONINSD/2*RFsigma,RFTRUNCATIONINSD/2*RFsigma,y_int_limit1,y_int_limit2,'method','iterated');
   % end
end
if any(isnan(S2(:)))
    error('Nan in S2');
end
S2 = S2./max(S2(:)); % Cov depends on doproduct between RF and if two RFs are identical cov = 1

mus = mus2./max(mus2);
weights = S2\mus; % these are the ideal weights. No need to normalize to max(weights)
tmp = (weights*weights').*S2;
population_scalefactor=(mus'*weights)./sqrt(sum(tmp(:)));

% Testing code below
% Need a quick way of calculating the inner product of two 2-D
% Gaussian PDFs (important for calculating RF overlap for S2 and for
% calculating stimulus/RF overlap for mus). Assuming cov = 0 but unequal variances. 

% Distribution 1 will be our reference distribution. We're looking at
% distribution 2 through the window of distribution 1.
mu1 = [0 0]; sigma1 = 1;
mu2 = [2 2]; sigma2 = 1.5;

% Way 1.
bpdf_vec=@(x,y,mu_x,mu_y,sigma)(exp(-((x-mu_x).^2)./(2*sigma^2)-((y-mu_y).^2)/(2*sigma^2))./(2*pi*sigma^2));
overlap_point = @(x,y) bpdf_vec(x,y,mu1(1),mu1(2),sigma1).*bpdf_vec(x,y,mu2(1),mu2(2),sigma2);
integral2(overlap_point,-4,4,-4,4)

% Another way of doing it
% assuming one Gaussian is at (0,0) and the other is at (0,d)
d = sqrt((mu1(1)-mu2(1))^2+(mu1(2)-mu2(2))^2);
d = d./sigma1; % Normalized distances (scaling so that STD = 1)
bpdf_vec1=@(x,y,mu_x,mu_y,sigma)(exp(-((x-mu_x).^2)./(sigma^2)-((y-mu_y).^2)/(sigma^2))./(2*pi*sigma^2));
overlap_point = @(x,y) bpdf_vec1(x,y,0,0,sigma1).*bpdf_vec1(x,y,0,d,sigma2);
integral2(overlap_point,-4,4,-4,4)

% Brute force confirmation
tmp = linspace(-4,4,1000);
x = normpdf(tmp,0,sigma1);
x = x./sum(x);
y1 = normpdf(tmp,0,sigma2);
y2 = normpdf(tmp,d,sigma2);
y1 = y1./sum(y1); % This normalization is a broken if the second Gaussian is too far off to the side
y2 = y2./sum(y2); % This normalization is a broken if the second Gaussian is too far off to the side

x*y';
xx = x'*x;
yy = y1'*y2;
xx(:)'*yy(:)./(tmp(2)-tmp(1))^2

% Discretized bivariate Gaussian weighting of white random variables. Brute
% force way to get the covariance between neighboring neurons
mu1=[0 0];
mu2=[2 0]; % displacement should be along the first dimension
sigma = 1;
x =linspace(-5*sigma+min([mu1(1), mu2(1)]),5*sigma+max([mu1(1), mu2(1)]),200);
y = x-mean(x)+mean([mu1(2),mu2(2)]);
w1 = normpdf(x,mu1(1),sigma)'*normpdf(y,mu1(2),sigma);
w1(w1<max(w1(:))*.1) = 0; % Truncating the Gaussian
%w1 = w1./sum(w1(:)); % Probably unnecessary
w2 = normpdf(x,mu2(1),sigma)'*normpdf(y,mu2(2),sigma);
w2(w2<max(w2(:))*.1) = 0; % Truncating the Gaussian

%w2 = w2./sum(w2(:)); % Probably unnecessary

figure; axes;
imagesc(w1+w2)

niter = 10000;
data = [];
for i =1:niter
    noise = normrnd(0,1,length(x),length(y));
    data = [data; w1(:)'*noise(:) w2(:)'*noise(:)];
end
S2 = cov(data);
S2 = S2./max(S2(:))
corrcoef(data)

bpdf_vec1=@(x,y,mu_x,mu_y,sigma)(exp(-((x-mu_x).^2./(2*sigma^2))-((y-mu_y).^2)/(2*sigma^2))./(2*pi*sigma^2));
overlap_point = @(x,y) bpdf_vec1(x,y,mu1(1),mu1(2),sigma).*bpdf_vec1(x,y,mu2(1),mu2(2),sigma);
a = integral2(overlap_point,min(x),max(x),min(y),max(y));
overlap_point = @(x,y) bpdf_vec1(x,y,mu1(1),mu1(2),sigma).*bpdf_vec1(x,y,mu1(1),mu1(2),sigma);
b = integral2(overlap_point,min(x),max(x),min(y),max(y)); % Integral of Gaussian with itself
a/b
%overlap_point = @(x,y) mvnpdf([x, y],mu1).*mvnpdf([x, y],mu2);
%integral2(overlap_point,min(x),max(x),min(y),max(y));

% RF center separation 2 sigma: covariance = 0.367
% RF center separation 4 sigma: covariance = 0.014

% Can I reduce this to a 1-D problem?
% seems like no.
f = @(x,mu) (1/sqrt(2*pi).*exp(-(x-mu).^2/2));
g = @(x) f(x,0).*f(x,0);
self = integral(g,-3, 3);
g = @(x) f(x,0).*f(x,2);
cross = integral(g,-2, 2);
cross/self



%%
% Section 14
% Hacking around with Victor and Purpura spike distance stuff
a = [1 2 3 1 2 3 7 7.9 9];
b = [1 4 7];
c = [3 6 9];

figure; axes; hold on;
for i = 1:length(b)
    spiketimes = a(b(i):c(i))
    plot(spiketimes, i,'o')
end
d = spkdl(a,b,c,.01); % 10 ms
d_mat = reshape(d,length(b),length(b));
d_mat = d_mat+diag(nan(length(b),1)); % filling in the diagonal with nan so a spike train isn't closest to itself
figure;
imagesc(d_mat);

%%
% Section 15
% Hacking around with maximum likelihood fitting of ROCs
% First I'm going to need a function that takes a parametrized ROC and a
% data point and returns the maximum likelihood of that point (finds the
% criterion that maximizes the likelihood of that point under the model and
% returns it).

% this is a two parameter fit. we can assume mu0 = 0 and sigma0 = 1.
mu0 = 0;
sigma0 = 1;
mu1 = 1.5;
sigma1 = 1;

% figure; axes; hold on;
% crits = linspace(-5,5,100);
% plot(1-normcdf(crits,mu0,sigma0),1-normcdf(crits,mu1,sigma1));
% axis square;
% xlabel('Prob FA'); ylabel('Prob HIT');

% simulating an couple of experiments
nreps = 10; % Number of simulation repeats
ntrials = 50; % Number of trials per condition
c = [1.5]; % criteria

pretend_data = [];
for i = 1:length(c)
    for j = 1:nreps
        pFA = 1-normcdf(c(i),mu0,sigma0);
        pHIT = 1-normcdf(c(i),mu1,sigma1);
        x = [binornd(ntrials, pFA) binornd(ntrials, pHIT)];
        pretend_data = [pretend_data; x ntrials ntrials];
    end
end

% For importing data from For SC_stimcueAnalyses.m
% Substituting FA_laser for FA and FA_nolaser for hits
%pretend_data = [choice_data(:,3) choice_data(:,4) choice_data(:,3)+choice_data(:,7) choice_data(:,4)+choice_data(:,8)];

% % Now creating a function that returns the (maximum over criteria) likelihood
% % for each point
% % Four terms: FA + CR + HIT + MISS
% generic_ll = @(crit,x,mu1,sigma1) x(1).*log(1-normcdf(crit,0,1)) + (x(3)-x(1)).*log(normcdf(crit,0,1))+x(2).*log(1-normcdf(crit,mu1,sigma1)) + (x(4)-x(2)).*log(normcdf(crit,mu1,sigma1));
% total_negllik = 0;
% for i = 1:size(pretend_data,1)
%     specific_negll = @(x) -1*generic_ll(x,pretend_data(i,:),mu1,sigma1);
%     [p_hat, negllik] = fminsearch(specific_negll,0);
%     total_negllik = total_negllik+negllik;
% end
% total_negllik

% Trying a range of mu1s and sigma1s.
% When I'm using real data, currently, the laser trials are N(0,1) and the NO
% laser trials are N(mu, sigma). When I do this for real I should reverse
% this.
data = [];
mus = linspace(0,.5,15);
sigmas = linspace(.5,1.5,15);
for i = 1:length(mus)
    mu_candidate = mus(i);
    for j = 1:length(sigmas)
        sigma_candidate = sigmas(j);
        
        % preparing for a lookup table for good initial criterion guesses
        tmpcrits1 = linspace(-6,6,10);
        tmpcrits2 = linspace(-6,6,10)*sigma_candidate+mu_candidate;

        generic_ll = @(crit,x,m,s) x(1).*log(1-normcdf(crit,0,1)) + (x(3)-x(1)).*log(normcdf(crit,0,1))+x(2).*log(1-normcdf(crit,m,s)) + (x(4)-x(2)).*log(normcdf(crit,m,s));
        total_negllik = 0;
        for k = 1:size(pretend_data,1)
            initialguess1 = interp1(1-normcdf(tmpcrits1,mu0,sigma0),tmpcrits1,pretend_data(k,1)/pretend_data(k,3));
            initialguess2 = interp1(1-normcdf(tmpcrits2,mu_candidate,sigma_candidate),tmpcrits2, pretend_data(k,2)/pretend_data(k,4));
            initialguess = nanmean([initialguess1,initialguess2]);

            specific_negll = @(x) -1*generic_ll(x,pretend_data(k,:),mu_candidate,sigma_candidate);
            [crit_hat, negllik] = fminsearch(specific_negll,initialguess);
            total_negllik = total_negllik+negllik;
        end
        data(i,j) = total_negllik;
    end
end

figure; 
subplot(2,1,1); hold on;
surf(mus,sigmas,data');
xlabel('mu'); ylabel('sigma');
[i,j] = ind2sub(size(data),find(data == min(data(:))));
plot3(mus(i),sigmas(j),min(data(:)),'m*');
axis square;
subplot(2,1,2); hold on;
crits = linspace(-5,5,100);
plot(crits,normpdf(crits,0,1),'b-'); % laser
plot(crits,normpdf(crits,mus(i),sigmas(j)),'k-'); % no laser

% Plotting "best fit"
figure; axes; hold on;
plot(1-normcdf(crits,mu0,sigma0),1-normcdf(crits,mus(i),sigmas(j)),'k-');
%plot(1-normcdf(crits,mu0,sigma0),1-normcdf(crits,mu1,sigma1),'b-'); % The truth
plot([0 1],[0 1],'k:'); % identity
axis square;
plot(pretend_data(:,1)./pretend_data(:,3),pretend_data(:,2)./pretend_data(:,4),'ko');
xlabel('Prob FA'); ylabel('Prob HIT');

% Trying to come up with a good estimate of mu1 and sigma1
prop_fa = pretend_data(:,1)./pretend_data(:,3);
prop_hit =  pretend_data(:,2)./pretend_data(:,4);
z_fa = 1-norminv(prop_fa); 
z_fa(z_fa < -4) = -4; z_fa(z_fa > 4) = 4;
z_hit = 1-norminv(prop_hit);
z_hit(z_hit < -4) = -4; z_hit(z_hit > 4) = 4;

figure; axes; hold on;
plot(z_fa,z_hit,'ko');
b = regress(z_hit, [z_fa ones(size(pretend_data,1),1)]);
plot([0 2],b(1)*[0 2]+b(2));

% estimate for sigma1 = 1/b(1)
% estimate for mu1 = 2/(1/b(1)+1)
[mu1 b(2) sigma1 1/b(1) ]

%%
% Making sure I am doing the right transformation to convert N(0,1) &
% N(mu1, sigma1) to N(mu2, sigma2) & N(0,1) while keeping the ROC curve
% fixed. I want to estimate the LASER noise distribution relative to the 
% NO LASER (standard) noise distribution, which is assumed to be N(0,1). 
% But FitROC.m assumes that distribution with the lower mean is N(0,1).

% Testing
muhat = 2;
sigmahat = 2;
tmp = linspace(-4,4,100);
figure; 
subplot(2,2,1); hold on;
plot(tmp, normpdf(tmp,0,1),'b-');
plot(tmp, normpdf(tmp,muhat,sigmahat),'k-');
subplot(2,2,2); hold on;
plot(tmp, normpdf(tmp,0,1),'k-');
plot(tmp, normpdf(tmp,-muhat/sigmahat,1/sigmahat),'b-');
subplot(2,2,3); hold on;
plot(1-normcdf(tmp,0,1),1-normcdf(tmp,muhat,sigmahat),'o-'); % These two ROCs are the same.
plot(1-normcdf(tmp,-muhat/sigmahat,1/sigmahat),1-normcdf(tmp,0,1),'ko-');
axis square;

%%
% Section 16
% Comparing difference in medians versus median of differences
niter = 2000;
n = 20;
data = [];
for i = 1:niter
   % a = normrnd(0,1,n,2);
    a = exprnd(1,n,2);
    data(i,:) = [median(a(:,1))-median(a(:,2)) median(a(:,1)-a(:,2))];
end
bins = linspace(min(data(:)), max(data(:)),30);
figure; axes; hold on;
hist(data,bins);
legend({'diff. of med.','med. of diff.'})


mean(data)
std(data)

% 