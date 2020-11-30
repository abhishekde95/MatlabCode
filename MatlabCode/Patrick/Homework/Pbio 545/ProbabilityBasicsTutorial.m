% probabilityBasicsTutorial.m
%
% This short tutorial provides an introduction to probability distributions and
% random variables. It covers the following concepts: mean, expectation,
% moments, standard deviation, sample mean and standard error. You may need to
% look up the definition of: probability function, probability density function
% and a cumulative probability function. You'll have a feel for these by the end
% even if you don't get formal.
%
% Requires the statistics toolbox.

% 2/12/05 mns updated
% 2/9/07  mns cleaning
% 2/10/07 mns adding sums and convolution
% 2/15/09 mns minor clean up
% 1/28/12 mns update

%% Matlab setup
clear all
rand('state',sum(100*clock))
format short g

%% Generating random numbers
% Fortunately, matlab has a number of tools for generating random numbers
% and plotting probability distributions. You must have the statistics
% toolbox. Try typing 'help stats' at the command prompt to see what's
% available. There are three types of functions that we will deal with. The
% ones ending in ...rnd.m generate random numbers.
%
% Let's generate 200 random numbers from the Poisson distribution
n = 200
mu = 10.3     % mean (please use a non-integer or later illustrations will suffer)
x = poissrnd(mu,n,1);
% Plot a Frequency histogram
figure(1), clf
[count vals] = hist(x,[0:1:40]);
freq = count/n;    % this is a critical step. Make sure you appreciate it
sum(freq)           % notice that the sum of all the frequencies come out to 1 (the sum of the counts is n)
bar(vals,freq) 
set(gca,'XLim',[-.5 max(vals)+.5])
xlabel('Value of random number (x)')
ylabel('Frequency')

%% Expectation
% The expectation of a random variable (RV) is what it sounds like it ought
% to be. Think of it as the mean of an infinitely large sample of random
% values drawn from the distribution. It's a deeper concept though, because
% it is a property that we can associate with a variable without ever
% measuring samples. It's so fundamental that we will actually think of the
% mean as something (that we calculate) and that will return the
% expectation of the RV.
% 
% The expecation of x, which I will denote E[X] is whatever you set mu to above.
% Let's try to estimate the expectation by calculating the mean. Here are three
% ways to do this in matlab. 
mean(x)     % 1. matlab's function
sum(x/n)    % 2. sum of the random numbers divided by n (which is what matlab's mean does)
sum(freq .* vals) % 3. Weighted sum of the values that x can take. 
% #3 is less obvious, but very important. Make sure you understand why it works.
% It is a key step toward understanding lots of expectations. This is the
% expectation of x. 
% 
% If you've convinced yourself why the 3rd approach works, we can calculate the
% expectation, E[X]. We can do this without actually generating any samples, so
% long as we know the probability of observing the values. The matlab function
% POISSPDF will tell us these probabilities. For the Poisson distribution, the
% only possible values are integers. The probability of seeing any particular
% value is
theoryFreq = poisspdf(vals,mu);
% Notice that the total probability is 1
sum(theoryFreq)
% The E[X] is the weighted sum of the possible values of X. The weights are the
% probabilities (i.e., the expected frequencies of observing each possible value
% of X.
expectation = sum(vals .* theoryFreq)
% You can see why the sample mean only approximated the expectation. Look at the
% theoretical frequencies superimposed on the frequency
% histogram we made with our n samples. 
figure(1), hold on
h = plot(vals, theoryFreq, 'rV','MarkerFaceColor','r')

% There's a 4th way to compute E[x]. It uses the moment generating function. But
% I'll hold off on this for now. 

%% moments
% The expectation of X is also known as the first moment of the random variable,
% x. 
firstMoment = expectation

% The second moment is the expectation of x^2. The 3rd moment is the
% expectation of x^3, and so on. The second moment is related to the variance.
% The 3rd moment is related to the skew of (the degree of asymmetry of the
% distribution about its mode). The 4th moment is related to the kurtosis of the
% distribution (the relative concentration of probability around the mode
% compared the tails -- in other words, the peakiness of the distribution). We
% rarely worry about moments greater than 2nd.

% Important point. You can calculate the expectation of any function of X by
% taking the weighted sum of all possible f(x), the weights given by the
% probability of observing each possible value of x.
%
% So, if for some crazy reason, we wanted to know the expected value of
% log(1+x), we could compute it as
weirdExpectation = sum(log(1+vals) .* theoryFreq)

% A more useful calculation is the E[X^2] (the 2nd moment). 
secondMoment = sum((vals.^2) .* theoryFreq)
% In class we will derive that E[(x-mu)^2], which is known as the variance or 2nd central moment is
secondMoment - mu^2
% or
variance = secondMoment - firstMoment^2


 
%% Standard error of the mean (SEM), Variance and Standard Deviation. 
%
% Let's start with a brief overview about what these terms mean and then we'll
% convince ourselves of the math. The standard deviation is a measure of the
% spread of the random values. It conforms to our intuitive sense of how
% variable the random numbers are, hence how wide a distribution is from which
% they are drawn. Standard deviation is also the number we would want to rely on
% when we talk about uncertainty. What I mean by that is that if I tell you that
% I measured the number 5.2 and you wanted to know how close this might be to a
% true value (what I would get if there were no souces of noise affecting my
% measurement), you would want to know the standard deviation. You couldn't know
% that from my one sample, but if you did know it, you would be able to guess
% how close 5.2 might be to other samples I might draw on another occasion. 
%
% If instead of a single value, I told you a number that I calculated, say a
% mean from 20 samples, you might want to know how close this sample mean might
% be to other sample means, were I to obtain 20 samples again. Again, you would
% want to know the standard deviation of this sample mean. This standard
% deviation of a statistic (in this case the sample mean) is called the standard
% error (in this case, the standard error of the mean). If we calculated the
% slope of a line, then we would want to know the standard deviation of such
% slope estimates in order to get a sense of how close our estimate was to other
% estimates we might obtain from another data set. Again, this standard
% deviation is called a standard error (but not standard error of the MEAN, of
% course).
%
% I have already defined the variance as the expectation[(x-mu)^2]. In words:
% take each X and subtract off its expectation; then square. Do this for an
% infinite sample size of X's. 
secondCentralMoment = sum((vals-mu).^2 .* theoryFreq)

% Notice that for the Poisson distribution, the variance equals the expectation
% (mu). 

% This is different than the variance you would compute from the n samples. Just
% as you don't get E[X] when you compute the mean from n samples, 
sampleMean = mean(x)
% you don't get the true variance. If you ask matlab to calculate the variance
% from the n samples, you get
var(x)
% For the sample mean, matlab takes the average value of the X's. But, for the
% sample variance, matlab multiplies the average of the (X-sampleMean)^2 terms by
% n/(n-1). Let's convince ourselves of this.
(n/(n-1)) * mean((x-sampleMean).^2)
% This is the same as
sum((x-sampleMean).^2) / (n-1)
% Notice that this number is a little bigger than the average of the squared
% departures from the mean
mean((x-sampleMean).^2)

% This business of dividing by n-1 instead of n drives a lot of people
% crazy. I don't want to get bogged down in the derivation (but see Brian
% Lundstrom's tutorial on this). Nonetheless, I think it's worth have an
% intuition about why we do this. Notice that the definition of variance is
% E[(X-mu)^2]. But we don't know mu from a sample of data. All we can do is
% estimate mu by calculating the mean. Now, let me say something that
% sounds overly wordy and almost obvious: the sample mean is the best
% estimate we can make of the true E[X]. What I mean by that is that if we
% could draw another n samples and compute another sampleMean and we did
% this a gzillion times, the average of the sampleMeans would be very near
% E[X]. More formally E[sampleMean] = E[X]. Technically, the sampleMean is
% an unbiased estimator of the E[X]. Better yet, it is the best unbiased
% estimator, meaning that is has the smallest dispersion (standard
% devation) around E[X]. Let's illustrate that by repeating the process of
% making sampleMeans several times
nreps = 1000     % the number of repetitions
Y = poissrnd(mu,n,nreps);
sampleMeanYs = mean(Y);
% Let's look at these sample means superimposed on the distribution of 200
% samples we saw before.
[count vals] = hist(sampleMeanYs,[0:.1:40]);
freq = count/nreps;    
sum(freq)           % notice that the sum of all the frequencies come out to 1 (the sum of the counts is n)
figure(1), hold on
h = bar(vals,freq)
set(h,'FaceColor', 'y','EdgeColor','none')
s = sprintf('Blue: sample n = %d; Yellow: distribution of means from %d repetitions', n, nreps)
title(s)
% Notice that the sample means fall near the true theoretical mean. I didn't
% prove that expectation of these sample means is the true mean, but it's true.
% The dispersion around this true mean can be calculated by taking the standard
% deviation. But since I have yet to develop standard deviation, let's use
% a more intuitive statistic. What is the average absolute difference between
% the sample means and mu (the E[X])?
meanAbsDev_sampleMeans = mean(abs(sampleMeanYs - mu))
% To illustrate my point about this being the best you can do, let's come up
% with an alternative measure of the E[X] based on a sample. How about the
% sample median?
sampleMedianYs = median(Y);
meanAbsDev_sampleMedianYs = mean(abs(sampleMedianYs - mu))
% Notice that the deviations are smaller on average for the sample means.
% (Note, this will fail in the special case that mu is an integer.) This is
% not a proof that the sample mean is the most efficient statistic for
% estimating the E[X], but it gives you a feel for what we're talking
% about.
%
% Now, let's return to the sample variance. Why the 1/(n-1) instead of 1/n
% in the formula for the sample variance? It's because we want an estimate
% of the variance from our n samples that is unbiased and maximally
% efficient. if we happen to know mu, then taking the average of the n
% terms, (X-mu)^2, would be an unbiased estimator of E[(X-mu)^2]. But we
% don't know mu from our n samples. All we know is some estimate of mu,
% like the sample mean. But the sample mean is computed from the n samples.
% This is kind of a cheat. The sample mean is the value we would choose
% that is nearest all the samples, on average. The actual expectation,
% which is similar to but not exactly the same as the sample mean, is going
% to be further away from the samples, on average. This means that
var1 = sum((x-sampleMean).^2) / n
% is a biased estimator of the true variance. It is biased to be smaller
% than the true variance: E[var1] < E[(X-mu)^2] We can see that from the
% simulated samples we took. Matlab will compute this biased version of the
% variance if you pass a second argument var(x,1)
var(x,1)
% If we use the standard matlab var(x) we get back the 1/(n-1) version of
% the sample variance. Confirm this:
sum((x-sampleMean).^2) / (n-1)
var(x)

% Let's use these two versions of var() to operate column-wise on the Y's we
% made above. Remember, Y contains nreps columns; each contains a sample of n
% random values.
sampleVarWithDenomNs = var(Y,1);
sampleVarEithDenomNminus1s = var(Y);
% Take the mean of these nrep samples
mean(sampleVarWithDenomNs)
mean(sampleVarEithDenomNminus1s)
% compare to the true variance
variance

% I'm going to come back to these examples in a moment. Let's save the key
% variables in a matfile.
save poissonRandNumbers.mat mu Y x

% I didn't show you that the sample variance using the n-1 denominator is
% unbiased and maximally efficient. But you can make nreps bigger to
% convince yourself that it's on track.

% So we know the definition of the mean and variance of a random variable
% and we know how to estimate the mean and variance from a finite sample of
% observations. Just to be clear, the mean of the RV is the RV's
% expectation, which is also known as the 1st moment. The variance of the
% RV is the 2nd central moment, that is the expectation of the square of
% the quantity, RV minus its expectation. The statistics we compute from
% samples are unbiased, maximally efficient estimators of these
% expectations.

% The standard deviation of a random variable is the square root of its
% variance. We could state that in 4 or 5 ways using all the terminology
% for the variance, but let's get a feel for standard deviation instead.
% This is the statistic that conforms best to the spread of the
% distribution. It's the square root of a squared difference between a RV
% and its mean. So it's a lot like taking an absolute value. You should
% think for a second about why sqrt(var(x)) is not the same thing as
% mean(abs(x-mu)) Matlab's
std(x)
% is the same as
sqrt(var(x))
% but is not the same as
mean(abs(x-mu))
% or even 
(n/(n-1)) * mean(abs(x-mu))
%
% The reason stdev is so useful is this. It captures the sense of spread of a
% distribution. Let's consider the normal distribution
mu = 50
sigma= 10
n = 100
x = normrnd(mu,sigma,n,1)
% You can try others with the stats toolbox


figure(2), clf, hold on
% make an x axis
dx = 5;
xax = [0:dx:100]';
[nn xx] = hist(x,xax);
hb1 = bar(xx,nn/n,1)
set(gca,'FontSize',18,'Tickdir','out','XLim',[0 100])
xlabel('x')
ylabel('Frequency of observation')
% Here's the theoretical distribution
% See if you can figure out what this next line achieves
pdfBinned = diff(normcdf([xax-dx/2; xax(end)+dx/2],mu,sigma));
sum(pdfBinned)  % total probability should be 1
% Superimpose the predicted probabilities on the histogram
h = plot(xx,pdfBinned,'rV-','MarkerFaceColor','r');

% Here are our sample mean, standard deviation and variance
xmean = mean(x)
xs = std(x)
xvar = var(x)

% The true standard deviation is related to the width of the distribution. At
% half height, the width of the Normal is 2.33 sigma.
figure(2), hold on
hw = line([mu-2.33*sigma/2 mu+2.33*sigma/2], max(pdfBinned)*[1 1]/2,'Color','r','LineWidth',3)

% There's another way to think about the spread of the distribution. Suppose we
% were to ask, What fraction of the samples are greater than 60? In other words,
% what fraction of the total probability is accounted for by values >60? This is
% also a question about spread. Here's how we answer it using the standard
% deviation.
% 
crit = 60
% How many standard deviations from the mean is our criterion? 
crit_in_stdevs = (crit-mu)/sigma
% What's the cumulative probability up to this value? To get that use the
% cumulative probability function.
normcdf(crit,mu,sigma)
% The remaining probability -- that is, the area under the probability density
% to the right of this 
1 - normcdf(crit,mu,sigma)
% Let's see often we really did see a value in our sample greater than 60?
sum(x > crit)/n

% Now, suppose we were to turn a gain dial and thereby change all the x values
% by some factor. Where would we place our criterion to extract the same
% fraction of the observations? 
gain = .5
x2 = gain * x;
figure(2), hold on
[nn2 xx] = hist(x2,xax);
hb2 = bar(xx,nn2/n,1)
set(hb2,'FaceColor','none','EdgeColor','g')
% Notice that the spread has changed. This is captured by applying the gain to
% the standard deviation.
mu2 = gain * mu
sigma2 = gain * sigma
pdfBinned2 = diff(normcdf([xax-dx/2; xax(end)+dx/2],mu2,sigma2));
sum(pdfBinned2)  % total probability should be 1
% Superimpose the predicted probabilities on the histogram
h2 = plot(xx,pdfBinned2,'m^--','MarkerFaceColor','m');
% Not surprisingly, our criterion should be the same distance from the new mean in
% units of the new standard deviation.
crit2 = mu2 + crit_in_stdevs * sigma2
sum(x2 > crit2)
% Think about what we did here. By scaling the observations, it's as if we
% recoded them in some new measurement unit (like changing from inches to
% centimeters). The basic characteristization of the spread in the
% measurements should not be affected by this. And that is captured by
% thinking about the observations in terms of distance from the mean in
% units of standard deviation. Importantly, muliplying the random values by
% a constant multiplies the mean and the standard deviation by this same
% amount. The ratio, standard deviation divided by the mean, gives us an
% indication of the spread of the distribution, regardless of what units we
% use to describe the observations. This ratio is called the coefficient of
% variation (CV) and its reciprocal is sometimes called the signal to noise
% ratio. Changing the units, or turning up the amplifier, or increasing the
% power on the microscope, etc. do not change the dispersion of the data,
% just the units we use to describe the data.
%
% Notice that the CV does not change. But, since the variance is the square
% of the stdev, it follows that another ratio, the variance divided by the
% mean, does change when we scale units differently. We'll discuss
% variance/mean in the stochastic processes tutorial. It's worth keeping
% this point in mind or coming back to it when we encounter this ratio.

%% Finally, let's turn to the standard error of the mean.
%
% Let's look at Figure 1 again
figure(1)
% and reload our Poisson random numbers we used to make it.
load poissonRandNumbers.mat

% Recall that the blue shows a single sample of n observations. This has a
% single sample mean. The yellow shows the sample means from nreps
% repetitions of this experiment. It's pretty obvious that the sample means
% fall near the true mean. And, as we noted earlier, the spread in the
% sample means about the true mean is a lot less than the spread of the
% observations in any one experiment. The standard error of the mean is
% simply the standard deviation of these sample means.
stdOfSampleMeans = std(sampleMeanYs)
% Compare this to the standard deviation from the first sample of
% observations
std(Y(:,1))
% This set of observations is not off the mark. Remember, we used a Poisson
% distribution whose mean = 10. So the variance should also be 10 and the
% standard deviation should be
sqrt(mu)
% Now, here's the really magical thing. In real life, we rarely get a chance to
% repeat a set of measurements a bunch of times. In real life, we usually get a
% single set of n observations. We only get one sample mean. How can we
% tell how near it might be to other sample means? In other words, how can we
% estimate the standard deviation of the sample mean from just one sample set?
% The answer is, take the sample standard deviation and divide by the square
% root of the number of observations that comprise it.     
n = size(Y,1)
std(Y(:,1)) / sqrt(n)
% It's very easy to derive this expression. It's based solely on two
% points. First, if you make a new random variable by adding two random
% numbers, A + B, the variance of the sum is the sum of the variances
% associated with A and B. This holds if A and B are independent, that is
% the joint probability that A=x and B=y is the product of the probability
% that A=x (regardless of B) times the probability that B=y (regardless of
% A). The second principle is that the sample variance is an unbiased
% estimator of the true variance. We've already made that point. It follows
% that the sum of n variables, each with mean mu and variance s^2, has
% E[sum] = n*mu
% and
% Var[sum] = n * s^2
% Stdev[sum] = sqrt(n) * s
% To change from sums to the mean, we could scale the sums by 1/n.
% The mean has expectation E[sum]/n = mu
% The standard deviation of the mean = sqrt(n) * s /n = s/sqrt(n)



% Excercise (optional): Graph three distributions, Poisson, Normal (aka,
% Gaussian), and Uniform. Do this so that they each have the identical mean
% = 9 and the same standard deviation = 3. How do these distributions
% differ from each other. [Hint: skewness & kurtosis].

%% Sums and convolution
% What is the distribution of a random number, z, formed as the sum of two
% random numbers, x and y?
%
% Let's fool around with matlab's random number generator to get a feel.
% Set the number of samples
n = 1000
% Suppose x is Poisson distributed with mean
xmean = 4
x = poissrnd(xmean,n,1);
% Suppose y is also a Poisson distributed with mean
ymean = 10
y = poissrnd(ymean, n, 1);
% z is the sum
z = x+y;
figure(1), clf, hold on
rvaxis = [0:30];
xf = hist(x,rvaxis) / n;
yf = hist(y,rvaxis) / n;
zf = hist(z,rvaxis) / n;
stem(rvaxis, xf,'b')
stem(rvaxis, yf,'g')
stem(rvaxis, zf,'r')

% Verify that the mean of the sum is the sum of the means.
mean(z)
mean(x) + mean(y)
% Of course, mean(x) is a sample mean. We would also say that the
% expectation of z is the sum of the expectations for x and y.
% 

% Same hold for for the variance
var(z,1)
var(x,1) + var(y,1)


% Notice that the distribution of z is shifted to the right (higher mean)
% and also spread out. We might say it is blurred. It so happens that z is
% also a Poisson distribution. It is not the case in general that the sum
% of two random variabls retains its orginal form. It can only hold
% if the variance is proportional to the mean. That is still not a
% guarantee.
% 
% It may be easier to see this if we overlay the theoretical distributions 
plot(rvaxis,poisspdf(rvaxis,xmean),'b--')
plot(rvaxis,poisspdf(rvaxis,ymean),'g--')
plot(rvaxis,poisspdf(rvaxis,xmean+ymean),'r--')

% Convince yourself that the red distribution looks like the blue one (x)
% blurred by the green one (y), or the green one blurred by the blue one. 
% This is very clear if we run through the steps above using a narrow
% gaussian for x and a broad distribution for y.
n = 10000
xmean = 0, sigma = .2
x = normrnd(xmean,sigma,n,1);
ymean = 1
y = 1 + exprnd(ymean, n, 1);
z = x+y; % z is the sum
figure(2), clf, hold on
drv = .1;
rvaxis = [-1:drv:10];
xf = hist(x,rvaxis) / n;
yf = hist(y,rvaxis) / n;
zf = hist(z,rvaxis) / n;

stairs(rvaxis, xf,'b')
stairs(rvaxis, yf,'g')
stairs(rvaxis, zf,'r')

% This blurring ought to remind you of the convolution operation we did in
% the Fourier analysis module. In fact the distribution of z is the
% convolution of the PDFs for x and y. To see this, consider what it takes
% to get any particular value of z. If you happen to pick the value tau
% from the x distribution, then in order to get a particular value of z,
% you need to pick z-tau from the y distribution. 
%
% Let's call the probability distribution (PDF) of x, b(x) [the blue
% function]
% Let's call the probability distribution (PDF) of y, g(y) [the green
% function]
% Let's call the probability distribution (PDF) of z, r(z) [the red
% function]
% Our claim is that r(x+y) = b(x) ** g(y), where ** denotes convolution.
% To get any value z=zeta, we can pick any combination of tau from the x
% distribution [with probabilty b(tau)] and z minus tau from the
% y-distribution [with probability g(z-tau)]. Therefore to get any value of
% z, we take the weighted sum of all the z-tau, across all possible values
% of tau. The weighting is b(tau). 
% But this is the integal with respect to tau of b(tau)*g(z-tau). That's
% convolution!
%
% So here's the convolution of the x and y distributions
conv_xy = real(ifft(fft(xf).*fft(yf)));
% We labeled a point about 11 in as 0. But the function doesn't realize we
% mean that to be 0 and "thinks" that both of the functions are offset by
% this amount. So we'll just delete these points by shifting the result to
% the left by this amount. 
nz = find(rvaxis == 0)
conv_xy(1:nz-1) = [];
conv_xy(end:end+nz-1) = 0;
h = plot(rvaxis,conv_xy,'k--')
set(h,'Color',[.7 .7 .7],'LineWidth',3)

% To cap off this section, you might be wondering whether there's some kind
% of parallel between the convolution/multiplication duality we discussed
% in the Fourier transform section. Is there something analogous to fourier
% transform that would render convolution of PDFs of random variables
% multiplication of some transform of the RVs. The answer is yes. These are
% called moment generating functions and characteristic functions. They
% look just like Laplace and Fourier transforms. The transformed variable
% is not a frequency, but it has some interesting properties. That's all
% I'll say about this topic. See momentGeneratingFuncTutorial.m 





