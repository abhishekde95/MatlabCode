% Contents
%
% Section 1: Hacking around with linear regression and the inv(X'*X) term.
% Seeing how this affects STAs of linear neurons and LN neurons.
%
% Section 2: Comparing two potential corrections to the STA (stimulus white
% in gun space, STA calculated in cone space). inv(M') and inv(M'*M).
%
% Section 3: 
%% Section 1
% Hacking around with stuff for the Vision Research paper
% In particular, how does linear regression work on an LN
% neuron when the stimuli are equally spaced on the perimeter of
% an ellipse? 
% Transforming to a white space, calculating the STA and transforming back 
% works surpringly well for points that are uniformly positioned on an ellipse.
% RING = 0 Gaussian
% RING = 1 Uniform on a 1 STD ellipse  
% RING = 2 Uniform inside a box with ~same cov matrix as Gaussian
% RING = 3 Non-uniform on a 1 STD ellipse (uniform on a circle)

LN = 2; % 0 = linear, 1 = half-wave rectified, 2= half-squaring
WHITESTIM = 0; % mixMat is the identity matrix
v = [10; 0]; % preferred direction on neuron (L,M)
n = 100; % n stimuli
if WHITESTIM
    mixMat = sqrt(2)*eye(2);
else
    mixMat = mkbasis(normrnd(0,1,2,2));
   % mixMat = [sqrt(2)-.01 sqrt(2)+.01; sqrt(2)+.01 sqrt(2)-.01];
end
figure;
plotcounter = 1;
WHICHRINGS = [0,1,2,3] % Which conditions to try
for RING = WHICHRINGS
    if (RING == 0)
        X = normrnd(0,1,n,2)*mixMat;
    elseif (RING == 1) % Note, this doesn't make the covariances identical between the RING and ~RING cases
        % But it does put the ellipse at 1 STD (?) of the ~RING distribution.
        [eigv,d] = eig(mixMat'*mixMat);
        [xy,tk] = equiellipse_roots(sqrt(d(2,2)), sqrt(d(1,1)), n); % first two arguments are major and minor axis lengths
        theta = atan2(eigv(2,2),eigv(1,2));
        rotmat = [cos(theta) sin(theta); -sin(theta) cos(theta)];
        X = xy*rotmat;
    elseif (RING == 2) % Uniform within a parallelogram
        [eigv,d] = eig(mixMat'*mixMat);
        theta = atan2(eigv(2,2),eigv(1,2));
        rotmat = [cos(theta) sin(theta); -sin(theta) cos(theta)];
        xy = [unifrnd(-sqrt(3*d(2,2)),sqrt(3*d(2,2)),n,1),...
            unifrnd(-sqrt(3*d(1,1)),sqrt(3*d(1,1)),n,1)];
        X = xy*rotmat;
    elseif (RING == 3) % Uniform on a circle
       [x,y] = pol2cart(linspace(0,2*pi*(n-1)/n,n)',ones(n,1));
       X = [x,y]*mixMat;
    end
    cov(X)
    subplot(length(WHICHRINGS),2,plotcounter); plotcounter=plotcounter+1;
    plot(X(:,1),X(:,2),'.'); axis equal;
    lingensig = X*v;
    if (LN > 0)
        lingensig(lingensig<0) = 0;
    end
    if (LN == 2)
        lingensig(lingensig>0) = lingensig(lingensig>0).^2;
    end
    
    % Identical to the correctedSTA calculated below
    %STS = sum(X.*repmat(lingensig,1,2));
    %STA = STS/n; % Ignoring the fact that variance of X ~= 0
    %correctedSTA = inv(X'*X)*STS'
    
    % Whitening the stimuli to see what they look like
    whtmat = sqrtm(inv(X'*X));
    Xprime = X*whtmat;
    cov(Xprime)
    subplot(length(WHICHRINGS),2,plotcounter); plotcounter=plotcounter+1;
    plot(Xprime(:,1),Xprime(:,2),'o'); axis equal;
    STS = sum(Xprime.*repmat(lingensig,1,2));
    correctedSTA = STS*whtmat';
    title(num2str(mkbasis(correctedSTA)));
end


%STS_prime = sum(Xprime.*repmat(lingensigprime,1,2));
%STA_prime = STS_prime/n
%correctedSTAprime = inv(Xprime'*Xprime)*STS_prime'

% correctedSTA and correctedSTAprime are identical as long as neuron is
% linear (Gaussian distn is unimportant).
% If the neuron is LN, corrected STA is pretty good (essentially just STA)
% but correctedSTAprime can be quite poor depending on the distribution of
% the stimuli (Gaussian or Ring appear to be equally bad)

% If the cell is purely linear, the two OLS estimates of the STA
% are the same up to round off error
%correctedSTA-correctedSTAprime

% sanity check to make sure I'm implementing the normal equations correctly
% (I am).
%bprime = regress(lingensigprime,Xprime);
%correctedSTAprime;

%b = regress(lingensig,X);
%correctedSTA;

% Can I get the same answer by whitening X, calculating the STA in the
% whitened space, and transforming back? You get the same answer as the
% regular STA.
% In this section, below:
% 1) whtmat is the LIGHT transformation from the CORRELATED space to the WHITE space
% 2) inv(whtmat) is the LIGHT transformation from the WHITE space to the CORRELATED space
% 3) inv(whtmat') is the MECHANISM transformation from the CORRELATED space to the WHITE space
% 4) whtmat' is the MECHANISM transformation from the WHITE space to the CORRELATED space

%[eigv,d] = eig(Xprime'*Xprime);
%whtmat = eigv*inv(sqrt(d));
%whitened_X= Xprime*whtmat; % I'm sphering by multiplying by whtmat
%whitened_STS = sum(whitened_X.*repmat(lingensigprime,1,2));
%whitened_STS*inv(whtmat)/n; % Same as STA_prime (just whitening and unwhitening)

% How about if I transform back by the "mechanism transform" instead of the
% inverse lights transform? That's less trivial.

%whitened_STS*whtmat' % <-- the mechanism transform from white space back to original (non-white) space
% which is the same at the (X'*X)^-1 corrected STA

% Here's an equation that describes the sequence of events of the
% whiten/unwhiten strategy. I'm using VDV' = X'X instead of cov(X) because
% I'm assuming that the stimuli (X) are zero-mean and I don't want to deal 
% with all the (1/n)'s (or all the 1/sqrt(n)s in the D^-.5 matrix). VD^-.5
% is still a whitening matrix even without the n's in the denominator.
%
% ((X(VD^-.5))'y)'(VD^-.5)'
%    1         2      3
% 1 = initial whitening of the stimuli
% 2 = spike triggered sum
% 3 = unwhitening the spike triggered sum
%
% (VD^-.5)(X(VD^-.5))'y (reversing the order of operations and losing some
% transposes)
% (VD^-.5)(VD^-.5)'X'y (again but on the second factor only)
% (VD^-.5)(D^-.5V')X'y
% VD^-1V'X'y
% (X'X)^-1X'y


% ----------------------------------------
% in linear regression, we're basically calculating an STA and then
% transforming it by a matrix multiplication (inv(X'*X)) to "undo" the
% correlation in the stimulus. How different is this matrix from the
% inv(M') matrix that we  multiply the RGB STA by to get to cone weights?


%% Section 2
% K033108004 is a white noise file
stro=nex2stro(findfile('K033108004'))
gammaTable = stro.sum.exptParams.gamma_table;
gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
gammaTable1 = interp1(linspace(0,255,256),gammaTable,linspace(0,255,65536), 'spline');
invgamma = InvertGamma(gammaTable, 0);

% Reconstructing the M matrix
fundamentals = stro.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = stro.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

% % Getting the background rgb/lms
% ridx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_r'));
% gidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_g'));
% bidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_b'));
% bkgndRGB = [mode(stro.trial(:,ridx)), mode(stro.trial(:,gidx)), mode(stro.trial(:,bidx))];
% bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
% bkgndlms = M*bkgndrgb;
% 
% The covariance of the gun noise in cone space should be M*I*M' = M*M'
% % empirical sanity checking
% M*M'
% x = normrnd(0,1,1000,3); % RGB gun noise
% lmsnoise = transpose(M*x');
% cov(lmsnoise) 

convertmat1 = inv(M*M') % linear regression
convertmat2 = inv(M')

mkbasis(convertmat1)
mkbasis(convertmat2)
% They're close but not identical. They're doing two slighty different
% things. inv(M*M') is compensating for non-whiteness in the stimulus
% distribution (no changes in space here). inv(M') is taking an estimate of
% a mechanism in one space and putting it into another space.
% ----------------------------------------
%%
% Trying to general uniformly distrtibuted random points on an ellipse
n = 250;
a = 10;
b = 1;
% Using the thing that Zack found on Mathematica (that I still don't
% understand) - Points are *not* uniform on the ellipse, they are are 
% uniform on the whitened ellipse. This is just taking a multivariate
% gaussian and projecting it on to an isoprobability ellipse. Lame.
R = mvnrnd(zeros(1, 2), diag([a^2 b^2]), n);
XY = bsxfun(@rdivide, a*b*R, sqrt(R.^2 * [b^2; a^2]));
x = XY(:,1);
y = XY(:,2);

plot(x,y,'.');
axis equal;

%%
% Confirming that whiteing, response-triggered average (dividing by total
% number of stimuli, not the number of spikes) and then
% unwhitening via the mechanism transform is the same as linear regression.
% Remember, whitening doesn't just mean diagonalizing the covariance matrix
% but making it the *identity* matrix. Stimulus distribution has to have
% mean 0 for this to work.

n = 20;
mu = 20;
sigma = 10;
mixmat = normrnd(0,1,2,2);
x = normrnd(mu,sigma,n,2)*mixmat;
x = x-repmat(mean(x),n,1); % Try to relax this 
%x = [x, ones(n,1)];
y = unidrnd(2,n,1)-1; % response is completely unrelated to stimulus, but that's OK.
b_reg1 = regress(y,x)

[v,d] = eig(cov(x,1));
sqrtd = diag(1./sqrt(diag(d)));
whtmat = v*sqrtd;
%whtmat = sqrtm(inv(cov(x)));
x_prime = x*whtmat;

% sanity check
%cov(x_prime)

STS = x_prime'*y;
STA = STS./n;
b_reg2 = whtmat*STA

% Alternatively, bet we could whiten using the under centered second moment
% and it would work for non-zero mean data. Hmmm.. this didn't work. Not
% worrying about it for the moment.

x = normrnd(mu,sigma,n,2)*mixmat;
[v,d] = eig(x'*x);
sqrtd = diag(1./sqrt(diag(d)));
whtmat = v*sqrtd;
x_prime = x*whtmat;
cov(x_prime);
(x_prime'*x_prime);
b_reg1 = regress(y,x)
STA = (x_prime'*y)./n;
b_reg2 = whtmat*STA