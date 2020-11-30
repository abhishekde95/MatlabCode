% J Patrick Weller
% PBio 545 Homework Assignment #4
% Created   2/13/12     JPW

clear all
close all

%% Question 1

% (1) play with the rate constant alpha and explain what happens.

% The update rule in the eqation x(n) = x(n-1) + timestep * (y(n-1) - alpha
% * (x(n-1) ) relies on the difference between two values: y at some point
% n and x at some point n.  Alpha, the rate constant, serves to scale x(n-1),
% which in turn enlarges or diminishes the difference between y(n-1) and
% x(n-1).  The greater the difference between x(n-1) and y(n-1), the greater 
% the effect on x(n).  Thus, the point where the difference between x and y 
% is the greatest (at n=100, where y(n) = 1 and x(n) = 0), the slope of
% x is also the greatest.


% (2) try different shaped inputs (other than step) and explain the results.

figure(1);
tx(1) = 0;				% initial condition
alpha = 20;				% rate constant for decay of x (in 1/sec)
TimeStep = 0.001;		% time step for difference equation (in sec)
PrePts = 200;			% points before step in y
StmPts = 400;			% number of time points that y is active
NumPts = 1000;			% total points to simulate

% Upward Ramp
ty(1:PrePts) = 0;
ty(PrePts+1:PrePts+StmPts) = linspace(0,1,StmPts);
ty(StmPts+PrePts+1:NumPts) = 0;
y(1,:) = ty;

% Downward Ramp
ty(1:PrePts) = 0;
ty(PrePts+1:PrePts+StmPts) = linspace(1,0,StmPts);
ty(StmPts+PrePts+1:NumPts) = 0;
y(2,:) = ty;

% Gaussian
ty(1:StmPts) = 0;
ty(PrePts+1:PrePts+StmPts) = normpdf([1:StmPts],(StmPts)/2,(StmPts)/8);
ty(PrePts+1:PrePts+StmPts) = (ty(PrePts+1:PrePts+StmPts)./max(ty));
ty(StmPts+PrePts+1:NumPts) = 0;
y(3,:) = ty;

% Comb Function
ty(1:PrePts+StmPts) = 0;
ty(PrePts+1:25:PrePts+StmPts) = 1;
ty(StmPts+PrePts+1:NumPts) = 0;
y(4,:) = ty;

% Sinusoid
ty(1:StmPts) = 0;
npi = 20;
ty(PrePts+1:PrePts+StmPts) = cos((npi*pi/StmPts:npi*pi/StmPts:npi*pi)+pi)./2 +.5;
ty(StmPts+PrePts+1:NumPts) = 0;
y(5,:) = ty;

% Plot y
q = 1;
tme = 1:NumPts;
tme = (tme - PrePts) * TimeStep;
for n = 1:5
    subplot(5, 2, q); q= q+2;
    plot(tme, y(n,:));
    ylim([-.1 1.1]);
    xlim([min(tme) max(tme)])
    if n == 1
        title('input y');
    end
end
xlabel('time (sec)');


% Plot x
q=2;
for neq = 1:5
    for pnt = 2:NumPts
        tx(pnt) = tx(pnt-1) + (y(neq,pnt-1) - alpha * tx(pnt-1)) * TimeStep;
    end
    x(neq,:) = tx;
    subplot(5, 2, q); q = q+2;
    plot(tme, x(neq,:));
    xlim([min(tme)-max(tme)*.1 max(tme)*1.1])
    ylim([min(x(neq,:))-max(x(neq,:))*.1 max(x(neq,:))*1.1])
    if neq == 1
        title('output x');
    end
end
xlabel('time (sec)');


% (3) play with TimeStep to see over what range of time bins the numerical solution
%	  is accurate (i.e. compare to an analytical solution to the same DEQ).

clear y x tme

figure(2);
rc = .02;			% rate constant for decay of x (in sec)
TimeStep = [0.0001 0.001 0.01 0.1];		% time step for difference equation (in sec)
Pret = .2;			% time before step in y (in s)
Stmt = .5;			% time that y is active (in s)
Totalt = 1;			% total time to simulate (in s)
alpha = rc ./ TimeStep; %rate constant for decay of x (in 1/sec)
PrePts = Pret ./ TimeStep;
StmPts = Stmt ./ TimeStep;
NumPts = Totalt ./ TimeStep;
tme = cell(numel(TimeStep),1);
y = cell(numel(TimeStep),1);

% Plot y
for n = 1:numel(TimeStep)
    tme{n} = 1:NumPts(n);
    tme{n} = (tme{n} - PrePts(n)) * TimeStep(n);
    y{n}(1:PrePts(n)) = 0;
    y{n}((PrePts(n)+1) : (PrePts(n)+StmPts(n))) = 1;
    y{n}(PrePts(n) + StmPts(n) + 1:NumPts(n)) = 0;
end
subplot(numel(TimeStep),2,1:2:numel(TimeStep)*2)
plot(tme{1},y{1})
ylim([-.1 1.1]);
xlim([min(tme{1}) max(tme{1})])
title('Input y')
xlabel('time (sec)');


% Plot x
q = 2;
for n = 1:numel(TimeStep)
    tx = nan(NumPts(n),1);
    tx(1) = 0;
    for pnt = 2:NumPts(n)
        tx(pnt) = tx(pnt-1) + (y{n}(pnt-1) - alpha(n) * tx(pnt-1)) * TimeStep(n);
    end
    x{n} = tx;
    subplot(numel(TimeStep), 2, q);q=q+2;
    plot(tme{n}, x{n});
    xlim([min(tme{n})-max(tme{n})*.1 max(tme{n})*1.1])
    ylim([min(x{n})-max(x{n})*.1 max(x{n})*1.1])
    if n == 1
        title('output x');
    end
    ylabel(['Time Step = ' num2str(TimeStep(n))])
end
xlabel('time (sec)');


%% Questions/Problems #2:
% (1) How does the solution differ in this case from that above?  Why?  

% The solution differs in this case only by some constant offset.  That is,
% the total modulation of the output determined by the input is equal (.05 
% units above baseline), but the baseline is different in each case.  In
% the first, there is no spontaneous activity, so the output is 0 during
% the pre-stimulus time.  In the present case, there is spontaneous
% activity from the output, so the stimulus modulates the output from some
% non-zero value.  

% (2) Can you explain the change from the differential equation?

% The change, I take it, is that we see some constant non-zero output when there is no
% input?  This happens for two reasons: 1. We see an additive term because we have a 
% linear equation to which some constant, eta, is being added (and scaled);
% and 2. This term, eta, keeps the output activity constant during the
% pre-stimulus phase because we begin at the value to which this equation
% tends (eta/alpha).  That is, if the input y is 0, the output x stays 
% constant when the previous output (the state of the output at one time
% step prior) is equal to eta/alpha (since we want eta - alpha * x(n-1) to
% be 0).

% (3) Why is eta/alpha is reasonable initial condition?  

% This is a resonable initial condition because it is the only one that
% does not evoke some immediate change in the output at x(2). In other 
% words, we want our output to be constant during the pre-stimulus phase so
% that we can confidently attribute any change in the output to be the
% result of the stimulus, and not the result of some prior activity.  If
% x(1) is not equal to eta/alpha, then the output will be modulating in the
% pre-stimulus phase toward the steady state eta/alpha.


%% Questions/Problems #3:

% (1) How does the behavior compare to the case without feedback?
%	  Compare both the amplitude and kinetics of x.

% Both the amplitude and kinetics of x have changed using this new
% equation. Because the output x is being scaled at each step by some
% term (1 + gx).^Power, and the Power term is less than 0 (in this case, it 
% is -4), the output x is always less than the original output (the output
% without the feedback mechanism).  If the Power term were greater than
% zero, the output would be greater than the original (without feedback).
% The feedback mechanism also changes the kinetics of the output: when the
% feedback is negative, a mechanism with power of a greater magnitude seems 
% to asymptote more quickly than a muchanism with power of a smaller 
% magnitude.  If the feedback mechanism raises the term to a positive
% power, the output quickly becomes enormous.  In my example below, I use
% only values between 0 and 1, for outputs with feedback mechanisms raised
% to a power greater than one cannot be seen on the same scale as those
% without a feedback mechanism.  Raised to a large enough power (somewhere
% between 1 and 2), the output does not even come close to its asymptotic
% value over this time range.


% (2) How does effect of feedback change as power changed?  Why? Try both negative
%	  and positive values for the power.

% The effect of the feedback is dependent on both the sign and magnitude of
% the power because it ultimately scales the input up or down.  That is, at
% each point in time, the input y is mulitplied by some term (1 + g * x(n-1))
% raised to some power.  Since neither the output x nor the feeback gain g
% can be negative by convention (to make it a more realistic output of a 
% cell), the smallest value that this term can take is 1.  If x has some
% positive value and g is non-zero, this whole term takes some value
% greater than 1.  When raised to a negative power, this whole term will
% become small (less than 1), and scale the input y down.  If the power to
% which this term is raised is positive, it quickly becomes very large and
% scales the input y up.  If this whole term is raised to the 0 power, it
% becomes 1 and the input y remains as it is.

clear x y tme

figure(3); clf
alpha = 10;
TimeStep = 0.001;		% time step for difference equation (in sec)
PrePts = 200;			% points before step in y
StmPts = 400;			% number of time points that y is active
NumPts = 1000;			% total points to simulate
NegPower = -4:0;						% feedback power
PosPower = 0:.2:.6;
xN = nan(numel(NegPower),NumPts);
xP = nan(numel(PosPower),NumPts);
xN(:,1) = 0;						% initial condition
xP(:,1) = 0;
g = 20;							% feedback gain
tme = 1:NumPts;
tme = (tme - PrePts) * TimeStep;

% initialize y; in this case y is a simple step
y(1:PrePts) = 0;
y(PrePts+1:PrePts + StmPts) = 1;
y(PrePts + StmPts + 1:NumPts) = 0;

% plot input y
subplot(4, 1, 1)
plot(tme, y);
ylim([-.1 1.1]);
xlabel('time (sec)');
ylabel('input');
title('Activation of x by y with feedback');

% No Feedback
xnfb = nan(1,NumPts);
xnfb(1) = 0;
for pnt = 2:NumPts
    xnfb(pnt) = xnfb(pnt-1) + (y(pnt-1) - alpha * xnfb(pnt-1)) * TimeStep;
end

% Negative Powers
for pnt = 2:NumPts
    xN(:,pnt) = xN(:,pnt-1) + (y(pnt-1) * (1 + g * xN(:,pnt-1))'.^NegPower' - alpha * xN(:,pnt-1)) .* TimeStep;
end
xN = cat(1,xN,xnfb);

% Positive Powers
for pnt = 2:NumPts
    xP(:,pnt) = xP(:,pnt-1) + (y(pnt-1) * (1 + g * xP(:,pnt-1))'.^PosPower' - alpha * xP(:,pnt-1)) .* TimeStep;
end
xP = cat(1,xP,xnfb);

% Plot Negative Powers
subplot(4, 1, 2); hold on;
plot(tme, xN(1:end-1,:));
plot(tme(1:20:end), xN(end,1:20:end),'ko')
for n = 1:numel(NegPower)
    legN{n} = ['Power = ' num2str(NegPower(n))];
end
legN{n+1} = 'No Feedback';
legend(legN);
ylim([0 1.1*max(max(xN))])
title('Negative Exponential Powers')
xlabel('time (sec)');
ylabel('output');

% Plot Positive Powers
subplot(4, 1, 3); hold on;
plot(tme, xP(1:end-1,:));
plot(tme(1:20:end), xP(end,1:20:end),'ko')
for n = 1:numel(PosPower)
    legP{n} = ['Power = ' num2str(PosPower(n))];
end
legP{n+1} = 'No Feedback';
legend(legP);
ylim([0 1.1*max(max(xP))])
title('Positive Exponential Powers')
xlabel('time (sec)');
ylabel('output');


% (3) How does changing g change things?  Why?

% The term g is the gain on the feedback - that is, it is a measure of how
% influential the preceding output is on the present output.  When g is
% large, the previous values of x will be heaily influential; when g is
% small, the previous values of x will be only mildly influential.  As I
% mentioned previously, g should remain positive if we wish to make use
% of it as a description of a real neuron.  Below is an example where
% larger values for g result in smaller outputs x.  This is because we are
% raising this whole term to a negative value (here, as in the original, to
% -4); thus, the larger the value inside the parentheses, the smaller the
% scaling value of input y.

g = 0:5:20;
xg = nan(numel(g),NumPts);
xg(:,1) = 0;
Power = -4;

for pnt = 2:NumPts
    xg(:,pnt) = xg(:,pnt-1) + (y(pnt-1) * (1 + g' .* xg(:,pnt-1)).^Power - alpha * xg(:,pnt-1)) .* TimeStep;
end

subplot(4,1,4); hold on
plot(tme,xg)
plot(tme(1:20:end), xnfb(end,1:20:end),'ko')
for n = 1:numel(g)
    legg{n} = ['Gain = ' num2str(g(n))];
end
legg{n+1} = 'No Feedback';
legend(legg);
ylim([0 1.1*max(max(xg))])
title('Various Feedback Gains')
xlabel('time (sec)');
ylabel('output');


% (4) What happens when you change TimeStep?  Is the numerical solution more 
%	  or less sensitive to changes in TimeStep than that above?  Why?

% Changing TimeStep should have the same influence over the output as
% above.  When you change the sampling of the equation - that is, how often
% the output is influenced by the input, the solution becomes more or less
% sensitive.  Above, the input was only y, and the larger the sampling
% rate of y, the more sensitive the output x became to y.  Here, the input
% is not only y, but the previous values of x.  The sensitivity still
% increases as we increase the sampling rate.


%% Questions/Problems #4:

% (1) why do we choose an exponential decay for the shape of rhodopsin's 
% activity?  What might change that?

% I assume that we've chosen an exponential decay for the shape of
% rhodopsin because this accurately reflects the actual activation of 
% rhodopsin in the eye by a single photon.  Several photons activating
% rhodopsin in close succession might change the shape of this activity.

% (2) Explore different combinations of sigma and phi and their impact.

% When sigma is large, the rhodopsin activity dissipates quickly.  The more
% quickly this activity approaches 0, the less effect it has on the
% phosphodiesterase.  When phi, the phosphodiesterase activity decay rate,
% is large, the effect of rhodopsin on phe is weakened, and less activity
% in a given time is generated in the phe.  Not only does does phe activity
% approach its baseline more quickly for larger values of phi, but both its
% baseline and peak activity is smaller.

clear p r

figure(4); clf
sigma = [1 5 15];          % rhodopsin activity decay rate constant (1/sec)
phi = 4.7:.3:5.3;		% phosphodiesterase activity decay rate constant (1/sec)
eta = 10;				% phosphodiesterase activation rate constant (1/sec)
p(:,:,1) = repmat(eta./phi,numel(sigma),1);			% initial condition for p
r = ones(numel(sigma),1);     % initial condition for r

NumPts = 1000;			% number of points to simulate
TimeStep = 0.001;		% time step 
tme = 1:NumPts;
tme = tme * TimeStep;

% solve difference equation
for pnt = 2:NumPts
    r(:,pnt) = r(:,pnt-1) + TimeStep * (-sigma' .* r(:,pnt-1));
    p(:,:,pnt) = p(:,:,pnt-1) + TimeStep * (repmat(r(:,pnt-1),1,numel(phi))...
        + eta - repmat(phi,numel(sigma),1) .* p(:,:,pnt-1));
end

% plot time course of rhodopsin activity
q = 1;
for n = 1:numel(sigma)
    subplot(numel(sigma),2,q);
    plot(tme, r(n,:));
    if n == ceil(numel(sigma)/2)
        ylabel('rhodopsin activity')
    end
    legend(['sigma = ' num2str(sigma(n))])
    q = q+2;
end
xlabel('time (sec)')

% plot time course of phosphodiesterase activity
q = 2;
for n = 1:numel(phi)
    leg{n} = ['phi = ' num2str(phi(n))];
end
for n = 1:numel(sigma)
    subplot(numel(sigma),2,q);
    pplot = shiftdim(p(n,:,:),1);
    plot(tme, pplot);
    ylim([min(min(pplot))*.9 max(max(pplot))*1.15])
    if n == ceil(numel(sigma)/2)
        ylabel('pde activity')
    end
    legend(leg,'Location','NorthEast')
    q = q+2;
end
xlabel('time (sec)')


%% Questions/Problems #5:

% (1) Explain the relation between the steady state conditions
%	  and the constants q and smax.

% We use the steady state conditions determine both q and smax because we
% are beginning our simulation in darkness, where the system has reached
% some equilibrium point - a state in which none of our variables (current,
% pde activity, cGMP, etc) are not changing: a steady state. q determines
% our initial quantity of current, which, in the dark, is determined by the
% amounts of calcium and cGMP present in the cell: cdark and gdark.  smax,
% the maximum rate of synthesis of cGMP, is also determined by the steady
% states of calcium and cGMP.  Since the rate of cGMP synthesis in the cell
% is inversely proportional to the amout of calcium in the cell, we can
% discover the their proprotionality constant by taking the state when it
% is at rest - when neither the rate of cGMP synthesis or calcium
% concentration is changing.

% (2) Play with the various parameters and see how they alter the calculated
%	  light response.  Explain why things change they way they do?

clear p r
figure(5)

for n = 1:5
    
    sigma = 5;				% rhodopsin activity decay rate constant (1/sec)
    phi = 5;				% phosphodiesterase activity decay rate constant (1/sec)
    eta = 10;				% phosphodiesterase activation rate constant (1/sec)
    gdark = 15;				% concentration of cGMP in darkness
    if n == 2
        gdark = 2 * gdark;
    elseif n == 3
        gdark = .5 * gdark;
    end
    cgmp2cur = 8e-3;		% constant relating cGMP to current
    cdark = 0.5;			% dark calcium concentration
    if n == 4
        cdark = 2 * cdark;
    elseif n == 5
        cdark = .5 * cdark;
    end
    beta = 2;				% rate constant for calcium removal in 1/sec
    hillcoef = 4;			% cooperativity
    hillaffinity = 0.3;		% affinity
    cur2ca = beta * cdark / (cgmp2cur * gdark^3);				% get q using steady state
    smax = eta/phi * gdark * (1 + (cdark / hillaffinity)^hillcoef);		% get smax using steady state
    
    NumPts = 1000;			% number of points to simulate
    TimeStep = 0.001;		% time step
    tme = 1:NumPts;
    tme = tme * TimeStep;
    
    
    % initial conditions
    g(1) = gdark;
    s(1) = gdark * eta/phi;
    c(1) = cdark;
    p(1) = eta/phi;
    r(1) = 1;				% initial condition for r
    
    
    % solve difference equations
    for pnt = 2:NumPts
        r(pnt) = r(pnt-1) + TimeStep * (-sigma * r(pnt-1));
        p(pnt) = p(pnt-1) + TimeStep * (r(pnt-1) + eta - phi * p(pnt-1));
        c(pnt) = c(pnt-1) + TimeStep * (cur2ca * cgmp2cur * g(pnt-1)^3 - beta * c(pnt-1));
        s(pnt) = smax / (1 + (c(pnt) / hillaffinity)^hillcoef);
        g(pnt) = g(pnt-1) + TimeStep * (s(pnt-1) - p(pnt-1) * g(pnt-1));
    end
    % determine current change
    cur = cgmp2cur * g.^3;
    
    % plot current, pde, synthesis, cGMP and calcium
    subplot(5, 5, 1 + 5*(n-1));
    plot(tme, cur);
    xlabel('time (sec)');
    ylabel('current');
    subplot(5, 5, 2 + 5*(n-1));
    plot(tme, p);
    xlabel('time (sec)');
    ylabel('pde activity');
    subplot(5, 5, 3 + 5*(n-1));
    plot(tme, s);
    xlabel('time (sec)');
    ylabel('synthesis rate');
    subplot(5, 5, 4 + 5*(n-1));
    plot(tme, g);
    xlabel('time (sec)');
    ylabel('[cGMP]');
    if n == 1 || n == 2 || n == 3
        legend(['gdark = ' num2str(gdark)])
    end
    subplot(5, 5, 5 + 5*(n-1))
    plot(tme, c)
    xlabel('time (sec)');
    ylabel('[calcium]');
    if n == 1 || n == 4 || n == 5
        legend(['cdark = ' num2str(cdark)])
    end
    
end

% In the first row, I have reproduced the simulation from the tutorial with
% normal values for darkc and darkg.  In each of the following rows, I have
% altered either darkc or darkg, first by doubling their values, then my
% halving them.  In the second row I have simulated the
% effect of the steady state amount of cGMP in the cell as twice its
% initial value from row 1.  Since darkg is twice its normal value, its 
% resting rate of synthesis must also be twice its original value.  And
% since the current in cell is also directly proportional to the amount of
% cGMP, the steady state of the current is also increased (by quite a
% bit!).  In the second row, the resting amount of cGMP is half of its
% original value, making the rate of snythesis also half of its original
% rate, and the current is reduced (again, by quite a bit).  In the fourth
% row, the amount of calcium present in the cell in its resting state is
% double its original value.  Changing the resting amount of calcium in the
% cell also changes the relationship between current and calicum (through
% smax), so when we increase the resting state of calcium in the cell, we
% also increase the effect of a change in current on that amount of calcium;
% thus, our amount of calcium drops by twice the level as it originally
% does, and, in turn, increases the maximum rate of cGMP synthesis.  When we
% half the resting amount of calcium from its original value, we also
% reduce the effect that a change in current has on the amount of calcium
% in the cell, which in turn reduces the rate of cGMP synthesis, and
% thereby the amount of cGMP in the cell.

%% (3) The model will generate damped oscillations for some values of beta.  Why?

% The model will generate damped oscillations for some values of beta
% because of the relationship between cGMP and calcium.  beta controls how
% quickly calcium leaves the cell, and, in turn, how quickly the amount of
% cGMP in the cell effects the amount of calcium.  If calcium leaves the
% cell slowly enough, it lingers long enough to raise the rate of cGMP
% synthesis to such a level that more cGMP is present in the cell than was
% present in the dark.  When this happens, the current surpasses its
% original level, as well as calcium.  When the amount of calcium suprasses
% its resting level, the rate of cGMP synthesis is slowed, and the amount 
% of current and calcium fall below their original level, and all three
% values act as a damped oscilator tending toward the resting values.  If
% beta is too high, calcium is removed from the cell at such a rate that it
% cannot bring the rate of cGMP synthesis high enough to surpass its
% resting value, and we get no oscillation.

figure(6)

betachoices = 1:2:10;

for n = 1:numel(betachoices)
    
    sigma = 5;				% rhodopsin activity decay rate constant (1/sec)
    phi = 5;				% phosphodiesterase activity decay rate constant (1/sec)
    eta = 10;				% phosphodiesterase activation rate constant (1/sec)
    gdark = 15;				% concentration of cGMP in darkness
    cgmp2cur = 8e-3;		% constant relating cGMP to current
    cdark = 0.5;			% dark calcium concentration
    beta = betachoices(n);	% rate constant for calcium removal in 1/sec
    hillcoef = 4;			% cooperativity
    hillaffinity = 0.3;		% affinity
    cur2ca = beta * cdark / (cgmp2cur * gdark^3);	% get q using steady state
    smax = eta/phi * gdark * (1 + (cdark / hillaffinity)^hillcoef);	% get smax using steady state
    
    NumPts = 2000;			% number of points to simulate
    TimeStep = 0.001;		% time step
    tme = 1:NumPts;
    tme = tme * TimeStep;
    
    
    % initial conditions
    g(1) = gdark;
    s(1) = gdark * eta/phi;
    c(1) = cdark;
    p(1) = eta/phi;
    r(1) = 1;				% initial condition for r
    
    
    % solve difference equations
    for pnt = 2:NumPts
        r(pnt) = r(pnt-1) + TimeStep * (-sigma * r(pnt-1));
        p(pnt) = p(pnt-1) + TimeStep * (r(pnt-1) + eta - phi * p(pnt-1));
        c(pnt) = c(pnt-1) + TimeStep * (cur2ca * cgmp2cur * g(pnt-1)^3 - beta * c(pnt-1));
        s(pnt) = smax / (1 + (c(pnt) / hillaffinity)^hillcoef);
        g(pnt) = g(pnt-1) + TimeStep * (s(pnt-1) - p(pnt-1) * g(pnt-1));
    end
    % determine current change
    cur = cgmp2cur * g.^3;
    
    % plot current, pde, synthesis, cGMP and calcium
    subplot(numel(betachoices), 5, 1 + 5*(n-1));
    plot(tme, cur);
    xlabel('time (sec)');
    ylabel('current');
    legend(['beta = ' num2str(beta)],'Location','SouthEast')
    subplot(numel(betachoices), 5, 2 + 5*(n-1));
    plot(tme, p);
    xlabel('time (sec)');
    ylabel('pde activity');
    subplot(numel(betachoices), 5, 3 + 5*(n-1));
    plot(tme, s);
    xlabel('time (sec)');
    ylabel('synthesis rate');
    subplot(numel(betachoices), 5, 4 + 5*(n-1));
    plot(tme, g);
    xlabel('time (sec)');
    ylabel('[cGMP]');
    subplot(numel(betachoices), 5, 5 + 5*(n-1))
    plot(tme, c)
    xlabel('time (sec)');
    ylabel('[calcium]');
    
end

