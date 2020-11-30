% J Patrick Weller
% PBio 545 Homework Assignment #3
% Created   1/31/12     JPW


%% Question #1

dt = .01;
t = 0:dt:2;
f = zeros(size(t));
f(find(t==.5):find(t==1.5)) = 1/numel(f(find(t==.5):find(t==1.5)));

cumf = nan(size(f));
for n = 1:numel(f)
    cumf(n) = sum(f(1:n));
end

hazf = f./(1-cumf);
hazf(sign(hazf)==-1)=1;

figure(1);clf;hold on;axis([0 2 0 1.1])
plot(t,f)
plot(t,cumf,'r')
plot(t,hazf,'g')
xlabel('Time(ms)')
ylabel('Probability')
legend('Probability Density Function','Cumulative Probability Function','Hazard Function','Location','NorthWest')

% If the runner has the knowledge of the hazard function that we have, it
% may very well affect her start time.  That is, if the condition arises in
% which the first gun goes off, then nearly 1.5 seconds goes by, the runner
% knows that the probability of the gun firing at t = 1.5 is equal to 1,
% and she should take off.  If t is less than 1.5, the runner has little
% extra knowledge that will be helpful.  Although she knows that the 
% probability of the gun going off steadily increases as t approaches t=1.5
% (and the gun has not already gone off), the runner cannot make use of
% this without risking a fasle start until t=1.5.


%% Question 2

% The general relationship between a probability function (f(t)) and its
% cumulative probability function (F(t)) is that the cumulative probability
% is the integral of the probability function (that is, the integral of the
% probability mass to the left on the graph, or prior in time or space).
% The probability function will tell you the probability of some event
% occuring at any particular time (or space, etc), but the cumulative 
% probability function tells you the probability of an event occuring at 
% OR BEFORE that point.  The probability function integrates to one, and
% the final point on the cumulative probability function should have a
% height of one.


%% Question 4

% The randomness in this Markov chain comes from the X = rand step and
% taking advantage of Matlab's serial execution of conditionals.  When the
% probability of transitioning into the different states is recalculated at
% each execution of the loop, the cumulative probability (p_s) is
% calculated.  Then, a random number is generated from a standard uniform
% distribution (from 0 to 1).  The probability that this random number
% falls in the range between 0 and the first value of p_s is equal to the
% probability of transitioning into the first state, just as the
% probability that the random number falls in the range between the first
% and second values of p_s is equal to the probability of transitioning
% into the second state, etc.


%% Question 5

p  = .5
q = 1-p
P = diag(p*ones(8,1),1) + diag(q*ones(8,1),-1); P(1,:) = 0; P(1,1) = 1; P(end,:)=0; P(end,end)=1
T = P'

s0 = [0 0 0 0 1 0 0 0 0]

% After 16 rounds, the probability of being any any state is:
P_16f = T^16 * s0'

% And just for fun, lets look out a few more rounds:

P_20f = T^25 * s0'
P_30f = T^30 * s0'
P_50 = T^50 * s0'

% We can see that the probability of being in the absorbing states is 
% increasing as we go futher and further out (project into more and more
% rounds).  We can also see that we have a 0 percent probability of being
% in some states after every round.  This makes sense: if we begin in state
% 5, we cannot still be in state 5 after one round, but we can return to
% that state after 2 rounds; cannot be in state 5 after 3 rounds, but
% again, can return to state 5 after 4 rounds, etc.  This goes for other
% states as well: we can be in states 4 or 6 after one round, but not after
% two; can enter them after three rounds, but not after four.  Hence, after
% 16 rounds, we cannot be in any of the even states: states 2, 4, 6, or 8.

% Now let us mess with the odds: we will use an unfair coin that favors
% Adrienne:
p  = .6
q = 1-p
P = diag(p*ones(8,1),1) + diag(q*ones(8,1),-1); P(1,:) = 0; P(1,1) = 1; P(end,:)=0; P(end,end)=1
T = P'

P_16uf = T^16 * s0'

P_16uf(1) + P_16uf(9)
% We see that the probability of being in one of the two absorbing states
% is rougly .75.  The probability of Adrienne being the winner is .62,
% while the probability of Fred being the winner is only .12!  As we
% increase the 
