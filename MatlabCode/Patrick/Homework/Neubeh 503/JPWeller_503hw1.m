%% This is the homework assigned by Mike Shadlen for 503 due 4/4/11
clear all
close all
load sptimesPS1_2011.mat

%% Question 1
disp 'Question1 = Left Hemisphere'


%% Question 2
disp 'Question2 = Between 2 and 4 degrees vertically, 1 and 3.5 degrees horizontally'


%% Question 3
% What is the average spike rate during the 1st 100 ms of the stimulus sweep? 
% Express your answer in spikes per second (spikes/s)

% Set Variables
ms = 100;
sprbins3 = nan(1000/ms,size(s,2));

% Calculations
for n=1:(1000/ms)
    sprbins3(n,:) = sum(s>=((n-1)*ms) & s<(ms*n));
end
meanspr3 = mean(sprbins3,2).*(1000/ms);
semspr3 = (std(sprbins3(1,:)).*(1000/ms))./sqrt(size(s,2));
Answer3_1 = meanspr3(1)
Answer3_2 = semspr3 


%% Question 4
% What is the average spike rate in the 50 ms epoch beginning at t=150 ms after the 
% stimulus was turned on?

% Set Variables
ms = 50;
sprbins4 = nan(1000/ms,size(s,2));

% Calculations
for n=1:(1000/ms)
    sprbins4(n,:) = sum(s>=(n-1)*ms & s<(ms*n));
end
meanspr4 = mean(sprbins4,2)*(1000/ms);
semspr4 = (std(sprbins4(4,:)).*(1000/ms))./sqrt(size(s,2));
Answer4_1 = meanspr4(4)
Answer4_2 = semspr4


%% Quesiton 5
% Make a graph of mean spike rate as a function of time. Use 50 msec bins so that 
% the spike rate is plotted for the epoch from 0-49.999 ms, then 50-99.999, then 
% 100-149.999, etc. Label the ordinate "Spikes/sec.

% User-Defined Variables
ms = 50;

% Set Variables
sprbins5 = nan(1000/ms,size(s,2));

% Calculations
for n=1:(1000/ms)
    sprbins5(n,:) = sum(s>=(n-1)*ms & s<(ms*n));
end
meanspr5 = mean(sprbins5,2)*(1000/ms);
sespr5  = std(sprbins5,0,2)./sqrt(10);

% Plot
figure(1); hold on; grid on;
axis([0 1000 0 max(meanspr5)])
title('Question 5')
ylabel('Spikes/s')
xlabel('Time From Onset of Movement (ms)')
plot(25:50:975,meanspr5)
errorbar(25:50:975,meanspr5,sespr5,'o')


%% Quesiton 6
% What is the average spike rate in the last 500 ms of the stimulus presentation?

% User-Defined Variables
ms = 500;

% Set Variables
sprbins6 = nan(1000/ms,size(s,2));

% Calculations
for n=1:(1000/ms)
    sprbins6(n,:) = sum(s>=(n-1)*ms & s<(ms*n));
end
meanspr6 = mean(sprbins6,2)*(1000/ms);
semspr6 = (std(sprbins6(2,:)).*(1000/ms))./sqrt(size(s,2));
Answer6_1 = meanspr6(2) 
Answer6_2 = semspr6


%% Question 7
% In the epoch from 100 to 400 ms, how many spikes were emitted on each of the 
% ten trials?

% User-Defined Variables
ms = 100;

% Set Variables
sprbins7 = nan(1000/ms,size(s,2));

% Calculations
for n=1:(1000/ms)
    sprbins7(n,:) = sum(s>=(n-1)*ms & s<(ms*n));
end

spikesinepoch7 = sum(sprbins3(2:4,:));
Answer7 = spikesinepoch7


%% Question 8
% What is the mean and the variance of the spike counts from these 10 trials. Use 
% the counts, not the spike rate.

spcount8 = sum(~isnan(s));
Answer8 = mean(spcount8),var(spcount8)


%% Question 9
% What is the mean and variance of the spike count in the last 500 ms of the trial

spcount9 = sum(s>500);
Answer9 = mean(spcount9),var(spcount9)


%% Question 10
disp 'Question10 = 9'

%% Question 11
disp 'Question11 = 3'

%% Question 12
disp 'Question12 = 0.6'

%% Question 13
disp 'Question13 = 0.3'

%% Question 14
disp 'Question14 = 2 trials'

%% Question 15

% Variables
lam = 1.9;
k   = 0;

% Calculations
poissprob15 = (lam^k * exp(-lam)) / factorial(k);
Question15 = poissprob15

