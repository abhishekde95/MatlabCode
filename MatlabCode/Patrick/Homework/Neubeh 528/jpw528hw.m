%% This code is for the completion of homework assignment No 1 in Pbio 528
close all;
clear all;
load c1p8.mat

%% Question 1 (Chapter 1, Exercise 1)

% User-Defined Variables
Hz   = 100;
time = 10;  % In seconds
bins = .1;  % In seconds

% Auto-Defined Variables
ms   = time*1000;
t    = zeros(1,ms);
sumn = nan(bins*1000,1);

% Generate Data
for i=1:length(t)-1
    t(i+1) = (t(i) - log(rand)/Hz);
end

% Specify window of spikes that we're interested in
t = t(t<=time & t>0);

% Bin Data to Determine Fano Factor
for i=1:bins*1000
    bin = t(t<=i*bins & t>(i-1)*bins);
    sumn(i) = length(bin);
end

%Results
co_of_var = std(diff(t))/mean(diff(t))
Fano_fac = var(sumn)/mean(sumn)

figure(1); hold on;
hist(diff(t),100)
title('Histogram of Interspike Intervals')
xlabel('Interspike Intervals (s)')
ylabel('Number of Occurances')


%% Question 2 (Chapter 1, Exercise 8)
clear all
load c1p8.mat

% User-Defined Variables
stim_range = 300; % In ms

% Auto-Defined Variables
range2 = stim_range/2; % Data is in 2 ms bins
spikes2 = find(rho==1);
spikes2 = spikes2(spikes2>range2);
stimmat2 = nan(length(spikes2),range2);

% Records the stimulus at set time intervals before a spike
for s=1:length(spikes2)
    for t=1:range2
        stimmat2(s,range2+1-t) = stim(spikes2(s)-t);
    end
end

% Calculate mean stimulus at each time point preceeding a spike
sta = mean(stimmat2);


% Results
figure(2); clf, hold on;
plot(-range2*2:2:-2,sta)
title('Single-Spike-Triggered Average')
xlabel('Time Preceeding Spike (ms)')
ylabel('Mean Stimulus Value')



%% Question 3 (Chapter 1, Exercise 9)

% User-Defined Variables
stim_range    = 300; % In ms
pairing_range = 100; % In ms

% Auto-Defined Variables
stimrange3 = stim_range/2; % Data is in 2 ms bins
pairrange3 = pairing_range/2; % Data is in 2 ms bins
pairedsta = nan(pairrange3,stimrange3);
twostimcomb = pairedsta;
diffs = nan(1,pairrange3);

% Record the stimulus at set time intervals before particular spike pairings
% Calculate the mean stimulus at each time point preceeding the spike pair
for p = 1:pairrange3
    spikes3 = find(conv(rho,[1 zeros(1,p-1) 1])==2)-(1+p); %1+p to identify the first spike of the pair
    spikes3 = spikes3(spikes3>stimrange3);
    stimmat3 = nan(length(spikes3),stimrange3);
    for s=1:length(spikes3)
        for t=1:stimrange3
            stimmat3(s,stimrange3+1-t) = stim(spikes3(s)-t);
        end
    end
    pairedsta(p,:) = mean(stimmat3);
end

% Plot Paired STAs
figure(3); clf, hold on;
plot(-stimrange3*2:2:-2,pairedsta)
title('Paired-Spike-Triggered Averages')
xlabel('Time Preceeding Spike (ms)')
ylabel('Mean Stimulus Value')

% Construct Single-STAs separated by various intervals
% Compare Added Single-STAs to Paired-STAs
for dt = 1:pairrange3
    twostimcomb(dt,:) = sta(1:end) + [sta(dt:end) zeros(1,dt-1)];
    diffs(dt) = sum(twostimcomb(dt,:)-pairedsta(dt,:));
end

% Plot Added Single-STAs and Paired-STAs 
figure(4);hold on;
title('Added Single- and Paired-Spike-Triggered Averages')
xlabel('Time Preceeding Spike (ms)')
ylabel('Average Stimulus')
plot(-stimrange3*2:2:-2,twostimcomb,'b')
plot(-stimrange3*2:2:-2,pairedsta,'k')
legend('Added Single-STA','Paired STA','Location','Northwest')

% Plot Difference between Added Single- and Paired-STA's
figure(5);clf;hold on;
title('Single- vs. Paired-Spike-Triggered Averages')
xlabel('Time Between Paired Spikes (ms)')
ylabel('Magnitude of Divergence Between Spike Triggered Averages')
plot(1:2:2*pairrange3,diffs)


%% Question 4
clear all
load c1p8.mat

% User-Defined Variables
stim_range = 300; % In ms
window     = 1:5000; % A range in ms, eg. 1:500

% Auto-Defined Variables
testrho = rho(length(rho)*.8:end);
teststim = stim(length(stim)*.8:end);
trainrho = rho(1:length(rho)*.8);
trainstim = stim(1:length(stim)*.8);
stimrange4 = stim_range/2; % Data is in 2 ms bins
spikes4 = find(trainrho==1);
spikes4 = spikes4(spikes4>stimrange4);
stimmat4 = nan(length(spikes4),stimrange4);
trainmodel4 = nan(1,length(trainrho));
testmodel4 = nan(1,length(testrho));

% Records the stimulus at set time intervals before a spike
for s=1:length(spikes4)
    for t=1:stimrange4
        stimmat4(s,stimrange4+1-t) = trainstim(spikes4(s)-t);
    end
end

% Calculate mean stimulus at each time point preceeding a spike
sta4 = mean(stimmat4);

% Results
figure(6); clf, hold on;
plot(-stimrange4*2:2:-2,sta4)
title('Training Spike-Triggered Average')
xlabel('Time Preceeding Spike (ms)')
ylabel('Mean Stimulus Value')

% Convolve STA with stimuli
avgfrtrain4 = sum(trainrho)/length(trainrho);
kernel = fliplr((avgfrtrain4*sta4)/var(trainstim));
resttrain = conv(kernel,trainstim);
avgfrest4 = sum(resttrain)/length(resttrain);
r0 = avgfrtrain4-avgfrest4;
resttrain = resttrain+r0;

% Poisson Spike Generator
for t=1:length(trainrho)
    if resttrain(t) > rand
        trainmodel4(t) = 1;
    else
        trainmodel4(t) = 0;
    end
end

% Calculate error or modeled spike train
training_error = 1-sum(trainmodel4==trainrho')/length(trainmodel4==trainrho')

% %Plot results - Sanity check to see if training data predicts itself
% x=window;
% figure(7);clf;
% subplot(2,1,1); hold on;
% title('Sample from Empirical Spike Train')
% xlabel('Time (ms)')
% ylabel('Spike')
% bar(2*x,trainrho(x))
% subplot(2,1,2);hold on;
% title('Sample from Modeled Spike Train')
% xlabel('Time (ms)')
% ylabel('Spike')
% bar(2*x,trainmodel4(x))


% Use kernel to predict test data
resttest = conv(kernel,teststim);
avgfrest4 = sum(resttest)/length(resttest);
r0 = avgfrtrain4-avgfrest4;
resttest = resttest+r0;

% Poisson Spike Generator
for t=1:length(testrho)
    if resttest(t) > rand
        testmodel4(t) = 1;
    else
        testmodel4(t) = 0;
    end
end

% Plot results
x=window;
figure(7);clf;
subplot(2,1,1); hold on;
title('Sample from Empirical Spike Train')
xlabel('Time (ms)')
ylabel('Spike')
bar(2*x,testrho(x))
subplot(2,1,2);hold on;
title('Sample from Modeled Spike Train')
xlabel('Time (ms)')
ylabel('Spike')
bar(2*x,testmodel4(x))

% Calculate error or modeled spike train
test_error = 1-sum(testmodel4==testrho')/length(testmodel4==testrho')

%% Question 5
clear all
load c1p8.mat

covar = cov(rho,stim);


disp 'My apologies, but I have no answer for problem #5.  I have spent'
disp 'many hours on it but have been unable to make much progress.  The problem,' 
disp 'I believe, is that I dont yet understand a few fundamental concepts, '
disp 'eg: Eigenvalues, eigenvectors, or eigenmodes.  Ive been reading'
disp 'supplementary sources to try to improve, but thus far, Im at a loss.'
disp 'Ill continue to work on this problem, however, and would be happy to'
disp 'share my work with you if you care to see.  Thanks, and, again, my apoligies.'



