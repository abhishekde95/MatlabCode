% Code for Mike Shadlen's 5/11/11 Neubeh 503 hw
clear all
close all

% Data
NoStim.coh = [-50 -40 -30 -20 -10 0 10 20 30 40 50];
NoStim.nupchoice = [4 3 3 6 14 20 30 39 45 46 47];
NoStim.n = [50 50 50 50 50 50 50 50 50 50 50];
Stim.coh = NoStim.coh
Stim.nupchoice = [6 11 16 25 38 40 47 46 50 48 50];
Stim.n = NoStim.n

%% Question 1
disp('Question 1')
conditions = 11;
trials_per_condition = 100;
total_trials = conditions * trials_per_condition;
view_prd = 1; %in s
rew_time = 1; %in s
ITI = 1.5; %in s
time_per_trial = rew_time+view_prd+ITI;
total_sec = total_trials*time_per_trial; %in s
total_min = total_sec/60


%% Question 2
disp('Question 2')
neu_per_mm3 = 100000;
mm3 = 1;
um3 = mm3 * 1000^3;
neu_per_um3 = neu_per_mm3/um3
r = 100;
sphere_vol = .75*pi*r^3;
neurons = neu_per_um3*sphere_vol


%% Question 3
disp('Question 3')
stim_total_upchoices = sum(Stim.nupchoice)
nostim_total_upchoices = sum(NoStim.nupchoice)


%% Question 4
disp('Question 4')
correct_up = NoStim.nupchoice(NoStim.coh==40);
correct_down = NoStim.n(NoStim.coh==-40) - NoStim.nupchoice(NoStim.coh==-40);
c40_ns_correct = correct_up + correct_down


%% Question 5
disp('Question 5')
correct_up = Stim.nupchoice(Stim.coh==40);
correct_down = Stim.n(Stim.coh==-40) - Stim.nupchoice(Stim.coh==-40);
c40_s_correct = correct_up + correct_down


%% Question 6
disp('Question 6')
correct_stim = fliplr(Stim.n(sign(Stim.coh)==-1) - Stim.nupchoice(sign(Stim.coh)==-1));
correct_stim = [50 (correct_stim + Stim.nupchoice(sign(Stim.coh)==1))]/100;
correct_nostim = fliplr(NoStim.n(sign(NoStim.coh)==-1) - NoStim.nupchoice(sign(NoStim.coh)==-1));
correct_nostim = [50 (correct_nostim + NoStim.nupchoice(sign(NoStim.coh)==1))]/100;

figure(1);clf;hold on;grid on;
title('Proportion Correct Choice (No Stimulation)')
xlabel('% Coherence')
ylabel('Proportion Correct')
plot(0:10:50,correct_stim,'ro--')
plot(0:10:50,correct_nostim,'bo--')
legend('Stimulated Trials','Non-Stimulated Trials','location','SouthEast')


%% Question 8
disp('Question 8')
p_ns = (NoStim.nupchoice./NoStim.n);
logit_pns = log(p_ns./(1-p_ns))
p_s = (Stim.nupchoice./Stim.n);
logit_ps = log(p_s./(1-p_s))

figure(2);clf;hold on;grid on;
title('Log Odds Ratios vs. Motion Coherence/Direction')
xlabel('Motion Coherence/Direction')
ylabel('Log Odds Ratio')
plot(Stim.coh,logit_pns,'bo')
plot(Stim.coh,logit_ps,'ro')
plot(Stim.coh,.07*Stim.coh,'b--')
plot(Stim.coh,1.5+.07*Stim.coh,'r--')
legend('No Stimulation','Stimulation','location','SouthEast')


%% Question 9
disp('Question 9')
y_ns = exp(logit_pns)./(1+exp(logit_pns))
y_s  = exp(logit_ps)./(1+exp(logit_ps))


%% Question 10
disp('Question 10')
p = .5;
n = 50;
k = 20;
prob_20 = binopdf(k,n,p)


%% Question 11
disp('Question 11')
prob_47 = binopdf(47,50,.9)
prob_40 = binopdf(40,50,.9)
prob_50 = binopdf(50,50,.9)


%% Question 12
disp('Question 12')
prob_48 = binopdf(48,50,.9)
prob_28 = binopdf(28,50,.5)
prob_8  = binopdf(8,50,.1)