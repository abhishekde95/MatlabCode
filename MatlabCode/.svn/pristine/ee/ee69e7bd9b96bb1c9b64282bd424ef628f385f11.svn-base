% This code is for Mike Shadlen's 3rd 503 assignment due 5/13/11

clear all
close all

load xposVals.mat
pos0 = y(:,1)
pos1 = y(:,2)

%% Question 1
disp('Question 1')

meansr0 = mean(pos0)*2
sem0    = std(pos0)/sqrt(length(pos0))
meansr1 = mean(pos1)*2
sem1    = std(pos1)/sqrt(length(pos1))


%% Question 2
disp('Question 2')

[passttest,p] = ttest2(pos0,pos1,.05,'both','unequal')

%% Question 3
disp('Question 3')

stdpos0 = std(pos0)
stdpos1 = std(pos1)


%% Question 4
disp('Question 4')

diffmeans = mean(pos0) - mean(pos1)

%% Quesation 5
disp('Question 5')

thresh   = 40;
probtrupos0 = sum(pos0>=thresh)/length(pos0)


%% Question 6
disp('Question 6')

probcorrej1 = sum(pos1<thresh)/length(pos1)
 
%% Question 7
disp('Question 7')

probfalpos = sum(pos1>=thresh)/length(pos1)


%% Question 8
disp('Question 8')

thresh = 35;
probtrupos35 = sum(pos0>thresh)/length(pos0)
probtruneg35 = sum(pos1<=thresh)/length(pos1)
probfalpos35 = sum(pos1>thresh)/length(pos1)
probfalneg35 = sum(pos0<=thresh)/length(pos0)


%% Question 9
disp('Question 9')

thresh = [min(pos1):1:max(pos0)+1];
probtrupos = nan(1,length(thresh));
probfalpos = probtrupos;

for i = 1:length(thresh)
    probtrupos(i) = sum(pos0>=thresh(i))/length(pos0);
    probfalpos(i) = sum(pos1>=thresh(i))/length(pos1);
end

figure(1);clf;hold on;grid on;
axis([0 1 0 1])
title('ROC Curve')
xlabel('False Alarm Rate')
ylabel('Hit Rate')
plot(probfalpos,probtrupos,'o--')


%% Question 10
disp('Question 10')

AUC = trapz(fliplr(probfalpos),fliplr(probtrupos))


%% Question 11
disp('Question 11')

Dprime = (mean(pos0) - mean(pos1))/sqrt(.5*(var(pos0)+var(pos1)))
