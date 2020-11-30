function [p_pred, prodLt] = chronometric(tg2,df)

%% Set variables
cohL=df.coh(1:4);
meanrtL=df.meanrt(1:4);
seL=df.sert(1:4);
cohR=df.coh(6:9);
meanrtR=df.meanrt(6:9);
seR=df.sert(6:9);
ML=1./seL;
MR=1./seR;
A=tg2(1);
B=tg2(2);
k=tg2(3);


%% Loop through predicted mean (pm) and likelihood (Lt) forumulas

% Adjust free parameters for left choice
pmL=(A./(k.*cohL)).*tanh(A*k.*cohL)+B;
pmL(isnan(pmL))=A^2+B;

% Calculate error for left parameters
NL=(pmL'-meanrtL).^2./(2*seL.^2);
logLtL=log(ML)-NL;

% Adjust free parameters for right choice
pmR=(A./(k.*cohR)).*tanh(A*k.*cohR)+B;
pmR(isnan(pmR))=A^2+B;

%Calculate error for right parameters
NR=(pmR'-meanrtR).^2./(2*seR.^2);
logLtR=log(MR)-NR;

% Sum errors
p_pred=-sum(logLtL)-sum(logLtR);
prodLt=-sum(logLtL)-sum(logLtR)

%% Plot results in realtime to observe progression of fitting proceedure
figure(5); clf; hold on; grid on;
title('Working Chronometric Fit')
xlabel('% Coherence / Direction')
ylabel('Response Time')
axis([min(df.coh) max(df.coh) (.9*min(df.meanrt)) (1.1*max(df.meanrt))])
plot(df.coh,df.meanrt,'bo')
errorbar(df.coh,df.meanrt,df.sert,'bo')
plot(cohL,pmL,'-g')
plot(cohR,pmR,'-g')


