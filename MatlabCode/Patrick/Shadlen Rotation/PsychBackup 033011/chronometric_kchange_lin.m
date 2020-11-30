function [p_pred, prodLt] = chronometric_kchange_lin(tg2,df)

%% Set variables
cohL=df.coh(1:4);
meanrtL=df.meanrt(1:4);
seL=df.sert(1:4);
cohR=df.coh(6:9);
meanrtR=df.meanrt(6:9);
seR=df.sert(6:9);
ML=1./seL;
MR=1./seR;
B=tg2(1);
kL=tg2(2);
kR=tg2(3);


%% Loop through predicted mean (pm) and likelihood (Lt) forumulas

% Adjust free parameters for left choice
pmL=B+kL*cohL;

% Calculate error for left parameters
NL=(pmL'-meanrtL).^2./(2*seL.^2);
logLtL=log(ML)-NL;

% Adjust free parameters for right choice
pmR=B+kR.*cohR;

% Calculate error for right parameters
NR=(pmR'-meanrtR).^2./(2*seR.^2);
logLtR=log(MR)-NR;

% Sum errors
p_pred=-sum(logLtL)-sum(logLtR);
prodLt=-sum(logLtL)-sum(logLtR);

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


