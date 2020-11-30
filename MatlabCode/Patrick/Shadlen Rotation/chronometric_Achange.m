function [p_pred, prodLt] = chronometric_Achange(tg2,df)

%% Set variables
cohl=df.coh(1:4);
meanrtl=df.meanrt(1:4);
sel=df.sert(1:4);
cohr=df.coh(6:9);
meanrtr=df.meanrt(6:9);
ser=df.sert(6:9);
Ml=1./sel;
Mr=1./ser;
AL=tg2(1);
AR=tg2(2);
B=tg2(3);
k=tg2(4);


%% Loop through predicted mean (pm) and likelihood (Lt) forumulas

% Adjust free parameters for left choice
pml=(AL./(k.*cohl)).*tanh(AL*k.*cohl)+B;
pml(isnan(pml))=AL^2+B;

% Calculate error for left parameters
Nl=(pml'-meanrtl).^2./(2*sel.^2);
logLtl=log(Ml)-Nl;

% Adjust free parameters for right choice
pmr=(AR./(k.*cohr)).*tanh(AR*k.*cohr)+B;
pmr(isnan(pmr))=AR^2+B;

% Calculate error for right parameters
Nr=(pmr'-meanrtr).^2./(2*ser.^2);
logLtr=log(Mr)-Nr;

% Sum errors
p_pred=-sum(logLtl)-sum(logLtr);
prodLt=-sum(logLtl)-sum(logLtr)

%% Plot results in realtime to observe progression of fitting proceedure
figure(5); clf; hold on; grid on;
title('Working Chronometric Fit')
xlabel('% Coherence / Direction')
ylabel('Response Time')
axis([min(df.coh) max(df.coh) (.9*min(df.meanrt)) (1.1*max(df.meanrt))])
plot(df.coh,df.meanrt,'bo')
errorbar(df.coh,df.meanrt,df.sert,'bo')
plot(cohl,pml,'-g')
plot(cohr,pmr,'-g')


