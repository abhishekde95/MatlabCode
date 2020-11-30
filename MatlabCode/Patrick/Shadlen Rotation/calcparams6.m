function [err] = calcparams6(tg,df)

%% Set Variables (currently an inelligant hack)
i=1;
if abs(df(i).coh(find(df(i).coh==0)+1))>abs(df(i).coh(find(df(i).coh==0)-1))
    refH1=find(df(i).coh==0)+1:length(df(i).coh);
    refL1=1:find(df(i).coh==0)-1;
    cohH1=df(i).coh(refH1)';
    cohL1=df(i).coh(refL1)';
    sumtH1 = df(i).sumtrials(refH1);
    sumtL1 = df(i).sumtrials(refL1);
    proprH1 = df(i).propright(refH1).*sumtH1;
    proprL1 = df(i).propright(refL1).*sumtL1;
    rtH1 = df(i).meanrt(refH1);
    rtL1 = df(i).meanrt(refL1);
    seH1 = df(i).sert(refH1);
    seL1 = df(i).sert(refL1);
else
    refH1=1:find(df(i).coh==0)-1;
    refL1=find(df(i).coh==0)+1:length(df(i).coh);
    cohH1=df(i).coh(refH1)';
    cohL1=df(i).coh(refL1)';
    sumtH1 = df(i).sumtrials(refH1);
    sumtL1 = df(i).sumtrials(refL1);
    proprH1 = df(i).propright(refH1).*sumtH1;
    proprL1 = df(i).propright(refL1).*sumtL1;
    rtH1 = df(i).meanrt(refH1);
    rtL1 = df(i).meanrt(refL1);
    seH1 = df(i).sert(refH1);
    seL1 = df(i).sert(refL1);
end

i=2;
if abs(df(i).coh(find(df(i).coh==0)+1))>abs(df(i).coh(find(df(i).coh==0)-1))
    refH2=find(df(i).coh==0)+1:length(df(i).coh);
    refL2=1:find(df(i).coh==0)-1;
    cohH2=df(i).coh(refH2)';
    cohL2=df(i).coh(refL2)';
    sumtH2 = df(i).sumtrials(refH2);
    sumtL2 = df(i).sumtrials(refL2);
    proprH2 = df(i).propright(refH2).*sumtH2;
    proprL2 = df(i).propright(refL2).*sumtL2;
    rtH2 = df(i).meanrt(refH2);
    rtL2 = df(i).meanrt(refL2);
    seH2 = df(i).sert(refH2);
    seL2 = df(i).sert(refL2);
else
    refH2=1:find(df(i).coh==0)-1;
    refL2=find(df(i).coh==0)+1:length(df(i).coh);
    cohH2=df(i).coh(refH2)';
    cohL2=df(i).coh(refL2)';
    sumtH2 = df(i).sumtrials(refH2);
    sumtL2 = df(i).sumtrials(refL2);
    proprH2 = df(i).propright(refH2).*sumtH2;
    proprL2 = df(i).propright(refL2).*sumtL2;
    rtH2 = df(i).meanrt(refH2);
    rtL2 = df(i).meanrt(refL2);
    seH2 = df(i).sert(refH2);
    seL2 = df(i).sert(refL2);
end

k   = tg(1);
AH  = tg(2);
AL  = tg(3);
Tnd = tg(4);

%% Predict values of Pc and tT with adjusted variables

% Chronometric function
% High coherences
tTH1=(AH./(k.*cohH1)).*tanh(k.*AH.*cohH1)+Tnd;
tTH2=(AH./(k.*cohH2)).*tanh(k.*AH.*cohH2)+Tnd;

% Low coherences
tTL1=(AL./(k.*cohL1)).*tanh(k.*AL.*cohL1)+Tnd;
tTL2=(AL./(k.*cohL2)).*tanh(k.*AL.*cohL2)+Tnd;

% % Psychometric function
% High coherences
PcH1=(1./(1+exp(-2*k*AH*cohH1)))';
PcH2=(1./(1+exp(-2*k*AH*cohH2)))';

% Low coherences
PcL1=(1./(1+exp(-2*k*AL*cohL1)))';
PcL2=(1./(1+exp(-2*k*AL*cohL2)))';


%% Compare predicted values with observed values

% Cost function for chronometric predictions
% High coherences
err = sum((rtH1-tTH1).^2./(2*seH1.^2));
err = err + sum((rtH2-tTH2).^2./(2*seH2.^2));

% Low coherences
err = err + sum((rtL1-tTL1).^2./(2*seL1.^2));
err = err + sum((rtL2-tTL2).^2./(2*seL2.^2));


% Cost function for psychometric predictions
% High coherences
err = err - sum(proprH1.*log(PcH1+eps)' + (sumtH1-proprH1).*log(1-PcH1+eps)');
err = err - sum(proprH2.*log(PcH2+eps)' + (sumtH2-proprH2).*log(1-PcH2+eps)');

% Low coherences
err = err - sum(proprL1.*log(PcL1+eps)' + (sumtL1-proprL1).*log(1-PcL1+eps)');
err = err - sum(proprL2.*log(PcL2+eps)' + (sumtL2-proprL2).*log(1-PcL2+eps)');
 
 
%% Plot psychometric results in realtime to observe progression of fitting proceedure
% figure(5); clf; hold on; grid on;
% title('Working Psychometric Fit')
% xlabel('% Coherence / Direction')
% ylabel('Probability of Right Choice')
% axis([-.4 .4 0 1])
% plot(cohH1,proprH1./sumtH1,'bo')
% plot(cohH2,proprH2./sumtH2,'bo')
% plot(cohL1,proprL1./sumtL1,'ro')
% plot(cohL2,proprL2./sumtL2,'ro')
% errorbar(cohH1,proprH1./sumtH1,seH1./1000,'bo')
% errorbar(cohH2,proprH2./sumtH2,seH2./1000,'bo')
% errorbar(cohL1,proprL1./sumtL1,seL1./1000,'ro')
% errorbar(cohL2,proprL2./sumtL2,seL2./1000,'ro')
% plot(cohH1,PcH1,'-c')
% plot(cohH2,PcH2,'-c')
% plot(cohL1,PcL1,'-g')
% plot(cohL2,PcL2,'-g')

% figure(6); clf; hold on; grid on;
% title('Working Chronometric Fit')
% xlabel('% Coherence / Direction')
% ylabel('Reaction Time')
% axis([-.4 .4 min(min(df(1).meanrt),min(df(2).meanrt)) max(max(df(1).meanrt),max(df(2).meanrt))])
% plot(cohH1,rtH1,'bo')
% plot(cohH2,rtH2,'bo')
% plot(cohL1,rtL1,'ro')
% plot(cohL2,rtL2,'ro')
% errorbar(cohH1,rtH1,seH1./1000,'bo')
% errorbar(cohH2,rtH2,seH2./1000,'bo')
% errorbar(cohL1,rtL1,seL1./1000,'ro')
% errorbar(cohL2,rtL2,seL2./1000,'ro')
% plot(cohH1,tTH1,'-c')
% plot(cohH2,tTH2,'-c')
% plot(cohL1,tTL1,'-m')
% plot(cohL2,tTL2,'-m')
