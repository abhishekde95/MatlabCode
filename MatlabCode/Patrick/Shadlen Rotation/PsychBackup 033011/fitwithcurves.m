function [th_fit, err] = fitwithcurves(df,tg,q)

%% Specify variables

xax = -.4:.001:.4;
opts = optimset('MaxFunEvals',10.^6,'MaxIter',10.^6);
tg(1) = (max(df.propright)-min(df.propright))/length(df.coh);
tg(3) = min(df.meanrt);
tg(2) = sqrt(max(df.meanrt)-tg(3));

%% Evaluate the original data

[th_fit, err] = fminsearch('calcparams',tg(1:4),opts,df);
k   = th_fit(1);
A   = th_fit(2);
Tnd = th_fit(3);
B   = th_fit(4);
Cb  = B/(2*A*k);

p_pred_psych = (1./(1+exp(-2*k*A*(xax-Cb))));
p_pred_chron(xax~=0) = A./(k.*(xax(xax~=0)-Cb)).*tanh(k.*A.*(xax(xax~=0)-Cb))+Tnd;
p_pred_chron(xax==0) = A.^2+Tnd;


%% Plot results
figure(q); clf;
subplot(2,1,1); hold on; grid on;
axis([-.4 .4 0 1])
xlabel('% Coherence/Direction')
ylabel('Proportion Right Choice')
title('Psychometric Fit to Data')
plot(df.coh,df.propright,'bo');
errorbar(df.coh,df.propright,df.serc,'og')
plot(xax,p_pred_psych,'b--')

subplot(2,1,2); hold on; grid on;
axis([-.4 .4 .9*min(df.meanrt) 1.1*max(df.meanrt)])
xlabel('% Coherence/Direction')
ylabel('Reaction Time')
title('Chronometric Fit to Data')
axis([-.4 .4 (.9*min(df.meanrt)) (1.1*max(df.meanrt))])
plot(df.coh,df.meanrt,'go')
errorbar(df.coh,df.meanrt,df.sert,'go')
plot(xax,p_pred_chron,'b--')

