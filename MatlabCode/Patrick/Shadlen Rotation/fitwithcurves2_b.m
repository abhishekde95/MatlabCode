% This code fits curves to psych and chron data for two k's, one A, one B, and one Tnd
% for for H and L coherences

function [th_fit2,err2] = fitwithcurves2_b(df,tg)

%% Set Variables

opts = optimset('MaxFunEvals',10.^6,'MaxIter',10.^6);
tg(2)=sqrt(max(df(1).meanrt));
tg(3)=min(df(1).meanrt);

%% Fit Variables

[th_fit2,err2]=fminsearch('calcparams2_b',tg(1:4),opts,df);

kH  = th_fit2(1);
kL  = kH;
A   = th_fit2(2);
Tnd = th_fit2(3);
B   = th_fit2(4);
CbH = B/(2*kH*A);
CbL = B/(2*kL*A);

xax  = -.4:.001:.4;
p_pred_psych2H=(1./(1+exp(-2*kH*A*(xax-CbH))));
p_pred_psych2L=(1./(1+exp(-2*kL*A*(xax-CbL))));
p_pred_chron2H=A./(kH.*(xax-CbH)).*tanh(kH.*A.*(xax-CbH))+Tnd;
p_pred_chron2L=A./(kL.*(xax-CbL)).*tanh(kL.*A.*(xax-CbL))+Tnd;

%  p_pred_psych2Ht = 1./(1+exp(-2*kH*A*xax+B));
%  p_pred_psych2Lt = 1./(1+exp(-2*kL*A*xax+B));
%  p_pred_chron2Ht = (A./(kH.*xax-(B/(2*A)))).*tanh(kH.*A.*xax-(B/2))+Tnd;
%  p_pred_chron2Lt = (A./(kL.*xax-(B/(2*A)))).*tanh(kL.*A.*xax-(B/2))+Tnd;


%% Plot results
figure(8); clf;
subplot(2,1,1); hold on; grid on;
axis([-.4 .4 0 1])
xlabel('% Coherence/Direction')
ylabel('Proportion Right Choice')
title('Psychometric Fit to High and Low Coherences (Different ks)')
errorbar(df(1).coh(df(1).coh~=0),df(1).propright(df(1).coh~=0),df(1).serc(df(1).coh~=0),'ob')
errorbar(df(2).coh(df(2).coh~=0),df(2).propright(df(2).coh~=0),df(2).serc(df(2).coh~=0),'*r')
plot(df(1).coh(df(1).coh~=0),df(1).propright(df(1).coh~=0),'bo')
plot(df(2).coh(df(2).coh~=0),df(2).propright(df(2).coh~=0),'r*')
plot(xax(1:find(xax==0)),p_pred_psych2H(1:find(xax==0)),'b--')
plot(xax(find(xax==0):end),p_pred_psych2H(find(xax==0):end),'r-')
plot(xax(1:find(xax==0)),p_pred_psych2L(1:find(xax==0)),'r-')
plot(xax(find(xax==0):end),p_pred_psych2L(find(xax==0):end),'b--')

subplot(2,1,2); hold on; grid on;
axis([-.4 .4 min(min(df(1).meanrt),min(df(2).meanrt)) max(max(df(1).meanrt),max(df(2).meanrt))])
xlabel('% Coherence/Direction')
ylabel('Reaction Time')
title('Chronometric Fit to High and Low Coherences (Different ks)')
errorbar(df(1).coh(df(1).coh~=0),df(1).meanrt(df(1).coh~=0),df(1).sert(df(1).coh~=0),'ob')
errorbar(df(2).coh(df(2).coh~=0),df(2).meanrt(df(2).coh~=0),df(2).sert(df(2).coh~=0),'*r')
plot(df(1).coh(df(1).coh~=0),df(1).meanrt(df(1).coh~=0),'bo')
plot(df(2).coh(df(2).coh~=0),df(2).meanrt(df(2).coh~=0),'r*')
plot(xax(1:find(xax==0)),p_pred_chron2H(1:find(xax==0)),'b--')
plot(xax(find(xax==0):end),p_pred_chron2H(find(xax==0):end),'r-')
plot(xax(1:find(xax==0)),p_pred_chron2L(1:find(xax==0)),'r-')
plot(xax(find(xax==0):end),p_pred_chron2L(find(xax==0):end),'b--')

legend('Dataset 1 --','Dataset 2 -','Location','SouthEastOutside')

