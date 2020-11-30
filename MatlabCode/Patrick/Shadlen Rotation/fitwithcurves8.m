% This code fits curves to psych and chron data for two A's, two k's, one Tnd,
% and no B for H and L coherences

function [th_fit8,err8] = fitwithcurves8(df,tg)

%% Set Variables

opts = optimset('MaxFunEvals',10.^6,'MaxIter',10.^6);
tg(3)=sqrt(max(df(1).meanrt));
tg(4)=sqrt(max(df(2).meanrt));
tg(5)=min(df(1).meanrt);

[th_fit8,err8]=fminsearch('calcparams8',tg(1:5),opts,df);

kH  = th_fit8(1);
kL  = th_fit8(2);
AH  = th_fit8(3);
AL  = th_fit8(4);
Tnd = th_fit8(5);

xax=-.4:.001:.4;
p_pred_psych2H=(1./(1+exp(-2*kH*AH*xax)));
p_pred_psych2L=(1./(1+exp(-2*kL*AL*xax)));
p_pred_chron2H=AH./(kH.*xax).*tanh(kH.*AH.*xax)+Tnd;
p_pred_chron2L=AL./(kL.*xax).*tanh(kL.*AL.*xax)+Tnd;

%% Plot results
figure(14); clf;
subplot(2,1,1); hold on; grid on;
axis([-.4 .4 0 1])
xlabel('% Coherence/Direction')
ylabel('Proportion Right Choice')
title('Psychometric Fit to High and Low Coherences (Different ks and As, No B)')
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
title('Chronometric Fit to High and Low Coherences (Different ks and As, No B)')
errorbar(df(1).coh(df(1).coh~=0),df(1).meanrt(df(1).coh~=0),df(1).sert(df(1).coh~=0),'ob')
errorbar(df(2).coh(df(2).coh~=0),df(2).meanrt(df(2).coh~=0),df(2).sert(df(2).coh~=0),'*r')
plot(df(1).coh(df(1).coh~=0),df(1).meanrt(df(1).coh~=0),'bo')
plot(df(2).coh(df(2).coh~=0),df(2).meanrt(df(2).coh~=0),'r*')
plot(xax(1:find(xax==0)),p_pred_chron2H(1:find(xax==0)),'b--')
plot(xax(find(xax==0):end),p_pred_chron2H(find(xax==0):end),'r-')
plot(xax(1:find(xax==0)),p_pred_chron2L(1:find(xax==0)),'r-')
plot(xax(find(xax==0):end),p_pred_chron2L(find(xax==0):end),'b--')

legend('Dataset 1 --','Dataset 2 -','Location','SouthEastOutside')