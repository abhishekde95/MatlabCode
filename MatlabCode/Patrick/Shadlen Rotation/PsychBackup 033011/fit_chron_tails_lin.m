% This script fits linear chronometric functions to each tail of the RT data
% function fit_chron_tails(df)

function [errs1 errs2 th_fit_chron1 th_fit_chron2]=fit_chron_tails_lin(df)

%% Set Variables

opts = optimset('MaxFunEvals',10.^6,'MaxIter',10.^6);
Bguess=min(df.meanrt);
kguess=(max(df.meanrt)-Bguess)/(max(df.coh));
tg2=[Bguess kguess kguess];
xaxneg=-.4:.01:0;
xaxpos=0:.01:.4;


% holdABkL=@(x)(x(1)./(x(3).*cohL)).*tanh(x(1)*x(3).*cohL)+x(2);
% holdABkR=@(x)(x(1)./(x(3).*cohR)).*tanh(x(1)*x(3).*df.coh(6:9))+x(2);
% ML=1./seL;
% MR=1./seR;
% NR=(holdABkR(vars)'-meanrtR).^2./(2*seR.^2);
% logLtR=log(MR)-NR;
% 
% fminsearch(holdABkL,vars,opts,df)


%% Find best fits for variables

%Hold B and k constant between tails
[th_fit_chron1, errs1]=fminsearch('chronometric_lin',poissrnd(tg2([1 2])),opts,df);
B1=th_fit_chron1(1);
k1=th_fit_chron1(2);
p_pred_chron1L=B1+k1.*xaxneg;
p_pred_chron1R=B1+(-k1)*xaxpos;

%Plot results
figure; hold on; grid on;
axis([-.4 .4 min(.9*df.meanrt) max(1.1*df.meanrt)])
title('Chronometric Fit (B and k Held Constant)')
xlabel('% Coherence / Direction')
ylabel('Reaction Time')
plot(df.coh,df.meanrt,'bo')
errorbar(df.coh,df.meanrt,df.sert,'bo')
plot(xaxneg,p_pred_chron1L,'g-')
plot(xaxpos,p_pred_chron1R,'g-')


%% Hold B constant, allow different k's

[th_fit_chron2, errs2]=fminsearch('chronometric_kchange_lin',poissrnd(tg2(1:3)),opts,df);
B2=th_fit_chron2(1);
k2L=th_fit_chron2(2);
k2R=th_fit_chron2(3);
p_pred_chron2L=B2+k2L.*xaxneg;
p_pred_chron2R=B2+k2R.*xaxpos;

%Plot results
figure; hold on; grid on;
axis([-.4 .4 min(.9*df.meanrt) max(1.1*df.meanrt)])
title('Chronometric Fit (B Constant, Different k)')
xlabel('% Coherence / Direction')
ylabel('Reaction Time')
plot(df.coh,df.meanrt,'bo')
errorbar(df.coh,df.meanrt,df.sert,'bo')
plot(xaxneg,p_pred_chron2L,'g-')
plot(xaxpos,p_pred_chron2R,'g-')




