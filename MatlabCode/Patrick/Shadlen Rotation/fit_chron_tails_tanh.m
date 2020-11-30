% This script fits chronometric functions to each tail of the RT data
% function fit_chron_tails(df)

function [errs1 errs2 errs3 errs4 th_fit_chron1 th_fit_chron2 th_fit_chron3 th_fit_chron4]=fit_chron_tails_tanh(df)

%% Set Variables

resamps=5;

opts = optimset('MaxFunEvals',10.^6,'MaxIter',10.^6);
Aguess=sqrt(max(df.meanrt));
Bguess=min(df.meanrt);
kguess=(Aguess^2-Bguess)/(max(df.coh));
tg2=[Aguess Aguess Bguess kguess kguess];
xaxneg=-.4:.01:0;
xaxpos=0:.01:.4;
xax=-.4:.01:.4;


% holdABkL=@(x)(x(1)./(x(3).*cohL)).*tanh(x(1)*x(3).*cohL)+x(2);
% holdABkR=@(x)(x(1)./(x(3).*cohR)).*tanh(x(1)*x(3).*df.coh(6:9))+x(2);
% ML=1./seL;
% MR=1./seR;
% NR=(holdABkR(vars)'-meanrtR).^2./(2*seR.^2);
% logLtR=log(MR)-NR;
% 
% fminsearch(holdABkL,vars,opts,df)


%% Find best fits for variables

%Hold A, B, and k constant between tails
for i=1:resamps
    [th_fit_chron1, prodLt1]=fminsearch('chronometric',poissrnd(tg2([1 3 5])),opts,df);
    errs1(i)=prodLt1;
    vars1(i,:)=th_fit_chron1;
end
th_fit_chron1=vars1((errs1==min(errs1)),:);
errs1=min(errs1);
A1=th_fit_chron1(1);
B1=th_fit_chron1(2);
k1=th_fit_chron1(3);
p_pred_chron1=(A1./(k1.*xax)).*tanh(A1*k1.*xax)+B1;
p_pred_chron1(isnan(p_pred_chron1))=A1^2+B1;

%Plot results
figure; hold on; grid on;
axis([-.4 .4 min(.9*df.meanrt) max(1.1*df.meanrt)])
title('Chronometric Fit (A, B, and k Held Constant)')
xlabel('% Coherence / Direction')
ylabel('Reaction Time')
plot(df.coh,df.meanrt,'bo')
errorbar(df.coh,df.meanrt,df.sert,'bo')
plot(xax,p_pred_chron1,'g-')


%% Hold A and B constant, allow different k's

for i=1:resamps
    [th_fit_chron2, prodLt2]=fminsearch('chronometric_kchange',poissrnd(tg2([1 3:5])),opts,df);
    errs2(i)=prodLt2;
    vars2(i,:)=th_fit_chron2;
end
th_fit_chron2=vars2((errs2==min(errs2)),:);
errs2=min(errs2);
A2=th_fit_chron2(1);
B2=th_fit_chron2(2);
k2l=th_fit_chron2(3);
k2r=th_fit_chron2(4);
p_pred_chron2l=(A2./(k2l.*xaxneg)).*tanh(A2*k2l.*xaxneg)+B2;
p_pred_chron2l(isnan(p_pred_chron2l))=A2^2+B2;
p_pred_chron2r=(A2./(k2r.*xaxpos)).*tanh(A2*k2r.*xaxpos)+B2;
p_pred_chron2r(isnan(p_pred_chron2r))=A2^2+B2;

%Plot results
figure; hold on; grid on;
axis([-.4 .4 min(.9*df.meanrt) max(1.1*df.meanrt)])
title('Chronometric Fit (A and B constant, different k)')
xlabel('% Coherence / Direction')
ylabel('Reaction Time')
plot(df.coh,df.meanrt,'bo')
errorbar(df.coh,df.meanrt,df.sert,'bo')
plot(xaxneg,p_pred_chron2l,'g-')
plot(xaxpos,p_pred_chron2r,'g-')


%% Hold B and k constant, allow different A's
for i=1:resamps
    [th_fit_chron3, prodLt3]=fminsearch('chronometric_Achange',poissrnd(tg2(1:4)),opts,df);
    errs3(i)=prodLt3;
    vars3(i,:)=th_fit_chron3;
end
th_fit_chron3=vars3((errs3==min(errs3)),:);
errs3=min(errs3);
A3L=th_fit_chron3(1);
A3R=th_fit_chron3(2);
B3=th_fit_chron3(3);
k3=th_fit_chron3(4);
p_pred_chron3L=(A3L./(k3.*xaxneg)).*tanh(A3L*k3.*xaxneg)+B3;
p_pred_chron3L(isnan(p_pred_chron3L))=A3L^2+B3;
p_pred_chron3R=(A3R./(k3.*xaxpos)).*tanh(A3R*k3.*xaxpos)+B3;
p_pred_chron3R(isnan(p_pred_chron3R))=A3R^2+B3;

%Plot results
figure; hold on; grid on;
axis([-.4 .4 min(.9*df.meanrt) max(1.1*df.meanrt)])
title('Chronometric Fit (B and k constant, different A)')
xlabel('% Coherence / Direction')
ylabel('Reaction Time')
plot(df.coh,df.meanrt,'bo')
errorbar(df.coh,df.meanrt,df.sert,'bo')
plot(xaxneg,p_pred_chron3L,'g-')
plot(xaxpos,p_pred_chron3R,'g-')


%% Hold B constant, allow different k's and A's
for i=1:resamps
    [th_fit_chron4, prodLt4]=fminsearch('chronometric_Akchange',poissrnd(tg2(1:5)),opts,df);
    errs4(i)=prodLt4;
    vars4(i,:)=th_fit_chron4;
end
th_fit_chron4=vars4((errs4==min(errs4)),:);
errs4=min(errs4);
A4L=th_fit_chron4(1);
A4R=th_fit_chron4(2);
B4=th_fit_chron4(3);
k4L=th_fit_chron4(4);
k4R=th_fit_chron4(5);

p_pred_chron4L=(A4L./(k4L.*xaxneg)).*tanh(A4L*k4L.*xaxneg)+B4;
p_pred_chron4L(isnan(p_pred_chron4L))=A4L^2+B4;
p_pred_chron4R=(A4R./(k4R.*xaxpos)).*tanh(A4R*k4R.*xaxpos)+B4;
p_pred_chron4R(isnan(p_pred_chron4R))=A4R^2+B4;

%Plot results
figure; hold on; grid on;
axis([-.4 .4 min(.9*df.meanrt) max(1.1*df.meanrt)])
title('Chronometric Fit (B and k constant, different A)')
xlabel('% Coherence / Direction')
ylabel('Reaction Time')
plot(df.coh,df.meanrt,'bo')
errorbar(df.coh,df.meanrt,df.sert,'bo')
plot(xaxneg,p_pred_chron4L,'g-')
plot(xaxpos,p_pred_chron4R,'g-')








