function [th_fit_san,err_san] = sanitycheck_t2(df,tg,n)


%% Set variables to recover
%n=50;

kL   = 0.6;
kH   = 0.5;
A    = 16;
Tnd  = 400;
%B    = 0.2;
Cb   = 0.02;

CohL = [-.32 -.16 -.08 -.04];
CohH = [.05 .1 .2 .4];

xaxL = -.4:.001:0;
xaxH = 0:.001:.4;
rtL = nan(n,4);
rtH = nan(n,4);
opts = optimset('MaxFunEvals',10.^6,'MaxIter',10.^6);


%% Functions to generate data

% PL = 1./(1 + exp(-2*kL*A*CohL+B));
% PH = 1./(1 + exp(-2*kH*A*CohH+B));

PL = 1./(1 + exp(-2*kL*A*(CohL+Cb)));
PH = 1./(1 + exp(-2*kH*A*(CohH+Cb)));

% RtL = A./(kL*CohL-(B/(2*A))).*tanh(kL*A*CohL-(B/2)) + Tnd;
% RtH = A./(kH*CohH-(B/(2*A))).*tanh(kH*A*CohH-(B/2)) + Tnd;

RtL = A./(kL*(CohL-Cb)).*tanh(kL*A*(CohL-Cb)) + Tnd;
RtH = A./(kH*(CohH-Cb)).*tanh(kH*A*(CohH-Cb)) + Tnd;


ObsL = binornd(n,PL);
ObsH = binornd(n,PH);
for i=1:n
    rtL(i,:) = RtL+(df.stdrt(1:4)).*(randn(1,4));
    rtH(i,:) = RtH+(df.stdrt(6:9)).*(randn(1,4));
end
meanrtL = mean(rtL);
meanrtH = mean(rtH);
sertL = df.sert(1:4);
sertH = df.sert(6:9);
sumt  = [n n n n];


[df_sanity] = [CohL CohH; ObsL ObsH; meanrtL meanrtH; sertL sertH; sumt sumt];

%% Try to recover orignal variable values

tg(1) = max(CohL)-min(CohL)/length(CohL);
tg(2) = max(CohH)-min(CohH)/length(CohH);
tg(3) = sqrt(max(meanrtH));
tg(4) = min(meanrtH);

[th_fit_san,err_san] = fminsearch('calcsanitycheck_t2',tg(1:5),opts,df_sanity);

kL  = th_fit_san(1);
kH  = th_fit_san(2);
A   = th_fit_san(3);
Tnd = th_fit_san(4);
Cb  = th_fit_san(5);

% p_pred_psychL = 1./(1+exp(-2*kL*A*xaxL+B));
% p_pred_psychH = 1./(1+exp(-2*kH*A*xaxH+B));
% p_pred_chronL = A./(kL.*xaxL-(B/(2*A))).*tanh(kL.*A.*xaxL-(B/2))+Tnd;
% p_pred_chronH = A./(kH.*xaxH-(B/(2*A))).*tanh(kH.*A.*xaxH-(B/2))+Tnd;

p_pred_psychL = 1./(1+exp(-2*kL*A*(xaxL+Cb)));
p_pred_psychH = 1./(1+exp(-2*kH*A*(xaxH+Cb)));
p_pred_chronL = (A./(kL.*(xaxL-Cb))).*tanh(kL.*A.*(xaxL-Cb))+Tnd;
p_pred_chronH = (A./(kH.*(xaxH-Cb))).*tanh(kH.*A.*(xaxH-Cb))+Tnd;


%% Plot results
% Plot the psychometric function against the original data
figure(30); clf; 
subplot(2,1,1); hold on; grid on;
axis([-.4 .4 0 1]);
title('Psychometric Fit to Sanity Check Data (Two ks)');
xlabel('%Coherence/Direction');
ylabel('Proportion Right Choice');
plot(CohL,ObsL./n,'bo');
plot(CohH,ObsH./n,'bo');
plot(xaxL,p_pred_psychL,'b--')
plot(xaxH,p_pred_psychH,'b--')

% Plot the chronometric function against the original data
subplot(2,1,2); hold on; grid on;
title('Chronometric Fit to Sanity Check Data (Two ks)')
xlabel('% Coherence / Direction')
ylabel('Response Time')
axis([-.4 .4 .9*min(min(meanrtL),min(meanrtH)) 1.1*max(max(meanrtH),max(meanrtL))])
plot(CohL,meanrtL,'go')
plot(CohH,meanrtH,'go')
errorbar(CohL,meanrtL,sertL,'go')
errorbar(CohH,meanrtH,sertH,'go')
plot(xaxL,p_pred_chronL,'b--')
plot(xaxH,p_pred_chronH,'b--')



