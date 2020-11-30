function [th_fit_san,err_san] = sanitycheck3(df,tg,n)


%% Set variables to recover

k    = 0.55;
A    = 16;
Tnd  = 400;
BL   = .6;
BH   = .1;

CohL = [-.32 -.16 -.08 -.04];
CohH = [.05 .1 .2 .4];

xaxL = -.4:.001:0;
xaxH = 0:.001:.4;
rtL = nan(n,4);
rtH = nan(n,4);
opts = optimset('MaxFunEvals',10.^6,'MaxIter',10.^6);


%% Functions to generate data

PL = 1./(1 + exp(-2*k*A*CohL+BL));
PH = 1./(1 + exp(-2*k*A*CohH+BH));

RtL = A./(k*CohL-(BL/(2*A))).*tanh(k*A*CohL-(BL/2)) + Tnd;
RtH = A./(k*CohH-(BH/(2*A))).*tanh(k*A*CohH-(BH/2)) + Tnd;


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

tg(1) = (max(meanrtL)-min(meanrtL))/length(CohL);
tg(2) = sqrt(max(max(meanrtH),max(meanrtL)));
tg(3) = min(min(meanrtH),min(meanrtL));
tg(4) = 1;
tg(5) = 1;

[th_fit_san,err_san] = fminsearch('calcsanitycheck3',tg(1:5),opts,df_sanity);

k   = th_fit_san(1);
A   = th_fit_san(2);
Tnd = th_fit_san(3);
BL  = th_fit_san(4);
BH  = th_fit_san(5);

p_pred_psychL = 1./(1+exp(-2*k*A*xaxL+BL));
p_pred_psychH = 1./(1+exp(-2*k*A*xaxH+BH));
p_pred_chronL = A./(k.*xaxL-(BL/(2*A))).*tanh(k.*A.*xaxL-(BL/2))+Tnd;
p_pred_chronH = A./(k.*xaxH-(BH/(2*A))).*tanh(k.*A.*xaxH-(BH/2))+Tnd;

% Plot the psychometric function against the original data
figure(32); clf; 
subplot(2,1,1); hold on; grid on;
axis([-.4 .4 0 1]);
title('Psychometric Fit to Sanity Check Data (Two Bs)');
xlabel('%Coherence/Direction');
ylabel('Proportion Right Choice');
plot(CohL,ObsL./n,'bo');
plot(CohH,ObsH./n,'bo');
plot(xaxL,p_pred_psychL,'b--')
plot(xaxH,p_pred_psychH,'b--')

% Plot the chronometric function against the original data
subplot(2,1,2); hold on; grid on;
title('Chronometric Fit to Sanity Check Data (Two Bs)')
xlabel('% Coherence / Direction')
ylabel('Response Time')
axis([-.4 .4 .9*min(min(meanrtL),min(meanrtH)) 1.1*max(max(meanrtH),max(meanrtL))])
plot(CohL,meanrtL,'go')
plot(CohH,meanrtH,'go')
errorbar(CohL,meanrtL,sertL,'go')
errorbar(CohH,meanrtH,sertH,'go')
plot(xaxL,p_pred_chronL,'b--')
plot(xaxH,p_pred_chronH,'b--')

