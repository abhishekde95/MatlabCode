function [err] = calcsanitycheck2(tg,df_sanity)

k   = tg(1);
AL  = tg(2);
AH  = tg(3);
Tnd = tg(4);
B   = tg(5);

sumt   = df_sanity(5,1:4);
cohL   = df_sanity(1,1:4);
cohH   = df_sanity(1,5:8);
proprL = df_sanity(2,1:4);
proprH = df_sanity(2,5:8);
rtL    = df_sanity(3,1:4);
rtH    = df_sanity(3,5:8);
seL    = df_sanity(4,1:4);
seH    = df_sanity(4,5:8);


%% Predict values of Pc and tT with adjusted variables

% Chronometric function
tTL = AL./(k.*cohL-(B/(2*AL))).*tanh(k.*AL.*cohL-(B/2))+Tnd;
tTH = AH./(k.*cohH-(B/(2*AH))).*tanh(k.*AH.*cohH-(B/2))+Tnd;

% Psychometric function
PcL = 1./(1+exp(-2*k*AL*cohL+B))';
PcH = 1./(1+exp(-2*k*AH*cohH+B))';


%% Compare predicted values with observed values

% Cost function for chronometric predictions
LT = sum((rtL-tTL).^2./(2*seL.^2))+sum(log(seL))+length(cohL)/2*log(2*pi);
LT = LT + sum((rtH-tTH).^2./(2*seH.^2))+sum(log(seH))+length(cohH)/2*log(2*pi);

% Cost function for psychometric predictions
LP = LT - sum(proprL.*log(PcL+eps)' + (sumt-proprL).*log(1-PcL+eps)');
LP = LP - sum(proprH.*log(PcH+eps)' + (sumt-proprH).*log(1-PcH+eps)');

err = LP;
 
% %% Plot psychometric results in realtime to observe progression of fitting proceedure
% figure(5); clf; hold on; grid on;
% title('Working Psychometric Fit')
% xlabel('% Coherence / Direction')
% ylabel('Probability of Right Choice')
% axis([-.4 .4 0 1])
% plot(cohL,proprL./sumt,'bo')
% plot(cohH,proprH./sumt,'bo')
% plot(cohL,PcL,'-g')
% plot(cohH,PcH,'-g')
