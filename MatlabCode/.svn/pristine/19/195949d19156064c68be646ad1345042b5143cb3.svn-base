%% This program is an attempt to fit BOTH the psychometric and chronometric functions

function [err] = calcparams(tg, df)

k   = tg(1);
A   = tg(2);
Tnd = tg(3);
%B   = tg(4);
%CB  = B/(2*k*A);
CB = 0;
coh = df.coh;
sumt = df.sumtrials;
propr = df.propright.*sumt;
rt = df.meanrt;
se = df.sert;

%% Predict values of Pc and tT with adjusted variables

% Chronometric function
tT(coh~=0)=A./(k.*(coh(coh~=0)-CB)).*tanh(k.*A.*(coh(coh~=0)-CB))+Tnd;
tT(coh==0)=A.^2+Tnd;

% Psychometric function
Pc=(1./(1+exp(-2*k*A*(coh-CB))))';


%% Compare predicted values with observed values

% Cost function for chronometric predictions
err = sum((rt-tT).^2./(2*se.^2));

% Cost function for psychometric predictions
err = err - sum(propr.*log(Pc+eps) + (sumt-propr).*log(1-Pc+eps));

 
%% Plot psychometric results in realtime to observe progression of fitting proceedure
% figure(5); clf; hold on; grid on;
% title('Working Psychometric Fit')
% xlabel('% Coherence / Direction')
% ylabel('Probability of Right Choice')
% axis([min(coh) max(coh) 0 1])
% plot(coh,df.propright,'bo')
% errorbar(coh,df.propright,se./1000,'bo')
% plot(coh,Pc,'-g')

