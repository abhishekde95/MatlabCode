%function p_pred=psychometric(tg,coh,n,m)
function p_pred=psychometric(tg,df)

%Variables
coh=df.coh;
propright=df.propright;
sumtrials=df.sumtrials;
se=(sqrt(propright.*sumtrials.*(1-propright))./sumtrials);
se(isnan(se))=eps;
M=1./se;
M(isinf(M))=eps;
A=tg(1);
k=tg(2);
B=tg(3);

% p=zeros(length(coh),1);
% 
% for i=1:length(coh)    
%     p(i)=-log(1./(1+exp((-coh(i)+tg(2))/tg(1)))^n(i)...
%         *(1-(1./(1+exp((-coh(i)+tg(2))/(tg(1))))))^m(i));
% end

%% Loop through Predicted Probability of Right Choice (Ppr) and Likelihood (err) formulas

%Adjust free parameters A and k
Pc=1./(1+exp(-2*k*A*coh));


% Calculate error of parameters
Lp=sum(gammaln(sumtrials+1))-sum(gammaln(propright+1))-sum(gammaln(sumtrials-propright+1))...
    +propright.*log(Pc+eps)+(sumtrials-propright).*log(1-Pc+eps);

%N=(Pm'-propright).^2./(2*se.^2);
%N(isnan(N))=eps;

%Sum errors
% err=log(M)-N;
[p_pred]=sum(Lp)


%% Plot results in realtime to observe progression of fitting proceedure
figure(6); clf; hold on; grid on;
title('Working Psychometric Fit')
xlabel('% Coherence / Direction')
ylabel('Probability of Right Choice')
axis([min(coh) max(coh) 0 1])
plot(coh,propright,'bo')
errorbar(coh,propright,se,'bo')
plot(coh,Pc,'-g')
