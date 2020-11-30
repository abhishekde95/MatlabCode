%% This program fits a psychometric function to model behavior from cones
% Created Mar_2011 for Shadlen Rotation
% Modified Jun_2011 for Horwitz/Rieke Collaboration

function [err,k,A] = psych_curve(tg,df)

k   = tg(1);
A   = tg(2);
cont = df.contrast;
sumt = df.sumtrials;
propr = df.propright.*sumt;

%% Predict values of Pc with adjusted variables

% Psychometric function
Pc = 1./(1+exp(-2*k*A*cont))';


%% Compare predicted values with observed values

% Cost function for psychometric predictions
err = -sum(propr.*log(Pc'+eps) + (sumt-propr).*log(1-Pc'+eps));


%% Plot psychometric results in realtime to observe progression of fitting proceedure
figure(5); clf; hold on; grid on;
title('Working Psychometric Fit')
xlabel('% conterence / Direction')
ylabel('Probability of Right Choice')
%axis([min(cont) max(cont) 0.5 1])
plot(cont,df.propright,'bo')
plot(cont,Pc,'-g')

