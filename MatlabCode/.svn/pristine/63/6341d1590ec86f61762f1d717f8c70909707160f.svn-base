%fitPrmtsTim gives different offsets to right(T1) & left(T2) choices. 

clear all

load TimData

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part calculates RT1, RT2, and # of T1 choices from Monkey1 data at
% each coherence.  And it sets a criteria to exclude certain coherences at
% which monkey made few choices (less than 10).


Complete_ind = logical(data(:,1)==1 | data(:,1)==0);
Break_ind = logical(data(:,1)==-1);
complete_trials = data(Complete_ind,:);
T1_ind = logical(complete_trials(:,2)==1); 
T2_ind = logical(complete_trials(:,2)==2); 
complete_T1 = complete_trials(T1_ind,:);
complete_T2 = complete_trials(T2_ind,:);

num_cohs = length(cohs);

for i=1:num_cohs
    % get indices for current coh
    coh1_ind = logical(complete_T1(:,3)==cohs(i));
    coh2_ind = logical(complete_T2(:,3)==cohs(i));
    % calculate values (note: only use corrects for RTs)
        RT1_mean(i) = mean(complete_T1(coh1_ind, 6));
        RT1_se(i)   = std(complete_T1(coh1_ind,6))./sqrt(sum(coh1_ind));
        RT2_mean(i) = mean(complete_T2(coh2_ind, 6));
        RT2_se(i)   = std(complete_T2(coh2_ind,6))./sqrt(sum(coh2_ind));
        RT_mean(i)  = mean([complete_T1(coh1_ind, 6);complete_T2(coh2_ind, 6)]);
        RT_se(i)    = std([complete_T1(coh1_ind, 6);complete_T2(coh2_ind, 6)])./sqrt(sum(coh1_ind)+sum(coh2_ind));
        nT1(i)      = sum(coh1_ind);
        n(i)        = sum(coh1_ind) + sum(coh2_ind);
        
        
        % Set a criteria for data.
        % If # of choices are less than 10, the data is excluded
        % because the variance will not be trustable.
        if nT1(i) <= 10 
            RT1_mean(i)=NaN; RT1_se(i)=NaN;
        elseif (n(i)-nT1(i)) <= 10
            RT2_mean(i)=NaN; RT2_se(i)=NaN;
        end
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%load Tim_data2
%----------------

data1=[cohs' RT1_mean' RT1_se' RT2_mean' RT2_se' nT1' n'];
%data2=[cohs' RT_mean' RT_se' nT1' n'];

p_mean = data1(:,6)./data1(:,7) ;
p_se   = sqrt(p_mean.*(1-p_mean)./data1(:,6));


thetaGuess = [.000457 25 300 300];
opts = optimset('fminsearch');
opts = optimset('MaxFunEvals',10.^5,'MaxIter',10.^5);
% [thetaFit,errFit,xflg,oput] = fminsearch('fitDiff', thetaGuess, opts, data);
[thetaFit,errFit,xflg,oput] = fminsearch('fitDiff5', thetaGuess, opts, data1);

thetaFit
errFit

xax=-600:20:600;
[t1_pred, t2_pred, p_pred] = calcDiff5(xax, thetaFit);

cohs_positive = logical(cohs > 0);
cohs_negative = logical(cohs < 0);
cohs_zero     = logical(cohs == 0);

figure
errorbar(cohs,RT1_mean,RT1_se,'ob');hold on
errorbar(cohs,RT2_mean,RT2_se,'sr');
plot(cohs(cohs_positive),RT1_mean(cohs_positive),'ob','MarkerFaceColor','b');
plot(cohs(cohs_negative),RT2_mean(cohs_negative),'sr','MarkerFaceColor','r');
plot(xax,t1_pred,'b--',xax,t2_pred,'r--');
title('Reaction Time');xlabel('Coherence');ylabel('Reaction Time');hold off

figure
errorbar(cohs,p_mean,p_se,'ob');hold on
plot(xax,p_pred,'k--');
title('Proportion Correct');xlabel('Coherence');ylabel('Proportion Correct');hold off