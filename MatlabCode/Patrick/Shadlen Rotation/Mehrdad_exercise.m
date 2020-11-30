%% 1. Generate function, 2. get probability at each point, 3. generate data

%Variables
a=.1;
mu=.5;
coh=(-1:.1:1)';
trials_per_coh=50;
alpha_est=1;
mu_est=1;
theta_guess=[alpha_est,mu_est];

%Function
myfunc = @(coh)1./(1+exp((-coh+mu)/a));
myfunc(coh);

%Generate fake dataset
q=1;
subj_resp=nan(length(coh)*trials_per_coh,2);
for i=1:size(coh)
    for j=1:trials_per_coh
        subj_resp(q,1)=coh(i);
        if rand>(1-myfunc(coh(i)))
            subj_resp(q,2)=1;
        else
            subj_resp(q,2)=0;
        end
        q=q+1;
    end
end


%% Impliment fminsearch


n=zeros(length(coh),1);
m=zeros(size(n));

for i=1:length(coh)
    for j=1:length(subj_resp)
        if subj_resp(j,1)==coh(i)
            n(i)=n(i)+subj_resp(j,2);
        end
    end
    m(i)=trials_per_coh-n(i);
end

opts = optimset('MaxFunEvals',10.^6,'MaxIter',10.^6);
[theta_fit]=fminsearch('myfunc2',theta_guess,opts,coh,n,m)



p_pred=zeros(length(coh),1);
for i=1:length(coh)
    p_pred(i)=1./(1+exp((-coh(i)+theta_fit(2))/theta_fit(1)));
end
p_mean=n./(n+m);
p_se=sqrt(p_mean.*(1-p_mean)./n);



% Plot the data, sucka!
figure
errorbar(coh,p_mean,p_se,'ob');hold on
plot(coh,p_pred,'k--')
title('Proportion Correct');xlabel('%Coherence/Direction');ylabel('Proportion Correct');
grid on
