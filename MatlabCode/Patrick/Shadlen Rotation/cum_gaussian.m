function p_pred=cum_gaussian(tg,coh,n,m)

p=zeros(length(coh),1);

for i=1:length(coh)    
    p(i)=-log(1./(1+exp((-coh(i)+tg(2))/tg(1)))^n(i)...
        *(1-(1./(1+exp((-coh(i)+tg(2))/(tg(1))))))^m(i));
end

[p_pred]=sum(p);