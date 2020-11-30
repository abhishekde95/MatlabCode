
freqs = 0:1:100;
channel_noise = .001*ones(1,length(freqs));
transducer_noise = .001*ones(1,length(freqs));
signal = 1./(1+freqs);
transducer = 705*normpdf(freqs,20,5);
Sp = signal.*transducer;
lambda = -30;  % needs to be adjusted to meet power constraint
term1 = -channel_noise.*(2*transducer_noise+Sp);
term2 = sqrt((channel_noise.^2.*Sp.^2)-(4*transducer_noise.*Sp.*channel_noise.*lambda^-1));
denom = 2*transducer_noise.*(transducer_noise+Sp);
pn = (term1+term2)./denom;
figure;
subplot(3,2,1);
plot(freqs,signal)
title('signal');
subplot(3,2,2);
plot(freqs,transducer)
title('transducer');
subplot(3,2,3);
plot(freqs,Sp)
title('Sp');
subplot(3,2,4);
plot(freqs,transducer_noise,'k-');
plot(freqs,channel_noise,'r-')
title('noise');
subplot(3,2,5);
plot(freqs,pn)
title('neural filter');
subplot(3,2,6);
plot(freqs,sqrt(Sp))
title('sqrt(Sp)');

% lets see what K is equal to
K = sqrt(sum(Sp.*pn+transducer_noise.*pn+channel_noise))


% Varying lambda
for lambda = -1000:10:0001
term1 = -channel_noise.*(2*transducer_noise+Sp);
term2 = sqrt((channel_noise.^2.*Sp.^2)-(4*transducer_noise.*Sp.*channel_noise.*lambda^-1));
denom = 2*transducer_noise.*(transducer_noise+Sp);
pn = (term1+term2)./denom;
K = sqrt(sum(Sp.*pn+transducer_noise.*pn+channel_noise));
[lambda K]
end

