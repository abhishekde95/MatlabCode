function [kHdist kLdist AHdist ALdist] = resample_th_fit7(rs,tg,num_resamples)

th_fit7_rs = nan(num_resamples,6);
err7_rs = nan(num_resamples,1);
for z=1:num_resamples
    rsdf(1) = rs(1,z,1).df;
    rsdf(2) = rs(1,z,2).df;
    [th_fit7_rs(z,:),err7_rs(z,1)] = fitwithcurves7(rsdf,tg(1:6))
end

kHdist = th_fit7_rs(:,1);
kLdist = th_fit7_rs(:,2);
AHdist = th_fit7_rs(:,3);
ALdist = th_fit7_rs(:,4);

figure(35); clf; hold on; grid on;
title('Distribution of k values')
xlabel('k Value')
ylabel('Frequency')
hist(kHdist);
hist(kLdist);

figure(45); clf; hold on; grid on;
title('Distribution of A Values')
xlabel('A Value')
ylabel('Frequency')
hist(AHdist);
hist(ALdist);