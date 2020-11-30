% Compare STA in two different distributions

stim1 = randn(100,2);
resps = rand(100,1);

STA1 = mean(stim1 .* repmat(resps,1,2));

M = [15 3; .5 .5];

stim2 = stim1 * M';
stim2 = stim1.^2;
STA2 = mean(stim2 .* repmat(resps,1,2));
STA2_trans = STA1 * M';
STA2_trans = STA1.^2;

figure(1); clf; hold on;
plot(stim1(:,1),stim1(:,2),'ko')
plot(STA1(1),STA1(2),'r*')

figure(2); clf; hold on;
plot(stim2(:,1),stim2(:,2),'ko')
plot(STA2(1),STA2(2),'g*')
plot(STA2_trans(1),STA2_trans(2),'ro')



