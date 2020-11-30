% For playing with RWA

x = linspace(0,1,50);
nsp = ComputeNakaRushtonJPW([1 .5 3 0],x,'symmetric');
stanorm = sum(nsp./sum(nsp).*x);
staraw = mean(nsp.*x);
b = regress(nsp',x');

figure(44); clf; hold on;
plot(x,nsp,'ko')
plot(sta,mean(nsp),'go')
plot(x,x*b,'r--')
plot(x,x*stanorm,'b--')
plot(x,x*staraw,'g--')