% adding 2 naka-rushtons


x = linspace(-1,1,101);
[xx,yy] = meshgrid(x,x);

A = 50;
bl = 0;
sig1 = .5;
sig2 = Inf;
orthsig = Inf;
exp = 3;
ang1 = 0;
ang2 = pi/2;
params1 = [A 1/sig1 1/sig2 1/orthsig exp bl ang1];
params2 = [A 1/sig1 1/sig2 1/orthsig exp bl ang2];

nkr1 = ComputeNakaRushtonJPW(params1,[xx(:) yy(:)],'conicsection_xy');
nkr1 = reshape(nkr1,size(xx));

nkr2 = ComputeNakaRushtonJPW(params2,[xx(:) yy(:)],'conicsection_xy');
nkr2 = reshape(nkr2,size(xx));

% plot
figure(1); clf; hold on;
surfc(xx,yy,nkr1)
figure(2); clf; hold on;
surfc(xx,yy,nkr2)

figure(3); clf; hold on
surfc(xx,yy,nkr1+nkr2)
