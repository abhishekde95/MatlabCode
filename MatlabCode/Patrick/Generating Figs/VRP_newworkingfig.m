
meshx = linspace(0,1,100);
[meshxx,meshyy] = meshgrid(meshx,meshx);

stimx = linspace(0,1,5);
[stimxx,stimyy] = meshgrid(meshx,meshx);

A = 50;
sig1 = .5;
sig2 = Inf;
orthsig = Inf;
exp = 3;
bl = 0;
rot = pi/4;
params = [A 1/sig1 1/sig2 1/orthsig exp bl rot];

linmesh = reshape(ComputeNakaRushtonJPW(params,[meshxx(:) meshyy(:)],'conicsection'),size(meshxx));
nonlinmesh = reshape(ComputeNakaRushtonJPW(params,[meshxx(:).^2 meshyy(:).^2],'conicsection'),size(meshxx));

figure(1); clf; hold on; box on;
h = surf(meshxx,meshyy,linmesh); alpha(.5)
set(h,'edgecolor','none')
colormap spring
contour3(meshxx,meshyy,linmesh);
xlabel('Le')
ylabel('Me')
zlabel('Spikes')

figure(2); clf; hold on; box on;
h = surf(meshxx.^2,meshyy.^2,nonlinmesh); alpha(.5)
set(h,'edgecolor','none')
colormap winter
contour3(meshxx.^2,meshyy.^2,nonlinmesh);
xlabel('Le.^2')
ylabel('Me.^2')
zlabel('Spikes')

figure(3); clf; 
p = get(gca,'position');
contour(meshxx,meshyy,linmesh);
colormap(gca,'spring')
axes('units','normalized','pos',p,'color','none'); hold on;
contour(meshxx,meshyy,nonlinmesh);
colormap(gca,'winter')
xlabel('Le')
ylabel('Me')
zlabel('Spikes')

figure(4); clf;
p = get(gca,'pos');
contour(meshxx.^2,meshyy.^2,linmesh);
colormap(gca,'spring')
axes('units','normalized','pos',p,'color','none'); hold on;
contour3(meshxx.^2,meshyy.^2,nonlinmesh);
colormap(gca,'winter')
xlabel('Le^2')
ylabel('Me^2')
zlabel('Spikes')

