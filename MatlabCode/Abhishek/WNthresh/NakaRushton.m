% Author- Abhishek De, 8/18
% Trying to understand Naka-Rushton function

close all; clearvars; 
x = linspace(0,10,51);
Rmax= 5;
bl = 0;
N = 2;
figure(1);set(gcf,'Name','Understanding effects of parameters');
color = ['r','g','b','m','k'];
Chalf = logspace(-1,1,4);
N1 = Chalf;
for ii = 1:numel(Chalf)
    R1 = (Rmax*(x.^N./(x.^N + Chalf(ii).^N))) + bl;
    R2 = (Rmax*(x.^N1(ii)./(x.^N1(ii) + Chalf(1).^N1(ii)))) + bl;
    subplot(121); plot(x,R1,color(ii),'Linewidth',2); hold on;
    subplot(122); plot(x,R2,color(ii),'Linewidth',2); hold on;
end
subplot(121);xlabel('x'); ylabel('Resp'); title('Naka-Rushton: C50'); axis square; hold off;
subplot(122);xlabel('x'); ylabel('Resp'); title('Naka-Rushton: N'); axis square; hold off;