% probable RFs with virus expression
% Author - Abhishek De, 10/18
clearvars; close all;
RFsite1 = [-80 71; -85 68; -100 70; -62 64]; % drive coordinates 0.5 RIGHT 6.5 POST
RFsite2 = [-90 9; -61 84; -72 73; -55 76; -32 74]; % drive coordinates 2 LEFT 6.5 POST
K1 = convhull(RFsite1(:,1),RFsite1(:,2));
K2 = convhull(RFsite2(:,1),RFsite2(:,2));
figure(1); plot(RFsite1(K1,1),RFsite1(K1,2),'g-','Linewidth',2); hold on;
plot(RFsite2(K2,1),RFsite2(K2,2),'r-','Linewidth',2);
set(gca,'TickDir','out','Xlim',[-100 100],'Ylim',[-100 100]); 
line([-100 100],[0 0]); line([0 0],[-100 100]); legend('0.5 RIGHT 6.5 POST','2 LEFT 6.5 POST','x axis','y axis'); 
xlabel('X'); ylabel('Y'); title('RFs at injection sites');axis square; grid on; hold off;
