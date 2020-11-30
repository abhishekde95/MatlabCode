% Plotting laser intensity as function of dial reading: obtained the data from Yasmine and Robi
% Author - Abhishek De, 2/19

close all; clearvars;
dial = 0:0.5:10;
x = 0:0.05:10;
laserpower = [0.8;12.8;33.0;55.3;77.2;98.8;119;140.5;161.6;184.0;206;229;252;275;296;318;339;360;380;397;413];
y = spline(dial,laserpower,x);
figure(1); set(gcf,'Name','Laser dial-power relationship: sinusoidal blue laser');
plot(dial,laserpower,'-o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on; plot(x,y,'k','Linewidth',2);
xlabel('dial'); ylabel('power(mW)'); legend({'Data','Spline'},'FontSize',8);

laserdetails.dial = dial;
laserdetails.laserpower = laserpower;
laserdetails.dialspline = x;
laserdetails.laserpowerspline = y;
save laserdetails laserdetails