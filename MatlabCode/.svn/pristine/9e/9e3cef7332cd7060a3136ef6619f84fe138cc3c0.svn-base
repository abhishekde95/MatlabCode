% Plotting the Vlambda data into the RGB space along with 1L+1M axis 
% Author - Abhishek De, 12/18
close all; clearvars;
load('T_vos1978_Y');
load mon_spd.mat
load fundamentals.mat
plot_counter = 1;
Vlambda = T_vos1978_Y';
wave = linspace(380,780,numel(Vlambda));
fundamentals = reshape(fundamentals,[numel(fundamentals)/3 3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals
M = inv(M');

figure(plot_counter); subplot(221); plot(wave,Vlambda,'k','Linewidth',2); hold on; plot(wave,fundamentals,'Linewidth',2); set(gca,'Xlim',[wave(1) wave(end)]); hold off;
subplot(222); plot(wave,mon_spd,'Linewidth',2); set(gca,'Xlim',[wave(1) wave(end)]);
LMSweighting = fundamentals'*Vlambda; 
LMSweighting = LMSweighting/norm(LMSweighting);

RGBVlambda = inv(M)*LMSweighting; RGBVlambda = RGBVlambda/norm(RGBVlambda);
RGBLM = inv(M)*[1; 1; 0]; RGBLM = RGBLM/norm(RGBLM);
anglediff = acos(RGBLM'*RGBVlambda)*180/pi;
subplot(223); plot3([RGBVlambda(1) -1*RGBVlambda(1)],[RGBVlambda(2) -1*RGBVlambda(2)],[RGBVlambda(3) -1*RGBVlambda(3)],'k','Linewidth',2); hold on; 
plot3([RGBLM(1) -1*RGBLM(1)],[RGBLM(2) -1*RGBLM(2)],[RGBLM(3) -1*RGBLM(3)],'r','Linewidth',2); grid on; axis square; 
xlabel('R'); ylabel('G'); zlabel('B'); hold off;




