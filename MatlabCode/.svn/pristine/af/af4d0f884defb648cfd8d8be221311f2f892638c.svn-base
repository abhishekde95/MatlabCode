% Understanding the effect of macular pigmentation
% Trying to see the differences between the 2 degree and 10 degree cone fundamentals
% Author - Abhishek De, 8/17
close all; clearvars;
load T_cones_smj.mat % Loading the 2 deg cone fundamentals
load T_cones_smj10.mat % Loading the 10 deg cone fundamentals
load('C:\Users\Abhishek\Desktop\MatlabCode\Slave\Monitor Calibration\Monitor data\ProPixx\ProPixx.mat');
Propixx = cals{end};
bkgnd_rgb = Propixx.bgColor;
mon_spd = Propixx.P_device;
mon_spd = spline([380:4:780], mon_spd', [380:5:780]);

% Color contrasts which I am interested in testing
max_rg_contrast = 0.045;
max_by_contrast = 0.45;
t=-pi:0.1:pi;
x=max_rg_contrast*cos(t);
y=max_by_contrast*sin(t);
color_contrasts = [x'/2 -x'/2 y'];
L = size(color_contrasts,1);
M10 = T_cones_smj10*mon_spd';
bkgnd_coneexcitations = M10*bkgnd_rgb;
color_excitations = (1+color_contrasts).*repmat(bkgnd_coneexcitations',[L 1]);
rgb_vals = inv(M10)*color_excitations';

% Now convert rgb values cone contrasts using 2 degree fundamentals
M2 = T_cones_smj*mon_spd';
cone_excitations2 = M2*rgb_vals;
bkgnd_coneexcitations2 = M2*bkgnd_rgb;
tmp = repmat(bkgnd_coneexcitations2,[1 L]);
color_contrasts2 = (cone_excitations2 - tmp)./tmp;
color_contrasts2 = color_contrasts2';

figure(1),plot([color_contrasts(:,1)-color_contrasts(:,2); color_contrasts(1,1)-color_contrasts(1,2)] , [color_contrasts(:,3);color_contrasts(1,3)],'g'); hold on;
plot([color_contrasts2(:,1)-color_contrasts2(:,2);color_contrasts2(1,1)-color_contrasts2(1,2)], [color_contrasts2(:,3);color_contrasts2(1,3)],'r');
plot(0, 0, 'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor','k'); hold off;

