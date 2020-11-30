% trogdor_injection.m
% 3D plot of injection sites
% CRF 11-23-14

%% cortex
clear all; close all;

figure; set(gcf, 'Color', 'w', 'Position', [400 100 700 900], 'PaperPositionMode', 'auto');
subplot(2,1,2);
txt = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n', ...
'1. rAAV1-hSyn-ChR2(H134R)-mCherry (undiluted)', ...
'2. rAAV1-hSyn-ChR2(H134R)-mCherry (diluted 1:2)', ...
'3. rAAV1-hSyn-ChR2(H134R)-mCherry (diluted 1:10)', ...
'4. Lenti-CamKII-ArchT-GFP', ...
'5. rAAV1-hSyn-ArchT-eYFP', ...
'6. rAAV5-hSyn-ChR2(H134R)-eYFP', ...
'7. rAAV9-hSyn-ChR2(H134R)-eYFP', ...
'8. rAAV1-hSyn-ChR2(H134R)-mCherry + rAAV1-hSyn-ArchT-eYFP (mixture)', ...
'9. rAAV1-hSyn-ChR2(H134R)-mCherry (thawed, inj only)', ...
'10 ("A"). rAAV8-DREADD(gpcr)-mCherry');
text(-0.1,0.5,txt,'FontSize',17); axis off;


subplot(2,1,1);
xlim([-7 7]); ylim([-7 7]); zlim([-1 9]); hold on; axis equal;
% xlabel('M-L (mm)','FontSize',15); ylabel('A-P (mm)','FontSize',15)
% zlabel('Depth (mm, wrt 16mm gt)','FontSize',15);
text(0,7,'A','FontSize',18,'HorizontalAlignment','center');
text(0,-7,'P','FontSize',18,'HorizontalAlignment','center');
text(-7,0,'M','FontSize',18,'HorizontalAlignment','center');
text(7,0,'L','FontSize',18,'HorizontalAlignment','center');
r = 0.3; % radius (complete guess) for spread of eaceh injection (0.5 µl)

% all depths (in mm) will be wrt tip of 16 mm guide tube

% first plot grid holes and cannulae inserted during perfusion
for X = -6:6
    for Y = -6:6
        circle(X,Y,0.32,'k-');
    end
end
plot3(0,0,0,'kx','MarkerSize',24);
[X,Y,Z] = cylinder(0.25,30);
surf(X-5,Y-5,Z*10-1);
surf(X-5,Y+5,Z*10-1);
surf(X+5,Y-5,Z*10-1);
surf(X+5,Y+5,Z*10-1);
set(gca,'ZDir','reverse'); grid on;

% track 1 [-5,3], 4 inj, 12/23/13
% rAAV1-hSyn-ChR2(H134R)-mCherry (undiluted)
x1 = [-5 -5 -5 -5]; % x grid position -5 mm
y1 = [-3 -3 -3 -3]; % y grid position -3 mm
z1 = 0.6:0.5:2.1; % 16 mm guide                     
plotSphere(x1,y1,z1,r,'r');
x=x1(1); y=y1(1); z=0; st='1';
text(x+0.1,y,z,st,'FontSize',15,'FontWeight','bold','HorizontalAlignment','center');

% track 2 ("new 2") [3,-1] (NOTE: differs from red dot in grid photo), 4 inj, 12/26/13
% rAAV1-hSyn-ChR2(H134R)-mCherry (diluted 1:2)
x2 = [3 3 3 3];
y2 = [-1 -1 -1 -1];
z2 = [1.3 1.8 4.3 4.8] - 0.75; % 15.25 mm guide, hence -0.75
plotSphere(x2,y2,z2,r,'r');
x=x2(1); y=y2(1); z=0; st='2';
text(x+0.1,y,z,st,'FontSize',15,'FontWeight','bold','HorizontalAlignment','center');

% track 3 [1,1], 5 inj, 12/27/13
% rAAV1-hSyn-ChR2(H134R)-mCherry (diluted 1:10)
x3 = [1 1 1 1 1];
y3 = [1 1 1 1 1];
z3 = (6.5:0.5:8.5) - 0.75; % 15.25 mm guide
plotSphere(x3,y3,z3,r,'r');
x=x3(1); y=y3(1); z=0; st='3';
text(x+0.1,y,z,st,'FontSize',15,'FontWeight','bold','HorizontalAlignment','center');

% track 4 [-2,-5], 5 inj, 12/20/13
% Lenti-CamKII-ArchT-GFP
x4 = [-2 -2 -2 -2 -2];
y4 = [-5 -5 -5 -5 -5];
z4 = 0.3:0.5:2.3; % 16 mm guide
plotSphere(x4,y4,z4,r,'g');
x=x4(1); y=y4(1); z=0; st='4';
text(x+0.1,y,z,st,'FontSize',15,'FontWeight','bold','HorizontalAlignment','center');

% track 5 [-2,-2], 3 inj, 12/19/13
% rAAV1-hSyn-ArchT-eYFP
x5 = [-2 -2 -2];
y5 = [-2 -2 -2];
z5 = [4.5 5 5.5]; % 16 mm guide
plotSphere(x5,y5,z5,r,'y');
x=x5(1); y=y5(1); z=0; st='5';
text(x+0.1,y,z,st,'FontSize',15,'FontWeight','bold','HorizontalAlignment','center');

% track 6 [-3,2], 4 inj, 12/18/13 -- 32 ga cannula
% rAAV5-hSyn-ChR2(H134R)-eYFP
x6 = [-3 -3 -3 -3];
y6 = [2 2 2 2];
z6 = 2.7:0.5:4.2; % 16 mm guide
plotSphere(x6,y6,z6,r,'y');
x=x6(1); y=y6(1); z=0; st='6';
text(x+0.1,y,z,st,'FontSize',15,'FontWeight','bold','HorizontalAlignment','center');

% track 7 [4,2], 8 inj, 12/19/13
% rAAV9-hSyn-ChR2(H134R)-eYFP
x7 = [4 4 4 4 4 4 4 4];
y7 = [2 2 2 2 2 2 2 2];
z7 = 3:0.5:6.5; % 16 mm guide
plotSphere(x7,y7,z7,r,'y');
x=x7(1); y=y7(1); z=0; st='7';
text(x+0.1,y,z,st,'FontSize',15,'FontWeight','bold','HorizontalAlignment','center');

% track 8 [1,-3], 4 inj, 12/20/13
% rAAV1-hSyn-ChR2(H134R)-mCherry + rAAV1-hSyn-ArchT-eYFP (mixture)
x8 = [1 1 1 1];
y8 = [-3 -3 -3 -3];
z8 = [0.5 1 3.5 4]; % 16 mm guide
plotSphere(x8,y8,z8,r,'m');
x=x8(1); y=y8(1); z=0; st='8';
text(x+0.1,y,z,st,'FontSize',15,'FontWeight','bold','HorizontalAlignment','center');

% track 9 [2,4], 1 inj, 12/18/13
% rAAV1-hSyn-ChR2(H134R)-mCherry (thawed, inj only)
x9 = 2;
y9 = 4;
z9 = 5-0.25; % 15.75 mm guide
plotSphere(x9,y9,z9,r,'r');
x=x9(1); y=y9(1); z=0; st='9';
text(x+0.1,y,z,st,'FontSize',15,'FontWeight','bold','HorizontalAlignment','center');

% track 10 ("A") [-2,5], 4 inj, 6/10/14 
% rAAV8-DREADD(gpcr)-mCherry
x10 = [-2 -2 -2 -2];
y10 = [5 5 5 5];
z10 = (5.5:0.5:7) + 0.5; % 16.5 mm guide
plotSphere(x10,y10,z10,r,'r');
x=x10(1); y=y10(1); z=0; st='10';
text(x+0.1,y,z,st,'FontSize',15,'FontWeight','bold','HorizontalAlignment','center');


%% pulvinar
% 4 tracks (labeled 11-14), all Lenti-CamKII-ArchT-GFP
% 10-degree angled grid, rotated 45-deg clockwise relative to standard grid,
% such that penetration angle is lateral and posterior wrt block orientation.
% gt depth is 17.7, or 1.7 longer than orig gt

figure; set(gcf, 'Color', 'w', 'Position', [400 100 640 820], 'PaperPositionMode', 'auto');
xlim([-7 7]); ylim([-7 7]); zlim([-1 22]); hold on; axis equal;
% xlabel('M-L (mm)','FontSize',15); ylabel('A-P (mm)','FontSize',15)
% zlabel('Depth (mm, wrt 16mm gt)','FontSize',15);
text(0,7,'A','FontSize',18,'HorizontalAlignment','center');
text(0,-7,'P','FontSize',18,'HorizontalAlignment','center');
text(-7,0,'M','FontSize',18,'HorizontalAlignment','center');
text(7,0,'L','FontSize',18,'HorizontalAlignment','center');
r = 0.3; % radius (complete guess) for spread of eaceh injection (0.5 µl)

% first plot grid holes and cannulae inserted during perfusion
for X = -6:6
    for Y = -6:6
        circle(X,Y,0.32,'k:');
    end
end
plot3(0,0,0,'kx','MarkerSize',24);
[X,Y,Z] = cylinder(0.25,30);
surf(X-5,Y-5,Z*23-1);
surf(X-5,Y+5,Z*23-1);
surf(X+5,Y-5,Z*23-1);
surf(X+5,Y+5,Z*23-1);
set(gca,'ZDir','reverse'); grid on;

% grid holes for pulvinar inj, in the coordinates of original grid
% (based on manual measurements)
xy14 = [sqrt(2) sqrt(2)/4];
xy13 = xy14 + [1/sqrt(2) 1/sqrt(2)];
xy12 = xy14 + [0 sqrt(2)];
xy11 = xy13 + [0 sqrt(2)];
z_gt = 1.7; % depth of guide tube tip
circle(xy11(1),xy11(2),z_gt,0.32,'k-');
circle(xy12(1),xy12(2),z_gt,0.32,'k-');
circle(xy13(1),xy13(2),z_gt,0.32,'k-');
circle(xy14(1),xy14(2),z_gt,0.32,'k-');

delta = 0.5*sin(10*pi/180); % grid angle is 10 deg...
dx = delta/sqrt(2); % and its x-y plane is rotated 45 degrees...
dy = -dx; % so this is the change in x,y for each 0.5 mm of depth (y is neg bc oriented posterior-lateral)
dz = 0.5*cos(10*pi/180); % also each z step is a bit less than 0.5

% track 11, 4 inj
z0 = 15.5; % microdrive depth of uppermost inj
x = xy11(1) + z0/0.5 * dx; %
y = xy11(2) + z0/0.5 * dy; % x,y,z pos of uppermost inj (offset from grid hole bc 10 deg angle)
z = z0/0.5 * dz;           %
% now calculate coords for all injs
x11 = [x x+dx x+2*dx x+3*dx] + z_gt/0.5 * dx;
y11 = [y y+dy y+2*dy y+3*dy] + z_gt/0.5 * dy;
z11 = [z z+dz z+2*dz z+3*dz] + z_gt/0.5 * dz;
plotSphere(x11,y11,z11,r,'g');
plot3([xy11(1) x+z_gt/0.5*dx], [xy11(2) y+z_gt/0.5*dy], [z_gt z+z_gt/0.5*dz], 'k-');
text(xy11(1),xy11(2),z_gt,'11','FontSize',15,'FontWeight','bold','HorizontalAlignment','center');

% track 12, 3 inj
z0 = 17.7;
x = xy12(1) + z0/0.5 * dx; 
y = xy12(2) + z0/0.5 * dy;
z = z0/0.5 * dz;
x12 = [x x+dx x+2*dx] + z_gt/0.5 * dx;
y12 = [y y+dy y+2*dy] + z_gt/0.5 * dy;
z12 = [z z+dz z+2*dz] + z_gt/0.5 * dz;
plotSphere(x12,y12,z12,r,'g');
plot3([xy12(1) x+z_gt/0.5*dx], [xy12(2) y+z_gt/0.5*dy], [z_gt z+z_gt/0.5*dz], 'k-');
text(xy12(1),xy12(2),z_gt,'12','FontSize',15,'FontWeight','bold','HorizontalAlignment','center');

% track 13, 7 inj
z0 = 16.3;
x = xy13(1) + z0/0.5 * dx; 
y = xy13(2) + z0/0.5 * dy;
z = z0/0.5 * dz;
x13 = [x x+dx x+2*dx x+3*dx x+4*dx x+5*dx x+6*dx] + z_gt/0.5 * dx;
y13 = [y y+dy y+2*dy y+3*dy y+4*dy y+5*dy y+6*dy] + z_gt/0.5 * dy;
z13 = [z z+dz z+2*dz z+3*dz z+4*dz z+5*dz z+6*dz] + z_gt/0.5 * dz;
plotSphere(x13,y13,z13,r,'g');
plot3([xy13(1) x+z_gt/0.5*dx], [xy13(2) y+z_gt/0.5*dy], [z_gt z+z_gt/0.5*dz], 'k-');
text(xy13(1),xy13(2),z_gt,'13','FontSize',15,'FontWeight','bold','HorizontalAlignment','center');

% track 14, 5 inj
z0 = 16.1;
x = xy14(1) + z0/0.5 * dx; 
y = xy14(2) + z0/0.5 * dy;
z = z0/0.5 * dz;
x14 = [x x+dx x+2*dx x+3*dx x+4*dx] + z_gt/0.5 * dx;
y14 = [y y+dy y+2*dy y+3*dy y+4*dy] + z_gt/0.5 * dy;
z14 = [z z+dz z+2*dz z+3*dz z+4*dz] + z_gt/0.5 * dz;
plotSphere(x14,y14,z14,r,'g');
plot3([xy14(1) x+z_gt/0.5*dx], [xy14(2) y+z_gt/0.5*dy], [z_gt z+z_gt/0.5*dz], 'k-');
text(xy14(1),xy14(2),z_gt,'14','FontSize',15,'FontWeight','bold','HorizontalAlignment','center');


% re-plot cortex inj from above
plotSphere(x1,y1,z1,r,'r');
plotSphere(x2,y2,z2,r,'r');
plotSphere(x3,y3,z3,r,'r');
plotSphere(x4,y4,z4,r,'g');
plotSphere(x5,y5,z5,r,'y');
plotSphere(x6,y6,z6,r,'y');
plotSphere(x7,y7,z7,r,'y');
plotSphere(x8,y8,z8,r,'m');
plotSphere(x9,y9,z9,r,'r');
plotSphere(x10,y10,z10,r,'r');

