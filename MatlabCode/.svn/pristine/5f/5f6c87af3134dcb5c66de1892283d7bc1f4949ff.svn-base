% trogdor_injection.m
% 3D plot of injection sites
% CRF 11-23-14

% cortex
% 1. rAAV1-hSyn-ChR2(H134R)-mCherry (undiluted)
% 2. rAAV1-hSyn-ChR2(H134R)-mCherry (diluted 1:2)
% 3. rAAV1-hSyn-ChR2(H134R)-mCherry (diluted 1:10)
% 4. Lenti-CamKII-ArchT-GFP
% 5. rAAV1-hSyn-ArchT-eYFP
% 6. rAAV5-hSyn-ChR2(H134R)-eYFP
% 7. rAAV9-hSyn-ChR2(H134R)-eYFP
% 8. rAAV1-hSyn-ChR2(H134R)-mCherry + rAAV1-hSyn-ArchT-eYFP (mixture)
% 9. rAAV1-hSyn-ChR2(H134R)-mCherry (thawed, inj only)
% 10 ("A"). rAAV8-DREADD(gpcr)-mCherry

%  pulvinar
% 11. 
% 12. 


clear all;

figure; set(gcf, 'Color', 'w', 'Position', [400 100 480 480], 'PaperPositionMode', 'auto');
xlim([-7 7]); ylim([-7 7]); zlim([-1 9]); hold on; axis equal;
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

% track 2 ("new 2") [3,-1] (NOTE: differs from red dot in grid photo), 4 inj, 12/26/13
% rAAV1-hSyn-ChR2(H134R)-mCherry (diluted 1:2)
x2 = [3 3 3 3];
y2 = [-1 -1 -1 -1];
z2 = [1.3 1.8 4.3 4.8] - 0.75; % 15.25 mm guide, hence -0.75
plotSphere(x2,y2,z2,r,'r');

% track 3 [1,1], 5 inj, 12/27/13
% rAAV1-hSyn-ChR2(H134R)-mCherry (diluted 1:10)
x3 = [1 1 1 1 1];
y3 = [1 1 1 1 1];
z3 = (6.5:0.5:8.5) - 0.75; % 15.25 mm guide
plotSphere(x3,y3,z3,r,'r');

% track 4 [-2,-5], 5 inj, 12/20/13
% Lenti-CamKII-ArchT-G(?)FP
x4 = [-2 -2 -2 -2 -2];
y4 = [-5 -5 -5 -5 -5];
z4 = 0.3:0.5:2.3; % 16 mm guide
plotSphere(x4,y4,z4,r,'g');

% track 5 [-2,-2], 3 inj, 12/19/13
% rAAV1-hSyn-ArchT-eYFP
x5 = [-2 -2 -2];
y5 = [-2 -2 -2];
z5 = [4.5 5 5.5]; % 16 mm guide
plotSphere(x5,y5,z5,r,'g');

% track 6 [-3,2], 4 inj, 12/18/13 -- 32 ga cannula
% rAAV5-hSyn-ChR2(H134R)-eYFP
x6 = [-3 -3 -3 -3];
y6 = [2 2 2 2];
z6 = 2.7:0.5:4.2; % 16 mm guide
plotSphere(x6,y6,z6,r,'g');

% track 7 [4,2], 8 inj, 12/19/13
% rAAV9-hSyn-ChR2(H134R)-eYFP
x7 = [4 4 4 4 4 4 4 4];
y7 = [2 2 2 2 2 2 2 2];
z7 = 3:0.5:6.5; % 16 mm guide
plotSphere(x7,y7,z7,r,'g');

% track 8 [1,-3], 4 inj, 12/20/13
% rAAV1-hSyn-ChR2(H134R)-mCherry + rAAV1-hSyn-ArchT-eYFP (mixture)
x8 = [1 1 1 1];
y8 = [-3 -3 -3 -3];
z8 = [0.5 1 3.5 4]; % 16 mm guide
plotSphere(x8,y8,z8,r,'m');

% track 9 [2,4], 1 inj, 12/18/13
% rAAV1-hSyn-ChR2(H134R)-mCherry (thawed, inj only)
x9 = 2;
y9 = 4;
z9 = 5-0.25; % 15.75 mm guide
plotSphere(x9,y9,z9,r,'r');

% track 10 ("A") [-2,5], 4 inj, 6/10/14 
% rAAV8-DREADD(gpcr)-mCherry
x10 = [-2 -2 -2 -2];
y10 = [5 5 5 5];
z10 = (5.5:0.5:7) + 0.5; % 16.5 mm guide
plotSphere(x10,y10,z10,r,'r');


