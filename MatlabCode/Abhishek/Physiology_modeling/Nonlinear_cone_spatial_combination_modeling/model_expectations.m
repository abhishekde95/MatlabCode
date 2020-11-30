% Started writing this script just to get a good hold on my concepts
% Author - Abhishek De

close all; clearvars;

% Model Expectations 3-D (More applicable for understanding the how cone signals are combined) 

[x y z] = meshgrid([-3:.1:3],[-3:.1:3],[-3:.1:3]);
fnvals = 0.1;

% Ellipsoid
fn = (x.^2/50)+(y.^2/100)+(z.^2/2);
v = isosurface(x,y,z,fn,fnvals/2);
figure(2),subplot(221);h = patch(v);set(h,'FaceColor',[1 0 0],'FaceAlpha',.5); xlabel('X'), xlabel('Y'), xlabel('Z'); title('Ellipse'); view(10,50);

% Hyperboloid of 1 sheet
fn = (x.^2/50)-(y.^2/100)+(z.^2/2);
v = isosurface(x,y,z,fn,fnvals/2);
figure(2),subplot(222);h = patch(v);set(h,'FaceColor',[0.5 0.5 0],'FaceAlpha',.5); xlabel('X'), xlabel('Y'), xlabel('Z'); title('Hyperboloid (1 sheet)'); view(40,70);

% Hyperboloid of 2 sheets
fn = (x.^2/50)-(y.^2/100)+(z.^2/2);
v = isosurface(x,y,z,fn,-1*fnvals/10);
figure(2),subplot(223);h = patch(v);set(h,'FaceColor',[0 1 0],'FaceAlpha',.5); xlabel('X'), xlabel('Y'), xlabel('Z');title('Hyperboloid (2 sheets)'); view(70,50);

% Parallel planes
fn = abs((x./50)+(y/100)+(z/2));
v = isosurface(x,y,z,fn,0.1);
figure(2),subplot(224);h = patch(v);set(h,'FaceColor',[0 0 1],'FaceAlpha',.5); xlabel('X'), xlabel('Y'), xlabel('Z'); title('Parallel Planes'); view(30,10);

  