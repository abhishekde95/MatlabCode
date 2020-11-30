% Analyzing cone weights 
% Author - Abhishek De, 12/18
close all; clearvars;
load S1LMS.mat
load S2LMS.mat
load S1RGB.mat
load S2RGB.mat
[OCidx, LMidx, LUMidx, SOidx, hardtoclassifyidx] = classifycells(S1LMS,S2LMS);
DOidx = [OCidx LMidx];

% Doing PCA on cone weights
LMS = [S1LMS;S2LMS];
[v,d] = eig(cov(LMS'));
k1 = v(:,end); k2 = v(:,end-1); k3 = v(:,end-2);
projLMS = [k1'*LMS; k2'*LMS; k3'*LMS];

% Doing PCA on gun weights
RGB = [S1RGB; S2RGB];
[v1,d1] = eig(cov(RGB'));
m1 = v1(:,end); m2 = v1(:,end-1); m3 = v1(:,end-2);
projRGB = [m1'*RGB; m2'*RGB; m3'*RGB];

figure(1),subplot(221);plot3(projLMS(1,DOidx),projLMS(2,DOidx),projLMS(3,DOidx),'o','MarkerFaceColor',[0 0 1]); hold on;
plot3(projLMS(1,LUMidx),projLMS(2,LUMidx),projLMS(3,LUMidx),'o','MarkerFaceColor',[0 1 0]);
plot3(projLMS(1,hardtoclassifyidx),projLMS(2,hardtoclassifyidx),projLMS(3,hardtoclassifyidx),'o','MarkerFaceColor',[1 0 0]); 
title('LMS 3D PCA'); hold off;
subplot(222); plot3(projRGB(1,DOidx),projRGB(2,DOidx),projRGB(3,DOidx),'o','MarkerFaceColor',[0 0 1]); hold on;
plot3(projRGB(1,LUMidx),projRGB(2,LUMidx),projRGB(3,LUMidx),'o','MarkerFaceColor',[0 1 0]);
plot3(projRGB(1,hardtoclassifyidx),projRGB(2,hardtoclassifyidx),projRGB(3,hardtoclassifyidx),'o','MarkerFaceColor',[1 0 0]); 
title('RGB 3D PCA'); hold off;
subplot(223);plot(projLMS(1,DOidx),projLMS(2,DOidx),'o','MarkerFaceColor',[0 0 1]); hold on;
plot(projLMS(1,LUMidx),projLMS(2,LUMidx),'o','MarkerFaceColor',[0 1 0]);
plot(projLMS(1,hardtoclassifyidx),projLMS(2,hardtoclassifyidx),'o','MarkerFaceColor',[1 0 0]); 
title('LMS 2D PCA'); hold off;
subplot(224); plot(projRGB(1,DOidx),projRGB(2,DOidx),'o','MarkerFaceColor',[0 0 1]); hold on;
plot(projRGB(1,LUMidx),projRGB(2,LUMidx),'o','MarkerFaceColor',[0 1 0]);
plot(projRGB(1,hardtoclassifyidx),projRGB(2,hardtoclassifyidx),'o','MarkerFaceColor',[1 0 0]); 
title('RGB 2D PCA'); hold off;
