% This is for generating figures relevant for rotation talk

close all; clear all;
nbins = 0.0:0.001:0.5;
[X,Y] = meshgrid(nbins',nbins);
lin_model = X + Y;

figure(1),
imagesc([min(nbins) max(nbins)],[min(nbins) max(nbins)], lin_model); 
xlabel('A'); ylabel('B');

nbins1 = nbins;
nbins1(nbins1<0) = 0;
% nbins1 = nbins1.^2;
quad_model = max(sigmoid(repmat(nbins1',[1 numel(nbins1)]),20,0.4),sigmoid(repmat(nbins1,[numel(nbins1) 1]),20,0.4));


figure(2),
imagesc([min(nbins) max(nbins)],[min(nbins) max(nbins)], quad_model); 
xlabel('A'); ylabel('B');

nbins2 = nbins;
nbins2(nbins2<0) = 0;
% nbins2 = nbins2.^(0.5);
% sqrt_model = repmat(nbins2',[1 numel(nbins2)]) + repmat(nbins2,[numel(nbins2) 1]);
sqrt_model = min(sigmoid(repmat(nbins2',[1 numel(nbins2)]),10,0.4),sigmoid(repmat(nbins2,[numel(nbins2) 1]),10,0.4));


figure(3),
imagesc([min(nbins) max(nbins)],[min(nbins) max(nbins)], sqrt_model); 
xlabel('A'); ylabel('B');

figure(4),
imagesc([min(nbins) max(nbins)],[min(nbins) max(nbins)], zeros(size(X))); colormap('gray');
xlabel('A'); ylabel('B');

figure(5),
imagesc([min(nbins) max(nbins)],[min(nbins) max(nbins)], lin_model.^3); 
xlabel('A'); ylabel('B');