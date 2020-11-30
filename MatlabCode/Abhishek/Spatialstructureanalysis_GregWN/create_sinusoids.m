% Creating sinusoids 
% Author - Abhishek De, 5/19
close all; clearvars;
[X, Y] = meshgrid(1:1:10,1:1:10);
phi = 0;
lambda = [2 4 8 12 16 20];
theta = linspace(0,360,17); theta(end) = [];
X = X-5; Y = Y+5; % So negative numbers mean down
count = 1;
figure(1)
for ii = 1:numel(lambda)
    for jj = 1:numel(theta)
        xprime = X.*cos(-theta(jj))+Y.*sin(-theta(jj));
        yprime = -X.*sin(-theta(jj))+Y.*cos(-theta(jj));
        sinusoid = cos((2.*pi.*yprime./lambda(ii))-phi);
        figure(1); subplot(numel(lambda),numel(theta),count); imagesc(sinusoid), set(gca,'XTick',[],'YTick',[]); axis square; colormap('gray');
        figure(2); subplot(numel(lambda),numel(theta),count); imagesc(abs(fftshift(fft2(sinusoid)))), set(gca,'XTick',[],'YTick',[]); axis square; colormap('gray')
        count = count + 1;
    end
end
