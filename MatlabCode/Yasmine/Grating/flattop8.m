%% Yasmine El-Shamayleh (6/2014)
% Vingetted grating
% Flattop 8 Gaussian profile

%% BASIC WINDOWED GRATING


clear all; close all;

fig1 =figure('position',[1400 1600 800 1300]);
annotation('textbox',[.4 .89 .2 .1],'string',['Grating apertures a la Movshon'],'HorizontalAlignment','Center','FitBoxToText','on','EdgeColor','none')
annotation('textbox',[.4 .88 .2 .1],'string',['Gaussian flattop w/ 25% dropoff'],'HorizontalAlignment','Center','FitBoxToText','on','EdgeColor','none')


freq = 8;
theta = -180:1:180;
sine1D = sind(freq*theta);
sine2D = ones(size(sine1D))'*sine1D;


npix = size(sine2D,1);                                          %number of pixels of image
R = 1;                                                          %radius of stimulus
[x,y] = meshgrid(linspace(-1,1,npix),linspace(-1,1,npix));
[~,r] = cart2pol(x,y);
cAp = zeros(npix,npix);                                         %circular aperture
cAp(r<=(R * 1)) = 1;

subplot(4,2,1)
imagesc(cAp); axis xy; axis square; colormap(gray);
set(gca,'xtick',[],'ytick',[],'tickDir','out','xlim',[0 360],'ylim',[0 360]); axis square;
title(sprintf('circular aperture'));

subplot(4,2,2);
gratC = sine2D.*cAp;
imagesc(gratC); axis xy; axis square;
set(gca,'xtick',[],'ytick',[],'tickDir','out','xlim',[0 360],'ylim',[0 360]); axis square;
title(sprintf('grating x circular aperture'));


%% FLATTOP

% Romesh Kumbhani 6/2014
% generic flattop, where img is your 2D image, and amount is
% the type of flattop. e.g. flattop(mypic,8) is a flattop 8 version.
% flattop(mypic,4) is a flattop4... etc..

img = gratC;         %circular grating
a = [20 8 4];       %amount of windowing

for i = 1:length(a)
    
    amount = a(i);
    
    thetype = class(img);
    img = double(img);
    [npixX,npixY] = size(img);
    [x,y] = meshgrid(linspace(-1,1,npixX),linspace(-1,1,npixY));
    [~,r] = cart2pol(x,y);
    gAp = exp(-0.5*((r-(1-(2/amount)))/(1/(3*(amount/2)))).^2);
    gAp(r<=(1-(2/amount))) = 1;
    gratG = eval(sprintf('%s(img.*gAp)',thetype));
    
    subplot(4,2,(2*i)-1+2)
    imagesc(gAp); axis square;
    set(gca,'xtick',[],'ytick',[],'tickDir','out','xlim',[0 360],'ylim',[0 360]); axis square;
    text =[];
    text{1} = (sprintf('flattop Gaussian'));
    text{2} = (sprintf('drop off at outer %.0f %%', 100*(2/amount)));
    title(text);

    subplot(4,2,(2*i)-1+3)
    imagesc(gratG); axis square;
    set(gca,'xtick',[],'ytick',[],'tickDir','out','xlim',[0 360],'ylim',[0 360]); axis square;
    text=[];
    text{1} = (sprintf('grating x flattop Gaussian'));
    text{2} = (sprintf('drop off at outer %.0f %%', 100*(2/amount)));
    title(text)
    
end
