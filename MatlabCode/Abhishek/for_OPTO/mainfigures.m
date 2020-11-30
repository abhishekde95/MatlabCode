% figures for OPTO paper
clearvars; close all;
% Creating the gabor
figure(1);
L = 1001; W = 1001;
theta = 0; sigma = 0.2; ggamma = 1; sf = 1.0 ; phase = pi/2; nsigmas = 1; 


for jj = 1:numel(sf)
    [x,y] = meshgrid(linspace(-1,1,L), linspace(-1,1,W));
    X = x*cos(-theta) + y*sin(-theta);
    Y =-x*sin(-theta) + y*cos(-theta);
    expterm = exp(-(X.^2 + ggamma^2*Y.^2)/2/sigma^2);
    costerm = cos(2*pi*Y*sf(jj) + phase);
    im = costerm.*expterm;
    normfactor = max(abs(im(:)));
    im = (0.5*im/(normfactor)) + 0.5;
    image(cat(3,im,im,im)), drawnow, set(gca,'Xtick',[],'Ytick',[]); colormap(gray); axis square;
end

%% 
%% 
% Section 24
% Stimulus movie
close all; 
writerObj = VideoWriter('Gaborstim.mp4','MPEG-4');
open(writerObj);

theta = 0; sigma = 0.2; ggamma = 1; sf = 1.0 ; phase = pi/2; nsigmas = 1; pixperdeg = 2000;
tf = 8; 
time = 0.5;% in seconds
refreshrate = 75; % in Hz

ramp = linspace(0, 1, nframesramp);
plateau = ones(1, nframesplateau);
temporalprofile = [ramp plateau fliplr(ramp)];
nframes = length(temporalprofile);
stimsizeindeg = sigma*nsigmas;
stimsizeinpix = round(stimsizeindeg*pixperdeg);
[x,y] = meshgrid(4*stimsizeindeg*linspace(-1, 1, stimsizeinpix), ...
    4*stimsizeindeg*linspace(-1, 1, stimsizeinpix));

X = x*cos(-theta) + y*sin(-theta);
Y =-x*sin(-theta) + y*cos(-theta);

deltaphase = 2*pi*tf/refreshrate;
phases = 0 + (0:nframes-1)*deltaphase;
phases = reshape(phases, [1 1 nframes]);
temporalprofile = reshape(temporalprofile, [1 1 nframes]);
expterm = bsxfun(@times, exp(-(X.^2 + ggamma^2*Y.^2)/2/sigma^2), temporalprofile);
costerm = cos(bsxfun(@plus, 2*pi*Y*sf, phases));
template = expterm .* costerm;
normfactor = max(abs(template(:)));
template = 0.5 + 0.5*template/(normfactor+0.001);

for i = 1:nframes
    im = squeeze(template(:,:,i));
    image(cat(3,im,im,im)),set(gca,'Xtick',[],'Ytick',[],'Visible','off'); colormap(gray); axis square;
    axis image;
    frame = getframe(gca);
    writeVideo(writerObj,frame);
end

close(writerObj);


