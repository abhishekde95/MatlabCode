% main for FVM 
clearvars; close all;
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S 
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = spline(380:4:780, mon_spd', 380:5:780); % fitting a cubic spline
M = mon_spd*fundamentals; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals
M = inv(M');

% Create the gabor
LMS = [0.06 -0.06 0.45; 0.06 -0.06 -0.45]; % L,M,S
figure(1);
L = 100; W = 100;
theta = 0; sigma = 0.2; ggamma = 1; sf = [0.5 3.0] ; phase = pi/2; 
for jj = 1:numel(sf)
    [x,y] = meshgrid(linspace(-1,1,L), linspace(-1,1,W));
    X = x*cos(-theta) + y*sin(-theta);
    Y =-x*sin(-theta) + y*cos(-theta);
    expterm = exp(-(X.^2 + ggamma^2*Y.^2)/2/sigma^2);
    costerm = cos(2*pi*Y*sf(jj) + phase);
    template = costerm.*expterm;
    for ii = 1:size(LMS,1)
        RGB = (LMS(ii,:)/norm(LMS(ii,:)))*M;
        im = cat(3,RGB(1)*template,RGB(2)*template,RGB(3)*template);
        normfactor = 1.01*max(abs(im(:)));
        im = (0.5*im/(normfactor)) + 0.5;
        subplot(2,2,2*(jj-1)+ii),image(im), drawnow, set(gca,'Xtick',[],'Ytick',[]); hold on;
    end
end
hold off;