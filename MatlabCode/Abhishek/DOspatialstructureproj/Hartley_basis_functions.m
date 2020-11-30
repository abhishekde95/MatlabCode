% Here I want simulate the V1 like neuronal responses to Hartley Basis functions
% Author - Abhishek De, 12/18
close all; clearvars;
wDeg = 1; % size of the image in degrees
nPix = 50; % Resoution of image (pixels)
plot_counter = 1;
[x,y] = meshgrid(linspace(-wDeg/2,wDeg/2,nPix+1));
x = x(1:end-1,1:end-1);
y = y(1:end-1,1:end-1);
ori = linspace(0,135,8);
sf = [0.5 1 2 4 8];
phases = linspace(0,135,8);
R = numel(ori);
C = numel(sf);
grating = cell(numel(phases)*R*C,1);
for kk = 1:numel(phases)
    count = 1;
    for ii = 1:numel(ori)
        ramp = sin(ori(ii)*pi/180)*x - cos(ori(ii)*pi/180)*y;
        for jj = 1:numel(sf)
            grating{(kk-1)*numel(ori)*numel(sf)+(ii-1)*numel(sf)+jj,1} = sin(2*pi*sf(jj)*ramp + (phases(kk)*pi/180));
            count = count + 1;
        end
    end
end

% Till this point I have created Hartley basis functions and now I am creating a RF and running it through a stimulus train
iter = 5;
N = 2000;
stim_train = randi(numel(phases)*R*C,[N 1]);
RefreshRate = 75;% Stim refresh rate (Hz)
dtbin = 0.1; % binsize for Poisson spike generation
sigma = 0.1;
for jj = 1:iter
    RF = grating{randi(numel(phases)*R*C),1}.*exp(-(x.^2 + y.^2)/2/sigma^2); % creating a gabor like RF
    resp = zeros(N,1);
    stim = [];
    for ii = 1:N
        k = grating{stim_train(ii),1};
        tmpresp = k(:)'*RF(:);
        r = max(0,tmpresp); % cumulative drive
        rbig = repmat(r/RefreshRate*dtbin,1./dtbin,1); % make Poisson spike train
        sp = sum(rand(size(rbig))<rbig)';
        resp(ii) = sp;
        stim = [stim k(:)];
    end
    spikeind = find(resp>0);
    nspikes = sum(resp);
    STS = sum(stim(:,spikeind).*repmat(resp(spikeind)',[nPix^2,1]),2); % spike triggered ensemble
    STA = STS./nspikes;
    STA = reshape(STA,[nPix nPix]);
    figure(plot_counter); subplot(iter,2,(jj-1)*2+1); imagesc(RF); colormap('gray'); axis square; set(gca,'XTick',[],'YTick',[]);
    subplot(iter,2,(jj-1)*2+2); imagesc(STA); colormap('gray'); axis square; set(gca,'XTick',[],'YTick',[]);
end

