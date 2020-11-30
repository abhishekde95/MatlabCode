%%
%define the constants
nPix = (stro.sum.exptParams.flash_size*2) * (unique(stro.trial(:, stro.idx.gaborSigma))/10) * stro.sum.exptParams.pixperdeg;
halfSize = round(nPix/2);
gaborTheta = unique(stro.trial(:,stro.idx.gaborTheta));
rgbTrial = [0.3 0.2 0.6]; %picking randomly for now... but we need to detemine this later.
bkgndrgb = [stro.sum.exptParams.bkgnd_r, stro.sum.exptParams.bkgnd_g, stro.sum.exptParams.bkgnd_b];
gaborGamma = unique(stro.trial(:, stro.idx.gaborGamma));
gaborLambda = unique(stro.trial(:,stro.idx.gaborLambda));
gaborSigma = unique(stro.trial(:,stro.idx.gaborSigma)) / 10 * stro.sum.exptParams.pixperdeg;
gaborPhase = 0;
flashDeltaPhase = unique(stro.trial(:,stro.idx.driftRate)) .* 2 .* pi .* (1 ./ stro.sum.exptParams.frame_rate); %the amount to advance each frame
flashLength = stro.sum.exptParams.flash_length;
invGamma = InvertGamma(reshape(stro.sum.exptParams.gamma_table,[],3),1);
maxdac = (2^16)-1;          %65536 dac values in colour mode

flashNumFrames = ceil(stro.sum.exptParams.frame_rate .* (flashLength./1000));
rampLength = ceil(flashNumFrames ./ 4);
ramp = linspace(0, 1, rampLength);  %ramp is 1/4th of the total duration on either side
plateau = ones(1,flashNumFrames - (rampLength .* 2));
flashTimeProf = [ramp, plateau, fliplr(ramp)];

row = -halfSize:halfSize-1; %subtract one so that you don't overlap the texture window;
col = -halfSize:halfSize-1;
[X, Y] = meshgrid(row, col);
xprime = X .* cos(-gaborTheta) + Y .* sin(-gaborTheta);
yprime = -X .* sin(-gaborTheta) + Y .* cos(-gaborTheta);
rgb_increment = rgbTrial - bkgndrgb;

for a = 1:flashNumFrames;
    
    %start by making the gabor, then multiplying by the increment from
    %background and by the temporal weighting fxn.
    gabor = exp(-(xprime.^2 + gaborGamma.^2 .* yprime.^2) ./ (2.* gaborSigma.^2)) .* cos(2 .* pi .* yprime ./ gaborLambda + gaborPhase);
    flashImg = repmat(gabor, [1, 1, 3]);
    flashImg(:,:,1) = (flashImg(:,:,1) .* flashTimeProf(a) .* rgb_increment(1)) + bkgndrgb(1);
    flashImg(:,:,2) = (flashImg(:,:,2) .* flashTimeProf(a) .* rgb_increment(2)) + bkgndrgb(2);
    flashImg(:,:,3) = (flashImg(:,:,3) .* flashTimeProf(a) .* rgb_increment(3)) + bkgndrgb(3);
    
    %now convert to DAC values
    imgSize = size(gabor);
    flashrgb = round(maxdac .* flashImg) + 1; %intensities b/w 1 & maxdac +1
    flashRGB = ones(size(flashImg)); %preallocate so that indexing works
    flashRGB(:,:,1) = reshape(round(maxdac .* invGamma(flashrgb(:,:,1), 1)), imgSize(1), imgSize(2));
    flashRGB(:,:,2) = reshape(round(maxdac .* invGamma(flashrgb(:,:,2), 2)), imgSize(1), imgSize(2));
    flashRGB(:,:,3) = reshape(round(maxdac .* invGamma(flashrgb(:,:,3), 3)), imgSize(1), imgSize(2));
    
    %update the phase
    gaborPhase = gaborPhase + flashDeltaPhase;
    unique(flashRGB)
    keyboard
    %plot for fun
    imagesc(flashImg)
    pause(0.05)
end