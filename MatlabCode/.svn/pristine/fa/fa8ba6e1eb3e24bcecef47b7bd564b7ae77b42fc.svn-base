% This code works in conjuction with Sig2Noise_RiekeHorwitzCollab.m 
% It should return the stim intensity (possibly also gun voltage) at each point in space and time
% Created 5/25/2011 by CH and JPW

function [flashrgb,gaussian] = stimDynamics(rgbTrial,stro,Stim)

% Define the constants
nPix = (stro.sum.exptParams.flash_size*2) * (unique(stro.trial(:, stro.idx.gaborSigma))/10) * stro.sum.exptParams.pixperdeg;
halfSize = round(nPix/2); %does rounding matter here?
gaborTheta = unique(stro.trial(:,stro.idx.gaborTheta));
gaborGamma = unique(stro.trial(:, stro.idx.gaborGamma));
gaborLambda = unique(stro.trial(:,stro.idx.gaborLambda));
gaborSigma = unique(stro.trial(:,stro.idx.gaborSigma)) / 10 * stro.sum.exptParams.pixperdeg;
gaborPhase = 0;
driftRate = unique(stro.trial(:,stro.idx.driftRate));
frameRate = stro.sum.exptParams.frame_rate;
flashDeltaPhase = driftRate * 2 * pi * (1 / frameRate); %the amount to advance each frame
flashLength = stro.sum.exptParams.flash_length;
invGamma = InvertGamma(reshape(stro.sum.exptParams.gamma_table,[],3),1);
maxdac = (2^16)-1;          %65536 dac values in colour mode


% Parameters for ramping on/off stimulus
flashNumFrames = ceil(stro.sum.exptParams.frame_rate * (flashLength/1000)); %Is this rounding ok?
rampLength = ceil(flashNumFrames / 4); %Is this rounding ok?
ramp = linspace(0, 1, rampLength);  %ramp is 1/4th of the total duration on either side
plateau = ones(1,flashNumFrames - (rampLength .* 2));
flashTimeProf = [ramp, plateau, fliplr(ramp)];

% Gabor coordinates
row = -halfSize:halfSize-1; %subtract one so that you don't overlap the texture window;
col = -halfSize:halfSize-1;
[X, Y] = meshgrid(row, col);
xprime = X .* cos(-gaborTheta) + Y .* sin(-gaborTheta);
yprime = -X .* sin(-gaborTheta) + Y .* cos(-gaborTheta);
rgb_increment = rgbTrial - Stim.bkgndrgb;

% Make gabor, then vary it with time
gaussian = exp(-(xprime.^2 + gaborGamma.^2 .* yprime.^2)./ (2.* gaborSigma.^2));
gabor = gaussian  .* cos(2 .* pi .* yprime ./ gaborLambda + gaborPhase);

% Preallocate variables
flashImg.R = nan(size(gabor,1),size(gabor,1),flashNumFrames);
flashImg.G = nan(size(gabor,1),size(gabor,1),flashNumFrames);
flashImg.B = nan(size(gabor,1),size(gabor,1),flashNumFrames);
flashImg.R(:,:,1) = gabor;
flashImg.G(:,:,1) = gabor;
flashImg.B(:,:,1) = gabor;
flashRGB.R = nan(size(flashImg.R)); 
flashRGB.G = nan(size(flashImg.G));
flashRGB.B = nan(size(flashImg.B));

% Loop through each frame refresh, accomodate for retina sampling rate
for t = 1:flashNumFrames;
   
    % Multiply the gabor by the increment from background and by the temporal weighting fxn.
    gabor = gaussian .* cos(2 .* pi .* yprime ./ gaborLambda + gaborPhase);

    flashImg.R(:,:,t) = (gabor .* flashTimeProf(t) .* rgb_increment(1)) + Stim.bkgndrgb(1);
    flashImg.G(:,:,t) = (gabor .* flashTimeProf(t) .* rgb_increment(2)) + Stim.bkgndrgb(2);
    flashImg.B(:,:,t) = (gabor .* flashTimeProf(t) .* rgb_increment(3)) + Stim.bkgndrgb(3);
   
    % Convert to DAC values
    imgSize = size(gabor);
    flashrgb.R(:,:,t) = round(maxdac .* flashImg.R(:,:,t)) + 1; %intensities b/w 1 & maxdac +1
    flashrgb.G(:,:,t) = round(maxdac .* flashImg.G(:,:,t)) + 1; %intensities b/w 1 & maxdac +1
    flashrgb.B(:,:,t) = round(maxdac .* flashImg.B(:,:,t)) + 1; %intensities b/w 1 & maxdac +1
    flashRGB.R(:,:,t) = reshape(round(maxdac .* invGamma(flashrgb.R(:,:,t), 1)), imgSize(1), imgSize(2));
    flashRGB.G(:,:,t) = reshape(round(maxdac .* invGamma(flashrgb.G(:,:,t), 2)), imgSize(1), imgSize(2));
    flashRGB.B(:,:,t) = reshape(round(maxdac .* invGamma(flashrgb.B(:,:,t), 3)), imgSize(1), imgSize(2));
    
    % Update the phase
    gaborPhase = gaborPhase + flashDeltaPhase;

    % Plot for fun
    figure(1)
    Image = cat(3,flashImg.R(:,:,t),flashImg.G(:,:,t),flashImg.B(:,:,t));
    imagesc(Image)
    pause(0.0013)
end

flashrgb.R = flashrgb.R ./ maxdac;
flashrgb.G = flashrgb.G ./ maxdac;
flashrgb.B = flashrgb.B ./ maxdac;