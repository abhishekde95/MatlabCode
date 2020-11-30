%% Parent Flash Signal-to-Noise Script
% Created Apr_2011 Angueyra
% Modified June_2011 JPW

close all
clear all

%% Inputs and Parameters

% User defined parameters
Stim.SamplingRate = 975; % in Hz - for simplicity, made comensurable with monitor refresh rate. JPW
SampLBookend = 200;
trim = 0; %Reduce gabor AOI (entry must be <1. This amount is taken off both gabor dimensions)
repetitions = 200;
savefile = 'Responses_NoSpatialGaussWeight';

% Input file(s)
stro = dtobj('./K020311002.nex');
%stro = dtobj('K020311002.nex');


%% Auto-Defined Parameters

h = waitbar(0, 'Loading Parameters...');

Stim.Area = (unique(stro.trial(:,stro.idx.gaborSigma))/10*stro.sum.exptParams.flash_size*2)^2; %in VA degrees^2, JPW
Stim.bkgndrgb = [stro.sum.exptParams.bkgnd_r, stro.sum.exptParams.bkgnd_g, stro.sum.exptParams.bkgnd_b]; %JPW
x = 0:255; %the normal range of the gamma look up table (CH)
xx = linspace(0, 255, 2^16); %the desired quantization of the gammaTable (CH)
g1 = reshape(stro.sum.exptParams.gamma_table, 256, 3); % CH
gammaTable = [spline(x, g1(:,1), xx)', spline(x, g1(:,2), xx)', spline(x, g1(:,3), xx)']; %CH/JPW
RGB = unique(stro.trial(:,[stro.idx.flashR,stro.idx.flashG,stro.idx.flashB]),'rows')+1; %CH/JPW (Voltage of stimuli(0:2^16);+1 for matlab indexing)
rgb = [gammaTable(RGB(1:end,1), 1), gammaTable(RGB(1:end,2), 2), gammaTable(RGB(1:end,3), 3)]; %CH/JPW (Convert voltage(0:2^16) to intensity(0:1); first RGB is bkgnd)
M = reshape(stro.sum.exptParams.m_mtx, 3, 3); %JPW/CH (converts from gun intensity to cone intensity
Stim.lms = rgb * M;
mon_spect = reshape(stro.sum.exptParams.mon_spect,[],3); %JPW mon_spect is already in Joules/s
Stim.FrameRate = stro.sum.exptParams.frame_rate;
Stim.Wavelength = (380:5:780)'; % From GH Readme in the dropbox
Stim.Spectra.R = mon_spect(:,1);
Stim.Spectra.G = mon_spect(:,2);
Stim.Spectra.B = mon_spect(:,3);
w = 0; % Progress of waitbar


%% Create Cone Grid

waitbar(0,h,'Performing Initial Calculations...');

StimLengthInPixels = round((stro.sum.exptParams.flash_size*2) * (unique(stro.trial(:, stro.idx.gaborSigma))/10) * stro.sum.exptParams.pixperdeg);
Stim.Pixel.Area = (sqrt(Stim.Area))/StimLengthInPixels; % deg va per pixel
center_x = stro.sum.exptParams.rf_x/10; % in deg va
center_y = stro.sum.exptParams.rf_y/10; % in deg va
loc_x = linspace(center_x-(sqrt(Stim.Area)/2),center_x+(sqrt(Stim.Area)/2),StimLengthInPixels*2+1);%finding half pixels
loc_x = loc_x(2:2:end); %Pixel centers
loc_y = linspace(center_y-(sqrt(Stim.Area)/2),center_y+(sqrt(Stim.Area)/2),StimLengthInPixels*2+1);
loc_y = loc_y(2:2:end);
Stim.Pixel.Eccentricity = nan(StimLengthInPixels);
for x = 1:StimLengthInPixels
    for y = 1:StimLengthInPixels
        Stim.Pixel.Eccentricity(x,y) = sqrt(loc_x(x)^2 + loc_y(y)^2);
        Cones.Pixels(x,y).Centers = ConeGrid(Stim.Pixel.Area,Stim.Pixel.Eccentricity(x,y)); %in mm of retina
    end
end
Cones.CollectingArea = 0.7 * 1e-6; %in mm^2 -JA
Cones.Size = 2e-3; %in mm (cell body diameter) -JA

%Trimming down periphery of stimulus by specified amount
toss = round(size(Cones.Pixels,1)*trim/2);
Cones.Pixels = Cones.Pixels(toss+1:end-toss,toss+1:end-toss);
Stim.Pixel.Eccentricity = Stim.Pixel.Eccentricity(toss+1:end-toss,toss+1:end-toss);


%% Calibration and Linear Responses

% Normalized Cone Spectral Sensitivities across spectral stimulus wavelength
SpectrumLCone = GeneratePhotoreceptorSpectrum('LCone',Stim.Wavelength);
SpectrumMCone = GeneratePhotoreceptorSpectrum('MCone',Stim.Wavelength);
SpectrumSCone = GeneratePhotoreceptorSpectrum('SCone',Stim.Wavelength);

StimEnergySpectraR = RadianceToRetIrradiance(Stim.Spectra.R,[],3,20)*.5;
StimEnergySpectraG = RadianceToRetIrradiance(Stim.Spectra.G,[],3,20)*.5;
StimEnergySpectraB = RadianceToRetIrradiance(Stim.Spectra.B,[],3,20)*.5;
% Dividing by 2 because cone size is .5 um^2
% Value for pupilAreaMM comes from Ostrin and Glasser, 'The Effects of Phenylephrine on Pupil Diameter and Accommodation in Rhesus Monkeys' Inves. Opth. and Vis. Sci., 2004
% Value for eyeSizeMM comes from Wilson and Gammon, 'Abnormol Development of the Axiol Length of Aphokic Monkey Eyes', Inves. Opth. and Vis. Sci., 1987

% Transforming to photons across wavelengths: Photons(lambda)=Energy(lambda)*lambda/hc
ScFact = 6.6e-34 * 3e8; %Planck's constant * Speed of Light
PhotonSpectraR = StimEnergySpectraR .* Stim.Wavelength * 1e-9 / ScFact; % (wavelength is in nm)
PhotonSpectraG = StimEnergySpectraG .* Stim.Wavelength * 1e-9 / ScFact;
PhotonSpectraB = StimEnergySpectraB .* Stim.Wavelength * 1e-9 / ScFact;

% Integrating across wavelength to get total number of photons
FluxLConeR = sum(PhotonSpectraR .* SpectrumLCone); %Should be ballpark 1000 (300 seems dim, 10,000 seems too bright)
FluxLConeG = sum(PhotonSpectraG .* SpectrumLCone);
FluxLConeB = sum(PhotonSpectraB .* SpectrumLCone);
FluxMConeR = sum(PhotonSpectraR .* SpectrumMCone);
FluxMConeG = sum(PhotonSpectraG .* SpectrumMCone);
FluxMConeB = sum(PhotonSpectraB .* SpectrumMCone);
FluxSConeR = sum(PhotonSpectraR .* SpectrumSCone);
FluxSConeG = sum(PhotonSpectraG .* SpectrumSCone);
FluxSConeB = sum(PhotonSpectraB .* SpectrumSCone);

clear Photon* ScFact StimEnergy* center* loc*


%% Grand Loop
%Preallocate Space and Loop Over Each Pixel, Then Frame Refresh, Then Resample

Stim.SampPerFrRef = ones(1,round(Stim.SamplingRate/stro.sum.exptParams.frame_rate));

% Loop through each stimulus
for s = 1:size(rgb,1)
        
    % Find power of each pixel at each frame refresh (rewrites for each stim)
    [Stim.Power,Stim.Pixel.Gaussian] = stimDynamics(rgb(s,:),stro,Stim); % Retuns a matrix at each frame refresh
    
    % Trim down stimulus to center pixels only (assumes stim is a square)
    Stim.Power.R = Stim.Power.R(toss+1:end-toss,toss+1:end-toss,:);
    Stim.Power.G = Stim.Power.G(toss+1:end-toss,toss+1:end-toss,:);
    Stim.Power.B = Stim.Power.B(toss+1:end-toss,toss+1:end-toss,:);
    Stim.Pixel.Gaussian = Stim.Pixel.Gaussian(toss+1:end-toss,toss+1:end-toss,:);
    
    % Reshape matrices to put all pixels in single frame refresh into one vector
    Stim.Power.R = shiftdim(reshape(Stim.Power.R,1,[],size(Stim.Power.R,3)));
    Stim.Power.G = shiftdim(reshape(Stim.Power.G,1,[],size(Stim.Power.G,3)));
    Stim.Power.B = shiftdim(reshape(Stim.Power.B,1,[],size(Stim.Power.B,3)));
    Stim.Pixel.Gaussian = reshape(Stim.Pixel.Gaussian,[],1);
    Stim.Pixel.Eccentricity = reshape(Stim.Pixel.Eccentricity,[],1);
    Cones.Pixels = reshape(Cones.Pixels,[],1);
    
    % Loop each pixel through time
    for p = 1:length(Cones.Pixels)
        
        % Stimulus area on retina
        RetinalStimArea = VisAngle2Retina_Area(Stim.Pixel.Eccentricity(p),Stim.Pixel.Area);% in mm^2
        
        % Convert From Frame Refresh Rate to Cone Sampling Rate
        LConeRStarR = Stim.SampPerFrRef .* FluxLConeR;
        LConeRStarG = Stim.SampPerFrRef .* FluxLConeG;
        LConeRStarB = Stim.SampPerFrRef .* FluxLConeB;
        MConeRStarR = Stim.SampPerFrRef .* FluxMConeR;
        MConeRStarG = Stim.SampPerFrRef .* FluxMConeG;
        MConeRStarB = Stim.SampPerFrRef .* FluxMConeB;
        SConeRStarR = Stim.SampPerFrRef .* FluxSConeR;
        SConeRStarG = Stim.SampPerFrRef .* FluxSConeG;
        SConeRStarB = Stim.SampPerFrRef .* FluxSConeB;
        
        TimeVectorL=nan(1,ceil(Stim.SamplingRate/Stim.FrameRate*size(Stim.Power.R,2)));
        TimeVectorM=nan(1,ceil(Stim.SamplingRate/Stim.FrameRate*size(Stim.Power.R,2)));
        TimeVectorS=nan(1,ceil(Stim.SamplingRate/Stim.FrameRate*size(Stim.Power.R,2)));
        
        for t = 1:size(Stim.Power.R,2)
            
            % Sum vectors by Cone Type.  Use as input for ConeVolution.
            LConeRStar = LConeRStarR.*Stim.Power.R(p,t) + LConeRStarG.*Stim.Power.G(p,t) + LConeRStarB.*Stim.Power.B(p,t);
            MConeRStar = MConeRStarR.*Stim.Power.R(p,t) + MConeRStarG.*Stim.Power.G(p,t) + MConeRStarB.*Stim.Power.B(p,t);
            SConeRStar = SConeRStarR.*Stim.Power.R(p,t) + SConeRStarG.*Stim.Power.G(p,t) + SConeRStarB.*Stim.Power.B(p,t);
            
            TimeVectorL(t*numel(Stim.SampPerFrRef)-(numel(Stim.SampPerFrRef)-1):t*numel(Stim.SampPerFrRef)) = LConeRStar;
            TimeVectorM(t*numel(Stim.SampPerFrRef)-(numel(Stim.SampPerFrRef)-1):t*numel(Stim.SampPerFrRef)) = MConeRStar;
            TimeVectorS(t*numel(Stim.SampPerFrRef)-(numel(Stim.SampPerFrRef)-1):t*numel(Stim.SampPerFrRef)) = SConeRStar;
        
        end
        if s == 1
            
            Stim.Bookend.L(p,:) = TimeVectorL(1:SampLBookend);
            Stim.Bookend.M(p,:) = TimeVectorM(1:SampLBookend);
            Stim.Bookend.S(p,:) = TimeVectorS(1:SampLBookend);
            
        end
        
        % Bookend each vector with sample points of background (for ConeVolution)
        PixelTVectL = [Stim.Bookend.L(p,:) TimeVectorL Stim.Bookend.L(p,:)];
        PixelTVectM = [Stim.Bookend.M(p,:) TimeVectorM Stim.Bookend.M(p,:)];
        PixelTVectS = [Stim.Bookend.S(p,:) TimeVectorS Stim.Bookend.S(p,:)];
        
        % Calculate Linear Responses for each Cone Type (modulated output should be less than a picoamp)
        LinRespLCone = ConeVolution(PixelTVectL,Stim.SamplingRate);
        LinRespMCone = ConeVolution(PixelTVectM,Stim.SamplingRate);
        LinRespSCone = ConeVolution(PixelTVectS,Stim.SamplingRate);
        
        % Dropping off bookends
        LinRespLCone = LinRespLCone(SampLBookend+1:end-SampLBookend);
        LinRespMCone = LinRespMCone(SampLBookend+1:end-SampLBookend);
        LinRespSCone = LinRespSCone(SampLBookend+1:end-SampLBookend);
        
        % Defining variables (only need defining once)
        if s == 1
            Stim.DataTime = size(LinRespLCone,2)/Stim.SamplingRate;
            scale = size(LinRespLCone,2) / 2;
            ncycles = ceil((stro.sum.exptParams.flash_length/1000)*unique(stro.trial(:,stro.idx.driftRate)));
            tsin = ncycles*2*pi/(2*scale):ncycles*2*pi/(2*scale):ncycles*2*pi;
        end
        if p == 1
            % Preallocate space for collapsing pixels over time
            SampSigs.LCones = nan(1,repetitions);
            SampSigs.MCones = nan(1,repetitions);
            SampSigs.SCones = nan(1,repetitions);
            SpaceCollapseL  = zeros(1,repetitions);
            SpaceCollapseM  = zeros(1,repetitions);
            SpaceCollapseS  = zeros(1,repetitions);

        end
        for z=1:repetitions
            
            w = w+1;
            waitbar(w/(repetitions*size(rgb,1)*size(Cones.Pixels,1)),h,'Looping Through Stimuli...')
            
            % Adding noise to each cone type by pixel, collapsing over time
            % L Cones
            Noise = ConeNoise(Stim.DataTime,Stim.SamplingRate);
            ncones = length(Cones.Pixels(p).Centers.LCones.X);
            LinRespLCone = LinRespLCone - LinRespLCone(1);
            SigPlusNoiseL = ncones*LinRespLCone + sqrt(ncones)*Noise;
            SigPlusNoiseL = abs(SigPlusNoiseL);
            dotsin = SigPlusNoiseL * sin(tsin)';
            dotcos = SigPlusNoiseL * cos(tsin)';
            scaledsin = dotsin/scale * sin(tsin);
            scaledcos = dotcos/scale * cos(tsin);
            freqphase = scaledsin + scaledcos;
            %TimeCollapseL = sum(LinRespLCone.*SigPlusNoiseL)*Stim.Pixel.Gaussian(p);
            TimeCollapseL = sum(freqphase.*SigPlusNoiseL);%*Stim.Pixel.Gaussian(p);
            %TimeCollapseL = sum(abs(SigPlusNoiseL))*Stim.Pixel.Gaussian(p);
            
            % M Cones
            Noise = ConeNoise(Stim.DataTime,Stim.SamplingRate);
            ncones = length(Cones.Pixels(p).Centers.MCones.X);
            LinRespMCone = LinRespMCone - LinRespMCone(1);
            SigPlusNoiseM = ncones*LinRespMCone + sqrt(ncones)*Noise;
            dotsin = SigPlusNoiseM * sin(tsin)';
            dotcos = SigPlusNoiseM * cos(tsin)';
            scaledsin = dotsin/scale * sin(tsin);
            scaledcos = dotcos/scale * cos(tsin);
            freqphase = scaledsin + scaledcos;
            %TimeCollapseM = sum(LinRespMCone.*SigPlusNoiseM)*Stim.Pixel.Gaussian(p);
            TimeCollapseM = sum(freqphase.*SigPlusNoiseM);%*Stim.Pixel.Gaussian(p);
            %TimeCollapseM = sum(abs(SigPlusNoiseM))*Stim.Pixel.Gaussian(p);
            
            % S Cones
            Noise = ConeNoise(Stim.DataTime,Stim.SamplingRate);
            ncones = length(Cones.Pixels(p).Centers.SCones.X);
            LinRespSCone = LinRespSCone - LinRespSCone(1);
            SigPlusNoiseS = ncones*LinRespSCone + sqrt(ncones)*Noise;
            dotsin = SigPlusNoiseS * sin(tsin)';
            dotcos = SigPlusNoiseS * cos(tsin)';
            scaledsin = dotsin/scale * sin(tsin);
            scaledcos = dotcos/scale * cos(tsin);
            freqphase = scaledsin + scaledcos;
            %TimeCollapseS = sum(LinRespSCone.*SigPlusNoiseS)*Stim.Pixel.Gaussian(p);
            TimeCollapseS = sum(freqphase.*SigPlusNoiseS);%*Stim.Pixel.Gaussian(p);
            %TimeCollapseS = sum(abs(SigPlusNoiseS))*Stim.Pixel.Gaussian(p);

            
            % Save sample responses
            SampSigs.LCones(z) = TimeCollapseL;
            SampSigs.MCones(z) = TimeCollapseM;
            SampSigs.SCones(z) = TimeCollapseS;
            
        end
        
        SpaceCollapseL = SpaceCollapseL + SampSigs.LCones;
        SpaceCollapseM = SpaceCollapseM + SampSigs.MCones;
        SpaceCollapseS = SpaceCollapseS + SampSigs.SCones;
        
    end
    
    %Collapse over pixels space, as well.
    %Final Product (3 numbers per sample. Plot in Cone Space)
    Responses.LCones(s,:) = SpaceCollapseL;
    Responses.MCones(s,:) = SpaceCollapseM;
    Responses.SCones(s,:) = SpaceCollapseS;
    
    clear SampSigs* Pixel* TimeVector* SpaceCollapse*
  
end

% Save Results
save(savefile, 'Responses')

% Plot Results
figure;hold on;grid on;
title('Response Clusters')
xlabel('L Cones')
ylabel('M Cones')
zlabel('S Cones')
for s=1:size(Responses.LCones,1)
    plot3(Responses.LCones(s,:),Responses.MCones(s,:),Responses.SCones(s,:),'o','MarkerFaceColor',unifrnd(0,1,1,3));
end


%% Compute ROC for each response cluster compared to the background

waitbar(1,h,'Performing Final Calculations...');
AUC = HorwitzRiekeROC(Responses);

ordercd1 = [1,3:2:size(AUC,1)];
ordercd2 = [1,2:2:size(AUC,1)];
AUCcd(1,:) = AUC(ordercd1);
AUCcd(2,:) = AUC(ordercd2);
rgbdir(:,:,1) = rgb(ordercd1,:);
rgbdir(:,:,2) = rgb(ordercd2,:);
contrast(1,:) = (rgbdir(:,:,1)-repmat(rgbdir(1,:,1),[8,1]))/rgbdir(1,:,1);
contrast(2,:) = (rgbdir(:,:,2)-repmat(rgbdir(1,:,2),[8,1]))/rgbdir(1,:,2);

% Display AUC
figure;hold on;grid on;
title('AUC')
plot(contrast(1,:),AUCcd(1,:),'-*')
plot(contrast(2,:),AUCcd(2,:),'g-*')
xlabel('Stimuli/Contrasts')
ylabel('AUC')


%% Attempting to fit Charlie's Psych Curve
% expt.norms = all stimuli of one color direction [1x8] (in lms units?)**No, units are strange.  Check with Charlie.
% monkFits.nTrials = # of trials at each contrast
% a = number of colors (one at a time, so a=1)
% b = number of spatial frequencies (b=1)
% nFigRows = Just what it sounds like? Use one for each color direction.
% nSptFreqs = same as b? =1 ** Practically, this is used for graphing columns, so renamed.
% monkFits = [2x8] first row: performance; second row: ntrials
global ModelVMonk 
ModelVMonk = figure; hold on;

[expt] = DTunpack_JPW(stro,1);

for s = 1:2
    
    norms(1,:) = expt.norms{1};
    norms(2,:) = expt.norms{2};
    nTrials = repetitions*ones(1,8);
    a = 1;
    b = 1;
    nFigRows = 2;
    nFigCols = 1;
    monkFits = [AUCcd(s,:);nTrials];
    
    monkFits = psychFun_JPW(norms(s,:),nTrials,a,b,nFigRows,nFigCols,monkFits,s);
    
end

% Display Model vs Monkey
figure(ModelVMonk);hold on;
subplot(2,1,1);hold on;
legend('Monkey','Model')
subplot(2,1,2);hold on;
legend('Monkey','Model')

% Close Waitbar
close(h)

% Taunt User
disp('Well...? Who won?')
