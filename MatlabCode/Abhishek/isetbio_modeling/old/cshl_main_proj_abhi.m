%% Initialize parameters
% code pretty much derived from t_colorGaborConeAbsorptionMovie.m which is
% writen by Dr. David Brainard
% Author - Abhishek De, CSHL project
close all; clear all;
ieInit;

% Define Gabor params 
gaborParams.fieldOfViewDegs = 1.0;
gaborParams.gaussianFWHMDegs = 0.942;
gaborParams.cyclesPerDegree = 1;
gaborParams.row = 128;
gaborParams.col = 128;
gaborParams.contrast = 1;
gaborParams.ang = 0;
gaborParams.ph = 0;
%gaborParams.coneContrasts = [0.06 -0.06 0.45]';% L-M+S, S-ON
gaborParams.coneContrasts = [0.06 -0.06 -0.45]';% L-M-S, S-OFF
%gaborParams.coneContrasts = [0.06 0.06 0.45]';% Luminance
%gaborParams.coneContrasts = [0.06 -0.06 0.0]';% L-M
%gaborParams.coneContrasts = [0.0 0.0 1.0]';% S
%gaborParams.coneContrasts = [0.0 0.0 -1.0]';% -S
gaborParams.coneContrasts = gaborParams.coneContrasts/norm(gaborParams.coneContrasts); % normalizing the cone contrast
gaborParams.backgroundxyY = [0.27 0.30 49.8]';
gaborParams.monitorFile = 'CRT-MODEL';
gaborParams.viewingDistance = 1.0;

% Define temporal parameters
frameRate = 60;
temporalParams.windowTauInSeconds = 0.1;
temporalParams.stimulusDurationInSeconds = 0.5; % 500 ms stimulus presentation
temporalParams.stimulusSamplingIntervalInSeconds = 1/frameRate;
temporalParams.eyesDoNotMove = false; 
temporalParams.millisecondsToInclude = 500;
temporalParams.millisecondsToIncludeOffset = 35;
% [sampleTimes,gaussianTemporalWindow] = gaussianTemporalWindowCreate(temporalParams);
[sampleTimes,TemporalWindow] = trapezoidalTemporalWindowCreate(temporalParams);

% define sample times 
nSampleTimes = length(sampleTimes);
gaborScene = cell(nSampleTimes,1);

% Setting the optical parameters
oiParams.fieldOfViewDegs = gaborParams.fieldOfViewDegs;
oiParams.offAxis = false;
oiParams.blur = false; % no blurring due to optics
oiParams.lens = false; % optical lens inactive

% Creating the optics structure
theBaseOI = colorDetectOpticalImageConstruct(oiParams);
theOI = cell(nSampleTimes,1);

% Defining the cone mosaic parameters
mosaicParams.fieldOfViewDegs = gaborParams.fieldOfViewDegs;
mosaicParams.macular = false;
mosaicParams.LMSRatio = [0.6 0.3 0.1];
mosaicParams.timeStepInSeconds = temporalParams.stimulusSamplingIntervalInSeconds;
mosaicParams.integrationTimeInSeconds = mosaicParams.timeStepInSeconds;
mosaicParams.photonNoise = false; % Not adding photon noise
mosaicParams.osNoise = false;% Adding noise, Fred's model
mosaicParams.osModel = 'Linear';
mosaicParams.append = false;
% Defining the cone mosaic
theMosaic = colorDetectConeMosaicConstruct(mosaicParams);
mosaicParams.pattern = theMosaic.pattern;


% Defining the contrast levels you would want to test for 
reference_contrast = 0; % No gabor presentation onto the cone mosaic
contrast_lattice = [reference_contrast 0.005 0.01 0.02 0.04 0.08 0.9]; 
save('contrast_lattice.mat','contrast_lattice');
NCONTRASTS = numel(contrast_lattice);
coneIsomerizationRate = cell(1,NCONTRASTS);
coneCurrent = cell(1,NCONTRASTS);

disp('Intialized the parameters');
%% Computing the cone currents for each contrasts

tic;
gaborScene = cell(NCONTRASTS,nSampleTimes);
for jj = 1:NCONTRASTS % Repeat for a number of trials to incorporate the effect of noise
    % As of now I am showing the stimulus movie to the cone mosaic just once
    disp(jj);
    for ii = 1:nSampleTimes
        gaborParams.contrast = contrast_lattice(jj) * TemporalWindow(ii); % This is where u kind of define the peak contrast value within the gaussian temporal window
        gaborScene{jj,ii} = colorGaborSceneCreate(gaborParams);
        theOI{jj,ii} = oiCompute(theBaseOI,gaborScene{jj,ii});
        gaborConeAbsorptions(:,:,ii) = theMosaic.compute(theOI{jj,ii},'currentFlag',false); % captures the cone photoisomerizations
%         current(:,:,ii) = theMosaic.os.coneCurrentSignal;
    end
    coneIsomerizationRate{jj} = gaborConeAbsorptions/theMosaic.integrationTime; % calculating the cone isomerization rate
    coneCurrent{jj} = theMosaic.os.compute(coneIsomerizationRate{jj},theMosaic.pattern);
%     coneCurrent{jj} = current;
end
matchedfilter = coneCurrent{end} - coneCurrent{1};
matchedfilter = RGB2XWFormat(matchedfilter);
toc;
disp('Obtained the cone currents for different contrast values');

%% Visualize the matched filter
k = XW2RGBFormat(matchedfilter,theMosaic.mosaicSize(1),theMosaic.mosaicSize(2));
for ii = 1:size(k,3)
    figure(2), imagesc(k(:,:,ii)); drawnow; pause(0.05);
end


%% Visualize the optical image
disp('Displaying the image falling on the retina');
for jj = 1: NCONTRASTS
    for ii = 1:nSampleTimes
        figure(2); oiShowImage(theOI{jj,ii}); drawnow;
    end
    pause(0.1); disp('Presenting the next optical image');
end
    
%% Plot the cone mosaic and Visualize the cone photocurrents
disp('Displaying the cone currents and the retinal mosaic');
%load coneCurrent_Siso.mat
%load coneCurrent_negSiso.mat

for jj = NCONTRASTS
    for ii = 1:nSampleTimes
        figure(3), imagesc(coneCurrent{jj}(:,:,ii)); drawnow;
    end
end
%S_plus = RGB2XWFormat(coneCurrent_Siso{end});
%S_minus = RGB2XWFormat(coneCurrent_negSiso{end});
%figure, subplot(121),plot(S_plus'); subplot(122),plot(S_minus');
theMosaic.plot('conemosaic','hf',figure(4));

%% Visualize the cone photoisomerizations
disp('Displaying the cone isomerizations');
for jj = 1: NCONTRASTS
    for ii = 1:nSampleTimes
        figure(5), imagesc(coneIsomerizationRate{jj}(:,:,ii)); drawnow;
    end
end
%% Next step is to generate pooled responses over noisy data
NTRIALS = 1000; % Number of trials 
sensor = theMosaic;
% Find coordinates of L, M and S cones.
cone_mosaic = theMosaic.pattern;
% The next step is to convolve the 1D filters with the 1D current data at each point in the cone mosaic.
[sz1, sz2, nSteps] = size(coneIsomerizationRate{1});
params = [];
pooledData = zeros(NCONTRASTS,NTRIALS,3); % for 3 different cone types
noiseparams.sampTime = theMosaic.sampleTime;  % Sec
theMosaic.os.noiseFlag = false;
tic;
for ii = 1:NCONTRASTS
    disp(ii);
    for iter = 1:NTRIALS 
        coneCurrentRS = osAddNoise(coneCurrent{ii},noiseparams);
        coneCurrentRS = RGB2XWFormat(coneCurrentRS-coneCurrent{1});
        coneCurrentRS = coneCurrentRS.*matchedfilter; % multiply conecurrents by the matched filter
        for cone_type = 2:4
            cone_locations = find(cone_mosaic==cone_type);
            coneCurrentSingleType = (coneCurrentRS(cone_locations,:));
            pooledData(ii,iter, cone_type-1) = mean(rms(coneCurrentSingleType,2)); % As of now don't know which one is better
            %pooledData(ii,iter, cone_type-1) = mean(mean(abs(coneCurrentSingleType),2));
        end
    end
end
toc;
disp('Pooled all the Data');

%% Visualize the 3-D data
color = ['r','g'];
figure(3);
for ii = 2: NCONTRASTS
    noisy_data = squeeze(pooledData(1,:,:));
    data = squeeze(pooledData(ii,:,:));
    figure(3), subplot(2,3,ii-1); scatter3(noisy_data(:,1),noisy_data(:,2),noisy_data(:,3),color(1)); hold on; 
    scatter3(data(:,1),data(:,2),data(:,3),color(2)); 
    xlabel('L'); ylabel('M'); zlabel('S'); title(strcat('Contrast=',num2str(contrast_lattice(ii)))); 
    set(gca,'Fontsize',13); hold off; drawnow;
end

%% Perform LDA on the data and find the Linear Discriminant vector that best classifies the data
Target = zeros(2*NTRIALS,1);
Target(NTRIALS+1:2*NTRIALS) = deal(1);
factor = 1e8;

neurometric_curve = zeros(1,NCONTRASTS-1);
figure(4);
for ii = 2:NCONTRASTS 
    LDA_vec = [];
    LDA_vec = cov(squeeze(pooledData(1,:,:))) * (mean(squeeze(pooledData(1,:,:))) - mean(squeeze(pooledData(ii,:,:))))'; % find the LDA vector
    LDA_vec = LDA_vec * factor;
    LDA_vec = LDA_vec/norm(LDA_vec);  
    % Project the noise and signal onto the Linear Discriminant vector
    noise = squeeze(pooledData(1,:,:))*LDA_vec;
    signal = squeeze(pooledData(ii,:,:))*LDA_vec;
    [counts,centers] = hist([noise signal],20);
    figure(4),subplot(2,3,ii-1),hist([noise signal],20); drawnow;
    xlabel('Projection onto LDA vector'); ylabel('frequency');title(strcat('Contrast=',num2str(contrast_lattice(ii)))); 
    set(gca,'Fontsize',13);hold off; drawnow;
    % ROC analysis and neurometric curve - continue
    % Ratio of false positives to true positives
    neurometric_curve(ii-1) = 1-(sum(min(counts,[],2))/(2*NTRIALS));
end

%%
[final_model, final_fval, final_success] = weibullFit_iset(contrast_lattice(2:end), neurometric_curve);
model_fit = (1 - 0.5.*exp(-((contrast_lattice(2:end)./final_model(1)).^final_model(2))));
% Plot the neurometric curve and fit a weibull distribution
figure; semilogx(contrast_lattice(2:end), 100*neurometric_curve,'-bo'); xlabel('Contrast levels'); ylabel('Performance in percentage'); title('Neurometric curve'); hold on;
semilogx(contrast_lattice(2:end), 100*model_fit, '-ro'); legend('data','fit');hold off;
threshold = final_model(1);
slope = final_model(2);

%% Compute bipolar cell responses from the cone outsegment currents

% Initialize the bipolar object

% Compute the bipolar mosaic response


%% Compare neurometric of lime-magenta and orange-cyan stimuli; The neurometric curve should be indistinguishable at the level of cones
load neurometric_SON2.mat;
load neurometric_SOFF2.mat;
load neurometric_red_green2.mat
load neurometric_green_red2.mat
load neurometric_Siso2.mat;
load contrast_lattice.mat
figure(5); semilogx(contrast_lattice(2:end),neurometric_SON2,'-bo'); ylabel('Performance'); title(' Neurometric Curve'); hold on;
semilogx(contrast_lattice(2:end),neurometric_SOFF2,'-ro'); hold on
semilogx(contrast_lattice(2:end),neurometric_Siso2,'-go'); hold on;
semilogx(contrast_lattice(2:end),neurometric_red_green2,'-mo'); hold on;
semilogx(contrast_lattice(2:end),neurometric_green_red2,'-ko'); hold on;
legend('lime-magenta','orange-cyan','Siso','L-M','M-L'); hold off;

%% Compare neurometric of brainard and anti-brainard approach 
load neurometric_brainard2.mat;
load neurometric_antibrainard2.mat;
figure(5); semilogx(contrast_lattice(2:end),neurometric_brainard2,'-bo', 'Linewidth',2); ylabel('Performance'); title(' Neurometric Curve'); hold on;
semilogx(contrast_lattice(2:end),neurometric_antibrainard2,'-ro','Linewidth',2); legend('In built photocurrent','Current calculated from cone photoisomerization'); hold off;

%% Compare neurometric of lime-magenta and orange-cyan stimuli; The neurometric curve should be indistinguishable at the level of cones
load neurometric_lm2.mat;
load neurometric_oc2.mat;
%load neurometric_lum.mat
load neurometric_RG.mat;
%load neurometric_GR.mat;
load neurometric_Siso.mat
%load neurometric_negSiso.mat

figure(5); semilogx(contrast_lattice(2:end),neurometric_lm2,'-bo','Linewidth',2); ylabel('Performance'); xlabel('contrast');title(' Neurometric Curve'); hold on;
semilogx(contrast_lattice(2:end),neurometric_oc2,'-ro','Linewidth',2); hold on
%semilogx(contrast_lattice(2:end),neurometric_lum,'-ko','Linewidth',2); hold on
semilogx(contrast_lattice(2:end),neurometric_RG,'-mo','Linewidth',2); hold on
%semilogx(contrast_lattice(2:end),neurometric_GR,'--m','Linewidth',2); hold on
semilogx(contrast_lattice(2:end),neurometric_Siso,'-go','Linewidth',2); hold on
%semilogx(contrast_lattice(2:end),neurometric_negSiso,'--g','Linewidth',2); hold on
legend('L-M+S (lime magenta)','L-M-S (orange cyan)','L-M','S'); 
set(gca,'Fontsize',20);hold off;

%% Compare the responses of lime magenta stimulus on different cone ratios
load neurometric_lm_coneratio1.mat;
load neurometric_lm_coneratio2.mat;

figure(5); semilogx(contrast_lattice(2:end),neurometric_lm_coneratio1,'-bo','Linewidth',2); ylabel('Performance'); title(' Neurometric Curve'); xlabel('contrast');hold on;
semilogx(contrast_lattice(2:end),neurometric_lm_coneratio2,'-ro','Linewidth',2); hold on
legend('L:M:S = 0.45:0.45:0.1','L:M:S = 0.6:0.3:0.1'); title('Effect of cone ratio'); set(gca,'Fontsize',20); hold off;

%% Compare the responses of lime magenta stimulus for different spatial phases
load neurometric_lm_coneratio2.mat;
load neurometric_lm_coneratio2_phase45.mat;
load neurometric_lm_coneratio2_phase90.mat;
load neurometric_lm_coneratio2_phase180.mat;

figure(5); semilogx(contrast_lattice(2:end),neurometric_lm_coneratio2,'-bo','Linewidth',2); ylabel('Performance'); title(' Neurometric Curve'); xlabel('contrast');hold on;
semilogx(contrast_lattice(2:end),neurometric_lm_coneratio2_phase45,'-ro','Linewidth',2); hold on
semilogx(contrast_lattice(2:end),neurometric_lm_coneratio2_phase90,'-ko','Linewidth',2); hold on
semilogx(contrast_lattice(2:end),neurometric_lm_coneratio2_phase180,'-go','Linewidth',2); hold on
legend('0','45','90','180'); title('Effect of phase'); set(gca,'Fontsize',20);hold off;

% Don't see any effect of the phase of the stimuli
%% Checking the same thing for L-M 180 and M-L
load neurometric_RG.mat
load neurometric_RG180.mat
load neurometric_GR.mat

figure(5); semilogx(contrast_lattice(2:end),neurometric_RG180,'-bo','Linewidth',2); ylabel('Performance'); title(' Neurometric Curve'); hold on;
semilogx(contrast_lattice(2:end),neurometric_GR,'-ro','Linewidth',2); hold on
semilogx(contrast_lattice(2:end),neurometric_RG,'-go','Linewidth',2); hold on
legend('L-M phase 180','L-M','M-L'); title('L-M: Effect of phase'); hold off;


% The code is working till this point
