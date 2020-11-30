%% Initialize parameters
% code pretty much derived from t_colorGaborConeAbsorptionMovie.m which is
% writen by Dr. David Brainard
% Author - Abhishek De, CSHL project
close all; clear all;
ieInit;

ecc = 0;
ecc_deg = RetinalEccentricityMMToDegrees(ecc); % Retinal eccentricities in degrees
fov = 1.0;
% Define Gabor params 
gaborParams.fieldOfViewDegs = fov;
gaborParams.gaussianFWHMDegs = 0.942;
gaborParams.cyclesPerDegree = 1;
gaborParams.row = 128;
gaborParams.col = 128;
gaborParams.contrast = 1;
gaborParams.ang = 0;
gaborParams.ph = 0;
%gaborParams.coneContrasts = [0.06 -0.06 0.45]';% L-M+S, S-ON
gaborParams.coneContrasts = [0.06 -0.06 -0.45]';% L-M-S, S-OFF
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
[sampleTimes,TemporalWindow] = trapezoidalTemporalWindowCreate(temporalParams);

% define sample times 
nSampleTimes = length(sampleTimes);
gaborScene = cell(nSampleTimes,1);

% Setting the optical parameters
oiParams.fieldOfViewDegs = fov;
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

% bp params
bpParams.rectifyType = 1;
bpParams.cellType = 'offMidget';
bpParams.filterType = 1;
bpParams.ecc = 0;
bpParams.coneMosaic = theMosaic.pattern;

% ir params
irparams.name = 'macaque phys';
irparams.eyeSide = 'left';
irparams.eyeRadius = 0;
irparams.fov = fov;       % set at top
irparams.cellType = 'offparasolrpe';
irparams.averageMosaic = 1;
irparams.eyeAngle = 0;

% Defining the contrast levels you would want to test for 
reference_contrast = 0; % No gabor presentation onto the cone mosaic
contrast_lattice = [reference_contrast 0.9]; 
NCONTRASTS = numel(contrast_lattice);
coneIsomerizationRate = cell(1,NCONTRASTS);
coneCurrent = cell(1,NCONTRASTS);

disp('Intialized the parameters');
%% Computing the cone currents for each contrasts

tic;
gaborScene = cell(NCONTRASTS,nSampleTimes);
noiseparams.sampTime = theMosaic.sampleTime;
BP = cell(NCONTRASTS,1);
IR = cell(NCONTRASTS,1);
for jj = 1:NCONTRASTS % Repeat for a number of trials to incorporate the effect of noise
    % As of now I am showing the stimulus movie to the cone mosaic just once
    clear bp ir
    disp(jj);
    for ii = 1:nSampleTimes
        gaborParams.contrast = contrast_lattice(jj) * TemporalWindow(ii); % This is where u kind of define the peak contrast value within the gaussian temporal window
        gaborScene{jj,ii} = colorGaborSceneCreate(gaborParams);
        theOI{jj,ii} = oiCompute(theBaseOI,gaborScene{jj,ii});
        gaborConeAbsorptions(:,:,ii) = theMosaic.compute(theOI{jj,ii},'currentFlag',false); % captures the cone photoisomerizations
    end
    coneIsomerizationRate{jj} = gaborConeAbsorptions/theMosaic.integrationTime; % calculating the cone isomerization rate
    coneCurrent{jj} = theMosaic.os.compute(coneIsomerizationRate{jj},theMosaic.pattern);
%     coneCurrent{jj} = osAddNoise(coneCurrent{jj},noiseparams); % Add Fred's noise
    osSet(theMosaic.os,'coneCurrentSignal',coneCurrent{jj})
   % bipolar computation
    bp = bipolar(theMosaic.os, bpParams);
    bp.set('sRFcenter',1);
    bp.set('sRFsurround',1);
    BP{jj} = bipolarCompute(bp, theMosaic.os);
    % RGC computation
    irparams.inputScale = size(bp.sRFcenter,1);
    irparams.inputSize = size(bp.responseCenter);
    ir = irPhys(bp, irparams);
    IR{jj}  = irCompute(ir, bp); % Storing the RGC responses
end

toc;
disp('Obtained the RGC responses for different contrast values');

%% Debugging stuff
curr_diff = coneCurrent{end} - coneCurrent{1};
curr_diff = RGB2XWFormat(curr_diff);
figure,plot(curr_diff');

bp_cur = BP{end}.rectificationCenter - BP{1}.responseCenter;
bp_cur = RGB2XWFormat(bp_cur);
figure,plot(bp_cur');

%% Visualize the gaborScene 
disp('Displaying the stimulus');
for jj = NCONTRASTS
    for ii = 1:nSampleTimes
        figure(1); sceneShowImage(gaborScene{jj,ii}); drawnow; pause(0.1);
    end
    pause(0.1); disp('Presenting the next stimulus');
end

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
[final_model, final_fval, final_success] = weibullFit(contrast_lattice(2:end), neurometric_curve);
model_fit = (1 - 0.5.*exp(-((contrast_lattice(2:end)./final_model(1)).^final_model(2))));
% Plot the neurometric curve and fit a weibull distribution
figure; semilogx(contrast_lattice(2:end), 100*neurometric_curve,'-bo'); xlabel('Contrast levels'); ylabel('Performance in percentage'); title('Neurometric curve'); hold on;
semilogx(contrast_lattice(2:end), 100*model_fit, '-ro'); legend('data','fit');hold off;
threshold = final_model(1);
slope = final_model(2);

