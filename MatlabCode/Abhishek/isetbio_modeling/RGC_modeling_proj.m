%% Initialize parameters
% code pretty much derived from t_colorGaborConeAbsorptionMovie.m which is
% writen by Dr. David Brainard 
% Some of the previous code was written during the CSHL project
% Author - Abhishek De, CSHL project
close all; clear all;
ieInit; 

ecc = 0; % For the time being, I am just looking at the foveal data
fov = 1.0; % in degrees of visual angle

% Define Gabor params 
gaborParams.fieldOfViewDegs = fov;
gaborParams.gaussianFWHMDegs = 0.942;
gaborParams.cyclesPerDegree = 3;
gaborParams.row = 128;
gaborParams.col = 128;
gaborParams.contrast = 1;
gaborParams.ang = 0;
gaborParams.ph = 0;
gaborParams.coneContrasts = [0.06 -0.06 0.45]';% L-M+S, Lime-Magenta
gaborParams.coneContrasts = gaborParams.coneContrasts/norm(gaborParams.coneContrasts); % normalizing the cone contrast
gaborParams.backgroundxyY = [0.27 0.30 49.8]';
gaborParams.monitorFile = 'CRT-MODEL';
gaborParams.viewingDistance = 1.0;

% Define temporal parameters
frameRate = 60;
temporalParams.windowTauInSeconds = 0.1;
temporalParams.stimulusDurationInSeconds = 0.5; % 500 ms stimulus presentation
temporalParams.stimulusSamplingIntervalInSeconds = 1/frameRate;
temporalParams.eyesDoNotMove = true; 
temporalParams.millisecondsToInclude = 500;
temporalParams.millisecondsToIncludeOffset = 35;
[sampleTimes,gaussianTemporalWindow] = gaussianTemporalWindowCreate(temporalParams);

% define sample times 
nSampleTimes = length(sampleTimes);

% Setting the optical parameters
oiParams.fieldOfViewDegs = gaborParams.fieldOfViewDegs;
oiParams.offAxis = false;
oiParams.blur = false; % no blurring due to optics
oiParams.lens = false; % optical lens inactive

% Creating the optics structure
oi = colorDetectOpticalImageConstruct(oiParams);

% Creating sensor 
absorptions = sensorCreate('human'); % foveal cone mosaic
absorptions = sensorSetSizeToFOV(absorptions, gaborParams.fieldOfViewDegs);
absorptions.noiseFlag = 0; % Don't add photoisomerization noise
cone_mosaic = sensorGet(absorptions,'cone type');

% Defining the contrast levels you would want to test for 
reference_contrast = 0; % No gabor presentation onto the cone mosaic
contrast_lattice = [0.5];
NCONTRASTS = numel(contrast_lattice);

disp('Intialized the parameters');
%% Computing the cone currents for each contrasts

tic;
gaborConeAbsorptions = zeros([sensorGet(absorptions, 'size') nSampleTimes]);
for jj = 1:NCONTRASTS % Repeat for a number of trials to incorporate the effect of noise
    disp(jj);
    for ii = 1:nSampleTimes
        gaborParams.contrast = contrast_lattice(jj) * gaussianTemporalWindow(ii); % This is where u kind of define the peak contrast value within the gaussian temporal window
        gaborScene = colorGaborSceneCreate(gaborParams);
        oi = oiCompute(oi,gaborScene);
        absorptions = sensorCompute(absorptions,oi );
        gaborConeAbsorptions(:,:,ii) = sensorGet(absorptions, 'volts');
    end
    absorptions = sensorSet(absorptions, 'volts', gaborConeAbsorptions);
    sensor = absorptions;
end

% Debugging stuff
% k = RGB2XWFormat(gaborConeAbsorptions);
% figure,plot(k');
% vcNewGraphWin;
% ieMovie(sensor.data.volts);

toc;
disp('Obtained the cone currents for different contrast values');

%% Outer segment calculation - biophysical model
% The next step is to get the outer segment current. 
osB = osCreate('biophys');
patchSize = sensorGet(sensor,'width','m');
osB = osSet(osB, 'patch size', patchSize);
timeStep = sensorGet(sensor,'time interval','sec');
osB = osSet(osB, 'time step', timeStep);
sensorVolts = sensorGet(sensor,'volts');
paramsOS.bgVolts = 10*mean(sensorVolts(:));
paramsOS.noiseFlag = 1; % Adds noise - Fred's model
clear sensorVolts
osB = osCompute(osB,sensor,paramsOS);

% Debugging stuff
% k = RGB2XWFormat(osB.coneCurrentSignal);
% figure,plot(k');

%% Find bipolar responses
clear bp
bpParams.rectifyType = 1;
bpParams.cellType = 'offMidget';
bpParams.filterType = 1; 
bpParams.ecc = ecc;
bpParams.coneMosaic = cone_mosaic;
bp = bipolar(osB, bpParams);
bp.set('sRFcenter',1);
bp.set('sRFsurround',1);
bp = bipolarCompute(bp, osB);

% Visualize the bipolar mosaic response
% bipolarPlot(bp,'response');
% vcNewGraphWin;
% ieMovie(bp.responseCenter);
%% Set other RGC mosaic parameters

clear irparams innerRetinaSU
irparams.name = 'macaque phys';
irparams.eyeSide = 'left'; 
irparams.eyeRadius = ecc; 
irparams.fov = fov;       % set at top
irparams.cellType = 'offparasolrpe';
irparams.averageMosaic = 1;
irparams.eyeAngle = 0;
irparams.inputScale = size(bp.sRFcenter,1);
irparams.inputSize = size(bp.responseCenter);
ir = irPhys(bp, irparams); 
%%
ir = irCompute(ir, bp); 
mosaicPlot(ir, bp, sensor, irparams, irparams.cellType, ecc)


