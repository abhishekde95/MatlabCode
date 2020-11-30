%% I am writing this script to see if the population cone responses have any bias towards the intermediate directions stimuli
% am using for my Psychophysics experiment
% A big chunk of this code was borrowed from my CSHL project
% Author - Abhishek De
clearvars; ieInit;

% Initialize basic parameters
ecc_inDeg = 5; % in visual degrees
tmp = RetinalEccentricityMMToDegrees(1); % Retinal eccentricities in degrees
ecc_inMM = ecc_inDeg/tmp;
fov = 1.0; % Field of View in degrees, keeping to 1x1 degrees

% getting the gabor parameters, see spatialParamsGenerate.m
spatialParams.type = 'Spatial';
spatialParams.spatialType = 'Gabor';
spatialParams.windowType = 'Gaussian';
spatialParams.fieldOfViewDegs = fov;
spatialParams.gaussianFWHMDegs = 0.942; % 2.3548 x standard deviation = 0.4
spatialParams.cyclesPerDegree = 1.0;
spatialParams.row = 128;
spatialParams.col = 128;
spatialParams.ang = 0;
spatialParams.ph = 0;
spatialParams.viewingDistance = 0.7; % in metres, check from t_colorGaborConeAbsorptionMovie.m

% getting the color modulation parameters, see colorModulationParamsGenerate.m
colorModulationParams.type = 'ColorModulation';
colorModulationParams.modulationType = 'monitor';
colorModulationParams.contrast = 1;
colorModulationParams.coneContrasts = [1 1 1]';
colorModulationParams.coneContrasts = colorModulationParams.coneContrasts./norm(colorModulationParams.coneContrasts); % normlizing the contrast
colorModulationParams.startWl = 380;
colorModulationParams.endWl = 780;
colorModulationParams.deltaWl = 4;

% getting the background parameters, see backgroundParamsGenerate.m
backgroundParams.type = 'Background';
backgroundParams.backgroundType = 'monitor';
backgroundParams.backgroundxyY = [0.27 0.30 9.8]';
backgroundParams.monitorFile = 'CRT-MODEL';
backgroundParams.leakageLum = 1.0;
backgroundParams.lumFactor = 5;

% Obtaining the temporal parameters, see temporalParamsGenerate.m
temporalParams.type = 'Temporal';
temporalParams.frameRate = 60;
temporalParams.windowTauInSeconds = 0.1;
temporalParams.stimulusDurationInSeconds = 0.5; % 1000 ms stimulus presentation
temporalParams.stimulusSamplingIntervalInSeconds = 1/temporalParams.frameRate;
temporalParams.eyesDoNotMove = false;
temporalParams.secondsToInclude = 1.0;
temporalParams.secondsToIncludeOffset = 0.01;
temporalParams.secondsForResponseStabilization = 0.00;
temporalParams.secondsForResponseExtinction = 0.00;
temporalParams.emPathType = 'none';
[temporalParams.sampleTimes,temporalParams.TemporalWindow] = trapezoidalTemporalWindowCreate(temporalParams);
temporalParams.nSampleTimes = length(temporalParams.sampleTimes);

% Obtaining optical image params, see oiParamsGenerate.m
oiParams.type = 'Optics';
oiParams.offAxis = false;
oiParams.blur = false; % no blurring
oiParams.lens = false; % no effect of lens transmittance
oiParams.pupilDiamMm = 3;
oiParams.opticsModel = 'WvfHuman';
oiParams.fieldOfViewDegs = spatialParams.fieldOfViewDegs;
BaseOI = colorDetectOpticalImageConstruct(oiParams);

% Obtaining the cone mosaic params, see mosaicParamsGenerate.m
mosaicParams.type = 'Mosaic';
mosaicParams.conePacking = 'rect';
mosaicParams.realisticSconeSubmosaic = false;       % if this is set to true, there will be a 0.3 deg S-cone free region and the S-cone lattice will be semiregular
mosaicParams.LMSRatio = [0.62 0.31 0.07];
mosaicParams.innerSegmentSizeMicrons = 1.4;
mosaicParams.apertureBlur = false;
mosaicParams.coneSpacingMicrons = 2.0;
mosaicParams.mosaicRotationDegs = 0;
mosaicParams.macular = false;
mosaicParams.eccentricityDegs = ecc_inDeg;
mosaicParams.integrationTimeInSeconds = 5/1000;
mosaicParams.osTimeStepInSeconds = 0.1/1000;
cm = coneMosaic('os',osLinear); % creating the cone mosaic, currently on hold
cm.setSizeToFOV(fov);
cm.integrationTime = mosaicParams.integrationTimeInSeconds;
cm.noiseFlag = 'none';
cm.os.noiseFlag = 'none';


theBaseColorModulationParams = colorModulationParams;
theBaseColorModulationParams.coneContrasts = [0 0 0]';
theBaseColorModulationParams.contrast = 0;
backgroundScene = colorSceneCreate(spatialParams, backgroundParams, theBaseColorModulationParams, []);
BaseOI = colorDetectOpticalImageConstruct(oiParams);
oiBackground = BaseOI;
oiBackground = oiCompute(oiBackground, backgroundScene);
oiModulated = oiBackground;
zeroContrastOIsequence = oiSequence(oiBackground, oiModulated, temporalParams.sampleTimes, temporalParams.TemporalWindow, 'composition', 'blend');
eyeMovementsNum = zeroContrastOIsequence.maxEyeMovementsNumGivenIntegrationTime(cm.integrationTime);
theZeroContrastEMpaths = colorDetectMultiTrialEMPathGenerate(cm, 1, eyeMovementsNum, 'none');
[isomerizationsAdaptingStim, photocurrentsAdaptingStim, osImpulseResponseFunctions, osMeanCurrents] = ...
    cm.computeForOISequence(zeroContrastOIsequence, ...
    'emPaths', theZeroContrastEMpaths, ...
    'interpFilters', [], ...
    'meanCur', [], ...
    'currentFlag', true);

%% Creating a battery of Gabor, the Scene and the OI are working correctly
NTRIALS = 11;
contrast_lattice = [0 0.1 0.5];
NCONTRASTS = numel(contrast_lattice);
isomerizations = cell(NCONTRASTS,NTRIALS);
photocurrents = cell(NCONTRASTS,NTRIALS);
for ii = 1: NCONTRASTS
    % Create the modulated scene (varied contrast)
    colorModulationParams.contrast = NCONTRASTS(ii);
    modulatedScene = colorSceneCreate(spatialParams, backgroundParams, colorModulationParams, []);
    oiModulated = BaseOI;
    oiModulated = oiCompute(oiModulated, modulatedScene);
    stimulusOIsequence = oiSequence(oiBackground, oiModulated, temporalParams.sampleTimes, temporalParams.TemporalWindow, 'composition', 'blend');
    % Generate eye movement paths for all response instances
    eyeMovementsNum = stimulusOIsequence.maxEyeMovementsNumGivenIntegrationTime(cm.integrationTime);
    theEMpaths = colorDetectMultiTrialEMPathGenerate(cm, 1, eyeMovementsNum, temporalParams.emPathType);
    
    for jj = 1:NTRIALS
        % Compute the noisy responses
        
        cm.noiseFlag = 'random';
        cm.os.noiseFlag = 'random';
        
        [isomerizations{ii,jj}, photocurrents{ii,jj}, ~, ~] = ...
            cm.computeForOISequence(stimulusOIsequence, ...
            'emPaths', theEMpaths, ...
            'interpFilters', osImpulseResponseFunctions, ...
            'meanCur', osMeanCurrents, ...
            'currentFlag', true);
    end
end

matchedfilter = photocurrents{end,1} - photocurrents{1,1};
matchedfilter = RGB2XWFormat(matchedfilter);













































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































