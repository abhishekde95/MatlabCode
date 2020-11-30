%% Initialize parameters
% Author - Fred Rieke & Abhishek De
clearvars; ieInit;

% Initialize basic parameters
ecc_inDeg = 5; % in visual degrees
tmp = RetinalEccentricityMMToDegrees(1); % Retinal eccentricities in degrees
ecc_inMM = ecc_inDeg/tmp;
fov = 1.2; % Field of View in degrees, from Greg and Fred's paper 

% Obtaining Gabor, color Modulation and backgroundparams 
[gaborParams, colorModulationParams, backgroundParams] = getSceneParams_abhi(fov);

% Obtaining temporal parameters
temporalParams = getTemporalParams_abhi();
[sampleTimes,TemporalWindow] = trapezoidalTemporalWindowCreate(temporalParams); 
nSampleTimes = length(sampleTimes);

% Setting the optical parameters
oiParams = getoiParams_abhi(fov);
theBaseOI = colorDetectOpticalImageConstruct(oiParams);

% Defining the cone mosaic parameters
mosaicParams = getMosaicParams_abhi( ecc_inMM, fov,temporalParams );
theMosaic = createConeMosaic_abhi(mosaicParams,'coneMosaicHex');
% theMosaic.visualizeGrid(); % Visualize the cone maps

% bp params
bpParams.rectifyType = 1;
bpParams.cellType = 'offMidget';
bpParams.filterType = 1;
bpParams.ecc = ecc_inMM;
bpParams.coneMosaic = theMosaic.pattern;

% ir params
irparams.name = 'macaque phys';
irparams.eyeSide = 'left';
irparams.eyeRadius = 0;
irparams.fov = fov;       % set at top
irparams.cellType = 'offmidget';
irparams.averageMosaic = 1;
irparams.eyeAngle = 0;

% Defining the contrast levels you would want to test for 
reference_contrast = 0; % No gabor presentation onto the cone mosaic
contrast_lattice = [reference_contrast 0.5];%logspace(-4,-1,12)]; 
NCONTRASTS = numel(contrast_lattice);
gaborScene = cell(NCONTRASTS, nSampleTimes);
theOI = cell(NCONTRASTS, nSampleTimes);
coneIsomerizationRate = cell(1,NCONTRASTS);
coneCurrent = cell(1,NCONTRASTS);
BP = cell(NCONTRASTS,1);
IR = cell(NCONTRASTS,1);

disp('Intialized the parameters');
%%
tic;
noiseparams.sampTime = theMosaic.integrationTime;
for jj = 1:NCONTRASTS % Repeat for a number of trials to incorporate the effect of noise
    % As of now I am showing the stimulus movie to the cone mosaic just on
    clear bp 
    disp(jj);
    for ii = 1:nSampleTimes
        colorModulationParams.contrast = contrast_lattice(jj) * TemporalWindow(ii); % This is where u kind of define the peak contrast value within the gaussian temporal window
        gaborScene{jj,ii} = colorGaborSceneCreate(gaborParams,backgroundParams,colorModulationParams);
        theOI{jj,ii} = oiCompute(theBaseOI,gaborScene{jj,ii});
        gaborConeAbsorptions(:,:,ii) = theMosaic.compute(theOI{jj,ii},'currentFlag',false); % captures the cone photoisomerizations
    end
    coneIsomerizationRate{jj} = gaborConeAbsorptions/theMosaic.integrationTime; % calculating the cone isomerization rate
    tmp_curr = theMosaic.os.compute(cat(3,coneIsomerizationRate{1},coneIsomerizationRate{jj}),theMosaic.pattern);
    coneCurrent{jj} = tmp_curr(:,:,nSampleTimes+1:end);
    tempCurrent = tmp_curr(:,:,nSampleTimes+1:end);
%     coneCurrent{jj} = osAddNoise(coneCurrent{jj},noiseparams); % Add Fred's noise
%     osSet(theMosaic.os,'coneCurrentSignal',coneCurrent{jj})
%     % bipolar computation
%     bp = bipolar(theMosaic.os, bpParams);
%     bp.set('sRFcenter',1);
%     bp.set('sRFsurround',1);
%     bp = bipolarCompute(bp, theMosaic.os);
%     % RGC computation
%     irparams.inputScale = size(bp.sRFcenter,1);
%     irparams.inputSize = size(bp.responseCenter);
%     innerRetinaSU = ir(bp, irparams);
%     innerRetinaSU.mosaicCreate('type',cellType,'model','GLM');
   
end

iso_matchedfilter = RGB2XWFormat(coneIsomerizationRate{end}-coneIsomerizationRate{1});
cur_matchedfilter = RGB2XWFormat(coneCurrent{end}-coneCurrent{1});
toc;
disp('Obtained the RGC responses for different contrast values');

%% For debugging purposes
% Background cone excitations in Hass et al; 2015: L=7131 R*/s, M=60 R*/s, S=1973 R*/s
L = numel([1 NCONTRASTS]);
count = 1;
noiseparams.sampTime = theMosaic.integrationTime;  % Sec
cone_mosaic = theMosaic.pattern;
c = ['r','g','b'];
for ii = [1 NCONTRASTS];
    coneCurrentwnoise = zeros(size(coneCurrent{1}));
    coneCurrentwonoise = RGB2XWFormat(coneCurrent{ii});
    coneIsoRatewnoise = zeros(size(coneIsomerizationRate{1}));
    coneIsoRatewonoise =  RGB2XWFormat(coneIsomerizationRate{ii});
    for jj = 1:1        
        coneCurrentwnoise = coneCurrentwnoise + osAddNoise(coneCurrent{ii},noiseparams); % Adds noise to cone photocurrents
        coneIsoRatewnoise = coneIsoRatewnoise + theMosaic.photonNoise(coneIsomerizationRate{ii});
    end
    coneCurrentwnoise = RGB2XWFormat(coneCurrentwnoise);
    coneIsoRatewnoise = RGB2XWFormat(coneIsoRatewnoise);   
%    coneCurrentwnoise = (((RGB2XWFormat(coneCurrentwnoise)))-RGB2XWFormat(coneCurrent{1})).*cur_matchedfilter;
%    coneIsoRatewnoise = (((RGB2XWFormat(coneIsoRatewnoise)))- RGB2XWFormat(coneIsomerizationRate{1})).*iso_matchedfilter;
    for cone_type = 2:4
        cone_locations = find(cone_mosaic==cone_type);
        coneCurrentSingleTypewnoise = (coneCurrentwnoise(cone_locations,:));
        coneCurrentSingleTypewonoise = (coneCurrentwonoise(cone_locations,:));
        coneIsoRateSingleTypewonoise = (coneIsoRatewonoise(cone_locations,:));
        coneIsoRateSingleTypewnoise = (coneIsoRatewnoise(cone_locations,:));
        figure(10),subplot(L,4,4*count-3),plot(sampleTimes,coneIsoRateSingleTypewonoise',c(cone_type-1)); xlabel('Time(s)'); ylabel('*R/s'); hold on;
        subplot(L,4,4*count-2),plot(sampleTimes,coneIsoRateSingleTypewnoise',c(cone_type-1)); xlabel('Time(s)');ylabel('*R/s'); hold on;
        subplot(L,4,4*count-1),plot(sampleTimes,coneCurrentSingleTypewonoise',c(cone_type-1)); xlabel('Time(s)');ylabel('pA'); hold on;
        subplot(L,4,4*count),plot(sampleTimes,coneCurrentSingleTypewnoise',c(cone_type-1)); xlabel('Time(s)');ylabel('pA^2'); hold on; drawnow;
    end
    count = count + 1;
end
figure(10); subplot(L,4,1); title('Photosiomerizations');
subplot(L,4,2); title('Noisy Photosiomerizations * matched filter');
subplot(L,4,3);title('Cone photocurrents');
subplot(L,4,4);title('Noisy photocurrent * matched filter');
hold off;

figure(11); plot(theMosaic.os.lmsConeFilter,'LineWidth',2); xlabel('frames'), ylabel('pA/*R'); title('Temporal Impulse Response Function');
legend('L','M','S');

%%
targetCone = 150;
figure(12);clf; subplot(1, 2, 1);hold on; 
plot(coneCurrentwonoise(targetCone, :));
plot(coneCurrentwnoise(targetCone, :));
subplot(1, 2, 2); hold on
plot(coneIsoRatewonoise(targetCone, :)); %-coneIsoRatewonoise(targetCone, 1));
plot(coneIsoRatewnoise(targetCone, :));

