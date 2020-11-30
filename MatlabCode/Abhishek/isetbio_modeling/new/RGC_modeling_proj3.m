%% Initialize parameters
% Author - Abhishek De, improvising RGC_modeling_proj2.m
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
mosaicParams = getMosaicParams_abhi(ecc_inMM, fov,temporalParams );
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
for ii = [1 NCONTRASTS]
    coneCurrentwnoise = zeros(size(coneCurrent{1}));
    coneCurrentwonoise = RGB2XWFormat(coneCurrent{ii});
    coneIsoRatewnoise = zeros(size(coneIsomerizationRate{1}));
    coneIsoRatewonoise =  RGB2XWFormat(coneIsomerizationRate{ii});
    for jj = 1:100
        coneCurrentwnoise = coneCurrentwnoise + osAddNoise(coneCurrent{ii},noiseparams); % Adds noise to cone photocurrents
        coneIsoRatewnoise = coneIsoRatewnoise + theMosaic.photonNoise(coneIsomerizationRate{ii});
    end
    coneCurrentwnoise = ((0.01*(RGB2XWFormat(coneCurrentwnoise)))-RGB2XWFormat(coneCurrent{1})).*cur_matchedfilter;
    coneIsoRatewnoise = ((0.01*(RGB2XWFormat(coneIsoRatewnoise)))- RGB2XWFormat(coneIsomerizationRate{1})).*iso_matchedfilter;
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
%% Visualize the gaborScene 
disp('Displaying the stimulus');
for jj = NCONTRASTS
    for ii = 1:nSampleTimes
        figure(1); sceneShowImage(gaborScene{jj,ii}); drawnow; pause(0.1);
    end
    pause(0.1); disp('Presenting the next stimulus');
end
%% Next step is to pool responses over noisy data
NTRIALS = 1000; % Number of trials
% Find coordinates of L, M and S cones.
cone_mosaic = theMosaic.pattern;
% The next step is to convolve the 1D filters with the 1D current data at each point in the cone mosaic.
[sz1, sz2, nSteps] = size(coneIsomerizationRate{1});
params = [];
noiseparams.sampTime = theMosaic.integrationTime;  % Sec
theMosaic.os.noiseFlag = false;
cur_fact = [1 1];
mode = [0 1]; % 0-isomerization, 1-photocurrents
neurometric_curve = zeros(numel(cur_fact),NCONTRASTS-1);
tic;
for kk = 1:numel(cur_fact)
    pooledData = zeros(NCONTRASTS,NTRIALS,3); % for 3 different cone types
    disp(kk);
    for ii = 1:NCONTRASTS
        for iter = 1:NTRIALS
            if mode(kk) == 0
                coneisoRS = theMosaic.photonNoise(coneIsomerizationRate{ii}); % Adds noise to the cone isomerizations
                coneisoRS = RGB2XWFormat(coneisoRS-coneIsomerizationRate{1});
                coneCurrentRS = coneisoRS.*iso_matchedfilter; % multiply photoisomerizations by the matched filter
            else
                coneCurrentRS = osAddNoise(cur_fact(kk)*coneCurrent{ii},noiseparams); % Adds noise to cone photocurrents
                coneCurrentRS = RGB2XWFormat(coneCurrentRS-cur_fact(kk)*coneCurrent{1});
                coneCurrentRS = coneCurrentRS.*cur_matchedfilter; % multiply conecurrents by the matched filter
            end
            
            for cone_type = 2:4
                cone_locations = find(cone_mosaic==cone_type);
                coneCurrentSingleType = (coneCurrentRS(cone_locations,:));
                pooledData(ii,iter, cone_type-1) = mean(mean(coneCurrentSingleType,2)); % As of now don't know which one is better
            end
        end
    end
       
    % Perform LDA on the data and find the Linear Discriminant vector that best classifies the data
    for ii = 2:NCONTRASTS
        LDA_vec = [];
        cov_mat = diag(diag(cov(squeeze(pooledData(1,:,:))))');
        LDA_vec = inv(cov_mat) * (mean(squeeze(pooledData(1,:,:))) - mean(squeeze(pooledData(ii,:,:))))'; % find the LDA vector
        LDA_vec = LDA_vec;
        LDA_vec = LDA_vec/norm(LDA_vec);
        % Project the noise and signal onto the Linear Discriminant vector
        noise = squeeze(pooledData(1,:,:))*LDA_vec;
        signal = squeeze(pooledData(ii,:,:))*LDA_vec;
        [counts,centers] = hist([noise signal],20);
        % ROC analysis and neurometric curve - continue
        % Ratio of false positives to true positives
        neurometric_curve(kk,ii-1) = 1-(sum(min(counts,[],2))/(2*NTRIALS));
    end
end
toc;
figure; semilogx(contrast_lattice(2:end),neurometric_curve,'LineWidth',2); 
xlabel('contrast'); title('Neurometric Curve'); legend('Isomerization','photocurrent');

%% Plot the neurometric curve and fit a weibull distribution
[final_model, final_fval, final_success] = weibullFit_iset(contrast_lattice(2:end), neurometric_curve);
model_fit = (1 - 0.5*exp(-((contrast_lattice(2:end)./final_model(1)).^final_model(2))));
figure; semilogx(contrast_lattice(2:end), 100*neurometric_curve,'-bo'); xlabel('Contrast levels'); ylabel('Performance in percentage'); title('Neurometric curve'); hold on;
semilogx(contrast_lattice(2:end), 100*model_fit, '-ro'); legend('data','fit');hold off;
threshold = final_model(1);
slope = final_model(2);

%% Visualize the 3-D data
color = ['r','g'];
figure(3);
for ii = 2: NCONTRASTS
    noisy_data = squeeze(pooledData(1,:,:));
    data = squeeze(pooledData(ii,:,:));
    figure(3), subplot(4,3,ii-1); scatter3(noisy_data(:,1),noisy_data(:,2),noisy_data(:,3),color(1),'filled'); hold on; 
    scatter3(data(:,1),data(:,2),data(:,3),color(2),'filled'); 
    xlabel('L'); ylabel('M'); zlabel('S'); title(strcat('Contrast=',num2str(contrast_lattice(ii)))); 
    set(gca,'Fontsize',13); hold off; drawnow;
end

%% Neurometric curve for different color directions (photocurrents)
load('neurometric_LM.mat');
load('neurometric_Siso.mat');
load('neurometric_negSiso.mat');
load('neurometric_SON.mat');
load('neurometric_SOFF.mat');
figure,semilogx(contrast_lattice(2:end),[neurometric_LM;neurometric_Siso;neurometric_negSiso;neurometric_SON;neurometric_SOFF],'LineWidth',2);
xlabel('contrast');title('Neurometric curve'); legend('LM','Siso','negSiso','+L-M+S','+L-M-S');

%% Neurometric curve for different color directions (photoisomerizations)
load('isometric_LM.mat');
load('isometric_Siso.mat');
load('isometric_negSiso.mat');
load('isometric_SON.mat');
load('isometric_SOFF.mat');

figure,semilogx(contrast_lattice(2:end),[isometric_LM;isometric_Siso;isometric_negSiso;isometric_SON;isometric_SOFF],'LineWidth',2);
xlabel('contrast');title('Isomerization based Neurometric curve');legend('L-M','Siso','negSiso','+L-M+S','L-M-S');
