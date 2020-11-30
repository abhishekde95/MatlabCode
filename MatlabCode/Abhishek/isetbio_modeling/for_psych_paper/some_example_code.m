% Demonstration of a pre-image effect in the bipolar model

% The following program returns a movie of the bipolar response when a

% vertical grating is shown for a brief period in front of a uniform

% background. When the composition mode of oiSequence is set to 'blend',

% the bipolar response "knows" what the grating looks like before the

% grating appears

if ~(exist('bpL'))
    
%% Parameters

close all; clearvars;
p.meanLuminance=200;
 
p.sceneSize=64;
    
p.fov=1.0;
    
p.sceneDistance=0.7;
    
p.cellType='onmidget';
    
p.frequency=2;
    
p.contrast=1;
    
p.phase=pi;
    
p.angle=0;
    
p.leadUp=1;
    
p.coolDown=0.5;
    
p.duration=0.2;
    
p.integrationTime=0.025;
    
p.sRFcenter=5;
    
p.sRFsurround=15;
    
%% Creating the background scene
    
fprintf('Generating scenes...\n');
    
uniformScene=sceneCreate('uniform equal photon',p.sceneSize);
    
uniformScene=sceneSet(uniformScene,'wAngular',p.fov);
    
uniformScene=sceneSet(uniformScene,'distance',p.sceneDistance);
    
uniformScene=sceneAdjustLuminance(uniformScene,p.meanLuminance);
    
%% Creating the grating
    
gratingParams.freq=p.frequency;
    
gratingParams.contrast=p.contrast;
    
gratingParams.ph=p.phase;
    
gratingParams.ang=p.angle;
    
gratingParams.row=p.sceneSize;
    
gratingParams.col=p.sceneSize;
    
gratingParams.GaborFlag=0;
    
gratingScene=sceneCreate('harmonic',gratingParams);
    
gratingScene=sceneSet(gratingScene,'distance',p.sceneDistance);
    
gratingScene=sceneAdjustLuminance(gratingScene,p.meanLuminance);
    
gratingScene=sceneSet(gratingScene,'wAngular',p.fov);
    
%% Stimulus modulation function
    
% Creating a time axis
    
leadUp=p.leadUp;
    
coolDown=p.coolDown;
    
fullTime=p.duration+leadUp+coolDown;
    
timeAxis=linspace(0,fullTime,fullTime/p.integrationTime);
    
% Modulation function?there's a uniform background, a grating flashes on,
    % and it goes away.
    
for ii=1:length(timeAxis)
        
modulationFunction(ii)=1*(timeAxis(ii)>=leadUp && timeAxis(ii)<=(fullTime-coolDown));
    
end
    
%% Optical Image
    
fprintf('Creating human optics model...\n');
    
% Create human optics
    
oi=oiCreate('human');
    
oi=oiCompute(oi, uniformScene);
    
oiBackground=oiCompute(oi,uniformScene);
    
oiModulated=oiBackground;
    
oiModulated=oiCompute(oi, gratingScene);
    
% Create image sequence
    
imageSequence=oiSequence(oiBackground,oiModulated,timeAxis,...
        modulationFunction,'composition','blend');

%% Cone mosaic response
    fprintf('Creating cone mosaic...\n');

cMosaic=coneMosaic('integrationTime',p.integrationTime);
    
cMosaic.setSizeToFOV(p.fov);
    
% Eye movement
    
cMosaic.emGenSequence(length(timeAxis));
    
% Set eye movement to zero
    
cMosaic.emPositions=cMosaic.emPositions*0;
    
% Compute currents
    
cMosaic.computeForOISequence(imageSequence);
    
cMosaic.computeCurrent;
    
    
%% Bipolar response
    
fprintf('Generating bipolar response...\n');
    
bpL=bipolarLayer(cMosaic);
    
bpL.mosaic{1}=bipolarMosaic(cMosaic,p.cellType,...
        'sRFcenter',p.sRFcenter,...
        'sRFsurround',p.sRFsurround);
    
bpL.mosaic{1}.compute;

end


%% Plots
% As these plots show, even though a grating has not been shown and the

% cone mosaic has not given the grating response, the bipolar layer gives a

% grating response.
figure
subplot(1,2,1)

contourf(cMosaic.current(:,:,1));

title('Cone mosaic photocurrent at t=0');
pbaspect([1,1,1]);
subplot(1,2,2)
contourf(bpL.mosaic{1}.responseCenter(:,:,1));

title('Bipolar center response at t=0');
pbaspect([1,1,1]);
