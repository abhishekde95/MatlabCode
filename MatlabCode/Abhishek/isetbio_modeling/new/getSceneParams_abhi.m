function [gaborParams, colorModulationParams, backgroundParams] = getSceneParams_abhi(fov)

gaborParams.type = 'Gabor';
gaborParams.windowType = 'Gaussian';
gaborParams.fieldOfViewDegs = fov;
gaborParams.gaussianFWHMDegs = 0.942; % 2.3548 x standard deviation 
gaborParams.cyclesPerDegree = 1;
gaborParams.row = 128;
gaborParams.col = 128;
gaborParams.ang = 0;
gaborParams.ph = 0;
gaborParams.viewingDistance = 1.0;

colorModulationParams.type = 'ColorModulation';
colorModulationParams.contrast = 1;
% colordir = [0.06 -0.06 0.45]';% S-ON direction
% colordir = [0.06 -0.06 -0.45]';% S-OFF direction
% colordir = [0.06 0.06 0.45]';% luminance direction
% colordir = [1 -1 0]';% L-M direction
% colordir = [0 0 1]';% Siso direction
colordir = [0 0 -1]';% negSiso direction
colorModulationParams.coneContrasts = colordir/norm(colordir);

backgroundParams.type = 'Background';
backgroundParams.backgroundxyY = [0.33 0.33 3.98]'; % [0.27 0.30 3.98]; [X Y intensity];
backgroundParams.monitorFile = 'CRT-MODEL';
backgroundParams.leakageLum = 1.0;
backgroundParams.lumFactor = 1;

end

