
function mosaicParams = getMosaicParams_abhi( ecc_inMM, fov, temporalParams )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
mosaicParams.resamplingFactor = 1; % Am not sure what exactly this variable means
mosaicParams.spatiallyVaryingConeDensity = logical(1);
mosaicParams.fieldOfViewDegs = fov;
mosaicParams.macular = false;
mosaicParams.LMSRatio = [0.45 0.45 0.1];
mosaicParams.center = [ecc_inMM 0]*1e-3;
mosaicParams.timeStepInSeconds = temporalParams.stimulusSamplingIntervalInSeconds;
mosaicParams.integrationTimeInSeconds = mosaicParams.timeStepInSeconds;
mosaicParams.osNoise = 'none';% Adding noise, Fred's model
mosaicParams.noiseFlag = 'none';
mosaicParams.osModel = 'biophys';
mosaicParams.append = false;

end

