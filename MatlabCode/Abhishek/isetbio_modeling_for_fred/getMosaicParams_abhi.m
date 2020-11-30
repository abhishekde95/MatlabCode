function mosaicParams = getMosaicParams_abhi( ecc_inMM, fov, temporalParams )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
mosaicParams.resamplingFactor = 1;
mosaicParams.spatiallyVaryingConeDensity = true;
mosaicParams.fieldOfViewDegs = fov;
mosaicParams.macular = false;
mosaicParams.LMSRatio = [0.45 0.45 0.1];
mosaicParams.center = [ecc_inMM 0]*1e-3;
mosaicParams.timeStepInSeconds = temporalParams.stimulusSamplingIntervalInSeconds;
mosaicParams.integrationTimeInSeconds = mosaicParams.timeStepInSeconds;
mosaicParams.osNoise = false;% Adding noise, Fred's model
mosaicParams.noiseFlag = false;
mosaicParams.osModel = 'biophys';
mosaicParams.append = false;

end

