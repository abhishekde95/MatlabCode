function theMosaic = createConeMosaic_abhi(mosaicParams,option)

% Cone mosaic of interest
if strcmp(option,'coneMosaicHex')
    theMosaic = coneMosaicHex(mosaicParams.resamplingFactor,...
        mosaicParams.spatiallyVaryingConeDensity,[],...
        'center',mosaicParams.center);
else
    theMosaic = coneMosaic( 'integrationTime',mosaicParams.integrationTimeInSeconds, ...
        'center',mosaicParams.center);
end

if (isfield(mosaicParams, 'fieldOfViewDegs'))
    theMosaic.setSizeToFOV(mosaicParams.fieldOfViewDegs);
end

if (isfield(mosaicParams, 'isomerizationNoise'))
    theMosaic.noiseFlag = mosaicParams.isomerizationNoise;
end
if (isfield(mosaicParams, 'osNoise'))
    theMosaic.os.noiseFlag = mosaicParams.osNoise;
end
if (isfield(mosaicParams, 'noiseFlag'))
    theMosaic.noiseFlag = mosaicParams.noiseFlag;
end

if (isfield(mosaicParams, 'timeStepInSeconds'))
    theMosaic.os.timeStep = mosaicParams.timeStepInSeconds;
end

if (isfield(mosaicParams, 'LMSRatio'))
    if (numel(mosaicParams.LMSRatio) == 3)
        theMosaic.spatialDensity = [0 mosaicParams.LMSRatio(1) mosaicParams.LMSRatio(2) mosaicParams.LMSRatio(3)]';
    else
        theMosaic.spatialDensity = mosaicParams.LMSRatio(:);
    end
end
end

