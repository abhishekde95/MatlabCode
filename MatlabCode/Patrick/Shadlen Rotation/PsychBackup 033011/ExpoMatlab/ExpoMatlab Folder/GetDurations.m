function durations = GetDurations(expoDataSet, passIDs, unit)
% function durations = GetDurations(expoDataSet, passIDs, unit)
%
% See also ReadExpoXML, GetSlots, GetPasses, GetEvents, GetSpikeTimes, GetAnalog, 
% GetPSTH, PlotPSTH, GetStartTimes, GetEndTimes, GetDuration,
% GetConversionFactor.
%
%   Author:      Jim Muller
%   Last updated:  2005-03-28
%	Version:     1.1
%   E-mail:      jim@monkeybiz.stanford.edu

    matlabImportVersion = '1.1';
    CheckExpoVersion(expoDataSet, matlabImportVersion);
    
    units = expoDataSet.environment.Conversion.units;
    
    if exist('unit') && length(unit) > 0
        unitNum = GetUnitNumFromName(units, unit, units.U_BASETIME); 
    else
        unitNum = units.U_BASETIME; 
    end

    [numOfPasses passIDs] = TransformToColumnVector(passIDs);

    if ~isnumeric(passIDs)
        error('passIDs should be a vector of numbers. It appears to be non numeric.');
    end

    durations = double(expoDataSet.passes.EndTimes(passIDs+1) - expoDataSet.passes.StartTimes(passIDs+1));

    if unitNum ~= units.U_BASETIME;
        durations = ConvertNum(units, durations, units.U_BASETIME, unitNum);
    end


