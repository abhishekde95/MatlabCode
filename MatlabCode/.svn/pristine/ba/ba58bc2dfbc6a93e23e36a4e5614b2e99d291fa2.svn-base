function endTimes = GetEndTimes(expoDataSet, passIDs, unit)
% function endTimes = GetEndTimes(expoDataSet, passIDs, unit)
%
% 
% See also ReadExpoXML, GetSlots, GetPasses, GetEvents, GetSpikeTimes, GetAnalog, 
% GetPSTH, PlotPSTH, GetStartTimes, GetEndTimes, GetDuration,
% GetConversionFactor.
%
%   Author:      Jim Muller
%	Version:     1.1
%   Last updated:  2005-03-28
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

    endTimes = double(expoDataSet.passes.EndTimes(passIDs+1));

    if unitNum ~= units.U_BASETIME;
        endTimes = ConvertNum(units, endTimes, units.U_BASETIME, unitNum);
    end


