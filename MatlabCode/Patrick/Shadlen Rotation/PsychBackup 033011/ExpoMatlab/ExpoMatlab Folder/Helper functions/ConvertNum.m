function convertedNum = ConvertNum(units, numToConvert, captureUnitNums, requestedUnitNums)
%function convertedNum = ConvertNum(units, numToConvert, captureUnitNums, requestedUnitNums)
%
% An Expo helper function used in conjunction with GetUnitNumFromName for converting values 
% according to the requested unit. captureUnitNums and requestedUnitNums may be vectors
% because some unit names are ambiguous and in these cases GetUnitNumFromName returns 
% a vector. If only one or the other of these variables is ambiguous
% ConverNum is able to resolve the ambiguity by looking to see which
% results are not zero.
%
% Used in GetEvents, GetAnalog, GetSpikeTimes, GetWaveforms etc.
%
%   Author:      Julian Brown
%   Last updated:  2004-12-28
%   E-mail:      julianb@stanford.edu
    

    % if -1 then no unit provided and no conversion needed
    if requestedUnitNums == -1
        convertedNum = numToConvert;
        return
    end
    
    if max(captureUnitNums) < 2
        error(sprintf('Incompatible unit requested: the requested parameter has a capture unit value of %d which suggests it cannot be converted.', captureUnitNum));
    end
    
    conversionFactor = full(units.ConversionMatrix(captureUnitNums, requestedUnitNums));
    
    if size(captureUnitNums, 1) > 1
        % get the only compatible unit in case of ambiguous numbers 
        captureUnitNum = (conversionFactor'>0) * captureUnitNums;
        conversionFactor = (conversionFactor'>0) * conversionFactor;
    else
        captureUnitNum = captureUnitNums;
    end
    
    if size(requestedUnitNums, 1) > 1
        % get the only compatible unit in case of ambiguous numbers 
        requestedUnitNum = (conversionFactor>0) * requestedUnitNums;
        conversionFactor = (conversionFactor>0) * conversionFactor';
    else
        requestedUnitNum = requestedUnitNums;
    end
    
    if conversionFactor ==0 
        RaiseErrorInUnitConversion(units, captureUnitNum, requestedUnitNums);
    end
   
    doRecipricol = units.ReciprocalMatrix(captureUnitNum, requestedUnitNum);
    
    if doRecipricol == 0
        convertedNum = conversionFactor .* numToConvert;
    else
        convertedNum = conversionFactor ./ numToConvert;
    end
return


function RaiseErrorInUnitConversion(units, captureUnitNum, requestedUnitNum)
    error(sprintf('Incompatible unit requested: the parameter`s capture unit is %s whereas you asked for %s.\nCompatible units are:%s', ...
        GetUnitNameFromNum(units, captureUnitNum), GetUnitNameFromNum(units, requestedUnitNum), GetListOfCompatibleUnitTypes(units, captureUnitNum)));
return


function unitNamesList = GetListOfCompatibleUnitTypes(units, unitNums)
    
    unitNamesList = '';
    
    % build a list of compatible types
    foundNames = 0;
    for i=unitNums'
        for j=1:size(units.Names,2);
            if units.ConversionMatrix(i,j+1)>0
                foundNames = 1;
                unitNamesList = [unitNamesList, sprintf('%20s', units.Names{j})];
            end
        end
        if foundNames ==1, break, end
    end
   
return