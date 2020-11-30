function unitNum = GetUnitNumFromName(units, unitName, compatibleUnit)
%function unitNum = GetUnitNumFromName(units, unitName)
%
% An Expo helper function used for returning the unitNum from the unitName
% as part of unit conversions
%
% Used in GetEvents
%
%   Author:      Julian Brown
%   Last updated:  2004-12-22
%   E-mail:      julianb@stanford.edu

            
    unitNum = strmatch (unitName, units.Names, 'exact') + 1;
    
    if size(unitNum, 1) == 0
        if exist('compatibleUnit')
            error(sprintf('Could not find the requested unit.\n\nCompatible unit types are:%s', GetListOfCompatibleUnitTypes(units, compatibleUnit)));
        else
            error(sprintf('Could not find the requested unit.\n\nAvailable unit types and groups are:%s', GetListOfTypes(units)));
        end
    end

return


function unitNamesList = GetListOfTypes(units)

    unitNamesList = '';
    
    % build a list of unit types and groups
    for i=1:size(units.GroupRepUnits,2);
        repUnitNum = units.GroupRepUnits(i);
        isFirstInRow = 1;
        for j=1:size(units.Names,2);
            if units.ConversionMatrix(repUnitNum,j+1)>0
                if isFirstInRow
                    unitNamesList = [unitNamesList, sprintf('\n%20s', units.Names{j})];
                    isFirstInRow = 0;
                else
                    unitNamesList = [unitNamesList, sprintf('%20s', units.Names{j})];
                end
            end
        end
    end
    
return

function unitNamesList = GetListOfCompatibleUnitTypes(units, unitNums)
    
    unitNamesList = '';
    
    % build a list of compatible types
    foundNames = 0;
    for i=unitNums'
        for j=1:size(units.Names,2);
            if units.ConversionMatrix(i,j+1)>0
                foundNames = 1;
                unitNamesList = [unitNamesList, sprintf('`%s` ', units.Names{j})];
            end
        end
        if foundNames ==1, break, end
    end
   
return
