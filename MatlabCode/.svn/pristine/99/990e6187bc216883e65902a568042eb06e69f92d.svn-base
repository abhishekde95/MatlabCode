function unitName = GetUnitNameFromNum(units, unitNums)
%function unitName = GetUnitNameFromNum(units, unitNums)
%
% An Expo helper function used for converting values according to the
% requested unit
%
% Used in GetEvents, ConvertNum
%
%   Author:      Julian Brown
%   Last updated:  2004-12-24
%   E-mail:      julianb@stanford.edu
       
    unitName = '';
    
    if size(unitNums,1)>1
        % if more than one unit num just pick the first name we find (the
        % others will most likely be the same as we should only get
        % multiple unitNums when the the user requested an ambiguous name in the first
        % place)
        unitNum = unitNums(1);
    else
        unitNum = unitNums;
    end
    
    if unitNum > size(units.Names, 2) + 1|| unitNum < 0
        error('unitNum is out of range');
    elseif unitNum==0
        unitName = 'None';
    elseif unitNum==1
        unitName = 'Any';
    else
        unitName = units.Names{1, unitNum - 1};
    end
return