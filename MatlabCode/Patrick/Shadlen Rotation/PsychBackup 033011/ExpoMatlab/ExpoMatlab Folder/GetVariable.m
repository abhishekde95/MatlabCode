function variableName = GetVariable(expoDataSet, paramName, routine, block, routineInstance)
%function variableName = GetVariable(expoDataSet, parameter, routine, block, routineInstance)
%
% This Expo utility returns the name of a variable specified by parameter,
% routine, block, and routineInstance.  The routine and block can be
% specified either by name or by ID.  The routineInstance is optional and 0 based.
% If there are several instances of the routine in the block you can specify the
% instance using 0 for the first instance, 1 for the second and so on. If you
% do not provide the routineInstance, a value of 0 is assumed.
%
% See also ReadExpoXML, GetPasses, GetEvents, GetSpikeTimes, GetAnalog,
% GetPSTH, PlotPSTH, GetWaveforms, GetStartTimes, GetEndTimes, GetDuration,
% GetConversionFactor, MergeExpoDataSets, GetTransitionProbabilities.
%
%   Author:      Julian Brown
%	Version:     1.1
%   Last updated:  2005-03-28
%   E-mail:      julian@monkeybiz.stanford.edu


variableName = '';

matlabImportVersion = '1.1';
CheckExpoVersion(expoDataSet, matlabImportVersion);

if isnumeric(block)
    blockID = block;
    if blockID > length(expoDataSet.blocks.IDs) - 1
        error('BlockID %i is out of range.  The maximum block ID is %i.', blockID, length(expoDataSet.blocks.IDs)-1)
    end
    blockName = expoDataSet.blocks.Names{1, blockID + 1};
else
    blockName = block;
    blockID = expoDataSet.blocks.IDs(strmatch(block, expoDataSet.blocks.Names, 'exact'));
    if length(blockID) < 1
        disp('Available block names are:')
        disp(expoDataSet.blocks.Names')
        error('No matching block name found. Aborting')
    end
end

routines = expoDataSet.blocks.routines{1, blockID + 1};

if ~exist('routineInstance');
    routineInstance = 0;
else
    if ~isnumeric(routineInstance)
        error('routineInstance must be a 0 based numeric value.')
    end
end


if isnumeric(routine)
    routineNums = find(routines.IDs == routine);
    if length(routineNums) < 1
        if length(routines.IDs)>0
            disp('Available routine IDs are:')
            disp(unique(routines.IDs))
            error('No matching routine ID found in block %s', blockName);
        else
            error('Block %s has no routines.', blockName);
        end
    end
    if length(routineNums) < routineInstance + 1
        error('Invalid routineInstance value.  There are only %i instances of routine ID %i in block %s and instance numbers start from 0.', length(routineNums), routine, blockName);
    end
    routineName = routines.Names{1, routineNums(routineInstance + 1)};
else
    routineName = routine;
    routineNums = strmatch(routine, routines.Names, 'exact');
    if length(routineNums) < 1
        disp('Available routine names are:')
        disp(routines.Names)
        error('No matching routine name found in block %s', blockName);
    end
    if length(routineNums) < routineInstance + 1
        error('Invalid routineInstance value.  There are only %i instances of routine %s in block %s and instance numbers start from 0.', length(routineNums), routine, blockName);
    end
end

params = routines.Params{1, routineNums(routineInstance + 1)};

paramNum = strmatch(paramName, params.Names, 'exact');
if length(paramNum) < 1
    disp('Available parameter names are:')
    disp(params.Names)
    error('No matching parameter name found in routine %s of block %s', routineName, blockName);
end

variableName = params.VarNames(paramNum);



