function slotIDs = GetSlots(expoDataSet, blockName)
%function slotIDs = GetSlots(expoDataSet, blockName)
%
% This Expo utility returns a vector of slotIDs belonging to the named
% block.  If the blockname is 'MATRIX' then the function returns the 
% the slot that has blockID -1 if there is such a slot.
%
% See also ReadExpoXML, GetPasses, GetEvents, GetSpikeTimes, GetAnalog, 
% GetPSTH, PlotPSTH, GetWaveforms, GetStartTimes, GetEndTimes, GetDuration,
% GetConversionFactor, MergeExpoDataSets, GetTransitionProbabilities.
%
%   Author:      Julian Brown
%	Version:     1.1
%   Last updated:  2005-03-28
%   E-mail:      julian@monkeybiz.stanford.edu


slotIDs = [];

matlabImportVersion = '1.1';
CheckExpoVersion(expoDataSet, matlabImportVersion);

numOfBlocks = size(expoDataSet.blocks.Names, 2);

if length(blockName) == 0
    if numOfBlocks<30
        disp('Available block names are:');
        disp(expoDataSet.blocks.Names')
    end
    error('blockName was empty. Aborting')
    return
end

if strcmp(blockName, 'MATRIX')
    slotNums = find(expoDataSet.slots.BlockIDs == -1);
else

    blockNum = strmatch(blockName, expoDataSet.blocks.Names, 'exact');

    if size(blockNum, 1) == 0
        if numOfBlocks<30
            disp('Available block names are:')
            disp(expoDataSet.blocks.Names')
        end
        error('No matching block name found. Aborting')
        return
    end

    if isfield(expoDataSet, 'matrix') 
    	if isempty(expoDataSet.matrix)
			expoDataSet = rmfield(expoDataSet,'matrix');
		end
	end

    if isfield(expoDataSet, 'matrix') 
        if expoDataSet.blocks.IDs(blockNum) >= expoDataSet.matrix.MatrixBaseID
            slotNums = find(expoDataSet.slots.BlockIDs == -1);
        else
            slotNums = find(expoDataSet.slots.BlockIDs == expoDataSet.blocks.IDs(blockNum));
        end
    elseif isfield(expoDataSet.blocks, 'MatrixBaseID') % for backward compatibility with Expo v1.0
        if expoDataSet.blocks.IDs(blockNum) >= expoDataSet.blocks.MatrixBaseID
            slotNums = find(expoDataSet.slots.BlockIDs == -1);
        else
            slotNums = find(expoDataSet.slots.BlockIDs == expoDataSet.blocks.IDs(blockNum));
        end
    else
       slotNums = find(expoDataSet.slots.BlockIDs == expoDataSet.blocks.IDs(blockNum)); 
    end
end
   
slotIDs = double(expoDataSet.slots.IDs(slotNums));

return