function passIDs = GetPasses(expoDataSet, slot, offset, matrixElement)
% function passIDs = GetPasses(expoDataSet, slot, offset, matrixElement)
%
% This Expo utility takes a slotID or slot label and returns a vector
% of passes that have that slot ID or label.  It then adds an offset (default 0) to all of the IDs
% returning NaN if any go out of range.
%
% The matrixElement is an optional parameter that allows you to select
% passes that belong to a particular Matrix block associated with that slot.
%
% See also ReadExpoXML, GetSlots, GetEvents, GetSpikeTimes, GetAnalog,
% GetPSTH, PlotPSTH, GetWaveforms, GetStartTimes, GetEndTimes, GetDuration,
% GetConversionFactor, MergeExpoDataSets, GetTransitionProbabilities.
%
%   Author:      Julian Brown
%	Version:     1.1
%   Last updated:  2005-03-29
%   E-mail:      julianb@stanford.edu


matlabImportVersion = '1.1';
CheckExpoVersion(expoDataSet, matlabImportVersion);

if ~exist('offset'), offset = 0;, end

if isnumeric(slot)
    slotID = slot;
else
    % retrieve slotID from slot array by matching with the slot label
    slotID = expoDataSet.slots.IDs(strmatch(slot, expoDataSet.slots.Labels, 'exact'));
    if size(slotID, 2) < 1
        error(sprintf('Could not find a match for the slot label %s', slot));
    end
end

passIDs = double(expoDataSet.passes.IDs(find (expoDataSet.passes.SlotIDs == slotID)));

if length(passIDs) == 0
    if isnumeric(slot) 
        warning(sprintf('There are no passes associated with slot %d', slot));
    else
        warning(sprintf('There are no passes associated with slot %s', slot));
    end
    return
end

passIDs = passIDs + offset;


% now check for boundary conditions
minPassID = min(expoDataSet.passes.IDs); % probably always 0 but let's not assume
maxPassID = max(expoDataSet.passes.IDs);

tooLow = find(passIDs<minPassID);
passIDs(tooLow) = NaN;

tooHigh = find(passIDs>maxPassID);
passIDs(tooHigh) = NaN;

if exist('matrixElement')
    
    if ~isfield(expoDataSet, 'matrix')
        error('You requested a matrixElement but the expoDataSet does not contain any matrix information');
    end
    
    blockIDs = expoDataSet.passes.BlockIDs(passIDs + 1);
    
    if length(blockIDs)==0 || max(blockIDs) < expoDataSet.matrix.MatrixBaseID || min(blockIDs) > expoDataSet.matrix.MatrixBaseID + expoDataSet.matrix.NumOfBlocks
        warning('The requested slot/offset did not yield any passes in the Matrix blocks');
    end
    
    if  ~isnumeric(matrixElement)
        error('matrixElement should be a vector of numbers.  It appears to be non numeric.');
    end

	numOfDimensions = expoDataSet.matrix.NumOfDimensions;

    if length(matrixElement) ~= numOfDimensions;
        error(sprintf('matrixElement is of length %d but the expoDataset''s matrix has %d dimensions.  These numbers should be the same.',length(matrixElement), expoDataSet.matrix.NumOfDimensions));
    end
    
    if numOfDimensions ~= 0
		pos = 0;
		for i=numOfDimensions:-1:1
			value = matrixElement(i);
			dimensionSize = expoDataSet.matrix.Dimensions{i}.Size;

			if value < 0
				error(sprintf('Value #%d=%d of matrixElement is negative. Values for this dimension should range from 0 to %d.', i, value, dimensionSize));
			end

			if value > dimensionSize - 1
                if dimensionSize ==1
                    error(sprintf('Value #%d=%d of matrixElement exceeds the dimension size. Values for this dimension can only be 0.', i, value));
                else
                    error(sprintf('Value #%d=%d of matrixElement exceeds the dimension size. Values for this dimension should range from 0 to %d', i, value, dimensionSize-1));
                end
			end

			if i==numOfDimensions
				pos = value;
			else
				pos = pos * dimensionSize + value;
			end
		end

		% pos now points to the block position within the matrix so add to matrixBaseID
		blockID = expoDataSet.matrix.MatrixBaseID + pos;
        
        i = find(~isnan(passIDs));

		% remove any NaNs in the passIDs vector (can arise from the offset)
        passIDs = passIDs(find(~isnan(passIDs)));
        
        % now filter out all passes not associated with this blockID
		passIDs = double(passIDs(find (expoDataSet.passes.BlockIDs(passIDs + 1) == blockID)));

	end
end

