function [slotAndBlockNames slotIDs slotPortions blockIDs] ...
    = GetTransitionProbabilities(expoDataSet, passIDs, offsets, pruneSlots, printme)
% function [slotAndBlockNames slotIDs slotPortions blockIDs] ...
%   = GetTransitionProbabilities(expoDataSet, passIDs, offsets, pruneSlots, printme)
%
% Takes an expoDataSet and a list of passIDs (e.g. from GetPasses) and
% computes the probability of transitioning into various state/blocks
% therefrom as a function of slot offset.
% [optional]
% offsets: a vector of offsets treated as per GetPasses [default 0].
% pruneSlots: Prune the tree if execution passes through one of these slots (a vector of slotIDs or a cell 
%   array containing a mixture of string slot labels and slotIDs).  This is useful for stopping computation 
%   if a trial ends, is aborted, etc.  Specifically, as the computation of transition probabilities 
%   progresses through increasing offsets above 0, if we encounter one of pruneSlots,
%   we stop tracking it.  Independently this is also done as decreasing offsets progress below 0.
% print: prints a table of the results before returning [default 1].
% Returns several arguments of the same form.  Each is a cell array, each row r of which
%   gives statistics about which unique slots are contained in passIDs+offsets(r), and blocks within them
%   (if a slot contains multiple blocks, i.e. a MATRIX, just 1 example block is considered):
% slotAndBlockNames contains a string for each unique slot of the form "<slotLabel>/<blockName> <slotPortion>%".
% slotIDs contains the corresponding list of unique slotIDs.
% slotPortions is the portion of all passes at this offset that fall within this unique slot (add to 1).
% blockIDs contains the blockIDs of the blocks within those slots.
%
% See also ReadExpoXML, GetEvents, GetSlots, GetPasses, GetSpikeTimes, GetAnalog, 
% GetPSTH, PlotPSTH, GetWaveforms, GetStartTimes, GetEndTimes, GetDuration,
% GetConversionFactor.
%
%   Author:      Jim Muller
%	Version:     1.1
%   Last updated:  2005-03-30
%   E-mail:      jim@monkeybiz.stanford.edu

if ~exist('offsets') offsets = 0; end
if ~exist('pruneSlots') pruneSlots = []; end
if ~exist('printme') printme = 1; end

matlabImportVersion = '1.1';
CheckExpoVersion(expoDataSet, matlabImportVersion);

[numOfPasses passIDs] = TransformToColumnVector(passIDs);

if ~isnumeric(passIDs)
    error('passIDs should be a vector of numbers.  It appears to be non numeric.');
end
if iscell(pruneSlots)
    pruneSlotIDs = [];
    for idx = 1:length(pruneSlots)
        if isnumeric(pruneSlots{idx})
            pruneSlotIDs(idx) = pruneSlots{idx};
        else
            % retrieve slotID from slot array by matching with the slot label
            pruneSlotIDs(idx) = expoDataSet.slots.IDs(strmatch(pruneSlots{idx}, expoDataSet.slots.Labels, 'exact'));
            if size(pruneSlotIDs(idx), 2) < 1 
                error(sprintf('Could not find a match for the label %s', slot));
            end
        end
    end
else
    pruneSlotIDs = pruneSlots;
end
% now check for boundary conditions
minPassID = min(expoDataSet.passes.IDs); % probably always 0 but let's not assume
maxPassID = max(expoDataSet.passes.IDs); % probably the last element, but let's not assume.

initialPassIDs = passIDs; % Save the full passIDs list
offsets = unique(offsets);
start = min(find(offsets >= 0));
% Proceed from tree's ground-level root to top, then from ground-level root to deep roots. 
for offset_idx =  [start:length(offsets) (start-1):-1:1]; 
    if offset_idx == (start-1) % We've returned to the ground-level root of the tree.  
        passIDs = initialPassIDs; % Reset to full passID list. No we will continue from ground-level root to deep roots. 
    end
    pIDs = passIDs + offsets(offset_idx);
    rangeError = pIDs<minPassID | pIDs>maxPassID;
    pIDs = pIDs(~rangeError);
    if isempty(rangeError)
		rangeErrorPortion = 0;
    else
		rangeErrorPortion = mean(rangeError);
    end
    % get unique set of slot/block IDs and Names corresponding to the slots
    [uniqueSlotIDs unique_pass_idxs] = unique(expoDataSet.passes.SlotIDs(pIDs + 1));
    slotIDs{offset_idx,1} = uniqueSlotIDs;
    slotNames = expoDataSet.slots.Labels(uniqueSlotIDs + 1);

    blockIDs{offset_idx,1} = expoDataSet.passes.BlockIDs(pIDs(unique_pass_idxs)+1);% unique(expoDataSet.passes.BlockIDs(pIDs + 1));
    blockNames = expoDataSet.blocks.Names(blockIDs{offset_idx,1} + 1);
   
    % For further progress up the tree, prune slots the user wants pruned, e.g. trial-abort/invalid slots
    passesToPrune = zeros(size(passIDs));
    for idx = 1:length(uniqueSlotIDs)
        if find(uniqueSlotIDs(idx) == pruneSlotIDs) % If this is a slot to prune.
            % ...then add it to the prune mask -- it's based on pIDs, so it only applies to the ~rangeError mask.
            passesToPrune(~rangeError) = passesToPrune(~rangeError) | (uniqueSlotIDs(idx) == expoDataSet.passes.SlotIDs(pIDs + 1)');
        end
    end
    passIDs = passIDs(~passesToPrune);
    
    slotAndBlockNames{offset_idx,1} = {};
    for idx = 1:length(slotNames)
        slotPortions{offset_idx,1}(idx) = mean(slotIDs{offset_idx,1}(idx) == expoDataSet.passes.SlotIDs(pIDs + 1));
        slotAndBlockNames{offset_idx,1}{idx} = sprintf('%-20.20s %s%%      ',[slotNames{idx} '/' blockNames{idx}], num2str(100 * slotPortions{offset_idx,1}(idx),3));
        slotAndBlockNames{offset_idx,1}{idx} = sprintf('%-20.20s %3.0f%%      ',[slotNames{idx} '/' blockNames{idx}], 100 * slotPortions{offset_idx,1}(idx));
   end
    if rangeErrorPortion > 0 
        slotAndBlockNames{offset_idx,1}{end+1} = ['Range errors ' num2str(100 * rangeErrorPortion) '%'];
    end
end
if printme
    disp('---------------------------');
    disp(sprintf('Off %-20s %%', 'Slot/Block'));
    disp('---------------------------');
	for offset_idx = 1:length(offsets)
    	%disp(sprintf('%3d%s', offsets(offset_idx), slotAndBlockNames{offset_idx,1}{:}))
    	disp([num2str(offsets(offset_idx)) ' ' slotAndBlockNames{offset_idx,1}{:}])
    end
end



