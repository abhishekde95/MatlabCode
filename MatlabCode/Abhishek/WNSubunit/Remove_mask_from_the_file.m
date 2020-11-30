function stro = Remove_mask_from_the_file(stro)
% Cuts the part of the structure which has the mask (subunit) and keeps the rest
% of the file. Useful when u want to analyse the cell's response to only
% the white noise pixelated stimuli. Returns back the modified structure.
% Author - Abhishek De, 06/26/2015
global maskidx ntrials
maskidx = strcmp(stro.sum.rasterCells(1,:), 'subunit_mask');
ntrials = size(stro.trial,1);

% Determine when the mask changed in order to determine from what trials we started our experiment to study 
% the computation of the subunits which we are interested in.
mask_changes = 1;
all_masks = stro.ras(:,maskidx);
for k = 2:ntrials
    if isequal(all_masks{k}, all_masks{k-1}) %|| all(all_masks{k} == 0) && any(isnan(all_masks{k-1}))
        continue
    else
        mask_changes = [mask_changes k-1 k]; %#ok<AGROW>
    end
end

% Each column of `mask_changes` holds the start and end trial index where each
% mask was used. For example, the subunits A and B relevant to my project were added on trial number 69.
mask_changes = reshape([mask_changes ntrials], 2, []);
temp = [];
for i = 1:size(mask_changes,2)
    if (mask_changes(2,i) - mask_changes(1,i) > 4)
        temp = [temp, mask_changes(:,i);];
    end
end
mask_changes = temp;

cut_trial = mask_changes(2,1);
stro.ras = stro.ras(1:cut_trial,:);
stro.trial = stro.trial(1:cut_trial,:);
stro.sum.absTrialNum = stro.sum.absTrialNum(1:cut_trial);
% filename = stro.sum.fileName(end-13:end);
% save(horzcat(filename(1:end-4),'_no_mask','.mat'),'-struct','stro');

end

