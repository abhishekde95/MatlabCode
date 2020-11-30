%Part of an overly complicated solution to a problem that didn't need any
%of this nonsense. Do not use.
function p_val_map = analysis(actual_raw_data_map, distinct_IDs)
%[actual_raw_data_map, distinct_IDs] = pullDataForAnalysis();
full_keys = keys(actual_raw_data_map);
false_error_ratios = {}; p_val_map = containers.Map(); num_runs = 100; %eventually 1000
for i = 1:length(distinct_IDs)
    key_IDX = full_keys(~cellfun(@isempty,strfind(full_keys, distinct_IDs{i})));
    curr_vals = values(actual_raw_data_map, key_IDX);
    sizes = cellfun('length', curr_vals);
    composite_IDX = find(sizes==max(sizes));
    other_IDX = curr_vals(sizes~=max(sizes));
    composite_subset_ID = key_IDX(composite_IDX);
    composite_val = curr_vals(composite_IDX);
    composite_err = error_fit(composite_val);
    actual_error = error_fit(other_IDX);
    actual_error_ratio = composite_err/sum(actual_error);
    n = 1;
    while n <= num_runs
        alt_vals = {}; temp_comp = composite_val{1}; %should I shuffle it? Does it matter?
        max_size = (size(temp_comp, 1))-2;
        for s = 1:size(other_IDX, 2)
            num_samps = randi([4 max_size], 1);
            [comp_sample, sample_IDX] = datasample(temp_comp, num_samps);                
            if s==size(other_IDX,2)
                alt_vals{end+1} = temp_comp;
            else
                alt_vals{end+1} = temp_comp(unique(sample_IDX), :);
                temp_comp(sample_IDX,:) = [];
                max_size = (size(temp_comp, 1)-1);
            end
        end
        try
            false_error = error_fit(alt_vals);
        catch
            keyboard;
        end
        false_error_ratios{end+1} = composite_err/sum(false_error);
        n = n+1;
    end
    figure('Name',composite_subset_ID{1},'NumberTitle','off'); hold on;
    h = histogram(cell2mat(false_error_ratios)); 
    xlabel('error ratios');
    ylabel('num occurrences');
    %BELOW IS WRONG, CHECK THAT IT WORKS FIRST THEN DETERMINE SIGNIFICANCE
    average_error = mean(h.Data);
    p_value = abs(actual_error_ratio - average_error);
    z = h.BinEdges(h.BinEdges > p_value);
    z = z(1:2);
    z = [z fliplr(z)];
    y = [0 0 repmat(p_value, 1, 2)];
    patch(z, y, 'red');
    sig = 'can''t tell automagically yet';
    title_string = sprintf('p-val is %d and therefore we %s', p_value, sig);
    title(title_string);
    keyset = {'composite_subset_ID', 'actual_error_ratio', 'average_error', 'p_value', 'sig'};
    valueset = {composite_subset_ID{1}, actual_error_ratio, average_error, p_value, sig};
    temp_map = containers.Map(keyset, valueset);
    p_val_map = [p_val_map; temp_map];
end

    function error = error_fit(tot_dat)
        error = []; %data = {}; modelParams = {};
        for t = 1:length(tot_dat)
            d = cell2mat(tot_dat(t));
            modelParams = LMTFfitting(d);
            %if we're not using more than one eccentricity (which we're not
            %atm) then we don't have to use Lecc for model params
            predr = LMTF_thresh_from_model(d(:,1),d(:,2),d(:,3), modelParams); % vector lengths in cc space
            r = sqrt(d(:,1).^2+d(:,2).^2);
            resids = r - predr;
            error = [error ; sum(abs(resids))];
        end
    end
end