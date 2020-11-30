%part of the overly complicated solution to a problem with a simple answer.
%Not useful for anything.E
function p_val_map = LMTF_comparison_analysis(actual_raw_data_map, distinct_IDs)
[actual_raw_data_map, distinct_IDs] = pullDataForAnalysis();
full_keys = keys(actual_raw_data_map);
num_runs = 2000; false_error_ratios = zeros(num_runs,1); p_val_map = containers.Map(); 
for i = 1:length(distinct_IDs)
    key_IDX = full_keys(~cellfun(@isempty,strfind(full_keys, distinct_IDs{i})));
    curr_vals = values(actual_raw_data_map, key_IDX);
    sizes = cellfun('length', curr_vals);
    composite_dataset = vertcat(curr_vals{:});
    composite_err = error_fit({composite_dataset});
    actual_error = error_fit(curr_vals);
    actual_error_ratio = composite_err/sum(actual_error);
    n = 1;
    while n <= num_runs
        alt_vals = {}; temp_comp = composite_dataset; a = 0;
        [temp_comp, temp_IDXs] = Shuffle(temp_comp, 2); %randomly sort rows of composite matrix
        for s = 1:length(sizes)
            alt_vals{s} = temp_comp((a+1):(sizes(s)+a), :);
            a = sizes(s);
        end
        try
            false_error = error_fit(alt_vals);
        catch
            keyboard;
        end
        false_error_ratios(n) = composite_err/sum(false_error);
        n = n+1;
        disp(n);
    end
    fig_string = sprintf('%s ', key_IDX{:});
    figure('Name',fig_string,'NumberTitle','off'); hold on;
    hist(false_error_ratios);
    plot(actual_error_ratio,0, 'm*');
    xlabel('error ratios');
    ylabel('num occurrences');
    p_val = sum(false_error_ratios >= actual_error_ratio)/length(false_error_ratios);
    if p_val <= 0.05
        sig = 'can''t pool across values';
    else
        sig = 'can pool across values';
    end
    title_string = sprintf('p-val is %s and therefore we %s', p_val, sig);
    title(title_string);
    drawnow;
    whereToSave = fullfile('H:/EmilyAnalysisImagesTemp', fig_string);
    saveas(figure(1), whereToSave);
    keyset = {'ID', 'actual_error_ratio', 'p_value', 'sig'};
    valueset = {fig_string, actual_error_ratio, p_val, sig};
    temp_map = containers.Map(keyset, valueset);
    p_val_map = [p_val_map; temp_map];
end

    function error = error_fit(tot_dat)
        error = [];
        for t = 1:length(tot_dat)
            d = cell2mat(tot_dat(t));
            [modelParams, curr_error] = LMTFfitting(d);
            error = [error ; curr_error];
        end
    end
end