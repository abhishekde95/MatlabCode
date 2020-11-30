textfiles_path = nexfilepath('nexfilelists', 'Greg', 'DTEM');
text_filenames = dir(fullfile(textfiles_path, '*DTscot.txt'));
text_files = strcat(textfiles_path, filesep, {text_filenames.name});
text_files(~cellfun('isempty', regexp(text_files, 'Zack|GregDT|Leah'))) = [];

for text_file = text_files
    nominal_filenames = fnamesFromTxt2(text_file{1});
    valid_names = isvalidnexfilename(flatten(nominal_filenames));
    nex_filepaths = findfile(nominal_filenames(valid_names));
    [~, titles] = cellfun(@fileparts, flatten(nex_filepaths), 'unif', 0);
    plot_fixation_positions(nex_filepaths, titles, [-5 5 -5 5], true);
end
