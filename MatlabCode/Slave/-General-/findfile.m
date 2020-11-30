function found_file = findfile(target_file, start_path)
%
%   This function will take a nex file name (e.g. 'K010309002')
% and find it in the data directory tree (and will append a ".nex"
% if necessary).  The output is the full name (with the path) which can
% be the passed to nex2stro.
%
%   findfile begins searching in `start_path` if it's provided.
%
% ZALB 12/21/11

if ~nargin, return; end

if nargin < 2 || isempty(start_path)
    start_path = {};
end

% if the user passes in a cell, recursively find all files and
% return the results as a cell.
if iscell(target_file)
    found_file = cellfun(@(x) {findfile(x,start_path)}, target_file);
    return
end

target_file = strtrim(target_file);
% We expect one, non-empty file name (i.e., 1xn char array)
if isempty(target_file) || size(target_file, 1) ~= 1
    error('findfile:badfilename', 'Ill-formed target file name')
end

suffix = '.nex';
if ~strncmpi(fliplr(target_file), fliplr(suffix), length(suffix))
    target_file = [target_file suffix];
end

persistent saved_match_path

% handle the trivial case if we've been handed a full/relative file path already
% and also the case where there's an erroneous path in front of the target filename
[folder_path,test_file_name,file_ext] = fileparts(target_file);
if exist(target_file, 'file')
    found_file = which(target_file);
    if isempty(found_file) % the file exists but 'which' doesn't see it
        % it must be relative to the pwd then
        old_dir = cd(folder_path);
        found_file = fullfile(pwd, [test_file_name file_ext]);
        saved_match_path = pwd;
        cd(old_dir);
    else
        saved_match_path = fileparts(found_file);
    end
    return
else
    target_file = [test_file_name file_ext];
end

path_guesses = {};
if ~isempty(saved_match_path)
    path_guesses = [path_guesses saved_match_path];
end

% try the user-defined path next
path_guesses = [path_guesses start_path];

% try the NexFiles path
path_guesses = [path_guesses nexfilepath];

if isempty(path_guesses) % last resort
    user_dir = uigetdir('', 'Select a path as the starting point');
    if ~user_dir, found_file = ''; return; end
    path_guesses = {user_dir};
end

path_counter = 1;

% Let's look through the likely places first
while true
    if path_counter > length(path_guesses), break; end
    
    current_path = path_guesses{path_counter};
    [files,directories] = process_path(current_path);
    found_file = find_target_file(target_file, current_path, files);
    
    if ~isempty(found_file)
        saved_match_path = current_path;
        break
    end
    
    if ~isempty(directories)
        more_paths = arrange_paths(target_file, current_path, ...
            {directories.name});
        path_guesses = [path_guesses(1:path_counter) more_paths ...
            path_guesses(path_counter+1:end)];
        path_guesses = remove_dupe_paths(path_guesses, path_counter);
    end
    path_counter = path_counter + 1;
end

% Remove paths in counter+1:end that are present in 1:counter
function paths = remove_dupe_paths(paths, counter)
[null,removal_idxs] = ismember(paths(1:counter), paths(counter+1:end)); %#ok<ASGLU>
if any(removal_idxs)
    removal_idxs = removal_idxs(removal_idxs > 0) + counter;
    paths(removal_idxs) = [];
end

% Arrange the directory queue in a smart way based on the target file
function dirs = arrange_paths(target_file, curr_path, dirs)
% find the first six consecutive digits in the filename
% note: this assumes the year is the last two digits of the six
digits = regexp(target_file, '\d{6}', 'match', 'once');

if ~isempty(digits)
    year = ['20' digits(end-1:end)];
    year_match = strcmp(year, dirs);
else
    year_match = [];
end

% put paths that match the year or first character in the target filename
% higher in the queue.
folder_match = strncmpi(target_file(1), dirs, 1);
matched_idxs = [find(year_match) find(folder_match)];
unmatched_idxs = setdiff(1:length(dirs), matched_idxs);
dirs = strcat(curr_path, filesep, dirs([matched_idxs unmatched_idxs]));

function rslt = find_target_file(target_file, curr_path, files)
rslt = '';
if ~isempty(files)
    file_match_idx = find(strcmpi(target_file, {files.name}), 1);
    if ~isempty(file_match_idx)
        rslt = fullfile(curr_path, files(file_match_idx).name);
    end
end

function [files,directories] = process_path(in_path)
dir_ignore_list = {'.' '..' '.svn' '.git' 'nexfilelists' ...
    'Text Files' 'BatchData' 'text files' 'Batch Data And Text Files'};

dirstruct = dir(in_path);
is_dir = [dirstruct.isdir];
files = dirstruct(~is_dir);
directories = dirstruct(is_dir);
directories = directories(~ismember({directories.name}, dir_ignore_list));
