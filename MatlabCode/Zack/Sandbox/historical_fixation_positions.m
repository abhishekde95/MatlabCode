%HISTORICAL_FIXATION_POSITIONS   Plot a subject's fixation positions over time
%   HISTORICAL_FIXATION_POSITIONS(SUBJECT, START, END, LIMITS, SAVE_PLOTS) loads
%   SUBJECT's nex files between START and END dates and plots the median eye
%   position for each day. If SAVE_PLOTS is true then this function saves the
%   plots in the current directory. The save feature requires `export_fig` to be
%   in the path. SUBJECT is the identifying string at the beginning of each nex
%   file name (e.g., 'S' for Sedna). START and END are dates represented in
%   `datenums` format (see help `datenums`). If either START or END are empty,
%   this function defaults to the beginning of time date and the present date,
%   respectively. LIMITS is a 4-element vector with minX, maxX, minY, and maxY.
%
%   See also DATENUM
%
% Zack L.-B. 2013
function historical_fixation_positions(subject_id, date_start, date_end, limits, save_plots)
if nargin == 0, return; end

if nargin < 5 || isempty(save_plots)
    save_plots = false;
end
if nargin < 4 || isempty(limits)
    limits = [-2 2 -2 2];
end
if nargin < 3 || isempty(date_end)
    date_end = floor(now);
end
if nargin < 2 || isempty(date_start)
    date_start = 0;
end

% get all files in `nexfilepath` with this pattern; get filename tokens
file_paths = getAllFiles(nexfilepath, [upper(subject_id) '\d{9}\.nex$']); % this filter excludes filenames with _# or .#
[valid_filename, filename_tokens] = isvalidnexfilename(file_paths);
file_paths = file_paths(valid_filename);
filename_tokens = filename_tokens(valid_filename);

% get rid of duplicate filenames
[~,filenames] = cellfun(@fileparts, file_paths, 'unif', 0);
[~,unique_idxs] = unique(filenames, 'stable');
file_paths = file_paths(unique_idxs);
filename_tokens = filename_tokens(unique_idxs);

% convert tokens to datenums; remove files outside the date range
datenums = datenums_from_tokens(filename_tokens);
valid_dates = datenums >= date_start & datenums <= date_end;
file_paths = file_paths(valid_dates);
datenums = datenums(valid_dates);

% sort the files chronologically
[~,sorted_pos] = sort(datenums);
datenums = datenums(sorted_pos);
file_paths = file_paths(sorted_pos);

% stack all file names with the same date into a cell array
[datenums, first_uniq_idx] = unique(datenums);
dupe_idxs = arrayfun(@colon, first_uniq_idx, ...
    [first_uniq_idx(2:end); length(file_paths)]-1, 'unif', 0);
file_paths_merged = cellfun(@(x) {file_paths(x)}, dupe_idxs);

plot_fixation_positions(file_paths_merged, cellstr(datestr(datenums, 26)), limits, save_plots);

function datenums = datenums_from_tokens(tokens)
iMONTH = 2; iDAY = 3; iYEAR = 4; % from isvalidnexfilename
dates = cellfun(@(x) {[str2double(['20' x{iYEAR}]) ...
    str2double(x{iMONTH}) str2double(x{iDAY})]}, tokens);
datenums = datenum(cell2mat(dates));
