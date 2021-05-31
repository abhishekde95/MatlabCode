function [filenames, mode, spikeidx] = get_WNthreshdata_from_csvfile(csv_filename, WNthresh_mode)
% This function extracts the filenames, mode and spike index from a given csv file
%
% csv_filename <string>: Name of the csv file
% WNthresh_mode <string> Name of the Whitenoise thresh subunit selection mode 

% Reading the table
C = readtable(csv_filename);

% Extracting the filenames, mode and spike idx
mode = C.Var5;
filenames = C.Var2(strcmp(mode, WNthresh_mode));
spikeidx = C.Var3(strcmp(mode, WNthresh_mode));
mode = C.Var5(strcmp(mode, WNthresh_mode));

end

