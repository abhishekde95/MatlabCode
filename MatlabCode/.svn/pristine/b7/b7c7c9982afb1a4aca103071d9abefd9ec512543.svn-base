function [fnames, spikenums] = fnamesFromTxt2(txtFile)
%
%   EXAMPLE: [fnames, spikenums] = fnamesFromTxt2(txtFile)
%
% opens a .txt file and returns each line of the text file to a row of a
% matrix.  This is helpful when batch procesing data from a list of data
% files contained in a text file. Each row should have the same number of
% character elements. 'txtFile' should be the entire path to the txt file.
% 'spikeID' is a flag, which when set, returns the variable 'spikeIdx'.
% spikeIdx is a vector of unit numbers to consider for each file.
%
% This program ignores any text following a percent sign so that the user
% can include comments in the text file. file names and spike IDs should be
% tab delimited in the text file.
%
% CAH 06/08
%
% Modified by GDLH to return cell arrays (which may contain more than one
% file name per cell) and getting rid of the second input argument.
%
% GDLH 03/10/09

fnames = {};
spikenums = [];
% GDLH added 5/1/09
if nargin == 0 || isempty(txtFile)
    [filename, pathname] = uigetfile(nexfilepath('nexfilelists','*.txt'), ...
        'Select a list of files');
    if isequal(filename,0), return; end
    txtFile = fullfile(pathname, filename);
end

fid = fopen(txtFile);
obj = onCleanup(@()fclose(fid));
cellcounter = 1;
if fid > 0
    while ~feof(fid)
        line = fgetl(fid);
        
        %skip over the commented lines and carriage returns
        if numel(line) > 0
            %strip out the comments
            if line(1) == '%'
                continue
            end
            
            while ischar(line)
                %strip off any leading space
                line = strtrim(line);
                firstspace = find(isspace(line),1);
                if ~isempty(firstspace)
                    tmp = line(1:firstspace-1);
                else
                    tmp = line(1:end);
                end
                if isempty(tmp) || tmp(1) == '%' % stripping out comments at the end of a line
                    line = [];
                    continue;
                end
                if length(tmp) > 1 % It's a filename
                    if length(fnames) < cellcounter
                        fnames{cellcounter,1} = {tmp};
                    else
                        fnames{cellcounter,1} = [fnames{cellcounter,1}, {tmp}];
                    end
                else
                    spikenums(cellcounter,1) = str2num(tmp);
                end
                line(1:length(tmp)) = [];
            end
            cellcounter = cellcounter + 1;
        end
    end
else
    error('Unable to open the specified file: \n %s \n', txtFile);
end
