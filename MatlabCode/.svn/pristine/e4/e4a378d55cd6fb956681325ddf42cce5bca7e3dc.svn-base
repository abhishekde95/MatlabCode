%STROCAT   Concatenate STRO files
%   OUTSTRO = STROCAT(stro1, stro2, ...) concatenates any number of STRO
%   structures into a single OUTSTRO. The first structure's (`stro1`)
%   experimental parameters are used as a template to assure subsequent
%   structures came from the same paradigm.
%
%   OUTSTRO = STROCAT(STROS) does the same as above but accepts a cell array
%   STROS containing any number of STRO structures.
%
% CAH 9.10.08

%#ok<*AGROW>
function out = strocat(varargin)

varargin = flatten(varargin);
varargin(cellfun('isempty', varargin)) = [];
nFiles = length(varargin);
TIMEMARGIN = 1; % assuming 1 sec dead time between files

if nFiles == 0
    out = [];
    return
end

% use the first file as the standard for comparison
out = varargin{1};
if nFiles == 1, return; end

id = out.sum.paradigmID;
fld_names = fieldnames(out.sum.exptParams);
[paths,filenames] = cellfun(@(x) fileparts(x.sum.fileName), varargin, 'unif', 0);

[filenames,uniq_location] = unique(filenames, 'stable');
if numel(uniq_location) < nFiles
    filenames = filenames(uniq_location);
    warning(['Found STROs with matching filenames. ' ...
        'There may be duplicate STROs in the concatenated result!']);
end
out.sum.fileName = [paths{1} filesep integrateFilenames(upper(filenames)) '.nex'];

% iterate over the subsequent files and concatenate
for ii = 2:nFiles
    tmp = varargin{ii};
    
    % break if the files don't have the same id
    if tmp.sum.paradigmID ~= id
        error('One or more of the files came from different paradigms');
    end
    
    % make sure that the experimental parameters are the same
    for b = 1:length(fld_names)
        if ~isequaln(out.sum.exptParams.(fld_names{b}), tmp.sum.exptParams.(fld_names{b}))
            warning('\n*** File <%d> has an inconsistent <%s> parameter ***', ii, fld_names{b})
        end
    end
    % If you've made it to this point everything should check out, so
    % concatenate away.
    
    % Finding the time of the last event and then adding on TIMEMARGIN 
    timefields = strcmp(out.sum.trialFields(2,:),'time');
    
    LASTTIME = nanmax(out.trial(end,timefields));
    STARTTIMECONST = LASTTIME+TIMEMARGIN;
    tmp.trial(:,timefields) = tmp.trial(:,timefields)+STARTTIMECONST;
    out.trial = [out.trial; tmp.trial];
    [out.ras, out.sum.rasterCells] = integrateRaster(out.ras, out.sum.rasterCells, tmp.ras, tmp.sum.rasterCells,STARTTIMECONST);
    out.other = [out.other; tmp.other];
    out.sum.absTrialNum = [out.sum.absTrialNum; out.sum.absTrialNum(end)+tmp.sum.absTrialNum];
    
   % Old Charlie stuff
   % if any([210 212] == tmp.sum.paradigmID)
   %     out.LMS = [out.LMS; tmp.LMS];
   % end
end

% Concatenates filenames. `filenames{1}` serves as the root and the rest are
% used to modify it. What is returned is: "filenames{1}&<non-matching part of
% filenames{2}>&<and so on>".
function combined_fname = integrateFilenames(filenames)
template = filenames{1};
filenames = filenames(2:end);
append_strs = cell(1, length(filenames));
for n_chars = length(template):-1:1
    unmatched_pos = ~strncmp(template, filenames, n_chars);
    if ~any(unmatched_pos), break; end
    append_strs(unmatched_pos) = cellfun(@(x) {x(n_chars:end)}, filenames(unmatched_pos));
end
combined_fname = [template sprintf('&%s', append_strs{:})];

% Concatenates raster matrices in the event that they are unequally sized. This
% happens when a neuron is picked up or lost between data files.
function [ras, cells] = integrateRaster(ras1, cells1, ras2, cells2, STARTTIMECONST)

spiketimefields = strncmp(cells1,'sig',3);
anlgStartTimefield = strcmp(cells1,'anlgStartTime');
timefields = spiketimefields | anlgStartTimefield;

ras = {};
cells = {};
% First iterating over the columns of the first file
for i = 1:length(cells1)
    cells(i) = cells1(i);
    idx = find(ismember(cells2, cells1(i)));
    if ~isempty(idx)
        if timefields(i)
            ras(:,i) = [ras1(:,i); cellfun(@(x)(x+STARTTIMECONST), ras2(:,idx), 'UniformOutput',0)];
        else
            ras(:,i) = [ras1(:,i); ras2(:,idx)];            
        end
    else
        ras(:,i) = [ras1(:,1); cell(size(ras2,1),1)];
    end
end
% Now taking care of any columns in the second file that didn't appear in the
% first file

spiketimefields = strncmp(cells2,'sig',3);
anlgStartTimefield = strcmp(cells2,'anlgStartTime');
timefields = spiketimefields | anlgStartTimefield;
for i = 1:length(cells2)
    idx = find(ismember(cells2(i), cells1), 1);
    if isempty(idx)
        tmpras = ras(:,i:end);
        tmpcell = cells(i:end);
        if timefields(i)
            ras(:,i) = [cell(size(ras1,1),1); cellfun(@(x)(x+STARTTIMECONST), ras2(:,i), 'UniformOutput',0)];
        else
            ras(:,i) = [cell(size(ras1,1),1); ras2(:,idx)];
        end
        cells(i) = cells2(i);
        ras(:,i+1:i+size(tmpras,2)) = tmpras;
        cells(i+1:i+length(tmpcell)) = tmpcell;
    end
end

% http://www.mathworks.com/matlabcentral/fileexchange/27009-flatten-nested-cell-arrays
% Manu Raghavan, modified by Zack L-B
function C = flatten(A)
C = {};
for i = 1:numel(A)
    if ~iscell(A{i})
        C = cat(2, C, A(i));
    else
        C = cat(2, C, flatten(A{i}));
    end
end
