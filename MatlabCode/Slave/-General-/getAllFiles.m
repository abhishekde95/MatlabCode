% from http://stackoverflow.com/a/2654459
function fileList = getAllFiles(dirName, pattern)
if nargin < 2 || isempty(pattern), pattern = []; end
dirData = dir(dirName);
dirIndex = [dirData.isdir];
fileList = {dirData(~dirIndex).name}';
if ~isempty(fileList)
    fileList = cellfun(@(x) {fullfile(dirName,x)}, fileList);
    if ~isempty(pattern)
        matchstart = regexp(fileList, pattern);
        fileList = fileList(~cellfun(@isempty, matchstart));
    end
end
subDirs = {dirData(dirIndex).name};
for iDir = 3:numel(subDirs)
    nextDir = fullfile(dirName,subDirs{iDir});
    fileList = [fileList; getAllFiles(nextDir, pattern)]; %#ok<AGROW>
end
