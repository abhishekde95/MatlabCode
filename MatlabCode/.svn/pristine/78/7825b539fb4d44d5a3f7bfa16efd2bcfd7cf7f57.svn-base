%called by at least nex2db, if not other functions. 
%Starting at directory dirName, moves through all subdirectories and stores/outputs all filenames of type fileType.
%Originally created for storing all nex files, but can store any arbitrary fileType.
%Also can store/output the full directory tree.
function [fileList, dirList] = scrapeFilesFromDir(dirName, fileType)
persistent FullFileList
persistent validDirectoryList
cd(dirName);
cdir = dir;
fullFileInfo = struct2cell(dir(fileType));
fileNames = fullFileInfo(1, :);
fileNamesWithDir = fullfile(dirName, fileNames);
FullFileList = unique([FullFileList, fileNamesWithDir]);
validDirectoryList{end+1} = {dirName};
fullDirNames = {cdir([cdir.isdir]).name};
validDirIdx = ~ismember(fullDirNames,{'.','..'});
validDirNames = fullDirNames(validDirIdx);
if ~isempty(validDirNames)
    for d = 1:length(validDirNames)
        nextDir = fullfile(dirName, validDirNames(d));
        [fileList, dirList] = scrapeFilesFromDir(nextDir{1}, fileType);        
    end
else
    FullFileList = [FullFileList, fileNamesWithDir];
    fileList = FullFileList;
    dirList = validDirectoryList;
end
keyboard;
end