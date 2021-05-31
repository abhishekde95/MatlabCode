% A new script to check if I can used the names from the CSV files and get
% the contents of the nex2stro files

close all;
clearvars;

C = readtable('WNthresh.csv');
mode = C.Var5; 
filenames = C.Var2(strcmp(mode,"subunit"));
spikidx = C.Var3(strcmp(mode,"subunit"));

for ii = 1: numel(filenames)
    fileofinterest = filenames{ii};
    stro = nex2stro(findfile(fileofinterest,'/Users/abhishekde/Google Drive/Data_Physiology'));
end
