function whichmon = FindMon(monspd)
% function whichmon = FindMon(monspd)
%
% Takes an nx3 matrix of monitor SPDs and looks through all of the
% calibration structures to find the one that matches. This will help up
% track down which monitor any given nex file came from (as long as the
% monitor spectra were dropped).

persistent data;
data = [];
dirpath = [fileparts(which('TestCal')),filesep,'Monitor data'];
d = dir(dirpath);
for i = 1:length(d)
    if d(i).name(1) == '.'
        continue
    end
    if d(i).isdir ~= 1
        continue
    end
    if strcmp(d(i).name, 'Rig1') | strcmp(d(i).name, 'Rig2')
        continue
    end
    mon_name = d(i).name;
    dircontents = dir([dirpath,filesep,d(i).name]);
    dataidx = length(data)+1;
    data(dataidx).name = mon_name;
    data(dataidx).spds = {};
    % Could make this recursive
    for j = 1:length(dircontents) % looping over the contents of the subfolders
        if ~dircontents(j).isdir
            filecontents = load([dirpath,filesep,d(i).name,filesep,dircontents(j).name]);
            if isfield(filecontents,'cals')
                cals = filecontents.cals;
                for k = 1:length(cals)
                    data(dataidx).spds{length(data(dataidx).spds)+1}=cals{k}.P_device;
                end
            end
        end
    end
end

whichmon = 'Cannot find';
% Now doing the checking
for i = 1:length(data)
    for j = 1:length(data(i).spds)
        tmp = data(i).spds{j};
        if all(size(tmp) == size(monspd))
            if mean(mean((monspd-tmp).^2)) == 0
                whichmon = data(i).name;
            end
        end
    end
end
end