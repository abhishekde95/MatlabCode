% This script should open txt files from iModel and put them into GLMS
% format.
clear all
close all
global GLMP DN


filename = 'zz.txt';
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/iModel/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\iModel\';
end

try
    load([library 'iModelDN.mat'])
    
catch
    
    fID = fopen([library filename]);
    tline = fgetl(fID);
    nframes = 292;
    lumcc = .5;
    colcc = .09;
    nstix = 11;
    
    % Set up DN structure
    while ischar(tline)
        temp = textscan(tline,'%s');
        if strcmp(temp{1}(1),'Channels')
            nchannels = str2double(tline(10:end));
            for c = 1:nchannels
                DN{c}.datafile = ['iModel ' num2str(c)];
                DN{c}.framerate = 75;
            end
            
        elseif strcmp(temp{1}(1),'T')
            t = str2double(temp{1}(2));
            seed = str2double(temp{1}(3));
            for c = 1:nchannels
                tline = fgetl(fID);
                DN{c}.seed(t,1) = seed;
                DN{c}.nFrames(t,1) = nframes;
                DN{c}.NStixGrid(t,1) = nstix;
                spiketimes = str2num(tline(3:end))';
                DN{c}.spiketimes{t} = spiketimes(2:end);
                DN{c}.lumCC(t,1) = lumcc;
                DN{c}.colCC(t,1) = colcc;
            end
        end
        tline = fgetl(fID);
    end
    
    
    % STA for all DN stimuli
    nframesback = 10;
    msperframe = 1000/75;
    colordirlms = [1 1 0; -1 -1 0; 1 -1 0; -1 1 0];
    idxnames = {'pLpML','mLmML','pLmML','mLpML'};
    fieldnames = {'pLpM','mLmM','pLmM','mLpM'};
    for c = 1:nchannels
        for f = 1:numel(fieldnames)
            DN{c}.stim{1}.stats.(fieldnames{f}).STS = zeros(DN{c}.NStixGrid(1).^2*3,nframesback);
            DN{c}.stim{1}.stats.(fieldnames{f}).nspikes = 0;
        end
        for t = 1:numel(DN{c}.seed)
            randnums = getEJrandnums(DN{c}.NStixGrid(t)^2*DN{c}.nFrames(t), DN{c}.seed(t));
            randnums = reshape(randnums, [DN{c}.NStixGrid(t)^2, DN{c}.nFrames(t)]);
            randnums = mod(randnums, size(colordirlms,1))+1;
            randnums_orig = colordirlms(randnums,:);
            %Seperate +/- lum/chrom conditions...
            pLpML = find(randnums_orig(:,1) > 0 & randnums_orig(:,2) > 0);
            mLmML = find(randnums_orig(:,1) < 0 & randnums_orig(:,2) < 0);
            pLmML = find(randnums_orig(:,1) > 0 & randnums_orig(:,2) < 0);
            mLpML = find(randnums_orig(:,1) < 0 & randnums_orig(:,2) > 0);
            for f = 1:numel(fieldnames)
                randnums = zeros(size(randnums_orig));
                randnums(eval(idxnames{f}),:) = randnums_orig(eval(idxnames{f}),:);
                randnums = reshape(randnums, [DN{c}.NStixGrid(t).^2, DN{c}.nFrames(t), 3]);
                randnums = permute(randnums,[1 3 2]);
                randnums = reshape(randnums, [DN{c}.NStixGrid(t)^2*3, DN{c}.nFrames(t)]);
                % Each column should be Ls followed by Ms followed by Ss for each frame.
                if any(DN{c}.spiketimes{t}) && DN{c}.nFrames(t) ~= 0
                    frametimes = linspace(0, DN{c}.nFrames(t) * msperframe, DN{c}.nFrames(t))+(msperframe/2)';
                    sptimes = DN{c}.spiketimes{t};
                    sptimes(sptimes<0) = [];
                    sptimes(sptimes>frametimes(end)) = [];
                    [n,x] = hist(sptimes, frametimes);
                    STCOVmex('init',{DN{c}.NStixGrid(t).^2,3,nframesback})
                    STCOVmex(randnums(:), n);
                    out = STCOVmex('return');
                    DN{c}.stim{1}.stats.(fieldnames{f}).STS = DN{c}.stim{1}.stats.(fieldnames{f}).STS+out{1};
                    DN{c}.stim{1}.stats.(fieldnames{f}).nspikes = DN{c}.stim{1}.stats.(fieldnames{f}).nspikes+out{3};
                    DN{c}.stim{1}.stats.(fieldnames{f}).STA = DN{c}.stim{1}.stats.(fieldnames{f}).STS ./ DN{c}.stim{1}.stats.(fieldnames{f}).nspikes;
                end
            end
        end
    end
    
    
    iModelDN = DN;
    save([library 'iModelDN'],'iModelDN')
    
end



DN = iModelDN{16};
GLMP = [];

GLMS_AnalysisGUI(GLMP,DN);


