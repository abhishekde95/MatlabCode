%%
%SMurray paradigm analysis
%plots based on analysis of a specified LFP frequency
clear all; clc



%%
%choose analysis parameters


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%choose which files to analyze
num_files = 3;
%'1' to analyze selected file
%'2' to analyze all files recorded on the same day as the selected file
%'3' to analyze all files in the folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%choose LFP frequency to analyze
freq = 40;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%choose LFP analysis method
analysis_method = 1;
%'1' for multitaper analysis
%'2' for filter analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%choose array channels to analyze
%channels = [3:16 19:32]; %Freya
channels = [1:13 17:29]; %Apollo 1st round of data
%channels = [1:32]; %Apollo 2nd round of data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%choose analysis window beginning and end times
aw = [0.05; 0.05]; %[seconds post-stim onset; seconds post-stim offset]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%
%store LFP data from file/s in struct 'days'

if num_files == 1
    [fname, pathname] = uigetfile('*.nex', 'Select a NeuroExplorer file');
    days{1}.dayname = fname(1:7);
    stro = nex2stro(strcat(pathname, fname));
    disp(stro);
    bkgnd = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd'));
    if bkgnd == 2 %gray background
        [LFPsamprate, nearG, farG, uniqueors, nearBaseG, farBaseG] = SMurray_freq_LFP(stro, aw);
        days{1}.LFPsamprate = LFPsamprate;
        days{1}.uniqueors = uniqueors;
        days{1}.nearG = nearG;
        days{1}.farG = farG;
        days{1}.nearBaseG = nearBaseG;
        days{1}.farBaseG = farBaseG;
    elseif bkgnd == 0 %corridor background
        [LFPsamprate, nearC, farC, uniqueors, nearBaseC, farBaseC] = SMurray_freq_LFP(stro, aw);
        days{1}.LFPsamprate = LFPsamprate;
        days{1}.uniqueors = uniqueors;
        days{1}.nearC = nearC;
        days{1}.farC = farC;
        days{1}.nearBaseC = nearBaseC;
        days{1}.farBaseC = farBaseC;
    elseif bkgnd == 1 %brick background
        [LFPsamprate, nearB, farB, uniqueors, nearBaseB, farBaseB] = SMurray_freq_LFP(stro, aw);
        days{1}.LFPsamprate = LFPsamprate;
        days{1}.uniqueors = uniqueors;
        days{1}.nearB = nearB;
        days{1}.farB = farB;
        days{1}.nearBaseB = nearBaseB;
        days{1}.farBaseB = farBaseB;
    end
elseif num_files == 2
    [fname, pathname] = uigetfile('*.nex', 'Select a NeuroExplorer file');
    days{1}.dayname = fname(1:7);
    day = fname([1:8 11:14]);
    allfiles = dir(pathname);
    rowfiles = 1;
    for i = 1:size(allfiles,1)
        if size(allfiles(i).name,2) == 14
            if allfiles(i).name([1:8 11:14]) == day
                allblocks(rowfiles,:) = allfiles(i).name;
                rowfiles = rowfiles + 1;
            end
        end
    end
    rowG = 1;
    rowC = 1;
    rowB = 1;
    for i = 1:size(allblocks,1)
        stro = nex2stro(strcat(pathname, allblocks(i,:)));
        bkgnd = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd'));
        %combine trials from all files, according to background
        if bkgnd == 2 %gray background
            stroG.sum = stro.sum;
            stroG.trial(rowG:rowG+length(stro.trial)-1,:) = stro.trial(1:length(stro.trial),:);
            stroG.ras(rowG:rowG+length(stro.ras)-1,:) = stro.ras(1:length(stro.ras),:);
            rowG = rowG+length(stro.trial);
        elseif bkgnd == 0 %corridor background
            stroC.sum = stro.sum;
            stroC.trial(rowC:rowC+length(stro.trial)-1,:) = stro.trial(1:length(stro.trial),:);
            stroC.ras(rowC:rowC+length(stro.ras)-1,:) = stro.ras(1:length(stro.ras),:);
            rowC = rowC+length(stro.trial);
        elseif bkgnd == 1 %brick background
            stroB.sum = stro.sum;
            stroB.trial(rowB:rowB+length(stro.trial)-1,:) = stro.trial(1:length(stro.trial),:);
            stroB.ras(rowB:rowB+length(stro.ras)-1,:) = stro.ras(1:length(stro.ras),:);
            rowB = rowB+length(stro.trial);
        end
    end
    if exist('stroG', 'var')
        [LFPsamprate, nearG, farG, uniqueors, nearBaseG, farBaseG] = SMurray_freq_LFP(stroG, aw);
        days{1}.LFPsamprate = LFPsamprate;
        days{1}.uniqueors = uniqueors;
        days{1}.nearG = nearG;
        days{1}.farG = farG;
        days{1}.nearBaseG = nearBaseG;
        days{1}.farBaseG = farBaseG;
    end
    if exist('stroC', 'var')
        [~,           nearC, farC, ~,         nearBaseC, farBaseC] = SMurray_freq_LFP(stroC, aw);
        days{1}.LFPsamprate = LFPsamprate;
        days{1}.uniqueors = uniqueors;
        days{1}.nearC = nearC;
        days{1}.farC = farC;
        days{1}.nearBaseC = nearBaseC;
        days{1}.farBaseC = farBaseC;
    end
    if exist('stroB', 'var')
        [~,           nearB, farB, ~,         nearBaseB, farBaseB] = SMurray_freq_LFP(stroB, aw);
        days{1}.LFPsamprate = LFPsamprate;
        days{1}.uniqueors = uniqueors;
        days{1}.nearB = nearB;
        days{1}.farB = farB;
        days{1}.nearBaseB = nearBaseB;
        days{1}.farBaseB = farBaseB;
    end
    
    
elseif num_files == 3
    [~, pathname] = uigetfile('*.nex', 'Select a NeuroExplorer file');
    folder = pathname;
    allfiles = dir(folder);
    for i = 1:length(allfiles)
        if allfiles(i).name(1) == '.'
            continue
        else
            break
        end
    end
    startfile = i;
    allnames = NaN(length(allfiles)-startfile+1, 1);
    for i = startfile:length(allfiles)
        allnames(i-startfile+1, 1:8) = allfiles(i).name(1:8);
    end
    daynames = unique(allnames, 'rows');
    for i = 1:size(daynames,1)
        day = char(daynames(i,:));
        days{i}.dayname = day(1:7);
        clear stroG stroC stroB
        rowG = 1;
        rowC = 1;
        rowB = 1;
        for j = startfile:size(allfiles,1)
            if allfiles(j).name(1:8) == day
                stro = nex2stro(strcat(folder, allfiles(j).name));
                bkgnd = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd'));
                %combine trials from all files, according to background
                if bkgnd == 2 %gray background
                    stroG.sum = stro.sum;
                    stroG.trial(rowG:rowG+length(stro.trial)-1,:) = stro.trial(1:length(stro.trial),:);
                    stroG.ras(rowG:rowG+length(stro.ras)-1,:) = stro.ras(1:length(stro.ras),:);
                    rowG = rowG+length(stro.trial);
                elseif bkgnd == 0 %corridor background
                    stroC.sum = stro.sum;
                    stroC.trial(rowC:rowC+length(stro.trial)-1,:) = stro.trial(1:length(stro.trial),:);
                    stroC.ras(rowC:rowC+length(stro.ras)-1,:) = stro.ras(1:length(stro.ras),:);
                    rowC = rowC+length(stro.trial);
                elseif bkgnd == 1 %brick background
                    stroB.sum = stro.sum;
                    stroB.trial(rowB:rowB+length(stro.trial)-1,:) = stro.trial(1:length(stro.trial),:);
                    stroB.ras(rowB:rowB+length(stro.ras)-1,:) = stro.ras(1:length(stro.ras),:);
                    rowB = rowB+length(stro.trial);
                end
            end
        end
        if exist('stroG', 'var')
            [LFPsamprate, nearG, farG, uniqueors, nearBaseG, farBaseG] = SMurray_freq_LFP(stroG, aw);
            days{i}.LFPsamprate = LFPsamprate;
            days{i}.uniqueors = uniqueors;
            days{i}.nearG = nearG;
            days{i}.farG = farG;
            days{i}.nearBaseG = nearBaseG;
            days{i}.farBaseG = farBaseG;
        end
        if exist('stroC', 'var')
            [~,           nearC, farC, ~,         nearBaseC, farBaseC] = SMurray_freq_LFP(stroC, aw);
            days{i}.LFPsamprate = LFPsamprate;
            days{i}.uniqueors = uniqueors;
            days{i}.nearC = nearC;
            days{i}.farC = farC;
            days{i}.nearBaseC = nearBaseC;
            days{i}.farBaseC = farBaseC;
        end
        if exist('stroB', 'var')
            [~,           nearB, farB, ~,         nearBaseB, farBaseB] = SMurray_freq_LFP(stroB, aw);
            days{i}.LFPsamprate = LFPsamprate;
            days{i}.uniqueors = uniqueors;
            days{i}.nearB = nearB;
            days{i}.farB = farB;
            days{i}.nearBaseB = nearBaseB;
            days{i}.farBaseB = farBaseB;
        end
    end
end






    
    
    
    
    

%%
%set LFP analysis parameters

if analysis_method == 1 %multitaper analysis (Chronux)
    %Preliminary parameter settings
    params.Fs=days{1}.LFPsamprate; % sampling frequency
    params.fpass=freq; % frequency of interest
    params.tapers=[3 5]; % tapers
    params.trialave=0; % 1=average over trials
    %params.err=0; % 1 = error computation
elseif analysis_method == 2 %filter analysis
    %filters
    cyclespersample = freq./LFPsamprate;
    ncyclesper6sigma = 5;  % Controls the bandwidth (higher numbers = tighter BW)
    nsamplesper6sigma = ceil(ncyclesper6sigma/cyclespersample);
    ncycles = nsamplesper6sigma*cyclespersample;
    filtkernel1 = normpdf(linspace(-3,3,nsamplesper6sigma),0,1).*cos(linspace(0,2*pi*ncycles,nsamplesper6sigma));
    filtkernel2 = normpdf(linspace(-3,3,nsamplesper6sigma),0,1).*sin(linspace(0,2*pi*ncycles,nsamplesper6sigma));
    filtkernel1 = filtkernel1./norm(filtkernel1);
    filtkernel2 = filtkernel2./norm(filtkernel2);
end



%%
%LFP power (frequency specified above) for individual channels

for i = 1:length(days)
    if isfield(days{i}, 'nearG')
        if analysis_method == 1 %multitaper analysis (Chronux)
            [LFPmeanG, LFPsemG, BASEmeanG, BASEsemG, RESPmeanG, RESPsemG] = SMurray_freq_spectrum(days{i}.nearG, days{i}.farG, ...
                                                                               days{i}.nearBaseG, days{i}.farBaseG, params);
        elseif analysis_method == 2 %filter analysis
            [LFPmeanG, LFPsemG, BASEmeanG, BASEsemG, RESPmeanG, RESPsemG] = SMurray_freq_filter(days{i}.nearG, days{i}.farG, ...
                                                                            days{i}.nearBaseG, days{i}.farBaseG, filtkernel1, filtkernel2);
        end
        days{i}.LFPmeanG = LFPmeanG;
        days{i}.LFPsemG = LFPsemG;
        days{i}.BASEmeanG = BASEmeanG;
        days{i}.BASEsemG = BASEsemG;
        days{i}.RESPmeanG = RESPmeanG;
        days{i}.RESPsemG = RESPsemG;
    end
    
    if isfield(days{i}, 'nearC')
        if analysis_method == 1 %multitaper analysis (Chronux)
            [LFPmeanC, LFPsemC, BASEmeanC, BASEsemC, RESPmeanC, RESPsemC] = SMurray_freq_spectrum(days{i}.nearC, days{i}.farC, ...
                                                                               days{i}.nearBaseC, days{i}.farBaseC, params);
        elseif analysis_method == 2 %filter analysis
            [LFPmeanC, LFPsemC, BASEmeanC, BASEsemC, RESPmeanC, RESPsemC] = SMurray_freq_filter(days{i}.nearC, days{i}.farC, ...
                                                                            days{i}.nearBaseC, days{i}.farBaseC, filtkernel1, filtkernel2);
        end
        days{i}.LFPmeanC = LFPmeanC;
        days{i}.LFPsemC = LFPsemC;
        days{i}.BASEmeanC = BASEmeanC;
        days{i}.BASEsemC = BASEsemC;
        days{i}.RESPmeanC = RESPmeanC;
        days{i}.RESPsemC = RESPsemC;
    end
    
    if isfield(days{i}, 'nearB')
        if analysis_method == 1 %multitaper analysis (Chronux)
            [LFPmeanB, LFPsemB, BASEmeanB, BASEsemB, RESPmeanB, RESPsemB] = SMurray_freq_spectrum(days{i}.nearB, days{i}.farB, ...
                                                                               days{i}.nearBaseB, days{i}.farBaseB, params);
        elseif analysis_method == 2 %filter analysis
            [LFPmeanB, LFPsemB, BASEmeanB, BASEsemB, RESPmeanB, RESPsemB] = SMurray_freq_filter(days{i}.nearB, days{i}.farB, ...
                                                                            days{i}.nearBaseB, days{i}.farBaseB, filtkernel1, filtkernel2);
        end
        days{i}.LFPmeanB = LFPmeanB;
        days{i}.LFPsemB = LFPsemB;
        days{i}.BASEmeanB = BASEmeanB;
        days{i}.BASEsemB = BASEsemB;
        days{i}.RESPmeanB = RESPmeanB;
        days{i}.RESPsemB = RESPsemB;
    end
end



%}
%%
%{
% load 'day' struct if skipping above cells

load days_LFP_FixInMove_Ap
days = days_LFP_FixInMove_Ap;

%}
%%
% plots of individual channels, LFP power

for i = 1:length(days)
    
    %Gray background, LFP power
    if isfield(days{i}, 'LFPmeanG')
        figure;
        numchans = size(days{i}.LFPmeanG,3);
        for j = 1:numchans
            subplot(ceil(sqrt(numchans)),ceil(sqrt(numchans)),j);
            set(gca,'yscale','log');
            hold on;
            errorbar(days{i}.uniqueors, days{i}.LFPmeanG(1,:,j), days{i}.LFPsemG(1,:,j), '.-m') %near
            errorbar(days{i}.uniqueors, days{i}.LFPmeanG(2,:,j), days{i}.LFPsemG(2,:,j), '.-k') %far
            xlabel('Radius (deg)');
            ylabel('Power (dB)');
        end
        set(gcf, 'Name', ['Day ' num2str(days{i}.dayname) ': Gray LFP power'])
        pause
        close(gcf)
    end
    
    %Gray background, LFP response
    if isfield(days{i}, 'RESPmeanG')
        figure;
        numchans = size(days{i}.RESPmeanG,3);
        for j = 1:numchans
            subplot(ceil(sqrt(numchans)),ceil(sqrt(numchans)),j);
            hold on;
            errorbar(days{i}.uniqueors, days{i}.RESPmeanG(1,:,j), days{i}.RESPsemG(1,:,j), '.-m') %near
            errorbar(days{i}.uniqueors, days{i}.RESPmeanG(2,:,j), days{i}.RESPsemG(2,:,j), '.-k') %far
            xlabel('Radius (deg)');
            ylabel('Response (dB)');
        end
        set(gcf, 'Name', ['Day ' num2str(days{i}.dayname) ': Gray LFP change-from-baseline'])
        pause
        close(gcf)
    end
    
    %Corridor background, LFP power
    if isfield(days{i}, 'LFPmeanC')
        figure;
        numchans = size(days{i}.LFPmeanC,3);
        for j = 1:numchans
            subplot(ceil(sqrt(numchans)),ceil(sqrt(numchans)),j);
            set(gca,'yscale','log');
            hold on;
            errorbar(days{i}.uniqueors, days{i}.LFPmeanC(1,:,j), days{i}.LFPsemC(1,:,j), '.-m') %near
            errorbar(days{i}.uniqueors, days{i}.LFPmeanC(2,:,j), days{i}.LFPsemC(2,:,j), '.-k') %far
            xlabel('Radius (deg)');
            ylabel('Power (dB)');
        end
        set(gcf, 'Name', ['Day ' num2str(days{i}.dayname) ': Corridor LFP power'])
        pause
        close(gcf)
    end
    
    %Corridor background, LFP respone
    if isfield(days{i}, 'RESPmeanC')
        figure;
        numchans = size(days{i}.RESPmeanC,3);
        for j = 1:numchans
            subplot(ceil(sqrt(numchans)),ceil(sqrt(numchans)),j);
            hold on;
            errorbar(days{i}.uniqueors, days{i}.RESPmeanC(1,:,j), days{i}.RESPsemC(1,:,j), '.-m') %near
            errorbar(days{i}.uniqueors, days{i}.RESPmeanC(2,:,j), days{i}.RESPsemC(2,:,j), '.-k') %far
            xlabel('Radius (deg)');
            ylabel('Response (dB)');
        end
        set(gcf, 'Name', ['Day ' num2str(days{i}.dayname) ': Corridor LFP change-from-baseline'])
        pause
        close(gcf)
    end
    
    %Brick background, LFP power
    if isfield(days{i}, 'LFPmeanB')
        figure;
        numchans = size(days{i}.LFPmeanB,3);
        for j = 1:numchans
            subplot(ceil(sqrt(numchans)),ceil(sqrt(numchans)),j);
            set(gca,'yscale','log');
            hold on;
            errorbar(days{i}.uniqueors, days{i}.LFPmeanB(1,:,j), days{i}.LFPsemB(1,:,j), '.-m') %near
            errorbar(days{i}.uniqueors, days{i}.LFPmeanB(2,:,j), days{i}.LFPsemB(2,:,j), '.-k') %far
            xlabel('Radius (deg)');
            ylabel('Power (dB)');
        end
        set(gcf, 'Name', ['Day ' num2str(days{i}.dayname) ': Brick LFP power'])
        pause
        close(gcf)
    end
    
    %Brick background, LFP response
    if isfield(days{i}, 'RESPmeanB')
        figure;
        numchans = size(days{i}.RESPmeanB,3);
        for j = 1:numchans
            subplot(ceil(sqrt(numchans)),ceil(sqrt(numchans)),j);
            hold on;
            errorbar(days{i}.uniqueors, days{i}.RESPmeanB(1,:,j), days{i}.RESPsemB(1,:,j), '.-m') %near
            errorbar(days{i}.uniqueors, days{i}.RESPmeanB(2,:,j), days{i}.RESPsemB(2,:,j), '.-k') %far
            xlabel('Radius (deg)');
            ylabel('Response (dB)');
        end
        set(gcf, 'Name', ['Day ' num2str(days{i}.dayname) ': Brick LFP change-from-baseline'])
        pause
        close(gcf)
    end
end



    

%%
% Mean of LFP power, across channels, per day
% Exclude any signals not passing the LFP power threshold
% Exclusion based on LFP power recorded on gray background

for i = 1:length(days)

    day_avg{1}.uniqueors(i,:) = days{i}.uniqueors';
    day_avg{1}.ex.uniqueors(i,:) = days{i}.uniqueors';
    
    if isfield(days{i}, 'LFPmeanG')
        day_avg{1}.LFPmeanG(:,:,i) = mean(days{i}.LFPmeanG,3);
        day_avg{1}.RESPmeanG(:,:,i) = mean(days{i}.RESPmeanG,3);
        day_avg{1}.LFPsemG(:,:,i)  = std(days{i}.LFPmeanG,0,3)  ./ sqrt(size(days{i}.LFPmeanG,3));
        day_avg{1}.RESPsemG(:,:,i) = std(days{i}.RESPmeanG,0,3) ./ sqrt(size(days{i}.RESPmeanG,3));
        
        %Exclusion threshold for signals of lower power, based on LFP power during Gray background
        avgmean = mean(reshape(day_avg{1}.LFPmeanG(:,:,i), 1, 14)); % day's average LFP power on gray background
        avgsem = mean(reshape(day_avg{1}.LFPsemG(:,:,i), 1, 14)); % day's SEM of LFP power on gray background
        thresh = avgmean - 5*avgsem; % minumum LFP power required
        excluded = find(days{i}.LFPmeanG < thresh); %signals to exclude for this day, based on LFP power on gray background
        
        %Only signals that have passed exclusion threshold
        days{i}.ex.LFPmeanG = days{i}.LFPmeanG;
        days{i}.ex.LFPmeanG(excluded) = NaN;
        days{i}.ex.RESPmeanG = days{i}.RESPmeanG;
        days{i}.ex.RESPmeanG(excluded) = NaN;
        
        %Mean of LFP power, across channels, per day (excluding signals of low power)
        day_avg{1}.ex.LFPmeanG(:,:,i) = nanmean(days{i}.ex.LFPmeanG,3);
        day_avg{1}.ex.RESPmeanG(:,:,i) = nanmean(days{i}.ex.RESPmeanG,3);
        day_avg{1}.ex.LFPsemG(:,:,i)  = nanstd(days{i}.ex.LFPmeanG,0,3)  ./ sqrt(sum(~isnan(days{i}.ex.LFPmeanG),3));   
        day_avg{1}.ex.RESPsemG(:,:,i) = nanstd(days{i}.ex.RESPmeanG,0,3) ./ sqrt(sum(~isnan(days{i}.ex.RESPmeanG),3));
    end
    
    if isfield(days{i}, 'LFPmeanC')
        day_avg{1}.LFPmeanC(:,:,i) = mean(days{i}.LFPmeanC,3);
        day_avg{1}.RESPmeanC(:,:,i) = mean(days{i}.RESPmeanC,3);
        day_avg{1}.LFPsemC(:,:,i)  = std(days{i}.LFPmeanC,0,3)  ./ sqrt(size(days{i}.LFPmeanC,3));
        day_avg{1}.RESPsemC(:,:,i) = std(days{i}.RESPmeanC,0,3) ./ sqrt(size(days{i}.RESPmeanC,3));
        
        %Only signals that have passed exclusion threshold
        days{i}.ex.LFPmeanC = days{i}.LFPmeanC;
        days{i}.ex.LFPmeanC(excluded) = NaN;
        days{i}.ex.RESPmeanC = days{i}.RESPmeanC;
        days{i}.ex.RESPmeanC(excluded) = NaN;
        
        %Mean of LFP power, across channels, per day (excluding signals of low power)
        day_avg{1}.ex.LFPmeanC(:,:,i) = nanmean(days{i}.ex.LFPmeanC,3);
        day_avg{1}.ex.RESPmeanC(:,:,i) = nanmean(days{i}.ex.RESPmeanC,3);
        day_avg{1}.ex.LFPsemC(:,:,i)  = nanstd(days{i}.ex.LFPmeanC,0,3)  ./ sqrt(sum(~isnan(days{i}.ex.LFPmeanC),3));   
        day_avg{1}.ex.RESPsemC(:,:,i) = nanstd(days{i}.ex.RESPmeanC,0,3) ./ sqrt(sum(~isnan(days{i}.ex.RESPmeanC),3));
    end
    
    if isfield(days{i}, 'LFPmeanB')
        day_avg{1}.LFPmeanB(:,:,i) = mean(days{i}.LFPmeanB,3);
        day_avg{1}.RESPmeanB(:,:,i) = mean(days{i}.RESPmeanB,3);
        day_avg{1}.LFPsemB(:,:,i)  = std(days{i}.LFPmeanB,0,3)  ./ sqrt(size(days{i}.LFPmeanB,3));
        day_avg{1}.RESPsemB(:,:,i) = std(days{i}.RESPmeanB,0,3) ./ sqrt(size(days{i}.RESPmeanB,3));
        
        %Only signals that have passed exclusion threshold
        days{i}.ex.LFPmeanB = days{i}.LFPmeanB;
        days{i}.ex.LFPmeanB(excluded) = NaN;
        days{i}.ex.RESPmeanB = days{i}.RESPmeanB;
        days{i}.ex.RESPmeanB(excluded) = NaN;
        
        %Mean of LFP power, across channels, per day (excluding signals of low power)
        day_avg{1}.ex.LFPmeanB(:,:,i) = nanmean(days{i}.ex.LFPmeanB,3);
        day_avg{1}.ex.RESPmeanB(:,:,i) = nanmean(days{i}.ex.RESPmeanB,3);
        day_avg{1}.ex.LFPsemB(:,:,i)  = nanstd(days{i}.ex.LFPmeanB,0,3)  ./ sqrt(sum(~isnan(days{i}.ex.LFPmeanB),3));   
        day_avg{1}.ex.RESPsemB(:,:,i) = nanstd(days{i}.ex.RESPmeanB,0,3) ./ sqrt(sum(~isnan(days{i}.ex.RESPmeanB),3));
    end
    
end




%% 
% Plots of means of LFP power, across channels, per day
% Based on signals that passed exclusion threshold above

for i = 1:size(day_avg{1}.ex.LFPmeanG, 3)
    
    %LFP power
    figure;
    set(gcf, 'Name', ['Day ' num2str(days{i}.dayname) ': Mean of LFP, across channels'])
    %Gray background
    if isfield(day_avg{1}.ex, 'LFPmeanG')
        subplot(2,3,1); hold on;
        errorbar(day_avg{1}.ex.uniqueors(i,:), day_avg{1}.ex.LFPmeanG(1,:,i), day_avg{1}.ex.LFPsemG(1,:,i), '.-m') %near
        errorbar(day_avg{1}.ex.uniqueors(i,:), day_avg{1}.ex.LFPmeanG(2,:,i), day_avg{1}.ex.LFPsemG(2,:,i), '.-k') %far
        xlabel('Radius (deg)');
        ylabel('Power (dB)');
        title('Gray: LFP power')
    end
    %Corridor background
    if isfield(day_avg{1}.ex, 'LFPmeanC')
        if size(day_avg{1}.ex.LFPmeanC,3) >= i 
            subplot(2,3,2); hold on;
            errorbar(day_avg{1}.ex.uniqueors(i,:), day_avg{1}.ex.LFPmeanC(1,:,i), day_avg{1}.ex.LFPsemC(1,:,i), '.-m') %near
            errorbar(day_avg{1}.ex.uniqueors(i,:), day_avg{1}.ex.LFPmeanC(2,:,i), day_avg{1}.ex.LFPsemC(2,:,i), '.-k') %far
            xlabel('Radius (deg)');
            ylabel('Power (dB)');
            title('Corridor: LFP power')
        end
    end
    %Brick background
    if isfield(day_avg{1}.ex, 'LFPmeanB')
        if size(day_avg{1}.ex.LFPmeanB,3) >= i
            subplot(2,3,3); hold on;
            errorbar(day_avg{1}.ex.uniqueors(i,:), day_avg{1}.ex.LFPmeanB(1,:,i), day_avg{1}.ex.LFPsemB(1,:,i), '.-m') %near
            errorbar(day_avg{1}.ex.uniqueors(i,:), day_avg{1}.ex.LFPmeanB(2,:,i), day_avg{1}.ex.LFPsemB(2,:,i), '.-k') %far
            xlabel('Radius (deg)');
            ylabel('Power (dB)');
            title('Brick: LFP power')
        end
    end

    
    %LFP response
    %Gray background
    if isfield(day_avg{1}.ex, 'RESPmeanG')
        subplot(2,3,4); hold on;
        errorbar(day_avg{1}.ex.uniqueors(i,:), day_avg{1}.ex.RESPmeanG(1,:,i), day_avg{1}.ex.RESPsemG(1,:,i), '.-m') %near
        errorbar(day_avg{1}.ex.uniqueors(i,:), day_avg{1}.ex.RESPmeanG(2,:,i), day_avg{1}.ex.RESPsemG(2,:,i), '.-k') %far
        xlabel('Radius (deg)');
        ylabel('Response (dB)');
        title('Gray: LFP response')
    end
    %Corridor background
    if isfield(day_avg{1}.ex, 'RESPmeanC')
        if size(day_avg{1}.ex.RESPmeanC,3) >= i
            subplot(2,3,5); hold on;
            errorbar(day_avg{1}.ex.uniqueors(i,:), day_avg{1}.ex.RESPmeanC(1,:,i), day_avg{1}.ex.RESPsemC(1,:,i), '.-m') %near
            errorbar(day_avg{1}.ex.uniqueors(i,:), day_avg{1}.ex.RESPmeanC(2,:,i), day_avg{1}.ex.RESPsemC(2,:,i), '.-k') %far
            xlabel('Radius (deg)');
            ylabel('Response (dB)');
            title('Corridor: LFP response')
        end
    end
    %Brick background
    if isfield(day_avg{1}.ex, 'RESPmeanB')
        if size(day_avg{1}.ex.RESPmeanB,3) >= i
            subplot(2,3,6); hold on;
            errorbar(day_avg{1}.ex.uniqueors(i,:), day_avg{1}.ex.RESPmeanB(1,:,i), day_avg{1}.ex.RESPsemB(1,:,i), '.-m') %near
            errorbar(day_avg{1}.ex.uniqueors(i,:), day_avg{1}.ex.RESPmeanB(2,:,i), day_avg{1}.ex.RESPsemB(2,:,i), '.-k') %far
            xlabel('Radius (deg)');
            ylabel('Response (dB)');
            title('Brick: LFP response')
        end
    end
    pause
    close(gcf)
end



%%
% Mean of LFP power, across all days (each day = average across channels above)

if isfield(day_avg{1}.ex, 'LFPmeanG')
    day_avg{1}.ex.alldays.LFPmeanG  = mean(day_avg{1}.ex.LFPmeanG,3);
    day_avg{1}.ex.alldays.LFPsemG   = std(day_avg{1}.ex.LFPmeanG,0,3)  ./ sqrt(size(day_avg{1}.ex.LFPmeanG,3));
    day_avg{1}.ex.alldays.RESPmeanG = mean(day_avg{1}.ex.RESPmeanG,3);
    day_avg{1}.ex.alldays.RESPsemG   = std(day_avg{1}.ex.RESPmeanG,0,3)  ./ sqrt(size(day_avg{1}.ex.RESPmeanG,3));  
end

if isfield(day_avg{1}.ex, 'LFPmeanC')
    day_avg{1}.ex.alldays.LFPmeanC  = mean(day_avg{1}.ex.LFPmeanC,3);
    day_avg{1}.ex.alldays.LFPsemC   = std(day_avg{1}.ex.LFPmeanC,0,3)  ./ sqrt(size(day_avg{1}.ex.LFPmeanC,3));
    day_avg{1}.ex.alldays.RESPmeanC = mean(day_avg{1}.ex.RESPmeanC,3);
    day_avg{1}.ex.alldays.RESPsemC   = std(day_avg{1}.ex.RESPmeanC,0,3)  ./ sqrt(size(day_avg{1}.ex.RESPmeanC,3));  
end

if isfield(day_avg{1}.ex, 'LFPmeanB')
    day_avg{1}.ex.alldays.LFPmeanB  = mean(day_avg{1}.ex.LFPmeanB,3);
    day_avg{1}.ex.alldays.LFPsemB   = std(day_avg{1}.ex.LFPmeanB,0,3)  ./ sqrt(size(day_avg{1}.ex.LFPmeanB,3));
    day_avg{1}.ex.alldays.RESPmeanB = mean(day_avg{1}.ex.RESPmeanB,3);
    day_avg{1}.ex.alldays.RESPsemB   = std(day_avg{1}.ex.RESPmeanB,0,3)  ./ sqrt(size(day_avg{1}.ex.RESPmeanB,3));  
end

%%
% Population plot, mean of LFP power across all days (each day = average across channels)

ordradii = 1:1:length(days{1}.uniqueors); %ordinals for increasing ring size
figure;

%Gray background, LFP power
if isfield(day_avg{1}.ex.alldays, 'LFPmeanG')
    subplot(2,3,1); hold on;
    errorbar(ordradii, day_avg{1}.ex.alldays.LFPmeanG(1,:), day_avg{1}.ex.alldays.LFPsemG(1,:), '.-m') %near
    errorbar(ordradii, day_avg{1}.ex.alldays.LFPmeanG(2,:), day_avg{1}.ex.alldays.LFPsemG(2,:), '.-k') %far
    xlabel('Radius (ord)');
    ylabel('Power (dB)');
    title('Gray: Mean of LFP power, across days')
    hold off
end

%Corridor background, LFP power
if isfield(day_avg{1}.ex.alldays, 'LFPmeanC')
    subplot(2,3,2); hold on;
    errorbar(ordradii, day_avg{1}.ex.alldays.LFPmeanC(1,:), day_avg{1}.ex.alldays.LFPsemC(1,:), '.-m') %near
    errorbar(ordradii, day_avg{1}.ex.alldays.LFPmeanC(2,:), day_avg{1}.ex.alldays.LFPsemC(2,:), '.-k') %far
    xlabel('Radius (ord)');
    ylabel('Power (dB)');
    title('Corridor: Mean of LFP power, across days')
    hold off
end

%Brick background, LFP power
if isfield(day_avg{1}.ex.alldays, 'LFPmeanB')
    subplot(2,3,3); hold on;
    errorbar(ordradii, day_avg{1}.ex.alldays.LFPmeanB(1,:), day_avg{1}.ex.alldays.LFPsemB(1,:), '.-m') %near
    errorbar(ordradii, day_avg{1}.ex.alldays.LFPmeanB(2,:), day_avg{1}.ex.alldays.LFPsemB(2,:), '.-k') %far
    xlabel('Radius (ord)');
    ylabel('Power (dB)');
    title('Brick: Mean of LFP power, across days')
    hold off
end

%Gray background, LFP response
if isfield(day_avg{1}.ex.alldays, 'RESPmeanG')
    subplot(2,3,4); hold on;
    errorbar(ordradii, day_avg{1}.ex.alldays.RESPmeanG(1,:), day_avg{1}.ex.alldays.RESPsemG(1,:), '.-m') %near
    errorbar(ordradii, day_avg{1}.ex.alldays.RESPmeanG(2,:), day_avg{1}.ex.alldays.RESPsemG(2,:), '.-k') %far
    xlabel('Radius (ord)');
    ylabel('Power (dB)');
    title('Gray: Mean of LFP response, across days')
    hold off
end

%Corridor background, LFP response
if isfield(day_avg{1}.ex.alldays, 'RESPmeanC')
    subplot(2,3,5); hold on;
    errorbar(ordradii, day_avg{1}.ex.alldays.RESPmeanC(1,:), day_avg{1}.ex.alldays.RESPsemC(1,:), '.-m') %near
    errorbar(ordradii, day_avg{1}.ex.alldays.RESPmeanC(2,:), day_avg{1}.ex.alldays.RESPsemC(2,:), '.-k') %far
    xlabel('Radius (ord)');
    ylabel('Power (dB)');
    title('Corridor: Mean of LFP response, across days')
    hold off
end

%Brick background, LFP response
if isfield(day_avg{1}.ex.alldays, 'RESPmeanB')
    subplot(2,3,6); hold on;
    errorbar(ordradii, day_avg{1}.ex.alldays.RESPmeanB(1,:), day_avg{1}.ex.alldays.RESPsemB(1,:), '.-m') %near
    errorbar(ordradii, day_avg{1}.ex.alldays.RESPmeanB(2,:), day_avg{1}.ex.alldays.RESPsemB(2,:), '.-k') %far
    xlabel('Radius (ord)');
    ylabel('Power (dB)');
    title('Brick: Mean of LFP response, across days')
    hold off
end


%%
% Mean of LFP power, across days, per channel

 LFPchannelsG = [];
RESPchannelsG = [];
 LFPchannelsC = [];
RESPchannelsC = [];
 LFPchannelsB = [];
RESPchannelsB = [];
for i = 1:length(day_avg)
    if isfield(day_avg{i}, 'LFPmeanG')
        LFPchannelsG = cat(4, LFPchannelsG, days{i}.LFPmeanG);
        RESPchannelsG = cat(4, RESPchannelsG, days{i}.RESPmeanG);
    end
    if isfield(days{i}, 'LFPmeanC')
        LFPchannelsC = cat(4, LFPchannelsC, days{i}.LFPmeanC);
        RESPchannelsC = cat(4, RESPchannelsC, days{i}.RESPmeanC);
    end
    if isfield(days{i}, 'LFPmeanB')
        LFPchannelsB = cat(4, LFPchannelsB, days{i}.LFPmeanB);
        RESPchannelsB = cat(4, RESPchannelsB, days{i}.RESPmeanB);
    end
end

if isfield(days{1}, 'LFPmeanG')
    chan_avg{1}.LFPmeanG  = mean(LFPchannelsG,4);
    chan_avg{1}.RESPmeanG = mean(RESPchannelsG,4);
    chan_avg{1}.LFPsemG  = std(LFPchannelsG,0,4)  ./ sqrt(size(LFPchannelsG,4));
    chan_avg{1}.RESPsemG = std(RESPchannelsG,0,4) ./ sqrt(size(RESPchannelsG,4));
end

if isfield(days{1}, 'LFPmeanC')
    chan_avg{1}.LFPmeanC  = mean(LFPchannelsC,4);
    chan_avg{1}.RESPmeanC = mean(RESPchannelsC,4);
    chan_avg{1}.LFPsemC  = std(LFPchannelsC,0,4)  ./ sqrt(size(LFPchannelsC,4));
    chan_avg{1}.RESPsemC = std(RESPchannelsC,0,4) ./ sqrt(size(RESPchannelsC,4));
end

if isfield(days{1}, 'LFPmeanB')
    chan_avg{1}.LFPmeanB  = mean(LFPchannelsB,4);
    chan_avg{1}.RESPmeanB = mean(RESPchannelsB,4);
    chan_avg{1}.LFPsemB  = std(LFPchannelsB,0,4)  ./ sqrt(size(LFPchannelsB,4));
    chan_avg{1}.RESPsemB = std(RESPchannelsB,0,4) ./ sqrt(size(RESPchannelsB,4));
end


%%
% Plots of individual channels, mean of LFP power across days

ordradii = 1:1:length(days{1}.uniqueors); %ordinals for increasing ring size

%Gray background, LFP power
if isfield(chan_avg{1}, 'LFPmeanG')
    figure;
    numchans = size(chan_avg{1}.LFPmeanG,3);
    for j = 1:numchans
        subplot(ceil(sqrt(numchans)),ceil(sqrt(numchans)),j);
        hold on;
        errorbar(ordradii, chan_avg{1}.LFPmeanG(1,:,j), chan_avg{1}.LFPsemG(1,:,j), '.-m') %near
        errorbar(ordradii, chan_avg{1}.LFPmeanG(2,:,j), chan_avg{1}.LFPsemG(2,:,j), '.-k') %far
        xlabel('Radius (ord)');
        ylabel('Power (dB)');
    end
    set(gcf, 'Name', 'Gray: Mean of LFP power, across days')
    pause
    close(gcf)
end

%Gray background, LFP response
if isfield(chan_avg{1}, 'RESPmeanG')
    figure;
    numchans = size(chan_avg{1}.RESPmeanG,3);
    for j = 1:numchans
        subplot(ceil(sqrt(numchans)),ceil(sqrt(numchans)),j);
        hold on;
        errorbar(ordradii, chan_avg{1}.RESPmeanG(1,:,j), chan_avg{1}.RESPsemG(1,:,j), '.-m') %near
        errorbar(ordradii, chan_avg{1}.RESPmeanG(2,:,j), chan_avg{1}.RESPsemG(2,:,j), '.-k') %far
        xlabel('Radius (ord)');
        ylabel('Response (dB)');
    end
    set(gcf, 'Name', 'Gray: Mean of LFP change-from-baseline, across days')
    pause
    close(gcf)
end

%Corridor background, LFP power
if isfield(chan_avg{1}, 'LFPmeanC')
    figure;
    numchans = size(chan_avg{1}.LFPmeanC,3);
    for j = 1:numchans
        subplot(ceil(sqrt(numchans)),ceil(sqrt(numchans)),j);
        hold on;
        errorbar(ordradii, chan_avg{1}.LFPmeanC(1,:,j), chan_avg{1}.LFPsemC(1,:,j), '.-m') %near
        errorbar(ordradii, chan_avg{1}.LFPmeanC(2,:,j), chan_avg{1}.LFPsemC(2,:,j), '.-k') %far
        xlabel('Radius (ord)');
        ylabel('Power (dB)');
    end
    set(gcf, 'Name', 'Corridor: Mean of LFP power, across days')
    pause
    close(gcf)
end

%Corridor background, LFP response
if isfield(chan_avg{1}, 'RESPmeanC')
    figure;
    numchans = size(chan_avg{1}.RESPmeanC,3);
    for j = 1:numchans
        subplot(ceil(sqrt(numchans)),ceil(sqrt(numchans)),j);
        hold on;
        errorbar(ordradii, chan_avg{1}.RESPmeanC(1,:,j), chan_avg{1}.RESPsemC(1,:,j), '.-m') %near
        errorbar(ordradii, chan_avg{1}.RESPmeanC(2,:,j), chan_avg{1}.RESPsemC(2,:,j), '.-k') %far
        xlabel('Radius (ord)');
        ylabel('Response (dB)');
    end
    set(gcf, 'Name', 'Corridor: Mean of LFP change-from-baseline, across days')
    pause
    close(gcf)
end

%Brick background, LFP power
if isfield(chan_avg{1}, 'LFPmeanB')
    figure;
    numchans = size(chan_avg{1}.LFPmeanB,3);
    for j = 1:numchans
        subplot(ceil(sqrt(numchans)),ceil(sqrt(numchans)),j);
        hold on;
        errorbar(ordradii, chan_avg{1}.LFPmeanB(1,:,j), chan_avg{1}.LFPsemB(1,:,j), '.-m') %near
        errorbar(ordradii, chan_avg{1}.LFPmeanB(2,:,j), chan_avg{1}.LFPsemB(2,:,j), '.-k') %far
        xlabel('Radius (ord)');
        ylabel('Power (dB)');
    end
    set(gcf, 'Name', 'Brick: Mean of LFP power, across days')
    pause
    close(gcf)
end

%Brick background, LFP response
if isfield(chan_avg{1}, 'RESPmeanB')
    figure;
    numchans = size(chan_avg{1}.RESPmeanB,3);
    for j = 1:numchans
        subplot(ceil(sqrt(numchans)),ceil(sqrt(numchans)),j);
        hold on;
        errorbar(ordradii, chan_avg{1}.RESPmeanB(1,:,j), chan_avg{1}.RESPsemB(1,:,j), '.-m') %near
        errorbar(ordradii, chan_avg{1}.RESPmeanB(2,:,j), chan_avg{1}.RESPsemB(2,:,j), '.-k') %far
        xlabel('Radius (ord)');
        ylabel('Response (dB)');
    end
    set(gcf, 'Name', 'Brick: Mean of LFP change-from-baseline, across days')
    pause
    close(gcf)
end







%%
% Z-scored LFP power for individual channels

for i = 1:length(days)
    if isfield(days{i}, 'LFPmeanG')
        [LFPzG] = zscorefunc(days{i}.LFPmeanG, days{i}.uniqueors);
        days{i}.LFPzG = LFPzG;
    end
    if isfield(days{i}, 'RESPmeanG')
        [RESPzG] = zscorefunc(days{i}.RESPmeanG, days{i}.uniqueors);
        days{i}.RESPzG = RESPzG;
    end
    if isfield(days{i}, 'LFPmeanC')
        [LFPzC] = zscorefunc(days{i}.LFPmeanC, days{i}.uniqueors);
        days{i}.LFPzC = LFPzC;
    end
    if isfield(days{i}, 'RESPmeanC')
        [RESPzC] = zscorefunc(days{i}.RESPmeanC, days{i}.uniqueors);
        days{i}.RESPzC = RESPzC;
    end
    if isfield(days{i}, 'LFPmeanB')
        [LFPzB] = zscorefunc(days{i}.LFPmeanB, days{i}.uniqueors);
        days{i}.LFPzB = LFPzB;
    end
    if isfield(days{i}, 'RESPmeanB')
        [RESPzB] = zscorefunc(days{i}.RESPmeanB, days{i}.uniqueors);
        days{i}.RESPzB = RESPzB;
    end
end


%%
% plots of individual channels, Z-scored LFP power 

for i = 1:length(days)
    
    %Gray background, LFP power
    if isfield(days{i}, 'LFPzG')
        figure;
        numchans = size(days{i}.LFPzG,3);
        for j = 1:numchans
            subplot(ceil(sqrt(numchans)),ceil(sqrt(numchans)),j);
            hold on;
            plot(days{i}.uniqueors, days{i}.LFPzG(1,:,j), '.-m') %near
            plot(days{i}.uniqueors, days{i}.LFPzG(2,:,j), '.-k') %far
            xlabel('Radius (deg)');
            ylabel('Z-score');
        end
        set(gcf, 'Name', ['Day ' num2str(days{i}.dayname) ': Gray Z-scored LFP power'])
        pause
        close(gcf)
    end
    
    %Gray background, LFP response
    if isfield(days{i}, 'RESPzG')
        figure;
        numchans = size(days{i}.RESPzG,3);
        for j = 1:numchans
            subplot(ceil(sqrt(numchans)),ceil(sqrt(numchans)),j);
            hold on;
            plot(days{i}.uniqueors, days{i}.RESPzG(1,:,j), '.-m') %near
            plot(days{i}.uniqueors, days{i}.RESPzG(2,:,j), '.-k') %far
            xlabel('Radius (deg)');
            ylabel('Z-score');
        end
        set(gcf, 'Name', ['Day ' num2str(days{i}.dayname) ': Gray Z-scored LFP change-from-baseline'])
        pause
        close(gcf)
    end
    
    %Corridor background, LFP power
    if isfield(days{i}, 'LFPzC')
        figure;
        numchans = size(days{i}.LFPzC,3);
        for j = 1:numchans
            subplot(ceil(sqrt(numchans)),ceil(sqrt(numchans)),j);
            hold on;
            plot(days{i}.uniqueors, days{i}.LFPzC(1,:,j), '.-m') %near
            plot(days{i}.uniqueors, days{i}.LFPzC(2,:,j), '.-k') %far
            xlabel('Radius (deg)');
            ylabel('Z-score');
        end
        set(gcf, 'Name', ['Day ' num2str(days{i}.dayname) ': Corridor Z-scored LFP power'])
        pause
        close(gcf)
    end
    
    %Corridor background, LFP response
    if isfield(days{i}, 'RESPzC')
        figure;
        numchans = size(days{i}.RESPzC,3);
        for j = 1:numchans
            subplot(ceil(sqrt(numchans)),ceil(sqrt(numchans)),j);
            hold on;
            plot(days{i}.uniqueors, days{i}.RESPzC(1,:,j), '.-m') %near
            plot(days{i}.uniqueors, days{i}.RESPzC(2,:,j), '.-k') %far
            xlabel('Radius (deg)');
            ylabel('Z-score');
        end
        set(gcf, 'Name', ['Day ' num2str(days{i}.dayname) ': Corridor Z-scored LFP change-from-baseline'])
        pause
        close(gcf)
    end
    
    %Brick background, LFP power
    if isfield(days{i}, 'LFPzB')
        figure;
        numchans = size(days{i}.LFPzB,3);
        for j = 1:numchans
            subplot(ceil(sqrt(numchans)),ceil(sqrt(numchans)),j);
            hold on;
            plot(days{i}.uniqueors, days{i}.LFPzB(1,:,j), '.-m') %near
            plot(days{i}.uniqueors, days{i}.LFPzB(2,:,j), '.-k') %far
            xlabel('Radius (deg)');
            ylabel('Z-score');
        end
        set(gcf, 'Name', ['Day ' num2str(days{i}.dayname) ': Brick Z-scored LFP power'])
        pause
        close(gcf)
    end
    
    %Brick background, LFP response
    if isfield(days{i}, 'RESPzB')
        figure;
        numchans = size(days{i}.RESPzB,3);
        for j = 1:numchans
            subplot(ceil(sqrt(numchans)),ceil(sqrt(numchans)),j);
            hold on;
            plot(days{i}.uniqueors, days{i}.RESPzB(1,:,j), '.-m') %near
            plot(days{i}.uniqueors, days{i}.RESPzB(2,:,j), '.-k') %far
            xlabel('Radius (deg)');
            ylabel('Z-score');
        end
        set(gcf, 'Name', ['Day ' num2str(days{i}.dayname) ': Brick Z-scored LFP change-from-baseline'])
        pause
        close(gcf)
    end
end




%%
% Mean of z-scored LFP power, across days, per channel

 LFPchannelsG = [];
RESPchannelsG = [];
 LFPchannelsC = [];
RESPchannelsC = [];
 LFPchannelsB = [];
RESPchannelsB = [];
for i = 1:length(days)
    if isfield(days{i}, 'LFPzG')
        LFPchannelsG = cat(4, LFPchannelsG, days{i}.LFPzG);
        RESPchannelsG = cat(4, RESPchannelsG, days{i}.RESPzG);
    end
    if isfield(days{i}, 'LFPzC')
        LFPchannelsC = cat(4, LFPchannelsC, days{i}.LFPzC);
        RESPchannelsC = cat(4, RESPchannelsC, days{i}.RESPzC);
    end
    if isfield(days{i}, 'LFPzB')
        LFPchannelsB = cat(4, LFPchannelsB, days{i}.LFPzB);
        RESPchannelsB = cat(4, RESPchannelsB, days{i}.RESPzB);
    end
end

if isfield(days{1}, 'LFPzG')
    chan_avg{1}.LFPmeanzG  = mean(LFPchannelsG,4);
    chan_avg{1}.RESPmeanzG = mean(RESPchannelsG,4);
    chan_avg{1}.LFPsemzG  = std(LFPchannelsG,0,4)  ./ sqrt(size(LFPchannelsG,4));
    chan_avg{1}.RESPsemzG = std(RESPchannelsG,0,4) ./ sqrt(size(RESPchannelsG,4));
end

if isfield(days{1}, 'LFPzC')
    chan_avg{1}.LFPmeanzC  = mean(LFPchannelsC,4);
    chan_avg{1}.RESPmeanzC = mean(RESPchannelsC,4);
    chan_avg{1}.LFPsemzC  = std(LFPchannelsC,0,4)  ./ sqrt(size(LFPchannelsC,4));
    chan_avg{1}.RESPsemzC = std(RESPchannelsC,0,4) ./ sqrt(size(RESPchannelsC,4));
end

if isfield(days{1}, 'LFPzB')
    chan_avg{1}.LFPmeanzB  = mean(LFPchannelsB,4);
    chan_avg{1}.RESPmeanzB = mean(RESPchannelsB,4);
    chan_avg{1}.LFPsemzB  = std(LFPchannelsB,0,4)  ./ sqrt(size(LFPchannelsB,4));
    chan_avg{1}.RESPsemzB = std(RESPchannelsB,0,4) ./ sqrt(size(RESPchannelsB,4));
end


%%
% Plots of individual channels, mean of z-scored LFP power across days

ordradii = 1:1:length(days{1}.uniqueors); %ordinals for increasing ring size

%Gray background, LFP power
if isfield(chan_avg{1}, 'LFPmnzG')
    figure;
    numchans = size(chan_avg{1}.LFPmnzG,3);
    for j = 1:numchans
        subplot(ceil(sqrt(numchans)),ceil(sqrt(numchans)),j);
        hold on;
        errorbar(ordradii, chan_avg{1}.LFPmnzG(1,:,j), chan_avg{1}.LFPsemG(1,:,j), '.-m') %near
        errorbar(ordradii, chan_avg{1}.LFPmnzG(2,:,j), chan_avg{1}.LFPsemG(2,:,j), '.-k') %far
        xlabel('Radius (ord)');
        ylabel('Z-score');
    end
    set(gcf, 'Name', 'Gray: Mean of z-scored LFP power, across days')
    pause
    close(gcf)
end

%Gray background, LFP response
if isfield(chan_avg{1}, 'RESPmnzG')
    figure;
    numchans = size(chan_avg{1}.RESPmnzG,3);
    for j = 1:numchans
        subplot(ceil(sqrt(numchans)),ceil(sqrt(numchans)),j);
        hold on;
        errorbar(ordradii, chan_avg{1}.RESPmnzG(1,:,j), chan_avg{1}.RESPsemG(1,:,j), '.-m') %near
        errorbar(ordradii, chan_avg{1}.RESPmnzG(2,:,j), chan_avg{1}.RESPsemG(2,:,j), '.-k') %far
        xlabel('Radius (ord)');
        ylabel('Z-score');
    end
    set(gcf, 'Name', 'Gray: Mean of z-scored LFP change-from-baseline, across days')
    pause
    close(gcf)
end

%Corridor background, LFP power
if isfield(chan_avg{1}, 'LFPmnzC')
    figure;
    numchans = size(chan_avg{1}.LFPmnzC,3);
    for j = 1:numchans
        subplot(ceil(sqrt(numchans)),ceil(sqrt(numchans)),j);
        hold on;
        errorbar(ordradii, chan_avg{1}.LFPmnzC(1,:,j), chan_avg{1}.LFPsemC(1,:,j), '.-m') %near
        errorbar(ordradii, chan_avg{1}.LFPmnzC(2,:,j), chan_avg{1}.LFPsemC(2,:,j), '.-k') %far
        xlabel('Radius (ord)');
        ylabel('Z-score');
    end
    set(gcf, 'Name', 'Corridor: Mean of z-scored LFP power, across days')
    pause
    close(gcf)
end

%Corridor background, LFP response
if isfield(chan_avg{1}, 'RESPmnzC')
    figure;
    numchans = size(chan_avg{1}.RESPmnzC,3);
    for j = 1:numchans
        subplot(ceil(sqrt(numchans)),ceil(sqrt(numchans)),j);
        hold on;
        errorbar(ordradii, chan_avg{1}.RESPmnzC(1,:,j), chan_avg{1}.RESPsemC(1,:,j), '.-m') %near
        errorbar(ordradii, chan_avg{1}.RESPmnzC(2,:,j), chan_avg{1}.RESPsemC(2,:,j), '.-k') %far
        xlabel('Radius (ord)');
        ylabel('Z-score');
    end
    set(gcf, 'Name', 'Corridor: Mean of z-scored LFP change-from-baseline, across days')
    pause
    close(gcf)
end

%Brick background, LFP power
if isfield(chan_avg{1}, 'LFPmnzB')
    figure;
    numchans = size(chan_avg{1}.LFPmnzB,3);
    for j = 1:numchans
        subplot(ceil(sqrt(numchans)),ceil(sqrt(numchans)),j);
        hold on;
        errorbar(ordradii, chan_avg{1}.LFPmnzB(1,:,j), chan_avg{1}.LFPsemB(1,:,j), '.-m') %near
        errorbar(ordradii, chan_avg{1}.LFPmnzB(2,:,j), chan_avg{1}.LFPsemB(2,:,j), '.-k') %far
        xlabel('Radius (ord)');
        ylabel('Z-score');
    end
    set(gcf, 'Name', 'Brick: Mean of z-scored LFP power, across days')
    pause
    close(gcf)
end

%Brick background, LFP response
if isfield(chan_avg{1}, 'RESPmnzB')
    figure;
    numchans = size(chan_avg{1}.RESPmnzB,3);
    for j = 1:numchans
        subplot(ceil(sqrt(numchans)),ceil(sqrt(numchans)),j);
        hold on;
        errorbar(ordradii, chan_avg{1}.RESPmnzB(1,:,j), chan_avg{1}.RESPsemB(1,:,j), '.-m') %near
        errorbar(ordradii, chan_avg{1}.RESPmnzB(2,:,j), chan_avg{1}.RESPsemB(2,:,j), '.-k') %far
        xlabel('Radius (ord)');
        ylabel('Z-score');
    end
    set(gcf, 'Name', 'Brick: Mean of z-scored LFP change-from-baseline, across days')
    pause
    close(gcf)
end

%%
% Mean of z-scored LFP power, across channels, per day

%Gray background
for i = 1:length(days)
    if isfield(days{i}, 'zscoresG')
        mnzG = mean(days{i}.zscoresG(:,:,channels),3);
        semG = std(days{i}.zscoresG(:,:,channels),0,3) ./ sqrt(size(zscoresG(:,:,channels),3));
        days{i}.mnzG = mnzG;
        days{i}.semG = semG;
    end

end

