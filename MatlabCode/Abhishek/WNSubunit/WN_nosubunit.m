% processing

tic;
clear all;
close all;
plot_counter = 1;
[stro,filename,no_subunit] = library(); % get one file from the library

global spikename maskidx spikeidx nstixperside ngammasteps seedidx nframesidx stimonidx muidxs sigmaidxs
global msperframe ntrials maxT xx yy M
spikename = getSpikenum(stro);
maskidx = strcmp(stro.sum.rasterCells(1,:), 'subunit_mask');
spikeidx = strcmp(stro.sum.rasterCells(1,:),spikename);
nstixperside = stro.sum.exptParams.nstixperside;
ngammasteps = 2^16; % 65536
seedidx = strcmp(stro.sum.trialFields(1,:),'seed');
nframesidx = strcmp(stro.sum.trialFields(1,:),'num_frames');
stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on');
muidxs = [find(strcmp(stro.sum.trialFields(1,:),'mu1')), ...
    find(strcmp(stro.sum.trialFields(1,:),'mu2')), ...
    find(strcmp(stro.sum.trialFields(1,:),'mu3'))];
sigmaidxs = [find(strcmp(stro.sum.trialFields(1,:),'sigma1')), ...
    find(strcmp(stro.sum.trialFields(1,:),'sigma2')), ...
    find(strcmp(stro.sum.trialFields(1,:),'sigma3'))];
msperframe = 1000/stro.sum.exptParams.framerate;
ntrials = size(stro.trial,1);
maxT = 15; % this represents the temporal part in the spatiotemporal receptive field
xx = linspace(stro.sum.exptParams.gauss_locut/1000, stro.sum.exptParams.gauss_hicut/1000,ngammasteps); % xx represents the probabilities. For more info, have a look at the MATLAB 'norminv' function.
yy = norminv(xx'); % defining norminv to extract the values for which the cdf values range between gauss_locut and gauss_hicut

%%
% Obtaining the M matrix, code extracted from Greg, fitting a cubic spline
% using the command 'spline'. 'SplineRaw' only availabe through
% psychtoolbox which I currently don't have now.
fundamentals = stro.sum.exptParams.fundamentals; % CONE FUNDAMENTALS: L,M,S
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S 
mon_spd = stro.sum.exptParams.mon_spd; % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals

%%
%%
% Determine when the mask changed in order to determine from what trials we started our experiment to study 
% the computation of the subunits which we are interested in.
if no_subunit ~= 1
    mask_changes = 1;
    all_masks = stro.ras(:,maskidx);
    for k = 2:ntrials
        if isequal(all_masks{k}, all_masks{k-1}) %|| all(all_masks{k} == 0) && any(isnan(all_masks{k-1}))
            continue
        else
            mask_changes = [mask_changes k-1 k]; %#ok<AGROW>
        end
    end
    % Each column of `mask_changes` holds the start and end trial index where each
    % mask was used. For example, the subunits A and B relevant to my project were added on trial number 69.
    mask_changes = reshape([mask_changes ntrials], 2, []);
    % disp(stro.ras{mask_changes(1,2),4}); % displays the time when the mask was selected
    temp = [];
    
    for i = 1:size(mask_changes,2)
        if strcmp('N072215002.nex',filename)
            temp = [temp, mask_changes(:,i);];
        elseif (mask_changes(2,i) - mask_changes(1,i) > 4)
            temp = [temp, mask_changes(:,i);];
        end
    end
    mask_changes = temp;
    
    % Ad-hoc chnages to some of the nex files
    if strcmp('N050415001.nex',filename)
        % file-specific changes of the masking trials
        mask_changes = user_defined_mask_changes(mask_changes,[2 3]);
    elseif strcmp('N061515001.nex',filename)
        mask_changes = user_defined_mask_changes(mask_changes,[3 4]);
    else
    end
    mask_changes = mask_changes(:,1);
    
else
    mask_changes = [1; ntrials];
end
%%

WNpixelplot(stro,plot_counter,mask_changes,'WN_nosubunit'); % calculates and plot the STAs and the PCs, performs significance testing also on them
toc;