% Analysis of eye position on some of the Whitenoise files
% Author - Abhishek De, 11/19
close all; clearvars

load fundamentals.mat 
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;
resize_fact = 1;
% Include Gun noise data for statistical tests in order to estimate the RF size
channels = 3;
NPOINTS = 65536;
resize_fact2 = 1;
CHI2CRIT = 0.95; % For flaging a stixel as significant (adding gun and cone noise z-scores, squared)
maxT = 9;
crit = chi2inv(CHI2CRIT,300); % 3 color channels
spikename_options = ['sig001a'; 'sig001b'];
[filename_Lum, spikeIdx_Lum] = fnamesFromTxt2('Lum.txt');
[filename_ColorOpponent, spikeIdx_ColorOpponent] = fnamesFromTxt2('ColorOpponent.txt');

% Loading the Neurothresh files
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
tmp_filename = fetch(conn,'SELECT filename FROM WNthresh');
NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
close(conn);
tmp_filename = tmp_filename(strcmp('subunit',NTmode));
spikeidx_NT = spikeidx_NT(strcmp('subunit',NTmode));
filename_NT = cell(size(tmp_filename));
for kk = 1:size(tmp_filename)
    filename_NT(kk) = {tmp_filename(kk)};
end

% Loading the just the subunit files 
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
tmp_filename = fetch(conn,'SELECT filename FROM WNSubunit');
subunit_mode = fetch(conn,'SELECT mode FROM WNSubunit');
spikeidx_subunit = cell2mat(fetch(conn,'SELECT spikeidx FROM WNSubunit'));
close(conn);
tmp_filename = tmp_filename(strcmp('STA',subunit_mode));
spikeidx_subunit = spikeidx_subunit(strcmp('STA',subunit_mode));
filename_subunit = cell(size(tmp_filename));
for kk = 1:size(tmp_filename)
    filename_subunit(kk) = {tmp_filename(kk)};
end

Input_List = [filename_Lum; filename_ColorOpponent; filename_NT; filename_subunit];
spikeIdx = [spikeIdx_Lum; spikeIdx_ColorOpponent; spikeidx_NT; spikeidx_subunit];

numcells = numel(Input_List);

% Loading a file
ii = 3;
flag = 0;
stro = {};
for jj = 1:size(Input_List{ii},2)
    try
        tmpstro = nex2stro(findfile(char(Input_List{ii}(jj))));
    catch
        flag = 1;
        break;
    end
    if (isempty(stro))
        stro = tmpstro;
    else
        stro = strocat(stro, tmpstro);
    end
    if ~any(strcmp(stro.sum.rasterCells,'sig001a_wf'))
        flag = 1;
    end
end
noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
basisvecidx = strcmp(stro.sum.rasterCells(1,:),'basis_vec');
maskidx = strcmp(stro.sum.rasterCells(1,:),'subunit_mask');
T = stro.trial(:,noisetypeidx)==1;
stro.ras(~T ,:) = []; % modiftying the WN structure
stro.trial(~T,:) = []; % modiftying the WN structure
mask_changes = [2 size(stro.trial,1)];
if any(basisvecidx)
    mask_changes = [2];
    all_masks = stro.ras(:,maskidx);
    Fx = @(xi) any(isnan(xi)); % function that finds 'NaN' in a cell array
    inds = find(cellfun(Fx,stro.ras(:,basisvecidx))==0);
    if isempty(inds)
        inds = size(stro.trial,1)-1;
    end
    last_wntrial =  inds(1)-1;
    for k = 3:last_wntrial
        if isequal(all_masks{k}, all_masks{k-1}) %|| all(all_masks{k} == 0) && any(isnan(all_masks{k-1}))
            continue
        else
            mask_changes = [mask_changes k-1 k]; %#ok<AGROW>
        end
    end
    if mask_changes(end) == last_wntrial
        mask_changes(end) = [];
    else
        mask_changes = [mask_changes  last_wntrial];
    end
    mask_changes = reshape(mask_changes , 2, []);
    mask_changes  = mask_changes(:,1);
    
    idxs = zeros(size(stro.trial,1),1);
    idxs(mask_changes(2,1)+1:end) = 1;
    idxs = logical(idxs);
    stro.ras(idxs,:) = []; % modiftying the WN structure
    stro.trial(idxs,:) = []; % modiftying the WN structure
end

samplingrate = stro.sum.analog.storeRates{1};
stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on');
fponidx = strcmp(stro.sum.trialFields(1,:),'fp_on');
fpacqidx = strcmp(stro.sum.trialFields(1,:),'fpacq');
spikeidx = strcmp(stro.sum.rasterCells(1,:),'sig001a');
anlgStartTimeidx = strcmp(stro.sum.rasterCells(1,:),'anlgStartTime');
eyexidx = strcmp(stro.sum.rasterCells(1,:),'AD11');
eyeyidx = strcmp(stro.sum.rasterCells(1,:),'AD12');

% Comparison with GetSacData
SacData = getSacData(stro,1,Inf,10); % Extracting get sac data from stro structure
timedur = [];
baselineFR = [];
count = 1;
for kk = mask_changes(1):mask_changes(2)
    Amp = SacData.amplitudes{kk};
    K = SacData.endtimes{kk};
    startime = K(find(K > stro.trial(kk,fponidx) & K < stro.trial(jj,stimonidx) & Amp > 1,1,'last'));
    if isempty(startime)
        startime = stro.trial(kk,fponidx);
    end
    stimontime = stro.trial(kk,stimonidx);
    timedur = [timedur; stimontime - startime];
    numspikes = sum(stro.ras{kk,spikeidx}>startime & stro.ras{kk,spikeidx}<stimontime);
    baselineFR = [baselineFR; numspikes/timedur(end)];
    
    % Just analyzing the eye traces
    analogstartime = stro.ras{kk,anlgStartTimeidx};
    eyex = stro.ras{kk,eyexidx}*4096/400; % Converting the eye traces into degrees of visual angle
    eyey = stro.ras{kk,eyeyidx}*4096/400;
    stimont = stro.trial(kk,stimonidx);
    timepts = (analogstartime:(1/samplingrate):analogstartime + (numel(eyey)-1)/samplingrate) - stimont;
    figure(1); subplot(231); plot(timepts,eyex,'r'); hold on;
    figure(1); subplot(232); plot(timepts,eyey,'g'); hold on;
    figure(1); subplot(233); plot(timepts,sqrt(eyex.^2+eyey.^2),'k'); hold on;
    figure(1); subplot(234); plot(stro.ras{kk,spikeidx} - stimontime,count*ones(size(stro.ras{kk,spikeidx})),'k.'); hold on;
    plot(-timedur(end),count,'rv'); 
    
    count = count + 1;
end
figure(1); subplot(231); axis square; set(gca,'Tickdir','out','Xlim',[-1 1]); xlabel('time'), ylabel('X pos'); hold off;
subplot(232); axis square; set(gca,'Tickdir','out','Xlim',[-1 1]); xlabel('time'), ylabel('X pos'); hold off;
subplot(233); axis square; set(gca,'Tickdir','out','Xlim',[-1 1]); xlabel('time'), ylabel('Amp'); hold off;
subplot(234); axis square; set(gca,'Tickdir','out','Xlim',[-1 1]); xlabel('time'), ylabel('Rasters'); hold off;
subplot(235); histogram(baselineFR,20,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); axis square; set(gca,'Tickdir','out'); ylabel('counts'), xlabel('Baseline FR'); hold off;
subplot(236); histogram(timedur,40,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); axis square; set(gca,'Tickdir','out','Xlim',[0 1]); ylabel('counts'), xlabel('Time dur'); hold off;



            