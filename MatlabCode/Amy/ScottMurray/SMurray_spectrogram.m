% LFP analysis using Chronux software
% from SMurray paradigm with two fixation points



%stro = nex2stro(findfile('A113011008'));

fp_acq = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_acq'));
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
fpx = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_x'));
fpy = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_y'));
stim_ir = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_ir'));
stim_or = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_or'));
bkgnd = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd'));
ntrials = size(stro.trial,1);
LFPsamprate = stro.sum.analog.storeRates{1};  % Assuming all channels are sampled at the same rate
LFPstarttimes = [stro.ras{:,strcmp(stro.sum.rasterCells,'anlgStartTime')}];

uniquefpxy = unique([fpx fpy],'rows');
uniquefpx = unique(fpx);
uniquefpy = unique(fpy);
uniqueors = unique(stim_or);
uniqueirs = unique(stim_ir);

%determine if using a electrode/array
electrode = ismember(stro.sum.rasterCells, 'AD01');
if  any(electrode) %electrode
    numelec = 1;
else %array
    numelec = 1:32;
end



%%
%pull LFP data for electrode/array
%first determine size of 'data' to be filled
PreStimTime = 0.05; %time analyzed prior to stimulus onset (seconds)
PostStimTime = 0.05; %time analyzed post stimulus offset (seconds)
StimTime = mean(stimoff_t - stimon_t); %stimulus time (seconds)
AW = PreStimTime + StimTime + PostStimTime; %total analysis window time period
data = nan*ones(ntrials, round(AW*LFPsamprate), numelec(end));

for whichADchan = numelec
    for i = 1:ntrials
        if size(numelec, 2) > 1 %array
            LFP = stro.ras{i,strcmp(stro.sum.rasterCells,['AD',num2str(16+whichADchan)])};
        else %electrode
            LFP = stro.ras{i,strcmp(stro.sum.rasterCells,'AD01')};
        end
        LFPtimes = LFPstarttimes(i)+(0:length(LFP)-1)/LFPsamprate; % sec
        LFPstim = LFP(LFPtimes >= stimon_t(i)-PreStimTime & LFPtimes <= stimoff_t(i)+PostStimTime);
        len = length(LFPstim);
        data(i,1:len,whichADchan) = LFPstim;
    end
end


%mean of LFP across channels
dataAllChannels = mean(data, 3);
nearfar = nan(size(uniquefpx,1), length(uniqueors), size(dataAllChannels,2));
for i = 1:size(uniquefpx,1)
    for j = 1:length(uniqueors)
        Lor = stim_or == uniqueors(j);
        Lfp = fpx == uniquefpx(i,1);
        L = Lor & Lfp;
        nearfar(i,j,:) = mean(dataAllChannels(L,:)); % near, then far
    end
end


%%

% LFP spectrogram (using Chronux software)


% Preliminary parameter settings
movingwin=[0.05 0.005]; % set the moving window dimensions
params.Fs=LFPsamprate; % sampling frequency
params.fpass=[0 200]; % frequency of interest
params.tapers=[5 9]; % tapers
params.trialave=0; % 1 = average over trials
params.err=0; % 1 = error computation

% Plot spectrograms
%figure
for i = 1:size(nearfar,1)
    for j = 1:size(nearfar,2)
        figure%subplot(size(nearfar,1),size(nearfar,2), j + (i-1)*7)
        hold on
        dataplot = nan((size(nearfar,3)-1), 1);
        for k = 1:size(dataplot,1)
            dataplot(k,1) = nearfar(i,j,k);
        end
        [S1,t,f]=mtspecgramc(dataplot,movingwin,params); % compute spectrogram
        plot_matrix(10*log10(S1),t,f); % plot spectrogram
        caxis([8 28]); colorbar;
        set(gca,'FontName','Times New Roman','Fontsize', 14);
        title({['LFP,  W=' num2str(params.tapers(1)/movingwin(1)) 'Hz']; ['moving window = ' num2str(movingwin(1)) 's, step = ' num2str(movingwin(2)) 's']});
        xlabel('time s');
        ylabel('frequency Hz');
        position_figure
    end
end