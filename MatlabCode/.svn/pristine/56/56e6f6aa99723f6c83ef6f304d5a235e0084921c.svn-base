function [LFPsamprate, near, far, uniqueors, nearBase, farBase] = SMurray_freq_LFP(stro, aw)

% pulls LFP data recorded during analysis window for SMurray paradigm
% 'aw' specifies analysis window in seconds: [time after stim onset; time after stim offset]

stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
fpx = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_x'));
stim_or = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_or'));
ntrials = size(stro.trial,1);
LFPsamprate = stro.sum.analog.storeRates{1};  % Assuming all channels are sampled at the same rate
LFPstarttimes = [stro.ras{:,strcmp(stro.sum.rasterCells,'anlgStartTime')}];
uniquefpx = unique(fpx);
uniqueors = unique(stim_or);

%determine if using a electrode/array
electrode = ismember(stro.sum.rasterCells, 'AD01');
if  any(electrode) %electrode
    numelec = 1;
else %array
    numelec = 1:32;
end

%pull LFP data for electrode/array
%first determine size of 'data' to be filled
awlength = ((mean(stimoff_t)+aw(2)) - (mean(stimon_t)+aw(1))); %total analysis window time period
data = NaN(ntrials, round(awlength*LFPsamprate)-1, numelec(end));
database = NaN(ntrials, round(awlength*LFPsamprate)-1, numelec(end));

for whichADchan = numelec
    for i = 1:ntrials
        if size(numelec, 2) > 1 %array
            LFPraw = stro.ras{i,strcmp(stro.sum.rasterCells,['AD',num2str(16+whichADchan)])};
        else %electrode
            LFPraw = stro.ras{i,strcmp(stro.sum.rasterCells,'AD01')};
        end
        LFPtimes = LFPstarttimes(i)+(0:length(LFPraw)-1)/LFPsamprate; % sec
        LFPstim = LFPraw(LFPtimes > stimon_t(i)+aw(1) & LFPtimes <= stimoff_t(i)+aw(2));
        LFPbase = LFPraw(LFPtimes > stimon_t(i)-awlength & LFPtimes <= stimon_t(i));
        if size(LFPstim,1) < size(data,2)
            LFPstim = LFPraw(LFPtimes > stimon_t(i)+aw(1) & LFPtimes <= stimoff_t(i)+aw(2)+.01);
            LFPbase = LFPraw(LFPtimes > stimon_t(i)-awlength & LFPtimes <= stimon_t(i)+.01);
        end
        data(i,1:round(awlength*LFPsamprate)-1,whichADchan) = LFPstim(1:round(awlength*LFPsamprate)-1);
        database(i,1:round(awlength*LFPsamprate)-1,whichADchan) = LFPbase(1:round(awlength*LFPsamprate)-1);
    end
end

%organize LFP data by stimulus type
trialsperstim = round(ntrials/size(uniquefpx,1)/length(uniqueors));
near = NaN(trialsperstim, size(data,2), size(data,3), length(uniqueors)); 
far =  NaN(trialsperstim, size(data,2), size(data,3), length(uniqueors));
nearBase = NaN(trialsperstim, size(data,2), size(data,3), length(uniqueors)); 
farBase =  NaN(trialsperstim, size(data,2), size(data,3), length(uniqueors));

for j = 1:length(uniqueors)
    Lor = stim_or == uniqueors(j);
    Lnear = fpx == uniquefpx(1,1);
    Lfar = fpx == uniquefpx(2,1);
    Ln = Lor & Lnear;
    Lf = Lor & Lfar;
    near(1:size(data(Ln),1),:,:,j) = data(Ln,:,:);
    far(1:size(data(Lf),1),:,:,j) = data(Lf,:,:);
    nearBase(1:size(database(Ln),1),:,:,j) = database(Ln,:,:);
    farBase(1:size(database(Lf),1),:,:,j) = database(Lf,:,:);
end


%removes trials if the number of trials/stim condition are not equal
if length(find(isnan(near))) + length(find(isnan(far))) > 0
    near = near(1:size(near,1)-1, :, :, :);
    far = far(1:size(far,1)-1, :, :, :);
    nearBase = nearBase(1:size(nearBase,1)-1, :, :, :);
    farBase = farBase(1:size(farBase,1)-1, :, :, :);
end

%radius size in degrees
uniqueors = uniqueors ./10;


end