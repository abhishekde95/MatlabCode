function [bkgnd, near, far, uniqueors, spikename, interstim, nearPre farPre] = SMurray_unit_spikes(stro, aw)

% spike data analysis of SMurray paradigm
% 'aw' specifies analysis window in seconds: [time after stim onset; time after stim offset]

stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
fix_acq = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_acq'));
fpx = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_x'));
stim_or = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_or'));
bkgnd = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd'));
ntrials = size(stro.trial,1);
uniquefpx = unique(fpx);
uniqueors = unique(stim_or);

%open each isolated unit
choices = {};
for i = 1:length(stro.sum.rasterCells)
    if (strncmp(stro.sum.rasterCells(i),'sig0',4))
        choices{length(choices)+1} = stro.sum.rasterCells{i};
    end
end
%weed out the lfp channels (if present)
lfpChannels = strfind(choices, '_wf');
lfpChannels = cellfun(@any, lfpChannels);
choices(lfpChannels) = [];
spikename = char(choices);

near = NaN(ceil(ntrials/size(uniquefpx,1)/length(uniqueors)), length(uniqueors), length(choices));
far = NaN(ceil(ntrials/size(uniquefpx,1)/length(uniqueors)), length(uniqueors), length(choices));  
spikerates = NaN(ntrials,length(choices));
interstim = NaN(ntrials,length(choices));
nearPre = NaN(ceil(ntrials/size(uniquefpx,1)/length(uniqueors)), length(uniqueors), length(choices));
farPre = NaN(ceil(ntrials/size(uniquefpx,1)/length(uniqueors)), length(uniqueors), length(choices));   

for k = 1:length(choices)
    spikeidx = strcmp(stro.sum.rasterCells(1,:),spikename(k,:));
    
    for i = 1:ntrials
        spiketimes = stro.ras{i,spikeidx};
        nspikes = sum(spiketimes > stimon_t(i)+aw(1) & spiketimes <= stimoff_t(i)+aw(2)); %analysis window
        awlength = ((stimoff_t(i)+aw(2)) - (stimon_t(i)+aw(1))); %total analysis window time period
        spikerates(i,k) = nspikes./awlength;
        interspikes = sum(spiketimes > fix_acq(i)+0.05 & spiketimes <= fix_acq(i)+0.05+awlength); %interstimulus analysis window
        interstim(i,k) = interspikes./awlength;
    end
    
    for j = 1:length(uniqueors)
        Lor = stim_or == uniqueors(j);
        Lnear = fpx == uniquefpx(1,1);
        Lfar = fpx == uniquefpx(2,1);
        Ln = Lor & Lnear;
        Lf = Lor & Lfar;
        near(1:size(spikerates(Ln,k),1),j,k) = spikerates(Ln,k);
        far(1:size(spikerates(Lf,k),1),j,k) = spikerates(Lf,k);
        nearPre(1:size(interstim(Ln,k),1),j,k) = interstim(Ln,k);
        farPre(1:size(interstim(Lf,k),1),j,k) = interstim(Lf,k);
    end
end

%removes trials if the number of trials/stim condition are not equal
if length(find(isnan(near))) + length(find(isnan(far))) > 0
    near = near(1:size(near,1)-1, :, :);
    far = far(1:size(far,1)-1, :, :);
    nearPre = nearPre(1:size(nearPre,1)-1, :, :);
    farPre = farPre(1:size(farPre,1)-1, :, :);
end

end