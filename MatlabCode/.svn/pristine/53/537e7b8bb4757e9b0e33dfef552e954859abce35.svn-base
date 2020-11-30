function [bkgnd, near, far, uniqueors, spikename] = SMurray_unit_spikes_all(stro)

% spike data analysis of SMurray paradigm
% average spike rate per smaller time windows
aw = .3; %large analysis window time, starting at stimulus onset
sw = 6; %number of smaller analysis windows
ws = aw/sw; % smaller analysis window length (s)

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

near = NaN(ceil(ntrials/size(uniquefpx,1)/length(uniqueors)), sw, length(uniqueors), length(choices)); 
far = NaN(ceil(ntrials/size(uniquefpx,1)/length(uniqueors)), sw, length(uniqueors), length(choices));  
spikerates = NaN(ntrials, sw, length(choices));  

for k = 1:length(choices)
    spikeidx = strcmp(stro.sum.rasterCells(1,:),spikename(k,:));

    for i = 1:ntrials
        spiketimes = stro.ras{i,spikeidx};
        nspikes = spiketimes > stimon_t(i) & spiketimes <= stimon_t(i)+aw; 
        allspikes = spiketimes(nspikes);
        allspikes = allspikes - stimon_t(i); %convert spike time to 'time post-stim onset'
        for m = 1:sw;
            windowend = m * ws;
            winsp = sum(allspikes > windowend - ws & allspikes <= windowend);
            spikerates(i,m,k) = winsp/ws; 
        end 
    end
    
    for j = 1:length(uniqueors)
        Lor = stim_or == uniqueors(j);
        Lnear = fpx == uniquefpx(1,1);
        Lfar = fpx == uniquefpx(2,1);
        Ln = Lor & Lnear;
        Lf = Lor & Lfar;
        near(1:size(spikerates(Ln,:,k),1),:,j,k) = spikerates(Ln,:,k);
        far(1:size(spikerates(Lf,:,k),1),:,j,k) = spikerates(Lf,:,k);
    end
end

%removes trials if the number of trials/stim condition are not equal
if length(find(isnan(near))) + length(find(isnan(far))) > 0
    near = near(1:size(near,1)-1, :, :, :);
    far = far(1:size(far,1)-1, :, :, :);
end

end