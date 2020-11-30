
Q = 10; % Cost for moving a spike by 1 s (adding or removing a spike costs '1')
signaloffset = 0.05;

DEBUG = 1;
if DEBUG % Randomizing the Lcc, Mcc, and TF values. This is destructive!
    idxs = randperm(length(Lcc));
    Lcc = Lcc(idxs);
    Mcc = Mcc(idxs);
    TF = TF(idxs);
end
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));

spikeidx = ismember(stro.sum.rasterCells(1,:),{'sig001a'});
% Setting up the spikes so that time 0 is stim on and cutting off
% spikes after stimoff
ntrials = length(stro.ras(:,spikeidx));
spikes = {};
startidxs = nan*ones(ntrials,1);
endidxs = nan*ones(ntrials,1);
for i = 1:ntrials
    tmpspikes = stro.ras{i,spikeidx}-stimon_t(i)-signaloffset;
    tmpspikes(tmpspikes<0) = [];
    tmpspikes(tmpspikes>stimoff_t(i)-stimon_t(i)) = [];
    if i == 1
        if isempty(tmpspikes)
            startidxs(i) = 0;
            endidxs(i) = -1;
        end
        startidxs(i) = 1;
        endidxs(i) = length(tmpspikes);
    else
        if isempty(tmpspikes)
            startidxs(i) = startidxs(i-1);
            endidxs(i) = startidxs(i)-1;
        end
        startidxs(i) = numel(vertcat(spikes{:}))+1;
        endidxs(i) = startidxs(i)+length(tmpspikes)-1;
    end
    spikes{i} = tmpspikes;
end
spiketimes = vertcat(spikes{:});

%figure; axes; hold on;
%for i = 1:length(startidxs)
%    plot(spiketimes(startidxs(i):endidxs(i)))
%end
% spkdl is blazingly fast so calculating *all* pairwise distances
d = spkdl(spiketimes,startidxs,endidxs,Q); % 10 ms
d_mat = reshape(d,length(startidxs),length(startidxs));
% imagesc(d_mat);

% ---------------------------------------------------------
% Doing the k-nearest neighbors calculation
% Getting the blanks in the correct format for spkdl
Lblank = Lcc == 0 & Mcc == 0;
nblanktrials = sum(Lblank);
data = [];
for i = 1:size(uniquestim,1)
    Lstim = Lcc == uniquestim(i,1) & Mcc == uniquestim(i,2) & TF == uniquestim(i,3);
    nsignaltrials = sum(Lstim);
    
    ntrials_tot = nblanktrials+nsignaltrials;
    sub_d_mat = [d_mat(Lblank,Lblank), d_mat(Lblank,Lstim);...
        d_mat(Lstim,Lblank), d_mat(Lstim,Lstim)];
    sub_d_mat = sub_d_mat+diag(nan(ntrials_tot,1));
    sorted_d_mat = sort(sub_d_mat); % sorts each column in ascending order (so first rows correspond to nearby spike trains)
    k = (nsignaltrials/2)+1-mod(nsignaltrials/2,2); % an odd number of nearest neighbors to avoid ties
    critical_values = sorted_d_mat(k,:); % K nearest neighbor
    k_close_points = sub_d_mat<=repmat(critical_values,size(sub_d_mat,1),1);
    if ~all(sum(k_close_points) == k)
        error('bug in nearest neighbor selection');
    end
    
    mask = [-1*ones(nblanktrials,ntrials_tot); ones(nsignaltrials,ntrials_tot)]; % -1 means blank, 1 means stim
    whichcategoryassigned = sum(k_close_points.*mask)>0; % false (0) means classified as a blank, true (1) means classified as a stim
    % For blank-against-blank classification *every* trial is
    % misclassified because every trial has a duplicate of itself with
    % distance 0
    if any(sum(k_close_points.*mask) == 0)
        disp('got a tie');
        if i ~= 1 & ~DEBUG
            keyboard
        end
        tieidxs = sum(k_close_points.*mask) == 0;
        whichcategoryassigned(tieidxs) = unidrnd(2,1,sum(tieidxs))-1; % flipping a coin
    end
    
    nmisclassified_blank = sum(whichcategoryassigned(1:nblanktrials) ~= 0);
    nmisclassified_signal = sum(whichcategoryassigned(nblanktrials+1:ntrials_tot) == 0);
    % Each term in the nmisclassified sum, above, is biased upward because there are more "unlike" trials than "like" trials
    ncorrectlyclassified_blank = sum(whichcategoryassigned(1:nblanktrials) == 0);
    ncorrectlyclassified_signal = sum(whichcategoryassigned(nblanktrials+1:ntrials_tot) ~= 0);
   
    % similarly, ncorrectlyclassified is biased downward
    
    % Prob. that a blank is correctly classified is:
    % P(# nearby blanks > # nearby stims) =
    % P(# nearby blanks > k/2) =
    % 1-hygecdf(k/2,ntrials_tot-1, nblanktrials-1, k)
    % probability that a blank is *incorrectly* classified is
    % hygecdf(k/2,ntrials_tot-1, nblanktrials-1, k)
    % So E(# of correctly classified blanks) = nblanktrials*P(a blank is correctly classified)
    E_nmisclassified_blank = nblanktrials*hygecdf(k/2,ntrials_tot-1, nblanktrials-1, k);
    E_nmisclassified_signal = nsignaltrials*hygecdf(k/2,ntrials_tot-1, nsignaltrials-1, k);
    E_correctlyclassified_blank = nblanktrials-E_nmisclassified_blank;
    E_correctlyclassified_signal = nsignaltrials-E_nmisclassified_signal;
    
    % Gives same answer as: nblanktrials*(1-hygecdf(k/2,ntrials_tot-1, nblanktrials-1, k)) + nsignaltrials*(1-hygecdf(k/2,ntrials_tot-1, nsignaltrials-1, k))
    
    data = [data; nmisclassified_blank nmisclassified_signal nblanktrials nsignaltrials E_nmisclassified_blank E_nmisclassified_signal];
    % figure; subplot(2,2,1); hist(data(2:end,2))
    % data is in # correct classified beyond what is expected by chance.
    
    %data = [data; nmisclassified ncorrectlyclassified];
    %data = [data; nmisclassified*(ntrials_tot/2)/E_nmis ncorrectlyclassified*(ntrials_tot/2)/E_cor];% Correcting for the bias
end

% Correction for zeros and ones in hits and FAs
% https://stats.stackexchange.com/questions/134779/d-prime-with-100-hit-rate-probability-and-0-false-alarm-probability
if any(data(:,1) == 0)
    data(data(:,1) == 0,1) = 0.5;
end
if any(data(:,2) == 0)
    data(data(:,2) == 0,2) = 0.5;
end
if any(data(:,1) == data(:,3))
    data(data(:,1) == data(:,3),1) = data(data(:,1) == data(:,3),3)-0.5;
end
if any(data(:,2) == data(:,4))
    data(data(:,2) == data(:,4),2) = data(data(:,2) == data(:,4),4)-0.5;
end


H = (data(:,4)-data(:,2))./data(:,4);
FA = data(:,1)./data(:,3); 
dprime_raw = (1-norminv(FA,0,1))-(1-norminv(H,0,1));

E_H = (data(:,4)-data(:,6))./data(:,4);
E_FA = data(:,5)./data(:,3); 
dprime_correction = (1-norminv(E_FA,0,1))-(1-norminv(E_H,0,1));

dprime_corrected = dprime_raw-dprime_correction;

hold on;
Lblank = all(uniquestim == 0,2);
Lrg = sign(uniquestim(:,1)) ~= sign(uniquestim(:,2)) & ~Lblank;
Llum = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & ~Lblank;
plot(uniquestim(Lrg,3), dprime_corrected(Lrg),'ro-');
plot(uniquestim(Llum,3), dprime_corrected(Llum),'ko-');
set(gca,'Xscale','log');
xlabel('TF'); ylabel('"non-parametric" d-prime');
