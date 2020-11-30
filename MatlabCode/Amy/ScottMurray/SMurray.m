% SMurray with two fixation points
% And backgrounds

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

uniquefpxy = unique([fpx fpy],'rows');
uniquefpx = unique([fpx]);
uniquefpy = unique([fpy]);
uniqueors = unique(stim_or);
uniqueirs = unique(stim_ir);

% Looking at LFPs
LFPsamprate = stro.sum.analog.storeRates{1};  % Assuming all channels are sampled at the same rate
centerfreq = 40; % Hz
cyclespersample = centerfreq./LFPsamprate;
ncyclesper6sigma = 5;  % Controls the bandwidth (higher numbers = tighter BW)
nsamplesper6sigma = ceil(ncyclesper6sigma/cyclespersample);
ncycles = nsamplesper6sigma*cyclespersample;
filtkernel1 = normpdf(linspace(-3,3,nsamplesper6sigma),0,1).*cos(linspace(0,2*pi*ncycles,nsamplesper6sigma));
filtkernel2 = normpdf(linspace(-3,3,nsamplesper6sigma),0,1).*sin(linspace(0,2*pi*ncycles,nsamplesper6sigma));
filtkernel1 = filtkernel1./norm(filtkernel1);
filtkernel2 = filtkernel2./norm(filtkernel2);
LFPstarttimes = [stro.ras{:,strcmp(stro.sum.rasterCells,'anlgStartTime')}];


%determine if using an array or a single electrode
electrode = ismember(stro.sum.rasterCells, 'AD01');
if  any(electrode) %electrode
    array = 0;
    data = nan*ones(ntrials,1);
    numelec = 1;
else %array
    array = 1;
    data = nan*ones(ntrials,32);
    numelec = [1:32];
end

for whichADchan = numelec
    for i = 1:ntrials
        if array == 1
            LFP = stro.ras{i,strcmp(stro.sum.rasterCells,['AD',num2str(16+whichADchan)])};
        else
            LFP = stro.ras{i,strcmp(stro.sum.rasterCells,['AD01'])};
        end
        filteredLFP = conv(LFP, filtkernel1,'same').^2+conv(LFP, filtkernel2,'same').^2;
        LFPtimes = LFPstarttimes(i)+[0:length(filteredLFP)-1]/LFPsamprate; % sec
        data(i,whichADchan) = sum(filteredLFP(LFPtimes > stimon_t(i)+0.05 & LFPtimes < stimoff_t(i)+0.05));
    end
end

tmp = nan(size(uniquefpx,1), length(uniqueors), size(data,2));
for i = 1:size(uniquefpx,1)
    for j = 1:length(uniqueors)
        Lor = stim_or == uniqueors(j);
        Lfp = fpx == uniquefpx(i,1);
        L = Lor & Lfp;
        %tmp(i,j,:) = median(data(L,:)); % near, then far
        %tmp(i,j,:) = mean(data(L,:)); % near, then far
        tmp(i,j,:) = trimmean(data(L,:),30); % near, then far
    end
end

% Plotting raw LFPs
figure;
for i = 1:size(tmp,3)
    subplot(ceil(sqrt(size(tmp,3))),ceil(sqrt(size(tmp,3))),i); hold on;
    plot(uniqueors,tmp(1,:,i),'m.-'); % near
    plot(uniqueors,tmp(2,:,i),'k.-'); % far
end
set(gcf,'Name','raw LFP');

% Plotting Z-scored LFPs
mns = mean(tmp,2);
stds = std(tmp,0,2);
zscores = (tmp-repmat(mns,[1,length(uniqueors),1]))./repmat(stds,[1,length(uniqueors),1]);
figure; axes; hold on;
mnz = mean(zscores,3);
plot(uniqueors,mnz(1,:),'m.-'); % near
plot(uniqueors,mnz(2,:),'k.-'); % far
title([num2str(centerfreq),' Hz']);
legend('near','far');
figure;
for i = 1:size(zscores,3)
    subplot(ceil(sqrt(size(zscores,3))),ceil(sqrt(size(zscores,3))),i); hold on;
    plot(uniqueors,zscores(1,:,i),'m.-'); % near
    plot(uniqueors,zscores(2,:,i),'k.-'); % far
end
set(gcf,'Name','Z scores');

% Plotting optimally scaled LFPs
figure;
for i = 1:size(tmp,3)
    subplot(ceil(sqrt(size(tmp,3))),ceil(sqrt(size(tmp,3))),i); hold on;
    scalefactor = tmp(2,:,i)'\tmp(1,:,i)';
    plot(uniqueors,tmp(1,:,i),'m.-'); % near
    plot(uniqueors,scalefactor*tmp(2,:,i),'k.-'); % far
end
set(gcf,'Name','LS scaling');


%%
%{
stro = nex2stro(findfile('A111811001'));
stimx = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_x'));
stimy = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_y'));
stim_ir = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_ir'));
stim_or = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_or'));
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
spikeidx = find(strcmp(stro.sum.rasterCells(1,:),getSpikenum(stro)));

spikerates = [];
for i = 1:size(stro.trial,1)
    spiketimes = stro.ras{i,spikeidx};
    nspikes = sum(spiketimes > stimon_t(i) & spiketimes < stimoff_t(i));
    spikerates = [spikerates; nspikes./(stimoff_t(i)-stimon_t(i))];
end

uniqueirs = unique(stim_ir);
data = [];
for i = 1:length(uniqueirs)
    L = stim_ir == uniqueirs(i);
    data = [data; uniqueirs(i)/10 mean(spikerates(L)) std(spikerates(L)) sum(L)];
end
%figure;
hold on;
h = errorbar(data(:,1),data(:,2),data(:,3)./sqrt(data(:,4)),'k.-','linewidth',2)
ylabel('response (sp/sec)');
xlabel('Annulus radius (deg)');
title(stro.sum.fileName)
set(gca,'Ylim',[min(data(:,2)), max(10*ceil(data(:,2)/10))]);
%}

%%

% --------------------------------------------------------------------
% From here down all scripts are for analyzingdata from SMURRAY1 
% which is GABORS WITH FLANKERS.
% --------------------------------------------------------------------

%%
%{
% Code for looking at data collected with the SMurray1 paradigm (flankers)
% Configurations: 
% 0 = target only
% 1 = flankers only (pref orientation), 
% 2 = flankers only (orthogonal to pref), 
% 3 = target and flankers (pref orient), 
% 4 = target and flankers (flankers orthogonal orient).
% (Conditions 5 - 10 contain flankers and "flankers of flankers")
% 5 = flankers only (pref orientation)
% 6 = flankers only (orthogonal to pref)
% 7 = flankers only (inner: orthogonal to pref, outer: pref orientation)
% 8 = target and flankers (pref orientation)
% 9 = target and flankers (orthogonal to pref)
%10 = target and flankers (inner: orthogonal to pref, outer: pref
%orientation)

% Best cell K040611004 (6/1/11)
% K051311003
stro = nex2stro;
config = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'config'));
configs = unique(config)';
flankerson_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flankers_on'));
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));

spikename = getSpikenum(stro);
spikeidx = find(strcmp(stro.sum.rasterCells(1,:),spikename));
spikes = stro.ras(:,spikeidx);
% rasters
figure; set(gcf,'Name',stro.sum.fileName(end-13:end));
for i = 1:length(configs)
    subplot(length(configs)+1,1,i+1); hold on;
    L = logical(config == configs(i));
    trlidxs = find(L);
    for counter = 1:sum(L)
        trlidx = trlidxs(counter);
        plot(0,counter,'g.');
        plot(stimoff_t(trlidx)-stimon_t(trlidx),counter,'r.');
%        plot(flankerson_t(trlidx)-stimon_t(trlidx),counter,'m*');
        
        Lplotspike = logical(spikes{trlidx}-stimon_t(trlidx) > -0.5);
        nspikestot = sum(Lplotspike);
        plot([spikes{trlidx}(Lplotspike) spikes{trlidx}(Lplotspike)]'-stimon_t(trlidx),[zeros(nspikestot,1) .5*ones(nspikestot,1)]'+counter,'k-');
    end
    set(gca,'Xlim',[-.1 .3]);
end

% psths
figure; set(gcf,'Name',stro.sum.fileName(end-13:end));
bins = linspace(-.05,.35,100);
binwidth = bins(2)-bins(1);
for i = 1:length(configs)
    subplot(length(configs)+1,1,i+1); hold on;
    L = logical(config == configs(i));
    trlidxs = find(L);
    psth = nan*ones(sum(L),length(bins));
    for counter = 1:sum(L)
        trlidx = trlidxs(counter);
        Lplotspike = logical(spikes{trlidx}-stimon_t(trlidx) > -0.5);
        psth(counter,:) = hist(spikes{trlidx}(Lplotspike)-stimon_t(trlidx),bins);
    end
    h = bar(bins(2:end-1),mean(psth(:,2:end-1))/binwidth);
    set(gca,'Xlim',[bins(2) bins(end-1)],'Ylim',[0 500]);
end

%barplot
figure; set(gcf,'Name',stro.sum.fileName(end-13:end)); axes; hold on;
for i = unique(config)'
    L = logical(config == i);
    trlidxs = find(L);
    tmp = nan(sum(L),1);
    for counter = 1:sum(L)
        st = spikes{trlidxs(counter)}-stimon_t(trlidxs(counter));
        tmp(counter) = sum(st > 0 & st < .2);
    end
    se = std(tmp)./sqrt(sum(L));
    bar(i,mean(tmp));
    plot([i i],[-se se]+mean(tmp),'k-')
end
%}
%% 
%{
% Mini-population plot for SMurray1 paradigm (Gabor with flankers ,but not
% flankers-of-flankers)
filenames = {'K030411004','K021711002','S020111002','S012111003','S011711005','S011711001'};
data = [];
for i = 1:length(filenames)
    stro = nex2stro(findfile(filenames{i}));
    config = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'config'));
    flankerson_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flankers_on'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
    
    spikename = getSpikenum(stro);
    spikeidx = find(strcmp(stro.sum.rasterCells(1,:),spikename));
    spikes = stro.ras(:,spikeidx);
    for j = 1:2
        L = logical(config == j+2);
        trlidxs = find(L);
        tmp = nan(sum(L),1);
        for counter = 1:sum(L)
            st = spikes{trlidxs(counter)}-stimon_t(trlidxs(counter));
            tmp(counter) = sum(st > 0 & st < .2);
        end
        if (j == 1)
            data(i,1:2) = [mean(tmp) std(tmp)./sqrt(sum(L))];
        else
            data(i,3:4) = [mean(tmp) std(tmp)./sqrt(sum(L))];
        end
    end
end

figure; axes; hold on;
for i = 1:length(filenames)
    plot(data(i,1),data(i,3),'ko','MarkerFaceColor','black','MarkerSize',5);
    plot(data(i,1)+data(i,2)*[-1 1],[data(i,3) data(i,3)],'k-');
    plot([data(i,1) data(i,1)],data(i,3)+data(i,4)*[-1 1],'k-');    
end
set(gca,'XScale','log','YScale','log');
plot([1 100],[1 100],'k-');
xlabel('Flankers/targets same');
ylabel('Flankers/targets orthogonal');
axis square;
[h,p] = ttest(data(:,1)-data(:,3))
%}
%%
%{
% Mini-population plot for SMurray1 paradigm (including
% flankers-of-flankers)
% Conditions
% 8:  VVVVV
% 9:  HHVHH
% 10: VHVHV
% Is the response to 9 greater than the response to 10?

% All relevant files

filenames = {'K033011003','K033011004','K033111004','K033111005','K033111006','K040411002','K040411003',...
    'K040511002','K040511003','K040511005','K040611002','K040611004','K040611005','K041111002',...
    'K041111003','K041111004'};

% One file per cell

ps = [];
data = [];
for i = 1:length(filenames)
    stro = nex2stro(findfile(filenames{i}));
    config = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'config'));
    flankerson_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flankers_on'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
    
    spikename = getSpikenum(stro);
    spikeidx = find(strcmp(stro.sum.rasterCells(1,:),spikename));
    spikes = stro.ras(:,spikeidx);
    tmp = [];
    for j = 1:3
        L = logical(config == j+7);
        trlidxs = find(L);
        for counter = 1:sum(L)
            st = spikes{trlidxs(counter)}-stimon_t(trlidxs(counter));
            tmp(j,counter) = sum(st > 0 & st < .2);
        end
        if (j == 1)
            data(i,1:2) = [mean(tmp(1,:)) std(tmp(1,:))./sqrt(sum(L))];
        elseif (j == 2)
            data(i,3:4) = [mean(tmp(2,:)) std(tmp(2,:))./sqrt(sum(L))];
        elseif (j == 3)
            data(i,5:6) = [mean(tmp(3,:)) std(tmp(3,:))./sqrt(sum(L))];
        end
    end
    
    % significance tests on individual cells
    p1 = ranksum(tmp(1,:), tmp(2,:))
    p2 = ranksum(tmp(2,:), tmp(3,:))
    ps = [ps; p1 p2]
end

% VVVVV vs HHVHH
L = data(:,1) > 0 & data(:,3) > 0;  % min spike count
Lsig = ps(:,1) < 0.05;
figure; axes; hold on;
plot(data(L&Lsig,3),data(L&Lsig,5),'k*');
plot(data(L,1),data(L,3),'k.');
plot([0 25],[0 25],'k-');
[h,p] = ttest(data(L,1)-data(L,3));
title(['p = ',num2str(p)])
xlabel('VVVVV'); ylabel('HHVHH');

% HHVHH vs VHVHV
L = data(:,3) > 0 & data(:,5) > 0;  % min spike count
Lsig = ps(:,2) < 0.05;
figure; axes; hold on;
plot(data(L&Lsig,3),data(L&Lsig,5),'k*');
plot(data(L,3),data(L,5),'m.');
plot([0 25],[0 25],'k-');
[h,p] = ttest(data(L,3)-data(L,5));
title(['p = ',num2str(p)])
xlabel('HHVHH'); ylabel('VHVHV');

%}
%%
%{
% making a .mat file for Scott so that he can make a figure for the grant


% All relevant files
filenames = {'K030411004','K021711002','S020111002','S012111003','S011711005','S011711001',...
    'K033011003','K033011004','K033111004','K033111005','K033111006','K040411002','K040411003',...
    'K040511002','K040511003','K040511005','K040611002','K040611004','K040611005','K041111002',...
    'K041111003','K041111004'};
countwin = [.06 .2]; % time in ms re stimon to count spikes
data = [];
for i = 1:length(filenames)
    stro = nex2stro(findfile(filenames{i}));
    config = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'config'));
    fix_t =  stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_acq'));
    flankerson_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flankers_on'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
    
    spikename = getSpikenum(stro);
    spikeidx = find(strcmp(stro.sum.rasterCells(1,:),spikename));
    spikes = stro.ras(:,spikeidx);

    baselines = [];
    for j = 1:size(stro.trial,1)
        spiketimes = stro.ras{j,spikeidx};
        nspikes = sum(spiketimes > fix_t(j) & spiketimes < flankerson_t(j));
        baselines = [baselines; nspikes./(flankerson_t(j)-fix_t(j))];
    end
    
    tmp = [];
    for j = unique(config)'
        L = logical(config == j);
        trlidxs = find(L);
        for counter = 1:sum(L)
            st = spikes{trlidxs(counter)}-stimon_t(trlidxs(counter));
            tmp(j==unique(config),counter) = sum(st > countwin(1) & st < countwin(2))./(sum(L)*countwin(2)-countwin(1));
            st = spikes{trlidxs(counter)}-fix_t(trlidxs(counter));
        end
    end
    
    data(i).baselinerate = mean(baselines);
    data(i).configs = unique(config);
    data(i).spikerate = tmp;
    data(i).filename = filenames(i);
    data(i).parameters = stro.sum.exptParams;
end
%}
%%
%{
% As above, but making a file containing the entire stro structures

[fnames, spikeIdx] = fnamesFromTxt2();
data = cell(length(fnames),3);
for cellcounter = 1:size(fnames,1)
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        paradigmID = getparadigmID(filename);
        if (paradigmID ~= 104)
            error('Wrong paradigm');
        end
        data{cellcounter,i} = nex2stro(filename);
    end
end
%}