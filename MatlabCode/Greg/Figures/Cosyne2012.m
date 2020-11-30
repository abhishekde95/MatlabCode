% DTSpot reviewer figure
WN = nex2stro(findfile('K092810007.nex'));


framerate = WN.sum.exptParams.framerate;
nstixperside = WN.sum.exptParams.nstixperside;
ntrials = length(WN.sum.absTrialNum);
stimonidx = find(strcmp(WN.sum.trialFields(1,:),'stim_on'));
stimoffidx = find(strcmp(WN.sum.trialFields(1,:),'all_off'));
nframesidx = find(strcmp(WN.sum.trialFields(1,:),'num_frames'));
noisetypeidx = find(strcmp(WN.sum.trialFields(1,:),'noise_type'));
sigmaidxs = strmatch('sigma',WN.sum.trialFields(1,:));

hepidx = find(strcmp(WN.sum.rasterCells(1,:),'AD11'));
vepidx = find(strcmp(WN.sum.rasterCells(1,:),'AD12'));
anlgStartTimeidx = find(strcmp(WN.sum.rasterCells(1,:),'anlgStartTime'));
eyestart_t = [WN.ras{:,anlgStartTimeidx}]';
eyesampperiod = 1/WN.sum.analog.storeRates{1};
gammaTable = WN.sum.exptParams.gamma_table;
gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
gammaTable1 = interp1(linspace(0,255,256),gammaTable,linspace(0,255,65536), 'spline');
invgamma = InvertGamma(gammaTable, 0);

% Reconstructing the M matrix
fundamentals = WN.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = WN.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

% Getting the background rgb/lms
ridx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_r'));
gidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_g'));
bidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_b'));
bkgndRGB = [mode(WN.trial(:,ridx)), mode(WN.trial(:,gidx)), mode(WN.trial(:,bidx))];
bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
bkgndlms = M*bkgndrgb;

spikename = getSpikenum(WN);
spikeidx = find(strcmp(WN.sum.rasterCells(1,:),spikename));
maxT = 9;
Lgunnoise = WN.trial(:,noisetypeidx) == 1;
Lconenoise = WN.trial(:,noisetypeidx) == 2;
if (sum(Lconenoise) == 0)
    whichnoisetype = 1;
    disp('Gun noise only');
elseif (sum(Lgunnoise) == 0)
    whichnoisetype = 2;
    disp('Cone noise only');
else
    querystr = ['Which noisetype? 1=gun (n=',num2str(sum(Lgunnoise)),') 2=cone (n=',num2str(sum(Lconenoise)),'): '];
    whichnoisetype = input(querystr);
end
tmpstro = WN;
if (whichnoisetype == 2) && (any(Lgunnoise))  % Removing gun noise
    disp(['Getting rid of ',num2str(sum(Lgunnoise)),' gun noise trials!']);
    tmpstro.ras(Lgunnoise,:) = [];
    tmpstro.trial(Lgunnoise,:) = [];
elseif (whichnoisetype == 1) && (any(Lconenoise))  % Removing cone noise
    disp(['Getting rid of ',num2str(sum(Lconenoise)),' cone noise trials!']);
    tmpstro.ras(Lconenoise,:) = [];
    tmpstro.trial(Lconenoise,:) = [];
end
out = getWhtnsStats(tmpstro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, spikename);
tmpstro = [];
STAs = out{1};
STCs = out{2};
nspikes = out{3};

% Plotting 
figure;
for i = 1:size(STAs,2)
    tmpSTA = reshape(STAs(:,i),nstixperside.^2,3);
    LMSSTA = inv(M')*tmpSTA';
    LMSSTA = permute(reshape(LMSSTA,[3, nstixperside nstixperside]),[2 3 1]);
    for j = 1:3
        normfact = max(abs(LMSSTA(:)));
        subplot(size(LMSSTA,2),3,(3*(i-1))+j);
        imagesc(LMSSTA(:,:,j)./(normfact*2.1)+.5);
        axis square;
        colormap(gray);
    end
end

% Ploting average
tmpSTA = reshape(mean(STAs(:,[4 5 6]),2),nstixperside.^2,3);
LMSSTA = inv(M')*tmpSTA';
LMSSTA = permute(reshape(LMSSTA,[3, nstixperside nstixperside]),[2 3 1]);
normfact = max(abs(LMSSTA(:)));
titles = {'L-cone weight','M-cone weight','S-cone weight'};
for j = 1:3
    subplot(3,1,j);
    image(255*(LMSSTA(:,:,j)./(normfact)+.5));
    axis square;
    colormap(gray(255))
    title(titles{j});
    set(gca,'xtick',[],'ytick',[]);
end

