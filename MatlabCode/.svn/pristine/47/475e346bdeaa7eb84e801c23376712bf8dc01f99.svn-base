% Some additional analysis (eye movement, RF size) after the JoV rejection
% Author - Abhishek De, 08/2020


close all; clearvars;
plot_counter = 1;

%% Description of analyses

% Section 1: RF sizes as a function of eccentricity

% Section 2: Pulling out the example STA filenames for which Greg might have Gratings data

% Section 3: Eye movement analyis: Goal is to check if the eye movement statistics are the same for spatially opponent and spatially-non-opponent cells

% Section 4: Eye movement analysis for population of cells

% Section 5: Some from Greg

% Section 6: R as a function of SNR 

% Section 7: Greg's approach for multiclass linear regression

% Section 8: Checking whether DO and simple cell RF locations differed systematically

%% Section 1: RF sizes as a function of eccentricity 
% Check how the sqrt(RF area) varies with eccentricity 

if ~exist('plot_counter')
    plot_counter = 1;
end
% load Output_ListWN.mat
load Output_ListWN2.mat
load Singleopponent.mat
load Ratio_of_power.mat

crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
NLI = cell2mat(Output_List(:,13));
simplecells = NLI<0;
Singleopponent = logical(Singleopponent);

% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

% calculating the cone weights
conewts_svd = cell2mat(Output_List(~Singleopponent & simplecells,23)');
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

% Next, I am plotting the SVD derived luminance, color-opponent and other cells
ind = find(~Singleopponent & simplecells);
Lumind = ind(LumIds_conewts);
COind = ind(ColorOpponentIds_conewts); 
Sconeind = ind(Sconedominated_conewts);

idxs = [Lumind; COind; Sconeind];

ecc = [];
RFarea = [];
for aa = 1:numel(idxs)
    ii = idxs(aa);
    % The signficant RFs are stored in column 18
    RFarea = [RFarea; sum(Output_List{ii,18}(:))*0.04];
    
    % RF eccentricity is stored in column 8
    tmpecc = cell2mat(Output_List(ii,8));
    ecc = [ecc; sqrt(sum(tmpecc.^2,2))/10];
end

% plotting the results 
figure(plot_counter);
plot(ecc(ecc<10), sqrt(RFarea(ecc<10)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out'); xlabel('Ecc'); ylabel('RF field size');
plot_counter = plot_counter + 1;

%% Section 2: Pulling out the example STA filenames for which Greg might have Gratings data

if ~exist('plot_counter')
    plot_counter = 1;
end
% load Output_ListWN.mat
load Output_ListWN2.mat
load Singleopponent.mat
load modelfits.mat
load Ratio_of_power.mat
crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
NLI = cell2mat(Output_List(:,13));
simplecells = NLI<0;
Singleopponent = logical(Singleopponent);
ind = find(~Singleopponent & simplecells);

% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

% calculating the cone weights
conewts_svd = cell2mat(Output_List(~Singleopponent & simplecells,23)');
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

% Next, I am plotting the SVD derived luminance, color-opponent and other cells
Lumind = ind(LumIds_conewts); 
COind = ind(ColorOpponentIds_conewts); 
Sconeind = ind(Sconedominated_conewts);

% idx = [Lumind(8); COind(30); Sconeind(20)];
idx = [Lumind(18); Lumind(2); COind(23); COind(18); Sconeind(22); Sconeind(11)];
filesofinterest = Output_List(idx,1);

%% Section 3: Eye movement analyis: Goal is to check if the eye movement statistics are the same for 
% spatially opponent and spatially-non-opponent cells
close all;
if ~exist('plot_counter')
    plot_counter = 1;
end
% load Output_ListWN.mat
load Output_ListWN2.mat
load Singleopponent.mat
load modelfits.mat
load Ratio_of_power.mat
crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
NLI = cell2mat(Output_List(:,13));
simplecells = NLI<0;
Singleopponent = logical(Singleopponent);
ind = find(~Singleopponent & simplecells);

% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

% calculating the cone weights
conewts_svd = cell2mat(Output_List(~Singleopponent & simplecells,23)');
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

% Next, I am plotting the SVD derived luminance, color-opponent and other cells
Lumind = ind(LumIds_conewts); 
COind = ind(ColorOpponentIds_conewts); 
Sconeind = ind(Sconedominated_conewts);

idx = [Lumind(18); Lumind(2); COind(23); COind(18); Sconeind(22); Sconeind(11)];

% Let's start with the analysis of eye movements for all these example files
% The good thing is that all these files were collected from monkey K

count = 1;
color_counter = 1;
c = [0 0 0; 1 0 0; 0 0.5 1.0];
for jj = 1:numel(idx)
    
    if ~contains(Output_List{idx(jj),1},'.nex')
        WN = nex2stro(findfile(strcat(Output_List{idx(jj),1},'.nex')));
    else
        WN = nex2stro(findfile(Output_List{idx(jj),1},'.nex'));
    end
    
    framerate = WN.sum.exptParams.framerate;
    stimonidx = find(strcmp(WN.sum.trialFields(1,:),'stim_on'));
    stimoffidx = find(strcmp(WN.sum.trialFields(1,:),'all_off'));
    fpacqidx = find(strcmp(WN.sum.trialFields(1,:),'fpacq'));
    
    noisetypeidx = find(strcmp(WN.sum.trialFields(1,:),'noise_type'));
    
    hepidx = find(strcmp(WN.sum.rasterCells(1,:),'AD11'));
    vepidx = find(strcmp(WN.sum.rasterCells(1,:),'AD12'));
    
    anlgStartTimeidx = find(strcmp(WN.sum.rasterCells(1,:),'anlgStartTime'));
    L = WN.trial(:,noisetypeidx)==1;
    WN.ras(~L ,:) = []; % modiftying the WN structure
    WN.trial(~L,:) = []; % modiftying the WN structure
    
    samplingrate = WN.sum.analog.storeRates{1};
    H = [];
    V = [];
    for ii = 1:size(WN.trial,1)
        
        anlgStarttime = WN.ras{ii,anlgStartTimeidx};
        stimont = WN.trial(ii,stimonidx);
        stimofft = WN.trial(ii,stimoffidx);
        
        h_eyepos = WN.ras{ii,hepidx}*4096/400;
        v_eyepos = WN.ras{ii,vepidx}*4096/400;
        t = anlgStarttime: (1/samplingrate): anlgStarttime+(numel(h_eyepos)-1)*(1/samplingrate);
        tmp = t>=stimont & t<=stimofft;
        
        % Storing the H and V position
        H = [H; h_eyepos(tmp)];
        V = [V; v_eyepos(tmp)];
        
    end
    
    % Storing the amplitude of the eye position
    Amp = sqrt((H-mean(H)).^2 + (V-mean(V)).^2);
    
    % Plotting the polar histogram
    [THETA,~] = cart2pol(H-mean(H),V-mean(V));
    figure(plot_counter); subplot(2,6,count);
    polarhistogram(THETA,0:2*pi/50:2*pi,'FaceColor',c(color_counter,:),'EdgeColor',[1 1 1]);
    hold off;
    
    % Histogram of amplitudes
    figure(plot_counter); subplot(2,6,6+count);
    histogram(Amp,0:0.05:2.0,'FaceColor',c(color_counter,:), 'EdgeColor', [1 1 1]);
    axis square; set(gca,'Tickdir','out','Xlim',[0 1], 'XTick',0:0.25:1.0); xlabel('Amp'); ylabel('Count'); hold off;
    
    count = count + 1;
    if mod(count,2)
        color_counter = color_counter + 1;
    end
    disp(median(Amp))
end
plot_counter = plot_counter + 1;

%% Section 4: Eye movement analysis for population of cells

close all; 
if ~exist('plot_counter')
    plot_counter = 1;
end
% load Output_ListWN.mat
load Output_ListWN2.mat
load Singleopponent.mat
load modelfits.mat
load Ratio_of_power.mat
crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
NLI = cell2mat(Output_List(:,13));
simplecells = NLI<0;
Singleopponent = logical(Singleopponent);
ind = find(~Singleopponent & simplecells);

% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

% calculating the cone weights
conewts_svd = cell2mat(Output_List(~Singleopponent & simplecells,23)');
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

% Next, I am plotting the SVD derived luminance, color-opponent and other cells
Lumind = ind(LumIds_conewts); 
COind = ind(ColorOpponentIds_conewts); 
Sconeind = ind(Sconedominated_conewts);

% Pulling the eye movements
eyeAmp = Output_List(:,29);
load Peraccuracy.mat
load SSE.mat
load Deviation.mat
meanperaccuracy = zeros(size(Peraccuracy));
meanSSE = zeros(size(SSE));
meanR = zeros(size(SSE));
for ii = 1:size(Peraccuracy,1)
    for jj = 1:size(Peraccuracy,2)
        meanperaccuracy(ii,jj) = mean(Peraccuracy{ii,jj});
        meanSSE(ii,jj) = mean(SSE{ii,jj});
        meanR(ii,jj) = mean(cos(Deviation{ii,jj}*pi/180));
    end
end
meanperaccuracy = meanperaccuracy(~Singleopponent & simplecells,:);
meanSSE = meanSSE(~Singleopponent & simplecells,:);
meanR = meanR(~Singleopponent & simplecells,:);

% Calculating the spatial opponency index  
SPATIALRF = Output_List(:,4);
SOI_index_svd = []; % Spatial opponency index using SVD spatial RF 
for ii=1:size(SPATIALRF,1)
    % Based on SVD derived maps 
    tmp_svd = SPATIALRF{ii};
    P_svd = abs(sum(sum(tmp_svd(tmp_svd>0)))); N_svd = abs(sum(sum(tmp_svd(tmp_svd<0))));
    SOI_index_svd = [SOI_index_svd; 1-abs((P_svd-N_svd)/(P_svd+N_svd))];
    
end

% Plotting the population results
figure(plot_counter);
Amp = cell2mat(Output_List(:,29));
subplot(311); 
plot(Amp(Lumind), meanR(LumIds_conewts,2)-meanR(LumIds_conewts,3), 'o', 'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(Amp(COind), meanR(ColorOpponentIds_conewts,2)-meanR(ColorOpponentIds_conewts,3), 'o', 'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(Amp(Sconeind), meanR(Sconedominated_conewts,2)-meanR(Sconedominated_conewts,3), 'o', 'MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Xlim',[0.0 0.40], 'Ylim', [-0.5 0.5]); 
xlabel('median eye displacement (in degrees)'); ylabel('R(Gabor)-R(DoG)');

subplot(312); 
plot(Amp(Lumind), meanR(LumIds_conewts,2)-meanR(LumIds_conewts,1), 'o', 'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(Amp(COind), meanR(ColorOpponentIds_conewts,2)-meanR(ColorOpponentIds_conewts,1), 'o', 'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(Amp(Sconeind), meanR(Sconedominated_conewts,2)-meanR(Sconedominated_conewts,1), 'o', 'MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Xlim',[0.0 0.40], 'Ylim', [-0.5 0.5]); 
xlabel('median eye displacement (in degrees)'); ylabel('R(Gabor)-R(modified DoG)');

% Plotting the relationship between eye movements and spatial opponency index 
subplot(313); 
plot(Amp(Lumind), SOI_index_svd(Lumind), 'o', 'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(Amp(COind), SOI_index_svd(COind), 'o', 'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(Amp(Sconeind), SOI_index_svd(Sconeind), 'o', 'MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Xlim',[0.0 0.40], 'Ylim', [0.3 1.0]); 
xlabel('median eye displacement (in degrees)'); ylabel('Spatial Opponency Index');
plot_counter = plot_counter + 1;


% Doing some stats for eye movement amplitude vs. R(Gabor)-R(DoG)
combined_ind = [LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts];
[r1,p1] = corr(meanR(LumIds_conewts,2)-meanR(LumIds_conewts,3), Amp(Lumind),'type','Spearman');
[r2,p2] = corr(meanR(ColorOpponentIds_conewts,2)-meanR(ColorOpponentIds_conewts,3), Amp(COind),'type','Spearman','rows','complete');
[r3,p3] = corr(meanR(Sconedominated_conewts,2)-meanR(Sconedominated_conewts,3), Amp(Sconeind),'type','Spearman','rows','complete');
[r_combined,p_combined] = corr(meanR(combined_ind,2)-meanR(combined_ind,3), Amp([Lumind; COind; Sconeind]),'type','Spearman','rows','complete');

% Eye movement amplitude vs. R(Gabor)-R(modified DoG)
[r4,p4] = corr(meanR(LumIds_conewts,2)-meanR(LumIds_conewts,1), Amp(Lumind),'type','Spearman');
[r5,p5] = corr(meanR(ColorOpponentIds_conewts,2)-meanR(ColorOpponentIds_conewts,1), Amp(COind),'type','Spearman','rows','complete');
[r6,p6] = corr(meanR(Sconedominated_conewts,2)-meanR(Sconedominated_conewts,1), Amp(Sconeind),'type','Spearman','rows','complete');
[r_combined2, p_combined2] = corr(meanR(combined_ind,2)-meanR(combined_ind,1), Amp([Lumind; COind; Sconeind]),'type','Spearman','rows','complete');

% Eye movement amplitude vs. SOI index
[r7,p7] = corr(SOI_index_svd(Lumind), Amp(Lumind),'type','Spearman');
[r8,p8] = corr(SOI_index_svd(COind), Amp(COind),'type','Spearman','rows','complete');
[r9,p9] = corr(SOI_index_svd(Sconeind), Amp(Sconeind),'type','Spearman','rows','complete');
[r_combined3, p_combined3] = corr(SOI_index_svd([Lumind; COind; Sconeind]), Amp([Lumind; COind; Sconeind]),'type','Spearman','rows','complete');


%% Section 5: Some from Greg 
% Eye movement sample 95% ellipses 

filenames = {'K040708003.nex','M121417002.nex','K110408003.nex','M082217003.nex','K052308001.nex','S071510004.nex'};
figure; axes; hold on;
for fileidx = 1:length(filenames)
    if filenames{fileidx}(1) == 'K'
        minipath = ['Greg',filesep,'Kali',filesep,'2008'];
    elseif filenames{fileidx}(1) == 'S'
        minipath = ['Greg',filesep,'Sedna',filesep,'2010'];
    elseif filenames{fileidx}(1) == 'M'
        minipath = ['Abhishek',filesep,'Maui'];
    else
        minipath = [];
    end
    WN=nex2stro([nexfilepath,filesep,minipath,filesep,filenames{fileidx}]);
    framerate = WN.sum.exptParams.framerate;
    hepidx =strcmp(WN.sum.rasterCells,'AD11');
    vepidx =strcmp(WN.sum.rasterCells,'AD12');
    stimonidx = find(strcmp(WN.sum.trialFields(1,:),'stim_on'));
    nframesidx = find(strcmp(WN.sum.trialFields(1,:),'num_frames'));
    stimon_t = WN.trial(:,stimonidx);
    numframes = WN.trial(:,nframesidx);
    stimoff_t = stimon_t+numframes/framerate;
    stimoff_t = stimon_t+(numframes-9)/framerate; % Eliminating saccades that leave fixwin
    ntrials = length(WN.sum.absTrialNum);
    anlgStartTimeidx = find(strcmp(WN.sum.rasterCells(1,:),'anlgStartTime'));
    eyestart_t = [WN.ras{:,anlgStartTimeidx}]';
    eyesampperiod = 1/WN.sum.analog.storeRates{1};
    data = [];
    for i = 1:ntrials
        hep = WN.ras{i,hepidx}*4096/400;
        vep = WN.ras{i,vepidx}*4096/400;
        nsamps = length(hep);
        t = [0:1:nsamps-1]*eyesampperiod+eyestart_t(i);
        data = [data; hep(t>stimon_t(i) & t < stimoff_t(i)) vep(t>stimon_t(i) & t < stimoff_t(i))];
        %plot(hep(t>stimon_t(i) & t < stimoff_t(i)),vep(t>stimon_t(i) & t < stimoff_t(i)),'r-')
    end
    mn = mean(data);
    sd = std(data);
    S2 = cov(data);
    [v,d] = eig(S2);
    [counts,xedges,yedges] = histcounts2(data(:,1),data(:,2),20);
    x = linspace(xedges(1),xedges(end),size(counts,2));
    y = linspace(yedges(1),yedges(end),size(counts,1));
    prob = counts./sum(counts(:)); % probability
    n = size(data,1);
    % crit = sqrt((2*(n-1)/(n-2))*finv(.95,2,n-2)); % <--- Hotelling's T^2 (since we're estimating S2)
    crit = sqrt(chi2inv(.95,2));
    tmp = [cos(linspace(0,2*pi,100))' sin(linspace(0,2*pi,100))'];
    xy = tmp*crit*sqrt(d)*v;
    %plot(xy(:,1)+mn(1),xy(:,2)+mn(2),'b-');
    plot(xy(:,1),xy(:,2),'r-'); % Assuming same mean for all data sets, and why not?
    %contourvals = chi2pdf(crit.^2,2)./chi2pdf(0,2)
    %countourvals = max(prob(:)./10);
    %contour(repmat(x',1,length(y)),repmat(y,length(x),1),prob,[countourvals countourvals]);
    set(gca,'Xlim',[-1 1]);
    set(gca,'Ylim',[-1 1]);
end

%% Section 6: Calculating the fraction of explained variance (Ringach, 2002)

close all; 
if ~exist('plot_counter')
    plot_counter = 1;
end
% load Output_ListWN.mat
load Output_ListWN2.mat
load Singleopponent.mat
load Ratio_of_power.mat
crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
NLI = cell2mat(Output_List(:,13));
simplecells = NLI<0;
Singleopponent = logical(Singleopponent);
ind = find(~Singleopponent & simplecells);

% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

% calculating the cone weights
conewts_svd = cell2mat(Output_List(~Singleopponent & simplecells,23)');
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

% Next, I am plotting the SVD derived luminance, color-opponent and other cells
Lumind = ind(LumIds_conewts); 
COind = ind(ColorOpponentIds_conewts); 
Sconeind = ind(Sconedominated_conewts);

% Loading R_square and Pearson R for the entire dataset
load R_square.mat
load Pearson_R.mat

% Calculating the fraction of unexplained variance for Gabor and DOG models
FUV_Gabor = 1-R_square([Lumind; COind; Sconeind],2);
FUV_DOG = 1-R_square([Lumind; COind; Sconeind],3);

% Calculating the Pearson correlation coefficient for Gabor and DOG models
Pearson_Gabor = abs(cos(Pearson_R([Lumind; COind; Sconeind],2)));
Pearson_DOG = abs(cos(Pearson_R([Lumind; COind; Sconeind],3)));

figure(plot_counter);
subplot(221); histogram(FUV_Gabor, 0:0.1:1.0,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
xlabel('Fraction of unexplained variance'); ylabel('Count'); axis square; title('Gabor model');
set(gca,'Tickdir','out','Xlim',[0 1],'XTick',0:0.2:1.0);
subplot(222); histogram(FUV_DOG, 0:0.1:1.0,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
xlabel('Fraction of unexplained variance'); ylabel('Count'); axis square; title('DOG model');
set(gca,'Tickdir','out','Xlim',[0 1],'XTick',0:0.2:1.0);
subplot(223); histogram(Pearson_Gabor, 0:0.1:1.0,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
xlabel('Pearson R'); ylabel('Count'); axis square; title('Gabor model');
set(gca,'Tickdir','out','Xlim',[0 1],'XTick',0:0.2:1.0);
subplot(224); histogram(Pearson_DOG, 0:0.1:1.0,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
xlabel('Pearson R'); ylabel('Count'); axis square; title('DOG model');
set(gca,'Tickdir','out','Xlim',[0 1],'XTick',0:0.2:1.0);
plot_counter = plot_counter + 1;

%% Section 6: R as a function of SNR and comparing that relationship for DO and simple cells
% Based on Greg's suggestion, I need to Z-transform the R and fit a regression model

if ~exist('plot_counter')
    plot_counter = 1;
end
% load Output_ListWN.mat
load Output_ListWN2.mat
load Singleopponent.mat
load modelfits.mat
load Ratio_of_power.mat
crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;

Output_List(~Z_cellsofinterest,:) = [];
NLI = cell2mat(Output_List(:,13));
simplecells = NLI<0;
Singleopponent = logical(Singleopponent);
% Singleopponent = logical(Ratio_of_power<1.2);

SNR = Zmax(Z_cellsofinterest)/300; % Dividing by the number of elements
SNR = SNR(~Singleopponent & simplecells);

% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

% calculating the cone weights
conewts_svd = cell2mat(Output_List(~Singleopponent & simplecells,23)');
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

load Peraccuracy.mat
load SSE.mat
load Deviation.mat
meanperaccuracy = zeros(size(Peraccuracy));
meanSSE = zeros(size(SSE));
meanR = zeros(size(SSE));
for ii = 1:size(Peraccuracy,1)
    for jj = 1:size(Peraccuracy,2)
        meanperaccuracy(ii,jj) = mean(Peraccuracy{ii,jj});
        meanSSE(ii,jj) = mean(SSE{ii,jj});
        meanR(ii,jj) = mean(cos(Deviation{ii,jj}*pi/180));
    end
end
meanperaccuracy = meanperaccuracy(~Singleopponent & simplecells,:);
meanSSE = meanSSE(~Singleopponent & simplecells,:);
meanR = meanR(~Singleopponent & simplecells,:);

bins = 0:0.1:1.2;
figure(plot_counter); set(gcf,'Name','R=f(SNR)');
subplot(121); hold on;
plot(SNR(ColorOpponentIds_conewts),atanh(meanR(ColorOpponentIds_conewts,2)),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(SNR(Sconedominated_conewts),atanh(meanR(Sconedominated_conewts,2)),'o','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]);
plot(SNR(LumIds_conewts),atanh(meanR(LumIds_conewts,2)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
set(gca,'Tickdir','out','Xlim',[1.0 16],'Ylim',[0 1.6],'XTick',[1 2 4 8 16],'YTick',0:0.8:1.6,'XScale','log'); axis square; xlabel('log10(SNR)'); ylabel('Z(R)'); title('Gabor'); hold off;
subplot(122); hold on;
plot(SNR(ColorOpponentIds_conewts),atanh(meanR(ColorOpponentIds_conewts,3)),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(SNR(Sconedominated_conewts),atanh(meanR(Sconedominated_conewts,3)),'o','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]);
plot(SNR(LumIds_conewts),atanh(meanR(LumIds_conewts,3)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[1.0 16],'Ylim',[0 1.6],'XTick',[1 2 4 8 16],'YTick',0:0.8:1.6,'XScale','log'); axis square; xlabel('log10(SNR)'); ylabel('Z(R)'); title('DoG'); hold off;
plot_counter = plot_counter + 1;

% Now have to use regression
DOIds_conewts = [ColorOpponentIds_conewts Sconedominated_conewts];

% Z-transforming the R values
data1 = [SNR(LumIds_conewts); SNR(ColorOpponentIds_conewts); SNR(Sconedominated_conewts)];
data1log = log10(data1);
data2 = atanh([meanR(LumIds_conewts,2); meanR(ColorOpponentIds_conewts,2); meanR(Sconedominated_conewts,2)]); % Gabor 
data3 = atanh([meanR(LumIds_conewts,3); meanR(ColorOpponentIds_conewts,3); meanR(Sconedominated_conewts,3)]); % DoG
group = [ones(size(SNR(LumIds_conewts))); 2*ones(size(SNR(ColorOpponentIds_conewts))) ; 3*ones(size(SNR(Sconedominated_conewts)))];

% Need to quantify the the differences in DoG as a function of SNR:
celltype_dv = zeros(numel(group),2);
for ii = 1:size(celltype_dv,1)
    celltype_dv(ii,group(ii)) = 1;
end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
colors = {'black','red','cyan'};

numcells = size(celltype_dv,2);
 
X = [repmat(data1log,1,numcells).*celltype_dv celltype_dv]; % Individual model
[b, bint, r_Gabor,~, stats] = regress(data2,X);
[b_full, bint_full, r_full_Gabor] = regress(data2,[data1log ones(size(data1log))]);

% Trying linear regression
figure(plot_counter);
subplot(121); plot(SNR(LumIds_conewts),atanh(meanR(LumIds_conewts,2)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(SNR(ColorOpponentIds_conewts),atanh(meanR(ColorOpponentIds_conewts,2)),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(SNR(Sconedominated_conewts),atanh(meanR(Sconedominated_conewts,2)),'o','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[1 16],'Ylim',[0 1.6],'XTick',[1 2 4 8 16],'YTick',0:0.8:1.6); 
lims = [1 2 4 8 16];

for i = 1:numcells
    tmp = [b(i) b(i+numcells)]*[log10(lims); ones(size(lims))];
    plot(lims,tmp,'color',colors{i}); hold on;
end
plot(lims,[b_full(1) b_full(2)]*[log10(lims); ones(size(lims))],'g', 'Linewidth', 2); hold on;
axis square; xlabel('SNR'); ylabel('Z(R)'); title('Gabor'); set(gca,'XScale','log'); 
hold off;

subplot(122); plot(SNR(LumIds_conewts),atanh(meanR(LumIds_conewts,3)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(SNR(ColorOpponentIds_conewts),atanh(meanR(ColorOpponentIds_conewts,3)),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(SNR(Sconedominated_conewts),atanh(meanR(Sconedominated_conewts,3)),'o','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[1 16],'Ylim',[0 1.6],'XTick',[1 2 4 8 16],'YTick',0:0.8:1.6,'XScale','log'); 

X = [repmat(data1log,1,numcells).*celltype_dv celltype_dv];
[b, bint, r_DoG,~, stats1_DoG] = regress(data3,X);
[b_full, bint_full, r_full_DoG,~,stats0_DoG] = regress(data3,[data1log ones(size(data1log))]);
for i = 1:numcells
    tmp = [b(i) b(i+numcells)]*[log10(lims); ones(size(lims))];
    plot(lims,tmp,'color',colors{i}); hold on;
end
plot(lims,[b_full(1) b_full(2)]*[log10(lims); ones(size(lims))],'g', 'Linewidth', 2); hold on;
axis square; xlabel('SNR'); ylabel('Z(R)'); title('DoG'); hold off;
plot_counter = plot_counter + 1;

% Doing some stats: F test
dp = numcells*2-2;
RSSE_Gabor = sum(r_Gabor.^2);
RSSE_full_Gabor = sum(r_full_Gabor.^2);
Fstat_Gabor = ((RSSE_full_Gabor - RSSE_Gabor)/dp)/(RSSE_Gabor/(numel(data1)-numcells*2));
p_Gabor = 1-fcdf(Fstat_Gabor,dp,(numel(data1)-numcells*2));

RSSE_DoG = sum(r_DoG.^2);
RSSE_full_DoG = sum(r_full_DoG.^2);
Fstat_DoG = ((RSSE_full_DoG - RSSE_DoG)/dp)/(RSSE_DoG/(numel(data1)-numcells*2));
p_DoG = 1-fcdf(Fstat_DoG,dp,(numel(data1)-numcells*2));

%% Combining the LM and S cells in a DO cell category 
% Now have to use regression
DOIds_conewts = [ColorOpponentIds_conewts Sconedominated_conewts];

% Z-transforming the R values
data1 = [SNR(LumIds_conewts); SNR(DOIds_conewts)];
data1log = log10(data1);
data2 = atanh([meanR(LumIds_conewts,2); meanR(DOIds_conewts,2)]); % Gabor 
data3 = atanh([meanR(LumIds_conewts,3); meanR(DOIds_conewts,3)]); % DoG
group = [ones(size(SNR(LumIds_conewts))); 2*ones(size(SNR(DOIds_conewts)))];

% Need to quantify the the differences in DoG as a function of SNR:
celltype_dv = zeros(numel(group),2);
for ii = 1:size(celltype_dv,1)
    celltype_dv(ii,group(ii)) = 1;
end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
colors = {'black','red','cyan'};

numcells = size(celltype_dv,2);
 
X = [repmat(data1log,1,numcells).*celltype_dv celltype_dv]; % Individual model
[b, bint, r_Gabor,~, stats] = regress(data2,X);
[b_full, bint_full, r_full_Gabor] = regress(data2,[data1log ones(size(data1log))]);

% Trying linear regression
figure(plot_counter);
subplot(121); plot(SNR(LumIds_conewts),atanh(meanR(LumIds_conewts,2)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(SNR(DOIds_conewts),atanh(meanR(DOIds_conewts,2)),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[1 16],'Ylim',[0 1.6],'XTick',[1 2 4 8 16],'YTick',0:0.8:1.6); 
lims = [1 2 4 8 16];

for i = 1:numcells
    tmp = [b(i) b(i+numcells)]*[log10(lims); ones(size(lims))];
    plot(lims,tmp,'color',colors{i}); hold on;
end
plot(lims,[b_full(1) b_full(2)]*[log10(lims); ones(size(lims))],'g', 'Linewidth', 2); hold on;
axis square; xlabel('SNR'); ylabel('Z(R)'); title('Gabor'); set(gca,'XScale','log'); 
hold off;

subplot(122); plot(SNR(LumIds_conewts),atanh(meanR(LumIds_conewts,3)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(SNR(DOIds_conewts),atanh(meanR(DOIds_conewts,3)),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[1 16],'Ylim',[0 1.6],'XTick',[1 2 4 8 16],'YTick',0:0.8:1.6,'XScale','log'); 

X = [repmat(data1log,1,numcells).*celltype_dv celltype_dv];
[b, bint, r_DoG,~, stats1_DoG] = regress(data3,X);
[b_full, bint_full, r_full_DoG,~,stats0_DoG] = regress(data3,[data1log ones(size(data1log))]);
for i = 1:numcells
    tmp = [b(i) b(i+numcells)]*[log10(lims); ones(size(lims))];
    plot(lims,tmp,'color',colors{i}); hold on;
end
plot(lims,[b_full(1) b_full(2)]*[log10(lims); ones(size(lims))],'g', 'Linewidth', 2); hold on;
axis square; xlabel('SNR'); ylabel('Z(R)'); title('DoG'); hold off;
plot_counter = plot_counter + 1;

% Doing some stats: F test
dp = numcells*2-2;
RSSE_Gabor = sum(r_Gabor.^2);
RSSE_full_Gabor = sum(r_full_Gabor.^2);
Fstat_Gabor = ((RSSE_full_Gabor - RSSE_Gabor)/dp)/(RSSE_Gabor/(numel(data1)-numcells*2));
p_Gabor = 1-fcdf(Fstat_Gabor,dp,(numel(data1)-numcells*2));

RSSE_DoG = sum(r_DoG.^2);
RSSE_full_DoG = sum(r_full_DoG.^2);
Fstat_DoG = ((RSSE_full_DoG - RSSE_DoG)/dp)/(RSSE_DoG/(numel(data1)-numcells*2));
p_DoG = 1-fcdf(Fstat_DoG,dp,(numel(data1)-numcells*2));

%% Section 7: Greg's approach for multiclass linear regression
% making some fake data
n = 10; % n per group
Z = unifrnd(0,1.6,n*3,1);
celltype = [repmat(1,n,1); repmat(2,n,1);repmat(3,n,1)];
SNR = exprnd(4,n*3,1);
colors = {'red','blue','green'};
 
% Making the design matrix
tmp = [];
for i = 1:2
  L = celltype == i;
  tmp = [tmp, L];
end
X1 = [ones(length(Z),1), SNR, tmp, tmp.*repmat(SNR,1,2)]; % Separate slopes and intercepts for each group
X0 = [ones(length(Z),1), SNR]; % Same slope and intercept for each group
 
% Fitting the models
[b1,bint1,r1,~,stats1] = regress(Z,X1); % full model
[b0,bint0,r0,~,stats0] = regress(Z,X0); % restricted model
 
% Plotting the model fits
figure; axes; hold on;
for celltype_idx = 1:3
  L = celltype == celltype_idx;
  plot(SNR(L),Z(L),'o','color',colors{celltype_idx});
  intercept = b1(1);
  slope = b1(2);
  if celltype_idx == 1
    intercept = intercept+b1(3);
    slope = slope+b1(5);
  elseif celltype_idx == 2
    intercept = intercept+b1(4);
    slope =slope+b1(6);
  elseif celltype_idx == 3
    % Do nothing
  else
    error('unknown cell type');
  end
  plot([1 12],[1 1; 1 12]*[intercept; slope],'-','color',colors{celltype_idx})
end
plot([1 12],[1 1; 1 12]*b0,'k:')
 
% Doing the nested F test
SSE1 = sum(r1.^2);
SSE0 = sum(r0.^2);
dp = (size(X1,2)-size(X0,2)); % difference in number of parameters
F = (SSE0-SSE1)/dp/stats1(4);
p = fcdf(F,dp,n*3-dp);

%% Section 8: Checking whether DO and simple cell RF locations differed systematically

if ~exist('plot_counter')
    plot_counter = 1;
end
% load Output_ListWN.mat
load Output_ListWN2.mat
load Singleopponent.mat
load modelfits.mat
load Ratio_of_power.mat
crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;

Output_List(~Z_cellsofinterest,:) = [];
NLI = cell2mat(Output_List(:,13));
simplecells = NLI<0;
Singleopponent = logical(Singleopponent);
% Singleopponent = logical(Ratio_of_power<1.2);

SNR = Zmax(Z_cellsofinterest)/300; % Dividing by the number of elements
SNR = SNR(~Singleopponent & simplecells);
ind = find(~Singleopponent & simplecells);

% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

% calculating the cone weights
conewts_svd = cell2mat(Output_List(~Singleopponent & simplecells,23)');
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

Lumind = ind(LumIds_conewts); 
COind = ind(ColorOpponentIds_conewts); 
Sconeind = ind(Sconedominated_conewts);
DOind = [COind; Sconeind];

% Extracting the RF locations
RF_loc = cell2mat(Output_List(:,8))/10;
ecc = sqrt(sum(RF_loc.^2,2));

% Mann-Whitney U test
[p,h] = ranksum(ecc(Lumind),ecc(DOind));


% Doing a bootstrap (randomization) test: Similar to what Greg performed in 1999 paper 
niter = 5000;
RF_interest = RF_loc([Lumind; DOind]);
indices = [zeros(size(Lumind)); ones(size(DOind))];
randomizedteststats = zeros(niter,1);
for ii = 1:niter
    L = logical(indices(randperm(length(indices))));
    
    centroid_Lum = mean(RF_interest(L,:),1);
    centroid_DO = mean(RF_interest(~L,:),1);
    dist = sqrt(sum((centroid_Lum - centroid_DO).^2));
    
    randomizedteststats(ii) = dist;
end

actualdist = sqrt(sum((mean(RF_loc(Lumind,:),1) - mean(RF_loc(DOind,:),1)).^2));
p_value = sum(randomizedteststats>actualdist)/niter;

% Plotting the eccentricities 
figure(plot_counter);
subplot(121); hold on;
plot(RF_loc(DOind,1),RF_loc(DOind,2),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(RF_loc(Lumind,1),RF_loc(Lumind,2),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Ylim',[-8 2]); axis square; xlabel('X'); ylabel('Y'); title('RF locations'); 
legend('Simple','DO'); hold off;
subplot(122); histogram(randomizedteststats,20,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); hold on;
plot([actualdist actualdist],[0 700],'color','k','Linewidth',2)
set(gca,'Tickdir','out'); axis square; title('Randomization test'); 
ylabel('Count'); xlabel('Distance between simple and DO centroids'); 
legend('Randomization data', 'Actual data'); hold off;
plot_counter = plot_counter + 1;


