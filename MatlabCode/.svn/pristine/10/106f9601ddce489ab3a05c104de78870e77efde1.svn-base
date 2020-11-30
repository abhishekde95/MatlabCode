% Analyzing direction selectivity of V1 data 
% Author - Abhishek De, 1/19
close all; clearvars;
plot_counter = 1;
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
maxT = 12;
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

% Merging all the files in the list
Input_List = [filename_Lum; filename_ColorOpponent; filename_NT; filename_subunit];
spikeIdx = [spikeIdx_Lum; spikeIdx_ColorOpponent; spikeidx_NT; spikeidx_subunit];

numcells = numel(Input_List);
files_not_working = [];
files_not_working_idxs = [];
count = 1;
CHI2CRIT = .95;
nrows = 10;
STOI = [];
STOI_euclidean = [];
for ii = 1:numcells
    flag = 0;
    disp(ii);
    filename = char(Input_List{ii}{1}); % acquiring the filename (1st column) from the List
    
    WN = {};
    for jj = 1:size(Input_List{ii},2)
        try 
            tmpstro = nex2stro(findfile(char(Input_List{ii}(jj))));
        catch 
            files_not_working = [files_not_working; Input_List{ii}];
            files_not_working_idxs = [files_not_working_idxs; ii];
            flag = 1;
            break;
        end
        if (isempty(WN))
            WN = tmpstro;
        else
            WN = strocat(WN, tmpstro);
        end
        if ~any(strcmp(WN.sum.rasterCells,'sig001a_wf'))
            files_not_working = [files_not_working; Input_List{ii}];
            files_not_working_idxs = [files_not_working_idxs; ii];
            flag = 1;
        end
    end
    if flag 
        continue;
    end
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
    maskidx = strcmp(WN.sum.rasterCells(1,:), 'subunit_mask');
    anlgStartTimeidx = find(strcmp(WN.sum.rasterCells(1,:),'anlgStartTime'));
    L = WN.trial(:,noisetypeidx)==1;
    mu1idx = find(strcmp(WN.sum.trialFields(1,:),'mu1'));
    mu2idx = find(strcmp(WN.sum.trialFields(1,:),'mu2'));
    mu3idx = find(strcmp(WN.sum.trialFields(1,:),'mu3'));
    sigma1idx = find(strcmp(WN.sum.trialFields(1,:),'sigma1'));
    sigma2idx = find(strcmp(WN.sum.trialFields(1,:),'sigma2'));
    sigma3idx = find(strcmp(WN.sum.trialFields(1,:),'sigma3'));
    maskidx = strcmp(WN.sum.rasterCells(1,:),'subunit_mask');
    basisvecidx = strcmp(WN.sum.rasterCells(1,:),'basis_vec');
    gammaTable = WN.sum.exptParams.gamma_table;
    gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
    gammaTable1 = interp1(linspace(0,255,256),gammaTable,linspace(0,255,65536), 'spline');
    invgamma = InvertGamma(gammaTable, 0);
    sigmavect = unique(WN.trial(L,[sigma1idx sigma2idx sigma3idx]),'rows')/1000;
    sigmavect(all(any(sigmavect == 0),2),:) = [];
    gausslims = [WN.sum.exptParams.gauss_locut WN.sum.exptParams.gauss_hicut]/1000;
    x = linspace(gausslims(1),gausslims(2),NPOINTS);
    Fx = norminv(x)*sigmavect(1);
    sigmacorrectionfactor = std(Fx)./sigmavect(1);
    muvar = (sigmavect(1)*sigmacorrectionfactor)^2;
    
    % Getting the background rgb/lms
    ridx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_r'));
    gidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_g'));
    bidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_b'));
    bkgndRGB = [mode(WN.trial(:,ridx)), mode(WN.trial(:,gidx)), mode(WN.trial(:,bidx))];
    bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
    bkgndlms = M*bkgndrgb;

    WN.ras(~L ,:) = []; % modiftying the WN structure
    WN.trial(~L,:) = []; % modiftying the WN structure
    mask_changes = [2 size(WN.trial,1)];
    if any(basisvecidx)
        mask_changes = [2];
        all_masks = WN.ras(:,maskidx);
        Fx = @(xi) any(isnan(xi)); % function that finds 'NaN' in a cell array
        inds = find(cellfun(Fx,WN.ras(:,basisvecidx))==0);
        if isempty(inds)
            inds = size(WN.trial,1)-1;
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
        
        idxs = zeros(size(WN.trial,1),1);
        idxs(mask_changes(2,1)+1:end) = 1;
        idxs = logical(idxs);
        WN.ras(idxs,:) = []; % modiftying the WN structure 
        WN.trial(idxs,:) = []; % modiftying the WN structure
    end
        
    spikeidx = spikeIdx(ii);
    spikename = spikename_options(spikeidx,:); 
    
    % Calculating STA and STC for frames which triggered spikes
    out_gun = getWhtnsStats(WN,maxT,'STCOVmex',{nstixperside^2,3,maxT},spikename);
    STS_gun = out_gun{1}; STCross_gun = out_gun{2}; nspikes_gun = out_gun{3}; clear out_gun;
    STAs_gun = STS_gun/nspikes_gun; 
    
    T = reshape(STAs_gun,[nstixperside^2 channels maxT]);
    T = RGB2XWFormat(permute(T,[1 3 2]));
    [u,s,v] = svd(T);
    
    % calculating cone weights
    Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
    Mrgbtocc = inv(Mrgbtocc');
    wts = Mrgbtocc * v(:,1);
    wts = wts./repmat(sum(abs(wts),1),[3 1]);
    wts = wts .* repmat(sign(wts(2,:)),[3 1]);
    
    % calculating space time map
    spacetimemap = reshape(u(:,1),[nstixperside nstixperside maxT]);
    maxt = sum(spacetimemap.^2,3);  
    [row,col] = find(maxt==max(maxt(:)));
    
    xtmap = permute(squeeze(spacetimemap(:,col,:)),[2 1]);
    Fxtmap = fftshift(fft2(xtmap));
    ytmap = permute(squeeze(spacetimemap(row,:,:)),[2 1]);
    Fytmap = fftshift(fft2(ytmap));
    figure(plot_counter); subplot(nrows,4,4*count-3); imagesc(xtmap); colormap('gray'); axis xy; axis square; set(gca,'Tickdir','out','XTick',[0 nstixperside],'YTick',[0 maxT]);
    subplot(nrows,4,4*count-2); imagesc(abs(Fxtmap)); colormap('gray'); axis xy; axis square; set(gca,'Tickdir','out','XTick',[0 nstixperside],'YTick',[0 maxT]);
    subplot(nrows,4,4*count-1); imagesc(ytmap); colormap('gray'); axis xy; axis square; set(gca,'Tickdir','out','XTick',[0 nstixperside],'YTick',[0 maxT]);
    subplot(nrows,4,4*count); imagesc(abs(Fytmap)); colormap('gray'); axis xy; axis square; set(gca,'Tickdir','out','XTick',[0 nstixperside],'YTick',[0 maxT]);
    
    xmid = maxT/2; ymid = nstixperside/2;
    Fxpower = abs(Fxtmap).^2; Fypower = abs(Fytmap).^2;
    
    p1x = Fxpower(1:xmid,ymid+1:nstixperside);
    p2x = Fxpower(1:xmid,1:ymid);
    p3x = Fxpower(xmid+1:maxT,1:ymid);
    p4x = Fxpower(xmid+1:maxT,ymid+1:nstixperside);
    STOI_x = abs(sum(p1x(:))+ sum(p3x(:)) - sum(p2x(:))-sum(p4x(:)))/abs(sum(Fxpower(:))); 
    
    p1y = Fypower(1:xmid,ymid+1:nstixperside);
    p2y = Fypower(1:xmid,1:ymid);
    p3y = Fypower(xmid+1:maxT,1:ymid);
    p4y = Fypower(xmid+1:maxT,ymid+1:nstixperside);
    STOI_y = abs(sum(p1y(:))+ sum(p3y(:)) - sum(p2y(:))-sum(p4y(:)))/abs(sum(Fypower(:)));
     
    STOI = [STOI; max([STOI_x STOI_y])];
    STOI_euclidean = [STOI_euclidean; sqrt(STOI_x^2 + STOI_y^2)];
    
    count = count + 1;
    if count == nrows+1
        count = 1;
        plot_counter = plot_counter + 1;
    end
end

savevariables = 0;
if savevariables
    save STOI STOI
    save STOI_euclidean STOI_euclidean
end
%% Further direction selectivty analyses  

load STOI.mat
load STOI_euclidean.mat
load Output_ListWN2.mat
load Singleopponent.mat
crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
STOI(~Z_cellsofinterest) = [];
STOI_euclidean(~Z_cellsofinterest) = [];
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

STOI_euclidean = STOI_euclidean(~Singleopponent & simplecells); 
STOI = STOI(~Singleopponent & simplecells);
RGB_svd = cell2mat(Output_List(~Singleopponent & simplecells,5)');
bkgndlms = cell2mat(Output_List(1,16)');
Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
Mrgbtocc = inv(Mrgbtocc');
conewts_svd = Mrgbtocc * RGB_svd;
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];

% plotting results from STOI analyses 
figure(plot_counter); set(gcf,'Name','Cone wts') 
plot(conewts_svd(1,LumIds_conewts),conewts_svd(2,LumIds_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(1,ColorOpponentIds_conewts),conewts_svd(2,ColorOpponentIds_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(1,Sconedominated_conewts),conewts_svd(2,Sconedominated_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]);
plot(conewts_svd(1,Other_conewts),conewts_svd(2,Other_conewts),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out'); plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');  xlabel('L'), ylabel('M');
plot_counter = plot_counter + 1;

figure(plot_counter); set(gcf,'Name','STOI') 
subplot(221); histogram(STOI(LumIds_conewts),0:0.05:1,'FaceColor',[0 1 0]); hold on; plot(median(STOI(LumIds_conewts)),0,'kv'); axis square; set(gca,'Tickdir','out','XTick',0:0.25:1.0,'YTick',0:5:15,'Ylim',[0 15]); ylabel('# cells'); xlabel('STOI'); title('LM non-opponent'); hold off;
subplot(222); histogram(STOI(ColorOpponentIds_conewts),0:0.05:1,'FaceColor',[1 0 0]);  hold on; plot(median(STOI(ColorOpponentIds_conewts)),0,'kv'); axis square; set(gca,'Tickdir','out','XTick',0:0.25:1.0,'YTick',0:10:30,'Ylim',[0 30]); ylabel('# cells'); xlabel('STOI'); title('LM opponent'); hold off;
subplot(223); histogram(STOI(Sconedominated_conewts),0:0.05:1,'FaceColor',[0 0 1]);  hold on; plot(median(STOI(Sconedominated_conewts)),0,'kv'); axis square; set(gca,'Tickdir','out','XTick',0:0.25:1.0,'YTick',0:5:15,'Ylim',[0 15]); ylabel('# cells'); xlabel('STOI'); title('S-cone sensitive'); hold off;
subplot(224); histogram(STOI(Other_conewts),0:0.05:1,'FaceColor',[0 0 0]);  hold on; plot(median(STOI(Other_conewts)),0,'kv'); axis square; set(gca,'Tickdir','out','XTick',0:0.25:1.0,'YTick',0:5:15,'Ylim',[0 15]); ylabel('# cells'); xlabel('STOI'); title('Others'); hold off;
plot_counter = plot_counter + 1;

figure(plot_counter); set(gcf,'Name','STOI euclidean') 
subplot(221); histogram(STOI_euclidean(LumIds_conewts),0:0.05:1,'FaceColor',[0 1 0]); hold on; plot(median(STOI_euclidean(LumIds_conewts)),0,'kv'); axis square; set(gca,'Tickdir','out','XTick',0:0.25:1.0,'YTick',0:5:15,'Ylim',[0 15]); ylabel('# cells'); xlabel('STOI'); title('LM non-opponent'); hold off;
subplot(222); histogram(STOI_euclidean(ColorOpponentIds_conewts),0:0.05:1,'FaceColor',[1 0 0]);  hold on; plot(median(STOI_euclidean(ColorOpponentIds_conewts)),0,'kv'); axis square; set(gca,'Tickdir','out','XTick',0:0.25:1.0,'YTick',0:10:30,'Ylim',[0 30]); ylabel('# cells'); xlabel('STOI'); title('LM opponent'); hold off;
subplot(223); histogram(STOI_euclidean(Sconedominated_conewts),0:0.05:1,'FaceColor',[0 0 1]);  hold on; plot(median(STOI_euclidean(Sconedominated_conewts)),0,'kv'); axis square; set(gca,'Tickdir','out','XTick',0:0.25:1.0,'YTick',0:5:15,'Ylim',[0 15]); ylabel('# cells'); xlabel('STOI'); title('S-cone sensitive'); hold off;
subplot(224); histogram(STOI_euclidean(Other_conewts),0:0.05:1,'FaceColor',[0 0 0]);  hold on; plot(median(STOI_euclidean(Other_conewts)),0,'kv'); axis square; set(gca,'Tickdir','out','XTick',0:0.25:1.0,'YTick',0:5:15,'Ylim',[0 15]); ylabel('# cells'); xlabel('STOI'); title('Others'); hold off;
plot_counter = plot_counter + 1;

% Trying out the Violin plot
figure(plot_counter); set(gcf,'Name','Violin plot: STOI'); 
subplot(121); violinPlot(STOI(LumIds_conewts), 'histOri','right','xyOri','flipped','widthDiv', [9 9], 'showMM', 3,'color',  [0 1 0]); hold on;
violinPlot(STOI(ColorOpponentIds_conewts), 'histOri', 'right','xyOri','flipped', 'widthDiv', [9 6], 'showMM', 3,'color',  [1 0 0]);
violinPlot(STOI(Sconedominated_conewts), 'histOri', 'right','xyOri','flipped','widthDiv', [9 3], 'showMM', 3,'color',  [0 0 1]);
violinPlot(STOI(Other_conewts), 'histOri', 'right','xyOri','flipped', 'widthDiv', [9 0], 'showMM', 3,'color',  [0 0 0]); 
set(gca,'Tickdir','out','XTick',0:0.25:1.0,'Xlim',[0 1.0]); xlabel('STOI'); ylabel('Data'); axis square; hold off;
subplot(122); violinPlot(STOI_euclidean(LumIds_conewts), 'histOri','right','xyOri','flipped','widthDiv', [9 9], 'showMM', 3,'color',  [0 1 0]); hold on;
violinPlot(STOI_euclidean(ColorOpponentIds_conewts), 'histOri', 'right','xyOri','flipped', 'widthDiv', [9 6], 'showMM', 3,'color',  [1 0 0]);
violinPlot(STOI_euclidean(Sconedominated_conewts), 'histOri', 'right','xyOri','flipped','widthDiv', [9 3], 'showMM', 3,'color',  [0 0 1]);
violinPlot(STOI_euclidean(Other_conewts), 'histOri', 'right','xyOri','flipped', 'widthDiv', [9 0], 'showMM', 3,'color',  [0 0 0]); 
set(gca,'Tickdir','out','XTick',0:0.25:1.0,'Xlim',[0 1.0]); xlabel('STOI euclidean'); ylabel('Data'); axis square; hold off;
plot_counter = plot_counter + 1;


% Need to do a kruskal-wallis test 

data1 = [STOI(LumIds_conewts); STOI(ColorOpponentIds_conewts); STOI(Sconedominated_conewts); STOI(Other_conewts)]; % STOI 
data2 = [STOI_euclidean(LumIds_conewts); STOI_euclidean(ColorOpponentIds_conewts); STOI_euclidean(Sconedominated_conewts); STOI_euclidean(Other_conewts)]; % STOI euclidean
group =[ones(size(STOI(LumIds_conewts))); 2*ones(size(STOI(ColorOpponentIds_conewts))); 3*ones(size(STOI(Sconedominated_conewts))); 4*ones(size(STOI(Other_conewts)))];
p1 = kruskalwallis(data1,group,'off');
p2 = kruskalwallis(data2,group,'off');






