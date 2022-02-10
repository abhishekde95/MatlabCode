% Post-rebutall script
% Author -Abhishek De
% Date - 02/22

close all; clearvars;

%% 1. Relationship between Isoresponse NLI and target firing rates

if ~exist('plot_counter')
    plot_counter = 1;
end

% Loading all the files 
try 
    % Using the JDBC connection
    conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
    filename = fetch(conn,'SELECT filename FROM WNthresh');
    NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
    spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
    close(conn);
    filename = filename(strcmp(string(NTmode),"subunit"));
    NTmode = NTmode(strcmp(string(NTmode),"subunit"));
    spikeidx_NT = spikeidx_NT(strcmp(string(NTmode),"subunit"));

catch
    csv_filename = '/Users/abhishekde/Desktop/MatlabCode/Abhishek/CSV_PHPmyadmin_files/WNthresh.csv';
    [filename, NTmode, spikeIdx] = get_WNthreshdata_from_csvfile(csv_filename, 'subunit');
    spikeidx_NT = str2num(cell2mat(spikeIdx));
end


load conewts_svd.mat
load vals.mat
load vals.mat
load S1RGB_svd.mat
load S2RGB_svd.mat
load angulardifferences_RGB.mat
anglebwvectors = angulardifference_RGB;
S1RGB = S1RGB_svd;
S2RGb = S2RGB_svd;
% SpatiallyOpponent = sum(sign(S1RGB).*sign(S2RGB),1)<3;
SpatiallyOpponent = anglebwvectors'>90;

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

LUMidx = LumIds_conewts;
DOidx = [ColorOpponentIds_conewts Sconedominated_conewts];
hardtoclassifyidx = [Other_conewts];
hardtoclassifyidx = [hardtoclassifyidx LUMidx(vals(LUMidx)>=95) DOidx(vals(DOidx)>=95)];
LUMidx = LUMidx(vals(LUMidx)<95);
DOidx = DOidx(vals(DOidx)<95);

% Considering only the spatially opponent subunits
LUMidx = LUMidx(SpatiallyOpponent(LUMidx));
DOidx = DOidx(SpatiallyOpponent(DOidx));
hardtoclassifyidx = hardtoclassifyidx(SpatiallyOpponent(hardtoclassifyidx));

% Checking how target firing rates of the DO, simple and unclassified cells
load baselineFRstats.mat % Baseline FR rates
load TFR.mat % Target firing rates


% Load the isoresponse data
load RSSE_linearmodel_CV.mat % Robust regression
load RSSE_quadmodel_CV.mat


% For storing the Isoresponse NLI
Isoresponse_NLI = [];
BaselineFR = [];

for ii = 1:numel(RSSE_linearmodel) 
   
    % Isoresponse NLI
    Isoresponse_NLI = [Isoresponse_NLI; log10(median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii}))];
    BaselineFR = [BaselineFR; mean(baselineFRstats{ii})];
end

% Measuring correlation between TFR and NLI
idx = [LUMidx, DOidx, hardtoclassifyidx];
[r, p] = corr(Isoresponse_NLI(idx), TFR(1,idx)', 'type', 'Spearman');
[r2, p2] = corr(Isoresponse_NLI(idx), TFR(1,idx)'-BaselineFR(idx), 'type', 'Spearman');  % Baseline subtracted


%% 2. z-scores of Target firing rates w.r.t Baseline firing rates

close all; clearvars;

if ~exist('plot_counter')
    plot_counter = 1;
end

% Loading all the files 
try 
    % Using the JDBC connection
    conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
    filename = fetch(conn,'SELECT filename FROM WNthresh');
    NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
    spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
    close(conn);
    filename = filename(strcmp(string(NTmode),"subunit"));
    NTmode = NTmode(strcmp(string(NTmode),"subunit"));
    spikeidx_NT = spikeidx_NT(strcmp(string(NTmode),"subunit"));

catch
    csv_filename = '/Users/abhishekde/Desktop/MatlabCode/Abhishek/CSV_PHPmyadmin_files/WNthresh.csv';
    [filename, NTmode, spikeIdx] = get_WNthreshdata_from_csvfile(csv_filename, 'subunit');
    spikeidx_NT = str2num(cell2mat(spikeIdx));
end


load conewts_svd.mat
load vals.mat
load vals.mat
load S1RGB_svd.mat
load S2RGB_svd.mat
load angulardifferences_RGB.mat
anglebwvectors = angulardifference_RGB;
S1RGB = S1RGB_svd;
S2RGb = S2RGB_svd;
% SpatiallyOpponent = sum(sign(S1RGB).*sign(S2RGB),1)<3;
SpatiallyOpponent = anglebwvectors'>90;

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

LUMidx = LumIds_conewts;
DOidx = [ColorOpponentIds_conewts Sconedominated_conewts];
hardtoclassifyidx = [Other_conewts];
hardtoclassifyidx = [hardtoclassifyidx LUMidx(vals(LUMidx)>=95) DOidx(vals(DOidx)>=95)];
LUMidx = LUMidx(vals(LUMidx)<95);
DOidx = DOidx(vals(DOidx)<95);

% Considering only the spatially opponent subunits
LUMidx = LUMidx(SpatiallyOpponent(LUMidx));
DOidx = DOidx(SpatiallyOpponent(DOidx));
hardtoclassifyidx = hardtoclassifyidx(SpatiallyOpponent(hardtoclassifyidx));

% Checking how target firing rates of the DO, simple and unclassified cells
load baselineFRstats.mat % Baseline FR rates
load TFR.mat % Target firing rates


% Error plot representation
figure(plot_counter);
cell_idx = [DOidx LUMidx hardtoclassifyidx];
[~, indices] = sort(TFR(1,cell_idx));
indices = cell_idx(indices);
TFRprctile = [];
TFRzscore = []; % for storing the z-scores
TFRprctile_phase3 = [];
TFRzscore_phase3 = []; % for storing the z-scores
for iter = 1:numel(indices)
    disp(iter);
    
    ii = indices(iter);
    
    % Pulling out the respective file
    fileofinterest = char(filename(ii,:));
    stro = nex2stro(findfile(fileofinterest));
    spikename = 'sig001a';
    maskidx = strcmp(stro.sum.rasterCells(1,:),'subunit_mask');
    spikeidx = strcmp(stro.sum.rasterCells(1,:),spikename);
    basisvecidx = strcmp(stro.sum.rasterCells(1,:),'basis_vec');
    stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on');
    fpacqidx = strcmp(stro.sum.trialFields(1,:),'fpacq');
    
    mask_changes = [2];
    all_masks = stro.ras(:,maskidx);
    Fx = @(xi) any(isnan(xi)); % function that finds 'NaN' in a cell array
    inds = find(cellfun(Fx,stro.ras(:,basisvecidx))==0);
    if isempty(inds)
        inds = size(stro.trial,1)-1;
    end
    last_wntrial =  inds(1)-1;
    
    % Separately plotting the baseline firing rates for Phases 1 & 2, and Phase 3
    num_spikes =[];
    num_dur = [];
    for jj = 1:size(stro.ras,1)
        tmp = stro.ras{jj,spikeidx} ;
        if ~isnan((stro.trial(jj,stimonidx)- stro.trial(jj,fpacqidx)))
            spikes = tmp(tmp<stro.trial(jj,stimonidx) & tmp>stro.trial(jj,fpacqidx));
            num_spikes = [num_spikes; numel(spikes)];
            num_dur = [num_dur; (stro.trial(jj,stimonidx)- stro.trial(jj,fpacqidx))];
        end
    end
    %baselineFRstats{aa} = num_spikes./num_dur;
    
    baselineFR = num_spikes./num_dur;
    baselineFR_phase12 = baselineFR(1:last_wntrial);
    baselineFR_phase3 = baselineFR(last_wntrial+1:end);
    
    if ismember(ii,DOidx)
        c = [1 0 0];
    elseif ismember(ii, LUMidx)
        c = [0 0 0];
    else
        c = [0.5 0.5 0.5];
    end
    
    % Prctile for target firing rate as a function of baseline FRs from all
    % phases 
    TFRprctile = [TFRprctile; invprctile(baselineFR, TFR(1,ii))];
    TFRzscore = [TFRzscore; (TFR(1,ii)-mean(baselineFR))/std(baselineFR)];
    
    % Prctile for target firing rate as a function of baseline FRs from all
    % phases 
    TFRprctile_phase3 = [TFRprctile_phase3; invprctile(baselineFR_phase3, TFR(1,ii))];
    TFRzscore_phase3 = [TFRzscore_phase3; (TFR(1,ii)-mean(baselineFR_phase3))/std(baselineFR_phase3)];
    
end

%%



%% 3. Data points and model fits for example cells in R-theta space

close all; clearvars;


if ~exist('plot_counter')
    plot_counter = 1;
end

load conewts_svd.mat
load vals.mat
thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

% Classifying cells 
LUMidx = LumIds_conewts;
DOidx = [ColorOpponentIds_conewts Sconedominated_conewts];
hardtoclassifyidx = [Other_conewts];
hardtoclassifyidx = [hardtoclassifyidx LUMidx(vals(LUMidx)>=95) DOidx(vals(DOidx)>=95)];
LUMidx = LUMidx(vals(LUMidx)<95);
DOidx = DOidx(vals(DOidx)<95);

load RSSE_linearmodel_CV.mat % Robust regression
load RSSE_quadmodel_CV.mat

% Calculating the median of differences/ratios
RSSEisoresp_medianofratios = [];
for ii = 1:numel(RSSE_linearmodel)   
    % computation for calculating median of differences/ratios
    RSSEisoresp_medianofratios = [RSSEisoresp_medianofratios; median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii})];   
end

% Loading all the files 
try 
    % Using the JDBC connection
    conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
    filename = fetch(conn,'SELECT filename FROM WNthresh');
    NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
    spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
    close(conn);
    filename = filename(strcmp(string(NTmode),"subunit"));
    NTmode = NTmode(strcmp(string(NTmode),"subunit"));
    spikeidx_NT = spikeidx_NT(strcmp(string(NTmode),"subunit"));

catch
    csv_filename = '/Users/abhishekde/Desktop/MatlabCode/Abhishek/CSV_PHPmyadmin_files/WNthresh.csv';
    [filename, NTmode, spikeIdx] = get_WNthreshdata_from_csvfile(csv_filename, 'subunit');
    spikeidx_NT = str2num(cell2mat(spikeIdx));
end

% These are files that I am concerned with
load RHO_all.mat
load THETA_all.mat
load not_oog_idx_all.mat
load oog_idx_all.mat
load linear_modelparams.mat
load quad_modelparams.mat
load subunitbasisvec.mat

indices = [109 24 74];% example LUM, DO and HTC cells
plot_counter = 1;
for zz = 1:numel(indices)
    
    ind = indices(zz);
    THETA = THETA_all{1,ind};
    THETA = THETA * pi/180; % converting to radians
    if any(THETA>(135*pi/180))
        allthetas = linspace(-pi,pi,100);
        newtheta = linspace(-pi,pi,101);
    else
        allthetas = linspace(-pi/4,3*pi/4,100);
        newtheta = linspace(-pi/4,3*pi/4,101);
    end
    RHO = RHO_all{1,ind};
    oog_idx = oog_idx_all{1,ind};
    not_oog_idx = not_oog_idx_all{1,ind};
    outofgamut = zeros(size(THETA));
    outofgamut(oog_idx) = 1;
    outofgamut = logical(outofgamut);
    [x_orig, y_orig] = pol2cart(THETA,RHO);

    % Linear model predictions
    rho1 = 1./(linear_modelparams(ind,:)*[cos(allthetas); sin(allthetas)]);
    LOOGtmp1= rho1<0;
    [x_lin,y_lin] = pol2cart(allthetas(~LOOGtmp1),rho1(~LOOGtmp1));
    
    % Quadratic model predictions
    [x_quad,y_quad,rho3] = calc_xyvalues(allthetas, quad_modelparams(ind,:));
    L = rho3>0 & rho3==real(rho3);
    [x_quad2,y_quad2] = pol2cart(newtheta(L),rho3(L)');
    
    
    % Plotting the data, linear and non-linear fits in the R-Theta space
    figure(plot_counter); subplot(1,numel(indices),zz);
    plot(allthetas(~LOOGtmp1),log10(rho1(~LOOGtmp1)),'g','Linewidth',2); hold on; 
    plot(newtheta(L),log10(rho3(L)),'r','LineWidth',2); hold on; 
    plot(THETA(not_oog_idx), log10(RHO(not_oog_idx)),'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'PickableParts','none','MarkerEdgeColor',[1 1 1]);
    plot(THETA(oog_idx), log10(RHO(oog_idx)),'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'PickableParts','none','MarkerEdgeColor',[1 1 1]); 
    xlabel('Theta'); ylabel('log10(R)');
    if zz == 1
        title('Simple');
        %set(gca,'Xlim',[-pi pi],'XTick',[-pi -pi/2 0 pi/2 pi],'XTicklabel',{'-pi','-pi/2','0','pi/2','pi'},'Tickdir','out'); axis square;
        set(gca,'Xlim',[-pi/4 0.75*pi],'XTick',[-pi/4 0 pi/4 pi/2 0.75*pi],'XTicklabel',{'-pi/4', '0', 'pi/4','pi/2','0.75pi'},'Tickdir','out'); axis square;
        
    elseif zz==2
        title('DO');
        set(gca,'Xlim',[-pi/4 0.75*pi],'XTick',[-pi/4 0 pi/4 pi/2 0.75*pi],'XTicklabel',{'-pi/4', '0', 'pi/4','pi/2','0.75pi'},'Tickdir','out'); axis square;

    else
        title('Unclassified')
        set(gca,'Xlim',[-pi/4 0.75*pi],'XTick',[-pi/4 0 pi/4 pi/2 0.75*pi],'XTicklabel',{'-pi/4', '0', 'pi/4','pi/2','0.75pi'},'Tickdir','out'); axis square;

    end
    axis square; hold off;
   
end


%% 4. RF as a function of cell types


if ~exist('plot_counter')
    plot_counter = 1;
end

load conewts_svd.mat
load vals.mat
load S1RGB_svd.mat
load S2RGB_svd.mat
load angulardifferences_RGB.mat
S1RGB = S1RGB_svd;
S2RGB = S2RGB_svd;
anglebwvectors = angulardifference_RGB;
% SpatiallyOpponent = sum(sign(S1RGB).*sign(S2RGB),1)<3;
SpatiallyOpponent = anglebwvectors'>90;


thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

LUMidx = LumIds_conewts;
DOidx = [ColorOpponentIds_conewts Sconedominated_conewts];
hardtoclassifyidx = [Other_conewts];
hardtoclassifyidx = [hardtoclassifyidx LUMidx(vals(LUMidx)>=95) DOidx(vals(DOidx)>=95)];
LUMidx = LUMidx(vals(LUMidx)<95);
DOidx = DOidx(vals(DOidx)<95);

% Considering only the spatially opponent subunits
LUMidx = LUMidx(SpatiallyOpponent(LUMidx));
DOidx = DOidx(SpatiallyOpponent(DOidx));
hardtoclassifyidx = hardtoclassifyidx(SpatiallyOpponent(hardtoclassifyidx));


% Loading all the files 
try 
    % Using the JDBC connection
    conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
    filename = fetch(conn,'SELECT filename FROM WNthresh');
    NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
    spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
    close(conn);
    filename = filename(strcmp(string(NTmode),"subunit"));
    NTmode = NTmode(strcmp(string(NTmode),"subunit"));
    spikeidx_NT = spikeidx_NT(strcmp(string(NTmode),"subunit"));

catch
    csv_filename = '/Users/abhishekde/Desktop/MatlabCode/Abhishek/CSV_PHPmyadmin_files/WNthresh.csv';
    [filename, NTmode, spikeIdx] = get_WNthreshdata_from_csvfile(csv_filename, 'subunit');
    spikeidx_NT = str2num(cell2mat(spikeIdx));
end



% Pulling out the RFs
RF = [];
for ii=1:numel(filename)
    WN = nex2stro(findfile(char(filename(ii))));
    x = WN.sum.exptParams.rf_x;
    y = WN.sum.exptParams.rf_y;
    
    RF = [RF; sqrt(x^2 + y^2)/10];
 
end

% Comparing RFs across cell types
group = [ones(size(LUMidx)) 2*ones(size(DOidx)) 3*ones(size(hardtoclassifyidx))];
data = RF([LUMidx'; DOidx'; hardtoclassifyidx']); 
p1 = kruskalwallis(data,group','off');