% A new script for reviewers eyes only for Elife rebuttal
% Author- Abhishek De, 6/27/21

close all;
clearvars; 

%% 1. Reviewer figure: Consistency of Isoresponse contours: 
% Robustness of the isoresponse curves based on the last point
% Plotting the average of the last 5 points as well for each probed direction: checking for the example cells only

if ~exist('plot_counter')
    plot_counter = 1;
end

load conewts_svd.mat
load vals.mat
load S1RGB_svd.mat
load S2RGB_svd.mat
load angulardifferences_RGB.mat
anglebwvectors = angulardifference_RGB;
S1RGB = S1RGB_svd;
S2RGb = S2RGB_svd;
SpatiallyOpponent = anglebwvectors'>90;

% Classifying cells based on cone-weights and PC1 signficance
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

load RSSE_linearmodel_CV.mat % Robust regression
load RSSE_quadmodel_CV.mat

% Calculating the median of differences/ratios
RSSEisoresp_medianofratios = [];
for ii = 1:numel(RSSE_linearmodel)   
    % computation for calculating median of differences/ratios
    RSSEisoresp_medianofratios = [RSSEisoresp_medianofratios; median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii})];   
end

% Loading all the files 

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
plot_counter = 3;
for zz = 1:numel(indices)
    stro = nex2stro(findfile(char(filename(indices(zz)))));
    global spikename maskidx spikeidx neurothreshidx nstixperside ngammasteps seedidx nframesidx correctidx
    global fponidx stimoffidx stimonidx muidxs sigmaidxs basisvecidx weightsidx fpacqidx targetspikerateidx basisvecdiridx latencyidx reversalflagidx parentverticesidx
    global msperframe ntrials maxT xx yy M linepredtol stepsizescale stepsize nreversals oogscale
    spikename = 'sig001a';
    maskidx = strcmp(stro.sum.rasterCells(1,:),'subunit_mask');
    spikeidx = strcmp(stro.sum.rasterCells(1,:),spikename);
    basisvecidx = strcmp(stro.sum.rasterCells(1,:),'basis_vec');
    weightsidx = strcmp(stro.sum.rasterCells(1,:),'weights');
    parentverticesidx = strcmp(stro.sum.rasterCells(1,:),'parentvertices');
    nstixperside = stro.sum.exptParams.nstixperside;
    ngammasteps = 2^16; % 65536
    linepredtol = stro.sum.exptParams.linepredtol;
    stepsizescale = stro.sum.exptParams.stepsizescale;
    stepsize = stro.sum.exptParams.stepsize;
    nreversals = stro.sum.exptParams.nreversals;
    oogscale = stro.sum.exptParams.oogscale;
    seedidx = strcmp(stro.sum.trialFields(1,:),'seed');
    nframesidx = strcmp(stro.sum.trialFields(1,:),'num_frames');
    stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on');
    stimoffidx = strcmp(stro.sum.trialFields(1,:),'stim_off');
    fponidx = strcmp(stro.sum.trialFields(1,:),'fp_on');
    fpacqidx = strcmp(stro.sum.trialFields(1,:),'fpacq');
    basisvecdiridx = strcmp(stro.sum.trialFields(1,:),'weights_idx');
    neurothreshidx = strcmp(stro.sum.trialFields(1,:),'neurothresh'); % when exactly the neurothresh trials started
    targetspikerateidx = strcmp(stro.sum.trialFields(1,:),'targetspikerate');
    correctidx = strcmp(stro.sum.trialFields(1,:),'correct');
    muidxs = [find(strcmp(stro.sum.trialFields(1,:),'mu1')), ...
        find(strcmp(stro.sum.trialFields(1,:),'mu2')), ...
        find(strcmp(stro.sum.trialFields(1,:),'mu3'))];
    sigmaidxs = [find(strcmp(stro.sum.trialFields(1,:),'sigma1')), ...
        find(strcmp(stro.sum.trialFields(1,:),'sigma2')), ...
        find(strcmp(stro.sum.trialFields(1,:),'sigma3'))];
    latencyidx = strcmp(stro.sum.trialFields(1,:),'latency');
    reversalflagidx = strcmp(stro.sum.trialFields(1,:),'reversalflag');
    msperframe = 1000/stro.sum.exptParams.framerate;
    ntrials = size(stro.trial,1);
    maxT = 10; % this represents the temporal part in the spatiotemporal receptive field
    xx = linspace(stro.sum.exptParams.gauss_locut/1000, stro.sum.exptParams.gauss_hicut/1000,ngammasteps); % xx represents the probabilities. For more info, have a look at the MATLAB 'norminv' function.
    yy = norminv(xx'); % defining norminv to extract the values for which the cdf values range between gauss_locut and gauss_hicut
    t_offset = stro.trial(end,latencyidx)/1000;
    
    % Obtaining the M matrix, code extracted from Greg, fitting a cubic spline
    % using the command 'spline'. 'SplineRaw' only availabe through
    % psychtoolbox which I currently don't have now.
    fundamentals = stro.sum.exptParams.fundamentals; % CONE FUNDAMENTALS: L,M,S
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
    mon_spd = stro.sum.exptParams.mon_spd; % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
    M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals
    M = inv(M');
    
    
    Fx = @(xi) any(isnan(xi)); % function that finds 'NaN' in a cell array
    inds = find(cellfun(Fx,stro.ras(:,basisvecidx))==0);
    if isempty(inds)
        inds = size(stro.trial,1)-1;
    end
    
    % Determining when Neurothresh mode was active, plotting the basis vector, working correctly
    neurothreshmode = stro.trial(:,neurothreshidx);
    basisvec_dropidx = inds(end);
    neurothresh_startidx = max(find(neurothreshmode==1,1),basisvec_dropidx+1); % There is possibility of pressing the Neurothresh button multiple times
    num_targetspikerates = unique(stro.trial(neurothresh_startidx:end,targetspikerateidx));
    vect = stro.ras{basisvec_dropidx,basisvecidx};
    basisvec_size = nstixperside*nstixperside*3;
    numvect = (numel(vect)/basisvec_size)-1;
    basisvec = cell(1,numvect);
    % Actual basis vec
    for ii = 1:numvect
        tmp_vec = vect((ii-1)*basisvec_size+1:basisvec_size*ii) + vect(numvect*basisvec_size+1:basisvec_size*(numvect+1));
        basisvec{ii} = reshape(tmp_vec,[nstixperside nstixperside 3]);
    end
    bkgnd_monitor = reshape(vect(numvect*basisvec_size+1:basisvec_size*(numvect+1)),[nstixperside nstixperside 3]);
    tmp1 = unique(basisvec{1}-bkgnd_monitor,'stable');
    tmp2 = unique(basisvec{2}-bkgnd_monitor,'stable');

    % plotting the rasterplots for neurothresh trials
    norms = cell(1,numel(num_targetspikerates));
    completed_search_alongdir = cell(1,numel(num_targetspikerates));
    for jj = 1: numel(num_targetspikerates)
        idxs = find(~isnan(stro.trial(:,correctidx)) & stro.trial(:,targetspikerateidx)==num_targetspikerates(jj));
        idxs(idxs<=neurothresh_startidx) = [];
        different_weights = unique(stro.trial(idxs,basisvecdiridx));
        tmp_norm = [];
        tmp_completed_search_alongdir = [];
        
        for kk = 1:numel(different_weights)
            idxs1 = find(stro.trial(:,basisvecdiridx) == different_weights(kk));
            idxs1(idxs1<neurothresh_startidx) = [];
            raster_data = stro.ras(idxs1,1);
            tmp_norm = [tmp_norm; stro.ras{idxs1(end),weightsidx}'];
            for ii = 1:size(raster_data,1)
                tmp = raster_data{ii} ;
                spikes = tmp(tmp>stro.trial(idxs1(ii),stimonidx)+t_offset & tmp < stro.trial(idxs1(ii),stimoffidx));
                spikes = spikes - stro.trial(idxs1(ii),stimonidx)-t_offset;
            end
            [flag, gamutViolation] = Check_ifcompletely_probed(stro,basisvec,bkgnd_monitor,idxs1,tmp_norm(end,:));
            % flag = 0, incompletely probed
            % flag = 1, completely probed
            % gamutViolation = 1, out of gamut point
            tmp_completed_search_alongdir = [tmp_completed_search_alongdir; flag gamutViolation];
        end
        norms{jj} = tmp_norm;
        completed_search_alongdir{jj} = tmp_completed_search_alongdir;
    end
    
    
    % Plotting the staircase termination and the geometric mean of last 3 staircase
    % termination points
    completed_dir = completed_search_alongdir{1};
    probed_dirs = logical(completed_dir(:,1)==1); % only including the directions that have been completely probed
    oog_idx = completed_dir(:,1)==1 & completed_dir(:,2)==1; % probed and out of gamut
    not_oog_idx = completed_dir(:,1)==1 & completed_dir(:,2)==0;
    figure(1);
    for ii = 1:size(norms{1},1)
        dir = ii;
        for kk = 1: numel(num_targetspikerates)
            tmp_n = [];
            tmp_wts = [];
            tmp_parentvertices = [];
            idxs1 = find(stro.trial(:,basisvecdiridx) == dir);
            if ~isempty(idxs1)
                for jj = 1:numel(idxs1)
                    tmp_parentvertices = [tmp_parentvertices; stro.ras{idxs1(jj),parentverticesidx}];
                    tmp_wts = [tmp_wts; stro.ras{idxs1(jj),weightsidx}'];
                    tmp_n = [tmp_n; norm(stro.ras{idxs1(jj),weightsidx})];
                end
                raster_data = stro.ras(idxs1,1);
                num_dur =[];
                firing_rate = [];
                if (~isempty(tmp_n) & not_oog_idx(ii)) 
                    figure(1), 
                    subplot(3,2,2*zz-1); plot(sign(tmp_wts(end,1))*geomean(abs(tmp_wts(end-4:end,1))),sign(tmp_wts(end,2))*geomean(abs(tmp_wts(end-4:end,2))),'o','MarkerFacecolor',[0 0 0], 'MarkerEdgeColor',[1 1 1]); hold on;
                    subplot(3,2,2*zz); plot(tmp_wts(end,1),tmp_wts(end,2),'o','MarkerFacecolor',[0 0 0], 'MarkerEdgeColor',[1 1 1]); hold on;     
                elseif (~isempty(tmp_n) & oog_idx(ii))
                    figure(1), 
                    subplot(3,2,2*zz-1); plot([0 tmp_wts(end,1)],[0 tmp_wts(end,2)],'k'); hold on;
                    subplot(3,2,2*zz); plot([0 tmp_wts(end,1)],[0 tmp_wts(end,2)],'k'); hold on;     
                end
            end
        end
    end
 
   
    xlabel('X'), ylabel('Y');
    if zz == 1
        title('Simple');
    elseif zz==2
        title('DO');
    else
        title('Unclassified')
    end
    
     
    if zz==1 
        subplot(3,2,2*zz-1); set(gca,'Xlim',[-1.6 1.6],'Ylim',[-1.6 1.6],'YTick',-1.6:0.8:1.6,'XTick',[-1.6:0.8:1.6]); set(gca,'Tickdir','out'); drawnow; axis square; hold off;
        subplot(3,2,2*zz); set(gca,'Xlim',[-1.6 1.6],'Ylim',[-1.6 1.6],'YTick',-1.6:0.8:1.6,'XTick',[-1.6:0.8:1.6]);set(gca,'Tickdir','out'); drawnow; axis square; hold off;
    elseif zz==2
        subplot(3,2,2*zz-1); set(gca,'Xlim',[-2.0 2.0],'Ylim',[-2.0 2.0],'YTick',-2.0:1.0:2.0,'XTick',[-2.0:1.0:2.0]); set(gca,'Tickdir','out'); drawnow; axis square; hold off;
        subplot(3,2,2*zz); set(gca,'Xlim',[-2.0 2.0],'Ylim',[-2.0 2.0],'YTick',-2.0:1.0:2.0,'XTick',[-2.0:1.0:2.0]); set(gca,'Tickdir','out'); drawnow; axis square; hold off;
    elseif zz==3
        subplot(3,2,2*zz-1); set(gca,'Xlim',[-0.2 0.2],'Ylim',[-0.2 0.2],'YTick',-0.2:0.1:0.2,'XTick',[-0.2:0.1:0.2]); set(gca,'Tickdir','out'); drawnow; axis square; hold off;
        subplot(3,2,2*zz); set(gca,'Xlim',[-0.2 0.2],'Ylim',[-0.2 0.2],'YTick',-0.2:0.1:0.2,'XTick',[-0.2:0.1:0.2]); set(gca,'Tickdir','out'); drawnow; axis square; hold off;
    end
    
    % Now, I am going to plot the
    wts  = [1 0;1 1;0 1;-1 1;-1 0;-1 -1;0 -1;1 -1];
    subplotidxs = [1;2;3;4;6;7;8;9];
    for ii = 1:numel(subplotidxs)
        vec = wts(ii,1)*(basisvec{1}-bkgnd_monitor) + wts(ii,2)*(basisvec{2}-bkgnd_monitor);
        normfactor = 0.5/(max(abs(vec(:)))+0.01);
        vec = normfactor*vec + 0.5;
        figure(2),subplot(numel(indices),size(wts,1),size(wts,1)*(zz-1)+ii); image(vec); set(gca,'XTick',[],'YTick',[]); axis square;
    end
    
    if zz == 1
        title('Simple');
    elseif zz==2
        title('DO');
    else
        title('Unclassified')
    end
    axis square; hold off;
    
    plot_counter = plot_counter + 1;
end


%% 2. Reviewer analysis: Rebound OFF responses-part 1
% Expected OFF responses from monophasic and biphasic temporal IR functions
if ~exist('plot_counter')
    plot_counter = 1;
end

t = 0:1:80;
mp = poisspdf(t,30);
bp1 = poisspdf(t,30)-poisspdf(t,50);

% Building an artifical stimulus
stim = zeros(500,1);
stim(150:350) = 1;

% response of temporal IRF
mp_resp = conv(stim, mp,'full');
bp_resp1 = conv(stim,bp1,'full');

% Plotting the figures
figure(plot_counter);
subplot(221), plot(mp, 'Color',[1 0 0], 'linewidth', 1); hold on;
plot(-mp, 'Color',[0 1 0], 'linewidth', 1);
plot(t,zeros(size(t)), 'Color',[0 0 0]);
legend('L+ IRF', 'L- IRF'); xlabel('Time (ms)');  axis square; 
title('Monophasic mechanism Temporal IRF')
set(gca, 'Tickdir', 'out', 'Ylim', [-0.1, 0.1], 'Xlim', [0, max(t)]); 

subplot(222), plot(stim, 'Color',[0 0 0], 'linewidth', 2); hold on;
plot(2*mp_resp, '-.', 'Color',[0 0 0], 'linewidth', 2);
plot(-stim, 'Color',[0.5 0.5 0.5], 'linewidth', 2);
plot(-2*mp_resp,'-.', 'Color',[0.5 0.5 0.5], 'linewidth', 2);
legend('L+ Stimulus', 'Response to L+', 'L- Stimulus', 'Response to L-');
title('Monophasic responses')
xlabel('Time (ms)'); set(gca, 'Tickdir', 'out', 'Ylim', [-3, 3]); axis square;

subplot(223), hold on;
plot(bp1, 'Color',[1 0 0],'linewidth', 1);
plot(-bp1, 'Color',[0 1 0], 'linewidth', 1);
plot(t,zeros(size(t)), 'Color',[0 0 0]);
legend('L+ IRF', ' L- IRF'); xlabel('Time (ms)');  axis square;
title('Biphasic mechanism Temporal IRF')
set(gca, 'Tickdir', 'out', 'Ylim', [-0.1, 0.1], 'Xlim', [0, max(t)]);

subplot(224), plot(stim, 'Color',[0 0 0], 'linewidth', 1); hold on;
plot(2*bp_resp1,'-.', 'Color',[0 0 0], 'linewidth', 2);
plot(-stim, 'Color',[0.5 0.5 0.5], 'linewidth', 1);
plot(-2*bp_resp1,  '-.', 'Color',[0.5 0.5 0.5],'linewidth', 2);
legend('L+ Stimulus', 'Response to L+', 'L- Stimulus', 'Response to L-');
title('Biphasic responses')
xlabel('Time (ms)'); set(gca, 'Tickdir', 'out', 'Ylim', [-3, 3]); axis square;

set(gcf, 'renderer', 'painters');
plot_counter = plot_counter + 1;

%% Reviewer analysis: Rebound OFF responses-part 2
if ~exist('plot_counter')
    plot_counter = 1;
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

% Load the isoresponse data
load RSSE_linearmodel_CV.mat % Robust regression
load RSSE_quadmodel_CV.mat

% For storing median of differences/ratios
RSSEisoresp_medianofratios = [];
Isoresponse_NLI = [];

for ii = 1:numel(RSSE_linearmodel)   
    % computation for calculating median of differences/ratios
    RSSEisoresp_medianofratios = [RSSEisoresp_medianofratios; median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii})];  
    Isoresponse_NLI = [Isoresponse_NLI; log10(RSSEisoresp_medianofratios(end))];
    
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


% Measuring the OFF-rebound responses from the Isoresponse phase
binwidth = 0.010;
N = numel(filename);
t1 = 0.15; t2 = 0.15;
PSTHbins = -t1:binwidth:t2;
PSTH = cell(1,N);
OFFrespindices= find(PSTHbins>=0 & PSTHbins<=0.15); % from STIMOFFCD to STIMOFFCD + 150 ms
afterOFFrespindices= find(PSTHbins>0.15);
eyepossamplingrate = 500; % in Hz
sampligrate = 1000; % in Hz for spikes 
baselineFR = [];
OFF_FR = [];
pOFF = [];
rOFF = [];
plot_counter = 1;
numsubplots = ceil(sqrt(N));
plotrasters = 0;
for ii = 1:N
    ind = ii;
    % Plotting the OFF responses
    fileofinterest = char(filename(ind,:));
    stro = nex2stro(findfile(fileofinterest));
    spikename = 'sig001a';%getSpikenum(stro);
    neurothreshidx = strcmp(stro.sum.trialFields(1,:),'neurothresh');
    maskidx = strcmp(stro.sum.rasterCells(1,:),'subunit_mask');
    basisvecidx = strcmp(stro.sum.rasterCells(1,:),'basis_vec');
    stimoffidx = strcmp(stro.sum.trialFields(1,:),'stim_off');
    latencyidx = strcmp(stro.sum.trialFields(1,:),'latency');
    stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on');
    fpacqidx = strcmp(stro.sum.trialFields(1,:),'fpacq');
    spikestampidx = find(strcmp(stro.sum.rasterCells(1,:),spikename));
    mask_changes = [2];
    all_masks = stro.ras(:,maskidx);
    Fx = @(xi) any(isnan(xi)); % function that finds 'NaN' in a cell array
    inds = find(cellfun(Fx,stro.ras(:,basisvecidx))==0);
    if isempty(inds)
        inds = size(stro.trial,1)-1;
    end
    t_offset = stro.trial(end,latencyidx)/1000;
    neurothreshmode = stro.trial(:,neurothreshidx);
    basisvec_dropidx = inds(end);
    neurothresh_startidx = max(find(neurothreshmode==1,1),basisvec_dropidx+1);
    baselinespikecounts = []; % from stimon - 150ms to stimon
    spikecountsduringstimpresent = []; % from stimon to stimoff
    spikecountsafter = []; % from stimoff to stimoff + 150ms
    idxs = neurothresh_startidx:1:size(stro.trial,1);
    %t_offset = stro.trial(end,latencyidx)/1000;
    t_offset = 0.05; % hardocding it to 50 ms
    
    num_dur = [];
    for jj = 1:numel(idxs)
        tmp = cell2mat(stro.ras(idxs(jj),spikestampidx));
        baselinespikecounts = [baselinespikecounts; numel(tmp(tmp<stro.trial(idxs(jj),stimonidx) & tmp>stro.trial(idxs(jj),stimonidx)-0.15))];
        spikecountsduringstimpresent = [spikecountsduringstimpresent; numel(tmp(tmp<stro.trial(idxs(jj),stimoffidx) & tmp>stro.trial(idxs(jj),stimonidx)+t_offset))];
        spikecountsafter = [spikecountsafter; numel(tmp(tmp<stro.trial(idxs(jj),stimoffidx)+0.15 & tmp>stro.trial(idxs(jj),stimoffidx)+t_offset))];
        num_dur = [num_dur; stro.trial(idxs(jj),stimonidx)-stro.trial(idxs(jj),fpacqidx)];    
    end
    [newspikecounts,newidxs] = sort(spikecountsduringstimpresent);
    
    [r,p] = corr(spikecountsduringstimpresent,spikecountsafter,'type','Spearman');
    baselineFR = [baselineFR; mean(baselinespikecounts(spikecountsduringstimpresent==0))/0.15];
    OFF_FR = [OFF_FR; mean(spikecountsafter(spikecountsduringstimpresent==0))/(0.15-t_offset)];
    rOFF = [rOFF; r];
    pOFF = [pOFF; p];
    
    count = 1; spcounts = [];
    % for plotting purposes
    if plotrasters
        analogstridx = find(strcmp(stro.sum.rasterCells(1,:),'anlgStartTime'));
        idxs = idxs(newidxs);
        spcount = [];
        for jj = 1:numel(idxs)
            tmp = cell2mat(stro.ras(idxs(jj),spikestampidx));
            spikes = tmp(tmp<stro.trial(idxs(jj),stimoffidx)+t2 & tmp>stro.trial(idxs(jj),stimonidx)-t1);
            spikes = spikes - stro.trial(idxs(jj),stimonidx);
            figure(plot_counter);subplot(numsubplots,numsubplots,ii); plot(spikes,(count)*ones(1,length(spikes)),'k.'); hold on;
            count = count + 1;
            spcount = [spcount; numel(tmp(tmp<stro.trial(idxs(jj),stimoffidx) & tmp>stro.trial(idxs(jj),stimonidx)))];
        end
        if any(DOidx == ind)
            c = 'r';
        elseif any(LUMidx == ind)
            c = 'g';
        else
            c = 'k';
        end
        set(gca,'XColor',c,'YColor',c,'Xlim',[-0.2 0.45],'Ylim',[0 count]);
        line([0 0],[0 numel(idxs)],'Color',[1 0 0]);
        line([0.3 0.3],[0 numel(idxs)],'Color',[1 0 0]);
        line([-0.15 -0.15],[0 numel(idxs)],'Color',[1 0 0]);
        line([0.45 0.45],[0 numel(idxs)],'Color',[1 0 0]); drawnow; hold off;
    end
    
end

newLUMidx = LUMidx(baselineFR(LUMidx)>0);
newDOidx = DOidx(baselineFR(DOidx)>0);
newhardtoclassifyidx = hardtoclassifyidx(baselineFR(hardtoclassifyidx)>0);

% Some stats on the OFF response and the Isoresponse NLI
[r2,p2] = corr(Isoresponse_NLI([LUMidx, DOidx, hardtoclassifyidx]),OFF_FR([LUMidx, DOidx, hardtoclassifyidx]),'type','Spearman', 'rows','complete');
[r3,p3] = corr(Isoresponse_NLI([LUMidx, DOidx, hardtoclassifyidx]),OFF_FR([LUMidx, DOidx, hardtoclassifyidx]),'type','Pearson', 'rows','complete');

% Only considering cells where the baseline FR is > 0
[r4,p4] = corr(Isoresponse_NLI([newLUMidx, newDOidx, newhardtoclassifyidx]),OFF_FR([newLUMidx, newDOidx, newhardtoclassifyidx]),'type','Spearman', 'rows','complete');


% Plotting the result
indices = [109 24 74];
figure(plot_counter);
plot(Isoresponse_NLI(LUMidx), OFF_FR(LUMidx), 'o', 'MarkerFaceColor', [0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(Isoresponse_NLI(indices(1)), OFF_FR(indices(1)), 'o', 'MarkerFaceColor', [0 0 0],'MarkerEdgeColor',[0 1 0]);
plot(Isoresponse_NLI(DOidx), OFF_FR(DOidx), 'o', 'MarkerFaceColor', [1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(Isoresponse_NLI(indices(2)), OFF_FR(indices(2)), 'o', 'MarkerFaceColor', [1 0 0],'MarkerEdgeColor',[0 1 0]);
plot(Isoresponse_NLI(hardtoclassifyidx), OFF_FR(hardtoclassifyidx), 'o', 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
plot(Isoresponse_NLI(indices(3)), OFF_FR(indices(3)), 'o', 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor',[0 1 0]);
set(gca,'Tickdir','out','Xlim',[-1 2],'XTick',[-1 0 1 2],'Ylim',[0 100],'YTick',[0:20:100]); 
xlabel('Isoresponse NLI'); ylabel('OFF response'); axis square; hold off;
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;


%% 3. Reviewer analysis: Eye movement (median eye displacement) vs. Isoresponse NLIs
% This part is similar to the Figure 9 of Jneurophysiology paper
% The code of computing the median eye displacement can be found in 
% Spatialstructureanalysis_GregWN/PopAnalysis_createOutputList.m

if ~exist('plot_counter')
    plot_counter = 1;
end

% Laoding data for cone weight calculation  
load conewts_svd.mat
load vals.mat
load angulardifferences_RGB.mat
anglebwvectors = angulardifference_RGB;
SpatiallyOpponent = anglebwvectors'>90;

% Classifying cells based on cone-weights and PC1 signficance
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

% Load the isoresponse data
load RSSE_linearmodel_CV.mat % Robust regression
load RSSE_quadmodel_CV.mat

% For storing median of differences/ratios
RSSEisoresp_medianofratios = [];
Isoresponse_NLI = [];
num_isoresponse_pts = [];

for ii = 1:numel(RSSE_linearmodel)   
    % computation for calculating median of differences/ratios
    RSSEisoresp_medianofratios = [RSSEisoresp_medianofratios; median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii})];  
    Isoresponse_NLI = [Isoresponse_NLI; log10(RSSEisoresp_medianofratios(end))];  
    num_isoresponse_pts = [num_isoresponse_pts; numel(RSSE_quadmodel{ii})];
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

% Measuring the OFF-rebound responses from the Isoresponse phase
N = numel(filename);
plotrasters = 0;
median_eye_displacement = [];
subject_ID  = [];
for ii = 1:N
    disp(ii);
    
    subject_ID = [subject_ID; filename{ii}(1)];
    % Extracting parameters from the stro structure
    fileofinterest = char(filename(ii,:));
    stro = nex2stro(findfile(fileofinterest));
    spikename = 'sig001a';
    neurothreshidx = strcmp(stro.sum.trialFields(1,:),'neurothresh');
    maskidx = strcmp(stro.sum.rasterCells(1,:),'subunit_mask');
    basisvecidx = strcmp(stro.sum.rasterCells(1,:),'basis_vec');
    stimoffidx = strcmp(stro.sum.trialFields(1,:),'stim_off');
    latencyidx = strcmp(stro.sum.trialFields(1,:),'latency');
    stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on');
    fpacqidx = strcmp(stro.sum.trialFields(1,:),'fpacq');
    spikestampidx = find(strcmp(stro.sum.rasterCells(1,:),spikename));
    anlgStartTimeidx = find(strcmp(stro.sum.rasterCells(1,:),'anlgStartTime'));
    hepidx = find(strcmp(stro.sum.rasterCells(1,:),'AD11'));
    vepidx = find(strcmp(stro.sum.rasterCells(1,:),'AD12'));
    Fx = @(xi) any(isnan(xi)); % function that finds 'NaN' in a cell array
    inds = find(cellfun(Fx,stro.ras(:,basisvecidx))==0);
    if isempty(inds)
        inds = size(stro.trial,1)-1;
    end
    
    if ~isempty(stro.sum.analog.storeRates)
        samplingrate = stro.sum.analog.storeRates{1};
        H = [];
        V = [];
        
        % Looping only the Isoresponse part
        for ll = inds+1:size(stro.trial,1)

            anlgStarttime = stro.ras{ll,anlgStartTimeidx};
            stimont = stro.trial(ll,stimonidx);
            stimofft = stro.trial(ll,stimoffidx);

            h_eyepos = stro.ras{ll,hepidx}*4096/400;
            v_eyepos = stro.ras{ll,vepidx}*4096/400;
            t = anlgStarttime: (1/samplingrate): anlgStarttime+(numel(h_eyepos)-1)*(1/samplingrate);
            tmp = t>=stimont & t<=stimofft;

            % Storing the H and V position
            H = [H; h_eyepos(tmp)];
            V = [V; v_eyepos(tmp)];
        end
        % Storing the median amplitude of the eye position
        Amp = sqrt((H-mean(H)).^2 + (V-mean(V)).^2);
        median_eye_displacement = [median_eye_displacement; median(Amp)];
    else
        median_eye_displacement = [median_eye_displacement; nan];
    end
 
end

% Median eye displacement vs. Isoresponse NLI
indices = [109 24 74];
figure(plot_counter); 
plot(median_eye_displacement(LUMidx),Isoresponse_NLI(LUMidx),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(median_eye_displacement(indices(1)),Isoresponse_NLI(indices(1)),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 1 0]); hold on;
plot(median_eye_displacement(DOidx),Isoresponse_NLI(DOidx),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(median_eye_displacement(indices(2)),Isoresponse_NLI(indices(2)),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 1 0]); hold on;
plot(median_eye_displacement(hardtoclassifyidx),Isoresponse_NLI(hardtoclassifyidx),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
plot(median_eye_displacement(indices(3)),Isoresponse_NLI(indices(3)),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0 1 0]); hold on;
axis square; set(gca,'Tickdir','out','Ylim',[-1 2],'YTick',-1:1:2,'Xlim',[0.05 0.20],'XTick',0.05:0.05:0.20); 
ylabel('Isoresponse NLI'); xlabel('Median eye displacement');
plot_counter = plot_counter + 1;


%  Spearman's correlation coefficient between Isoresponse NLI and median eye displacement
[r,p] = corr(Isoresponse_NLI([LUMi dx, DOidx, hardtoclassifyidx]), median_eye_displacement([LUMidx, DOidx, hardtoclassifyidx]), 'type', 'Spearman');
[r1,p1] = corr(Isoresponse_NLI([LUMidx, DOidx, hardtoclassifyidx]), num_isoresponse_pts([LUMidx, DOidx, hardtoclassifyidx]), 'type', 'Spearman');
[r2,p2] = corr(median_eye_displacement([LUMidx, DOidx, hardtoclassifyidx]), num_isoresponse_pts([LUMidx, DOidx, hardtoclassifyidx]), 'type', 'Spearman');


% Files for Pangu
ind_interest = logical(zeros(numel(Isoresponse_NLI),1));
ind_interest([LUMidx, DOidx, hardtoclassifyidx]) = 1;
[r3,p3] = corr(Isoresponse_NLI(ind_interest & subject_ID=='P'), median_eye_displacement(ind_interest & subject_ID=='P'), 'type', 'Spearman');
[r4,p4] = corr(Isoresponse_NLI(ind_interest & subject_ID=='M'), median_eye_displacement(ind_interest & subject_ID=='M'), 'type', 'Spearman');



%% 4. Reviewer analysis: Reliability of NLI, using Jacknife procedure

if ~exist('plot_counter')
    plot_counter = 1;
end

% Laoding data for cone weight calculation  
load conewts_svd.mat
load vals.mat
load angulardifferences_RGB.mat
anglebwvectors = angulardifference_RGB;
SpatiallyOpponent = anglebwvectors'>90;

% Classifying cells based on cone-weights and PC1 signficance
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

% Load the integration within the subunit data
load AUROClinS1_CV.mat
load AUROCquadS1_CV.mat
load AUROClinS2_CV.mat
load AUROCquadS2_CV.mat

% Load the integration within the subunit whitenoise analysis data
load AUROClinsubunits_CV.mat
load AUROCquadsubunits_CV.mat

% Load the isoresponse data
load RSSE_linearmodel_CV.mat % Robust regression
load RSSE_quadmodel_CV.mat


Within_subunit_NLI = []; % Storing the within subunit NLI
Isoresponse_NLI = []; % For storing the Isoresponse NLI
Whitenoise_NLI = []; % For storing the white noise NLI

JK_error_whitenoise = [];
JK_error_isoresponse = [];

for ii = 1:numel(AUROClin1) 
    
    % GLM error for both the subunits
    Error_lin_subunit1 = 1-AUROClin1{ii};
    Error_lin_subunit2 = 1-AUROClin2{ii};
    
    % GQM error for both the subunits
    Error_quad_subunit1 = 1-AUROCquad1{ii};
    Error_quad_subunit2 = 1-AUROCquad2{ii};
    
    % Calculating the median ratio of errors both the subunits in log scale
    median_NLI_subunit1 = log10(median(Error_lin_subunit1./ Error_quad_subunit1));
    median_NLI_subunit2 = log10(median(Error_lin_subunit2./ Error_quad_subunit2));
    
    % Storing the NLIs
    Within_subunit_NLI = [Within_subunit_NLI; median([median_NLI_subunit1 median_NLI_subunit2])];
    
    % White noise NLI 
    Error_quad = 1-(AUROCquadsubunits{ii});
    Error_lin = 1-(AUROClinsubunits{ii}); 
    Whitenoise_NLI = [Whitenoise_NLI; log10(median(Error_lin./Error_quad))];
    
    % Jackknife White noise error
    X1 = log10(jackknife(@median, Error_lin./Error_quad));
    JK_error_whitenoise = [JK_error_whitenoise; std(X1)];
    
    % Isoresponse NLI
    Isoresponse_NLI = [Isoresponse_NLI; log10(median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii}))];
    
    % Jackknife Isoresponse error
    X2 = log10(jackknife(@median, RSSE_linearmodel{ii}./RSSE_quadmodel{ii}));
    JK_error_isoresponse = [JK_error_isoresponse; std(X2)];
end


%% 5. A control analysis: Some more analyses that Greg suggested (a continuation of Figure 3 but not for the paper)
%  1) To make the isoprobability NLI and isoresponse NLI definition more consistent 
%  2) To split the other cell category into a) significant PC1 b) non-significant PC1

if ~exist('plot_counter')
    plot_counter = 1;
end

load conewts_svd.mat
load vals.mat
load S1RGB_svd.mat
load S2RGB_svd.mat
load angulardifferences_RGB.mat
anglebwvectors = angulardifference_RGB;
S1RGB = S1RGB_svd;
S2RGB = S2RGB_svd;
SpatiallyOpponent = anglebwvectors'>90;

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); 
Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

LUMidx = LumIds_conewts;
DOidx = [ColorOpponentIds_conewts Sconedominated_conewts];
LUMidx = LUMidx(vals(LUMidx)<95);
DOidx = DOidx(vals(DOidx)<95);

% Introducing new categories based on Greg's suggestion
hardtoclassifyidx = 1:size(conewts_svd,2);
hardtoclassifyidx([LUMidx DOidx]) = [];
hardtoclassifyidx_woPC1 = [hardtoclassifyidx(vals(hardtoclassifyidx)<95)];
hardtoclassifyidx_wPC1 = [hardtoclassifyidx(vals(hardtoclassifyidx)>=95)]; 

LUMidx = LUMidx(SpatiallyOpponent(LUMidx));
DOidx = DOidx(SpatiallyOpponent(DOidx));
hardtoclassifyidx_woPC1 = hardtoclassifyidx_woPC1(SpatiallyOpponent(hardtoclassifyidx_woPC1));
hardtoclassifyidx_wPC1 = hardtoclassifyidx_wPC1(SpatiallyOpponent(hardtoclassifyidx_wPC1)); 

% Checking the correlation with non-linearity indices 
% Load the isoresponse data
load RSSE_linearmodel_CV.mat % Robust regression
load RSSE_quadmodel_CV.mat

% Load the integration within the subunit data
load AUROClinsubunits_CV.mat
load AUROCquadsubunits_CV.mat

% For storing median of differences/ratios
Acrosssubunits_medianofdifferences = [];
Acrosssubunits_lin_median = []; Acrosssubunits_quad_median = [];

% For storing median of differences/ratios
RSSEisoresp_medianofratios = [];
RSSEisoresp_lin_median = []; RSSEisoresp_quad_median = []; % Isoresponse data
Isoresponse_NLI = [];

for ii = 1:numel(RSSE_linearmodel)   
    % Isoresponse data - computation for calculating median of differences/ratios
    RSSEisoresp_medianofratios = [RSSEisoresp_medianofratios; median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii})]; 
    Isoresponse_NLI = [Isoresponse_NLI; log10(RSSEisoresp_medianofratios(end))];
    RSSEisoresp_lin_median = [RSSEisoresp_lin_median; median(RSSE_linearmodel{ii})];
    RSSEisoresp_quad_median = [RSSEisoresp_quad_median; median(RSSE_quadmodel{ii})];
    
    % Storing the WN subunit spatial interaction data 
    Error_quad = 1-(AUROCquadsubunits{ii});
    Error_lin = 1-(AUROClinsubunits{ii});
    
    % A new definition of Acrosssubunits_medianofdifferences
    Acrosssubunits_medianofdifferences = [Acrosssubunits_medianofdifferences; median(Error_lin./Error_quad)];
    Acrosssubunits_lin_median = [Acrosssubunits_lin_median; median(Error_lin)];
    Acrosssubunits_quad_median = [Acrosssubunits_quad_median; median(Error_quad)];
    
end
RSSEisoresp_medianofratios(RSSEisoresp_medianofratios<0.1) = 0.1;
indices = [109 31 74];

% Isoresponse data: Plotting the results for SVD based cone weight classification including the PC1 z-scores 
figure(plot_counter);

subplot(121); histogram(log10(Acrosssubunits_medianofdifferences(hardtoclassifyidx_woPC1)),-0.02:0.005:0.1,'DisplayStyle','stairs','EdgeColor',[1 0 0.5],'Linewidth',2); hold on;
histogram(log10(Acrosssubunits_medianofdifferences(hardtoclassifyidx_wPC1)),-0.02:0.005:0.1,'DisplayStyle','stairs','EdgeColor',[1 0.5 0],'Linewidth',2); hold on;
plot(log10(median(Acrosssubunits_medianofdifferences(hardtoclassifyidx_woPC1))),13,'v','MarkerSize',8,'MarkerFaceColor',[1.0 0 0.5],'MarkerEdgeColor',[1 1 1]); 
plot(log10(median(Acrosssubunits_medianofdifferences(hardtoclassifyidx_wPC1))),14,'v','MarkerSize',8,'MarkerFaceColor',[1.0 0.5 0],'MarkerEdgeColor',[1 1 1]); 
set(gca,'Tickdir','out','Xlim',[-0.02 0.1],'XTick',-0.02:0.02:0.1,'Ylim',[0 15],'YTick',[0 5 10 15]); 
ylabel('# cells'); title('White noise spatial integration'); xlabel('White noise NLI'); legend('OSO woPC1','OSO wPC1'); axis square; hold off;

subplot(122); histogram(log10(RSSEisoresp_medianofratios(hardtoclassifyidx_woPC1)),-1:0.1:3,'DisplayStyle','stairs','EdgeColor',[1 0 0.5],'Linewidth',2); hold on;
histogram(log10(RSSEisoresp_medianofratios(hardtoclassifyidx_wPC1)),-1:0.1:3,'DisplayStyle','stairs','EdgeColor',[1 0.5 0],'Linewidth',2); hold on;
plot(log10(median(RSSEisoresp_medianofratios(hardtoclassifyidx_woPC1))),13,'v','MarkerSize',8,'MarkerFaceColor',[1.0 0 0.5],'MarkerEdgeColor',[1 1 1]); 
plot(log10(median(RSSEisoresp_medianofratios(hardtoclassifyidx_wPC1))),14,'v','MarkerSize',8,'MarkerFaceColor',[1.0 0.5 0],'MarkerEdgeColor',[1 1 1]); 
set(gca,'Tickdir','out','Xlim',[-1 3],'XTick',-1.0:0.5:3.0,'Ylim',[0 15],'YTick',[0 5 10 15]); 
ylabel('# cells'); title('Isoresponse spatial integration'); xlabel('Isoresponse NLI'); legend('OSO woPC1','OSO wPC1'); axis square; hold off;
set(gcf, 'renderer', 'painters')
plot_counter = plot_counter + 1;



% Some more stats
[p1,~] = ranksum(Isoresponse_NLI([hardtoclassifyidx_woPC1]),Isoresponse_NLI([hardtoclassifyidx_wPC1])); % Isoresponse NLI
[p2,~] = ranksum(log10(Acrosssubunits_medianofdifferences([hardtoclassifyidx_woPC1])),log10(Acrosssubunits_medianofdifferences([hardtoclassifyidx_wPC1]))); % Isoresponse NLI


%% 6. A nx7 matrix for Greg 
% L-cone weight subunit 1
% M-cone weight subunit 1
% S-cone weight subunit 1
% L-cone weight subunit 2
% M-cone weight subunit 2
% S-cone weight subunit 2
% Cell type (simple, DO, NSNDO)

load S1LMS_svd.mat
load S2LMS_svd.mat
load conewts_svd.mat
load vals.mat
load angulardifferences_RGB.mat
anglebwvectors = angulardifference_RGB;
SpatiallyOpponent = anglebwvectors'>90;

% Classifying cells based on cone-weights and PC1 signficance
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


% Creating the nx7 matrix
ID = [ones(numel(LUMidx),1); 2*ones(numel(DOidx),1); 3*ones(numel(hardtoclassifyidx),1)];
CONEWTS = [S1LMS_svd(:, [LUMidx, DOidx, hardtoclassifyidx])' S2LMS_svd(:, [LUMidx, DOidx, hardtoclassifyidx])' ID];
save('CONEWTS.mat', 'CONEWTS')


%% 7. Inferring linearity of spatial summation from hyperpixel STA

% Showing that it is possible to make a model neuron that does not have 
% a significant PC1 but appears nonlinear in the white noise NLI analysis.
 
% Simulation 1
% Two linear subunits, half rectified and multiplied.
% PC1 criterion will not be triggered, but white noise NLI will show the
% nonlinearity. This was a kind of nonlinearity we were concerned about.
 
n = 1000000;
nbins = 20;
stim = normrnd(0,1,n,6);
kernel1 = [.5 -.5 .2];
kernel2 = -kernel1;
 
% Simulation 1
lingen = stim*[[kernel1';0;0;0], [0;0;0;kernel2']];
resp = prod(max(lingen+1,0),2);
resp = (resp-min(resp))./(max(resp)-min(resp));
resp = resp>unifrnd(zeros(n,1),ones(n,1));
 
STA = resp'*stim;
[v,d] = eig(cov(stim(resp,:)));
figure; subplot(2,1,1);
plot(flipud(diag(d)),'ko','MarkerFaceColor','black')
 
clear allstim spikestim
allstim(:,1)=stim*[STA(1:3),0 0 0]';
allstim(:,2)=stim*[0 0 0 STA(4:6)]';
spikestim(:,1)=stim(resp,:)*[STA(1:3),0 0 0]';
spikestim(:,2)=stim(resp,:)*[0 0 0 STA(4:6)]';
 
bins = [nbins prctile(allstim(:),5) prctile(allstim(:),95)]';
 
[allstim_hist,~]=hist2(allstim,[bins bins]);
[spikestim_hist,~]=hist2(spikestim,[bins bins]);
subplot(2,1,2);
imagesc(spikestim_hist./allstim_hist*255); colormap(gray(255))
axis square;
set(gca,'Xtick',[],'Ytick',[]);

% Fit a GLM & GQM after the frames have been projected onto the subunits
mdllin1 =  fitglm(lingen,resp,'linear','Distribution','binomial','Link','logit');
mdlquad1 =  fitglm(lingen,resp,'quadratic','Distribution','binomial','Link','logit');
predlin1 = predict(mdllin1,lingen); % perdiction from GLM
predquad1 = predict(mdlquad1,lingen); % perdiction from GQM
Error_lin1 = 1-rocN(predlin1(resp),predlin1(~resp));
Error_quad1 = 1- rocN(predquad1(resp),predquad1(~resp)); 
WhiteNoise_NLI_1 = log10(Error_lin1/Error_quad1);


 
% Simulation 2
% Conversely: Two subunits with a halfwave rectified response to one color
% channel and a fullwave rectified response to another. Added linearly.
% This cell will have a significant PC1 but the GQM will not fit better
% than the GLM (because the fitting is in 2D).
lingen = [];
lingen(:,1) = stim(:,1).*kernel1(1)+(stim(:,2).*kernel1(2)).^2+stim(:,3).*kernel1(3); % nonlinear subfield
lingen(:,2) = stim(:,[4:6])*kernel2'; % linear subfield
resp = sum(lingen,2);
resp = (resp-min(resp))./(max(resp)-min(resp));
resp = resp>unifrnd(zeros(n,1),ones(n,1));
 
STA = resp'*stim;
[v,d] = eig(cov(stim(resp,:)));
figure; subplot(2,1,1);
plot(flipud(diag(d)),'ko','MarkerFaceColor','black')
 
clear allstim spikestim
allstim(:,1)=stim*[STA(1:3),0 0 0]';
allstim(:,2)=stim*[0 0 0 STA(4:6)]';
spikestim(:,1)=stim(resp,:)*[STA(1:3),0 0 0]';
spikestim(:,2)=stim(resp,:)*[0 0 0 STA(4:6)]';
 
bins = [nbins prctile(allstim(:),5) prctile(allstim(:),95)]';
 
[allstim_hist,~]=hist2(allstim,[bins bins]);
[spikestim_hist,~]=hist2(spikestim,[bins bins]);
subplot(2,1,2);
imagesc(spikestim_hist./allstim_hist*255); colormap(gray(255))
axis square;
set(gca,'Xtick',[],'Ytick',[]);

% Fit a GLM & GQM after the frames have been projected onto the subunits
mdllin2 =  fitglm(lingen,resp,'linear','Distribution','binomial','Link','logit');
mdlquad2 =  fitglm(lingen,resp,'quadratic','Distribution','binomial','Link','logit');
predlin2 = predict(mdllin2,lingen); % perdiction from GLM
predquad2 = predict(mdlquad2,lingen); % perdiction from GQM
Error_lin2 = 1-rocN(predlin2(resp),predlin2(~resp));
Error_quad2 = 1- rocN(predquad2(resp),predquad2(~resp)); 
WhiteNoise_NLI_2 = log10(Error_lin2/Error_quad2);


%% 8. Scatterplot of firing rates between the pixel WN and hyperpixel WN noise

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


% Classifying cells into simple, DO and hardtoclassify cells
load conewts_svd.mat
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


% Matrix for storing the firing rates (mean and std)
mean_firing_rates = zeros(2, numel(filename));
std_firig_rates = zeros(2, numel(filename));
for ii= 1:numel(filename)
    fileofinterest = char(filename(ii,:));
    disp(fileofinterest);
    stro = nex2stro(findfile(fileofinterest));
    spikename = 'sig001a';%getSpikenum(stro);
    maskidx = strcmp(stro.sum.rasterCells(1,:),'subunit_mask');
    spikeidx = strcmp(stro.sum.rasterCells(1,:),spikename);
    basisvecidx = strcmp(stro.sum.rasterCells(1,:),'basis_vec');
    weightsidx = strcmp(stro.sum.rasterCells(1,:),'weights');
    parentverticesidx = strcmp(stro.sum.rasterCells(1,:),'parentvertices');
    nstixperside = stro.sum.exptParams.nstixperside;
    ngammasteps = 2^16; % 65536
    linepredtol = stro.sum.exptParams.linepredtol;
    stepsizescale = stro.sum.exptParams.stepsizescale;
    stepsize = stro.sum.exptParams.stepsize;
    nreversals = stro.sum.exptParams.nreversals;
    oogscale = stro.sum.exptParams.oogscale;
    seedidx = strcmp(stro.sum.trialFields(1,:),'seed');
    nframesidx = strcmp(stro.sum.trialFields(1,:),'num_frames');
    stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on');
    stimoffidx = strcmp(stro.sum.trialFields(1,:),'stim_off');
    fponidx = strcmp(stro.sum.trialFields(1,:),'fp_on');
    fpacqidx = strcmp(stro.sum.trialFields(1,:),'fpacq');
    basisvecdiridx = strcmp(stro.sum.trialFields(1,:),'weights_idx');
    neurothreshidx = strcmp(stro.sum.trialFields(1,:),'neurothresh'); % when exactly the neurothresh trials started
    targetspikerateidx = strcmp(stro.sum.trialFields(1,:),'targetspikerate');
    correctidx = strcmp(stro.sum.trialFields(1,:),'correct');
    muidxs = [find(strcmp(stro.sum.trialFields(1,:),'mu1')), ...
        find(strcmp(stro.sum.trialFields(1,:),'mu2')), ...
        find(strcmp(stro.sum.trialFields(1,:),'mu3'))];
    sigmaidxs = [find(strcmp(stro.sum.trialFields(1,:),'sigma1')), ...
        find(strcmp(stro.sum.trialFields(1,:),'sigma2')), ...
        find(strcmp(stro.sum.trialFields(1,:),'sigma3'))];
    latencyidx = strcmp(stro.sum.trialFields(1,:),'latency');
    reversalflagidx = strcmp(stro.sum.trialFields(1,:),'reversalflag');
    msperframe = 1000/stro.sum.exptParams.framerate;
    ntrials = size(stro.trial,1);
    maxT = 15; % this represents the temporal part in the spatiotemporal receptive field
    xx = linspace(stro.sum.exptParams.gauss_locut/1000, stro.sum.exptParams.gauss_hicut/1000,ngammasteps); % xx represents the probabilities. For more info, have a look at the MATLAB 'norminv' function.
    yy = norminv(xx'); % defining norminv to extract the values for which the cdf values range between gauss_locut and gauss_hicut

    
    mask_changes = [2];
    all_masks = stro.ras(:,maskidx);
    Fx = @(xi) any(isnan(xi)); % function that finds 'NaN' in a cell array
    inds = find(cellfun(Fx,stro.ras(:,basisvecidx))==0);
    if isempty(inds)
        inds = size(stro.trial,1)-1;
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
    
    % Calculate Firing rates in pixel WN and hyperpixel WN
    for index = 1:2
        num_spikes =[];
        num_dur = [];
        for jj = mask_changes(1,index):mask_changes(2,index)
            nframes = stro.trial(jj,nframesidx);
            
            if nframes >0
                t_stimon = stro.trial(jj, stimonidx);
                num_spikes = [num_spikes; numel((stro.ras{jj,spikeidx}-t_stimon)*1000)];

                frametimes = linspace(0, nframes*msperframe, nframes)+(msperframe/2)';
                num_dur = [num_dur; range(frametimes)/1000];
            end

        end
        
        % Storing the mean and standard deviation of firing rates
        fr = num_spikes./num_dur;
        mean_firing_rates(index,ii) = mean(fr);
        %std_firig_rates(index,ii) = sem(fr);
        
        M = mean(fr);
        V = var(fr);
        n = numel(fr);
        std_firig_rates(index,ii) = sqrt(V./n./((log(10)*M).^2));
        
    end
    
end

% Comparing the Firing rates from pixel and hyperpixel WN 
FR_max = 1000;
FR_min = 0.3;

figure(plot_counter), hold on;
errorbar(log(mean_firing_rates(1,LUMidx)), log(mean_firing_rates(2,LUMidx)),std_firig_rates(2,LUMidx), std_firig_rates(2,LUMidx), std_firig_rates(1,LUMidx), std_firig_rates(1,LUMidx), 'o', 'MarkerFaceColor', [0 0 0],'MarkerEdgeColor',[1 1 1], 'color', [0 0 0]); 
errorbar(log(mean_firing_rates(1,DOidx)), log(mean_firing_rates(2,DOidx)),std_firig_rates(2,DOidx), std_firig_rates(2,DOidx),std_firig_rates(1,DOidx), std_firig_rates(1,DOidx), 'o', 'MarkerFaceColor', [1 0 0],'MarkerEdgeColor',[1 1 1], 'color', [1 0 0]); 
errorbar(log(mean_firing_rates(1,hardtoclassifyidx)), log(mean_firing_rates(2,hardtoclassifyidx)), std_firig_rates(2,hardtoclassifyidx), std_firig_rates(2,hardtoclassifyidx),  std_firig_rates(1,hardtoclassifyidx), std_firig_rates(1,hardtoclassifyidx),'o', 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1], 'color', [0.5 0.5 0.5]); 
plot(log([FR_min FR_max]), log([FR_min FR_max]), 'color', 'k')
set(gca,'Tickdir','out','Xlim',[log(FR_min) log(FR_max)],'Ylim',[log(FR_min) log(FR_max)], 'XTick', log([0.3 1 10 100 1000]),'YTick', log([0.3 1 10 100 1000]), 'XTickLabel', [0.3 1 10 100 1000],'YTickLabel', [0.3 1 10 100 1000]); 
xlabel('Pixel WN firing rates'); ylabel('Hyperpixel WN firing rates'); axis square; hold off;
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;


%% 9. Testing Jacknife estimation procedure
close all; clearvars;

T = round(logspace(1,4,10)); % # of samples to be tested
sem_original = [];
jackknife_std = [];

for j=1:numel(T)
    
    % selecting N random values from a normal distribution
    N = T(j);
    nums = randn(N,1); 
    
    % comuting the standard error of mean
    sem_original = [sem_original; std(nums)/sqrt(N)]; 
    
    % Computing the jackknife estimate
    jackknife_mean = jackknife(@mean, nums);
    n = numel(jackknife_mean);
    tmp_se = sqrt(sum((jackknife_mean - mean(nums)).^2)*(n-1)/n);
    jackknife_std = [jackknife_std; tmp_se];

end

figure(1); 
scatter(sem_original, jackknife_std, 50, 1:numel(sem_original), 'filled'); 
hold on;
plot([0 0.35], [0 0.35], 'k');
xlabel('SEM');
ylabel('Jacknife error estimate');
set(gca, 'Xlim', [0 0.35], 'Ylim', [0 0.35]);
hold off;

% Some more info on the sample size
disp('Printing the Jacknife error for the following sample size');
disp(T)

%%
close all; clearvars;

T = round(logspace(1,4,10)); % # of samples to be tested
sem_original = [];
jackknife_std = [];
jackknife_std2 = [];

for j=1:numel(T)
    
    % selecting N random values from a normal distribution
    N = T(j);
    nums = randn(N,1); 
    
    % comuting the standard error of mean
    sem_original = [sem_original; std(nums)/sqrt(N)]; 
    
    % Computing the jackknife estimate
    jackknife_mean = jackknife(@mean, nums);
    jackknife_std = [jackknife_std; std(jackknife_mean)];
    tmp = sqrt(sum((jackknife_mean-mean(nums)).^2)*(N-1)/N);
    jackknife_std2 = [jackknife_std2; tmp];

end

figure(1); 
scatter(sem_original, jackknife_std2, 50, 1:numel(sem_original), 'filled'); 
hold on;
plot([0 0.35], [0 0.35], 'k');
xlabel('SEM');
ylabel('Jacknife error estimate' );
set(gca, 'Xlim', [0 0.35], 'Ylim', [0 0.35]);
hold off;

% Some more info on the sample size
disp('Printing the Jacknife error for the following sample size');
disp(T)



%% 10. Send the staircases to Greg for Naka-Rushton fitting
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

global reversalflagidx stepsizescale stepsize nreversals
stro = nex2stro(findfile(char(filename(74))));
spikename = 'sig001a';
maskidx = strcmp(stro.sum.rasterCells(1,:),'subunit_mask');
spikeidx = strcmp(stro.sum.rasterCells(1,:),spikename);
basisvecidx = strcmp(stro.sum.rasterCells(1,:),'basis_vec');
weightsidx = strcmp(stro.sum.rasterCells(1,:),'weights');
parentverticesidx = strcmp(stro.sum.rasterCells(1,:),'parentvertices');
nstixperside = stro.sum.exptParams.nstixperside;
ngammasteps = 2^16; % 65536
linepredtol = stro.sum.exptParams.linepredtol;
stepsizescale = stro.sum.exptParams.stepsizescale;
stepsize = stro.sum.exptParams.stepsize;
nreversals = stro.sum.exptParams.nreversals;
oogscale = stro.sum.exptParams.oogscale;
seedidx = strcmp(stro.sum.trialFields(1,:),'seed');
nframesidx = strcmp(stro.sum.trialFields(1,:),'num_frames');
stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on');
stimoffidx = strcmp(stro.sum.trialFields(1,:),'stim_off');
fponidx = strcmp(stro.sum.trialFields(1,:),'fp_on');
fpacqidx = strcmp(stro.sum.trialFields(1,:),'fpacq');
basisvecdiridx = strcmp(stro.sum.trialFields(1,:),'weights_idx');
neurothreshidx = strcmp(stro.sum.trialFields(1,:),'neurothresh'); % when exactly the neurothresh trials started
targetspikerateidx = strcmp(stro.sum.trialFields(1,:),'targetspikerate');
correctidx = strcmp(stro.sum.trialFields(1,:),'correct');
muidxs = [find(strcmp(stro.sum.trialFields(1,:),'mu1')), ...
    find(strcmp(stro.sum.trialFields(1,:),'mu2')), ...
    find(strcmp(stro.sum.trialFields(1,:),'mu3'))];
sigmaidxs = [find(strcmp(stro.sum.trialFields(1,:),'sigma1')), ...
    find(strcmp(stro.sum.trialFields(1,:),'sigma2')), ...
    find(strcmp(stro.sum.trialFields(1,:),'sigma3'))];
latencyidx = strcmp(stro.sum.trialFields(1,:),'latency');
reversalflagidx = strcmp(stro.sum.trialFields(1,:),'reversalflag');
msperframe = 1000/stro.sum.exptParams.framerate;
ntrials = size(stro.trial,1);
maxT = 15; % this represents the temporal part in the spatiotemporal receptive field
xx = linspace(stro.sum.exptParams.gauss_locut/1000, stro.sum.exptParams.gauss_hicut/1000,ngammasteps); % xx represents the probabilities. For more info, have a look at the MATLAB 'norminv' function.
yy = norminv(xx'); % defining norminv to extract the values for which the cdf values range between gauss_locut and gauss_hicut

mask_changes = [2];
all_masks = stro.ras(:,maskidx);
Fx = @(xi) any(isnan(xi)); % function that finds 'NaN' in a cell array
inds = find(cellfun(Fx,stro.ras(:,basisvecidx))==0);
if isempty(inds)
    inds = size(stro.trial,1)-1;
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

% Determining when Neurothresh mode was active, plotting the basis vector, working correctly
t_offset = stro.trial(end,latencyidx)/1000;
neurothreshmode = stro.trial(:,neurothreshidx);
basisvec_dropidx = inds(end);
neurothresh_startidx = max(find(neurothreshmode==1,1),basisvec_dropidx+1); % There is possibility of pressing the Neurothresh button multiple times
num_targetspikerates = unique(stro.trial(neurothresh_startidx:end,targetspikerateidx));
vect = stro.ras{basisvec_dropidx,basisvecidx};
basisvec_size = nstixperside*nstixperside*3;
numvect = (numel(vect)/basisvec_size)-1;
basisvec = cell(1,numvect);

% Actual basis vec
bkgnd_monitor = reshape(vect(numvect*basisvec_size+1:basisvec_size*(numvect+1)),[nstixperside nstixperside 3]);
for ii = 1:numvect
    tmp_vec = vect((ii-1)*basisvec_size+1:basisvec_size*ii) + vect(numvect*basisvec_size+1:basisvec_size*(numvect+1));
    basisvec{ii} = reshape(tmp_vec,[nstixperside nstixperside 3]);  
end

% Extracting and storing the staircases
cum_flag = []; % to check if I am looking at all the probed directions
dir = 24; % Good directions: 14, 24, 25, 26
dirs = unique(stro.trial(:,basisvecdiridx));
dirs = dirs(2:end); % Because the direction 0 is just the whitenoise phase
staircase.contrast_1D = cell(1,numel(dirs));
staircase.firing_rate = cell(1,numel(dirs));
staircase.contrast_2D = cell(1,numel(dirs));
staircase.target_firing_rate = num_targetspikerates;
staircase.gamut_Violation = cell(1,numel(dirs));
staircase.completely_probed = cell(1,numel(dirs));

for dir_idx = 1:numel(dirs) 
    dir = dirs(dir_idx);
    
    for kk = 1: numel(num_targetspikerates)
        tmp_n = [];
        tmp_wts = [];
        tmp_parentvertices = [];
        idxs1 = find(stro.trial(:,basisvecdiridx) == dir);

        for jj = 1:numel(idxs1)
            tmp_wts = [tmp_wts; stro.ras{idxs1(jj),weightsidx}'];
            tmp_n = [tmp_n; norm(stro.ras{idxs1(jj),weightsidx})];
        end
        [flag, gamutViolation] = Check_ifcompletely_probed(stro,basisvec,bkgnd_monitor,idxs1,tmp_wts(end,:));
        cum_flag = [cum_flag; flag];
        
        % flag = 0, incompletely probed
        % flag = 1, completely probed
        % gamutViolation = 1, out of gamut point
        if gamutViolation == 1
            c = [0 1 0];
        else
            c = [0 0 1];
        end
                
        raster_data = stro.ras(idxs1,spikeidx);
        num_dur =[];
        firing_rate = [];

        for ii = 1:size(raster_data,1)
            tmp = raster_data{ii} ;
            spikes = tmp(tmp<stro.trial(idxs1(ii),stimoffidx) & tmp>stro.trial(idxs1(ii),stimonidx)+t_offset);
            num_dur = [num_dur; (stro.trial(idxs1(ii), stimoffidx)- stro.trial(idxs1(ii),stimonidx))-t_offset];
            firing_rate = [firing_rate; numel(spikes)/num_dur(end)];
        end
        
    end
    
    staircase.contrast_1D{dir_idx} = tmp_n;
    staircase.firing_rate{dir_idx} = firing_rate;
    staircase.contrast_2D{dir_idx} = tmp_wts;
    staircase.gamut_Violation{dir_idx} = gamutViolation;
    staircase.completely_probed{dir_idx} = flag;
    
end

plot_counter = plot_counter + 1;


%% 11. Compare the firing rates between Phase 2 (hyperpixel WN) and Phase 3 (isoresponse)

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


% Classifying cells into simple, DO and hardtoclassify cells
load conewts_svd.mat
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


% Matrix for storing the firing rates (mean and std)
mean_firing_rates_Phase2 = zeros(numel(filename),1);
std_firing_rates_Phase2 = zeros(numel(filename),1);
mean_firing_rates_Phase3 = zeros(numel(filename),1);
std_firing_rates_Phase3 = zeros(numel(filename),1);

for ii= 1:numel(filename)
    fileofinterest = char(filename(ii,:));
    disp(fileofinterest);
    stro = nex2stro(findfile(fileofinterest));
    spikename = 'sig001a';%getSpikenum(stro);
    maskidx = strcmp(stro.sum.rasterCells(1,:),'subunit_mask');
    spikeidx = strcmp(stro.sum.rasterCells(1,:),spikename);
    basisvecidx = strcmp(stro.sum.rasterCells(1,:),'basis_vec');
    weightsidx = strcmp(stro.sum.rasterCells(1,:),'weights');
    parentverticesidx = strcmp(stro.sum.rasterCells(1,:),'parentvertices');
    nstixperside = stro.sum.exptParams.nstixperside;
    ngammasteps = 2^16; % 65536
    linepredtol = stro.sum.exptParams.linepredtol;
    stepsizescale = stro.sum.exptParams.stepsizescale;
    stepsize = stro.sum.exptParams.stepsize;
    nreversals = stro.sum.exptParams.nreversals;
    oogscale = stro.sum.exptParams.oogscale;
    seedidx = strcmp(stro.sum.trialFields(1,:),'seed');
    nframesidx = strcmp(stro.sum.trialFields(1,:),'num_frames');
    stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on');
    stimoffidx = strcmp(stro.sum.trialFields(1,:),'stim_off');
    fponidx = strcmp(stro.sum.trialFields(1,:),'fp_on');
    fpacqidx = strcmp(stro.sum.trialFields(1,:),'fpacq');
    basisvecdiridx = strcmp(stro.sum.trialFields(1,:),'weights_idx');
    neurothreshidx = strcmp(stro.sum.trialFields(1,:),'neurothresh'); % when exactly the neurothresh trials started
    targetspikerateidx = strcmp(stro.sum.trialFields(1,:),'targetspikerate');
    correctidx = strcmp(stro.sum.trialFields(1,:),'correct');
    muidxs = [find(strcmp(stro.sum.trialFields(1,:),'mu1')), ...
        find(strcmp(stro.sum.trialFields(1,:),'mu2')), ...
        find(strcmp(stro.sum.trialFields(1,:),'mu3'))];
    sigmaidxs = [find(strcmp(stro.sum.trialFields(1,:),'sigma1')), ...
        find(strcmp(stro.sum.trialFields(1,:),'sigma2')), ...
        find(strcmp(stro.sum.trialFields(1,:),'sigma3'))];
    latencyidx = strcmp(stro.sum.trialFields(1,:),'latency');
    reversalflagidx = strcmp(stro.sum.trialFields(1,:),'reversalflag');
    msperframe = 1000/stro.sum.exptParams.framerate;
    ntrials = size(stro.trial,1);
    maxT = 15; % this represents the temporal part in the spatiotemporal receptive field
    xx = linspace(stro.sum.exptParams.gauss_locut/1000, stro.sum.exptParams.gauss_hicut/1000,ngammasteps); % xx represents the probabilities. For more info, have a look at the MATLAB 'norminv' function.
    yy = norminv(xx'); % defining norminv to extract the values for which the cdf values range between gauss_locut and gauss_hicut

    
    fundamentals = stro.sum.exptParams.fundamentals; % CONE FUNDAMENTALS: L,M,S
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
    mon_spd = stro.sum.exptParams.mon_spd; % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
    M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals
    M = inv(M');
    mask_changes = [2];
    all_masks = stro.ras(:,maskidx);
    Fx = @(xi) any(isnan(xi)); % function that finds 'NaN' in a cell array
    inds = find(cellfun(Fx,stro.ras(:,basisvecidx))==0);
    if isempty(inds)
        inds = size(stro.trial,1)-1;
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
    
    
    
    % Phase 2: Calculate Firing rates hyperpixel WN
    index = 2;
    num_spikes =[];
    num_dur = [];
    for jj = mask_changes(1,index):mask_changes(2,index)
        nframes = stro.trial(jj,nframesidx);

        if nframes >0
            t_stimon = stro.trial(jj, stimonidx);
            num_spikes = [num_spikes; numel((stro.ras{jj,spikeidx}-t_stimon)*1000)];

            frametimes = linspace(0, nframes*msperframe, nframes)+(msperframe/2)';
            num_dur = [num_dur; range(frametimes)/1000];
        end

    end

    % Storing the mean and standard deviation of firing rates
    fr = num_spikes./num_dur;
    mean_firing_rates_Phase2(ii) = mean(fr);
    M1 = mean(fr);
    V1 = var(fr);
    n1 = numel(fr);
    std_firing_rates_Phase2(ii) = sqrt(V1./n1./((log(10)*M1).^2));
    
    
    % Phase 3: Calculating firing rates for Isoresponse phase 
    spikerates = [];
    for jj = last_wntrial+2:size(stro.trial,1)
        stimont = stro.trial(jj,stimonidx);
        stimofft = stro.trial(jj,stimoffidx);
        t = stimofft-stimont;
        spiketimes = stro.ras{jj,spikeidx};
        nspikes = sum(spiketimes > stimont & spiketimes < stimofft);
        spikerates = [spikerates; nspikes/t];
    end
    
    mean_firing_rates_Phase3(ii) = mean(spikerates);
    M2 = mean(spikerates);
    V2 = var(spikerates);
    n2 = numel(spikerates);
    std_firing_rates_Phase3(ii) = sqrt(V2./n2./((log(10)*M2).^2));
    
end


% Plotting the data
FR_min = 0.3;
FR_max = 300;

figure(plot_counter), hold on;
errorbar(log10(mean_firing_rates_Phase2(LUMidx)), log10(mean_firing_rates_Phase3(LUMidx)),std_firing_rates_Phase2(LUMidx), std_firing_rates_Phase2(LUMidx), std_firing_rates_Phase3(LUMidx), std_firing_rates_Phase3(LUMidx), 'o', 'MarkerFaceColor', [0 0 0],'MarkerEdgeColor',[1 1 1], 'color', [0 0 0]); 
errorbar(log10(mean_firing_rates_Phase2(DOidx)), log10(mean_firing_rates_Phase3(DOidx)),std_firing_rates_Phase2(DOidx), std_firing_rates_Phase2(DOidx),std_firing_rates_Phase3(DOidx), std_firing_rates_Phase3(DOidx), 'o', 'MarkerFaceColor', [1 0 0],'MarkerEdgeColor',[1 1 1], 'color', [1 0 0]); 
errorbar(log10(mean_firing_rates_Phase2(hardtoclassifyidx)), log10(mean_firing_rates_Phase3(hardtoclassifyidx)), std_firing_rates_Phase2(hardtoclassifyidx), std_firing_rates_Phase2(hardtoclassifyidx),  std_firing_rates_Phase3(hardtoclassifyidx), std_firing_rates_Phase3(hardtoclassifyidx),'o', 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1], 'color', [0.5 0.5 0.5]); 
plot(log10([FR_min FR_max]), log10([FR_min FR_max]), 'color', 'k')
set(gca,'Tickdir','out','Xlim',log10([FR_min FR_max]),'Ylim',log10([FR_min FR_max]), 'XTick', log10([0.3 1 10 100 1000]),'YTick', log10([0.3 1 10 100 1000]), 'XTickLabel', [0.3 1 10 100 1000],'YTickLabel', [0.3 1 10 100 1000]); 
xlabel('Phase 2 firing rates'); ylabel('Phase 3 firing rates'); title('Phase 2 vs. Phase 3 FR'); 
axis square; hold off;
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;


%% 12. Extracting the phosphor emission spectra, gamma function, cone fundamentals, and background RGB

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


for ii= 1:numel(1)
    
    fileofinterest = char(filename(ii,:));
    disp(fileofinterest);
    
    stro = nex2stro(findfile(fileofinterest));
    spikename = 'sig001a';%getSpikenum(stro);
    spikeidx = strcmp(stro.sum.rasterCells(1,:),spikename);
    
    % Calculating the cone fundamentals 
    fundamentals = stro.sum.exptParams.fundamentals; % CONE FUNDAMENTALS: L,M,S
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
    mon_spd = stro.sum.exptParams.mon_spd; % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
    M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals
    
    % Gamma table
    gammaTable = stro.sum.exptParams.gamma_table;
    gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
    gammaTable1 = interp1(linspace(0,255,256),gammaTable,linspace(0,255,65536), 'spline');
    invgamma = InvertGamma(gammaTable1, 0);
    
    % background RGB
    ridx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_r'));
    gidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_g'));
    bidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_b'));
    bkgndRGB = [mode(stro.trial(:,ridx)), mode(stro.trial(:,gidx)), mode(stro.trial(:,bidx))];
    bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
    bkgndlms = M*bkgndrgb;
    Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
    Mrgbtocc = inv(Mrgbtocc');
    
end

% Plotting the emission spectra
wavelengths = [380:5:780];
figure(plot_counter);
plot(wavelengths,mon_spd(1,:)'/5, 'r'); hold on;
plot(wavelengths,mon_spd(2,:)'/5, 'g')
plot(wavelengths,mon_spd(3,:)'/5, 'b') 
set(gca, 'Tickdir', 'out'); axis square;
xlabel('Wavelength (nm)'); ylabel('watts/sr?m^2?nm');
plot_counter = plot_counter + 1;


monitor_cone_details.background_rgb = bkgndrgb;
monitor_cone_details.wavelengths = wavelengths;
monitor_cone_details.blue_phosphor_spd = mon_spd(3,:)'/5;
monitor_cone_details.green_phosphor_spd = mon_spd(2,:)'/5;
monitor_cone_details.red_phosphor_spd = mon_spd(1,:)'/5;
monitor_cone_details.Lcone_fundamentals = fundamentals(:,1);
monitor_cone_details.Mcone_fundamentals = fundamentals(:,2);
monitor_cone_details.Scone_fundamentals = fundamentals(:,3);



%% 13.  Comparing reliability between Whitenoise and Isoresponse NLI (Not sure about this analysis)
close all; clearvars

if ~exist('plot_counter')
    plot_counter = 1;
end

load conewts_svd.mat
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

% Checking the correlation with non-linearity indices 
% Load the isoresponse data
load RSSE_linearmodel_CV.mat % Robust regression
load RSSE_quadmodel_CV.mat
% Load the integration within the subunit data
load AUROClinsubunits_CV.mat
load AUROCquadsubunits_CV.mat


Isoresponse_NLI = []; % For storing the Isoresponse NLI
Whitenoise_NLI = []; % For storing the white noise NLI
JK_error_whitenoise = []; % Storing Jackknife based error estimates
JK_error_isoresponse = []; % Storing Jackknife based error estimates

for ii = 1:numel(AUROClinsubunits) 
  
    % White noise NLI 
    Error_quad = 1-(AUROCquadsubunits{ii});
    Error_lin = 1-(AUROClinsubunits{ii});
    Whitenoise_NLI = [Whitenoise_NLI; log10(median(Error_lin./Error_quad))];
    
    % Isoresponse NLI
    Isoresponse_NLI = [Isoresponse_NLI; log10(median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii}))];
    
    % Computing error via Jackknife resampling
    % Jackknife White noise error
    N1 = numel(Error_quad);
    X1 = log10(jackknife(@median, Error_lin./Error_quad));
    var_WN = sum((X1-mean(X1)).^2)*(N1-1);
    JK_error_whitenoise = [JK_error_whitenoise; var_WN];
    
    % Jackknife Isoresponse error
    N2 = numel(RSSE_linearmodel{ii});
    X2 = log10(jackknife(@median, RSSE_linearmodel{ii}./RSSE_quadmodel{ii}));
    var_Iso = sum((X2-mean(X2)).^2)*(N2-1);
    JK_error_isoresponse = [JK_error_isoresponse; var_Iso];
end


%% 14. Comparison of Naka-Rushton fitted isoresponse points to the actual isoresponse points
close all; clearvars;

if ~exist('plot_counter')
    plot_counter = 1;
end

% DO cells
load DO1.fig
load DO2.fig

% Simple cells
load Simple1.fig
load Simple2.fig

% Other spatially opponent cells
load OSO1.fig
load OSO2.fig



%% 15. Robustness of the isoresponse curves based on the last point
% Plotting the average of the last 5 points as well for each probed direction: checking for the example cells only

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
plot_counter = 3;
for zz = 1:numel(indices)
    stro = nex2stro(findfile(char(filename(indices(zz)))));
    global spikename maskidx spikeidx neurothreshidx nstixperside ngammasteps seedidx nframesidx correctidx
    global fponidx stimoffidx stimonidx muidxs sigmaidxs basisvecidx weightsidx fpacqidx targetspikerateidx basisvecdiridx latencyidx reversalflagidx parentverticesidx
    global msperframe ntrials maxT xx yy M linepredtol stepsizescale stepsize nreversals oogscale
    spikename = 'sig001a';
    maskidx = strcmp(stro.sum.rasterCells(1,:),'subunit_mask');
    spikeidx = strcmp(stro.sum.rasterCells(1,:),spikename);
    basisvecidx = strcmp(stro.sum.rasterCells(1,:),'basis_vec');
    weightsidx = strcmp(stro.sum.rasterCells(1,:),'weights');
    parentverticesidx = strcmp(stro.sum.rasterCells(1,:),'parentvertices');
    nstixperside = stro.sum.exptParams.nstixperside;
    ngammasteps = 2^16; % 65536
    linepredtol = stro.sum.exptParams.linepredtol;
    stepsizescale = stro.sum.exptParams.stepsizescale;
    stepsize = stro.sum.exptParams.stepsize;
    nreversals = stro.sum.exptParams.nreversals;
    oogscale = stro.sum.exptParams.oogscale;
    seedidx = strcmp(stro.sum.trialFields(1,:),'seed');
    nframesidx = strcmp(stro.sum.trialFields(1,:),'num_frames');
    stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on');
    stimoffidx = strcmp(stro.sum.trialFields(1,:),'stim_off');
    fponidx = strcmp(stro.sum.trialFields(1,:),'fp_on');
    fpacqidx = strcmp(stro.sum.trialFields(1,:),'fpacq');
    basisvecdiridx = strcmp(stro.sum.trialFields(1,:),'weights_idx');
    neurothreshidx = strcmp(stro.sum.trialFields(1,:),'neurothresh'); % when exactly the neurothresh trials started
    targetspikerateidx = strcmp(stro.sum.trialFields(1,:),'targetspikerate');
    correctidx = strcmp(stro.sum.trialFields(1,:),'correct');
    muidxs = [find(strcmp(stro.sum.trialFields(1,:),'mu1')), ...
        find(strcmp(stro.sum.trialFields(1,:),'mu2')), ...
        find(strcmp(stro.sum.trialFields(1,:),'mu3'))];
    sigmaidxs = [find(strcmp(stro.sum.trialFields(1,:),'sigma1')), ...
        find(strcmp(stro.sum.trialFields(1,:),'sigma2')), ...
        find(strcmp(stro.sum.trialFields(1,:),'sigma3'))];
    latencyidx = strcmp(stro.sum.trialFields(1,:),'latency');
    reversalflagidx = strcmp(stro.sum.trialFields(1,:),'reversalflag');
    msperframe = 1000/stro.sum.exptParams.framerate;
    ntrials = size(stro.trial,1);
    maxT = 10; % this represents the temporal part in the spatiotemporal receptive field
    xx = linspace(stro.sum.exptParams.gauss_locut/1000, stro.sum.exptParams.gauss_hicut/1000,ngammasteps); % xx represents the probabilities. For more info, have a look at the MATLAB 'norminv' function.
    yy = norminv(xx'); % defining norminv to extract the values for which the cdf values range between gauss_locut and gauss_hicut
    t_offset = stro.trial(end,latencyidx)/1000;
    
    % Obtaining the M matrix, code extracted from Greg, fitting a cubic spline
    % using the command 'spline'. 'SplineRaw' only availabe through
    % psychtoolbox which I currently don't have now.
    fundamentals = stro.sum.exptParams.fundamentals; % CONE FUNDAMENTALS: L,M,S
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
    mon_spd = stro.sum.exptParams.mon_spd; % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
    M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals
    M = inv(M');
    
    
    Fx = @(xi) any(isnan(xi)); % function that finds 'NaN' in a cell array
    inds = find(cellfun(Fx,stro.ras(:,basisvecidx))==0);
    if isempty(inds)
        inds = size(stro.trial,1)-1;
    end
    
    % Determining when Neurothresh mode was active, plotting the basis vector, working correctly
    neurothreshmode = stro.trial(:,neurothreshidx);
    basisvec_dropidx = inds(end);
    neurothresh_startidx = max(find(neurothreshmode==1,1),basisvec_dropidx+1); % There is possibility of pressing the Neurothresh button multiple times
    num_targetspikerates = unique(stro.trial(neurothresh_startidx:end,targetspikerateidx));
    vect = stro.ras{basisvec_dropidx,basisvecidx};
    basisvec_size = nstixperside*nstixperside*3;
    numvect = (numel(vect)/basisvec_size)-1;
    basisvec = cell(1,numvect);
    % Actual basis vec
    for ii = 1:numvect
        tmp_vec = vect((ii-1)*basisvec_size+1:basisvec_size*ii) + vect(numvect*basisvec_size+1:basisvec_size*(numvect+1));
        basisvec{ii} = reshape(tmp_vec,[nstixperside nstixperside 3]);
    end
    bkgnd_monitor = reshape(vect(numvect*basisvec_size+1:basisvec_size*(numvect+1)),[nstixperside nstixperside 3]);
    tmp1 = unique(basisvec{1}-bkgnd_monitor,'stable');
    tmp2 = unique(basisvec{2}-bkgnd_monitor,'stable');

    % plotting the rasterplots for neurothresh trials
    norms = cell(1,numel(num_targetspikerates));
    completed_search_alongdir = cell(1,numel(num_targetspikerates));
    for jj = 1: numel(num_targetspikerates)
        idxs = find(~isnan(stro.trial(:,correctidx)) & stro.trial(:,targetspikerateidx)==num_targetspikerates(jj));
        idxs(idxs<=neurothresh_startidx) = [];
        different_weights = unique(stro.trial(idxs,basisvecdiridx));
        tmp_norm = [];
        tmp_completed_search_alongdir = [];
        
        for kk = 1:numel(different_weights)
            idxs1 = find(stro.trial(:,basisvecdiridx) == different_weights(kk));
            idxs1(idxs1<neurothresh_startidx) = [];
            raster_data = stro.ras(idxs1,1);
            tmp_norm = [tmp_norm; stro.ras{idxs1(end),weightsidx}'];
            for ii = 1:size(raster_data,1)
                tmp = raster_data{ii} ;
                spikes = tmp(tmp>stro.trial(idxs1(ii),stimonidx)+t_offset & tmp < stro.trial(idxs1(ii),stimoffidx));
                spikes = spikes - stro.trial(idxs1(ii),stimonidx)-t_offset;
            end
            [flag, gamutViolation] = Check_ifcompletely_probed(stro,basisvec,bkgnd_monitor,idxs1,tmp_norm(end,:));
            % flag = 0, incompletely probed
            % flag = 1, completely probed
            % gamutViolation = 1, out of gamut point
            tmp_completed_search_alongdir = [tmp_completed_search_alongdir; flag gamutViolation];
        end
        norms{jj} = tmp_norm;
        completed_search_alongdir{jj} = tmp_completed_search_alongdir;
    end
    
    
    % Plotting the staircase termination and the geometric mean of last 3 staircase
    % termination points
    completed_dir = completed_search_alongdir{1};
    probed_dirs = logical(completed_dir(:,1)==1); % only including the directions that have been completely probed
    oog_idx = completed_dir(:,1)==1 & completed_dir(:,2)==1; % probed and out of gamut
    not_oog_idx = completed_dir(:,1)==1 & completed_dir(:,2)==0;
    
    % Important for storing and comparing the staircase termination points
    % to the last 5 geometric mean points
    staircase_termination = [];
    geometricmean_last5points = [];
    
    figure(1);
    for ii = 1:size(norms{1},1)
        dir = ii;
        for kk = 1: numel(num_targetspikerates)
            tmp_n = [];
            tmp_wts = [];
            tmp_parentvertices = [];
            idxs1 = find(stro.trial(:,basisvecdiridx) == dir);
            if ~isempty(idxs1)
                for jj = 1:numel(idxs1)
                    tmp_parentvertices = [tmp_parentvertices; stro.ras{idxs1(jj),parentverticesidx}];
                    tmp_wts = [tmp_wts; stro.ras{idxs1(jj),weightsidx}'];
                    tmp_n = [tmp_n; norm(stro.ras{idxs1(jj),weightsidx})];
                end
                raster_data = stro.ras(idxs1,1);
                num_dur =[];
                firing_rate = [];
                if (~isempty(tmp_n) & not_oog_idx(ii)) 
                    figure(1), hold on;
                    subplot(3,2,2*zz-1); plot(sign(tmp_wts(end,1))*geomean(abs(tmp_wts(end-4:end,1))),sign(tmp_wts(end,2))*geomean(abs(tmp_wts(end-4:end,2))),'o','MarkerFacecolor',[0 0 0], 'MarkerEdgeColor',[1 1 1]); hold on;
                    subplot(3,2,2*zz); plot(tmp_wts(end,1),tmp_wts(end,2),'o','MarkerFacecolor',[0 0 0], 'MarkerEdgeColor',[1 1 1]); hold on;   
                  
                    
                    tmp_termination_pt = norm([tmp_wts(end,1), tmp_wts(end,2)]);
                    staircase_termination = [staircase_termination; tmp_termination_pt];
                    
                    tmp_geometricmean_last5points = norm([sign(tmp_wts(end,1))*geomean(abs(tmp_wts(end-4:end,1))), sign(tmp_wts(end,2))*geomean(abs(tmp_wts(end-4:end,2)))]);
                    geometricmean_last5points = [geometricmean_last5points; tmp_geometricmean_last5points];
                    
                elseif (~isempty(tmp_n) & oog_idx(ii))
                    figure(1), hold on;
                    subplot(3,2,2*zz-1); plot([0 tmp_wts(end,1)],[0 tmp_wts(end,2)],'k'); hold on;
                    subplot(3,2,2*zz); plot([0 tmp_wts(end,1)],[0 tmp_wts(end,2)],'k'); hold on;     
                end
            end
        end
    end
 
   
    xlabel('X'), ylabel('Y');
    if zz == 1
        title('Simple');
    elseif zz==2
        title('DO');
    else
        title('Unclassified')
    end
    
     
    if zz==1
        figure(1); hold on;
        subplot(3,2,2*zz-1); set(gca,'Xlim',[-1.6 1.6],'Ylim',[-1.6 1.6],'YTick',-1.6:0.8:1.6,'XTick',[-1.6:0.8:1.6]); set(gca,'Tickdir','out'); drawnow; axis square; hold off;
        subplot(3,2,2*zz); set(gca,'Xlim',[-1.6 1.6],'Ylim',[-1.6 1.6],'YTick',-1.6:0.8:1.6,'XTick',[-1.6:0.8:1.6]);set(gca,'Tickdir','out'); drawnow; axis square; hold off;
    elseif zz==2
        figure(1); hold on;
        subplot(3,2,2*zz-1); set(gca,'Xlim',[-2.0 2.0],'Ylim',[-2.0 2.0],'YTick',-2.0:1.0:2.0,'XTick',[-2.0:1.0:2.0]); set(gca,'Tickdir','out'); drawnow; axis square; hold off;
        subplot(3,2,2*zz); set(gca,'Xlim',[-2.0 2.0],'Ylim',[-2.0 2.0],'YTick',-2.0:1.0:2.0,'XTick',[-2.0:1.0:2.0]); set(gca,'Tickdir','out'); drawnow; axis square; hold off;
    elseif zz==3
        figure(1); hold on;
        subplot(3,2,2*zz-1); set(gca,'Xlim',[-0.2 0.2],'Ylim',[-0.2 0.2],'YTick',-0.2:0.1:0.2,'XTick',[-0.2:0.1:0.2]); set(gca,'Tickdir','out'); drawnow; axis square; hold off;
        subplot(3,2,2*zz); set(gca,'Xlim',[-0.2 0.2],'Ylim',[-0.2 0.2],'YTick',-0.2:0.1:0.2,'XTick',[-0.2:0.1:0.2]); set(gca,'Tickdir','out'); drawnow; axis square; hold off;
    end
    
    
    figure(2); hold on;
    subplot(3,1,zz); plot(staircase_termination, geometricmean_last5points, 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [1 1 1]); 
    axis square; set(gca, 'Tickdir', 'out')
    
    if zz==1
        figure(2); hold on;
        subplot(3,1,zz); set(gca,'Xlim',[0 3],'Ylim',[0 3], 'Tickdir', 'out'); 
        plot([0 3], [0 3], 'color', [0 0 0]); axis square; hold off;
    elseif zz==2
        figure(2); hold on;
        subplot(3,1,zz); set(gca,'Xlim',[0, 3],'Ylim',[0 3],'Tickdir','out'); 
        plot([0 3], [0 3], 'color', [0 0 0]); drawnow; axis square; hold off;
    elseif zz==3
        figure(2); hold on;
        subplot(3,1,zz); set(gca,'Xlim',[0 0.3],'Ylim',[0 0.3],'Tickdir','out'); 
        plot([0 0.3], [0 0.3], 'color', [0 0 0]); drawnow; axis square; hold off;
    end
    
    
    % Now, I am going to plot the basis vectors
    wts  = [1 0;1 1;0 1;-1 1;-1 0;-1 -1;0 -1;1 -1];
    subplotidxs = [1;2;3;4;6;7;8;9];
    for ii = 1:numel(subplotidxs)
        vec = wts(ii,1)*(basisvec{1}-bkgnd_monitor) + wts(ii,2)*(basisvec{2}-bkgnd_monitor);
        normfactor = 0.5/(max(abs(vec(:)))+0.01);
        vec = normfactor*vec + 0.5;
        figure(3),subplot(numel(indices),size(wts,1),size(wts,1)*(zz-1)+ii); image(vec); set(gca,'XTick',[],'YTick',[]); axis square;
    end
    
    if zz == 1
        title('Simple');
    elseif zz==2
        title('DO');
    else
        title('Unclassified')
    end
    axis square; hold off;
    
    plot_counter = plot_counter + 1;
end



%% 16. Computing the number of pixels within each subregion 

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


% Acquiring the total pixels in subregions 1 and 2
subregion_pixels = [];
for zz = 1:numel(filename)
    disp(zz)
    stro = nex2stro(findfile(char(filename(zz))));
    
    spikename = 'sig001a';
    maskidx = strcmp(stro.sum.rasterCells(1,:),'subunit_mask');
    mask = stro.ras{end,maskidx};
    
    subregion_pixels = [subregion_pixels; sum(mask==1) sum(mask==2)];

end


% STATS

pixels = subregion_pixels([LUMidx, DOidx, hardtoclassifyidx], :);
mean_pixels = mean(pixels(:));
SD_pixels = std(pixels(:));

