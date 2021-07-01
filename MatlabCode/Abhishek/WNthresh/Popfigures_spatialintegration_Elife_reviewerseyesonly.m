% A new script for reviewers eyes only for Elife rebuttal
% Author- Abhishek De, 6/27/21

close all;
clearvars; 

%% Reviewer figure: Consistency of Isoresponse contours: 
% Robustness of the isoresponse curves based on the last point
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


%% Reviewer analysis: Rebound OFF responses
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


% Loading the files 
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filename = fetch(conn,'SELECT filename FROM WNthresh');
NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
close(conn);
filename = filename(strcmp(string(NTmode),"subunit"));
NTmode = NTmode(strcmp(string(NTmode),"subunit"));
spikeidx_NT = spikeidx_NT(strcmp(string(NTmode),"subunit"));

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
    t_offset = stro.trial(end,latencyidx)/1000;
    
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
    baselineFR = [baselineFR; mean(baselinespikecounts)/0.15];
    OFF_FR = [OFF_FR; mean(spikecountsafter)/(0.15-t_offset)];
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



%% Reviewer analysis: Eye movement (median eye displacement) vs. Isoresponse NLIs


%% Reviewer analysis: Reliability of NLI

if ~exist('plot_counter')
    plot_counter = 1;
end

% Laoding data for cone weight calculation  
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

% Storing the within subunit NLI
Within_subunit_NLI = [];

% For storing the Isoresponse NLI
Isoresponse_NLI = [];

% For storing the white noise NLI
Whitenoise_NLI = [];

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
    
    % Isoresponse NLI
    Isoresponse_NLI = [Isoresponse_NLI; log10(median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii}))];
end


%% A control analysis: Some more analyses that Greg suggested (a continuation of Figure 3 but not for the paper)
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
% SpatiallyOpponent = sum(sign(S1RGB).*sign(S2RGB),1)<3;
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
subplot(221); plot(RSSEisoresp_lin_median(LUMidx),RSSEisoresp_quad_median(LUMidx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(RSSEisoresp_lin_median(DOidx),RSSEisoresp_quad_median(DOidx),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); plot([0.0001 10],[0.0001 10],'k');
axis square; set(gca,'Tickdir','out','Xlim',[0.0001 10],'Ylim',[0.0001 10],'YScale','log','XScale','log','XTick',[0.0001 0.001 0.01 0.1 1 10],'YTick',[0.0001 0.001 0.01 0.1 1 10]); 
xlabel('median Linear error'); ylabel('median Quadratic error'); legend ('LUM','DO'); title('Isoresponse'); hold off;

subplot(222); plot(RSSEisoresp_lin_median(hardtoclassifyidx_woPC1),RSSEisoresp_quad_median(hardtoclassifyidx_woPC1),'o','MarkerFaceColor',[1 0 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(RSSEisoresp_lin_median(hardtoclassifyidx_wPC1),RSSEisoresp_quad_median(hardtoclassifyidx_wPC1),'o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]); hold on;  
axis square; set(gca,'Tickdir','out','Xlim',[0.0001 10],'Ylim',[0.0001 10],'YScale','log','XScale','log','XTick',[0.0001 0.001 0.01 0.1 1 10],'YTick',[0.0001 0.001 0.01 0.1 1 10]); plot([0.0001 10],[0.0001 10],'k'); 
xlabel('median Linear error'); ylabel('median Quadratic error'); legend('NSNDO woPC1','NSNDO wPC1'); title('Isoresponse'); hold off;

subplot(223); histogram(log10(RSSEisoresp_medianofratios(LUMidx)),-1:0.1:3,'DisplayStyle','stairs','EdgeColor',[0 0 0],'Linewidth',2); hold on;
histogram(log10(RSSEisoresp_medianofratios(DOidx)),-1:0.1:3,'DisplayStyle','stairs','EdgeColor',[1 0 0],'Linewidth',2); hold on;
plot(log10(median(RSSEisoresp_medianofratios(LUMidx))),14,'v','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
plot(log10(median(RSSEisoresp_medianofratios(DOidx))),15,'v','MarkerSize',8,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-1 3],'XTick',-1.0:0.5:3.0,'Ylim',[0 15],'YTick',[0 5 10 15]); 
ylabel('Count'); title('Isoresponse'); xlabel('Isoresponse NLI'); legend ('LUM','DO'); axis square; hold off;

subplot(224); histogram(log10(RSSEisoresp_medianofratios(hardtoclassifyidx_woPC1)),-1:0.1:3,'DisplayStyle','stairs','EdgeColor',[1 0 0.5],'Linewidth',2); hold on;
histogram(log10(RSSEisoresp_medianofratios(hardtoclassifyidx_wPC1)),-1:0.1:3,'DisplayStyle','stairs','EdgeColor',[1 0.5 0],'Linewidth',2); hold on;
plot(log10(median(RSSEisoresp_medianofratios(hardtoclassifyidx_woPC1))),13,'v','MarkerSize',8,'MarkerFaceColor',[1.0 0 0.5],'MarkerEdgeColor',[1 1 1]); 
plot(log10(median(RSSEisoresp_medianofratios(hardtoclassifyidx_wPC1))),14,'v','MarkerSize',8,'MarkerFaceColor',[1.0 0.5 0],'MarkerEdgeColor',[1 1 1]); 
set(gca,'Tickdir','out','Xlim',[-1 3],'XTick',-1.0:0.5:3.0,'Ylim',[0 15],'YTick',[0 5 10 15]); 
ylabel('Count'); title('Isoresponse'); xlabel('Isoresponse NLI'); legend('NSNDO woPC1','NSNDO wPC1'); axis square; hold off;
plot_counter = plot_counter + 1;


% WN subunit spatial interaction data 
figure(plot_counter);
subplot(221); plot(Acrosssubunits_lin_median(LUMidx),Acrosssubunits_quad_median(LUMidx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(Acrosssubunits_lin_median(DOidx),Acrosssubunits_quad_median(DOidx),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); plot([0.08 1],[0.08 1],'k');
axis square; set(gca,'Tickdir','out','Xlim',[0.08 1],'Ylim',[0.08 1],'YScale','log','XScale','log','XTick',[0.08 0.1 1],'YTick',[0.08 0.1 1]); 
xlabel('median GLM error'); ylabel('median GQM error'); legend ('LUM','DO'); title('Isoprobability'); hold off;

subplot(222); plot(Acrosssubunits_lin_median(hardtoclassifyidx_woPC1),Acrosssubunits_quad_median(hardtoclassifyidx_woPC1),'o','MarkerFaceColor',[1 0 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(Acrosssubunits_lin_median(hardtoclassifyidx_wPC1),Acrosssubunits_quad_median(hardtoclassifyidx_wPC1),'o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]); hold on;  
axis square; set(gca,'Tickdir','out','Xlim',[0.08 1],'Ylim',[0.08 1],'YScale','log','XScale','log','XTick',[0.08 0.1 1],'YTick',[0.08 0.1 1]); plot([0.08 1],[0.08 1],'k'); 
xlabel('median GLM error'); ylabel('median GQM error'); legend('HTC woPC1','HTC wPC1'); title('Isoprobability'); hold off;

subplot(223); histogram(log10(Acrosssubunits_medianofdifferences(LUMidx)),-0.02:0.005:0.1,'DisplayStyle','stairs','EdgeColor',[0 0 0],'Linewidth',2); hold on;
histogram(log10(Acrosssubunits_medianofdifferences(DOidx)),-0.02:0.005:0.1,'DisplayStyle','stairs','EdgeColor',[1 0 0],'Linewidth',2); hold on;
plot(log10(median(Acrosssubunits_medianofdifferences(LUMidx))),14,'v','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
plot(log10(median(Acrosssubunits_medianofdifferences(DOidx))),15,'v','MarkerSize',8,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-0.02 0.1],'XTick',-0.02:0.02:0.1,'Ylim',[0 15],'YTick',[0 5 10 15]); 
ylabel('Count'); title('Isoprobability'); xlabel('Isoprobability NLI'); legend ('LUM','DO'); axis square; hold off;

subplot(224); histogram(log10(Acrosssubunits_medianofdifferences(hardtoclassifyidx_woPC1)),-0.02:0.005:0.1,'DisplayStyle','stairs','EdgeColor',[1 0 0.5],'Linewidth',2); hold on;
histogram(log10(Acrosssubunits_medianofdifferences(hardtoclassifyidx_wPC1)),-0.02:0.005:0.1,'DisplayStyle','stairs','EdgeColor',[1 0.5 0],'Linewidth',2); hold on;
plot(log10(median(Acrosssubunits_medianofdifferences(hardtoclassifyidx_woPC1))),13,'v','MarkerSize',8,'MarkerFaceColor',[1.0 0 0.5],'MarkerEdgeColor',[1 1 1]); 
plot(log10(median(Acrosssubunits_medianofdifferences(hardtoclassifyidx_wPC1))),14,'v','MarkerSize',8,'MarkerFaceColor',[1.0 0.5 0],'MarkerEdgeColor',[1 1 1]); 
set(gca,'Tickdir','out','Xlim',[-0.02 0.1],'XTick',-0.02:0.02:0.1,'Ylim',[0 15],'YTick',[0 5 10 15]); 
ylabel('Count'); title('Isoprobability'); xlabel('Isoprobability NLI'); legend('HTC woPC1','HTC wPC1'); axis square; hold off;

plot_counter = plot_counter + 1;

% Some more stats
[p1,~] = ranksum(Isoresponse_NLI([LUMidx DOidx]),Isoresponse_NLI(hardtoclassifyidx_woPC1));
[p2,~] = ranksum(Isoresponse_NLI([LUMidx DOidx]),Isoresponse_NLI(hardtoclassifyidx_wPC1));
[p3,~] = ranksum(Isoresponse_NLI([LUMidx DOidx]),Isoresponse_NLI([hardtoclassifyidx_woPC1 hardtoclassifyidx_wPC1]));
