% Author - Abhishek De, 2/18, modifying a bit on 10/19
% I have used certain parts of it for CoSYNE
close all; clearvars;

plot_counter = 1;
% Loading the files 
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filename = fetch(conn,'SELECT filename FROM WNthresh');
NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
close(conn);
filename = filename(strcmp(string(NTmode),"subunit"));
NTmode = NTmode(strcmp(string(NTmode),"subunit"));
spikeidx_NT = spikeidx_NT(strcmp(string(NTmode),"subunit"));

THETA_all = cell(2,numel(filename));
RHO_all = cell(2,numel(filename));
oog_idx_all = cell(2,numel(filename));
not_oog_idx_all = cell(2,numel(filename));
bkgndRGB = cell(1,numel(filename)); % for storing the background RGBs 
subunitbasisvec = cell(1,numel(filename));
baselineFRstats = cell(1,numel(filename));
TFR = zeros(2,numel(filename));
baselineexceedsTFR = zeros(1,numel(filename));
confidence_score = zeros(1,numel(filename)); % a measure of how often the baselineFR exceeded the TFR
RF_loc = [];
Pangu = [];
Maui = [];
anglebwvectorsRGB = [];
WTS = cell(2,numel(filename));
SpikeCounts = cell(2,numel(filename));
S1LMS = []; S2LMS = [];
S1RGB = []; S2RGB = [];
pixelswithinsubunit1 = [];
pixelswithinsubunit2 = [];
count = 0;
for aa= 1:numel(filename)
    global reversalflagidx stepsizescale stepsize nreversals
    fileofinterest = char(filename(aa,:));
    disp(fileofinterest);
    if strcmp(fileofinterest(1),'M') == 1
        Maui = [Maui; aa];
    else
        Pangu = [Pangu; aa];
    end
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
    pixelswithinsubunit1 = [pixelswithinsubunit1; sum(stro.ras{inds(end)-1,maskidx}==1)];
    pixelswithinsubunit2 = [pixelswithinsubunit2; sum(stro.ras{inds(end)-1,maskidx}==2)];
    neurothreshmode = stro.trial(:,neurothreshidx);
    basisvec_dropidx = inds(end);
    neurothresh_startidx = max(find(neurothreshmode==1,1),basisvec_dropidx+1); % There is possibility of pressing the Neurothresh button multiple times
    num_targetspikerates = unique(stro.trial(neurothresh_startidx:end,targetspikerateidx));
    t_offset = stro.trial(end,latencyidx)/1000;
    vect = stro.ras{basisvec_dropidx,basisvecidx};
    basisvec_size = nstixperside*nstixperside*3;
    numvect = (numel(vect)/basisvec_size)-1;
    basisvec = cell(1,numvect);
    for ii = 1:numvect
        tmp_vec = vect((ii-1)*basisvec_size+1:basisvec_size*ii) + vect(numvect*basisvec_size+1:basisvec_size*(numvect+1));
        basisvec{ii} = reshape(tmp_vec,[nstixperside nstixperside 3]);
    end
    bkgnd_monitor = reshape(vect(numvect*basisvec_size+1:basisvec_size*(numvect+1)),[nstixperside nstixperside 3]);
    % Calculating cone and gun weights
    tmp1 = unique(basisvec{1}-bkgnd_monitor,'stable');
    tmp2 = unique(basisvec{2}-bkgnd_monitor,'stable');
    S1RGB = [S1RGB tmp1(tmp1~=0)];
    S2RGB = [S2RGB tmp2(tmp2~=0)];
    tmpLMS1 = M*tmp1(tmp1~=0);
    tmpLMS1 = tmpLMS1/sum(abs(tmpLMS1));
    S1LMS = [S1LMS tmpLMS1]; % subunit 1 cone weights
    tmpLMS2 = M*tmp2(tmp2~=0);
    tmpLMS2 = tmpLMS2/sum(abs(tmpLMS2));
    S2LMS = [S2LMS tmpLMS2]; % subunit 2 cone weights

    norms = cell(1,numel(num_targetspikerates));
    completed_search_alongdir = cell(1,numel(num_targetspikerates));
    for jj = 1:numel(num_targetspikerates) % as of now just focusing on 1 target firing rates
        tmp_wts = [];
        spikecounts = [];
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
            tmp_wts = [tmp_wts;cell2mat(stro.ras(idxs1,weightsidx)')'];
            for ii = 1:size(raster_data,1)
                tmp = raster_data{ii} ;
                spikes = tmp(tmp>stro.trial(idxs1(ii),stimonidx)+t_offset & tmp < stro.trial(idxs1(ii),stimoffidx));
                spikes = spikes - stro.trial(idxs1(ii),stimonidx)-t_offset;
                spikecounts = [spikecounts; numel(spikes)];
            end
            [flag, gamutViolation] = Check_ifcompletely_probed(stro,basisvec,bkgnd_monitor,idxs1,tmp_norm(end,:));
            % flag = 0, incompletely probed
            % flag = 1, completely probed
            % gamutViolation = 1, out of gamut point
            tmp_completed_search_alongdir = [tmp_completed_search_alongdir; flag gamutViolation];
        end
        norms{jj} = tmp_norm;
        completed_search_alongdir{jj} = tmp_completed_search_alongdir;
        WTS{jj,aa} = tmp_wts;
        SpikeCounts{jj,aa} = spikecounts;
        
    end
    basisvec1 = basisvec{1}-bkgnd_monitor;
    basisvec2 = basisvec{2}-bkgnd_monitor;
    oog1 = min(abs((1-bkgnd_monitor(:))./basisvec1(:)));
    oog2 = min(abs((1-bkgnd_monitor(:))./basisvec2(:)));
    tmp = basisvec1 + basisvec2;
    subunitbasisvec{aa} = (0.5 * tmp./(max(abs(tmp(:))) + 0.01)) + 0.5;
    
    % Calculate baselineFR
    num_spikes =[];
    num_dur = [];
    for ii = 1:size(stro.ras,1)
        tmp = stro.ras{ii,1} ;
        if ~isnan((stro.trial(ii,stimonidx)- stro.trial(ii,fpacqidx)))
            spikes = tmp(tmp<stro.trial(ii,stimonidx) & tmp>stro.trial(ii,fpacqidx));
            num_spikes = [num_spikes; numel(spikes)];
            num_dur = [num_dur; (stro.trial(ii,stimonidx)- stro.trial(ii,fpacqidx))];
        end
    end
    baselineFRstats{aa} = num_spikes./num_dur;
    TFR(1,aa) = num_targetspikerates(1);
    if numel(num_targetspikerates)>1
        TFR(2,aa) = num_targetspikerates(2);
    end
    baselineexceedsTFR(aa) = any(baselineFRstats{aa}>TFR(aa));
    confidence_score(aa) = 1-(numel(find(baselineFRstats{aa}>TFR(aa)))/numel(baselineFRstats{aa}));
        
    for ii = 1:numel(num_targetspikerates)
        tmp = norms{ii};
        completed_dir = completed_search_alongdir{ii};
        probed_dirs = logical(completed_dir(:,1)==1); % only including the directions that have been completely probed
        oog_idx = find(completed_dir(:,1)==1 & completed_dir(:,2)==1); % probed and out of gamut
        not_oog_idx = find(completed_dir(:,1)==1 & completed_dir(:,2)==0);
        fact = 0.5./sqrt(tmp(probed_dirs,1).^2 + tmp(probed_dirs,2).^2); % factor needed to extract unit vector
        [THETA1,RHO1] = cart2pol(tmp(:,1),tmp(:,2));
        ind = (1:numel(THETA1))';
        r = fliplr(linspace(0,1,numel(ind)));
        b = fliplr(r);
        THETA1 = THETA1 * (180/pi);
        THETA_all{ii,aa} = THETA1;
        RHO_all{ii,aa} = RHO1;
        oog_idx_all{ii,aa} = oog_idx;
        not_oog_idx_all{ii,aa} = not_oog_idx;
    end
%     if aa == numel(filename_c)
%         plot_counter = plot_counter + 1;
%         count = 0;
%     end
    count = count + 1;
    
    % Acquiring the receptive field location
    RF_loc = [RF_loc; stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y]; 
    
    % Calculating angle between vectors
    vec3 = S1RGB(:,aa);
    vec4 = S2RGB(:,aa);
    anglebwvectorsRGB = [anglebwvectorsRGB; 180*acos(dot(vec3,vec4)/(norm(vec3)*norm(vec4)))/pi];
end

savevariables = 0;
if savevariables
    save S1LMS S1LMS
    save S2LMS S2LMS
    save S1RGB S1RGB
    save S2RGB S2RGB
    save pixelswithinsubunit1 pixelswithinsubunit1
    save pixelswithinsubunit2 pixelswithinsubunit2 
    save WTS WTS 
    save SpikeCounts SpikeCounts
    save TFR TFR
    save baselineFRstats baselineFRstats
    save RF_loc RF_loc
    save Pangu Pangu
    save Maui Maui
    save anglebwvectorsRGB anglebwvectorsRGB
end

%% Next task is to fit the points and see which curve/ line describes the data well
% I am using linefit_AD2.m to fit a line and quadfit_AD.m to fit a quadratic equation which is currently doing a good job.
% I am storing the model parameters in linear_modelparams.mat and quad_modelparams.mat

Ssignal = abs(S1LMS(3,:))+abs(S2LMS(3,:)); % absolute S cone signal from 2 subunits
Ssignal = Ssignal';
count = 1;
GAMUTEDGE = 10;
num_rows = 4;
linear_modelparams = [];
quad_modelparams = [];
SSE_linearmodel = []; % residuals from the linear model 
SSE_quadmodel = []; % residuals from the quadratic model
RSSE_linearmodel = []; % residuals from the linear model from robust regression 
RSSE_quadmodel = []; % residuals from the quadratic model from robust regression
hruns_linearmodel = []; 
hruns_quadmodel = [];
pruns_linearmodel = []; % p value from runs test: linear model
pruns_quadmodel = []; % p value from runs test: quadratic model
rho_spearman = []; % spearman correlation coefficient for x and y
p_spearman = []; % p value for spearmann correlation coefficient
Fprob = [];
plotmode = 0; % 0-don't plot, 1-plot
ratioeig = [];
numInFs = [];
garbageval = 5;
Monkeyidentity = zeros(numel(filename),1);
Monkeyidentity(Maui) = 1;
Monkeyidentity(Pangu) = 2;
linear_FRsurfacefitparams = [];
quad_FRsurfacefitparams = [];
Rlin = [];
Rquad = [];
adjR2lin = [];
adjR2quad = [];
colormap('gray');
for ii = 1:numel(filename)
    ind = ii;
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
    [~,eigval] = eig(cov([x_orig(not_oog_idx) y_orig(not_oog_idx)]));
    ratioeig = [ratioeig; max(diag(eigval))/sum(abs(diag(eigval)))];
    
    % Fitting the linear model
    initguess1 = [x_orig(not_oog_idx) y_orig(not_oog_idx)]\ones(numel(x_orig(not_oog_idx)),1);
    [final_model1] = linefit_AD2(RHO, THETA,not_oog_idx,outofgamut,initguess1'); 
    rho1 = 1./(final_model1*[cos(allthetas); sin(allthetas)]);
    LOOGtmp1= rho1<0;
    [x_lin,y_lin] = pol2cart(allthetas(~LOOGtmp1),rho1(~LOOGtmp1));
    [h_runlin,p_runlin,nInF,fval1R,fval1OLS,tmpRlin,tmpadjR2lin] = calclinSSE(final_model1,RHO,THETA,not_oog_idx);
    
    initguess3 = [0 0 0 final_model1];
    [final_model3] = quadfit_AD2(RHO, THETA, not_oog_idx,outofgamut,initguess3);
    [final_model3] = conicsectionfit(x_orig(not_oog_idx),y_orig(not_oog_idx),initguess3); % passing on the values to conicsection fit which fits in x-y plane
    [final_model3] = quadfit_AD2(RHO, THETA, not_oog_idx,outofgamut,final_model3); % repassing the values to quadfit which fits in r-theta plane
    [x_quad,y_quad,rho3] = calc_xyvalues(allthetas, final_model3);
    L = rho3>0 & rho3==real(rho3);
    [x_quad2,y_quad2] = pol2cart(newtheta(L),rho3(L)');
    [h_runquad,p_runquad,fval3R,fval3OLS,tmpRquad,tmpadjR2quad] = calcquadSSE(final_model3,RHO,THETA,not_oog_idx);
            
    % Storing the linear and quadratic model parameters
    linear_modelparams = [linear_modelparams; final_model1];
    quad_modelparams = [quad_modelparams; final_model3];
    RSSE_linearmodel = [RSSE_linearmodel; fval1R];
    RSSE_quadmodel = [RSSE_quadmodel; fval3R];
    SSE_linearmodel = [SSE_linearmodel; fval1OLS];
    SSE_quadmodel = [SSE_quadmodel; fval3OLS];
    adjR2lin = [adjR2lin; tmpadjR2lin];
    adjR2quad = [adjR2quad; tmpadjR2quad];
    
    % Calculating F_statistic and its associated p value 
    F_stat = ((fval1R - fval3R)/(3))/(fval1R/(numel(not_oog_idx)-5));
    Fprob = [Fprob; 1-fcdf(F_stat,3,numel(not_oog_idx)-5)];
    
    % Storing the Pearson's correlation coefficient from the linear and the non-linear models 
    Rlin = [Rlin; tmpRlin];
    Rquad = [Rquad; tmpRquad];
    
    [rho,p] = corr(x_orig(not_oog_idx),y_orig(not_oog_idx),'Type','Spearman');
    rho_spearman = [rho_spearman; rho];
    p_spearman = [p_spearman; p];
    numInFs = [numInFs; nInF];
   
    if p_runlin<0.05
        c = [1 0 1]; % magenta
        disp(filename(ind));
    else 
        c = [0 1 1]; % sky blue 
    end
    
    %Plotting the figures
    figure(plot_counter), subplot(num_rows,7,7*count-6),plot(x_lin,y_lin,'g','Linewidth',2); hold on; 
    plot(x_orig(not_oog_idx), y_orig(not_oog_idx),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'PickableParts','none','MarkerEdgeColor',[0 0 1]);
    set(gca,'XLim',[-4 4],'YLim',[-4 4]); drawnow; axis square; hold off;
    subplot(num_rows,7,7*count-5),plot(allthetas(~LOOGtmp1),log10(rho1(~LOOGtmp1)),'g','Linewidth',2); hold on; 
    plot(THETA(not_oog_idx), log10(RHO(not_oog_idx)),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'PickableParts','none','MarkerEdgeColor',[1 0 0]);
    axis equal; drawnow; title(filename{ind});hold off;
    subplot(num_rows,7,7*count-4), bar([S1LMS(:,ind) S2LMS(:,ind)]); set(gca,'XTick',[1 2 3],'XTickLabel',{'L','M','S'}); axis square;% plotting the cone weights 
    subplot(num_rows,7,7*count-3), bar([S1RGB(:,ind) S2RGB(:,ind)]); set(gca,'XTick',[1 2 3],'XTickLabel',{'R','G','B'}); axis square;% plotting the gun weights 
    subplot(num_rows,7,7*count-2), image(subunitbasisvec{ind}); set(gca,'XTick',[],'YTick',[]);axis square;
    subplot(num_rows,7,7*count-1), plot(x_quad2,y_quad2,'k','Linewidth',2); hold on; 
    plot(x_orig(not_oog_idx), y_orig(not_oog_idx),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'PickableParts','none','MarkerEdgeColor',[0 0 1]);
    plot(x_orig(oog_idx), y_orig(oog_idx),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'PickableParts','none','MarkerEdgeColor',[0 1 0]);
    set(gca,'XLim',[-4 4],'YLim',[-4 4]); axis square; drawnow;hold off;
    subplot(num_rows,7,7*count),plot(newtheta(L),log10(rho3(L)),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'PickableParts','none','MarkerEdgeColor',[0 0 0]); hold on; 
    plot(THETA(not_oog_idx), log10(RHO(not_oog_idx)),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',c,'PickableParts','none','MarkerEdgeColor',c);
    plot(THETA(oog_idx), log10(RHO(oog_idx)),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'PickableParts','none','MarkerEdgeColor',[0 1 0]); axis equal; 
    drawnow;

count  = count + 1;
if count == (num_rows + 1)
    count = 1;
    plot_counter = plot_counter + 2;
end
end

savevariables = 0;
if savevariables
    save linear_modelparams linear_modelparams
    save quad_modelparams quad_modelparams
    save SSE_linearmodel SSE_linearmodel 
    save SSE_quadmodel SSE_quadmodel
    save RSSE_linearmodel RSSE_linearmodel
    save RSSE_quadmodel RSSE_quadmodel 
    save ratioeig ratioeig
    save anglebwvectorsRGB anglebwvectorsRGB
end

%% Now I want to compare metrics across different cells types 
load newOCidx.mat
load newLMidx.mat
load newLUMidx.mat 
load newSOidx.mat
load newhardtoclassifyidx.mat
OCidx = newOCidx;
LMidx = newLMidx;
DOidx = [OCidx LMidx];
LUMidx = newLUMidx';
SOidx = newSOidx;
hardtoclassifyidx = [SOidx newhardtoclassifyidx];
idx = [DOidx LUMidx hardtoclassifyidx];

load numsubunitspikes.mat
Q_rob = log10(RSSE_linearmodel./RSSE_quadmodel);
RFsize = pixelswithinsubunit1+pixelswithinsubunit2;
ecc = sqrt(sum(RF_loc.^2,2))/10;
num_pts_oog = []; num_pts_not_oog = [];
for ii = 1:numel(filename)
    num_pts_oog = [num_pts_oog; numel(oog_idx_all{1,ii})];
    num_pts_not_oog = [num_pts_not_oog; numel(not_oog_idx_all{1,ii})];
end

group = [ones(numel(DOidx),1); 2*ones(numel(LUMidx),1); 3*ones(numel(hardtoclassifyidx),1)];

[r1,p1] = corr(Q_rob,RFsize,'type','Spearman');
[r2,p2] = corr(Q_rob,anglebwvectorsRGB,'type','Spearman');
[r3,p3] = corr(Q_rob,Ssignal,'type','Spearman');
[r4,p4] = corr(Q_rob,ecc,'type','Spearman');
[r5,p5] = corr(Q_rob,confidence_score','type','Spearman');
[r6,p6] = corr(Q_rob,numsubunitspikes,'type','Spearman');

figure(plot_counter); % First looking at if there is any correlation of any parameter with spatial non-linearity
subplot(321); plot(Q_rob,RFsize,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); axis square;
xlabel('Spatial non-linearity'); ylabel('RF size'); set(gca,'Tickdir','out','Xlim',[0 5]); title(strcat('p=',num2str(p1,2))); 
subplot(322); plot(Q_rob,anglebwvectorsRGB,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); axis square;
xlabel('Spatial non-linearity'); ylabel('angle bw subunits'); set(gca,'Tickdir','out','Xlim',[0 5]); title(strcat('p=',num2str(p2,2)));  
subplot(323); plot(Q_rob,Ssignal,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); axis square;
xlabel('Spatial non-linearity'); ylabel('S-cone signal'); set(gca,'Tickdir','out','Xlim',[0 5]); title(strcat('p=',num2str(p3,2)));
subplot(324); plot(Q_rob,ecc,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); axis square;
xlabel('Spatial non-linearity'); ylabel('Eccentricity'); set(gca,'Tickdir','out','Xlim',[0 5]); title(strcat('p=',num2str(p4,2)));
subplot(325); plot(Q_rob,confidence_score,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); axis square;
xlabel('Spatial non-linearity'); ylabel('Confidence score'); set(gca,'Tickdir','out','Xlim',[0 5]); title(strcat('p=',num2str(p5,2)));
subplot(326); plot(Q_rob,numsubunitspikes,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); axis square;
xlabel('Spatial non-linearity'); ylabel('num spikes'); set(gca,'Tickdir','out','Xlim',[0 5]); title(strcat('p=',num2str(p6,2)));
plot_counter = plot_counter + 1;

p11 = kruskalwallis(Q_rob(idx),group,'off'); % ANOVA on residuals from robust regression
p12 = kruskalwallis(Ssignal(idx),group,'off'); % ANOVA on absolute S signal
p13 = kruskalwallis(RFsize(idx),group,'off'); % ANOVA on RFsize
p14 = kruskalwallis(num_pts_not_oog(idx),group,'off'); % ANOVA on number of not out-of-gamut points
p15 = kruskalwallis(num_pts_oog(idx),group,'off'); % ANOVA on number of out-of-gamut` points
p16 = kruskalwallis(num_pts_oog(idx)+num_pts_not_oog(idx),group,'off'); % ANOVA on number of out-of-gamut` points
p17 = kruskalwallis(numsubunitspikes(idx),group,'off'); % ANOVA on number of spikes collected during the subunit whitenoise phase

figure(plot_counter); % Secondly looking whether there is any pattern within the subclasses of the cells
subplot(331); boxplot(Q_rob(idx),group); set(gca,'Tickdir','out','XTicklabel',{'DO','Lum','oth'}); ylabel('Spatial non-linearity'); 
title(strcat('p=',num2str(p11,2))); axis square;
subplot(332); boxplot(Ssignal(idx),group); set(gca,'Tickdir','out','XTicklabel',{'DO','Lum','oth'}); ylabel('S-cone signal'); 
title(strcat('p=',num2str(p12,2))); axis square;
subplot(333); boxplot(RFsize(idx),group); set(gca,'Tickdir','out','XTicklabel',{'DO','Lum','oth'}); ylabel('RF size'); 
title(strcat('p=',num2str(p13,2))); axis square;
subplot(334); boxplot(num_pts_not_oog(idx),group); set(gca,'Tickdir','out','XTicklabel',{'DO','Lum','oth'}); ylabel('Not oog pts'); 
title(strcat('p=',num2str(p14,2))); axis square;
subplot(335); boxplot(num_pts_oog(idx),group); set(gca,'Tickdir','out','XTicklabel',{'DO','Lum','oth'}); ylabel('OOG pts'); 
title(strcat('p=',num2str(p15,2))); axis square;
subplot(336); boxplot(num_pts_not_oog(idx)+num_pts_oog(idx),group); set(gca,'Tickdir','out','XTicklabel',{'DO','Lum','oth'}); ylabel('total pts probed'); 
title(strcat('p=',num2str(p16,2))); axis square;
subplot(337); boxplot(numsubunitspikes(idx),group); set(gca,'Tickdir','out','XTicklabel',{'DO','Lum','oth'}); ylabel('subunit spikes'); 
title(strcat('p=',num2str(p17,2))); axis square;
plot_counter = plot_counter + 1;

% Analyzing whether there are monkey specific differences
load Maui.mat
load Pangu.mat
idx2 = [Maui; Pangu];
group2 = [ones(numel(Maui),1); 2*ones(numel(Pangu),1)];
p21 = kruskalwallis(Q_rob(idx2),group2,'off'); % ANOVA on residuals from robust regression

load RFstructure_c.mat
load RFstructure_l.mat
RFstructure = [RFstructure_c; RFstructure_l];
[~,p22] = fishertest([sum(RFstructure(Maui)==1) sum(RFstructure(Pangu)==1); sum(RFstructure(Maui)==2) sum(RFstructure(Pangu)==2)]);
figure(plot_counter); % Secondly looking whether there is any pattern within the subclasses of the cells
boxplot(Q_rob(idx2),group2); set(gca,'Tickdir','out','XTicklabel',{'Maui','Pangu'}); ylabel('Spatial non-linearity'); 
title(strcat('p=',num2str(p21,2))); axis square;
