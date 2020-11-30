% Starting a new script to visualize the isoresponse contours and do further analysis
% Obsolete script
% Author - Abhishek De, 09/17
close all; clearvars;
plot_counter = 1;
load filename_l.mat
load filename_c.mat
% Figure 1 - color cells, Figure 2 - luminance cells, 
filename = [filename_c; filename_l];
count = 1;
THETA_all = cell(1,numel(filename));
RHO_all = cell(1,numel(filename));
oog_idx_all = cell(1,numel(filename));
not_oog_idx_all = cell(1,numel(filename));
S1RGB = [];
S2RGB = [];
S1LMS = []; % for storing cone weights
S2LMS = [];
S1xyY = cell(1,numel(filename)); % storing xy chromaticities of subunit 1 
S2xyY = cell(1,numel(filename)); % storing xy chromaticities of subunit 2 
bkgndRGB = cell(1,numel(filename)); % for storing the background RGBs 
coneexcitationssub1 = cell(1,numel(filename)); % for storing cone excitations of subunit1
coneexcitationssub2 = cell(1,numel(filename)); % for storing cone excitations of subunit2
subunitbasisvec = cell(1,numel(filename));
RF_loc = [];
Pangu = [];
Maui = [];
for aa= 1:numel(filename)
    global reversalflagidx stepsizescale stepsize nreversals
    fileofinterest = char(filename(aa,:));
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
    basisvec1 = basisvec{1}-bkgnd_monitor;
    basisvec2 = basisvec{2}-bkgnd_monitor;
    oog1 = min(abs((1-bkgnd_monitor(:))./basisvec1(:)));
    oog2 = min(abs((1-bkgnd_monitor(:))./basisvec2(:)));
    tmp = basisvec1 + basisvec2;
    subunitbasisvec{aa} = (0.5 * tmp./(max(abs(tmp(:))) + 0.01)) + 0.5;

    for ii = 1:1
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
    end
    if aa == numel(filename_c)
        plot_counter = plot_counter + 1;
        count = 0;
    end
    count = count + 1;
    THETA_all{aa} = THETA1;
    RHO_all{aa} = RHO1;
    oog_idx_all{aa} = oog_idx;
    not_oog_idx_all{aa} = not_oog_idx;
    
    % Write the code for calculating the CIE chromaticities of subunit 1 and subunit 2
    subunit1RGB = S1RGB(:,end);
    subunit2RGB = S2RGB(:,end);
    bkgndRGB{aa} = squeeze(bkgnd_monitor(1,1,:));
    contrast1  = tmp(not_oog_idx,1)';
    RGBtripletssub1 = repmat(subunit1RGB,[1 numel(contrast1)]).*repmat(contrast1,[3 1]) + repmat(bkgndRGB{aa},[1 numel(contrast1)]);
    S1xyY{aa} = XYZToxyY(SRGBPrimaryToXYZ(RGBtripletssub1));
    contrast2  = tmp(not_oog_idx,2)';
    RGBtripletssub2 = repmat(subunit2RGB,[1 numel(contrast2)]).*repmat(contrast2,[3 1]) + repmat(bkgndRGB{aa},[1 numel(contrast2)]);
    S2xyY{aa} = XYZToxyY(SRGBPrimaryToXYZ(RGBtripletssub2));
    
    coneexcitationssub1{aa} = inv(M')*RGBtripletssub1;
    coneexcitationssub2{aa} = inv(M')*RGBtripletssub2;
    
    % Acquiring the receptive field location
    RF_loc = [RF_loc; stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y]; 
   
end
plot_counter = plot_counter + 1;


% Next task is to fit the points and see which curve/ line describes the data well
% I am using linefit_AD2.m to fit a line and quadfit_AD.m to fit a quadratic equation which is currently doing a good job.
% I am storing the model parameters in linear_modelparams.mat and quad_modelparams.mat
count = 1;
GAMUTEDGE = 5;
num_rows = 6;
linear_modelparams = [];
quad_modelparams = [];
SSE_linearmodel = []; % residuals from the linear model 
SSE_quadmodel = []; % residuals from the quadratic model
F_statistic = [];
F_testresult = [];
plotmode = 1; % 0-don't plot, 1-plot
for ind = 1:numel(filename)
    THETA = THETA_all{ind};
    THETA = THETA * pi/180;
    if any(THETA>3*pi/4)
        allthetas = linspace(-pi,pi,100);
    else
        allthetas = linspace(-pi/4,3*pi/4,100);
    end
    RHO = RHO_all{ind};
    oog_idx = oog_idx_all{ind};
    not_oog_idx = not_oog_idx_all{ind};
    outofgamut = zeros(size(THETA));
    outofgamut(oog_idx) = 1;
    outofgamut = logical(outofgamut);
    [x_orig, y_orig] = pol2cart(THETA,RHO);
    
    % Fitting the linear model
    [final_model1,fval1] = linefit_AD2(RHO, THETA, outofgamut,[100 100]);
    rho1 = 1./(final_model1*[cos(allthetas); sin(allthetas)]);
    LOOGtmp1= rho1>GAMUTEDGE|rho1<0;
    [x_lin,y_lin] = pol2cart(allthetas(~LOOGtmp1),rho1(~LOOGtmp1));
    
    %Fitting the quadratic model
    [final_model2,fval2] = quadfit_AD(RHO, THETA, outofgamut,[100 100 100]);
    rho2 = 1./(final_model2*[cos(allthetas).^2; sin(allthetas).^2; cos(allthetas).*sin(allthetas)]);
    LOOGtmp2= rho2>GAMUTEDGE|rho2<0;
    [x_quad,y_quad] = pol2cart(allthetas(~LOOGtmp2),rho2(~LOOGtmp2));
    
    % Storing the linear and quadratic model parameters
    linear_modelparams = [linear_modelparams; final_model1];
    quad_modelparams = [quad_modelparams; final_model2];
    
    % Now I am calculating the SSE of the linear and quad model
    THETA_sub = THETA(not_oog_idx);
    RHO_sub = RHO(not_oog_idx);
    linmodel_pred = 1./(final_model1*[cos(THETA_sub) sin(THETA_sub)]');
    L_lin = linmodel_pred<0;
    quadmodel_pred = 1./(final_model2*[cos(THETA_sub).^2 sin(THETA_sub).^2 cos(THETA_sub).*sin(THETA_sub)]');
    L_quad = quadmodel_pred<0;
    SSE_linearmodel = [SSE_linearmodel; sum((log(linmodel_pred(~L_lin)) - log(RHO_sub(~L_lin)')).^2) + sum(L_lin)];
    SSE_quadmodel = [SSE_quadmodel; sum((log(quadmodel_pred(~L_quad)) - log(RHO_sub(~L_quad)')).^2) + sum(L_quad)];
    F_statistic = [F_statistic; (SSE_linearmodel(end)-SSE_quadmodel(end))/(SSE_quadmodel(end)/(numel(not_oog_idx)-3))];
    F_testresult = [F_testresult; F_statistic(end)>finv(0.99,1,numel(not_oog_idx)-3)]; % i.e p<0.01
    
    if plotmode
        %Plotting the figures
        figure(plot_counter), subplot(num_rows,6,6*count-5),plot(x_lin,y_lin,'g','Linewidth',2); hold on;
        plot(x_orig(not_oog_idx), y_orig(not_oog_idx),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'PickableParts','none','MarkerEdgeColor',[0 0 1]);
        axis equal; drawnow; hold off;
        subplot(num_rows,6,6*count-4),plot(allthetas(~LOOGtmp1),log10(rho1(~LOOGtmp1)),'g','Linewidth',2); hold on;
        plot(THETA(not_oog_idx), log10(RHO(not_oog_idx)),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'PickableParts','none','MarkerEdgeColor',[1 0 0]);
        axis equal; drawnow; hold off;
        subplot(num_rows,6,6*count-3), plot(x_quad,y_quad,'g','Linewidth',2); hold on;
        plot(x_orig(not_oog_idx), y_orig(not_oog_idx),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'PickableParts','none','MarkerEdgeColor',[0 0 1]);
        axis equal; drawnow;hold off;
        
        if F_testresult(end) == 1
            c = [1 0 1];
        else
            c = [0 1 1];
        end
        subplot(num_rows,6,6*count-2),plot(allthetas(~LOOGtmp2),log10(rho2(~LOOGtmp2)),'k','Linewidth',2); hold on;
        plot(THETA(not_oog_idx), log10(RHO(not_oog_idx)),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',c,'PickableParts','none','MarkerEdgeColor',c);
        axis equal; drawnow; hold off;
        subplot(num_rows,6,6*count-1), bar([S1LMS(:,ind) S2LMS(:,ind)]); set(gca,'XTick',[1 2 3],'XTickLabel',{'L','M','S'}); % plotting the cone weights
        subplot(num_rows,6,6*count), image(subunitbasisvec{ind}); set(gca,'XTick',[],'YTick',[]);
        count  = count + 1;
        if count == (num_rows + 1)
            count = 1;
            plot_counter = plot_counter + 1;
        end
    end
    
end
plot_counter = plot_counter + 1;
% save('linear_modelparams.mat','linear_modelparams');
% save('quad_modelparams.mat','quad_modelparams');
% save('SSE_linearmodel.mat','SSE_linearmodel');
% save('SSE_quadmodel.mat','SSE_quadmodel');
% save('F_statistic.mat','F_statistic');
% save('F_testresult.mat','F_testresult');

% Plotting the cone weights in 3D LMS space: This might not be the only way to look at the data therefore I am plotting the RGB weights
N = numel(filename_c);
figure(plot_counter),subplot(121);plot3(S1LMS(1,1:N),S1LMS(2,1:N),S1LMS(3,1:N),'o','LineWidth',1.0,'MarkerFaceColor',[1 0 0],'MarkerSize',6); hold on;
plot3(S2LMS(1,1:N),S2LMS(2,1:N),S2LMS(3,1:N),'o','LineWidth',1.0,'MarkerFaceColor',[1 0 0],'MarkerSize',6);
plot3(S1LMS(1,N+1:end),S1LMS(2,N+1:end),S1LMS(3,N+1:end),'o','LineWidth',1.0,'MarkerFaceColor',[0 0 1],'MarkerSize',6); 
plot3(S2LMS(1,N+1:end),S2LMS(2,N+1:end),S2LMS(3,N+1:end),'o','LineWidth',1.0,'MarkerFaceColor',[0 0 1],'MarkerSize',6);
plot3([-1 1],[0 0],[0 0],'-k','Linewidth',2); % x- axis 
plot3([0 0],[-1 1],[0 0],'-k','Linewidth',2); % y- axis
plot3([0 0],[0 0],[-1 1],'-k','Linewidth',2); % z- axis
xlabel('L'),ylabel('M'),zlabel('S'), title('Cone Weights'); hold off;
% plotting the RGB weights now 
subplot(122);plot3(S1RGB(1,1:N),S1RGB(2,1:N),S1RGB(3,1:N),'o','LineWidth',1.0,'MarkerFaceColor',[1 0 0],'MarkerSize',6); hold on;
plot3(S2RGB(1,1:N),S2RGB(2,1:N),S2RGB(3,1:N),'o','LineWidth',1.0,'MarkerFaceColor',[1 0 0],'MarkerSize',6);
plot3(S1RGB(1,N+1:end),S1RGB(2,N+1:end),S1RGB(3,N+1:end),'o','LineWidth',1.0,'MarkerFaceColor',[0 0 1],'MarkerSize',6); 
plot3(S2RGB(1,N+1:end),S2RGB(2,N+1:end),S2RGB(3,N+1:end),'o','LineWidth',1.0,'MarkerFaceColor',[0 0 1],'MarkerSize',6);
plot3([-1 1],[0 0],[0 0],'-k','Linewidth',2); % x- axis 
plot3([0 0],[-1 1],[0 0],'-k','Linewidth',2); % y- axis
plot3([0 0],[0 0],[-1 1],'-k','Linewidth',2); % z- axis
xlabel('R'),ylabel('G'),zlabel('B'), title('Gun Weights'); hold off;
plot_counter = plot_counter + 1;

% Doing some further analysis
N1 = numel(filename_c);
absL = abs(S1LMS(1,:)) + abs(S2LMS(1,:));
absM = abs(S1LMS(2,:)) + abs(S2LMS(2,:));
absS = abs(S1LMS(3,:)) + abs(S2LMS(3,:));
absR = abs(S1RGB(1,:)) + abs(S2RGB(1,:));
absG = abs(S1RGB(2,:)) + abs(S2RGB(2,:));
absB = abs(S1RGB(3,:)) + abs(S2RGB(3,:));
absLMdiff = abs(S1LMS(1,:)-S1LMS(2,:)-S2LMS(1,:)+S2LMS(2,:));
logresidualratio = log(SSE_linearmodel./SSE_quadmodel);
corrL = corrcoef(logresidualratio(1:N1),absL(1:N1)');
corrM = corrcoef(logresidualratio(1:N1),absM(1:N1)');
corrS = corrcoef(logresidualratio(1:N1),absS(1:N1)');
corrR = corrcoef(logresidualratio(1:N1),absR(1:N1)');
corrG = corrcoef(logresidualratio(1:N1),absG(1:N1)');
corrB = corrcoef(logresidualratio(1:N1),absB(1:N1)');
corrLMdiff = corrcoef(logresidualratio(1:N1),absLMdiff(1:N1)');
figure(plot_counter); subplot(331),hist(logresidualratio(1:N1),20); xlabel('Ratio(Lin/Quad)'); title('Ratio of residuals');
subplot(332),hist(F_statistic,20); xlabel('F statistic'); title('Dist of F statistic');
subplot(333),plot(logresidualratio(1:N1),absL((1:N1)),'o','LineWidth',1.0,'MarkerFaceColor',[1 0 0],'MarkerSize',6); hold on; lsline; ylabel('L'); xlabel('Ratio(Lin/Quad)'); 
text(-4.5,1.25,num2str(corrL(1,2))); set(gca,'Xlim',[-5 5]); hold off;
subplot(334),plot(logresidualratio(1:N1),absM((1:N1)),'o','LineWidth',1.0,'MarkerFaceColor',[0 1 0],'MarkerSize',6); hold on; lsline; ylabel('M'); xlabel('Ratio(Lin/Quad)');
text(-4.5,1.25,num2str(corrM(1,2))); set(gca,'Xlim',[-5 5]); hold off;
subplot(335),plot(logresidualratio(1:N1),absS((1:N1)),'o','LineWidth',1.0,'MarkerFaceColor',[0 0 1],'MarkerSize',6); hold on; lsline; ylabel('S'); xlabel('Ratio(Lin/Quad)'); 
text(-4.5,1.25,num2str(corrS(1,2))); set(gca,'Xlim',[-5 5]); hold off;
subplot(336),plot(logresidualratio(1:N1),absLMdiff((1:N1)),'o','LineWidth',1.0,'MarkerFaceColor',[0 1 1],'MarkerSize',6); hold on; lsline; ylabel('L-M'); xlabel('Ratio(Lin/Quad)'); 
text(-4.5,1.25,num2str(corrLMdiff(1,2))); set(gca,'Xlim',[-5 5]); hold off;
subplot(337),plot(logresidualratio(1:N1),absR((1:N1)),'o','LineWidth',1.0,'MarkerFaceColor',[1 0 0],'MarkerSize',6); hold on; lsline; ylabel('R'); xlabel('Ratio(Lin/Quad)'); 
text(-4.5,0.4,num2str(corrR(1,2))); set(gca,'Xlim',[-5 5]); hold off;
subplot(338),plot(logresidualratio(1:N1),absG((1:N1)),'o','LineWidth',1.0,'MarkerFaceColor',[0 1 0],'MarkerSize',6); hold on; lsline; ylabel('G'); xlabel('Ratio(Lin/Quad)');
text(-4.5,0.4,num2str(corrG(1,2))); set(gca,'Xlim',[-5 5]); hold off;
subplot(339),plot(logresidualratio(1:N1),absB((1:N1)),'o','LineWidth',1.0,'MarkerFaceColor',[0 0 1],'MarkerSize',6); hold on; lsline; ylabel('B'); xlabel('Ratio(Lin/Quad)'); 
text(-4.5,0.4,num2str(corrB(1,2))); set(gca,'Xlim',[-5 5]); hold off;
plot_counter = plot_counter + 1;

% Selecting those neurons where the residual ratio of linear model is higher than the quadratic model by a factor of 2
cellsofinterest = find(logresidualratio>2);
count = 1;
linear_modelparams = [];
quad_modelparams = [];
SSE_linearmodel = []; % residuals from the linear model 
SSE_quadmodel = []; % residuals from the quadratic model
F_statistic = [];
F_testresult = [];

for ii = 1: numel(cellsofinterest)
    ind = cellsofinterest(ii);
    THETA = THETA_all{ind};
    THETA = THETA * pi/180;
    if any(THETA>3*pi/4)
        allthetas = linspace(-pi,pi,100);
    else
        allthetas = linspace(-pi/4,3*pi/4,100);
    end
    RHO = RHO_all{ind};
    oog_idx = oog_idx_all{ind};
    not_oog_idx = not_oog_idx_all{ind};
    outofgamut = zeros(size(THETA));
    outofgamut(oog_idx) = 1;
    outofgamut = logical(outofgamut);
    [x_orig, y_orig] = pol2cart(THETA,RHO);
    
    % Fitting the linear model
    [final_model1,fval1] = linefit_AD2(RHO, THETA, outofgamut,[100 100]);
    rho1 = 1./(final_model1*[cos(allthetas); sin(allthetas)]);
    LOOGtmp1= rho1>GAMUTEDGE|rho1<0;
    [x_lin,y_lin] = pol2cart(allthetas(~LOOGtmp1),rho1(~LOOGtmp1));
    
    %Fitting the quadratic model
    [final_model2,fval2] = quadfit_AD(RHO, THETA, outofgamut,[100 100 100]);
    rho2 = 1./(final_model2*[cos(allthetas).^2; sin(allthetas).^2; cos(allthetas).*sin(allthetas)]);
    LOOGtmp2= rho2>GAMUTEDGE|rho2<0;
    [x_quad,y_quad] = pol2cart(allthetas(~LOOGtmp2),rho2(~LOOGtmp2));
    
    % Storing the linear and quadratic model parameters
    linear_modelparams = [linear_modelparams; final_model1];
    quad_modelparams = [quad_modelparams; final_model2];
    
    % Now I am calculating the SSE of the linear and quad model
    THETA_sub = THETA(not_oog_idx);
    RHO_sub = RHO(not_oog_idx);
    linmodel_pred = 1./(final_model1*[cos(THETA_sub) sin(THETA_sub)]');
    L_lin = linmodel_pred<0;
    quadmodel_pred = 1./(final_model2*[cos(THETA_sub).^2 sin(THETA_sub).^2 cos(THETA_sub).*sin(THETA_sub)]');
    L_quad = quadmodel_pred<0;
    SSE_linearmodel = [SSE_linearmodel; sum((log(linmodel_pred(~L_lin)) - log(RHO_sub(~L_lin)')).^2) + sum(L_lin)];
    SSE_quadmodel = [SSE_quadmodel; sum((log(quadmodel_pred(~L_quad)) - log(RHO_sub(~L_quad)')).^2) + sum(L_quad)];
    F_statistic = [F_statistic; (SSE_linearmodel(end)-SSE_quadmodel(end))/(SSE_quadmodel(end)/(numel(not_oog_idx)-3))];
    F_testresult = [F_testresult; F_statistic(end)>finv(0.99,1,numel(not_oog_idx)-3)]; % i.e p<0.01
    
    %Plotting the figures
    figure(plot_counter), subplot(num_rows,7,7*count-6),plot(x_lin,y_lin,'g','Linewidth',2); hold on; 
    plot(x_orig(not_oog_idx), y_orig(not_oog_idx),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'PickableParts','none','MarkerEdgeColor',[0 0 1]);
    axis equal; drawnow; hold off;
    subplot(num_rows,7,7*count-5),plot(allthetas(~LOOGtmp1),log10(rho1(~LOOGtmp1)),'g','Linewidth',2); hold on; 
    plot(THETA(not_oog_idx), log10(RHO(not_oog_idx)),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'PickableParts','none','MarkerEdgeColor',[1 0 0]);
    axis equal; drawnow; hold off;
    subplot(num_rows,7,7*count-4), plot(x_quad,y_quad,'g','Linewidth',2); hold on; 
    plot(x_orig(not_oog_idx), y_orig(not_oog_idx),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'PickableParts','none','MarkerEdgeColor',[0 0 1]);
    axis equal; drawnow;hold off;
        
    if F_testresult(end) == 1
        c = [1 0 1];
    else 
        c = [0 1 1];
    end
    subplot(num_rows,7,7*count-3),plot(allthetas(~LOOGtmp2),log10(rho2(~LOOGtmp2)),'k','Linewidth',2); hold on; 
    plot(THETA(not_oog_idx), log10(RHO(not_oog_idx)),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',c,'PickableParts','none','MarkerEdgeColor',c);
    axis equal; drawnow; hold off;
    subplot(num_rows,7,7*count-2), bar([S1LMS(:,ind) S2LMS(:,ind)]); set(gca,'XTick',[1 2 3],'XTickLabel',{'L','M','S'}); % plotting the cone weights 
    subplot(num_rows,7,7*count-1), image(subunitbasisvec{ind}); set(gca,'XTick',[],'YTick',[]);
    subplot(num_rows,7,7*count), bar([S1RGB(:,ind)/sum(abs(S1RGB(:,ind))) S2RGB(:,ind)/sum(abs(S2RGB(:,ind)))]); set(gca,'XTick',[1 2 3],'XTickLabel',{'R','G','B'}); % plotting the gun weights
    count  = count + 1;
    if count == (num_rows + 1)
        count = 1;
        plot_counter = plot_counter + 1;
    end
end
plot_counter = plot_counter + 1;

% I think the next logical step is to calculate the SSE (sum of squared errors and do a hypothesis testing using F test)
% Plotting the difference in the chromaticities between the subunits in CIE space and some additional analysis
load T_xyz1964.mat;
T_xyY = T_xyz1964./(repmat(sum(T_xyz1964),[3 1])+0.0001);
rgb = reshape([T_xyY(1,1:end-20)',T_xyY(2,1:end-20)',T_xyY(3,1:end-20)'],[size(T_xyY(:,1:end-20),2) 1 3]);
figure(plot_counter),patch(T_xyY(1,1:end-20),T_xyY(2,1:end-20),rgb); hold on; plot(0.33,0.33,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
xlabel('x'), ylabel('y'); title('CIE');hold off;
plot_counter = plot_counter + 1;
count = 1;
num_rows = 6;
anglebwvectors = [];
anglebwvectorsRGB = [];
chrom2lumvar = [];
plotmode = 0; % 0-don't plot, 1-plot
for ii = 1:numel(filename)
    diffxyY = S1xyY{ii} - S2xyY{ii};
    k1 = S1xyY{ii};
    k2 = S2xyY{ii};
    bkgndxyY = XYZToxyY(SRGBPrimaryToXYZ(bkgndRGB{ii}));
    basisvec1xyY = XYZToxyY(SRGBPrimaryToXYZ(S1RGB(:,ii)+bkgndRGB{ii}));
    basisvec2xyY = XYZToxyY(SRGBPrimaryToXYZ(S2RGB(:,ii)+bkgndRGB{ii}));
    
    THETA = THETA_all{ii}*pi/180;
    RHO = RHO_all{ii};
    not_oog_idx = not_oog_idx_all{ii};
    [x_orig, y_orig] = pol2cart(THETA,RHO);
    
    LMSsub1 = coneexcitationssub1{ii};
    LMSsub2 = coneexcitationssub2{ii};
    vec1 = basisvec1xyY - bkgndxyY;
    vec2 = basisvec2xyY - bkgndxyY;
    vec3 = S1RGB(:,ii);
    vec4 = S2RGB(:,ii);
    anglebwvectors = [anglebwvectors; 180*acos(dot(vec1,vec2)/(norm(vec1)*norm(vec2)))/pi];
    anglebwvectorsRGB = [anglebwvectorsRGB; 180*acos(dot(vec3,vec4)/(norm(vec3)*norm(vec4)))/pi];
    
    % calculating chromatic to luminance variance
    chrom2lumvar = [chrom2lumvar; max(var([k1(1,:) k2(1,:)]),var([k1(2,:) k2(2,:)]))/var([k1(3,:) k2(3,:)])];     
    
    if plotmode
        % plotting the xy values in CIE coordinates for both the subunits value
        figure(plot_counter); subplot(num_rows,6,6*count-5), plot(k1(1,:), k1(2,:),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1]); hold on;
        plot(k2(1,:), k2(2,:),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 1 0]); hold on;
        %     plot(T_xyY(1,:), T_xyY(2,:),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0]);
        plot(bkgndxyY(1,:), bkgndxyY(2,:),'*'); drawnow; hold off;
        
        % want to check if the chromaticities of the 2 subunits are opposite to each other
        figure(plot_counter), subplot(num_rows,6,6*count-4), plot([bkgndxyY(1) basisvec1xyY(1)],[bkgndxyY(2) basisvec1xyY(2)],'r','Linewidth',2); hold on;
        plot([bkgndxyY(1) basisvec2xyY(1)],[bkgndxyY(2) basisvec2xyY(2)],'g','Linewidth',2); drawnow; hold off;
        
        % plotting the delta_x and delta_y of the chromaticities of the 2 subunits
        figure(plot_counter); subplot(num_rows,6,6*count-3), plot(diffxyY(1,:), diffxyY(2,:),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 1 1]);
        drawnow;
        
        % plotting the isoresponse contour
        figure(plot_counter); subplot(num_rows,6,6*count-2), plot(x_orig(not_oog_idx), y_orig(not_oog_idx),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0]);
        drawnow;
        
        % plotting the subunit STA
        subplot(num_rows,6,6*count-1), image(subunitbasisvec{ii}); set(gca,'XTick',[],'YTick',[]); drawnow;
        
        % plotting the ratio of cone excitations between the subunits
        subplot(num_rows,6,6*count),plot(sum((k1(1:2,:)-bkgndxyY(1:2)).*(k2(1:2,:)-bkgndxyY(1:2)),1),RHO(not_oog_idx),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1]);
        
        count  = count + 1;
        if count == num_rows + 1
            plot_counter = plot_counter + 1;
            count = 1;
        end
    end
    
    
end
plot_counter = plot_counter + 1;

figure(plot_counter), subplot(321),hist(anglebwvectors(1:N1)); xlabel('Angle between subunits CIE'); ylabel('Frequency'); title('DO neurons'); hold on;
plot(mean(anglebwvectors(1:N1)),0,'kv','MarkerFacecolor','g'); hold off;
subplot(322),hist(anglebwvectors(N1+1:end)); xlabel('Angle between subunits CIE'); ylabel('Frequency'); title('Luminance neurons'); hold on;
plot(mean(anglebwvectors(N1+1:end)),0,'kv','MarkerFacecolor','g'); hold off;
subplot(323),hist(anglebwvectorsRGB(1:N1)); xlabel('Angle between subunits RGB'); ylabel('Frequency'); title('DO neurons'); hold on;
plot(mean(anglebwvectorsRGB(1:N1)),0,'kv','MarkerFacecolor','g'); hold off;
subplot(324),hist(anglebwvectorsRGB(N1+1:end)); xlabel('Angle between subunits RGB'); ylabel('Frequency'); title('Luminance neurons'); hold on;
plot(mean(anglebwvectorsRGB(N1+1:end)),0,'kv','MarkerFacecolor','g'); hold off;
subplot(325),hist(log10(chrom2lumvar(1:N1))); xlabel('Ratio(Chrom\Lum)'); ylabel('Frequency'); title('DO neurons'); hold on;
plot(mean(log10(chrom2lumvar(1:N1))),0,'kv','MarkerFacecolor','g'); hold off;
subplot(326),hist(log10(chrom2lumvar(N1+1:end))); xlabel('Ratio(Chrom\Lum)'); ylabel('Frequency'); title('Luminance neurons'); hold on;
plot(mean(log10(chrom2lumvar(N1+1:end))),0,'kv','MarkerFacecolor','g'); hold off;
plot_counter = plot_counter + 1;

corrangleresidualratio = corrcoef(logresidualratio(1:N1),anglebwvectorsRGB(1:N1)');
figure(plot_counter); subplot(121),plot(logresidualratio(1:N1),anglebwvectorsRGB(1:N1),'o','LineWidth',1.0,'MarkerFaceColor',[1 0 0],'MarkerSize',6); hold on; lsline;;
xlabel('Ratio(Lin/Quad)'); ylabel('Angle between subunits RGB'); text(-0.5,20,num2str(corrangleresidualratio(1,2))); hold off; 
subplot(122), plot(RF_loc(Pangu,1), RF_loc(Pangu,2),'o','LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerSize',5); hold on;
plot(RF_loc(Maui,1), RF_loc(Maui,2),'o','LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerSize',5); 
xlabel('X'), ylabel('Y'); title('RF locations'); hold off;
plot_counter = plot_counter + 1;

% save logresidualratio logresidualratio
% save absS absS
% save anglebwvectorsRGB anglebwvectorsRGB 

%%

load latencySTAcheck.mat
load latencySTAsubunit.mat
load latencybetweensubunits.
load latencydiffWNchecksubunit.mat
load anglediffENchecksubunit.mat
% Here I am using a function suggested by Greg : prettypaircomp to see what variables have a strong correlation coefficient among themlves
% For color cells 
prettycorr([logresidualratio(1:N1),anglebwvectorsRGB(1:N1),anglebwvectors(1:N1),absL(1:N1)',absM(1:N1)',absS(1:N1)',absR(1:N1)',absG(1:N1)',absB(1:N1)'],{'LogResRat','AngleRGB','AngleCIE','L','M','S','R','G','B'});
prettycorr([logresidualratio(1:N1),anglediffENchecksubunit(1:N1),absS(1:N1)',latencybetweensubunits(1:N1),latencySTAsubunit(1:N1),latencySTAcheck(1:N1),latencydiffWNchecksubunit(1:N1),anglebwvectorsRGB(1:N1)],{'LogResRat','Anglediffchecksubunit','S','Latencybwsubunits','LatencySTAsubunit','LatencySTAcheck','Latencydiffchecksubunits','AngleRGB'});

% For luminance cells
prettycorr([logresidualratio(N1+1:end),anglebwvectorsRGB(N1+1:end),anglebwvectors(N1+1:end),absL(N1+1:end)',absM(N1+1:end)',absS(N1+1:end)',absR(N1+1:end)',absG(N1+1:end)',absB(N1+1:end)'],{'LogResRat','AngleRGB','AngleCIE','L','M','S','R','G','B'});
prettycorr([logresidualratio(N1+1:end),anglediffENchecksubunit(N1+1:end),absS(N1+1:end)',latencybetweensubunits(N1+1:end),latencySTAsubunit(N1+1:end),latencySTAcheck(N1+1:end),latencydiffWNchecksubunit(N1+1:end),anglebwvectorsRGB(N1+1:end)],{'LogResRat','Anglediffchecksubunit','S','Latencybwsubunits','LatencySTAsubunit','LatencySTAcheck','Latencydiffchecksubunits','AngleRGB'});


%%
%******************************************************************************************
%    Below this section are the preliminary codes/scripts for fitting a
%    linear/quadratic equation
%******************************************************************************************

% Am trying grid search to hunt for the best set of parameters that could minimize the fitting error, Currently searching the paramater in 2-D polar space. 
% Does a pretty good job of fitting the data. (Grid Search in r-theta space)
count = 1;
model_theta = -pi/4:pi/36:3*pi/4;
model_r = logspace(-1,0,20);
GAMUTEDGE = 5;
for ind = 1:numel(filename)
    data = zeros(length(model_theta),length(model_r));
    thetas = THETA_all{ind};
    thetas = thetas*pi/180;
    RHO = RHO_all{ind};
    not_oog_idx = not_oog_idx_all{ind};
    for i = 1:numel(model_theta)
        for j = 1:length(model_r)
            pred_staircase_terminations = model_r(j)./(cos(thetas(not_oog_idx)-model_theta(i)));
            L = pred_staircase_terminations > GAMUTEDGE | pred_staircase_terminations < 0;
            pred_staircase_terminations(L) = GAMUTEDGE;
            % Ignoring OOG points in the calculation of error for now.
            if sum(L)./length(L) > 0.7 % if > 0.5 of the seaches go out of gamut, this is probably not a good fit (e.g. weights are too small)
                err = nan;
            else
                err = mean((log(pred_staircase_terminations)-log(RHO(not_oog_idx))).^2);
            end
            data(i,j) = err;
        end
    end
    [tmp_i, tmp_j] = ind2sub([length(model_theta) length(model_r)],find(data==min(min(data))));
    bestparams = [model_theta(tmp_i) model_r(tmp_j)];
    preds = bestparams(2)./(cos(model_theta-bestparams(1)));
    LOOGtmp= preds>GAMUTEDGE|preds<0;
    [x1,y1] = pol2cart(model_theta(~LOOGtmp)', preds(~LOOGtmp)');
    [x2,y2] = pol2cart(thetas, RHO);
    figure(plot_counter); subplot(5,6,count); 
    plot(x1,y1,'g','Linewidth',2); hold on; 
    plot(x2(not_oog_idx), y2(not_oog_idx),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'PickableParts','none','MarkerEdgeColor',[0 0 1]);
    drawnow; axis equal; hold off;
    % plotting in log R - theta domain
    figure(plot_counter+2); subplot(5,6,count); 
    plot(model_theta(~LOOGtmp),log10(preds(~LOOGtmp)),'g','Linewidth',2); hold on; 
    plot(thetas(not_oog_idx), log10(RHO(not_oog_idx)),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'PickableParts','none','MarkerEdgeColor',[1 0 0]);
    drawnow; axis equal; hold off;
    count  = count + 1;
    if count == 31
        count = 1;
        plot_counter = plot_counter + 1;
    end
end
plot_counter = plot_counter + 2;

%%  This part of the code is inspired by GregAnalysis. Currently I am just fitting a line. Line works pretty well. (Grid search in x-y space)
count = 1;
A = linspace(-10,10,100);
B = A;
GAMUTEDGE = 5;
SSE_linearmodel = zeros(numel(filename),1);
for ind = 1:numel(filename)
    data = zeros(length(A),length(B));
    thetas = THETA_all{ind};
    thetas = thetas*pi/180;
    RHO = RHO_all{ind};
    not_oog_idx = not_oog_idx_all{ind};
    if any(thetas>3*pi/4)
        allthetas = linspace(-pi,pi,100);
    else
        allthetas = linspace(-pi/4,3*pi/4,100);
    end
    for i = 1:numel(A)
        for j = 1:numel(B)
            modelparams = [A(i) B(j)];
            pred_staircase_terminations = 1./(modelparams*[cos(thetas(not_oog_idx)) sin(thetas(not_oog_idx))]');
            L = pred_staircase_terminations > GAMUTEDGE | pred_staircase_terminations < 0;
            pred_staircase_terminations(L) = GAMUTEDGE; % A penalty is when prediction falls short of OOG
            % Ignoring OOG points in the calculation of error for now.
            if sum(L)/length(L) > 0.7 % if > 0.5 of the seaches go out of gamut, this is probably not a good fit (e.g. weights are too small)
                err = nan;
            else
                err = mean((log(pred_staircase_terminations)-log(RHO(not_oog_idx)')).^2);
            end
            data(i,j) = err;
        end
    end
    [tmp_i, tmp_j] = ind2sub([length(A) length(B)],find(data==min(min(data))));
    bestparams = [A(tmp_i) B(tmp_j)];
    preds = 1./(bestparams*[cos(allthetas); sin(allthetas)]);
    LOOGtmp= preds>GAMUTEDGE|preds<0;
    [x1,y1] = pol2cart(allthetas(~LOOGtmp)', preds(~LOOGtmp)');
    [x2,y2] = pol2cart(thetas, RHO);
    figure(plot_counter); subplot(5,6,count); 
    plot(x1,y1,'g','Linewidth',2); hold on; 
    plot(x2(not_oog_idx), y2(not_oog_idx),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'PickableParts','none','MarkerEdgeColor',[0 0 1]);
    drawnow; axis equal; hold off;
    % plotting in log R - theta domain
    figure(plot_counter+2); subplot(5,6,count); 
    plot(allthetas(~LOOGtmp),log10(preds(~LOOGtmp)),'g','Linewidth',2); hold on; 
    plot(thetas(not_oog_idx), log10(RHO(not_oog_idx)),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'PickableParts','none','MarkerEdgeColor',[1 0 0]);
    drawnow; axis equal; hold off;
    count  = count + 1;
    if count == 31
        count = 1;
        plot_counter = plot_counter + 1;
    end
    
   % Calculating the sum of squared errors for each neuron for a linear model
   preds_pts = 1./(bestparams*[cos(thetas) sin(thetas)]');
   L1 = preds_pts>GAMUTEDGE | preds_pts<0;
   preds_pts(L1) = GAMUTEDGE;
   SSE_linearmodel(ind) = sum((log(preds_pts)- log(RHO')).^2);
end
plot_counter = plot_counter + 2;

%% Now I am gonna implement fitting a conic section in 2-D plane (Grid search in x,y,z space)
count = 1;
A = logspace(0,3,40);
B = A;
C = A;
GAMUTEDGE = 5;
SSE_quadmodel = zeros(numel(filename),1);
for ind = 1:numel(filename)
    data = zeros(length(A),length(B),length(C));
    thetas = THETA_all{ind};
    thetas = thetas*pi/180;
    RHO = RHO_all{ind};
    if any(thetas>3*pi/4)
        allthetas = linspace(-pi,pi,100);
    else
        allthetas = linspace(-pi/4,3*pi/4,100);
    end
    not_oog_idx = not_oog_idx_all{ind};
    for i = 1:numel(A)
        for j = 1:numel(B)
            for k = 1:numel(C)
                modelparams = [A(i) B(j) C(k)];
                pred_staircase_terminations = 1./(modelparams*[cos(thetas(not_oog_idx)).^2 sin(thetas(not_oog_idx)).^2 cos(thetas(not_oog_idx)).*sin(thetas(not_oog_idx))]');
                if any(pred_staircase_terminations<0)
                    data(i,j,k) = nan;
                else
                    pred_staircase_terminations = sqrt(pred_staircase_terminations);
                    L = pred_staircase_terminations > GAMUTEDGE;
                    pred_staircase_terminations(L) = GAMUTEDGE;
                    % Ignoring OOG points in the calculation of error for now.
                    if sum(L)/length(L) > 0.7 % if > 0.5 of the seaches go out of gamut, this is probably not a good fit (e.g. weights are too small)
                        err = nan;
                    else
                        err = mean((log(pred_staircase_terminations)-log(RHO(not_oog_idx)')).^2);
                    end
                    data(i,j,k) = err;
                    if ~isreal(err)
                        keyboard;
                    end
                end
            end
        end
    end
    [tmp_i, tmp_j, tmp_k] = ind2sub([length(A) length(B) length(C)],find(data==min(min(min(data)))));
    bestparams = [A(tmp_i) B(tmp_j) C(tmp_k)];
    preds = 1./(bestparams*[cos(allthetas).^2; sin(allthetas).^2; cos(allthetas).*sin(allthetas)]);
    LOOGtmp= preds>GAMUTEDGE.^2|preds<0;
    preds = sqrt(preds);
    [x1,y1] = pol2cart(allthetas(~LOOGtmp)', preds(~LOOGtmp)');
    [x2,y2] = pol2cart(thetas, RHO);
    figure(plot_counter); subplot(5,6,count); 
    plot(x1,y1,'g','Linewidth',2); hold on; 
    plot(x2(not_oog_idx), y2(not_oog_idx),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'PickableParts','none','MarkerEdgeColor',[0 0 1]);
    drawnow; axis equal; hold off;
    % plotting in log R - theta domain
    figure(plot_counter+2); subplot(5,6,count); 
    plot(allthetas(~LOOGtmp),log10(preds(~LOOGtmp)),'g','Linewidth',2); hold on; 
    plot(thetas(not_oog_idx), log10(RHO(not_oog_idx)),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'PickableParts','none','MarkerEdgeColor',[1 0 0]);
    drawnow; axis equal; hold off;
    count  = count + 1;
    if count == 31
        count = 1;
        plot_counter = plot_counter + 1;
    end
    
   preds_pts = 1./(bestparams*[cos(thetas).^2 sin(thetas).^2 cos(thetas).*sin(thetas)]');
   L1 = preds_pts > GAMUTEDGE.^2 | preds_pts < 0;
   preds_pts(L1) = GAMUTEDGE.^2;
   preds_pts = sqrt(preds_pts);
   SSE_quadmodel(ind) = sum((log(preds_pts)- log(RHO')).^2);
end
plot_counter = plot_counter + 2;