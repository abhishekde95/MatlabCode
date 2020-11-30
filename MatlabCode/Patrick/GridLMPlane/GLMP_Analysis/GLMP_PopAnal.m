%GLMP_PopAnal

% Created   5/13/14     JPW

% This code is designed to pull out some key features of the GLMP datasets
% and perform some analyses.

% This code was written in preparation for a potential abstract.

clear all
close all

datafiles = {
    % GLMP v1.0 (asymmetric cartesian grid)
    %'S102111002.nex'; % Pan color, mostly lum (#1)
    'S102811002.nex'; % +/- Luminance (#2)
    'S103111002.nex'; % -L-M (# assign number)
    'S103111003.nex'; % -L-M, with some +L+M (#3)
    'S110211007.nex'; % +/- Luminance  (Most Data) (#4)
    'S110911003.nex'; % +/- Luminance (#5)
    'S110911007.nex'; % +/- Luminance (#6)
    'S111011003.nex'; % -L-M, with some +L+M (Second to Most Data) (#7)
    'S112311006.nex'; % -L-M, with some +L+M (#8)
    'S120211008.nex'; % -L+M (# assign number)
    'S120511008.nex'; % +/- Luminance (#9)

    % GLMP v2.0 (asymmetric Radial grid)
    'S121211003.nex';% +/- Luminance (#10)
    'S121211006.nex';% +/- Luminance (#11)
    'S121511004.nex';% +/- Luminance (#12)

    % GLMP v2.1 (S-cone Plateau)
    %'S010212003.nex';% Possible Elipse (mostly unmodulated fr) (plat 1 is #13)
    'S011312002.nex';% +/- Luminance - maybe just off lum axis (plat 1 is #14)
    'S020112002.nex';% +/- Luminance (plat 1 is #15)
    'S020212007.nex';% -L-M, with some +L+M (plat 1 is #16)

    % GLMP v2.2 (2 Directions of Motion)
    'S021012007.nex';% +L+M, with some -L-M (plat 1 #17)
    'S021712006.nex';% Clean up - both plats look luminance

    % GLMP v3.2 (symmetric radial grid) (2 directions of motion)
    'S032912004.nex';% +/- Luminance (plat 2 is #18)
    'S040412002.nex';% R-G Chromatic (plat 2 is #19)
    %'S041012002.nex';% Ellipsoidal Cell!! (#20)
    
    % GLMP v3.3
    % Logrythmic spacing
    'S041912004.nex'; % Luminance (# assign number)
    'S042012002.nex'; % -L+M (# assign number)

    % GLMP v3.0
    %'A060112003.nex';% Ellispoidal Cell! (Apollo's First Cell!) (#21)
    %'A062612003.nex';% Ellispoidal Cell! (#22)
    'A070312003.nex';% +/- Luminance (#23)
    %'A070312005.nex';% Ellipsoidal Cell (#24)
    %'A071012004.nex';% Ellipsoidal (10 repeats) (#25)
    'A071212005.nex';% R-G Chromatic (#26)
    %'A071612002.nex';% Ellipsoidal Cell (#27)
    %'A071612005.nex';% Ellipsoidal Cell (Same cell as A071612002) (#28)
    
    % GLMP v4.0 (Adaptive + Nonadaptive)
    'S101812003.nex'; % Luminance! (plat 1 is # assign number)
    'S102312007.nex'; % +L-M Chromatic (plat 1 is # assign number)
    %datafile = 'S102412003.nex'; % Ellipse! (plat 1 is # assign number)
    %datafile = 'S102412004.nex'; % Maybe ellipse... sorting is strange. Work on this.
    'S110212005.nex'; % -L+M Chromatic (plat 1 is # assign number)
    
    };

% Set up large structure to hold values
popVals.plat = ones(numel(datafiles),1);
popVals.plat(end-9:end-8) = 2;
popVals.monkey = nan(numel(datafiles),1);
popVals.rf_x = nan(numel(datafiles),1);
popVals.rf_y = nan(numel(datafiles),1);


%% Receptive Field Location

for c = 1:numel(datafiles)
    
    datafile = datafiles{c};
    %library = 'Users/jpatrickweller/Documents/Matlab_Working_Directory/Patrick/GridLMPlane/GLMPDatafiles/';
    %library = '/Users/jpatrickweller/Dropbox/Patrick/';
    library = 'C:\Users\jpweller\Dropbox\Patrick\DatafilesAll\';

    disp(['Now processing ' datafile])
    rawdata = nex2stro([char(library) datafile]);
    
    % Pull out some key features of the cell
    popVals.monkey(c) = datafile(1);
    popVals.rf_x(c) = rawdata.sum.exptParams.rf_x; 
    popVals.rf_y(c) = rawdata.sum.exptParams.rf_y;
    
    % Plotting retinal eccentricity
    figure(3); hold on; grid on;
    if strcmp(datafile(1),'S')
        plot(rawdata.sum.exptParams.rf_x,rawdata.sum.exptParams.rf_y,'r*')
    elseif strcmp(datafile(1),'A')
        plot(rawdata.sum.exptParams.rf_x,rawdata.sum.exptParams.rf_y,'b*')
    end
    
end


%% 1D Fits

% Define some fields
popVals.angs = cell(numel(datafiles),1);
popVals.GOF = cell(numel(datafiles),1);
popVals.params = cell(numel(datafiles),1);
popVals.bestIdx = nan(numel(datafiles),1);


% Set up figure
figure(60); clf;
set(gcf,'units','pixels','pos',[200 200 1100 550],'NumberTitle','off',...
    'Name',['Surface Analyses (' datafiles{1} ')']);
surffig.disppanel = uipanel('Pos',[.025 .025 .6 .95],'Parent',gcf);
surffig.poppanel = uipanel('pos',[.65 .025 .325 .95],'parent',gcf);
poppanel.axes.direction = axes('parent',surffig.poppanel,'units','normalized',...
    'pos',[.1 .69 .8 .25],'visible','off');
poppanel.axes.sigma = axes('parent',surffig.poppanel,'units','normalized',...
    'pos',[.1 .37 .8 .25],'box','on');
poppanel.axes.exp = axes('parent',surffig.poppanel,'units','normalized',...
    'pos',[.1 .05 .8 .25],'box','on');
poppanel.labels.direction = uicontrol('style','text','parent',surffig.poppanel,...
    'units','normalized','pos',[.1 .95 .8 .025],'string','Best Directions','fontsize',12);
poppanel.labels.sigma = uicontrol('style','text','parent',surffig.poppanel,...
    'units','normalized','pos',[.1 .63 .8 .025],'string','Sigma Values','fontsize',12);
poppanel.labels.exp = uicontrol('style','text','parent',surffig.poppanel,...
    'units','normalized','pos',[.1 .31 .8 .025],'string','Exponent Values','fontsize',12);
disppanel.axes.surf3d = axes('parent',surffig.disppanel,'units','pixels',...
    'pos',[50 200 250 250],'box','on','XTick',[],'YTick',[]);
disppanel.labels.surf3d = uicontrol('style','text','parent',surffig.disppanel,...
    'units','pixels','pos',[50 460 250 20],'string','Surface Fit','fontsize',12);
disppanel.axes.proj2d = axes('parent',surffig.disppanel,'units','pixels',...
    'pos',[375 200 250 250],'box','on','XTick',[],'YTick',[]);
disppanel.labels.proj2d = uicontrol('style','text','parent',surffig.disppanel,...
    'units','pixels','pos',[375 460 250 20],'string','2D Projection','fontsize',12);
disppanel.axes.GOF = axes('parent',surffig.disppanel,'units','pixels',...
    'pos',[80 40 555 100],'box','on','XTick',[],'YTick',[]);
disppanel.labels.GOF = uicontrol('style','text','parent',surffig.disppanel,...
    'units','pixels','pos',[80 140 555 20],'string','Goodness of Fit','fontsize',12);

for c = 1:numel(datafiles)
    
    datafile = datafiles{c};
    %library = 'Users/jpatrickweller/Documents/Matlab_Working_Directory/Patrick/GridLMPlane/GLMPDatafiles/';
    %library = '/Users/jpatrickweller/Dropbox/Patrick/';
    library = 'C:\Users\jpweller\Dropbox\Patrick\DatafilesAll\';

    disp(['Now processing ' datafile])
    rawdata = nex2stro([char(library) datafile]);
    [trial par plat] = OrganizeRawGLMPData(rawdata);
    p = popVals.plat(c);
    set(gcf,'Name',['Surface Analyses (' datafiles{c} ')'])
        
    % Set up some variables
    %startang = randi(180)/180*pi;
    startang = 0;
    angs = linspace(startang,startang+pi,181)';
    nrots = numel(angs);
    GOF = nan(nrots,1);
    vlb = [0     0     .1   .1   1    1  .0001];
    vub = [1000 1000    2    2  10   10     50];
    options = optimset('MaxFunEvals',5000,'MaxIter',5000,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
    sigmaguess = 0.4;
    expguess = 2;
    baselineguess = 5;
    params = nan(nrots,7);
    
    % Generating an initial guess
    maxx = [];
    topfr = max(plat{p}.trial.nspikes);
    if topfr == 0
        topfr = baselineguess;
    end
    params0 = [topfr, topfr, sigmaguess, sigmaguess, expguess, expguess, baselineguess];  % need to constrain B across color directions
    
    % Rotate through angles
    for rot = 1:nrots
        
        % Rotate L and M coordinates to make the 'fitting axis' the x-axis
        rotMat = [cos(angs(rot)) -sin(angs(rot)); sin(angs(rot)) cos(angs(rot))];
        tempRotPts = [rotMat * [plat{p}.trial.Lcc_orig plat{p}.trial.Mcc_orig]']';
        projs = tempRotPts(:,1);
        spikes = plat{p}.trial.nspikes;
        
        % Use previous parameters
        if rot>1
            params1 = mean([params(rot-1,:);params0]);
        else
            params1 = params0;
        end
        
        % Fit all rotated points
        [f1,fval] = fmincon('FitNakaRushtonFunJPW',params1,[],[],[],[],vlb,vub,[],options,projs,spikes,'asymmetric');
        params(rot,:) = f1;
        GOF(rot) = -fval;
        
        %Variables for plotting
        maxx = cat(1,maxx,max(projs));
        pts = min(projs):.001:max(projs);
        fitPts = ComputeNakaRushtonJPW(params(rot,:),pts,'asymmetric');
        axlim = max(maxx);
        projline = [rotMat\[-axlim 0; 0 0; axlim 0]']';
        [x,y] = meshgrid(linspace(-axlim,axlim,50));
        rotxy = [rotMat * [x(:) y(:)]']';
        surface = ComputeNakaRushtonJPW(params(rot,:),rotxy(:,1),'asymmetric');
        surface = reshape(surface,size(x));

        
        % Plot figure
        axes(disppanel.axes.surf3d); cla; hold on; grid on;
        surf(x,y,surface)
        plot3(plat{p}.trial.Lcc_orig,plat{p}.trial.Mcc_orig,plat{p}.trial.nspikes,'k*')
        plot3(projline(:,1),projline(:,2),[topfr 0 topfr]','m')
        set(gca,'XLim',[-axlim axlim],'XTickMode','auto','YLim',[-axlim axlim],'YTickMode','auto')
        xlabel('L Cone Contrast')
        ylabel('M Cone Contrast')
        zlabel('Number of Spikes')
        
        % Plot Projections
        axes(disppanel.axes.proj2d);
        cla; hold on; grid on;
        xlim([-max(maxx) max(maxx)])
        set(gca,'XTickMode','auto','YTickMode','auto')
        plot(projs,plat{p}.trial.nspikes,'k*');
        plot(pts,fitPts,'r--')
        xlabel('Cone Contrast')
        ylabel('Number of Spikes')
        
        % Plot Goodness of Fit
        axes(disppanel.axes.GOF)
        if rot == 1
            cla; hold on; grid on;
            xlim([min(angs)/pi*180 max(angs)/pi*180])
            set(gca,'XTick',[linspace(min(angs),max(angs),5)/pi*180],'YTickMode','auto')
            xlabel('Rotation of Fitting Axis')
            ylabel('Negative Log Likelihood')
        else
            plot(angs(2:rot)/pi*180,GOF(2:rot),'ko-')
        end
    end
    
    %Variables for best fitting surface
    [~,bestIdx] = max(GOF);
    rotMat = [cos(angs(bestIdx)) -sin(angs(bestIdx)); sin(angs(bestIdx)) cos(angs(bestIdx))];
    projline = [inv(rotMat) * [-axlim 0; 0 0; axlim 0]']';
    [x,y] = meshgrid(linspace(-axlim,axlim,50));
    rotxy = [rotMat * [x(:) y(:)]']';
    surface = ComputeNakaRushtonJPW(params(bestIdx,:),rotxy(:,1),'asymmetric');
    surface = reshape(surface,size(x));
    fitPts = ComputeNakaRushtonJPW(params(bestIdx,:),pts,'asymmetric');
    tempRotPts = [rotMat * [plat{p}.trial.Lcc_orig plat{p}.trial.Mcc_orig]']';
    projs = tempRotPts(:,1);
    
    % Plot best 3d surface
    axes(disppanel.axes.surf3d); cla; hold on; grid on;
    surf(x,y,surface)
    plot3(plat{p}.trial.Lcc_orig,plat{p}.trial.Mcc_orig,plat{p}.trial.nspikes,'k*')
    plot3(projline(:,1),projline(:,2),[topfr 0 topfr]','m')
    set(gca,'XLim',[-axlim axlim],'XTickMode','auto','YLim',[-axlim axlim],'YTickMode','auto')
    xlabel('L Cone Contrast')
    ylabel('M Cone Contrast')
    zlabel('Number of Spikes')
    
    % Plot Projections
    axes(disppanel.axes.proj2d)
    cla; hold on; grid on;
    xlim([-max(maxx) max(maxx)])
    plot(projs,plat{p}.trial.nspikes,'k*');
    plot(pts,fitPts,'r--')
    set(gca,'XTickMode','auto','YTickMode','auto');
    xlabel('Cone Contrast')
    ylabel('Number of Spikes')
    
    % Indicate best fitting surface
    axes(disppanel.axes.GOF); cla; hold on; grid on;
    h = plot(angs(2:end)/pi*180,GOF(2:end),'ko-');
    plot(angs(bestIdx)/pi*180,GOF(bestIdx),'r*')
    %set(h,'ButtonDownFcn',@SurfaceSelectionFit1Ax);
    
    % Save results
    popVals.angs{c} = angs';
    popVals.GOF{c} = GOF;
    popVals.params{c} = params(bestIdx,:);
    popVals.bestIdx(c) = bestIdx;
    
    % Plot population results
    axes(poppanel.axes.direction);
    %hist(angs(popVals.bestIdx(1:c)),angs)
    %xlim([0 pi])
    %keyboard
    rose(pi-angs(popVals.bestIdx(1:c)))
    axis equal tight
    
    axes(poppanel.axes.sigma); cla; hold on; grid on;
    paramsMat = cell2mat(popVals.params);
    hist(paramsMat(:,3:4))

    axes(poppanel.axes.exp); cla; hold on; grid on;
    hist(paramsMat(:,5:6))
    
end







