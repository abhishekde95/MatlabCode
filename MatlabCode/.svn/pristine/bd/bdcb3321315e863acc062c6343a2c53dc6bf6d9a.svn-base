function [] = GLMSGUI_Surface_highlum(~,~)
global GLMP

% This code fits surfaces to GLMP datasets.

% Set up figure
figure(160); clf;
set(gcf,'units','pixels','pos',[200 200 1100 550],'NumberTitle','off',...
    'Name',['Surface Analyses -high lum- (' GLMP.datafile ')']);
surffig = get(gcf,'UserData');
surffig.disppanel = uipanel('Pos',[.025 .025 .6 .95],'Parent',gcf);
surffig.conpanel = uipanel('pos',[.65 .025 .325 .95],'parent',gcf); 

disppanel = get(surffig.disppanel,'UserData');
conpanel = get(surffig.conpanel,'UserData');


% Set up control panel
conpanel.analyses.buttons = uibuttongroup('parent',surffig.conpanel,...
    'units','normalized','pos',[.025 .825 .55 .15]);
conpanel.analyses.fit = uicontrol('style','radiobutton',...
    'parent',conpanel.analyses.buttons,'string','Fit to All Data',...
    'fontsize',12,'units','normalized','pos',[.05 .6 .9 .3]);
conpanel.analyses.predict = uicontrol('style','radiobutton',...
    'parent',conpanel.analyses.buttons,'string','Predict from Subset of Data',...
    'fontsize',12,'units','normalized','pos',[.05 .2 .9 .3]);
conpanel.analprops.naxes = uibuttongroup('parent',surffig.conpanel,...
    'units','normalized','pos',[.025 .65 .55 .15],...
    'SelectionChangeFcn',@OrthogOnOff);
conpanel.analprops.nax1 = uicontrol('style','radiobutton',...
    'parent',conpanel.analprops.naxes,'string','One Axis',...
    'fontsize',12,'units','normalized','pos',[.05 .6 .9 .3]);
conpanel.analprops.nax2 = uicontrol('style','radiobutton',...
    'parent',conpanel.analprops.naxes,'string','Two Axes',...
    'fontsize',12,'units','normalized','pos',[.05 .2 .9 .3]);
conpanel.axisprops.orthog = uicontrol('style','checkbox',...
    'parent',surffig.conpanel,'units','normalized','pos',[.05 .5 .9 .1],...
    'string','Constrain Axes to be Orthogonal','fontsize',12,'enable','off');
conpanel.startanalysis = uicontrol('style','pushbutton',...
    'parent',surffig.conpanel,'units','pixels','pos',[50 100 100 50],...
    'string','Start Analysis','fontsize',12,...
    'callback',@startanalysis);
conpanel.subbuttons = uibuttongroup('Parent',surffig.conpanel,...
    'units','normalized','pos',[.6 .775 .35 .2]);
conpanel.sub1 = uicontrol('style','radiobutton',...
    'parent',conpanel.subbuttons,'pos',[15 65 100 25],...
    'string','Subunit #1','fontsize',12);
conpanel.sub2 = uicontrol('style','radiobutton',...
    'parent',conpanel.subbuttons,'pos',[15 37 100 25],...
    'string','Subunit #2','fontsize',12);
conpanel.sub3 = uicontrol('style','radiobutton',...
    'parent',conpanel.subbuttons,'pos',[15 10 100 25],...
    'string','Subunit #3','fontsize',12);
if numel(GLMP.subunit) == 1
    set(conpanel.sub2,'enable','off')
    set(conpanel.sub3,'enable','off')
end


% Set up display panel
disppanel.axes.surf3d = axes('parent',surffig.disppanel,'units','pixels',...
    'pos',[50 200 250 250],'box','on','XTick',[],'YTick',[]);
disppanel.labels.surf3d = uicontrol('style','text','parent',surffig.disppanel,...
    'units','pixels','pos',[50 460 250 20],'string','Surface Fit','fontsize',12);
disppanel.axes.proj2d = axes('parent',surffig.disppanel,'units','pixels',...
    'pos',[375 200 250 250],'box','on','XTick',[],'YTick',[]);
disppanel.labels.proj2d = uicontrol('style','text','parent',surffig.disppanel,...
    'units','pixels','pos',[375 460 250 20],'string','2D Projection','fontsize',12);
disppanel.axes.GOF = axes('parent',surffig.disppanel,'units','pixels',...
    'pos',[75 30 555 100],'box','on','XTick',[],'YTick',[]);
disppanel.labels.GOF = uicontrol('style','text','parent',surffig.disppanel,...
    'units','pixels','pos',[75 130 555 20],'string','Goodness of Fit','fontsize',12);


%Save variables
set(surffig.conpanel,'UserData',conpanel);
set(surffig.disppanel,'UserData',disppanel);
set(gcf,'UserData',surffig);


end

function OrthogOnOff(~,~)

% Load Figure Variables
surffig = get(gcf,'UserData');
conpanel = get(surffig.conpanel,'UserData');
disppanel = get(surffig.disppanel,'UserData');

if get(conpanel.analprops.naxes,'SelectedObject') == conpanel.analprops.nax2
    set(conpanel.axisprops.orthog,'Enable','on');
else
    set(conpanel.axisprops.orthog,'Enable','off');
end

%Save variables
set(surffig.conpanel,'UserData',conpanel);
set(surffig.disppanel,'UserData',disppanel);
set(gcf,'UserData',surffig);

end


function startanalysis(~,~)

% I think this should route the analysis to different functions
surffig = get(gcf,'UserData');
conpanel = get(surffig.conpanel,'UserData');
disppanel = get(surffig.disppanel,'UserData');

fit = 0;
pred = 0;
nax = 0;

if get(conpanel.analyses.buttons,'selectedobject') == conpanel.analyses.fit
    fit = 1;
elseif get(conpanel.analyses.buttons,'selectedobject') == conpanel.analyses.predict
    pred = 1;
end

if get(conpanel.analprops.naxes,'selectedobject') == conpanel.analprops.nax1
    nax = 1;
elseif get(conpanel.analprops.naxes,'selectedobject') == conpanel.analprops.nax2
    nax = 2;
end

orthog = get(conpanel.axisprops.orthog,'value');

if fit == 1 && nax == 1
    Fit1Ax
elseif pred == 1 && nax == 1
    Pred1Ax
elseif fit == 1 && nax == 2 && orthog == 1
    Fit2AxOrthog
else
    disp('Analysis does not yet exist.')
end

    
    
end


function Fit1Ax
global GLMP

% Load Figure Variables
surffig = get(gcf,'UserData');
conpanel = get(surffig.conpanel,'UserData');
disppanel = get(surffig.disppanel,'UserData');

% Set up some variables
angs = linspace(0,pi,181)';
nrots = numel(angs);
GOF = nan(nrots,7,numel(GLMP.subunit));
vlb = [0     0   0.0001 0.0001  0  0  0];
vub = [1000 1000  100    100   10 10 50];
options = optimset('MaxFunEvals',300,'MaxIter',300,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
sigmaguess = 0.1;
baselineguess = 0;
params = nan(nrots,7);

% Find the selected subunit
if get(conpanel.subbuttons,'SelectedObject') == conpanel.sub1
    p = 1;
elseif get(conpanel.subbuttons,'SelectedObject') == conpanel.sub2
    p = 2;
elseif get(conpanel.subbuttons,'SelectedObject') == conpanel.sub3
    p = 3;
end

% Generating an initial guess
maxx = [];
topfr = max(GLMP.subunit{p}.specialcases.highlum.nspikes);
if (topfr == 0)
    params0 = [0 0 1 1 2 2 0];
elseif isempty(topfr)
    params0 = [0 0 1 1 2 2 0];
else
    params0 = [topfr, topfr, sigmaguess, sigmaguess, 2, 2, baselineguess];  % need to constrain B across color directions
end

% Rotate through angles
for rot = 1:nrots
    
    % Rotate L and M coordinates to make the 'fitting axis' the x-axis
    rotMat = [cos(angs(rot)) -sin(angs(rot)); sin(angs(rot)) cos(angs(rot))];
    tempRotPts = [rotMat * [GLMP.subunit{p}.specialcases.highlum.Lcc GLMP.subunit{p}.specialcases.highlum.Mcc]']';
    projs = tempRotPts(:,1);
    
%     figure(100000); clf; hold on; grid on; axis equal
%     plot(GLMP.subunit{p}.Lcc(1:10),GLMP.subunit{p}.Mcc(1:10),'bo')
%     figure(100001); clf; hold on; axis equal
%     polar(GLMP.subunit{p}.theta(1:10),GLMP.subunit{p}.rho(1:10),'ko')
%     plot(tempRotPts(1:10,1),tempRotPts(1:10,2),'r*')
    
    % Use previous parameters
    if rot>1
        params1 = mean([params(rot-1,:,p);params0]);
    else
        params1 = params0;
    end
    
    % Fit all rotated points
    [f1,fval] = fmincon('FitNakaRushtonFunJPW',params1,[],[],[],[],vlb,vub,[],options,projs,GLMP.subunit{p}.specialcases.highlum.nspikes,'asymmetric');
    params(rot,:,p) = f1;
    GOF(rot,p) = -fval;
    
    %Variables for plotting
    maxx = cat(1,maxx,max(projs));
    pts = min(projs):.001:max(projs);
    fitPts = ComputeNakaRushtonJPW(params(rot,:,p),pts,'asymmetric');
    axlim = max(maxx);
    projline = [rotMat\[-axlim 0; 0 0; axlim 0]']';
    [x,y] = meshgrid(linspace(-axlim,axlim,50));
    rotxy = [rotMat * [x(:) y(:)]']';
    surface = ComputeNakaRushtonJPW(params(rot,:),rotxy(:,1),'asymmetric');
    surface = reshape(surface,size(x));
    
    % Plot figure
    axes(disppanel.axes.surf3d); cla; hold on; grid on;
    surf(x,y,surface)
    plot3(GLMP.subunit{p}.specialcases.highlum.Lcc,GLMP.subunit{p}.specialcases.highlum.Mcc,GLMP.subunit{p}.specialcases.highlum.nspikes,'k*')
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
    plot(projs,GLMP.subunit{p}.specialcases.highlum.nspikes,'k*');
    plot(pts,fitPts,'r--')
    xlabel('Cone Contrast')
    ylabel('Number of Spikes')
    
    % Plot Goodness of Fit
    axes(disppanel.axes.GOF)
    if rot == 1
        cla; hold on; grid on;
        xlim([min(angs)/pi*180 max(angs)/pi*180])
        set(gca,'XTick',[0 45 90 135 180],'YTickMode','auto')
        xlabel('Rotation of Fitting Axis')
        ylabel('Negative Log Likelihood')
    end
    plot(angs(1:rot)/pi*180,GOF(1:rot,p),'ko-')
end

%Variables for best fitting surface
[~,bestIdx] = max(GOF(:,p));
rotMat = [cos(angs(bestIdx)) -sin(angs(bestIdx)); sin(angs(bestIdx)) cos(angs(bestIdx))];
projline = [inv(rotMat) * [-axlim 0; 0 0; axlim 0]']';
[x,y] = meshgrid(linspace(-axlim,axlim,50));
rotxy = [rotMat * [x(:) y(:)]']';
surface = ComputeNakaRushtonJPW(params(bestIdx,:),rotxy(:,1),'asymmetric');
surface = reshape(surface,size(x));
fitPts = ComputeNakaRushtonJPW(params(bestIdx,:),pts,'asymmetric');
tempRotPts = [rotMat * [GLMP.subunit{p}.specialcases.highlum.Lcc GLMP.subunit{p}.specialcases.highlum.Mcc]']';
projs = tempRotPts(:,1);

% Plot best 3d surface
axes(disppanel.axes.surf3d); cla; hold on; grid on;
surf(x,y,surface)
plot3(GLMP.subunit{p}.specialcases.highlum.Lcc,GLMP.subunit{p}.specialcases.highlum.Mcc,GLMP.subunit{p}.specialcases.highlum.nspikes,'k*')
plot3(projline(:,1),projline(:,2),[topfr 0 topfr]','m')
set(gca,'XLim',[-axlim axlim],'XTickMode','auto','YLim',[-axlim axlim],'YTickMode','auto')
xlabel('L Cone Contrast')
ylabel('M Cone Contrast')
zlabel('Number of Spikes')
    
% Plot Projections
axes(disppanel.axes.proj2d)
cla; hold on; grid on;
xlim([-max(maxx) max(maxx)])
plot(projs,GLMP.subunit{p}.specialcases.highlum.nspikes,'k*');
plot(pts,fitPts,'r--')
set(gca,'XTickMode','auto','YTickMode','auto');
xlabel('Cone Contrast')
ylabel('Number of Spikes')

% Indicate best fitting surface
axes(disppanel.axes.GOF); cla; hold on; grid on;
h = plot(angs/pi*180,GOF(:,p),'ko-');
plot(angs(bestIdx)/pi*180,GOF(bestIdx,p),'r*')
set(h,'ButtonDownFcn',@SurfaceSelectionFit1Ax);

% Save results
disppanel.results.fit1ax.angs{p} = angs;
disppanel.results.fit1ax.GOF{p} = GOF;
disppanel.results.fit1ax.fitparams{p} = params;
disppanel.results.fit1ax.bestIdx(p) = bestIdx;

%Save variables
set(surffig.conpanel,'UserData',conpanel);
set(surffig.disppanel,'UserData',disppanel);
set(gcf,'UserData',surffig);
    
end

function SurfaceSelectionFit1Ax(~,~)
global GLMP

h = gca;
whichpt = get(h,'CurrentPoint');
whichpt = whichpt(1,[1 2]);

% Load variables
surffig = get(gcf,'UserData');
conpanel = get(surffig.conpanel,'UserData');
disppanel = get(surffig.disppanel,'UserData');

% Find the selected subunit
if get(conpanel.subbuttons,'SelectedObject') == conpanel.sub1
    whichsub = 1;
elseif get(conpanel.subbuttons,'SelectedObject') == conpanel.sub2
    whichsub = 2;
elseif get(conpanel.subbuttons,'SelectedObject') == conpanel.sub3
    whichsub = 3;
end

% Plot selected point
GOF = disppanel.results.fit1ax.GOF{whichsub};
angs = disppanel.results.fit1ax.angs{whichsub};
bestIdx = disppanel.results.fit1ax.bestIdx(whichsub);
params = disppanel.results.fit1ax.fitparams{whichsub};

[~,idx] = min(sum(abs(angs - repmat(whichpt(1)/180*pi,size(angs,1),1)),2));
axlim = max(GLMP.subunit{whichsub}.specialcases.highlum.rho);
rotMat = [cos(angs(idx)) -sin(angs(idx)); sin(angs(idx)) cos(angs(idx))];
projline = [rotMat\[-axlim 0; 0 0; axlim 0]']';
[x,y] = meshgrid(linspace(-axlim,axlim,50));
rotxy = [rotMat*[x(:) y(:)]']';
surface = ComputeNakaRushtonJPW(params(idx,:),rotxy(:,1),'asymmetric');
surface = reshape(surface,size(x));
tempRotPts = [rotMat*[GLMP.subunit{whichsub}.specialcases.highlum.Lcc GLMP.subunit{whichsub}.specialcases.highlum.Mcc]']';
projs = tempRotPts(:,1);
pts = min(projs):.001:max(projs);
fitPts = ComputeNakaRushtonJPW(params(idx,:),pts,'asymmetric');
topnsp = max(GLMP.subunit{whichsub}.specialcases.highlum.nspikes);

% Plot GOF selection
axes(disppanel.axes.GOF)
cla; hold on;
h0 = plot(angs/pi*180,GOF,'ko-');
h1 = plot(angs(bestIdx)/pi*180,GOF(bestIdx),'r*');
disppanel.surfsel = plot(angs(idx)/pi*180,GOF(idx),'c*');
set(h0,'ButtonDownFcn',@SurfaceSelectionFit1Ax);
set(h1,'ButtonDownFcn',@SurfaceSelectionFit1Ax);

% Plot selected surface
axes(disppanel.axes.surf3d)
cla; hold on; grid on;
surf(x,y,surface)
plot3(GLMP.subunit{whichsub}.specialcases.highlum.Lcc,GLMP.subunit{whichsub}.specialcases.highlum.Mcc,GLMP.subunit{whichsub}.specialcases.highlum.nspikes,'k*')
plot3(projline(:,1),projline(:,2),[topnsp 0 topnsp]','m')
set(gca,'XLim',[-axlim axlim],'XTickMode','auto','YLim',[-axlim axlim],'YTickMode','auto')
xlabel('L Cone Contrast')
ylabel('M Cone Contrast')
zlabel('Number of Spikes')

% Plot selected CRF
axes(disppanel.axes.proj2d); cla; hold on; grid on;
plot(projs,GLMP.subunit{whichsub}.specialcases.highlum.nspikes,'k*')
plot(pts,fitPts,'r--')

%Save variables
set(surffig.conpanel,'UserData',conpanel);
set(surffig.disppanel,'UserData',disppanel);
set(gcf,'UserData',surffig);

end


function Pred1Ax
global GLMP

% Load Figure Variables
surffig = get(gcf,'UserData');
conpanel = get(surffig.conpanel,'UserData');
disppanel = get(surffig.disppanel,'UserData');


% Find the selected subunit
if get(conpanel.subbuttons,'SelectedObject') == conpanel.sub1
    p = 1;
elseif get(conpanel.subbuttons,'SelectedObject') == conpanel.sub2
    p = 2;
elseif get(conpanel.subbuttons,'SelectedObject') == conpanel.sub3
    p = 3;
end

% Set up some variables
angs = unique(GLMP.subunit{p}.specialcases.highlum.theta);
angs = angs(angs<=0);
angs = sort(angs,'descend');
nrots = numel(angs);
GOF = nan(nrots,numel(GLMP.subunit));
vlb = [0     0   0.0001 0.0001  0  0  0];
vub = [1000 1000  100    100   10 10 50];
options = optimset('MaxFunEvals',5000,'MaxIter',5000,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
sigmaguess = 0.1;
baselineguess = 0;
params = nan(nrots,7,numel(GLMP.subunit));

% Generating an initial guess
topfr = max(GLMP.subunit{p}.specialcases.highlum.nspikes);
if (topfr == 0)
    params0 = [0 0 1 1 2 2 0];
elseif isempty(topfr)
    params0 = [0 0 1 1 2 2 0];
else
    params0 = [topfr, topfr, sigmaguess, sigmaguess, 2, 2, baselineguess];  % need to constrain B across color directions
end

% Rotate through angles
for rot = 1:nrots
    
    % Rotate L and M coordinates to make the 'fitting axis' the x-axis
    rotL = GLMP.subunit{p}.specialcases.highlum.theta == angs(rot)...
        | softEq(GLMP.subunit{p}.specialcases.highlum.theta,angs(rot)+pi)...
        | softEq(GLMP.subunit{p}.specialcases.highlum.theta,angs(rot)-pi);
    rotMat = [cos(angs(rot)) -sin(angs(rot)); sin(angs(rot)) cos(angs(rot))];
    tempRotPts = [rotMat\[GLMP.subunit{p}.specialcases.highlum.Lcc GLMP.subunit{p}.specialcases.highlum.Mcc]']';
    projs = tempRotPts(rotL,1);

%     figure(100001); clf; hold on; axis equal
%     temprotL = find(rotL)
%     polar(GLMP.subunit{p}.theta(temprotL),GLMP.subunit{p}.rho(temprotL),'ko')
%     plot(tempRotPts(temprotL,1),tempRotPts(temprotL,2),'r*')
%     polar(GLMP.subunit{p}.theta(temprotL(1)),GLMP.subunit{p}.rho(temprotL(1)),'k*')
%     plot(tempRotPts(temprotL(1),1),tempRotPts(temprotL(1),2),'ro')
     
    % Use previous parameters
    if rot>1
        params1 = mean([params(rot-1,:);params0]);
    else
        params1 = params0;
    end
    
    % Fit all rotated points
    [f1,fval] = fmincon('FitNakaRushtonFunJPW',params1,[],[],[],[],vlb,vub,[],options,projs,GLMP.subunit{p}.specialcases.highlum.nspikes(rotL),'asymmetric');
    params(rot,:,p) = f1;
    
    % Calculate SSE as a metric for GOF
    rotstim = [rotMat\[GLMP.subunit{p}.specialcases.highlum.Lcc GLMP.subunit{p}.specialcases.highlum.Mcc]']';
    predresp = ComputeNakaRushtonJPW(params(rot,:),rotstim(:,1),'asymmetric');
    GOF(rot,p) = sum((predresp - GLMP.subunit{p}.specialcases.highlum.nspikes).^2);
    
    %Variables for plotting
    axlim = max(GLMP.subunit{p}.specialcases.highlum.rho);
    pts = -axlim:.001:axlim;
    fitPts = ComputeNakaRushtonJPW(params(rot,:),pts,'asymmetric');
    axlim = max(GLMP.subunit{p}.specialcases.highlum.rho);
    projline = [rotMat*[-axlim 0; 0 0; axlim 0]']';
    [x,y] = meshgrid(linspace(-axlim,axlim,50));
    rotxy = [rotMat\[x(:) y(:)]']';
    surface = ComputeNakaRushtonJPW(params(rot,:),rotxy(:,1),'asymmetric');
    surface = reshape(surface,size(x));
    
    % Plot figure
    axes(disppanel.axes.surf3d); cla; hold on; grid on;
    surf(x,y,surface)
    plot3(GLMP.subunit{p}.specialcases.highlum.Lcc,GLMP.subunit{p}.specialcases.highlum.Mcc,GLMP.subunit{p}.specialcases.highlum.nspikes,'k*')
    plot3(projline(:,1),projline(:,2),[topfr 0 topfr]','m')
    set(gca,'XLim',[-axlim axlim],'XTickMode','auto','YLim',[-axlim axlim],'YTickMode','auto')
    xlabel('L Cone Contrast')
    ylabel('M Cone Contrast')
    zlabel('Number of Spikes')
    

    % Plot Projections
    axes(disppanel.axes.proj2d);
    cla; hold on; grid on;
    xlim([-axlim axlim])
    set(gca,'XTickMode','auto','YTickMode','auto');
    plot(projs,GLMP.subunit{p}.specialcases.highlum.nspikes(rotL),'k*');
    plot(pts,fitPts,'r--')
    ylabel('Number of Spikes')
    xlabel('Cone Cotrast')
    
    % Plot Goodness of Fit
    axes(disppanel.axes.GOF)
    if rot == 1
        cla; hold on; grid on;
        xlim([0 180])
        set(gca,'XTick',linspace(min(xlim),max(xlim),10),'YTickMode','auto')
        xlabel('Rotation of Fitting Line (deg)')
        ylabel('Sum of Squarred Error')
    end
    plot(abs((angs(1:rot)))/pi*180,GOF(1:rot,p),'o-')

end

%Variables for best fitting surface
[~,bestIdx] = min(GOF(:,p));
rotMat = [cos(angs(bestIdx)) -sin(angs(bestIdx)); sin(angs(bestIdx)) cos(angs(bestIdx))];
projline = [rotMat*[-axlim 0; 0 0; axlim 0]']';
[x,y] = meshgrid(linspace(-axlim,axlim,50));
rotxy = [rotMat\[x(:) y(:)]']';
surface = ComputeNakaRushtonJPW(params(bestIdx,:),rotxy(:,1),'asymmetric');
surface = reshape(surface,size(x));
fitPts = ComputeNakaRushtonJPW(params(bestIdx,:),pts,'asymmetric');
tempRotPts = [rotMat\[GLMP.subunit{p}.specialcases.highlum.Lcc GLMP.subunit{p}.specialcases.highlum.Mcc]']';
projs = tempRotPts(:,1);

% Plot best 3d surface
axes(disppanel.axes.surf3d); cla; hold on; grid on;
surf(x,y,surface)
plot3(GLMP.subunit{p}.specialcases.highlum.Lcc,GLMP.subunit{p}.specialcases.highlum.Mcc,GLMP.subunit{p}.specialcases.highlum.nspikes,'k*')
plot3(projline(:,1),projline(:,2),[topfr 0 topfr]','m')
set(gca,'XLim',[-axlim axlim],'XTickMode','auto','YLim',[-axlim axlim],'YTickMode','auto')
xlabel('L Cone Contrast')
ylabel('M Cone Contrast')
zlabel('Number of Spikes')
    
% Plot Projections
axes(disppanel.axes.proj2d)
cla; hold on; grid on
xlim([-axlim axlim])
set(gca,'XTickMode','auto','YTickMode','auto')
plot(projs,GLMP.subunit{p}.specialcases.highlum.nspikes,'k*');
plot(pts,fitPts,'r--')
ylabel('Number of Spikes')
xlabel('Cone Cotrast')

% Indicate best fitting surface
axes(disppanel.axes.GOF); cla; hold on; grid on;
plot(abs(angs(bestIdx)/pi*180),GOF(bestIdx,p),'r*')
h = plot(abs(angs/pi*180),GOF(:,p),'ko-');
set(h,'ButtonDownFcn',@SurfaceSelectionPred1Ax);

% Save results
disppanel.results.pred1ax.GOF{p} = GOF;
disppanel.results.pred1ax.fitparams{p} = params;
disppanel.results.pred1ax.angs{p} = angs;
disppanel.results.pred1ax.bestIdx(p) = bestIdx;

%Save variables
set(surffig.conpanel,'UserData',conpanel);
set(surffig.disppanel,'UserData',disppanel);
set(gcf,'UserData',surffig);
    
end


function SurfaceSelectionPred1Ax(~,~)
global GLMP

h = gca;
whichpt = get(h,'CurrentPoint');
whichpt = whichpt(1,[1 2]);

% Load variables
surffig = get(gcf,'UserData');
conpanel = get(surffig.conpanel,'UserData');
disppanel = get(surffig.disppanel,'UserData');

% Find the selected subunit
if get(conpanel.subbuttons,'SelectedObject') == conpanel.sub1
    whichsub = 1;
elseif get(conpanel.subbuttons,'SelectedObject') == conpanel.sub2
    whichsub = 2;
elseif get(conpanel.subbuttons,'SelectedObject') == conpanel.sub3
    whichsub = 3;
end

% Plot selected point
GOF = disppanel.results.pred1ax.GOF{whichsub};
angs = disppanel.results.pred1ax.angs{whichsub};
bestIdx = disppanel.results.pred1ax.bestIdx(whichsub);
params = disppanel.results.pred1ax.fitparams{whichsub};

[~,idx] = min(sum(abs(abs(angs) - repmat(whichpt(1)/180*pi,size(angs,1),1)),2));
axlim = max(GLMP.subunit{whichsub}.specialcases.highlum.rho);
rotMat = [cos(angs(idx)) -sin(angs(idx)); sin(angs(idx)) cos(angs(idx))];
projline = [rotMat*[-axlim 0; 0 0; axlim 0]']';
[x,y] = meshgrid(linspace(-axlim,axlim,50));
rotxy = [rotMat\[x(:) y(:)]']';
surface = ComputeNakaRushtonJPW(params(idx,:),rotxy(:,1),'asymmetric');
surface = reshape(surface,size(x));
tempRotPts = [rotMat\[GLMP.subunit{whichsub}.specialcases.highlum.Lcc GLMP.subunit{whichsub}.specialcases.highlum.Mcc]']';
rotL = GLMP.subunit{whichsub}.specialcases.highlum.theta == angs(idx)...
    | softEq(GLMP.subunit{whichsub}.specialcases.highlum.theta,angs(idx)+pi)...
    | softEq(GLMP.subunit{whichsub}.specialcases.highlum.theta,angs(idx)-pi);
projs = tempRotPts(rotL,1);
pts = min(projs):.001:max(projs);
fitPts = ComputeNakaRushtonJPW(params(idx,:),pts,'asymmetric');
topnsp = max(GLMP.subunit{whichsub}.specialcases.highlum.nspikes);

% Plot GOF selection
axes(disppanel.axes.GOF)
cla; hold on;
h0 = plot(abs(angs/pi*180),GOF,'ko-');
h1 = plot(abs(angs(bestIdx)/pi*180),GOF(bestIdx),'r*');
disppanel.surfsel = plot(abs(angs(idx)/pi*180),GOF(idx),'c*');
set(h0,'ButtonDownFcn',@SurfaceSelectionPred1Ax);
set(h1,'ButtonDownFcn',@SurfaceSelectionPred1Ax);

% Plot selected surface
axes(disppanel.axes.surf3d)
cla; hold on; grid on;
surf(x,y,surface)
plot3(GLMP.subunit{whichsub}.specialcases.highlum.Lcc,GLMP.subunit{whichsub}.specialcases.highlum.Mcc,GLMP.subunit{whichsub}.specialcases.highlum.nspikes,'k*')
plot3(projline(:,1),projline(:,2),[topnsp 0 topnsp]','m')
set(gca,'XLim',[-axlim axlim],'XTickMode','auto','YLim',[-axlim axlim],'YTickMode','auto')
xlabel('L Cone Contrast')
ylabel('M Cone Contrast')
zlabel('Number of Spikes')

% Plot selected CRF
axes(disppanel.axes.proj2d); cla; hold on; grid on;
plot(projs,GLMP.subunit{whichsub}.specialcases.highlum.nspikes(rotL),'k*')
plot(pts,fitPts,'r--')

%Save variables
set(surffig.conpanel,'UserData',conpanel);
set(surffig.disppanel,'UserData',disppanel);
set(gcf,'UserData',surffig);

end


function Fit2AxOrthog

% Load Figure Variables
surffig = get(gcf,'UserData');
conpanel = get(surffig.conpanel,'UserData');
disppanel = get(surffig.disppanel,'UserData');


% Find the selected subunit
if get(conpanel.subbuttons,'SelectedObject') == conpanel.sub1
    p = 1;
elseif get(conpanel.subbuttons,'SelectedObject') == conpanel.sub2
    p = 2;
elseif get(conpanel.subbuttons,'SelectedObject') == conpanel.sub3
    p = 3;
end

keyboard

end
