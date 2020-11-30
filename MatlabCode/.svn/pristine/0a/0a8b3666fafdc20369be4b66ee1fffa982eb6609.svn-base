function [] = GLMSGUI_Surface_CSsym(~,~,varargin)
% This code fits surfaces to GLMP datasets.

SetUpFig()


end


%%% Interactive functions %%%

function reanalall(~,~)
global ConSecSymFit GLMP library

% Load variables
surffig = get(gcf,'UserData');
conpanel = get(surffig.conpanel,'UserData');
fitspanel = get(surffig.fitspanel,'UserData');
parpanel = get(surffig.parpanel,'UserData');


for n = 1:numel(ConSecSymFit.glmps)

    GLMP = ConSecSymFit.glmps{n};
    CalcNormLogLLBounds()
    LN_HR()
    LN_FR()
    oLN_HR()
    oLN_FR()
    
    fitspanel = get(surffig.fitspanel,'UserData');
    ConSecSymFit.params{n} = fitspanel.params;
    ConSecSymFit.LL{n} = fitspanel.LL;
    
    save([library 'ConSecSymFit'],'ConSecSymFit')
    
end
    


% Save variables
set(surffig.conpanel,'UserData',conpanel);
set(surffig.fitspanel,'UserData',fitspanel);
set(surffig.parpanel,'UserData',parpanel);
set(gcf,'UserData',surffig);


end

function startanalysis(~,~)
global GLMP
% This is the heart of the analysis. Here, we will route the analysis to
% different subfunctions.

% Load variables
surffig = get(gcf,'UserData');
conpanel = get(surffig.conpanel,'UserData');
fitspanel = get(surffig.fitspanel,'UserData');
parpanel = get(surffig.parpanel,'UserData');

% idx = 31;
% GLMP = GLMSPopData{idx+1,11};

disp('')
disp(['Beginning Surface Analysis on ' GLMP.datafile...
    '(sub # ' num2str(conpanel.subselect) ')']);

% Run through fits and compare
CalcNormLogLLBounds()
LN_HR()
LN_FR()
oLN_HR()
oLN_FR()




disp('Surface Analysis Completed.')
    
% Save variables
set(surffig.conpanel,'UserData',conpanel);
set(surffig.fitspanel,'UserData',fitspanel);
set(surffig.parpanel,'UserData',parpanel);
set(gcf,'UserData',surffig);

end


%%% Analyses %%%

function LN_HR()
global GLMP
disp('Fitting a half-rectified LN model to the data...')

% Load Figure Variables
surffig = get(gcf,'UserData');
conpanel = get(surffig.conpanel,'UserData');
fitspanel = get(surffig.fitspanel,'UserData');
parpanel = get(surffig.parpanel,'UserData');

% Load data parameters
sub = conpanel.subselect;
Lcc = cat(1,GLMP.subunit{sub}.Lcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
Mcc = cat(1,GLMP.subunit{sub}.Mcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
nsp = cat(1,GLMP.subunit{sub}.nspikes,GLMP.subunit{sub}.blnspikes);

% Call bounds
vub = fitspanel.ub.LN_hr;
vlb = fitspanel.lb.LN_hr;

% Set initial guesses and options for 1D fit
angs = unique(GLMP.subunit{sub}.theta);
Aguess = max(GLMP.subunit{sub}.nspikes)*.8;
expguess = 3;
blguess = mean(GLMP.subunit{sub}.blnspikes);
kappaguess = regress(GLMP.subunit{sub}.varnspikes-GLMP.subunit{sub}.meannspikes,...
    GLMP.subunit{sub}.meannspikes.^2);% Variance for negative binomial fit
if kappaguess < vlb(end)
    kappaguess = vlb(end);
end

% Rotate through axes for which we have data and fit each one. This serves
% as the initial guess for the larger surface.
nLL = Inf;
for rot = 1:size(angs)
    
    % New initial guess
    [x,y] = pol2cart(angs(rot),1);
    tempx = [Lcc Mcc] * [x y]';
    sigguess = max(tempx)/2;
    paramsGuess = [Aguess sigguess 0 0 expguess blguess angs(rot) kappaguess];
    [f1,fval] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
        fitspanel.a,fitspanel.b,[],[],vlb,vub,[],fitspanel.options,[Lcc Mcc],nsp,...
        fitspanel.surftype,fitspanel.errortype);
        
    % If better than previous best fit, replace parmeter values
    if fval < nLL
        params = f1;
        nLL = fval;
    end
end


% Display surface and pts
cmap = repmat(linspace(.7,0,32)',1,3);
colormap(cmap)
x = linspace(-max(GLMP.subunit{sub}.rho),max(GLMP.subunit{sub}.rho),50);
[xx,yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(params,[xx(:) yy(:)],fitspanel.surftype);
surface = reshape(surface,size(xx));
axes(fitspanel.axes.LN_HR); cla; hold on; grid on;
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(GLMP.subunit{sub}.meannspikes);
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
h = contour(xx,yy,surface);
set(gca,'view',[0 90])

% Fill in parameter values
set(parpanel.paramvals.LN_HR.A,'string',params(1));
set(parpanel.paramvals.LN_HR.sig1,'string',params(2));
set(parpanel.paramvals.LN_HR.sig2,'string',params(3)*params(2));
set(parpanel.paramvals.LN_HR.osig,'string',1/params(4));
set(parpanel.paramvals.LN_HR.exp,'string',params(5));
set(parpanel.paramvals.LN_HR.bl,'string',params(6));
set(parpanel.paramvals.LN_HR.rot,'string',params(7)/pi*180);
set(parpanel.paramvals.LN_HR.kappa,'string',params(8));
set(parpanel.paramvals.LN_HR.LL,'string',-nLL);

drawnow;

fitspanel.params.LN_HR = params;
fitspanel.LL.LN_HR = -nLL;

% Save variables
set(surffig.conpanel,'UserData',conpanel);
set(surffig.fitspanel,'UserData',fitspanel);
set(surffig.parpanel,'UserData',parpanel);
set(gcf,'UserData',surffig);


end

function LN_FR()
global GLMP
disp('Fitting a fully rectified LN model to the data...')

% Load Figure Variables
surffig = get(gcf,'UserData');
conpanel = get(surffig.conpanel,'UserData');
fitspanel = get(surffig.fitspanel,'UserData');
parpanel = get(surffig.parpanel,'UserData');

% Load data parameters
sub = conpanel.subselect;
Lcc = cat(1,GLMP.subunit{sub}.Lcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
Mcc = cat(1,GLMP.subunit{sub}.Mcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
nsp = cat(1,GLMP.subunit{sub}.nspikes,GLMP.subunit{sub}.blnspikes);

% Call bounds
vub = fitspanel.ub.LN_fr;
vlb = fitspanel.lb.LN_fr;

% Set initial guesses and options for 1D fit
angs = unique(GLMP.subunit{sub}.theta);
Aguess = max(GLMP.subunit{sub}.nspikes)*.8;
expguess = 3;
blguess = mean(GLMP.subunit{sub}.blnspikes);
kappaguess = regress(GLMP.subunit{sub}.varnspikes-GLMP.subunit{sub}.meannspikes,...
    GLMP.subunit{sub}.meannspikes.^2);% Variance for negative binomial fit
if kappaguess < vlb(end)
    kappaguess = vlb(end);
end

% Rotate through axes for which we have data and fit each one. This serves
% as the initial guess for the larger surface.
nLL = Inf;
for rot = 1:size(angs)
    
    % New initial guess
    [x,y] = pol2cart(angs(rot),1);
    tempx = [Lcc Mcc] * [x y]';
    sigguess = max(tempx)/2;
    paramsGuess = [Aguess sigguess 1 0 expguess blguess angs(rot) kappaguess];
    [f1,fval] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
        fitspanel.a,fitspanel.b,[],[],vlb,vub,[],fitspanel.options,[Lcc Mcc],nsp,...
        fitspanel.surftype,fitspanel.errortype);
        
    % If better than previous best fit, replace parmeter values
    if fval < nLL
        params = f1;
        nLL = fval;
    end
end

% Display surface and pts
cmap = repmat(linspace(.7,0,32)',1,3);
colormap(cmap)
x = linspace(-max(GLMP.subunit{sub}.rho),max(GLMP.subunit{sub}.rho),50);
[xx,yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(params,[xx(:) yy(:)],fitspanel.surftype);
surface = reshape(surface,size(xx));
axes(fitspanel.axes.LN_FR); cla; hold on; grid on;
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(GLMP.subunit{sub}.meannspikes);
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
h = contour(xx,yy,surface);
%h = surf(xx,yy,surface);
set(gca,'view',[0 90])

% Fill in parameter values
set(parpanel.paramvals.LN_FR.A,'string',params(1));
set(parpanel.paramvals.LN_FR.sig1,'string',params(2));
set(parpanel.paramvals.LN_FR.sig2,'string',params(3)*params(2));
set(parpanel.paramvals.LN_FR.osig,'string',1/params(4));
set(parpanel.paramvals.LN_FR.exp,'string',params(5));
set(parpanel.paramvals.LN_FR.bl,'string',params(6));
set(parpanel.paramvals.LN_FR.rot,'string',params(7)/pi*180);
set(parpanel.paramvals.LN_FR.kappa,'string',params(8));
set(parpanel.paramvals.LN_FR.LL,'string',-nLL);

drawnow;

fitspanel.params.LN_FR = params;
fitspanel.LL.LN_FR = -nLL;

% Save variables
set(surffig.conpanel,'UserData',conpanel);
set(surffig.fitspanel,'UserData',fitspanel);
set(surffig.parpanel,'UserData',parpanel);
set(gcf,'UserData',surffig);


end

function oLN_HR()
global GLMP
disp('Fitting a fully rectified oLN model to the data...')

% Load Figure Variables
surffig = get(gcf,'UserData');
conpanel = get(surffig.conpanel,'UserData');
fitspanel = get(surffig.fitspanel,'UserData');
parpanel = get(surffig.parpanel,'UserData');

% Load data parameters
sub = conpanel.subselect;
Lcc = cat(1,GLMP.subunit{sub}.Lcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
Mcc = cat(1,GLMP.subunit{sub}.Mcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
nsp = cat(1,GLMP.subunit{sub}.nspikes,GLMP.subunit{sub}.blnspikes);

% Call bounds
vub = fitspanel.ub.oLN_hr;
vlb = fitspanel.lb.oLN_hr;

% Set initial guesses and options for 1D fit
angs = unique(GLMP.subunit{sub}.theta);
Aguess = max(GLMP.subunit{sub}.nspikes)*.8;
expguess = 3;
blguess = mean(GLMP.subunit{sub}.blnspikes);
kappaguess = regress(GLMP.subunit{sub}.varnspikes-GLMP.subunit{sub}.meannspikes,...
    GLMP.subunit{sub}.meannspikes.^2);% Variance for negative binomial fit
if kappaguess < vlb(end)
    kappaguess = vlb(end);
end

% Rotate through axes for which we have data and fit each one. This serves
% as the initial guess for the larger surface.
nLL = Inf;
for rot = 1:size(angs)
    
    % New initial guess
    [x,y] = pol2cart(angs(rot),1);
    tempx = [Lcc Mcc] * [x y]';
    sigguess = max(tempx)/2;
    [x,y] = pol2cart(angs(rot)+pi/2,1);
    tempx = [Lcc Mcc] * [x y]';
    osigguess = max(tempx)/2;
    paramsGuess = [Aguess sigguess 0 osigguess expguess blguess angs(rot) kappaguess];
    [f1,fval] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
        fitspanel.a,fitspanel.b,[],[],vlb,vub,[],fitspanel.options,[Lcc Mcc],nsp,...
        fitspanel.surftype,fitspanel.errortype);
        
    % If better than previous best fit, replace parmeter values
    if fval < nLL
        params = f1;
        nLL = fval;
    end
end

% Display surface and pts
cmap = repmat(linspace(.7,0,32)',1,3);
colormap(cmap)
x = linspace(-max(GLMP.subunit{sub}.rho),max(GLMP.subunit{sub}.rho),50);
[xx,yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(params,[xx(:) yy(:)],fitspanel.surftype);
surface = reshape(surface,size(xx));
axes(fitspanel.axes.oLN_HR); cla; hold on; grid on;
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(GLMP.subunit{sub}.meannspikes);
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
h = contour(xx,yy,surface);
%h = surf(xx,yy,surface);
set(gca,'view',[0 90])

% Fill in parameter values
set(parpanel.paramvals.oLN_HR.A,'string',params(1));
set(parpanel.paramvals.oLN_HR.sig1,'string',params(2));
set(parpanel.paramvals.oLN_HR.sig2,'string',params(3)*params(2));
set(parpanel.paramvals.oLN_HR.osig,'string',1/params(4));
set(parpanel.paramvals.oLN_HR.exp,'string',params(5));
set(parpanel.paramvals.oLN_HR.bl,'string',params(6));
set(parpanel.paramvals.oLN_HR.rot,'string',params(7)/pi*180);
set(parpanel.paramvals.oLN_HR.kappa,'string',params(8));
set(parpanel.paramvals.oLN_HR.LL,'string',-nLL);

drawnow;

fitspanel.params.oLN_HR = params;
fitspanel.LL.oLN_HR = -nLL;

% Save variables
set(surffig.conpanel,'UserData',conpanel);
set(surffig.fitspanel,'UserData',fitspanel);
set(surffig.parpanel,'UserData',parpanel);
set(gcf,'UserData',surffig);


end

function oLN_FR()
global GLMP
disp('Fitting a fully rectified oLN model to the data...')

% Load Figure Variables
surffig = get(gcf,'UserData');
conpanel = get(surffig.conpanel,'UserData');
fitspanel = get(surffig.fitspanel,'UserData');
parpanel = get(surffig.parpanel,'UserData');

% Load data parameters
sub = conpanel.subselect;
Lcc = cat(1,GLMP.subunit{sub}.Lcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
Mcc = cat(1,GLMP.subunit{sub}.Mcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
nsp = cat(1,GLMP.subunit{sub}.nspikes,GLMP.subunit{sub}.blnspikes);

% Call bounds
vub = fitspanel.ub.oLN_fr;
vlb = fitspanel.lb.oLN_fr;

% Set initial guesses and options for 1D fit
angs = unique(GLMP.subunit{sub}.theta);
Aguess = max(GLMP.subunit{sub}.nspikes)*.8;
expguess = 3;
blguess = mean(GLMP.subunit{sub}.blnspikes);
kappaguess = regress(GLMP.subunit{sub}.varnspikes-GLMP.subunit{sub}.meannspikes,...
    GLMP.subunit{sub}.meannspikes.^2);% Variance for negative binomial fit
if kappaguess < vlb(end)
    kappaguess = vlb(end);
end

% Rotate through axes for which we have data and fit each one. This serves
% as the initial guess for the larger surface.
nLL = Inf;
for rot = 1:size(angs)
    
    % New initial guess
    [x,y] = pol2cart(angs(rot),1);
    tempx = [Lcc Mcc] * [x y]';
    sigguess = max(tempx)/2;
    [x,y] = pol2cart(angs(rot)+pi/2,1);
    tempx = [Lcc Mcc] * [x y]';
    osigguess = max(tempx)/2;
    paramsGuess = [Aguess sigguess 1 osigguess expguess blguess angs(rot) kappaguess];
    [f1,fval] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
        fitspanel.a,fitspanel.b,[],[],vlb,vub,[],fitspanel.options,[Lcc Mcc],nsp,...
        fitspanel.surftype,fitspanel.errortype);
        
    % If better than previous best fit, replace parmeter values
    if fval < nLL
        params = f1;
        nLL = fval;
    end
end

% Display surface and pts
cmap = repmat(linspace(.7,0,32)',1,3);
colormap(cmap)
x = linspace(-max(GLMP.subunit{sub}.rho),max(GLMP.subunit{sub}.rho),50);
[xx,yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(params,[xx(:) yy(:)],fitspanel.surftype);
surface = reshape(surface,size(xx));
axes(fitspanel.axes.oLN_FR); cla; hold on; grid on;
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(GLMP.subunit{sub}.meannspikes);
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
h = contour(xx,yy,surface);
%h = surf(xx,yy,surface);
set(gca,'view',[0 90])

% Fill in parameter values
set(parpanel.paramvals.oLN_FR.A,'string',params(1));
set(parpanel.paramvals.oLN_FR.sig1,'string',params(2));
set(parpanel.paramvals.oLN_FR.sig2,'string',params(3)*params(2));
set(parpanel.paramvals.oLN_FR.osig,'string',1/params(4));
set(parpanel.paramvals.oLN_FR.exp,'string',params(5));
set(parpanel.paramvals.oLN_FR.bl,'string',params(6));
set(parpanel.paramvals.oLN_FR.rot,'string',params(7)/pi*180);
set(parpanel.paramvals.oLN_FR.kappa,'string',params(8));
set(parpanel.paramvals.oLN_FR.LL,'string',-nLL);

drawnow;

fitspanel.params.oLN_FR = params;
fitspanel.LL.oLN_FR = -nLL;

% Save variables
set(surffig.conpanel,'UserData',conpanel);
set(surffig.fitspanel,'UserData',fitspanel);
set(surffig.parpanel,'UserData',parpanel);
set(gcf,'UserData',surffig);


end

function CalcNormLogLLBounds()

global GLMP
% Normalized log likelhood between mean and exact fit

disp('Calculating Normalized Log Likelihood...')

% Load Figure Variables
surffig = get(gcf,'UserData');
conpanel = get(surffig.conpanel,'UserData');
fitspanel = get(surffig.fitspanel,'UserData');
parpanel = get(surffig.parpanel,'UserData');

% Load data variables
sub = conpanel.subselect;
nsp = cat(1,GLMP.subunit{sub}.nspikes,GLMP.subunit{sub}.blnspikes);
kappaguess = .5;

% Find minimum LL (using Negative Binomial)
mu = repmat(mean(nsp),size(nsp));
x = nsp;
[kappa,f] = fmincon('minimizeKappa',kappaguess,[],[],[],[],0,[],[],fitspanel.options,mu,x);
fitspanel.normLL.lbLL = -f;
fitspanel.normLL.lbkappa = kappa;

% Find maximum LL (using negative binomial)
kappaguess = .5;
nsp = [];
prednsp = [];
for t = 1:numel(GLMP.subunit{sub}.uniqueIdx)
    nsp = cat(1,nsp,GLMP.subunit{sub}.nspikes(GLMP.subunit{sub}.uniqueIdx{t}));
    prednsp = cat(1,prednsp,repmat(GLMP.subunit{sub}.meannspikes(t),numel(GLMP.subunit{sub}.uniqueIdx{t}),1));
end
meanbl = repmat(mean(GLMP.subunit{sub}.blnspikes),numel(GLMP.subunit{sub}.blnspikes),1);
prednsp = cat(1,prednsp,meanbl);
nsp = cat(1,nsp,GLMP.subunit{sub}.blnspikes);
mu = prednsp;
x = nsp;
[kappa,f] = fmincon('minimizeKappa',kappaguess,[],[],[],[],0,[],[],fitspanel.options,mu,x);
fitspanel.normLL.ubLL = -f;
fitspanel.normLL.ubkappa = kappa;

% Calculate normalized LL values
% ub = fitspanel.normLL.ubLL;
% lb = fitspanel.normLL.lbLL;
% fitspanel.fitcomps.oneD.normLL = (-fitspanel.fitcomps.oneD.nLLs - lb)/(ub-lb);
% fitspanel.fitcomps.twoD.normLL = (-fitspanel.fitcomps.twoD.nLLs - lb)/(ub-lb);
% oneDLL = fitspanel.LL.oneD;
% twoDLL = fitspanel.LL.twoD;
% fitspanel.normLL.oneDnormLL = (oneDLL-lb)/(ub-lb);
% fitspanel.normLL.twoDnormLL = (twoDLL-lb)/(ub-lb);

% Fill in fields
% set(conpanel.paramvals.oneD.normLL,'string',round(fitspanel.normLL.oneDnormLL,2))
% set(conpanel.paramvals.twoD.normLL,'string',round(fitspanel.normLL.twoDnormLL,2))
% set(conpanel.paramvals.normLL,'string',round(fitspanel.normLL.twoDnormLL-fitspanel.normLL.oneDnormLL,4))

% Save variables
set(surffig.conpanel,'UserData',conpanel);
set(surffig.fitspanel,'UserData',fitspanel);
set(surffig.parpanel,'UserData',parpanel);
set(gcf,'UserData',surffig);


end


%%% Setup %%%

function SetUpFig()
global GLMP ConSecSymFit library

% Set up figure
figure(60); clf;
set(gcf,'units','pixels','pos',[100 150 1200 700],'NumberTitle','off',...
    'Name','Surface Analyses');

% Grab saved population data 
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
load([library 'ConSecSymFit.mat'])

% Set up panels
surffig = get(gcf,'UserData');
surffig.fitspanel = uipanel('Pos',[.01 .5 .98 .49],'Parent',gcf,'title','Fits to Data');
surffig.conpanel = uipanel('pos',[.01 .01 .98 .48],'parent',gcf);
surffig.parpanel = uipanel('parent',surffig.conpanel,'units','normalized',...
    'pos',[.01 .01 .78 .98],'title','Parameters');

% Set up control panel
conpanel.startanalysis = uicontrol('style','pushbutton',...
    'parent',surffig.conpanel,'units','normalized','pos',[.8 .05 .09 .1],...
    'string','Analyze Cell','fontsize',12,'callback',@startanalysis,...
    'backgroundcolor',[.2 1 0]);
conpanel.analyzeall  = uicontrol('style','pushbutton',...
    'parent',surffig.conpanel,'units','normalized','pos',[.9 .05 .09 .1],...
    'string','Reanalyze All','fontsize',12,'callback',@reanalall,...
    'backgroundcolor',[.2 1 0]);

% Table
data = ConSecSymFit.datafiles;
% data(:,2) = ConSecSymFit(1:end,strcmp(datatypes,'Subunit'));
% data(:,3) = ConSecSymFit(1:end,strcmp(datatypes,'normLLdiff'));
%data(:,4) = num2cell(conpanel.conf);

% Repackage data and name columns
colname = {'Datafile'};
colformat = cell(1,size(data,2));
colformat(:) = {'char'};

% Display table 
conpanel.table = uitable('parent',surffig.conpanel,...
    'units','normalized','pos',[.8 .2  .19 .75],...
    'data',data,'columnname',colname,...
    'columnformat',colformat,...
    'CellSelectionCallback',@cellselect);

% Parameter labels
parpanel.paramslabels.A = uicontrol('style','text','parent',surffig.parpanel,...
    'units','normalized','pos',[.1 .9 .08 .08],'fontsize',12,'string','A =');
parpanel.paramslabels.sig1 = uicontrol('style','text','parent',surffig.parpanel,...
    'units','normalized','pos',[.2 .9 .08 .08],'fontsize',12,'string','Sig1 =');
parpanel.paramslabels.sig2 = uicontrol('style','text','parent',surffig.parpanel,...
    'units','normalized','pos',[.3 .9 .08 .08],'fontsize',12,'string','Sig2 =');
parpanel.paramslabels.osig = uicontrol('style','text','parent',surffig.parpanel,...
    'units','normalized','pos',[.4 .9 .08 .08],'fontsize',12,'string','oSig =');
parpanel.paramslabels.exp = uicontrol('style','text','parent',surffig.parpanel,...
    'units','normalized','pos',[.5 .9 .08 .08],'fontsize',12,'string','Exponent =');
parpanel.paramslabels.bl = uicontrol('style','text','parent',surffig.parpanel,...
    'units','normalized','pos',[.6 .9 .08 .08],'fontsize',12,'string','Baseline =');
parpanel.paramslabels.rot = uicontrol('style','text','parent',surffig.parpanel,...
    'units','normalized','pos',[.7 .9 .08 .08],'fontsize',12,'string','Rotation =');
parpanel.paramslabels.kappa = uicontrol(surffig.parpanel,'style','text',...
    'units','normalized','pos',[.8 .9 .08 .08],'fontsize',12,'string','Kappa =');
parpanel.paramslabels.nLL = uicontrol(surffig.parpanel,'style','text',...
    'units','normalized','pos',[.9 .9 .08 .08],'fontsize',12,'string','norm LL =');


%%% Parameter values
% LN half-rectified
parpanel.paramvals.LN_HR.A = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.1 .7 .08 .08],'fontsize',12);
parpanel.paramvals.LN_HR.sig1 = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.2 .7 .08 .08],'fontsize',12);
parpanel.paramvals.LN_HR.sig2 = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.3 .7 .08 .08],'fontsize',12);
parpanel.paramvals.LN_HR.osig = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.4 .7 .08 .08],'fontsize',12);
parpanel.paramvals.LN_HR.exp = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.5 .7 .08 .08],'fontsize',12);
parpanel.paramvals.LN_HR.bl = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.6 .7 .08 .08],'fontsize',12);
parpanel.paramvals.LN_HR.rot = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.7 .7 .08 .08],'fontsize',12);
parpanel.paramvals.LN_HR.kappa = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.8 .7 .08 .08],'fontsize',12);
parpanel.paramvals.LN_HR.LL = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.9 .7 .08 .08],'fontsize',12);

% LN fully-rectified
parpanel.paramvals.LN_FR.A = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.1 .5 .08 .08],'fontsize',12);
parpanel.paramvals.LN_FR.sig1 = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.2 .5 .08 .08],'fontsize',12);
parpanel.paramvals.LN_FR.sig2 = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.3 .5 .08 .08],'fontsize',12);
parpanel.paramvals.LN_FR.osig = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.4 .5 .08 .08],'fontsize',12);
parpanel.paramvals.LN_FR.exp = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.5 .5 .08 .08],'fontsize',12);
parpanel.paramvals.LN_FR.bl = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.6 .5 .08 .08],'fontsize',12);
parpanel.paramvals.LN_FR.rot = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.7 .5 .08 .08],'fontsize',12);
parpanel.paramvals.LN_FR.kappa = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.8 .5 .08 .08],'fontsize',12);
parpanel.paramvals.LN_FR.LL = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.9 .5 .08 .08],'fontsize',12);


% oLN fully-rectified
parpanel.paramvals.oLN_HR.A = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.1 .3 .08 .08],'fontsize',12);
parpanel.paramvals.oLN_HR.sig1 = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.2 .3 .08 .08],'fontsize',12);
parpanel.paramvals.oLN_HR.sig2 = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.3 .3 .08 .08],'fontsize',12);
parpanel.paramvals.oLN_HR.osig = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.4 .3 .08 .08],'fontsize',12);
parpanel.paramvals.oLN_HR.exp = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.5 .3 .08 .08],'fontsize',12);
parpanel.paramvals.oLN_HR.bl = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.6 .3 .08 .08],'fontsize',12);
parpanel.paramvals.oLN_HR.rot = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.7 .3 .08 .08],'fontsize',12);
parpanel.paramvals.oLN_HR.kappa = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.8 .3 .08 .08],'fontsize',12);
parpanel.paramvals.oLN_HR.LL = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.9 .3 .08 .08],'fontsize',12);


% oLN fully-rectified
parpanel.paramvals.oLN_FR.A = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.1 .1 .08 .08],'fontsize',12);
parpanel.paramvals.oLN_FR.sig1 = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.2 .1 .08 .08],'fontsize',12);
parpanel.paramvals.oLN_FR.sig2 = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.3 .1 .08 .08],'fontsize',12);
parpanel.paramvals.oLN_FR.osig = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.4 .1 .08 .08],'fontsize',12);
parpanel.paramvals.oLN_FR.exp = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.5 .1 .08 .08],'fontsize',12);
parpanel.paramvals.oLN_FR.bl = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.6 .1 .08 .08],'fontsize',12);
parpanel.paramvals.oLN_FR.rot = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.7 .1 .08 .08],'fontsize',12);
parpanel.paramvals.oLN_FR.kappa = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.8 .1 .08 .08],'fontsize',12);
parpanel.paramvals.oLN_FR.LL = uicontrol('style','edit','parent',surffig.parpanel,...
    'units','normalized','pos',[.9 .1 .08 .08],'fontsize',12);


% Set up fits panel (for displaying surfaces)
fitspanel.axes.LN_HR = axes('parent',surffig.fitspanel,'units','normalized',...
    'pos',[.04 .15 .2 .75],'box','on','cameraposition',[5 -5 4]);
title('LN Half-Rectified'); xlabel('Lcc'); ylabel('Mcc'); zlabel('Response'); grid on; hold on;
fitspanel.axes.LN_FR = axes('parent',surffig.fitspanel,'units','normalized',...
    'pos',[.28 .15 .2 .75],'box','on','cameraposition',[5 -5 4]);
title('LN Fully-Rectified'); xlabel('Lcc'); ylabel('Mcc'); zlabel('Response'); grid on; hold on;
fitspanel.axes.oLN_HR = axes('parent',surffig.fitspanel,'units','normalized',...
    'pos',[.52 .15 .2 .75],'box','on','cameraposition',[5 -5 4]);
title('oLN Fully-Rectified'); xlabel('Lcc'); ylabel('Mcc'); zlabel('Response'); grid on; hold on;
fitspanel.axes.oLN_FR = axes('parent',surffig.fitspanel,'units','normalized',...
    'pos',[.76 .15 .2 .75],'box','on','cameraposition',[5 -5 4]);
title('oLN Fully-Rectified'); xlabel('Lcc'); ylabel('Mcc'); zlabel('Response'); grid on; hold on;

% Default subunit is 1
conpanel.subselect = 1;

% Set up some fitting specifics
fitspanel.options = optimset('Algorithm','interior-point','MaxFunEvals',5000,...
    'MaxIter',5000,'FinDiffType','central','Hessian','bfgs','display','off',...
    'TolFun',0,'TolCon',0,'FunValCheck','on','AlwaysHonorConstraints','bounds');

% Bounds for oLN half-rectified
fitspanel.ub.oLN_hr = [max(GLMP.nspikes)  Inf 0  300 10 max(GLMP.nspikes)  pi 5];
fitspanel.lb.oLN_hr = [min(GLMP.nspikes) .003 0 -300  1              .001 -pi 0];

% Bounds for oLN fully-rectified
fitspanel.ub.oLN_fr = fitspanel.ub.oLN_hr;
fitspanel.lb.oLN_fr = fitspanel.lb.oLN_hr;
fitspanel.ub.oLN_fr(3) = 1;
fitspanel.lb.oLN_fr(3) = 1;

% Bounds for LN fully-rectified
fitspanel.ub.LN_fr = fitspanel.ub.oLN_fr;
fitspanel.lb.LN_fr = fitspanel.lb.oLN_fr;
fitspanel.ub.LN_fr(4) = 0;
fitspanel.lb.LN_fr(4) = 0;

% Bounds for LN half-rectified
fitspanel.ub.LN_hr = fitspanel.ub.oLN_hr;
fitspanel.lb.LN_hr = fitspanel.lb.oLN_hr;
fitspanel.ub.LN_hr(4) = 0;
fitspanel.lb.LN_hr(4) = 0;

fitspanel.surftype = 'conicsection_sym';
fitspanel.errortype = 'NegativeBinomial';
fitspanel.thresh = .05; % For comparing types of 1D fits
fitspanel.params.LN_HR = [];
fitspanel.params.LN_FR = [];
fitspanel.params.oLN_HR = [];
fitspanel.params.oLN_FR = [];
fitspanel.LL.LN_HR = [];
fitspanel.LL.LN_FR = [];
fitspanel.LL.oLN_HR = [];
fitspanel.LL.oLN_FR = [];
fitspanel.a = [-1 0 0 0 0 1 0 0];
fitspanel.b = zeros(size(fitspanel.a,1),1);

%Save variables
set(surffig.conpanel,'UserData',conpanel);
set(surffig.fitspanel,'UserData',fitspanel);
set(surffig.parpanel,'UserData',parpanel);
set(gcf,'UserData',surffig);


end

function cellSelect(a,b)






%Save variables
set(surffig.conpanel,'UserData',conpanel);
set(surffig.fitspanel,'UserData',fitspanel);
set(surffig.parpanel,'UserData',parpanel);
set(gcf,'UserData',surffig);

end