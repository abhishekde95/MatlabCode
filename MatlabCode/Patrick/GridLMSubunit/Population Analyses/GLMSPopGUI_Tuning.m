function GLMSPopGUI_Tuning(varargin)

SetUpFig()
UnpackPopulationData()
resetAll()


if ~isempty(varargin)
    for n = 1:numel(varargin)
        eval(varargin{n})
    end
end

end


%%% Analyses %%%

function ReanalAll(~,~,~)
global GLMSPopData

% Load figure and pop variables
TunFig = get(150,'userdata');
conpanel = get(TunFig.conpanel,'userdata');

% Rotate through selected index
for n = 1:(size(GLMSPopData,1)-1)
    conpanel.selectedidx = n;
    set(TunFig.conpanel,'userdata',conpanel)
    set(150,'userdata',TunFig)
    Reanal()
end

disp('Finished reanalyzing population surface fits.')

end

function Reanal(~,~,~)
global GLMSPopData GLMP DN

% Load figure and pop variables
TunFig = get(150,'userdata');
conpanel = get(TunFig.conpanel,'userdata');

% Which dataset
if isempty(conpanel.selectedidx)
    disp('Must Select Cell First')
    return
end

% Reload nex file, unpack, analyze, and save
datatypes = GLMSPopData(1,:);
if ismac
    library1 = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/nex files/';
    library2 = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library1 = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\nex files\';
    library2 = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
    %library1 = 'C:\Users\Patty\Dropbox\Patrick\GLMS Data\nex files\';
    %library2 = 'C:\Users\Patty\Dropbox\Patrick\GLMS Data\';
end
datafile = conpanel.table.Data(conpanel.selectedidx,1);
sub = GLMSPopData{conpanel.selectedidx+1,strcmp(datatypes,'Subunit')};
disp(['Analyzing datafile No. ' num2str(conpanel.selectedidx) ': ' datafile{1} ' sub # ' num2str(sub)])

% Unpack data and begin surface analysis
popidx =  conpanel.selectedidx+1;

% Unpack raw data, reorganize, and display overview
%rawdata = nex2stro([library1 char(datafile) '.nex']);
%[GLMP,DN] = OrganizeRawGLMSData(rawdata);
%GLMSPopData{popidx,strcmp(datatypes,'GLMP')} = GLMP;
%GLMSPopData{popidx,strcmp(datatypes,'DN')} = DN;
%GLMSGUI_Overview();

% To load previously unpacked raw data
GLMP = GLMSPopData{popidx,strcmp(datatypes,'GLMP')};
DN = GLMSPopData{popidx,strcmp(datatypes,'DN')};

% Analyze
GLMSGUI_Surface([],[],'Analyze',sub)

% Get parameter values from figure variables
SurfFig = get(60,'userdata');
fitspanel = get(SurfFig.fitspanel,'userdata');

% Organize parameters into population structure and save
oneD.surftype = fitspanel.surftype;
twoD.surftype = fitspanel.surftype;
oneD.parvals = fitspanel.params.oneD;
twoD.parvals = fitspanel.params.twoD;
oneD.LL = fitspanel.LL.oneD;
twoD.LL = fitspanel.LL.twoD;
oneD.Hessian = fitspanel.Hessian.oneD;
twoD.Hessian = fitspanel.Hessian.twoD;
oneD.normLL = fitspanel.normLL.oneDnormLL;
twoD.normLL = fitspanel.normLL.twoDnormLL;
oneD.rsq = fitspanel.rsq.oneD;
twoD.rsq = fitspanel.rsq.twoD;
surfparams = struct('oneD',oneD,'twoD',twoD);
surfparams.diffnormLL = twoD.normLL - oneD.normLL;
surfparams.LLRTpval = fitspanel.LLRTpval;
surfparams.fitcomps = fitspanel.fitcomps;

Lcc = GLMP.subunit{sub}.uniqueLcc;
Mcc = GLMP.subunit{sub}.uniqueMcc;
pred1d = ComputeNakaRushtonJPW(oneD.parvals,[Lcc Mcc],fitspanel.surftype);
pred2d = ComputeNakaRushtonJPW(twoD.parvals,[Lcc Mcc],fitspanel.surftype);
diff1d2d = diff([pred2d pred1d],[],2);
surfparams.maxdiffnsp = max(abs(diff1d2d));

GLMSPopData{popidx,strcmp(datatypes,'Surface Parameters')} = surfparams;

% Organize into 1D and 2D cells
GLMSPopData{popidx,strcmp(datatypes,'normLLdiff')} = surfparams.diffnormLL;
if surfparams.diffnormLL < conpanel.normLLthresh
    GLMSPopData{popidx,strcmp(datatypes,'nD')} = '1D';
    nd = 'oneD';
else
    GLMSPopData{popidx,strcmp(datatypes,'nD')} = '2D';
    nd = 'twoD';
end
GLMSPopData{popidx,strcmp(datatypes,'Tuning')} = surfparams.(nd).parvals(end-1)/pi*180;

% Save user data
disp('Saving surface data into the population...')
save([library2 'GLMSPopData'],'GLMSPopData')
set(TunFig.conpanel,'userdata',conpanel)
set(150,'userdata',TunFig)
disp('Data saved.')

% Now that new data is in the pop structure, reload whole analysis
UnpackPopulationData()
resetAll()
b.Indices = conpanel.selectedidx;
cellselect([],b)

disp(['Datafile No. ' num2str(conpanel.selectedidx) ': ' datafile{1} ' sub # ' num2str(sub) ' Reanalyzed.'])

end


%%% Interactive %%%

function cellselect(a,b)
global GLMSPopData GLMP

tic
% Load figure and pop variables
TunFig = get(150,'userdata');
conpanel = get(TunFig.conpanel,'userdata');
popanel = get(TunFig.popanel,'userdata');
surfpanel = get(TunFig.surfpanel,'userdata');
cellpanel = get(TunFig.cellpanel,'userdata');

% Set aside the index (+1 for referencing GLMSPopData)
if isempty(b)
    conpanel.selectedidx = a;
elseif strcmp(get(a,'type'),'uitable')
    if isempty(b.Indices)
        return
    elseif conpanel.selectedidx == b.Indices(1)
        return
    end
    conpanel.selectedidx = b.Indices(1);
end
set(conpanel.table,'CellSelectionCallback',[])
popidx = conpanel.selectedidx + 1;

% Unpack variables
datatypes = GLMSPopData(1,:);
surfparams = GLMSPopData{popidx,strcmp(datatypes,'Surface Parameters')};
GLMP = GLMSPopData{popidx,strcmp(datatypes,'GLMP')};
DN = GLMSPopData{popidx,strcmp(datatypes,'DN')};
sub = GLMSPopData{popidx,strcmp(datatypes,'Subunit')};
Lcc = GLMP.subunit{sub}.Lcc;
%Mcc = GLMP.subunit{sub}.Mcc;
%nsp = GLMP.subunit{sub}.nspikes;

% Pull out specifics of selected datafile
if conpanel.diffnormLL(conpanel.selectedidx) < conpanel.normLLthresh
    params1 = conpanel.oneD.params(conpanel.selectedidx,:);
    ndL = conpanel.oneD.L;
    ax = popanel.axes.oneD;    
    if popanel.uicontrols.which1Dsurf.Value == 1
        surfL = conpanel.oneD.L;
    elseif popanel.uicontrols.which1Dsurf.Value == 2
        surfL = conpanel.oneD.unichrom.L;
    elseif popanel.uicontrols.which1Dsurf.Value == 3
        surfL = conpanel.oneD.bichrom.L;
    end
else
    params1 = conpanel.twoD.params(conpanel.selectedidx,:);
    ndL = ~conpanel.oneD.L;
    ax = popanel.axes.twoD;
    if popanel.uicontrols.which2Dsurf.Value == 1
        surfL = ones(size(conpanel.oneD.L));
    elseif popanel.uicontrols.which2Dsurf.Value == 2
        surfL = conpanel.twoD.unichrom.hypL;
    elseif popanel.uicontrols.which2Dsurf.Value == 3
        surfL = conpanel.twoD.bichrom.hypL;
    elseif popanel.uicontrols.which2Dsurf.Value == 4
        surfL = conpanel.twoD.unichrom.eliL;
    elseif popanel.uicontrols.which2Dsurf.Value == 5
        surfL = conpanel.twoD.bichrom.eliL;
    end
end

% Which monkey
if popanel.uicontrols.whichMonk.Value == 1
    monkL = ones(size(conpanel.monk.M));
elseif popanel.uicontrols.whichMonk.Value == 2
    monkL = conpanel.monk.N;
elseif popanel.uicontrols.whichMonk.Value == 3
    monkL = conpanel.monk.M;
end

%%% Plot whitenoise data %%%
filters = GLMSPopData{popidx,strcmp(datatypes,'DN Filters')};
dev = filters.subfilter(sub).deviation;
[~,idx] = max(dev)
stas = filters.norm4filterSTAs(:,idx);
linsta = filters.norm1filterSTA(:,:,:,idx);

axes(cellpanel.pLpM.axes); cla; box on; hold on;
cellpanel.pLpM.STA = image(repmat(stas{1},[1 1 3])); axis tight
set(cellpanel.pLpM.STA,'ButtonDownFcn',@subonoff)

axes(cellpanel.mLmM.axes); cla; box on; hold on;
cellpanel.mLmM.STA = image(repmat(stas{2},[1 1 3])); axis tight
set(cellpanel.mLmM.STA,'ButtonDownFcn',@subonoff)

axes(cellpanel.pLmM.axes); cla; box on; hold on;
cellpanel.pLmM.STA = image(repmat(stas{3},[1 1 3])); axis tight
set(cellpanel.pLmM.STA,'ButtonDownFcn',@subonoff)

axes(cellpanel.mLpM.axes); cla; box on; hold on;
cellpanel.mLpM.STA = image(repmat(stas{4},[1 1 3])); axis tight
set(cellpanel.mLpM.STA,'ButtonDownFcn',@subonoff)

axes(cellpanel.linfil.axes); cla; box on; hold on;
cellpanel.linfil.STA = image(linsta); axis tight;
set(cellpanel.linfil.STA,'ButtonDownFcn',@subonoff)

% Draw subunits
cellpanel.pLpM.subs = [];
cellpanel.mLmM.subs = [];
cellpanel.pLmM.subs = [];
cellpanel.mLpM.subs = [];
cellpanel.linfil.subs = [];
if sub == 1
    col = 'r';
elseif sub == 2
    col = 'g';
else
    col = 'y';
end

offset = ceil(DN.NStixGrid(1)/2);
axes(cellpanel.pLpM.axes); hold on;
cellpanel.pLpM.sub = plot(GLMP.subunit{sub}.gridX{1}+offset,-GLMP.subunit{sub}.gridY{1}+offset,[col 'o']);
axes(cellpanel.mLmM.axes);
cellpanel.mLmM.sub = plot(GLMP.subunit{sub}.gridX{1}+offset,-GLMP.subunit{sub}.gridY{1}+offset,[col 'o']);
axes(cellpanel.pLmM.axes);
cellpanel.pLmM.sub = plot(GLMP.subunit{sub}.gridX{1}+offset,-GLMP.subunit{sub}.gridY{1}+offset,[col 'o']);
axes(cellpanel.mLpM.axes);
cellpanel.mLpM.sub = plot(GLMP.subunit{sub}.gridX{1}+offset,-GLMP.subunit{sub}.gridY{1}+offset,[col 'o']);
axes(cellpanel.linfil.axes);
cellpanel.linfil.sub = plot(GLMP.subunit{sub}.gridX{1}+offset,-GLMP.subunit{sub}.gridY{1}+offset,[col 'o']);


%%% Plot surfaces %%%
PlotPSTH()

% Plot 1D Surface
colormap('cool')
params1d = conpanel.oneD.params(conpanel.selectedidx,:);
x = linspace(-max(GLMP.subunit{sub}.Lcc),max(GLMP.subunit{sub}.Mcc),50);
[xx,yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(params1d,[xx(:) yy(:)],surfparams.oneD.surftype);
surface = reshape(surface,size(xx));
axes(surfpanel.axes.surface1d); cla; hold on; grid on;
ticks = rndofferr(linspace(min(Lcc),max(Lcc),5),2);
set(gca,'XTick',ticks,'YTick',ticks,'xlim',[min(xx(:)) max(xx(:))],'ylim',[min(yy(:)) max(yy(:))]);
p = surf(gca,xx,yy,surface);
set(p,'edgecolor','none')
alpha(.3);
contour(xx,yy,surface,'linewidth',2);
uniqLcc = GLMP.subunit{sub}.uniqueLcc;
uniqMcc = GLMP.subunit{sub}.uniqueMcc;
meannsp = GLMP.subunit{sub}.meannspikes;
for i = 1:numel(uniqLcc)
    mn = meannsp(i)/max(meannsp)*10+3;
    h = plot3(uniqLcc(i),uniqMcc(i),meannsp(i),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',mn,'MarkerEdgeColor','white')
end
xlabel('Lcc');
ylabel('Mcc');
zlabel('# of spikes')
title('1D Fit')

% Plot 2D Surface
params2d = conpanel.twoD.params(conpanel.selectedidx,:);
surface = ComputeNakaRushtonJPW(params2d,[xx(:) yy(:)],surfparams.twoD.surftype);
surface = reshape(surface,size(xx));
axes(surfpanel.axes.surface2d); cla; hold on; grid on;
set(gca,'XTick',ticks,'YTick',ticks,'xlim',[min(xx(:)) max(xx(:))],'ylim',[min(yy(:)) max(yy(:))]);
p = surf(xx,yy,surface);
set(p,'edgecolor','none')
alpha(.3);
contour(xx,yy,surface,'linewidth',2);
for i = 1:numel(uniqLcc)
    mn = meannsp(i)/max(meannsp)*10+3;
    h = plot3(uniqLcc(i),uniqMcc(i),meannsp(i),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',mn,'MarkerEdgeColor','white')
end
xlabel('Lcc');
ylabel('Mcc');
zlabel('# of spikes')
title('2D Fit')

% Link plots! This is cool!
% linkax = [surfpanel.axes.surface1d, surfpanel.axes.surface2d];
% hlink = linkprop(linkax,...
%     {'cameraposition','cameraupvector','cameratarget','cameraviewangle'});
%keyboard

% Fill in 1d parameter values
set(surfpanel.uicontrols.params.oneD.A,'string',round(params1d(1),1));
set(surfpanel.uicontrols.params.oneD.sig1,'string',round(1/params1d(2),3));
set(surfpanel.uicontrols.params.oneD.sig2,'string',round(1/params1d(3),3));
set(surfpanel.uicontrols.params.oneD.exp,'string',round(params1d(5),1));
set(surfpanel.uicontrols.params.oneD.bl,'string',round(params1d(6),2));
set(surfpanel.uicontrols.params.oneD.rot,'string',round(params1d(7)/pi*180,0));
set(surfpanel.uicontrols.params.oneD.kappa,'string',round(params1d(8),1));

% Fill in 2d parameter values
set(surfpanel.uicontrols.params.twoD.A,'string',round(params2d(1),1));
set(surfpanel.uicontrols.params.twoD.sig1,'string',round(1/params2d(2),3));
set(surfpanel.uicontrols.params.twoD.sig2,'string',round(1/params2d(3),3));
set(surfpanel.uicontrols.params.twoD.orthsig,'string',round(1/params2d(4),3));
set(surfpanel.uicontrols.params.twoD.exp,'string',round(params2d(5),1));
set(surfpanel.uicontrols.params.twoD.bl,'string',round(params2d(6),2));
set(surfpanel.uicontrols.params.twoD.rot,'string',round(params2d(7)/pi*180,0));
set(surfpanel.uicontrols.params.twoD.kappa,'string',round(params2d(8),1));

% Fill in Norm LL vals
nLL1d = rndofferr(conpanel.oneD.bichrom.normLLs(conpanel.selectedidx),4);
nLL2d = rndofferr(conpanel.twoD.bichrom.normLLs(conpanel.selectedidx),4);
nLLdiff = rndofferr(nLL2d - nLL1d,4);
set(conpanel.uicontrols.normLL1D,'string',nLL1d);
set(conpanel.uicontrols.normLL2D,'string',nLL2d);
set(conpanel.uicontrols.normLLdiff,'string',nLLdiff);

[conpanel.twoD.bichrom.normLLs(conpanel.selectedidx) conpanel.oneD.bichrom.normLLs(conpanel.selectedidx);...
    conpanel.twoD.unichrom.normLLs(conpanel.selectedidx) conpanel.oneD.unichrom.normLLs(conpanel.selectedidx)]

%%% Replot everything %%%
set(TunFig.conpanel,'userdata',conpanel)
set(gcf,'userdata',TunFig)
resetAll()
TunFig = get(gcf,'userdata');
conpanel = get(TunFig.conpanel,'userdata');


%%% Highlight bin in normLL hist %%%
if ~conpanel.excludeL(conpanel.selectedidx)
    axes(conpanel.axes.normLLhist)
    [~,idx] = min(abs(nLLdiff - conpanel.hist.centers));
    counts = zeros(size(conpanel.hist.counts));
    counts(idx) = conpanel.hist.counts(idx);
    bar(conpanel.hist.centers,counts,'r')
end

%%% Highlight bin in pop tuning %%%
if ~conpanel.excludeL(conpanel.selectedidx)
        
    % Find the collection of cells in same tuning bin
    [tout,~] = rose(conpanel.tuning,popanel.nangs);
    edges = unique(tout-pi);
    lowerTheta = edges(find(sign(params1(end-1)-edges)==1,1,'last'));
    upperTheta = edges(find(sign(params1(end-1)-edges)~=1,1,'first'));
    angL = conpanel.tuning > lowerTheta & conpanel.tuning < upperTheta;
    
    % Find and plot selected angles
    col = [1 0 0];
    L = ndL & angL & monkL & surfL & ~conpanel.excludeL;
    if any(find(L) == conpanel.selectedidx)
        L = ndL & angL & monkL & surfL & conpanel.uniqL & ~conpanel.excludeL;
        axes(ax);
        h = rose(conpanel.tuning(L),popanel.nangs);
        set(h,'Color',col,'LineWidth',2)
    end
    
end

%%% Highlight cell in table %%%
col = [1 0 0];
cols = get(conpanel.table,'BackgroundColor');
cols(conpanel.selectedidx,:) = col;
set(conpanel.table,'BackgroundColor',cols)


% Save user data
set(conpanel.table,'CellSelectionCallback',@cellselect)
set(TunFig.conpanel,'userdata',conpanel)
set(TunFig.surfpanel,'userdata',surfpanel)
set(TunFig.popanel,'userdata',popanel)
set(TunFig.cellpanel,'userdata',cellpanel)
set(150,'userdata',TunFig)
toc
end

function SelectTuningDirs(a,~)

% Grab current point before anything else
whichpt = get(gca,'CurrentPoint');
whichpt = whichpt(1,[1 2]);
[theta,~] = cart2pol(whichpt(1),whichpt(2));

% Load Variables
TunFig = get(150,'UserData');
conpanel = get(TunFig.conpanel,'UserData');
popanel = get(TunFig.popanel,'UserData');
set(conpanel.table,'CellSelectionCallback',[])

% Which nd
if exist('a','var') && strcmp(a.Parent.Title.String,'1D Tuning')
    whichnd = 'oneD';
    ndL = conpanel.oneD.L;
    if popanel.uicontrols.which1Dsurf.Value == 1
        surfL = conpanel.oneD.L;
    elseif popanel.uicontrols.which1Dsurf.Value == 2
        surfL = conpanel.oneD.unichrom.L;
    elseif popanel.uicontrols.which1Dsurf.Value == 3
        surfL = conpanel.oneD.bichrom.L;
    end
elseif exist('a','var') && strcmp(a.Parent.Title.String,'2D Tuning')
    whichnd = 'twoD';
    ndL = ~conpanel.oneD.L;
    if popanel.uicontrols.which2Dsurf.Value == 1
        surfL = ~conpanel.oneD.L;
    elseif popanel.uicontrols.which2Dsurf.Value == 2
        surfL = conpanel.twoD.unichrom.hypL;
    elseif popanel.uicontrols.which2Dsurf.Value == 3
        surfL = conpanel.twoD.bichrom.hypL;
    elseif popanel.uicontrols.which2Dsurf.Value == 4
        surfL = conpanel.twoD.unichrom.eliL;
    elseif popanel.uicontrols.which2Dsurf.Value == 5
        surfL = conpanel.twoD.bichrom.eliL;
    end
end

% Which monkey
if popanel.uicontrols.whichMonk.Value == 1
    monkL = ones(size(conpanel.monk.M));
elseif popanel.uicontrols.whichMonk.Value == 2
    monkL = conpanel.monk.N;
elseif popanel.uicontrols.whichMonk.Value == 3
    monkL = conpanel.monk.M;
end

% Find the selected collection of cells
[tout,~] = rose(conpanel.tuning,popanel.nangs);
edges = unique(tout-pi);
lowerTheta = edges(find(sign(theta-edges)==1,1,'last'));
upperTheta = edges(find(sign(theta-edges)==-1,1,'first'));
angL = conpanel.tuning >= lowerTheta & conpanel.tuning <= upperTheta;

% FInd intersection of all criteria
L = angL & ndL & surfL & monkL & ~conpanel.excludeL;
if all(L == conpanel.(whichnd).selectedangs)
    conpanel.(whichnd).selectedangs = zeros(size(conpanel.tuning));
else
    conpanel.(whichnd).selectedangs = L;
end

% Save user data
set(TunFig.conpanel,'userdata',conpanel)
set(TunFig.popanel,'userdata',popanel)
set(150,'userdata',TunFig)

% Replot both pop axes
resetAll()

end

function resetAll(~,~)

exclusioncriteria()
DispTable()
ShowPopTun()
PlotNLLs()

end

function SelectNLLthresh_setval(~,~)

% Load figure and pop variables
TunFig = get(150,'userdata');
conpanel = get(TunFig.conpanel,'userdata');

% Assign Thresh and adjust graphs
conpanel.normLLthresh = str2double(get(conpanel.uicontrols.normLLthresh,'string'));

% Distinguish 1D from 2D using bichromatic fits
conpanel.oneD.L = (conpanel.twoD.bichrom.normLLs - conpanel.oneD.bichrom.normLLs) < conpanel.normLLthresh;

% Distinguish unichrom from bichrom
conpanel.oneD.unichrom.L = (conpanel.oneD.bichrom.normLLs - conpanel.oneD.unichrom.normLLs) < conpanel.normLLthresh...
    & conpanel.oneD.L;
conpanel.twoD.unichrom.L = (conpanel.twoD.bichrom.normLLs - conpanel.twoD.unichrom.normLLs) < conpanel.normLLthresh...
    & ~conpanel.oneD.L;

% Parvals
conpanel.oneD.params = nan(size(conpanel.oneD.unichrom.parvals));
conpanel.oneD.params(conpanel.oneD.unichrom.L,:) = conpanel.oneD.unichrom.parvals(conpanel.oneD.unichrom.L,:);
conpanel.oneD.params(conpanel.oneD.bichrom.L,:) = conpanel.oneD.bichrom.parvals(conpanel.oneD.bichrom.L,:);
conpanel.twoD.params = nan(size(conpanel.oneD.unichrom.parvals));
conpanel.twoD.params(conpanel.twoD.unichrom.L,:) = conpanel.twoD.unichrom.parvals(conpanel.twoD.unichrom.L,:);
conpanel.twoD.params(~conpanel.twoD.unichrom.L,:) = conpanel.twoD.bichrom.parvals(~conpanel.twoD.unichrom.L,:);

% Sort 2D surfs by type
orax = conpanel.twoD.params(:,4);
conpanel.twoD.unichrom.hypL = orax < 0 & conpanel.twoD.unichrom.L & ~conpanel.oneD.L;
conpanel.twoD.unichrom.eliL = orax > 0 & conpanel.twoD.unichrom.L & ~conpanel.oneD.L;
conpanel.twoD.bichrom.hypL = orax < 0 & ~conpanel.twoD.unichrom.L & ~conpanel.oneD.L;
conpanel.twoD.bichrom.eliL = orax > 0 & ~conpanel.twoD.unichrom.L & ~conpanel.oneD.L;

% Tuning
conpanel.tuning = nan(size(conpanel.oneD.L));
conpanel.tuning(conpanel.oneD.L) = conpanel.oneD.params(conpanel.oneD.L,end-1);
conpanel.tuning(~conpanel.oneD.L) = conpanel.twoD.params(~conpanel.oneD.L,end-1);
if isempty(conpanel.oneD.selectedangs)
    conpanel.oneD.selectedangs = zeros(size(conpanel.oneD.L));
    conpanel.twoD.selectedangs = zeros(size(conpanel.oneD.L));
end
conpanel.diffnormLL = conpanel.twoD.bichrom.normLLs - conpanel.oneD.bichrom.normLLs;

% Conf
conpanel.conf(conpanel.oneD.L) = conpanel.oneD.conf(conpanel.oneD.L); 
conpanel.conf(~conpanel.oneD.L) = conpanel.twoD.conf(~conpanel.oneD.L);

% Save user data
set(TunFig.conpanel,'userdata',conpanel)
set(150,'userdata',TunFig)

resetAll()

end

function SelectNLLthresh_clickhist(~,~)

% Grab current point before anything else
h = gca;
whichpt = get(h,'CurrentPoint');

% Load figure and pop variables
TunFig = get(150,'userdata');
conpanel = get(TunFig.conpanel,'userdata');

% Assign Thresh and adjust graphs
conpanel.normLLthresh = rndofferr(whichpt(1,1),3);
set(conpanel.uicontrols.normLLthresh,'string',conpanel.normLLthresh);

% Distinguish 1D from 2D using bichromatic fits
conpanel.oneD.L = (conpanel.twoD.bichrom.normLLs - conpanel.oneD.bichrom.normLLs) < conpanel.normLLthresh;

% Distinguish unichrom from bichrom
conpanel.oneD.unichrom.L = (conpanel.oneD.bichrom.normLLs - conpanel.oneD.unichrom.normLLs) < conpanel.normLLthresh...
    & conpanel.oneD.L;
conpanel.twoD.unichrom.L = (conpanel.twoD.bichrom.normLLs - conpanel.twoD.unichrom.normLLs) < conpanel.normLLthresh...
    & ~conpanel.oneD.L;

% Parvals
conpanel.oneD.params = nan(size(conpanel.oneD.unichrom.parvals));
conpanel.oneD.params(conpanel.oneD.unichrom.L,:) = conpanel.oneD.unichrom.parvals(conpanel.oneD.unichrom.L,:);
conpanel.oneD.params(conpanel.oneD.bichrom.L,:) = conpanel.oneD.bichrom.parvals(conpanel.oneD.bichrom.L,:);
conpanel.twoD.params = nan(size(conpanel.oneD.unichrom.parvals));
conpanel.twoD.params(conpanel.twoD.unichrom.L,:) = conpanel.twoD.unichrom.parvals(conpanel.twoD.unichrom.L,:);
conpanel.twoD.params(~conpanel.twoD.unichrom.L,:) = conpanel.twoD.bichrom.parvals(~conpanel.twoD.unichrom.L,:);

% Sort by surface type
orax = conpanel.twoD.params(:,4);
conpanel.twoD.unichrom.hypL = orax < 0 & conpanel.twoD.unichrom.L & ~conpanel.oneD.L;
conpanel.twoD.unichrom.eliL = orax > 0 & conpanel.twoD.unichrom.L & ~conpanel.oneD.L;
conpanel.twoD.bichrom.hypL = orax < 0 & ~conpanel.twoD.unichrom.L & ~conpanel.oneD.L;
conpanel.twoD.bichrom.eliL = orax > 0 & ~conpanel.twoD.unichrom.L & ~conpanel.oneD.L;

% Tuning
conpanel.tuning = nan(size(conpanel.oneD.L));
conpanel.tuning(conpanel.oneD.L) = conpanel.oneD.params(conpanel.oneD.L,end-1);
conpanel.tuning(~conpanel.oneD.L) = conpanel.twoD.params(~conpanel.oneD.L,end-1);
if isempty(conpanel.oneD.selectedangs)
    conpanel.oneD.selectedangs = zeros(size(conpanel.oneD.L));
    conpanel.twoD.selectedangs = zeros(size(conpanel.oneD.L));
end
conpanel.diffnormLL = conpanel.twoD.bichrom.normLLs - conpanel.oneD.bichrom.normLLs;

% Conf
% conpanel.conf(conpanel.oneD.L) = conpanel.oneD.conf(conpanel.oneD.L); 
% conpanel.conf(~conpanel.oneD.L) = conpanel.twoD.conf(~conpanel.oneD.L);

% Save user data
set(TunFig.conpanel,'userdata',conpanel)
set(150,'userdata',TunFig)

resetAll()


end


%%% Setup %%%

function SetUpFig()
global hlink

figure(150); clf;
set(150,'pos',[50 75 1200 750],'numbertitle','off','name','Population Tuning');

% Set up panels
TunFig.conpanel = uipanel(gcf,'units','normalized','pos',[.01 .01 .3 .98]);
TunFig.popanel = uipanel(gcf,'units','normalized','pos',[.32 .01 .22 .98]);
TunFig.surfpanel = uipanel(gcf,'units','normalized','pos',[.55 .01 .22 .98]);
TunFig.cellpanel = uipanel(gcf,'units','normalized','pos',[.78 .01 .21 .98]);

%%% Set up conpanel %%%
conpanel.uicontrols.reanalyzeall = uicontrol('parent',TunFig.conpanel,...
    'style','pushbutton','string','Reanalyze All Cells','units','normalized',...
    'pos',[.05 .01 .425 .05],'backgroundcolor',[1 0 0],'callback',@ReanalAll);
conpanel.uicontrols.reanalyze = uicontrol('parent',TunFig.conpanel,...
    'style','pushbutton','string','Reanalyze Cell','units','normalized',...
    'pos',[.525 .01 .425 .05],'backgroundcolor',[.2 1 0],'callback',@Reanal);
conpanel.uicontrols.overview = uicontrol('parent',TunFig.conpanel,...
    'style','pushbutton','string','Overview','units','normalized',...
    'pos',[.025 .07 .3 .05],'callback',@LoadOverview);
conpanel.uicontrols.analindiv = uicontrol('parent',TunFig.conpanel,...
    'style','pushbutton','string','Surface Analysis','units','normalized',...
    'pos',[.35 .07 .3 .05],'callback',@callanalysis);
conpanel.uicontrols.analGUI = uicontrol('parent',TunFig.conpanel,...
    'style','pushbutton','string','Load Control Panel','units','normalized',...
    'pos',[.675 .07 .3 .05],'callback',@LoadControlPanel);

% Norm LL Histogram
conpanel.axes.normLLhist = axes('parent',TunFig.conpanel,'units','normalized',...
    'pos',[.15 .45 .75 .1],'box','on');
conpanel.labels.normLL2D = uicontrol(TunFig.conpanel,'units','normalized',...
    'style','text','pos',[.05 .38 .25 .03],'string','Norm LL 2D','fontsize',10);
conpanel.labels.normLL1D = uicontrol(TunFig.conpanel,'units','normalized',...
    'style','text','pos',[.35 .38 .25 .03],'string','Norm LL 1D','fontsize',10);
conpanel.labels.normLLdiff = uicontrol(TunFig.conpanel,'units','normalized',...
    'style','text','pos',[.7 .38 .25 .03],'string','Diff Norm LL','fontsize',10);
conpanel.uicontrols.normLL2D = uicontrol(TunFig.conpanel,'units','normalized',...
    'style','edit','pos',[.05 .35 .25 .03]);
conpanel.uicontrols.normLL1D = uicontrol(TunFig.conpanel,'units','normalized',...
    'style','edit','pos',[.35 .35 .25 .03]);
conpanel.uicontrols.normLLdiff = uicontrol(TunFig.conpanel,'units','normalized',...
    'style','edit','pos',[.7 .35 .25 .03]);
conpanel.normLLthresh = .07;
conpanel.labels.normLLthresh = uicontrol(TunFig.conpanel,'units','normalized',...
    'style','text','pos',[.2 .3 .4 .03],'string','Norm LL Thresh','fontsize',12);
conpanel.uicontrols.normLLthresh = uicontrol(TunFig.conpanel,'units','normalized',...
    'style','edit','pos',[.6 .3 .25 .03],'string',conpanel.normLLthresh,...
    'Callback',@SelectNLLthresh_setval);
conpanel.oneDcol = [1 .5 0];
conpanel.twoDcol = [0 0 1];
conpanel.oneD.selectedangs = [];
conpanel.twoD.selectedangs = [];
conpanel.selectedidx = [];

% Exclusion criteria
conpanel.uicontrols.exclusioncrit = uipanel('parent',TunFig.conpanel,...
    'units','normalized','pos',[.05 .15 .9 .1],'title','Exclusion Criteria');
conpanel.uicontrols.modmin = uicontrol('parent',conpanel.uicontrols.exclusioncrit,...
    'style','edit','units','normalized','pos',[.1 .1 .2 .5],...
    'string','10','callback',@resetAll);
conpanel.labels.modmin = uicontrol('parent',conpanel.uicontrols.exclusioncrit,...
    'style','text','units','normalized','pos',[.1 .6 .2 .3],...
    'string','Mod min');
conpanel.uicontrols.kappamax = uicontrol('parent',conpanel.uicontrols.exclusioncrit,...
    'style','edit','units','normalized','pos',[.4 .1 .2 .5],...
    'string','4','callback',@resetAll);
conpanel.labels.kappamax = uicontrol('parent',conpanel.uicontrols.exclusioncrit,...
    'style','text','units','normalized','pos',[.4 .6 .2 .3],...
    'string','Kappa max');
conpanel.uicontrols.confmax = uicontrol('parent',conpanel.uicontrols.exclusioncrit,...
    'style','edit','units','normalized','pos',[.7 .1 .2 .5],...
    'string','360','callback',@resetAll);
conpanel.labels.confmax = uicontrol('parent',conpanel.uicontrols.exclusioncrit,...
    'style','text','units','normalized','pos',[.7 .6 .2 .3],...
    'string','1D Conf min');


%%% Set up popanel %%%
popanel.nangs = 20;
popanel.axes.oneD = axes('parent',TunFig.popanel,'units','normalized',...
    'pos',[.15 .625 .8 .35]); box on;
popanel.axes.twoD = axes('parent',TunFig.popanel,'units','normalized',...
    'pos',[.15 .2 .8 .35]); box on;
strs = {'All 1D Surfs' 'Unichromatic' 'Bichromatic'};
popanel.uicontrols.which1Dsurf = uicontrol(TunFig.popanel,...
    'style','popup','units','normalized','pos',[.2 .575 .8 .02],...
    'string',strs,'value',1,'callback',@resetAll);
strs = {'All 2D Surfs' 'Unichromatic hyperbola' 'Bichromatic hyperbola' 'Unichromatic elipse' 'Bichromatic elipse'};
popanel.uicontrols.which2Dsurf = uicontrol(TunFig.popanel,...
    'style','popup','units','normalized','pos',[.2 .15 .8 .02],...
    'string',strs,'value',1,'callback',@resetAll);
strs = {'Both Monkeys' 'Nut Data' 'Maui Data'};
popanel.uicontrols.whichMonk = uicontrol('parent',TunFig.popanel,...
    'style','popup','units','normalized','pos',[.1 .05 .8 .05],...
    'string',strs,'callback',@resetAll);
% strs = {'0-350 ms free fit' '50-250 ms fixed' '80-280 ms fixed'};
% popanel.uicontrols.whichcntwin = uicontrol('parent',TunFig.popanel,...
%     'style','popup','units','normalized','pos',[.1 .01 .8 .05],...
%     'string',strs,'callback',@resetAll);

%%% Set up surfpanel %%%
surfpanel.axes.surface1d = axes('parent',TunFig.surfpanel,'units','normalized',...
    'pos',[.2 .7 .7 .25],'box','on'); axis square;
surfpanel.axes.surface2d = axes('parent',TunFig.surfpanel,'units','normalized',...
    'pos',[.2 .2 .7 .25],'box','on'); axis square;
% linkax = [surfpanel.axes.surface1d, surfpanel.axes.surface2d];
% hlink = linkprop(linkax,...
%     {'cameraposition','cameraupvector','cameratarget','cameraviewangle'});

% 1D params
surfpanel.labels.params.oneD.A = uicontrol(TunFig.surfpanel,'units','normalized',...
    'style','text','pos',[.05 .62 .2 .03],'string','A','fontsize',12);
surfpanel.uicontrols.params.oneD.A = uicontrol(TunFig.surfpanel,'units','normalized',...
    'style','edit','pos',[.25 .62 .2 .03],'fontsize',12);
surfpanel.labels.params.oneD.bl = uicontrol(TunFig.surfpanel,'units','normalized',...
    'style','text','pos',[.5 .62 .2 .03],'string','Bl','fontsize',12);
surfpanel.uicontrols.params.oneD.bl = uicontrol(TunFig.surfpanel,'units','normalized',...
    'style','edit','pos',[.75 .62 .2 .03],'fontsize',12);
surfpanel.labels.params.oneD.sig = uicontrol(TunFig.surfpanel,'units','normalized',...
    'style','text','pos',[.05 .58 .2 .03],'string','Sigmas','fontsize',12);
surfpanel.uicontrols.params.oneD.sig1 = uicontrol(TunFig.surfpanel,'units','normalized',...
    'style','edit','pos',[.25 .58 .2 .03],'fontsize',12);
surfpanel.uicontrols.params.oneD.sig2 = uicontrol(TunFig.surfpanel,'units','normalized',...
    'style','edit','pos',[.5 .58 .2 .03],'fontsize',12);
surfpanel.labels.params.oneD.exp = uicontrol(TunFig.surfpanel,'units','normalized',...
    'style','text','pos',[.25 .51 .2 .03],'string','Exp','fontsize',12);
surfpanel.uicontrols.params.oneD.exp = uicontrol(TunFig.surfpanel,'units','normalized',...
    'style','edit','pos',[.25 .54 .2 .03],'fontsize',12);
surfpanel.labels.params.oneD.rot = uicontrol(TunFig.surfpanel,'units','normalized',...
    'style','text','pos',[.5 .51 .2 .03],'string','Rot','fontsize',12);
surfpanel.uicontrols.params.oneD.rot = uicontrol(TunFig.surfpanel,'units','normalized',...
    'style','edit','pos',[.5 .54 .2 .03],'fontsize',12);
surfpanel.labels.params.oneD.kappa = uicontrol(TunFig.surfpanel,'units','normalized',...
    'style','text','pos',[.75 .51 .2 .03],'string','Kappa','fontsize',12);
surfpanel.uicontrols.params.oneD.kappa = uicontrol(TunFig.surfpanel,'units','normalized',...
    'style','edit','pos',[.75 .54 .2 .03],'fontsize',12);

% 2D params
surfpanel.labels.params.twoD.A = uicontrol(TunFig.surfpanel,'units','normalized',...
    'style','text','pos',[.05 .12 .2 .03],'string','A','fontsize',12);
surfpanel.uicontrols.params.twoD.A = uicontrol(TunFig.surfpanel,'units','normalized',...
    'style','edit','pos',[.25 .12 .2 .03],'fontsize',12);
surfpanel.labels.params.twoD.bl = uicontrol(TunFig.surfpanel,'units','normalized',...
    'style','text','pos',[.5 .12 .2 .03],'string','Bl','fontsize',12);
surfpanel.uicontrols.params.twoD.bl = uicontrol(TunFig.surfpanel,'units','normalized',...
    'style','edit','pos',[.75 .12 .2 .03],'fontsize',12);
surfpanel.labels.params.twoD.sig = uicontrol(TunFig.surfpanel,'units','normalized',...
    'style','text','pos',[.05 .08 .2 .03],'string','Sigmas','fontsize',12);
surfpanel.uicontrols.params.twoD.sig1 = uicontrol(TunFig.surfpanel,'units','normalized',...
    'style','edit','pos',[.25 .08 .2 .03],'fontsize',12);
surfpanel.uicontrols.params.twoD.sig2 = uicontrol(TunFig.surfpanel,'units','normalized',...
    'style','edit','pos',[.5 .08 .2 .03],'fontsize',12);
surfpanel.uicontrols.params.twoD.orthsig = uicontrol(TunFig.surfpanel,'units','normalized',...
   'style','edit','pos',[.75 .08 .2 .03],'fontsize',12);
surfpanel.labels.params.twoD.exp = uicontrol(TunFig.surfpanel,'units','normalized',...
    'style','text','pos',[.25 .01 .2 .02],'string','Exp','fontsize',12);
surfpanel.uicontrols.params.twoD.exp = uicontrol(TunFig.surfpanel,'units','normalized',...
    'style','edit','pos',[.25 .04 .2 .03],'fontsize',12);
surfpanel.labels.params.twoD.rot = uicontrol(TunFig.surfpanel,'units','normalized',...
    'style','text','pos',[.5 .01 .2 .02],'string','Rot','fontsize',12);
surfpanel.uicontrols.params.twoD.rot = uicontrol(TunFig.surfpanel,'units','normalized',...
    'style','edit','pos',[.5 .04 .2 .03],'fontsize',12);
surfpanel.labels.params.twoD.kappa = uicontrol(TunFig.surfpanel,'units','normalized',...
    'style','text','pos',[.75 .01 .2 .02],'string','Kappa','fontsize',12);
surfpanel.uicontrols.params.twoD.kappa = uicontrol(TunFig.surfpanel,'units','normalized',...
    'style','edit','pos',[.75 .04 .2 .03],'fontsize',12);

%%% Set up stats panel %%%
cellpanel.psth.axes = axes('parent',TunFig.cellpanel,'units','normalized',...
    'pos',[.1 .7 .8 .25]); cla; hold on; box on;
cellpanel.pLpM.axes = axes('parent',TunFig.cellpanel,'units','normalized',...
    'pos',[.55 .5 .35 .1],'xtick',[],'ytick',[]); box on;
cellpanel.mLmM.axes = axes('parent',TunFig.cellpanel,'units','normalized',...
    'pos',[.1 .35 .35 .1],'xtick',[],'ytick',[]); box on;
cellpanel.pLmM.axes = axes('parent',TunFig.cellpanel,'units','normalized',...
    'pos',[.55 .35 .35 .1],'xtick',[],'ytick',[]); box on;
cellpanel.mLpM.axes = axes('parent',TunFig.cellpanel,'units','normalized',...
    'pos',[.1 .5 .35 .1],'xtick',[],'ytick',[]); box on;
cellpanel.linfil.axes = axes('parent',TunFig.cellpanel,'units','normalized',...
    'pos',[.1 .05 .8 .25],'xtick',[],'ytick',[]); box on;

% Save user data
set(TunFig.conpanel,'userdata',conpanel)
set(TunFig.popanel,'userdata',popanel)
set(TunFig.surfpanel,'userdata',surfpanel)
set(TunFig.cellpanel,'userdata',cellpanel)
set(150,'userdata',TunFig)

end

function UnpackPopulationData()
global GLMSPopData

% Load figure variables
TunFig = get(150,'userdata');
conpanel = get(TunFig.conpanel,'userdata');
popanel = get(TunFig.popanel,'userdata');

% Grab saved population data
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
    %library = 'C:\Users\Patty\Dropbox\Patrick\GLMS Data\';
end

load([library 'GLMSPopData.mat'])
datatypes = GLMSPopData(1,:);

% Monkey
filenames = cat(1,GLMSPopData{2:end,strcmp(datatypes,'Datafile')});
monk = filenames(:,1);
conpanel.monk.M = monk == 'M';
conpanel.monk.N = monk == 'N';

% Surf parameters
surfparams = [GLMSPopData{2:end,strcmp(datatypes,'Surface Parameters')}]';
fitcomps = [surfparams.fitcomps]';
oneD = [fitcomps.oneD]';
parvals = cat(1,oneD.parvals);
conpanel.oneD.unichrom.parvals = parvals(1:2:end,:);
conpanel.oneD.bichrom.parvals = parvals(2:2:end,:);
nLLs = [oneD.nLLs]';
conpanel.oneD.unichrom.LLs = -nLLs(:,1);
conpanel.oneD.bichrom.LLs = -nLLs(:,2);
normLLs = cat(1,oneD.normLL);
conpanel.oneD.unichrom.normLLs = normLLs(1:2:end);
conpanel.oneD.bichrom.normLLs = normLLs(2:2:end);

twoD = [fitcomps.twoD]';
parvals = cat(1,twoD.parvals);
conpanel.twoD.unichrom.parvals = parvals(1:2:end,:);
conpanel.twoD.bichrom.parvals = parvals(2:2:end,:);
nLLs = [twoD.nLLs]';
conpanel.twoD.unichrom.LLs = -nLLs(:,1);
conpanel.twoD.bichrom.LLs = -nLLs(:,2);
normLLs = cat(1,twoD.normLL);
conpanel.twoD.unichrom.normLLs = normLLs(1:2:end);
conpanel.twoD.bichrom.normLLs = normLLs(2:2:end);

% Distinguish 1D from 2D using bichromatic fits
conpanel.oneD.L = (conpanel.twoD.bichrom.normLLs - conpanel.oneD.bichrom.normLLs) < conpanel.normLLthresh;

% Distinguish unichrom from bichrom
conpanel.oneD.unichrom.L = (conpanel.oneD.bichrom.normLLs - conpanel.oneD.unichrom.normLLs) < conpanel.normLLthresh;
conpanel.oneD.bichrom.L = (conpanel.oneD.bichrom.normLLs - conpanel.oneD.unichrom.normLLs) >= conpanel.normLLthresh;
conpanel.twoD.unichrom.L = (conpanel.twoD.bichrom.normLLs - conpanel.twoD.unichrom.normLLs) < conpanel.normLLthresh;
conpanel.twoD.bichrom.L = (conpanel.twoD.bichrom.normLLs - conpanel.twoD.unichrom.normLLs) >= conpanel.normLLthresh;


% Parvals
conpanel.oneD.params = nan(size(conpanel.oneD.unichrom.parvals));
conpanel.oneD.params(conpanel.oneD.unichrom.L,:) = conpanel.oneD.unichrom.parvals(conpanel.oneD.unichrom.L,:);
conpanel.oneD.params(conpanel.oneD.bichrom.L,:) = conpanel.oneD.bichrom.parvals(conpanel.oneD.bichrom.L,:);
conpanel.twoD.params = nan(size(conpanel.oneD.unichrom.parvals));
conpanel.twoD.params(conpanel.twoD.unichrom.L,:) = conpanel.twoD.unichrom.parvals(conpanel.twoD.unichrom.L,:);
conpanel.twoD.params(conpanel.twoD.bichrom.L,:) = conpanel.twoD.bichrom.parvals(conpanel.twoD.bichrom.L,:);

% Sort by surface type
orax = conpanel.twoD.params(:,4);
conpanel.twoD.unichrom.hypL = orax < 0 & conpanel.twoD.unichrom.L & ~conpanel.oneD.L;
conpanel.twoD.unichrom.eliL = orax > 0 & conpanel.twoD.unichrom.L & ~conpanel.oneD.L;
conpanel.twoD.bichrom.hypL = orax < 0 & ~conpanel.twoD.unichrom.L & ~conpanel.oneD.L;
conpanel.twoD.bichrom.eliL = orax > 0 & ~conpanel.twoD.unichrom.L & ~conpanel.oneD.L;

% Tuning
conpanel.tuning = nan(size(conpanel.oneD.L));
conpanel.tuning(conpanel.oneD.L) = conpanel.oneD.params(conpanel.oneD.L,end-1);
conpanel.tuning(~conpanel.oneD.L) = conpanel.twoD.params(~conpanel.oneD.L,end-1);
if isempty(conpanel.oneD.selectedangs)
    conpanel.oneD.selectedangs = zeros(size(conpanel.oneD.L));
    conpanel.twoD.selectedangs = zeros(size(conpanel.oneD.L));
end
conpanel.diffnormLL = conpanel.twoD.bichrom.normLLs - conpanel.oneD.bichrom.normLLs; %using bichrom fits

% % Tuning Std
% conpanel.conf = nan(size(GLMSPopData,1)-1,1);
% confint = [GLMSPopData{2:end,strcmp(datatypes,'Confidence Intervals')}]';
% t = [confint.bfconfdeg]';
% oneD = [t.oneD]';
% twoD = [t.twoD]';
% conpanel.oneD.unichrom.conf = [oneD.uc]';
% conpanel.oneD.bichrom.conf = [oneD.uc]';
% conpanel.twoD.unichrom.conf = [twoD.uc]';
% conpanel.twoD.bichrom.conf = [twoD.uc]';
% conpanel.oneD.conf(conpanel.oneD.unichrom.L) = conpanel.oneD.unichrom.conf(conpanel.oneD.unichrom.L); 
% conpanel.oneD.conf(conpanel.oneD.bichrom.L) = conpanel.oneD.bichrom.conf(conpanel.oneD.bichrom.L); 
% conpanel.twoD.conf(conpanel.twoD.unichrom.L) = conpanel.twoD.unichrom.conf(conpanel.twoD.unichrom.L); 
% conpanel.twoD.conf(conpanel.twoD.bichrom.L) = conpanel.twoD.bichrom.conf(conpanel.twoD.bichrom.L); 
% conpanel.conf(conpanel.oneD.L) = conpanel.oneD.conf(conpanel.oneD.L);
% conpanel.conf(~conpanel.oneD.L) = conpanel.twoD.conf(~conpanel.oneD.L);

% Note: I should make one field called conpanel.params where all the
% correct parameters are loaded (2D vs 1D, and bi vs uni), like in
% GLMSPopGUI_Params.


% Save user data
set(TunFig.conpanel,'userdata',conpanel)
set(TunFig.popanel,'userdata',popanel)
set(150,'userdata',TunFig)

end

function DispTable()
global GLMSPopData

% Load figure and pop variables
TunFig = get(150,'userdata');
conpanel = get(TunFig.conpanel,'userdata');
popanel = get(TunFig.popanel,'userdata');

% Unpack data
datatypes = GLMSPopData(1,:);
data = cell(size(GLMSPopData,1)-1,3);
data(:,1) = GLMSPopData(2:end,strcmp(datatypes,'Datafile'));
data(:,2) = GLMSPopData(2:end,strcmp(datatypes,'Subunit'));
data(:,3) = GLMSPopData(2:end,strcmp(datatypes,'normLLdiff'));
%data(:,4) = num2cell(conpanel.conf);

% Repackage data and name columns
colname = {'Datafile' 'Sub' 'Norm LL Diff' 'Conf'};
colformat = cell(1,size(data,2));
colformat(:) = {'char'};

% For highlighting selected cells
cols = ones(size(data,1),3);
if sum(conpanel.oneD.selectedangs) > 0
    cols(conpanel.oneD.selectedangs,:) = repmat(conpanel.oneDcol,sum(conpanel.oneD.selectedangs),1);
end
if sum(conpanel.twoD.selectedangs) > 0
    cols(conpanel.twoD.selectedangs,:) = repmat(conpanel.twoDcol,sum(conpanel.twoD.selectedangs),1);
end

% Mark excluded data
cols(conpanel.excludeL,:) = repmat([.5 .5 .5],sum(conpanel.excludeL),1);

% Display table 
conpanel.table = uitable('parent',TunFig.conpanel,...
    'units','normalized','pos',[.01 .6  .98 .39],...
    'data',data,'columnname',colname,...
    'columnformat',colformat,...
    'BackgroundColor',cols,...
    'CellSelectionCallback',@cellselect);

% Save user data
set(TunFig.conpanel,'userdata',conpanel)
set(TunFig.popanel,'userdata',popanel)
set(150,'userdata',TunFig)

end

function PlotNLLs()

% Load figure and pop variables
TunFig = get(150,'userdata');
conpanel = get(TunFig.conpanel,'userdata');
popanel = get(TunFig.popanel,'userdata');
surfpanel = get(TunFig.surfpanel,'userdata');
cellpanel = get(TunFig.cellpanel,'userdata');

% Which monkey
if popanel.uicontrols.whichMonk.Value == 1
    monkL = true(size(conpanel.monk.M));
elseif popanel.uicontrols.whichMonk.Value == 2
    monkL = logical(conpanel.monk.N);
elseif popanel.uicontrols.whichMonk.Value == 3
    monkL = logical(conpanel.monk.M);
end

% Which 1D surface
if popanel.uicontrols.which1Dsurf.Value == 1
    surf1dL = conpanel.oneD.L;
elseif popanel.uicontrols.which1Dsurf.Value == 2
    surf1dL = conpanel.oneD.unichrom.L;
elseif popanel.uicontrols.which1Dsurf.Value == 3
    surf1dL = conpanel.oneD.bichrom.L;
end

% Which 2D Surface
if popanel.uicontrols.which2Dsurf.Value == 1
    surf2dL = ones(size(conpanel.oneD.L));
elseif popanel.uicontrols.which2Dsurf.Value == 2
    surf2dL = conpanel.twoD.unichrom.hypL;
elseif popanel.uicontrols.which2Dsurf.Value == 3
    surf2dL = conpanel.twoD.bichrom.hypL;
elseif popanel.uicontrols.which2Dsurf.Value == 4
    surf2dL = conpanel.twoD.unichrom.eliL;
elseif popanel.uicontrols.which2Dsurf.Value == 5
    surf2dL = conpanel.twoD.bichrom.eliL;
end

L = monkL & (surf1dL | surf2dL) & ~conpanel.excludeL;

% Plot histogram of deviations from unity line
conpanel.diffnormLL = conpanel.twoD.bichrom.normLLs - conpanel.oneD.bichrom.normLLs;
maxval = max(conpanel.diffnormLL);
axes(conpanel.axes.normLLhist); cla; hold on;
[conpanel.hist.counts,conpanel.hist.centers] = hist(conpanel.diffnormLL(L),[0:.005:maxval]);
conpanel.bar = bar(conpanel.hist.centers,conpanel.hist.counts,...
    'FaceColor','w','EdgeColor','k');
xlim([min(conpanel.diffnormLL(L)) max(conpanel.diffnormLL(L))])
xlabel('2D - 1D Norm LL')
ylabel('# of Subunits')
title('Difference b/w Normalized LL Values')
plot([conpanel.normLLthresh conpanel.normLLthresh],...
    [0 max(conpanel.hist.counts)],'r--')
set(gca,'ButtonDownFcn',@SelectNLLthresh_clickhist)
set(conpanel.bar,'ButtonDownFcn',@SelectNLLthresh_clickhist)

% Save user data
set(TunFig.conpanel,'userdata',conpanel)
set(TunFig.popanel,'userdata',popanel)
set(TunFig.surfpanel,'userdata',surfpanel)
set(TunFig.cellpanel,'userdata',cellpanel)
set(150,'userdata',TunFig)


end

function ShowPopTun()

% Load figure variables again
TunFig = get(150,'userdata');
conpanel = get(TunFig.conpanel,'userdata');
popanel = get(TunFig.popanel,'userdata');

% Which monkey
if popanel.uicontrols.whichMonk.Value == 1
    monkL = ones(size(conpanel.monk.M));
elseif popanel.uicontrols.whichMonk.Value == 2
    monkL = conpanel.monk.N;
elseif popanel.uicontrols.whichMonk.Value == 3
    monkL = conpanel.monk.M;
end

% % For plotting 1 subunit per cell
% datafiles = conpanel.table.Data(:,1);
% [uniquefile,a,b] = unique(datafiles); 
% uniqL = zeros(size(datafiles)); 
% uniqL(a) = 1; %for 1 sub/cell
% repeatsidx = find(~uniqL);
% for n = 1:numel(repeatsidx)
%   idx = repeatsidx(n);
%   
% end

% Which 1D surface
if popanel.uicontrols.which1Dsurf.Value == 1
    surf1dL = conpanel.oneD.L;
elseif popanel.uicontrols.which1Dsurf.Value == 2
    surf1dL = conpanel.oneD.unichrom.L;
elseif popanel.uicontrols.which1Dsurf.Value == 3
    surf1dL = conpanel.oneD.bichrom.L;
end

% Which 2D Surface
if popanel.uicontrols.which2Dsurf.Value == 1
    surf2dL = ones(size(conpanel.oneD.L));
elseif popanel.uicontrols.which2Dsurf.Value == 2
    surf2dL = conpanel.twoD.unichrom.hypL;
elseif popanel.uicontrols.which2Dsurf.Value == 3
    surf2dL = conpanel.twoD.bichrom.hypL;
elseif popanel.uicontrols.which2Dsurf.Value == 4
    surf2dL = conpanel.twoD.unichrom.eliL;
elseif popanel.uicontrols.which2Dsurf.Value == 5
    surf2dL = conpanel.twoD.bichrom.eliL;
end

% Plot preferred axes of 1D cells
angs = conpanel.tuning(conpanel.oneD.L & monkL & surf1dL & conpanel.uniqL & ~conpanel.excludeL);
axes(popanel.axes.oneD); cla;
popanel.hist.oneD = rose(angs,popanel.nangs); hold on;
t = findall(gca,'type','text');
if size(t,1) > 13
    delete(t(1:13)); % Deleting angle labels and outermost rho label
end
axis tight
set(popanel.hist.oneD,'color',conpanel.oneDcol,'ButtonDownFcn',@SelectTuningDirs)
title('1D Tuning')
xlabel('L-Cone Contrast')
ylabel('M-Cone Contrast')
t = get(gca,'xlim');
text(t(2)*.6,t(2)*.6,['n = ' num2str(numel(angs))])
angs = conpanel.tuning(conpanel.oneD.L & monkL & surf1dL & conpanel.uniqL & ~conpanel.excludeL & conpanel.oneD.selectedangs);
if sum(conpanel.oneD.selectedangs > 0)
    popanel.hist.oneDselected = rose2(angs,popanel.nangs,conpanel.oneDcol,...
        'edgecolor',conpanel.oneDcol,'ButtonDownFcn',@SelectTuningDirs);
end

% Plot axes of 2D cells
angs = conpanel.tuning(~conpanel.oneD.L & monkL & surf2dL & conpanel.uniqL & ~conpanel.excludeL);
axes(popanel.axes.twoD); cla;
popanel.hist.twoD = rose(angs,popanel.nangs); hold on; axis tight
t = findall(gca,'type','text');
if size(t,1) > 13
    delete(t(1:13));% Deleting angle labels and outermost rho label
end
set(popanel.hist.twoD,'color',conpanel.twoDcol,'ButtonDownFcn',@SelectTuningDirs)
title('2D Tuning')
xlabel('L-Cone Contrast')
ylabel('M-Cone Contrast')
t = get(gca,'xlim');
text(t(2)*.6,t(2)*.6,['n = ' num2str(numel(angs))])
angs = conpanel.tuning(~conpanel.oneD.L & monkL & surf2dL & conpanel.uniqL & ~conpanel.excludeL & conpanel.twoD.selectedangs);
if sum(conpanel.twoD.selectedangs > 0)
    popanel.hist.twoDselected = rose2(angs,popanel.nangs,conpanel.twoDcol);
    set(popanel.hist.twoDselected,'ButtonDownFcn',@SelectTuningDirs);
end


% Save user data
set(TunFig.conpanel,'userdata',conpanel)
set(TunFig.popanel,'userdata',popanel)
set(150,'userdata',TunFig)

end

function PlotPSTH()
global GLMP

% Load figure variables
TunFig = get(150,'userdata');
%conpanel = get(TunFig.conpanel,'userdata');
cellpanel = get(TunFig.cellpanel,'userdata');

% Bin and smooth baseline 
binsize = .001; % in seconds
bins = -.2:binsize:.5;
gaussSize = .01;
gaussSize = ceil(gaussSize/binsize);
gaussfilter = gausswin(gaussSize,3);
gaussfilter = gaussfilter./sum(gaussfilter); %normalizing
PSTH = histc(cat(1,GLMP.normspiketimes{GLMP.flashL}),bins);
PSTH = (PSTH./sum(GLMP.flashL)) ./ binsize;
smoothPSTH = conv(PSTH,gaussfilter,'same');
nrows = sum(GLMP.flashL);
rowcoords = linspace(max(smoothPSTH),0,nrows*2+1);
rowcents = rowcoords(2:2:end);
rowspace = rowcoords(2)-rowcoords(1);

% plot hist
axes(cellpanel.psth.axes); cla;
plot(bins,smoothPSTH,'linestyle','-','color',[.5 .5 .5])
plot([bins(1) bins(end)],[GLMP.blfrthresh GLMP.blfrthresh],'r--')
plot([GLMP.countingwin(1) GLMP.countingwin(1)],[0 max(smoothPSTH)],'r--')
plot([GLMP.countingwin(2) GLMP.countingwin(2)],[0 max(smoothPSTH)],'r--')
x = [];
y = [];
for n = 1:nrows
    tpts = GLMP.normspiketimes{n};
    if ~isempty(tpts)
        x = cat(1,x,repmat(tpts,1,2));
        y = cat(1,y,repmat([rowcents(n)-rowspace rowcents(n)+rowspace],numel(tpts),1));
    end
end
plot(x',y','k-')
xlim([bins(1) bins(end)])
ylim([0 max(smoothPSTH)])
set(gca,'xlim',[min(bins) max(bins)])
title('PSTH')
xlabel('Time from Stim Onset (ms)')
ylabel('Mean Firing Rate (sp/s)')

end

function exclusioncriteria()
global GLMSPopData

% Load figure variables
TunFig = get(150,'userdata');
conpanel = get(TunFig.conpanel,'userdata');

% parmeter bounds
modmin = str2double(get(conpanel.uicontrols.modmin,'string')); % in sp/s
kappamax = str2double(get(conpanel.uicontrols.kappamax,'string'));
%confmax = str2double(get(conpanel.uicontrols.confmax,'string'));

% Load data
datatypes = GLMSPopData(1,:);

% Modulation
modvals = nan(size(GLMSPopData,1)-1,1);
for n = 2:size(GLMSPopData,1)
    GLMP = GLMSPopData{n,strcmp(datatypes,'GLMP')};
    sub = GLMSPopData{n,strcmp(datatypes,'Subunit')};
    modvals(n-1) = max(GLMP.subunit{sub}.meannspikes)...
        ./ mean(GLMP.subunit{sub}.stimDur) - mean(GLMP.subunit{sub}.blfr);
end

% Kappa
kappavals = nan(size(GLMSPopData,1)-1,1);
kappavals(conpanel.oneD.L) = conpanel.oneD.params(conpanel.oneD.L,end);
kappavals(~conpanel.oneD.L) = conpanel.twoD.params(~conpanel.oneD.L,end);

% Tuning Std
% confvals = nan(size(GLMSPopData,1)-1,1);
% confint = [GLMSPopData{2:end,strcmp(datatypes,'Confidence Intervals')}]';
% t = [confint.bfconfdeg]';
% onedconf = [t.oneD]';
% confvals(conpanel.oneD.L) = onedconf(conpanel.oneD.L); % including just 1D confs

% Which don't meet criteria?
conpanel.critL = ~(modvals < modmin | kappavals > kappamax);
disp([num2str(sum(~conpanel.critL)) ' datafiles excluded based on modulation and variance.'])

% For selecting one subunit per datafile
filenames = GLMSPopData(2:end,strcmp(datatypes,'Datafile'));
[~,fnidx,ufnidx] = unique(filenames); 
repsL = ones(size(filenames,1),1); 
repsL(fnidx) = 0; %for 1 sub/cell
repeatsidx = find(repsL);
conpanel.uniqL = ones(size(filenames,1),1);
for n = 1:numel(repeatsidx)
    idx = find(ufnidx == ufnidx(repeatsidx(n)));
    glmp1 = GLMSPopData{idx(1)+1,strcmp(datatypes,'GLMP')};
    sub1 = GLMSPopData{idx(1)+1,strcmp(datatypes,'Subunit')};
    glmp2 = GLMSPopData{idx(2)+1,strcmp(datatypes,'GLMP')};
    sub2 = GLMSPopData{idx(2)+1,strcmp(datatypes,'Subunit')};
    maxnsp1 = max(glmp1.subunit{sub1}.meannspikes);
    maxnsp2 = max(glmp2.subunit{sub2}.meannspikes);
    if any(~conpanel.critL(idx))
        conpanel.uniqL(idx) = conpanel.critL(idx);
    elseif maxnsp1 > maxnsp2
        conpanel.uniqL(idx(2)) = 0;
    else
        conpanel.uniqL(idx(1)) = 0;
    end
end
disp([num2str(sum(~conpanel.uniqL & conpanel.critL)) ' secondary subunits excluded.'])

% Exclude RS datasets
conpanel.RSL = zeros(size(modvals)); 
conpanel.RSL(1:28) = 1; % hard coded, last RS dataset
disp([num2str(sum(conpanel.RSL & conpanel.uniqL & conpanel.critL)) ' RS datafiles excluded.'])

% All criteria together
conpanel.excludeL = ~conpanel.critL | ~conpanel.uniqL | conpanel.RSL;
disp([num2str(sum(conpanel.excludeL)) ' total datafiles excluded.'])
disp([num2str(sum(~conpanel.excludeL)) ' datafiles remaining.'])

% Save user data
set(TunFig.conpanel,'userdata',conpanel)
set(150,'userdata',TunFig)

end

function subonoff(~,~)

% Load figure variables
TunFig = get(150,'userdata');
cellpanel = get(TunFig.cellpanel,'userdata');

if strcmp(get(cellpanel.pLpM.sub,'Visible'),'on')

    set(cellpanel.pLpM.sub,'Visible','off')
    set(cellpanel.mLmM.sub,'Visible','off')
    set(cellpanel.pLmM.sub,'Visible','off')
    set(cellpanel.mLpM.sub,'Visible','off')
    set(cellpanel.linfil.sub,'Visible','off')
    
else
    
    set(cellpanel.pLpM.sub,'Visible','on')
    set(cellpanel.mLmM.sub,'Visible','on')
    set(cellpanel.pLmM.sub,'Visible','on')
    set(cellpanel.mLpM.sub,'Visible','on')
    set(cellpanel.linfil.sub,'Visible','on')

end

% Save user data
set(TunFig.cellpanel,'userdata',cellpanel)
set(150,'userdata',TunFig)

end


%%% Call to other functions %%%

function LoadOverview(~,~)
global GLMSPopData GLMP DN

% Load user data
TunFig = get(150,'userdata');
conpanel = get(TunFig.conpanel,'userdata');

% Load selected data
datanames = GLMSPopData(1,:);
idx = conpanel.selectedidx+1;
if isempty(idx)
    disp('Must select a datafile before overview can be shown.')
    return
end
DN = GLMSPopData{idx,strcmp(datanames,'DN')};
GLMP = GLMSPopData{idx,strcmp(datanames,'GLMP')};
GLMSGUI_Overview();

end

function LoadControlPanel(~,~)
global GLMSPopData GLMP DN

% Load user data
TunFig = get(150,'userdata');
conpanel = get(TunFig.conpanel,'userdata');

% Load selected data
datanames = GLMSPopData(1,:);
idx = conpanel.selectedidx+1;
if isempty(idx)
    disp('Must select a datafile before Control Panel can be loaded.')
    return
end
DN = GLMSPopData{idx,strcmp(datanames,'DN')};
GLMP = GLMSPopData{idx,strcmp(datanames,'GLMP')};
GLMS_AnalysisGUI();

end

function callanalysis(~,~)
global GLMSPopData GLMP

% Load figure and pop variables
TunFig = get(150,'userdata');
conpanel = get(TunFig.conpanel,'userdata');

% If no cell selected, return prompt. If selected, call analysis.
if isempty(conpanel.selectedidx)
    disp('Must select a cell first')
else
    datatype = GLMSPopData(1,:);
    GLMP = GLMSPopData{conpanel.selectedidx+1,strcmp(datatype,'GLMP')};
    sub = GLMSPopData{conpanel.selectedidx+1,strcmp(datatype,'Subunit')};
    GLMSGUI_Surface([],[],[],sub)
end

% Save user data
set(TunFig.conpanel,'userdata',conpanel)
set(150,'userdata',TunFig)

b.Indices = conpanel.selectedidx; 
cellselect([],b)

end

